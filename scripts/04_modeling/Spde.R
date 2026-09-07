# ==============================================================================
# Step 3: SPDE Model Setup
# ==============================================================================

library(INLA)
library(data.table)
library(here)
library(sf)

# 0. Configuración: año y escala se declaran UNA vez y se reutilizan tanto para
# leer el maestro como para construir las rutas de salida. Repetir los literales
# en cada saveRDS es lo que provoca erratas del tipo "EspaciA" / "EspaciAl".
ANIO   <- 2020
ESCALA <- "DIARIO" # "DIARIO" | "HORARIO"
MALLA  <- "media" # "gruesa" | "media" | "fina"

# 1. Cargar el dataset maestro y la malla (si estás en una sesión nueva)
dt_maestro <- readRDS(here(
  "data", "processed", "Maestro", as.character(ANIO),
  paste0("dataset_maestro_inla_", ANIO, "_", ESCALA, ".rds")
))
malla_madrid <- readRDS(here(
  "data", "processed", paste0("malla_spde_madrid_", MALLA, ".rds")
))
setDT(dt_maestro)

# 2. Definir las dimensiones temporales
# ID_TIEMPO no viene en el dataset maestro: se construye aquí, que es donde se
# monta la dimensión temporal del campo. Numera los instantes de forma
# consecutiva —días en escala diaria, pares fecha-hora en horaria— que es lo que
# el AR1 necesita para saber qué grupo sigue a cuál.
llaves_tiempo <- if ("HORA" %in% names(dt_maestro)) c("FECHA", "HORA") else "FECHA"

# Comprobación de que el fichero cargado es de la escala que se ha declarado.
escala_detectada <- if ("HORA" %in% names(dt_maestro)) "HORARIO" else "DIARIO"
if (escala_detectada != ESCALA) {
  stop("ESCALA declarada '", ESCALA, "' pero el maestro es '", escala_detectada, "'.")
}

instantes <- unique(dt_maestro[, llaves_tiempo, with = FALSE])
setorderv(instantes, llaves_tiempo)
instantes[, ID_TIEMPO := .I]

dt_maestro <- merge(dt_maestro, instantes, by = llaves_tiempo, all.x = TRUE)
stopifnot(!any(is.na(dt_maestro$ID_TIEMPO)))

# ¿Cuántos instantes únicos hay en tu dataset? (Esto formará nuestro modelo AR1)
ndays <- nrow(instantes)

# 3. Coordenadas en kilómetros (UTM 30N)
# El maestro ya las trae desde el Paso 2, en las mismas unidades que la malla
# (EPSG:25830 dividido por 1000). No se reproyectan aquí: mantener dos cálculos
# del mismo dato es la vía por la que los datos y la malla acaban descuadrados.
if (!all(c("X_km", "Y_km") %in% names(dt_maestro))) {
  stop("El maestro no trae X_km/Y_km: revisa el Bloque 1 del Paso 2.")
}

# MUY IMPORTANTE: se ordena por tiempo y luego por estación, y a partir de aquí
# NO se vuelve a reordenar ni a hacer ningún merge. La matriz A se construye
# sobre este orden exacto; si las filas se movieran después, A quedaría
# desalineada con los datos sin que R lanzara ningún error.
setorder(dt_maestro, ID_TIEMPO, ESTACION)

coords_puntos <- as.matrix(dt_maestro[, .(X_km, Y_km)])

# ==============================================================================
# 3b. ESCALA EMPÍRICA DE LA VARIABILIDAD ESPACIAL -> sigma_0 del PC prior
# ==============================================================================
# Idea: ajustar un modelo con un efecto fijo por instante. Ese efecto absorbe
# TODO lo temporal —ciclo anual, episodios, meteorología común a la ciudad— y lo
# que queda en los residuos es la variabilidad ENTRE estaciones dentro de un
# mismo instante, que es justo lo que el campo espacial tiene que explicar. Su
# desviación típica es, por tanto, una cota superior razonable para sigma.
#
#   Equivale a:  m0 <- lm(y ~ factor(FECHA)); sigma_0 <- summary(m0)$sigma
#   pero centrando por grupo, que da el MISMO número sin construir una matriz de
#   diseño con un dummy por instante (inviable en escala horaria: 8760 niveles).
#
# OJO: la respuesta YA está transformada — LOG_NO2_* es log1p(DATO_*), no log().
# No hay que volver a aplicar ningún logaritmo.
col_resp <- grep("^LOG_NO2_", names(dt_maestro), value = TRUE)[1]
if (is.na(col_resp)) stop("No encuentro la columna de respuesta LOG_NO2_*")

dt_sig <- dt_maestro[!is.na(get(col_resp))]
dt_sig[, .resid_esp := get(col_resp) - mean(get(col_resp)), by = ID_TIEMPO]
gl <- nrow(dt_sig) - uniqueN(dt_sig$ID_TIEMPO) # grados de libertad
sigma_0 <- sqrt(sum(dt_sig$.resid_esp^2) / gl)

cat("\n--- Escala de la variabilidad espacial (", col_resp, ") ---\n", sep = "")
cat(sprintf("  sd total                        : %.4f\n", sd(dt_sig[[col_resp]])))
cat(sprintf("  sd entre estaciones (medias)    : %.4f\n",
  sd(dt_sig[, mean(get(col_resp)), by = ESTACION]$V1)))
cat(sprintf("  sigma_0 = sd residual intra-instante: %.4f  <- referencia del prior\n", sigma_0))

# Qué implica cada elección de prior.sigma = c(sigma_0, p), con
# P(sigma > sigma_0) = p  =>  lambda = -log(p)/sigma_0, mediana = log(2)/lambda.
cat("  prior.sigma candidatos:\n")
for (p in c(0.50, 0.10, 0.05, 0.01)) {
  lambda_s <- -log(p) / sigma_0
  cat(sprintf(
    "    c(%.2f, %.2f) -> mediana a priori de sigma = %.3f\n",
    sigma_0, p, log(2) / lambda_s
  ))
}
dt_sig[, .resid_esp := NULL]

# 4. Crear el objeto SPDE (El modelo matemático espacial de Matérn)
# alpha = 2 es el estándar en R-INLA para superficies en 2D (equivale a nu = 1).
#
# pcmatern, NO matern: los argumentos prior.range/prior.sigma solo existen en
# inla.spde2.pcmatern. inla.spde2.matern los interpreta como su propio
# prior.range.nominal (un escalar) y falla en qr.coef.
#
#   prior.range = c(r0, p) -> P(rango < r0) = p
#   prior.sigma = c(s0, p) -> P(sigma > s0) = p   (referencia: sigma_0 de 3b)
spde <- inla.spde2.pcmatern(
  mesh = malla_madrid,
  alpha = 2,
  prior.range = c(9.3, 0.5), # mediana a priori del rango = 9.3 km = D/3
  prior.sigma = c(0.6, 0.01), # mediana a priori de sigma = 0.151
  constr = FALSE
)

# 5. Crear el Índice Espacio-Temporal
# Esto le dice a INLA: "Crea una variable latente espacial y multiplícala por los 'ndays'"

indice_s <- inla.spde.make.index(name = "campo_espacial", 
                                 n.spde = spde$n.spde, 
                                 n.group = ndays)

# 6. Crear la Matriz de Proyección (Matriz A) para los Puntos
# - loc: las coordenadas de cada fila
# - group: a qué día (1 a 365) pertenece esa coordenada
A_espacial <- inla.spde.make.A(mesh = malla_madrid, # Mi malla 
                               loc = coords_puntos, # Dónde están las estaciones 
                               group = dt_maestro$ID_TIEMPO, # En qué dia ocurrió la medicion 
                               n.group = ndays) # Total de dias 

# ==============================================================================
# MODELO SOLO ESPACIAL: Índice y Matriz A sin dimensión temporal
# Sin 'n.group' → un único campo latente espacial compartido por todos los días.
# Sin 'group'   → cada observación se proyecta solo en el espacio (sin AR1).

indice_s_solo <- inla.spde.make.index(
  name   = "campo_espacial_s",
  n.spde = spde$n.spde
  # n.group = 1 es el defecto → sin replicación temporal
)

A_espacial_s <- inla.spde.make.A(
  mesh = malla_madrid,
  loc  = coords_puntos   # Sin 'group': proyección puramente espacial
)

# 7. GUARDAR LOS OBJETOS PARA EL PASO 4
# Estructura:  data/processed/SPDE/<ANIO>/<escala>/
#                ├─ spde_madrid.rds            (común a los dos modelos)
#                ├─ Espacio_temporal/          (índice y A con grupo temporal)
#                └─ Espacial/                  (índice y A sin dimensión temporal)
# Las carpetas se crean si no existen: saveRDS no las crea y fallaría con
# "cannot open the connection".
dir_base <- here("data", "processed", "SPDE", as.character(ANIO), tolower(ESCALA))
dir_st   <- file.path(dir_base, "Espacio_temporal")
dir_s    <- file.path(dir_base, "Espacial")
for (d in c(dir_base, dir_st, dir_s)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# El objeto SPDE es el mismo para ambos modelos, así que se guarda una sola vez
# en la raíz en lugar de duplicarlo en las dos subcarpetas.
saveRDS(spde, file.path(dir_base, "spde_madrid.rds"))

# Objetos del modelo espacio-temporal
saveRDS(indice_s,   file.path(dir_st, "indice_s_madrid.rds"))
saveRDS(A_espacial, file.path(dir_st, "A_espacial_madrid.rds"))

# Objetos del modelo solo espacial (baseline)
saveRDS(indice_s_solo, file.path(dir_s, "indice_s_solo_madrid.rds"))
saveRDS(A_espacial_s,  file.path(dir_s, "A_espacial_s_madrid.rds"))

cat("Objetos guardados en:", dir_base, "\n")

# ------------------------------------------------------------------------------
# VALIDACIÓN (el script se detiene si algo no cuadra)
# ------------------------------------------------------------------------------
stopifnot("Filas de A_espacial != filas de dt_maestro" =
            nrow(A_espacial) == nrow(dt_maestro))
stopifnot("Columnas de A_espacial != nodos * días" =
            ncol(A_espacial) == malla_madrid$n * ndays)
stopifnot("Filas de A_espacial_s != filas de dt_maestro" =
            nrow(A_espacial_s) == nrow(dt_maestro))
stopifnot("Columnas de A_espacial_s != nodos de la malla" =
            ncol(A_espacial_s) == malla_madrid$n)

cat("Paso 3 completado.\n")
cat("- Nodos de la malla:", malla_madrid$n, "\n")
cat("- Días (grupos):", ndays, "\n")
cat("- A_espacial:", nrow(A_espacial), "x", ncol(A_espacial), "\n")
cat("- A_espacial_s:", nrow(A_espacial_s), "x", ncol(A_espacial_s), "\n")

