# ==============================================================================
# COMPARACION DE DOS MODELOS ESPACIALES INLA-SPDE (2019-2021)
# ==============================================================================
#
# Objetivo metodologico
# ---------------------
# Comparar dos procedimientos sobre el mismo conjunto de datos:
#
#   M1. Modelo espacial base
#       Incluye efectos fijos de temperatura, velocidad del viento, lluvia y
#       presion barometrica, ademas de un campo espacial Matérn para el NO2.
#
#   M2. Modelo con coeficiente espacialmente variable
#       Conserva exactamente la estructura de M1 y añade un segundo campo SPDE
#       que permite que el efecto de la temperatura cambie espacialmente.
#
# La ecuacion conceptual del segundo modelo es:
#
#   y(s,t) = beta_0 + beta' x(s,t) + w_NO2(s)
#            + Temperatura(s,t) * w_Temperatura(s) + error(s,t)
#
# De este modo, beta_Temperatura representa el efecto medio de la temperatura
# y w_Temperatura(s) representa la desviacion espacial respecto a ese efecto.
# Todos los demas elementos se mantienen iguales para que la comparacion sea
# atribuible al campo adicional.
#
# Salidas principales
# -------------------
#   1. tabla_comparacion_modelos_2019_2021.csv
#   2. tabla_coeficientes_modelos_2019_2021.csv
#
# Nota: RMSE_ajuste se calcula con los mismos datos empleados para ajustar el
# modelo. Sirve como diagnostico descriptivo, no sustituye una
# validacion cruzada o un conjunto de prueba independiente.
# ==============================================================================


# ==============================================================================
# BLOQUE 1. CONFIGURACION DEL ANALISIS
# ==============================================================================

library(INLA)
library(data.table)
library(here)
library(Matrix)


source(here("R", "modeling", "spde_config.R"))



# Periodo y escala de los datos.
ANIOS <- 2019:2021
ESCALA <- "DIARIO"
RESPUESTA <- "LOG_NO2_DIARIO"

# Nombre completo del maestro conjunto 2019-2021.
NOMBRE_ARCHIVO_DATOS <-
    "dataset_maestro_inla_20190101_20211231_DIARIO2.rds"
RUTA_DATOS <- here(
    "data", "processed", "Maestro", "diario", NOMBRE_ARCHIVO_DATOS
)

# Las cuatro covariables solicitadas. Llueve es binaria:
#   0 = no llueve; 1 = precipitacion diaria igual o superior a 1 mm.
COVARIABLES <- c(
    "Temperatura",
    "Velocidad_Viento",
    "Llueve",
    "Presion_Barometrica",
    "intensidad"
)

# Covariable cuyo coeficiente se deja variar espacialmente en M2.
# Para estudiar otra covariable basta con cambiar este unico nombre.
COVARIABLE_CAMPO <- "Velocidad_Viento"

# Malla y opciones comunes de INLA.
NOMBRE_MALLA <- "media"
RUTA_MALLA <- here(
    "data", "processed", "Malla", "NO2",
    sprintf("malla_spde_madrid_%s.rds", NOMBRE_MALLA)
)
FAMILIA <- "gaussian"
NUM_THREADS <- 5L
ESTRATEGIA_INTEGRACION <- "eb"
VERBOSE_INLA <- TRUE
GUARDAR_MODELOS <- TRUE

DIR_SALIDA <- here(
    "outputs", "modelo", "Modelo_2019_2021", "comparacion_campos_latentes"
)
DIR_MODELOS <- here(
    "data", "processed", "Modelos", "Modelo_2019_2021"
)
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
if (GUARDAR_MODELOS) {
    dir.create(DIR_MODELOS, recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================
# BLOQUE 2. CARGA Y PREPARACION DE LOS DATOS
# ==============================================================================

if (!file.exists(RUTA_DATOS)) {
    stop("No se encuentra el dataset maestro: ", RUTA_DATOS)
}
if (!file.exists(RUTA_MALLA)) {
    stop("No se encuentra la malla: ", RUTA_MALLA)
}
if (!COVARIABLE_CAMPO %in% COVARIABLES) {
    stop("COVARIABLE_CAMPO debe pertenecer al vector COVARIABLES.")
}

datos <- as.data.table(readRDS(RUTA_DATOS))

columnas_necesarias <- c(
    "ESTACION", "FECHA", "X_km", "Y_km", RESPUESTA, COVARIABLES
)
columnas_ausentes <- setdiff(columnas_necesarias, names(datos))
if (length(columnas_ausentes) > 0L) {
    stop(
        "Faltan columnas necesarias en el maestro: ",
        paste(columnas_ausentes, collapse = ", ")
    )
}

datos[, FECHA := as.Date(FECHA)]
datos <- datos[as.integer(format(FECHA, "%Y")) %in% ANIOS]

# Se emplean casos completos para que M1 y M2 utilicen exactamente las mismas
# observaciones. Esto es imprescindible para comparar WAIC y DIC.
n_antes <- nrow(datos)
datos <- datos[complete.cases(datos[, ..columnas_necesarias])]
n_eliminadas <- n_antes - nrow(datos)

if (nrow(datos) == 0L) {
    stop("No quedan observaciones completas para ajustar los modelos.")
}
if (!all(datos$Llueve %in% c(0, 1))) {
    stop("La variable Llueve debe estar codificada exclusivamente como 0 y 1.")
}

datos[, y := get(RESPUESTA)]
setorder(datos, FECHA, ESTACION)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("DATOS DEL ANALISIS\n")
cat(strrep("=", 72), "\n", sep = "")
cat("Archivo:       ", NOMBRE_ARCHIVO_DATOS, "\n", sep = "")
cat(
    "Periodo:       ", as.character(min(datos$FECHA)), " a ",
    as.character(max(datos$FECHA)), "\n",
    sep = ""
)
cat("Observaciones: ", nrow(datos), "\n", sep = "")
cat("Estaciones:    ", uniqueN(datos$ESTACION), "\n", sep = "")
cat("Filas omitidas por datos incompletos: ", n_eliminadas, "\n", sep = "")


# ==============================================================================
# BLOQUE 3. MALLA Y CAMPOS LATENTES SPDE
# ==============================================================================

malla <- readRDS(RUTA_MALLA)

# Campo espacial comun que representa la variacion espacial residual del NO2.
spde_no2 <- crear_spde(malla, escala = ESCALA)
indice_no2 <- inla.spde.make.index(
    name = "campo_no2",
    n.spde = spde_no2$n.spde
)

# Campo adicional e independiente para el coeficiente espacial de la
# covariable elegida. Utiliza los mismos PC priors que el campo del NO2.
# constr = TRUE centra el campo y separa su desviacion espacial del efecto fijo
# medio de la covariable.
prior_spde <- PRIORS_SPDE[[ESCALA]]
spde_covariable <- inla.spde2.pcmatern(
    mesh = malla,
    alpha = 2,
    prior.range = prior_spde$prior.range,
    prior.sigma = prior_spde$prior.sigma,
    constr = TRUE
)
indice_covariable <- inla.spde.make.index(
    name = "campo_covariable",
    n.spde = spde_covariable$n.spde
)

coordenadas <- as.matrix(datos[, .(X_km, Y_km)])
A_no2 <- inla.spde.make.A(mesh = malla, loc = coordenadas)

# Matriz del coeficiente espacialmente variable:
# cada fila de A se multiplica por el valor observado de la covariable.
A_covariable <- Matrix::Diagonal(x = datos[[COVARIABLE_CAMPO]]) %*% A_no2

stopifnot(
    "La matriz A del NO2 tiene un numero de filas incorrecto" =
        nrow(A_no2) == nrow(datos),
    "La matriz A de la covariable tiene un numero de filas incorrecto" =
        nrow(A_covariable) == nrow(datos),
    "La matriz A y la malla tienen dimensiones incompatibles" =
        ncol(A_no2) == malla$n
)

# Efectos fijos compartidos por ambos modelos.
efectos_fijos <- data.frame(
    Intercept = rep(1, nrow(datos)),
    Llueve = datos$Llueve,
    Presion_Barometrica = datos$Presion_Barometrica,
    intendidad = datos$intensidad
)


# ==============================================================================
# BLOQUE 4. FORMULACION Y AJUSTE DE LOS MODELOS
# ==============================================================================

# M1: campo espacial del NO2 y cuatro covariables como efectos fijos.
formula_M1 <- y_response ~ -1 + Intercept +
    Temperatura + Velocidad_Viento + Llueve + Presion_Barometrica +
    f(campo_no2, model = spde_no2)

# M2: M1 mas un campo que permite que el coeficiente de la temperatura cambie
# entre localizaciones. A_covariable introduce el producto
# Temperatura(s,t) * w_Temperatura(s) en el predictor.
formula_M2 <- y_response ~ -1 + Intercept +
    Temperatura + Velocidad_Viento + Llueve + Presion_Barometrica +
    f(campo_no2, model = spde_no2) +
    f(campo_covariable, model = spde_covariable)

catalogo_modelos <- data.table(
    Modelo = c("M1", "M2"),
    Procedimiento = c(
        "Campo espacial del NO2",
        paste0(
            "Campo espacial del NO2 + coeficiente espacial de ",
            COVARIABLE_CAMPO
        )
    ),
    Formula = c(
        paste(deparse(formula_M1), collapse = " "),
        paste(deparse(formula_M2), collapse = " ")
    )
)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("MODELOS COMPARADOS\n")
cat(strrep("=", 72), "\n", sep = "")
print(catalogo_modelos[, .(Modelo, Procedimiento)], nrows = Inf)

crear_stack <- function(incluir_campo_covariable = FALSE) {
    if (incluir_campo_covariable) {
        inla.stack(
            data = list(y_response = datos$y),
            A = list(A_no2, A_covariable, 1),
            effects = list(
                indice_no2,
                indice_covariable,
                efectos_fijos
            ),
            tag = "estimacion"
        )
    } else {
        inla.stack(
            data = list(y_response = datos$y),
            A = list(A_no2, 1),
            effects = list(indice_no2, efectos_fijos),
            tag = "estimacion"
        )
    }
}

ajustar_modelo <- function(formula, stack, id_modelo) {
    cat("\nAjustando ", id_modelo, "...\n", sep = "")
    inicio <- Sys.time()

    modelo <- inla(
        formula = formula,
        data = inla.stack.data(stack),
        family = FAMILIA,
        num.threads = NUM_THREADS,
        verbose = VERBOSE_INLA,
        control.predictor = list(
            A = inla.stack.A(stack),
            compute = TRUE
        ),
        control.inla = list(
            strategy = "gaussian",
            int.strategy = ESTRATEGIA_INTEGRACION
        ),
        control.compute = list(
            dic = TRUE,
            waic = TRUE,
            cpo = TRUE,
            openmp.strategy = "huge"
        )
    )

    minutos <- as.numeric(difftime(Sys.time(), inicio, units = "mins"))
    cat(sprintf("%s completado en %.2f minutos.\n", id_modelo, minutos))

    list(modelo = modelo, stack = stack, minutos = minutos)
}

stack_M1 <- crear_stack(incluir_campo_covariable = FALSE)
stack_M2 <- crear_stack(incluir_campo_covariable = TRUE)

ajuste_M1 <- ajustar_modelo(formula_M1, stack_M1, "M1")
ajuste_M2 <- ajustar_modelo(formula_M2, stack_M2, "M2")


# ==============================================================================
# BLOQUE 5. COMPARACION DE LOS DOS PROCEDIMIENTOS
# ==============================================================================

obtener_varianza_residual <- function(modelo) {
    nombre_precision <- grep(
        "^Precision for the Gaussian observations",
        names(modelo$marginals.hyperpar),
        value = TRUE
    )
    if (length(nombre_precision) != 1L) {
        stop("No se ha podido identificar la precision residual gaussiana.")
    }

    # E(1 / precision) es la varianza residual posterior media.
    inla.emarginal(
        function(precision) 1 / precision,
        modelo$marginals.hyperpar[[nombre_precision]]
    )
}

resumir_modelo <- function(ajuste, id_modelo, procedimiento) {
    modelo <- ajuste$modelo
    indices <- inla.stack.index(ajuste$stack, tag = "estimacion")$data
    resumen_ajuste <- modelo$summary.fitted.values[indices, , drop = FALSE]

    media <- resumen_ajuste[, "mean"]
    sd_media <- resumen_ajuste[, "sd"]
    error <- media - datos$y

    # Para la cobertura de observaciones se suma la varianza residual. El sd de
    # fitted.values por si solo describe la incertidumbre de la media, no toda
    # la incertidumbre predictiva de una observacion nueva.
    varianza_residual <- obtener_varianza_residual(modelo)
    sd_predictiva <- sqrt(sd_media^2 + varianza_residual)
    dentro_95 <- datos$y >= media - 1.96 * sd_predictiva &
        datos$y <= media + 1.96 * sd_predictiva

    data.table(
        Modelo = id_modelo,
        Procedimiento = procedimiento,
        WAIC = modelo$waic$waic,
        WAIC_p_eff = modelo$waic$p.eff,
        DIC = modelo$dic$dic,
        DIC_p_eff = modelo$dic$p.eff,
        RMSE_ajuste = sqrt(mean(error^2, na.rm = TRUE)),
        COV95_predictiva_ajuste = 100 * mean(dentro_95, na.rm = TRUE)
    )
}

tabla_comparacion <- rbindlist(list(
    resumir_modelo(
        ajuste_M1,
        "M1",
        "Campo espacial del NO2"
    ),
    resumir_modelo(
        ajuste_M2,
        "M2",
        paste0(
            "Campo espacial del NO2 + coeficiente espacial de ",
            COVARIABLE_CAMPO
        )
    )
))

# Los deltas se calculan respecto a M1. Un valor negativo indica que M2 reduce
# la metrica correspondiente. Para WAIC, DIC y RMSE, menor es mejor.
referencia_M1 <- tabla_comparacion[Modelo == "M1"]
tabla_comparacion[, `:=`(
    Delta_WAIC_vs_M1 = WAIC - referencia_M1$WAIC,
    Delta_DIC_vs_M1 = DIC - referencia_M1$DIC,
    Delta_RMSE_vs_M1 = RMSE_ajuste - referencia_M1$RMSE_ajuste,
    Orden_WAIC = frank(WAIC, ties.method = "min")
)]

# Tabla conjunta de coeficientes fijos. Las covariables continuas ya proceden
# estandarizadas del maestro; el coeficiente de Llueve compara dias con lluvia
# frente a dias sin lluvia. Un efecto es significativo al 95 % cuando su
# intervalo posterior no contiene cero.
extraer_coeficientes <- function(modelo, id_modelo) {
    resumen <- as.data.table(
        modelo$summary.fixed,
        keep.rownames = "Variable"
    )
    resumen[Variable != "Intercept", .(
        Modelo = id_modelo,
        Variable,
        Coeficiente = mean,
        IC95 = sprintf(
            "[%.4f, %.4f]",
            `0.025quant`,
            `0.975quant`
        ),
        Significativa_95 = fifelse(
            `0.025quant` > 0 | `0.975quant` < 0,
            "Si",
            "No"
        )
    )]
}

tabla_coeficientes <- rbindlist(list(
    extraer_coeficientes(
        ajuste_M1$modelo,
        "M1"
    ),
    extraer_coeficientes(
        ajuste_M2$modelo,
        "M2"
    )
))


# ==============================================================================
# BLOQUE 6. EXPORTACION Y PRESENTACION DE RESULTADOS
# ==============================================================================

fwrite(
    tabla_comparacion,
    file.path(DIR_SALIDA, "tabla_comparacion_modelos_2019_2021.csv")
)
fwrite(
    tabla_coeficientes,
    file.path(DIR_SALIDA, "tabla_coeficientes_modelos_2019_2021.csv")
)
if (GUARDAR_MODELOS) {
    saveRDS(
        ajuste_M1$modelo,
        file.path(DIR_MODELOS, "modelo_M1_espacial_NO2.rds")
    )
    saveRDS(
        ajuste_M2$modelo,
        file.path(DIR_MODELOS, "modelo_M2_campo_temperatura.rds")
    )
}

cat("\n", strrep("=", 72), "\n", sep = "")
cat("TABLA DE COMPARACION\n")
cat(strrep("=", 72), "\n", sep = "")
print(tabla_comparacion, nrows = Inf)
cat("\nInterpretacion de los deltas:\n")
cat("  Delta < 0: M2 mejora la metrica respecto a M1.\n")
cat("  Delta > 0: M2 empeora la metrica respecto a M1.\n")
cat("  COV95 debe aproximarse a 95 %, pero aqui sigue siendo una medida de ajuste.\n")

cat("\n", strrep("=", 72), "\n", sep = "")
cat("TABLA DE COEFICIENTES FIJOS\n")
cat(strrep("=", 72), "\n", sep = "")
print(tabla_coeficientes, nrows = Inf)
cat(
    "\nEl campo de ", COVARIABLE_CAMPO,
    " no tiene un unico coeficiente: es una superficie espacial.\n",
    sep = ""
)

if (interactive()) {
    View(tabla_comparacion, title = "Comparacion de modelos 2019-2021")
    View(tabla_coeficientes, title = "Coeficientes de los modelos 2019-2021")
}

cat("\nResultados guardados en: ", DIR_SALIDA, "\n", sep = "")
