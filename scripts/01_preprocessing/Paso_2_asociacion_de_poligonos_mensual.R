# ==============================================================================
# STEP 2 (MONTHLY): CREATION OF THE SPATIO-TEMPORAL MASTER DATASET, 2019-2025
# ==============================================================================
# Versión mensual de Paso_2_asociacion_de_poligonos.R. Mismo pipeline —cruce de
# polígonos, interpolación del clima y join maestro— pero partiendo de los
# ficheros mensuales ya preprocesados, y recorriendo varios años.
#
# Se construye de forma NATIVA a partir del mensual, no agregando los maestros
# diarios: los tres orígenes mensuales existen para los 7 años con el mismo
# esquema, y así el clima se interpola directamente sobre el agregado mensual en
# vez de interpolarse en diario y promediarse después.
#
# La clave temporal es MES, character "YYYY-MM" en las tres fuentes.
# ==============================================================================

library(data.table)
library(sf)
library(gstat)
library(here)
source(here("R", "cleaning", "cleaning_functions.R"))
source(here("R", "interpolation", "FUNCIONES_INTERPOLACION.R"))

# Deactivate the S2 engine to avoid issues with spatial operations
sf_use_s2(FALSE)
on.exit(sf_use_s2(TRUE), add = TRUE)

# ==============================================================================
# BLOCK 0: CONFIGURATION (Only change this here!)
# ==============================================================================

ANIOS <- 2019:2025

ruta_no2_anio <- function(a) {
  here("data", "processed", "Contaminacion", "mensual",
       paste0("aire_madrid_", a, "_log_No2_mensuales1.rds"))
}
ruta_trafico_anio <- function(a) {
  here("data", "processed", "Trafico", "Mensual_Barrio", a,
       paste0("trafico_madrid_", a, "_mensual_barrio1.rds"))
}
ruta_meteo_anio <- function(a) {
  here("data", "processed", "Clima", "mensual",
       paste0("meteo_madrid_", a, "_mensual5.rds"))
}

ruta_out_clima <- here(
  "data", "processed", "Clima", "mensual",
  "clima_interpolado_mensual_2019_2025.rds"
)
ruta_out_maestro <- here(
  "data", "processed", "Maestro", "mensual",
  "dataset_maestro_inla_2019_2025_MENSUAL.rds"
)

# ------------------------------------------------------------------------------
# CONFIGURACIÓN DE LA INTERPOLACIÓN — es lo único que hay que tocar aquí.
# ------------------------------------------------------------------------------
#   metodo : "Media" -> media de TODAS las estaciones con dato ese mes; el mismo
#                       valor para las 24 estaciones (ignora k y p)
#            "kNN"   -> media simple de los k vecinos       (ignora p)
#            "IDW"   -> ponderación 1/d^p sobre k vecinos
#            "Vecino Cercano" -> el más próximo             (fuerza k = 1)
#            "Ensemble"       -> mezcla 1/RMSE de 1-NN, IDW p=1 y kNN
#   k      : número de vecinos usados.
#   p      : exponente de la distancia, solo aplica a "IDW".
#
# A escala mensual la interpolación ES la media de todas las estaciones: con la
# variable ya promediada sobre un mes entero, la variabilidad espacial dentro de
# Madrid es despreciable frente a la temporal, así que cada variable toma un
# único valor por mes, común a las 24 estaciones.
config_clima <- data.table(
  variable = c(
    "Temperatura", "Humedad_Relativa", "Precipitaciones",
    "Presion_Barometrica", "Radiacion_Solar", "Velocidad_Viento"
  ),
  metodo = rep("Media", 6),
  k      = NA_integer_,
  p      = NA_real_
)

LLAVE_TIEMPO <- "MES"

cat("\n", strrep("=", 70), "\n", sep = "")
cat("MAESTRO MENSUAL ", min(ANIOS), "-", max(ANIOS), "\n", sep = "")
cat(strrep("=", 70), "\n", sep = "")

# ==============================================================================
# BLOCK 1: GEOMETRIES (Districts and Neighborhoods)
# ==============================================================================
# Fuera del bucle: las geometrías no cambian de un año a otro.
mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(25830)
mapa_distritos$distrito <- limpiar_nombres(mapa_distritos$name)

mapa_barrios <- st_read(here("data", "raw", "geometrias", "BARRIOS.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(25830)
mapa_barrios$barrio <- limpiar_nombres(mapa_barrios$NOMBRE)

# ==============================================================================
# BLOCK 2: PER-YEAR PROCESSING
# ==============================================================================
# Por cada año: cargar, validar, cruzar polígonos, interpolar el clima y unir.
# Los resultados se apilan y la estandarización se hace al final, sobre el
# conjunto completo, para que los coeficientes sean comparables entre años.

# Comprobación de que MES viene como "YYYY-MM" del año que toca.
validar_mes <- function(dt, anio, nombre) {
  if (!LLAVE_TIEMPO %in% names(dt)) {
    stop("Falta la columna ", LLAVE_TIEMPO, " en ", nombre, " (", anio, ")")
  }
  mes <- as.character(dt[[LLAVE_TIEMPO]])
  if (!all(grepl("^[0-9]{4}-[0-9]{2}$", mes))) {
    stop("MES no tiene formato 'YYYY-MM' en ", nombre, " (", anio, ")")
  }
  if (!all(substr(mes, 1, 4) == as.character(anio))) {
    stop("Hay meses fuera del año ", anio, " en ", nombre)
  }
  invisible(NULL)
}

lista_anios <- list()
lista_clima <- list()

for (anio in ANIOS) {
  cat(sprintf("\n--- %d ---\n", anio))

  rutas <- c(ruta_no2_anio(anio), ruta_trafico_anio(anio), ruta_meteo_anio(anio))
  no_existen <- rutas[!file.exists(rutas)]
  if (length(no_existen) > 0) {
    stop("Faltan ficheros de entrada:\n  ", paste(no_existen, collapse = "\n  "))
  }

  dt_no2 <- setDT(readRDS(rutas[1]))
  dt_trafico <- setDT(readRDS(rutas[2]))
  dt_meteo <- setDT(readRDS(rutas[3]))

  validar_mes(dt_no2, anio, "NO2")
  validar_mes(dt_trafico, anio, "tráfico")
  validar_mes(dt_meteo, anio, "clima")

  # MES como character: garantiza que la clave del join sea del mismo tipo en
  # las tres fuentes aunque alguna llegue como factor.
  dt_no2[, (LLAVE_TIEMPO) := as.character(get(LLAVE_TIEMPO))]
  dt_trafico[, (LLAVE_TIEMPO) := as.character(get(LLAVE_TIEMPO))]
  dt_meteo[, (LLAVE_TIEMPO) := as.character(get(LLAVE_TIEMPO))]

  dt_trafico[, barrio := limpiar_nombres(barrio)]

  stopifnot(
    !any(duplicated(dt_no2[, c("ESTACION", LLAVE_TIEMPO), with = FALSE])),
    !any(duplicated(dt_meteo[, c("ESTACION", LLAVE_TIEMPO), with = FALSE])),
    !any(duplicated(dt_trafico[, c("barrio", LLAVE_TIEMPO), with = FALSE]))
  )

  # --- Columnas que no entran al maestro --------------------------------------
  # Indicadores de calidad del preprocesamiento (incluido num_dias del tráfico,
  # equivalente mensual de num_medidores) y columnas redundantes.
  cols_fuera <- c(
    "ESTADO", "fila_ausente", "num_dias", "num_medidores", "MAGNITUD", "ANO",
    grep("_estado$", names(dt_meteo), value = TRUE)
  )
  for (nm in c("dt_no2", "dt_trafico", "dt_meteo")) {
    d <- get(nm)
    sobran <- intersect(cols_fuera, names(d))
    if (length(sobran) > 0) d[, (sobran) := NULL]
  }

  if (!"NOM_TIPO" %in% names(dt_no2)) {
    stop("Falta NOM_TIPO en el NO2 mensual de ", anio)
  }

  # --- Coordenadas UTM 30N en kilómetros (las de la malla SPDE) ---------------
  coords_km <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])
  xy_km <- st_coordinates(st_transform(
    st_as_sf(coords_km, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830
  )) / 1000
  coords_km[, `:=`(X_km = xy_km[, 1], Y_km = xy_km[, 2])]
  dt_no2 <- merge(dt_no2, coords_km[, .(ESTACION, X_km, Y_km)],
    by = "ESTACION", all.x = TRUE
  )
  stopifnot(!any(is.na(dt_no2$X_km)), !any(is.na(dt_no2$Y_km)))

  # --- Asociación espacial de las estaciones ----------------------------------
  estaciones_sf <- st_as_sf(
    unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)]),
    coords = c("LONGITUD", "LATITUD"), crs = 4326
  ) |> st_transform(25830)

  est_distrito <- st_join(estaciones_sf, mapa_distritos[, "distrito"], join = st_intersects)
  est_barrio <- st_join(estaciones_sf, mapa_barrios[, "barrio"], join = st_intersects)
  dt_est_geo <- merge(
    as.data.table(est_distrito)[, .(ESTACION, distrito)],
    as.data.table(est_barrio)[, .(ESTACION, barrio)],
    by = "ESTACION"
  )
  if (any(is.na(dt_est_geo$barrio)) || any(is.na(dt_est_geo$distrito))) {
    stop("Hay estaciones sin barrio o distrito en ", anio)
  }
  dt_no2 <- merge(dt_no2, dt_est_geo, by = "ESTACION", all.x = TRUE)

  # El nº de barrios con tráfico crece con los años (127 en 2019 -> 130 en 2025),
  # así que se avisa si algún barrio con estación se queda sin cobertura.
  sin_trafico <- setdiff(unique(dt_no2$barrio), unique(dt_trafico$barrio))
  if (length(sin_trafico) > 0) {
    cat("  AVISO: barrios con estación y sin tráfico:",
      paste(sin_trafico, collapse = ", "), "\n"
    )
  }

  # --- Interpolación del clima a las estaciones de NO2 ------------------------
  # min_estaciones = 1: con el método "Media" basta con que haya una estación
  # con dato. El umbral de 7 existe para poder resolver un campo espacial, que
  # aquí no se estima. Entre 2022-09 y 2023-01 la red baja hasta 3 estaciones y
  # así esos meses siguen teniendo valor en vez de quedarse en NA.
  dt_clima_interp <- interpolar_clima_por_metodo(
    dt_meteo = dt_meteo,
    dt_objetivo = unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)]),
    config_variables = config_clima,
    llaves_tiempo = LLAVE_TIEMPO,
    min_estaciones = 1L
  )
  lista_clima[[as.character(anio)]] <- dt_clima_interp

  # --- Join maestro -----------------------------------------------------------
  cols_trafico <- intersect(
    c(LLAVE_TIEMPO, "barrio", "intensidad", "carga"), names(dt_trafico)
  )
  dt_anio <- merge(dt_no2, dt_trafico[, ..cols_trafico],
    by = c(LLAVE_TIEMPO, "barrio"), all.x = TRUE
  )
  dt_anio <- merge(dt_anio, dt_clima_interp,
    by = c("ESTACION", LLAVE_TIEMPO), all.x = TRUE
  )
  dt_anio[, ANIO := anio]

  faltan_maestro <- setdiff(config_clima$variable, names(dt_anio))
  if (length(faltan_maestro) > 0) {
    stop("Faltan variables climáticas en ", anio, ": ", paste(faltan_maestro, collapse = ", "))
  }

  cat(sprintf(
    "  %d filas | %d estaciones | %d meses\n",
    nrow(dt_anio), uniqueN(dt_anio$ESTACION), uniqueN(dt_anio[[LLAVE_TIEMPO]])
  ))
  lista_anios[[as.character(anio)]] <- dt_anio
}

dt_maestro <- rbindlist(lista_anios, use.names = TRUE, fill = TRUE)
saveRDS(rbindlist(lista_clima, use.names = TRUE, fill = TRUE), ruta_out_clima)

# ==============================================================================
# BLOCK 3: DISTRICT ID
# ==============================================================================
# Codificación determinista a partir del mapa: si dependiera del orden de
# aparición, un mismo distrito podría recibir IDs distintos según el año.
distritos_ref <- sort(unique(mapa_distritos$distrito))
dt_maestro[, ID_DISTRITO := as.integer(factor(distrito, levels = distritos_ref))]
stopifnot(!any(is.na(dt_maestro$ID_DISTRITO)))

# ==============================================================================
# BLOCK 4: COVARIATE STANDARDIZATION (Z-SCORE)
# ==============================================================================
# Sobre el apilado completo de los 7 años, NO año a año: si cada año usara su
# propia media y desviación, un mismo valor significaría cosas distintas en 2019
# y en 2025, y se perdería justo la variabilidad interanual que busca la escala
# mensual.
cols_std <- intersect(c("intensidad", "carga", config_clima$variable), names(dt_maestro))

for (v in cols_std) {
  dt_maestro[, (paste0(v, "_raw")) := get(v)]
  dt_maestro[, (v) := scale(get(v))[, 1]]
}

cat("\n--- Estandarización (sobre los 7 años juntos) ---\n")
for (v in cols_std) {
  cat(sprintf(
    "%-25s mean = %7.4f | sd = %6.4f\n",
    v, mean(dt_maestro[[v]], na.rm = TRUE), sd(dt_maestro[[v]], na.rm = TRUE)
  ))
}

setorderv(dt_maestro, c(LLAVE_TIEMPO, "ESTACION"))

# ==============================================================================
# BLOCK 5: SAVING AND QUALITY CONTROL
# ==============================================================================
dir.create(dirname(ruta_out_maestro), recursive = TRUE, showWarnings = FALSE)
saveRDS(dt_maestro, ruta_out_maestro)

cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf(
  "Maestro mensual: %d filas | %d estaciones | %d meses | %d años\n",
  nrow(dt_maestro), uniqueN(dt_maestro$ESTACION),
  uniqueN(dt_maestro[[LLAVE_TIEMPO]]), uniqueN(dt_maestro$ANIO)
))
cat(sprintf(
  "Periodo: %s -> %s\n",
  min(dt_maestro[[LLAVE_TIEMPO]]), max(dt_maestro[[LLAVE_TIEMPO]])
))
cat("Guardado en:", ruta_out_maestro, "\n")

# Meses con dato por variable: deja ver de un vistazo qué covariable cojea.
cat("\n--- Cobertura por variable (meses con al menos un dato) ---\n")
n_meses <- uniqueN(dt_maestro[[LLAVE_TIEMPO]])
cobertura <- data.table(
  variable = c("DATO_MENSUAL", cols_std),
  meses_con_dato = sapply(c("DATO_MENSUAL", cols_std), function(v) {
    uniqueN(dt_maestro[!is.na(get(v))][[LLAVE_TIEMPO]])
  }),
  pct_NA = sapply(c("DATO_MENSUAL", cols_std), function(v) {
    round(100 * mean(is.na(dt_maestro[[v]])), 2)
  })
)
cobertura[, meses_totales := n_meses]
print(cobertura)

cat("\n✅ Maestro mensual completado.\n")
