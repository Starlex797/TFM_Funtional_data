# ==============================================================================
# PREPARA LOS DATOS PARA VALIDAR LA INTERPOLACIÓN SOLO SOBRE MEDIDAS REALES
# ==============================================================================
# Problema que resuelve: los ficheros procesados ya llevan imputación (lineal +
# vecino más cercano). Validar la interpolación sobre ellos es circular y sesga
# el resultado hacia 1-NN/IDW. Aquí regeneramos, desde el crudo:
#
#   1. HORARIO CON FLAGS  -> meteo_madrid_<anio>_horario_flag.rds
#      Igual que el horario procesado, pero con una columna lógica "imp_<var>"
#      que marca qué celdas fueron imputadas. Permite enmascararlas al validar.
#
#   2. DIARIO SOLO REAL   -> meteo_madrid_<anio>_diario_REAL.rds
#      Agregado diario calculado ÚNICAMENTE con las horas medidas (sin imputar),
#      con el mismo umbral del 30 % de NA. Un dato diario aquí es un promedio de
#      medidas reales, no de imputaciones.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("R", "utilities", "dictionaries.R"))
source(here("R", "cleaning", "cleaning_functions.R"))

ANIO <- 2025L
carpeta_base_meteo <- here("data", "raw", "Datos metereologicos")
ruta_estaciones    <- here("data", "raw", "Datos metereologicos",
                           "Estaciones_2019", "estaciones.csv")

# --- 1. Cargar crudo y limpiar (SIN imputar) ----------------------------------
cat("Cargando y limpiando datos crudos de", ANIO, "...\n")
dt_ubicaciones <- as.data.table(
  read.csv(ruta_estaciones, sep = ";", fileEncoding = "latin1"))
data_raw <- cargar_datos_metereo_anual(ANIO, carpeta_base_meteo)

dt_real <- limpiar_datos_metereo(data_raw, dt_ubicaciones)   # medidas reales
dt_real <- dt_real[year(FECHA) == ANIO]

id_cols <- c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "FECHA", "HORA")
cols_clima <- setdiff(names(dt_real), id_cols)
cols_clima <- cols_clima[sapply(dt_real[, ..cols_clima], is.numeric)]
cat("Variables climáticas detectadas:", paste(cols_clima, collapse = ", "), "\n")

# --- 2. HORARIO CON FLAGS (imputación marcada) --------------------------------
dt_flag <- imputar_na_horario(dt_real, cols_clima = cols_clima, maxgap = 3L,
                              marcar_imputados = TRUE)
dt_flag <- dt_flag[year(FECHA) == ANIO]

ruta_flag <- here("data", "processed", "Clima", "horario",
                  sprintf("meteo_madrid_%d_horario_flag.rds", ANIO))
saveRDS(dt_flag, ruta_flag)

# --- 3. DIARIO SOLO REAL (agregado de horas medidas, sin imputar) -------------
by_diario  <- c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "FECHA")
cols_suma  <- intersect("Precipitaciones", cols_clima)
cols_media <- setdiff(cols_clima, cols_suma)

dt_media <- dt_real[, lapply(.SD, function(x) {
  if (sum(is.na(x)) / .N >= 0.3) NA_real_ else mean(x, na.rm = TRUE)
}), by = by_diario, .SDcols = cols_media]

if (length(cols_suma) > 0) {
  dt_suma <- dt_real[, lapply(.SD, function(x) {
    if (sum(is.na(x)) / .N >= 0.3) NA_real_ else sum(x, na.rm = TRUE)
  }), by = by_diario, .SDcols = cols_suma]
  dt_diario_real <- merge(dt_media, dt_suma, by = by_diario)
} else {
  dt_diario_real <- dt_media
}
dt_diario_real[, ANO := ANIO]
setorder(dt_diario_real, FECHA)

ruta_diario <- here("data", "processed", "Clima", "diario",
                    sprintf("meteo_madrid_%d_diario_REAL.rds", ANIO))
saveRDS(dt_diario_real, ruta_diario)

# --- 4. Estadísticas: cuánto se imputó (la magnitud del problema #1) ----------
cat("\n", strrep("=", 66), "\n", sep = "")
cat("PROPORCIÓN DE VALORES IMPUTADOS (horario) por variable\n")
cat(strrep("=", 66), "\n", sep = "")
resumen <- rbindlist(lapply(cols_clima, function(v) {
  flag <- paste0("imp_", v)
  if (!flag %in% names(dt_flag)) return(NULL)
  n_real <- sum(!is.na(dt_real[[v]]))
  n_imp  <- sum(dt_flag[[flag]], na.rm = TRUE)
  data.table(Variable = v, N_real = n_real, N_imputado = n_imp,
             Pct_imputado = round(100 * n_imp / (n_real + n_imp), 1))
}))
print(resumen)

cat("\nGuardado:\n")
cat("  Horario con flags :", ruta_flag, "\n")
cat("  Diario solo real  :", ruta_diario, "\n")
cat("\nEn la validación: enmascarar las celdas con imp_<var> == TRUE (horario) y\n")
cat("usar el diario_REAL en lugar del diario procesado.\n")
