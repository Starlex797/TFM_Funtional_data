# ==============================================================================
# CARGA, LIMPIEZA Y AGREGACIÓN DIARIA DE DATOS METEOROLÓGICOS
# ==============================================================================

library(tidyverse)
library(data.table)
library(here)

# 1. Cargar diccionarios y funciones
source(here("R", "dictionaries.R"))
source(here("R", "cleaning_functions.R"))

# 2. Rutas y lectura del archivo bruto
ruta_datos_metereo <- here("data", "raw", "Datos metereologicos", "2025", "2025_datos_metereo.csv")
data_raw_metereo <- fread(ruta_datos_metereo, sep = ";")

# 3. Cargar estaciones meteorológicas
# Nota: fread falla en Windows con espacios en la ruta, usamos read.csv
ruta_est_metereo <- here("data", "raw", "Datos metereologicos", "Estaciones_2019", "estaciones_mete.csv")
dt_ubicaciones_metereo <- as.data.table(read.csv(ruta_est_metereo, sep = ";", fileEncoding = "latin1"))

# 4. Limpieza y reestructuración (Llamada a tu función modular)
cat("⏳ Limpiando datos meteorológicos...\n")
datos_metereo_horarios <- limpiar_datos_metereo(data_raw_metereo, dt_ubicaciones_metereo)
head(datos_metereo_horarios)
# 4. Agregación diaria múltiple (Controlando NAs)
cat("⏳ Agregando a escala diaria (Umbral 20% NA)...\n")

# Extraemos qué columnas pertenecen al clima
cols_clima <- setdiff(names(datos_metereo_horarios), c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "FECHA", "HORA"))

datos_metereo_diarios <- datos_metereo_horarios[, lapply(.SD, function(x) {
  if (sum(is.na(x)) / .N >= 0.2) {
    NA_real_
  } else {
    mean(x, na.rm = TRUE)
  }
}), by = .(ESTACION, LONGITUD, LATITUD, X_km, Y_km, FECHA), .SDcols = cols_clima]

# 5. Creación del Índice Temporal (ID_TIEMPO)
setorder(datos_metereo_diarios, FECHA)
datos_metereo_diarios[, ID_TIEMPO := .GRP, by = FECHA]
view(datos_metereo_diarios)
# 6. Verificación y Guardado
cat("\n✅ Procesamiento completado. Resumen diario:\n")
cat("Filas totales (Estaciones × Días):", nrow(datos_metereo_diarios), "\n\n")
print(head(datos_metereo_diarios[, 1:6], 10))

saveRDS(datos_metereo_diarios, here("data", "processed", "meteo_madrid_2025_diario.rds"))
cat("\n💾 Guardado en data/processed/meteo_madrid_2025_diario.rds\n")
