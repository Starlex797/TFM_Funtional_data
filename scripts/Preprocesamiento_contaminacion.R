# scripts/01_load_clean_data.R
# Objetivo: Cargar los datos crudos de contaminación, limpiarlos y guardarlos por año.

library(here)
library(data.table)
library(lubridate)

# 2. Cargar herramientas --------------------------------------------------
source(here("R", "dictionaries.R"))
source(here("R", "cleaning_functions.R"))

# 3. Cargar ubicaciones espaciales ----------------------------------------
dt_ubicaciones_aire <- fread(here("data", "raw", "Datos_contaminacion", "Estaciones", "datos.csv"))
setDT(dt_ubicaciones_aire)

# 4. Bucle de procesamiento y limpieza ------------------------------------
anios_analisis <- c(2019, 2023, 2025)
lista_historico_aire <- list()

for (anio in anios_analisis) {
  
  # Construimos la ruta de forma segura y limpia delegando en 'here'
  ruta_archivo <- here("data", "raw", "Datos_contaminacion", "Datos", paste0("calidad_", anio, ".csv"))
  
  if(file.exists(ruta_archivo)){
    cat("✅ Procesando contaminación del año:", anio, "\n")
    
    # Leemos y limpiamos
    dt_temporal <- fread(ruta_archivo, sep = ";")
    dt_limpio   <- limpiar_aire_madrid(dt_temporal, dt_ubicaciones_aire)
    
    # Guardamos en la lista indexando por el año
    lista_historico_aire[[as.character(anio)]] <- dt_limpio
  } else {
    cat("❌ Archivo no encontrado para el año:", anio, "\n")
    cat("   Ruta buscada:", ruta_archivo, "\n\n")
  }
}

# 5. Guardado seguro en disco ---------------------------------------------
# Solo entramos a guardar si logramos limpiar algún año
if (length(lista_historico_aire) > 0) {
  for (anio in names(lista_historico_aire)) {
    dt_limpio <- lista_historico_aire[[anio]]
    
    # Guardamos cada año limpio de forma independiente en 'data/processed'
    saveRDS(dt_limpio, here("data", "processed", paste0("aire_madrid_", anio, "_limpio.rds")))
    cat("💾 Guardado con éxito el año:", anio, "en data/processed/\n")
  }
} else {
  cat("⚠️ No se ha guardado nada porque la lista de datos está vacía. Revisa las rutas de arriba.\n")
}


# ==============================================================================
# FILTRO DE DATOS A NO2, AGREGACIÓN DIARIA Y LOG-TRANSFORMACIÓN
# ==============================================================================

# 1. Cargar el año de análisis limpio (Formato Largo)
aire_madrid_2025_limpio <- readRDS(here("data", "processed", "aire_madrid_2025_limpio.rds"))
setDT(aire_madrid_2025_limpio)

# 2. Filtrar estrictamente por la magnitud de interés (NO2)
dt_no2 <- aire_madrid_2025_limpio[MAGNITUD == "NO2"]
head(dt_no2)
cat("📊 Registros horarios brutos de NO2 para 2025: ", nrow(dt_no2), "\n")

# [CORREGIDO] Guardado del dataset horario filtrado de NO2 antes de promediar
saveRDS(dt_no2, here("data", "processed", "aire_madrid_2025_No2_horarios.rds"))


# 3. PASAR DE DATOS HORARIOS A DIARIOS CON CONTROL DE NAs
# Calculamos la media aritmética de la concentración real (valor_NO2)
# Si el día tiene más del 20% de NAs (umbral_na = 0.2), el día se descarta como NA
dt_no2_diario <- agregar_a_diario(dt_no2, 
                                  col_grupo = "ESTACION", 
                                  col_fecha = "FECHA", 
                                  col_valor = "DATO", 
                                  col_longitud ="LONGITUD", 
                                  col_latitud = "LATITUD",
                                  umbral_na = 0.2)

View(dt_no2_diario)


# 4. APLICAR TRANSFORMACIÓN LOGARÍTMICA AL PROMEDIO DIARIO
# Ahora aplicamos log(x + 1) sobre el valor diario para aproximar la normalidad en INLA
# R de forma nativa mantendrá los NA intactos si DATO_DIARIO es NA.


dt_no2_diario[!is.na(DATO_DIARIO), LOG_NO2_DIARIO := log(DATO_DIARIO + 1)]


# 5. CREACIÓN DEL ÍNDICE TEMPORAL SECUENCIAL (ID_TIEMPO)
# Con el objetivo que INLA pueda procesar el tiempo en el caso que se quiera poner un 
#Autoregresivo. En caso de que no se ponga un autoregresivo, solo estaríamos ante una predicción 
#espacia. Con el autoregresivo es espacio-temporal y el ID_TIEMPO es necesario para que INLA sepa el orden temporal de los datos.
# Ordenamos cronológicamente por fecha (indispensable para el AR1)
setorder(dt_no2_diario, FECHA)

# Asignamos el índice entero continuo (1, 2, 3... T) tal como hace Wong en prep.R
dt_no2_diario[, ID_TIEMPO := .GRP, by = FECHA]


# 6. VERIFICACIONES DE CONTROL DE CALIDAD Y NAs
cat("📊 Resumen del Dataset Diario generado para INLA:\n")
cat("Total filas (Estaciones × Días): ", nrow(dt_no2_diario), "\n")
cat("Días válidos (Sin NA):            ", sum(!is.na(dt_no2_diario$LOG_NO2_DIARIO)), "\n")
cat("Días imputables por INLA (NA):    ", sum(is.na(dt_no2_diario$LOG_NO2_DIARIO)), "\n\n")

# Visualizar la estructura en la consola para confirmar el ID_TIEMPO
print(head(dt_no2_diario[, .(ESTACION, FECHA, DATO_DIARIO, LOG_NO2_DIARIO, ID_TIEMPO,LONGITUD,LATITUD)], 15))


# 7. GUARDADO FINAL SEGURO DEL DATASET DIARIO MODELIZABLE
# [CORREGIDO] Guardamos el objeto 'dt_no2_diario' con nombre explícito del año
saveRDS(dt_no2_diario, here("data", "processed", "aire_madrid_2025_No2_trans_diarios.rds"))
cat("💾 Dataset diario de NO2 guardado con éxito en data/processed/\n")


