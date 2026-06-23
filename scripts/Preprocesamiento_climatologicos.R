# ==============================================================================
# Preprocesing climatological data 
# ==============================================================================

library(tidyverse)
library(data.table)
library(here)

source(here("R", "dictionaries.R"))
source(here("R", "cleaning_functions.R"))

# ==============================================================================
# Years to be processed
# ==============================================================================

anios_procesar <- c(2019)



# ==============================================================================
# Principal loop to process each year
# ==============================================================================

# Select the base folder and the stations file
carpeta_base_meteo <- here("data", "raw", "Datos metereologicos")
ruta_estaciones    <- here("data", "raw", "Datos metereologicos",
                           "Estaciones_2019", "estaciones.csv")

anios_ok <- character(0)

# Loop in order to process each year and save the three levels of resolution
# (hourly, daily, monthly). TryCatch is used to handle errors and continue
# processing the next year if an error occurs.

for (anio in anios_procesar) {
  
  resultado <- tryCatch(
    procesar_anio_meteo(anio, carpeta_base_meteo, ruta_estaciones),
    error = function(e) {
      warning("❌ Error procesando el año ", anio, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(resultado)) next
  
  # Save the processed data to RDS files
  cat("\n💾 Guardando archivos del año", anio, "...\n")
  
  ruta_h <- here("data", "processed","Clima","horario", paste0("meteo_madrid_", anio, "_horario.rds"))
  ruta_d <- here("data", "processed","Clima","diario" , paste0("meteo_madrid_", anio, "_diario.rds"))
  ruta_m <- here("data", "processed", "Clima","mensual", paste0("meteo_madrid_", anio, "_mensual.rds"))
  
  saveRDS(resultado$horario,  ruta_h)
  saveRDS(resultado$diario,   ruta_d)
  saveRDS(resultado$mensual,  ruta_m)
  
  cat("   Horario  →", basename(ruta_h), "\n")
  cat("   Diario   →", basename(ruta_d), "\n")
  cat("   Mensual  →", basename(ruta_m), "\n")
  
  anios_ok <- c(anios_ok, as.character(anio))
}

# ==============================================================================
# Summary of the processing results
# ==============================================================================
cat("\n", strrep("=", 60), "\n")
if (length(anios_ok) == 0) {
  stop("No se pudo procesar ningún año correctamente.")
} else {
  cat("✅ Procesamiento completado.\n")
  cat("   Años guardados:", paste(anios_ok, collapse = ", "), "\n")
  cat("   Archivos generados por año: horario / diario / mensual\n")
  cat(strrep("=", 60), "\n")
}


horario_2025<-readRDS(here("data", "processed","Clima","diario","meteo_madrid_2019_diario.rds"))
view(horario_2025)
#===============================================================================

