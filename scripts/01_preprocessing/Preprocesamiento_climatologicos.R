# ==============================================================================
# Preprocesing climatological data
# ==============================================================================

library(tidyverse)
library(data.table)
library(here)

source(here("R", "utilities", "dictionaries.R"))
source(here("R", "cleaning", "cleaning_functions.R"))

# ==============================================================================
# Years to be processed
# ==============================================================================

anios_procesar <- 2019:2025

# Sufijo de la salida. Se escribe en ficheros NUEVOS para no pisar los que ya
# existen: los originales (sin sufijo) siguen siendo los que consumen los
# scripts actuales hasta que la nueva versión esté validada.
# Poner "" para volver a escribir sobre los originales.
SUFIJO <- "1"

# Con TRUE el script se detiene si el fichero de destino ya existe, en lugar de
# sobrescribirlo sin avisar.
PROTEGER_EXISTENTES <- TRUE


# ==============================================================================
# Principal loop to process each year
# ==============================================================================

# Select the base folder and the stations file
carpeta_base_meteo <- here("data", "raw", "Datos metereologicos")
ruta_estaciones <- here(
  "data", "raw", "Datos metereologicos",
  "Estaciones_2019", "estaciones.csv"
)

anios_ok <- character(0)

# Loop in order to process each year and save the three levels of resolution
# (hourly, daily, monthly). TryCatch is used to handle errors and continue
# processing the next year if an error occurs.

for (anio in anios_procesar) {
  ruta_h <- here("data", "processed", "Clima", "horario", paste0("meteo_madrid_", anio, "_horario", SUFIJO, ".rds"))
  ruta_d <- here("data", "processed", "Clima", "diario", paste0("meteo_madrid_", anio, "_diario", SUFIJO, ".rds"))
  ruta_m <- here("data", "processed", "Clima", "mensual", paste0("meteo_madrid_", anio, "_mensual", SUFIJO, ".rds"))

  # Una version anual completa no se recalcula. Esto permite reanudar los siete
  # anos sin sobrescribir resultados ya validados (por ejemplo, 2019 v4).
  destinos <- c(ruta_h, ruta_d, ruta_m)
  if (all(file.exists(destinos))) {
    cat("\nSaltando", anio, ": los tres ficheros ya existen.\n")
    anios_ok <- c(anios_ok, as.character(anio))
    next
  }
  if (PROTEGER_EXISTENTES && any(file.exists(destinos))) {
    stop(
      "Salida parcial para ", anio, ": ",
      paste(basename(destinos[file.exists(destinos)]), collapse = ", "),
      "\nElimina o renombra esos ficheros, o usa otro SUFIJO."
    )
  }

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

  if (PROTEGER_EXISTENTES) {
    ya_existen <- Filter(file.exists, c(ruta_h, ruta_d, ruta_m))
    if (length(ya_existen) > 0) {
      stop(
        "Ya existen ficheros de destino para ", anio, ": ",
        paste(basename(ya_existen), collapse = ", "),
        "\nCambia SUFIJO o pon PROTEGER_EXISTENTES <- FALSE para sobrescribir."
      )
    }
  }

  saveRDS(resultado$horario, ruta_h)
  saveRDS(resultado$diario, ruta_d)
  saveRDS(resultado$mensual, ruta_m)

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
Check_1 <- readRDS(here("data", "processed", "Clima", "diario", "meteo_madrid_2025_diario5.rds"))
View(Check_1)
