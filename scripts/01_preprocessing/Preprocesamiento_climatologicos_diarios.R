# ==============================================================================
# Preprocessing of the official daily meteorological files
# ==============================================================================
#
# Independent pipeline for data/raw/Datos metereologicos/diarios.
# The daily values published by the provider are retained as supplied. This
# script does NOT read hourly files and does NOT create monthly aggregates.

library(data.table)
library(here)

source(here("R", "utilities", "dictionaries.R"))
source(here("R", "cleaning", "cleaning_climate_daily_functions.R"))


# ==============================================================================
# Configuration
# ==============================================================================

ANIOS_PROCESAR <- 2019:2025

# This is expressed in DAYS. It is not the three-hour rule from the hourly
# pipeline. Only invalid observations (FALLO) may be interpolated; days absent
# from the source (AUSENTE) always remain NA.
IMPUTAR_FALLOS <- TRUE
MAXGAP_DIAS <- 3L

# Never overwrite validated results silently.
PROTEGER_EXISTENTES <- TRUE

CARPETA_DIARIOS <- here(
  "data", "raw", "Datos metereologicos", "diarios"
)
RUTA_ESTACIONES <- here(
  "data", "raw", "Datos metereologicos",
  "Estaciones_2019", "estaciones.csv"
)
CARPETA_SALIDA <- here("data", "processed", "Clima", "diario")


# ==============================================================================
# Annual processing
# ==============================================================================

if (!dir.exists(CARPETA_DIARIOS)) {
  stop("No existe la carpeta de datos diarios: ", CARPETA_DIARIOS)
}
if (!file.exists(RUTA_ESTACIONES)) {
  stop("No existe el catalogo de estaciones: ", RUTA_ESTACIONES)
}
dir.create(CARPETA_SALIDA, recursive = TRUE, showWarnings = FALSE)

anios_ok <- integer(0)
anios_error <- integer(0)

for (anio in ANIOS_PROCESAR) {
  cat("\n", strrep("=", 68), "\n", sep = "")
  cat("Procesando datos meteorologicos diarios de ", anio, "\n", sep = "")
  cat(strrep("=", 68), "\n", sep = "")

  ruta_salida <- here(
    "data", "processed", "Clima", "diario", "nuevo",
    paste0("meteo_madrid_", anio, "_diario_fuente.rds")
  )

  if (file.exists(ruta_salida)) {
    if (PROTEGER_EXISTENTES) {
      cat("   Saltando: ya existe ", basename(ruta_salida), "\n", sep = "")
      anios_ok <- c(anios_ok, anio)
      next
    }
    warning("Se sobrescribira ", basename(ruta_salida), ".")
  }

  resultado <- tryCatch(
    procesar_anio_meteo_diario_fuente(
      anio = anio,
      carpeta_diarios = CARPETA_DIARIOS,
      ruta_estaciones = RUTA_ESTACIONES,
      maxgap_dias = MAXGAP_DIAS,
      imputar_fallos = IMPUTAR_FALLOS
    ),
    error = function(e) {
      warning("Error procesando ", anio, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(resultado)) {
    anios_error <- c(anios_error, anio)
    next
  }

  # Validation has already completed inside procesar_anio_meteo_diario_fuente().
  # saveRDS is only reached for a coherent, complete station x day grid.
  saveRDS(resultado, ruta_salida)
  cat("   Guardado: ", basename(ruta_salida), "\n", sep = "")
  anios_ok <- c(anios_ok, anio)
}


# ==============================================================================
# Final summary
# ==============================================================================

cat("\n", strrep("=", 68), "\n", sep = "")
cat("Procesamiento diario finalizado.\n")
cat(
  "   Anios correctos: ",
  if (length(anios_ok)) paste(anios_ok, collapse = ", ") else "ninguno",
  "\n",
  sep = ""
)
if (length(anios_error)) {
  cat("   Anios con error: ", paste(anios_error, collapse = ", "), "\n", sep = "")
}
cat("   No se ha generado ninguna escala mensual.\n")
cat(strrep("=", 68), "\n", sep = "")

if (length(anios_error)) {
  stop("El procesamiento termino con errores en algun anio.")
}

Check2 <- readRDS(here("data", "processed", "Clima", "diario", "meteo_madrid_2025_diario_fuente.rds"))
View(Check2)
