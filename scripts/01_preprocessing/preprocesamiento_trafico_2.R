library(data.table)
library(sf)
library(here)
here::i_am("scripts/01_preprocessing/preprocesamiento_trafico_2.R")
source(here("R", "cleaning", "cleaning_functions.R"))

# --- Global Parameetrs ---

# ===============================================================================
# Global Parameter in order to process the traffic data of Madrid for the year
# Traffic data is only available for the whole year, but I had to download it month by month
# ===============================================================================
anios_procesar <- 2021:2025

meses <- c(
  "Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio",
  "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre"
)

# ==============================================================================
# 0. Load and charge the official maps of districts and neighborhoods of Madrid (shapefiles)
# ==============================================================================
cat("Preparando mapas oficiales de distritos y barrios de Madrid...\n")
carpeta_geo <- here("data", "raw", "Geometrias")

# We define the paths for the zipped shapefiles of districts and neighborhoods.
# If they don't exist, we download them from the official Madrid Geoportal.
# After downloading, we unzip them to the specified folder.
ruta_zip <- here(carpeta_geo, "Distritos.zip")
ruta_zip_barrios <- here(carpeta_geo, "Barrios.zip")

if (!dir.exists(carpeta_geo)) dir.create(carpeta_geo, recursive = TRUE)

if (!file.exists(ruta_zip)) {
  url_distritos <- "https://geoportal.madrid.es/fsdescargas/IDEAM_WBGEOPORTAL/LIMITES_ADMINISTRATIVOS/Distritos/Distritos.zip"
  download.file(url_distritos, destfile = ruta_zip, mode = "wb", quiet = TRUE)
}
unzip(ruta_zip, exdir = carpeta_geo)

if (!file.exists(ruta_zip_barrios)) {
  url_barrios <- "https://geoportal.madrid.es/fsdescargas/IDEAM_WBGEOPORTAL/LIMITES_ADMINISTRATIVOS/Barrios/Barrios.zip"
  download.file(url_barrios, destfile = ruta_zip_barrios, mode = "wb", quiet = TRUE)
}
unzip(ruta_zip_barrios, exdir = carpeta_geo)

# ignore.case = TRUE to avoid problems with uppercase/lowercase in the filenames.
# We take the first match just in case there are multiple (e.g., with different cases).

archivo_shp_distritos <- list.files(carpeta_geo,
  pattern = "Distritos.*\\.shp$",
  full.names = TRUE, ignore.case = TRUE
)[1]
archivo_shp_barrios <- list.files(carpeta_geo,
  pattern = "Barrios.*\\.shp$",
  full.names = TRUE, ignore.case = TRUE
)[1]

if (is.na(archivo_shp_distritos)) stop("No se encontró el shapefile de distritos en: ", carpeta_geo)
if (is.na(archivo_shp_barrios)) stop("No se encontró el shapefile de barrios en: ", carpeta_geo)

mapa_distritos <- st_read(archivo_shp_distritos, quiet = TRUE)
mapa_barrios <- st_read(archivo_shp_barrios, quiet = TRUE)

# Standarize the names to Title Case (e.g., "CENTRO" -> "centro") for better
# matching later
mapa_distritos$distrito <- tolower(mapa_distritos$NOMBRE)
mapa_barrios$barrio <- tolower(mapa_barrios$NOMBRE)

mapa_distritos <- mapa_distritos[, c("distrito", "geometry")]
mapa_barrios <- mapa_barrios[, c("barrio", "geometry")]

cat("Mapas cargados:", nrow(mapa_distritos), "distritos |", nrow(mapa_barrios), "barrios\n")

# ==============================================================================
# 2. LOOP MONTH PROCESS - CLEANING + SPATIAL JOIN + IMPUTATION + AGGREGATION
# ==============================================================================

# We create three output folders for the three levels of resolution (hourly,
# daily by district, daily by neighborhood)

# Estructura de salida: <base>/<escala>/<anio>/. Agrupar primero por escala y
# despues por anio permite leer una escala completa de varios anios con un solo
# list.files(), que es como la consumen los analisis.
carpeta_base <- here("data", "processed", "Trafico")

# Sufijo de version: se escribe en ficheros NUEVOS para no pisar los existentes.
# Poner "" para volver a escribir sobre los originales.
SUFIJO <- "1"

# Con TRUE el script se detiene si el fichero de destino ya existe.
PROTEGER_EXISTENTES <- TRUE

escalas <- c(
  horario = "Horario",
  hora_distrito = "Horario_Distrito",
  hora_barrio = "Horario_Barrio",
  distrito = "Diario_Distrito",
  barrio = "Diario_Barrio",
  mensual_distrito = "Mensual_Distrito",
  mensual_barrio = "Mensual_Barrio"
)

# Helper unico de guardado: con siete llamadas repartidas, el sufijo y la
# comprobacion no pueden quedar a criterio de cada una.
guardar <- function(objeto, carpeta, nombre) {
  ruta <- file.path(carpeta, paste0(nombre, SUFIJO, ".rds"))
  if (PROTEGER_EXISTENTES && file.exists(ruta)) {
    stop(
      "Ya existe el fichero de destino:
  ", ruta,
      "
Cambia SUFIJO o pon PROTEGER_EXISTENTES <- FALSE."
    )
  }
  dir.create(dirname(ruta), recursive = TRUE, showWarnings = FALSE)
  saveRDS(objeto, ruta)
}

procesar_anio_trafico <- function(anio) {
  anio <- as.integer(anio)

  cat("\n==================================================\n")
  cat("PROCESANDO ANIO:", anio, "\n")

  carpetas_salida <- lapply(escalas, function(e) file.path(carpeta_base, e, anio))
  invisible(lapply(carpetas_salida, dir.create, recursive = TRUE, showWarnings = FALSE))

  # Loop through each month, process the data, and save the three resolutions
  # (hourly, daily by district, daily by neighborhood)

  for (mes in meses) {
    cat("\n--------------------------------------------------\n")
    cat("Procesando MES:", mes, "\n")

    ruta_trafico_mes <- here(
      "data", "raw", "Datos_trafico",
      paste0("Datos_trafico_", anio),
      paste0(mes, "_", anio),
      paste0(mes, "_", anio, ".csv")
    )

    ruta_ubica_mes <- here(
      "data", "raw", "Datos_trafico",
      paste0("Detectores_", anio),
      paste0(mes, ".csv")
    )

    # Check if both files exist before processing
    if (file.exists(ruta_trafico_mes) & file.exists(ruta_ubica_mes)) {
      # Load the raw traffic data and the detector locations for the month
      dt_raw <- fread(ruta_trafico_mes, sep = ";")
      dt_ubi <- fread(ruta_ubica_mes, sep = ";")

      # 1) Clean and impute hourly data with spatial joins to assign district and neighborhood to each sensor reading
      cat("Empieza el proceso de limpieza de los datos horarios")
      dt_horario <- limpiar_trafico_espacial_horario(dt_raw,
        dt_ubicaciones = dt_ubi,
        mapa_distritos = mapa_distritos,
        mapa_barrios   = mapa_barrios
      )


      # 2) Aggregation to hourly by zone -> list(distrito = data.table, barrio = data.table)
      cat("Agregacion Horario zona")
      lista_horario_zona <- agregar_trafico_horario_zona(dt_horario)

      # Rejilla horaria completa. El "by" de la agregacion solo crea grupos para
      # las combinaciones observadas: si ningun sensor de la zona reporta una
      # hora, esa fila no existe y desaparece del denominador. Las zonas son
      # fijas (21 distritos, 127 barrios) y el calendario tambien, asi que la
      # rejilla es inequivoca.
      dias_mes <- sort(unique(dt_horario$FECHA))
      lista_horario_zona$distrito <- completar_rejilla_trafico(
        lista_horario_zona$distrito, "distrito", dias_mes
      )
      lista_horario_zona$barrio <- completar_rejilla_trafico(
        lista_horario_zona$barrio, "barrio", dias_mes
      )

      cat(sprintf(
        "  Rejilla horaria: +%d filas distrito, +%d barrio
",
        sum(lista_horario_zona$distrito$fila_ausente),
        sum(lista_horario_zona$barrio$fila_ausente)
      ))

      # 3) Aggregation to daily -> list(distrito = data.table, barrio = data.table)

      lista_diario <- agregar_trafico_diario(dt_horario)

      # 4) Aggregation to monthly -> list(distrito = data.table, barrio = data.table))

      lista_mensual <- agregar_trafico_mensual(
        dt_diario_distrito = lista_diario$distrito,
        dt_diario_barrio   = lista_diario$barrio,
        umbral_na          = 0.2
      )

      # Guardado de las siete resoluciones. El nombre ya no lleva el anio: lo
      # aporta la carpeta.
      guardar(
        dt_horario, carpetas_salida$horario,
        sprintf("Trafico_Horario_%s_", mes)
      )
      guardar(
        lista_horario_zona$distrito, carpetas_salida$hora_distrito,
        sprintf("Trafico_Horario_Distrito_%s_", mes)
      )
      guardar(
        lista_horario_zona$barrio, carpetas_salida$hora_barrio,
        sprintf("Trafico_Horario_Barrio_%s_", mes)
      )
      guardar(
        lista_diario$distrito, carpetas_salida$distrito,
        sprintf("Trafico_Distrito_%s_", mes)
      )
      guardar(
        lista_diario$barrio, carpetas_salida$barrio,
        sprintf("Trafico_Barrio_%s_", mes)
      )
      guardar(
        lista_mensual$distrito, carpetas_salida$mensual_distrito,
        sprintf("Trafico_Mensual_Distrito_%s_", mes)
      )
      guardar(
        lista_mensual$barrio, carpetas_salida$mensual_barrio,
        sprintf("Trafico_Mensual_Barrio_%s_", mes)
      )

      cat("ÉXITO:", mes, "guardado (horario raw, horario por zona, diario y mensual).\n")

      rm(
        dt_raw, dt_ubi, dt_horario, lista_horario_zona,
        lista_diario, lista_mensual, dias_mes
      )
      gc()
    } else {
      cat("AVISO: Faltan archivos para el mes de", mes, ". Saltando...\n")
    }
  }

  # ==============================================================================
  # 3.Final task: save the datasets and create the anual data set
  # ==============================================================================
  cat("\n==================================================\n")
  cat("INICIANDO ...\n")

  # 1.Define the route
  configuraciones <- list(
    list(
      nombre      = "Horario_Distrito",
      carpeta     = carpetas_salida$hora_distrito,
      orden       = c("FECHA", "HORA", "distrito"),
      archivo_out = paste0("trafico_madrid_", anio, "_horario_distrito")
    ),
    list(
      nombre      = "Horario_Barrio",
      carpeta     = carpetas_salida$hora_barrio,
      orden       = c("FECHA", "HORA", "barrio"),
      archivo_out = paste0("trafico_madrid_", anio, "_horario_barrio")
    ),
    list(
      nombre      = "Diario_Distrito",
      carpeta     = carpetas_salida$distrito,
      orden       = c("FECHA", "distrito"),
      archivo_out = paste0("trafico_madrid_", anio, "_diario_distrito")
    ),
    list(
      nombre      = "Diario_Barrio",
      carpeta     = carpetas_salida$barrio,
      orden       = c("FECHA", "barrio"),
      archivo_out = paste0("trafico_madrid_", anio, "_diario_barrio")
    ),
    list(
      nombre      = "Mensual_Distrito",
      carpeta     = carpetas_salida$mensual_distrito,
      orden       = c("MES", "distrito"),
      archivo_out = paste0("trafico_madrid_", anio, "_mensual_distrito")
    ),
    list(
      nombre      = "Mensual_Barrio",
      carpeta     = carpetas_salida$mensual_barrio,
      orden       = c("MES", "barrio"),
      archivo_out = paste0("trafico_madrid_", anio, "_mensual_barrio")
    )
  )

  # 2.Loop
  for (conf in configuraciones) {
    # Search for all RDS files in the specified folder
    archivos <- list.files(conf$carpeta, pattern = "\\.rds$", full.names = TRUE)
    archivos <- archivos[
      startsWith(basename(archivos), "Trafico_") &
        endsWith(basename(archivos), paste0(SUFIJO, ".rds"))
    ]

    if (length(archivos) > 0) {
      cat(sprintf("\n--- Compilando Dataset %s (%d meses encontrados) ---\n", conf$nombre, length(archivos)))

      # Combine all RDS files into a single data.table and sort by the specified columns
      dt_anual <- rbindlist(lapply(sort(archivos), readRDS))
      setorderv(dt_anual, conf$orden) # setorderv permite pasar un vector de nombres de columnas

      cat("Filas totales:", nrow(dt_anual), "\n")
      print(head(dt_anual, 3))

      # El compilado anual se guarda junto a los meses que lo componen, en la
      # misma carpeta <escala>/<anio>/, y respeta el sufijo de version.
      ruta_guardado <- file.path(
        conf$carpeta,
        paste0(conf$archivo_out, SUFIJO, ".rds")
      )
      if (PROTEGER_EXISTENTES && file.exists(ruta_guardado)) {
        stop(
          "Ya existe el fichero de destino:\n  ", ruta_guardado,
          "\nCambia SUFIJO o pon PROTEGER_EXISTENTES <- FALSE."
        )
      }
      saveRDS(dt_anual, ruta_guardado)
      cat("Guardado con éxito en:", ruta_guardado, "\n")

      # Clean up memory
      rm(dt_anual)
      gc()
    } else {
      cat(sprintf("\nAVISO: No se encontraron archivos para %s en %s\n", conf$nombre, conf$carpeta))
    }
  }

  cat("\nPREPROCESAMIENTO FINALIZADO PARA", anio, ".\n")
}

for (anio in anios_procesar) {
  procesar_anio_trafico(anio)
}

cat("\nPREPROCESAMIENTO FINALIZADO.\n")

check1 <- readRDS(here("data", "processed", "Trafico", "Diario_Barrio", "2021", "trafico_madrid_2021_diario_barrio1.rds"))
View(check1)
