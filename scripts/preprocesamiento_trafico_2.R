
library(data.table)
library(sf)
library(here)
here::i_am("scripts/preprocesamiento_trafico_2.R")
source(here("R", "cleaning_functions.R"))

# --- Global Parameetrs ---

#===============================================================================
# Global Parameter in order to process the traffic data of Madrid for the year
# Traffic data is only available for the whole year, but I had to download it month by month
#===============================================================================
anio <- 2020 # Year of the traffic data to process
meses <- c("Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio",
           "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre")

# ==============================================================================
# 0. Load and charge the official maps of districts and neighborhoods of Madrid (shapefiles)
# ==============================================================================
cat("Preparando mapas oficiales de distritos y barrios de Madrid...\n")
carpeta_geo <- here("data", "raw", "Geometrias")

# We define the paths for the zipped shapefiles of districts and neighborhoods. 
# If they don't exist, we download them from the official Madrid Geoportal. 
# After downloading, we unzip them to the specified folder.
ruta_zip         <- here(carpeta_geo, "Distritos.zip")
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

archivo_shp_distritos <- list.files(carpeta_geo, pattern = "Distritos.*\\.shp$",
                                    full.names = TRUE, ignore.case = TRUE)[1]
archivo_shp_barrios   <- list.files(carpeta_geo, pattern = "Barrios.*\\.shp$",
                                    full.names = TRUE, ignore.case = TRUE)[1]

if (is.na(archivo_shp_distritos)) stop("No se encontró el shapefile de distritos en: ", carpeta_geo)
if (is.na(archivo_shp_barrios))   stop("No se encontró el shapefile de barrios en: ",   carpeta_geo)

mapa_distritos <- st_read(archivo_shp_distritos, quiet = TRUE)
mapa_barrios   <- st_read(archivo_shp_barrios,   quiet = TRUE)

# Standarize the names to Title Case (e.g., "CENTRO" -> "centro") for better
# matching later
mapa_distritos$distrito <- tolower(mapa_distritos$NOMBRE)
mapa_barrios$barrio     <- tolower(mapa_barrios$NOMBRE)

mapa_distritos <- mapa_distritos[, c("distrito", "geometry")]
mapa_barrios   <- mapa_barrios[,   c("barrio",   "geometry")]

cat("Mapas cargados:", nrow(mapa_distritos), "distritos |", nrow(mapa_barrios), "barrios\n")

# ==============================================================================
# 2. LOOP MONTH PROCESS - CLEANING + SPATIAL JOIN + IMPUTATION + AGGREGATION
# ==============================================================================

# We create three output folders for the three levels of resolution (hourly, 
# daily by district, daily by neighborhood)

carpeta_base <- here("data", "raw", "Datos_trafico", paste0("Datos_limpios_", anio))
carpetas_salida <- list(
  horario          = file.path(carpeta_base, "Horario"),
  hora_distrito    = file.path(carpeta_base, "Horario_Distrito"),
  hora_barrio      = file.path(carpeta_base, "Horario_Barrio"),
  distrito         = file.path(carpeta_base, "Diario_Distrito"),
  barrio           = file.path(carpeta_base, "Diario_Barrio"),
  mensual_distrito = file.path(carpeta_base, "Mensual_Distrito"),
  mensual_barrio   = file.path(carpeta_base, "Mensual_Barrio")
)

lapply(carpetas_salida, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

# Loop through each month, process the data, and save the three resolutions 
# (hourly, daily by district, daily by neighborhood)

for (mes in meses) {

  cat("\n--------------------------------------------------\n")
  cat("Procesando MES:", mes, "\n")

  ruta_trafico_mes <- here("data", "raw", "Datos_trafico", 
                           paste0("Datos_trafico_", anio),
                           paste0(mes, "_", anio), 
                           paste0(mes, "_", anio, ".csv"))
  
  ruta_ubica_mes <- here("data", "raw", "Datos_trafico", 
                         paste0("Detectores_", anio), 
                         paste0(mes, ".csv"))

  # Check if both files exist before processing
  if (file.exists(ruta_trafico_mes) & file.exists(ruta_ubica_mes)) {
    
    # Load the raw traffic data and the detector locations for the month
    dt_raw <- fread(ruta_trafico_mes, sep = ";")
    dt_ubi <- fread(ruta_ubica_mes,   sep = ";")

    # 1) Clean and impute hourly data with spatial joins to assign district and neighborhood to each sensor reading
    cat("Empieza el proceso de limpieza de los datos horarios")
    dt_horario <- limpiar_trafico_espacial_horario(dt_raw,
                                                   dt_ubicaciones = dt_ubi,
                                                   mapa_distritos = mapa_distritos,
                                                   mapa_barrios   = mapa_barrios)
    
   

    # 2) Aggregation to hourly by zone -> list(distrito = data.table, barrio = data.table)
    cat("Agregacion Horario zona")
    lista_horario_zona <- agregar_trafico_horario_zona(dt_horario)

    # 3) Aggregation to daily -> list(distrito = data.table, barrio = data.table)
    
    lista_diario <- agregar_trafico_diario(dt_horario)
    
    # 4) Aggregation to monthly -> list(distrito = data.table, barrio = data.table))
    
    lista_mensual <- agregar_trafico_mensual(
      dt_diario_distrito = lista_diario$distrito,
      dt_diario_barrio   = lista_diario$barrio,
      umbral_na          = 0.2
    )
 
    # Guardado de las cinco resoluciones
    saveRDS(dt_horario,                    file.path(carpetas_salida$horario,          sprintf("Trafico_Horario_%s_%s.rds",                  mes, anio)))
    saveRDS(lista_horario_zona$distrito,   file.path(carpetas_salida$hora_distrito,    sprintf("Trafico_Horario_Distrito_%s_%s.rds",         mes, anio)))
    saveRDS(lista_horario_zona$barrio,     file.path(carpetas_salida$hora_barrio,      sprintf("Trafico_Horario_Barrio_%s_%s.rds",           mes, anio)))
    saveRDS(lista_diario$distrito,         file.path(carpetas_salida$distrito,         sprintf("Trafico_Distrito_%s_%s.rds",                 mes, anio)))
    saveRDS(lista_diario$barrio,           file.path(carpetas_salida$barrio,           sprintf("Trafico_Barrio_%s_%s.rds",                   mes, anio)))
    saveRDS(lista_mensual$distrito,        file.path(carpetas_salida$mensual_distrito, sprintf("Trafico_Mensual_Distrito_%s_%s.rds",         mes, anio)))
    saveRDS(lista_mensual$barrio,          file.path(carpetas_salida$mensual_barrio,   sprintf("Trafico_Mensual_Barrio_%s_%s.rds",           mes, anio)))

    cat("ÉXITO:", mes, "guardado (horario raw, horario por zona, diario y mensual).\n")

    rm(dt_raw, dt_ubi, dt_horario_barrio)
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

carpeta_salida_distrito <- file.path(carpeta_base, "Diario_Distrito")
carpeta_salida_barrio   <- file.path(carpeta_base, "Diario_Barrio")
# 1.Define the route 
configuraciones <- list(

  list(
    nombre      = "Horario_Distrito",
    carpeta     = file.path(carpeta_base, "Horario_Distrito"),
    orden       = c("FECHA", "HORA", "distrito"),
    archivo_out = paste0("trafico_madrid_", anio, "_horario_distrito.rds")
  ),
  list(
    nombre      = "Horario_Barrio",
    carpeta     = file.path(carpeta_base, "Horario_Barrio"),
    orden       = c("FECHA", "HORA", "barrio"),
    archivo_out = paste0("trafico_madrid_", anio, "_horario_barrio.rds")
  ),
  list(
    nombre      = "Diario_Distrito",
    carpeta     = carpeta_salida_distrito,
    orden       = c("FECHA", "distrito"),
    archivo_out = paste0("trafico_madrid_", anio, "_diario_distrito.rds")
  ),
  list(
    nombre      = "Diario_Barrio",
    carpeta     = carpeta_salida_barrio,
    orden       = c("FECHA", "barrio"),
    archivo_out = paste0("trafico_madrid_", anio, "_diario_barrio.rds")
  ),
  list(
    nombre      = "Mensual_Distrito",
    carpeta     = file.path(carpeta_base, "Mensual_Distrito"),
    orden       = c("MES", "distrito"),
    archivo_out = paste0("trafico_madrid_", anio, "_mensual_distrito.rds")
  ),
  list(
    nombre      = "Mensual_Barrio",
    carpeta     = file.path(carpeta_base, "Mensual_Barrio"),
    orden       = c("MES", "barrio"),
    archivo_out = paste0("trafico_madrid_", anio, "_mensual_barrio.rds")
  )
)

# 2.Loop 
for (conf in configuraciones) {
  
  # Search for all RDS files in the specified folder
  archivos <- list.files(conf$carpeta, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(archivos) > 0) {
    cat(sprintf("\n--- Compilando Dataset %s (%d meses encontrados) ---\n", conf$nombre, length(archivos)))
    
    # Combine all RDS files into a single data.table and sort by the specified columns
    dt_anual <- rbindlist(lapply(archivos, readRDS))
    setorderv(dt_anual, conf$orden) # setorderv permite pasar un vector de nombres de columnas
    
    cat("Filas totales:", nrow(dt_anual), "\n")
    print(head(dt_anual, 3))
    
    # Save the final annual dataset in the processed folder
    ruta_guardado <- here("data", "processed", conf$archivo_out)
    saveRDS(dt_anual, ruta_guardado)
    cat("Guardado con éxito en:", ruta_guardado, "\n")
    
    # Clean up memory
    rm(dt_anual)
    gc()
    
  } else {
    cat(sprintf("\nAVISO: No se encontraron archivos para %s en %s\n", conf$nombre, conf$carpeta))
  }
}

cat("\nPREPROCESAMIENTO FINALIZADO.\n")

# ==============================================================================
# 4. Generate a map of the sensor network for visual inspection
# ==============================================================================

library(ggplot2)
library(sf)
library(here)

cat("\nGenerando cartografía de la red de sensores...\n")

ruta_ubica_mes   <- here("data", "raw", "Datos_trafico", paste0("Detectores_", anio), "Enero.csv")
dt_ubi <- fread(ruta_ubica_mes,   sep = ";")
# 1. Create the output folder for figures if it doesn't exist
carpeta_figuras <- here("output", "figures")
if (!dir.exists(carpeta_figuras)) dir.create(carpeta_figuras, recursive = TRUE)

# 2. Convert the sensor locations to an sf object for plotting
# We filter out any rows with NA coordinates to avoid errors in the spatial conversion

sf_medidores <- st_as_sf(dt_ubi[!is.na(utm_x)],
                         coords = c("utm_x", "utm_y"), 
                         crs = 25830)

# 3. Create a ggplot map with the districts and sensor points.
mapa_sensores <- ggplot() +
  # Background: Districts are filled with a light gray color and outlined with a slightly darker gray
  geom_sf(data = mapa_distritos, fill = "gray95", color = "gray60", linewidth = 0.3) +
  # Sensor points are plotted with a semi-transparent dark blue color to avoid cluttering the map
  geom_sf(data = sf_medidores, size = 0.6, alpha = 0.5, color = "#2c3e50") +
  # Add titles, subtitles, captions, and axis labels for clarity
  labs(
    title    = "Distribución Espacial de los Detectores de Tráfico",
    subtitle = sprintf("Ciudad de Madrid (Año %s) | N = %d sensores operativos", anio, nrow(sf_medidores)),
    caption  = "Fuente: Elaboración propia a partir del Portal de Datos Abiertos del Ayto. de Madrid",
    x = "Longitud (UTM)",
    y = "Latitud (UTM)"
  ) +
  
  theme_minimal(base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 14, margin = margin(b = 5)),
    plot.subtitle    = element_text(color = "gray40", size = 11, margin = margin(b = 15)),
    plot.caption     = element_text(color = "gray50", size = 8, hjust = 1, margin = margin(t = 10)),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    axis.title       = element_text(size = 9, color = "gray30"),
    axis.text        = element_text(size = 8, color = "gray50")
  )

print(mapa_sensores)

# 4. Save the map to a PNG. 
ruta_guardado_mapa <- file.path(carpeta_figuras, "mapa_distribucion_sensores",anio,".png")

ggsave(filename = ruta_guardado_mapa, 
       plot   = mapa_sensores, 
       width  = 10,    # Ancho en pulgadas
       height = 8,     # Alto en pulgadas
       dpi    = 300,   # 300 DPI es el estándar de publicación universitaria
       bg     = "white") # Fuerza el fondo blanco para que no salga transparente

cat("ÉXITO: Mapa cartográfico guardado en:", ruta_guardado_mapa, "\n")

