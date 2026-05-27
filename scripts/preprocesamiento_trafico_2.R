# ==============================================================================
# PREPROCESAMIENTO DE TRÁFICO: AGREGACIÓN DIARIA MEDIANTE CRUCE ESPACIAL
# ==============================================================================

library(data.table)
library(sf)
library(here)
source(here("R", "cleaning_functions.R"))

# ----------------- 0. DESCARGA Y CARGA DEL MAPA OFICIAL DE DISTRITOS -----------------
cat("⏳ Preparando mapa oficial de distritos de Madrid...\n")
carpeta_geo <- here("data", "raw", "Geometrias")
ruta_zip    <- here(carpeta_geo, "Distritos.zip")

if (!dir.exists(carpeta_geo)) dir.create(carpeta_geo, recursive = TRUE)

if (!file.exists(ruta_zip)) {
  url_distritos <- "https://geoportal.madrid.es/fsdescargas/IDEAM_WBGEOPORTAL/LIMITES_ADMINISTRATIVOS/Distritos/Distritos.zip"
  download.file(url_distritos, destfile = ruta_zip, mode = "wb", quiet = TRUE)
}

unzip(ruta_zip, exdir = carpeta_geo)
archivo_shp <- list.files(carpeta_geo, pattern = "\\.shp$", full.names = TRUE)[1]

mapa_distritos <- st_read(archivo_shp, quiet = TRUE)
# Estandarizamos el nombre a Título (Ej: "CENTRO" -> "Centro") para que cruce bien luego
mapa_distritos$distrito <- tools::toTitleCase(tolower(mapa_distritos$NOMBRE))


# ----------------- 2. BUCLE DE PROCESAMIENTO MENSUAL -----------------
# (El bucle de ejecución se mantiene igual, pero ahora ejecuta la función horaria)

meses <- c("Enero","Febrero", "Marzo", "Abril", "Mayo", "Junio", 
           "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre")

carpeta_salida <- here("data", "raw", "Datos_trafico", "Datos_limpios_2025")
if(!dir.exists(carpeta_salida)) dir.create(carpeta_salida)

for (mes in meses) {
  
  cat("\n--------------------------------------------------\n")
  print(paste("Procesando MES:", mes, "(Imputación Horaria)"))
  
  ruta_trafico_mes <- here("data", "raw", "Datos_trafico", "Datos_trafico_2025", paste0(mes, "_2025"), paste0(mes, "_2025.csv"))
  ruta_ubica_mes   <- here("data", "raw", "Datos_trafico", "Detectores_2025", paste0(mes, ".csv"))
  
  if (file.exists(ruta_trafico_mes) & file.exists(ruta_ubica_mes)) {
    
    dt_raw <- fread(ruta_trafico_mes, sep = ";")
    dt_ubi <- fread(ruta_ubica_mes, sep = ";")
    
    # Ejecutamos la función pesada
    dt_limpio_diario <- limpiar_trafico_espacial_horario(dt_raw, dt_ubicaciones = dt_ubi, mapa_poligonos = mapa_distritos)
    
    ruta_guardado <- file.path(carpeta_salida, paste0("Trafico_Areas_Diario_", mes, "_2025.rds"))
    saveRDS(dt_limpio_diario, ruta_guardado)
    
    print(paste("ÉXITO: Tráfico diario del mes de", mes, "guardado con éxito."))
    
    rm(dt_raw, dt_ubi, dt_limpio_diario)
    gc() # Vaciamos la memoria tras el esfuerzo intensivo
    
  } else {
    print(paste("⚠️ AVISO: Faltan archivos para el mes de", mes, ". Saltando..."))
  }
}
# ==============================================================================
# 3. APILACIÓN FINAL EN UN ÚNICO DATASET ANUAL
# ==============================================================================

# Buscamos automáticamente todos los archivos procesados
archivos_rds <- list.files(path = carpeta_salida, 
                           pattern = "^Trafico_Areas_Diario_.*\\.rds$", 
                           full.names = TRUE)

cat("\n📂 Archivos de tráfico encontrados:", length(archivos_rds), "\n")

# Leemos todos los archivos y los apilamos
lista_meses <- lapply(archivos_rds, readRDS)
dt_trafico_anual_2025 <- rbindlist(lista_meses)

# Ordenamos cronológicamente
setorder(dt_trafico_anual_2025, FECHA, distrito)

cat("📊 Filas totales del tráfico anual por distrito:", nrow(dt_trafico_anual_2025), "\n")
print(head(dt_trafico_anual_2025))

# Guardado final unificado
saveRDS(dt_trafico_anual_2025, here("data", "processed", "trafico_madrid_2025_anual.rds"))
cat("💾 Dataset anual guardado en: data/processed/trafico_madrid_2025_anual.rds\n")

view

library(sf)
library(ggplot2)

# Convertir detectores a sf con el mismo CRS que el mapa (UTM zona 30N)
sf_medidores <- st_as_sf(dt_ubi, coords = c("utm_x", "utm_y"), crs = 25830)

ggplot() +
  geom_sf(data = mapa_distritos, aes(fill = distrito), alpha = 0.3, color = "gray30", linewidth = 0.4) +
  geom_sf(data = sf_medidores, size = 0.8, alpha = 0.6, color = "black") +
  labs(
    title    = "Medidores de tráfico por distrito — Madrid 2025",
    subtitle = paste("N =", nrow(sf_medidores), "detectores"),
    fill     = "Distrito",
    caption  = "Fuente: Detectores_2025/Enero.csv"
  ) +
  theme_minimal() +
  theme(legend.position  = "right",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 7))
