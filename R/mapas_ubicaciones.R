# ==============================================================================
# MAPA DE UBICACIONES: ESTACIONES METEOROLÓGICAS Y DE CONTAMINACIÓN — 2025
# ==============================================================================

library(data.table)
library(sf)
library(ggplot2)
library(here)

source(here("R", "dictionaries.R"))      # nombres_estaciones_aire, nombres_estaciones_clima
source(here("R", "cleaning_functions.R"))

# ==============================================================================
# 1. MAPAS BASE (distritos y barrios)
# ==============================================================================
carpeta_geo <- here("data", "raw", "Geometrias")

archivo_shp_distritos <- list.files(carpeta_geo, pattern = "Distritos.*\\.shp$",
                                    full.names = TRUE, ignore.case = TRUE)[1]
archivo_shp_barrios   <- list.files(carpeta_geo, pattern = "Barrios.*\\.shp$",
                                    full.names = TRUE, ignore.case = TRUE)[1]

mapa_distritos <- st_read(archivo_shp_distritos, quiet = TRUE)
mapa_barrios   <- st_read(archivo_shp_barrios,   quiet = TRUE)

mapa_distritos$distrito <- tools::toTitleCase(tolower(mapa_distritos$NOMBRE))
mapa_barrios$barrio     <- tools::toTitleCase(tolower(mapa_barrios$NOMBRE))

# Proyectar a WGS84
mapa_distritos_wgs <- st_transform(mapa_distritos, crs = 4326)
mapa_barrios_wgs   <- st_transform(mapa_barrios,   crs = 4326)

# Colorear barrios por distrito: usamos el centroide de cada barrio para garantizar
# exactamente una asignación (st_intersects sobre polígonos puede dar duplicados en bordes)
centroides_barrios <- st_centroid(mapa_barrios_wgs["barrio"])
barrio_distrito    <- st_join(centroides_barrios, mapa_distritos_wgs["distrito"],
                              join = st_intersects)
mapa_barrios_wgs$distrito <- barrio_distrito$distrito

# ==============================================================================
# 2. ESTACIONES METEOROLÓGICAS (2025)
# ==============================================================================
# fread falla en Windows con espacios en la ruta → usamos read.csv
ruta_est_metereo <- here("data", "raw", "Datos metereologicos", "Estaciones_2019", "estaciones.csv")
dt_ubi_metereo <- as.data.table(read.csv(ruta_est_metereo, sep = ";", fileEncoding = "latin1"))

# Filtrar estaciones del diccionario 2025
dt_ubi_metereo_filtrado <- dt_ubi_metereo[CODIGO_CORTO %in% as.integer(names(nombres_estaciones_clima))]

# Usamos COORDENADA_X/Y_ETRS89 (UTM zona 30N) porque LONGITUD/LATITUD tiene
# formato inconsistente entre estaciones antiguas y nuevas
dt_ubi_metereo_filtrado[, X_utm := as.numeric(gsub(",", ".", COORDENADA_X_ETRS89))]
dt_ubi_metereo_filtrado[, Y_utm := as.numeric(gsub(",", ".", COORDENADA_Y_ETRS89))]

sf_metereo <- st_as_sf(dt_ubi_metereo_filtrado[!is.na(X_utm)],
                       coords = c("X_utm", "Y_utm"),
                       crs = 25830) |>
  st_transform(crs = 4326)

cat("Estaciones meteorológicas 2025:", nrow(sf_metereo), "\n")

# ==============================================================================
# 3. ESTACIONES DE CONTAMINACIÓN DEL AIRE (2025)
# ==============================================================================
dt_ubi_aire <- fread(here("data", "raw", "Datos_contaminacion", "Estaciones", "datos.csv"))

dt_ubi_aire_filtrado <- dt_ubi_aire[CODIGO_CORTO %in% as.integer(names(nombres_estaciones_aire))]

sf_aire <- st_as_sf(dt_ubi_aire_filtrado,
                    coords = c("LONGITUD", "LATITUD"),
                    crs = 4326)

cat("Estaciones contaminación 2025:", nrow(sf_aire), "\n")

# ==============================================================================
# 4. MAPA: BARRIOS (coloreados por distrito) + ESTACIONES
# ==============================================================================
mapa_estaciones <- ggplot() +
  # Barrios coloreados por distrito
  geom_sf(data = mapa_barrios_wgs,
          aes(fill = distrito),
          color = "white", linewidth = 0.15, alpha = 0.6) +
  # Bordes de distritos encima (más gruesos)
  geom_sf(data = mapa_distritos_wgs,
          fill = NA, color = "gray20", linewidth = 0.6) +
  # Estaciones de contaminación
  geom_sf(data  = sf_aire,
          aes(shape = "Contaminación"),
          color = "firebrick", fill = "firebrick", size = 3) +
  # Estaciones meteorológicas
  geom_sf(data  = sf_metereo,
          aes(shape = "Meteorología"),
          color = "steelblue", fill = "steelblue", size = 3) +
  scale_shape_manual(
    name   = "Tipo de estación",
    values = c("Contaminación" = 21, "Meteorología" = 24)
  ) +
  guides(
    fill  = guide_legend(title = "Distrito", ncol = 2,
                         override.aes = list(alpha = 0.7, color = NA)),
    shape = guide_legend(title = "Tipo de estación")
  ) +
  labs(
    title    = "Estaciones meteorológicas y de contaminación — Madrid 2025",
    subtitle = paste0(nrow(sf_aire), " estaciones de contaminación  |  ",
                      nrow(sf_metereo), " estaciones meteorológicas"),
    caption  = "Barrios coloreados por distrito. Bordes negros: límites de distrito.",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text     = element_text(size = 8),
    axis.text       = element_blank(),
    panel.grid      = element_blank()
  )


# ==============================================================================
# 5. GUARDADO
# ==============================================================================
carpeta_graficos <- here("results", "figures")
if (!dir.exists(carpeta_graficos)) dir.create(carpeta_graficos, recursive = TRUE)

ggsave(filename = file.path(carpeta_graficos, "mapa_estaciones_madrid_2025.png"),
       plot     = mapa_estaciones,
       width    = 12, height = 10, dpi = 300, bg = "white")

cat("Mapa guardado en:", file.path(carpeta_graficos, "mapa_estaciones_madrid_2025.png"), "\n")

