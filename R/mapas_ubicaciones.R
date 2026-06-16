# ==============================================================================
# GENERACIÓN DE 3 MAPAS: CONTAMINACIÓN, METEOROLOGÍA Y FUSIÓN DUAL
# ==============================================================================

library(data.table)
library(sf)
library(ggplot2)
library(ggrepel)
library(here)

source(here("R", "dictionaries.R"))      
source(here("R", "cleaning_functions.R"))

# ==============================================================================
# 1. CARGA Y PROCESAMIENTO ESPACIAL DE DISTRITOS Y BARRIOS
# ==============================================================================
carpeta_geo <- here("data", "raw", "Geometrias")

mapa_distritos <- st_read(list.files(carpeta_geo, "Distritos.*\\.shp$", full.names = TRUE)[1], quiet = TRUE) |> st_transform(4326)
mapa_barrios   <- st_read(list.files(carpeta_geo, "Barrios.*\\.shp$", full.names = TRUE)[1], quiet = TRUE) |> st_transform(4326)

mapa_distritos$distrito <- tools::toTitleCase(tolower(mapa_distritos$NOMBRE))
mapa_barrios$barrio     <- tools::toTitleCase(tolower(mapa_barrios$NOMBRE))
mapa_barrios$distrito   <- st_join(st_centroid(mapa_barrios["barrio"]), mapa_distritos["distrito"])$distrito

# ==============================================================================
# 2. CARGA Y LIMPIEZA DE ESTACIONES
# ==============================================================================
# Meteorología
sf_metereo <- as.data.table(read.csv(here("data", "raw", "Datos metereologicos", "Estaciones_2019", "estaciones.csv"), sep = ";", fileEncoding = "latin1"))[
  CODIGO_CORTO %in% as.integer(names(nombres_estaciones_clima))
][, `:=`(
  Nombre_Estacion = nombres_estaciones_clima[as.character(CODIGO_CORTO)],
  X_utm = as.numeric(gsub(",", ".", COORDENADA_X_ETRS89)),
  Y_utm = as.numeric(gsub(",", ".", COORDENADA_Y_ETRS89))
)][!is.na(X_utm)] |> 
  st_as_sf(coords = c("X_utm", "Y_utm"), crs = 25830) |> st_transform(4326)

# Contaminación
sf_aire <- fread(here("data", "raw", "Datos_contaminacion", "Estaciones", "datos.csv"))[
  CODIGO_CORTO %in% as.integer(names(nombres_estaciones_aire))
][, Nombre_Estacion := nombres_estaciones_aire[as.character(CODIGO_CORTO)]] |> 
  st_as_sf(coords = c("LONGITUD", "LATITUD"), crs = 4326)

cat("✅ Datos listos: ", nrow(sf_aire), " Contaminación | ", nrow(sf_metereo), " Clima\n")

# ==============================================================================
# 3. CONSTRUCCIÓN DE LA PLANTILLA BASE DEL MAPA
# ==============================================================================
mapa_base <- ggplot() +
  geom_sf(data = mapa_barrios, aes(fill = distrito), color = "white", linewidth = 0.15, alpha = 0.6) +
  geom_sf(data = mapa_distritos, fill = NA, color = "gray20", linewidth = 0.6) +
  guides(fill = guide_legend(title = "Distrito", ncol = 2, override.aes = list(alpha = 0.7, color = NA))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right", axis.text = element_blank(), panel.grid = element_blank())

# ==============================================================================
# 4. MAPA 1: SOLO CONTAMINACIÓN
# ==============================================================================
mapa_1_aire <- mapa_base +
  geom_sf(data = sf_aire, color = "#c0392b", fill = "#e74c3c", size = 2.5, shape = 21) +
  geom_text_repel(
    data = sf_aire, aes(label = Nombre_Estacion, geometry = geometry), stat = "sf_coordinates",
    size = 3.2, fontface = "bold", color = "#c0392b", bg.color = "white", bg.r = 0.15, 
    box.padding = 0.5, point.padding = 0.4, min.segment.length = 0
  ) +
  labs(title = "Red de Estaciones de Contaminación Atmosférica",
       subtitle = paste("Ciudad de Madrid 2025 |", nrow(sf_aire), "estaciones activas"), x = NULL, y = NULL)

# ==============================================================================
# 5. MAPA 2: SOLO METEOROLOGÍA
# ==============================================================================
mapa_2_clima <- mapa_base +
  geom_sf(data = sf_metereo, color = "#2980b9", fill = "#3498db", size = 2.5, shape = 24) +
  geom_text_repel(
    data = sf_metereo, aes(label = Nombre_Estacion, geometry = geometry), stat = "sf_coordinates",
    size = 3.2, fontface = "bold", color = "#2980b9", bg.color = "white", bg.r = 0.15, 
    box.padding = 0.5, point.padding = 0.4, min.segment.length = 0
  ) +
  labs(title = "Red de Estaciones Meteorológicas",
       subtitle = paste("Ciudad de Madrid 2025 |", nrow(sf_metereo), "estaciones activas"), x = NULL, y = NULL)

# ==============================================================================
# 6. MAPA 3: FUSIÓN DE AMBAS REDES (ETIQUETAS POR COINCIDENCIA DE NOMBRE)
# ==============================================================================
# Identificamos qué nombres exactos están en ambas listas
nombres_comunes <- intersect(sf_aire$Nombre_Estacion, sf_metereo$Nombre_Estacion)

# Dividimos las estaciones en 3 grupos lógicos para etiquetar
sf_etiquetas_aire  <- sf_aire[!sf_aire$Nombre_Estacion %in% nombres_comunes, ]
sf_etiquetas_clima <- sf_metereo[!sf_metereo$Nombre_Estacion %in% nombres_comunes, ]
sf_etiquetas_dual  <- sf_aire[sf_aire$Nombre_Estacion %in% nombres_comunes, ] # Usamos coordenadas de aire como ancla

mapa_3_dual <- mapa_base +
  # Dibujamos TODOS los puntos geométricos superpuestos
  geom_sf(data = sf_aire, aes(shape = "Contaminación"), color = "#c0392b", fill = "#e74c3c", size = 2.5) +
  geom_sf(data = sf_metereo, aes(shape = "Meteorología"), color = "#2980b9", fill = "#3498db", size = 2.5) +
  
  # 1. Etiquetas EXCLUSIVAS de Contaminación (Rojas)
  geom_text_repel(data = sf_etiquetas_aire, aes(label = Nombre_Estacion, geometry = geometry), stat = "sf_coordinates",
                  size = 2.8, fontface = "bold", color = "#c0392b", bg.color = "white", bg.r = 0.15, box.padding = 0.4) +
  
  # 2. Etiquetas EXCLUSIVAS de Meteorología (Azules)
  geom_text_repel(data = sf_etiquetas_clima, aes(label = Nombre_Estacion, geometry = geometry), stat = "sf_coordinates",
                  size = 2.8, fontface = "bold", color = "#2980b9", bg.color = "white", bg.r = 0.15, box.padding = 0.4) +
  
  # 3. Etiquetas DOBLES/COMPARTIDAS (Amarillo Dorado)
  geom_text_repel(data = sf_etiquetas_dual, aes(label = Nombre_Estacion, geometry = geometry), stat = "sf_coordinates",
                  size = 2.8, fontface = "bold", color = "#d4ac0d", bg.color = "white", bg.r = 0.15, box.padding = 0.4) +
  
  scale_shape_manual(name = "Tipo de estación", values = c("Contaminación" = 21, "Meteorología" = 24)) +
  labs(title = "Red Dual de Estaciones: Meteorología y Contaminación",
       subtitle = "Ciudad de Madrid 2025 | En amarillo estaciones que miden ambos parámetros", x = NULL, y = NULL)

# ==============================================================================
# 7. GUARDADO MÚLTIPLE
# ==============================================================================
carpeta_graficos <- here("results", "figures")
if (!dir.exists(carpeta_graficos)) dir.create(carpeta_graficos, recursive = TRUE)

ggsave(file.path(carpeta_graficos, "mapa_01_solo_aire.png"), plot = mapa_1_aire, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(file.path(carpeta_graficos, "mapa_02_solo_clima.png"), plot = mapa_2_clima, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(file.path(carpeta_graficos, "mapa_03_dual.png"), plot = mapa_3_dual, width = 15, height = 11, dpi = 300, bg = "white")

cat("✅ Los 3 mapas se han guardado correctamente en la carpeta:", carpeta_graficos, "\n")

# ==============================================================================
# 8. MAPA 4: COBERTURA POR VARIABLE CLIMATOLÓGICA (FACETADO)
# ==============================================================================
vars_clima_cols <- c("Dir.Viento", "Humedad_Relativa", "Precipitaciones",
                     "Presion Barométrica", "Temperatura", "Velocidad Viento")

etiquetas_vars <- c(
  "Dir.Viento"          = "Dirección del Viento",
  "Humedad_Relativa"    = "Humedad Relativa",
  "Precipitaciones"     = "Precipitaciones",
  "Presion Barométrica" = "Presión Barométrica",
  "Temperatura"         = "Temperatura",
  "Velocidad Viento"    = "Velocidad del Viento"
)

# Presencia de cada variable por estación (TRUE/FALSE)
presencia <- dt_meteo[, lapply(.SD, function(x) sum(!is.na(x)) > 0),
                      by = ESTACION, .SDcols = vars_clima_cols]

# Coordenadas UTM (X_km/Y_km están en km → convertir a metros)
coords_estaciones <- dt_meteo[, .(X_m = mean(X_km, na.rm = TRUE) * 1000,
                                   Y_m = mean(Y_km, na.rm = TRUE) * 1000),
                               by = ESTACION]
presencia <- merge(presencia, coords_estaciones, by = "ESTACION")

# Todas las estaciones (puntos de referencia en gris)
sf_todas <- st_as_sf(presencia, coords = c("X_m", "Y_m"), crs = 25830) |>
  st_transform(4326)

# Formato largo: solo estaciones que miden cada variable
dt_long_vars <- melt(presencia,
                     id.vars       = c("ESTACION", "X_m", "Y_m"),
                     measure.vars  = vars_clima_cols,
                     variable.name = "Variable",
                     value.name    = "Mide")[Mide == TRUE]

dt_long_vars[, Variable_label := factor(
  etiquetas_vars[as.character(Variable)],
  levels = c("Temperatura", "Humedad Relativa", "Precipitaciones",
             "Dirección del Viento", "Velocidad del Viento", "Presión Barométrica")
)]

sf_long_vars <- st_as_sf(dt_long_vars, coords = c("X_m", "Y_m"), crs = 25830) |>
  st_transform(4326)

# Paleta de colores
colores_var <- c(
  "Temperatura"          = "#e74c3c",
  "Humedad Relativa"     = "#3498db",
  "Precipitaciones"      = "#2ecc71",
  "Dirección del Viento" = "#9b59b6",
  "Velocidad del Viento" = "#e67e22",
  "Presión Barométrica"  = "#1abc9c"
)

mapa_4_variables <- ggplot() +
  geom_sf(data = mapa_barrios,   fill = "#f0f0f0", color = "white", linewidth = 0.1) +
  geom_sf(data = mapa_distritos, fill = NA,        color = "gray40", linewidth = 0.5) +
  geom_sf(data = sf_todas,       color = "gray70", size = 1.5, shape = 16) +
  geom_sf(data = sf_long_vars,   aes(color = Variable_label), size = 3, shape = 16) +
  geom_text_repel(
    data              = sf_long_vars,
    aes(label = ESTACION, geometry = geometry, color = Variable_label),
    stat              = "sf_coordinates",
    size              = 2.2,
    bg.color          = "white",
    bg.r              = 0.12,
    box.padding       = 0.3,
    point.padding     = 0.2,
    min.segment.length = 0.3,
    show.legend       = FALSE
  ) +
  scale_color_manual(values = colores_var, guide = "none") +
  facet_wrap(~ Variable_label, ncol = 3) +
  labs(
    title    = "Cobertura de la Red de Estaciones Meteorológicas de Madrid",
    subtitle = "Puntos grises: estaciones sin datos para esa variable | Puntos coloreados: estaciones con datos",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text        = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray95", color = NA),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 8, color = "gray40")
  )

ggsave(file.path(carpeta_graficos, "mapa_04_variables_climatologicas.png"),
       plot = mapa_4_variables, width = 16, height = 11, dpi = 300, bg = "white")

cat("✅ Mapa 4 guardado:", file.path(carpeta_graficos, "mapa_04_variables_climatologicas.png"), "\n")

