# ==============================================================================
# MAPA DE ESTACIONES DE MONITOREO DE NO₂ — Madrid 2025
# Por distritos · Diferenciadas por tipología · Con etiquetas de estación
# ==============================================================================

library(data.table)
library(sf)
library(ggplot2)
library(ggrepel)
library(here)


# ------------------------------------------------------------------------------
# 1. DATOS
# ------------------------------------------------------------------------------

# Distritos de Madrid (ETRS89 / UTM zona 30N — EPSG:25830)
distritos <- st_read(here("data", "raw", "geometrias", "DISTRITOS.shp"),
                     quiet = TRUE)

# Estaciones: coordenadas + tipología (una fila por estación)
dt_mens <- readRDS(here("data", "processed", "contaminacion", "mensual",
                         "aire_madrid_2025_log_No2_mensuales.rds"))

estaciones <- as.data.frame(
  dt_mens[, .(
    LONGITUD = unique(LONGITUD),
    LATITUD  = unique(LATITUD),
    NOM_TIPO = unique(NOM_TIPO)
  ), by = ESTACION]
)

# Convertir a sf (WGS84) → reprojectar al mismo CRS que los distritos
sf_estaciones <- st_as_sf(estaciones,
                           coords = c("LONGITUD", "LATITUD"),
                           crs    = 4326) |>
  st_transform(st_crs(distritos))

# Coordenadas XY explícitas para ggrepel
sf_estaciones$X <- st_coordinates(sf_estaciones)[, 1]
sf_estaciones$Y <- st_coordinates(sf_estaciones)[, 2]

# Orden fijo de tipologías en leyenda
sf_estaciones$NOM_TIPO <- factor(sf_estaciones$NOM_TIPO,
                                  levels = c("Urbana tr\u00e1fico",
                                             "Urbana fondo",
                                             "Suburbana"))

cat("Distritos cargados :", nrow(distritos), "\n")
cat("Estaciones cargadas:", nrow(sf_estaciones), "\n")


# ------------------------------------------------------------------------------
# 2. PALETAS
# ------------------------------------------------------------------------------

# 21 colores pastel — uno por distrito (orden alfabético)
colores_distritos <- setNames(
  hcl.colors(21, palette = "Pastel 1"),
  sort(distritos$NOMBRE)
)

# Estaciones: colour + shape sólidos (sin fill, para no competir con fill de distritos)
colores_estacion <- c(
  "Urbana tr\u00e1fico" = "#c0392b",
  "Urbana fondo"        = "#1a5276",
  "Suburbana"           = "#1e8449"
)
formas_estacion <- c(
  "Urbana tr\u00e1fico" = 16,   # círculo sólido
  "Urbana fondo"        = 17,   # triángulo sólido
  "Suburbana"           = 15    # cuadrado sólido
)


# ------------------------------------------------------------------------------
# 3. MAPA
# ------------------------------------------------------------------------------

mapa_estaciones <- ggplot() +

  # — Polígonos de distritos coloreados —
  geom_sf(data = distritos,
          aes(fill = NOMBRE),
          colour    = "white",
          linewidth = 0.55) +
  scale_fill_manual(
    values = colores_distritos,
    name   = "Distritos de Madrid",
    guide  = guide_legend(ncol = 2,
                          keywidth  = unit(0.45, "cm"),
                          keyheight = unit(0.38, "cm"),
                          override.aes = list(colour    = "grey70",
                                              linewidth = 0.3))
  ) +

  # — Puntos de estaciones (colour + shape, sin fill) —
  geom_point(data  = sf_estaciones,
             aes(x      = X,
                 y      = Y,
                 colour = NOM_TIPO,
                 shape  = NOM_TIPO),
             size   = 3.8,
             stroke = 0.5) +
  scale_colour_manual(
    values = colores_estacion,
    name   = "Tipolog\u00eda de estaci\u00f3n",
    guide  = guide_legend(keywidth  = unit(0.5, "cm"),
                          keyheight = unit(0.5, "cm"),
                          override.aes = list(size = 4))
  ) +
  scale_shape_manual(
    values = formas_estacion,
    name   = "Tipolog\u00eda de estaci\u00f3n"
  ) +

  # — Etiquetas de estaciones con fondo blanco (ggrepel) —
  geom_text_repel(data  = sf_estaciones,
                  aes(x = X, y = Y, label = ESTACION),
                  size          = 2.3,
                  colour        = "grey10",
                  fontface      = "bold",
                  bg.color      = "white",
                  bg.r          = 0.12,
                  box.padding   = 0.45,
                  point.padding = 0.3,
                  segment.color = "grey40",
                  segment.size  = 0.35,
                  max.overlaps  = 40,
                  seed          = 7231) +

  labs(
    title    = "Estaciones de Monitoreo de NO\u2082 \u2014 Madrid 2025",
    subtitle = "24 estaciones \u00b7 Red Municipal de Calidad del Aire",
    caption  = "Fuente: Ayuntamiento de Madrid \u00b7 Cartograf\u00eda: DISTRITOS.shp (ETRS89 / UTM 30N)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle    = element_text(colour = "grey40", size = 10, hjust = 0),
    plot.caption     = element_text(colour = "grey55", size = 8),
    axis.title       = element_blank(),
    axis.text        = element_text(size = 7, colour = "grey60"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    legend.spacing.y = unit(0.15, "cm"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  coord_sf(crs = st_crs(distritos))


# ------------------------------------------------------------------------------
# 4. GUARDAR
# ------------------------------------------------------------------------------

carpeta_out <- here("outputs", "mapas")
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(carpeta_out, "mapa_estaciones_NO2_Madrid_2025.png"),
       plot   = mapa_estaciones,
       width  = 13,
       height = 10,
       dpi    = 300,
       bg     = "white")

cat("\n\u2705 Mapa guardado en:\n  ",
    file.path(carpeta_out, "mapa_estaciones_NO2_Madrid_2025.png"), "\n")

print(mapa_estaciones)
