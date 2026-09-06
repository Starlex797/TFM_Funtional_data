# ==============================================================================
# MAPA DE ESTACIONES DE MONITOREO DE NO₂ — Madrid 2025
# Por distritos · Diferenciadas por tipología · Con etiquetas de estación
# ==============================================================================
# Aporta el contexto espacia del estudio. Permite conocer dónde se encuentran las estaciones,
# qué tipo de entorno representan y cómo se distribuyen territorialmente los niveles medios de No2.
# ==============================================================================
# Resumen del Script:
# El script genera dos tipos de mapas.
# Mapa de localización de estacionesLee la cartografía de los distritos de Madrid.
# Extrae las coordenadas y la tipología de cada estación:urbana de tráfico;
# urbana de fondo;
# suburbana.

# Convierte las estaciones a objetos espaciales y las reproyecta al sistema de referencia del mapa.
# Representa los distritos y coloca las estaciones con colores, símbolos y etiquetas diferentes según su tipología.

# Mapas de contaminación anualDefine una función aplicable a distintos años.
# Calcula el NO₂ medio anual de cada estación.
# Colorea y dimensiona las estaciones según su concentración media.
# Destaca las estaciones con una media igual o superior a 40 µg/m³.
# Utiliza una escala fija de 0 a 60 para que los mapas de distintos años sean comparables.
# Genera los mapas correspondientes a 2025 y 2019.
# =============================================================================================




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
  quiet = TRUE
)

# Estaciones: coordenadas + tipología (una fila por estación)
dt_mens <- readRDS(here(
  "data", "processed", "contaminacion", "mensual",
  "aire_madrid_2025_log_No2_mensuales1.rds"
))

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
  crs    = 4326
) |>
  st_transform(st_crs(distritos))

# Coordenadas XY explícitas para ggrepel
sf_estaciones$X <- st_coordinates(sf_estaciones)[, 1]
sf_estaciones$Y <- st_coordinates(sf_estaciones)[, 2]

# Orden fijo de tipologías en leyenda
sf_estaciones$NOM_TIPO <- factor(sf_estaciones$NOM_TIPO,
  levels = c(
    "Urbana tr\u00e1fico",
    "Urbana fondo",
    "Suburbana"
  )
)

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
  "Urbana tr\u00e1fico" = 16, # círculo sólido
  "Urbana fondo"        = 17, # triángulo sólido
  "Suburbana"           = 15 # cuadrado sólido
)


# ------------------------------------------------------------------------------
# 3. MAPA
# ------------------------------------------------------------------------------

mapa_estaciones <- ggplot() +

  # — Polígonos de distritos coloreados —
  geom_sf(
    data = distritos,
    aes(fill = NOMBRE),
    colour = "white",
    linewidth = 0.55
  ) +
  scale_fill_manual(
    values = colores_distritos,
    name = "Distritos de Madrid",
    guide = guide_legend(
      ncol = 2,
      keywidth = unit(0.45, "cm"),
      keyheight = unit(0.38, "cm"),
      override.aes = list(
        colour = "grey70",
        linewidth = 0.3
      )
    )
  ) +

  # — Puntos de estaciones (colour + shape, sin fill) —
  geom_point(
    data = sf_estaciones,
    aes(
      x = X,
      y = Y,
      colour = NOM_TIPO,
      shape = NOM_TIPO
    ),
    size = 3.8,
    stroke = 0.5
  ) +
  scale_colour_manual(
    values = colores_estacion,
    name = "Tipolog\u00eda de estaci\u00f3n",
    guide = guide_legend(
      keywidth = unit(0.5, "cm"),
      keyheight = unit(0.5, "cm"),
      override.aes = list(size = 4)
    )
  ) +
  scale_shape_manual(
    values = formas_estacion,
    name   = "Tipolog\u00eda de estaci\u00f3n"
  ) +

  # — Etiquetas de estaciones con fondo blanco (ggrepel) —
  geom_text_repel(
    data = sf_estaciones,
    aes(x = X, y = Y, label = ESTACION),
    size = 2.3,
    colour = "grey10",
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.12,
    box.padding = 0.45,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size = 0.35,
    max.overlaps = 40,
    seed = 7231
  ) +
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

carpeta_out <- here("outputs", "figures", "EDA", "NO2", "mapa_estaciones_NO2")
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(carpeta_out, "mapa_estaciones_NO2_Madrid_2025.png"),
  plot   = mapa_estaciones,
  width  = 13,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

cat(
  "\n\u2705 Mapa guardado en:\n  ",
  file.path(carpeta_out, "mapa_estaciones_NO2_Madrid_2025.png"), "\n"
)

print(mapa_estaciones)


# ==============================================================================
# BLOQUE 5: MAPA POR NIVEL DE CONTAMINACI\u00d3N (parametrizado por a\u00f1o)
# Estaciones coloreadas por NO\u2082 medio anual \u00b7 Umbral regulatorio 40 \u00b5g/m\u00b3
# ==============================================================================

# Umbral regulatorio (valor l\u00edmite anual UE \u2014 Directiva 2008/50/CE)
UMBRAL_NO2 <- 40

# Escala de color FIJA para todos los a\u00f1os \u2192 mapas comparables entre s\u00ed
LIMITES_NO2 <- c(0, 60)

mapa_contaminacion_anual <- function(anio) {
  # -- NO\u2082 medio anual por estaci\u00f3n (datos crudos diarios, \u00b5g/m\u00b3) --
  # El archivo diario ya trae LONGITUD/LATITUD, no hace falta cruzar con nada.
  dt_diario <- readRDS(here(
    "data", "processed", "contaminacion", "diario",
    sprintf("aire_madrid_%d_No2_trans_diarios.rds", anio)
  ))

  media_est <- dt_diario[!is.na(DATO_DIARIO), .(
    NO2_medio = mean(DATO_DIARIO),
    LONGITUD  = LONGITUD[1],
    LATITUD   = LATITUD[1]
  ), by = ESTACION]

  media_est[, GRAVE := NO2_medio >= UMBRAL_NO2]

  # sf (WGS84) \u2192 mismo CRS que los distritos
  sf_conta <- st_as_sf(media_est,
    coords = c("LONGITUD", "LATITUD"),
    crs    = 4326
  ) |>
    st_transform(st_crs(distritos))
  sf_conta$X <- st_coordinates(sf_conta)[, 1]
  sf_conta$Y <- st_coordinates(sf_conta)[, 2]
  sf_conta$ETIQUETA <- sprintf("%s\n%.1f", sf_conta$ESTACION, sf_conta$NO2_medio)

  n_graves <- sum(sf_conta$GRAVE)
  cat(sprintf(
    "\n[%d] Estaciones \u2265 %d \u00b5g/m\u00b3 (grave): %d de %d\n",
    anio, UMBRAL_NO2, n_graves, nrow(sf_conta)
  ))

  # -- MAPA --
  p <- ggplot() +

    # Distritos en gris neutro: no compiten con el color de NO\u2082
    geom_sf(
      data = distritos,
      fill = "grey96", colour = "grey80", linewidth = 0.4
    ) +

    # Halo rojo bajo las estaciones graves (\u2265 40) para destacarlas
    geom_point(
      data = subset(sf_conta, GRAVE),
      aes(x = X, y = Y),
      colour = "#7b241c", size = 8.5, shape = 21,
      fill = NA, stroke = 1.3
    ) +

    # Estaciones: color = NO\u2082 medio (divergente centrado en 40), tama\u00f1o = NO\u2082
    geom_point(
      data = sf_conta,
      aes(x = X, y = Y, fill = NO2_medio, size = NO2_medio),
      shape = 21, colour = "grey25", stroke = 0.6
    ) +
    scale_fill_gradient2(
      low = "#1a5276", mid = "#f0f0f0", high = "#7b241c",
      midpoint = UMBRAL_NO2,
      limits = LIMITES_NO2,
      oob = scales::squish,
      name = "NO\u2082 medio\n(\u00b5g/m\u00b3)",
      guide = guide_colourbar(
        barheight = unit(3.5, "cm"),
        barwidth = unit(0.45, "cm"), order = 1
      )
    ) +
    scale_size_continuous(range = c(3, 9), limits = LIMITES_NO2, guide = "none") +

    # Etiqueta con nombre + valor (los graves en rojo y negrita)
    geom_text_repel(
      data = sf_conta,
      aes(
        x = X, y = Y, label = ETIQUETA,
        colour = GRAVE, fontface = ifelse(GRAVE, "bold", "plain")
      ),
      size = 2.4, bg.color = "white", bg.r = 0.12,
      box.padding = 0.5, point.padding = 0.4,
      segment.color = "grey45", segment.size = 0.35,
      max.overlaps = 40, seed = 7231
    ) +
    scale_colour_manual(
      values = c(`TRUE` = "#7b241c", `FALSE` = "grey25"),
      guide = "none"
    ) +
    labs(
      title = sprintf("Nivel de contaminaci\u00f3n por estaci\u00f3n \u2014 NO\u2082 medio %d", anio),
      subtitle = sprintf(
        "%d de %d estaciones superan el valor l\u00edmite anual de 40 \u00b5g/m\u00b3 (rojo)",
        n_graves, nrow(sf_conta)
      ),
      caption = "Fuente: Ayuntamiento de Madrid \u00b7 Umbral: Directiva 2008/50/CE (valor l\u00edmite anual NO\u2082)"
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
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    coord_sf(crs = st_crs(distritos))

  # -- Guardar --
  archivo <- file.path(
    carpeta_out,
    sprintf("mapa_contaminacion_NO2_Madrid_%d.png", anio)
  )
  ggsave(archivo, plot = p, width = 13, height = 10, dpi = 300, bg = "white")
  cat("\u2705 Mapa de contaminaci\u00f3n guardado en:\n  ", archivo, "\n")

  p
}

# Generar los mapas (a\u00f1ade aqu\u00ed m\u00e1s a\u00f1os si lo necesitas)
mapa_conta_2025 <- mapa_contaminacion_anual(2025)
mapa_conta_2019 <- mapa_contaminacion_anual(2019)

print(mapa_conta_2019)
