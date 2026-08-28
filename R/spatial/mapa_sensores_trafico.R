
# ==============================================================================
# MAPA DE UBICACIÓN DE LOS SENSORES DE TRÁFICO DE MADRID (2025)
# ==============================================================================
# Genera un único mapa con la localización de todos los puntos de medida de
# tráfico activos durante 2025, diferenciados por tipología (URB / M-30 / Otros).
# Sistema de referencia: ETRS89 / UTM zona 30N (EPSG:25830), coherente con el
# resto de mapas del TFM (véase R/spatial/mapas_ubicaciones_clima.R).
# ==============================================================================

library(data.table)
library(sf)
library(ggplot2)
library(ggspatial)
library(here)

# ==============================================================================
# 1. GEOMETRÍAS ADMINISTRATIVAS (DISTRITOS Y BARRIOS)
# ==============================================================================
carpeta_geo <- here("data", "raw", "geometrias")

mapa_distritos <- st_read(
  file.path(carpeta_geo, "DISTRITOS.shp"),
  quiet = TRUE
) |>
  st_transform(25830)

mapa_barrios <- st_read(
  file.path(carpeta_geo, "BARRIOS.shp"),
  quiet = TRUE
) |>
  st_transform(25830)

mapa_distritos$codigo_distrito <- as.character(mapa_distritos$COD_DIS_TX)
mapa_distritos$distrito <- trimws(mapa_distritos$NOMBRE)

orden_distritos <- order(as.integer(mapa_distritos$codigo_distrito))
niveles_distritos <- mapa_distritos$distrito[orden_distritos]
mapa_distritos$distrito <- factor(mapa_distritos$distrito, levels = niveles_distritos)

# Contorno del municipio: se usa como marco del mapa y para validar que los
# sensores caen dentro del término municipal.
limite_municipio <- st_union(st_geometry(mapa_distritos))

# Puntos interiores para etiquetar los distritos con su código.
puntos_distritos <- st_sf(
  codigo_distrito = mapa_distritos$codigo_distrito,
  distrito = mapa_distritos$distrito,
  geometry = st_point_on_surface(st_geometry(mapa_distritos))
)

# ==============================================================================
# 2. INVENTARIO DE SENSORES DE TRÁFICO 2025
# ==============================================================================
# El Ayuntamiento publica un fichero mensual de puntos de medida. Se combinan
# los 12 meses de 2025 y se consolida un registro por sensor: la coordenada es
# la mediana mensual (robusta frente a pequeños reposicionamientos) y se guarda
# el número de meses en que el sensor aparece en el inventario.
carpeta_detectores <- here("data", "raw", "Datos_trafico", "Detectores_2025")

archivos_detectores <- list.files(
  carpeta_detectores,
  pattern = "[.]csv$",
  full.names = TRUE
)

if (length(archivos_detectores) == 0) {
  stop("No se han encontrado ficheros de detectores en: ", carpeta_detectores)
}

detectores <- rbindlist(
  lapply(
    archivos_detectores,
    function(ruta) fread(ruta, sep = ";", dec = ".", encoding = "Latin-1")
  ),
  fill = TRUE
)

columnas_requeridas <- c("tipo_elem", "id", "utm_x", "utm_y")
columnas_faltantes <- setdiff(columnas_requeridas, names(detectores))
if (length(columnas_faltantes) > 0) {
  stop("Faltan columnas en los ficheros de detectores: ",
       paste(columnas_faltantes, collapse = ", "))
}

detectores <- detectores[is.finite(utm_x) & is.finite(utm_y) & utm_x > 0 & utm_y > 0]

sensores <- detectores[, .(
  tipo_elem = tipo_elem[1],
  n_meses   = .N,
  utm_x     = median(utm_x),
  utm_y     = median(utm_y)
), by = id]

# La tipología es estable dentro del año; se verifica antes de usar el primer valor.
tipos_por_sensor <- detectores[, .(n_tipos = uniqueN(tipo_elem)), by = id]
if (any(tipos_por_sensor$n_tipos > 1)) {
  warning("Hay sensores con más de una tipología en 2025: ",
          sum(tipos_por_sensor$n_tipos > 1))
}

# Los rótulos que aparecen en la figura están en inglés (requisito de la memoria).
etiquetas_tipo <- c(
  "URB"   = "Urban (municipal road network)",
  "M30"   = "M-30 ring road and accesses",
  "other" = "Other measurement points"
)

sensores[, Tipologia := factor(
  etiquetas_tipo[tipo_elem],
  levels = unname(etiquetas_tipo)
)]

if (any(is.na(sensores$Tipologia))) {
  stop("Tipologías no contempladas: ",
       paste(unique(sensores[is.na(Tipologia)]$tipo_elem), collapse = ", "))
}

sf_sensores <- st_as_sf(sensores, coords = c("utm_x", "utm_y"), crs = 25830, remove = FALSE)

# ==============================================================================
# 3. CONTROL DE CALIDAD ESPACIAL
# ==============================================================================
dentro_municipio <- lengths(st_intersects(sf_sensores, limite_municipio)) > 0
sf_sensores$dentro_municipio <- dentro_municipio

if (any(!dentro_municipio)) {
  warning(sum(!dentro_municipio),
          " sensor(es) fuera del término municipal; se excluyen del mapa")
}

sf_sensores <- sf_sensores[dentro_municipio, ]

# Asignación de cada sensor a su distrito (para la tabla resumen).
sf_sensores <- st_join(
  sf_sensores,
  mapa_distritos["distrito"],
  join = st_within,
  left = TRUE
)

cat("Sensores de tráfico 2025 cartografiados:", nrow(sf_sensores), "\n")
print(table(sf_sensores$Tipologia))

# ==============================================================================
# 4. LEYENDA CON RECUENTOS Y PALETA
# ==============================================================================
recuento_tipo <- sf_sensores |>
  st_drop_geometry() |>
  as.data.table() |>
  (\(dt) dt[, .N, by = Tipologia])()

etiquetas_leyenda <- setNames(
  sprintf("%s  (n = %s)",
          recuento_tipo$Tipologia,
          format(recuento_tipo$N, big.mark = ",", trim = TRUE)),
  as.character(recuento_tipo$Tipologia)
)

# Paleta segura para daltonismo (Okabe-Ito).
colores_tipo <- setNames(
  c("#0072B2", "#D55E00", "#009E73"),
  unname(etiquetas_tipo)
)

formas_tipo <- setNames(c(21, 24, 22), unname(etiquetas_tipo))

# Los sensores urbanos son mayoría y se dibujan primero; M-30 y "otros" se
# superponen después para que no queden ocultos por la nube de puntos urbanos.
sf_urbanos <- sf_sensores[sf_sensores$Tipologia == etiquetas_tipo[["URB"]], ]
sf_resto   <- sf_sensores[sf_sensores$Tipologia != etiquetas_tipo[["URB"]], ]

# ==============================================================================
# 5. MAPA
# ==============================================================================
bb <- st_bbox(limite_municipio)
margen <- 900  # metros de aire alrededor del municipio

mapa_sensores <- ggplot() +
  geom_sf(data = mapa_barrios, fill = "#F2F3F5", color = "white", linewidth = 0.18) +
  geom_sf(data = mapa_distritos, fill = NA, color = "gray55", linewidth = 0.35) +
  geom_sf(data = limite_municipio, fill = NA, color = "gray20", linewidth = 0.7) +
  geom_sf(
    data = sf_urbanos,
    aes(fill = Tipologia, shape = Tipologia),
    color = "white", size = 1.5, stroke = 0.12, alpha = 0.85
  ) +
  geom_sf(
    data = sf_resto,
    aes(fill = Tipologia, shape = Tipologia),
    color = "white", size = 2.1, stroke = 0.25, alpha = 0.95
  ) +
  geom_sf_text(
    data = puntos_distritos,
    aes(label = codigo_distrito),
    size = 2.5, fontface = "bold", color = "gray15"
  ) +
  scale_fill_manual(
    values = colores_tipo,
    labels = etiquetas_leyenda,
    drop = FALSE,
    name = "Sensor type"
  ) +
  scale_shape_manual(
    values = formas_tipo,
    labels = etiquetas_leyenda,
    drop = FALSE,
    name = "Sensor type"
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.22, height = unit(0.18, "cm"),
    text_cex = 0.65, line_width = 0.5, pad_y = unit(0.35, "cm")
  ) +
  annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(0.95, "cm"), width = unit(0.75, "cm"),
    style = north_arrow_minimal(line_width = 0.8, text_size = 7)
  ) +
  coord_sf(
    crs = 25830,
    xlim = c(bb["xmin"] - margen, bb["xmax"] + margen),
    ylim = c(bb["ymin"] - margen, bb["ymax"] + margen),
    datum = 25830,
    expand = FALSE
  ) +
  scale_x_continuous(labels = function(x) sprintf("%.0f", x / 1000)) +
  scale_y_continuous(labels = function(y) sprintf("%.0f", y / 1000)) +
  labs(
    title = "Traffic sensor network of the city of Madrid (2025)",
    subtitle = sprintf(
      "%s active measurement points, classified by road type",
      format(nrow(sf_sensores), big.mark = ",", trim = TRUE)
    ),
    x = "UTM Easting (km)",
    y = "UTM Northing (km)",
    caption = paste0(
      "Numbers on the map identify the 21 districts of Madrid. ",
      "Coordinates: median of the 2025 monthly sensor inventories.\n",
      "Source: Madrid City Council Open Data Portal; own elaboration. ",
      "Coordinate reference system: ETRS89 / UTM zone 30N (EPSG:25830)."
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3.4, alpha = 1, stroke = 0.3),
      nrow = 1
    ),
    shape = guide_legend(
      override.aes = list(size = 3.4, alpha = 1, stroke = 0.3),
      nrow = 1
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.border      = element_rect(fill = NA, color = "gray30", linewidth = 0.4),
    panel.grid.major  = element_line(color = "gray92", linewidth = 0.2),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 7.5, color = "gray35"),
    axis.title        = element_text(size = 8.5, color = "gray25"),
    plot.title        = element_text(face = "bold", size = 14.5, color = "gray10"),
    plot.subtitle     = element_text(size = 10, color = "gray35", margin = margin(b = 9)),
    plot.caption      = element_text(size = 7, color = "gray40", hjust = 0,
                                     margin = margin(t = 9)),
    plot.margin       = margin(12, 14, 10, 14),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8.5),
    legend.key.spacing.x = unit(0.35, "cm"),
    legend.box.margin = margin(t = 2)
  )

# ==============================================================================
# 6. EXPORTACIÓN
# ==============================================================================
carpeta_graficos <- here("outputs", "figures", "Mapas")
if (!dir.exists(carpeta_graficos)) dir.create(carpeta_graficos, recursive = TRUE)

ruta_png <- file.path(carpeta_graficos, "mapa_sensores_trafico_2025.png")
ruta_pdf <- file.path(carpeta_graficos, "mapa_sensores_trafico_2025.pdf")

ggsave(ruta_png, plot = mapa_sensores, width = 9, height = 10, dpi = 320, bg = "white")
ggsave(ruta_pdf, plot = mapa_sensores, width = 9, height = 10,
       device = grDevices::cairo_pdf, bg = "white")

# Tabla reproducible con el inventario cartografiado.
carpeta_tablas <- here("outputs", "tables")
if (!dir.exists(carpeta_tablas)) dir.create(carpeta_tablas, recursive = TRUE)

tabla_sensores <- as.data.table(st_drop_geometry(sf_sensores))[, .(
  id,
  tipo_elem,
  Tipologia = as.character(Tipologia),
  distrito = as.character(distrito),
  n_meses,
  utm_x,
  utm_y
)][order(id)]

fwrite(tabla_sensores, file.path(carpeta_tablas, "sensores_trafico_2025.csv"))

cat("Mapa guardado:", ruta_png, "\n")
cat("Versión vectorial:", ruta_pdf, "\n")
cat("Inventario guardado:", file.path(carpeta_tablas, "sensores_trafico_2025.csv"), "\n")
