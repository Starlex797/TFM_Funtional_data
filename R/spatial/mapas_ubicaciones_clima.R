
# ==============================================================================
# GENERACIÓN DE Mapa METEOROLOGÍA 
# ==============================================================================

library(data.table)
library(sf)
library(ggplot2)
library(ggrepel)
library(here)

source(here("R", "utilities", "dictionaries.R"))      
source(here("R", "cleaning", "cleaning_functions.R"))

# ==============================================================================
# 1. CARGA Y PROCESAMIENTO ESPACIAL DE DISTRITOS Y BARRIOS
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
mapa_barrios$barrio     <- trimws(mapa_barrios$NOMBRE)

orden_distritos <- order(as.integer(mapa_distritos$codigo_distrito))
niveles_distritos <- mapa_distritos$distrito[orden_distritos]
mapa_distritos$distrito <- factor(
  mapa_distritos$distrito,
  levels = niveles_distritos
)

# Los puntos interiores se calculan en UTM para evitar operaciones métricas
# sobre coordenadas geográficas y garantizar que el punto quede en el barrio.
puntos_barrios <- st_sf(
  barrio = mapa_barrios$barrio,
  geometry = st_point_on_surface(st_geometry(mapa_barrios))
)
mapa_barrios$distrito <- st_join(
  puntos_barrios,
  mapa_distritos["distrito"],
  left = TRUE
)$distrito

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
  st_as_sf(coords = c("X_utm", "Y_utm"), crs = 25830)

# Contaminación
sf_aire <- fread(here("data", "raw", "Datos_contaminacion", "Estaciones", "datos.csv"))[
  CODIGO_CORTO %in% as.integer(names(nombres_estaciones_aire))
][, Nombre_Estacion := nombres_estaciones_aire[as.character(CODIGO_CORTO)]] |> 
  st_as_sf(coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)

# Datos horarios empleados para cuantificar la cobertura de cada variable.
ruta_meteo <- here(
  "data", "processed", "Clima", "horario",
  "meteo_madrid_2025_horario.rds"
)

if (!file.exists(ruta_meteo)) {
  stop("No se encuentra el archivo meteorológico: ", ruta_meteo)
}

dt_meteo <- as.data.table(readRDS(ruta_meteo))

vars_clima_cols <- c(
  "Temperatura",
  "Humedad_Relativa",
  "Precipitaciones",
  "Presion Barométrica",
  "Radiación Solar",
  "Velocidad Viento"
)

columnas_faltantes <- setdiff(
  c("ESTACION", "X_km", "Y_km", vars_clima_cols),
  names(dt_meteo)
)

if (length(columnas_faltantes) > 0) {
  stop(
    "Faltan columnas meteorológicas: ",
    paste(columnas_faltantes, collapse = ", ")
  )
}

cat("✅ Datos listos: ", nrow(sf_aire), " Contaminación | ", nrow(sf_metereo), " Clima\n")

carpeta_graficos <- here("outputs", "figures")
if (!dir.exists(carpeta_graficos)) dir.create(carpeta_graficos, recursive = TRUE)

# Los mapas 1–3 se conservan en el script, pero dejan de generarse por defecto.
# Cambiar a TRUE únicamente si se desean volver a exportar.
generar_mapas_red <- FALSE

if (generar_mapas_red) {

# ==============================================================================
# 3. CONSTRUCCIÓN DE LA PLANTILLA BASE DEL MAPA
# ==============================================================================
mapa_base <- ggplot() +
  geom_sf(data = mapa_barrios, fill = "#f4f4f4", color = "white", linewidth = 0.15) +
  geom_sf(data = mapa_distritos, fill = NA, color = "gray25", linewidth = 0.65) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )


# ==============================================================================
# 8. MAPA 4: COBERTURA POR VARIABLE CLIMATOLÓGICA (FACETADO)
# ==============================================================================
# Los rótulos que aparecen en la figura están en inglés (requisito de la memoria).
etiquetas_vars <- c(
  "Temperatura"          = "Temperature",
  "Humedad_Relativa"    = "Relative Humidity",
  "Precipitaciones"     = "Precipitation",
  "Presion Barométrica" = "Barometric Pressure",
  "Radiación Solar"     = "Solar Radiation",
  "Velocidad Viento"    = "Wind Speed"
)

# Porcentaje de observaciones disponibles por estación y variable.
cobertura <- dt_meteo[, lapply(
  .SD,
  function(x) round(mean(!is.na(x)) * 100, 1)
), by = ESTACION, .SDcols = vars_clima_cols]

# Coordenadas UTM (X_km/Y_km están en km → convertir a metros)
coords_estaciones <- dt_meteo[, .(X_m = mean(X_km, na.rm = TRUE) * 1000,
                                   Y_m = mean(Y_km, na.rm = TRUE) * 1000),
                               by = ESTACION]
cobertura <- merge(cobertura, coords_estaciones, by = "ESTACION")

if (any(!is.finite(cobertura$X_m)) || any(!is.finite(cobertura$Y_m))) {
  stop("Existen estaciones sin coordenadas UTM válidas")
}

# Formato largo: una fila por estación y variable.
dt_long_vars <- melt(cobertura,
                     id.vars       = c("ESTACION", "X_m", "Y_m"),
                     measure.vars  = vars_clima_cols,
                     variable.name = "Variable",
                     value.name    = "pct_disponible")

dt_long_vars[, Variable_label := factor(
  etiquetas_vars[as.character(Variable)],
  levels = c(
    "Temperature",
    "Relative Humidity",
    "Barometric Pressure",
    "Precipitation",
    "Solar Radiation",
    "Wind Speed"
  )
)]

dt_long_vars[, Categoria := factor(
  fcase(
    pct_disponible == 0, "No data",
    pct_disponible < 80, "<80%",
    pct_disponible < 95, "80–<95%",
    default = "≥95%"
  ),
  levels = c("No data", "<80%", "80–<95%", "≥95%")
)]

# Objeto espacial y asignación de cada estación a su distrito.
sf_long_vars <- st_as_sf(
  copy(dt_long_vars),
  coords = c("X_m", "Y_m"),
  crs = 25830,
  remove = FALSE
)

sf_long_vars <- st_join(
  sf_long_vars,
  mapa_distritos["distrito"],
  join = st_within,
  left = TRUE
)

if (nrow(sf_long_vars) != nrow(dt_long_vars)) {
  stop("La unión espacial con distritos ha duplicado estaciones")
}

if (any(is.na(sf_long_vars$distrito))) {
  warning("Hay estaciones que no se han podido asociar a un distrito")
}

# Puntos interiores empleados exclusivamente para etiquetar los distritos.
puntos_distritos <- st_sf(
  distrito = mapa_distritos$distrito,
  codigo_distrito = mapa_distritos$codigo_distrito,
  geometry = st_point_on_surface(st_geometry(mapa_distritos))
)

# Tabla reproducible de cobertura.
carpeta_tablas <- here("outputs", "tables")
if (!dir.exists(carpeta_tablas)) dir.create(carpeta_tablas, recursive = TRUE)

tabla_cobertura <- as.data.table(st_drop_geometry(sf_long_vars))[, .(
  ESTACION,
  distrito,
  Variable = as.character(Variable_label),
  pct_disponible,
  Categoria = as.character(Categoria),
  X_m,
  Y_m
)]

fwrite(
  tabla_cobertura,
  file.path(carpeta_tablas, "cobertura_estaciones_clima.csv")
)

resumen_cobertura <- tabla_cobertura[, .(
  estaciones_con_datos = sum(pct_disponible > 0),
  estaciones_sin_datos = sum(pct_disponible == 0),
  disponibilidad_mediana = median(pct_disponible)
), by = Variable]
print(resumen_cobertura)

colores_variables <- c(
  "Temperature"          = "#D55E00",
  "Relative Humidity"    = "#0072B2",
  "Precipitation"        = "#009E73",
  "Barometric Pressure"  = "#8C564B",
  "Solar Radiation"      = "#E69F00",
  "Wind Speed"           = "#6A3D9A"
)

colores_distritos <- setNames(
  grDevices::hcl.colors(
    n = nrow(mapa_distritos),
    palette = "Set 3"
  ),
  niveles_distritos
)

etiquetas_leyenda_distritos <- setNames(
  paste(
    mapa_distritos$codigo_distrito[orden_distritos],
    niveles_distritos,
    sep = " · "
  ),
  niveles_distritos
)

mapa_4_variables <- ggplot() +
  geom_sf(
    data = mapa_distritos,
    aes(fill = distrito),
    color = "white",
    linewidth = 0.35,
    alpha = 0.78
  ) +
  geom_sf(
    data = mapa_barrios,
    fill = NA,
    color = "white",
    linewidth = 0.08,
    alpha = 0.55
  ) +
  geom_sf(
    data = mapa_distritos,
    fill = NA,
    color = "gray25",
    linewidth = 0.5
  ) +
  geom_sf(
    data = sf_long_vars[sf_long_vars$pct_disponible == 0, ],
    color = "gray25",
    size = 2.4,
    shape = 4,
    stroke = 0.75
  ) +
  geom_sf(
    data = sf_long_vars[sf_long_vars$pct_disponible > 0, ],
    color = "white",
    size = 4,
    shape = 16
  ) +
  geom_sf(
    data = sf_long_vars[sf_long_vars$pct_disponible > 0, ],
    aes(color = Variable_label),
    size = 2.9,
    shape = 16
  ) +
  geom_sf_label(
    data = puntos_distritos,
    aes(label = codigo_distrito),
    color = "gray15",
    fill = "white",
    size = 1.7,
    fontface = "bold",
    linewidth = 0.12,
    label.padding = grid::unit(0.08, "lines"),
    label.r = grid::unit(0.04, "lines")
  ) +
  scale_fill_manual(
    values = colores_distritos,
    breaks = niveles_distritos,
    labels = etiquetas_leyenda_distritos,
    drop = FALSE,
    name = "District"
  ) +
  scale_color_manual(
    values = colores_variables,
    guide = "none"
  ) +
  facet_wrap(~ Variable_label, ncol = 3) +
  coord_sf(datum = NA, expand = FALSE, clip = "off") +
  labs(
    title = "Spatial coverage of the meteorological network by variable",
    subtitle = paste0(
      "Stations with valid hourly observations in 2025. ",
      "Colored dot: data available; grey cross: no data."
    ),
    caption = paste0(
      "Numbers on the map identify the 21 districts of Madrid. ",
      "Source: Madrid City Council Open Data Portal; own elaboration. ",
      "Coordinate reference system: ETRS89 / UTM zone 30N (EPSG:25830)."
    ),
    x = NULL,
    y = NULL
  ) +
  guides(
    fill = guide_legend(
      title = "District",
      ncol = 3,
      byrow = TRUE,
      override.aes = list(alpha = 0.9, color = "white")
    )
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.spacing = grid::unit(8, "pt"),
    strip.text = element_text(face = "bold", size = 9.5, color = "gray15"),
    strip.background = element_rect(fill = "#EEF1F4", color = NA),
    plot.title = element_text(face = "bold", size = 16, color = "gray10"),
    plot.subtitle = element_text(size = 9.5, color = "gray35", margin = margin(b = 8)),
    plot.caption = element_text(size = 7.5, color = "gray40", hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(12, 14, 10, 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 7.5),
    legend.key.height = grid::unit(0.42, "cm"),
    legend.key.width = grid::unit(0.42, "cm"),
    legend.spacing.x = grid::unit(0.12, "cm"),
    legend.box.margin = margin(t = 4)
  )

ggsave(file.path(carpeta_graficos, "mapa_04_variables_climatologicas.png"),
       plot = mapa_4_variables, width = 16, height = 11, dpi = 300, bg = "white")
ggsave(file.path(carpeta_graficos, "mapa_04_variables_climatologicas.pdf"),
       plot = mapa_4_variables, width = 16, height = 11,
       device = grDevices::cairo_pdf, bg = "white")

cat("✅ Mapa 4 guardado:", file.path(carpeta_graficos, "mapa_04_variables_climatologicas.png"), "\n")
cat("✅ Versión vectorial guardada:", file.path(carpeta_graficos, "mapa_04_variables_climatologicas.pdf"), "\n")
cat("✅ Tabla de cobertura guardada:", file.path(carpeta_tablas, "cobertura_estaciones_clima.csv"), "\n")

}
