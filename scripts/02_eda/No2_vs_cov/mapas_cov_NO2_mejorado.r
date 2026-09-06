# ==============================================================================
# Mapa mejorado: anomalía climática (ensemble) + NO2 dentro
# ==============================================================================
# Versión cartográficamente más legible del mapa covariables + NO2:
#   - La superficie se muestra como ANOMALÍA respecto a la media de Madrid, con
#     escala DIVERGENTE (azul = por debajo, rojo = por encima) -> los focos
#     cálidos/secos/etc. resaltan mucho más que con una escala secuencial.
#   - El NO2 se dibuja como BURBUJAS con halo, coloreadas en verde-amarillo
#     (viridis) para contrastar con el fondo rojo/azul.
#   - Las estaciones METEO se marcan con triángulos: dónde hay dato real frente
#     a interpolación (no sobre-interpretar el NW despoblado).
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(gstat)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})
source(here("R", "interpolation", "FUNCIONES_INTERPOLACION.R"))
sf_use_s2(FALSE)

ANIO <- 2025L
K_VECINOS <- 3L
DIR_SALIDA <- here("outputs", "figures", "mapas_cov_NO2")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Utilidades
# ------------------------------------------------------------------------------

normalizar_coordenadas_madrid <- function(dt) {
  filas_lat <- !is.na(dt$LATITUD) & dt$LATITUD < 35
  filas_lon <- !is.na(dt$LONGITUD) & dt$LONGITUD > -1
  dt[filas_lat, LATITUD := LATITUD * 10^round(log10(40.4 / abs(LATITUD)))]
  dt[filas_lon, LONGITUD := LONGITUD * 10]
  invisible(dt)
}

sf_puntos <- function(dt) {
  st_transform(st_as_sf(dt, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830)
}

superficie_ensemble <- function(sf_est, sf_target, k, pesos) {
  f <- Z ~ 1
  k_real <- min(k, nrow(sf_est))
  p_nn <- idw(f, sf_est, sf_target, nmax = 1, debug.level = 0)$var1.pred
  p_idw <- idw(f, sf_est, sf_target, nmax = k_real, idp = 1, debug.level = 0)$var1.pred
  p_knn <- idw(f, sf_est, sf_target, nmax = k_real, idp = 0, debug.level = 0)$var1.pred
  pesos[1] * p_nn + pesos[2] * p_idw + pesos[3] * p_knn
}

# ------------------------------------------------------------------------------
# 2. Datos
# ------------------------------------------------------------------------------

dt_meteo <- as.data.table(readRDS(
  here("data", "processed", "Clima", "diario",
       sprintf("meteo_madrid_%d_diario.rds", ANIO))
))
normalizar_coordenadas_madrid(dt_meteo)

v_pres <- grep("^Presion Bar", names(dt_meteo), value = TRUE)
v_rad <- grep("Solar$", names(dt_meteo), value = TRUE)
VARIABLES <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  v_pres, v_rad, "Velocidad Viento"
)
ETIQUETAS <- c(
  "Temperatura" = "Temperatura (°C)",
  "Humedad_Relativa" = "Humedad relativa (%)",
  "Precipitaciones" = "Precipitación (mm/día)",
  "Velocidad Viento" = "Viento (m/s)"
)
ETIQUETAS[v_pres] <- "Presión (hPa)"
ETIQUETAS[v_rad] <- "Radiación (W/m²)"

med_est <- dt_meteo[
  , lapply(.SD, mean, na.rm = TRUE),
  by = .(ESTACION, LONGITUD, LATITUD), .SDcols = VARIABLES
]
sf_meteo <- sf_puntos(unique(med_est[, .(ESTACION, LONGITUD, LATITUD)]))
df_meteo <- as.data.frame(st_coordinates(sf_meteo))

dt_no2 <- as.data.table(readRDS(
  here("data", "processed", "Contaminacion", "horario",
       sprintf("aire_madrid_%d_No2_horarios.rds", ANIO))
))
no2_est <- dt_no2[
  , .(NO2 = mean(DATO, na.rm = TRUE)),
  by = .(ESTACION, LONGITUD, LATITUD)
][!is.na(NO2)]
xy_no2 <- st_coordinates(sf_puntos(no2_est))
df_no2 <- data.frame(X = xy_no2[, 1], Y = xy_no2[, 2], NO2 = no2_est$NO2)

# ------------------------------------------------------------------------------
# 3. Geometría
# ------------------------------------------------------------------------------

mapa_utm <- st_transform(st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE
)), 25830)
rejilla <- st_as_sf(st_make_grid(mapa_utm, n = c(120, 120), what = "centers"))
rejilla_madrid <- st_intersection(rejilla, st_union(mapa_utm))
coords_rejilla <- st_coordinates(rejilla_madrid)
dx <- diff(sort(unique(coords_rejilla[, 1])))[1]
dy <- diff(sort(unique(coords_rejilla[, 2])))[1]
bordes_coords <- as.data.frame(st_coordinates(
  st_cast(st_cast(st_geometry(mapa_utm), "MULTILINESTRING"), "LINESTRING")
))

# ------------------------------------------------------------------------------
# 4. Pesos ensemble
# ------------------------------------------------------------------------------

pesos_ens <- lapply(VARIABLES, function(v) pesos_ensemble_loocv(dt_meteo, v, K_VECINOS))
names(pesos_ens) <- VARIABLES

# ------------------------------------------------------------------------------
# 5. Paneles de anomalía
# ------------------------------------------------------------------------------

lim_no2 <- range(df_no2$NO2, na.rm = TRUE)

# Zoom a la almendra central: bounding box de las estaciones de NO2 + margen.
bb <- st_bbox(sf_puntos(no2_est))
mx <- 0.10 * (bb$xmax - bb$xmin)
my <- 0.10 * (bb$ymax - bb$ymin)
ZOOM_X <- c(bb$xmin - mx, bb$xmax + mx)
ZOOM_Y <- c(bb$ymin - my, bb$ymax + my)

panel_anomalia <- function(v, mostrar_no2 = FALSE) {
  me <- med_est[!is.na(get(v)), .(ESTACION, LONGITUD, LATITUD, Z = get(v))]
  surf <- superficie_ensemble(sf_puntos(me), rejilla_madrid, K_VECINOS, pesos_ens[[v]])
  media_madrid <- mean(surf, na.rm = TRUE)
  anom <- surf - media_madrid
  lim <- max(abs(anom), na.rm = TRUE)
  df_s <- data.frame(X = coords_rejilla[, 1], Y = coords_rejilla[, 2], Anom = anom)

  ggplot() +
    geom_tile(data = df_s, aes(X, Y, fill = Anom), width = dx, height = dy) +
    geom_path(data = bordes_coords, aes(X, Y, group = L1),
              color = "grey30", linewidth = 0.2, inherit.aes = FALSE) +
    # Estaciones meteo (anclas de la interpolación).
    geom_point(data = df_meteo, aes(X, Y), shape = 2,
               color = "grey25", size = 1.2, stroke = 0.5) +
    # NO2: halo + burbuja coloreada.
    geom_point(data = df_no2, aes(X, Y), color = "black", size = 3.3) +
    geom_point(data = df_no2, aes(X, Y, color = NO2), size = 2.6) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "grey96", high = "#b2182b", midpoint = 0,
      limits = c(-lim, lim),
      name = paste0("Anomalía\n", unname(ETIQUETAS[v]))
    ) +
    scale_color_viridis_c(
      option = "viridis", limits = lim_no2, name = "NO₂ medio\n(µg/m³)",
      guide = if (mostrar_no2) "colourbar" else "none"
    ) +
    coord_equal(xlim = ZOOM_X, ylim = ZOOM_Y) +
    labs(title = unname(ETIQUETAS[v])) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.32, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
}

paneles <- lapply(seq_along(VARIABLES), function(i) {
  panel_anomalia(VARIABLES[i], mostrar_no2 = (i == length(VARIABLES)))
})

fig <- wrap_plots(paneles, ncol = 3) +
  plot_annotation(
    title = "Anomalía climática (ensemble) y NO₂ — media anual 2025",
    subtitle = paste(
      "Relleno = desviación de cada zona respecto a la media de Madrid",
      "(rojo = por encima, azul = por debajo).",
      "\nBurbujas = NO₂ medio por estación. Triángulos = estaciones meteo (dato real)."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

ggsave(
  file.path(DIR_SALIDA, "mapa_anomalias_NO2.png"),
  plot = fig, width = 16, height = 9.5, dpi = 200, bg = "white"
)

cat("Mapa mejorado guardado en:\n", file.path(DIR_SALIDA, "mapa_anomalias_NO2.png"), "\n")
