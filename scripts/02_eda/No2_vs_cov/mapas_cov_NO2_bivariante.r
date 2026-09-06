# ==============================================================================
# Mapa BIVARIANTE 3x3: clima (ensemble) x NO2 (ensemble)
# ==============================================================================
# Cada celda de Madrid se colorea según la COMBINACIÓN de su nivel de clima y su
# nivel de NO2 (ambos en terciles: bajo/medio/alto). La paleta 2D permite ver de
# un vistazo las zonas donde coinciden clima alto y NO2 alto (esquina oscura).
#   - Eje "Clima": verde a la derecha (clima alto).
#   - Eje "NO2":   azul hacia arriba (NO2 alto).
#   - Clima alto + NO2 alto = verde-azulado oscuro.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(gstat)
  library(ggplot2)
  library(patchwork)
})
source(here("R", "interpolation", "FUNCIONES_INTERPOLACION.R"))
sf_use_s2(FALSE)

ANIO <- 2025L
K_VECINOS <- 3L
DIR_SALIDA <- here("outputs", "figures", "mapas_cov_NO2")

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
tercil <- function(x) {
  q <- quantile(x, c(1 / 3, 2 / 3), na.rm = TRUE)
  cut(x, c(-Inf, q, Inf), labels = c("1", "2", "3"))
}

# ------------------------------------------------------------------------------
# Datos
# ------------------------------------------------------------------------------

dt_meteo <- as.data.table(readRDS(
  here("data", "processed", "Clima", "diario",
       sprintf("meteo_madrid_%d_diario.rds", ANIO))
))
normalizar_coordenadas_madrid(dt_meteo)
v_pres <- grep("^Presion Bar", names(dt_meteo), value = TRUE)
v_rad <- grep("Solar$", names(dt_meteo), value = TRUE)
VARIABLES <- c("Temperatura", "Humedad_Relativa", "Precipitaciones",
               v_pres, v_rad, "Velocidad Viento")
ETIQUETAS <- c(
  "Temperatura" = "Temperatura", "Humedad_Relativa" = "Humedad relativa",
  "Precipitaciones" = "Precipitación", "Velocidad Viento" = "Viento"
)
ETIQUETAS[v_pres] <- "Presión"
ETIQUETAS[v_rad] <- "Radiación"

med_est <- dt_meteo[, lapply(.SD, mean, na.rm = TRUE),
                    by = .(ESTACION, LONGITUD, LATITUD), .SDcols = VARIABLES]

dt_no2 <- as.data.table(readRDS(
  here("data", "processed", "Contaminacion", "horario",
       sprintf("aire_madrid_%d_No2_horarios.rds", ANIO))
))
no2_est <- dt_no2[, .(NO2 = mean(DATO, na.rm = TRUE)),
                  by = .(ESTACION, LONGITUD, LATITUD)][!is.na(NO2)]

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
bb <- st_bbox(sf_puntos(no2_est))
mx <- 0.10 * (bb$xmax - bb$xmin)
my <- 0.10 * (bb$ymax - bb$ymin)
ZOOM_X <- c(bb$xmin - mx, bb$xmax + mx)
ZOOM_Y <- c(bb$ymin - my, bb$ymax + my)

# ------------------------------------------------------------------------------
# Superficies: NO2 (una vez) y cada covariable
# ------------------------------------------------------------------------------

pesos_ens <- lapply(VARIABLES, function(v) pesos_ensemble_loocv(dt_meteo, v, K_VECINOS))
names(pesos_ens) <- VARIABLES

sf_no2e <- sf_puntos(no2_est[, .(ESTACION, LONGITUD, LATITUD, Z = NO2)])
no2_surf <- superficie_ensemble(sf_no2e, rejilla_madrid, K_VECINOS, c(1, 1, 1) / 3)
no2_cl <- tercil(no2_surf)

# Paleta bivariante 3x3: filas = NO2 (1..3 abajo->arriba), col = clima (1..3).
bipal <- c(
  "1-1" = "#e8e8e8", "2-1" = "#b8d6be", "3-1" = "#73ae80",
  "1-2" = "#b5c0da", "2-2" = "#90b2b3", "3-2" = "#5a9178",
  "1-3" = "#6c83b5", "2-3" = "#567994", "3-3" = "#2a5a5b"
)

panel_bivar <- function(v) {
  me <- med_est[!is.na(get(v)), .(ESTACION, LONGITUD, LATITUD, Z = get(v))]
  cs <- superficie_ensemble(sf_puntos(me), rejilla_madrid, K_VECINOS, pesos_ens[[v]])
  df <- data.frame(
    X = coords_rejilla[, 1], Y = coords_rejilla[, 2],
    biv = paste0(as.character(tercil(cs)), "-", as.character(no2_cl))
  )
  ggplot() +
    geom_tile(data = df, aes(X, Y, fill = biv), width = dx, height = dy) +
    geom_path(data = bordes_coords, aes(X, Y, group = L1),
              color = "grey20", linewidth = 0.2, inherit.aes = FALSE) +
    scale_fill_manual(values = bipal, guide = "none") +
    coord_equal(xlim = ZOOM_X, ylim = ZOOM_Y) +
    labs(title = unname(ETIQUETAS[v])) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
}

paneles <- lapply(VARIABLES, panel_bivar)

# Leyenda 3x3.
leg_df <- data.table(expand.grid(clima = 1:3, no2 = 1:3))
leg_df[, biv := paste0(clima, "-", no2)]
p_leg <- ggplot(leg_df, aes(clima, no2, fill = biv)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = bipal, guide = "none") +
  labs(x = "Clima →", y = "NO₂ →") +
  coord_equal() +
  theme_minimal(base_size = 9) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold", angle = 90)
  )

fig <- (wrap_plots(paneles, ncol = 3) / p_leg) +
  plot_layout(heights = c(9, 1.4)) +
  plot_annotation(
    title = "Mapa bivariante: clima (ensemble) × NO₂ — media anual 2025",
    subtitle = paste(
      "Cada zona se colorea por su combinación de nivel de clima y de NO₂",
      "(terciles). Esquina verde-azulada oscura = clima alto y NO₂ alto."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

ggsave(file.path(DIR_SALIDA, "mapa_bivariante_NO2.png"),
       plot = fig, width = 15, height = 11, dpi = 200, bg = "white")
cat("Mapa bivariante guardado en:\n",
    file.path(DIR_SALIDA, "mapa_bivariante_NO2.png"), "\n")
