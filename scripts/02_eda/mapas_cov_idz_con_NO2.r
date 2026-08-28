# ==============================================================================
# Mapas exploratorios: covariables climáticas (ensemble) + NO2 (respuesta)
# ==============================================================================
# Objetivo: estudiar de forma exploratoria el patrón espacial de las covariables
# climáticas y su relación con el NO2. Para cada covariable se interpola su
# MEDIA ANUAL 2025 sobre Madrid con el método ENSEMBLE (pesos 1/RMSE de 1-NN,
# IDW beta=1 y kNN) y se superponen las estaciones de NO2 como burbujas según su
# concentración media anual. Además se cuantifica la relación interpolando cada
# covariable a las estaciones de NO2 y calculando correlaciones + dispersión.
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

# Predicción ensemble en las ubicaciones objetivo (rejilla o puntos).
superficie_ensemble <- function(sf_est, sf_target, k, pesos) {
  f <- Z ~ 1
  k_real <- min(k, nrow(sf_est))
  p_nn <- idw(f, sf_est, sf_target, nmax = 1, debug.level = 0)$var1.pred
  p_idw <- idw(f, sf_est, sf_target, nmax = k_real, idp = 1, debug.level = 0)$var1.pred
  p_knn <- idw(f, sf_est, sf_target, nmax = k_real, idp = 0, debug.level = 0)$var1.pred
  pesos[1] * p_nn + pesos[2] * p_idw + pesos[3] * p_knn
}

# ------------------------------------------------------------------------------
# 2. Clima: media anual por estación meteorológica
# ------------------------------------------------------------------------------

dt_meteo <- as.data.table(readRDS(
  here("data", "processed", "Clima", "diario",
       sprintf("meteo_madrid_%d_diario.rds", ANIO))
))
normalizar_coordenadas_madrid(dt_meteo)

# Nombres reales de las columnas con acentos (robusto ante el locale).
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

# ------------------------------------------------------------------------------
# 3. NO2: media anual por estación
# ------------------------------------------------------------------------------

dt_no2 <- as.data.table(readRDS(
  here("data", "processed", "Contaminacion", "horario",
       sprintf("aire_madrid_%d_No2_horarios.rds", ANIO))
))
no2_est <- dt_no2[
  , .(NO2 = mean(DATO, na.rm = TRUE)),
  by = .(ESTACION, LONGITUD, LATITUD)
][!is.na(NO2)]

sf_no2 <- st_transform(
  st_as_sf(no2_est, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830
)
xy_no2 <- st_coordinates(sf_no2)
df_no2 <- data.frame(X = xy_no2[, 1], Y = xy_no2[, 2], NO2 = no2_est$NO2)

cat(sprintf(
  "Estaciones clima: %d | Estaciones NO2: %d | NO2 medio: %.1f ug/m3\n",
  nrow(med_est), nrow(no2_est), mean(no2_est$NO2)
))

# ------------------------------------------------------------------------------
# 4. Geometría de Madrid: rejilla y bordes
# ------------------------------------------------------------------------------

mapa <- st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE
))
mapa_utm <- st_transform(mapa, 25830)
rejilla <- st_as_sf(st_make_grid(mapa_utm, n = c(100, 100), what = "centers"))
rejilla_madrid <- st_intersection(rejilla, st_union(mapa_utm))
coords_rejilla <- st_coordinates(rejilla_madrid)
dx <- diff(sort(unique(coords_rejilla[, 1])))[1]
dy <- diff(sort(unique(coords_rejilla[, 2])))[1]

bordes <- st_cast(st_cast(st_geometry(mapa_utm), "MULTILINESTRING"), "LINESTRING")
bordes_coords <- as.data.frame(st_coordinates(bordes))

# ------------------------------------------------------------------------------
# 5. Pesos del ensemble por variable (LOOCV sobre el clima diario)
# ------------------------------------------------------------------------------

pesos_ens <- lapply(VARIABLES, function(v) {
  w <- pesos_ensemble_loocv(dt_meteo, v, K_VECINOS)
  cat(sprintf(
    "Pesos %-22s 1-NN=%.3f IDW b1=%.3f kNN=%.3f\n",
    v, w[1], w[2], w[3]
  ))
  w
})
names(pesos_ens) <- VARIABLES

# ------------------------------------------------------------------------------
# 6. Superficies ensemble + paneles, e interpolación a estaciones NO2
# ------------------------------------------------------------------------------

paneles <- list()
clima_en_no2 <- list()

for (v in VARIABLES) {
  me <- med_est[!is.na(get(v)), .(ESTACION, LONGITUD, LATITUD, Z = get(v))]
  sf_e <- st_transform(
    st_as_sf(me, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830
  )

  # Superficie sobre la rejilla.
  surf <- superficie_ensemble(sf_e, rejilla_madrid, K_VECINOS, pesos_ens[[v]])
  df_s <- data.frame(
    X = coords_rejilla[, 1], Y = coords_rejilla[, 2], Valor = surf
  )

  # Interpolación a las estaciones de NO2 (para las correlaciones).
  clima_en_no2[[v]] <- superficie_ensemble(sf_e, sf_no2, K_VECINOS, pesos_ens[[v]])

  paneles[[v]] <- ggplot() +
    geom_tile(data = df_s, aes(X, Y, fill = Valor), width = dx, height = dy) +
    geom_path(
      data = bordes_coords, aes(X, Y, group = L1),
      color = "white", linewidth = 0.25, inherit.aes = FALSE
    ) +
    geom_point(
      data = df_no2, aes(X, Y, size = NO2),
      shape = 21, fill = "#d73027", color = "white", stroke = 0.5, alpha = 0.9
    ) +
    scale_fill_viridis_c(option = "plasma", name = unname(ETIQUETAS[v])) +
    scale_size_continuous(range = c(1.3, 6), name = "NO₂ medio\n(µg/m³)") +
    coord_equal() +
    labs(title = unname(ETIQUETAS[v])) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.35, "cm")
    )
}

fig <- wrap_plots(paneles, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Covariables climáticas (ensemble) y NO₂ — media anual 2025",
    subtitle = paste(
      "Superficie = interpolación ensemble de la media anual.",
      "Burbujas rojas = NO₂ medio por estación (tamaño ∝ concentración)."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

ggsave(
  file.path(DIR_SALIDA, "mapa_covariables_NO2.png"),
  plot = fig, width = 15, height = 9, dpi = 200, bg = "white"
)

# ------------------------------------------------------------------------------
# 7. Relación cuantitativa: correlaciones y dispersión NO2 vs covariable
# ------------------------------------------------------------------------------

dt_corr <- data.table(ESTACION = no2_est$ESTACION, NO2 = no2_est$NO2)
for (v in VARIABLES) dt_corr[, (v) := clima_en_no2[[v]]]

tabla_cor <- rbindlist(lapply(VARIABLES, function(v) {
  data.table(
    Variable = unname(ETIQUETAS[v]),
    Pearson = round(cor(dt_corr$NO2, dt_corr[[v]], use = "complete.obs"), 3),
    Spearman = round(
      cor(dt_corr$NO2, dt_corr[[v]], method = "spearman", use = "complete.obs"), 3
    )
  )
}))
cat("\nCorrelaciones NO2 vs covariable (interpoladas a estaciones NO2):\n")
print(tabla_cor)
fwrite(tabla_cor, file.path(DIR_SALIDA, "correlaciones_NO2_covariables.csv"))

# Dispersión facetada.
dlong <- melt(
  dt_corr, id.vars = c("ESTACION", "NO2"),
  measure.vars = VARIABLES, variable.name = "Cov", value.name = "Valor"
)
dlong[, Cov := unname(ETIQUETAS[as.character(Cov)])]
etiquetas_cor <- tabla_cor[, .(
  Cov = Variable,
  label = sprintf("r = %.2f | rho = %.2f", Pearson, Spearman)
)]

p_scatter <- ggplot(dlong, aes(Valor, NO2)) +
  geom_point(color = "#2166ac", size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#d73027", linewidth = 0.7) +
  geom_text(
    data = etiquetas_cor, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.4, size = 3.2, inherit.aes = FALSE, fontface = "bold"
  ) +
  facet_wrap(~Cov, scales = "free_x", ncol = 3) +
  labs(
    title = "Relación NO₂ vs covariables climáticas (media anual 2025)",
    subtitle = "Covariables interpoladas (ensemble) a las estaciones de NO₂",
    x = "Valor de la covariable", y = "NO₂ medio (µg/m³)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(DIR_SALIDA, "dispersion_NO2_covariables.png"),
  plot = p_scatter, width = 12, height = 7, dpi = 200, bg = "white"
)

cat("\nSalidas guardadas en:\n", DIR_SALIDA, "\n")
