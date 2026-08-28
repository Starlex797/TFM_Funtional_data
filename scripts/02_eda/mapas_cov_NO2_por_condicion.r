# ==============================================================================
# EDA espacio-temporal: covariables (ensemble) vs NO2 por CONDICIÓN
# ==============================================================================
# La media anual mezcla condiciones opuestas y borra la señal. Aquí se estratifica
# para poder sacar conclusiones:
#   - Por estación del año (invierno, primavera, verano, otoño)
#   - Hora punta (franja de mañana 7-9 h, pico de tráfico) a escala horaria
#   - Un día concreto con lluvia y viento (para precipitación y viento)
# Para cada condición se interpola cada covariable con el ENSEMBLE (media/estación
# de la condición), se superpone el NO2 de esa misma condición y se calcula la
# correlación espacial NO2 vs covariable. Todo se resume en un heatmap.
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
HORAS_PUNTA <- 7:9 # franja de mañana (pico de tráfico)
DIR_SALIDA <- here("outputs", "figures", "mapas_cov_NO2", "por_condicion")
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

superficie_ensemble <- function(sf_est, sf_target, k, pesos) {
  f <- Z ~ 1
  k_real <- min(k, nrow(sf_est))
  p_nn <- idw(f, sf_est, sf_target, nmax = 1, debug.level = 0)$var1.pred
  p_idw <- idw(f, sf_est, sf_target, nmax = k_real, idp = 1, debug.level = 0)$var1.pred
  p_knn <- idw(f, sf_est, sf_target, nmax = k_real, idp = 0, debug.level = 0)$var1.pred
  pesos[1] * p_nn + pesos[2] * p_idw + pesos[3] * p_knn
}

estacion_anio <- function(mes) {
  fifelse(mes %in% c(12, 1, 2), "Invierno",
    fifelse(mes %in% 3:5, "Primavera",
      fifelse(mes %in% 6:8, "Verano", "Otoño")
    )
  )
}

sf_puntos <- function(dt) {
  st_transform(st_as_sf(dt, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830)
}

# ------------------------------------------------------------------------------
# 2. Datos
# ------------------------------------------------------------------------------

dt_meteo_d <- as.data.table(readRDS(
  here("data", "processed", "Clima", "diario",
       sprintf("meteo_madrid_%d_diario.rds", ANIO))
))
normalizar_coordenadas_madrid(dt_meteo_d)

v_pres <- grep("^Presion Bar", names(dt_meteo_d), value = TRUE)
v_rad <- grep("Solar$", names(dt_meteo_d), value = TRUE)
VARIABLES <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  v_pres, v_rad, "Velocidad Viento"
)
ETIQUETAS <- c(
  "Temperatura" = "Temperatura (°C)",
  "Humedad_Relativa" = "Humedad rel. (%)",
  "Precipitaciones" = "Precipitación (mm)",
  "Velocidad Viento" = "Viento (m/s)"
)
ETIQUETAS[v_pres] <- "Presión (hPa)"
ETIQUETAS[v_rad] <- "Radiación (W/m²)"

# NO2 horario (para agregarlo por condición)
dt_no2 <- as.data.table(readRDS(
  here("data", "processed", "Contaminacion", "horario",
       sprintf("aire_madrid_%d_No2_horarios.rds", ANIO))
))
dt_no2[, HORA_INT := as.integer(gsub("[^0-9]", "", as.character(HORA)))]

# Meteo horaria (para la hora punta)
dt_meteo_h <- as.data.table(readRDS(
  here("data", "processed", "Clima", "horario",
       sprintf("meteo_madrid_%d_horario.rds", ANIO))
))
normalizar_coordenadas_madrid(dt_meteo_h)
if (is.factor(dt_meteo_h$HORA) || is.character(dt_meteo_h$HORA)) {
  dt_meteo_h[, HORA := as.integer(gsub("[^0-9]", "", as.character(HORA)))]
}

# ------------------------------------------------------------------------------
# 3. Geometría
# ------------------------------------------------------------------------------

mapa_utm <- st_transform(st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE
)), 25830)
rejilla <- st_as_sf(st_make_grid(mapa_utm, n = c(100, 100), what = "centers"))
rejilla_madrid <- st_intersection(rejilla, st_union(mapa_utm))
coords_rejilla <- st_coordinates(rejilla_madrid)
dx <- diff(sort(unique(coords_rejilla[, 1])))[1]
dy <- diff(sort(unique(coords_rejilla[, 2])))[1]
bordes_coords <- as.data.frame(st_coordinates(
  st_cast(st_cast(st_geometry(mapa_utm), "MULTILINESTRING"), "LINESTRING")
))

# ------------------------------------------------------------------------------
# 4. Pesos del ensemble (una vez, sobre el clima diario; método = propiedad fija)
# ------------------------------------------------------------------------------

pesos_ens <- lapply(VARIABLES, function(v) pesos_ensemble_loocv(dt_meteo_d, v, K_VECINOS))
names(pesos_ens) <- VARIABLES

# ------------------------------------------------------------------------------
# 5. Funciones de análisis por condición
# ------------------------------------------------------------------------------

# Correlación espacial NO2 vs cada covariable, interpolando la covariable a las
# estaciones de NO2 de esa condición.
correlaciones_condicion <- function(med_clima, no2_cond, condicion,
                                    variables = VARIABLES) {
  sf_no2c <- sf_puntos(no2_cond)
  res <- lapply(variables, function(v) {
    me <- med_clima[!is.na(get(v)), .(ESTACION, LONGITUD, LATITUD, Z = get(v))]
    if (nrow(me) < 4L) return(NULL)
    pred <- superficie_ensemble(sf_puntos(me), sf_no2c, K_VECINOS, pesos_ens[[v]])
    data.table(
      Condicion = condicion,
      Variable = unname(ETIQUETAS[v]),
      Pearson = cor(no2_cond$NO2, pred, use = "complete.obs"),
      Spearman = cor(no2_cond$NO2, pred, method = "spearman", use = "complete.obs")
    )
  })
  rbindlist(res)
}

# Panel de mapa de una covariable con burbujas de NO2 (escalas opcionales
# compartidas entre paneles comparables).
panel_mapa <- function(med_clima, v, no2_cond, titulo,
                       lim_fill = NULL, lim_no2 = NULL, mostrar_leyenda = TRUE) {
  me <- med_clima[!is.na(get(v)), .(ESTACION, LONGITUD, LATITUD, Z = get(v))]
  surf <- superficie_ensemble(sf_puntos(me), rejilla_madrid, K_VECINOS, pesos_ens[[v]])
  df_s <- data.frame(X = coords_rejilla[, 1], Y = coords_rejilla[, 2], Valor = surf)
  sfn <- sf_puntos(no2_cond)
  xy <- st_coordinates(sfn)
  df_n <- data.frame(X = xy[, 1], Y = xy[, 2], NO2 = no2_cond$NO2)

  ggplot() +
    geom_tile(data = df_s, aes(X, Y, fill = Valor), width = dx, height = dy) +
    geom_path(data = bordes_coords, aes(X, Y, group = L1),
              color = "white", linewidth = 0.2, inherit.aes = FALSE) +
    geom_point(data = df_n, aes(X, Y, size = NO2), shape = 21,
               fill = "#d73027", color = "white", stroke = 0.4, alpha = 0.9) +
    scale_fill_viridis_c(option = "plasma", limits = lim_fill,
                         name = unname(ETIQUETAS[v])) +
    scale_size_continuous(range = c(1, 5), limits = lim_no2,
                          name = "NO₂\n(µg/m³)") +
    coord_equal() +
    labs(title = titulo) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = if (mostrar_leyenda) "right" else "none",
      legend.key.size = unit(0.3, "cm")
    )
}

# Figura de 4 paneles (una covariable, las 4 estaciones del año) con escalas
# compartidas para que sean comparables.
figura_estacional <- function(v, med_por_est, no2_por_est) {
  ests <- c("Invierno", "Primavera", "Verano", "Otoño")
  # Rango común de la covariable y del NO2 entre estaciones.
  vals <- unlist(lapply(ests, function(s) {
    me <- med_por_est[EST == s & !is.na(get(v)), .(ESTACION, LONGITUD, LATITUD, Z = get(v))]
    superficie_ensemble(sf_puntos(me), rejilla_madrid, K_VECINOS, pesos_ens[[v]])
  }))
  lim_fill <- range(vals, na.rm = TRUE)
  lim_no2 <- range(unlist(lapply(ests, function(s) no2_por_est[[s]]$NO2)), na.rm = TRUE)

  paneles <- lapply(ests, function(s) {
    panel_mapa(
      med_por_est[EST == s], v, no2_por_est[[s]], titulo = s,
      lim_fill = lim_fill, lim_no2 = lim_no2
    )
  })
  wrap_plots(paneles, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("NO₂ y ", unname(ETIQUETAS[v]), " por estación del año — ", ANIO),
      subtitle = "Superficie = ensemble de la media estacional | burbujas = NO₂ medio estacional",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
}

# ------------------------------------------------------------------------------
# 6. Estaciones del año
# ------------------------------------------------------------------------------

dt_meteo_d[, EST := estacion_anio(month(FECHA))]
dt_no2[, EST := estacion_anio(month(FECHA))]

med_por_est <- dt_meteo_d[
  , lapply(.SD, mean, na.rm = TRUE),
  by = .(EST, ESTACION, LONGITUD, LATITUD), .SDcols = VARIABLES
]
no2_por_est_dt <- dt_no2[
  , .(NO2 = mean(DATO, na.rm = TRUE)),
  by = .(EST, ESTACION, LONGITUD, LATITUD)
][!is.na(NO2)]
no2_por_est <- split(no2_por_est_dt, no2_por_est_dt$EST)

# Correlaciones por estación (todas las covariables).
corr_estaciones <- rbindlist(lapply(
  c("Invierno", "Primavera", "Verano", "Otoño"),
  function(s) correlaciones_condicion(med_por_est[EST == s], no2_por_est[[s]], s)
))

# Mapas estacionales (covariables continuas más informativas).
for (v in c("Temperatura", "Humedad_Relativa", v_rad)) {
  fig <- figura_estacional(v, med_por_est, no2_por_est)
  nombre <- gsub("[^A-Za-z0-9]", "_", v)
  ggsave(file.path(DIR_SALIDA, paste0("mapa_estacional_", nombre, ".png")),
         plot = fig, width = 11, height = 9, dpi = 200, bg = "white")
  cat("Mapa estacional guardado:", v, "\n")
}

# ------------------------------------------------------------------------------
# 7. Hora punta (franja de mañana, escala horaria)
# ------------------------------------------------------------------------------

med_punta <- dt_meteo_h[HORA %in% HORAS_PUNTA][
  , lapply(.SD, mean, na.rm = TRUE),
  by = .(ESTACION, LONGITUD, LATITUD), .SDcols = VARIABLES
]
no2_punta <- dt_no2[HORA_INT %in% HORAS_PUNTA][
  , .(NO2 = mean(DATO, na.rm = TRUE)),
  by = .(ESTACION, LONGITUD, LATITUD)
][!is.na(NO2)]

corr_punta <- correlaciones_condicion(med_punta, no2_punta, "Hora punta 7-9h")

# Mapa hora punta: temperatura y humedad.
fig_punta <- (panel_mapa(med_punta, "Temperatura", no2_punta, "Temperatura") |
  panel_mapa(med_punta, "Humedad_Relativa", no2_punta, "Humedad relativa")) +
  plot_annotation(
    title = sprintf("NO₂ y clima en hora punta (%d-%d h) — %d",
                    min(HORAS_PUNTA), max(HORAS_PUNTA), ANIO),
    subtitle = "Media en la franja de mañana | burbujas = NO₂ medio en esa franja",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
ggsave(file.path(DIR_SALIDA, "mapa_hora_punta.png"),
       plot = fig_punta, width = 12, height = 6, dpi = 200, bg = "white")
cat("Mapa hora punta guardado.\n")

# ------------------------------------------------------------------------------
# 8. Día con lluvia y viento (para precipitación y viento)
# ------------------------------------------------------------------------------

resumen_dia <- dt_meteo_d[
  , .(
    precip = mean(Precipitaciones, na.rm = TRUE),
    viento = median(`Velocidad Viento`, na.rm = TRUE),
    n_precip = sum(!is.na(Precipitaciones))
  ),
  by = FECHA
][n_precip >= 7]
# Día que combina lluvia y viento: máximo del producto de rangos normalizados.
resumen_dia[, score := scale(precip)[, 1] + scale(viento)[, 1]]
dia_evento <- resumen_dia[which.max(score), FECHA]
cat(sprintf(
  "Día evento (lluvia+viento): %s | precip media=%.2f mm | viento mediana=%.2f m/s\n",
  dia_evento, resumen_dia[FECHA == dia_evento, precip],
  resumen_dia[FECHA == dia_evento, viento]
))

med_evento <- dt_meteo_d[FECHA == dia_evento,
  c("ESTACION", "LONGITUD", "LATITUD", VARIABLES), with = FALSE
]
no2_evento <- dt_no2[FECHA == dia_evento,
  .(NO2 = mean(DATO, na.rm = TRUE)),
  by = .(ESTACION, LONGITUD, LATITUD)
][!is.na(NO2)]

corr_evento <- correlaciones_condicion(
  med_evento, no2_evento, paste0("Día evento ", format(dia_evento, "%d-%b"))
)

fig_evento <- (panel_mapa(med_evento, "Precipitaciones", no2_evento, "Precipitación") |
  panel_mapa(med_evento, "Velocidad Viento", no2_evento, "Viento")) +
  plot_annotation(
    title = sprintf("NO₂, lluvia y viento — %s", format(dia_evento, "%d %b %Y")),
    subtitle = "Día seleccionado por combinar lluvia y viento apreciables",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
ggsave(file.path(DIR_SALIDA, "mapa_evento_lluvia_viento.png"),
       plot = fig_evento, width = 12, height = 6, dpi = 200, bg = "white")
cat("Mapa día evento guardado.\n")

# ------------------------------------------------------------------------------
# 9. Heatmap resumen: correlación Pearson covariable x condición
# ------------------------------------------------------------------------------

corr_total <- rbindlist(list(corr_estaciones, corr_punta, corr_evento))
orden_cond <- c(
  "Invierno", "Primavera", "Verano", "Otoño", "Hora punta 7-9h",
  paste0("Día evento ", format(dia_evento, "%d-%b"))
)
orden_var <- unname(ETIQUETAS[VARIABLES])
corr_total[, Condicion := factor(Condicion, levels = orden_cond)]
corr_total[, Variable := factor(Variable, levels = rev(orden_var))]

fwrite(corr_total, file.path(DIR_SALIDA, "correlaciones_por_condicion.csv"))

p_heat <- ggplot(corr_total, aes(Condicion, Variable, fill = Pearson)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Pearson)), size = 3.3) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-1, 1), name = "Pearson\nNO₂ vs cov."
  ) +
  labs(
    title = "Correlación espacial NO₂ vs covariables por condición",
    subtitle = "Correlación entre estaciones (Pearson). Rojo = positiva, azul = negativa.",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(file.path(DIR_SALIDA, "heatmap_correlaciones_condicion.png"),
       plot = p_heat, width = 10, height = 6, dpi = 200, bg = "white")

cat("\nCorrelaciones por condición:\n")
print(dcast(corr_total, Variable ~ Condicion, value.var = "Pearson"))
cat("\nSalidas en:\n", DIR_SALIDA, "\n")
