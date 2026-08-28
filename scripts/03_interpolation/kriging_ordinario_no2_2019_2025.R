# ==============================================================================
# Kriging Ordinario (escalar) del NO2 en Madrid, 2019-2025
# ==============================================================================
# Baseline "solo-espacio" para comparar con el Kriging Funcional Ordinario (OFK)
# de kriging_funcional_no2.R. En el paper de Montero-Lorenzo et al. (2013) es la
# columna "Spatial" de la Tabla 3: kriging ordinario de un ESCALAR por estacion.
#
# Aqui el escalar es la MEDIA ANUAL de log-NO2 por estacion (misma variable y
# misma escala logaritmica que usa el OFK), de modo que el RMSE de la validacion
# LOOCV es directamente comparable entre ambos metodos y entre años.
#
# Para cada año 2019..2025:
#   1. Media anual de Ln(NO2) por estacion.
#   2. Variograma empirico + ajuste esferico.
#   3. Kriging ordinario sobre la rejilla de Madrid (misma rejilla que el OFK).
#   4. LOOCV espacial -> RMSE (log), global y por NOM_TIPO.
# Salidas: mapas por año (comparables), evolucion del RMSE y de la media, y CSV.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(gstat)
  library(ggplot2)
  library(viridis)
})
sf_use_s2(FALSE)

ANIOS      <- 2019:2025
DIR_SALIDA <- here("outputs", "figures", "kriging_ordinario_no2")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

col_tipo <- c("Urbana tráfico" = "#d73027", "Urbana fondo" = "#4575b4",
              "Suburbana" = "#1a9850")

# ------------------------------------------------------------------------------
# 1. Rejilla y bordes de Madrid (idénticos al script OFK, para comparabilidad)
# ------------------------------------------------------------------------------

mapa_utm <- st_transform(st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE
)), 25830)
rej <- st_as_sf(st_make_grid(mapa_utm, n = c(90, 90), what = "centers"))
rej <- st_intersection(rej, st_union(mapa_utm))
GR  <- st_coordinates(rej)
dx  <- diff(sort(unique(GR[, 1])))[1]
dy  <- diff(sort(unique(GR[, 2])))[1]
bordes <- as.data.frame(st_coordinates(
  st_cast(st_cast(st_geometry(mapa_utm), "MULTILINESTRING"), "LINESTRING")))

# ------------------------------------------------------------------------------
# 2. Bucle por año: media anual de log-NO2, variograma, kriging y LOOCV
# ------------------------------------------------------------------------------

pred_all <- list(); est_all <- list(); rmse_all <- list(); vgm_all <- list()

for (a in ANIOS) {
  dt <- as.data.table(readRDS(here(
    "data", "processed", "Contaminacion", "horario",
    sprintf("aire_madrid_%d_No2_horarios.rds", a)
  )))

  # Media anual de Ln(NO2) por estacion (mismo objetivo, en log, que el OFK).
  est <- dt[!is.na(LOG_NO2_HORARIO), .(
    z = mean(LOG_NO2_HORARIO),
    media = mean(DATO, na.rm = TRUE)
  ), by = .(ESTACION, NOM_TIPO, LONGITUD, LATITUD)]

  sf_est <- st_transform(
    st_as_sf(est, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830)

  # Variograma empirico + ajuste esferico (con arranque robusto y fallback).
  v <- variogram(z ~ 1, sf_est)
  ini <- vgm(psill = 0.8 * var(est$z), model = "Sph",
             range = 8000, nugget = 0.2 * var(est$z))
  fit <- tryCatch(suppressWarnings(fit.variogram(v, ini)),
                  error = function(e) ini)
  if (any(fit$range < 0) || sum(fit$psill) <= 0) fit <- ini
  vgm_all[[as.character(a)]] <-
    data.table(v)[, .(anio = a, h = dist, gamma, np)]

  # Kriging ordinario sobre la rejilla.
  k <- krige(z ~ 1, sf_est, rej, model = fit, debug.level = 0)
  pred_all[[as.character(a)]] <- data.table(
    anio = a, X = GR[, 1], Y = GR[, 2],
    NO2 = exp(k$var1.pred),          # back-transf. (media geometrica)
    sd_log = sqrt(k$var1.var)        # incertidumbre de kriging (log)
  )

  # Puntos de estacion (para las burbujas del mapa).
  xy <- st_coordinates(sf_est)
  est_all[[as.character(a)]] <- data.table(
    anio = a, X = xy[, 1], Y = xy[, 2],
    NO2 = est$media, NOM_TIPO = est$NOM_TIPO)

  # Validacion LOOCV espacial -> residuo = observado - predicho (en log).
  cv <- suppressWarnings(krige.cv(z ~ 1, sf_est, model = fit, verbose = FALSE))
  rcv <- data.table(NOM_TIPO = est$NOM_TIPO, res = cv$residual)
  rmse_all[[as.character(a)]] <- rbind(
    data.table(anio = a, grupo = "GLOBAL",
               RMSE_log = sqrt(mean(rcv$res^2)), n = nrow(rcv),
               media_anual = mean(est$media)),
    rcv[, .(anio = a, grupo = NOM_TIPO,
            RMSE_log = sqrt(mean(res^2)), n = .N,
            media_anual = mean(est$media)), by = NOM_TIPO][, -1]
  )

  cat(sprintf("%d | n=%2d | media=%.1f ug/m3 | rango vgm=%.0f m | RMSE_LOOCV(log)=%.4f\n",
              a, nrow(est), mean(est$media),
              fit$range[nrow(fit)], sqrt(mean(rcv$res^2))))
}

pred_dt <- rbindlist(pred_all)
est_dt  <- rbindlist(est_all)
rmse_dt <- rbindlist(rmse_all)
fwrite(rmse_dt, file.path(DIR_SALIDA, "loocv_rmse_por_anio.csv"))

# ------------------------------------------------------------------------------
# 3. Mapa: superficie de NO2 kriged por año (escala común -> comparable)
# ------------------------------------------------------------------------------

p_maps <- ggplot() +
  geom_tile(data = pred_dt, aes(X, Y, fill = NO2), width = dx, height = dy) +
  geom_path(data = bordes, aes(X, Y, group = L1),
            color = "white", linewidth = 0.15, inherit.aes = FALSE) +
  geom_point(data = est_dt, aes(X, Y, size = NO2), shape = 21,
             fill = "white", color = "black", stroke = 0.3, alpha = 0.8) +
  facet_wrap(~anio, ncol = 4) +
  scale_fill_viridis_c(option = "plasma", name = "NO₂\n(µg/m³)") +
  scale_size_continuous(range = c(0.5, 3.5), guide = "none") +
  coord_equal() +
  labs(title = "Kriging ordinario del NO₂ en Madrid (media anual, 2019–2025)",
       subtitle = "Escala común entre años. Superficie = kriging de Ln(NO₂); puntos = estaciones.") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey30"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "right")
ggsave(file.path(DIR_SALIDA, "mapas_kriging_ordinario_por_anio.png"), p_maps,
       width = 15, height = 8, dpi = 200, bg = "white")

# ------------------------------------------------------------------------------
# 4. Mapa: incertidumbre de kriging (sd en log) por año
# ------------------------------------------------------------------------------

p_sd <- ggplot() +
  geom_tile(data = pred_dt, aes(X, Y, fill = sd_log), width = dx, height = dy) +
  geom_path(data = bordes, aes(X, Y, group = L1),
            color = "white", linewidth = 0.15, inherit.aes = FALSE) +
  geom_point(data = est_dt, aes(X, Y), shape = 21, size = 1,
             fill = "white", color = "black", stroke = 0.3) +
  facet_wrap(~anio, ncol = 4) +
  scale_fill_viridis_c(option = "viridis", name = "sd\n(log)") +
  coord_equal() +
  labs(title = "Incertidumbre del kriging ordinario (desv. típica, log-NO₂)",
       subtitle = "Crece al alejarse de las estaciones — comparar con la trace-variance del OFK") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey30"),
        strip.text = element_text(face = "bold", size = 12))
ggsave(file.path(DIR_SALIDA, "mapas_incertidumbre_por_anio.png"), p_sd,
       width = 15, height = 8, dpi = 200, bg = "white")

# ------------------------------------------------------------------------------
# 5. Evolución del RMSE LOOCV y de la media anual (para la comparación)
# ------------------------------------------------------------------------------

rmse_glob <- rmse_dt[grupo == "GLOBAL"]
p_rmse <- ggplot() +
  geom_line(data = rmse_dt[grupo != "GLOBAL"],
            aes(anio, RMSE_log, color = grupo), linewidth = 0.7, alpha = 0.7) +
  geom_point(data = rmse_dt[grupo != "GLOBAL"],
             aes(anio, RMSE_log, color = grupo), size = 2) +
  geom_line(data = rmse_glob, aes(anio, RMSE_log), color = "black",
            linewidth = 1.2) +
  geom_point(data = rmse_glob, aes(anio, RMSE_log), color = "black", size = 2.8) +
  scale_color_manual(values = col_tipo, name = "Tipología") +
  scale_x_continuous(breaks = ANIOS) +
  labs(title = "RMSE de la validación LOOCV del kriging ordinario (log-NO₂)",
       subtitle = "Línea negra = global. Métrica comparable con el OFK y entre años.",
       x = NULL, y = "RMSE integrado (log)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_SALIDA, "evolucion_rmse_loocv.png"), p_rmse,
       width = 10, height = 6, dpi = 200, bg = "white")

cat("\nResumen RMSE LOOCV (log) por año:\n")
print(dcast(rmse_dt, anio + media_anual ~ grupo, value.var = "RMSE_log"))
cat("\nSalidas guardadas en:\n", DIR_SALIDA, "\n")
