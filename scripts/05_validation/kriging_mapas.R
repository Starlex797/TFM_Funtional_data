# ==============================================================================
# MAPAS DE KRIGING ORDINARIO (predicción + incertidumbre) - 2025
# ==============================================================================
# Para las variables donde el kriging tiene más sentido (Temperatura, Humedad,
# Radiación solar). Para cada una, en la fecha de referencia:
#   1) ajusta un variograma diario,
#   2) krigea sobre la rejilla 100x100 de Madrid -> superficie de predicción,
#   3) genera además el MAPA DE INCERTIDUMBRE (sd de kriging), que IDW no da.
# ==============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(sf); library(gstat)
  library(ggplot2); library(viridis); library(gridExtra)
})

ANIO       <- 2025L
FECHA_MAPA <- as.Date("2025-04-03")
VARIABLES  <- c("Temperatura", "Humedad_Relativa", "Radiacion_Solar")
ETIQUETAS  <- c(Temperatura = "Temperatura (°C)",
                Humedad_Relativa = "Humedad relativa (%)",
                Radiacion_Solar = "Radiación solar (W/m²)")
DIR_OUT <- here("outputs", "figures", "kriging")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# --- Datos y normalización de coordenadas (idéntica al script de interpolación)-
d <- readRDS(here("data","processed","Clima","diario",
                  sprintf("meteo_madrid_%d_diario.rds", ANIO)))
setDT(d)
col_rad <- grep("Solar$", names(d), value = TRUE)
if (length(col_rad) == 1L) setnames(d, col_rad, "Radiacion_Solar")
fl <- !is.na(d$LATITUD)  & d$LATITUD  < 35
fo <- !is.na(d$LONGITUD) & d$LONGITUD > -1
d[fl, LATITUD  := LATITUD * 10^round(log10(40.4 / abs(LATITUD)))]
d[fo, LONGITUD := LONGITUD * 10]

# --- Rejilla de Madrid (igual que simulacion_interpolacion_clima.R) ------------
mapa <- st_make_valid(st_read(here("data","raw","geometrias","madrid_distritos.geojson"),
                              quiet = TRUE))
mapa_utm <- st_transform(mapa, 25830)
rejilla <- st_as_sf(st_make_grid(mapa_utm, n = c(100,100), what = "centers"))
rejilla_madrid <- st_intersection(rejilla, st_union(mapa_utm))
bordes <- as.data.frame(st_coordinates(
  st_cast(st_cast(st_geometry(mapa_utm), "MULTILINESTRING"), "LINESTRING")))

# --- Bucle por variable --------------------------------------------------------
for (v in VARIABLES) {
  dd <- d[FECHA == FECHA_MAPA & !is.na(get(v)) &
            LATITUD > 40.2 & LATITUD < 40.7 & LONGITUD > -4 & LONGITUD < -3.3]
  sfp <- st_transform(st_as_sf(dd, coords = c("LONGITUD","LATITUD"), crs = 4326), 25830)
  sfp$VALOR <- sfp[[v]]

  # Colapsar estaciones CO-LOCALIZADAS (misma coordenada -> matriz singular en
  # kriging). Se promedia su valor y se conserva un único punto.
  key <- do.call(paste, as.data.frame(round(st_coordinates(sfp))))
  sfp$VALOR <- ave(sfp$VALOR, key, FUN = function(x) mean(x, na.rm = TRUE))
  sfp <- sfp[!duplicated(key), ]

  # Variograma + ajuste (con reserva exponencial si falla)
  vemp <- variogram(VALOR ~ 1, sfp)
  psill0 <- var(sfp$VALOR); rng0 <- as.numeric(median(dist(st_coordinates(sfp))))
  m <- tryCatch(suppressWarnings(fit.variogram(vemp, vgm(psill0,"Sph",rng0,0.1*psill0))),
                error = function(e) NULL)
  if (is.null(m) || any(m$psill < 0)) m <- vgm(psill0,"Exp",rng0,0.1*psill0)

  # Kriging sobre la rejilla -> predicción + varianza
  kr <- krige(VALOR ~ 1, sfp, newdata = rejilla_madrid, model = m, debug.level = 0)
  co <- st_coordinates(rejilla_madrid)
  df <- data.frame(X = co[,1], Y = co[,2],
                   pred = kr$var1.pred, sd = sqrt(kr$var1.var))
  est <- as.data.frame(st_coordinates(sfp)); est$VALOR <- sfp$VALOR
  etq <- ETIQUETAS[[v]]

  base_theme <- theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.title = element_blank(), plot.title = element_text(face="bold", size=12))

  # (1) Mapa de PREDICCIÓN
  p_pred <- ggplot(df, aes(X, Y, fill = pred)) +
    geom_raster() +
    geom_path(data = bordes, aes(X, Y, group = L1), color = "white",
              linewidth = 0.3, inherit.aes = FALSE) +
    geom_point(data = est, aes(X, Y, fill = VALOR), shape = 21, size = 2.6,
               color = "black", stroke = 0.6, inherit.aes = FALSE) +
    scale_fill_viridis_c(option = "plasma", name = etq) +
    coord_equal() +
    labs(title = paste("Kriging ordinario —", v),
         subtitle = sprintf("%s | predicción | %d estaciones",
                            format(FECHA_MAPA, "%d %b %Y"), nrow(sfp))) +
    base_theme

  # (2) Mapa de INCERTIDUMBRE (sd de kriging)
  p_sd <- ggplot(df, aes(X, Y, fill = sd)) +
    geom_raster() +
    geom_path(data = bordes, aes(X, Y, group = L1), color = "white",
              linewidth = 0.3, inherit.aes = FALSE) +
    geom_point(data = est, aes(X, Y), shape = 21, size = 2.2, fill = "white",
               color = "black", stroke = 0.6, inherit.aes = FALSE) +
    scale_fill_viridis_c(option = "viridis", direction = -1, name = "sd kriging") +
    coord_equal() +
    labs(title = paste("Incertidumbre del kriging —", v),
         subtitle = "Desviación típica de predicción (mayor lejos de estaciones)") +
    base_theme

  ggsave(file.path(DIR_OUT, sprintf("kriging_prediccion_%s.png", v)),
         p_pred, width = 8, height = 6.5, dpi = 200, bg = "white")
  ggsave(file.path(DIR_OUT, sprintf("kriging_incertidumbre_%s.png", v)),
         p_sd, width = 8, height = 6.5, dpi = 200, bg = "white")
  # Combinado lado a lado
  ggsave(file.path(DIR_OUT, sprintf("kriging_panel_%s.png", v)),
         arrangeGrob(p_pred, p_sd, ncol = 2), width = 15, height = 6.5, dpi = 170, bg = "white")

  cat(sprintf("%-18s | variograma: %s | rango=%.0f m | pred [%.1f, %.1f] | sd [%.2f, %.2f]\n",
      v, m$model[nrow(m)], m$range[nrow(m)],
      min(df$pred), max(df$pred), min(df$sd), max(df$sd)))
}
cat("\nMapas guardados en:", DIR_OUT, "\n")
