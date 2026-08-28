# ==============================================================================
# KRIGING PARA UN SOLO DÍA (rápido): comparación de métodos + mapas
# ==============================================================================
# Para cada variable se elige UN día representativo y en él se hace:
#   (1) Comparación LOOCV de 6 métodos (Media, 1-NN, kNN, IDW b1, IDW b2, Kriging)
#   (2) Mapa de predicción por kriging + mapa de incertidumbre (sd)
#
# Días:
#   - Temperatura, Humedad, Presión, Radiación : día de referencia 2025-04-03
#   - Precipitaciones : el día de LLUVIA MÁS ABUNDANTE (mayor precipitación media)
#   - Velocidad viento: el día MÁS VENTOSO (mayor mediana de viento)
#
# Nota: al validar en un solo día la muestra es pequeña (pocas estaciones), así
# que el RMSE es orientativo; sirve para ilustrar, no para decidir el método
# definitivo (eso requiere la validación anual).
# ==============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(sf); library(gstat)
  library(ggplot2); library(viridis); library(gridExtra)
})

ANIO <- 2025L
K <- 3L
FECHA_REF <- as.Date("2025-04-03")
MIN_EST <- 7L                       # mínimo de estaciones para elegir día
DIR_FIG <- here("outputs", "figures", "kriging")
DIR_TAB <- here("outputs", "tables", "calidad_datos")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

ETI <- c(Temperatura = "Temperatura (°C)", Humedad_Relativa = "Humedad relativa (%)",
         Precipitaciones = "Precipitación (mm)", Presion_Barometrica = "Presión (hPa)",
         Radiacion_Solar = "Radiación solar (W/m²)", `Velocidad Viento` = "Viento (m/s)")
# Cotas físicas para recortar la predicción en el mapa (kriging puede salirse)
COTA_MIN <- c(Temperatura = -Inf, Humedad_Relativa = 0, Precipitaciones = 0,
              Presion_Barometrica = -Inf, Radiacion_Solar = 0, `Velocidad Viento` = 0)

rmse <- function(r, p) sqrt(mean((r - p)^2))
mae  <- function(r, p) mean(abs(r - p))

# --- Datos + normalización de coordenadas -------------------------------------
d <- readRDS(here("data","processed","Clima","diario",
                  sprintf("meteo_madrid_%d_diario.rds", ANIO)))
setDT(d)
cp <- grep("^Presion Barom", names(d), value = TRUE)
cr <- grep("Solar$", names(d), value = TRUE)
if (length(cp) == 1L) setnames(d, cp, "Presion_Barometrica")
if (length(cr) == 1L) setnames(d, cr, "Radiacion_Solar")
fl <- !is.na(d$LATITUD)  & d$LATITUD  < 35
fo <- !is.na(d$LONGITUD) & d$LONGITUD > -1
d[fl, LATITUD  := LATITUD * 10^round(log10(40.4 / abs(LATITUD)))]
d[fo, LONGITUD := LONGITUD * 10]
d <- d[LATITUD > 40.2 & LATITUD < 40.7 & LONGITUD > -4 & LONGITUD < -3.3]

# --- Selección del día de lluvia abundante y del día más ventoso --------------
res_dia <- d[, .(
  n_precip = sum(!is.na(Precipitaciones)),
  precip_media = mean(Precipitaciones, na.rm = TRUE),
  n_viento = sum(!is.na(`Velocidad Viento`)),
  viento_mediana = median(`Velocidad Viento`, na.rm = TRUE)
), by = FECHA]

fecha_lluvia <- res_dia[n_precip >= MIN_EST][which.max(precip_media), FECHA]
fecha_viento <- res_dia[n_viento >= MIN_EST][which.max(viento_mediana), FECHA]

cat(sprintf("Día de lluvia más abundante: %s (%.1f mm de media)\n",
            fecha_lluvia, res_dia[FECHA == fecha_lluvia, precip_media]))
cat(sprintf("Día más ventoso: %s (mediana %.2f m/s)\n\n",
            fecha_viento, res_dia[FECHA == fecha_viento, viento_mediana]))

VARIABLES <- c("Temperatura", "Humedad_Relativa", "Precipitaciones",
               "Presion_Barometrica", "Radiacion_Solar", "Velocidad Viento")
DIA_VAR <- c(Temperatura = FECHA_REF, Humedad_Relativa = FECHA_REF,
             Precipitaciones = fecha_lluvia, Presion_Barometrica = FECHA_REF,
             Radiacion_Solar = FECHA_REF, `Velocidad Viento` = fecha_viento)

# --- Rejilla de Madrid --------------------------------------------------------
mapa <- st_make_valid(st_read(here("data","raw","geometrias","madrid_distritos.geojson"),
                              quiet = TRUE))
mapa_utm <- st_transform(mapa, 25830)
rejilla_madrid <- st_intersection(
  st_as_sf(st_make_grid(mapa_utm, n = c(100,100), what = "centers")),
  st_union(mapa_utm))
bordes <- as.data.frame(st_coordinates(
  st_cast(st_cast(st_geometry(mapa_utm), "MULTILINESTRING"), "LINESTRING")))
co_rej <- st_coordinates(rejilla_madrid)

# --- Función: prepara los puntos del día (normaliza + colapsa co-localizadas) --
puntos_dia <- function(v, fecha) {
  dd <- d[FECHA == fecha & !is.na(get(v))]
  sfp <- st_transform(st_as_sf(dd, coords = c("LONGITUD","LATITUD"), crs = 4326), 25830)
  xy <- st_coordinates(sfp); val <- sfp[[v]]
  key <- do.call(paste, as.data.frame(round(xy)))
  val <- ave(val, key, FUN = function(x) mean(x, na.rm = TRUE))
  keep <- !duplicated(key)
  list(xy = xy[keep, , drop = FALSE], val = val[keep])
}

filas <- list()

for (v in VARIABLES) {
  fecha <- as.Date(DIA_VAR[[v]])
  P <- puntos_dia(v, fecha)
  xy <- P$xy; val <- P$val; n <- length(val)
  if (n < 5) { cat("  ", v, ": muy pocas estaciones, se omite.\n"); next }

  # ---------- (1) LOOCV de los 6 métodos en ese día ----------
  D <- as.matrix(dist(xy)); diag(D) <- Inf
  pred_media <- (sum(val) - val) / (n - 1)
  pred_nn <- pred_knn <- pred_b1 <- pred_b2 <- numeric(n)
  for (i in seq_len(n)) {
    idx <- order(D[i, ])[seq_len(K)]; vv <- val[idx]; dv2 <- D[i, idx]
    pred_nn[i] <- vv[1]; pred_knn[i] <- mean(vv)
    w1 <- 1/pmax(dv2,1)^1; pred_b1[i] <- sum(w1*vv)/sum(w1)
    w2 <- 1/pmax(dv2,1)^2; pred_b2[i] <- sum(w2*vv)/sum(w2)
  }
  sfd <- st_as_sf(data.frame(x=xy[,1], y=xy[,2], VALOR=val), coords=c("x","y"), crs=25830)
  vemp <- variogram(VALOR ~ 1, sfd)
  ps <- var(val); rg <- as.numeric(median(dist(xy)))
  m <- tryCatch(suppressWarnings(fit.variogram(vemp, vgm(ps,"Sph",rg,0.1*ps))),
                error=function(e) NULL)
  if (is.null(m) || any(m$psill < 0)) m <- vgm(ps,"Exp",rg,0.1*ps)
  cv <- suppressWarnings(krige.cv(VALOR ~ 1, sfd, model=m, verbose=FALSE))

  rmses <- c(Media=rmse(val,pred_media), `1-NN`=rmse(val,pred_nn), kNN=rmse(val,pred_knn),
             `IDW b=1`=rmse(val,pred_b1), `IDW b=2`=rmse(val,pred_b2),
             Kriging=rmse(cv$observed, cv$var1.pred))
  filas[[v]] <- data.table(
    Variable=v, Fecha=as.character(fecha), N_estaciones=n,
    RMSE_Media=round(rmses["Media"],3), RMSE_NN=round(rmses["1-NN"],3),
    RMSE_kNN=round(rmses["kNN"],3), RMSE_IDW_b1=round(rmses["IDW b=1"],3),
    RMSE_IDW_b2=round(rmses["IDW b=2"],3), RMSE_Kriging=round(rmses["Kriging"],3),
    Mejor_Metodo=names(which.min(rmses)))

  # ---------- (2) Mapas de kriging (predicción + incertidumbre) ----------
  kr <- krige(VALOR ~ 1, sfd, newdata=rejilla_madrid, model=m, debug.level=0)
  pred <- pmax(kr$var1.pred, COTA_MIN[[v]])
  if (v == "Humedad_Relativa") pred <- pmin(pred, 100)
  df <- data.frame(X=co_rej[,1], Y=co_rej[,2], pred=pred, sd=sqrt(kr$var1.var))
  est <- data.frame(X=xy[,1], Y=xy[,2], VALOR=val); etq <- ETI[[v]]
  bt <- theme_minimal(base_size=11) + theme(panel.grid=element_blank(),
        axis.text=element_blank(), axis.title=element_blank(),
        plot.title=element_text(face="bold", size=12))

  p_pred <- ggplot(df, aes(X,Y,fill=pred)) + geom_raster() +
    geom_path(data=bordes, aes(X,Y,group=L1), color="white", linewidth=0.3, inherit.aes=FALSE) +
    geom_point(data=est, aes(X,Y,fill=VALOR), shape=21, size=2.6, color="black", stroke=0.6, inherit.aes=FALSE) +
    scale_fill_viridis_c(option="plasma", name=etq) + coord_equal() +
    labs(title=paste("Kriging ordinario —", v),
         subtitle=sprintf("%s | predicción | %d estaciones", format(fecha,"%d %b %Y"), n)) + bt
  p_sd <- ggplot(df, aes(X,Y,fill=sd)) + geom_raster() +
    geom_path(data=bordes, aes(X,Y,group=L1), color="white", linewidth=0.3, inherit.aes=FALSE) +
    geom_point(data=est, aes(X,Y), shape=21, size=2.2, fill="white", color="black", stroke=0.6, inherit.aes=FALSE) +
    scale_fill_viridis_c(option="viridis", direction=-1, name="sd kriging") + coord_equal() +
    labs(title=paste("Incertidumbre del kriging —", v),
         subtitle="Desviación típica (mayor lejos de estaciones)") + bt

  vn <- gsub("[^A-Za-z0-9]", "_", v)
  ggsave(file.path(DIR_FIG, sprintf("kriging_prediccion_%s.png", vn)), p_pred, width=8, height=6.5, dpi=200, bg="white")
  ggsave(file.path(DIR_FIG, sprintf("kriging_incertidumbre_%s.png", vn)), p_sd, width=8, height=6.5, dpi=200, bg="white")
  ggsave(file.path(DIR_FIG, sprintf("kriging_panel_%s.png", vn)), arrangeGrob(p_pred,p_sd,ncol=2), width=15, height=6.5, dpi=170, bg="white")

  cat(sprintf("  %-20s | %s | n=%2d | RMSE Kriging=%.3f | mejor: %s\n",
      v, format(fecha,"%Y-%m-%d"), n, rmses["Kriging"], names(which.min(rmses))))
}

tabla <- rbindlist(filas)
fwrite(tabla, file.path(DIR_TAB, "comparacion_kriging_1dia.csv"))
cat("\n================ TABLA COMPARATIVA (un día por variable) ================\n")
print(tabla)
cat("\nTabla:", file.path(DIR_TAB, "comparacion_kriging_1dia.csv"), "\n")
cat("Mapas en:", DIR_FIG, "\n")
