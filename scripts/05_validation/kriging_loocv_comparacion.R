# ==============================================================================
# COMPARACIÓN LOOCV: KRIGING ORDINARIO vs Media / 1-NN / kNN / IDW  (2025)
# ==============================================================================
# Calcula, sobre EXACTAMENTE los mismos datos (coordenadas normalizadas y
# estaciones co-localizadas promediadas), el RMSE y MAE por LOOCV de los 6
# métodos para Temperatura, Humedad relativa y Radiación solar.
#
# Escribe el resultado de CADA variable en cuanto termina, para no perder el
# progreso si la sesión se corta.
# Salida: outputs/tables/calidad_datos/comparacion_kriging.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(sf); library(gstat)
})

ANIO <- 2025L
K    <- 3L
VARIABLES <- c("Temperatura", "Humedad_Relativa", "Radiacion_Solar")
OUT <- here("outputs", "tables", "calidad_datos", "comparacion_kriging.csv")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

rmse <- function(r, p) sqrt(mean((r - p)^2))
mae  <- function(r, p) mean(abs(r - p))

# --- Datos + normalización de coordenadas (igual que el script oficial) --------
d <- readRDS(here("data","processed","Clima","diario",
                  sprintf("meteo_madrid_%d_diario.rds", ANIO)))
setDT(d)
col_rad <- grep("Solar$", names(d), value = TRUE)
if (length(col_rad) == 1L) setnames(d, col_rad, "Radiacion_Solar")
fl <- !is.na(d$LATITUD)  & d$LATITUD  < 35
fo <- !is.na(d$LONGITUD) & d$LONGITUD > -1
d[fl, LATITUD  := LATITUD * 10^round(log10(40.4 / abs(LATITUD)))]
d[fo, LONGITUD := LONGITUD * 10]
d <- d[LATITUD > 40.2 & LATITUD < 40.7 & LONGITUD > -4 & LONGITUD < -3.3]

resultados <- list()

for (v in VARIABLES) {
  cat("\n===== ", v, " =====\n", sep = "")
  dv <- d[!is.na(get(v)), .(ESTACION, LONGITUD, LATITUD, FECHA, VALOR = get(v))]
  fechas <- sort(unique(dv$FECHA))

  # Acumuladores de (real, pred) por método
  ac <- list(real=c(), Media=c(), NN=c(), kNN=c(), IDWb1=c(), IDWb2=c(), Krig=c())
  krig_real <- c(); krig_pred <- c()
  n_dias <- 0L; t0 <- Sys.time()

  for (f in fechas) {
    dd <- dv[FECHA == f]
    # Proyección a UTM (metros)
    sfp <- st_transform(st_as_sf(dd, coords = c("LONGITUD","LATITUD"), crs = 4326), 25830)
    xy  <- st_coordinates(sfp)
    val <- sfp$VALOR
    # Colapsar co-localizadas: promediar valor, un punto por coordenada
    key <- do.call(paste, as.data.frame(round(xy)))
    val <- ave(val, key, FUN = function(x) mean(x, na.rm = TRUE))
    keep <- !duplicated(key)
    xy <- xy[keep, , drop = FALSE]; val <- val[keep]
    n <- length(val)
    if (n < 5L) next
    n_dias <- n_dias + 1L

    # --- Métodos manuales (k vecinos) ---
    D <- as.matrix(dist(xy)); diag(D) <- Inf
    pred_media <- (sum(val) - val) / (n - 1)
    pred_nn <- pred_knn <- pred_b1 <- pred_b2 <- numeric(n)
    for (i in seq_len(n)) {
      idx <- order(D[i, ])[seq_len(K)]
      vv <- val[idx]; dv2 <- D[i, idx]
      pred_nn[i]  <- vv[1]
      pred_knn[i] <- mean(vv)
      w1 <- 1 / pmax(dv2, 1)^1; pred_b1[i] <- sum(w1 * vv) / sum(w1)
      w2 <- 1 / pmax(dv2, 1)^2; pred_b2[i] <- sum(w2 * vv) / sum(w2)
    }
    ac$real  <- c(ac$real,  val)
    ac$Media <- c(ac$Media, pred_media)
    ac$NN    <- c(ac$NN,    pred_nn)
    ac$kNN   <- c(ac$kNN,   pred_knn)
    ac$IDWb1 <- c(ac$IDWb1, pred_b1)
    ac$IDWb2 <- c(ac$IDWb2, pred_b2)

    # --- Kriging ordinario (variograma diario + LOOCV) ---
    sfd <- st_as_sf(data.frame(x = xy[,1], y = xy[,2], VALOR = val),
                    coords = c("x","y"), crs = 25830)
    vemp <- tryCatch(variogram(VALOR ~ 1, sfd), error = function(e) NULL)
    if (is.null(vemp) || nrow(vemp) < 3) next
    ps <- var(val); rg <- as.numeric(median(dist(xy)))
    m <- tryCatch(suppressWarnings(fit.variogram(vemp, vgm(ps,"Sph",rg,0.1*ps))),
                  error = function(e) NULL)
    if (is.null(m) || any(m$psill < 0)) m <- vgm(ps,"Exp",rg,0.1*ps)
    cv <- tryCatch(suppressWarnings(krige.cv(VALOR ~ 1, sfd, model = m, verbose = FALSE)),
                   error = function(e) NULL)
    if (is.null(cv)) next
    krig_real <- c(krig_real, cv$observed)
    krig_pred <- c(krig_pred, cv$var1.pred)
  }

  fila <- data.table(
    Variable = v,
    N_Dias = n_dias,
    N_Predicciones = length(ac$real),
    RMSE_Media = round(rmse(ac$real, ac$Media), 3),
    RMSE_NN    = round(rmse(ac$real, ac$NN), 3),
    RMSE_kNN   = round(rmse(ac$real, ac$kNN), 3),
    RMSE_IDW_b1= round(rmse(ac$real, ac$IDWb1), 3),
    RMSE_IDW_b2= round(rmse(ac$real, ac$IDWb2), 3),
    RMSE_Kriging = round(rmse(krig_real, krig_pred), 3),
    MAE_Media = round(mae(ac$real, ac$Media), 3),
    MAE_kNN   = round(mae(ac$real, ac$kNN), 3),
    MAE_IDW_b1= round(mae(ac$real, ac$IDWb1), 3),
    MAE_Kriging = round(mae(krig_real, krig_pred), 3)
  )
  # Mejor método por RMSE
  rmses <- unlist(fila[, .(RMSE_Media, RMSE_NN, RMSE_kNN, RMSE_IDW_b1, RMSE_IDW_b2, RMSE_Kriging)])
  fila[, Mejor_Metodo := c("Media","1-NN","kNN","IDW b=1","IDW b=2","Kriging")[which.min(rmses)]]

  resultados[[v]] <- fila
  # Guardado incremental (persiste aunque se corte)
  fwrite(rbindlist(resultados), OUT)
  cat(sprintf("  %s | dias=%d | pred=%d | RMSE Kriging=%.3f (mejor: %s) | %.0f s\n",
      v, n_dias, length(ac$real), fila$RMSE_Kriging, fila$Mejor_Metodo,
      as.numeric(difftime(Sys.time(), t0, units="secs"))))
  print(fila)
}

cat("\n================ TABLA FINAL ================\n")
print(rbindlist(resultados))
cat("\nGuardado en:", OUT, "\n")
