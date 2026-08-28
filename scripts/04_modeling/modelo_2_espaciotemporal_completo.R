# ==============================================================================
# MODELO 2 — ESPACIAL vs ESPACIO-TEMPORAL (comparación justa) — 2024-2025 DIARIO
# ==============================================================================
# Ajusta DOS modelos sobre los MISMOS datos, malla y periodo, cambiando SOLO la
# estructura espacio-temporal, para comparar de forma limpia:
#
#   Modelo A — campo espacial ESTÁTICO   : f(campo, model=spde)
#   Modelo B — campo ESPACIO-TEMPORAL    : f(campo, model=spde, group=día, ar1)
#                                          (la "otra forma", como Cameletti 2013)
#
# Reutiliza las anomalías del Modelo 1 (media mensual restada con los 3 años ->
# tendencia y estacionalidad ya quitadas). NO toca el Modelo 1.
#
# Periodo: 2024-2025 (2 años). Malla BASTA: imprescindible para que el campo
# espacio-temporal (nodos × días) quepa en memoria.
#
# Objetivo: ver si la interacción espacio-temporal MEJORA (DIC/WAIC/RMSE) y deja
# mejor ruido blanco que el campo estático. verbose=TRUE para seguir el avance.
# ==============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(INLA)
  library(ggplot2); library(gridExtra); library(grid)
})
set.seed(4827)

# ---- Configuración ----
PERIODO_INI <- as.Date("2024-01-01")   # 2 años: 2024-2025
FECHA_CORTE_TEST <- as.Date("2025-12-01")
# Malla basta (menos nodos = factible con mapa por día). Sube max.edge/cutoff
# si crashea; bájalos si quieres más detalle espacial (más lento).
MESH_MAXEDGE <- c(12, 24)
MESH_CUTOFF  <- 2

DIR_OUT <- here("outputs", "figures", "modelos", "modelo_2_espaciotemporal")
DIR_DAT <- here("data", "processed", "Modelos")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

guardar_tabla_png <- function(df, titulo, subtitulo = NULL, ruta) {
  df <- as.data.frame(df)
  th <- ttheme_minimal(
    core = list(fg_params = list(fontsize = 9),
                bg_params = list(fill = c("white", "#F5F5F5"), col = NA)),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"),
                   bg_params = list(fill = "#DDEEFF", col = NA)))
  comp <- arrangeGrob(textGrob(titulo, gp = gpar(fontsize = 13, fontface = "bold")),
    if (!is.null(subtitulo)) textGrob(subtitulo, gp = gpar(fontsize = 9, col = "grey40")) else nullGrob(),
    tableGrob(df, rows = NULL, theme = th), nrow = 3,
    heights = unit(c(0.5, 0.3, nrow(df) * 0.3 + 0.5), "in"))
  ggsave(ruta, comp, width = max(6, ncol(df) * 1.5), height = max(3, nrow(df) * 0.35 + 1.3),
         dpi = 150, bg = "white")
  cat("  PNG:", basename(ruta), "\n")
}

# ==============================================================================
# 1. Datos (anomalías del Modelo 1) restringidos a 2024-2025
# ==============================================================================
dt <- readRDS(file.path(DIR_DAT, "modelo_1_2023_2025_diario_preparado.rds"))
setDT(dt)
dt <- dt[FECHA >= PERIODO_INI]
setorder(dt, FECHA, ESTACION)

sel_f <- file.path(DIR_DAT, "modelo_1_paso2_seleccion.rds")
vars_cov <- if (file.exists(sel_f)) readRDS(sel_f)$vars_finales else
  c("temperatura", "humedad", "precipitacion", "presion", "radiacion", "viento", "trafico_int")

dt[, DIA_IDX := match(FECHA, sort(unique(FECHA)))]   # índice de día para el AR1
n_dias <- max(dt$DIA_IDX)
y_masked <- ifelse(dt$es_train, dt$anom_NO2, NA_real_)
idx_test <- which(!dt$es_train)
cat(sprintf("Filas: %d | días (grupos AR1): %d | estaciones: %d | test: %d filas\n",
            nrow(dt), n_dias, uniqueN(dt$ESTACION), length(idx_test)))
cat("Covariables:", paste(vars_cov, collapse = ", "), "\n")

# ==============================================================================
# 2. Malla BASTA + SPDE (priors PC débiles) — la MISMA para ambos modelos
# ==============================================================================
coords_est <- as.matrix(unique(dt[, .(X_km, Y_km)]))
bnd_i <- inla.nonconvex.hull(coords_est, convex = -0.08, resolution = 40)
bnd_o <- inla.nonconvex.hull(coords_est, convex = -0.25)
malla <- inla.mesh.2d(loc = coords_est, boundary = list(bnd_i, bnd_o),
                      max.edge = MESH_MAXEDGE, cutoff = MESH_CUTOFF)
spde <- inla.spde2.pcmatern(malla, alpha = 2,
                            prior.range = c(5, 0.5), prior.sigma = c(1, 0.5))
coords_all <- as.matrix(dt[, .(X_km, Y_km)])
FML_SEMANAL <- "f(dow, model='rw1', cyclic=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))"
cat(sprintf("Malla: %d nodos | campo estático: %d incógnitas | espacio-temporal: %d\n",
            malla$n, malla$n, malla$n * n_dias))

# ==============================================================================
# 3. Helpers de métricas y diagnóstico de residuos
# ==============================================================================
metricas <- function(mod, stk, etiqueta) {
  idd <- inla.stack.index(stk, tag = "est")$data
  pred <- mod$summary.fitted.values$mean[idd]; psd <- mod$summary.fitted.values$sd[idd]
  yt <- dt$anom_NO2[idx_test]; pt <- pred[idx_test]
  hp <- mod$summary.hyperpar
  fila <- grep("Precision for the Gaussian", rownames(hp))
  sdt <- sqrt(psd[idx_test]^2 + 1 / hp[fila, "mean"])
  gr <- function(pat) { i <- grep(pat, rownames(hp)); if (length(i)) hp[i[1], "mean"] else NA_real_ }
  data.table(Modelo = etiqueta, DIC = round(mod$dic$dic, 1), WAIC = round(mod$waic$waic, 1),
    RMSE = round(sqrt(mean((pt - yt)^2, na.rm = TRUE)), 4),
    Cob95 = round(100 * mean(yt >= pt - 1.96 * sdt & yt <= pt + 1.96 * sdt, na.rm = TRUE), 1),
    Rango_km = round(gr("^Range for campo"), 1),
    Sigma = round(gr("^Stdev for campo"), 3),
    AR1 = round(gr("GroupRho"), 3))
}
residuos_diag <- function(mod, stk, etiqueta, archivo, color) {
  idd <- inla.stack.index(stk, tag = "est")$data
  pred <- mod$summary.fitted.values$mean[idd]
  res_tr <- dt$anom_NO2[dt$es_train] - pred[dt$es_train]
  serie <- data.table(ID = dt[es_train == TRUE, DIA_IDX], r = res_tr)[
    , .(r = mean(r, na.rm = TRUE)), by = ID][order(ID)]
  png(file.path(DIR_OUT, archivo), width = 900, height = 480)
  acf(serie$r, lag.max = 30, main = sprintf("ACF residuos — %s", etiqueta), col = color, lwd = 2)
  dev.off()
  Box.test(serie$r, lag = 30, type = "Ljung-Box")$p.value
}

# ==============================================================================
# 4. MODELO A — campo espacial ESTÁTICO
# ==============================================================================
cat("\n== MODELO A: campo espacial ESTÁTICO ==\n")
idx_A <- inla.spde.make.index("campo", n.spde = spde$n.spde)
A_A   <- inla.spde.make.A(mesh = malla, loc = coords_all)
ef_A  <- c(list(dow = dt$dow), setNames(lapply(vars_cov, function(v) dt[[v]]), vars_cov))
stk_A <- inla.stack(tag = "est", data = list(y = y_masked), A = list(A_A, 1),
  effects = list(c(idx_A, list(intercept = 1)), ef_A), compress = FALSE)
fml_A <- as.formula(paste("y ~ 0 + intercept +", paste(vars_cov, collapse = " + "),
  "+", FML_SEMANAL, "+ f(campo, model = spde)"))
t0 <- Sys.time()
mod_A <- inla(fml_A, data = inla.stack.data(stk_A, spde = spde), family = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_A), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE),
  control.inla = list(strategy = "gaussian", int.strategy = "eb"),
  verbose = TRUE)
cat(sprintf("  Modelo A terminado en %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ==============================================================================
# 5. MODELO B — campo ESPACIO-TEMPORAL (interacción, AR1)
# ==============================================================================
cat("\n== MODELO B: campo ESPACIO-TEMPORAL (group = día, ar1) ==\n")
idx_B <- inla.spde.make.index("campo", n.spde = spde$n.spde, n.group = n_dias)
A_B   <- inla.spde.make.A(mesh = malla, loc = coords_all, group = dt$DIA_IDX, n.group = n_dias)
ef_B  <- c(list(dow = dt$dow), setNames(lapply(vars_cov, function(v) dt[[v]]), vars_cov))
stk_B <- inla.stack(tag = "est", data = list(y = y_masked), A = list(A_B, 1),
  effects = list(c(idx_B, list(intercept = 1)), ef_B), compress = FALSE)
fml_B <- as.formula(paste("y ~ 0 + intercept +", paste(vars_cov, collapse = " + "),
  "+", FML_SEMANAL,
  "+ f(campo, model = spde, group = campo.group, control.group = list(model = 'ar1'))"))
cat("  (puede tardar; sigue el avance con verbose)\n")
t0 <- Sys.time()
mod_B <- inla(fml_B, data = inla.stack.data(stk_B, spde = spde), family = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_B), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE),
  control.inla = list(strategy = "gaussian", int.strategy = "eb"),
  verbose = TRUE)
cat(sprintf("  Modelo B terminado en %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ==============================================================================
# 6. COMPARACIÓN + RESIDUOS
# ==============================================================================
comp <- rbind(metricas(mod_A, stk_A, "A · Espacial estático"),
              metricas(mod_B, stk_B, "B · Espacio-temporal (AR1)"))
guardar_tabla_png(comp, "Modelo 2 — Espacial vs Espacio-temporal",
  "Mismos datos, malla y periodo (2024-2025) | menor DIC/WAIC/RMSE = mejor",
  file.path(DIR_OUT, "comparacion.png"))
fwrite(comp, file.path(DIR_OUT, "comparacion.csv"))

pA <- residuos_diag(mod_A, stk_A, "Espacial estático", "acf_A_estatico.png", "#2166AC")
pB <- residuos_diag(mod_B, stk_B, "Espacio-temporal", "acf_B_espaciotemporal.png", "#1A9850")
guardar_tabla_png(
  data.table(Modelo = c("A · Estático", "B · Espacio-temporal"),
             LjungBox_p = signif(c(pA, pB), 4),
             Ruido_blanco = ifelse(c(pA, pB) >= 0.05, "Sí", "No")),
  "Modelo 2 — ¿Ruido blanco? (residuos)", NULL, file.path(DIR_OUT, "ljung_box.png"))

saveRDS(list(A = mod_A, B = mod_B), file.path(DIR_DAT, "modelo_2_comparacion.rds"))

cat("\n=== COMPARACIÓN ===\n"); print(comp)
cat(sprintf("Ljung-Box  A(estático)=%.4f  B(esp-temp)=%.4f\n", pA, pB))
cat("Guardado en:", DIR_OUT, "\n")
cat("--- Modelo 2 completado (NO se ha tocado el Modelo 1) ---\n")
