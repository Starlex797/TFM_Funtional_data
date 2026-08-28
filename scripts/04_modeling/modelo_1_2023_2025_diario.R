# ==============================================================================
# MODELO 1 — 2023-2025 DIARIO  (INLA-SPDE, bayesiano jerárquico)
# ==============================================================================
# Estrategia:
#   - Quitar ESTACIONALIDAD y TENDENCIA con ANOMALÍAS: a la respuesta y a cada
#     covariable se les resta la media mensual GLOBAL (por mes-año, toda la
#     ciudad), conservando las diferencias entre estaciones para el campo espacial.
#   - Mantener el CICLO SEMANAL (la media mensual no lo elimina): se modela con
#     f(dow, model="rw1", cyclic=TRUE) en el modelo.
#   - Modelo primero SÓLO ESPACIAL (SPDE Matérn); diagnosticar residuos (ACF /
#     Ljung-Box). Objetivo: RUIDO BLANCO. Si queda autocorrelación temporal ->
#     efecto temporal compartido f(ID_TIEMPO, "ar1"/"rw1").
#   - Covariables colineales se descartan por VIF (Paso 2).
#
# Desalineamiento: las covariables ya vienen INTERPOLADAS por IDW a las 24
# estaciones de NO2 (columnas *_raw del maestro). No hay que re-alinear.
#
# PASOS:  0-1 Preparación + anomalías  <-- ESTE BLOQUE
#         2   Selección (VIF/signif/DIC) | 3 Espacial SPDE | 4 Diagnóstico |
#         5   AR(1)/temporal si procede  | 6 Comparación
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(ggplot2)
})
set.seed(4827)

# --- Utilidad de verificación ---
verif <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "OK " else "FALLO", msg))
  if (!isTRUE(cond)) warning("Verificación fallida: ", msg)
  invisible(cond)
}

# ==============================================================================
# CONFIGURACIÓN
# ==============================================================================
ANIOS <- 2023:2025
FECHA_CORTE_TEST <- as.Date("2025-12-01")   # test = diciembre 2025

COVS_RAW <- c("Temperatura_raw", "Humedad_Relativa_raw", "Precipitaciones_raw",
              "Presion Barométrica_raw", "Radiación Solar_raw",
              "Velocidad Viento_raw", "intensidad_raw", "carga_raw")
COVS_ALIAS <- c("temperatura", "humedad", "precipitacion",
                "presion", "radiacion", "viento",
                "trafico_int", "trafico_carga")

DIR_OUT <- here("outputs", "figures", "modelos", "modelo_1_2023_2025_diario")
DIR_DAT <- here("data", "processed", "Modelos")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_DAT, recursive = TRUE, showWarnings = FALSE)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("MODELO 1 — 2023-2025 DIARIO | PASO 0-1: preparación + anomalías\n")
cat(strrep("=", 72), "\n", sep = "")

# ==============================================================================
# PASO 0 — CARGA Y COMBINACIÓN DE LOS MAESTROS DIARIOS 2023-2025
# ==============================================================================
COLS_CARGAR <- c("ESTACION", "FECHA", "LONGITUD", "LATITUD",
                 "LOG_NO2_DIARIO", "DATO_DIARIO", COVS_RAW)
dt <- rbindlist(lapply(ANIOS, function(a) {
  f <- here("data", "processed", "Maestro", "diario",
            sprintf("dataset_maestro_inla_%d_DIARIO.rds", a))
  d <- readRDS(f); setDT(d)
  d[, ..COLS_CARGAR]
}), use.names = TRUE)
setnames(dt, "LOG_NO2_DIARIO", "LOG_NO2")
dt[, FECHA := as.Date(FECHA)]
dt <- dt[!is.na(LOG_NO2)]
# Quitar filas con alguna covariable cruda NA ANTES de calcular anomalías
# (así las medias mensuales y la estandarización se computan sobre datos limpios)
n_na_raw <- sum(rowSums(is.na(dt[, ..COVS_RAW])) > 0)
if (n_na_raw > 0) {
  cat(sprintf("Filas eliminadas por covariable cruda NA: %d\n", n_na_raw))
  dt <- dt[rowSums(is.na(dt[, ..COVS_RAW])) == 0]
}
setorder(dt, FECHA, ESTACION)

cat(sprintf("\nFilas: %d | estaciones: %d | días: %d | rango: %s a %s\n",
            nrow(dt), uniqueN(dt$ESTACION), uniqueN(dt$FECHA),
            as.character(min(dt$FECHA)), as.character(max(dt$FECHA))))

# --- Coordenadas UTM en km ---
coords <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
sf_c <- st_transform(st_as_sf(coords, coords = c("LONGITUD", "LATITUD"),
                              crs = 4326), 25830)
coords[, `:=`(X_km = st_coordinates(sf_c)[, 1] / 1000,
              Y_km = st_coordinates(sf_c)[, 2] / 1000)]
dt <- merge(dt, coords[, .(ESTACION, X_km, Y_km)], by = "ESTACION")
setorder(dt, FECHA, ESTACION)

# --- Índice temporal + día de la semana + train/test ---
fechas <- sort(unique(dt$FECHA))
dt[, ID_TIEMPO := match(FECHA, fechas)]
dt[, dow := as.integer(format(FECHA, "%u"))]     # 1=lunes ... 7=domingo
dt[, finde := as.integer(dow >= 6)]
dt[, MES_ANIO := format(FECHA, "%Y-%m")]
dt[, es_train := FECHA < FECHA_CORTE_TEST]

# ------------------ VERIFICACIÓN PASO 0 ------------------
cat("\n--- Verificación Paso 0 ---\n")
verif(uniqueN(dt$ESTACION) == 24, "24 estaciones de NO2")
verif(!any(duplicated(dt[, .(ESTACION, FECHA)])), "sin duplicados estación-día")
verif(all(!is.na(dt$X_km)) && all(!is.na(dt$Y_km)), "coordenadas UTM completas")
verif(all(dt$dow %in% 1:7), "día de la semana en 1..7")
verif(sum(dt$es_train) > 0 && sum(!dt$es_train) > 0, "train y test no vacíos")
cat(sprintf("  Train: %d filas (hasta %s) | Test: %d filas (%s a %s)\n",
            sum(dt$es_train), as.character(max(dt[es_train == TRUE, FECHA])),
            sum(!dt$es_train), as.character(FECHA_CORTE_TEST),
            as.character(max(dt$FECHA))))

# ==============================================================================
# PASO 1 — ANOMALÍAS: restar la media mensual GLOBAL (por mes-año, toda la ciudad)
# ==============================================================================
# Quita estacionalidad Y tendencia de una vez (el nivel medio de la ciudad cada
# mes-año las contiene). GLOBAL (no por estación) para CONSERVAR las diferencias
# de nivel entre estaciones -> así el campo espacial estático SÍ tiene señal que
# capturar. Se aplica a la respuesta y a cada covariable.

dt[, anom_NO2 := LOG_NO2 - mean(LOG_NO2, na.rm = TRUE), by = .(MES_ANIO)]
for (i in seq_along(COVS_RAW)) {
  v <- COVS_RAW[i]; a <- paste0("anom_", COVS_ALIAS[i])
  dt[, (a) := get(v) - mean(get(v), na.rm = TRUE), by = .(MES_ANIO)]
}

# --- Estandarizar covariables-anomalía (estadísticos de TRAIN) ---
for (alias in COVS_ALIAS) {
  a <- paste0("anom_", alias)
  mu <- mean(dt[es_train == TRUE][[a]], na.rm = TRUE)
  sg <- sd(dt[es_train == TRUE][[a]], na.rm = TRUE)
  dt[, (alias) := (get(a) - mu) / sg]
}

# ------------------ VERIFICACIÓN PASO 1 ------------------
cat("\n--- Verificación Paso 1 (anomalías) ---\n")
# (a) La media GLOBAL por mes-año de la anomalía debe ser ~0 (quitamos el nivel mensual)
chk_mean <- dt[, .(m = mean(anom_NO2, na.rm = TRUE)), by = .(MES_ANIO)]
verif(max(abs(chk_mean$m), na.rm = TRUE) < 1e-8,
      sprintf("anom_NO2 con media ~0 por mes-año global (máx |media|=%.1e)",
              max(abs(chk_mean$m), na.rm = TRUE)))
# (b) Covariables estandarizadas: media~0, sd~1 en train
for (alias in COVS_ALIAS) {
  mu <- mean(dt[es_train == TRUE][[alias]], na.rm = TRUE)
  sg <- sd(dt[es_train == TRUE][[alias]], na.rm = TRUE)
  verif(abs(mu) < 1e-6 && abs(sg - 1) < 1e-6,
        sprintf("%-13s estandarizada (media=%.0e, sd=%.4f)", alias, mu, sg))
}
# (c) No se han introducido NAs nuevos en las covariables del modelo
n_na_cov <- sum(sapply(COVS_ALIAS, function(a) sum(is.na(dt[[a]]))))
verif(n_na_cov == 0, sprintf("sin NAs en covariables-anomalía (NAs=%d)", n_na_cov))

# ------------------ EVIDENCIA: ciclo semanal ------------------
cat("\n--- Evidencia del ciclo semanal (NO2 medio por día de la semana) ---\n")
dias_lab <- c("Lun", "Mar", "Mié", "Jue", "Vie", "Sáb", "Dom")
perfil_dow <- dt[, .(NO2_ugm3 = mean(DATO_DIARIO, na.rm = TRUE),
                     anom = mean(anom_NO2, na.rm = TRUE), n = .N), by = dow][order(dow)]
perfil_dow[, dia := factor(dias_lab[dow], levels = dias_lab)]
print(perfil_dow[, .(dia, NO2_ugm3 = round(NO2_ugm3, 1), anom = round(anom, 4), n)])
fwrite(perfil_dow, file.path(DIR_OUT, "01_perfil_dia_semana.csv"))

# ------------------ CORRELACIÓN entre covariables (colinealidad) ------------------
cat("\n--- Matriz de correlación de covariables-anomalía (para VIF del Paso 2) ---\n")
mat_cor <- cor(dt[es_train == TRUE, ..COVS_ALIAS], use = "pairwise.complete.obs")
print(round(mat_cor, 2))
cat(sprintf("  Pares con |cor| > 0.7: %d (candidatos a descarte por VIF)\n",
            (sum(abs(mat_cor) > 0.7) - length(COVS_ALIAS)) / 2))

# ==============================================================================
# GRÁFICOS: antes vs después (se elimina estacionalidad+tendencia) y ciclo semanal
# ==============================================================================
serie <- dt[, .(original = mean(LOG_NO2, na.rm = TRUE),
                anomalia = mean(anom_NO2, na.rm = TRUE)), by = FECHA]
serie_l <- melt(serie, id.vars = "FECHA", variable.name = "tipo", value.name = "valor")
g1 <- ggplot(serie_l, aes(FECHA, valor)) +
  geom_line(color = "grey65", linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, color = "#2166AC") +
  facet_wrap(~ tipo, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(original = "LOG_NO2 original (con estacionalidad + tendencia)",
                                      anomalia = "Anomalía (media mensual restada → debe ser plana)"))) +
  labs(title = "Efecto de restar la media mensual", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_OUT, "00_antes_vs_despues.png"), g1, width = 11, height = 6, dpi = 150)

g2 <- ggplot(perfil_dow, aes(dia, anom)) +
  geom_col(fill = "#2166AC", width = 0.65) +
  geom_hline(yintercept = 0, color = "grey40") +
  labs(title = "Ciclo semanal en la anomalía de NO2",
       subtitle = "La media mensual NO elimina el ciclo semanal → se modelará con f(dow, cyclic=TRUE)",
       x = NULL, y = "Anomalía media de LOG_NO2") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_OUT, "01_ciclo_semanal.png"), g2, width = 8, height = 4.5, dpi = 150)

# ==============================================================================
# GUARDAR DATASET PREPARADO
# ==============================================================================
saveRDS(dt, file.path(DIR_DAT, "modelo_1_2023_2025_diario_preparado.rds"))
cat("\n--- PASO 0-1 completado ---\n")
cat("Dataset preparado -> data/processed/Modelos/modelo_1_2023_2025_diario_preparado.rds\n")
cat("Figuras -> outputs/figures/modelos/modelo_1_2023_2025_diario/\n")

# ==============================================================================
# LIBRERÍAS DEL MODELO + HELPERS (tablas PNG, ajuste INLA, métricas)
# ==============================================================================
suppressPackageStartupMessages({
  library(INLA); library(car); library(gridExtra); library(grid)
})
if (!exists("dt")) dt <- readRDS(file.path(DIR_DAT, "modelo_1_2023_2025_diario_preparado.rds"))

# --- Guardar una tabla como PNG (justificación gráfica de cada decisión) ---
guardar_tabla_png <- function(df, titulo, subtitulo = NULL, ruta) {
  df <- as.data.frame(df)
  th <- ttheme_minimal(
    core = list(fg_params = list(fontsize = 9),
                bg_params = list(fill = c("white", "#F5F5F5"), col = NA)),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"),
                   bg_params = list(fill = "#DDEEFF", col = NA)))
  tg <- tableGrob(df, rows = NULL, theme = th)
  tit <- textGrob(titulo, gp = gpar(fontsize = 13, fontface = "bold"))
  comp <- if (!is.null(subtitulo)) {
    sg <- textGrob(subtitulo, gp = gpar(fontsize = 9, col = "grey40"))
    arrangeGrob(tit, sg, tg, nrow = 3,
                heights = unit(c(0.5, 0.35, nrow(df) * 0.3 + 0.5), "in"))
  } else {
    arrangeGrob(tit, tg, nrow = 2,
                heights = unit(c(0.5, nrow(df) * 0.3 + 0.5), "in"))
  }
  ggsave(ruta, comp, width = max(6, ncol(df) * 1.5),
         height = max(3, nrow(df) * 0.35 + 1.2), dpi = 150, bg = "white")
  cat("  PNG:", basename(ruta), "\n")
}

# ==============================================================================
# PRELIMINARES DEL MODELO: malla 8 km, SPDE (priors PC), matrices, helpers
# ==============================================================================
dt_tr <- dt[es_train == TRUE]

# Malla ÚNICA de 8 km (adecuada a solo 24 estaciones)
coords_est <- as.matrix(unique(dt[, .(X_km, Y_km)]))
bnd_i <- inla.nonconvex.hull(coords_est, convex = -0.05, resolution = 50)
bnd_o <- inla.nonconvex.hull(coords_est, convex = -0.2)
malla <- inla.mesh.2d(loc = coords_est, boundary = list(bnd_i, bnd_o),
                      max.edge = c(8, 12), cutoff = 0.5)

# SPDE con priors PC DÉBILMENTE informativos (que NO maten el campo):
#   P(rango < 5 km) = 0.5  -> mediana ~5 km (escala inter-estación), sin forzar
#                             rangos enormes
#   P(sigma > 1)    = 0.5  -> permite una SD espacial moderada; NO la aplasta a 0
# (Los priors anteriores c(5,0.1)/c(1,0.01) apagaban el campo: SD~0, rango enorme.)
crear_spde <- function(mesh) inla.spde2.pcmatern(
  mesh, alpha = 2, prior.range = c(5, 0.5), prior.sigma = c(1, 0.5))
spde <- crear_spde(malla)
idx  <- inla.spde.make.index("campo_espacial", n.spde = spde$n.spde)

y_masked <- ifelse(dt$es_train, dt$anom_NO2, NA_real_)   # test enmascarado (finales)
idx_test <- which(!dt$es_train)

# Subconjunto REPRESENTATIVO para la SELECCIÓN (por memoria): 24 estaciones x
# ~1 de cada 15 días de train (repartidos por estación del año y por año).
# La selección de covariables es estable y no necesita las ~24.800 filas; los
# modelos FINALES (Paso 3/5) sí se ajustan sobre TODOS los datos. Es la misma
# estrategia que simulacion.R (que reduce el diario "para memoria").
fechas_tr  <- sort(unique(dt[es_train == TRUE, FECHA]))
fechas_sel <- fechas_tr[seq(1, length(fechas_tr), by = 15)]
dt_sel <- dt[es_train == TRUE & FECHA %in% fechas_sel]

# Término del ciclo semanal (prior PC en su precisión)
FML_SEMANAL <- paste0("f(dow, model = 'rw1', cyclic = TRUE, ",
  "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))")

# --- Ajustar un modelo INLA-SPDE sobre un conjunto de datos 'd' ---
# d: data.table con las columnas del modelo | y_vec: respuesta (puede llevar NA)
# compute_fitted = TRUE solo en los modelos finales (necesitan los ajustados);
# en la SELECCIÓN se deja FALSE (solo DIC) para no agotar la memoria.
ajustar_spde <- function(d, y_vec, vars_cov,
                         temporal = c("ninguno", "ar1", "rw1"),
                         compute_fitted = FALSE, int_strategy = "eb",
                         verbose_inla = FALSE) {
  temporal <- match.arg(temporal)
  A_d <- inla.spde.make.A(mesh = malla, loc = as.matrix(d[, .(X_km, Y_km)]))
  ef2 <- c(list(dow = d$dow),
           setNames(lapply(vars_cov, function(v) d[[v]]), vars_cov))
  if (temporal != "ninguno") ef2$tiempo <- d$ID_TIEMPO
  stk <- inla.stack(tag = "est", data = list(y = y_vec), A = list(A_d, 1),
    effects = list(c(idx, list(intercept = 1)), ef2), compress = FALSE)
  ter <- c("0", "intercept", vars_cov, FML_SEMANAL,
           "f(campo_espacial, model = spde)")
  if (temporal == "ar1") ter <- c(ter, "f(tiempo, model = 'ar1')")
  if (temporal == "rw1") ter <- c(ter, "f(tiempo, model = 'rw1')")
  fml <- as.formula(paste("y ~", paste(ter, collapse = " + ")))
  mod <- inla(fml, data = inla.stack.data(stk, spde = spde), family = "gaussian",
    control.predictor = list(A = inla.stack.A(stk), compute = compute_fitted),
    control.compute = list(dic = TRUE, waic = compute_fitted, cpo = FALSE),
    control.inla = list(strategy = "gaussian", int.strategy = int_strategy),
    verbose = verbose_inla)
  list(mod = mod, stk = stk)
}

# --- Métricas de un modelo (DIC/WAIC + RMSE/MAE/cobertura en test) ---
metricas_modelo <- function(res, etiqueta) {
  mod <- res$mod
  idd <- inla.stack.index(res$stk, tag = "est")$data
  pred <- mod$summary.fitted.values$mean[idd]
  psd  <- mod$summary.fitted.values$sd[idd]
  yt <- dt$anom_NO2[idx_test]; pt <- pred[idx_test]
  fila <- grep("Precision for the Gaussian", rownames(mod$summary.hyperpar), value = TRUE)
  var_obs <- 1 / mod$summary.hyperpar[fila[1], "mean"]
  sdt <- sqrt(psd[idx_test]^2 + var_obs)
  data.table(Modelo = etiqueta, DIC = round(mod$dic$dic, 1),
    WAIC = if (!is.null(mod$waic)) round(mod$waic$waic, 1) else NA_real_,
    RMSE = round(sqrt(mean((pt - yt)^2, na.rm = TRUE)), 4),
    MAE = round(mean(abs(pt - yt), na.rm = TRUE), 4),
    Cob95 = round(100 * mean(yt >= pt - 1.96 * sdt & yt <= pt + 1.96 * sdt, na.rm = TRUE), 1))
}

# --- Residuos de entrenamiento (real - ajustado) ---
residuos_train <- function(res) {
  idd <- inla.stack.index(res$stk, tag = "est")$data
  pred <- res$mod$summary.fitted.values$mean[idd]
  dt$anom_NO2[dt$es_train] - pred[dt$es_train]
}

# ==============================================================================
# PASO 2 — SELECCIÓN DE COVARIABLES (cada decisión justificada en PNG)
# ==============================================================================
# Reglas (como en simulacion.R): NO se meten al modelo covariables que estén
# correlacionadas entre sí ni que no sean significativas.
#   2a) Correlación   -> descartar una de cada par muy correlacionado (|r|>0.8)
#   2b) VIF (<5)      -> descartar colinealidad residual
#   2c) Significancia -> descartar las que su IC 95% incluya el 0
#   2d) Stepwise DIC  -> conjunto final (mejora de DIC >= 2)
# ==============================================================================
cat("\n", strrep("=", 72), "\n", sep = "")
cat("PASO 2 — Selección de covariables (justificación PNG por paso)\n")
cat(strrep("=", 72), "\n", sep = "")

# --- 2a. Correlación: descartar una de cada par muy correlacionado ---
CORR_UMBRAL <- 0.8
mat_cor <- cor(dt_tr[, ..COVS_ALIAS], use = "pairwise.complete.obs")
# PNG del heatmap de correlación (justificación gráfica)
cor_long <- melt(data.table(v1 = rownames(mat_cor), mat_cor),
                 id.vars = "v1", variable.name = "v2", value.name = "r")
g_cor <- ggplot(cor_long, aes(v1, v2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlación entre covariables-anomalía",
       subtitle = sprintf("Se descarta una de cada par con |r| > %.1f", CORR_UMBRAL),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_OUT, "02a_correlacion.png"), g_cor, width = 8, height = 6.5, dpi = 150)

vars_cor <- COVS_ALIAS
descartes_cor <- character(0)
pares <- which(abs(mat_cor) > CORR_UMBRAL & upper.tri(mat_cor), arr.ind = TRUE)
for (k in seq_len(nrow(pares))) {
  v1 <- rownames(mat_cor)[pares[k, 1]]; v2 <- colnames(mat_cor)[pares[k, 2]]
  # se conserva la primera (más interpretable) y se descarta la segunda
  if (v1 %in% vars_cor && v2 %in% vars_cor) {
    vars_cor <- setdiff(vars_cor, v2); descartes_cor <- c(descartes_cor, v2)
    cat(sprintf("  Correlación: '%s' ~ '%s' (r=%.2f) -> descarto '%s'\n",
                v1, v2, mat_cor[v1, v2], v2))
  }
}
guardar_tabla_png(
  data.table(Par = if (length(descartes_cor)) "trafico_int ~ trafico_carga" else "-",
             r = if (length(descartes_cor)) round(mat_cor["trafico_int", "trafico_carga"], 2) else NA,
             Descartada = if (length(descartes_cor)) paste(descartes_cor, collapse = ", ") else "ninguna"),
  "Paso 2a — Descartes por correlación",
  sprintf("Umbral |r| > %.1f | conservadas: %s", CORR_UMBRAL, paste(vars_cor, collapse = ", ")),
  file.path(DIR_OUT, "02a_tabla_correlacion.png"))

# --- 2b. VIF hacia atrás (umbral 5) sobre las supervivientes a correlación ---
VIF_UMBRAL <- 5
vars_vif <- vars_cor
repeat {
  if (length(vars_vif) < 2) break
  vifs <- car::vif(lm(as.formula(paste("anom_NO2 ~",
              paste(vars_vif, collapse = " + "))), data = dt_tr))
  if (max(vifs) <= VIF_UMBRAL) break
  drop <- names(which.max(vifs))
  cat(sprintf("  VIF: elimino '%s' (VIF=%.2f)\n", drop, max(vifs)))
  vars_vif <- setdiff(vars_vif, drop)
}
vif_final <- car::vif(lm(as.formula(paste("anom_NO2 ~",
                paste(vars_vif, collapse = " + "))), data = dt_tr))
guardar_tabla_png(
  data.table(Variable = names(vif_final), VIF = round(as.numeric(vif_final), 3))[order(-VIF)],
  "Paso 2b — VIF (colinealidad)",
  sprintf("Umbral VIF = %g | %d variables retenidas", VIF_UMBRAL, length(vars_vif)),
  file.path(DIR_OUT, "02b_tabla_vif.png"))

# --- 2c. Significancia bayesiana (INLA-SPDE, sobre subconjunto) ---
cat(sprintf("  Selección sobre subconjunto: %d filas (%d días × %d est.)\n",
            nrow(dt_sel), uniqueN(dt_sel$FECHA), uniqueN(dt_sel$ESTACION)))
cat("  Ajustando modelo de significancia (INLA-SPDE 8 km)...\n")
res_sig <- ajustar_spde(dt_sel, dt_sel$anom_NO2, vars_vif)
sf <- res_sig$mod$summary.fixed[vars_vif, , drop = FALSE]
sig <- sf[, "0.025quant"] > 0 | sf[, "0.975quant"] < 0
tabla_sig <- data.table(
  Variable = rownames(sf), Media = round(sf$mean, 4),
  Q2.5 = round(sf[, "0.025quant"], 4), Q97.5 = round(sf[, "0.975quant"], 4),
  Significativa = ifelse(sig, "Sí", "No"))
guardar_tabla_png(tabla_sig, "Paso 2c — Significancia bayesiana (IC 95%)",
  sprintf("%d significativas de %d (se descartan las que incluyen el 0)",
          sum(sig), length(sig)),
  file.path(DIR_OUT, "02c_tabla_significancia.png"))
vars_sig <- rownames(sf)[sig]

# --- 2d. Stepwise DIC (conjunto final: mejora de DIC >= 2) ---
DIC_MIN <- 2
dic_sel <- function(vc) {  # DIC de un modelo en el subconjunto (+ liberar memoria)
  d <- ajustar_spde(dt_sel, dt_sel$anom_NO2, vc)$mod$dic$dic
  gc(verbose = FALSE); d
}
vars_sw <- if (length(vars_sig)) vars_sig else vars_vif
dic_actual <- dic_sel(vars_sw)
hist_sw <- list(data.table(Iter = 0L, Accion = "inicial", Variable = "-",
                           Modelo = paste(vars_sw, collapse = "+"), DIC = round(dic_actual, 1)))
for (it in seq_len(length(vars_vif) + 1L)) {
  cand <- data.table()
  for (v in setdiff(vars_vif, vars_sw))
    cand <- rbind(cand, data.table(accion = "añadir", v = v, DIC = dic_sel(c(vars_sw, v))))
  if (length(vars_sw) >= 2) for (v in vars_sw)
    cand <- rbind(cand, data.table(accion = "quitar", v = v, DIC = dic_sel(setdiff(vars_sw, v))))
  if (!nrow(cand)) break
  mejor <- cand[which.min(DIC)]
  if (dic_actual - mejor$DIC < DIC_MIN) break
  vars_sw <- if (mejor$accion == "añadir") c(vars_sw, mejor$v) else setdiff(vars_sw, mejor$v)
  dic_actual <- mejor$DIC
  hist_sw[[length(hist_sw) + 1]] <- data.table(Iter = it, Accion = mejor$accion,
    Variable = mejor$v, Modelo = paste(vars_sw, collapse = "+"), DIC = round(dic_actual, 1))
}
guardar_tabla_png(rbindlist(hist_sw), "Paso 2d — Stepwise DIC",
  sprintf("Mejora mínima de DIC = %g | finales: %s", DIC_MIN, paste(vars_sw, collapse = ", ")),
  file.path(DIR_OUT, "02d_tabla_stepwise_dic.png"))

vars_finales <- vars_sw
saveRDS(list(vars_cor = vars_cor, vars_vif = vars_vif, vars_sig = vars_sig,
             vars_finales = vars_finales, mat_cor = mat_cor),
        file.path(DIR_DAT, "modelo_1_paso2_seleccion.rds"))

# ------------------ VERIFICACIÓN PASO 2 ------------------
cat("\n--- Verificación Paso 2 ---\n")
# Tras el descarte por correlación, ningún par debe superar el umbral
mc2 <- cor(dt_tr[, ..vars_cor], use = "pairwise.complete.obs")
verif(max(abs(mc2[upper.tri(mc2)])) <= CORR_UMBRAL,
      sprintf("sin pares |r| > %.1f tras el descarte (máx=%.2f)",
              CORR_UMBRAL, max(abs(mc2[upper.tri(mc2)]))))
verif(max(vif_final) <= VIF_UMBRAL, sprintf("VIF final ≤ %g (máx=%.2f)", VIF_UMBRAL, max(vif_final)))
verif(length(vars_finales) >= 1, sprintf("hay covariables finales (%d)", length(vars_finales)))
cat(sprintf("\n>> FINALES tras correlación+VIF+signif+DIC: %s\n", paste(vars_finales, collapse = ", ")))
cat("--- PASO 2 completado ---\n")

# ==============================================================================
# PASO 3 — MODELO SÓLO ESPACIAL (SPDE) + ciclo semanal
# ==============================================================================
cat("\n", strrep("=", 72), "\n", sep = "")
cat("PASO 3 — Modelo espacial SPDE\n")
cat(strrep("=", 72), "\n", sep = "")

res_esp <- ajustar_spde(dt, y_masked, vars_finales, "ninguno",
                        compute_fitted = TRUE, int_strategy = "ccd",
                        verbose_inla = TRUE)
mod_esp <- res_esp$mod

# Tabla de efectos fijos (PNG)
ef_fix <- mod_esp$summary.fixed
tab_fix <- data.table(Efecto = rownames(ef_fix), Media = round(ef_fix$mean, 4),
  Q2.5 = round(ef_fix[, "0.025quant"], 4), Q97.5 = round(ef_fix[, "0.975quant"], 4))
guardar_tabla_png(tab_fix, "Paso 3 — Efectos fijos (modelo espacial)", NULL,
  file.path(DIR_OUT, "03_efectos_fijos.png"))

# Tabla de hiperparámetros del campo (rango, sigma) (PNG)
spde_res <- inla.spde.result(mod_esp, "campo_espacial", spde)
tab_hyp <- data.table(
  Parametro = c("Rango (km)", "Sigma (sd espacial)", "SD del ruido (obs.)"),
  Media = round(c(exp(spde_res$summary.log.range.nominal$mean),
                  exp(spde_res$summary.log.variance.nominal$mean)^0.5,
                  sqrt(1 / mod_esp$summary.hyperpar[
                    grep("Precision for the Gaussian", rownames(mod_esp$summary.hyperpar)), "mean"])), 3))
guardar_tabla_png(tab_hyp, "Paso 3 — Hiperparámetros del campo Matérn", NULL,
  file.path(DIR_OUT, "03_hiperparametros.png"))

# Métricas (PNG)
met_esp <- metricas_modelo(res_esp, "Espacial")
guardar_tabla_png(met_esp, "Paso 3 — Métricas del modelo espacial",
  "RMSE/MAE/Cob95 sobre el test (diciembre 2025)", file.path(DIR_OUT, "03_metricas.png"))

saveRDS(mod_esp, file.path(DIR_DAT, "modelo_1_espacial.rds"))
cat("\n--- Verificación Paso 3 ---\n")
verif(!is.null(mod_esp$summary.fixed), "modelo espacial ajustó")
verif(is.finite(mod_esp$dic$dic), "DIC finito")
cat("--- PASO 3 completado ---\n")

# ==============================================================================
# PASO 4 — DIAGNÓSTICO DE RESIDUOS (¿RUIDO BLANCO?)
# ==============================================================================
cat("\n", strrep("=", 72), "\n", sep = "")
cat("PASO 4 — Diagnóstico de residuos del modelo espacial\n")
cat(strrep("=", 72), "\n", sep = "")

res_tr <- residuos_train(res_esp)
# Serie temporal: media diaria de los residuos (señal temporal común)
dt_tr[, .res := res_tr]
serie_res <- dt_tr[, .(res = mean(.res, na.rm = TRUE)), by = ID_TIEMPO][order(ID_TIEMPO)]

# ACF (PNG)
lag_max <- 30
png(file.path(DIR_OUT, "04_acf_espacial.png"), width = 900, height = 500)
acf(serie_res$res, lag.max = lag_max, main = "ACF de residuos (modelo espacial)",
    col = "#2166AC", lwd = 2)
dev.off()

# Ljung-Box: ¿queda autocorrelación temporal?
lb <- Box.test(serie_res$res, lag = lag_max, type = "Ljung-Box")
autocorr <- lb$p.value < 0.05

# QQ e histograma (PNG)
png(file.path(DIR_OUT, "04_qq_espacial.png"), width = 600, height = 600)
qqnorm(res_tr, main = "Q-Q residuos (modelo espacial)", pch = 16, cex = 0.4, col = "grey40")
qqline(res_tr, col = "#B2182B", lwd = 2)
dev.off()

g_hist <- ggplot(data.table(r = res_tr), aes(r)) +
  geom_histogram(bins = 60, fill = "#2166AC", color = "white") +
  labs(title = "Histograma de residuos (modelo espacial)",
       subtitle = "Objetivo: distribución simétrica centrada en 0 (ruido blanco)",
       x = "Residuo", y = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(DIR_OUT, "04_hist_espacial.png"), g_hist, width = 8, height = 4.5, dpi = 150)

guardar_tabla_png(
  data.table(Prueba = "Ljung-Box (lag 30)", p_valor = signif(lb$p.value, 4),
             Conclusion = ifelse(autocorr, "Autocorrelación -> AR(1) justificado",
                                 "Sin autocorrelación (ruido blanco)")),
  "Paso 4 — ¿Ruido blanco?", NULL, file.path(DIR_OUT, "04_ljung_box.png"))

cat(sprintf("  Ljung-Box p=%.4f -> %s\n", lb$p.value,
            ifelse(autocorr, "AR(1) justificado", "ruido blanco")))
cat("--- PASO 4 completado ---\n")

# ==============================================================================
# PASO 5 — MODELO ESPACIO-TEMPORAL (AR1 COMPARTIDO) SI EL PASO 4 LO JUSTIFICA
# ==============================================================================
if (autocorr) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("PASO 5 — Modelo con efecto temporal compartido AR(1)\n")
  cat(strrep("=", 72), "\n", sep = "")

  res_st <- ajustar_spde(dt, y_masked, vars_finales, "ar1",
                         compute_fitted = TRUE, int_strategy = "ccd",
                         verbose_inla = TRUE)
  mod_st <- res_st$mod

  # Re-diagnóstico: ACF de residuos del modelo con AR1 (PNG)
  res_tr2 <- residuos_train(res_st)
  dt_tr[, .res2 := res_tr2]
  serie_res2 <- dt_tr[, .(res = mean(.res2, na.rm = TRUE)), by = ID_TIEMPO][order(ID_TIEMPO)]
  png(file.path(DIR_OUT, "05_acf_ar1.png"), width = 900, height = 500)
  acf(serie_res2$res, lag.max = lag_max, main = "ACF de residuos (modelo AR1)",
      col = "#1A9850", lwd = 2)
  dev.off()
  lb2 <- Box.test(serie_res2$res, lag = lag_max, type = "Ljung-Box")
  guardar_tabla_png(
    data.table(Prueba = "Ljung-Box tras AR1 (lag 30)", p_valor = signif(lb2$p.value, 4),
               Conclusion = ifelse(lb2$p.value < 0.05, "Aún autocorrelación",
                                   "Ruido blanco conseguido")),
    "Paso 5 — Residuos tras AR(1)", NULL, file.path(DIR_OUT, "05_ljung_box_ar1.png"))

  met_st <- metricas_modelo(res_st, "Espacio-temporal AR1")
  saveRDS(mod_st, file.path(DIR_DAT, "modelo_1_espaciotemporal.rds"))
  cat(sprintf("  Ljung-Box tras AR1 p=%.4f\n", lb2$p.value))
  cat("--- PASO 5 completado ---\n")
} else {
  met_st <- NULL
  cat("\n[PASO 5 OMITIDO] El modelo espacial ya deja ruido blanco: no se añade AR(1).\n")
}

# ==============================================================================
# PASO 6 — COMPARACIÓN Y MODELO FINAL
# ==============================================================================
cat("\n", strrep("=", 72), "\n", sep = "")
cat("PASO 6 — Comparación de modelos\n")
cat(strrep("=", 72), "\n", sep = "")

tab_comp <- if (!is.null(met_st)) rbind(met_esp, met_st) else met_esp
guardar_tabla_png(tab_comp, "Paso 6 — Comparación de modelos",
  "Menor DIC/WAIC/RMSE y cobertura ~95% = mejor", file.path(DIR_OUT, "06_comparacion.png"))
fwrite(tab_comp, file.path(DIR_OUT, "06_comparacion.csv"))

cat("\n--- Modelo 1 completado (Pasos 0-6) ---\n")
print(tab_comp)
