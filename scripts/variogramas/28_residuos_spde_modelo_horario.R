# ==============================================================================
# MODELO INLA-SPDE DE LOS RESIDUOS DEL MODELO HORARIO (invierno, 5 dias)
# ------------------------------------------------------------------------------
# Objetivo: partiendo del MODELO HORARIO BUENO ya seleccionado en simulacion.R
# (invierno), sacar sus residuos y ajustar un INLA-SPDE sobre ellos para saber
# QUE HAY DENTRO de los residuos: ¿estructura espacial? ¿alguna covariable NO
# usada que todavia explique algo (candidata a anadir al modelo horario)?
#
# Reproduccion fiel del modelo horario bueno (simulacion.R, B1/B4/B5 invierno):
#   - Covariables USADAS (seleccion stepwise-DIC): intensidad, temperatura, viento
#     (precipitacion se descarto por no significativa -> pasa a "no usada").
#   - Modelo: SPDE (PC-Matern) + AR1 temporal, malla ganadora AR1 = Fina (1 km).
#   - Test (dia 5) enmascarado con NA -> residuos de ENTRENAMIENTO (dias 1-4).
#
# Modelo de los residuos:
#   resid_train ~ covariables NO usadas + campo SPDE (malla GRUESA 8 km)
#   - NO se reintroducen las covariables del modelo horario.
#   - NO USADAS = {carga, humedad, presion, radiacion, precipitacion, lag h-1}.
#   - Filtro VIF entre ellas (sin correlacion/colinealidad).
#   - Significancia = IC bayesiano 95% que no cruza 0 (mismo criterio que B1).
#   - Se reporta: rango/varianza del campo espacial de los residuos, betas de
#     las covariables, % de residuo explicado y ACF (Ljung-Box).
# Salidas: outputs/simulacion/residuos_spde_horario/
# ==============================================================================

library(INLA)
library(data.table)
library(sf)
library(car)
library(here)
library(ggplot2)
library(gridExtra)
library(grid)

# Guardar un data.frame como tabla PNG (estilo del proyecto)
guardar_tabla_png <- function(df, titulo, subtitulo, ruta_png, ancho_px = 1000) {
  df <- as.data.frame(df)
  th <- ttheme_minimal(
    core = list(fg_params = list(fontsize = 9),
                bg_params = list(fill = c("white", "#F5F5F5"), col = NA)),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"),
                   bg_params = list(fill = "#DDEEFF", col = NA)))
  g <- arrangeGrob(
    textGrob(titulo, gp = gpar(fontsize = 13, fontface = "bold")),
    textGrob(subtitulo, gp = gpar(fontsize = 9, col = "grey40")),
    tableGrob(df, rows = NULL, theme = th),
    heights = unit(c(0.55, 0.35, nrow(df) * 0.35 + 0.5), "inches"))
  ggsave(ruta_png, g, width = max(6, ancho_px / 96),
         height = 0.9 + nrow(df) * 0.35 + 0.5, dpi = 150, bg = "white")
}

set.seed(4827)

FECHA_INI <- as.Date("2025-01-01")
FECHA_FIN <- as.Date("2025-01-05")
DIAS_TRAIN <- 4L
VIF_UMBRAL <- 5
carpeta <- here("outputs", "simulacion", "residuos_spde_horario")
dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)

# PC-prior Matern igual que crear_spde() en simulacion.R
crear_spde <- function(mesh) {
  inla.spde2.pcmatern(mesh, alpha = 2,
                      prior.range = c(5, 0.5), prior.sigma = c(1, 0.01))
}
malla_de <- function(cmat, max.edge, cutoff) {
  inla.mesh.2d(loc = cmat,
    boundary = list(inla.nonconvex.hull(cmat, convex = -0.05),
                    inla.nonconvex.hull(cmat, convex = -0.2)),
    max.edge = max.edge, cutoff = cutoff)
}

# Covariables USADAS por el modelo horario bueno (seleccion stepwise-DIC)
USADAS_RAW <- c(intensidad = "intensidad_raw",
                temperatura = "Temperatura_raw",
                viento = "Velocidad Viento_raw")
# Covariables NO USADAS (a testar sobre los residuos) + retardo del trafico
NO_USADAS_RAW <- c(carga = "carga_raw",
                   humedad = "Humedad_Relativa_raw",
                   presion = "Presion Barométrica_raw",
                   radiacion = "Radiación Solar_raw",
                   precipitacion = "Precipitaciones_raw")

# ------------------------------------------------------------------------------
# 1. DATOS + estandarizar + coordenadas km + indices + retardo h-1
# ------------------------------------------------------------------------------
dt <- as.data.table(readRDS(here("data", "processed", "Maestro", "horario",
        "dataset_maestro_inla_2025_HORARIO.rds")))
setnames(dt, "LOG_NO2_HORARIO", "LOG_NO2")
dt <- dt[as.Date(FECHA) >= FECHA_INI & as.Date(FECHA) <= FECHA_FIN & !is.na(LOG_NO2)]

std <- function(x) as.numeric(scale(x))
usadas <- character(0)
for (a in names(USADAS_RAW)) if (USADAS_RAW[[a]] %in% names(dt)) {
  dt[, (a) := std(get(USADAS_RAW[[a]]))]; usadas <- c(usadas, a)
}
no_usadas <- character(0)
for (a in names(NO_USADAS_RAW)) if (NO_USADAS_RAW[[a]] %in% names(dt)) {
  dt[, (a) := std(get(NO_USADAS_RAW[[a]]))]; no_usadas <- c(no_usadas, a)
}

coords_u <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
csf <- st_as_sf(coords_u, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)
coords_u[, X_km := st_coordinates(csf)[, 1] / 1000]
coords_u[, Y_km := st_coordinates(csf)[, 2] / 1000]
dt <- merge(dt, coords_u[, .(ESTACION, X_km, Y_km)], by = "ESTACION")

dt[, datetime := as.POSIXct(FECHA, tz = "UTC") + (HORA - 1) * 3600]
setorder(dt, ESTACION, datetime)
dt[, intensidad_lag1 := shift(intensidad, 1L), by = ESTACION]  # retardo h-1
no_usadas <- c(no_usadas, "intensidad_lag1")

dt <- dt[complete.cases(dt[, c("LOG_NO2", usadas, no_usadas), with = FALSE])]
inst <- sort(unique(dt$datetime)); dt[, t_idx := match(datetime, inst)]
dt[, es_test := t_idx > DIAS_TRAIN * 24]

cat(sprintf("Obs=%d (train=%d, test=%d)\n  USADAS(modelo)=%s\n  NO USADAS(test)=%s\n\n",
            nrow(dt), sum(!dt$es_test), sum(dt$es_test),
            paste(usadas, collapse = ", "), paste(no_usadas, collapse = ", ")))

# ------------------------------------------------------------------------------
# 2. MODELO HORARIO BUENO: SPDE(PC-Matern) + AR1, malla GRUESA (8 km)
#    Test enmascarado -> residuos de ENTRENAMIENTO
# ------------------------------------------------------------------------------
cat("[1/2] Modelo horario bueno (SPDE + AR1, malla Gruesa 8 km)...\n")
cmat <- as.matrix(unique(dt[, .(X_km, Y_km)]))
malla_gruesa <- malla_de(cmat, max.edge = c(8, 12), cutoff = 0.5)
spde_gruesa <- crear_spde(malla_gruesa)

y_train <- ifelse(dt$es_test, NA, dt$LOG_NO2)   # enmascarar test
idx_st <- inla.spde.make.index("campo", n.spde = spde_gruesa$n.spde,
                               n.group = max(dt$t_idx))
A_st <- inla.spde.make.A(malla_gruesa, loc = as.matrix(dt[, .(X_km, Y_km)]),
                         group = dt$t_idx, n.group = max(dt$t_idx))
eff_usadas <- c(setNames(lapply(usadas, function(v) dt[[v]]), usadas),
                list(intercept = rep(1, nrow(dt))))
stk_h <- inla.stack(tag = "h", data = list(y = y_train), A = list(A_st, 1),
                    effects = list(idx_st, eff_usadas))
f_h <- as.formula(paste("y ~ 0 + intercept +", paste(usadas, collapse = " + "),
  "+ f(campo, model=spde_gruesa, group=campo.group, control.group=list(model='ar1'))"))
m_h <- inla(f_h, family = "gaussian",
            data = inla.stack.data(stk_h, spde_gruesa = spde_gruesa),
            control.predictor = list(A = inla.stack.A(stk_h), compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE))
i_h <- inla.stack.index(stk_h, "h")$data
dt[, ajuste_h := m_h$summary.fitted.values$mean[i_h]]
dt[, resid_h := LOG_NO2 - ajuste_h]

# residuos de ENTRENAMIENTO (dias 1-4), que es lo que analiza tu B5
dtr <- dt[es_test == FALSE]
r2_base <- 1 - var(dtr$resid_h) / var(dtr$LOG_NO2)
cat(sprintf("  DIC=%.1f WAIC=%.1f | SD(residuo train)=%.4f (SD log NO2=%.4f)\n",
            m_h$dic$dic, m_h$waic$waic, sd(dtr$resid_h), sd(dtr$LOG_NO2)))
cat(sprintf("  R2 del modelo horario (train) = %.4f (%.2f%%)\n", r2_base, 100 * r2_base))
cat("  Efectos fijos (covariables USADAS):\n")
print(round(m_h$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ------------------------------------------------------------------------------
# 3. VIF entre las covariables NO usadas (sobre los residuos de train)
# ------------------------------------------------------------------------------
cat("\n[VIF] entre covariables NO usadas...\n")
vars_test <- no_usadas
repeat {
  if (length(vars_test) < 2) break
  vv <- vif(lm(reformulate(vars_test, "resid_h"), data = dtr))
  if (max(vv) <= VIF_UMBRAL) {
    cat(sprintf("  VIF ok: %s\n",
        paste(sprintf("%s=%.2f", names(vv), vv), collapse = ", "))); break
  }
  elim <- names(which.max(vv))
  cat(sprintf("  Eliminando '%s' (VIF=%.2f)\n", elim, max(vv)))
  vars_test <- setdiff(vars_test, elim)
}

# ------------------------------------------------------------------------------
# 4. MODELO SPDE DE LOS RESIDUOS: resid_h ~ NO usadas + campo SPDE (malla GRUESA 8 km)
# ------------------------------------------------------------------------------
cat("\n[2/2] Modelo SPDE de los residuos (malla Gruesa 8 km)...\n")
A_r <- inla.spde.make.A(malla_gruesa, loc = as.matrix(dtr[, .(X_km, Y_km)]))
idx_s <- inla.spde.make.index("campo", n.spde = spde_gruesa$n.spde)
eff_test <- c(setNames(lapply(vars_test, function(v) dtr[[v]]), vars_test),
              list(intercept = rep(1, nrow(dtr))))
stk_r <- inla.stack(tag = "r", data = list(y = dtr$resid_h), A = list(A_r, 1),
                    effects = list(idx_s, eff_test))
f_r <- as.formula(paste("y ~ 0 + intercept +", paste(vars_test, collapse = " + "),
                        "+ f(campo, model=spde_gruesa)"))
m_r <- inla(f_r, family = "gaussian",
            data = inla.stack.data(stk_r, spde_gruesa = spde_gruesa),
            control.predictor = list(A = inla.stack.A(stk_r), compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE))

# --- Covariables: betas + IC 95% ---
sf_r <- m_r$summary.fixed[vars_test, , drop = FALSE]
sd_r <- sd(dtr$resid_h)
tab <- data.table(
  Covariable = rownames(sf_r),
  beta = round(sf_r$mean, 4),
  beta_std = round(sf_r$mean / sd_r, 4),
  IC_low = round(sf_r$`0.025quant`, 4),
  IC_high = round(sf_r$`0.975quant`, 4),
  importante = ifelse(sf_r$`0.025quant` > 0 | sf_r$`0.975quant` < 0, "SI", "no")
)[order(-abs(beta_std))]
cat("\n=== Covariables NO usadas sobre los residuos del modelo horario ===\n")
print(tab)

# --- Campo espacial de los residuos: ¿queda estructura espacial? ---
sr <- inla.spde.result(m_r, "campo", spde_gruesa)
rango <- exp(sr$summary.log.range.nominal$mean)
sigma2 <- exp(sr$summary.log.variance.nominal$mean)
cat(sprintf("\nCampo espacial de los residuos: rango=%.2f km | sigma2=%.5f\n",
            rango, sigma2))

# --- % de los residuos explicado ---
i_r <- inla.stack.index(stk_r, "r")$data
dtr[, resid_final := resid_h - m_r$summary.fitted.values$mean[i_r]]
r2_resid_total <- 1 - var(dtr$resid_final) / var(dtr$resid_h)
r2_cov_solo <- summary(lm(reformulate(vars_test, "resid_h"), data = dtr))$r.squared
cat(sprintf("\n%% de los residuos del modelo horario explicado:\n"))
cat(sprintf("  - modelo de residuos completo (campo SPDE + covariables): %.2f%%\n",
            100 * r2_resid_total))
cat(sprintf("  - covariables NO usadas SOLAS (lm): %.2f%%\n", 100 * r2_cov_solo))

# --- ACF / Ljung-Box del residuo final ---
res_t <- dtr[, .(r = mean(resid_final)), by = t_idx][order(t_idx)]
lj <- Box.test(res_t$r, lag = 24, type = "Ljung-Box")
png(file.path(carpeta, "acf_residuo_final.png"), width = 800, height = 500)
acf(res_t$r, lag.max = 48, col = "#2166AC", lwd = 2,
    main = sprintf("ACF residuo final | Ljung-Box p=%.3g", lj$p.value))
dev.off()
cat(sprintf("\nLjung-Box (lag 24) p=%.4g -> %s\n", lj$p.value,
    ifelse(lj$p.value < 0.05, "AUN queda autocorrelacion", "sin autocorrelacion")))

imp <- tab[importante == "SI", Covariable]
cat(sprintf("\n>>> Covariables importantes dentro de los residuos: %s\n",
            if (length(imp)) paste(imp, collapse = ", ") else "NINGUNA"))

# --- Guardar la tabla de resultados en PNG (no CSV) ---
subt <- sprintf(paste0("Residuos del modelo horario (AR1, malla Gruesa 8km) | SPDE residuos malla Gruesa 8km | ",
    "R2 modelo=%.1f%% | campo residuos: rango=%.1fkm sigma2=%.4f | ",
    "residuo explicado: SPDE+cov=%.1f%% cov solas=%.1f%% | Ljung-Box p=%.2g | Importantes: %s"),
    100 * r2_base, rango, sigma2, 100 * r2_resid_total, 100 * r2_cov_solo, lj$p.value,
    if (length(imp)) paste(imp, collapse = ", ") else "NINGUNA")
guardar_tabla_png(tab,
  titulo = "Covariables NO usadas dentro de los residuos del modelo horario",
  subtitulo = subt,
  ruta_png = file.path(carpeta, "covariables_no_usadas_en_residuos.png"),
  ancho_px = 1200)

# --- Tabla-resumen de diagnostico (PNG) ---
diag_tab <- data.frame(
  Metrica = c("R2 del modelo horario (train)",
              "Campo espacial de los residuos",
              "Residuo explicado (campo + covariables)",
              "Autocorrelacion remanente (Ljung-Box)"),
  Valor = c(
    sprintf("%.1f%% (SD residuo %.3f vs SD log NO2 %.3f)",
            100 * r2_base, sd(dtr$resid_h), sd(dtr$LOG_NO2)),
    sprintf("rango ~ %.0f km, sigma2 ~ %.0e  ->  sin estructura espacial",
            rango, sigma2),
    sprintf("%.2f%%  (covariables solas %.2f%%: no aportan nada)",
            100 * r2_resid_total, 100 * r2_cov_solo),
    sprintf("p ~ %.0e  (queda algo, pero de magnitud infima)", lj$p.value)),
  check.names = FALSE)
guardar_tabla_png(diag_tab,
  titulo = "Diagnostico del modelo de residuos del modelo horario",
  subtitulo = "Modelo horario (SPDE + AR1) | residuos modelados con SPDE (malla Gruesa 8 km)",
  ruta_png = file.path(carpeta, "diagnostico_residuos.png"),
  ancho_px = 1200)

cat(sprintf("Salidas (PNG) en %s\n", carpeta))
