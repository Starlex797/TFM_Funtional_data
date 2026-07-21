# ==============================================================================
# MODELO DE RESIDUO COMBINADO — NIVEL HORARIO
# ¿Qué queda por explicar de las variables DESCARTADAS (+ su retardo)?
# ==============================================================================
# A diferencia del diagnóstico en dos etapas (script 18), que probaba por
# SEPARADO las covariables débiles (Etapa 2) y el retardo de tráfico (Etapa 3),
# aquí se ajusta UN ÚNICO modelo de segunda etapa que mete TODAS las variables
# no significativas de golpe, junto con el tráfico de la hora anterior (h-1).
#
#   ETAPA 1 — Modelo INLA-SPDE con las covariables candidatas (depuradas por VIF).
#             Se identifican las SIGNIFICATIVAS (por tamaño de efecto, porque a
#             nivel horario N es enorme y todo sale sig. por IC95%) y se construye
#             el residuo parcial restando SOLO su contribución:
#
#                 y*_sh = y_sh − Σ_j  β̂_j · X_j,sh      (j = covariables sig.)
#
#   ETAPA 2 — Una SOLA regresión SPDE del residuo y* sobre TODAS las variables
#             que NO salieron significativas antes, en su hora (sh) y — para el
#             tráfico — también en la hora anterior s(h−1):
#
#                 y*_sh = γ_0 + γ_1·intensidad_sh + γ_11·intensidad_s(h-1)
#                             + γ_1c·carga_sh      + γ_11c·carga_s(h-1)
#                             + γ_2·viento_sh + γ_3·humedad_sh + …   (resto débiles)
#                             + w(s)  (campo espacial Matérn) + ε
#
# Interpretación: si alguna γ es significativa, esa variable "descartada"
#   (o su retardo) todavía explica parte del residuo → nos la estábamos dejando.
#   Lo que quede sin explicar es estructura espacial (campo SPDE) + ruido.
#
# Nota metodológica: residuo de dos etapas (no Frisch–Waugh–Lovell completo);
#   los β̂ de la Etapa 1 se tratan como conocidos, así que la inferencia de la
#   Etapa 2 es orientativa (problema del regresor generado).
# ==============================================================================

library(INLA)
library(data.table)
library(sf)
library(car)
library(here)
library(ggplot2)
library(gridExtra)
library(grid)

set.seed(4827)

# ------------------------------------------------------------------------------
# PARÁMETROS
# ------------------------------------------------------------------------------
MES_ANALISIS <- 1L     # mes a analizar (1 = enero ... 9 = septiembre) · 2025
VIF_UMBRAL   <- 5      # umbral de multicolinealidad
UMBRAL_BETA  <- 0.10   # |beta estandarizado| mínimo para ser "fuerte" (Etapa 1)
# Variables de tráfico a las que se añade el retardo h-1 en la Etapa 2:
TRAFICO_LAG  <- c("intensidad", "carga")

# Covariables candidatas (columnas crudas → se estandarizan dentro del mes)
COVS_RAW   <- c("intensidad_raw", "carga_raw", "Temperatura_raw",
                "Humedad_Relativa_raw", "Precipitaciones_raw",
                "Presion Barométrica_raw", "Radiación Solar_raw",
                "Velocidad Viento_raw")
COVS_ALIAS <- c("intensidad", "carga", "temperatura", "humedad",
                "precipitacion", "presion", "radiacion", "viento")

meses_nombre <- c("enero","febrero","marzo","abril","mayo","junio",
                  "julio","agosto","septiembre","octubre","noviembre","diciembre")
mes_lab <- meses_nombre[MES_ANALISIS]

carpeta_out <- here("outputs", "residuo_combinado_horario",
                    sprintf("horario_%02d_%s", MES_ANALISIS, mes_lab))
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Guardar tabla como PNG (sin navegador)
# ------------------------------------------------------------------------------
guardar_tabla_png <- function(df, titulo, subtitulo = NULL, ruta_png, ancho_px = 950) {
  df <- as.data.frame(df)
  th <- gridExtra::ttheme_minimal(
    core    = list(fg_params = list(fontsize = 9),
                   bg_params = list(fill = c("white", "#F5F5F5"), col = NA)),
    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"),
                   bg_params = list(fill = "#DDEEFF", col = NA)))
  tg <- gridExtra::tableGrob(df, rows = NULL, theme = th)
  tit <- grid::textGrob(titulo, gp = grid::gpar(fontsize = 13, fontface = "bold"))
  if (!is.null(subtitulo)) {
    sub <- grid::textGrob(subtitulo, gp = grid::gpar(fontsize = 9, col = "grey40"))
    comb <- gridExtra::arrangeGrob(tit, sub, tg, nrow = 3,
              heights = grid::unit(c(0.55, 0.35, nrow(df) * 0.32 + 0.5), "inches"))
    alto <- 0.55 + 0.35 + nrow(df) * 0.32 + 0.5
  } else {
    comb <- gridExtra::arrangeGrob(tit, tg, nrow = 2,
              heights = grid::unit(c(0.55, nrow(df) * 0.32 + 0.5), "inches"))
    alto <- 0.55 + nrow(df) * 0.32 + 0.5
  }
  ggplot2::ggsave(ruta_png, comb, width = max(6, ancho_px / 96),
                  height = max(3, alto), dpi = 96, bg = "white")
  cat("  Tabla guardada:", basename(ruta_png), "\n")
}

# ==============================================================================
# 1. CARGA Y PREPARACIÓN (mes horario)
# ==============================================================================
cat(sprintf("\n== Modelo de residuo combinado — %s 2025 (horario) ==\n\n",
            toupper(mes_lab)))

dt <- readRDS(here("data", "processed", "Maestro", "horario",
                   "dataset_maestro_inla_2025_HORARIO.rds"))
setDT(dt)
setnames(dt, "LOG_NO2_HORARIO", "LOG_NO2")
dt[, MES := as.integer(format(FECHA, "%m"))]
dt <- dt[MES == MES_ANALISIS & !is.na(LOG_NO2)]

# Coordenadas UTM en km
coords_u <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
csf <- st_as_sf(coords_u, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)
coords_u[, X_km := st_coordinates(csf)[, 1] / 1000]
coords_u[, Y_km := st_coordinates(csf)[, 2] / 1000]
dt <- merge(dt, coords_u[, .(ESTACION, X_km, Y_km)], by = "ESTACION")

# Estandarizar covariables DENTRO del mes (alias sin caracteres especiales)
covs_disp <- character(0)
for (i in seq_along(COVS_RAW)) {
  v <- COVS_RAW[i]; a <- COVS_ALIAS[i]
  if (v %in% names(dt)) {
    dt[, (a) := as.numeric(scale(get(v)))]
    covs_disp <- c(covs_disp, a)
  }
}

dt <- dt[complete.cases(dt[, c("LOG_NO2", covs_disp), with = FALSE])]
cat(sprintf("Observaciones: %d | estaciones: %d | horas: %d\n",
            nrow(dt), uniqueN(dt$ESTACION), uniqueN(dt$HORA)))

# ==============================================================================
# 2. CONTROL DE MULTICOLINEALIDAD (VIF hacia atrás)
# ==============================================================================
cat("\n--- VIF hacia atrás (umbral =", VIF_UMBRAL, ") ---\n")
vars_vif <- covs_disp
repeat {
  if (length(vars_vif) < 2) break
  fml  <- as.formula(paste("LOG_NO2 ~", paste(vars_vif, collapse = " + ")))
  vifs <- car::vif(lm(fml, data = dt))
  if (max(vifs) <= VIF_UMBRAL) { cat("  -> Sin multicolinealidad.\n"); break }
  drop <- names(which.max(vifs))
  cat(sprintf("  -> Elimino '%s' (VIF = %.2f)\n", drop, max(vifs)))
  vars_vif <- setdiff(vars_vif, drop)
}
vars_excluidas_vif <- setdiff(covs_disp, vars_vif)
cat(sprintf("  Covariables tras VIF: %s\n", paste(vars_vif, collapse = ", ")))
if (length(vars_excluidas_vif) > 0)
  cat(sprintf("  Excluidas por colinealidad: %s\n",
              paste(vars_excluidas_vif, collapse = ", ")))

# ==============================================================================
# 3. MALLA + SPDE (compartida por las dos etapas)
# ==============================================================================
cmat  <- as.matrix(unique(dt[, .(X_km, Y_km)]))
bnd_i <- inla.nonconvex.hull(cmat, convex = -0.05, resolution = 50)
bnd_o <- inla.nonconvex.hull(cmat, convex = -0.2)
malla <- inla.mesh.2d(loc = cmat, boundary = list(bnd_i, bnd_o),
                      max.edge = c(4, 8), cutoff = 0.5)
spde  <- inla.spde2.matern(malla, alpha = 2)
idx   <- inla.spde.make.index("campo", n.spde = spde$n.spde)

es_signif <- function(fixed_row) {
  fixed_row[["0.025quant"]] > 0 | fixed_row[["0.975quant"]] < 0
}

# Ajusta un SPDE gaussiano de y ~ 0 + intercept + <vars> + campo, sobre un
# subconjunto de filas dado (por defecto todo dt). `datos` es un data.table.
ajustar_spde <- function(y, vars, datos, loc) {
  A_mat     <- inla.spde.make.A(malla, loc = as.matrix(loc))
  covs_list <- setNames(lapply(vars, function(v) datos[[v]]), vars)
  stk <- inla.stack(tag = "d", data = list(y = y), A = list(A_mat, 1),
                    effects = list(c(idx, list(intercept = 1)), covs_list))
  fml <- as.formula(paste("y ~ 0 + intercept +", paste(vars, collapse = " + "),
                          "+ f(campo, model = spde)"))
  inla(fml, data = inla.stack.data(stk, spde = spde), family = "gaussian",
       control.predictor = list(A = inla.stack.A(stk), compute = FALSE),
       control.compute = list(dic = TRUE, waic = TRUE),
       control.inla = list(int.strategy = "eb"))
}

# ==============================================================================
# 4. ETAPA 1 — modelo completo, identificar covariables significativas/fuertes
# ==============================================================================
cat("\n[ETAPA 1] Ajustando modelo SPDE con las covariables candidatas...\n")
m1     <- ajustar_spde(dt$LOG_NO2, vars_vif, dt, dt[, .(X_km, Y_km)])
fixed1 <- m1$summary.fixed
fixed1$sig <- es_signif(fixed1)

# Partición por TAMAÑO DE EFECTO (con N enorme, el IC95% siempre excluye 0)
betas_cov  <- setNames(fixed1[vars_vif, "mean"], vars_vif)
vars_sig   <- vars_vif[abs(betas_cov) >= UMBRAL_BETA]   # fuertes -> se restan
vars_nosig <- vars_vif[abs(betas_cov) <  UMBRAL_BETA]   # débiles -> Etapa 2

cat(sprintf("  FUERTES |β|>=%.2f (se restan en Etapa 1): %s\n", UMBRAL_BETA,
            ifelse(length(vars_sig) > 0, paste(vars_sig, collapse = ", "), "ninguna")))
cat(sprintf("  DÉBILES |β|<%.2f  (van a la Etapa 2): %s\n", UMBRAL_BETA,
            ifelse(length(vars_nosig) > 0, paste(vars_nosig, collapse = ", "), "ninguna")))

# --- Residuo parcial: y* = y − Σ β̂_sig · X_sig ---
dt[, y_estrella := LOG_NO2]
for (v in vars_sig) dt[, y_estrella := y_estrella - fixed1[v, "mean"] * get(v)]

sd_y  <- sd(dt$LOG_NO2)
sd_ye <- sd(dt$y_estrella)
cat(sprintf("  SD(y) = %.4f  ->  SD(y*) tras quitar fuertes = %.4f (%.1f%% de la var. original)\n",
            sd_y, sd_ye, 100 * (sd_ye / sd_y)^2))

# ==============================================================================
# 5. CONSTRUIR RETARDO DE TRÁFICO h-1 (dentro de cada estación, gap = 1 hora)
# ==============================================================================
dt[, datetime := as.POSIXct(FECHA, tz = "UTC") + (HORA - 1) * 3600]
setorder(dt, ESTACION, datetime)
for (v in TRAFICO_LAG) {
  if (v %in% covs_disp) dt[, (paste0(v, "_lag1")) := shift(get(v), 1), by = ESTACION]
}
dt[, dtprev := shift(datetime, 1), by = ESTACION]
dt[, gap_h  := as.numeric(difftime(datetime, dtprev, units = "hours"))]

# Términos de la Etapa 2: TODAS las débiles contemporáneas + retardo de tráfico
lag_terms  <- paste0(intersect(TRAFICO_LAG, covs_disp), "_lag1")
vars_e2    <- c(vars_nosig, lag_terms)

# Subconjunto con retardo válido (gap exacto de 1 h y sin NAs en los términos E2)
dt_e2 <- dt[gap_h == 1]
dt_e2 <- dt_e2[complete.cases(dt_e2[, vars_e2, with = FALSE])]
sd_ye_sub <- sd(dt_e2$y_estrella)   # baseline del residuo en el subconjunto E2
cat(sprintf("\n  Observaciones con retardo h-1 válido: %d (de %d)\n",
            nrow(dt_e2), nrow(dt)))
cat(sprintf("  Términos de la Etapa 2: %s\n", paste(vars_e2, collapse = ", ")))

# ==============================================================================
# 6. ETAPA 2 — UNA ÚNICA regresión de y* sobre TODAS las no-significativas + h-1
# ==============================================================================
cat("\n[ETAPA 2] Regresión combinada del residuo y* (débiles + tráfico h-1)...\n")

# VIF de control sobre el conjunto combinado (debe seguir limpio)
if (length(vars_e2) >= 2) {
  vif2 <- car::vif(lm(as.formula(paste("y_estrella ~",
                      paste(vars_e2, collapse = " + "))), data = dt_e2))
  cat(sprintf("  VIF Etapa 2 (máx = %.2f): %s\n", max(vif2),
              paste(sprintf("%s=%.2f", names(vif2), vif2), collapse = ", ")))
}

m2     <- ajustar_spde(dt_e2$y_estrella, vars_e2, dt_e2, dt_e2[, .(X_km, Y_km)])
fixed2 <- m2$summary.fixed
fixed2$sig <- es_signif(fixed2)

# Residuo final tras quitar la contribución conjunta de la Etapa 2
dt_e2[, resid2 := y_estrella]
for (v in vars_e2) dt_e2[, resid2 := resid2 - fixed2[v, "mean"] * get(v)]
sd_resid2 <- sd(dt_e2$resid2)

vars_recuperadas <- setdiff(rownames(fixed2)[fixed2$sig], "intercept")
cat(sprintf("  Variables que TODAVÍA explican el residuo (γ sig.): %s\n",
            ifelse(length(vars_recuperadas) > 0,
                   paste(vars_recuperadas, collapse = ", "),
                   "ninguna (residuo limpio de covariables)")))
cat(sprintf("  SD(y*) = %.4f  ->  SD(resid) tras Etapa 2 = %.4f (reduce %.2f%%)\n",
            sd_ye_sub, sd_resid2, 100 * (1 - sd_resid2 / sd_ye_sub)))

# ==============================================================================
# 7. SALIDAS
# ==============================================================================

# --- Tabla 1: estado de cada covariable candidata (Etapa 1) ---
estado <- data.frame(
  Variable = covs_disp,
  Estado = ifelse(covs_disp %in% vars_excluidas_vif, "Excluida (colinealidad)",
           ifelse(covs_disp %in% vars_sig, "Fuerte (Etapa 1, se resta)",
                  "Débil (Etapa 2)")),
  Beta_Etapa1 = sapply(covs_disp, function(v)
                 if (v %in% rownames(fixed1)) round(fixed1[v, "mean"], 4) else NA),
  row.names = NULL
)
estado <- estado[order(-abs(estado$Beta_Etapa1)), ]
guardar_tabla_png(estado,
  titulo    = sprintf("Etapa 1 — estado de las covariables — horario %s 2025", mes_lab),
  subtitulo = sprintf("Se restan las FUERTES |β|>=%.2f para formar y*. Las DÉBILES van a la Etapa 2.",
                      UMBRAL_BETA),
  ruta_png  = file.path(carpeta_out, "01_estado_covariables.png"))
fwrite(estado, file.path(carpeta_out, "01_estado_covariables.csv"))

# --- Tabla 2: la regresión combinada de la Etapa 2 (todas las γ) ---
etiqueta_var <- function(v) {
  if (grepl("_lag1$", v)) sprintf("%s (h-1)", sub("_lag1$", "", v)) else sprintf("%s (h)", v)
}
tabla_e2 <- data.frame(
  Termino       = c("intercept (γ0)", sapply(vars_e2, etiqueta_var)),
  Gamma         = round(fixed2[c("intercept", vars_e2), "mean"], 4),
  IC025         = round(fixed2[c("intercept", vars_e2), "0.025quant"], 4),
  IC975         = round(fixed2[c("intercept", vars_e2), "0.975quant"], 4),
  Significativa = ifelse(fixed2[c("intercept", vars_e2), "sig"], "Sí", "no"),
  row.names     = NULL
)
guardar_tabla_png(tabla_e2,
  titulo    = "Etapa 2 — regresión combinada del residuo y*",
  subtitulo = sprintf("y* ~ γ0 + Σ débiles(h) + tráfico(h-1) + campo SPDE · N=%d · SD(y*)=%.4f -> %.4f",
                      nrow(dt_e2), sd_ye_sub, sd_resid2),
  ruta_png  = file.path(carpeta_out, "02_etapa2_combinada.png"))
fwrite(tabla_e2, file.path(carpeta_out, "02_etapa2_combinada.csv"))

# --- Gráfico: reducción de la desviación típica por etapa ---
df_sd <- data.frame(
  Fase = factor(c("y (original)", "y* (tras Etapa 1)", "resid (tras Etapa 2)"),
                levels = c("y (original)", "y* (tras Etapa 1)", "resid (tras Etapa 2)")),
  SD   = c(sd_y, sd_ye_sub, sd_resid2)
)
ggplot(df_sd, aes(x = Fase, y = SD, fill = Fase)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("%.4f", SD)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("#1a5276", "#7fb3d5", "#c0392b"), guide = "none") +
  coord_cartesian(ylim = c(0, max(df_sd$SD) * 1.08)) +
  labs(title = sprintf("Reducción de variabilidad por etapa — horario %s 2025", mes_lab),
       subtitle = "y -> y* (quitando covariables fuertes) -> residuo (quitando débiles + tráfico h-1)",
       x = NULL, y = "Desviación típica de log(NO₂)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 9))
ggsave(file.path(carpeta_out, "03_reduccion_sd_por_etapa.png"),
       width = 9, height = 6, dpi = 300, bg = "white")

cat(sprintf("\n✅ Modelo de residuo combinado completado. Salidas en:\n   %s\n", carpeta_out))
