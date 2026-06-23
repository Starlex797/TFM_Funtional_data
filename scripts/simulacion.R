# ==============================================================================
# SIMULACIÓN INLA-SPDE: TAMAÑO DE MALLA Y COMPARACIÓN DE MODELOS
# ==============================================================================
# Objetivo:
#   - Seleccionar 19 días CONSECUTIVOS de noviembre (train: 1-15, test: 16-19)
#   - Comparar 3 tamaños de malla (gruesa, media, fina)
#   - Comparar 2 modelos: solo espacial vs. espacio-temporal (AR1)
# Métricas sobre días de test: RMSE, MAE, Coverage 95%
#
# Nota: Seguimos la estructura del Paso 3 (Spde.R) para la construcción del
# SPDE, índices y matrices A.
# ==============================================================================

library(INLA)
library(fmesher)
library(data.table)
library(sf)
library(here)
library(ggplot2)
library(car)

set.seed(8314)

carpeta_sim <- here("outputs", "simulacion")
dir.create(carpeta_sim, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. CARGA DE DATOS Y SELECCIÓN DE 19 DÍAS CONSECUTIVOS DE NOVIEMBRE
# ==============================================================================

dt_maestro <- readRDS(here("data", "processed",
                           "dataset_maestro_inla_2025_DIARIO.rds"))
setDT(dt_maestro)

# Fechas disponibles en noviembre
fechas_noviembre <- sort(unique(dt_maestro$FECHA[
  format(dt_maestro$FECHA, "%m") == "11"
]))

if (length(fechas_noviembre) < 19) {
  stop(sprintf("Solo hay %d días de noviembre; se necesitan al menos 19.",
               length(fechas_noviembre)))
}

# Tomamos los primeros 19 días consecutivos de noviembre
fechas_sim   <- fechas_noviembre[1:19]
fechas_train <- fechas_sim[1:15]
fechas_test  <- fechas_sim[16:19]

cat("Período de simulación:\n")
cat(sprintf("  Train: %s  →  %s  (%d días)\n",
            min(fechas_train), max(fechas_train), length(fechas_train)))
cat(sprintf("  Test:  %s  →  %s  (%d días)\n",
            min(fechas_test),  max(fechas_test),  length(fechas_test)))

# Filtrar y crear ID_TIEMPO_SIM (1 a 19, ordenado por fecha)
dt_sim <- dt_maestro[FECHA %in% fechas_sim]
dt_sim[, ID_TIEMPO_SIM := match(FECHA, fechas_sim)]
setorder(dt_sim, ID_TIEMPO_SIM, ESTACION)

n_train <- 15L
n_test  <- 4L
ndays   <- 19L

cat(sprintf("\nFilas: %d | Estaciones: %d | Días: %d (train=%d, test=%d)\n",
            nrow(dt_sim), uniqueN(dt_sim$ESTACION), ndays, n_train, n_test))

# ==============================================================================
# 2. COORDENADAS UTM EN KM (igual que en Paso 3 / Spde.R)
# ==============================================================================

coords_unicas <- unique(dt_sim[, .(ESTACION, LONGITUD, LATITUD)])
coords_sf <- st_as_sf(coords_unicas,
                      coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)

coords_unicas[, X_km := st_coordinates(coords_sf)[, 1] / 1000]
coords_unicas[, Y_km := st_coordinates(coords_sf)[, 2] / 1000]

dt_sim <- merge(dt_sim, coords_unicas[, .(ESTACION, X_km, Y_km)],
                by = "ESTACION", all.x = TRUE)
setorder(dt_sim, ID_TIEMPO_SIM, ESTACION)

# Matriz de coordenadas (una fila por observación, replicadas por día)
coords_puntos <- as.matrix(dt_sim[, .(X_km, Y_km)])
# Coordenadas únicas para definir los límites de la malla
coords_matriz <- as.matrix(unique(dt_sim[, .(X_km, Y_km)]))

# Fronteras para la construcción de mallas
bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.05, resolution = 50)
bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.2)

# ==============================================================================
# 3. COVARIABLES Y VECTOR RESPUESTA
# ==============================================================================

# Variable respuesta: NA en test para que INLA prediga esas filas
y_train_test <- ifelse(dt_sim$FECHA %in% fechas_train,
                       dt_sim$LOG_NO2_DIARIO, NA_real_)

covariables_candidatas <- list(
  trafico_intensidad  = dt_sim$intensidad,
  temperatura         = dt_sim$Temperatura,
  humedad             = dt_sim$Humedad_Relativa,
  precipitacion       = dt_sim$Precipitaciones,
  presion_barometrica = dt_sim$`Presion Barométrica`,
  radiacion_solar     = dt_sim$`Radiación Solar`,
  velocidad_viento    = dt_sim$`Velocidad Viento`
)

# ------------------------------------------------------------------------------
# 3b. SELECCIÓN DE COVARIABLES POR VIF HACIA ATRÁS (umbral = 5)
#     Se usa solo el conjunto de entrenamiento para evitar data leakage.
# ------------------------------------------------------------------------------
VIF_UMBRAL <- 5

# Construir data.frame de entrenamiento con respuesta + covariables candidatas
df_train_vif <- as.data.frame(
  lapply(covariables_candidatas, function(x) x[dt_sim$FECHA %in% fechas_train])
)
df_train_vif$LOG_NO2 <- dt_sim$LOG_NO2_DIARIO[dt_sim$FECHA %in% fechas_train]

# Eliminar filas con NA en cualquier columna
df_train_vif <- na.omit(df_train_vif)

vars_vif <- names(covariables_candidatas)

cat("\n--- VIF hacia atrás (umbral =", VIF_UMBRAL, ") ---\n")

repeat {
  formula_vif <- as.formula(
    paste("LOG_NO2 ~", paste(vars_vif, collapse = " + "))
  )
  lm_vif <- lm(formula_vif, data = df_train_vif)
  vif_vals <- car::vif(lm_vif)

  cat(sprintf("  Variables: %s\n", paste(vars_vif, collapse = ", ")))
  cat("  VIF actuales:\n")
  print(round(vif_vals, 3))

  vif_max <- max(vif_vals)
  if (vif_max <= VIF_UMBRAL) {
    cat(sprintf("  → Todas las variables tienen VIF ≤ %.1f. Selección finalizada.\n",
                VIF_UMBRAL))
    break
  }

  var_eliminar <- names(which.max(vif_vals))
  cat(sprintf("  → Eliminando '%s' (VIF = %.3f)\n\n", var_eliminar, vif_max))
  vars_vif <- setdiff(vars_vif, var_eliminar)

  if (length(vars_vif) < 2) {
    warning("VIF hacia atrás dejó menos de 2 variables; se revisa el umbral.")
    break
  }
}

cat(sprintf("\nVariables retenidas (%d): %s\n",
            length(vars_vif), paste(vars_vif, collapse = ", ")))

# Guardar tabla resumen VIF final
tabla_vif <- data.frame(
  Variable = names(vif_vals),
  VIF      = round(as.numeric(vif_vals), 4),
  row.names = NULL
)
tabla_vif <- tabla_vif[order(tabla_vif$VIF, decreasing = TRUE), ]

write.csv(tabla_vif,
          file.path(carpeta_sim, "seleccion_vif_covariables.csv"),
          row.names = FALSE)

cat("Tabla VIF guardada en: seleccion_vif_covariables.csv\n")

# Filtrar la lista de covariables a las retenidas por VIF
covariables <- covariables_candidatas[vars_vif]

# Índices y valores reales de las filas de test
idx_test_filas <- which(dt_sim$FECHA %in% fechas_test)
y_test_real    <- dt_sim$LOG_NO2_DIARIO[idx_test_filas]

# ------------------------------------------------------------------------------
# 3c. SIGNIFICANCIA BAYESIANA (IC95% excluye 0 en modelo INLA SPDE completo)
#     Se construye una malla/SPDE de referencia (media, 4 km) compartida
#     con el stepwise del bloque 3d.
# ------------------------------------------------------------------------------

malla_sw <- inla.mesh.2d(
  loc      = coords_matriz,
  boundary = list(bnd_inner, bnd_outer),
  max.edge = c(4, 8),
  cutoff   = 0.5
)
spde_sw <- inla.spde2.matern(mesh = malla_sw, alpha = 2)

indice_sw <- inla.spde.make.index(
  name    = "campo_espacial",
  n.spde  = spde_sw$n.spde,
  n.group = ndays
)
A_sw <- inla.spde.make.A(
  mesh    = malla_sw,
  loc     = coords_puntos,
  group   = dt_sim$ID_TIEMPO_SIM,
  n.group = ndays
)
cat(sprintf("\n[Selección] Malla de referencia: %d nodos\n", malla_sw$n))

# Modelo INLA SPDE AR1 con todas las variables retenidas por VIF
formula_full <- as.formula(paste(
  "y ~ 0 + intercept +",
  paste(vars_vif, collapse = " + "),
  "+ f(campo_espacial, model = spde_sw,",
  "    group = campo_espacial.group,",
  "    control.group = list(model = 'ar1'))"
))

stk_full <- inla.stack(
  tag      = "full",
  data     = list(y = y_train_test),
  A        = list(A_sw, 1),
  effects  = list(
    c(indice_sw, list(intercept = 1)),
    covariables_candidatas[vars_vif]
  ),
  compress = FALSE
)

cat("  Ajustando modelo INLA completo (significancia)...\n")
modelo_full <- inla(
  formula           = formula_full,
  data              = inla.stack.data(stk_full, spde = spde_sw),
  family            = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_full), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = FALSE, cpo = FALSE),
  control.inla      = list(strategy = "laplace"),
  verbose           = FALSE
)

# Efectos fijos (excluir intercepto)
sf_full <- modelo_full$summary.fixed
sf_vars <- sf_full[rownames(sf_full) != "intercept", ]

# Variable significativa ↔ IC95% no contiene el 0
sig_mask    <- sf_vars[, "0.025quant"] > 0 | sf_vars[, "0.975quant"] < 0
vars_sig    <- rownames(sf_vars)[ sig_mask]
vars_no_sig <- rownames(sf_vars)[!sig_mask]

tabla_sig <- data.frame(
  Variable = rownames(sf_vars),
  Media    = round(sf_vars$mean,            4),
  SD       = round(sf_vars$sd,              4),
  Q2.5     = round(sf_vars[, "0.025quant"], 4),
  Q97.5    = round(sf_vars[, "0.975quant"], 4),
  Sig_95   = ifelse(sig_mask, "Sí", "No"),
  row.names = NULL
)

cat("\n--- Significancia bayesiana (IC95% excluye 0) ---\n")
print(tabla_sig)
cat(sprintf("\n  Significativas   (%d): %s\n", length(vars_sig),
            paste(vars_sig,    collapse = ", ")))
cat(sprintf("  No significativas(%d): %s\n", length(vars_no_sig),
            paste(vars_no_sig, collapse = ", ")))

write.csv(tabla_sig,
          file.path(carpeta_sim, "significancia_bayesiana_covariables.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 3d. STEPWISE BASADO EN DIC (INLA SPDE AR1, misma malla de referencia)
#     Punto de partida : variables significativas (vars_sig)
#     Pool completo    : variables retenidas por VIF (vars_vif)
#     Criterio de parada: mejora en DIC < DIC_MEJORA_MIN
# ------------------------------------------------------------------------------

DIC_MEJORA_MIN <- 2

# Helper: ajusta INLA SPDE AR1 para un subconjunto de variables;
# devuelve DIC, RMSE y MAE sobre los días de test.
# Usa spde_sw, indice_sw, A_sw del entorno global.
ajustar_sw <- function(vars_modelo) {
  formula_sw <- as.formula(paste(
    "y ~ 0 + intercept +",
    paste(vars_modelo, collapse = " + "),
    "+ f(campo_espacial, model = spde_sw,",
    "    group = campo_espacial.group,",
    "    control.group = list(model = 'ar1'))"
  ))

  stk_sw <- inla.stack(
    tag      = "sw",
    data     = list(y = y_train_test),
    A        = list(A_sw, 1),
    effects  = list(
      c(indice_sw, list(intercept = 1)),
      covariables_candidatas[vars_modelo]
    ),
    compress = FALSE
  )

  mod <- tryCatch(
    inla(
      formula           = formula_sw,
      data              = inla.stack.data(stk_sw, spde = spde_sw),
      family            = "gaussian",
      control.predictor = list(A = inla.stack.A(stk_sw), compute = TRUE),
      control.compute   = list(dic = TRUE, waic = FALSE, cpo = FALSE),
      control.inla      = list(strategy = "laplace"),
      verbose           = FALSE
    ),
    error = function(e) { message("  ERROR ajustar_sw: ", e$message); NULL }
  )

  if (is.null(mod)) return(list(DIC = Inf, RMSE = NA_real_, MAE = NA_real_))

  idx_data  <- inla.stack.index(stk_sw, tag = "sw")$data
  pred_test <- mod$summary.fitted.values$mean[idx_data][idx_test_filas]
  rmse <- sqrt(mean((pred_test - y_test_real)^2, na.rm = TRUE))
  mae  <- mean(abs(pred_test  - y_test_real),    na.rm = TRUE)

  list(DIC = mod$dic$dic, RMSE = rmse, MAE = mae)
}

# ---- Bucle stepwise bidireccional ----
vars_actuales <- if (length(vars_sig) > 0) vars_sig else vars_vif
tabla_sw      <- list()

cat(sprintf("\n--- Stepwise DIC (ΔDIC mín. = %g) ---\n", DIC_MEJORA_MIN))

res_actual <- ajustar_sw(vars_actuales)
cat(sprintf("Modelo inicial : %s\n  DIC=%.2f | RMSE=%.4f | MAE=%.4f\n",
            paste(vars_actuales, collapse = " + "),
            res_actual$DIC, res_actual$RMSE, res_actual$MAE))

tabla_sw[[1]] <- data.frame(
  Iteracion = 0L, Accion = "inicial", Variable = NA_character_,
  Variables = paste(vars_actuales, collapse = " + "),
  DIC  = round(res_actual$DIC,  2),
  RMSE = round(res_actual$RMSE, 4),
  MAE  = round(res_actual$MAE,  4),
  stringsAsFactors = FALSE
)

for (iter in seq_len(length(vars_vif))) {

  cat(sprintf("\n[Iter %d] Actual: %s  (DIC=%.2f)\n",
              iter, paste(vars_actuales, collapse = " + "), res_actual$DIC))

  cands <- data.frame(accion = character(), variable = character(),
                      DIC = numeric(), RMSE = numeric(), MAE = numeric(),
                      stringsAsFactors = FALSE)

  # Probar añadir cada variable del pool no incluida aún
  for (v in setdiff(vars_vif, vars_actuales)) {
    r <- ajustar_sw(c(vars_actuales, v))
    cat(sprintf("  + %-22s DIC=%8.2f  RMSE=%.4f  MAE=%.4f\n",
                v, r$DIC, r$RMSE, r$MAE))
    cands <- rbind(cands, data.frame(accion = "añadir",   variable = v,
                                     DIC = r$DIC, RMSE = r$RMSE, MAE = r$MAE,
                                     stringsAsFactors = FALSE))
  }

  # Probar eliminar cada variable activa (mínimo 1 variable en el modelo)
  if (length(vars_actuales) >= 2) {
    for (v in vars_actuales) {
      r <- ajustar_sw(setdiff(vars_actuales, v))
      cat(sprintf("  - %-22s DIC=%8.2f  RMSE=%.4f  MAE=%.4f\n",
                  v, r$DIC, r$RMSE, r$MAE))
      cands <- rbind(cands, data.frame(accion = "eliminar", variable = v,
                                       DIC = r$DIC, RMSE = r$RMSE, MAE = r$MAE,
                                       stringsAsFactors = FALSE))
    }
  }

  mejor  <- cands[which.min(cands$DIC), ]
  mejora <- res_actual$DIC - mejor$DIC

  if (mejora < DIC_MEJORA_MIN) {
    cat(sprintf("  → Sin mejora suficiente (ΔDIC=%.2f). Stepwise convergido.\n",
                mejora))
    break
  }

  vars_actuales <- if (mejor$accion == "añadir")
    c(vars_actuales, mejor$variable) else setdiff(vars_actuales, mejor$variable)

  res_actual <- list(DIC = mejor$DIC, RMSE = mejor$RMSE, MAE = mejor$MAE)

  cat(sprintf("  → %s '%s'  ΔDIC=%.2f | RMSE=%.4f | MAE=%.4f\n",
              mejor$accion, mejor$variable, mejora, mejor$RMSE, mejor$MAE))

  tabla_sw[[iter + 1]] <- data.frame(
    Iteracion = iter,
    Accion    = mejor$accion,
    Variable  = mejor$variable,
    Variables = paste(vars_actuales, collapse = " + "),
    DIC  = round(mejor$DIC,  2),
    RMSE = round(mejor$RMSE, 4),
    MAE  = round(mejor$MAE,  4),
    stringsAsFactors = FALSE
  )
}

tabla_stepwise <- do.call(rbind, tabla_sw)

cat(sprintf("\nVariables finales stepwise (%d): %s\n",
            length(vars_actuales), paste(vars_actuales, collapse = ", ")))

write.csv(tabla_stepwise,
          file.path(carpeta_sim, "stepwise_dic_covariables.csv"),
          row.names = FALSE)
cat("Tabla stepwise guardada: stepwise_dic_covariables.csv\n")

# Actualizar vars_vif y covariables con la selección final del stepwise
vars_vif    <- vars_actuales
covariables <- covariables_candidatas[vars_vif]

rm(malla_sw, spde_sw, indice_sw, A_sw, stk_full, modelo_full)

# ==============================================================================
# 4. CONFIGURACIÓN DE LAS 3 MALLAS
# ==============================================================================

config_mallas <- list(
  gruesa = list(max.edge = c(8, 12), cutoff = 0.5,  label = "Gruesa (8 km)"),
  media  = list(max.edge = c(4,  8), cutoff = 0.5,  label = "Media (4 km)"),
  fina   = list(max.edge = c(1,  4), cutoff = 0.25, label = "Fina (1 km)")
)

# ==============================================================================
# 5. FÓRMULAS COMPARTIDAS (consistente con Paso 5)
# ==============================================================================

efectos_fijos <- as.formula(
  paste("y ~ 0 + intercept +", paste(vars_vif, collapse = " + "))
)

# Modelo espacio-temporal: campo SPDE con correlación AR1 entre días
formula_st <- update(efectos_fijos,
  . ~ . + f(campo_espacial,
            model         = spde,
            group         = campo_espacial.group,
            control.group = list(model = "ar1")))

# Modelo solo espacial: campo SPDE sin dimensión temporal
formula_s <- update(efectos_fijos,
  . ~ . + f(campo_espacial_s, model = spde))

# ==============================================================================
# 6. BUCLE PRINCIPAL: 3 MALLAS × 2 MODELOS
# ==============================================================================

resultados <- list()

for (malla_nombre in names(config_mallas)) {

  cfg <- config_mallas[[malla_nombre]]
  cat("\n", strrep("=", 70), "\n")
  cat(sprintf(" MALLA: %s\n", cfg$label))
  cat(strrep("=", 70), "\n")

  # --- Construir malla ---
  malla <- inla.mesh.2d(
    loc      = coords_matriz,
    boundary = list(bnd_inner, bnd_outer),
    max.edge = cfg$max.edge,
    cutoff   = cfg$cutoff
  )
  cat(sprintf("  Nodos en la malla: %d\n", malla$n))

  # --- SPDE (Matérn, alpha = 2, igual que Spde.R) ---
  spde <- inla.spde2.matern(mesh = malla, alpha = 2)

  # ==========================================================================
  # 6a. MODELO ESPACIO-TEMPORAL (SPDE + AR1)
  #     Siguiendo Spde.R: n.group = ndays (19), group = ID_TIEMPO_SIM
  # ==========================================================================
  cat("\n  [ST] Construyendo índice y Matriz A espacio-temporal...\n")

  indice_st <- inla.spde.make.index(
    name    = "campo_espacial",
    n.spde  = spde$n.spde,
    n.group = ndays
  )

  A_st <- inla.spde.make.A(
    mesh    = malla,
    loc     = coords_puntos,
    group   = dt_sim$ID_TIEMPO_SIM,
    n.group = ndays
  )

  # Validación de dimensiones (igual que Spde.R)
  stopifnot(nrow(A_st) == nrow(dt_sim))
  stopifnot(ncol(A_st) == malla$n * ndays)

  stk_st <- inla.stack(
    tag      = "sim_st",
    data     = list(y = y_train_test),
    A        = list(A_st, 1),
    effects  = list(
      c(indice_st, list(intercept = 1)),
      covariables
    ),
    compress = FALSE
  )

  cat("  [ST] Ajustando modelo...\n")
  t0 <- Sys.time()

  modelo_st <- tryCatch(
    inla(
      formula           = formula_st,
      data              = inla.stack.data(stk_st, spde = spde),
      family            = "gaussian",
      control.predictor = list(A = inla.stack.A(stk_st), compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE),
      control.inla      = list(strategy = "laplace"),
      verbose           = FALSE
    ),
    error = function(e) { message("  ERROR ST: ", e$message); NULL }
  )
  t_st <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (!is.null(modelo_st)) {
    # Extraer predicciones: idx_st_data mapea filas del stack → fitted.values
    # Las predicciones están en el mismo orden que las filas de dt_sim
    idx_st_data <- inla.stack.index(stk_st, tag = "sim_st")$data

    pred_st_mean <- modelo_st$summary.fitted.values$mean[idx_st_data]
    pred_st_sd   <- modelo_st$summary.fitted.values$sd[idx_st_data]

    # Seleccionamos solo los días de test
    pred_st_mean_test <- pred_st_mean[idx_test_filas]
    pred_st_sd_test   <- pred_st_sd[idx_test_filas]

    rmse_st  <- sqrt(mean((pred_st_mean_test - y_test_real)^2, na.rm = TRUE))
    mae_st   <- mean(abs(pred_st_mean_test  - y_test_real),    na.rm = TRUE)
    lo95_st  <- pred_st_mean_test - 1.96 * pred_st_sd_test
    hi95_st  <- pred_st_mean_test + 1.96 * pred_st_sd_test
    cov95_st <- mean(y_test_real >= lo95_st & y_test_real <= hi95_st, na.rm = TRUE)

    cat(sprintf("  [ST] RMSE=%.4f | MAE=%.4f | Cov95=%.1f%% | %.2f min\n",
                rmse_st, mae_st, cov95_st * 100, t_st))

    resultados[[length(resultados) + 1]] <- data.table(
      Malla      = cfg$label,
      Modelo     = "Espacio-temporal (AR1)",
      n_nodos    = malla$n,
      DIC        = modelo_st$dic$dic,
      WAIC       = modelo_st$waic$waic,
      RMSE       = rmse_st,
      MAE        = mae_st,
      Cov95      = round(cov95_st * 100, 1),
      Tiempo_min = round(t_st, 2)
    )
  }

  # ==========================================================================
  # 6b. MODELO SOLO ESPACIAL
  #     Sin n.group ni group: campo SPDE estático, igual que en Spde.R
  # ==========================================================================
  cat("\n  [S]  Construyendo índice y Matriz A solo espacial...\n")

  indice_s_solo <- inla.spde.make.index(
    name   = "campo_espacial_s",
    n.spde = spde$n.spde
  )

  A_s <- inla.spde.make.A(mesh = malla, loc = coords_puntos)

  # Validación de dimensiones
  stopifnot(nrow(A_s) == nrow(dt_sim))
  stopifnot(ncol(A_s) == malla$n)

  stk_s <- inla.stack(
    tag      = "sim_s",
    data     = list(y = y_train_test),
    A        = list(A_s, 1),
    effects  = list(
      c(indice_s_solo, list(intercept = 1)),
      covariables
    ),
    compress = FALSE
  )

  cat("  [S]  Ajustando modelo...\n")
  t0 <- Sys.time()

  modelo_s <- tryCatch(
    inla(
      formula           = formula_s,
      data              = inla.stack.data(stk_s, spde = spde),
      family            = "gaussian",
      control.predictor = list(A = inla.stack.A(stk_s), compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE),
      control.inla      = list(strategy = "laplace"),
      verbose           = FALSE
    ),
    error = function(e) { message("  ERROR S: ", e$message); NULL }
  )
  t_s <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (!is.null(modelo_s)) {
    idx_s_data <- inla.stack.index(stk_s, tag = "sim_s")$data

    pred_s_mean <- modelo_s$summary.fitted.values$mean[idx_s_data]
    pred_s_sd   <- modelo_s$summary.fitted.values$sd[idx_s_data]

    pred_s_mean_test <- pred_s_mean[idx_test_filas]
    pred_s_sd_test   <- pred_s_sd[idx_test_filas]

    rmse_s  <- sqrt(mean((pred_s_mean_test - y_test_real)^2, na.rm = TRUE))
    mae_s   <- mean(abs(pred_s_mean_test  - y_test_real),    na.rm = TRUE)
    lo95_s  <- pred_s_mean_test - 1.96 * pred_s_sd_test
    hi95_s  <- pred_s_mean_test + 1.96 * pred_s_sd_test
    cov95_s <- mean(y_test_real >= lo95_s & y_test_real <= hi95_s, na.rm = TRUE)

    cat(sprintf("  [S]  RMSE=%.4f | MAE=%.4f | Cov95=%.1f%% | %.2f min\n",
                rmse_s, mae_s, cov95_s * 100, t_s))

    resultados[[length(resultados) + 1]] <- data.table(
      Malla      = cfg$label,
      Modelo     = "Solo espacial",
      n_nodos    = malla$n,
      DIC        = modelo_s$dic$dic,
      WAIC       = modelo_s$waic$waic,
      RMSE       = rmse_s,
      MAE        = mae_s,
      Cov95      = round(cov95_s * 100, 1),
      Tiempo_min = round(t_s, 2)
    )
  }

  rm(malla, spde, indice_st, A_st, stk_st, modelo_st,
     indice_s_solo, A_s, stk_s, modelo_s)
  gc()
}

# ==============================================================================
# 7. TABLA COMPARATIVA FINAL
# ==============================================================================

tabla_comparativa <- rbindlist(resultados)
setorder(tabla_comparativa, RMSE)

cat("\n", strrep("=", 80), "\n")
cat("   RESULTADOS DE LA SIMULACIÓN: 3 MALLAS × 2 MODELOS\n")
cat(strrep("=", 80), "\n")
cat(sprintf("   Período: noviembre 2025 | Train: %d días | Test: %d días\n",
            n_train, n_test))
cat(sprintf("   Train: %s → %s\n", min(fechas_train), max(fechas_train)))
cat(sprintf("   Test:  %s → %s\n", min(fechas_test),  max(fechas_test)))
cat(strrep("-", 80), "\n")
print(tabla_comparativa, row.names = FALSE)

write.csv(tabla_comparativa,
          file.path(carpeta_sim, "comparacion_mallas_modelos.csv"),
          row.names = FALSE)

# ==============================================================================
# 8. VISUALIZACIÓN
# ==============================================================================

tabla_comparativa[, Malla := factor(Malla,
  levels = c("Gruesa (8 km)", "Media (4 km)", "Fina (1 km)"))]

# -- Gráfico RMSE --
ggplot(tabla_comparativa, aes(x = Malla, y = RMSE, fill = Modelo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(RMSE, 4)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Espacio-temporal (AR1)" = "#2166AC",
                                "Solo espacial"          = "#B2182B")) +
  labs(
    title    = "Comparación RMSE: Tamaño de Malla × Tipo de Modelo",
    subtitle = sprintf(
      "Test: %s → %s (%d días) | Noviembre 2025 | Madrid",
      min(fechas_test), max(fechas_test), n_test),
    x = "Resolución de malla", y = "RMSE (log NO₂)", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(carpeta_sim, "comparacion_rmse_malla_modelo.png"),
       width = 10, height = 6, dpi = 300)

# -- Gráfico tiempo de cómputo --
ggplot(tabla_comparativa, aes(x = Malla, y = Tiempo_min, fill = Modelo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = paste0(round(Tiempo_min, 1), " min")),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Espacio-temporal (AR1)" = "#2166AC",
                                "Solo espacial"          = "#B2182B")) +
  labs(
    title = "Tiempo de Cómputo: Tamaño de Malla × Tipo de Modelo",
    x     = "Resolución de malla",
    y     = "Tiempo (minutos)",
    fill  = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(carpeta_sim, "comparacion_tiempo_malla_modelo.png"),
       width = 10, height = 6, dpi = 300)

cat("\nSimulación completada. Resultados guardados en:", carpeta_sim, "\n")
