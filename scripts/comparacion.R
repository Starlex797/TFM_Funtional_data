# ==============================================================================
# COMPARACIÓN Y VALIDACIÓN DE MODELOS INLA
# ==============================================================================
# Contenido:
#   FASE 1 — Tablas de efectos fijos (summary.fixed) de ambos modelos.
#   FASE 2 — Validación cruzada espacial por bloques (k = 5 folds de estaciones).
#             Métricas: RMSE, MAE, Cobertura 95% y CRPS.
#   FASE 3 — Tabla comparativa final: DIC, WAIC + métricas de CV.
#
# Referencia metodológica: He & Wong (2021), validación sobre bloques espaciales
# para evaluar predicción en puntos no observados (interpolación espacial).
# ==============================================================================

library(INLA)
library(data.table)
library(sf)
library(here)
library(scoringRules)   # install.packages("scoringRules") si no está instalado

# ------------------------------------------------------------------------------
# 0. CARGAR DATOS Y RECONSTRUIR OBJETOS SPDE
# ------------------------------------------------------------------------------
dt_maestro   <- readRDS(here("data", "processed", "dataset_maestro_inla_2025.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid.rds"))
setDT(dt_maestro)
setorder(dt_maestro, ID_TIEMPO, ESTACION)

# Modelos ya ajustados
modelo_st <- readRDS(here("data", "processed", "modelo_final_no2_madrid.rds"))
modelo_s  <- readRDS(here("data", "processed", "modelo_espacial_no2_madrid.rds"))

# Coordenadas en km (UTM 30N) — necesarias para reconstruir las matrices A en el CV
coords_temp   <- st_as_sf(dt_maestro[, .(LONGITUD, LATITUD)],
                           coords = c("LONGITUD", "LATITUD"), crs = 4326)
coords_utm    <- st_transform(coords_temp, 25830)
coords_puntos <- st_coordinates(coords_utm) / 1000

ndays <- length(unique(dt_maestro$ID_TIEMPO))
spde  <- inla.spde2.matern(mesh = malla_madrid, alpha = 2)

# Objetos espacio-temporal
indice_s   <- inla.spde.make.index("campo_espacial",
                                    n.spde = spde$n.spde, n.group = ndays)
A_espacial <- inla.spde.make.A(mesh  = malla_madrid,
                                loc   = coords_puntos,
                                group = dt_maestro$ID_TIEMPO,
                                n.group = ndays)

# Objetos solo espacial
indice_s_solo <- inla.spde.make.index("campo_espacial_s", n.spde = spde$n.spde)
A_espacial_s  <- inla.spde.make.A(mesh = malla_madrid, loc = coords_puntos)

# Lista de covariables (compartida por ambos modelos)
covariables <- list(
  trafico_intensidad = dt_maestro$intensidad,
  trafico_carga      = dt_maestro$carga,
  temperatura        = dt_maestro$Temperatura,
  viento             = dt_maestro$`Velocidad Viento`,
  precipitacion      = dt_maestro$Precipitaciones,
  humedad            = dt_maestro$Humedad_Relativa,
  id_distrito        = dt_maestro$ID_DISTRITO
)

# Fórmulas de cada modelo
formula_st <- y ~ 0 + intercept +
  trafico_intensidad + trafico_carga +
  temperatura + viento + precipitacion + humedad +
  f(id_distrito, model = "iid") +
  f(campo_espacial, model = spde,
    group = campo_espacial.group,
    control.group = list(model = "ar1"))

formula_s <- y ~ 0 + intercept +
  trafico_intensidad + trafico_carga +
  temperatura + viento + precipitacion + humedad +
  f(id_distrito, model = "iid") +
  f(campo_espacial_s, model = spde)


# ==============================================================================
# FASE 1: TABLAS DE EFECTOS FIJOS (summary.fixed)
# ==============================================================================
cat("\n", strrep("=", 60), "\n")
cat("FASE 1: RESUMEN DE EFECTOS FIJOS\n")
cat(strrep("=", 60), "\n\n")

# Nombres legibles para las covariables
nombres_legibles <- c(
  intercept          = "Intercepto",
  trafico_intensidad = "Intensidad de tráfico (std)",
  trafico_carga      = "Carga de tráfico (std)",
  temperatura        = "Temperatura (ºC)",
  viento             = "Velocidad del viento (m/s)",
  precipitacion      = "Precipitación (mm)",
  humedad            = "Humedad relativa (%)"
)

formatear_tabla_fixed <- function(modelo, nombre_modelo) {
  sf <- as.data.frame(modelo$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")])
  sf$Covariable <- nombres_legibles[rownames(sf)]
  sf$Modelo     <- nombre_modelo
  sf <- sf[, c("Modelo", "Covariable", "mean", "sd", "0.025quant", "0.975quant")]
  colnames(sf) <- c("Modelo", "Covariable", "Media post.", "Desv. típ.", "IC 2.5%", "IC 97.5%")
  rownames(sf) <- NULL
  sf
}

tabla_st <- formatear_tabla_fixed(modelo_st, "Espacio-temporal (AR1)")
tabla_s  <- formatear_tabla_fixed(modelo_s,  "Solo espacial")
tabla_efectos_fijos <- rbind(tabla_st, tabla_s)

cat("Efectos fijos — Espacio-temporal (AR1):\n")
print(round(tabla_st[, 3:6], 4))
cat("\nEfectos fijos — Solo espacial:\n")
print(round(tabla_s[, 3:6], 4))

# Guardar tablas como CSV
fwrite(tabla_efectos_fijos,
       here("output", "tables", "tabla_efectos_fijos.csv"))
cat("\n✅ Tabla de efectos fijos guardada en output/tables/tabla_efectos_fijos.csv\n")


# ==============================================================================
# FASE 2: VALIDACIÓN CRUZADA ESPACIAL POR BLOQUES (k = 5 folds)
# ==============================================================================
cat("\n", strrep("=", 60), "\n")
cat("FASE 2: CROSS-VALIDATION ESPACIAL POR BLOQUES\n")
cat(strrep("=", 60), "\n")
cat("Estrategia: Cada fold oculta todas las observaciones de un subconjunto\n")
cat("de estaciones completas. Evalúa capacidad de interpolación espacial.\n\n")

K_FOLDS <- 5
set.seed(8314)
estaciones  <- unique(dt_maestro$ESTACION)
fold_asig   <- sample(rep(1:K_FOLDS, length.out = length(estaciones)))
fold_df     <- data.frame(ESTACION = estaciones, fold = fold_asig)

# Función para calcular métricas de un fold dado pred_mean, pred_sd (predictiva completa) y obs
calcular_metricas <- function(y_obs, pred_mean, pred_sd) {
  residuos <- y_obs - pred_mean
  rmse <- sqrt(mean(residuos^2, na.rm = TRUE))
  mae  <- mean(abs(residuos),   na.rm = TRUE)
  # Cobertura del IC 95% predictivo
  lower <- pred_mean - 1.96 * pred_sd
  upper <- pred_mean + 1.96 * pred_sd
  cp95  <- mean(y_obs >= lower & y_obs <= upper, na.rm = TRUE) * 100
  # CRPS — evalúa la distribución predictiva completa frente al valor observado
  crps_val <- mean(crps_norm(y = y_obs, mean = pred_mean, sd = pred_sd), na.rm = TRUE)
  c(RMSE = rmse, MAE = mae, CP95 = cp95, CRPS = crps_val)
}

# Función interna: ajusta un modelo INLA con y_cv (NAs en posiciones test)
ajustar_fold <- function(y_cv, A_mat, indice_ef, formula, spde_obj) {
  stk <- inla.stack(
    tag  = "cv",
    data = list(y = y_cv),
    A    = list(A_mat, 1),
    effects = list(
      c(indice_ef, list(intercept = 1)),
      covariables
    ),
    compress = FALSE
  )
  inla(
    formula  = formula,
    data     = inla.stack.data(stk, spde = spde_obj),
    family   = "gaussian",
    control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
    control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE),
    verbose = FALSE
  )
}

# Acumuladores de resultados
resultados_st <- vector("list", K_FOLDS)
resultados_s  <- vector("list", K_FOLDS)

for (k in 1:K_FOLDS) {
  cat(sprintf("  Fold %d / %d ", k, K_FOLDS))

  est_test <- fold_df$ESTACION[fold_df$fold == k]
  idx_test <- which(dt_maestro$ESTACION %in% est_test)
  y_obs    <- dt_maestro$LOG_NO2_DIARIO[idx_test]

  # --- Modelo espacio-temporal ---
  cat("[ST")
  y_cv_st        <- dt_maestro$LOG_NO2_DIARIO
  y_cv_st[idx_test] <- NA
  mod_cv_st      <- ajustar_fold(y_cv_st, A_espacial, indice_s, formula_st, spde)

  prec_st        <- mod_cv_st$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  fv_st          <- mod_cv_st$summary.fitted.values
  pred_mean_st   <- fv_st[idx_test, "mean"]
  pred_sd_st     <- sqrt(fv_st[idx_test, "sd"]^2 + 1 / prec_st)  # sd predictiva completa

  resultados_st[[k]] <- calcular_metricas(y_obs, pred_mean_st, pred_sd_st)
  cat(" ✓]")

  # --- Modelo solo espacial ---
  cat(" [S")
  y_cv_s        <- dt_maestro$LOG_NO2_DIARIO
  y_cv_s[idx_test] <- NA
  mod_cv_s      <- ajustar_fold(y_cv_s, A_espacial_s, indice_s_solo, formula_s, spde)

  prec_s        <- mod_cv_s$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  fv_s          <- mod_cv_s$summary.fitted.values
  pred_mean_s   <- fv_s[idx_test, "mean"]
  pred_sd_s     <- sqrt(fv_s[idx_test, "sd"]^2 + 1 / prec_s)

  resultados_s[[k]]  <- calcular_metricas(y_obs, pred_mean_s, pred_sd_s)
  cat(" ✓]\n")
}

# Promediar métricas sobre los K folds
metricas_st <- colMeans(do.call(rbind, resultados_st))
metricas_s  <- colMeans(do.call(rbind, resultados_s))

cat("\nMétricas CV — Espacio-temporal:\n")
print(round(metricas_st, 4))
cat("Métricas CV — Solo espacial:\n")
print(round(metricas_s, 4))


# ==============================================================================
# FASE 3: TABLA COMPARATIVA FINAL
# ==============================================================================
cat("\n", strrep("=", 60), "\n")
cat("FASE 3: TABLA COMPARATIVA FINAL\n")
cat(strrep("=", 60), "\n\n")

tabla_comparativa <- data.frame(
  Modelo   = c("Solo espacial", "Espacio-temporal (AR1)"),
  # Criterios de información (sobre datos de entrenamiento completos)
  DIC      = c(modelo_s$dic$dic,    modelo_st$dic$dic),
  WAIC     = c(modelo_s$waic$waic,  modelo_st$waic$waic),
  p_eff    = c(modelo_s$dic$p.eff,  modelo_st$dic$p.eff),
  # Métricas de validación cruzada (sobre estaciones no observadas)
  RMSE_cv  = c(metricas_s["RMSE"],  metricas_st["RMSE"]),
  MAE_cv   = c(metricas_s["MAE"],   metricas_st["MAE"]),
  CP95_cv  = c(metricas_s["CP95"],  metricas_st["CP95"]),
  CRPS_cv  = c(metricas_s["CRPS"],  metricas_st["CRPS"])
)

cat("Comparativa de modelos (menor DIC/WAIC/RMSE/CRPS y CP95 ≈ 95% es mejor):\n\n")
print(round(tabla_comparativa[, -1], 3), row.names = FALSE)
cat(sprintf("\n  %-30s %s\n", tabla_comparativa$Modelo[1], "← fila 1"))
cat(sprintf("  %-30s %s\n", tabla_comparativa$Modelo[2], "← fila 2"))

# Guardar tabla comparativa
fwrite(tabla_comparativa,
       here("output", "tables", "tabla_comparativa_modelos.csv"))
cat("\n✅ Tabla comparativa guardada en output/tables/tabla_comparativa_modelos.csv\n")

# Guardar también los resultados por fold (útil para el TFM)
resultados_fold_st <- as.data.frame(do.call(rbind, resultados_st))
resultados_fold_s  <- as.data.frame(do.call(rbind, resultados_s))
resultados_fold_st$fold  <- 1:K_FOLDS
resultados_fold_s$fold   <- 1:K_FOLDS
resultados_fold_st$Modelo <- "Espacio-temporal (AR1)"
resultados_fold_s$Modelo  <- "Solo espacial"

fwrite(rbind(resultados_fold_st, resultados_fold_s),
       here("output", "tables", "cv_metricas_por_fold.csv"))
cat("✅ Métricas por fold guardadas en output/tables/cv_metricas_por_fold.csv\n")
