# ==============================================================================
# PASO 5: LA ECUACIÓN ESPACIOTEMPORAL Y EL MOTOR INLA
# ==============================================================================

library(INLA)
library(data.table)
library(here)

# 1. Cargar dependencias de los pasos anteriores
dt_maestro   <- readRDS(here("data", "processed", "dataset_maestro_inla_2025_DIARIO.rds"))
stk_madrid   <- readRDS(here("data", "processed", "inla_stack_madrid_2025.rds"))
stk_madrid_s <- readRDS(here("data", "processed", "inla_stack_espacial_2025.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid_gruesa.rds"))
setDT(dt_maestro)
setorder(dt_maestro, ID_TIEMPO, ESTACION)

spde <- inla.spde2.matern(mesh = malla_madrid, alpha = 2)

# ==============================================================================
# 3. FÓRMULA COMPARTIDA
# Definida una sola vez para garantizar que ambos modelos son comparables.
# La única diferencia entre ellos es el campo latente (con o sin AR1).
# ==============================================================================

efectos_fijos <- y ~ 0 + intercept +
  trafico_intensidad  +
  temperatura + humedad + precipitacion +
  presion_barometrica + radiacion_solar 

formula_st <- update(efectos_fijos,
  . ~ . + f(campo_espacial,
             model = spde,
             group = campo_espacial.group,
             control.group = list(model = "ar1")))

formula_s <- update(efectos_fijos,
  . ~ . + f(campo_espacial_s, model = spde))

# ==============================================================================
# 4. MODELO ESPACIO-TEMPORAL (SPDE + AR1)
# ==============================================================================

cat("\nAjustando modelo espacio-temporal...\n")
t0 <- Sys.time()

modelo_st <- inla(
  formula = formula_st,
  data    = inla.stack.data(stk_madrid, spde = spde),
  family  = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_madrid), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE, config = TRUE),
  control.inla      = list(strategy = "laplace"),
  verbose = FALSE
)

cat("Completado en", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "minutos.\n")
saveRDS(modelo_st, here("data", "processed", "modelo_st_no2_madrid.rds"))

cat("\n--- Efectos fijos (espacio-temporal) ---\n")
print(round(modelo_st$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ==============================================================================
# 5. MODELO SOLO ESPACIAL (baseline)
# ==============================================================================

cat("\nAjustando modelo solo espacial...\n")
t0_s <- Sys.time()

modelo_s <- inla(
  formula = formula_s,
  data    = inla.stack.data(stk_madrid_s, spde = spde),
  family  = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_madrid_s), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE, config = TRUE),
  control.inla      = list(strategy = "laplace"),
  verbose = FALSE
)

cat("Completado en", round(as.numeric(difftime(Sys.time(), t0_s, units = "mins")), 2), "minutos.\n")
saveRDS(modelo_s, here("data", "processed", "modelo_espacial_no2_madrid.rds"))

cat("\n--- Efectos fijos (solo espacial) ---\n")
print(round(modelo_s$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ==============================================================================
# 6. COMPARACIÓN DE MODELOS
# ==============================================================================

# Índices de observación reales dentro del stack (excluye filas internas de la Matriz A)
idx_st <- inla.stack.index(stk_madrid,   tag = "estimacion_no2")$data
idx_s  <- inla.stack.index(stk_madrid_s, tag = "estimacion_no2_espacial")$data

y_obs <- dt_maestro$LOG_NO2_DIARIO

comparacion <- data.frame(
  Modelo = c("Solo espacial", "Espacio-temporal (AR1)"),
  DIC    = c(modelo_s$dic$dic,   modelo_st$dic$dic),
  WAIC   = c(modelo_s$waic$waic, modelo_st$waic$waic),
  p.eff  = c(modelo_s$dic$p.eff, modelo_st$dic$p.eff),
  RMSE   = c(
    sqrt(mean((modelo_s$summary.fitted.values$mean[idx_s]   - y_obs)^2, na.rm = TRUE)),
    sqrt(mean((modelo_st$summary.fitted.values$mean[idx_st] - y_obs)^2, na.rm = TRUE))
  )
)

cat("\n===== COMPARACIÓN DE MODELOS =====\n")
print(round(comparacion[, -1], 3))
cat("Menor DIC/WAIC indica mejor balance entre ajuste y complejidad.\n")

saveRDS(comparacion, here("data", "processed", "comparacion_modelos.rds"))
write.csv(comparacion, here("output", "comparacion_modelos.csv"), row.names = FALSE)
