# ==============================================================================
# PASO 5: LA ECUACIÓN ESPACIOTEMPORAL Y EL MOTOR INLA
# ==============================================================================

library(INLA)
library(here)

# 1. Cargar dependencias de los pasos anteriores
stk_madrid   <- readRDS(here("data", "processed", "inla_stack_madrid_2025.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid.rds"))

# Recreamos el SPDE el campo aleatorio gaussiano (es instantáneo y nos asegura tenerlo en el entorno)
spde <- inla.spde2.matern(mesh = malla_madrid, alpha = 2)

# 2. LA FÓRMULA MATEMÁTICA
# La 'y' representa tu LOG_NO2_DIARIO
# '0 + intercept' apaga el intercepto automático de R y usa el nuestro del stack
formula_no2 <- y ~ 0 + intercept + 
  
  # Efectos Fijos (Covariables)
  trafico_intensidad + trafico_carga + 
  temperatura + viento + precipitacion + humedad + 
  
  # Efecto Aleatorio 1: Ruido no estructurado por distrito (IID)
  f(id_distrito, model = "iid") + 
  
  # Efecto Aleatorio 2: El Campo Espaciotemporal Continuo (SPDE + AR1)
  # 'group' enlaza con la dimensión temporal y 'model="ar1"' le da la persistencia
  f(campo_espacial, 
    model = spde, 
    group = campo_espacial.group, 
    control.group = list(model = "ar1"))

# 3. EJECUCIÓN DEL MODELO INLA
cat("Arrancando el motor bayesiano de INLA...\n")
cat("Ve a por un café, esto puede tardar varios minutos dependiendo de tu PC ☕\n")

t0 <- Sys.time()

modelo_final <- inla(
  formula = formula_no2,
  data = inla.stack.data(stk_madrid, spde = spde), # Extrae los datos del stack
  family = "gaussian",                             # Asumimos distribución normal para el Log(NO2)
  control.predictor = list(A = inla.stack.A(stk_madrid), compute = TRUE), # Proyecta las predicciones
  control.compute = list(dic = TRUE, waic = TRUE, cpo = FALSE, config = TRUE), # Métricas de bondad de ajuste
  verbose = TRUE # Activa los mensajes en consola para ver cómo avanza el algoritmo
)

t1 <- Sys.time()

# 4. GUARDADO Y RESUMEN
tiempo_total <- difftime(t1, t0, units = "mins")
cat("\n✅ ¡Modelo espacio-temporal completado en", round(as.numeric(tiempo_total), 2), "minutos!\n")

saveRDS(modelo_final, here("data", "processed", "modelo_final_no2_madrid.rds"))

cat("\n--- Efectos fijos (espacio-temporal) ---\n")
print(round(modelo_final$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ==============================================================================
# MODELO SOLO ESPACIAL

stk_madrid_s <- readRDS(here("data", "processed", "inla_stack_espacial_2025.rds"))

# Fórmula idéntica salvo que el campo latente no lleva 'group' ni AR1
formula_espacial <- y ~ 0 + intercept +
  trafico_intensidad +
  temperatura + viento + precipitacion  +
  f(campo_espacial_s, model = spde)   # Campo Matérn estático, sin AR1

cat("\nAjustando modelo SOLO ESPACIAL...\n")
t0_s <- Sys.time()

modelo_espacial <- inla(
  formula  = formula_espacial,
  data     = inla.stack.data(stk_madrid_s, spde = spde),
  family   = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_madrid_s), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE, config = TRUE),
  verbose = TRUE
)

t1_s <- Sys.time()
cat("\n✅ Modelo espacial completado en",
    round(as.numeric(difftime(t1_s, t0_s, units = "mins")), 2), "minutos!\n")

saveRDS(modelo_espacial, here("data", "processed", "modelo_espacial_no2_madrid.rds"))

cat("\n--- Efectos fijos (solo espacial) ---\n")
print(round(modelo_espacial$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ==============================================================================
# COMPARACIÓN DE MODELOS: DIC y WAIC

comparacion <- data.frame(
  Modelo = c("Solo espacial", "Espacio-temporal (AR1)"),
  DIC    = c(modelo_espacial$dic$dic,   modelo_final$dic$dic),
  WAIC   = c(modelo_espacial$waic$waic, modelo_final$waic$waic),
  p.eff  = c(modelo_espacial$dic$p.eff, modelo_final$dic$p.eff), 
  RMSE   = c(
    sqrt(mean((modelo_espacial$summary.fitted.values$mean - inla.stack.data(stk_madrid_s)$y)^2)),
    sqrt(mean((modelo_final$summary.fitted.values$mean - inla.stack.data(stk_madrid)$y)^2))
  )
)

cat("\n===== COMPARACIÓN DE MODELOS =====\n")
print(comparacion)
cat("→ Menor DIC/WAIC indica mejor balance entre ajuste y complejidad.\n")

