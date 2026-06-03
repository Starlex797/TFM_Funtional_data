# ==============================================================================
# PASO 5b: MODELO SOLO ESPACIAL (baseline para comparación con espacio-temporal)
# ==============================================================================
# Fórmula matemática:
#   log(NO2) = intercepto
#             + β₁·intensidad + β₂·carga                    (tráfico)
#             + β₃·temp + β₄·viento + β₅·precip + β₆·hum   (clima)
#             + u_i   (IID por distrito)
#             + w(s)  (campo espacial Matérn, SIN AR1 temporal)
#             + ε
#
# Comparación con el modelo espacio-temporal:
#   - Espacio-temporal: w(s,t) con AR1 en t  →  captura persistencia temporal
#   - Solo espacial:    w(s)                 →  campo estático, mismas obs.
#
# Métrica de comparación: DIC, WAIC (ambos modelos los computan)
# ==============================================================================

library(INLA)
library(here)

# 1. Cargar stack y SPDE del Paso 4b
stk_madrid_s <- readRDS(here("data", "processed", "inla_stack_espacial_2025.rds"))
spde_s       <- readRDS(here("data", "processed", "spde_espacial.rds"))

# 2. Fórmula del modelo solo espacial
formula_espacial <- y ~ 0 + intercept +

  # Efectos Fijos
  trafico_intensidad + trafico_carga +
  temperatura + viento + precipitacion + humedad +

  # Efecto Aleatorio IID por distrito
  f(id_distrito, model = "iid") +

  # Campo espacial puro Matérn (SIN group ni AR1)
  f(campo_espacial_s, model = spde_s)

# 3. Ajustar el modelo
cat("Ajustando modelo SOLO ESPACIAL...\n")
t0 <- Sys.time()

modelo_espacial <- inla(
  formula  = formula_espacial,
  data     = inla.stack.data(stk_madrid_s, spde = spde_s),
  family   = "gaussian",
  control.predictor = list(
    A       = inla.stack.A(stk_madrid_s),
    compute = TRUE
  ),
  control.compute = list(
    dic    = TRUE,
    waic   = TRUE,
    cpo    = FALSE,
    config = TRUE
  ),
  verbose = TRUE
)

t1 <- Sys.time()
tiempo_total <- difftime(t1, t0, units = "mins")
cat("\n✅ Modelo espacial completado en", round(as.numeric(tiempo_total), 2), "minutos.\n")

# 4. Guardar
saveRDS(modelo_espacial, here("data", "processed", "modelo_espacial_no2_madrid.rds"))

# 5. Resumen rápido
cat("\n--- Efectos fijos ---\n")
print(round(modelo_espacial$summary.fixed[, c("mean","0.025quant","0.975quant")], 4))

# ------------------------------------------------------------------------------
# COMPARACIÓN DIRECTA DE DIC / WAIC
# ------------------------------------------------------------------------------
# Carga el modelo espacio-temporal si ya está guardado y compara en una tabla
modelo_st_path <- here("data", "processed", "modelo_final_no2_madrid.rds")

if (file.exists(modelo_st_path)) {
  modelo_st <- readRDS(modelo_st_path)
  
  comparacion <- data.frame(
    Modelo = c("Solo espacial", "Espacio-temporal (AR1)"),
    DIC    = c(modelo_espacial$dic$dic,  modelo_st$dic$dic),
    WAIC   = c(modelo_espacial$waic$waic, modelo_st$waic$waic),
    p.eff  = c(modelo_espacial$dic$p.eff, modelo_st$dic$p.eff)
  )
  
  cat("\n===== COMPARACIÓN DE MODELOS =====\n")
  print(comparacion)
  cat("  → Menor DIC/WAIC indica mejor balance ajuste-complejidad.\n")
} else {
  cat("\nℹ️  Modelo espacio-temporal no encontrado en disco. Ejecuta Paso_5 primero para comparar.\n")
}
