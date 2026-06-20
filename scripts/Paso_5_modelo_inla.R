# ==============================================================================
# PASO 5: LA ECUACIÓN ESPACIOTEMPORAL Y EL MOTOR INLA
# ==============================================================================
# NOTA: Al inicio se construye un dataset de "día tipo" (7 grupos, uno por día
# de la semana) para reducir observaciones y medir el tiempo de cómputo del
# modelo INLA-SPDE antes de lanzar el año completo.
# ==============================================================================

library(INLA)
library(data.table)
library(here)
library(sf)

# ==============================================================================
# 0. CONSTRUCCIÓN DEL DATASET "DÍA TIPO"
#    Se promedian todas las observaciones de cada lunes, martes, ... para cada
#    estación. Resultado: 24 estaciones × 7 días = 168 filas.
# ==============================================================================

cat("Construyendo dataset de día tipo...\n")

dt_diario <- readRDS(here("data", "processed", "dataset_maestro_inla_2025_DIARIO.rds"))
setDT(dt_diario)

# Día de la semana ISO: 1 = lunes ... 7 = domingo (base R, sin lubridate)
dt_diario[, dia_semana := as.integer(format(FECHA, "%u"))]

# Columnas numéricas a promediar
cols_num <- c("DATO_DIARIO", "LOG_NO2_DIARIO",
              "intensidad", "ocupacion", "carga",
              "Temperatura", "Humedad_Relativa", "Precipitaciones",
              "Presion Barométrica", "Radiación Solar",
              "intensidad_raw", "carga_raw",
              "Temperatura_raw", "Humedad_Relativa_raw", "Precipitaciones_raw",
              "Presion Barométrica_raw", "Radiación Solar_raw")

# Columnas fijas por estación (constantes entre días)
cols_fijas <- c("LONGITUD", "LATITUD", "barrio", "distrito", "ID_DISTRITO")

# Agregar: media numérica + primera ocurrencia de fijas
dt_tipo <- dt_diario[,
  c(lapply(.SD, mean, na.rm = TRUE),
    lapply(mget(cols_fijas), function(x) x[1])),
  by  = .(ESTACION, dia_semana),
  .SDcols = cols_num
]

# ID_TIEMPO = día de semana (1 = lun, 7 = dom)
dt_tipo[, ID_TIEMPO := dia_semana]
setorder(dt_tipo, ID_TIEMPO, ESTACION)

n_estaciones <- uniqueN(dt_tipo$ESTACION)
n_dias_tipo  <- uniqueN(dt_tipo$ID_TIEMPO)
cat(sprintf("  Estaciones: %d  |  Días tipo: %d  |  Total filas: %d\n",
            n_estaciones, n_dias_tipo, nrow(dt_tipo)))

# ==============================================================================
# 1. MALLA + SPDE + MATRICES A  (reconstruidas con ndays = 7)
# ==============================================================================

malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid_gruesa.rds"))
spde         <- inla.spde2.matern(mesh = malla_madrid, alpha = 2)

# Proyectar coordenadas a km (UTM 30N)
coords_unicas <- unique(dt_tipo[, .(ESTACION, LONGITUD, LATITUD)])
coords_sf     <- st_as_sf(coords_unicas,
                          coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)
coords_unicas[, X_km := st_coordinates(coords_sf)[, 1] / 1000]
coords_unicas[, Y_km := st_coordinates(coords_sf)[, 2] / 1000]

dt_tipo <- merge(dt_tipo, coords_unicas[, .(ESTACION, X_km, Y_km)],
                 by = "ESTACION", all.x = TRUE)
setorder(dt_tipo, ID_TIEMPO, ESTACION)

coords_puntos <- as.matrix(dt_tipo[, .(X_km, Y_km)])
ndays <- n_dias_tipo  # 7

# --- Espacio-temporal (AR1 sobre los 7 días tipo) ---
indice_s <- inla.spde.make.index(name    = "campo_espacial",
                                 n.spde  = spde$n.spde,
                                 n.group = ndays)

A_espacial <- inla.spde.make.A(mesh    = malla_madrid,
                               loc     = coords_puntos,
                               group   = dt_tipo$ID_TIEMPO,
                               n.group = ndays)

# --- Solo espacial (baseline, sin AR1) ---
indice_s_solo <- inla.spde.make.index(name   = "campo_espacial_s",
                                      n.spde = spde$n.spde)

A_espacial_s  <- inla.spde.make.A(mesh = malla_madrid, loc = coords_puntos)

cat(sprintf("  A espacio-temporal: %d x %d\n", nrow(A_espacial), ncol(A_espacial)))
cat(sprintf("  A solo espacial:    %d x %d\n", nrow(A_espacial_s), ncol(A_espacial_s)))

# ==============================================================================
# 2. INLA STACKS CON EL DATASET DÍA TIPO
# ==============================================================================

covariables <- list(
  trafico_intensidad  = dt_tipo$intensidad,
  trafico_carga       = dt_tipo$carga,
  temperatura         = dt_tipo$Temperatura,
  humedad             = dt_tipo$Humedad_Relativa,
  precipitacion       = dt_tipo$Precipitaciones,
  presion_barometrica = dt_tipo$`Presion Barométrica`,
  radiacion_solar     = dt_tipo$`Radiación Solar`,
  id_distrito         = dt_tipo$ID_DISTRITO
)

stk_madrid <- inla.stack(
  tag  = "estimacion_no2",
  data = list(y = dt_tipo$LOG_NO2_DIARIO),
  A    = list(A_espacial, 1),
  effects = list(
    c(indice_s, list(intercept = 1)),
    covariables
  ),
  compress = FALSE
)

stk_madrid_s <- inla.stack(
  tag  = "estimacion_no2_espacial",
  data = list(y = dt_tipo$LOG_NO2_DIARIO),
  A    = list(A_espacial_s, 1),
  effects = list(
    c(indice_s_solo, list(intercept = 1)),
    covariables
  ),
  compress = FALSE
)
stk_madrid_s
cat("Stacks construidos con datos de día tipo.\n")

# ==============================================================================
# 3. FÓRMULA COMPARTIDA
# ==============================================================================

efectos_fijos <- y ~ 0 + intercept +
  trafico_intensidad  +
  temperatura  + precipitacion 

formula_st <- update(efectos_fijos,
  . ~ . + f(campo_espacial,
             model = spde,
             group = campo_espacial.group,
             control.group = list(model = "ar1")))

formula_s <- update(efectos_fijos,
  . ~ . + f(campo_espacial_s, model = spde))

# ==============================================================================
# 4. MODELO ESPACIO-TEMPORAL (SPDE + AR1) — DÍA TIPO
# ==============================================================================

cat("\nAjustando modelo espacio-temporal (día tipo)...\n")
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
saveRDS(modelo_st, here("data", "processed", "modelo_st_no2_madrid_dia_tipo.rds"))

cat("\n--- Efectos fijos (espacio-temporal, día tipo) ---\n")
print(round(modelo_st$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ==============================================================================
# 5. MODELO SOLO ESPACIAL (baseline) — DÍA TIPO
# ==============================================================================

cat("\nAjustando modelo solo espacial (día tipo)...\n")
t0_s <- Sys.time()

modelo_s <- inla(
  formula = formula_s,
  data    = inla.stack.data(stk_madrid_s, spde = spde),
  family  = "gaussian ",
  control.predictor = list(A = inla.stack.A(stk_madrid_s), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE, config = TRUE),
  control.inla      = list(strategy = "laplace"),
  verbose = FALSE
)

cat("Completado en", round(as.numeric(difftime(Sys.time(), t0_s, units = "mins")), 2), "minutos.\n")
saveRDS(modelo_s, here("data", "processed", "modelo_espacial_no2_madrid_dia_tipo.rds"))

cat("\n--- Efectos fijos (solo espacial, día tipo) ---\n")
print(round(modelo_s$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ==============================================================================
# 6. COMPARACIÓN DE MODELOS
# ==============================================================================

idx_st <- inla.stack.index(stk_madrid,   tag = "estimacion_no2")$data
idx_s  <- inla.stack.index(stk_madrid_s, tag = "estimacion_no2_espacial")$data

y_obs <- dt_tipo$LOG_NO2_DIARIO

comparacion <- data.frame(
  Modelo = c("Solo espacial (día tipo)", "Espacio-temporal AR1 (día tipo)"),
  DIC    = c(modelo_s$dic$dic,   modelo_st$dic$dic),
  WAIC   = c(modelo_s$waic$waic, modelo_st$waic$waic),
  p.eff  = c(modelo_s$dic$p.eff, modelo_st$dic$p.eff),
  RMSE   = c(
    sqrt(mean((modelo_s$summary.fitted.values$mean[idx_s]   - y_obs)^2, na.rm = TRUE)),
    sqrt(mean((modelo_st$summary.fitted.values$mean[idx_st] - y_obs)^2, na.rm = TRUE))
  )
)

cat("\n===== COMPARACIÓN DE MODELOS (DÍA TIPO) =====\n")
print(round(comparacion[, -1], 3))
cat("Menor DIC/WAIC indica mejor balance entre ajuste y complejidad.\n")

saveRDS(comparacion, here("data", "processed", "comparacion_modelos_dia_tipo.rds"))
write.csv(comparacion, here("output", "comparacion_modelos_dia_tipo.csv"), row.names = FALSE)






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
  presion_barometrica + velocidad_viento 

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

modelo_st2 <- inla(
  formula = formula_st,
  data    = inla.stack.data(stk_madrid, spde = spde),
  family  = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_madrid), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE, config = TRUE),
  control.inla      = list(strategy = "laplace"),
  verbose = FALSE
)

cat("Completado en", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "minutos.\n")
saveRDS(modelo_st2, here("data", "processed", "modelo_st2_no2_madrid.rds"))

cat("\n--- Efectos fijos (espacio-temporal) ---\n")
print(round(modelo_st$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ==============================================================================
# 5. MODELO SOLO ESPACIAL (baseline)
# ==============================================================================

cat("\nAjustando modelo solo espacial...\n")
t0_s <- Sys.time()

modelo_s2 <- inla(
  formula = formula_s,
  data    = inla.stack.data(stk_madrid_s, spde = spde),
  family  = "gaussian",
  control.predictor = list(A = inla.stack.A(stk_madrid_s), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE, config = TRUE),
  control.inla      = list(strategy = "laplace"),
  verbose = FALSE
)

cat("Completado en", round(as.numeric(difftime(Sys.time(), t0_s, units = "mins")), 2), "minutos.\n")
saveRDS(modelo_s2, here("data", "processed", "modelo_espacial2_no2_madrid.rds"))

cat("\n--- Efectos fijos (solo espacial) ---\n")
print(round(modelo_s2$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 4))

# ==============================================================================
# 6. COMPARACIÓN DE MODELOS
# ==============================================================================

# Índices de observación reales dentro del stack (excluye filas internas de la Matriz A)
idx_st2 <- inla.stack.index(stk_madrid,   tag = "estimacion_no2")$data
idx_s2  <- inla.stack.index(stk_madrid_s, tag = "estimacion_no2_espacial")$data

y_obs <- dt_maestro$LOG_NO2_DIARIO

comparacion <- data.frame(
  Modelo = c("Solo espacial", "Espacio-temporal (AR1)"),
  DIC    = c(modelo_s2$dic$dic,   modelo_st2$dic$dic),
  WAIC   = c(modelo_s2$waic$waic, modelo_st2$waic$waic),
  p.eff  = c(modelo_s2$dic$p.eff, modelo_st2$dic$p.eff),
  RMSE   = c(
    sqrt(mean((modelo_s2$summary.fitted.values$mean[idx_s2]   - y_obs)^2, na.rm = TRUE)),
    sqrt(mean((modelo_st2$summary.fitted.values$mean[idx_st2] - y_obs)^2, na.rm = TRUE))
  )
)

cat("\n===== COMPARACIÓN DE MODELOS =====\n")
print(round(comparacion[, -1], 3))
cat("Menor DIC/WAIC indica mejor balance entre ajuste y complejidad.\n")

saveRDS(comparacion, here("data", "processed", "comparacion_modelos.rds"))
write.csv(comparacion, here("output", "comparacion_modelos.csv"), row.names = FALSE)
º