# ==============================================================================
# PASO 4: EL SÁNDWICH DE DATOS (CONSTRUCCIÓN DEL INLA STACK)
# ==============================================================================

library(INLA)
library(data.table)
library(here)

# 1. Cargar los objetos de los pasos anteriores
dt_maestro   <- readRDS(here("data", "processed", "dataset_maestro_inla_2025_DIARIO.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid_gruesa.rds"))
setDT(dt_maestro)

# Objetos SPDE del Paso 3
spde           <- readRDS(here("data", "processed", "spde_madrid.rds"))
indice_s       <- readRDS(here("data", "processed", "indice_s_madrid.rds"))
A_espacial     <- readRDS(here("data", "processed", "A_espacial_madrid.rds"))
indice_s_solo  <- readRDS(here("data", "processed", "indice_s_solo_madrid.rds"))
A_espacial_s   <- readRDS(here("data", "processed", "A_espacial_s_madrid.rds"))

# Ordenación estricta para mantener sincronía con la Matriz A
setorder(dt_maestro, ID_TIEMPO, ESTACION)

# 2. Definir covariables en un solo lugar (evita duplicación entre stacks)
covariables <- list(
  trafico_intensidad  = dt_maestro$intensidad,
  trafico_carga       = dt_maestro$carga,
  temperatura         = dt_maestro$Temperatura,
  humedad             = dt_maestro$Humedad_Relativa,
  precipitacion       = dt_maestro$Precipitaciones,
  presion_barometrica = dt_maestro$`Presion Barométrica`,
  radiacion_solar     = dt_maestro$`Radiación Solar`,
  id_distrito         = dt_maestro$ID_DISTRITO
)

# 3. Stack espacio-temporal (campo SPDE + AR1)
stk_madrid <- inla.stack(
  tag  = "estimacion_no2",
  data = list(y = dt_maestro$LOG_NO2_DIARIO),
  A    = list(A_espacial, 1),
  effects = list(
    c(indice_s, list(intercept = 1)),
    covariables
  ),
  compress = FALSE
)

# 4. Stack solo espacial (campo SPDE sin dimensión temporal)
stk_madrid_s <- inla.stack(
  tag  = "estimacion_no2_espacial",
  data = list(y = dt_maestro$LOG_NO2_DIARIO),
  A    = list(A_espacial_s, 1),
  effects = list(
    c(indice_s_solo, list(intercept = 1)),
    covariables
  ),
  compress = FALSE
)

# 5. Guardar
saveRDS(stk_madrid,   here("data", "processed", "inla_stack_madrid_2025.rds"))
saveRDS(stk_madrid_s, here("data", "processed", "inla_stack_espacial_2025.rds"))

# ------------------------------------------------------------------------------
# VALIDACIÓN
# ------------------------------------------------------------------------------
n_st <- length(inla.stack.data(stk_madrid)$y)
n_s  <- length(inla.stack.data(stk_madrid_s)$y)

stopifnot("Stack espacio-temporal: nº obs != nrow(dt_maestro)" = n_st == nrow(dt_maestro))
stopifnot("Stack solo espacial: nº obs != nrow(dt_maestro)"    = n_s  == nrow(dt_maestro))

cat("Paso 4 completado.\n")
cat("- Stack espacio-temporal:", n_st, "obs\n")
cat("- Stack solo espacial:   ", n_s,  "obs\n")

