# ==============================================================================
# PASO 4b: STACK PARA EL MODELO SOLO ESPACIAL (baseline de comparación)
# ==============================================================================
# Diferencia clave respecto al Paso 3-4 (espacio-temporal):
#   - La Matriz A NO lleva argumento 'group' → sin dimensión temporal en el SPDE.
#   - El índice SPDE tiene n.group = 1 → un único campo espacial estático.
#   - Las observaciones son las mismas (todos los días × estaciones), pero el 
#     campo latente no varía en el tiempo: solo captura correlación espacial.
# ==============================================================================

library(INLA)
library(data.table)
library(sf)
library(here)

# 1. Cargar datos y malla
dt_maestro   <- readRDS(here("data", "processed", "dataset_maestro_inla_2025.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid.rds"))
setDT(dt_maestro)

# Orden estricto (igual que en el modelo espacio-temporal)
setorder(dt_maestro, ID_TIEMPO, ESTACION)

# 2. Coordenadas en kilómetros (UTM 30N)
coords_temp  <- st_as_sf(dt_maestro[, .(LONGITUD, LATITUD)],
                         coords = c("LONGITUD", "LATITUD"), crs = 4326)
coords_utm   <- st_transform(coords_temp, 25830)
coords_puntos <- st_coordinates(coords_utm) / 1000

# 3. SPDE Matérn (idéntico al espacio-temporal)
spde_s <- inla.spde2.matern(mesh = malla_madrid, alpha = 2)

# 4. Índice ESPACIAL PURO: sin n.group → un único campo latente compartido por todos los días
indice_s_solo <- inla.spde.make.index(
  name  = "campo_espacial_s",
  n.spde = spde_s$n.spde
  # n.group = 1 es el valor por defecto → sin replicación temporal
)

# 5. Matriz A ESPACIAL PURA: sin 'group' → solo proyección espacial
A_espacial_s <- inla.spde.make.A(
  mesh = malla_madrid,
  loc  = coords_puntos   # Sin group: cada observación se proyecta solo en el espacio
)

# 6. Construir el Stack
stk_madrid_s <- inla.stack(
  tag  = "estimacion_no2_espacial",
  data = list(y = dt_maestro$LOG_NO2_DIARIO),
  A    = list(A_espacial_s, 1),
  effects = list(
    # Efecto 1: Campo espacial puro + intercepto
    c(indice_s_solo, list(intercept = 1)),
    # Efecto 2: Covariables fijas y efecto IID de distrito
    list(
      trafico_intensidad = dt_maestro$intensidad,
      trafico_carga      = dt_maestro$carga,
      temperatura        = dt_maestro$Temperatura,
      viento             = dt_maestro$`Velocidad Viento`,
      precipitacion      = dt_maestro$Precipitaciones,
      humedad            = dt_maestro$Humedad_Relativa,
      id_distrito        = dt_maestro$ID_DISTRITO
    )
  ),
  compress = FALSE
)

# 7. Guardar stack y objetos necesarios para el Paso 5b
saveRDS(stk_madrid_s,  here("data", "processed", "inla_stack_espacial_2025.rds"))
saveRDS(spde_s,        here("data", "processed", "spde_espacial.rds"))

# ------------------------------------------------------------------------------
# CONTROL DE CALIDAD
# ------------------------------------------------------------------------------
cat("✅ ¡Paso 4b completado!\n")
cat("- Filas en el stack:", inla.stack.ndata(stk_madrid_s), "\n")
cat("- Columnas de A_espacial_s:", ncol(A_espacial_s),
    "== nodos de la malla:", malla_madrid$n, "\n")
cat("  (Debe ser igual; sin group no hay replicación temporal)\n")
