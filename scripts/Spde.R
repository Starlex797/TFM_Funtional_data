# ==============================================================================
# Step 3: SPDE Model Setup
# ==============================================================================

library(INLA)
library(data.table)
library(here)
library(sf)

# 1. Cargar el dataset maestro y la malla (si estás en una sesión nueva)
dt_maestro   <- readRDS(here("data", "processed", "dataset_maestro_inla_2025_DIARIO.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid_gruesa.rds"))
setDT(dt_maestro)

# MUY IMPORTANTE: Ordenamos estrictamente por tiempo y luego por estación.
# La Matriz A exige que el orden de las filas coincida exactamente con tus datos.
setorder(dt_maestro, ID_TIEMPO, ESTACION)

# 2. Definir las dimensiones temporales
# ¿Cuántos días únicos hay en tu dataset? (Esto formará nuestro modelo AR1)
ndays <- length(unique(dt_maestro$ID_TIEMPO))

# 3. Proyectar coordenadas a kilómetros (UTM 30N)
# Solo proyectamos las 24 estaciones únicas y luego pegamos al maestro
coords_unicas <- unique(dt_maestro[, .(ESTACION, LONGITUD, LATITUD)])
coords_sf <- st_as_sf(coords_unicas, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)
coords_unicas[, X_km := st_coordinates(coords_sf)[, 1] / 1000]
coords_unicas[, Y_km := st_coordinates(coords_sf)[, 2] / 1000]

dt_maestro <- merge(dt_maestro, coords_unicas[, .(ESTACION, X_km, Y_km)],
                    by = "ESTACION", all.x = TRUE)
setorder(dt_maestro, ID_TIEMPO, ESTACION)

coords_puntos <- as.matrix(dt_maestro[, .(X_km, Y_km)])

# 4. Crear el objeto SPDE (El modelo matemático espacial de Matérn)
# alpha = 2 es el estándar en R-INLA para superficies en 2D
spde <- inla.spde2.matern(mesh = malla_madrid, alpha = 2)

# 5. Crear el Índice Espacio-Temporal
# Esto le dice a INLA: "Crea una variable latente espacial y multiplícala por los 'ndays'"

indice_s <- inla.spde.make.index(name = "campo_espacial", 
                                 n.spde = spde$n.spde, 
                                 n.group = ndays)

# 6. Crear la Matriz de Proyección (Matriz A) para los Puntos
# - loc: las coordenadas de cada fila
# - group: a qué día (1 a 365) pertenece esa coordenada
A_espacial <- inla.spde.make.A(mesh = malla_madrid, # Mi malla 
                               loc = coords_puntos, # Dónde están las estaciones 
                               group = dt_maestro$ID_TIEMPO, # En qué dia ocurrió la medicion 
                               n.group = ndays) # Total de dias 

# ==============================================================================
# MODELO SOLO ESPACIAL: Índice y Matriz A sin dimensión temporal
# Sin 'n.group' → un único campo latente espacial compartido por todos los días.
# Sin 'group'   → cada observación se proyecta solo en el espacio (sin AR1).

indice_s_solo <- inla.spde.make.index(
  name   = "campo_espacial_s",
  n.spde = spde$n.spde
  # n.group = 1 es el defecto → sin replicación temporal
)

A_espacial_s <- inla.spde.make.A(
  mesh = malla_madrid,
  loc  = coords_puntos   # Sin 'group': proyección puramente espacial
)

# 7. GUARDAR LOS OBJETOS PARA EL PASO 4
# Guardamos los objetos del modelo espacio-temporal
saveRDS(spde,       here("data", "processed", "spde_madrid.rds"))
saveRDS(indice_s,   here("data", "processed", "indice_s_madrid.rds"))
saveRDS(A_espacial, here("data", "processed", "A_espacial_madrid.rds"))

# Guardamos los objetos del modelo solo espacial (baseline)
saveRDS(indice_s_solo, here("data", "processed", "indice_s_solo_madrid.rds"))
saveRDS(A_espacial_s,  here("data", "processed", "A_espacial_s_madrid.rds"))

# ------------------------------------------------------------------------------
# VALIDACIÓN (el script se detiene si algo no cuadra)
# ------------------------------------------------------------------------------
stopifnot("Filas de A_espacial != filas de dt_maestro" =
            nrow(A_espacial) == nrow(dt_maestro))
stopifnot("Columnas de A_espacial != nodos * días" =
            ncol(A_espacial) == malla_madrid$n * ndays)
stopifnot("Filas de A_espacial_s != filas de dt_maestro" =
            nrow(A_espacial_s) == nrow(dt_maestro))
stopifnot("Columnas de A_espacial_s != nodos de la malla" =
            ncol(A_espacial_s) == malla_madrid$n)

cat("Paso 3 completado.\n")
cat("- Nodos de la malla:", malla_madrid$n, "\n")
cat("- Días (grupos):", ndays, "\n")
cat("- A_espacial:", nrow(A_espacial), "x", ncol(A_espacial), "\n")
cat("- A_espacial_s:", nrow(A_espacial_s), "x", ncol(A_espacial_s), "\n")

