# ==============================================================================
# PASO 3: EL PUENTE MATEMÁTICO (SPDE, ÍNDICE Y MATRIZ A)
# ==============================================================================

library(INLA)
library(data.table)
library(here)
library(sf)

# 1. Cargar el dataset maestro y la malla (si estás en una sesión nueva)
dt_maestro   <- readRDS(here("data", "processed", "dataset_maestro_inla_2025.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid.rds"))
setDT(dt_maestro)

# MUY IMPORTANTE: Ordenamos estrictamente por tiempo y luego por estación.
# La Matriz A exige que el orden de las filas coincida exactamente con tus datos.
setorder(dt_maestro, ID_TIEMPO, ESTACION)

# 2. Definir las dimensiones temporales
# ¿Cuántos días únicos hay en tu dataset? (Esto formará nuestro modelo AR1)
ndays <- length(unique(dt_maestro$ID_TIEMPO))

# 3. Extraer las coordenadas de todas las observaciones y pasarlas a Kilómetros
# Primero, le decimos a R cuáles son nuestras coordenadas actuales (en grados WGS84)
coords_temp <- st_as_sf(dt_maestro[, .(LONGITUD, LATITUD)], 
                        coords = c("LONGITUD", "LATITUD"), 
                        crs = 4326)

# Segundo, las proyectamos a metros (UTM 30N)
coords_utm <- st_transform(coords_temp, 25830)

# Tercero, extraemos la matriz de números y la dividimos entre 1000 para tener kilómetros
coords_puntos <- st_coordinates(coords_utm) / 1000 

# (Opcional) Si quieres guardar X_km e Y_km en tu tabla para el futuro:
dt_maestro[, X_km := coords_puntos[, 1]]
dt_maestro[, Y_km := coords_puntos[, 2]]

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

# ------------------------------------------------------------------------------
# CONTROL DE CALIDAD Y COMPROBACIÓN
# ------------------------------------------------------------------------------
cat("✅ ¡Paso 3 completado con éxito!\n")
cat("- Nodos de la malla (Vértices):", malla_madrid$n, "\n")
cat("- Días totales del modelo (Grupos):", ndays, "\n")
cat("- Dimensiones de A_espacial (espacio-temporal):\n")
cat("  -> Filas:", nrow(A_espacial), " | debe coincidir con nrow(dt_maestro):", nrow(dt_maestro), "\n")
cat("  -> Columnas:", ncol(A_espacial), " | debe ser malla_madrid$n * ndays\n")
cat("- Dimensiones de A_espacial_s (solo espacial):\n")
cat("  -> Filas:", nrow(A_espacial_s), " | Columnas:", ncol(A_espacial_s),
    " | debe ser malla_madrid$n:", malla_madrid$n, "\n")
