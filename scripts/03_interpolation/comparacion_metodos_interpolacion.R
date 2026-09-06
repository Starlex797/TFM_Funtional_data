# Comparacion de interpoladores: primero mirar las curvas, despues elegir k.
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(ggplot2)
  library(gridExtra)
})

# Paso 1: leer configuracion y funciones.
source(here("scripts", "03_interpolation", "configuracion_comparacion_clima.R"))
source(here("R", "interpolation", "preparacion_comparacion_clima.R"))
source(here("R", "interpolation", "comparacion_clima_multiescala.R"))
source(here("R", "interpolation", "resumen_comparacion_clima.R"))

# Paso 2: calcular las curvas RMSE-k para las tres escalas.
DIRECTORIO_RESULTADOS <- ejecutar_sensibilidad_manual(
  raiz = here(), escalas = c("horario", "diario", "mensual"),
  anio = ANIO, cobertura_minima = COBERTURA_MINIMA,
  umbral_lluvia = UMBRAL_LLUVIA, k_elegido = K_ELEGIDO
)

# Paso 3: una tabla por escala, con los k que TU hayas indicado.
# Si k es NA, la fila queda pendiente; no se elige un minimo automaticamente.
exportar_tablas_manuales(here(), K_ELEGIDO, DIRECTORIO_RESULTADOS)

# Despues de ver los graficos:
# 1. Editar K_ELEGIDO en configuracion_comparacion_clima.R.
# 2. Ejecutar scripts/05_validation/tabla_comparacion_escalas.R.
# Ese segundo script actualiza las tablas sin repetir la interpolacion.
