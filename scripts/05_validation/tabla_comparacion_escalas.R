# Actualizar tablas separadas despues de editar K_ELEGIDO.
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(ggplot2)
  library(gridExtra)
})
source(here('scripts', '03_interpolation', 'configuracion_comparacion_clima.R'))
source(here('R', 'interpolation', 'preparacion_comparacion_clima.R'))
source(here('R', 'interpolation', 'comparacion_clima_multiescala.R'))
source(here('R', 'interpolation', 'resumen_comparacion_clima.R'))

# Lee la ultima ejecucion manual completa de las tres escalas.
# No recalcula interpolaciones y no mezcla escalas en una tabla.
exportar_tablas_manuales(here(), K_ELEGIDO)
