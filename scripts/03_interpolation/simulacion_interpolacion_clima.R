# Acceso a la sensibilidad RMSE-k de escala diario.
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

DIRECTORIO_RESULTADOS <- ejecutar_sensibilidad_manual(
  raiz = here(), escalas = 'diario', anio = ANIO,
  cobertura_minima = COBERTURA_MINIMA, umbral_lluvia = UMBRAL_LLUVIA,
  k_elegido = K_ELEGIDO
)
exportar_tablas_manuales(here(), K_ELEGIDO, DIRECTORIO_RESULTADOS)
