# ==============================================================================
# ANALISIS HORARIO: DIAS AUSENTES Y HORAS FALTANTES POR MES Y AÑO (2019-2025)
# ==============================================================================
# Trabaja sobre los datos HORARIOS ya preprocesados y responde:
#
#   (1) ¿Cuantos dias (a nivel estacion-dia) NO tienen ningun valor en el mes?
#       -> "dias ausentes": las 24 horas de esa estacion ese dia estan vacias.
#
#   (2) De los dias que si tienen datos, ¿en cuantos faltan algunas horas y en
#       que porcentaje? -> "dias parciales" (1-23 horas) y el % medio de horas
#       que faltan en ellos.
#
# Unidad de analisis: ESTACION-DIA. Un "dia" aqui es una estacion concreta en
# una fecha concreta (24 horas esperadas). Se agrega despues por AÑO x MES.
#
# Fuentes (post-preprocesamiento):
#   - NO2:   data/processed/Contaminacion/horario/aire_madrid_<a>_No2_horarios.rds
#   - Clima: data/processed/Clima/horario/meteo_madrid_<a>_horario.rds
#
# Nota sobre el clima: cada variable se evalua SOLO sobre las estaciones que la
# miden (las que no tienen sensor se excluyen para no confundir la ausencia
# estructural con un hueco temporal).
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

ANIOS <- 2019:2025
DIR_TABLAS <- here("outputs", "tables", "calidad_datos")
dir.create(DIR_TABLAS, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Funcion nucleo: dado un data.table largo con (ESTACION, FECHA, valor),
# resume por AÑO x MES los dias ausentes, parciales y completos.
# ------------------------------------------------------------------------------
resumir_horario <- function(dt_largo, anio, estaciones_validas = NULL) {

  dt <- copy(dt_largo)
  setDT(dt)
  dt[, FECHA := as.Date(FECHA)]

  # Estaciones a considerar (para clima, solo las que miden la variable)
  if (is.null(estaciones_validas)) {
    estaciones_validas <- dt[!is.na(valor), unique(ESTACION)]
  }
  dt <- dt[ESTACION %in% estaciones_validas]

  # Horas presentes por estacion-dia (sobre las filas existentes)
  por_ed <- dt[, .(h_present = sum(!is.na(valor))), by = .(ESTACION, FECHA)]

  # Rejilla COMPLETA estacion x dia (para que los dias sin ninguna fila cuenten)
  dias_anio <- seq(as.Date(sprintf("%d-01-01", anio)),
                   as.Date(sprintf("%d-12-31", anio)), by = "day")
  grid <- CJ(ESTACION = estaciones_validas, FECHA = dias_anio)
  ed <- por_ed[grid, on = c("ESTACION", "FECHA")]
  ed[is.na(h_present), h_present := 0L]

  # Clasificacion de cada estacion-dia
  ed[, `:=`(
    ANIO = anio,
    MES  = month(FECHA),
    h_faltan = 24L - h_present,
    tipo = fifelse(h_present == 0L, "ausente",
            fifelse(h_present == 24L, "completo", "parcial"))
  )]

  # Resumen por AÑO x MES
  resumen <- ed[, .(
    n_estaciones          = uniqueN(ESTACION),
    dias_mes              = uniqueN(FECHA),
    est_dias_esperados    = .N,
    dias_ausentes         = sum(tipo == "ausente"),
    dias_parciales        = sum(tipo == "parcial"),
    dias_completos        = sum(tipo == "completo"),
    pct_dias_ausentes     = round(100 * mean(tipo == "ausente"), 2),
    # % medio de horas que faltan SOLO en los dias parciales
    pct_horas_faltan_en_parciales = round(
      100 * mean(h_faltan[tipo == "parcial"] / 24), 2),
    # % global de horas faltantes en el mes (incluye ausentes y parciales)
    pct_horas_faltantes_mes = round(100 * sum(h_faltan) / (.N * 24), 2)
  ), by = .(ANIO, MES)]

  setorder(resumen, ANIO, MES)
  return(resumen[])
}

# ==============================================================================
# 1. NO2 (variable respuesta, NO imputada)
# ==============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("NO2 HORARIO - dias ausentes y horas faltantes por mes\n")
cat(strrep("=", 70), "\n", sep = "")

res_no2 <- list()
for (a in ANIOS) {
  f <- here("data", "processed", "Contaminacion", "horario",
            sprintf("aire_madrid_%d_No2_horarios.rds", a))
  if (!file.exists(f)) next
  d <- readRDS(f); setDT(d)
  d_largo <- d[, .(ESTACION, FECHA, valor = DATO)]
  res_no2[[as.character(a)]] <- resumir_horario(d_largo, a)
}
tabla_no2 <- rbindlist(res_no2)
fwrite(tabla_no2, file.path(DIR_TABLAS, "no2_horario_dias_por_mes.csv"))

# Resumen anual del NO2
cat("\nResumen ANUAL NO2 (suma de estacion-dias):\n")
print(tabla_no2[, .(
  est_dias_esperados = sum(est_dias_esperados),
  dias_ausentes      = sum(dias_ausentes),
  dias_parciales     = sum(dias_parciales),
  dias_completos     = sum(dias_completos),
  pct_horas_faltantes = round(100 * sum(pct_horas_faltantes_mes * est_dias_esperados) /
                              sum(est_dias_esperados) / 100 * 100, 2)
), by = ANIO])

cat("\nMeses NO2 con MAS dias ausentes o parciales (top 12):\n")
print(head(tabla_no2[order(-(dias_ausentes + dias_parciales)),
      .(ANIO, MES, dias_ausentes, dias_parciales, dias_completos,
        pct_horas_faltan_en_parciales, pct_horas_faltantes_mes)], 12))

# ==============================================================================
# 2. COVARIABLES CLIMATICAS (post-imputacion, solo estaciones que miden)
# ==============================================================================
VARS_CLIMA <- c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                "Presion Barométrica", "Radiación Solar", "Velocidad Viento")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("CLIMA HORARIO - dias ausentes y horas faltantes por mes (por variable)\n")
cat(strrep("=", 70), "\n", sep = "")

res_clima <- list()
for (a in ANIOS) {
  f <- here("data", "processed", "Clima", "horario",
            sprintf("meteo_madrid_%d_horario.rds", a))
  if (!file.exists(f)) next
  d <- readRDS(f); setDT(d)
  for (v in VARS_CLIMA) {
    if (!v %in% names(d)) next
    d_largo <- d[, .(ESTACION, FECHA, valor = get(v))]
    r <- resumir_horario(d_largo, a)     # estaciones que miden = detectadas dentro
    r[, Variable := v]
    res_clima[[paste(a, v)]] <- r
  }
}
tabla_clima <- rbindlist(res_clima, use.names = TRUE)
setcolorder(tabla_clima, c("Variable", "ANIO", "MES"))
fwrite(tabla_clima, file.path(DIR_TABLAS, "clima_horario_dias_por_mes.csv"))

# Resumen anual del clima por variable
cat("\nResumen ANUAL CLIMA por variable (suma estacion-dias):\n")
print(tabla_clima[, .(
  n_estaciones       = max(n_estaciones),
  est_dias_esperados = sum(est_dias_esperados),
  dias_ausentes      = sum(dias_ausentes),
  dias_parciales     = sum(dias_parciales),
  dias_completos     = sum(dias_completos),
  pct_dias_ausentes  = round(100 * sum(dias_ausentes) / sum(est_dias_esperados), 2)
), by = .(Variable, ANIO)][order(Variable, ANIO)])

cat("\nMeses CLIMA con MAS dias ausentes (top 15):\n")
print(head(tabla_clima[order(-dias_ausentes),
      .(Variable, ANIO, MES, n_estaciones, dias_ausentes, dias_parciales,
        pct_horas_faltan_en_parciales)], 15))

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Tablas guardadas:\n")
cat("  ", file.path(DIR_TABLAS, "no2_horario_dias_por_mes.csv"), "\n")
cat("  ", file.path(DIR_TABLAS, "clima_horario_dias_por_mes.csv"), "\n")
cat(strrep("=", 70), "\n", sep = "")
