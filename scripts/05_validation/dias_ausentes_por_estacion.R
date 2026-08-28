# ==============================================================================
# DIAS AUSENTES POR ESTACION, VARIABLE Y AÑO (2019-2025) - vista interpretativa
# ==============================================================================
# Un "dia ausente" = una estacion concreta que, un dia concreto, tiene sus 24
# horas vacias (offline / sin ninguna fila). Aqui NO se suman entre estaciones:
# se reporta cuantos dias falto CADA estacion (de los 365/366 del año).
#
# Genera:
#   (A) Tabla ANUAL resumen: variable x año (total dias-estacion ausentes,
#       media por estacion y nº de estaciones afectadas).
#   (B) Una tabla por AÑO: estacion x variable, con los dias ausentes de cada
#       estacion. Las estaciones que NO miden una variable se marcan "-".
#
# Se hace para las 6 covariables climaticas y para el NO2 (por separado, porque
# las redes de estaciones son distintas).
#
# Fuentes (post-preprocesamiento, horario):
#   Clima: data/processed/Clima/horario/meteo_madrid_<a>_horario.rds
#   NO2:   data/processed/Contaminacion/horario/aire_madrid_<a>_No2_horarios.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

ANIOS <- 2019:2025
VARS_CLIMA <- c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                "Presion Barométrica", "Radiación Solar", "Velocidad Viento")

DIR_OUT <- here("outputs", "tables", "calidad_datos", "dias_ausentes_por_estacion")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Helper: dias ausentes por estacion para un data.table largo (ESTACION,FECHA,valor)
#   Devuelve data.table con ESTACION y dias_ausentes (de n_dias del año).
#   Solo estaciones que MIDEN la variable (>=1 valor en el año).
# ------------------------------------------------------------------------------
dias_ausentes_estacion <- function(dt_largo, anio) {
  dt <- copy(dt_largo); setDT(dt)
  dt[, FECHA := as.Date(FECHA)]
  est_miden <- dt[!is.na(valor), unique(ESTACION)]
  if (length(est_miden) == 0) return(data.table(ESTACION = character(), dias_ausentes = integer()))

  dias <- seq(as.Date(sprintf("%d-01-01", anio)),
              as.Date(sprintf("%d-12-31", anio)), by = "day")
  por_ed <- dt[ESTACION %in% est_miden, .(h = sum(!is.na(valor))), by = .(ESTACION, FECHA)]
  grid <- CJ(ESTACION = est_miden, FECHA = dias)
  ed <- por_ed[grid, on = c("ESTACION", "FECHA")]
  ed[is.na(h), h := 0L]
  ed[, .(dias_ausentes = sum(h == 0L)), by = ESTACION]
}

# ==============================================================================
# Recolectar dias ausentes por estacion x variable x año (CLIMA)
# ==============================================================================
lista_clima <- list()
n_dias_anio <- integer()

for (a in ANIOS) {
  f <- here("data", "processed", "Clima", "horario",
            sprintf("meteo_madrid_%d_horario.rds", a))
  if (!file.exists(f)) next
  d <- readRDS(f); setDT(d)
  # El guardado horario no filtra por año: el fichero de 2025 arrastra ~90 dias
  # de 2026. Nos quedamos SOLO con las fechas del año objetivo.
  d <- d[year(as.Date(FECHA)) == a]
  # nº real de dias del año (365 o 366), no lo que arrastre el fichero
  n_dias_anio[as.character(a)] <- as.integer(
    as.Date(sprintf("%d-12-31", a)) - as.Date(sprintf("%d-01-01", a))) + 1L
  for (v in VARS_CLIMA) {
    if (!v %in% names(d)) next
    r <- dias_ausentes_estacion(d[, .(ESTACION, FECHA, valor = get(v))], a)
    if (nrow(r) == 0) next
    r[, `:=`(ANIO = a, Variable = v)]
    lista_clima[[paste(a, v)]] <- r
  }
}
clima_long <- rbindlist(lista_clima, use.names = TRUE)

# ==============================================================================
# Recolectar dias ausentes por estacion x año (NO2)
# ==============================================================================
lista_no2 <- list()
for (a in ANIOS) {
  f <- here("data", "processed", "Contaminacion", "horario",
            sprintf("aire_madrid_%d_No2_horarios.rds", a))
  if (!file.exists(f)) next
  d <- readRDS(f); setDT(d)
  d <- d[year(as.Date(FECHA)) == a]   # por si arrastra fechas de otro año
  r <- dias_ausentes_estacion(d[, .(ESTACION, FECHA, valor = DATO)], a)
  r[, `:=`(ANIO = a, Variable = "NO2")]
  lista_no2[[as.character(a)]] <- r
}
no2_long <- rbindlist(lista_no2, use.names = TRUE)

# ==============================================================================
# (A) TABLA ANUAL RESUMEN (variable x año)
# ==============================================================================
resumen_anual <- function(long) {
  long[, .(
    n_estaciones        = uniqueN(ESTACION),
    total_dias_ausentes = sum(dias_ausentes),
    media_por_estacion  = round(mean(dias_ausentes), 1),
    max_una_estacion    = max(dias_ausentes),
    est_afectadas       = sum(dias_ausentes > 0)
  ), by = .(Variable, ANIO)]
}

anual <- rbind(resumen_anual(clima_long), resumen_anual(no2_long))
setorder(anual, Variable, ANIO)

# Matriz compacta: total dias ausentes (Variable x Año)
matriz_total <- dcast(anual, Variable ~ ANIO, value.var = "total_dias_ausentes")
matriz_media <- dcast(anual, Variable ~ ANIO, value.var = "media_por_estacion")

cat("\n===== (A) TABLA ANUAL: total dias-estacion ausentes (Variable x Año) =====\n")
cat("(suma de todas las estaciones; nº de dias del año entre parentesis)\n")
cat("Dias por año:", paste(sprintf("%s=%d", names(n_dias_anio), n_dias_anio), collapse=" "), "\n\n")
print(matriz_total)
cat("\n===== (A') TABLA ANUAL: media de dias ausentes POR ESTACION (Variable x Año) =====\n\n")
print(matriz_media)

fwrite(anual,        file.path(DIR_OUT, "resumen_anual_dias_ausentes.csv"))
fwrite(matriz_total, file.path(DIR_OUT, "matriz_total_dias_ausentes.csv"))
fwrite(matriz_media, file.path(DIR_OUT, "matriz_media_por_estacion.csv"))

# ==============================================================================
# (B) UNA TABLA POR AÑO: estacion x variable (dias ausentes)
# ==============================================================================
# Combinamos clima + NO2 en un mismo formato largo y pivotamos por año.
todo_long <- rbind(clima_long, no2_long)
orden_vars <- c(VARS_CLIMA, "NO2")

for (a in ANIOS) {
  sub <- todo_long[ANIO == a]
  if (nrow(sub) == 0) next
  w <- dcast(sub, ESTACION ~ Variable, value.var = "dias_ausentes")
  # Ordenar columnas segun orden_vars (las que existan)
  cols <- c("ESTACION", intersect(orden_vars, names(w)))
  w <- w[, ..cols]
  # Ordenar filas por total de ausencias (mas problematicas arriba)
  vcols <- setdiff(names(w), "ESTACION")
  w[, total := rowSums(.SD, na.rm = TRUE), .SDcols = vcols]
  setorder(w, -total)

  # Guardar version numerica (NA = no mide la variable)
  fwrite(w, file.path(DIR_OUT, sprintf("dias_ausentes_por_estacion_%d.csv", a)))

  cat(sprintf("\n===== (B) AÑO %d — dias ausentes por estacion (de %d dias) =====\n",
              a, n_dias_anio[as.character(a)]))
  cat("   ('-' = la estacion no mide esa variable)\n\n")
  # Version legible: NA -> "-"
  w_print <- copy(w)
  for (cc in vcols) w_print[, (cc) := ifelse(is.na(get(cc)), "-", as.character(get(cc)))]
  print(w_print, nrows = 40)
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Guardado en:", DIR_OUT, "\n")
cat("  - resumen_anual_dias_ausentes.csv     (detalle variable x año)\n")
cat("  - matriz_total_dias_ausentes.csv      (matriz total)\n")
cat("  - matriz_media_por_estacion.csv       (matriz media/estacion)\n")
cat("  - dias_ausentes_por_estacion_<año>.csv (una por año)\n")
cat(strrep("=", 70), "\n", sep = "")
