# ==============================================================================
# STEP 3: CONCATENATE THE PER-YEAR DAILY MASTER DATASETS INTO ONE
# ==============================================================================
# Paso_2_asociacion_de_poligonos.R produce un maestro DIARIO por año en
# data/processed/Maestro/<ANIO>/dataset_maestro_inla_<ANIO>_DIARIO.rds. Este
# script los apila en un único dataset multi-año.
#
# Dos cosas NO se pueden apilar tal cual porque Paso_2 las calcula año a año:
#
#  - ID_DISTRITO: se asigna con .GRP por distrito DENTRO de cada fichero anual,
#    así que el mismo distrito puede tener un código distinto en cada año. Se
#    recalcula aquí sobre el conjunto combinado para que el mapeo sea único y
#    estable en todo el periodo.
#  - Las covariables estandarizadas (intensidad, carga y las seis variables
#    climáticas): Paso_2 las estandariza (z-score) con la media/sd de ESE año,
#    así que la misma covariable no es comparable de un año a otro tal como
#    llega. Se recalculan aquí sobre el periodo completo a partir de las
#    columnas *_raw, que Paso_2 conserva justamente para poder re-estandarizar
#    más adelante sin volver a la fuente.
#
# Todo lo demás (coordenadas, tipología de estación, distrito/barrio, NO2,
# tráfico y clima en crudo) viaja tal cual.
#
# FECHA_INICIO / FECHA_FIN (Bloque 0) permiten además recortar el resultado a
# un periodo concreto -p.ej. solo el primer semestre de 2020- en vez de
# quedarte con el/los año(s) completo(s). Con periodo definido, solo se leen
# los ficheros anuales que se solapan con él (no hace falta que ANIOS cubra
# nada más), y basta con un único año si el periodo cae dentro de uno solo.
# ==============================================================================

library(data.table)
library(here)

# ==============================================================================
# BLOCK 0: CONFIGURATION
# ==============================================================================

# NULL = usa todos los años que tengan maestro diario en
# data/processed/Maestro/<ANIO>/. Para limitar el rango, p.ej. c(2024, 2025).
ANIOS <- c(2020)

# Periodo a conservar en el maestro de salida. NULL en ambos = todo lo cargado
# en ANIOS. Acepta Date o texto "YYYY-MM-DD". Ejemplo: primer semestre de 2020:
#   FECHA_INICIO <- "2020-01-01"
#   FECHA_FIN    <- "2020-06-30"
FECHA_INICIO <- "2020-01-01"
FECHA_FIN <- "2020-06-30"

dir_maestro <- here("data", "processed", "Maestro")
ruta_maestro_anio <- function(a) {
  here(
    "data", "processed", "Maestro", as.character(a),
    paste0("dataset_maestro_inla_", a, "_DIARIO.rds")
  )
}

ruta_out <- function(fecha_min, fecha_max) {
  here(
    "data", "processed", "Maestro", "diario",
    sprintf(
      "dataset_maestro_inla_%s_%s_DIARIO.rds",
      format(fecha_min, "%Y%m%d"), format(fecha_max, "%Y%m%d")
    )
  )
}

# ==============================================================================
# BLOCK 1: DETECT AVAILABLE YEARS
# ==============================================================================
if (is.null(ANIOS)) {
  carpetas <- list.dirs(dir_maestro, recursive = FALSE, full.names = FALSE)
  anios_candidatos <- sort(suppressWarnings(as.integer(carpetas)))
  anios_candidatos <- anios_candidatos[!is.na(anios_candidatos)]
  ANIOS <- anios_candidatos[file.exists(vapply(anios_candidatos, ruta_maestro_anio, character(1)))]
}

# Si hay periodo, no hace falta leer años que quedan fuera de él: se restringe
# ANIOS a los que se solapan con [FECHA_INICIO, FECHA_FIN] antes de leer nada.
if (!is.null(FECHA_INICIO)) FECHA_INICIO <- as.Date(FECHA_INICIO)
if (!is.null(FECHA_FIN)) FECHA_FIN <- as.Date(FECHA_FIN)
if (!is.null(FECHA_INICIO) || !is.null(FECHA_FIN)) {
  anio_desde <- if (!is.null(FECHA_INICIO)) as.integer(format(FECHA_INICIO, "%Y")) else min(ANIOS)
  anio_hasta <- if (!is.null(FECHA_FIN)) as.integer(format(FECHA_FIN, "%Y")) else max(ANIOS)
  ANIOS <- intersect(ANIOS, anio_desde:anio_hasta)
}

if (length(ANIOS) < 1) {
  stop(
    "No hay ningún año con maestro diario disponible para el periodo/rango pedido."
  )
}
cat("Años a cargar:", paste(ANIOS, collapse = ", "), "\n")

# ==============================================================================
# BLOCK 2: READ EACH YEAR AND VALIDATE SCHEMA
# ==============================================================================
leer_maestro_anio <- function(a) {
  ruta <- ruta_maestro_anio(a)
  if (!file.exists(ruta)) stop("No existe el maestro diario de ", a, ": ", ruta)
  dt <- readRDS(ruta)
  setDT(dt)
  if ("HORA" %in% names(dt)) {
    stop(
      "El maestro de ", a, " tiene columna HORA (parece horario, no diario): ", ruta
    )
  }
  dt
}

lista_maestros <- lapply(ANIOS, leer_maestro_anio)
names(lista_maestros) <- as.character(ANIOS)

# Aviso (no interrumpe) si algún año trae columnas distintas a las del resto:
# rbindlist(fill = TRUE) las rellena con NA, pero es mejor saberlo antes de que
# aparezcan como NA "misteriosos" aguas abajo.
cols_por_anio <- lapply(lista_maestros, names)
cols_ref <- cols_por_anio[[1]]
for (a in names(cols_por_anio)[-1]) {
  faltan <- setdiff(cols_ref, cols_por_anio[[a]])
  sobran <- setdiff(cols_por_anio[[a]], cols_ref)
  if (length(faltan) || length(sobran)) {
    cat(sprintf(
      "⚠️  %s difiere en columnas de %s -> le faltan: [%s] | de más: [%s]\n",
      a, names(cols_por_anio)[1],
      paste(faltan, collapse = ", "), paste(sobran, collapse = ", ")
    ))
  }
}

# ==============================================================================
# BLOCK 3: CONCATENATE
# ==============================================================================
dt_maestro <- rbindlist(lista_maestros, fill = TRUE, use.names = TRUE)

stopifnot(!any(duplicated(dt_maestro[, .(ESTACION, FECHA)])))

# --- Recorte al periodo pedido -------------------------------------------------
# Se aplica antes de recalcular ID_DISTRITO y la estandarización para que
# ambos se calculen ya sobre el periodo final, no sobre los años completos
# leídos en el Bloque 1.
if (!is.null(FECHA_INICIO)) dt_maestro <- dt_maestro[FECHA >= FECHA_INICIO]
if (!is.null(FECHA_FIN)) dt_maestro <- dt_maestro[FECHA <= FECHA_FIN]
if (!nrow(dt_maestro)) {
  stop("El periodo FECHA_INICIO/FECHA_FIN no tiene ninguna fila en los años cargados (ANIOS).")
}

# --- ID_DISTRITO: recalculado sobre el conjunto combinado ---------------------
if ("distrito" %in% names(dt_maestro)) {
  dt_maestro[, ID_DISTRITO := .GRP, by = distrito]
}

# --- Re-estandarización de covariables sobre el periodo completo --------------
cols_raw <- grep("_raw$", names(dt_maestro), value = TRUE)
cols_std <- sub("_raw$", "", cols_raw)
cols_std_recalculadas <- character()
for (i in seq_along(cols_std)) {
  v <- cols_std[i]
  craw <- cols_raw[i]
  if (v %in% names(dt_maestro)) {
    dt_maestro[, (v) := scale(get(craw))[, 1]]
    cols_std_recalculadas <- c(cols_std_recalculadas, v)
  }
}
cat(
  "Covariables re-estandarizadas sobre el periodo final :",
  paste(cols_std_recalculadas, collapse = ", "), "\n"
)

# --- Orden cronológico y espacial ----------------------------------------------
setorderv(dt_maestro, intersect(c("FECHA", "ESTACION"), names(dt_maestro)))

# ==============================================================================
# BLOCK 4: SAVE AND SUMMARY
# ==============================================================================
fecha_min <- min(dt_maestro$FECHA)
fecha_max <- max(dt_maestro$FECHA)
ruta_out_maestro <- ruta_out(fecha_min, fecha_max)
dir.create(dirname(ruta_out_maestro), recursive = TRUE, showWarnings = FALSE)
saveRDS(dt_maestro, ruta_out_maestro)

cat("\n✅ Maestro guardado en:", ruta_out_maestro, "\n")
cat("Total de filas:", nrow(dt_maestro), "\n")
cat("Estaciones:", uniqueN(dt_maestro$ESTACION), "\n")
cat("Rango de fechas:", format(fecha_min), "a", format(fecha_max), "\n")
cat("\nFilas por año:\n")
print(dt_maestro[, .N, by = .(anio = as.integer(format(FECHA, "%Y")))][order(anio)])
