# ==============================================================================
# VALORES AUSENTES DEL MAESTRO DIARIO POR VARIABLE Y ANO
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

source(here("R", "utilities", "academic_quality_tables.R"))

ANIOS <- 2019:2025

# Variables que apareceran en la tabla. Editar esta lista manualmente.
# Usar VARIABLES_ANALIZAR <- NULL para incluir todas las variables del maestro.
VARIABLES_ANALIZAR <- c(
  "DATO_DIARIO",
  "LOG_NO2_DIARIO",
  "intensidad",
  "Temperatura",
  "Humedad_Relativa",
  "Precipitaciones",
  "Presion_Barometrica",
  "Radiacion_Solar",
  "Velocidad_Viento"
)

RUTA_DATOS <- here(
  "data", "processed", "Maestro", "diario",
  "dataset_maestro_inla_2019_2025_DIARIO.rds"
)
DIR_OUT <- here("outputs", "Analisis de la calidad de datos", "Maestro")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(RUTA_DATOS)) stop("No existe el maestro: ", RUTA_DATOS)

datos <- readRDS(RUTA_DATOS)
setDT(datos)

if (!"FECHA" %in% names(datos)) stop("El maestro no contiene la variable FECHA.")
if (anyNA(datos$FECHA)) {
  stop("FECHA contiene NA y esas filas no pueden asignarse a un ano.")
}

datos[, ANIO := as.integer(format(as.Date(FECHA), "%Y"))]
anios_ausentes <- setdiff(ANIOS, unique(datos$ANIO))
if (length(anios_ausentes) > 0L) {
  stop("No hay datos para los anos: ", paste(anios_ausentes, collapse = ", "))
}

variables_disponibles <- setdiff(names(datos), "ANIO")
variables <- if (is.null(VARIABLES_ANALIZAR)) {
  variables_disponibles
} else {
  variables_ausentes <- setdiff(VARIABLES_ANALIZAR, variables_disponibles)
  if (length(variables_ausentes) > 0L) {
    stop(
      "Las siguientes variables no existen en el maestro: ",
      paste(variables_ausentes, collapse = ", ")
    )
  }
  VARIABLES_ANALIZAR
}

# Formato largo: una fila por variable y ano.
tabla_na <- rbindlist(lapply(ANIOS, function(anio) {
  datos_anio <- datos[ANIO == anio]
  rbindlist(lapply(variables, function(variable) {
    observada <- !is.na(datos_anio[[variable]])
    data.table(
      Year = anio,
      Variable = variable,
      Observed = sum(observada),
      Missing = sum(!observada),
      Total = nrow(datos_anio),
      `Days with data` = uniqueN(datos_anio$FECHA[observada]),
      `Total days` = uniqueN(datos_anio$FECHA)
    )
  }))
}))
tabla_na[, `Missing (%)` := 100 * Missing / Total]
setorder(tabla_na, Variable, Year)

fwrite(
  tabla_na,
  file.path(DIR_OUT, "missing_values_by_variable_and_year_daily.csv")
)

# Formato ancho para la tabla academica: cada celda contiene n (%).
tabla_na[, Cell := sprintf("%s (%.2f%%)", format(Missing, big.mark = ","), `Missing (%)`)]
tabla_png <- dcast(tabla_na, Variable ~ Year, value.var = "Cell")

ruta_png <- file.path(
  DIR_OUT,
  "table_missing_values_by_variable_and_year_daily.png"
)

booktabs_png(
  tabla_png,
  ruta_png,
  title = "Missing values by variable and year",
  subtitle = "Daily master dataset, Madrid, 2019-2025",
  note = paste0(
    "Each cell reports the number of missing observations and, in parentheses, ",
    "the percentage of all rows available in that year. NA denotes an explicitly ",
    "missing value in the processed daily master dataset."
  ),
  widths = c(2.15, rep(1.05, length(ANIOS))),
  align = c("left", rep("right", length(ANIOS))),
  font_size = 7.5,
  row_height = 0.22,
  resolution = 220
)

# Segunda tabla: numero de observaciones no ausentes y cobertura temporal. Los
# dias con datos son fechas con al menos una observacion valida de la variable.
tabla_na[, Availability := sprintf(
  "%s | %d/%d",
  format(Observed, big.mark = ","), `Days with data`, `Total days`
)]
tabla_disponibilidad <- dcast(tabla_na, Variable ~ Year, value.var = "Availability")

ruta_disponibilidad <- file.path(
  DIR_OUT,
  "table_observations_and_days_by_variable_and_year_daily.png"
)

booktabs_png(
  tabla_disponibilidad,
  ruta_disponibilidad,
  title = "Available observations and days by variable and year",
  subtitle = "Daily master dataset, Madrid, 2019-2025",
  note = paste0(
    "Each cell reports valid station-day observations | days with at least one ",
    "valid observation / total days represented in the master dataset."
  ),
  widths = c(2.15, rep(1.25, length(ANIOS))),
  align = c("left", rep("right", length(ANIOS))),
  font_size = 7.5,
  row_height = 0.22,
  resolution = 220
)

# Diagnostico de 2022-2023: distingue dias sin datos en ninguna estacion de
# aquellos en los que los NA solo afectan a parte de la red.
variables_red <- c(
  "Precipitaciones",
  "Presion_Barometrica",
  "Radiacion_Solar",
  "Velocidad_Viento"
)
tabla_red <- rbindlist(lapply(2022:2023, function(anio) {
  rbindlist(lapply(variables_red, function(variable) {
    por_dia <- datos[ANIO == anio, .(
      rows = .N,
      missing = sum(is.na(get(variable)))
    ), by = FECHA]
    data.table(
      Year = anio,
      Covariate = variable,
      `Missing observations` = sum(por_dia$missing),
      `Days with NA` = sum(por_dia$missing > 0),
      `All-station days` = sum(por_dia$missing == por_dia$rows & por_dia$missing > 0),
      `Partial days` = sum(por_dia$missing > 0 & por_dia$missing < por_dia$rows)
    )
  }))
}))
setorder(tabla_red, Year, Covariate)

fwrite(
  tabla_red,
  file.path(DIR_OUT, "network_wide_missingness_2022_2023.csv")
)

tabla_red_png <- copy(tabla_red)
tabla_red_png[, Covariate := c(
  Precipitaciones = "Precipitation",
  Presion_Barometrica = "Barometric pressure",
  Radiacion_Solar = "Solar radiation",
  Velocidad_Viento = "Wind speed"
)[Covariate]]
inicios_anio <- tabla_red_png[, .I[1L], by = Year]$V1
tabla_red_png[, Year := as.character(Year)]
tabla_red_png[duplicated(Year), Year := ""]

ruta_red <- file.path(
  DIR_OUT,
  "table_network_wide_missingness_2022_2023.png"
)
booktabs_png(
  tabla_red_png,
  ruta_red,
  title = "Network-wide missingness in climate covariates",
  subtitle = "Daily master dataset, Madrid, 2022-2023",
  note = paste0(
    "An all-station day has NA in all 24 station rows. A partial day has NA in ",
    "only part of the network. All missing observations in these four covariates ",
    "occur on all-station days; no partial station failures were identified."
  ),
  widths = c(0.65, 1.70, 1.35, 1.05, 1.25, 0.95),
  align = c("right", "left", "right", "right", "right", "right"),
  group_starts = inicios_anio,
  font_size = 8.5,
  row_height = 0.26,
  resolution = 220
)

cat("\nFilas del maestro:", nrow(datos), "\n")
cat("Variables analizadas:", length(variables), "\n")
cat("CSV:", file.path(DIR_OUT, "missing_values_by_variable_and_year_daily.csv"), "\n")
cat("Tabla PNG:", ruta_png, "\n")
cat("Tabla de observaciones y dias:", ruta_disponibilidad, "\n")
cat("Tabla de ausencias generales:", ruta_red, "\n")
