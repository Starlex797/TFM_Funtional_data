# ==============================================================================
# CLIMATE DATA QUALITY FROM THE NATIVE DAILY SOURCE, 2019-2025
# ==============================================================================
# Same metric as analisis_valores_faltantes.R, but reading the daily files
# published by the provider (meteo_madrid_<year>_diario_fuente.rds) instead of
# the daily scale aggregated from hourly data.
#
# The point of the comparison: the aggregated pipeline discards a whole day when
# >= 30% of its hours are missing, which turns partial losses into total ones.
# The native daily source has no such threshold, so any drop in NA (%) measures
# exactly that amplification, while the absent component should stay put.

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})
source(here("R", "utilities", "climate_quality_functions.R"))

# Configuration
YEARS <- 2019:2025
SCALES <- "diario"
SCALE_LABELS <- c(diario = "Daily")
VALID_STATES <- c("OK", "IMPUTADO", "FALLO", "AUSENTE", "SIN_SENSOR")
SUFFIX <- "_fuente"
BASELINE_DIR <- here("outputs", "Analisis de la calidad de datos", "Clima")
OUTPUT_DIR <- here(
  "outputs", "Analisis de la calidad de datos", "Clima_diario_fuente"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Read and summarise the seven year files sequentially.
inputs <- quality_inputs(YEARS, SCALES, SUFFIX)
raw <- process_quality_files(inputs, SCALE_LABELS, VALID_STATES)

covariate_table <- wide_quality_table(
  raw$covariates, "Covariate", YEARS, SCALE_LABELS
)
station_summary <- station_year_summary(raw$stations)
station_table <- wide_quality_table(
  station_summary, "Station", YEARS, SCALE_LABELS
)
station_ranking <- rank_station_days(raw$station_days)
longest_periods <- longest_missing_periods(raw$missing_cells, 15L)

# ------------------------------------------------------------------------------
# Validation, aware of the single scale
# ------------------------------------------------------------------------------
# validate_quality_results() cannot be reused: it hard-codes 18 covariate rows
# (6 x 3 scales) and 78 station rows (26 x 3). Here there is one scale only, and
# the station roster of the daily source is not assumed to match the hourly one.

expected_covariates <- c(
  "Temperature", "Relative humidity", "Precipitation",
  "Barometric pressure", "Solar radiation", "Wind speed"
)
n_scales <- length(SCALES)
n_stations <- uniqueN(station_table$Station)

absent_exceeds_na <- function(table) {
  absent_columns <- grep("_Absent_pct$", names(table), value = TRUE)
  na_columns <- grep("_NA_pct$", names(table), value = TRUE)
  any(mapply(
    function(absent, na) table[[absent]] > table[[na]] + 1e-9,
    absent_columns, na_columns
  ))
}

checks <- c(
  bad_counts = nrow(raw$covariates[
    n_ok + n_imputed + n_failure + n_absent + n_no_sensor != n_cells
  ]) > 0L,
  bad_sensors = nrow(raw$stations[n_no_sensor > 0 & n_evaluable > 0]) > 0L,
  covariates = !setequal(covariate_table$Covariate, expected_covariates),
  covariate_rows = nrow(covariate_table) != length(expected_covariates) * n_scales,
  station_rows = nrow(station_table) != n_stations * n_scales,
  absent_covariates = absent_exceeds_na(covariate_table),
  absent_stations = absent_exceeds_na(station_table)
)
if (any(checks)) {
  stop(
    "Consistency checks failed: ",
    paste(names(checks)[checks], collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# Difference against the hourly-aggregated daily scale
# ------------------------------------------------------------------------------
# Positive Delta_NA means the aggregated pipeline reports MORE missing data than
# the native source, i.e. the 30% threshold discarded days the provider had
# published.

baseline_file <- file.path(BASELINE_DIR, "covariates_by_scale_and_year.csv")
comparison <- NULL

if (file.exists(baseline_file)) {
  metric_columns <- unlist(lapply(
    YEARS, function(year) paste0(year, c("_NA_pct", "_Imputed_pct", "_Absent_pct"))
  ))
  to_long <- function(table, source) {
    melt(
      table[Scale == "Daily"],
      id.vars = "Covariate", measure.vars = metric_columns,
      variable.name = "Column", value.name = source
    )
  }
  aggregated <- to_long(
    fread(baseline_file)[, Scale := as.character(Scale)], "Aggregated"
  )
  native <- to_long(copy(covariate_table)[, Scale := as.character(Scale)], "Native")

  comparison <- merge(aggregated, native, by = c("Covariate", "Column"))
  comparison[, `:=`(
    Year = as.integer(sub("_.*$", "", Column)),
    Metric = sub("^[0-9]{4}_", "", Column),
    Delta = round(Aggregated - Native, 2)
  )]
  comparison <- dcast(
    comparison, Covariate + Year ~ Metric,
    value.var = c("Aggregated", "Native", "Delta")
  )
  setorder(comparison, Covariate, Year)
  fwrite(comparison, file.path(OUTPUT_DIR, "daily_source_vs_aggregated.csv"))
} else {
  warning(
    "Baseline not found, the comparison table is skipped: ", baseline_file
  )
}

# ------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------

fwrite(covariate_table, file.path(OUTPUT_DIR, "covariates_by_year.csv"))
fwrite(station_table, file.path(OUTPUT_DIR, "stations_by_year.csv"))
fwrite(longest_periods, file.path(OUTPUT_DIR, "longest_missing_periods.csv"))
fwrite(inputs, file.path(OUTPUT_DIR, "input_files_used.csv"))
plot_station_ranking(
  station_ranking,
  file.path(OUTPUT_DIR, "stations_with_most_missing_days.png")
)

cat("\n=== Validation ===\n")
cat("Files read                    :", nrow(inputs), "\n")
cat("Stations in the daily source  :", n_stations, "\n")
cat("Rows in covariate table       :", nrow(covariate_table), "\n")
cat("Rows in station table         :", nrow(station_table), "\n")
cat("Rows in longest-period table  :", nrow(longest_periods), "\n")

if (!is.null(comparison)) {
  cat("\n=== Native daily source vs hourly aggregation ===\n")
  cat("Positive delta = the aggregated pipeline reports more NA.\n\n")
  print(comparison[, .(
    Covariate, Year,
    NA_aggregated = Aggregated_NA_pct, NA_native = Native_NA_pct,
    NA_delta = Delta_NA_pct,
    Absent_delta = Delta_Absent_pct,
    Imputed_native = Native_Imputed_pct
  )])
  cat("\nMean NA delta by year:\n")
  print(comparison[, .(NA_delta = round(mean(Delta_NA_pct), 2)), by = Year])
}

cat("\nOutputs ->", OUTPUT_DIR, "\n")
