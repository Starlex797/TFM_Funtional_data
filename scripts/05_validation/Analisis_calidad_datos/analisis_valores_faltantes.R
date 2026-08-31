# ==============================================================================
# CLIMATE DATA QUALITY, 2019-2025
# ==============================================================================
# NA (%) = measurement failures + periods absent from the source.
# Absent (%) is a subset of NA (%), not an additional category.
# Stations without a sensor are excluded from every denominator.

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})
source(here("R", "utilities", "climate_quality_functions.R"))

# Configuration
YEARS <- 2019:2025
SCALES <- c("horario", "diario", "mensual")
SCALE_LABELS <- c(horario = "Hourly", diario = "Daily", mensual = "Monthly")
VALID_STATES <- c("OK", "IMPUTADO", "FALLO", "AUSENTE", "SIN_SENSOR")
SUFFIX <- "4"
OUTPUT_DIR <- here("outputs", "Analisis de la calidad de datos", "Clima")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Read and summarise the 21 year-scale files sequentially.
inputs <- quality_inputs(YEARS, SCALES, SUFFIX)
raw <- process_quality_files(inputs, SCALE_LABELS, VALID_STATES)

# Table 1: covariates x temporal scale x year.
covariate_table <- wide_quality_table(
  raw$covariates, "Covariate", YEARS, SCALE_LABELS
)

# Table 2: stations x temporal scale x year. Covariates without a sensor have
# already been excluded through n_evaluable.
station_summary <- station_year_summary(raw$stations)
station_table <- wide_quality_table(
  station_summary, "Station", YEARS, SCALE_LABELS
)

# Additional diagnostics at daily scale.
station_ranking <- rank_station_days(raw$station_days)
longest_periods <- longest_missing_periods(raw$missing_cells, 15L)

validate_quality_results(raw, covariate_table, station_table)

fwrite(
  covariate_table,
  file.path(OUTPUT_DIR, "covariates_by_scale_and_year.csv")
)
fwrite(
  station_table,
  file.path(OUTPUT_DIR, "stations_by_scale_and_year.csv")
)
fwrite(
  longest_periods,
  file.path(OUTPUT_DIR, "longest_missing_periods.csv")
)
plot_station_ranking(
  station_ranking,
  file.path(OUTPUT_DIR, "stations_with_most_missing_days.png")
)

cat("\n=== Validation ===\n")
cat("Files read                    :", nrow(inputs), "\n")
cat("Rows in covariate table       :", nrow(covariate_table), "\n")
cat("Rows in station table         :", nrow(station_table), "\n")
cat("Columns per main table        :", ncol(covariate_table), "\n")
cat("Rows in longest-period table  :", nrow(longest_periods), "\n")
cat("Absent is a subset of NA      : yes\n")
cat("\nOutputs ->", OUTPUT_DIR, "\n")
