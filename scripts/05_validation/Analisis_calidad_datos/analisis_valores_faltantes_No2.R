# ==============================================================================
# NO2 DATA QUALITY, 2019-2025
# ==============================================================================
# NA (%) = measurement failures + periods absent from the source.
# Absent (%) is a subset of NA (%). The current NO2 pipeline does not impute.

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})
source(here("R", "utilities", "academic_quality_tables.R"))

# ------------------------------------------------------------------------------
# Files to analyse. Change only this block to select another pollution version.
# Each entry must contain one file per year, in the same order as YEARS.
# ------------------------------------------------------------------------------
YEARS <- 2019:2025
INPUT_FILES <- list(
  Hourly = here(
    "data", "processed", "Contaminacion", "horario",
    sprintf("aire_madrid_%d_No2_horarios1.rds", YEARS)
  ),
  Daily = here(
    "data", "processed", "Contaminacion", "diario",
    sprintf("aire_madrid_%d_No2_trans_diarios1.rds", YEARS)
  ),
  Monthly = here(
    "data", "processed", "Contaminacion", "mensual",
    sprintf("aire_madrid_%d_log_No2_mensuales1.rds", YEARS)
  )
)

OUTPUT_DIR <- here(
  "outputs", "Analisis de la calidad de datos", "Contaminacion", "NO2"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

canonical_no2_station <- function(x) {
  x <- enc2utf8(as.character(x))
  fcase(
    grepl("^Av\\. Ram", x), paste0("Av. Ram", "\u00f3", "n y Cajal"),
    grepl("^M.*ndez .*lvaro$", x), paste0("M", "\u00e9", "ndez ", "\u00c1", "lvaro"),
    grepl("^Plaza El.*ptica$", x), paste0("Plaza El", "\u00ed", "ptica"),
    grepl("^Pza\\. de Espa", x), paste0("Pza. de Espa", "\u00f1", "a"),
    grepl("^P.* Castellana$", x), paste0("P", "\u00ba", " Castellana"),
    default = x
  )
}

input_manifest <- rbindlist(lapply(names(INPUT_FILES), function(scale) {
  files <- INPUT_FILES[[scale]]
  if (length(files) != length(YEARS)) {
    stop(scale, " must contain exactly one file per year.")
  }
  data.table(Year = YEARS, Scale = scale, File = files)
}))
input_manifest[, Detected_year := as.integer(sub(
  ".*(20[0-9]{2}).*", "\\1", basename(File)
))]
if (anyNA(input_manifest$Detected_year) ||
    any(input_manifest$Detected_year != input_manifest$Year)) {
  stop("A selected filename does not match its assigned year.")
}
missing_files <- input_manifest[!file.exists(File)]
if (nrow(missing_files)) {
  stop("Missing selected files:\n  ", paste(missing_files$File, collapse = "\n  "))
}
fwrite(input_manifest, file.path(OUTPUT_DIR, "input_files_used.csv"))

read_selected <- function(year, scale) {
  path <- input_manifest[Year == year & Scale == scale, File]
  as.data.table(readRDS(path))
}

summarise_states <- function(cells, groups) {
  result <- cells[, .(
    n_periods = .N,
    n_ok = sum(State == "OK"),
    n_failure = sum(State == "FALLO"),
    n_absent = sum(State == "AUSENTE")
  ), by = groups]
  result[, `:=`(
    pct_na = 100 * (n_failure + n_absent) / n_periods,
    pct_imputed = 0,
    pct_absent = 100 * n_absent / n_periods
  )]
  result[]
}

wide_quality_table <- function(data, id) {
  selected <- data[, c(
    id, "Scale", "Year", "pct_na", "pct_imputed", "pct_absent"
  ), with = FALSE]
  long <- melt(
    selected, id.vars = c(id, "Scale", "Year"),
    variable.name = "Metric", value.name = "Value"
  )
  labels <- c(pct_na = "NA_pct", pct_imputed = "Imputed_pct", pct_absent = "Absent_pct")
  long[, Column := paste(Year, labels[as.character(Metric)], sep = "_")]
  result <- dcast(
    long, as.formula(paste(id, "+ Scale ~ Column")), value.var = "Value"
  )
  result[, Scale := factor(
    Scale, levels = c("Hourly", "Daily", "Monthly")
  )]
  columns <- unlist(lapply(
    YEARS, function(year) paste0(year, c("_NA_pct", "_Imputed_pct", "_Absent_pct"))
  ))
  setcolorder(result, c(id, "Scale", columns))
  setorderv(result, c(id, "Scale"))
  result[, (columns) := lapply(.SD, round, 2), .SDcols = columns]
  result[]
}

covariate_results <- station_results <- daily_cells <- vector("list", length(YEARS))

for (i in seq_along(YEARS)) {
  year <- YEARS[i]
  cat(sprintf("[%d/%d] Processing NO2 %d\n", i, length(YEARS), year))
  hourly <- read_selected(year, "Hourly")
  daily <- read_selected(year, "Daily")
  monthly <- read_selected(year, "Monthly")

  required <- list(
    Hourly = c("ESTACION", "MAGNITUD", "FECHA", "DATO", "ESTADO"),
    Daily = c("ESTACION", "FECHA", "DATO_DIARIO"),
    Monthly = c("ESTACION", "MES", "DATO_MENSUAL")
  )
  selected <- list(Hourly = hourly, Daily = daily, Monthly = monthly)
  for (scale in names(required)) {
    missing <- setdiff(required[[scale]], names(selected[[scale]]))
    if (length(missing)) stop(scale, " file for ", year, " lacks: ", paste(missing, collapse = ", "))
  }
  if (!all(as.character(hourly$MAGNITUD) == "NO2")) {
    stop("The hourly file for ", year, " contains pollutants other than NO2.")
  }

  hourly[, `:=`(
    Station = canonical_no2_station(ESTACION),
    Date = as.IDate(FECHA),
    State = as.character(ESTADO)
  )]
  daily[, `:=`(
    Station = canonical_no2_station(ESTACION),
    Date = as.IDate(FECHA)
  )]
  monthly[, `:=`(
    Station = canonical_no2_station(ESTACION),
    Month = as.character(MES)
  )]

  unknown_states <- setdiff(unique(hourly$State), c("OK", "FALLO", "AUSENTE"))
  if (length(unknown_states) || anyNA(hourly$State)) {
    stop("Unknown hourly states in ", year, ".")
  }
  if (any(is.na(hourly$DATO) != hourly$State %chin% c("FALLO", "AUSENTE"))) {
    stop("DATO and ESTADO are inconsistent in the hourly file for ", year, ".")
  }

  station_sets <- lapply(selected, function(x) canonical_no2_station(x$ESTACION))
  if (!setequal(unique(station_sets$Hourly), unique(station_sets$Daily)) ||
      !setequal(unique(station_sets$Hourly), unique(station_sets$Monthly))) {
    stop("The three selected files use different NO2 stations in ", year, ".")
  }

  hourly_cells <- hourly[, .(
    Year = year, Scale = "Hourly", Station, Date, State
  )]
  hourly_day <- hourly[, .(
    n_hours = .N,
    n_absent = sum(State == "AUSENTE")
  ), by = .(Station, Date)]
  day_states <- merge(
    daily[, .(Station, Date, Is_NA = is.na(DATO_DIARIO))],
    hourly_day, by = c("Station", "Date"), all.x = TRUE
  )
  if (anyNA(day_states$n_hours)) stop("Daily and hourly dates do not match in ", year, ".")
  day_states[, State := fcase(
    !Is_NA, "OK",
    n_absent == n_hours, "AUSENTE",
    default = "FALLO"
  )]
  day_cells <- day_states[, .(
    Year = year, Scale = "Daily", Station, Date, State
  )]

  hourly[, Month := format(Date, "%Y-%m")]
  hourly_month <- hourly[, .(
    n_hours = .N,
    n_absent = sum(State == "AUSENTE")
  ), by = .(Station, Month)]
  month_states <- merge(
    monthly[, .(Station, Month, Is_NA = is.na(DATO_MENSUAL))],
    hourly_month, by = c("Station", "Month"), all.x = TRUE
  )
  if (anyNA(month_states$n_hours)) stop("Monthly and hourly periods do not match in ", year, ".")
  month_states[, State := fcase(
    !Is_NA, "OK",
    n_absent == n_hours, "AUSENTE",
    default = "FALLO"
  )]
  month_cells <- month_states[, .(
    Year = year, Scale = "Monthly", Station,
    Date = as.IDate(paste0(Month, "-01")), State
  )]

  cells <- rbindlist(list(hourly_cells, day_cells, month_cells))
  covariate_results[[i]] <- summarise_states(cells, c("Year", "Scale"))[
    , `:=`(Pollutant = "NO2", Stations = uniqueN(cells$Station))
  ]
  station_results[[i]] <- summarise_states(cells, c("Year", "Scale", "Station"))
  daily_cells[[i]] <- day_cells
  rm(hourly, daily, monthly, cells)
  invisible(gc(FALSE))
}

covariate_summary <- rbindlist(covariate_results)
station_summary <- rbindlist(station_results)
daily_data <- rbindlist(daily_cells)

covariate_table <- wide_quality_table(covariate_summary, "Pollutant")
station_table <- wide_quality_table(station_summary, "Station")

# Daily ranking and consecutive missing episodes.
station_ranking <- daily_data[, .(
  Expected_days = .N,
  Failure_days = sum(State == "FALLO"),
  Absent_days = sum(State == "AUSENTE"),
  Affected_days = sum(State != "OK"),
  Affected_pct = round(100 * mean(State != "OK"), 2)
), by = Station]
setorder(station_ranking, -Affected_days, -Affected_pct, Station)

missing_days <- copy(daily_data[State != "OK"])
missing_days[, Cause := fifelse(
  State == "AUSENTE", "Station absent from source", "Measurement failure"
)]
setorder(missing_days, Station, Cause, Date)
missing_days[, New_period :=
  is.na(shift(Date)) | as.integer(Date - shift(Date)) != 1L,
by = .(Station, Cause)]
missing_days[, Period_id := cumsum(New_period), by = .(Station, Cause)]
longest_periods <- missing_days[, .(
  Start = min(Date), End = max(Date), Duration_days = .N
), by = .(Station, Cause, Period_id)]
setorder(longest_periods, -Duration_days, Start, Station)
longest_periods <- head(longest_periods, 15L)
longest_periods[, Rank := seq_len(.N)]
longest_periods[, Period_id := NULL]
setcolorder(longest_periods, c(
  "Rank", "Station", "Cause", "Start", "End", "Duration_days"
))

# Consistency checks.
if (!setequal(covariate_summary$Year, YEARS) ||
    !setequal(covariate_summary$Scale, c("Hourly", "Daily", "Monthly")) ||
    any(covariate_summary$pct_absent > covariate_summary$pct_na + 1e-9) ||
    any(station_summary$pct_absent > station_summary$pct_na + 1e-9) ||
    any(covariate_summary$pct_imputed != 0) ||
    any(station_summary$pct_imputed != 0) ||
    nrow(covariate_table) != 3L ||
    nrow(station_table) != uniqueN(station_summary$Station) * 3L ||
    nrow(longest_periods) != 15L) {
  stop("NO2 data-quality consistency checks failed.")
}

fwrite(covariate_table, file.path(OUTPUT_DIR, "no2_by_scale_and_year.csv"))
fwrite(station_table, file.path(OUTPUT_DIR, "stations_by_scale_and_year.csv"))
fwrite(longest_periods, file.path(OUTPUT_DIR, "longest_missing_periods.csv"))

# Top ten stations with the most affected days.
top <- head(station_ranking, 10L)
top_long <- melt(
  top,
  id.vars = c("Station", "Expected_days", "Affected_days", "Affected_pct"),
  measure.vars = c("Failure_days", "Absent_days"),
  variable.name = "Cause", value.name = "Days"
)
top_long[, `:=`(
  Station = factor(Station, levels = rev(top$Station)),
  Cause = factor(
    Cause, levels = c("Failure_days", "Absent_days"),
    labels = c("Measurement failure", "Station absent from source")
  )
)]
ranking_plot <- ggplot(top_long, aes(Station, Days, fill = Cause)) +
  geom_col(width = 0.72) +
  geom_text(
    data = top,
    aes(
      factor(Station, levels = rev(top$Station)), Affected_days,
      label = sprintf("%d days (%.1f%%)", Affected_days, Affected_pct)
    ),
    inherit.aes = FALSE, hjust = -0.08, size = 3.1
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c("Measurement failure" = "#C43C39", "Station absent from source" = "#6C757D"),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = "NO2 stations with the most days without usable data",
    subtitle = "2019-2025; each calendar day is counted once per station",
    x = NULL, y = "Affected days",
    caption = paste(
      "Total NA days are divided into measurement failures and source-absent days.",
      "The current NO2 preprocessing does not impute missing values."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey30"),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "top", plot.margin = margin(5.5, 35, 5.5, 5.5)
  )
ggsave(
  file.path(OUTPUT_DIR, "stations_with_most_missing_days.png"),
  ranking_plot, width = 10, height = 6.5, dpi = 300, bg = "white"
)

# Academic table panels.
subtitle <- "Madrid NO2 monitoring network, 2019-2025"
panel_specs <- list(
  A = list(
    suffix = "NA_pct", label = "Total NA percentage",
    note = paste(
      "Total NA includes measurement failures and source-absent periods.",
      "Panel C reports the absent component included in this total."
    )
  ),
  B = list(
    suffix = "Imputed_pct", label = "Imputed percentage",
    note = "The selected NO2 preprocessing does not impute missing values; therefore all percentages are zero."
  ),
  C = list(
    suffix = "Absent_pct", label = "Absent component of total NA",
    note = paste(
      "Absent is a subset of total NA and denotes a complete period not published in the source.",
      "Measurement failure equals total NA minus Absent."
    )
  )
)
for (panel_name in names(panel_specs)) {
  spec <- panel_specs[[panel_name]]
  quality_metric_panel(
    covariate_table, "Pollutant", spec$suffix, YEARS,
    file.path(OUTPUT_DIR, sprintf("Table_1%s_NO2_%s.png", panel_name, spec$suffix)),
    sprintf("Table 1%s. %s by pollutant, scale and year", panel_name, spec$label),
    subtitle, spec$note, font_size = 9, row_height = 0.27, first_width = 1.4
  )
  quality_metric_panel(
    station_table, "Station", spec$suffix, YEARS,
    file.path(OUTPUT_DIR, sprintf("Table_2%s_stations_%s.png", panel_name, spec$suffix)),
    sprintf("Appendix Table 2%s. %s by station, scale and year", panel_name, spec$label),
    subtitle, spec$note, font_size = 7.5, row_height = 0.18, first_width = 2.2
  )
}

period_table <- copy(longest_periods)
period_table[, `:=`(
  Start = format(as.Date(Start), "%d/%m/%Y"),
  End = format(as.Date(End), "%d/%m/%Y")
)]
setnames(period_table, "Duration_days", "Duration (days)")
booktabs_png(
  period_table,
  file.path(OUTPUT_DIR, "Table_3_longest_missing_periods.png"),
  "Longest consecutive periods without NO2 data",
  "Fifteen longest daily episodes in the Madrid network, 2019-2025",
  "Consecutive missing days at the same station and with the same cause form one episode.",
  widths = c(0.5, 2.2, 2.4, 1.1, 1.1, 1.2),
  align = c("right", "left", "left", "center", "center", "right"),
  font_size = 8.5, row_height = 0.28
)

cat("\n=== NO2 validation ===\n")
cat("Selected files                 :", nrow(input_manifest), "\n")
cat("Years                          :", paste(YEARS, collapse = "-"), "\n")
cat("Scales                         : Hourly, Daily, Monthly\n")
cat("NO2 stations                   :", uniqueN(station_summary$Station), "\n")
cat("Rows in pollutant table        :", nrow(covariate_table), "\n")
cat("Rows in station table          :", nrow(station_table), "\n")
cat("Absent is a subset of NA       : yes\n")
cat("Imputation performed           : no\n")
cat("\nOutputs ->", OUTPUT_DIR, "\n")
