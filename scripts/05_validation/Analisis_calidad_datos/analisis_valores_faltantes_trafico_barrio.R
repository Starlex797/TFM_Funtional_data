# ==============================================================================
# NEIGHBORHOOD TRAFFIC-INTENSITY DATA QUALITY, 2019-2025
# ==============================================================================
# Periods outside a neighborhood's active detector coverage are structural and
# excluded. Imputation occurred before aggregation and is not traceable here.

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})
source(here("R", "utilities", "academic_quality_tables.R"))

# Files to analyse. Change only this block to select another traffic version.
YEARS <- 2019:2025
TRAFFIC_ROOT <- here("data", "processed", "Trafico")
INPUT_FILES <- list(
  Hourly = file.path(
    TRAFFIC_ROOT, "Horario_Barrio", YEARS,
    sprintf("trafico_madrid_%d_horario_barrio1.rds", YEARS)
  ),
  Daily = file.path(
    TRAFFIC_ROOT, "Diario_Barrio", YEARS,
    sprintf("trafico_madrid_%d_diario_barrio1.rds", YEARS)
  ),
  Monthly = file.path(
    TRAFFIC_ROOT, "Mensual_Barrio", YEARS,
    sprintf("trafico_madrid_%d_mensual_barrio1.rds", YEARS)
  )
)

OUTPUT_DIR <- here(
  "outputs", "Analisis de la calidad de datos", "Trafico", "Barrio"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

input_manifest <- rbindlist(lapply(names(INPUT_FILES), function(scale) {
  files <- INPUT_FILES[[scale]]
  if (length(files) != length(YEARS)) stop(scale, " must contain one file per year.")
  data.table(Year = YEARS, Scale = scale, File = files)
}))
input_manifest[, Detected_year := as.integer(sub(
  ".*(20[0-9]{2}).*", "\\1", basename(File)
))]
if (anyNA(input_manifest$Detected_year) ||
    any(input_manifest$Detected_year != input_manifest$Year)) {
  stop("A selected traffic filename does not match its assigned year.")
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

normalise_neighborhood <- function(x) {
  trimws(tolower(enc2utf8(as.character(x))))
}

display_neighborhood <- function(x) {
  x <- as.character(x)
  paste0(toupper(substr(x, 1L, 1L)), substring(x, 2L))
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
    pct_imputed = NA_real_,
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
  result[, Scale := factor(Scale, levels = c("Hourly", "Daily", "Monthly"))]
  columns <- unlist(lapply(
    YEARS, function(year) paste0(year, c("_NA_pct", "_Imputed_pct", "_Absent_pct"))
  ))
  setcolorder(result, c(id, "Scale", columns))
  setorderv(result, c(id, "Scale"))
  result[, (columns) := lapply(.SD, round, 2), .SDcols = columns]
  result[]
}

variable_results <- neighborhood_results <- daily_cells <- vector("list", length(YEARS))

for (i in seq_along(YEARS)) {
  year <- YEARS[i]
  cat(sprintf("[%d/%d] Processing neighborhood traffic %d\n", i, length(YEARS), year))
  hourly <- read_selected(year, "Hourly")
  daily <- read_selected(year, "Daily")
  monthly <- read_selected(year, "Monthly")

  required <- list(
    Hourly = c("barrio", "FECHA", "HORA", "intensidad", "ESTADO"),
    Daily = c("barrio", "FECHA", "intensidad"),
    Monthly = c("barrio", "MES", "intensidad")
  )
  selected <- list(Hourly = hourly, Daily = daily, Monthly = monthly)
  for (scale in names(required)) {
    missing <- setdiff(required[[scale]], names(selected[[scale]]))
    if (length(missing)) stop(scale, " file for ", year, " lacks: ", paste(missing, collapse = ", "))
  }

  hourly[, `:=`(
    Neighborhood = normalise_neighborhood(barrio),
    Date = as.IDate(FECHA),
    State = as.character(ESTADO)
  )]
  daily[, `:=`(
    Neighborhood = normalise_neighborhood(barrio),
    Date = as.IDate(FECHA), Present = TRUE
  )]
  monthly[, `:=`(
    Neighborhood = normalise_neighborhood(barrio),
    Month = as.character(MES), Present = TRUE
  )]

  if (any(!hourly$State %chin% c("OK", "AUSENTE")) || anyNA(hourly$State) ||
      any(is.na(hourly$intensidad) != (hourly$State == "AUSENTE"))) {
    stop("Hourly intensity and ESTADO are inconsistent in ", year, ".")
  }
  if (hourly[, anyDuplicated(.SD), .SDcols = c("Neighborhood", "Date", "HORA")] ||
      daily[, anyDuplicated(.SD), .SDcols = c("Neighborhood", "Date")] ||
      monthly[, anyDuplicated(.SD), .SDcols = c("Neighborhood", "Month")]) {
    stop("Duplicated neighborhood-period keys in ", year, ".")
  }

  # Hourly rows already cover every hour inside active neighborhood-months.
  hourly_cells <- hourly[, .(
    Year = year, Scale = "Hourly", Neighborhood, Date, State
  )]

  # Expected daily/monthly periods come from hourly active coverage. Entire
  # months without detectors are therefore structural and never penalised.
  expected_days <- unique(hourly[, .(Neighborhood, Date)])
  day_states <- merge(
    expected_days,
    daily[, .(Neighborhood, Date, Present, Value = intensidad)],
    by = c("Neighborhood", "Date"), all.x = TRUE
  )
  day_states[, State := fcase(
    is.na(Present), "AUSENTE",
    is.na(Value) | !is.finite(Value), "FALLO",
    default = "OK"
  )]
  day_cells <- day_states[, .(
    Year = year, Scale = "Daily", Neighborhood, Date, State
  )]

  hourly[, Month := format(Date, "%Y-%m")]
  expected_months <- unique(hourly[, .(Neighborhood, Month)])
  month_states <- merge(
    expected_months,
    monthly[, .(Neighborhood, Month, Present, Value = intensidad)],
    by = c("Neighborhood", "Month"), all.x = TRUE
  )
  month_states[, State := fcase(
    is.na(Present), "AUSENTE",
    is.na(Value) | !is.finite(Value), "FALLO",
    default = "OK"
  )]
  month_cells <- month_states[, .(
    Year = year, Scale = "Monthly", Neighborhood,
    Date = as.IDate(paste0(Month, "-01")), State
  )]

  cells <- rbindlist(list(hourly_cells, day_cells, month_cells))
  variable_results[[i]] <- summarise_states(cells, c("Year", "Scale"))[
    , `:=`(
      Variable = "Traffic intensity",
      Neighborhoods = uniqueN(cells$Neighborhood)
    )
  ]
  neighborhood_results[[i]] <- summarise_states(
    cells, c("Year", "Scale", "Neighborhood")
  )
  daily_cells[[i]] <- day_cells
  rm(hourly, daily, monthly, cells)
  invisible(gc(FALSE))
}

variable_summary <- rbindlist(variable_results)
neighborhood_summary <- rbindlist(neighborhood_results)
daily_data <- rbindlist(daily_cells)
neighborhood_summary[, Neighborhood := display_neighborhood(Neighborhood)]
daily_data[, Neighborhood := display_neighborhood(Neighborhood)]

variable_table <- wide_quality_table(variable_summary, "Variable")
neighborhood_table <- wide_quality_table(neighborhood_summary, "Neighborhood")

# Daily neighborhood ranking and consecutive missing episodes.
neighborhood_ranking <- daily_data[, .(
  Expected_days = .N,
  Failure_days = sum(State == "FALLO"),
  Absent_days = sum(State == "AUSENTE"),
  Affected_days = sum(State != "OK"),
  Affected_pct = round(100 * mean(State != "OK"), 2)
), by = Neighborhood]
setorder(neighborhood_ranking, -Affected_days, -Affected_pct, Neighborhood)

missing_days <- copy(daily_data[State != "OK"])
missing_days[, Cause := fifelse(
  State == "AUSENTE", "Neighborhood absent from source", "Measurement failure"
)]
setorder(missing_days, Neighborhood, Cause, Date)
missing_days[, New_period :=
  is.na(shift(Date)) | as.integer(Date - shift(Date)) != 1L,
by = .(Neighborhood, Cause)]
missing_days[, Period_id := cumsum(New_period), by = .(Neighborhood, Cause)]
longest_periods <- missing_days[, .(
  Start = min(Date), End = max(Date), Duration_days = .N
), by = .(Neighborhood, Cause, Period_id)]
setorder(longest_periods, -Duration_days, Start, Neighborhood)
longest_periods <- head(longest_periods, 15L)
longest_periods[, Rank := seq_len(.N)]
longest_periods[, Period_id := NULL]
setcolorder(longest_periods, c(
  "Rank", "Neighborhood", "Cause", "Start", "End", "Duration_days"
))

if (!setequal(variable_summary$Year, YEARS) ||
    !setequal(variable_summary$Scale, c("Hourly", "Daily", "Monthly")) ||
    any(variable_summary$pct_absent > variable_summary$pct_na + 1e-9) ||
    any(neighborhood_summary$pct_absent > neighborhood_summary$pct_na + 1e-9) ||
    any(!is.na(variable_summary$pct_imputed)) ||
    any(!is.na(neighborhood_summary$pct_imputed)) ||
    nrow(variable_table) != 3L ||
    nrow(neighborhood_table) != uniqueN(neighborhood_summary$Neighborhood) * 3L) {
  stop("Traffic data-quality consistency checks failed.")
}

fwrite(variable_table, file.path(OUTPUT_DIR, "traffic_intensity_by_scale_and_year.csv"))
fwrite(neighborhood_table, file.path(OUTPUT_DIR, "neighborhoods_by_scale_and_year.csv"))
fwrite(longest_periods, file.path(OUTPUT_DIR, "longest_missing_periods.csv"))

# Top ten neighborhoods with the most missing days.
top <- head(neighborhood_ranking, 10L)
top_long <- melt(
  top,
  id.vars = c("Neighborhood", "Expected_days", "Affected_days", "Affected_pct"),
  measure.vars = c("Failure_days", "Absent_days"),
  variable.name = "Cause", value.name = "Days"
)
top_long[, `:=`(
  Neighborhood = factor(Neighborhood, levels = rev(top$Neighborhood)),
  Cause = factor(
    Cause, levels = c("Failure_days", "Absent_days"),
    labels = c("Measurement failure", "Neighborhood absent from source")
  )
)]
ranking_plot <- ggplot(top_long, aes(Neighborhood, Days, fill = Cause)) +
  geom_col(width = 0.72) +
  geom_text(
    data = top,
    aes(
      factor(Neighborhood, levels = rev(top$Neighborhood)), Affected_days,
      label = sprintf("%d days (%.2f%%)", Affected_days, Affected_pct)
    ),
    inherit.aes = FALSE, hjust = -0.08, size = 3.1
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c("Measurement failure" = "#C43C39", "Neighborhood absent from source" = "#6C757D"),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(
    title = "Neighborhoods with the most days without traffic-intensity data",
    subtitle = "2019-2025; only periods with active detector coverage are evaluated",
    x = NULL, y = "Affected days",
    caption = paste(
      "Months without neighborhood detector coverage are structural and excluded.",
      "The imputed share is not identifiable in the aggregated files."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey30"),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "top", plot.margin = margin(5.5, 45, 5.5, 5.5)
  )
ggsave(
  file.path(OUTPUT_DIR, "neighborhoods_with_most_missing_days.png"),
  ranking_plot, width = 10, height = 6.5, dpi = 300, bg = "white"
)

# Academic outputs. The large neighborhood table is paginated alphabetically.
subtitle <- "Madrid neighborhood traffic network, 2019-2025"
panel_specs <- list(
  A = list(
    suffix = "NA_pct", label = "Total NA percentage",
    note = paste(
      "Total NA includes source-absent periods and failures in active coverage.",
      "Periods without neighborhood detector coverage are excluded."
    )
  ),
  B = list(
    suffix = "Imputed_pct", label = "Imputed percentage",
    note = paste(
      "Not identifiable from these neighborhood aggregates: sensor-level imputation occurred",
      "before aggregation and its provenance was not retained. Dashes do not denote missing output."
    )
  ),
  C = list(
    suffix = "Absent_pct", label = "Absent component of total NA",
    note = paste(
      "Absent is a subset of total NA and refers only to missing periods during active coverage.",
      "It must not be added again to total NA."
    )
  )
)
for (panel_name in names(panel_specs)) {
  spec <- panel_specs[[panel_name]]
  quality_metric_panel(
    variable_table, "Variable", spec$suffix, YEARS,
    file.path(OUTPUT_DIR, sprintf("Table_1%s_traffic_%s.png", panel_name, spec$suffix)),
    sprintf("Table 1%s. %s by variable, scale and year", panel_name, spec$label),
    subtitle, spec$note, font_size = 9, row_height = 0.27, first_width = 2.1
  )
}

neighborhoods <- sort(unique(neighborhood_table$Neighborhood))
pages <- split(neighborhoods, ceiling(seq_along(neighborhoods) / 45L))
for (panel_name in names(panel_specs)) {
  spec <- panel_specs[[panel_name]]
  for (page in seq_along(pages)) {
    page_data <- neighborhood_table[Neighborhood %chin% pages[[page]]]
    quality_metric_panel(
      page_data, "Neighborhood", spec$suffix, YEARS,
      file.path(OUTPUT_DIR, sprintf(
        "Table_2%s_neighborhoods_%s_page_%02d.png", panel_name, spec$suffix, page
      )),
      sprintf(
        "Appendix Table 2%s. %s by neighborhood, scale and year (%d/%d)",
        panel_name, spec$label, page, length(pages)
      ),
      subtitle, spec$note, font_size = 7.2, row_height = 0.17, first_width = 2.4
    )
  }
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
  "Longest consecutive periods without neighborhood traffic intensity",
  "Fifteen longest daily episodes in Madrid, 2019-2025",
  "Consecutive missing days in the same neighborhood and with the same cause form one episode.",
  widths = c(0.5, 2.4, 2.6, 1.1, 1.1, 1.2),
  align = c("right", "left", "left", "center", "center", "right"),
  font_size = 8.2, row_height = 0.28
)

cat("\n=== Traffic validation ===\n")
cat("Selected files                    :", nrow(input_manifest), "\n")
cat("Years                             :", paste(YEARS, collapse = "-"), "\n")
cat("Scales                            : Hourly, Daily, Monthly\n")
cat("Neighborhoods across seven years :", uniqueN(neighborhood_summary$Neighborhood), "\n")
cat("Rows in variable table            :", nrow(variable_table), "\n")
cat("Rows in neighborhood table        :", nrow(neighborhood_table), "\n")
cat("Imputed percentage identifiable   : no\n")
cat("\nOutputs ->", OUTPUT_DIR, "\n")
