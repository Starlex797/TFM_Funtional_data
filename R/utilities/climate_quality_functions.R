# Reusable functions for the 2019-2025 climate data-quality analysis.

canonical_covariate <- function(x) {
  key <- tolower(as.character(x))
  fcase(
    grepl("temperatura|temperature", key), "Temperature",
    grepl("humedad|humidity", key), "Relative humidity",
    grepl("precipit", key), "Precipitation",
    grepl("presion|barom", key), "Barometric pressure",
    grepl("solar", key) & !grepl("ultrav", key), "Solar radiation",
    grepl("viento|wind", key), "Wind speed",
    default = as.character(x)
  )
}

canonical_station <- function(x) {
  x <- as.character(x)
  fcase(
    grepl("^Centro Mpal\\. De Ac", x),
    paste0("Centro Mpal. De Ac", "\u00fa", "stica"),
    grepl("^J\\.M\\.D Chamart", x),
    paste0("J.M.D Chamart", "\u00ed", "n"),
    grepl("^J\\.M\\.D Chamber", x),
    paste0("J.M.D Chamber", "\u00ed"),
    grepl("^Pe.*agrande$", x), paste0("Pe", "\u00f1", "agrande"),
    grepl("^Plaza El.*ptica$", x), paste0("Plaza El", "\u00ed", "ptica"),
    grepl("^Plaza Espa.*a$", x), paste0("Plaza Espa", "\u00f1", "a"),
    default = x
  )
}

quality_inputs <- function(years, scales, suffix) {
  inputs <- CJ(Year = years, Scale = scales)
  inputs[, Path := mapply(
    function(year, scale) here(
      "data", "processed", "Clima", scale,
      sprintf("meteo_madrid_%d_%s%s.rds", year, scale, suffix)
    ),
    Year, Scale, USE.NAMES = FALSE
  )]
  missing <- inputs[!file.exists(Path)]
  if (nrow(missing)) {
    stop(
      "Missing ", nrow(missing), " expected files:\n  ",
      paste(basename(missing$Path), collapse = "\n  ")
    )
  }
  inputs[]
}

read_climate_states <- function(path, year, scale, valid_states) {
  data <- as.data.table(readRDS(path))
  state_columns <- grep("_estado$", names(data), value = TRUE)
  if (!length(state_columns) || !"ESTACION" %in% names(data)) {
    stop(basename(path), " has no station or state columns.")
  }

  if (scale == "mensual") {
    if (!"MES" %in% names(data)) stop(basename(path), " has no MES column.")
    data[, Analysis_date := as.IDate(paste0(MES, "-01"))]
  } else {
    if (!"FECHA" %in% names(data)) stop(basename(path), " has no FECHA column.")
    data[, Analysis_date := as.IDate(FECHA)]
  }

  long <- melt(
    data[, c("ESTACION", "Analysis_date", state_columns), with = FALSE],
    id.vars = c("ESTACION", "Analysis_date"),
    variable.name = "Covariate", value.name = "State"
  )
  long[, `:=`(
    Covariate = canonical_covariate(sub("_estado$", "", Covariate)),
    State = as.character(State),
    Station = canonical_station(ESTACION),
    Year = as.integer(year), Scale = scale
  )]
  long <- long[!grepl("ultrav", tolower(Covariate))]
  unknown <- setdiff(unique(long$State), valid_states)
  if (length(unknown) || anyNA(long$State)) {
    stop(basename(path), " contains unknown or missing states.")
  }
  long[, ESTACION := NULL][]
}

summarise_states <- function(data, groups) {
  result <- data[, .(
    n_cells = .N,
    n_ok = sum(State == "OK"),
    n_imputed = sum(State == "IMPUTADO"),
    n_failure = sum(State == "FALLO"),
    n_absent = sum(State == "AUSENTE"),
    n_no_sensor = sum(State == "SIN_SENSOR")
  ), by = groups]
  result[, n_evaluable := n_ok + n_imputed + n_failure + n_absent]
  result[, `:=`(
    pct_observed = 100 * n_ok / n_evaluable,
    pct_imputed = 100 * n_imputed / n_evaluable,
    pct_failure = 100 * n_failure / n_evaluable,
    pct_absent = 100 * n_absent / n_evaluable,
    pct_available = 100 * (n_ok + n_imputed) / n_evaluable,
    pct_na = 100 * (n_failure + n_absent) / n_evaluable
  )]
  result[]
}

daily_diagnostics <- function(data) {
  list(
    station_days = data[State != "SIN_SENSOR", .(
      Day_status = fcase(
        any(State == "AUSENTE"), "ABSENT",
        any(State == "FALLO"), "FAILURE",
        default = "AVAILABLE"
      ),
      Covariates_measured = uniqueN(Covariate)
    ), by = .(Station, Date = Analysis_date)],
    missing_cells = data[State %chin% c("FALLO", "AUSENTE"), .(
      Station, Covariate, Date = Analysis_date,
      Cause = fifelse(State == "FALLO", "FAILURE", "ABSENT")
    )]
  )
}

process_quality_files <- function(inputs, scale_labels, valid_states) {
  covariates <- stations <- station_days <- missing_cells <- vector("list", nrow(inputs))

  for (i in seq_len(nrow(inputs))) {
    input <- inputs[i]
    cat(sprintf(
      "[%02d/%02d] %d - %s\n", i, nrow(inputs), input$Year,
      scale_labels[[input$Scale]]
    ))
    long <- read_climate_states(input$Path, input$Year, input$Scale, valid_states)

    covariates[[i]] <- merge(
      summarise_states(long, c("Year", "Scale", "Covariate")),
      long[, .(
        Stations_total = uniqueN(Station),
        Stations_with_sensor = uniqueN(Station[State != "SIN_SENSOR"])
      ), by = .(Year, Scale, Covariate)],
      by = c("Year", "Scale", "Covariate")
    )
    stations[[i]] <- summarise_states(
      long, c("Year", "Scale", "Covariate", "Station")
    )
    if (input$Scale == "diario") {
      daily <- daily_diagnostics(long)
      station_days[[i]] <- daily$station_days
      missing_cells[[i]] <- daily$missing_cells
    }
    rm(long)
    invisible(gc(FALSE))
  }

  list(
    covariates = rbindlist(covariates, use.names = TRUE),
    stations = rbindlist(stations, use.names = TRUE),
    station_days = rbindlist(station_days, use.names = TRUE),
    missing_cells = rbindlist(missing_cells, use.names = TRUE)
  )
}

station_year_summary <- function(data) {
  result <- data[, .(
    n_ok = sum(n_ok), n_imputed = sum(n_imputed),
    n_failure = sum(n_failure), n_absent = sum(n_absent),
    n_evaluable = sum(n_evaluable)
  ), by = .(Year, Scale, Station)]
  result[, `:=`(
    pct_na = 100 * (n_failure + n_absent) / n_evaluable,
    pct_imputed = 100 * n_imputed / n_evaluable,
    pct_absent = 100 * n_absent / n_evaluable
  )]
  result[]
}

wide_quality_table <- function(data, id, years, scale_labels) {
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
    Scale, levels = names(scale_labels), labels = unname(scale_labels)
  )]
  columns <- unlist(lapply(
    years, function(year) paste0(year, c("_NA_pct", "_Imputed_pct", "_Absent_pct"))
  ))
  setcolorder(result, c(id, "Scale", columns))
  setorderv(result, c(id, "Scale"))
  numeric_columns <- setdiff(names(result), c(id, "Scale"))
  result[, (numeric_columns) := lapply(.SD, round, 2), .SDcols = numeric_columns]
  result[]
}

rank_station_days <- function(station_days) {
  result <- station_days[, .(
    Expected_days = .N,
    Failure_days = sum(Day_status == "FAILURE"),
    Absent_days = sum(Day_status == "ABSENT"),
    Affected_days = sum(Day_status != "AVAILABLE"),
    Affected_pct = 100 * mean(Day_status != "AVAILABLE"),
    Covariates_measured = max(Covariates_measured)
  ), by = Station]
  setorder(result, -Affected_days, -Affected_pct, Station)
  result[, Affected_pct := round(Affected_pct, 2)][]
}

plot_station_ranking <- function(ranking, output_file, top_n = 10L) {
  top <- head(ranking, top_n)
  plot_data <- melt(
    top,
    id.vars = c("Station", "Expected_days", "Affected_days", "Affected_pct"),
    measure.vars = c("Failure_days", "Absent_days"),
    variable.name = "Cause", value.name = "Days"
  )
  plot_data[, `:=`(
    Station = factor(Station, levels = rev(top$Station)),
    Cause = factor(
      Cause, levels = c("Failure_days", "Absent_days"),
      labels = c("Measurement failure", "Station absent from source")
    )
  )]

  plot <- ggplot(plot_data, aes(Station, Days, fill = Cause)) +
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
      title = "Stations with the most days without usable climate data",
      subtitle = "2019-2025; each calendar day is counted once per station",
      x = NULL, y = "Affected days",
      caption = paste(
        "Total NA days are divided into measurement failures and source-absent days.",
        "Imputed values and covariates without a sensor are excluded."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey30"),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      legend.position = "top", plot.margin = margin(5.5, 35, 5.5, 5.5)
    )

  ggsave(output_file, plot, width = 10, height = 6.5, dpi = 300, bg = "white")
  invisible(plot)
}

longest_missing_periods <- function(missing_cells, n = 15L) {
  setorder(missing_cells, Station, Covariate, Cause, Date)
  missing_cells[, New_period :=
    is.na(shift(Date)) | as.integer(Date - shift(Date)) != 1L,
  by = .(Station, Covariate, Cause)]
  missing_cells[, Period_id := cumsum(New_period),
    by = .(Station, Covariate, Cause)]

  periods <- missing_cells[, .(
    Start = min(Date), End = max(Date), Duration_days = .N
  ), by = .(Station, Covariate, Cause, Period_id)]
  periods <- periods[, .(
    Covariates_affected = paste(sort(unique(Covariate)), collapse = ", ")
  ), by = .(Station, Cause, Start, End, Duration_days)]
  periods[, Cause := fifelse(
    Cause == "FAILURE", "Data/sensor failure", "Station absent from source"
  )]
  setorder(periods, -Duration_days, Start, Station)
  result <- head(periods, n)
  result[, Rank := seq_len(.N)]
  setcolorder(result, c(
    "Rank", "Station", "Covariates_affected", "Cause",
    "Start", "End", "Duration_days"
  ))
  result[]
}

validate_quality_results <- function(raw, covariate_table, station_table) {
  bad_counts <- raw$covariates[
    n_ok + n_imputed + n_failure + n_absent + n_no_sensor != n_cells
  ]
  bad_sensors <- raw$stations[n_no_sensor > 0 & n_evaluable > 0]
  expected_covariates <- c(
    "Temperature", "Relative humidity", "Precipitation",
    "Barometric pressure", "Solar radiation", "Wind speed"
  )
  absent_exceeds_na <- function(table) {
    absent_columns <- grep("_Absent_pct$", names(table), value = TRUE)
    na_columns <- grep("_NA_pct$", names(table), value = TRUE)
    any(mapply(
      function(absent, na) table[[absent]] > table[[na]] + 1e-9,
      absent_columns, na_columns
    ))
  }
  errors <- c(
    nrow(bad_counts) > 0L, nrow(bad_sensors) > 0L,
    !setequal(covariate_table$Covariate, expected_covariates),
    nrow(covariate_table) != 18L,
    uniqueN(station_table$Station) != 26L,
    nrow(station_table) != 78L,
    absent_exceeds_na(covariate_table),
    absent_exceeds_na(station_table)
  )
  if (any(errors)) stop("Data-quality consistency checks failed.")
  invisible(list(bad_counts = bad_counts, bad_sensors = bad_sensors))
}
