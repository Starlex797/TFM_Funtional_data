# Publication-ready climate data-quality tables without a browser dependency.

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(grid)
  library(ragg)
})

YEARS <- 2019:2025
OUTPUT_DIR <- here("outputs", "Analisis de la calidad de datos", "Clima")

booktabs_png <- function(data, filename, title, subtitle, note, widths, align,
                         group_starts = integer(), font_size = 9,
                         row_height = 0.25, resolution = 200) {
  data <- as.data.frame(data, check.names = FALSE)
  note_lines <- strwrap(note, width = max(90, floor(sum(widths) * 14)))
  body_bottom <- 0.34 + 0.17 * length(note_lines)
  data_top <- body_bottom + nrow(data) * row_height
  header_y <- data_top + 0.25
  subtitle_y <- header_y + 0.38
  title_y <- subtitle_y + 0.28
  top_y <- title_y + 0.22
  height <- top_y + 0.12
  width <- sum(widths) + 0.36
  left <- 0.18
  boundaries <- left + c(0, cumsum(widths))

  agg_png(
    file.path(OUTPUT_DIR, filename), width = width, height = height,
    units = "in", res = resolution, background = "white"
  )
  on.exit(dev.off(), add = TRUE)
  grid.newpage()

  draw_line <- function(y, colour = "#111111", size = 1) {
    grid.lines(
      x = unit(c(left, width - left), "in"), y = unit(c(y, y), "in"),
      gp = gpar(col = colour, lwd = size)
    )
  }
  cell_x <- function(column, justification) {
    switch(
      justification,
      left = boundaries[column] + 0.06,
      right = boundaries[column + 1L] - 0.06,
      center = mean(boundaries[c(column, column + 1L)])
    )
  }

  draw_line(top_y, size = 2)
  grid.text(
    title, x = unit(left + 0.06, "in"), y = unit(title_y, "in"),
    just = "left", gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = font_size + 3)
  )
  grid.text(
    subtitle, x = unit(left + 0.06, "in"), y = unit(subtitle_y, "in"),
    just = "left", gp = gpar(fontfamily = "Arial", fontsize = font_size + 1)
  )

  for (column in seq_along(data)) {
    header_align <- if (column == 1L) "left" else "center"
    grid.text(
      names(data)[column],
      x = unit(cell_x(column, header_align), "in"), y = unit(header_y, "in"),
      just = header_align,
      gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = font_size)
    )
  }
  draw_line(data_top + 0.08, colour = "#B0B0B0", size = 1)

  for (row in seq_len(nrow(data))) {
    y <- data_top - (row - 0.5) * row_height
    for (column in seq_along(data)) {
      value <- as.character(data[row, column])
      grid.text(
        value, x = unit(cell_x(column, align[column]), "in"), y = unit(y, "in"),
        just = align[column],
        gp = gpar(fontfamily = "Arial", fontsize = font_size)
      )
    }
  }

  for (row in group_starts[group_starts > 1L]) {
    draw_line(data_top - (row - 1L) * row_height, colour = "#B0B0B0", size = 0.7)
  }
  draw_line(body_bottom, size = 1.5)
  for (line in seq_along(note_lines)) {
    grid.text(
      note_lines[line], x = unit(left + 0.06, "in"),
      y = unit(body_bottom - 0.16 * line, "in"), just = "left",
      gp = gpar(fontfamily = "Arial", fontsize = max(font_size - 1, 7))
    )
  }
  draw_line(0.12, size = 2)
  invisible(file.path(OUTPUT_DIR, filename))
}

metric_panel <- function(data, id, suffix, filename, title, note,
                         font_size, row_height, first_width) {
  metric_columns <- paste0(YEARS, "_", suffix)
  panel <- copy(data[, c(id, "Scale", metric_columns), with = FALSE])
  panel[, Scale_order := match(Scale, c("Hourly", "Daily", "Monthly"))]
  setorderv(panel, c(id, "Scale_order"))
  panel[, Scale_order := NULL]
  group_starts <- panel[, .I[1L], by = id]$V1
  panel[[id]][duplicated(panel[[id]])] <- ""
  for (column in metric_columns) {
    set(panel, j = column, value = ifelse(
      is.na(panel[[column]]), "--", sprintf("%.2f", panel[[column]])
    ))
  }
  setnames(panel, metric_columns, as.character(YEARS))

  booktabs_png(
    panel, filename, title,
    "Madrid climate-monitoring network, 2019-2025", note,
    widths = c(first_width, 0.9, rep(0.72, length(YEARS))),
    align = c("left", "center", rep("right", length(YEARS))),
    group_starts = group_starts, font_size = font_size,
    row_height = row_height
  )
}

required <- c(
  "covariates_by_scale_and_year.csv", "stations_by_scale_and_year.csv",
  "longest_missing_periods.csv"
)
if (any(!file.exists(file.path(OUTPUT_DIR, required)))) {
  stop("Run analisis_valores_faltantes.R before rendering the tables.")
}

covariates <- fread(
  file.path(OUTPUT_DIR, "covariates_by_scale_and_year.csv"), encoding = "UTF-8"
)
stations <- fread(
  file.path(OUTPUT_DIR, "stations_by_scale_and_year.csv"), encoding = "UTF-8"
)

covariate_panels <- list(
  A = list(
    suffix = "NA_pct", file = "Table_1A_covariates_NA_pct.png",
    title = "Table 1A. Total NA percentage by covariate, temporal scale and year",
    note = paste(
      "Total NA includes measurement failures and source-absent periods.",
      "Table 1C reports the absent component included in this total; it must not be added again.",
      "Stations without the sensor are excluded."
    )
  ),
  B = list(
    suffix = "Imputed_pct", file = "Table_1B_covariates_Imputed_pct.png",
    title = "Table 1B. Imputed percentage by covariate, temporal scale and year",
    note = paste(
      "Imputed denotes a value recovered by preprocessing.",
      "Stations without the corresponding sensor are excluded from the denominator."
    )
  ),
  C = list(
    suffix = "Absent_pct", file = "Table_1C_covariates_Absent_pct.png",
    title = "Table 1C. Absent component of total NA by covariate, scale and year",
    note = paste(
      "Absent is a subset of total NA in Table 1A and denotes a period not published in the source.",
      "Measurement failure equals total NA minus Absent."
    )
  )
)

station_panels <- list(
  A = list(
    suffix = "NA_pct", file = "Table_2A_stations_NA_pct.png",
    title = "Appendix Table 2A. Total NA percentage by station, scale and year",
    note = paste(
      "Total NA includes measurement failures and source-absent periods.",
      "Table 2C reports the absent component included in this total.",
      "Only covariates measured by each station enter the denominator."
    )
  ),
  B = list(
    suffix = "Imputed_pct", file = "Table_2B_stations_Imputed_pct.png",
    title = "Appendix Table 2B. Imputed percentage by station, scale and year",
    note = paste(
      "Imputed values are reported separately from observed and missing values.",
      "Only covariates measured by each station enter the denominator."
    )
  ),
  C = list(
    suffix = "Absent_pct", file = "Table_2C_stations_Absent_pct.png",
    title = "Appendix Table 2C. Absent component of total NA by station, scale and year",
    note = paste(
      "Absent is a subset of total NA in Table 2A and must not be added again.",
      "At the three scales it refers respectively to absent hours, days and months."
    )
  )
)

for (panel in covariate_panels) {
  metric_panel(
    covariates, "Covariate", panel$suffix, panel$file, panel$title, panel$note,
    font_size = 9, row_height = 0.25, first_width = 2.35
  )
  cat("  ->", panel$file, "\n")
}
for (panel in station_panels) {
  metric_panel(
    stations, "Station", panel$suffix, panel$file, panel$title, panel$note,
    font_size = 7.5, row_height = 0.18, first_width = 2.2
  )
  cat("  ->", panel$file, "\n")
}

periods <- fread(
  file.path(OUTPUT_DIR, "longest_missing_periods.csv"), encoding = "UTF-8"
)
periods[, `:=`(
  Start = format(as.Date(Start), "%d/%m/%Y"),
  End = format(as.Date(End), "%d/%m/%Y")
)]
setnames(
  periods, c("Covariates_affected", "Duration_days"),
  c("Covariates affected", "Duration (days)")
)
booktabs_png(
  periods, "Table_3_longest_missing_periods.png",
  "Longest consecutive periods without climate data",
  "Fifteen longest daily episodes in the Madrid network, 2019-2025",
  paste(
    "Consecutive days form one episode. Covariates with the same station, dates and cause",
    "are combined to avoid repeating station-wide outages."
  ),
  widths = c(0.5, 1.7, 6.4, 2.25, 1.05, 1.05, 1.15),
  align = c("right", "left", "left", "left", "center", "center", "right"),
  font_size = 7.8, row_height = 0.28
)
cat("  -> Table_3_longest_missing_periods.png\n")
