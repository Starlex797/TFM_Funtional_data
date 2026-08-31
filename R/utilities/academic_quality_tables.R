# Minimal booktabs-style PNG tables rendered directly with R.

booktabs_png <- function(data, output_file, title, subtitle, note, widths, align,
                         group_starts = integer(), font_size = 9,
                         row_height = 0.25, resolution = 200) {
  data <- as.data.frame(data, check.names = FALSE)
  if (length(widths) != ncol(data) || length(align) != ncol(data)) {
    stop("Table widths and alignments must match the number of columns.")
  }
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

  ragg::agg_png(
    output_file, width = width, height = height, units = "in",
    res = resolution, background = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()

  draw_line <- function(y, colour = "#111111", size = 1) {
    grid::grid.lines(
      x = grid::unit(c(left, width - left), "in"),
      y = grid::unit(c(y, y), "in"),
      gp = grid::gpar(col = colour, lwd = size)
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
  grid::grid.text(
    title, x = grid::unit(left + 0.06, "in"), y = grid::unit(title_y, "in"),
    just = "left",
    gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = font_size + 3)
  )
  grid::grid.text(
    subtitle, x = grid::unit(left + 0.06, "in"), y = grid::unit(subtitle_y, "in"),
    just = "left",
    gp = grid::gpar(fontfamily = "Arial", fontsize = font_size + 1)
  )

  for (column in seq_along(data)) {
    header_align <- if (column == 1L) "left" else "center"
    grid::grid.text(
      names(data)[column],
      x = grid::unit(cell_x(column, header_align), "in"),
      y = grid::unit(header_y, "in"), just = header_align,
      gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = font_size)
    )
  }
  draw_line(data_top + 0.08, colour = "#B0B0B0")

  for (row in seq_len(nrow(data))) {
    y <- data_top - (row - 0.5) * row_height
    for (column in seq_along(data)) {
      grid::grid.text(
        as.character(data[row, column]),
        x = grid::unit(cell_x(column, align[column]), "in"),
        y = grid::unit(y, "in"), just = align[column],
        gp = grid::gpar(fontfamily = "Arial", fontsize = font_size)
      )
    }
  }
  for (row in group_starts[group_starts > 1L]) {
    draw_line(data_top - (row - 1L) * row_height, colour = "#B0B0B0", size = 0.7)
  }

  draw_line(body_bottom, size = 1.5)
  for (line in seq_along(note_lines)) {
    grid::grid.text(
      note_lines[line], x = grid::unit(left + 0.06, "in"),
      y = grid::unit(body_bottom - 0.16 * line, "in"), just = "left",
      gp = grid::gpar(fontfamily = "Arial", fontsize = max(font_size - 1, 7))
    )
  }
  draw_line(0.12, size = 2)
  invisible(output_file)
}

quality_metric_panel <- function(data, id, suffix, years, output_file,
                                 title, subtitle, note, font_size = 9,
                                 row_height = 0.25, first_width = 2.35) {
  metric_columns <- paste0(years, "_", suffix)
  panel <- data.table::copy(data[, c(id, "Scale", metric_columns), with = FALSE])
  panel[, Scale_order := match(Scale, c("Hourly", "Daily", "Monthly"))]
  data.table::setorderv(panel, c(id, "Scale_order"))
  panel[, Scale_order := NULL]
  group_starts <- panel[, .I[1L], by = id]$V1
  panel[[id]][duplicated(panel[[id]])] <- ""
  for (column in metric_columns) {
    data.table::set(panel, j = column, value = ifelse(
      is.na(panel[[column]]), "--", sprintf("%.2f", panel[[column]])
    ))
  }
  data.table::setnames(panel, metric_columns, as.character(years))

  booktabs_png(
    panel, output_file, title, subtitle, note,
    widths = c(first_width, 0.9, rep(0.72, length(years))),
    align = c("left", "center", rep("right", length(years))),
    group_starts = group_starts, font_size = font_size,
    row_height = row_height
  )
}
