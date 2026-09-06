# ==============================================================================
# CONTINUOUS DESCRIPTIVE MAPS OF ANNUAL NO2, 2019-2025
# ==============================================================================
# The script:
#   1. calculates the annual NO2 mean at each monitoring station;
#   2. estimates an optimal common number of IDW neighbours using pooled
#      leave-one-out cross-validation (LOOCV);
#   3. produces one continuous descriptive map per year and a comparative atlas;
#   4. exports the sensitivity analysis, LOOCV predictions and station summaries.
#
# IDW is used only for exploratory visualisation. It does not quantify prediction
# uncertainty and the resulting surfaces are not regulatory compliance maps.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(ggrepel)
  library(here)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

YEARS <- 2019:2025

# Editable IDW parameters -------------------------------------------------------
# Change IDW_POWER to rerun the sensitivity analysis and maps with another
# inverse-distance exponent. Set K_OVERRIDE to an integer to force that number
# of neighbours in the maps while still calculating the LOOCV optimum. Leave it
# as NA_integer_ to use the LOOCV optimum automatically.
IDW_POWER <- 2
K_OVERRIDE <- 5

GRID_RESOLUTION_M <- 250

if (length(IDW_POWER) != 1L || !is.finite(IDW_POWER) || IDW_POWER <= 0) {
  stop("IDW_POWER must be one positive finite number.")
}

INPUT_FILES <- here(
  "data", "processed", "Contaminacion", "diario",
  sprintf("aire_madrid_%d_No2_trans_diarios1.rds", YEARS)
)
STATION_TYPES_FILE <- here(
  "data", "processed", "Contaminacion", "mensual",
  "aire_madrid_2025_log_No2_mensuales1.rds"
)
DISTRICTS_FILE <- here("data", "raw", "geometrias", "DISTRITOS.shp")
OUTPUT_DIR <- here(
  "outputs", "figures", "EDA", "mapas", "Mapas_NO2_EDA"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

missing_files <- c(
  INPUT_FILES[!file.exists(INPUT_FILES)],
  STATION_TYPES_FILE[!file.exists(STATION_TYPES_FILE)],
  DISTRICTS_FILE[!file.exists(DISTRICTS_FILE)]
)
if (length(missing_files)) {
  stop("Missing input files:\n  ", paste(missing_files, collapse = "\n  "))
}

# ------------------------------------------------------------------------------
# 2. Annual station means
# ------------------------------------------------------------------------------

annual_station_list <- Map(function(path, year) {
  x <- as.data.table(readRDS(path))
  required <- c("ESTACION", "FECHA", "DATO_DIARIO", "LONGITUD", "LATITUD")
  missing_columns <- setdiff(required, names(x))
  if (length(missing_columns)) {
    stop(
      "The daily NO2 file for ", year, " lacks: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  x[, FECHA := as.IDate(FECHA)]
  if (anyNA(x$FECHA)) stop("Invalid dates in the daily NO2 file for ", year, ".")

  coordinates_per_station <- x[, .(
    n_longitude = uniqueN(LONGITUD[!is.na(LONGITUD)]),
    n_latitude = uniqueN(LATITUD[!is.na(LATITUD)])
  ), by = ESTACION]
  if (any(coordinates_per_station$n_longitude != 1L) ||
    any(coordinates_per_station$n_latitude != 1L)) {
    stop("Station coordinates are missing or not constant in ", year, ".")
  }

  x[, .(
    Year = as.integer(year),
    Valid_days = sum(!is.na(DATO_DIARIO)),
    Expected_days = .N,
    Coverage_pct = 100 * mean(!is.na(DATO_DIARIO)),
    Annual_NO2 = if (all(is.na(DATO_DIARIO))) {
      NA_real_
    } else {
      mean(DATO_DIARIO, na.rm = TRUE)
    },
    Longitude = first(LONGITUD),
    Latitude = first(LATITUD)
  ), by = .(Station = ESTACION)]
}, INPUT_FILES, YEARS)

annual_stations <- rbindlist(annual_station_list, use.names = TRUE)
if (anyNA(annual_stations$Annual_NO2)) {
  stop("At least one station-year has no valid daily NO2 observations.")
}
if (!all(annual_stations[, .N, by = Year]$N == 24L)) {
  stop("The analysis expects 24 evaluable NO2 stations in every year.")
}

# Station type is stable metadata and is taken from the 2025 monthly product,
# where NOM_TIPO is stored explicitly for every monitoring station.
station_types_raw <- as.data.table(readRDS(STATION_TYPES_FILE))
required_type_columns <- c("ESTACION", "NOM_TIPO")
missing_type_columns <- setdiff(required_type_columns, names(station_types_raw))
if (length(missing_type_columns)) {
  stop(
    "The station-type file lacks: ",
    paste(missing_type_columns, collapse = ", ")
  )
}

station_types <- unique(station_types_raw[, .(
  Station = as.character(ESTACION),
  Station_type_original = as.character(NOM_TIPO)
)])
if (anyDuplicated(station_types$Station)) {
  stop("At least one station has more than one monitoring-site type.")
}

station_types[, Station_type := fcase(
  grepl("^Urbana tr", Station_type_original), "Urban traffic",
  Station_type_original == "Urbana fondo", "Urban background",
  Station_type_original == "Suburbana", "Suburban",
  default = NA_character_
)]
if (anyNA(station_types$Station_type)) {
  stop("Unknown monitoring-site type in the station metadata.")
}

annual_stations <- merge(
  annual_stations,
  station_types[, .(Station, Station_type)],
  by = "Station", all.x = TRUE, sort = FALSE
)
if (anyNA(annual_stations$Station_type)) {
  stop("At least one annual NO2 station has no monitoring-site type.")
}
annual_stations[, Station_type := factor(
  Station_type,
  levels = c("Urban traffic", "Urban background", "Suburban")
)]
setorder(annual_stations, Year, Station)

districts <- st_read(DISTRICTS_FILE, quiet = TRUE)
if (is.na(st_crs(districts))) stop("The Madrid district layer has no CRS.")
districts <- st_make_valid(districts)
madrid_boundary <- st_union(districts)

annual_stations_sf <- st_as_sf(
  annual_stations,
  coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE
) |>
  st_transform(st_crs(districts))

if (any(lengths(st_intersects(annual_stations_sf, madrid_boundary)) == 0L)) {
  stop("At least one NO2 station lies outside the Madrid municipal boundary.")
}

# Three reference stations are labelled to make their location easy to identify.
# The prefix is used for Plaza Eliptica to avoid depending on accent encoding.
reference_station_mask <-
  annual_stations_sf$Station == "El Pardo" |
    grepl("^Plaza El", annual_stations_sf$Station) |
    annual_stations_sf$Station == "Escuelas Aguirre"
reference_stations_sf <- annual_stations_sf[reference_station_mask, ]

reference_counts <- table(reference_stations_sf$Year)
if (length(reference_counts) != length(YEARS) || any(reference_counts != 3L)) {
  stop(
    "El Pardo, Plaza Eliptica and Escuelas Aguirre must each occur once per year."
  )
}

reference_xy <- st_coordinates(reference_stations_sf)
reference_stations_sf$Label_X <- reference_xy[, 1]
reference_stations_sf$Label_Y <- reference_xy[, 2]

fwrite(
  annual_stations,
  file.path(OUTPUT_DIR, "annual_station_means_NO2_2019_2025.csv")
)

# ------------------------------------------------------------------------------
# 3. IDW functions and neighbour sensitivity analysis
# ------------------------------------------------------------------------------

idw_predict <- function(train_xy, train_values, target_xy, k, power = 2) {
  train_xy <- as.matrix(train_xy)
  target_xy <- as.matrix(target_xy)
  k <- min(as.integer(k), nrow(train_xy))
  if (k < 1L) stop("k must be at least 1.")

  predictions <- numeric(nrow(target_xy))
  for (i in seq_len(nrow(target_xy))) {
    distances <- sqrt(
      (train_xy[, 1] - target_xy[i, 1])^2 +
        (train_xy[, 2] - target_xy[i, 2])^2
    )
    zero_distance <- which(distances == 0)
    if (length(zero_distance)) {
      predictions[i] <- mean(train_values[zero_distance])
      next
    }
    neighbours <- order(distances)[seq_len(k)]
    weights <- distances[neighbours]^(-power)
    predictions[i] <- sum(weights * train_values[neighbours]) / sum(weights)
  }
  predictions
}

loocv_idw <- function(stations_sf, k, power = 2) {
  stations_sf <- stations_sf[order(stations_sf$Station), ]
  xy <- st_coordinates(stations_sf)
  observed <- stations_sf$Annual_NO2

  predicted <- vapply(seq_along(observed), function(i) {
    idw_predict(
      train_xy = xy[-i, , drop = FALSE],
      train_values = observed[-i],
      target_xy = xy[i, , drop = FALSE],
      k = k,
      power = power
    )
  }, numeric(1))

  data.table(
    Year = stations_sf$Year,
    Station = stations_sf$Station,
    K = as.integer(k),
    Observed = observed,
    Predicted = predicted,
    Error = predicted - observed
  )
}

n_by_year <- annual_stations[, .N, by = Year]
K_CANDIDATES <- seq_len(min(n_by_year$N) - 1L)

if (!is.na(K_OVERRIDE) &&
  (length(K_OVERRIDE) != 1L || !is.finite(K_OVERRIDE) ||
    K_OVERRIDE != as.integer(K_OVERRIDE) || !K_OVERRIDE %in% K_CANDIDATES)) {
  stop(
    "K_OVERRIDE must be NA_integer_ or one integer between ",
    min(K_CANDIDATES), " and ", max(K_CANDIDATES), "."
  )
}

loocv_predictions <- rbindlist(lapply(YEARS, function(year) {
  stations_year <- annual_stations_sf[annual_stations_sf$Year == year, ]
  rbindlist(lapply(K_CANDIDATES, function(k) {
    loocv_idw(stations_year, k = k, power = IDW_POWER)
  }))
}))

metric_summary <- function(observed, predicted) {
  error <- predicted - observed
  sst <- sum((observed - mean(observed))^2)
  list(
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    Bias = mean(error),
    R2_predictive = if (sst == 0) NA_real_ else 1 - sum(error^2) / sst
  )
}

sensitivity_overall <- loocv_predictions[, metric_summary(Observed, Predicted),
  by = K
]
sensitivity_by_year <- loocv_predictions[, metric_summary(Observed, Predicted),
  by = .(Year, K)
]
setorder(sensitivity_overall, K)
setorder(sensitivity_by_year, Year, K)

BEST_K_CV <- sensitivity_overall[which.min(RMSE), K]
K_USED <- if (is.na(K_OVERRIDE)) BEST_K_CV else as.integer(K_OVERRIDE)
K_SELECTION_MODE <- if (is.na(K_OVERRIDE)) {
  "Automatic LOOCV optimum"
} else {
  "Manual override"
}
sensitivity_overall[, `:=`(
  LOOCV_optimum = K == BEST_K_CV,
  Used_for_maps = K == K_USED
)]

selected_parameters <- data.table(
  Parameter = c(
    "IDW inverse-distance power",
    "LOOCV-optimal neighbours",
    "Neighbours used for maps",
    "Map neighbour selection",
    "Selection criterion",
    "Validation design",
    "Grid resolution (m)"
  ),
  Value = c(
    IDW_POWER,
    BEST_K_CV,
    K_USED,
    K_SELECTION_MODE,
    "Minimum pooled LOOCV RMSE",
    "Leave one station out; all station-years pooled",
    GRID_RESOLUTION_M
  )
)

fwrite(
  sensitivity_overall,
  file.path(OUTPUT_DIR, "idw_neighbour_sensitivity_overall.csv")
)
fwrite(
  sensitivity_by_year,
  file.path(OUTPUT_DIR, "idw_neighbour_sensitivity_by_year.csv")
)
fwrite(
  loocv_predictions,
  file.path(OUTPUT_DIR, "idw_loocv_predictions.csv")
)
fwrite(
  selected_parameters,
  file.path(OUTPUT_DIR, "idw_selected_parameters.csv")
)

sensitivity_plot <- ggplot() +
  geom_line(
    data = sensitivity_by_year,
    aes(K, RMSE, group = factor(Year)),
    colour = "grey72", linewidth = 0.55, alpha = 0.8
  ) +
  geom_line(
    data = sensitivity_overall,
    aes(K, RMSE), colour = "#1F4E79", linewidth = 1.15
  ) +
  geom_point(
    data = sensitivity_overall,
    aes(K, RMSE), colour = "#1F4E79", size = 2
  ) +
  geom_vline(
    xintercept = BEST_K_CV, colour = "#B2182B", linewidth = 0.7,
    linetype = "dashed"
  ) +
  geom_point(
    data = sensitivity_overall[K == BEST_K_CV],
    aes(K, RMSE), colour = "#B2182B", size = 3.3
  ) +
  annotate(
    "label",
    x = BEST_K_CV, y = Inf,
    label = sprintf("LOOCV optimum: k = %d", BEST_K_CV),
    vjust = 1.4,
    hjust = if (BEST_K_CV <= median(K_CANDIDATES)) -0.08 else 1.08,
    colour = "#B2182B", fill = "white", linewidth = 0.2, size = 3.5
  ) +
  scale_x_continuous(breaks = K_CANDIDATES) +
  labs(
    title = expression(
      "Sensitivity of annual NO"[2] * " IDW interpolation to neighbour count"
    ),
    subtitle = paste0(
      "Pooled station-level LOOCV, 2019-2025; inverse-distance power p = ",
      IDW_POWER
    ),
    x = "Number of nearest neighbours (k)",
    y = expression("LOOCV RMSE (" * mu * "g/m"^3 * ")"),
    caption = paste(
      "Grey lines show individual years; the blue line pools all station-years.",
      "The red line marks the pooled-RMSE optimum."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey35"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(colour = "grey40", hjust = 0)
  )

if (K_USED != BEST_K_CV) {
  sensitivity_plot <- sensitivity_plot +
    geom_vline(
      xintercept = K_USED, colour = "#E08214", linewidth = 0.7,
      linetype = "dotdash"
    ) +
    geom_point(
      data = sensitivity_overall[K == K_USED],
      aes(K, RMSE), colour = "#E08214", size = 3.3
    ) +
    annotate(
      "label",
      x = K_USED, y = Inf,
      label = sprintf("k = %d", K_USED),
      vjust = 3.5,
      hjust = if (K_USED <= median(K_CANDIDATES)) -0.08 else 1.08,
      colour = "#E08214", fill = "white", linewidth = 0.2, size = 3.5
    )
}

ggsave(
  file.path(OUTPUT_DIR, "idw_neighbour_sensitivity.png"),
  sensitivity_plot,
  width = 10, height = 6.5, dpi = 300, bg = "white"
)

# ------------------------------------------------------------------------------
# 4. Prediction grid and annual IDW surfaces
# ------------------------------------------------------------------------------

grid_points <- st_make_grid(
  madrid_boundary,
  cellsize = GRID_RESOLUTION_M, what = "centers"
)
grid_points <- st_sf(Grid_ID = seq_along(grid_points), geometry = grid_points)
grid_points <- grid_points[
  lengths(st_intersects(grid_points, madrid_boundary)) > 0L,
]
grid_xy <- st_coordinates(grid_points)

colour_limits <- range(annual_stations$Annual_NO2, finite = TRUE)
colour_limits <- c(
  5 * floor(colour_limits[1] / 5),
  5 * ceiling(colour_limits[2] / 5)
)
colour_breaks <- seq(colour_limits[1], colour_limits[2], by = 5)
colour_values <- c(
  "#2166AC", "#67A9CF", "#D1E5F0", "#FFFFBF",
  "#FDAE61", "#EF8A62", "#B2182B"
)
station_shapes <- c(
  "Urban traffic" = 24,
  "Urban background" = 21,
  "Suburban" = 22
)

surface_list <- vector("list", length(YEARS))
names(surface_list) <- as.character(YEARS)

map_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey30", size = 10),
    plot.caption = element_text(colour = "grey40", size = 8, hjust = 0),
    axis.title = element_blank(),
    axis.text = element_text(size = 7, colour = "grey45"),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

for (i in seq_along(YEARS)) {
  year <- YEARS[i]
  stations_year <- annual_stations_sf[annual_stations_sf$Year == year, ]
  reference_stations_year <-
    reference_stations_sf[reference_stations_sf$Year == year, ]
  station_xy <- st_coordinates(stations_year)
  prediction <- idw_predict(
    train_xy = station_xy,
    train_values = stations_year$Annual_NO2,
    target_xy = grid_xy,
    k = K_USED,
    power = IDW_POWER
  )

  surface_year <- data.table(
    X = grid_xy[, 1], Y = grid_xy[, 2],
    Annual_NO2_IDW = prediction, Year = factor(year, levels = YEARS)
  )
  surface_list[[i]] <- surface_year

  annual_map <- ggplot() +
    geom_raster(
      data = surface_year,
      aes(X, Y, fill = Annual_NO2_IDW)
    ) +
    geom_sf(
      data = districts, fill = NA, colour = "white", linewidth = 0.35
    ) +
    geom_sf(
      data = st_as_sf(stations_year),
      aes(shape = Station_type),
      size = 2.7, stroke = 0.65,
      fill = "white", colour = "#1A1A1A"
    ) +
    geom_text_repel(
      data = st_drop_geometry(reference_stations_year),
      aes(Label_X, Label_Y, label = Station),
      size = 3, fontface = "bold", colour = "#1A1A1A",
      bg.color = "white", bg.r = 0.12,
      box.padding = 0.45, point.padding = 0.4,
      min.segment.length = 0,
      segment.color = "grey35", segment.size = 0.3,
      seed = 7231, max.overlaps = Inf
    ) +
    scale_fill_gradientn(
      colours = colour_values,
      limits = colour_limits,
      breaks = colour_breaks,
      oob = scales::squish,
      name = expression(atop("Annual NO"[2] * " mean", mu * "g/m"^3))
    ) +
    scale_shape_manual(
      values = station_shapes, drop = FALSE, name = "Station type"
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 3.2, fill = "white"))
    ) +
    coord_sf(crs = st_crs(districts), expand = FALSE) +
    labs(
      title = bquote("Annual mean NO"[2] * " concentration - Madrid, " * .(year)),
      subtitle = sprintf(
        "Descriptive IDW surface (p = %g, k = %d; %s)",
        IDW_POWER, K_USED, tolower(K_SELECTION_MODE)
      ),
      caption = paste(
        "Symbols show monitoring-site type; labels identify three reference stations. Source: Madrid City Council Air Quality Monitoring Network.",
        "Descriptive interpolation; not an official regulatory or uncertainty map."
      )
    ) +
    map_theme

  ggsave(
    file.path(OUTPUT_DIR, sprintf("NO2_IDW_annual_mean_Madrid_%d.png", year)),
    annual_map,
    width = 9.5, height = 8.5, dpi = 300, bg = "white"
  )
  cat(sprintf("[%d/%d] Map saved for %d\n", i, length(YEARS), year))
}

# ------------------------------------------------------------------------------
# 5. Comparative atlas with a common scale
# ------------------------------------------------------------------------------

all_surfaces <- rbindlist(surface_list)
annual_stations_atlas <- annual_stations_sf
annual_stations_atlas$Year <- factor(
  annual_stations_atlas$Year,
  levels = YEARS
)
reference_stations_atlas <- reference_stations_sf
reference_stations_atlas$Year <- factor(
  reference_stations_atlas$Year,
  levels = YEARS
)

atlas <- ggplot() +
  geom_raster(
    data = all_surfaces,
    aes(X, Y, fill = Annual_NO2_IDW)
  ) +
  geom_sf(
    data = districts, fill = NA, colour = "white", linewidth = 0.18
  ) +
  geom_sf(
    data = annual_stations_atlas,
    aes(shape = Station_type),
    size = 1.25, stroke = 0.35,
    fill = "white", colour = "#1A1A1A"
  ) +
  geom_text_repel(
    data = st_drop_geometry(reference_stations_atlas),
    aes(Label_X, Label_Y, label = Station),
    size = 1.8, fontface = "bold", colour = "#1A1A1A",
    bg.color = "white", bg.r = 0.09,
    box.padding = 0.25, point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "grey35", segment.size = 0.18,
    seed = 7231, max.overlaps = Inf
  ) +
  facet_wrap(~Year, ncol = 4) +
  scale_fill_gradientn(
    colours = colour_values,
    limits = colour_limits,
    breaks = colour_breaks,
    oob = scales::squish,
    name = expression(atop("Annual NO"[2] * " mean", mu * "g/m"^3))
  ) +
  scale_shape_manual(
    values = station_shapes, drop = FALSE, name = "Station type"
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 2.8, fill = "white"))
  ) +
  coord_sf(crs = st_crs(districts), expand = FALSE) +
  labs(
    title = expression("Annual mean NO"[2] * " concentration in Madrid, 2019-2025"),
    subtitle = sprintf(
      "Descriptive IDW surfaces with a common colour scale (p = %g, k = %d; %s)",
      IDW_POWER, K_USED, tolower(K_SELECTION_MODE)
    ),
    caption = paste(
      "Labels identify El Pardo, Plaza Eliptica and Escuelas Aguirre.",
      "Descriptive interpolation; not an official regulatory or uncertainty map."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "grey30"),
    plot.caption = element_text(colour = "grey40", hjust = 0),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  file.path(OUTPUT_DIR, "NO2_IDW_annual_mean_Madrid_2019_2025_atlas.png"),
  atlas,
  width = 15, height = 8.5, dpi = 300, bg = "white"
)

cat("\n=== Annual NO2 IDW maps completed ===\n")
cat("Years                  :", paste(YEARS, collapse = "-"), "\n")
cat("Stations per year      : 24\n")
cat("Inverse-distance power :", IDW_POWER, "\n")
cat("Candidate neighbours   :", paste(range(K_CANDIDATES), collapse = "-"), "\n")
cat("LOOCV-optimal neighbours:", BEST_K_CV, "\n")
cat("Neighbours used in maps:", K_USED, "(", K_SELECTION_MODE, ")\n")
cat("Grid resolution        :", GRID_RESOLUTION_M, "m\n")
cat("Common colour limits   :", paste(colour_limits, collapse = "-"), "ug/m3\n")
cat("Outputs                :", OUTPUT_DIR, "\n")
