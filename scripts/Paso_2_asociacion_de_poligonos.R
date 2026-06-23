# ==============================================================================
# STEP 2: CREATION OF THE SPATIO-TEMPORAL MASTER DATASET (DAILY OR HOURLY)
# ==============================================================================

library(data.table)
library(sf)
library(gstat)
library(here)
source(here("R", "cleaning_functions.R"))
source(here("R", "FUNCIONES_INTERPOLACION.R"))

# Deactivate the S2 engine to avoid issues with spatial operations
sf_use_s2(FALSE) 

# ==============================================================================
# BLOCK 0: CONFIGURATION (Only change the paths here!)
# ==============================================================================

ANIO <- 2019  

ruta_no2     <- here("data", "processed", "Contaminacion", "diario", paste0("aire_madrid_", ANIO, "_No2_trans_diarios.rds"))
ruta_trafico <- here("data", "processed", paste0("trafico_madrid_", ANIO, "_diario_barrio.rds"))
ruta_meteo   <- here("data", "processed", "Clima", "diario", paste0("meteo_madrid_", ANIO, "_diario.rds"))

# 2. OUTPUT paths
# Output paths are set after Block 1 detects the temporal scale

# ==============================================================================
# BLOCK 1: DATA LOADING AND TEMPORAL SCALE DETECTION
# ==============================================================================
dt_no2     <- readRDS(ruta_no2)
dt_trafico <- readRDS(ruta_trafico)
dt_meteo   <- readRDS(ruta_meteo)


setDT(dt_no2)
setDT(dt_trafico)
setDT(dt_meteo)

# HORA comes as factor "H01".."H24" — convert to integer in all tables
if ("HORA" %in% names(dt_no2) && is.factor(dt_no2$HORA)) {
  dt_no2[, HORA := as.integer(gsub("H", "", as.character(HORA)))]
}
if ("HORA" %in% names(dt_meteo) && is.factor(dt_meteo$HORA)) {
  dt_meteo[, HORA := as.integer(gsub("H", "", as.character(HORA)))]
}
# Automatic detection: Is it daily or hourly?
llaves_tiempo <- "FECHA"
if ("HORA" %in% names(dt_no2)) {
  llaves_tiempo <- c("FECHA", "HORA")
  cat("✅ HOURLY mode detected automatically.\n")
} else {
  cat("✅ DAILY mode detected automatically.\n")
}

# 2. OUTPUT paths (set after temporal detection)
escala_temporal <- if ("HORA" %in% llaves_tiempo) "horario" else "diario"
ruta_out_clima   <- here("data", "processed", "Clima", escala_temporal,
                         paste0("clima_interpolado_", escala_temporal, "_", ANIO, ".rds"))
ruta_out_maestro <- here("data", "processed",
                         paste0("dataset_maestro_inla_", ANIO, "_", toupper(escala_temporal), ".rds"))

# Normalization of neighborhood names in the traffic dataset
dt_trafico[, barrio := limpiar_nombres(barrio)]

# ==============================================================================
# BLOCK 2: GEOMETRIES (Districts and Neighborhoods)
# ==============================================================================
mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE) |>
  st_make_valid() |> st_transform(25830)
mapa_distritos$distrito <- limpiar_nombres(mapa_distritos$name)

mapa_barrios <- st_read(here("data", "raw", "Geometrias", "BARRIOS.shp"), quiet = TRUE) |>
  st_make_valid() |> st_transform(25830)
mapa_barrios$barrio <- limpiar_nombres(mapa_barrios$NOMBRE)

# ==============================================================================
# BLOCK 3: SPATIAL ASSOCIATION OF STATIONS
# ==============================================================================
estaciones_coords <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])
estaciones_sf <- st_as_sf(estaciones_coords, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)

est_distrito <- st_join(estaciones_sf, mapa_distritos[, "distrito"], join = st_intersects)
est_barrio   <- st_join(estaciones_sf, mapa_barrios[, "barrio"], join = st_intersects)

dt_est_geo <- merge(
  as.data.table(est_distrito)[, .(ESTACION, distrito)],
  as.data.table(est_barrio)[, .(ESTACION, barrio)],
  by = "ESTACION"
)

dt_no2 <- merge(dt_no2, dt_est_geo, by = "ESTACION", all.x = TRUE)

# ==============================================================================
# BLOCK 4: METEOROLOGY INTERPOLATION (IDW)
# ==============================================================================
coords_no2 <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])

dt_clima_interp <- interpolar_idw_clima(
  dt_meteo    = dt_meteo,
  dt_objetivo = coords_no2,
  variables   = c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                  "Presion Barométrica", "Radiación Solar","Velocidad Viento")
)

saveRDS(dt_clima_interp, ruta_out_clima)

# ==============================================================================
# BLOCK 5: THE MASTER JOIN (Traffic and Climate)
# ==============================================================================
# Join Traffic using dynamic time keys + neighborhood
llaves_trafico <- c(llaves_tiempo, "barrio")
cols_trafico   <- c(llaves_trafico, "intensidad", "ocupacion", "carga")# Ojo que en el 2020 no está ocupacion 

dt_maestro <- merge(dt_no2, dt_trafico[, ..cols_trafico], by = llaves_trafico, all.x = TRUE)

# Join interpolated Meteorology using same temporal keys as NO2
llaves_meteo <- c("ESTACION", llaves_tiempo)
dt_maestro <- merge(dt_maestro, dt_clima_interp, by = llaves_meteo, all.x = TRUE)

# Create a numeric ID for the districts
dt_maestro[, ID_DISTRITO := .GRP, by = distrito]

# ==============================================================================
# BLOCK 5b: REMOVE DAYS WITH NO COVARIATE DATA (e.g. city-wide blackout)
# A date is dropped if ALL climate variables are NA for ALL stations that day.
# ==============================================================================
cols_clima_check <- intersect(
  c("Temperatura", "Humedad_Relativa", "Precipitaciones",
    "Presion Barométrica", "Radiación Solar", "Velocidad Viento"),
  names(dt_maestro)
)

umbral_na_clima <- 0.10  # Drop day if any climate variable exceeds this NA rate

# For each date, compute the NA rate per variable; flag the date if any exceeds the threshold
fechas_exceden <- dt_maestro[,
  lapply(.SD, function(x) mean(is.na(x))),
  .SDcols = cols_clima_check,
  by = FECHA
][, excede := apply(.SD, 1, function(r) any(r > umbral_na_clima)),
  .SDcols = cols_clima_check
][excede == TRUE]

fechas_eliminar <- fechas_exceden$FECHA

if (length(fechas_eliminar) > 0) {
  cat(sprintf(
    "\n⚠️  Eliminando %d día(s) donde alguna variable climática supera el %.0f%% de NAs:\n",
    length(fechas_eliminar), umbral_na_clima * 100
  ))
  print(fechas_exceden[, c("FECHA", cols_clima_check), with = FALSE])
  dt_maestro <- dt_maestro[!FECHA %in% fechas_eliminar]
} else {
  cat(sprintf("\n✅ Ningún día supera el %.0f%% de NAs en covariables climáticas.\n",
              umbral_na_clima * 100))
}

# ==============================================================================
# BLOCK 6: COVARIATE STANDARDIZATION (Z-SCORE)
# ==============================================================================
# Traffic
dt_maestro[, intensidad_raw := intensidad]
dt_maestro[, carga_raw      := carga]
dt_maestro[, intensidad := scale(intensidad)[, 1]]
dt_maestro[, carga      := scale(carga)[, 1]]

# Climate (Wind included correctly)
cols_clima_std <- c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                    "Presion Barométrica", "Radiación Solar", "Velocidad Viento")
cols_clima_std <- intersect(cols_clima_std, names(dt_maestro))

for (v in cols_clima_std) {
  raw_name <- paste0(v, "_raw")
  dt_maestro[, (raw_name) := get(v)]
  dt_maestro[, (v) := scale(get(v))[, 1]]
}

cat("\n--- Covariate standardization ---\n")
cols_std <- c("intensidad", "carga", cols_clima_std)
for (v in cols_std) {
  cat(sprintf("%-25s mean = %7.4f | sd = %6.4f\n",
              v, mean(dt_maestro[[v]], na.rm = TRUE), sd(dt_maestro[[v]], na.rm = TRUE)))
}

# Strict chronological and spatial ordering
setorderv(dt_maestro, intersect(c("FECHA", "HORA", "ID_TIEMPO", "ESTACION"), names(dt_maestro)))

# ==============================================================================
# BLOCK 7: SAVING AND QUALITY CONTROL
# ==============================================================================
saveRDS(dt_maestro, ruta_out_maestro)
View(dt_maestro)
sf_use_s2(TRUE) 

cat("\n✅ Unification completed successfully!\n")
cat("Total rows in the Master Dataset:", nrow(dt_maestro), "\n")
if ("ID_TIEMPO" %in% names(dt_maestro))
  cat("Are there any NAs in ID_TIEMPO?:", any(is.na(dt_maestro$ID_TIEMPO)), "\n")

# Show a summary adapted to the existing columns
cols_print <- intersect(c("FECHA", "HORA", "ESTACION", "distrito", "ID_TIEMPO", "ID_DISTRITO"), names(dt_maestro))
print(head(dt_maestro[, ..cols_print]))

# ==============================================================================
# BLOCK 8: NA DIAGNOSTICS ON THE MASTER DATASET
# ==============================================================================

# 1. Count of NA/NaN per column (only columns with at least one NA)
cat("\n--- NA summary by column ---\n")
na_resumen <- dt_maestro[, lapply(.SD, function(x) sum(is.na(x)))] |>
  as.data.frame() |> t() |> as.data.frame() |>
  setNames("n_NA") |>
  tibble::rownames_to_column("columna") |>
  dplyr::filter(n_NA > 0) |>
  dplyr::mutate(pct = round(n_NA / nrow(dt_maestro) * 100, 2)) |>
  dplyr::arrange(dplyr::desc(n_NA))
print(na_resumen)

# 2. NO2: missing observations per station
col_no2 <- if ("DATO" %in% names(dt_maestro)) "DATO" else "DATO_DIARIO"
cat("\n--- NO2 NAs by station ---\n")
print(dt_maestro[is.na(get(col_no2)), .N, by = ESTACION][order(-N)])

# 3. Climate: which date(s) are fully missing?
cat("\n--- Climate NAs by date ---\n")
print(dt_maestro[is.na(Temperatura), .N, by = FECHA][order(-N)])

# 4. Traffic: which station-barrio-date combinations are missing?
cat("\n--- Traffic NAs (station / barrio / date) ---\n")
print(dt_maestro[is.na(intensidad), .(ESTACION, barrio, FECHA)])

# ==============================================================================
# BLOCK 9: NA DIAGNOSTICS IN THE MASTER DATASET — CONTAMINACIÓN AND TRÁFICO
# Same random day (dia) as the traffic section above.
# ==============================================================================

col_no2 <- if ("DATO" %in% names(dt_maestro)) "DATO" else "DATO_DIARIO"
maestro_dia <- dt_maestro[FECHA == dia]

# --- 9a. Contaminación (NO2): which stations are missing in dt_maestro? -------
no2_dia_sf <- as.data.frame(maestro_dia[, .(ESTACION, LONGITUD, LATITUD,
                                            no2 = get(col_no2))]) |>
  st_as_sf(coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830) |>
  mutate(tiene_dato = !is.na(no2))

cat("\n--- NO2 coverage on", format(dia), "---\n")
cat("Estaciones con dato NO2:", sum(no2_dia_sf$tiene_dato), "/", nrow(no2_dia_sf), "\n")
cat("Estaciones SIN dato NO2:", sum(!no2_dia_sf$tiene_dato), "\n")

p_no2 <- ggplot() +
  geom_sf(data = mapa_barrios_lower, fill = "grey95", color = "white", linewidth = 0.2) +
  geom_sf(data = no2_dia_sf, aes(color = tiene_dato), size = 3.5, alpha = 0.85) +
  scale_color_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#d73027"),
    labels = c("TRUE" = "Con dato", "FALSE" = "Sin dato (NA)"),
    name = "Estado"
  ) +
  labs(
    title = paste("Cobertura NO2 (maestro) —", format(dia, "%d %b %Y")),
    subtitle = paste0(
      sum(no2_dia_sf$tiene_dato), " estaciones con datos  |  ",
      sum(!no2_dia_sf$tiene_dato), " sin datos (rojo)"
    ),
    caption = paste("Fuente: dt_maestro | Año:", ANIO)
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(color = "grey40", margin = margin(b = 8)),
    legend.position = "right"
  )

# Guardado dinámico del mapa de NO2
ruta_plot_no2 <- here("outputs", "plots", paste0("cobertura_no2_maestro_", ANIO, "_", format(dia, "%Y%m%d"), ".png"))
ggsave(filename = ruta_plot_no2, plot = p_no2, width = 8, height = 6, dpi = 200, bg = "white")
cat("✅ Mapa de cobertura NO2 guardado en:", ruta_plot_no2, "\n")

# --- 9b. Tráfico: which barrios are missing in dt_maestro? --------------------
trafico_maestro_dia <- unique(maestro_dia[, .(barrio, intensidad)])

mapa_trafico_maestro <- mapa_barrios_lower |>
  left_join(as.data.frame(trafico_maestro_dia), by = "barrio")

cat("\n--- Traffic coverage in dt_maestro on", format(dia), "---\n")
cat("Barrios con dato de tráfico (maestro):",
    sum(!is.na(mapa_trafico_maestro$intensidad)), "/", nrow(mapa_trafico_maestro), "\n")
cat("Barrios SIN dato de tráfico (maestro):",
    sum(is.na(mapa_trafico_maestro$intensidad)), "\n")

p_trafico_maestro <- ggplot(mapa_trafico_maestro) +
  geom_sf(aes(fill = intensidad), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "grey80",
    name = "Intensidad\n(veh/h)"
  ) +
  labs(
    title = paste("Tráfico (maestro) —", format(dia, "%d %b %Y")),
    subtitle = paste0(
      sum(!is.na(mapa_trafico_maestro$intensidad)), " barrios con datos  |  ",
      sum(is.na(mapa_trafico_maestro$intensidad)), " sin datos (gris)"
    ),
    caption = paste("Fuente: dt_maestro | Año:", ANIO)
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(color = "grey40", margin = margin(b = 8)),
    legend.position = "right"
  )

# Guardado dinámico del mapa de tráfico maestro
ruta_plot_trafico <- here("outputs", "plots", paste0("cobertura_trafico_maestro_", ANIO, "_", format(dia, "%Y%m%d"), ".png"))
ggsave(filename = ruta_plot_trafico, plot = p_trafico_maestro, width = 8, height = 6, dpi = 200, bg = "white")
cat("✅ Mapa de tráfico maestro guardado en:", ruta_plot_trafico, "\n")
