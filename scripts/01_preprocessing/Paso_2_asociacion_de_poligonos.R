# ==============================================================================
# STEP 2: CREATION OF THE SPATIO-TEMPORAL MASTER DATASET (DAILY OR HOURLY)
# ==============================================================================

library(data.table)
library(sf)
library(gstat)
library(here)
source(here("R", "cleaning", "cleaning_functions.R"))
source(here("R", "interpolation", "FUNCIONES_INTERPOLACION.R"))
source(here("R", "utilities", "dictionaries.R"))

# Deactivate the S2 engine to avoid issues with spatial operations
sf_use_s2(FALSE)

# ==============================================================================
# BLOCK 0: CONFIGURATION (Only change the paths here!)
# ==============================================================================

ANIO <- 2025

ruta_no2 <- here("data", "processed", "Contaminacion", "diario", paste0("aire_madrid_", ANIO, "_No2_trans_diarios1.rds"))
ruta_trafico <-here("data", "processed", "Trafico", "Diario_Barrio", ANIO, paste0("trafico_madrid_", ANIO, "_diario_barrio1.rds"))
ruta_meteo <- here("data", "processed", "Clima", "diario", paste0("meteo_madrid_", ANIO, "_diario5.rds"))


# 2. OUTPUT paths
# Output paths are set after Block 1 detects the temporal scale

# ==============================================================================
# BLOCK 1: DATA LOADING, VALIDATION AND TEMPORAL SCALE DETECTION
# ==============================================================================
dt_no2 <- readRDS(ruta_no2)
dt_trafico <- readRDS(ruta_trafico)
dt_meteo <- readRDS(ruta_meteo)

setDT(dt_no2)
setDT(dt_trafico)
setDT(dt_meteo)

# --- HORA: normalización uniforme de las tres fuentes -------------------------
# Convenio verificado en el preprocesamiento: HORA = k <=> intervalo [k-1, k).
# El aire y el clima llegan como H01..H24 (formato ancho del Ayuntamiento) y el
# tráfico como hour(fecha_hora) + 1, así que las tres escalas ya coinciden y no
# hay desfase que corregir. Se normalizan las tres por igual: si una llegara
# como texto en vez de factor, la clave quedaría de otro tipo y el join daría
# 100% de NA sin lanzar ningún error.
normalizar_hora <- function(dt, nombre) {
  if (!"HORA" %in% names(dt)) {
    return(invisible(NULL))
  }
  if (!is.integer(dt$HORA)) {
    dt[, HORA := as.integer(sub("^H", "", as.character(HORA)))]
  }
  if (!all(dt$HORA %in% 1:24)) stop("HORA fuera del rango 1..24 en ", nombre)
  invisible(NULL)
}
normalizar_hora(dt_no2, "NO2")
normalizar_hora(dt_meteo, "clima")
normalizar_hora(dt_trafico, "tráfico")

# --- Validación de entradas ---------------------------------------------------
# El preprocesamiento ha cambiado nombres de columnas y de estaciones más de una
# vez. Estas comprobaciones convierten ese tipo de cambio en un error inmediato
# en lugar de en columnas que desaparecen calladamente aguas abajo.
for (nm in c("dt_no2", "dt_trafico", "dt_meteo")) {
  d <- get(nm)
  if (!inherits(d$FECHA, "Date")) stop("FECHA no es de clase Date en ", nm)
  if (!all(year(d$FECHA) == ANIO)) stop("Hay fechas fuera del año ", ANIO, " en ", nm)
}

# Automatic detection: Is it daily or hourly?
llaves_tiempo <- "FECHA"
if ("HORA" %in% names(dt_no2)) {
  llaves_tiempo <- c("FECHA", "HORA")
  cat("✅ HOURLY mode detected automatically.\n")
} else {
  cat("✅ DAILY mode detected automatically.\n")
}

# Unicidad de las claves de unión. Se comprueba después de detectar la escala
# porque las claves son (FECHA) en diario y (FECHA, HORA) en horario: fijar HORA
# aquí rompería el modo diario, donde esa columna no existe.
stopifnot(
  !any(duplicated(dt_no2[, c("ESTACION", llaves_tiempo), with = FALSE])),
  !any(duplicated(dt_meteo[, c("ESTACION", llaves_tiempo), with = FALSE])),
  !any(duplicated(dt_trafico[, c("barrio", llaves_tiempo), with = FALSE]))
)

# 2. OUTPUT paths (set after temporal detection)
escala_temporal <- if ("HORA" %in% llaves_tiempo) "horario" else "diario"
ruta_out_clima <- here(
  "data", "processed", "Clima", escala_temporal,
  paste0("clima_interpolado_", escala_temporal, "_", ANIO, ".rds")
)
# El maestro se guarda en una carpeta por año. El nombre del fichero ya lleva la
# escala, así que el diario y el horario de un mismo año conviven sin pisarse.
ruta_out_maestro <- here(
  "data", "processed", "Maestro", as.character(ANIO),
  paste0("dataset_maestro_inla_", ANIO, "_", toupper(escala_temporal), ".rds")
)

# --- Nombres de barrios -------------------------------------------------------
# El tráfico ya llega en minúsculas pero conserva tildes ("águilas", "chamartín")
# y el shapefile viene en formato propio, así que hay que normalizar ambos lados
# con la misma función antes de cruzarlos (Bloque 2 hace lo propio con el mapa).
dt_trafico[, barrio := limpiar_nombres(barrio)]

# --- Columnas que ya no aportan nada al maestro -------------------------------
# MAGNITUD es constante ("NO2") desde que el preprocesamiento aplica el
# diccionario de magnitudes: peso muerto que viajaría hasta el dataset final.
if ("MAGNITUD" %in% names(dt_no2)) dt_no2[, MAGNITUD := NULL]

# Los indicadores de calidad del preprocesamiento (ESTADO, fila_ausente) no
# entran en el maestro: se modela sobre el dato tal cual y esas columnas se
# quedan en los ficheros de origen para el análisis de calidad.
cols_calidad <- c("ESTADO", "fila_ausente")
for (nm in c("dt_no2", "dt_trafico", "dt_meteo")) {
  d <- get(nm)
  sobran <- intersect(cols_calidad, names(d))
  if (length(sobran) > 0) d[, (sobran) := NULL]
}

# --- Tipología de la estación de contaminación --------------------------------
# NOM_TIPO (Suburbana / Urbana fondo / Urbana tráfico) clasifica las estaciones
# por entorno de medida y explica parte de la variabilidad que, si no, acabaría
# absorbida por el campo espacial.
# El fichero horario la trae; el diario no. Como es un atributo fijo de cada
# estación, en ese caso se recupera del fichero de estaciones en lugar de
# renunciar a la covariable.
if (!"NOM_TIPO" %in% names(dt_no2)) {
  ubic_estaciones <- fread(
    here("data", "raw", "Datos_contaminacion", "Estaciones", "datos.csv")
  )
  tipos_estacion <- unique(
    ubic_estaciones[, .(
      ESTACION = unname(nombres_estaciones_aire[as.character(CODIGO_CORTO)]),
      NOM_TIPO
    )][!is.na(ESTACION)],
    by = "ESTACION"
  )
  sin_tipo <- setdiff(unique(dt_no2$ESTACION), tipos_estacion$ESTACION)
  if (length(sin_tipo) > 0) {
    stop(
      "Estaciones sin NOM_TIPO en el fichero de estaciones: ",
      paste(sin_tipo, collapse = ", ")
    )
  }
  dt_no2 <- merge(dt_no2, tipos_estacion, by = "ESTACION", all.x = TRUE)
  cat("NOM_TIPO recuperado del fichero de estaciones.\n")
}
dt_no2[, NOM_TIPO := factor(NOM_TIPO)]
cat("\nTipología de estaciones (NOM_TIPO):\n")
print(dt_no2[, .(n_estaciones = uniqueN(ESTACION)), by = NOM_TIPO])

# --- Coordenadas UTM 30N en kilómetros ----------------------------------------
# Las mismas unidades que la malla SPDE (EPSG:25830 dividido por 1000). Se
# calculan aquí para que el maestro las lleve y el modelo no tenga que
# reproyectar por su cuenta, que es donde se descuadran datos y malla.
coords_km <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])
xy_km <- st_coordinates(st_transform(
  st_as_sf(coords_km, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830
)) / 1000
coords_km[, `:=`(X_km = xy_km[, 1], Y_km = xy_km[, 2])]
dt_no2 <- merge(dt_no2, coords_km[, .(ESTACION, X_km, Y_km)],
  by = "ESTACION", all.x = TRUE
)
stopifnot(!any(is.na(dt_no2$X_km)), !any(is.na(dt_no2$Y_km)))

# ==============================================================================
# BLOCK 2: GEOMETRIES (Districts and Neighborhoods)
# ==============================================================================
mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(25830)
mapa_distritos$distrito <- limpiar_nombres(mapa_distritos$name)

mapa_barrios <- st_read(here("data", "raw", "geometrias", "BARRIOS.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(25830)
mapa_barrios$barrio <- limpiar_nombres(mapa_barrios$NOMBRE)

# ==============================================================================
# BLOCK 3: SPATIAL ASSOCIATION OF STATIONS
# ==============================================================================
estaciones_coords <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])
estaciones_sf <- st_as_sf(estaciones_coords, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)

est_distrito <- st_join(estaciones_sf, mapa_distritos[, "distrito"], join = st_intersects)
est_barrio <- st_join(estaciones_sf, mapa_barrios[, "barrio"], join = st_intersects)

dt_est_geo <- merge(
  as.data.table(est_distrito)[, .(ESTACION, distrito)],
  as.data.table(est_barrio)[, .(ESTACION, barrio)],
  by = "ESTACION"
)

dt_no2 <- merge(dt_no2, dt_est_geo, by = "ESTACION", all.x = TRUE)

# ==============================================================================
# BLOCK 4: METEOROLOGY INTERPOLATION (winning method per variable)
# ==============================================================================
coords_no2 <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])

# ------------------------------------------------------------------------------
# CONFIGURACIÓN DE LA INTERPOLACIÓN — es lo único que hay que tocar aquí.
# ------------------------------------------------------------------------------
# Una fila por variable, y una tabla por escala temporal: el script elige la que
# toca según la escala detectada en el Bloque 1. Métodos y vecinos validados por
# LOOCV (ver R/interpolation/ y outputs/figures/interpolacion_clima/).
#
#   metodo : "kNN"  -> media simple de los k vecinos       (ignora p)
#            "IDW"  -> ponderación 1/d^p sobre k vecinos
#            "Vecino Cercano" -> el más próximo            (fuerza k = 1)
#            "Ensemble"       -> mezcla 1/RMSE de 1-NN, IDW p=1 y kNN
#   k      : número de vecinos usados.
#   p      : exponente de la distancia, solo aplica a "IDW".
#
# Para cambiar un método o una k, edita la celda correspondiente y ya está: el
# resto del script (validaciones, filtro de calidad y estandarización) se deriva
# de esta tabla.
vars_clima <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  "Presion_Barometrica", "Radiacion_Solar", "Velocidad_Viento"
)

config_clima_por_escala <- list(
  horario = data.table(
    variable = vars_clima,
    metodo   = c("IDW", "IDW", "kNN", "kNN", "kNN", "kNN"),
    k        = c(3L, 3L, 4L, 3L, 3L, 4L),
    p        = c(1, 1, NA, NA, NA, NA)
  ),
  diario = data.table(
    variable = vars_clima,
    metodo   = c("IDW", "IDW", "kNN", "kNN", "kNN", "kNN"),
    k        = c(3L, 3L, 3L, 3L, 3L, 2L),
    p        = c(1, 1, NA, NA, NA, NA)
  )
)

if (!escala_temporal %in% names(config_clima_por_escala)) {
  stop("No hay configuración de interpolación para la escala '", escala_temporal, "'.")
}
config_clima <- config_clima_por_escala[[escala_temporal]]
cat(sprintf("\nInterpolación climática en escala %s.\n", toupper(escala_temporal)))

dt_clima_interp <- interpolar_clima_por_metodo(
  dt_meteo = dt_meteo,
  dt_objetivo = coords_no2,
  config_variables = config_clima
)

saveRDS(dt_clima_interp, ruta_out_clima)

# ==============================================================================
# BLOCK 5: THE MASTER JOIN (Traffic and Climate)
# ==============================================================================
# Join Traffic using dynamic time keys + neighborhood
llaves_trafico <- c(llaves_tiempo, "barrio")
cols_trafico_candidatas <- c(llaves_trafico, "intensidad", "ocupacion", "carga")
cols_trafico <- intersect(cols_trafico_candidatas, names(dt_trafico))

dt_maestro <- merge(dt_no2, dt_trafico[, ..cols_trafico], by = llaves_trafico, all.x = TRUE)

# Join interpolated Meteorology using same temporal keys as NO2
llaves_meteo <- c("ESTACION", llaves_tiempo)
dt_maestro <- merge(dt_maestro, dt_clima_interp, by = llaves_meteo, all.x = TRUE)

# Create a numeric ID for the districts
dt_maestro[, ID_DISTRITO := .GRP, by = distrito]

# Comprobación de integridad del join: si alguna variable interpolada no llegó
# al maestro es un fallo, no algo que deba pasar desapercibido. No filtra nada.
cols_clima <- config_clima$variable
faltan_maestro <- setdiff(cols_clima, names(dt_maestro))
if (length(faltan_maestro) > 0) {
  stop("Faltan variables climáticas en el maestro: ", paste(faltan_maestro, collapse = ", "))
}

# No se descarta ninguna fecha por cobertura de covariables: el maestro conserva
# todos los días del año y los NA se dejan tal cual para tratarlos en el modelo.

# ==============================================================================
# BLOCK 6: COVARIATE STANDARDIZATION (Z-SCORE)
# ==============================================================================
# Traffic
dt_maestro[, intensidad_raw := intensidad]
dt_maestro[, carga_raw := carga]
dt_maestro[, intensidad := scale(intensidad)[, 1]]
dt_maestro[, carga := scale(carga)[, 1]]

# Climate: las seis variables interpoladas, ya validadas más arriba.
cols_clima_std <- cols_clima

for (v in cols_clima_std) {
  raw_name <- paste0(v, "_raw")
  dt_maestro[, (raw_name) := get(v)]
  dt_maestro[, (v) := scale(get(v))[, 1]]
}

cat("\n--- Covariate standardization ---\n")
cols_std <- c("intensidad", "carga", cols_clima_std)
for (v in cols_std) {
  cat(sprintf(
    "%-25s mean = %7.4f | sd = %6.4f\n",
    v, mean(dt_maestro[[v]], na.rm = TRUE), sd(dt_maestro[[v]], na.rm = TRUE)
  ))
}

# Strict chronological and spatial ordering
setorderv(dt_maestro, intersect(c("FECHA", "HORA", "ID_TIEMPO", "ESTACION"), names(dt_maestro)))

# ==============================================================================
# BLOCK 7: SAVING AND QUALITY CONTROL
# ==============================================================================
dir.create(dirname(ruta_out_maestro), recursive = TRUE, showWarnings = FALSE)
saveRDS(dt_maestro, ruta_out_maestro)
View(dt_maestro)
sf_use_s2(TRUE)

cat("\n✅ Unification completed successfully!\n")
cat("Total rows in the Master Dataset:", nrow(dt_maestro), "\n")
if ("ID_TIEMPO" %in% names(dt_maestro)) {
  cat("Are there any NAs in ID_TIEMPO?:", any(is.na(dt_maestro$ID_TIEMPO)), "\n")
}

# Show a summary adapted to the existing columns
cols_print <- intersect(c("FECHA", "HORA", "ESTACION", "distrito", "ID_TIEMPO", "ID_DISTRITO"), names(dt_maestro))
print(head(dt_maestro[, ..cols_print]))

# ==============================================================================
# BLOCK 8: NA DIAGNOSTICS ON THE MASTER DATASET
# ==============================================================================

# 1. NA por columna, TODAS las variables (también las completas: saber que una
#    columna está a cero es tan informativo como saber que le faltan datos).
#
#    pct                 = % sobre el total de filas del maestro, donde una fila
#                          es estación x fecha (x hora en escala horaria).
#    instantes_afectados = fechas (o fecha+hora) con al menos un NA. Distingue un
#                          hueco concentrado en pocos días de otro repartido por
#                          todo el año, que en pct son indistinguibles.
cat("\n--- NA summary by column ---\n")
n_filas <- nrow(dt_maestro)
n_instantes <- uniqueN(dt_maestro[, llaves_tiempo, with = FALSE])

na_resumen <- data.table(
  columna = names(dt_maestro),
  n_NA = vapply(dt_maestro, function(x) sum(is.na(x)), integer(1))
)
na_resumen[, pct := round(100 * n_NA / n_filas, 2)]
na_resumen[, instantes_afectados := vapply(columna, function(v) {
  filas_na <- dt_maestro[is.na(get(v))]
  if (nrow(filas_na) == 0L) 0L else uniqueN(filas_na[, llaves_tiempo, with = FALSE])
}, integer(1))]
setorder(na_resumen, -n_NA, columna)

cat(sprintf(
  "Filas: %d (%d estaciones x %d instantes)\n",
  n_filas, uniqueN(dt_maestro$ESTACION), n_instantes
))
print(na_resumen, nrows = Inf)

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
  no2 = get(col_no2)
)]) |>
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
    plot.title = element_text(face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(color = "grey40", margin = margin(b = 8)),
    legend.position = "right"
  )

# Guardado dinámico del mapa de NO2
ruta_plot_no2 <- here("outputs", "analysis", "plots", paste0("cobertura_no2_maestro_", ANIO, "_", format(dia, "%Y%m%d"), ".png"))
ggsave(filename = ruta_plot_no2, plot = p_no2, width = 8, height = 6, dpi = 200, bg = "white")
cat("✅ Mapa de cobertura NO2 guardado en:", ruta_plot_no2, "\n")

# --- 9b. Tráfico: which barrios are missing in dt_maestro? --------------------
trafico_maestro_dia <- unique(maestro_dia[, .(barrio, intensidad)])

mapa_trafico_maestro <- mapa_barrios_lower |>
  left_join(as.data.frame(trafico_maestro_dia), by = "barrio")

cat("\n--- Traffic coverage in dt_maestro on", format(dia), "---\n")
cat(
  "Barrios con dato de tráfico (maestro):",
  sum(!is.na(mapa_trafico_maestro$intensidad)), "/", nrow(mapa_trafico_maestro), "\n"
)
cat(
  "Barrios SIN dato de tráfico (maestro):",
  sum(is.na(mapa_trafico_maestro$intensidad)), "\n"
)

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
    plot.title = element_text(face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(color = "grey40", margin = margin(b = 8)),
    legend.position = "right"
  )

# Guardado dinámico del mapa de tráfico maestro
ruta_plot_trafico <- here("outputs", "analysis", "plots", paste0("cobertura_trafico_maestro_", ANIO, "_", format(dia, "%Y%m%d"), ".png"))
ggsave(filename = ruta_plot_trafico, plot = p_trafico_maestro, width = 8, height = 6, dpi = 200, bg = "white")
cat("✅ Mapa de tráfico maestro guardado en:", ruta_plot_trafico, "\n")

