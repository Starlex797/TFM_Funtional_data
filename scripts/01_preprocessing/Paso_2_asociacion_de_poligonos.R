# ==============================================================================
# STEP 2: CREATION OF THE SPATIO-TEMPORAL MASTER DATASET (DAILY OR HOURLY)
# ==============================================================================

library(data.table)
library(sf)
library(gstat)
library(ggplot2)
library(here)
library(lubridate)
source(here("R", "cleaning", "cleaning_functions.R"))
source(here("R", "interpolation", "FUNCIONES_INTERPOLACION.R"))
source(here("R", "utilities", "dictionaries.R"))
source(here("R", "utilities", "academic_quality_tables.R"))

# Deactivate the S2 engine to avoid issues with spatial operations
sf_use_s2(FALSE)

# ==============================================================================
# BLOCK 0: CONFIGURATION (Only change the paths here!)
# ==============================================================================

ANIO <- 2025
DIA_DIAGNOSTICO <- NULL
HORA_DIAGNOSTICO <- NULL

ruta_no2 <- here("data", "processed", "Contaminacion", "diario", paste0("aire_madrid_", ANIO, "_No2_trans_diarios1.rds"))
ruta_trafico <- here("data", "processed", "Trafico", "Diario_Barrio", ANIO, paste0("trafico_madrid_", ANIO, "_diario_barrio1.rds"))
ruta_meteo <- here("data", "processed", "Clima", "diario", paste0("meteo_madrid_", ANIO, "_diario5.rds"))

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
escala_temporal_en <- c(horario = "hourly", diario = "daily")[[escala_temporal]]
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
dir_fig_calidad <- here(
  "outputs", "figures", "EDA", "calidad de los datos",
  as.character(ANIO), escala_temporal
)
dir.create(dir_fig_calidad, recursive = TRUE, showWarnings = FALSE)

tabla_vacia <- function(mensaje) {
  data.table(Message = mensaje)
}

academic_name_en <- function(x) {
  labels <- c(
    FECHA = "Date",
    HORA = "Hour",
    ESTACION = "Station",
    NOM_TIPO = "Station type",
    DATO = "NO2",
    DATO_DIARIO = "Daily NO2",
    LOG_NO2_HORARIO = "Log NO2",
    LOG_NO2_DIARIO = "Log daily NO2",
    LONGITUD = "Longitude",
    LATITUD = "Latitude",
    distrito = "District",
    barrio = "Neighborhood",
    ID_TIEMPO = "Time ID",
    ID_DISTRITO = "District ID",
    X_km = "X coordinate (km)",
    Y_km = "Y coordinate (km)",
    intensidad = "Traffic intensity",
    ocupacion = "Traffic occupancy",
    carga = "Traffic load",
    Temperatura = "Temperature",
    Temperatura_log = "Log-transformed temperature",
    Humedad_Relativa = "Relative humidity",
    Precipitaciones = "Precipitation",
    Llueve = "Rain indicator",
    Presion_Barometrica = "Barometric pressure",
    Radiacion_Solar = "Solar radiation",
    Radiacion_Solar_log = "Log-transformed solar radiation",
    Velocidad_Viento = "Wind speed",
    Velocidad_Viento_sqrt = "Square-root-transformed wind speed"
  )
  values <- as.character(x)
  is_raw <- grepl("_raw$", values)
  base_values <- sub("_raw$", "", values)
  translated <- unname(labels[base_values])
  translated[is.na(translated)] <- gsub("_", " ", base_values[is.na(translated)])
  translated[is_raw] <- paste0(translated[is_raw], " (raw)")
  translated
}

interpolation_method_en <- function(x) {
  methods <- c(
    "Vecino Cercano" = "Nearest neighbour",
    "kNN" = "k-nearest neighbours",
    "IDW" = "Inverse distance weighting",
    "Ensemble" = "Ensemble"
  )
  translated <- unname(methods[as.character(x)])
  translated[is.na(translated)] <- as.character(x)[is.na(translated)]
  translated
}

formatear_tabla_academica <- function(tabla) {
  salida <- as.data.table(copy(tabla))
  for (col in names(salida)) {
    valores <- salida[[col]]
    if (inherits(valores, "Date")) {
      nuevo <- format(valores, "%Y-%m-%d")
    } else if (is.numeric(valores)) {
      # El formato se decide por COLUMNA, no celda a celda. Redondear cada
      # valor por su cuenta mezclaba "3.83" con "0" en la misma columna y, peor,
      # colapsaba a "0"/"-0" columnas de valores diminutos: la media de una
      # z-score (~1e-17) se imprimía como "-0" y dejaba el chequeo de
      # estandarización de la Tabla 3 sin ninguna información utilizable.
      finitos <- valores[is.finite(valores)]
      nuevo <- if (!length(finitos)) {
        rep("--", length(valores))
      } else if (all(finitos == round(finitos))) {
        sprintf("%.0f", valores)
      } else if (max(abs(finitos)) < 0.01) {
        sprintf("%.2e", valores)
      } else {
        sprintf("%.2f", valores)
      }
      nuevo[is.na(valores)] <- "--"
    } else {
      nuevo <- as.character(valores)
      nuevo[is.na(nuevo)] <- "--"
    }
    set(salida, j = col, value = nuevo)
  }
  salida
}

guardar_tabla_academica <- function(tabla, nombre_base, titulo, subtitulo, nota,
                                    widths, align, rows_per_page = 36L,
                                    font_size = 8.5, row_height = 0.24) {
  tabla <- formatear_tabla_academica(tabla)
  if (!nrow(tabla)) {
    tabla <- tabla_vacia("No records.")
    widths <- 3.0
    align <- "left"
  }
  paginas <- split(seq_len(nrow(tabla)), ceiling(seq_len(nrow(tabla)) / rows_per_page))
  n_paginas <- length(paginas)
  rutas <- character(n_paginas)

  for (i in seq_along(paginas)) {
    sufijo <- if (n_paginas > 1L) sprintf("_page_%02d", i) else ""
    ruta <- file.path(dir_fig_calidad, paste0(nombre_base, sufijo, ".png"))
    titulo_i <- if (n_paginas > 1L) {
      sprintf("%s (%d/%d)", titulo, i, n_paginas)
    } else {
      titulo
    }
    nota_i <- if (n_paginas > 1L) {
      paste(nota, "The table is paginated because it exceeds one printable panel.")
    } else {
      nota
    }
    booktabs_png(
      tabla[paginas[[i]]], ruta, titulo_i, subtitulo, nota_i,
      widths = widths, align = align, font_size = font_size,
      row_height = row_height
    )
    rutas[i] <- ruta
  }
  cat(sprintf("Academic table saved: %s (%d page(s))\n", nombre_base, n_paginas))
  invisible(rutas)
}

station_type_en <- function(x) {
  # El fichero de estaciones trae "Urbana tráfico" con tilde, así que buscar
  # "traf" no casaba y esa categoría se quedaba sin traducir en la tabla. Se
  # normaliza con la misma función que el resto de nombres (minúsculas y sin
  # tildes) antes de comparar.
  y <- limpiar_nombres(x)
  fifelse(
    grepl("suburbana", y), "Suburban",
    fifelse(
      grepl("fondo", y), "Urban background",
      fifelse(grepl("trafico", y), "Urban traffic", as.character(x))
    )
  )
}

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
tabla_tipologia <- dt_no2[, .(Stations = uniqueN(ESTACION)), by = NOM_TIPO]
tabla_tipologia[, NOM_TIPO := station_type_en(NOM_TIPO)]
setnames(tabla_tipologia, "NOM_TIPO", "Station type")
print(tabla_tipologia)
guardar_tabla_academica(
  tabla_tipologia,
  "Table_1_station_typology",
  "Table 1. Air-quality monitoring stations by station type",
  sprintf("Madrid NO2 network, %d, %s scale", ANIO, escala_temporal_en),
  "Station type is a fixed monitoring-site attribute used as a modelling covariate in the master dataset.",
  widths = c(2.4, 1.0),
  align = c("left", "right"),
  rows_per_page = 30L
)

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
tabla_config_clima <- copy(config_clima)
tabla_config_clima[, variable := academic_name_en(variable)]
tabla_config_clima[, metodo := interpolation_method_en(metodo)]
setnames(
  tabla_config_clima,
  c("variable", "metodo", "k", "p"),
  c("Climate variable", "Method", "Neighbours (k)", "Distance power (p)")
)
guardar_tabla_academica(
  tabla_config_clima,
  "Table_2_climate_interpolation_setup",
  "Table 2. Climate interpolation setup used in the master dataset",
  sprintf("Madrid, %d, %s scale", ANIO, escala_temporal_en),
  "The selected method and hyperparameters are applied variable by variable before joining climate covariates to NO2 stations.",
  widths = c(2.2, 1.3, 1.1, 1.2),
  align = c("left", "center", "right", "right"),
  rows_per_page = 30L
)
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
# Indicador binario calculado sobre la precipitacion original en milimetros,
# antes de estandarizarla: 0 = no llueve (< 1 mm), 1 = llueve (>= 1 mm).
# Si falta la precipitacion, el indicador tambien queda como NA.
dt_maestro[, Llueve := fifelse(
  is.na(Precipitaciones),
  NA_integer_,
  as.integer(Precipitaciones >= 1)
)]

# Transformaciones no lineales calculadas en las unidades originales.
# Se conservan las variables originales y se crean columnas nuevas.
TEMP_LOG_OFFSET_C <- 20
if (any(dt_maestro$Temperatura + TEMP_LOG_OFFSET_C <= 0, na.rm = TRUE)) {
  stop(
    "No se puede calcular Temperatura_log: existen temperaturas <= -",
    TEMP_LOG_OFFSET_C, " grados C."
  )
}

n_radiacion_negativa <- sum(dt_maestro$Radiacion_Solar < 0, na.rm = TRUE)
n_viento_negativo <- sum(dt_maestro$Velocidad_Viento < 0, na.rm = TRUE)

if (n_radiacion_negativa > 0L) {
  warning(
    n_radiacion_negativa,
    " valores negativos de Radiacion_Solar se sustituyen por 0 antes del logaritmo."
  )
}
if (n_viento_negativo > 0L) {
  warning(
    n_viento_negativo,
    " valores negativos de Velocidad_Viento se sustituyen por 0 antes de la raiz cuadrada."
  )
}

dt_maestro[, Temperatura_log := log(Temperatura + TEMP_LOG_OFFSET_C)]
dt_maestro[, Radiacion_Solar_log := log1p(pmax(Radiacion_Solar, 0))]
dt_maestro[, Velocidad_Viento_sqrt := sqrt(pmax(Velocidad_Viento, 0))]

# BLOCK 6: COVARIATE STANDARDIZATION (Z-SCORE)
# ==============================================================================
# Traffic
dt_maestro[, intensidad_raw := intensidad]
dt_maestro[, carga_raw := carga]
dt_maestro[, intensidad := scale(intensidad)[, 1]]
dt_maestro[, carga := scale(carga)[, 1]]

# Climate: variables interpoladas y transformaciones candidatas.
cols_clima_transformadas <- c(
  "Temperatura_log",
  "Radiacion_Solar_log",
  "Velocidad_Viento_sqrt"
)
cols_clima_std <- c(cols_clima, cols_clima_transformadas)

for (v in cols_clima_std) {
  raw_name <- paste0(v, "_raw")
  dt_maestro[, (raw_name) := get(v)]
  dt_maestro[, (v) := scale(get(v))[, 1]]
}

cat("\n--- Covariate standardization ---\n")
cols_std <- c("intensidad", "carga", cols_clima_std)
tabla_estandarizacion <- rbindlist(lapply(cols_std, function(v) {
  data.table(
    Covariate = academic_name_en(v),
    Mean = mean(dt_maestro[[v]], na.rm = TRUE),
    SD = sd(dt_maestro[[v]], na.rm = TRUE)
  )
}))
for (v in cols_std) {
  cat(sprintf(
    "%-25s mean = %7.4f | sd = %6.4f\n",
    v, mean(dt_maestro[[v]], na.rm = TRUE), sd(dt_maestro[[v]], na.rm = TRUE)
  ))
}
guardar_tabla_academica(
  tabla_estandarizacion,
  "Table_3_covariate_standardisation",
  "Table 3. Standardisation check for model covariates",
  sprintf("Master dataset, Madrid, %d, %s scale", ANIO, escala_temporal_en),
  "Traffic and interpolated climate covariates are standardised as z-scores before modelling.",
  widths = c(2.3, 1.1, 1.1),
  align = c("left", "right", "right"),
  rows_per_page = 30L
)

# Strict chronological and spatial ordering
setorderv(dt_maestro, intersect(c("FECHA", "HORA", "ID_TIEMPO", "ESTACION"), names(dt_maestro)))

# ==============================================================================
# BLOCK 7: SAVING AND QUALITY CONTROL
# ==============================================================================
dir.create(dirname(ruta_out_maestro), recursive = TRUE, showWarnings = FALSE)
saveRDS(dt_maestro, ruta_out_maestro)
if (interactive()) View(dt_maestro)
sf_use_s2(TRUE)

cat("\n✅ Unification completed successfully!\n")
cat("Total rows in the Master Dataset:", nrow(dt_maestro), "\n")
if ("ID_TIEMPO" %in% names(dt_maestro)) {
  cat("Are there any NAs in ID_TIEMPO?:", any(is.na(dt_maestro$ID_TIEMPO)), "\n")
}

# Show a summary adapted to the existing columns
cols_print <- intersect(c("FECHA", "HORA", "ESTACION", "distrito", "ID_TIEMPO", "ID_DISTRITO"), names(dt_maestro))
tabla_preview_maestro <- head(dt_maestro[, ..cols_print], 10L)
setnames(
  tabla_preview_maestro,
  old = intersect(c("FECHA", "HORA", "ESTACION", "distrito", "ID_TIEMPO", "ID_DISTRITO"), names(tabla_preview_maestro)),
  new = c("Date", "Hour", "Station", "District", "Time ID", "District ID")[
    c("FECHA", "HORA", "ESTACION", "distrito", "ID_TIEMPO", "ID_DISTRITO") %in% names(tabla_preview_maestro)
  ]
)
print(tabla_preview_maestro)
guardar_tabla_academica(
  tabla_preview_maestro,
  "Table_4_master_dataset_preview",
  "Table 4. Preview of the spatio-temporal master dataset",
  sprintf("First ten rows after sorting, Madrid, %d, %s scale", ANIO, escala_temporal_en),
  "The preview verifies the temporal, spatial and administrative keys retained in the master dataset.",
  widths = c(
    1.0, if ("Hour" %in% names(tabla_preview_maestro)) 0.7 else NULL,
    2.2, 1.7,
    if ("Time ID" %in% names(tabla_preview_maestro)) 0.8 else NULL,
    0.9
  ),
  align = c(
    "center", if ("Hour" %in% names(tabla_preview_maestro)) "right" else NULL,
    "left", "left",
    if ("Time ID" %in% names(tabla_preview_maestro)) "right" else NULL,
    "right"
  ),
  rows_per_page = 20L,
  font_size = 8.0
)

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

tabla_na_resumen <- copy(na_resumen)
setnames(
  tabla_na_resumen,
  c("columna", "n_NA", "pct", "instantes_afectados"),
  c("Column", "Missing rows", "Missing (%)", "Affected time points")
)
tabla_na_resumen[, Column := academic_name_en(Column)]
guardar_tabla_academica(
  tabla_na_resumen,
  "Table_5_missing_values_by_column",
  "Table 5. Missing values by column in the master dataset",
  sprintf(
    "%d rows; %d stations; %d time points; %s scale",
    n_filas, uniqueN(dt_maestro$ESTACION), n_instantes, escala_temporal_en
  ),
  "Percentages use the full master dataset as denominator. Affected time points count dates, or date-hour pairs in hourly scale, with at least one missing value.",
  widths = c(2.4, 1.1, 1.0, 1.4),
  align = c("left", "right", "right", "right"),
  rows_per_page = 34L,
  font_size = 8.0,
  row_height = 0.22
)

# 2. NO2: missing observations per station
col_no2 <- if ("DATO" %in% names(dt_maestro)) "DATO" else "DATO_DIARIO"
cat("\n--- NO2 NAs by station ---\n")
tabla_no2_na <- dt_maestro[is.na(get(col_no2)), .N, by = ESTACION][order(-N)]
setnames(tabla_no2_na, c("ESTACION", "N"), c("Station", "Missing rows"))
print(tabla_no2_na)
guardar_tabla_academica(
  tabla_no2_na,
  "Table_6_NO2_missing_by_station",
  "Table 6. Missing NO2 observations by monitoring station",
  sprintf("Master dataset, Madrid, %d, %s scale", ANIO, escala_temporal_en),
  "Rows are station-date observations in daily scale, or station-date-hour observations in hourly scale.",
  widths = c(2.6, 1.2),
  align = c("left", "right"),
  rows_per_page = 34L,
  font_size = 8.0,
  row_height = 0.22
)

# 3. Climate: which time points are missing, and in which variable?
#    Se recorren las SEIS variables interpoladas, no solo la temperatura. Con
#    una sola variable de diagnóstico la tabla salía vacía siempre que el hueco
#    estuviera en otra (presión, radiación, viento...), que es justo lo que pasa
#    en la mayoría de años: parecía que no faltaba nada de clima cuando la
#    Tabla 5 sí contaba cientos de filas sin dato.
cat("\n--- Climate NAs by time point and variable ---\n")
tabla_clima_na_fecha <- rbindlist(lapply(cols_clima, function(v) {
  filas_na <- dt_maestro[is.na(get(v))]
  if (!nrow(filas_na)) {
    return(NULL)
  }
  resumen <- filas_na[, .(N = .N), by = llaves_tiempo]
  resumen[, variable := v]
  resumen[]
}), fill = TRUE)

hay_hora <- "HORA" %in% llaves_tiempo
if (nrow(tabla_clima_na_fecha)) {
  setcolorder(tabla_clima_na_fecha, c("variable", llaves_tiempo, "N"))
  setorderv(tabla_clima_na_fecha, c("variable", llaves_tiempo))
  tabla_clima_na_fecha[, variable := academic_name_en(variable)]
  setnames(
    tabla_clima_na_fecha,
    c("variable", llaves_tiempo, "N"),
    c("Climate variable", academic_name_en(llaves_tiempo), "Rows missing")
  )
}
print(tabla_clima_na_fecha, nrows = Inf)
guardar_tabla_academica(
  tabla_clima_na_fecha,
  "Table_7_climate_missing_dates",
  "Table 7. Time points with missing interpolated climate values",
  sprintf(
    "All six interpolated climate covariates, Madrid, %d, %s scale",
    ANIO, escala_temporal_en
  ),
  "One row per climate covariate and time point without an interpolated value in the master dataset; the count is the number of station rows affected.",
  widths = c(2.0, 1.2, if (hay_hora) 0.7 else NULL, 1.2),
  align = c("left", "center", if (hay_hora) "right" else NULL, "right"),
  rows_per_page = 34L,
  font_size = 8.0,
  row_height = 0.22
)

# 4. Traffic: which master-dataset observations are missing?
cat("\n--- Traffic NAs by master-dataset observation ---\n")
cols_trafico_na <- c("ESTACION", "barrio", llaves_tiempo)
tabla_trafico_na <- dt_maestro[is.na(intensidad), ..cols_trafico_na]
setnames(tabla_trafico_na, cols_trafico_na, academic_name_en(cols_trafico_na))
print(tabla_trafico_na)
guardar_tabla_academica(
  tabla_trafico_na,
  "Table_8_traffic_missing_by_observation",
  "Table 8. Missing traffic covariates by master-dataset observation",
  sprintf("Master dataset, Madrid, %d, %s scale", ANIO, escala_temporal_en),
  "Each row is a master-dataset observation without traffic intensity after joining neighborhood traffic to NO2 stations.",
  widths = c(2.4, 2.4, 1.1, if ("Hour" %in% names(tabla_trafico_na)) 0.7 else NULL),
  align = c("left", "left", "center", if ("Hour" %in% names(tabla_trafico_na)) "right" else NULL),
  rows_per_page = 34L,
  font_size = 7.2,
  row_height = 0.19
)

# ==============================================================================
# BLOCK 9: NA DIAGNOSTICS IN THE MASTER DATASET — CONTAMINACIÓN AND TRÁFICO
# Same random day (dia) as the traffic section above.
# ==============================================================================

col_no2 <- if ("DATO" %in% names(dt_maestro)) "DATO" else "DATO_DIARIO"
mapa_barrios_lower <- mapa_barrios
dia <- if (is.null(DIA_DIAGNOSTICO)) {
  fechas_disponibles <- sort(unique(dt_maestro$FECHA))
  fechas_disponibles[ceiling(length(fechas_disponibles) / 2)]
} else {
  as.Date(DIA_DIAGNOSTICO)
}
hora <- if (escala_temporal == "horario") {
  if (is.null(HORA_DIAGNOSTICO)) 12L else as.integer(HORA_DIAGNOSTICO)
} else {
  NULL
}
if (!is.null(hora) && !hora %in% 1:24) {
  stop("HORA_DIAGNOSTICO must be an integer between 1 and 24.")
}
maestro_instante <- if (is.null(hora)) {
  dt_maestro[FECHA == dia]
} else {
  dt_maestro[FECHA == dia & HORA == hora]
}
if (!nrow(maestro_instante)) {
  stop("No master-dataset observations are available for the selected diagnostic time point.")
}
etiqueta_instante <- if (is.null(hora)) {
  format(dia, "%Y-%m-%d")
} else {
  sprintf("%s, hour %02d", format(dia, "%Y-%m-%d"), hora)
}
sufijo_instante <- if (is.null(hora)) "" else sprintf("_H%02d", hora)

# --- 9a. Contaminación (NO2): which stations are missing in dt_maestro? -------
no2_dia_sf <- as.data.frame(maestro_instante[, .(ESTACION, LONGITUD, LATITUD,
  no2 = get(col_no2)
)]) |>
  st_as_sf(coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)
no2_dia_sf$tiene_dato <- !is.na(no2_dia_sf$no2)

cat("\n--- NO2 coverage on", format(dia), "---\n")
cat("Estaciones con dato NO2:", sum(no2_dia_sf$tiene_dato), "/", nrow(no2_dia_sf), "\n")
cat("Estaciones SIN dato NO2:", sum(!no2_dia_sf$tiene_dato), "\n")

p_no2 <- ggplot() +
  geom_sf(data = mapa_barrios_lower, fill = "grey95", color = "white", linewidth = 0.2) +
  geom_sf(data = no2_dia_sf, aes(color = tiene_dato), size = 3.5, alpha = 0.85) +
  scale_color_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#d73027"),
    labels = c("TRUE" = "Observed", "FALSE" = "Missing (NA)"),
    name = "Status"
  ) +
  labs(
    title = paste("NO2 coverage in the master dataset -", etiqueta_instante),
    subtitle = paste0(
      sum(no2_dia_sf$tiene_dato), " stations observed  |  ",
      sum(!no2_dia_sf$tiene_dato), " stations missing (red)"
    ),
    caption = paste("Source: dt_maestro | Year:", ANIO, "| Scale:", escala_temporal_en)
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(color = "grey40", margin = margin(b = 8)),
    legend.position = "right"
  )

# Guardado dinámico del mapa de NO2
ruta_plot_no2 <- file.path(
  dir_fig_calidad,
  paste0("coverage_NO2_master_", ANIO, "_", escala_temporal, "_", format(dia, "%Y%m%d"), sufijo_instante, ".png")
)
ggsave(filename = ruta_plot_no2, plot = p_no2, width = 8, height = 6, dpi = 200, bg = "white")
cat("✅ Mapa de cobertura NO2 guardado en:", ruta_plot_no2, "\n")

# --- 9b. Tráfico: which barrios are missing in dt_maestro? --------------------
trafico_maestro_dia <- unique(maestro_instante[, .(barrio, intensidad = intensidad_raw)])

mapa_trafico_maestro <- merge(
  mapa_barrios_lower, as.data.frame(trafico_maestro_dia),
  by = "barrio", all.x = TRUE
)

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
    name = "Traffic intensity\n(veh/h)"
  ) +
  labs(
    title = paste("Traffic coverage in the master dataset -", etiqueta_instante),
    subtitle = paste0(
      sum(!is.na(mapa_trafico_maestro$intensidad)), " neighborhoods observed  |  ",
      sum(is.na(mapa_trafico_maestro$intensidad)), " neighborhoods missing (grey)"
    ),
    caption = paste("Source: dt_maestro | Year:", ANIO, "| Scale:", escala_temporal_en)
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(color = "grey40", margin = margin(b = 8)),
    legend.position = "right"
  )

# Guardado dinámico del mapa de tráfico maestro
ruta_plot_trafico <- file.path(
  dir_fig_calidad,
  paste0("coverage_traffic_master_", ANIO, "_", escala_temporal, "_", format(dia, "%Y%m%d"), sufijo_instante, ".png")
)
ggsave(filename = ruta_plot_trafico, plot = p_trafico_maestro, width = 8, height = 6, dpi = 200, bg = "white")
cat("✅ Mapa de tráfico maestro guardado en:", ruta_plot_trafico, "\n")
