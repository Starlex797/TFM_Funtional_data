# ==============================================================================
# Comparación diaria y mapas ganadores de las variables climáticas
# ==============================================================================

library(here)
source(here("R", "interpolation", "FUNCIONES_INTERPOLACION.R"))
library(ggplot2)
library(viridis)
library(tidyr)
library(gridExtra)

# ==============================================================================
# 1. Configuración
# ==============================================================================

ANIO <- 2025L
K_VECINOS <- 3L
BETAS_IDW <- c(1, 2)
FECHA_MAPA_REFERENCIA <- as.Date("2025-04-03")

# Criterio de días de lluvia.
MIN_ESTACIONES_PRECIP <- 7L
MIN_ESTACIONES_CON_LLUVIA <- 3L
UMBRAL_LLUVIA_MM <- 0.1

# Criterio de días de viento.
MIN_ESTACIONES_VIENTO <- 7L
PERCENTIL_VIENTO_VALIDACION <- 0.75
PERCENTIL_VIENTO_MAPA <- 0.90

VARIABLES_CLIMA <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  "Presion_Barometrica", "Radiacion_Solar", "Velocidad Viento"
)

ETIQUETAS_VARIABLES <- c(
  Temperatura = "Temperatura",
  Humedad_Relativa = "Humedad relativa",
  Precipitaciones = "Precipitaciones",
  Presion_Barometrica = "Presión barométrica",
  Radiacion_Solar = "Radiación solar",
  `Velocidad Viento` = "Velocidad del viento"
)

DIR_SALIDA_BASE <- here("outputs", "figures", "interpolacion_clima")

# ==============================================================================
# 2. Funciones auxiliares
# ==============================================================================

crear_directorio_versionado <- function(directorio_base) {
  dir.create(dirname(directorio_base), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(directorio_base)) {
    dir.create(directorio_base, recursive = TRUE)
    return(normalizePath(directorio_base, winslash = "/", mustWork = TRUE))
  }

  marca_temporal <- format(Sys.time(), "%Y%m%d_%H%M%S")
  candidato_base <- paste0(directorio_base, "_", marca_temporal)
  candidato <- candidato_base
  contador <- 2L

  while (dir.exists(candidato)) {
    candidato <- sprintf("%s_%02d", candidato_base, contador)
    contador <- contador + 1L
  }

  dir.create(candidato, recursive = TRUE)
  normalizePath(candidato, winslash = "/", mustWork = TRUE)
}

normalizar_coordenadas_madrid <- function(dt) {
  filas_latitud <- !is.na(dt$LATITUD) & dt$LATITUD < 35
  filas_longitud <- !is.na(dt$LONGITUD) & dt$LONGITUD > -1

  n_latitud <- sum(filas_latitud)
  n_longitud <- sum(filas_longitud)

  dt[filas_latitud, LATITUD :=
    LATITUD * 10^round(log10(40.4 / abs(LATITUD)))]
  dt[filas_longitud, LONGITUD := LONGITUD * 10]

  coordenadas_invalidas <- dt[
    is.na(LATITUD) | is.na(LONGITUD) |
      LATITUD < 40.2 | LATITUD > 40.7 |
      LONGITUD < -4.0 | LONGITUD > -3.3
  ]

  if (nrow(coordenadas_invalidas) > 0L) {
    stop(
      "Hay ", nrow(coordenadas_invalidas),
      " filas con coordenadas fuera del rango esperado de Madrid."
    )
  }

  cat(sprintf(
    "Coordenadas normalizadas: %d latitudes y %d longitudes.\n",
    n_latitud, n_longitud
  ))

  invisible(dt)
}

crear_tabla_png <- function(tabla, ruta, subtitulo) {
  tabla_grob <- tableGrob(
    as.data.frame(tabla),
    rows = NULL,
    theme = ttheme_minimal(base_size = 6.5)
  )
  titulo_grob <- grid::textGrob(
    "Comparación diaria de métodos de interpolación — LOOCV",
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  subtitulo_grob <- grid::textGrob(
    subtitulo,
    gp = grid::gpar(fontsize = 9)
  )
  composicion <- arrangeGrob(
    titulo_grob, subtitulo_grob, tabla_grob,
    ncol = 1,
    heights = grid::unit(c(0.4, 0.35, 4.8), "in")
  )

  ggsave(
    ruta,
    plot = composicion,
    width = 18,
    height = 5.6,
    dpi = 170,
    bg = "white"
  )
}

predecir_superficie <- function(
    metodo, sf_estaciones, rejilla_madrid, k_vecinos,
    pesos_ensemble = NULL) {
  formula_mapa <- Z ~ 1
  k_real <- min(k_vecinos, nrow(sf_estaciones))

  if (metodo == "Media") {
    return(rep(mean(sf_estaciones$Z), nrow(rejilla_madrid)))
  }

  if (metodo %in% c("Vecino Cercano", "Closest observation")) {
    modelo <- idw(
      formula_mapa, sf_estaciones, rejilla_madrid,
      nmax = 1, debug.level = 0
    )
    return(modelo$var1.pred)
  }

  if (metodo == "kNN") {
    modelo <- idw(
      formula_mapa, sf_estaciones, rejilla_madrid,
      nmax = k_real, idp = 0, debug.level = 0
    )
    return(modelo$var1.pred)
  }

  if (metodo %in% c("IDW beta=1", "IDW beta=2")) {
    beta <- if (metodo == "IDW beta=1") 1 else 2
    modelo <- idw(
      formula_mapa, sf_estaciones, rejilla_madrid,
      nmax = k_real, idp = beta, debug.level = 0
    )
    return(modelo$var1.pred)
  }

  # Ensemble: combinación ponderada de los tres interpoladores locales
  # (closest observation / 1-NN, IDW beta=1 y kNN) con pesos proporcionales a
  # 1/RMSE (pasados en pesos_ensemble: closest, IDW, kNN). No incluye la Media.
  if (metodo == "Ensemble") {
    if (is.null(pesos_ensemble) || length(pesos_ensemble) != 3L) {
      stop(
        "El método Ensemble requiere 'pesos_ensemble' de longitud 3 ",
        "(closest, IDW beta=1, kNN)."
      )
    }
    p_closest <- predecir_superficie(
      "Closest observation", sf_estaciones, rejilla_madrid, k_vecinos
    )
    p_idw <- predecir_superficie(
      "IDW beta=1", sf_estaciones, rejilla_madrid, k_vecinos
    )
    p_knn <- predecir_superficie(
      "kNN", sf_estaciones, rejilla_madrid, k_vecinos
    )
    return(
      pesos_ensemble[1] * p_closest +
        pesos_ensemble[2] * p_idw +
        pesos_ensemble[3] * p_knn
    )
  }

  stop("Método no reconocido: ", metodo)
}

# Pesos del ensemble para una variable, proporcionales a 1/RMSE de closest
# observation (1-NN), IDW beta=1 y kNN, tomados de la tabla LOOCV. Devuelve un
# vector de longitud 3 que suma 1, en el orden (closest, IDW beta=1, kNN).
pesos_ensemble_variable <- function(variable, tabla) {
  fila <- tabla[Variable == variable]
  rmses <- c(fila$RMSE_NN, fila$RMSE_IDW_b1, fila$RMSE_kNN)
  pesos <- 1 / pmax(rmses, 1e-9)
  as.numeric(pesos / sum(pesos))
}

descripcion_metodo <- function(metodo) {
  switch(
    metodo,
    "Media" = "Media de todas las estaciones",
    "Vecino Cercano" = "1-NN",
    "Closest observation" = "Closest observation (Voronoi)",
    "kNN" = sprintf("kNN, k=%d", K_VECINOS),
    "IDW beta=1" = sprintf("IDW, k=%d, beta=1", K_VECINOS),
    "IDW beta=2" = sprintf("IDW, k=%d, beta=2", K_VECINOS),
    "Ensemble" = sprintf(
      "Ensemble (pesos 1/RMSE de closest + IDW b=1 + kNN k=%d)",
      K_VECINOS
    ),
    metodo
  )
}

# ==============================================================================
# 3. Lectura y control de calidad
# ==============================================================================

ruta_datos <- here(
  "data", "processed", "Clima", "diario",
  sprintf("meteo_madrid_%d_diario_REAL.rds", ANIO)
)
dt_meteo <- readRDS(ruta_datos)
setDT(dt_meteo)

# Normalizar nombres con acentos para no depender del locale de Windows/R.
col_presion <- grep("^Presion Barom", names(dt_meteo), value = TRUE)
col_radiacion <- grep("Solar$", names(dt_meteo), value = TRUE)
if (length(col_presion) == 1L) {
  setnames(dt_meteo, col_presion, "Presion_Barometrica")
}
if (length(col_radiacion) == 1L) {
  setnames(dt_meteo, col_radiacion, "Radiacion_Solar")
}

normalizar_coordenadas_madrid(dt_meteo)

# ==============================================================================
# 4. Selección de días relevantes para precipitación y viento
# ==============================================================================

resumen_dias <- dt_meteo[, {
  precip_valida <- Precipitaciones[!is.na(Precipitaciones)]
  viento_valido <- `Velocidad Viento`[!is.na(`Velocidad Viento`)]

  list(
    N_Estaciones_Precip = length(precip_valida),
    N_Estaciones_Lluvia = sum(precip_valida > UMBRAL_LLUVIA_MM),
    Precipitacion_Media = if (length(precip_valida)) {
      mean(precip_valida)
    } else {
      NA_real_
    },
    N_Estaciones_Viento = length(viento_valido),
    Viento_Mediana = if (length(viento_valido)) {
      median(viento_valido)
    } else {
      NA_real_
    }
  )
}, by = FECHA]

fechas_lluvia <- resumen_dias[
  N_Estaciones_Precip >= MIN_ESTACIONES_PRECIP &
    N_Estaciones_Lluvia >= MIN_ESTACIONES_CON_LLUVIA,
  FECHA
]

if (!length(fechas_lluvia)) {
  stop("No se encontraron días de lluvia con el criterio configurado.")
}

dias_con_viento_suficiente <- resumen_dias[
  N_Estaciones_Viento >= MIN_ESTACIONES_VIENTO & !is.na(Viento_Mediana)
]
umbral_viento_p75 <- quantile(
  dias_con_viento_suficiente$Viento_Mediana,
  probs = PERCENTIL_VIENTO_VALIDACION,
  names = FALSE
)
umbral_viento_p90 <- quantile(
  dias_con_viento_suficiente$Viento_Mediana,
  probs = PERCENTIL_VIENTO_MAPA,
  names = FALSE
)
fechas_viento <- dias_con_viento_suficiente[
  Viento_Mediana >= umbral_viento_p75,
  FECHA
]

if (!length(fechas_viento)) {
  stop("No se encontraron días de viento con el criterio configurado.")
}

# Mantener el 3 de abril para precipitación si cumple el criterio acordado.
if (FECHA_MAPA_REFERENCIA %in% fechas_lluvia) {
  fecha_mapa_precipitacion <- FECHA_MAPA_REFERENCIA
} else {
  fecha_mapa_precipitacion <- resumen_dias[
    FECHA %in% fechas_lluvia
  ][which.max(Precipitacion_Media), FECHA]
}

# Para viento, seleccionar el día cuya mediana esté más cerca del percentil 90.
fecha_mapa_viento <- dias_con_viento_suficiente[
  which.min(abs(Viento_Mediana - umbral_viento_p90)),
  FECHA
]

cat(sprintf(
  "Días de lluvia: %d | Días de viento >= P75: %d.\n",
  length(fechas_lluvia), length(fechas_viento)
))
cat(sprintf(
  "Fechas de mapa: precipitación %s | viento %s.\n",
  fecha_mapa_precipitacion, fecha_mapa_viento
))

# ==============================================================================
# 5. Validación diaria
# ==============================================================================

# Resultado anual completo: se conserva para comparar con el procedimiento
# anterior y para seleccionar temperatura, humedad, presión y radiación.
tabla_anual <- comparar_interpolaciones_loocv(
  dt_meteo = dt_meteo,
  variables = VARIABLES_CLIMA,
  k_vecinos = K_VECINOS,
  betas_idw = BETAS_IDW
)

# Para precipitación se evalúan días de lluvia completos, conservando los ceros
# observados dentro de esos días.
tabla_precipitacion <- comparar_interpolaciones_loocv(
  dt_meteo = dt_meteo[FECHA %in% fechas_lluvia],
  variables = "Precipitaciones",
  k_vecinos = K_VECINOS,
  betas_idw = BETAS_IDW
)

# Para viento se evalúan los días con mediana diaria por encima del percentil 75.
tabla_viento <- comparar_interpolaciones_loocv(
  dt_meteo = dt_meteo[FECHA %in% fechas_viento],
  variables = "Velocidad Viento",
  k_vecinos = K_VECINOS,
  betas_idw = BETAS_IDW
)

tabla_decision <- rbindlist(list(
  tabla_anual[
    !Variable %chin% c("Precipitaciones", "Velocidad Viento")
  ],
  tabla_precipitacion,
  tabla_viento
), use.names = TRUE)
tabla_decision <- tabla_decision[match(VARIABLES_CLIMA, Variable)]
tabla_decision[, Muestra := fifelse(
  Variable == "Precipitaciones",
  "Días de lluvia",
  fifelse(
    Variable == "Velocidad Viento",
    "Días de viento >= P75",
    "Todo 2025"
  )
)]
setcolorder(
  tabla_decision,
  c(
    "Variable", "Muestra", "N_Estaciones", "N_Dias", "N_Predicciones",
    grep("^RMSE_", names(tabla_decision), value = TRUE),
    grep("^MAE_", names(tabla_decision), value = TRUE),
    "Mejor_Metodo"
  )
)

cat("\nTabla utilizada para decidir el método de cada variable:\n")
print(tabla_decision)

# ==============================================================================
# 5b. Validación de UN ÚNICO DÍA por variable
# ==============================================================================
# La tabla anterior promedia el RMSE sobre muchos días (todo 2025, o todos los
# días de lluvia/viento). Aquí evaluamos la LOOCV en un SOLO día por variable
# para ver la capacidad predictiva puntual, sin promediar entre fechas. El día
# elegido es el mismo que se usa en los mapas: día de referencia para las
# variables continuas, un día con lluvia para precipitación y un día con viento
# para la velocidad del viento.

dias_un_dia <- data.table(
  Variable = VARIABLES_CLIMA,
  FECHA = as.Date(c(
    FECHA_MAPA_REFERENCIA,     # Temperatura
    FECHA_MAPA_REFERENCIA,     # Humedad_Relativa
    fecha_mapa_precipitacion,  # Precipitaciones (día con lluvia)
    FECHA_MAPA_REFERENCIA,     # Presion_Barometrica
    FECHA_MAPA_REFERENCIA,     # Radiacion_Solar
    fecha_mapa_viento          # Velocidad Viento (día con viento)
  ))
)

tabla_un_dia <- rbindlist(
  lapply(seq_len(nrow(dias_un_dia)), function(i) {
    variable <- dias_un_dia$Variable[i]
    dia <- dias_un_dia$FECHA[i]
    res <- comparar_interpolaciones_loocv(
      dt_meteo = dt_meteo[FECHA == dia],
      variables = variable,
      k_vecinos = K_VECINOS,
      betas_idw = BETAS_IDW
    )
    res[, Dia := dia]
    res
  }),
  use.names = TRUE
)
tabla_un_dia <- tabla_un_dia[match(VARIABLES_CLIMA, Variable)]
# N_Dias es siempre 1 en esta tabla; se elimina por redundante.
tabla_un_dia[, N_Dias := NULL]
setcolorder(
  tabla_un_dia,
  c(
    "Variable", "Dia", "N_Estaciones", "N_Predicciones",
    grep("^RMSE_", names(tabla_un_dia), value = TRUE),
    grep("^MAE_", names(tabla_un_dia), value = TRUE),
    "Mejor_Metodo"
  )
)

cat("\nValidación de un único día por variable (sin promediar fechas):\n")
print(tabla_un_dia)

# ==============================================================================
# 6. Directorio versionado y salidas tabulares
# ==============================================================================

DIR_SALIDA <- crear_directorio_versionado(DIR_SALIDA_BASE)
cat("\nDirectorio de esta ejecución:\n", DIR_SALIDA, "\n")

fwrite(
  tabla_decision,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos.csv")
)
saveRDS(
  tabla_decision,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos.rds")
)
fwrite(
  tabla_anual,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos_anual.csv")
)
saveRDS(
  tabla_anual,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos_anual.rds")
)
fwrite(
  tabla_un_dia,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos_un_dia.csv")
)
saveRDS(
  tabla_un_dia,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos_un_dia.rds")
)
fwrite(
  resumen_dias[FECHA %in% fechas_lluvia],
  file.path(DIR_SALIDA, "dias_seleccionados_precipitacion.csv")
)
fwrite(
  resumen_dias[FECHA %in% fechas_viento],
  file.path(DIR_SALIDA, "dias_seleccionados_viento.csv")
)

subtitulo_tabla <- paste0(
  "k=", K_VECINOS, " | IDW beta=1 y 2 | ",
  "precipitación: días de lluvia | viento: días >= P75"
)
crear_tabla_png(
  tabla_decision,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos.png"),
  subtitulo_tabla
)

# Gráfico RMSE de la tabla utilizada para seleccionar los métodos.
df_rmse <- tabla_decision[, .(
  Variable,
  Media = RMSE_Media,
  `1-NN` = RMSE_NN,
  `Voronoi (closest)` = RMSE_Voronoi,
  `kNN (k=3)` = RMSE_kNN,
  `IDW (k=3, beta=1)` = RMSE_IDW_b1,
  `IDW (k=3, beta=2)` = RMSE_IDW_b2,
  Ensemble = RMSE_Ensemble
)]
df_rmse_long <- pivot_longer(
  df_rmse,
  cols = -Variable,
  names_to = "Metodo",
  values_to = "RMSE"
)
df_rmse_long$Metodo <- factor(
  df_rmse_long$Metodo,
  levels = c(
    "Media", "1-NN", "Voronoi (closest)", "kNN (k=3)",
    "IDW (k=3, beta=1)", "IDW (k=3, beta=2)", "Ensemble"
  )
)

p_rmse <- ggplot(df_rmse_long, aes(Variable, RMSE, fill = Metodo)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", RMSE)),
    position = position_dodge(width = 0.78),
    vjust = -0.35,
    size = 2.7
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "RMSE por variable y método de interpolación",
    subtitle = paste(
      "Todo 2025 excepto precipitación (días de lluvia)",
      "y viento (días >= P75)"
    ),
    x = NULL,
    y = "RMSE",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(
  file.path(DIR_SALIDA, "grafico_rmse_metodos.png"),
  plot = p_rmse,
  width = 11,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# ==============================================================================
# 7. Rejilla espacial y mapas de los métodos ganadores
# ==============================================================================

mapa_distritos <- st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"),
  quiet = TRUE
) |>
  st_make_valid()
mapa_distritos_utm <- st_transform(mapa_distritos, 25830)

rejilla <- st_make_grid(
  mapa_distritos_utm,
  n = c(100, 100),
  what = "centers"
) |>
  st_as_sf()
rejilla_madrid <- st_intersection(
  rejilla,
  st_union(mapa_distritos_utm)
)
coords_rejilla <- st_coordinates(rejilla_madrid)
dx <- diff(sort(unique(coords_rejilla[, 1])))[1]
dy <- diff(sort(unique(coords_rejilla[, 2])))[1]

bordes_utm <- st_geometry(mapa_distritos_utm) |>
  st_cast("MULTILINESTRING") |>
  st_cast("LINESTRING")
bordes_coords <- as.data.frame(st_coordinates(bordes_utm))

tabla_mapas <- data.table(
  Variable = VARIABLES_CLIMA,
  FECHA = as.Date(c(
    FECHA_MAPA_REFERENCIA,
    FECHA_MAPA_REFERENCIA,
    fecha_mapa_precipitacion,
    FECHA_MAPA_REFERENCIA,
    FECHA_MAPA_REFERENCIA,
    fecha_mapa_viento
  )),
  Metodo = tabla_decision$Mejor_Metodo
)
tabla_mapas[, `:=`(
  N_Estaciones = NA_integer_,
  Min_Observado = NA_real_,
  Max_Observado = NA_real_
)]

for (i in seq_len(nrow(tabla_mapas))) {
  variable <- tabla_mapas$Variable[i]
  fecha <- tabla_mapas$FECHA[i]
  metodo <- tabla_mapas$Metodo[i]
  etiqueta <- unname(ETIQUETAS_VARIABLES[variable])

  dt_instante <- dt_meteo[
    FECHA == fecha & !is.na(get(variable))
  ]

  if (nrow(dt_instante) <= K_VECINOS) {
    warning(
      "No hay estaciones suficientes para ", variable,
      " en ", fecha, ". Se omite el mapa."
    )
    next
  }

  sf_estaciones <- st_as_sf(
    dt_instante,
    coords = c("LONGITUD", "LATITUD"),
    crs = 4326
  ) |>
    st_transform(25830)
  sf_estaciones$Z <- sf_estaciones[[variable]]

  prediccion <- predecir_superficie(
    metodo = metodo,
    sf_estaciones = sf_estaciones,
    rejilla_madrid = rejilla_madrid,
    k_vecinos = K_VECINOS,
    pesos_ensemble = pesos_ensemble_variable(variable, tabla_decision)
  )

  df_superficie <- data.frame(
    X = coords_rejilla[, 1],
    Y = coords_rejilla[, 2],
    Valor = prediccion
  )
  df_estaciones <- as.data.frame(st_coordinates(sf_estaciones))
  df_estaciones$Valor <- sf_estaciones$Z

  p_mapa <- ggplot(df_superficie, aes(X, Y, fill = Valor)) +
    geom_tile(width = dx, height = dy) +
    geom_path(
      data = bordes_coords,
      aes(x = X, y = Y, group = L1),
      color = "white",
      linewidth = 0.35,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = df_estaciones,
      aes(X, Y, fill = Valor),
      color = "black",
      shape = 21,
      size = 2.7,
      stroke = 0.7,
      inherit.aes = FALSE
    ) +
    scale_fill_viridis_c(
      option = "plasma",
      begin = 0.05,
      end = 0.95,
      name = etiqueta
    ) +
    coord_equal() +
    labs(
      title = paste("Mapa del método ganador:", etiqueta),
      subtitle = sprintf(
        "%s | %s | %d estaciones",
        format(fecha, "%d %b %Y"),
        descripcion_metodo(metodo),
        nrow(sf_estaciones)
      ),
      caption = paste(
        "Superficie = predicción del método ganador.",
        "Puntos con borde negro = observaciones reales."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "right"
    )

  nombre_archivo <- gsub("[^a-zA-Z0-9]", "_", variable)
  ggsave(
    file.path(
      DIR_SALIDA,
      paste0("mapa_ganador_", nombre_archivo, ".png")
    ),
    plot = p_mapa,
    width = 9,
    height = 7,
    dpi = 300,
    bg = "white"
  )

  tabla_mapas[i, `:=`(
    N_Estaciones = nrow(sf_estaciones),
    Min_Observado = min(sf_estaciones$Z),
    Max_Observado = max(sf_estaciones$Z)
  )]

  cat(sprintf(
    "Mapa guardado: %s | %s | %s.\n",
    variable, metodo, fecha
  ))
}

fwrite(
  tabla_mapas,
  file.path(DIR_SALIDA, "metodos_y_fechas_mapas.csv")
)

# ==============================================================================
# 7b. Mapas comparativos de 4 paneles (closest, IDW, kNN, ensemble)
# ==============================================================================
# Para cada variable se genera una figura con las predicciones de los tres
# interpoladores locales que componen el ensemble (closest observation, IDW
# beta=1 y nearest neighbors) más el propio ensemble, con pesos proporcionales
# a 1/RMSE. Reproduce el formato de comparación de métodos (4 paneles).

# Métodos mostrados y etiquetas de los paneles (orden = closest, IDW, kNN, ens).
METODOS_ENSEMBLE <- c(
  "Closest observation",
  "IDW beta=1",
  "kNN",
  "Ensemble"
)
ETIQUETAS_PANEL_ENSEMBLE <- c(
  "Closest observation" = "Closest observation (Voronoi / 1-NN)",
  "IDW beta=1" = "IDW (k=3, beta=1)",
  "kNN" = sprintf("Nearest neighbors (kNN, k=%d)", K_VECINOS),
  "Ensemble" = "Ensemble (pesos 1/RMSE)"
)

for (i in seq_len(nrow(tabla_mapas))) {
  variable <- tabla_mapas$Variable[i]
  fecha <- tabla_mapas$FECHA[i]
  etiqueta <- unname(ETIQUETAS_VARIABLES[variable])

  dt_instante <- dt_meteo[
    FECHA == fecha & !is.na(get(variable))
  ]

  if (nrow(dt_instante) <= K_VECINOS) {
    warning(
      "No hay estaciones suficientes para el ensemble de ", variable,
      " en ", fecha, ". Se omite la figura de 4 paneles."
    )
    next
  }

  sf_estaciones <- st_as_sf(
    dt_instante,
    coords = c("LONGITUD", "LATITUD"),
    crs = 4326
  ) |>
    st_transform(25830)
  sf_estaciones$Z <- sf_estaciones[[variable]]

  df_estaciones <- as.data.frame(st_coordinates(sf_estaciones))
  df_estaciones$Valor <- sf_estaciones$Z

  # Un data.frame largo con las cuatro superficies apiladas para el facet.
  df_paneles <- rbindlist(lapply(METODOS_ENSEMBLE, function(metodo) {
    prediccion <- predecir_superficie(
      metodo = metodo,
      sf_estaciones = sf_estaciones,
      rejilla_madrid = rejilla_madrid,
      k_vecinos = K_VECINOS,
      pesos_ensemble = pesos_ensemble_variable(variable, tabla_decision)
    )
    data.table(
      X = coords_rejilla[, 1],
      Y = coords_rejilla[, 2],
      Valor = prediccion,
      Metodo = factor(
        unname(ETIQUETAS_PANEL_ENSEMBLE[metodo]),
        levels = unname(ETIQUETAS_PANEL_ENSEMBLE[METODOS_ENSEMBLE])
      )
    )
  }))

  p_ensemble <- ggplot(df_paneles, aes(X, Y, fill = Valor)) +
    geom_tile(width = dx, height = dy) +
    geom_path(
      data = bordes_coords,
      aes(x = X, y = Y, group = L1),
      color = "white",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = df_estaciones,
      aes(X, Y, fill = Valor),
      color = "black",
      shape = 21,
      size = 2,
      stroke = 0.6,
      inherit.aes = FALSE
    ) +
    scale_fill_viridis_c(
      option = "plasma",
      begin = 0.05,
      end = 0.95,
      name = etiqueta
    ) +
    facet_wrap(~Metodo, ncol = 2) +
    coord_equal() +
    labs(
      title = paste("Comparación de métodos de interpolación:", etiqueta),
      subtitle = sprintf(
        "%s | %d estaciones | escala común entre paneles",
        format(fecha, "%d %b %Y"),
        nrow(sf_estaciones)
      ),
      caption = paste(
        "Superficies = predicción de cada método sobre la misma rejilla.",
        "Puntos con borde negro = observaciones reales."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )

  nombre_archivo <- gsub("[^a-zA-Z0-9]", "_", variable)
  ggsave(
    file.path(
      DIR_SALIDA,
      paste0("mapa_ensemble_", nombre_archivo, ".png")
    ),
    plot = p_ensemble,
    width = 11,
    height = 9,
    dpi = 300,
    bg = "white"
  )

  cat(sprintf(
    "Mapa de 4 paneles (ensemble) guardado: %s | %s.\n",
    variable, fecha
  ))
}

# ==============================================================================
# 8. Metadatos de la ejecución
# ==============================================================================

lineas_metadatos <- c(
  paste("Fecha y hora de ejecución:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Directorio de salida:", DIR_SALIDA),
  paste("Directorio base histórico:", DIR_SALIDA_BASE),
  paste("Periodo anual:", min(dt_meteo$FECHA), "a", max(dt_meteo$FECHA)),
  paste("Vecinos comunes para kNN e IDW:", K_VECINOS),
  paste("Betas IDW comparadas:", paste(BETAS_IDW, collapse = ", ")),
  paste(
    "Closest observation (Voronoi): equivalente numéricamente a 1-NN;",
    "se informa como formulación geométrica separada."
  ),
  paste(
    "Criterio precipitación: al menos", MIN_ESTACIONES_PRECIP,
    "estaciones disponibles y", MIN_ESTACIONES_CON_LLUVIA,
    "estaciones por encima de", UMBRAL_LLUVIA_MM, "mm"
  ),
  paste("Número de días de lluvia:", length(fechas_lluvia)),
  paste(
    "Criterio viento: mediana diaria >= percentil",
    100 * PERCENTIL_VIENTO_VALIDACION,
    sprintf("(umbral %.4f)", umbral_viento_p75)
  ),
  paste("Número de días de viento seleccionados:", length(fechas_viento)),
  paste("Fecha del mapa de precipitación:", fecha_mapa_precipitacion),
  paste("Fecha del mapa de viento:", fecha_mapa_viento),
  "",
  "Métodos ganadores y muestra de validación:",
  paste0(
    "- ", tabla_decision$Variable,
    ": ", tabla_decision$Mejor_Metodo,
    " [", tabla_decision$Muestra, "]"
  ),
  "",
  paste(
    "Nota: esta versión no genera mapas de errores de interpolación.",
    "Los resultados de carpetas anteriores no se han modificado."
  )
)

writeLines(
  lineas_metadatos,
  file.path(DIR_SALIDA, "METADATOS_EJECUCION.txt"),
  useBytes = TRUE
)

cat("\nEjecución terminada correctamente.\n")
cat("Resultados guardados en:\n", DIR_SALIDA, "\n")
