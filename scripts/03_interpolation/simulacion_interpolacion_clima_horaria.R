# ==============================================================================
# Comparacion horaria de los metodos de interpolacion climatica
# ==============================================================================
#
# Este script reproduce a escala horaria la comparacion realizada a escala
# diaria. Cada observacion se valida usando exclusivamente las otras estaciones
# disponibles en la misma fecha y hora.
#
# Metodos:
#   - Media de las estaciones restantes
#   - Vecino mas cercano (1-NN)
#   - Closest observation mediante Voronoi (equivalente a 1-NN)
#   - kNN con k = 3
#   - IDW con k = 3 y beta = 1
#   - IDW con k = 3 y beta = 2
#
# Importante:
#   Si el archivo contiene indicadores imp_<variable>, los valores imputados se
#   enmascaran y la validacion utiliza solamente observaciones reales.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(ggplot2)
  library(gridExtra)
})

# ==============================================================================
# 1. Configuracion
# ==============================================================================

ANIO <- 2025L
K_VECINOS <- 3L
BETAS_IDW <- c(1, 2)
N_BOOTSTRAP_BLOQUES <- 3000L
SEMILLA_BOOTSTRAP <- 2025L

# Criterio para seleccionar horas de lluvia.
MIN_ESTACIONES_PRECIP <- 7L
MIN_ESTACIONES_CON_LLUVIA <- 3L
UMBRAL_LLUVIA_MM <- 0.1

# Criterio para seleccionar horas con radiación solar informativa.
MIN_ESTACIONES_RADIACION <- 7L
UMBRAL_RADIACION_WM2 <- 10

# Criterio para seleccionar horas de viento.
MIN_ESTACIONES_VIENTO <- 7L
PERCENTIL_VIENTO_VALIDACION <- 0.75

VARIABLES_CLIMA <- c(
  "Temperatura",
  "Humedad_Relativa",
  "Precipitaciones",
  "Presion_Barometrica",
  "Radiacion_Solar",
  "Velocidad Viento"
)

DIR_SALIDA_BASE <- here(
  "output", "figures", "interpolacion_clima_horaria"
)

# ==============================================================================
# 2. Funciones auxiliares
# ==============================================================================

crear_directorio_versionado <- function(directorio_base) {
  dir.create(
    dirname(directorio_base),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(directorio_base)) {
    dir.create(directorio_base, recursive = TRUE)
    return(normalizePath(
      directorio_base,
      winslash = "/",
      mustWork = TRUE
    ))
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

  dt[
    filas_latitud,
    LATITUD := LATITUD * 10^round(log10(40.4 / abs(LATITUD)))
  ]
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
    n_latitud,
    n_longitud
  ))

  invisible(dt)
}

crear_tabla_png <- function(tabla, ruta, subtitulo) {
  tabla_grob <- tableGrob(
    as.data.frame(tabla),
    rows = NULL,
    theme = ttheme_minimal(base_size = 6.3)
  )
  titulo_grob <- grid::textGrob(
    "Comparación horaria de métodos de interpolación — LOOCV",
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  subtitulo_grob <- grid::textGrob(
    subtitulo,
    gp = grid::gpar(fontsize = 9)
  )
  composicion <- arrangeGrob(
    titulo_grob,
    subtitulo_grob,
    tabla_grob,
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

# ==============================================================================
# 3. Validacion LOOCV horaria
# ==============================================================================

comparar_interpolaciones_loocv_horaria <- function(
    dt_meteo,
    variables,
    k_vecinos = 3L,
    betas_idw = c(1, 2)) {
  if (k_vecinos < 1L) {
    stop("k_vecinos debe ser al menos 1.")
  }

  betas_idw <- sort(unique(as.numeric(betas_idw)))
  if (!identical(betas_idw, c(1, 2))) {
    stop("betas_idw debe contener exactamente los valores 1 y 2.")
  }

  # Métodos base con acumulación incremental de SSE. El Ensemble no está aquí:
  # se calcula al final con pesos 1/RMSE (que requieren el RMSE global de cada
  # método base) a partir de las predicciones locales guardadas.
  metodos <- c(
    "Media",
    "Vecino Cercano",
    "Closest observation",
    "kNN",
    "IDW beta=1",
    "IDW beta=2"
  )
  resultados_globales <- list()
  detalles_diarios_globales <- list()

  cat("\nVALIDACIÓN LOOCV HORARIA\n")
  cat(
    "kNN e IDW: k =", k_vecinos,
    "| IDW beta = 1 y 2\n"
  )

  for (var_name in variables) {
    tiempo_inicio <- Sys.time()
    cat("\nVariable:", var_name, "\n")

    dt_var <- copy(dt_meteo)[
      !is.na(get(var_name)),
      .(
        FECHA,
        HORA,
        ESTACION,
        LONGITUD,
        LATITUD,
        VALOR = as.numeric(get(var_name))
      )
    ]

    # Estaciones co-localizadas: si dos comparten la ubicación exacta se
    # conserva solo una (la de más observaciones de esta variable), para que la
    # LOOCV no prediga una estación con su gemela a distancia cero.
    conteo_estaciones <- dt_var[
      ,
      .(n = .N, LON = first(LONGITUD), LAT = first(LATITUD)),
      by = ESTACION
    ]
    setorder(conteo_estaciones, LON, LAT, -n)
    estaciones_conservadas <- conteo_estaciones[
      !duplicated(paste(LON, LAT)), ESTACION
    ]
    if (length(estaciones_conservadas) < uniqueN(dt_var$ESTACION)) {
      dt_var <- dt_var[ESTACION %in% estaciones_conservadas]
    }

    n_estaciones <- uniqueN(dt_var$ESTACION)
    if (n_estaciones <= k_vecinos) {
      warning(
        "Se omite ", var_name, ": solo hay ", n_estaciones,
        " estaciones."
      )
      next
    }

    # Las coordenadas son fijas por estación. Se proyectan una única vez y la
    # matriz de distancias se reutiliza en todos los instantes horarios.
    dt_coords <- dt_var[
      ,
      .(
        LONGITUD = first(LONGITUD),
        LATITUD = first(LATITUD)
      ),
      by = ESTACION
    ]
    sf_coords <- st_as_sf(
      dt_coords,
      coords = c("LONGITUD", "LATITUD"),
      crs = 4326
    ) |>
      st_transform(25830)
    xy <- st_coordinates(sf_coords)
    dt_coords[, `:=`(X = xy[, 1], Y = xy[, 2])]
    matriz_distancias_estaciones <- as.matrix(
      dist(as.matrix(dt_coords[, .(X, Y)]))
    )
    rownames(matriz_distancias_estaciones) <- dt_coords$ESTACION
    colnames(matriz_distancias_estaciones) <- dt_coords$ESTACION

    dt_var <- dt_coords[dt_var, on = "ESTACION"]
    setorder(dt_var, FECHA, HORA, ESTACION)

    grupos <- dt_var[
      ,
      .(
        INICIO = first(.I),
        FIN = last(.I)
      ),
      by = .(FECHA, HORA)
    ]

    suma_error_cuadrado <- setNames(
      numeric(length(metodos)),
      metodos
    )
    suma_error_absoluto <- setNames(
      numeric(length(metodos)),
      metodos
    )
    fechas_variable <- sort(unique(dt_var$FECHA))
    n_predicciones_dia <- integer(length(fechas_variable))
    sse_dia <- matrix(
      0,
      nrow = length(fechas_variable),
      ncol = length(metodos),
      dimnames = list(NULL, metodos)
    )
    # Predicciones locales por instante para construir el ensemble al final.
    predicciones_ens <- vector("list", nrow(grupos))
    n_predicciones <- 0L
    n_instantes_validos <- 0L

    for (j in seq_len(nrow(grupos))) {
      dt_instante <- dt_var[
        grupos$INICIO[j]:grupos$FIN[j]
      ]
      n_instante <- nrow(dt_instante)

      # Tras ocultar una estación deben quedar al menos k vecinos.
      if (n_instante <= k_vecinos) {
        next
      }

      valores <- dt_instante$VALOR
      estaciones_instante <- dt_instante$ESTACION
      distancias <- matriz_distancias_estaciones[
        estaciones_instante,
        estaciones_instante,
        drop = FALSE
      ]
      diag(distancias) <- Inf

      pred_media <- (sum(valores) - valores) / (n_instante - 1L)
      pred_nn <- numeric(n_instante)
      pred_knn <- numeric(n_instante)
      pred_idw_b1 <- numeric(n_instante)
      pred_idw_b2 <- numeric(n_instante)

      for (i in seq_len(n_instante)) {
        idx_vecinos <- order(distancias[i, ])[seq_len(k_vecinos)]
        valores_vecinos <- valores[idx_vecinos]
        distancias_vecinos <- distancias[i, idx_vecinos]

        pred_nn[i] <- valores_vecinos[1]
        pred_knn[i] <- mean(valores_vecinos)

        pesos_b1 <- 1 / pmax(
          distancias_vecinos,
          1
        )^betas_idw[1]
        pesos_b2 <- 1 / pmax(
          distancias_vecinos,
          1
        )^betas_idw[2]
        pred_idw_b1[i] <- sum(
          pesos_b1 * valores_vecinos
        ) / sum(pesos_b1)
        pred_idw_b2[i] <- sum(
          pesos_b2 * valores_vecinos
        ) / sum(pesos_b2)
      }

      predicciones <- list(
        "Media" = pred_media,
        "Vecino Cercano" = pred_nn,
        "Closest observation" = pred_nn,
        "kNN" = pred_knn,
        "IDW beta=1" = pred_idw_b1,
        "IDW beta=2" = pred_idw_b2
      )

      indice_dia <- match(
        grupos$FECHA[j],
        fechas_variable
      )
      for (metodo in metodos) {
        error <- valores - predicciones[[metodo]]
        error_cuadrado_instante <- sum(error^2)
        suma_error_cuadrado[metodo] <-
          suma_error_cuadrado[metodo] + error_cuadrado_instante
        suma_error_absoluto[metodo] <-
          suma_error_absoluto[metodo] + sum(abs(error))
        sse_dia[indice_dia, metodo] <-
          sse_dia[indice_dia, metodo] + error_cuadrado_instante
      }

      # Se guardan las predicciones de los tres interpoladores locales para
      # construir el ensemble ponderado 1/RMSE al final del bucle.
      predicciones_ens[[j]] <- data.table(
        indice_dia = indice_dia,
        Real = valores,
        NN = pred_nn,
        IDW = pred_idw_b1,
        kNN = pred_knn
      )

      n_predicciones <- n_predicciones + n_instante
      n_predicciones_dia[indice_dia] <-
        n_predicciones_dia[indice_dia] + n_instante
      n_instantes_validos <- n_instantes_validos + 1L
    }

    if (n_predicciones == 0L) {
      warning(
        "Se omite ", var_name,
        ": no hay instantes con estaciones suficientes."
      )
      next
    }

    rmses <- sqrt(suma_error_cuadrado / n_predicciones)
    maes <- suma_error_absoluto / n_predicciones

    # Ensemble ponderado: pesos proporcionales a 1/RMSE de los tres
    # interpoladores locales (closest observation / 1-NN, IDW beta=1 y kNN).
    # No incluye la Media por no ser un interpolador espacial.
    dt_ens <- rbindlist(predicciones_ens)
    pesos_ens <- 1 / pmax(
      c(rmses["Vecino Cercano"], rmses["IDW beta=1"], rmses["kNN"]),
      1e-9
    )
    pesos_ens <- pesos_ens / sum(pesos_ens)
    dt_ens[, Pred_Ensemble :=
      pesos_ens[1] * NN + pesos_ens[2] * IDW + pesos_ens[3] * kNN]
    dt_ens[, err2 := (Real - Pred_Ensemble)^2]
    rmses["Ensemble"] <- sqrt(mean(dt_ens$err2))
    maes["Ensemble"] <- mean(abs(dt_ens$Real - dt_ens$Pred_Ensemble))
    mejor_metodo <- names(rmses)[which.min(rmses)]

    # SSE del ensemble por día, alineado con fechas_variable (para el bootstrap).
    sse_ens_por_dia <- dt_ens[, .(SSE = sum(err2)), by = indice_dia]
    sse_ens_vec <- numeric(length(fechas_variable))
    sse_ens_vec[sse_ens_por_dia$indice_dia] <- sse_ens_por_dia$SSE

    cat(sprintf(
      "  Pesos ensemble %s (1/RMSE): 1-NN=%.3f | IDW b1=%.3f | kNN=%.3f\n",
      var_name, pesos_ens[1], pesos_ens[2], pesos_ens[3]
    ))

    resultados_globales[[var_name]] <- data.table(
      Variable = var_name,
      N_Estaciones = n_estaciones,
      N_Instantes = n_instantes_validos,
      N_Predicciones = n_predicciones,
      RMSE_Media = round(rmses["Media"], 3),
      RMSE_NN = round(rmses["Vecino Cercano"], 3),
      RMSE_Voronoi = round(rmses["Closest observation"], 3),
      RMSE_kNN = round(rmses["kNN"], 3),
      RMSE_IDW_b1 = round(rmses["IDW beta=1"], 3),
      RMSE_IDW_b2 = round(rmses["IDW beta=2"], 3),
      RMSE_Ensemble = round(rmses["Ensemble"], 3),
      MAE_Media = round(maes["Media"], 3),
      MAE_NN = round(maes["Vecino Cercano"], 3),
      MAE_Voronoi = round(maes["Closest observation"], 3),
      MAE_kNN = round(maes["kNN"], 3),
      MAE_IDW_b1 = round(maes["IDW beta=1"], 3),
      MAE_IDW_b2 = round(maes["IDW beta=2"], 3),
      MAE_Ensemble = round(maes["Ensemble"], 3),
      Mejor_Metodo = mejor_metodo
    )

    detalle_variable <- data.table(
      Variable = var_name,
      FECHA = fechas_variable,
      N_Predicciones = n_predicciones_dia,
      SSE_Media = sse_dia[, "Media"],
      SSE_NN = sse_dia[, "Vecino Cercano"],
      SSE_Voronoi = sse_dia[, "Closest observation"],
      SSE_kNN = sse_dia[, "kNN"],
      SSE_IDW_b1 = sse_dia[, "IDW beta=1"],
      SSE_IDW_b2 = sse_dia[, "IDW beta=2"],
      SSE_Ensemble = sse_ens_vec
    )[N_Predicciones > 0L]
    detalles_diarios_globales[[var_name]] <- detalle_variable

    cat(sprintf(
      "%s: mejor = %s | RMSE = %.3f | %d predicciones | %.1f s\n",
      var_name,
      mejor_metodo,
      rmses[mejor_metodo],
      n_predicciones,
      as.numeric(difftime(
        Sys.time(),
        tiempo_inicio,
        units = "secs"
      ))
    ))
  }

  tabla_final <- rbindlist(
    resultados_globales,
    use.names = TRUE
  )
  print(tabla_final)
  list(
    resumen = tabla_final,
    detalle_diario = rbindlist(
      detalles_diarios_globales,
      use.names = TRUE
    )
  )
}

# Evalúa si la pequeña diferencia entre el mejor método y el segundo se
# mantiene al respetar la estructura temporal. El bootstrap remuestrea semanas
# completas, no observaciones horarias individuales.
analizar_estabilidad_temporal <- function(
    detalle_diario,
    n_bootstrap = 3000L,
    semilla = 2025L) {
  columnas_sse <- c(
    "SSE_Media",
    "SSE_NN",
    "SSE_Voronoi",
    "SSE_kNN",
    "SSE_IDW_b1",
    "SSE_IDW_b2",
    "SSE_Ensemble"
  )
  nombres_metodos <- c(
    "Media",
    "Vecino Cercano",
    "Closest observation",
    "kNN",
    "IDW beta=1",
    "IDW beta=2",
    "Ensemble"
  )
  # El ensemble es el más complejo: va al final del ranking de simplicidad para
  # que, en un empate técnico, se prefiera un método más simple.
  ranking_simplicidad <- nombres_metodos

  resumen_estabilidad <- list()

  for (indice_variable in seq_along(unique(detalle_diario$Variable))) {
    variable <- unique(detalle_diario$Variable)[indice_variable]
    dt_var <- copy(
      detalle_diario[Variable == variable]
    )

    sse_total <- colSums(
      as.matrix(dt_var[, ..columnas_sse])
    )
    n_total <- sum(dt_var$N_Predicciones)
    rmse_total <- sqrt(sse_total / n_total)
    names(rmse_total) <- nombres_metodos

    orden <- order(rmse_total)
    ganador <- names(rmse_total)[orden[1]]
    segundo <- names(rmse_total)[orden[2]]
    col_ganador <- columnas_sse[
      match(ganador, nombres_metodos)
    ]
    col_segundo <- columnas_sse[
      match(segundo, nombres_metodos)
    ]

    # Estabilidad mensual: el método ganador global puede no ganar en todos los
    # meses, incluso si su RMSE anual es el menor.
    dt_var[, MES := format(FECHA, "%Y-%m")]
    tabla_mensual <- dt_var[
      ,
      c(
        list(N_Predicciones = sum(N_Predicciones)),
        lapply(.SD, sum)
      ),
      by = MES,
      .SDcols = columnas_sse
    ]

    matriz_rmse_mensual <- sqrt(
      as.matrix(tabla_mensual[, ..columnas_sse]) /
        tabla_mensual$N_Predicciones
    )
    colnames(matriz_rmse_mensual) <- nombres_metodos
    ganador_mensual <- nombres_metodos[
      max.col(
        -matriz_rmse_mensual,
        ties.method = "first"
      )
    ]

    # Porcentaje de días en los que gana el método anual.
    matriz_rmse_diaria <- sqrt(
      as.matrix(dt_var[, ..columnas_sse]) /
        dt_var$N_Predicciones
    )
    ganador_diario <- nombres_metodos[
      max.col(
        -matriz_rmse_diaria,
        ties.method = "first"
      )
    ]

    # Bootstrap por bloques semanales. Cada bloque conserva juntas las horas y
    # los días consecutivos de esa semana.
    dt_var[, BLOQUE_SEMANAL := as.character(
      cut(FECHA, breaks = "week")
    )]
    bloques <- dt_var[
      ,
      .(
        N_Predicciones = sum(N_Predicciones),
        SSE_Ganador = sum(get(col_ganador)),
        SSE_Segundo = sum(get(col_segundo))
      ),
      by = BLOQUE_SEMANAL
    ]

    set.seed(semilla + indice_variable)
    diferencias_bootstrap <- replicate(
      n_bootstrap,
      {
        seleccion <- sample.int(
          nrow(bloques),
          nrow(bloques),
          replace = TRUE
        )
        n_b <- sum(bloques$N_Predicciones[seleccion])
        rmse_ganador_b <- sqrt(
          sum(bloques$SSE_Ganador[seleccion]) / n_b
        )
        rmse_segundo_b <- sqrt(
          sum(bloques$SSE_Segundo[seleccion]) / n_b
        )
        rmse_segundo_b - rmse_ganador_b
      }
    )

    intervalo <- quantile(
      diferencias_bootstrap,
      probs = c(0.025, 0.975),
      names = FALSE
    )
    diferencia_observada <- rmse_total[segundo] -
      rmse_total[ganador]
    decision <- if (intervalo[1] > 0) {
      "Ganador estable"
    } else {
      "Empate técnico"
    }

    candidatos_empate <- c(ganador, segundo)
    metodo_simple <- ranking_simplicidad[
      ranking_simplicidad %in% candidatos_empate
    ][1]
    recomendacion <- if (decision == "Ganador estable") {
      ganador
    } else if (variable == "Radiacion_Solar") {
      paste(
        "Sin método único:",
        "evaluar por mes o franja de radiación"
      )
    } else {
      paste0(
        metodo_simple,
        " (preferido por simplicidad en el empate)"
      )
    }

    resumen_estabilidad[[variable]] <- data.table(
      Variable = variable,
      Mejor_RMSE = ganador,
      Segundo_Metodo = segundo,
      RMSE_Mejor = rmse_total[ganador],
      RMSE_Segundo = rmse_total[segundo],
      Diferencia_RMSE = diferencia_observada,
      Diferencia_Relativa_Pct = 100 *
        diferencia_observada / rmse_total[ganador],
      IC95_Diferencia_Inferior = intervalo[1],
      IC95_Diferencia_Superior = intervalo[2],
      Probabilidad_Mejor = mean(
        diferencias_bootstrap > 0
      ),
      Meses_Ganados = sum(
        ganador_mensual == ganador
      ),
      Meses_Evaluados = length(ganador_mensual),
      Pct_Dias_Ganados = 100 * mean(
        ganador_diario == ganador
      ),
      Decision = decision,
      Metodo_Recomendado = recomendacion
    )
  }

  rbindlist(
    resumen_estabilidad,
    use.names = TRUE
  )
}

# ==============================================================================
# 4. Lectura y control de calidad
# ==============================================================================

ruta_datos <- here(
  "data",
  "processed",
  "Clima",
  "horario",
  sprintf("meteo_madrid_%d_horario_flag.rds", ANIO)
)

dt_meteo <- readRDS(ruta_datos)
setDT(dt_meteo)

# El archivo denominado 2025 contiene también fechas de 2026. Se filtra el
# intervalo explícitamente antes de cualquier validación.
fecha_inicio <- as.Date(sprintf("%d-01-01", ANIO))
fecha_fin <- as.Date(sprintf("%d-12-31", ANIO))
dt_meteo <- dt_meteo[
  FECHA >= fecha_inicio & FECHA <= fecha_fin
]

# --- Enmascarar valores IMPUTADOS: la validación usa SOLO medidas reales ------
# El fichero *_flag.rds trae columnas lógicas imp_<var> que marcan qué celdas se
# rellenaron en el preprocesamiento (interpolación lineal + vecino más cercano).
# Las ponemos a NA para que la LOOCV no valide -ni use como vecinos- valores
# imputados; así se evita el sesgo circular hacia 1-NN/IDW.
cols_flag <- grep("^imp_", names(dt_meteo), value = TRUE)
n_enmascarados <- 0L
for (columna_flag in cols_flag) {
  variable_valor <- sub("^imp_", "", columna_flag)
  if (variable_valor %in% names(dt_meteo)) {
    filas_imputadas <- which(dt_meteo[[columna_flag]])
    if (length(filas_imputadas) > 0L) {
      set(dt_meteo, i = filas_imputadas, j = variable_valor, value = NA_real_)
      n_enmascarados <- n_enmascarados + length(filas_imputadas)
    }
  }
}
if (length(cols_flag) > 0L) dt_meteo[, (cols_flag) := NULL]
cat(sprintf(
  "Valores imputados enmascarados (excluidos de la validación): %d\n",
  n_enmascarados
))

# Convertir H01,...,H24 a 1,...,24.
dt_meteo[, HORA := as.integer(
  gsub("[^0-9]", "", as.character(HORA))
)]

if (anyNA(dt_meteo$HORA) ||
    any(!dt_meteo$HORA %between% c(1L, 24L))) {
  stop("Se encontraron horas fuera del intervalo H01-H24.")
}

# Normalizar los nombres con acentos para evitar dependencias del locale.
col_presion <- grep(
  "^Presion Barom",
  names(dt_meteo),
  value = TRUE
)
col_radiacion <- grep(
  "Solar$",
  names(dt_meteo),
  value = TRUE
)
if (length(col_presion) == 1L) {
  setnames(
    dt_meteo,
    col_presion,
    "Presion_Barometrica"
  )
}
if (length(col_radiacion) == 1L) {
  setnames(
    dt_meteo,
    col_radiacion,
    "Radiacion_Solar"
  )
}

normalizar_coordenadas_madrid(dt_meteo)

columnas_faltantes <- setdiff(
  c(
    "ESTACION",
    "LONGITUD",
    "LATITUD",
    "FECHA",
    "HORA",
    VARIABLES_CLIMA
  ),
  names(dt_meteo)
)
if (length(columnas_faltantes) > 0L) {
  stop(
    "Faltan columnas necesarias: ",
    paste(columnas_faltantes, collapse = ", ")
  )
}

duplicados <- dt_meteo[
  ,
  .N,
  by = .(ESTACION, FECHA, HORA)
][N > 1L]
if (nrow(duplicados) > 0L) {
  stop(
    "Hay filas duplicadas para una misma estación, fecha y hora."
  )
}

ubicaciones_estaciones <- unique(
  dt_meteo[, .(ESTACION, LONGITUD, LATITUD)]
)
ubicaciones_coincidentes <- ubicaciones_estaciones[
  ,
  .(
    N_Estaciones = .N,
    Estaciones = paste(sort(ESTACION), collapse = " | ")
  ),
  by = .(LONGITUD, LATITUD)
][N_Estaciones > 1L]

cat(sprintf(
  "Datos horarios: %d filas | %d estaciones | %d instantes.\n",
  nrow(dt_meteo),
  uniqueN(dt_meteo$ESTACION),
  uniqueN(dt_meteo[, .(FECHA, HORA)])
))

if (nrow(ubicaciones_coincidentes) > 0L) {
  cat(
    "Advertencia:",
    nrow(ubicaciones_coincidentes),
    "ubicaciones contienen más de una estación.\n"
  )
}

# ==============================================================================
# 5. Seleccion de horas relevantes de precipitacion, radiacion y viento
# ==============================================================================

resumen_instantes <- dt_meteo[, {
  precip_valida <- Precipitaciones[
    !is.na(Precipitaciones)
  ]
  radiacion_valida <- Radiacion_Solar[
    !is.na(Radiacion_Solar)
  ]
  viento_valido <- get("Velocidad Viento")[
    !is.na(get("Velocidad Viento"))
  ]

  list(
    N_Estaciones_Precip = length(precip_valida),
    N_Estaciones_Lluvia = sum(
      precip_valida > UMBRAL_LLUVIA_MM
    ),
    Precipitacion_Media = if (length(precip_valida)) {
      mean(precip_valida)
    } else {
      NA_real_
    },
    N_Estaciones_Radiacion = length(radiacion_valida),
    Radiacion_Mediana = if (length(radiacion_valida)) {
      median(radiacion_valida)
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
}, by = .(FECHA, HORA)]

instantes_lluvia <- resumen_instantes[
  N_Estaciones_Precip >= MIN_ESTACIONES_PRECIP &
    N_Estaciones_Lluvia >= MIN_ESTACIONES_CON_LLUVIA
]

if (nrow(instantes_lluvia) == 0L) {
  stop(
    "No se encontraron horas de lluvia con el criterio configurado."
  )
}

instantes_radiacion <- resumen_instantes[
  N_Estaciones_Radiacion >= MIN_ESTACIONES_RADIACION &
    Radiacion_Mediana > UMBRAL_RADIACION_WM2
]

if (nrow(instantes_radiacion) == 0L) {
  stop(
    "No se encontraron horas de radiación solar informativa."
  )
}

instantes_viento_disponibles <- resumen_instantes[
  N_Estaciones_Viento >= MIN_ESTACIONES_VIENTO &
    !is.na(Viento_Mediana)
]

umbral_viento_p75 <- quantile(
  instantes_viento_disponibles$Viento_Mediana,
  probs = PERCENTIL_VIENTO_VALIDACION,
  names = FALSE
)

instantes_viento <- instantes_viento_disponibles[
  Viento_Mediana >= umbral_viento_p75
]

if (nrow(instantes_viento) == 0L) {
  stop(
    "No se encontraron horas de viento con el criterio configurado."
  )
}

cat(sprintf(
  paste0(
    "Horas de lluvia: %d | Horas con radiación > %.1f W/m²: %d | ",
    "Horas de viento >= P75: %d.\n"
  ),
  nrow(instantes_lluvia),
  UMBRAL_RADIACION_WM2,
  nrow(instantes_radiacion),
  nrow(instantes_viento)
))

# ==============================================================================
# 6. Validacion horaria
# ==============================================================================

# Resultado de todas las horas: se conserva como referencia completa y se usa
# para temperatura, humedad, presión y radiación.
resultado_anual <- comparar_interpolaciones_loocv_horaria(
  dt_meteo = dt_meteo,
  variables = VARIABLES_CLIMA,
  k_vecinos = K_VECINOS,
  betas_idw = BETAS_IDW
)
tabla_anual <- resultado_anual$resumen

claves_lluvia <- instantes_lluvia[, .(FECHA, HORA)]
dt_lluvia <- dt_meteo[
  claves_lluvia,
  on = .(FECHA, HORA),
  nomatch = 0L
]

resultado_precipitacion <- comparar_interpolaciones_loocv_horaria(
  dt_meteo = dt_lluvia,
  variables = "Precipitaciones",
  k_vecinos = K_VECINOS,
  betas_idw = BETAS_IDW
)
tabla_precipitacion <- resultado_precipitacion$resumen

claves_radiacion <- instantes_radiacion[, .(FECHA, HORA)]
dt_radiacion <- dt_meteo[
  claves_radiacion,
  on = .(FECHA, HORA),
  nomatch = 0L
]

resultado_radiacion <- comparar_interpolaciones_loocv_horaria(
  dt_meteo = dt_radiacion,
  variables = "Radiacion_Solar",
  k_vecinos = K_VECINOS,
  betas_idw = BETAS_IDW
)
tabla_radiacion <- resultado_radiacion$resumen

claves_viento <- instantes_viento[, .(FECHA, HORA)]
dt_viento <- dt_meteo[
  claves_viento,
  on = .(FECHA, HORA),
  nomatch = 0L
]

resultado_viento <- comparar_interpolaciones_loocv_horaria(
  dt_meteo = dt_viento,
  variables = "Velocidad Viento",
  k_vecinos = K_VECINOS,
  betas_idw = BETAS_IDW
)
tabla_viento <- resultado_viento$resumen

tabla_decision <- rbindlist(
  list(
    tabla_anual[
      !Variable %chin% c(
        "Precipitaciones",
        "Radiacion_Solar",
        "Velocidad Viento"
      )
    ],
    tabla_precipitacion,
    tabla_radiacion,
    tabla_viento
  ),
  use.names = TRUE
)
tabla_decision <- tabla_decision[
  match(VARIABLES_CLIMA, Variable)
]
tabla_decision[, Muestra := fifelse(
  Variable == "Precipitaciones",
  "Horas de lluvia",
  fifelse(
    Variable == "Radiacion_Solar",
    "Horas con radiación > 10 W/m²",
    fifelse(
      Variable == "Velocidad Viento",
      "Horas de viento >= P75",
      "Todas las horas de 2025"
    )
  )
)]
setcolorder(
  tabla_decision,
  c(
    "Variable",
    "Muestra",
    "N_Estaciones",
    "N_Instantes",
    "N_Predicciones",
    grep(
      "^RMSE_",
      names(tabla_decision),
      value = TRUE
    ),
    grep(
      "^MAE_",
      names(tabla_decision),
      value = TRUE
    ),
    "Mejor_Metodo"
  )
)

cat("\nTabla utilizada para decidir el método horario:\n")
print(tabla_decision)

detalle_decision <- rbindlist(
  list(
    resultado_anual$detalle_diario[
      !Variable %chin% c(
        "Precipitaciones",
        "Radiacion_Solar",
        "Velocidad Viento"
      )
    ],
    resultado_precipitacion$detalle_diario,
    resultado_radiacion$detalle_diario,
    resultado_viento$detalle_diario
  ),
  use.names = TRUE
)

tabla_estabilidad <- analizar_estabilidad_temporal(
  detalle_diario = detalle_decision,
  n_bootstrap = N_BOOTSTRAP_BLOQUES,
  semilla = SEMILLA_BOOTSTRAP
)[
  match(VARIABLES_CLIMA, Variable)
]

cat("\nEstabilidad temporal y decisión final:\n")
print(tabla_estabilidad)

# ==============================================================================
# 6b. Validacion de UN UNICO DIA por variable (todas sus horas)
# ==============================================================================
# Las tablas anteriores promedian el RMSE sobre miles de instantes (todo 2025,
# u horas de lluvia/radiacion/viento). Aqui evaluamos la LOOCV horaria en UN
# SOLO DIA por variable (todas sus horas) para ver la capacidad predictiva sin
# promediar entre fechas. Dia de referencia para las variables continuas, el
# dia con mas horas de lluvia para precipitacion y el dia con mas horas de
# viento (>= P75) para la velocidad del viento.

FECHA_REFERENCIA <- as.Date("2025-04-03")
dia_lluvia <- instantes_lluvia[, .N, by = FECHA][which.max(N), FECHA]
dia_viento <- instantes_viento[, .N, by = FECHA][which.max(N), FECHA]

dias_un_dia_h <- data.table(
  Variable = VARIABLES_CLIMA,
  FECHA = as.Date(c(
    FECHA_REFERENCIA,  # Temperatura
    FECHA_REFERENCIA,  # Humedad_Relativa
    dia_lluvia,        # Precipitaciones (dia con mas horas de lluvia)
    FECHA_REFERENCIA,  # Presion_Barometrica
    FECHA_REFERENCIA,  # Radiacion_Solar
    dia_viento         # Velocidad Viento (dia con mas horas de viento)
  ))
)

cat(sprintf(
  "\nDias del analisis de un dia: referencia %s | lluvia %s | viento %s.\n",
  FECHA_REFERENCIA, dia_lluvia, dia_viento
))

tabla_un_dia_horaria <- rbindlist(
  lapply(seq_len(nrow(dias_un_dia_h)), function(i) {
    variable <- dias_un_dia_h$Variable[i]
    dia <- dias_un_dia_h$FECHA[i]
    res <- comparar_interpolaciones_loocv_horaria(
      dt_meteo = dt_meteo[FECHA == dia],
      variables = variable,
      k_vecinos = K_VECINOS,
      betas_idw = BETAS_IDW
    )$resumen
    res[, Dia := dia]
    res
  }),
  use.names = TRUE
)
tabla_un_dia_horaria <- tabla_un_dia_horaria[
  match(VARIABLES_CLIMA, Variable)
]
setcolorder(
  tabla_un_dia_horaria,
  c(
    "Variable", "Dia", "N_Estaciones", "N_Instantes", "N_Predicciones",
    grep("^RMSE_", names(tabla_un_dia_horaria), value = TRUE),
    grep("^MAE_", names(tabla_un_dia_horaria), value = TRUE),
    "Mejor_Metodo"
  )
)

cat("\nValidacion horaria de un unico dia por variable (sin promediar):\n")
print(tabla_un_dia_horaria)

# ==============================================================================
# 7. Salida PNG
# ==============================================================================

DIR_SALIDA <- crear_directorio_versionado(
  DIR_SALIDA_BASE
)
cat("\nDirectorio de esta ejecución:\n", DIR_SALIDA, "\n")

subtitulo_tabla <- paste0(
  "Escala horaria | k=", K_VECINOS,
  " | IDW beta=1 y 2 | precipitación: horas de lluvia",
  " | radiación: horas > 10 W/m² | viento: horas >= P75"
)
crear_tabla_png(
  tabla_decision,
  file.path(
    DIR_SALIDA,
    "tabla_comparacion_metodos.png"
  ),
  subtitulo_tabla
)

fwrite(
  tabla_decision,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos.csv")
)
fwrite(
  tabla_un_dia_horaria,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos_un_dia.csv")
)
crear_tabla_png(
  tabla_un_dia_horaria,
  file.path(DIR_SALIDA, "tabla_comparacion_metodos_un_dia.png"),
  paste0(
    "Escala horaria | UN DIA por variable | referencia ",
    FECHA_REFERENCIA, " | lluvia ", dia_lluvia,
    " | viento ", dia_viento
  )
)

cat("\nEjecución horaria terminada correctamente.\n")
cat(
  "Tabla PNG guardada en:\n",
  file.path(DIR_SALIDA, "tabla_comparacion_metodos.png"),
  "\n"
)
