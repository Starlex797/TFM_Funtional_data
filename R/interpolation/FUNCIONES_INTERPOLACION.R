# ==============================================================================
# FUNCIONES DE INTERPOLACIÓN ESPACIAL
# ==============================================================================
# Funciones para interpolar variables climáticas mediante IDW, kNN y Vecino
# Más Cercano desde las estaciones meteorológicas a ubicaciones objetivo.
# Validación mediante Leave-One-Out Cross-Validation (LOOCV).
# ==============================================================================

library(sf)
library(gstat)
library(data.table)
library(here)

#' Interpolación IDW diaria de variables climáticas a ubicaciones objetivo
#'
#' Para cada fecha presente en dt_meteo, interpola espacialmente las variables
#' climáticas desde las estaciones meteorológicas con datos válidos ese día
#' a las ubicaciones objetivo (e.g. estaciones de NO2) usando IDW.
#'
#' @param dt_meteo data.table con columnas ESTACION, LONGITUD, LATITUD, FECHA
#'   y las variables climáticas a interpolar.
#' @param dt_objetivo data.table con columnas ESTACION, LONGITUD, LATITUD de
#'   las ubicaciones donde se quieren estimar los valores.
#' @param variables character vector con los nombres de las columnas climáticas.
#' @param idp numeric. Potencia de la distancia inversa (por defecto 2).
#' @param min_estaciones integer. Número mínimo de estaciones requerido.
#' @param crs_orig integer. CRS de entrada (por defecto 4326 WGS84).
#' @param crs_proj integer. CRS proyectado (por defecto 25830 UTM 30N).
#'
#' @return data.table con columnas ESTACION, FECHA y variables interpoladas.
interpolar_idw_clima <- function(dt_meteo,
                                 dt_objetivo,
                                 variables = c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                                               "Presion Barométrica", "Radiación Solar", "Velocidad Viento"),
                                 idp       = 2, # Potencia de la distancia inversa
                                 min_estaciones = 7,# Número mínimo de estaciones requeridas para interpolar
                                 crs_orig  = 4326,
                                 crs_proj  = 25830) {

  vars_disponibles <- intersect(variables, names(dt_meteo))
  if (length(vars_disponibles) == 0) {
    stop("Ninguna de las variables solicitadas existe en dt_meteo: ", paste(variables, collapse = ", "))
  }

  coords_objetivo <- unique(dt_objetivo[, .(ESTACION, LONGITUD, LATITUD)])
  sf_objetivo <- st_as_sf(coords_objetivo,
                          coords = c("LONGITUD", "LATITUD"),
                          crs = crs_orig) |>
    st_transform(crs = crs_proj)

  coords_meteo <- unique(dt_meteo[, .(ESTACION, LONGITUD, LATITUD)])
  sf_meteo_base <- st_as_sf(coords_meteo,
                            coords = c("LONGITUD", "LATITUD"),
                            crs = crs_orig) |>
    st_transform(crs = crs_proj)

  if ("HORA" %in% names(dt_meteo)) {
    instantes <- unique(dt_meteo[, .(FECHA, HORA)])
    setorder(instantes, FECHA, HORA)
    llaves_tiempo <- c("FECHA", "HORA")
  } else {
    instantes <- unique(dt_meteo[, .(FECHA)])
    setorder(instantes, FECHA)
    llaves_tiempo <- "FECHA"
  }

  n_instantes <- nrow(instantes)
  resultados <- vector("list", n_instantes)

  for (i in seq_len(n_instantes)) {
    instante_actual <- instantes[i]
    dt_momento <- dt_meteo[instante_actual, on = llaves_tiempo]

    dt_out <- data.table(ESTACION = coords_objetivo$ESTACION)
    for (col in llaves_tiempo) {
      dt_out[, (col) := instante_actual[[col]]]
    }

    for (var_name in vars_disponibles) {
      dt_var <- merge(dt_momento[, .(ESTACION, valor = get(var_name))],
                      coords_meteo, by = "ESTACION")
      dt_var <- dt_var[!is.na(valor)]

      if (nrow(dt_var) < min_estaciones) {
        dt_out[, (var_name) := NA_real_]
        next
      }

      sf_obs <- sf_meteo_base[sf_meteo_base$ESTACION %in% dt_var$ESTACION, ]
      sf_obs[["IDW_VAR"]] <- dt_var$valor[match(sf_obs$ESTACION, dt_var$ESTACION)]

      idw_result <- idw(IDW_VAR ~ 1,
                        locations = sf_obs,
                        newdata   = sf_objetivo,
                        idp       = idp,
                        debug.level = 0)

      dt_out[, (var_name) := idw_result$var1.pred]
    }

    resultados[[i]] <- dt_out

    if (i %% 50 == 0 || i == n_instantes) {
      message(sprintf("Interpolación IDW: %d / %d instantes completados", i, n_instantes))
    }
  }

  dt_final <- rbindlist(resultados)

  for (var_name in vars_disponibles) {
    n_ok <- sum(!is.na(dt_final[[var_name]]))
    n_total <- nrow(dt_final)
    message(sprintf("'%s': %d / %d valores interpolados (%.1f%%)",
                    var_name, n_ok, n_total, 100 * n_ok / n_total))
  }

  return(dt_final)
}


# ==============================================================================
# VALIDACIÓN LEAVE-ONE-OUT CROSS-VALIDATION (LOOCV)
# ==============================================================================
#' Compara métodos de interpolación usando LOOCV
#'
#' Para cada estación, la oculta y predice su valor con las restantes,
#' repitiendo para cada día y cada variable. Incluye: Media, Vecino Más
#' Cercano (1-NN), Closest observation mediante Voronoi, kNN (k vecinos,
#' promedio no ponderado), IDW y Ensemble. Closest observation y 1-NN son
#' numéricamente equivalentes; se conservan ambos nombres para documentar
#' las dos formulaciones del mismo estimador. El Ensemble combina los tres
#' interpoladores locales (closest observation / 1-NN, IDW beta=1 y kNN) con
#' pesos proporcionales a 1/RMSE de cada método, de modo que los más precisos
#' pesan más. No incluye la Media por no ser un interpolador espacial.
#'
#' @param dt_meteo data.table con ESTACION, LONGITUD, LATITUD, FECHA y variables.
#' @param variables character vector de variables a evaluar.
#' @param k_vecinos integer. Número de vecinos para kNN e IDW (por defecto 3).
#' @param betas_idw numeric. Potencias de IDW a comparar (beta = 1 y 2).
#'
#' @return data.table con RMSE y MAE por variable y método.
comparar_interpolaciones_loocv <- function(
    dt_meteo, variables, k_vecinos = 3L, betas_idw = c(1, 2)) {

  library(cli)

  if (k_vecinos < 1L) stop("k_vecinos debe ser al menos 1.")
  betas_idw <- sort(unique(as.numeric(betas_idw)))
  if (!identical(betas_idw, c(1, 2))) {
    stop("betas_idw debe contener exactamente los valores 1 y 2.")
  }

  rmse <- function(real, pred) sqrt(mean((real - pred)^2, na.rm = TRUE))
  mae <- function(real, pred) mean(abs(real - pred), na.rm = TRUE)
  resultados_globales <- list()

  cli_h1("LOOCV diaria - Validación Leave-One-Out")
  cli_alert_info(
    "Variables: {length(variables)} | kNN e IDW: k = {k_vecinos} | IDW: beta = 1 y 2"
  )

  for (var_name in variables) {
    dt_var <- copy(dt_meteo)[
      !is.na(get(var_name)),
      .(
        FECHA,
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
      , .(n = .N, LON = first(LONGITUD), LAT = first(LATITUD)),
      by = ESTACION
    ]
    setorder(conteo_estaciones, LON, LAT, -n)
    estaciones_conservadas <- conteo_estaciones[
      !duplicated(paste(LON, LAT)), ESTACION
    ]
    if (length(estaciones_conservadas) < uniqueN(dt_var$ESTACION)) {
      dt_var <- dt_var[ESTACION %in% estaciones_conservadas]
    }

    n_est <- uniqueN(dt_var$ESTACION)
    if (n_est <= k_vecinos) {
      cli_alert_warning(
        "Omitiendo '{var_name}': solo {n_est} estaciones (se necesitan al menos {k_vecinos + 1L})."
      )
      next
    }

    # Las coordenadas son fijas por estación. Se transforman una sola vez a
    # UTM 25830 para que todas las distancias estén expresadas en metros.
    dt_coords <- dt_var[
      , .(LONGITUD = first(LONGITUD), LATITUD = first(LATITUD)),
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
    dt_var <- dt_coords[dt_var, on = "ESTACION"]

    fechas <- sort(unique(dt_var$FECHA))
    resultados_dia <- vector("list", length(fechas))

    pb <- cli_progress_bar(
      name = var_name,
      total = length(fechas),
      format = "{cli::pb_name} | {cli::pb_bar} {cli::pb_percent} | Día {cli::pb_current}/{cli::pb_total}"
    )

    for (j in seq_along(fechas)) {
      cli_progress_update(id = pb)
      dt_dia <- dt_var[FECHA == fechas[j]]
      n_dia <- nrow(dt_dia)

      # Al ocultar una estación deben quedar exactamente k vecinos disponibles.
      if (n_dia <= k_vecinos) next

      valores <- dt_dia$VALOR
      matriz_xy <- as.matrix(dt_dia[, .(X, Y)])
      distancias <- as.matrix(dist(matriz_xy))
      diag(distancias) <- Inf

      pred_media <- (sum(valores) - valores) / (n_dia - 1L)
      pred_nn <- numeric(n_dia)
      pred_knn <- numeric(n_dia)
      pred_idw_b1 <- numeric(n_dia)
      pred_idw_b2 <- numeric(n_dia)

      for (i in seq_len(n_dia)) {
        idx_vecinos <- order(distancias[i, ])[seq_len(k_vecinos)]
        valores_vecinos <- valores[idx_vecinos]
        distancias_vecinos <- distancias[i, idx_vecinos]

        pred_nn[i] <- valores_vecinos[1]
        pred_knn[i] <- mean(valores_vecinos)

        pesos_b1 <- 1 / pmax(distancias_vecinos, 1)^betas_idw[1]
        pesos_b2 <- 1 / pmax(distancias_vecinos, 1)^betas_idw[2]
        pred_idw_b1[i] <- sum(pesos_b1 * valores_vecinos) / sum(pesos_b1)
        pred_idw_b2[i] <- sum(pesos_b2 * valores_vecinos) / sum(pesos_b2)
      }

      resultados_dia[[j]] <- data.table(
        FECHA = dt_dia$FECHA,
        ESTACION = dt_dia$ESTACION,
        Real = valores,
        Pred_Media = pred_media,
        Pred_NN = pred_nn,
        Pred_Voronoi = pred_nn,
        Pred_kNN = pred_knn,
        Pred_IDW_b1 = pred_idw_b1,
        Pred_IDW_b2 = pred_idw_b2
      )
    }

    cli_progress_done(id = pb)
    dt_res <- rbindlist(resultados_dia, use.names = TRUE)

    if (nrow(dt_res) == 0L) {
      cli_alert_warning("Omitiendo '{var_name}': no hay días con suficientes estaciones.")
      next
    }

    # Ensemble: combinación ponderada de los tres interpoladores locales
    # (closest observation / 1-NN, IDW beta=1 y kNN) con pesos proporcionales a
    # 1/RMSE de cada método (los más precisos pesan más; el peor queda casi
    # anulado). No incluye la Media porque no es un interpolador espacial.
    rmse_ens_nn <- rmse(dt_res$Real, dt_res$Pred_NN)
    rmse_ens_idw <- rmse(dt_res$Real, dt_res$Pred_IDW_b1)
    rmse_ens_knn <- rmse(dt_res$Real, dt_res$Pred_kNN)
    pesos_ens <- 1 / pmax(c(rmse_ens_nn, rmse_ens_idw, rmse_ens_knn), 1e-9)
    pesos_ens <- pesos_ens / sum(pesos_ens)
    dt_res[, Pred_Ensemble :=
      pesos_ens[1] * Pred_NN +
        pesos_ens[2] * Pred_IDW_b1 +
        pesos_ens[3] * Pred_kNN]
    cli_alert_info(
      "Pesos ensemble {var_name} (1/RMSE): 1-NN={round(pesos_ens[1], 3)} | IDW b1={round(pesos_ens[2], 3)} | kNN={round(pesos_ens[3], 3)}"
    )

    rmses <- c(
      "Media" = rmse(dt_res$Real, dt_res$Pred_Media),
      "Vecino Cercano" = rmse(dt_res$Real, dt_res$Pred_NN),
      "Closest observation" = rmse(dt_res$Real, dt_res$Pred_Voronoi),
      "kNN" = rmse(dt_res$Real, dt_res$Pred_kNN),
      "IDW beta=1" = rmse(dt_res$Real, dt_res$Pred_IDW_b1),
      "IDW beta=2" = rmse(dt_res$Real, dt_res$Pred_IDW_b2),
      "Ensemble" = rmse(dt_res$Real, dt_res$Pred_Ensemble)
    )
    maes <- c(
      "Media" = mae(dt_res$Real, dt_res$Pred_Media),
      "Vecino Cercano" = mae(dt_res$Real, dt_res$Pred_NN),
      "Closest observation" = mae(dt_res$Real, dt_res$Pred_Voronoi),
      "kNN" = mae(dt_res$Real, dt_res$Pred_kNN),
      "IDW beta=1" = mae(dt_res$Real, dt_res$Pred_IDW_b1),
      "IDW beta=2" = mae(dt_res$Real, dt_res$Pred_IDW_b2),
      "Ensemble" = mae(dt_res$Real, dt_res$Pred_Ensemble)
    )
    mejor_metodo <- names(rmses)[which.min(rmses)]

    resultados_globales[[var_name]] <- data.table(
      Variable = var_name,
      N_Estaciones = n_est,
      N_Dias = uniqueN(dt_res$FECHA),
      N_Predicciones = nrow(dt_res),
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

    cli_alert_success(
      "{var_name}: mejor = {mejor_metodo} (RMSE {round(rmses[mejor_metodo], 3)}) | {nrow(dt_res)} predicciones"
    )
  }

  tabla_final <- rbindlist(resultados_globales, use.names = TRUE)
  cli_h2("Resultados finales - LOOCV diaria")
  print(tabla_final)

  return(tabla_final)
}


# ==============================================================================
# INTERPOLACIÓN CON EL MÉTODO GANADOR POR VARIABLE
# ==============================================================================

#' Pesos del ensemble (1/RMSE) para una variable mediante LOOCV
#'
#' Calcula, sobre dt_meteo (diaria o horaria), el RMSE agregado por LOOCV de los
#' tres interpoladores locales (closest observation / 1-NN, IDW beta=1 y kNN) y
#' devuelve los pesos proporcionales a 1/RMSE, normalizados a suma 1. El orden
#' del vector devuelto es (1-NN, IDW beta=1, kNN).
#'
#' @param dt_meteo data.table con ESTACION, LONGITUD, LATITUD, FECHA (y HORA si
#'   es horaria) y la variable.
#' @param variable character. Nombre de la columna a evaluar.
#' @param k_vecinos integer. Vecinos para kNN e IDW (por defecto 3).
pesos_ensemble_loocv <- function(dt_meteo, variable, k_vecinos = 3L,
                                 crs_orig = 4326, crs_proj = 25830) {
  tiene_hora <- "HORA" %in% names(dt_meteo)
  llaves <- if (tiene_hora) c("FECHA", "HORA") else "FECHA"

  dt <- copy(dt_meteo)[!is.na(get(variable))]
  dt[, VALOR := as.numeric(get(variable))]

  # Distancias en metros entre estaciones (coordenadas fijas por estación).
  dt_coords <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
  sf_c <- st_transform(
    st_as_sf(dt_coords, coords = c("LONGITUD", "LATITUD"), crs = crs_orig),
    crs_proj
  )
  xy <- st_coordinates(sf_c)
  dt_coords[, `:=`(X = xy[, 1], Y = xy[, 2])]
  d_mat <- as.matrix(dist(as.matrix(dt_coords[, .(X, Y)])))
  rownames(d_mat) <- dt_coords$ESTACION
  colnames(d_mat) <- dt_coords$ESTACION

  setorderv(dt, c(llaves, "ESTACION"))
  grupos <- dt[, .(INI = first(.I), FIN = last(.I)), by = llaves]

  sse_nn <- 0
  sse_idw <- 0
  sse_knn <- 0
  n_tot <- 0L

  for (g in seq_len(nrow(grupos))) {
    idx <- grupos$INI[g]:grupos$FIN[g]
    est <- dt$ESTACION[idx]
    val <- dt$VALOR[idx]
    n <- length(val)
    if (n <= k_vecinos) next

    d_loc <- d_mat[est, est, drop = FALSE]
    diag(d_loc) <- Inf
    for (i in seq_len(n)) {
      ord <- order(d_loc[i, ])[seq_len(k_vecinos)]
      vv <- val[ord]
      dd <- d_loc[i, ord]
      p_nn <- vv[1]
      p_knn <- mean(vv)
      w1 <- 1 / pmax(dd, 1)^1
      p_idw <- sum(w1 * vv) / sum(w1)
      sse_nn <- sse_nn + (val[i] - p_nn)^2
      sse_idw <- sse_idw + (val[i] - p_idw)^2
      sse_knn <- sse_knn + (val[i] - p_knn)^2
    }
    n_tot <- n_tot + n
  }

  if (n_tot == 0L) {
    stop("No hay instantes con estaciones suficientes para '", variable, "'.")
  }
  rmses <- sqrt(c(sse_nn, sse_idw, sse_knn) / n_tot)
  pesos <- 1 / pmax(rmses, 1e-9)
  as.numeric(pesos / sum(pesos))
}


#' Interpola las variables climáticas a ubicaciones objetivo usando, para cada
#' variable, su método ganador.
#'
#' Para cada instante (día, o día+hora si hay HORA) interpola desde las
#' estaciones meteorológicas con datos válidos a las ubicaciones objetivo
#' (e.g. estaciones de NO2). El método se elige por variable mediante el mapa
#' metodos_por_variable. Métodos soportados: "IDW beta=1", "IDW beta=2", "kNN",
#' "Vecino Cercano"/"Closest observation" (1-NN) y "Ensemble" (combinación
#' ponderada 1/RMSE de 1-NN, IDW beta=1 y kNN).
#'
#' @param dt_meteo data.table con ESTACION, LONGITUD, LATITUD, FECHA (y HORA) y
#'   las variables.
#' @param dt_objetivo data.table con ESTACION, LONGITUD, LATITUD objetivo.
#' @param metodos_por_variable named character. variable -> método.
#' @param k_vecinos integer. Vecinos para kNN e IDW (por defecto 3).
#' @param min_estaciones integer. Estaciones mínimas para interpolar un instante.
#' @param pesos_ensemble named list opcional. variable -> vector de pesos
#'   (1-NN, IDW beta=1, kNN). Si falta para una variable Ensemble, se estima por
#'   LOOCV con pesos_ensemble_loocv().
#'
#' @return data.table con ESTACION, claves temporales y variables interpoladas.
interpolar_clima_por_metodo <- function(dt_meteo, dt_objetivo,
                                        metodos_por_variable,
                                        k_vecinos = 3L,
                                        min_estaciones = 7L,
                                        crs_orig = 4326, crs_proj = 25830,
                                        pesos_ensemble = NULL) {
  vars <- intersect(names(metodos_por_variable), names(dt_meteo))
  if (length(vars) == 0L) {
    stop("Ninguna variable de 'metodos_por_variable' existe en dt_meteo.")
  }

  coords_obj <- unique(dt_objetivo[, .(ESTACION, LONGITUD, LATITUD)])
  sf_obj <- st_transform(
    st_as_sf(coords_obj, coords = c("LONGITUD", "LATITUD"), crs = crs_orig),
    crs_proj
  )

  coords_meteo <- unique(dt_meteo[, .(ESTACION, LONGITUD, LATITUD)])
  sf_meteo_base <- st_transform(
    st_as_sf(coords_meteo, coords = c("LONGITUD", "LATITUD"), crs = crs_orig),
    crs_proj
  )

  # Pesos del ensemble (solo para variables cuyo método es Ensemble).
  if (is.null(pesos_ensemble)) pesos_ensemble <- list()
  for (v in vars) {
    if (metodos_por_variable[[v]] == "Ensemble" && is.null(pesos_ensemble[[v]])) {
      pesos_ensemble[[v]] <- pesos_ensemble_loocv(
        dt_meteo, v, k_vecinos, crs_orig, crs_proj
      )
      message(sprintf(
        "Pesos ensemble %s (1/RMSE): 1-NN=%.3f | IDW b1=%.3f | kNN=%.3f",
        v, pesos_ensemble[[v]][1], pesos_ensemble[[v]][2], pesos_ensemble[[v]][3]
      ))
    }
  }

  tiene_hora <- "HORA" %in% names(dt_meteo)
  llaves <- if (tiene_hora) c("FECHA", "HORA") else "FECHA"
  instantes <- unique(dt_meteo[, ..llaves])
  setorderv(instantes, llaves)
  n_inst <- nrow(instantes)
  resultados <- vector("list", n_inst)

  # Predice una variable en las ubicaciones objetivo según el método.
  predecir <- function(sf_obs, metodo, pesos = NULL) {
    k_real <- min(k_vecinos, nrow(sf_obs))
    f <- IDW_VAR ~ 1
    p <- function(nmax, idp) {
      idw(f, sf_obs, sf_obj, nmax = nmax, idp = idp, debug.level = 0)$var1.pred
    }
    switch(metodo,
      "IDW beta=1" = p(k_real, 1),
      "IDW beta=2" = p(k_real, 2),
      "kNN" = p(k_real, 0),
      "Vecino Cercano" = p(1, 2),
      "Closest observation" = p(1, 2),
      "Ensemble" = pesos[1] * p(1, 2) +
        pesos[2] * p(k_real, 1) +
        pesos[3] * p(k_real, 0),
      stop("Método no reconocido: ", metodo)
    )
  }

  for (t in seq_len(n_inst)) {
    inst <- instantes[t]
    dt_mom <- dt_meteo[inst, on = llaves]

    dt_out <- data.table(ESTACION = coords_obj$ESTACION)
    for (col in llaves) dt_out[, (col) := inst[[col]]]

    for (v in vars) {
      dv <- merge(
        dt_mom[, .(ESTACION, valor = get(v))], coords_meteo, by = "ESTACION"
      )
      dv <- dv[!is.na(valor)]
      if (nrow(dv) < min_estaciones) {
        dt_out[, (v) := NA_real_]
        next
      }
      sf_obs <- sf_meteo_base[sf_meteo_base$ESTACION %in% dv$ESTACION, ]
      sf_obs[["IDW_VAR"]] <- dv$valor[match(sf_obs$ESTACION, dv$ESTACION)]
      dt_out[, (v) := predecir(
        sf_obs, metodos_por_variable[[v]], pesos_ensemble[[v]]
      )]
    }

    resultados[[t]] <- dt_out
    if (t %% 200 == 0 || t == n_inst) {
      message(sprintf("Interpolación por método: %d / %d instantes", t, n_inst))
    }
  }

  rbindlist(resultados)
}
