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
#' Cercano (1-NN), kNN (k vecinos, promedio no ponderado) e IDW.
#'
#' @param dt_meteo data.table con ESTACION, LONGITUD, LATITUD, FECHA y variables.
#' @param variables character vector de variables a evaluar.
#' @param k_vecinos integer. Número de vecinos para kNN (por defecto 5).
#'
#' @return data.table con RMSE y MAE por variable y método.
comparar_interpolaciones_loocv <- function(dt_meteo, variables, k_vecinos = 5) {

  library(cli)
  options(gstat.messages = FALSE)
  on.exit(options(gstat.messages = TRUE))

  resultados_globales <- list()
  f_formula <- var_value ~ 1

  cli_h1("LOOCV - Validación Leave-One-Out")
  cli_alert_info("Variables a evaluar: {length(variables)} | k vecinos = {k_vecinos}")

  for (var_name in variables) {

    dt_var <- copy(dt_meteo)[!is.na(get(var_name))]

    estaciones_disponibles <- unique(dt_var$ESTACION)
    n_est <- length(estaciones_disponibles)

    if (n_est < 4) {
      cli_alert_warning("Omitiendo '{var_name}': solo {n_est} estaciones (insuficiente).")
      next
    }

    res_var <- list()

    pb <- cli_progress_bar(
      name   = var_name,
      total  = n_est,
      format = "{cli::pb_name} | {cli::pb_bar} {cli::pb_percent} | Est. {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
    )

    for (est_oculta in estaciones_disponibles) {
      cli_progress_update(id = pb)

      dt_train <- dt_var[ESTACION != est_oculta]
      dt_test  <- dt_var[ESTACION == est_oculta]
      fechas_test <- unique(dt_test$FECHA)

      for (f in fechas_test) {
        train_dia <- dt_train[FECHA == f]
        test_dia  <- dt_test[FECHA == f]

        if (nrow(train_dia) < 3) next

        train_dia[, var_value := get(var_name)]
        test_dia[, var_value  := get(var_name)]

        sf_train <- st_as_sf(train_dia, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
          st_transform(25830)
        sf_test  <- st_as_sf(test_dia, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
          st_transform(25830)

        valor_real <- sf_test[["var_value"]]

        # 1. Media global
        pred_media <- mean(sf_train[["var_value"]], na.rm = TRUE)

        # 2. Vecino más cercano (1-NN)
        mod_nn  <- idw(f_formula, sf_train, sf_test, nmax = 1, debug.level = 0)
        pred_nn <- mod_nn$var1.pred

        # 3. kNN (promedio no ponderado de los k vecinos más cercanos)
        k_real <- min(k_vecinos, nrow(sf_train))
        mod_knn  <- idw(f_formula, sf_train, sf_test, nmax = k_real, idp = 0, debug.level = 0)
        pred_knn <- mod_knn$var1.pred

        # 4. IDW (p = 2)
        mod_idw  <- idw(f_formula, sf_train, sf_test, idp = 2, debug.level = 0)
        pred_idw <- mod_idw$var1.pred

        res_var[[length(res_var) + 1]] <- data.table(
          Real      = valor_real,
          Pred_Media = pred_media,
          Pred_NN    = pred_nn,
          Pred_kNN   = pred_knn,
          Pred_IDW   = pred_idw
        )
      }
    }

    cli_progress_done(id = pb)

    dt_res <- rbindlist(res_var)

    rmse <- function(real, pred) sqrt(mean((real - pred)^2, na.rm = TRUE))
    mae  <- function(real, pred) mean(abs(real - pred), na.rm = TRUE)

    rmses <- c(
      "Media"           = rmse(dt_res$Real, dt_res$Pred_Media),
      "Vecino Cercano"  = rmse(dt_res$Real, dt_res$Pred_NN),
      "kNN"             = rmse(dt_res$Real, dt_res$Pred_kNN),
      "IDW"             = rmse(dt_res$Real, dt_res$Pred_IDW)
    )
    maes <- c(
      "Media"           = mae(dt_res$Real, dt_res$Pred_Media),
      "Vecino Cercano"  = mae(dt_res$Real, dt_res$Pred_NN),
      "kNN"             = mae(dt_res$Real, dt_res$Pred_kNN),
      "IDW"             = mae(dt_res$Real, dt_res$Pred_IDW)
    )
    mejor_metodo <- names(rmses)[which.min(rmses)]

    resultados_globales[[var_name]] <- data.table(
      Variable       = var_name,
      N_Estaciones   = n_est,
      N_Predicciones = nrow(dt_res),
      RMSE_Media     = round(rmses["Media"], 3),
      RMSE_NN        = round(rmses["Vecino Cercano"], 3),
      RMSE_kNN       = round(rmses["kNN"], 3),
      RMSE_IDW       = round(rmses["IDW"], 3),
      MAE_Media      = round(maes["Media"], 3),
      MAE_NN         = round(maes["Vecino Cercano"], 3),
      MAE_kNN        = round(maes["kNN"], 3),
      MAE_IDW        = round(maes["IDW"], 3),
      Mejor_Metodo   = mejor_metodo
    )

    cli_alert_success(
      "{var_name}: mejor = {mejor_metodo} (RMSE {round(rmses[mejor_metodo], 3)}) | {nrow(dt_res)} predicciones"
    )
  }

  tabla_final <- rbindlist(resultados_globales)

  cli_h2("Resultados finales - LOOCV")
  print(tabla_final)

  return(tabla_final)
}

