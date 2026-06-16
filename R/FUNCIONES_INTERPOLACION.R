# ==============================================================================
# FUNCIONES DE INTERPOLACIÓN ESPACIAL
# ==============================================================================
# Funciones para interpolar variables climáticas mediante IDW (Inverse Distance
# Weighting) desde las estaciones meteorológicas a las ubicaciones de las
# estaciones de calidad del aire, a nivel DIARIO.
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
#'   las ubicaciones donde se quieren estimar los valores (sin FECHA; se
#'   interpola para cada fecha disponible en dt_meteo).
#' @param variables character vector con los nombres de las columnas climáticas
#'   a interpolar.
#' @param idp numeric. Potencia de la distancia inversa (por defecto 2).
#' @param min_estaciones integer. Número mínimo de estaciones con datos
#'   requerido para interpolar (por defecto 3). Si hay menos, se devuelve NA.
#' @param crs_orig integer. CRS de las coordenadas de entrada (por defecto
#'   4326 WGS84).
#' @param crs_proj integer. CRS proyectado para el cálculo de distancias
#'   (por defecto 25830 UTM 30N).
#'
#' @return data.table con columnas ESTACION, FECHA y las variables climáticas
#'   interpoladas.
#'   
interpolar_idw_clima <- function(dt_meteo,
                                 dt_objetivo,
                                 variables = c("Temperatura", "Humedad_Relativa","Precipitaciones",
                                               "Presion_Barometrica", "Radiacion_Solar"),
                                 idp       = 2,
                                 min_estaciones = 3,
                                 crs_orig  = 4326,
                                 crs_proj  = 25830) {

  # --- 0. Validar variables disponibles ---
  vars_disponibles <- intersect(variables, names(dt_meteo))
  if (length(vars_disponibles) == 0) {
    stop("Ninguna de las variables solicitadas existe en dt_meteo: ",
         paste(variables, collapse = ", "))
  }
  if (length(vars_disponibles) < length(variables)) {
    warning("Variables no encontradas en dt_meteo: ",
            paste(setdiff(variables, vars_disponibles), collapse = ", "))
  }

  # --- 1. Preparar ubicaciones objetivo (sf proyectado, fijas) ---
  coords_objetivo <- unique(dt_objetivo[, .(ESTACION, LONGITUD, LATITUD)])
  sf_objetivo <- st_as_sf(coords_objetivo,
                           coords = c("LONGITUD", "LATITUD"),
                           crs = crs_orig) |>
    st_transform(crs = crs_proj)

  # --- 2. Coordenadas únicas de estaciones meteo (sf proyectado) ---
  coords_meteo <- unique(dt_meteo[, .(ESTACION, LONGITUD, LATITUD)])
  sf_meteo_base <- st_as_sf(coords_meteo,
                              coords = c("LONGITUD", "LATITUD"),
                              crs = crs_orig) |>
    st_transform(crs = crs_proj)

  # --- 3. Iterar sobre cada fecha ---
  fechas <- sort(unique(dt_meteo$FECHA))
  n_fechas <- length(fechas)
  n_obj <- nrow(sf_objetivo)

  # Pre-alocar lista de resultados
  resultados <- vector("list", n_fechas)

  for (i in seq_along(fechas)) {
    f <- fechas[i]

    # Datos meteo de este día
    dt_dia <- dt_meteo[FECHA == f, c("ESTACION", vars_disponibles), with = FALSE]

    # Plantilla de salida para este día
    dt_out <- data.table(ESTACION = coords_objetivo$ESTACION, FECHA = f)

    for (var_name in vars_disponibles) {
      # Unir valores del día con coordenadas meteo
      dt_var <- merge(dt_dia[, .(ESTACION, valor = get(var_name))],
                      coords_meteo, by = "ESTACION")
      dt_var <- dt_var[!is.na(valor)]

      if (nrow(dt_var) < min_estaciones) {
        dt_out[, (var_name) := NA_real_]
        next
      }

      # Construir sf con valores del día usando nombre seguro
      sf_obs <- sf_meteo_base[sf_meteo_base$ESTACION %in% dt_var$ESTACION, ]
      sf_obs[["IDW_VAR"]] <- dt_var$valor[match(sf_obs$ESTACION, dt_var$ESTACION)]

      # IDW con nombre temporal para evitar problemas con espacios/tildes
      idw_result <- idw(IDW_VAR ~ 1,
                        locations = sf_obs,
                        newdata   = sf_objetivo,
                        idp       = idp)
      dt_out[, (var_name) := idw_result$var1.pred]
    }

    resultados[[i]] <- dt_out

    # Progreso cada 50 días
    if (i %% 50 == 0 || i == n_fechas) {
      message(sprintf("Interpolación IDW: %d / %d fechas completadas", i, n_fechas))
    }
  }

  # --- 4. Combinar y reportar cobertura ---
  dt_final <- rbindlist(resultados)

  for (var_name in vars_disponibles) {
    n_ok <- sum(!is.na(dt_final[[var_name]]))
    n_total <- nrow(dt_final)
    message(sprintf("'%s': %d / %d valores interpolados (%.1f%%)",
                    var_name, n_ok, n_total, 100 * n_ok / n_total))
  }

  dt_final
}

