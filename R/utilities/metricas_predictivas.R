# ==============================================================================
# RMSE Y COV95 PARA MODELOS INLA GAUSSIANOS
# ==============================================================================
# RMSE y COV95 se calculaban de forma repetida y ligeramente distinta en cada
# script (01_modelo_2019_2021.R, simulacion.R, comparacion_mallas.R,
# prueba_transformaciones.R...). Esta es la version unica.
#
# Punto importante sobre COV95: con family = "gaussian", la sd de
# summary.fitted.values (o de modelo$summary.linear.predictor) solo recoge la
# incertidumbre de la MEDIA ajustada (intercepto + efectos fijos + RW2 +
# campo espacial). No incluye el ruido de observacion, es decir, la varianza
# residual gaussiana (1 / "Precision for the Gaussian observations"). Como
# y = media_ajustada + ruido, el intervalo para predecir una observacion
# individual necesita sd_total = sqrt(sd_media^2 + varianza_residual); sin ese
# termino el intervalo queda demasiado estrecho y COV95 sale muy por debajo de
# 95 aunque el modelo este bien especificado (ver 01_modelo_2019_2021.R).
#
# inla.group.cv() (validacion LOSO) ya devuelve una sd predictiva completa, asi
# que ahi varianza_residual debe quedarse en 0 (el valor por defecto).
# ==============================================================================

#' Varianza residual gaussiana de un modelo INLA.
#'
#' Media posterior de 1 / precision, calculada sobre la marginal completa con
#' inla.emarginal(). Es mas correcto que 1 / mean(precision), que ignora la
#' asimetria de la posterior de la precision.
#'
#' @param modelo objeto devuelto por inla() con family = "gaussian".
#' @return numero: varianza residual, en la escala de la respuesta modelada.
varianza_residual_gaussiana <- function(modelo) {
    nombre <- grep(
        "^Precision for the Gaussian observations",
        names(modelo$marginals.hyperpar),
        value = TRUE
    )
    if (length(nombre) != 1L) {
        stop("No se ha podido identificar la precision residual gaussiana.")
    }
    inla.emarginal(
        function(precision) 1 / precision,
        modelo$marginals.hyperpar[[nombre]]
    )
}

#' RMSE y cobertura del 95% de un conjunto de predicciones gaussianas.
#'
#' @param y_obs vector de valores observados.
#' @param media_pred vector de medias predichas (misma longitud que y_obs).
#' @param sd_pred vector de desviaciones tipicas de la prediccion.
#'   - Con summary.fitted.values$sd (ajuste en muestra): solo recoge la
#'     incertidumbre de la media ajustada. Pasa varianza_residual (ver
#'     varianza_residual_gaussiana()) para que COV95 sea comparable a una
#'     cobertura predictiva real.
#'   - Con inla.group.cv()$sd (validacion LOSO): ya es predictiva completa,
#'     deja varianza_residual = 0.
#' @param varianza_residual varianza residual gaussiana que se suma a
#'   sd_pred^2 antes de construir el intervalo. Por defecto 0.
#' @param nivel nivel de confianza del intervalo. Por defecto 0.95.
#'
#' @return data.table con RMSE, COV95 (en %) y N_validas (observaciones
#'   finitas usadas en el calculo).
calcular_rmse_cov95 <- function(y_obs, media_pred, sd_pred,
                                 varianza_residual = 0, nivel = 0.95) {
    validas <- is.finite(y_obs) & is.finite(media_pred) & is.finite(sd_pred)
    if (!any(validas)) {
        stop("No hay observaciones finitas para calcular RMSE/COV95.")
    }

    y_obs <- y_obs[validas]
    media_pred <- media_pred[validas]
    sd_total <- sqrt(sd_pred[validas]^2 + varianza_residual)

    z <- qnorm(1 - (1 - nivel) / 2)
    dentro_intervalo <- y_obs >= media_pred - z * sd_total &
        y_obs <= media_pred + z * sd_total

    data.table(
        RMSE = sqrt(mean((media_pred - y_obs)^2)),
        COV95 = 100 * mean(dentro_intervalo),
        N_validas = sum(validas)
    )
}
