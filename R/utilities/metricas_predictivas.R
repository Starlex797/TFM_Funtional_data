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
# La funcion que extrae la varianza residual de INLA se encuentra en
# R/modeling/inla_modeling.R. Este archivo conserva solamente metricas que se
# pueden aplicar a cualquier conjunto de predicciones gaussianas.
# ==============================================================================

#' RMSE y cobertura del 95% de un conjunto de predicciones gaussianas.
#'
#' @param y_obs vector de valores observados.
#' @param media_pred vector de medias predichas (misma longitud que y_obs).
#' @param sd_pred vector de desviaciones tipicas de la prediccion.
#'   - Con summary.fitted.values$sd (ajuste en muestra): solo recoge la
#'     incertidumbre de la media ajustada. Pasa varianza_residual (ver
#'     varianza_residual_gaussiana()) para que COV95 sea comparable a una
#'     cobertura predictiva real.
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
