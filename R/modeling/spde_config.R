# ==============================================================================
# DEFINICIÓN ÚNICA DEL CAMPO ESPACIAL MATÉRN (SPDE) PARA TODO EL TFM
# ==============================================================================
# Todos los scripts que ajusten un modelo INLA-SPDE deben construir el campo con
# crear_spde(). El objetivo es que el prior sea una decisión documentada en UN
# solo sitio: si cada script define el suyo, dos modelos dejan de ser comparables
# porque las diferencias de DIC/WAIC pueden venir del prior y no del modelo.
#
# PC priors (Fuglstad et al. 2019). Semántica de los argumentos:
#   prior.range = c(r0, p) -> P(rango < r0) = p
#   prior.sigma = c(s0, p) -> P(sigma > s0) = p
#
# Justificación de los valores por defecto (escala diaria, ver script
# 02_eda/29_distancias_estaciones_prior_spde.R y el diagnóstico de Spde.R):
#   - rango: el diámetro del dominio de estudio es D = 28 km y se toma la
#     mediana a priori en D/3 = 9.3 km, la regla habitual cuando no se conoce
#     el rango real. Con c(9.3, 0.5) la mediana a priori es exactamente 9.3 km.
#   - sigma: la escala de referencia es la sd de LOG_NO2. Del diagnóstico de
#     Spde.R: 0.59 la sd total, 0.32 la residual intra-instante y 0.26 la de
#     las medias por estación —esta última es la parte que el campo espacial
#     puede explicar—. Se toma s0 = 0.6 con p = 0.5, es decir MEDIANA a priori
#     = 0.6: por encima de la variabilidad espacial observada, para que el
#     prior no arrastre la sigma hacia abajo.
#
#     CUIDADO CON LA p. La mediana a priori NO es s0, es s0*log(2)/(-log(p)):
#         c(0.6, 0.50) -> mediana 0.600   <- el que se usa
#         c(0.6, 0.10) -> mediana 0.181
#         c(0.6, 0.01) -> mediana 0.090   <- aplasta el campo contra cero
#     Poner c(0.6, 0.01) equivale a afirmar "casi seguro que sigma vale ~0.09",
#     una décima parte de lo que se pretendía. Es el mismo error que documenta
#     modelo_1_2023_2025_diario.R, donde c(1, 0.01) "apagaba el campo: SD~0,
#     rango enorme" y hubo que relajarlo a c(1, 0.5).
# ==============================================================================

PRIORS_SPDE <- list(
  DIARIO = list(prior.range = c(9.3, 0.5), prior.sigma = c(0.6, 0.5)),
  HORARIO = list(prior.range = c(9.3, 0.5), prior.sigma = c(0.6, 0.5)),
  MENSUAL = list(prior.range = c(9.3, 0.5), prior.sigma = c(0.6, 0.5))
)

#' Construye el objeto SPDE con los PC priors del proyecto.
#'
#' @param malla objeto de malla (fm_mesh_2d / inla.mesh).
#' @param escala "DIARIO", "HORARIO" o "MENSUAL". Selecciona el juego de priors.
#' @param prior.range,prior.sigma opcionales. Solo para análisis de sensibilidad:
#'   sobrescriben los priors del proyecto. En los modelos finales NO se usan.
#'
#' @return objeto inla.spde2.
crear_spde <- function(malla, escala = "DIARIO",
                       prior.range = NULL, prior.sigma = NULL) {
  if (!escala %in% names(PRIORS_SPDE)) {
    stop("Escala '", escala, "' sin priors definidos. Válidas: ",
         paste(names(PRIORS_SPDE), collapse = ", "))
  }
  p <- PRIORS_SPDE[[escala]]
  if (!is.null(prior.range)) p$prior.range <- prior.range
  if (!is.null(prior.sigma)) p$prior.sigma <- prior.sigma

  inla.spde2.pcmatern(
    mesh = malla,
    alpha = 2, # nu = 1 en 2D
    prior.range = p$prior.range,
    prior.sigma = p$prior.sigma,
    constr = FALSE
  )
}

#' Mediana a priori implícita en un PC prior, para documentar la elección.
#'
#' Rango (d = 2): P(rango < r) = exp(-lambda / r), lambda = -r0 * log(p)
#' Sigma        : P(sigma > s) = exp(-lambda * s), lambda = -log(p) / s0
resumen_priors <- function(escala = "DIARIO") {
  p <- PRIORS_SPDE[[escala]]
  lam_r <- -p$prior.range[1] * log(p$prior.range[2])
  lam_s <- -log(p$prior.sigma[2]) / p$prior.sigma[1]
  data.frame(
    parametro = c("rango", "sigma"),
    prior = c(
      sprintf("c(%.2f, %.2f)", p$prior.range[1], p$prior.range[2]),
      sprintf("c(%.2f, %.2f)", p$prior.sigma[1], p$prior.sigma[2])
    ),
    interpretacion = c(
      sprintf("P(rango < %.2f km) = %.2f", p$prior.range[1], p$prior.range[2]),
      sprintf("P(sigma > %.2f) = %.2f", p$prior.sigma[1], p$prior.sigma[2])
    ),
    mediana_a_priori = c(lam_r / log(2), log(2) / lam_s)
  )
}
