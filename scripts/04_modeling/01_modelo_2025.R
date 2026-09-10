# ==============================================================================
# MODELO ESPACIAL INLA-SPDE - 2025
# ==============================================================================
# Un unico campo espacial Matern compartido por todas las observaciones.
# Las decisiones del modelo se modifican solo en el bloque CONFIGURACION.

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(here)
})

# ==============================================================================
# CONFIGURACION EDITABLE
# ==============================================================================
ANIO <- 2025
ESCALA <- "DIARIO"

# Periodo usado. Escribir NULL en ambos para utilizar todo el fichero.
FECHA_INICIO <- as.Date("2025-01-01")
FECHA_FIN <- as.Date("2025-12-31")

# Respuesta y covariables usadas para seleccionar y limpiar los datos. Si se
# cambia este vector, debe hacerse el mismo cambio en la formula del apartado 4.
RESPUESTA <- "LOG_NO2_DIARIO"
COVARIABLES <- c(
  "intensidad",
  "Temperatura",
  "Velocidad_Viento"
)

# Malla elegida en el analisis de sensibilidad.
MALLA <- "media" # "gruesa" | "media" | "fina"

# Hiperparametros del campo Matern (PC priors).
# prior.range = c(r0, p0): P(rango < r0) = p0, con el rango en km.
# prior.sigma = c(s0, p0): P(sigma > s0) = p0, en la escala de RESPUESTA.
PRIOR_RANGE <- c(9.3, 0.5)
PRIOR_SIGMA <- c(0.6, 0.2)
SPDE_ALPHA <- 2
SPDE_CONSTR <- FALSE

# Verosimilitud e integracion. Para fijar un prior de la precision gaussiana,
# sustituir list() por, por ejemplo:
# list(hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
FAMILIA <- "gaussian"
CONTROL_FAMILY <- list()
ESTRATEGIA_INTEGRACION <- "eb"
NUM_THREADS <- 4L

# Diagnosticos y ficheros de salida.
CALCULAR_CPO <- TRUE
CALCULAR_LOSO <- TRUE
GUARDAR_MODELO <- TRUE

DIR_OUT <- here("outputs", "modelo", paste0("modelo_", ANIO))
DIR_MODELO <- here("data", "processed", "Modelos", paste0("modelo_", ANIO))
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MODELO, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. DATOS
# ==============================================================================
ruta_datos <- here(
  "data", "processed", "Maestro", as.character(ANIO),
  sprintf("dataset_maestro_inla_%d_%s.rds", ANIO, ESCALA)
)
ruta_malla <- here(
  "data", "processed", "Malla", "NO2",
  sprintf("malla_spde_madrid_%s.rds", MALLA)
)

if (!file.exists(ruta_datos)) stop("No existe el maestro: ", ruta_datos)
if (!file.exists(ruta_malla)) stop("No existe la malla: ", ruta_malla)

df <- readRDS(ruta_datos)
setDT(df)

columnas_necesarias <- unique(c(
  "ESTACION", "FECHA", "X_km", "Y_km", RESPUESTA, COVARIABLES
))
columnas_ausentes <- setdiff(columnas_necesarias, names(df))
if (length(columnas_ausentes) > 0L) {
  stop("Faltan columnas en el maestro: ", paste(columnas_ausentes, collapse = ", "))
}

df[, FECHA := as.Date(FECHA)]
if (!is.null(FECHA_INICIO)) df <- df[FECHA >= FECHA_INICIO]
if (!is.null(FECHA_FIN)) df <- df[FECHA <= FECHA_FIN]

n_antes <- nrow(df)
df <- df[complete.cases(df[, c(RESPUESTA, COVARIABLES, "X_km", "Y_km"), with = FALSE])]
if (nrow(df) == 0L) stop("No quedan observaciones completas para ajustar el modelo.")

setorder(df, FECHA, ESTACION)
df[, y := get(RESPUESTA)]
coords <- as.matrix(df[, .(X_km, Y_km)])

cat("\n--- Datos del modelo ---\n")
cat("Periodo:", format(min(df$FECHA)), "a", format(max(df$FECHA)), "\n")
cat("Filas:", nrow(df), "| eliminadas por NA:", n_antes - nrow(df), "\n")
cat("Estaciones:", uniqueN(df$ESTACION), "| respuesta:", RESPUESTA, "\n")
cat("Covariables:", paste(COVARIABLES, collapse = ", "), "\n")

# ==============================================================================
# 2. MALLA Y CAMPO ESPACIAL
# ==============================================================================
mesh <- readRDS(ruta_malla)
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = SPDE_ALPHA,
  prior.range = PRIOR_RANGE,
  prior.sigma = PRIOR_SIGMA,
  constr = SPDE_CONSTR
)

indice_campo <- inla.spde.make.index(
  name = "campo_espacial",
  n.spde = spde$n.spde
)
A_campo <- inla.spde.make.A(mesh = mesh, loc = coords)

stopifnot(
  "Las filas de A no coinciden con los datos" = nrow(A_campo) == nrow(df),
  "Las columnas de A no coinciden con la malla" = ncol(A_campo) == mesh$n
)

# ==============================================================================
# 3. STACK
# ==============================================================================
efectos_fijos <- cbind(
  data.frame(Intercept = rep(1, nrow(df))),
  as.data.frame(df[, COVARIABLES, with = FALSE], check.names = FALSE)
)

stack_modelo <- inla.stack(
  data = list(y_response = df$y),
  A = list(A_campo, 1),
  effects = list(indice_campo, efectos_fijos),
  tag = "estimacion"
)

cat("Malla:", MALLA, "| vertices:", mesh$n, "\n")
cat(sprintf(
  "Prior Matern: P(range < %.2f km) = %.3f; P(sigma > %.3f) = %.3f\n\n",
  PRIOR_RANGE[1], PRIOR_RANGE[2], PRIOR_SIGMA[1], PRIOR_SIGMA[2]
))

# ==============================================================================
# 4. AJUSTE
# ==============================================================================
t0 <- Sys.time()
modelo_2025 <- inla(
  y_response ~ -1 + Intercept +
    intensidad + Temperatura + Velocidad_Viento +
    f(campo_espacial, model = spde),
  data = inla.stack.data(stack_modelo, spde = spde),
  family = FAMILIA,
  control.family = CONTROL_FAMILY,
  control.predictor = list(
    A = inla.stack.A(stack_modelo),
    compute = TRUE
  ),
  control.compute = list(
    cpo = CALCULAR_CPO,
    dic = TRUE,
    waic = TRUE,
    return.marginals.predictor = TRUE
  ),
  control.inla = list(int.strategy = ESTRATEGIA_INTEGRACION),
  num.threads = NUM_THREADS
)
minutos <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# ==============================================================================
# 5. RESUMEN POSTERIOR Y VALIDACION
# ==============================================================================
spde_result <- inla.spde2.result(
  inla = modelo_2025,
  name = "campo_espacial",
  spde = spde,
  do.transf = TRUE
)

marginal_rango <- spde_result$marginals.range.nominal[[1]]
marginal_sigma <- inla.tmarginal(
  sqrt,
  spde_result$marginals.variance.nominal[[1]]
)
resumen_rango <- inla.zmarginal(marginal_rango, silent = TRUE)
resumen_sigma <- inla.zmarginal(marginal_sigma, silent = TRUE)

metricas <- data.table(
  year = ANIO,
  scale = ESCALA,
  mesh = MALLA,
  vertices = mesh$n,
  observations = nrow(df),
  stations = uniqueN(df$ESTACION),
  DIC = modelo_2025$dic$dic,
  WAIC = modelo_2025$waic$waic,
  range_mean_km = resumen_rango$mean,
  range_q025_km = resumen_rango$quant0.025,
  range_q975_km = resumen_rango$quant0.975,
  spatial_sd_mean = resumen_sigma$mean,
  spatial_sd_q025 = resumen_sigma$quant0.025,
  spatial_sd_q975 = resumen_sigma$quant0.975,
  minutes = minutos
)

predicciones_loso <- NULL
if (CALCULAR_LOSO) {
  filas_por_estacion <- split(seq_len(nrow(df)), df$ESTACION)
  grupos_loso <- filas_por_estacion[as.character(df$ESTACION)]
  gcv <- inla.group.cv(result = modelo_2025, groups = grupos_loso)

  error_loso <- gcv$mean - df$y
  dentro_95 <- gcv$mean - 1.96 * gcv$sd <= df$y &
    df$y <= gcv$mean + 1.96 * gcv$sd

  metricas[, `:=`(
    RMSE_loso = sqrt(mean(error_loso^2, na.rm = TRUE)),
    MAE_loso = mean(abs(error_loso), na.rm = TRUE),
    Cov95_loso = 100 * mean(dentro_95, na.rm = TRUE)
  )]

  predicciones_loso <- df[, .(
    ESTACION, FECHA, observed = y
  )]
  predicciones_loso[, `:=`(
    predicted_mean = gcv$mean,
    predicted_sd = gcv$sd,
    lower_95 = gcv$mean - 1.96 * gcv$sd,
    upper_95 = gcv$mean + 1.96 * gcv$sd,
    error = error_loso,
    covered_95 = dentro_95
  )]
}

if (CALCULAR_CPO) {
  metricas[, `:=`(
    LCPO = mean(-log(modelo_2025$cpo$cpo), na.rm = TRUE),
    PIT_mean = mean(modelo_2025$cpo$pit, na.rm = TRUE),
    PIT_sd = sd(modelo_2025$cpo$pit, na.rm = TRUE)
  )]
}

# ==============================================================================
# 6. SALIDAS
# ==============================================================================
fwrite(metricas, file.path(DIR_OUT, "metricas_modelo_2025.csv"))
fwrite(
  as.data.table(modelo_2025$summary.fixed, keep.rownames = "term"),
  file.path(DIR_OUT, "efectos_fijos_modelo_2025.csv")
)
fwrite(
  as.data.table(modelo_2025$summary.hyperpar, keep.rownames = "hyperparameter"),
  file.path(DIR_OUT, "hiperparametros_modelo_2025.csv")
)
if (!is.null(predicciones_loso)) {
  fwrite(predicciones_loso, file.path(DIR_OUT, "predicciones_loso_modelo_2025.csv"))
}

configuracion <- list(
  year = ANIO,
  scale = ESCALA,
  start_date = FECHA_INICIO,
  end_date = FECHA_FIN,
  response = RESPUESTA,
  covariates = COVARIABLES,
  mesh = MALLA,
  prior_range = PRIOR_RANGE,
  prior_sigma = PRIOR_SIGMA,
  spde_alpha = SPDE_ALPHA,
  spde_constr = SPDE_CONSTR,
  family = FAMILIA,
  control_family = CONTROL_FAMILY,
  integration_strategy = ESTRATEGIA_INTEGRACION,
  num_threads = NUM_THREADS
)

if (GUARDAR_MODELO) {
  saveRDS(
    list(
      model = modelo_2025,
      mesh = mesh,
      spde = spde,
      stack = stack_modelo,
      data = df,
      configuration = configuracion,
      metrics = metricas
    ),
    file.path(DIR_MODELO, sprintf("modelo_espacial_%d_%s_malla_%s.rds", ANIO, ESCALA, MALLA))
  )
}

cat("\n--- Modelo completado ---\n")
print(metricas)
cat("Resultados:", DIR_OUT, "\n")
if (GUARDAR_MODELO) cat("Modelo:", DIR_MODELO, "\n")
