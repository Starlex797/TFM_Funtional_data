# ==============================================================================
# COMPARACION DE MALLAS SPDE - 2025 DIARIO
# ==============================================================================
# Objetivo: decidir con que malla trabajar y ver el rango y la sigma espacial.
#
# CRITERIO. La malla NO se elige por el mejor DIC/WAIC: al refinarla aumenta la
# dimension del campo latente y casi siempre "mejora" el ajuste. Se elige por
# ESTABILIDAD: la malla mas gruesa a partir de la cual el rango, la sigma y el
# error dejan de moverse.
#
# ==============================================================================

library(INLA)
library(data.table)
library(here)

# ------------------------------------------------------------------------------
# CONFIGURACION
# ------------------------------------------------------------------------------
ANIO <- 2025
ESCALA <- "DIARIO"
MALLAS <- c("gruesa", "media", "fina")

# Ventana de dias CONSECUTIVOS. Tiene que ser un bloque seguido: el AR1 modela
# la correlacion entre un dia y el SIGUIENTE. Si se espaciaran los dias, rho
# pasaria a medir la correlacion entre dias no contiguos y el efecto
# autorregresivo dejaria de ser interpretable.
# No se usa el ano entero porque el AR1 replica el campo espacial una vez por
# dia: con 365 dias y la malla fina serian ~213.000 nodos latentes.
# Se coge invierno, que es cuando hay episodios de NO2 y la estructura espacial
# esta mas marcada, que es justo lo que tiene que resolver la malla.
FECHA_INICIO <- as.Date("2025-01-01")
N_DIAS <- 60

# Estas son las que van en la formula. Si se anade o quita una aqui, hay que
# tocarla tambien en las dos formulas del bucle (estn escritas a mano, a
# proposito, para que se lea que entra en cada modelo).
COVARIABLES <- c("intensidad", "Temperatura", "Velocidad_Viento")

# PC priors fijos para todas las mallas.
PRIOR_RANGE <- c(9.3, 0.5) # P(rango < 9.3 km) = 0.5
PRIOR_SIGMA <- c(0.6, 0.2) # P(sigma > 0.6) = 0.05

# El analisis de sensibilidad ajusta el mismo modelo espacial en cada malla.
# La comparacion M1-M2 se hace despues contra un unico modelo sin campo espacial.
TIPOS <- c("espacial")

# VALIDACION: leave-one-station-out sobre las 24 estaciones.
#
# Se deja fuera una ESTACION entera cada vez, no una observacion suelta. Si se
# quitara solo la medicion de la estacion X del dia 12, la estacion X seguiria
# presente los otros 59 dias y el campo la reconstruiria casi perfectamente: se
# estaria midiendo la suavidad temporal, no la capacidad de interpolar donde no
# hay sensor, que es lo unico que justifica el campo espacial.
#
# No se hace reajustando 24 veces por malla (72 ajustes). Se usa
# inla.group.cv(), que da la predictiva leave-group-out a partir de UN solo
# ajuste. Es una aproximacion: condiciona en la posterior de los hiperparametros
# del ajuste completo en vez de reestimarlos en cada fold. Quitando 1 de 24
# estaciones eso apenas mueve theta, asi que la aproximacion es buena.

DIR_OUT <- here("outputs", "modelo", "comparacion_mallas")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PASO 1: DATOS
# ==============================================================================
df <- readRDS(here(
  "data", "processed", "Maestro", as.character(ANIO),
  paste0("dataset_maestro_inla_", ANIO, "_", ESCALA, ".rds")
))
setDT(df)

# La respuesta ya viene transformada: LOG_NO2_DIARIO es log1p(DATO_DIARIO).
# No hay que volver a aplicar ningun logaritmo.
df[, y := LOG_NO2_DIARIO]
df[, NO2 := DATO_DIARIO]

# Ventana consecutiva [FECHA_INICIO, FECHA_INICIO + N_DIAS - 1]
FECHA_FIN <- FECHA_INICIO + N_DIAS - 1L
df <- df[FECHA >= FECHA_INICIO & FECHA <= FECHA_FIN]
if (nrow(df) == 0) stop("La ventana de fechas no tiene datos en el maestro.")

# Filas sin respuesta o sin covariables fuera: todas las mallas se ajustan
# exactamente sobre el mismo conjunto de filas.
df <- df[complete.cases(df[, c("y", COVARIABLES), with = FALSE])]

# Indice temporal del AR1: se calcula desde la FECHA, NO con match() sobre los
# dias presentes. Asi, si algun dia se quedara sin observaciones, deja un hueco
# en la numeracion en vez de "cerrar filas" y hacer pasar por consecutivos dos
# dias que en el calendario no lo son.
df[, ID_TIEMPO := as.integer(FECHA - FECHA_INICIO) + 1L]
setorder(df, ID_TIEMPO, ESTACION) # a partir de aqui NO se reordena mas
n_dias <- N_DIAS # grupos del AR1 = dias del calendario

# ------------------------------------------------------------------------------
# Grupos del leave-one-station-out
# ------------------------------------------------------------------------------
# grupos[[i]] = las filas que se quitan para predecir la fila i. Como es la
# estacion entera, la prediccion de cada fila usa solo las otras 23 estaciones.
estaciones <- sort(unique(df$ESTACION))
filas_por_estacion <- split(seq_len(nrow(df)), df$ESTACION)
grupos <- filas_por_estacion[as.character(df$ESTACION)]

coords <- as.matrix(df[, .(X_km, Y_km)]) # UTM 30N en km, ya vienen del Paso 2

cat(
  "\nVentana:", format(FECHA_INICIO), "a", format(FECHA_FIN),
  "|", n_dias, "dias consecutivos\n"
)
cat(
  "Filas:", nrow(df), "| Estaciones:", length(estaciones),
  "| Dias con datos:", uniqueN(df$ID_TIEMPO), "de", n_dias, "\n\n"
)
cat(sprintf(
  "SD NO2 original: %.3f | SD log1p(NO2): %.3f\n\n",
  sd(df$NO2, na.rm = TRUE), sd(df$y, na.rm = TRUE)
))

# ==============================================================================
# PASO 2-8: UN MODELO POR (MALLA, TIPO)
# ==============================================================================
resultados <- list()
por_estacion <- list()

for (tipo in TIPOS) {
  for (nombre_malla in MALLAS) {
    clave <- paste(tipo, nombre_malla, sep = " | ")
    cat("--- ", clave, " ---\n", sep = "")

    # --- PASO 2: malla y objeto SPDE -------------------------------------------
    mesh <- readRDS(here(
      "data", "processed", "Malla", "NO2", paste0("malla_spde_madrid_", nombre_malla, ".rds")
    ))
    spde <- inla.spde2.pcmatern(
      mesh = mesh,
      prior.range = PRIOR_RANGE,
      prior.sigma = PRIOR_SIGMA
    )

    # --- PASOS 3, 4 y 6: indice, matriz A y formula -----------------------------
    # Los tres tienen que concordar en nombre y dimension, asi que se escriben
    # juntos: si el indice lleva n.group, la matriz A tiene que llevar group y la
    # formula control.group. Separarlos es la via para que se desincronicen.
    if (tipo == "espacial") {
      # Un unico campo espacial compartido por los 60 dias.
      s.index <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)
      A.estud <- inla.spde.make.A(mesh = mesh, loc = coords)

      formula <- y_response ~ -1 + Intercept +
        intensidad + Temperatura + Velocidad_Viento +
        f(spatial.field, model = spde)
    } else {
      # Una replica del campo espacial por dia (n.group = n_dias), encadenadas
      # por un AR1. Incognitas del campo: nodos x dias.
      s.index <- inla.spde.make.index("spatial.field",
        n.spde = spde$n.spde,
        n.group = n_dias
      )
      A.estud <- inla.spde.make.A(
        mesh = mesh, loc = coords,
        group = df$ID_TIEMPO, n.group = n_dias
      )

      formula <- y_response ~ -1 + Intercept + # Añadimos el -1 para quitar el intercepto que introduce inla automaticamente
        intensidad + Temperatura + Velocidad_Viento +
        f(spatial.field,
          model = spde, group = spatial.field.group,
          control.group = list(model = "ar1")
        )
    }

    # --- PASO 5: inla.stack ----------------------------------------------------
    # A altera la dimension del predictor lineal, asi que no se puede pasar un
    # data.frame clasico a inla(). Se empaquetan respuesta, covariables y campo
    # espacial en el stack:
    #   A = list(A.estud, 1) -> A.estud proyecta el campo, 1 (identidad) el resto
    #   El Intercept va DENTRO del stack, como efecto manual.
    stack.est <- inla.stack(
      data = list(y_response = df$y), # se ajusta con TODOS los datos
      A = list(A.estud, 1),
      effects = list(
        c(s.index, list(Intercept = 1)),
        as.data.frame(df[, COVARIABLES, with = FALSE])
      ),
      tag = "est"
    )

    # --- PASO 6b: ajuste --------------------------------------------------------
    # La formula (arriba) lleva -1 para quitar el intercepto automatico de R y
    # usar el Intercept que viaja dentro del stack.
    t0 <- Sys.time()
    m1 <- inla(
      formula,
      data = inla.stack.data(stack.est, spde = spde),
      family = "gaussian", # la respuesta ya esta en log: gaussiana, no gamma
      control.predictor = list(A = inla.stack.A(stack.est), compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, return.marginals.predictor = TRUE),
      control.inla = list(int.strategy = "eb") # Empirical Bayes: mas rapido
    )
    minutos <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

    # --- PASO 7: rango y varianza en escala fisica (km) ------------------------
    # INLA estima el campo en escala logaritmica interna (Theta1, Theta2).
    # inla.spde2.result lo devuelve a rango (km) y varianza marginal.
    spde.result <- inla.spde2.result(
      inla = m1, name = "spatial.field",
      spde = spde, do.transf = TRUE
    )

    rango <- inla.zmarginal(spde.result$marginals.range.nominal[[1]], silent = TRUE)
    varia <- inla.zmarginal(spde.result$marginals.variance.nominal[[1]], silent = TRUE)

    # --- PASO 8: validacion leave-one-station-out ------------------------------
    # gcv$mean y gcv$sd son la media y la sd de la distribucion PREDICTIVA de
    # cada observacion habiendo quitado su estacion entera. La sd ya incluye el
    # ruido de observacion (a diferencia de summary.fitted.values$sd, que solo
    # recoge la incertidumbre de la media ajustada), asi que Cov95 se calcula
    # directamente con ella.
    gcv <- inla.group.cv(result = m1, groups = grupos)

    err <- gcv$mean - df$y

    # RMSE y MAE en la escala log de la respuesta (LOG_NO2_DIARIO). No se deshace
    # el log1p: en la escala original el error lo dominarian los dias de episodio
    # y dejaria de medir la calidad de la interpolacion espacial.
    rmse <- sqrt(mean(err^2))
    mae <- mean(abs(err))

    # Cov95: % de observaciones dentro de su intervalo predictivo del 95%.
    # El RMSE mide si el punto acierta; Cov95, si la incertidumbre esta bien
    # calibrada. Muy por debajo de 95 = modelo demasiado confiado; muy por
    # encima = intervalos inutilmente anchos.
    cov95 <- 100 * mean(gcv$mean - 1.96 * gcv$sd <= df$y & df$y <= gcv$mean + 1.96 * gcv$sd)
    lcpo <- mean(-log(m1$cpo$cpo), na.rm = TRUE)
    pit_mean <- mean(m1$cpo$pit, na.rm = TRUE)
    pit_sd <- sd(m1$cpo$pit, na.rm = TRUE)

    # Rho del AR1. Solo existe en el modelo espacio-temporal; en el espacial
    # queda NA. Al ser dias consecutivos SI es la autocorrelacion de un dia al
    # siguiente. Si sale pegado a 1, el AR1 ha degenerado en paseo aleatorio y el
    # modelo esta absorbiendo los datos en vez de explicarlos.
    fila_rho <- grep("GroupRho", rownames(m1$summary.hyperpar))[1]
    rho_ar1 <- if (is.na(fila_rho)) NA_real_ else m1$summary.hyperpar$mean[fila_rho]

    # Error por estacion: dice DONDE falla el campo. Se espera que las peores
    # sean las aisladas (El Pardo) y las mejores las rodeadas de vecinas.
    por_estacion[[clave]] <- data.table(
      tipo = tipo, malla = nombre_malla, ESTACION = df$ESTACION, err = err
    )[, .(RMSE = sqrt(mean(err^2)), MAE = mean(abs(err)), n = .N),
      by = .(tipo, malla, ESTACION)
    ]

    resultados[[clave]] <- data.table(
      tipo       = tipo,
      malla      = nombre_malla,
      n_nodos    = mesh$n,
      n_latente  = ncol(A.estud), # incognitas del campo: nodos, o nodos x dias
      DIC        = m1$dic$dic,
      WAIC       = m1$waic$waic,
      LCPO       = lcpo,
      PIT_mean   = pit_mean,
      PIT_sd     = pit_sd,
      rango_km   = rango$mean,
      rango_q025 = rango$quant0.025,
      rango_q975 = rango$quant0.975,
      sigma      = sqrt(varia$mean),
      rho_ar1    = rho_ar1,
      RMSE_loso  = rmse,
      MAE_loso   = mae,
      Cov95_loso = cov95,
      minutos    = round(minutos, 1)
    )

    cat(sprintf(
      "  nodos=%d  rango=%.2f km  sigma=%.3f  rho=%s  RMSE=%.4f  Cov95=%.1f%%  (%.1f min)\n\n",
      mesh$n, rango$mean, sqrt(varia$mean),
      if (is.na(rho_ar1)) "-" else sprintf("%.3f", rho_ar1), rmse, cov95, minutos
    ))
  }
}

# ==============================================================================
# MODELO M2: SOLO COVARIABLES, SIN CAMPO ESPACIAL
# ==============================================================================
# Se ajusta una sola vez porque no depende de la malla. Usa exactamente las
# mismas filas, respuesta, covariables y grupos LOSO que M1.
datos_m2 <- cbind(
  data.frame(y_response = df$y, Intercept = 1),
  as.data.frame(df[, COVARIABLES, with = FALSE])
)

m2 <- inla(
  y_response ~ -1 + Intercept +
    intensidad + Temperatura + Velocidad_Viento,
  data = datos_m2,
  family = "gaussian",
  control.predictor = list(compute = TRUE),
  control.compute = list(
    cpo = TRUE, dic = TRUE, waic = TRUE,
    return.marginals.predictor = TRUE
  ),
  control.inla = list(int.strategy = "eb")
)

gcv_m2 <- inla.group.cv(result = m2, groups = grupos)
err_m2 <- gcv_m2$mean - df$y
rmse_m2 <- sqrt(mean(err_m2^2, na.rm = TRUE))
mae_m2 <- mean(abs(err_m2), na.rm = TRUE)
cov95_m2 <- 100 * mean(
  gcv_m2$mean - 1.96 * gcv_m2$sd <= df$y &
    df$y <= gcv_m2$mean + 1.96 * gcv_m2$sd,
  na.rm = TRUE
)

resultado_m2 <- data.table(
  modelo = "M2: Covariates only",
  DIC = m2$dic$dic,
  WAIC = m2$waic$waic,
  LCPO = mean(-log(m2$cpo$cpo), na.rm = TRUE),
  PIT_mean = mean(m2$cpo$pit, na.rm = TRUE),
  PIT_sd = sd(m2$cpo$pit, na.rm = TRUE),
  RMSE_loso = rmse_m2,
  MAE_loso = mae_m2,
  Cov95_loso = cov95_m2
)


# ==============================================================================
# TABLA FINAL
# ==============================================================================
tabla <- rbindlist(resultados)
tabla[, malla := factor(malla, levels = MALLAS)]
tabla[, tipo := factor(tipo, levels = TIPOS)]
setorder(tabla, tipo, malla)

# Estabilidad: cuanto cambian rango y RMSE al pasar a la malla siguiente, DENTRO
# de cada tipo de modelo (by = tipo; si no, la primera malla del segundo tipo se
# compararia con la ultima del primero, que no tiene sentido).
# Si apenas se mueven, la malla mas gruesa de las dos ya vale.
tabla[, d_rango_pct := round(100 * (rango_km - shift(rango_km)) / shift(rango_km), 1),
  by = tipo
]
tabla[, d_RMSE_pct := round(100 * (RMSE_loso - shift(RMSE_loso)) / shift(RMSE_loso), 1),
  by = tipo
]

cat("--- Resultados ---\n")
print(tabla)

# Aviso de degeneracion: un rho pegado a 1 significa que el AR1 ha dejado de ser
# estacionario y se ha convertido en un paseo aleatorio que absorbe los datos.
# Cuando pasa, el DIC se hunde pero el error fuera de muestra empeora, asi que
# ese modelo NO puede presentarse como "el mejor por DIC".
degenerados <- tabla[!is.na(rho_ar1) & rho_ar1 > 0.99]
if (nrow(degenerados) > 0) {
  cat(
    "\n*** AVISO: rho del AR1 pegado a 1 en:",
    paste(degenerados$malla, collapse = ", "),
    "\n    El campo espacio-temporal ha degenerado. No usar su DIC.\n"
  )
}

tabla_est <- dcast(rbindlist(por_estacion), ESTACION ~ tipo + malla, value.var = "RMSE")
# Se ordena por el RMSE medio de la estacion en todas las configuraciones: deja
# arriba las estaciones que le cuestan al campo espacial, sea cual sea la malla.
tabla_est[, RMSE_medio := rowMeans(.SD), .SDcols = setdiff(names(tabla_est), "ESTACION")]
setorder(tabla_est, -RMSE_medio)
cat("\n--- RMSE leave-one-station-out por estacion ---\n")
print(tabla_est)

fwrite(tabla, file.path(DIR_OUT, sprintf("comparacion_mallas_%d_%s.csv", ANIO, ESCALA)))
fwrite(tabla_est, file.path(DIR_OUT, sprintf("rmse_por_estacion_%d_%s.csv", ANIO, ESCALA)))
cat("\nGuardado en:", DIR_OUT, "\n")

# ==============================================================================
# TABLAS ACADEMICAS PARA EL TFM
# ==============================================================================
source(here("R", "utilities", "academic_quality_tables.R"))

DIR_FIG_MALLAS <- here("outputs", "figures", "modelo", "mallas")
dir.create(DIR_FIG_MALLAS, recursive = TRUE, showWarnings = FALSE)

# La primera tabla responde a la pregunta de la malla: ademas de los valores
# absolutos muestra cuanto cambian el error y el rango al refinarla.
tabla_mallas_tfm <- copy(tabla)
tabla_mallas_tfm[, Mesh := tools::toTitleCase(as.character(malla))]
tabla_mallas_tfm[, Vertices := n_nodos]
tabla_mallas_tfm[, DIC := sprintf("%.1f", DIC)]
tabla_mallas_tfm[, WAIC := sprintf("%.1f", WAIC)]
tabla_mallas_tfm[, RMSE := sprintf("%.3f", RMSE_loso)]
tabla_mallas_tfm[, `RMSE change (%)` := fifelse(
  is.na(d_RMSE_pct), "--", sprintf("%+.1f", d_RMSE_pct)
)]
tabla_mallas_tfm[, `Cov95 (%)` := sprintf("%.1f", Cov95_loso)]
tabla_mallas_tfm[, `Range (km)` := sprintf(
  "%.2f [%.2f, %.2f]", rango_km, rango_q025, rango_q975
)]
tabla_mallas_tfm[, `Range change (%)` := fifelse(
  is.na(d_rango_pct), "--", sprintf("%+.1f", d_rango_pct)
)]
tabla_mallas_tfm[, `Spatial SD` := sprintf("%.3f", sigma)]
tabla_mallas_tfm <- tabla_mallas_tfm[, .(
  Mesh, Vertices, DIC, WAIC, RMSE, `RMSE change (%)`, `Cov95 (%)`,
  `Range (km)`, `Range change (%)`, `Spatial SD`
)]

ruta_tabla_mallas_tfm <- file.path(
  DIR_FIG_MALLAS,
  sprintf("tabla_comparacion_mallas_%d_%s.png", ANIO, ESCALA)
)
booktabs_png(
  tabla_mallas_tfm,
  ruta_tabla_mallas_tfm,
  title = "Mesh sensitivity analysis",
  subtitle = "M1: Gaussian model with fixed covariates and an SPDE spatial field",
  note = paste0(
    "RMSE and Cov95 are based on leave-one-station-out prediction. Range is ",
    "posterior mean [95% credible interval]. Changes are relative to the ",
    "previous, coarser mesh; a small range change indicates mesh stability."
  ),
  widths = c(0.85, 0.65, 0.68, 0.68, 0.68, 1.18, 0.80, 1.55, 1.18, 0.85),
  align = c("left", rep("right", 9)),
  font_size = 8.5,
  row_height = 0.28
)
cat("Tabla academica de mallas guardada en:", ruta_tabla_mallas_tfm, "\n")

# La segunda tabla responde a la pregunta del efecto espacial. M2 aparece una
# sola vez porque no usa malla; M1 aparece en las tres resoluciones para mostrar
# si su ventaja frente a M2 es robusta a la discretizacion.
tabla_modelos <- rbindlist(list(
  data.table(
    Model = resultado_m2$modelo,
    Mesh = "--",
    DIC = resultado_m2$DIC,
    WAIC = resultado_m2$WAIC,
    RMSE = resultado_m2$RMSE_loso,
    Cov95 = resultado_m2$Cov95_loso
  ),
  tabla[, .(
    Model = "M1: Covariates + spatial field",
    Mesh = tools::toTitleCase(as.character(malla)),
    DIC,
    WAIC,
    RMSE = RMSE_loso,
    Cov95 = Cov95_loso
  )]
), use.names = TRUE)

tabla_modelos[, `DIC difference` := DIC - resultado_m2$DIC]
tabla_modelos[, `WAIC difference` := WAIC - resultado_m2$WAIC]
tabla_modelos[, `RMSE gain (%)` := 100 * (resultado_m2$RMSE_loso - RMSE) /
  resultado_m2$RMSE_loso]

fwrite(
  tabla_modelos,
  file.path(DIR_OUT, sprintf("comparacion_modelos_%d_%s.csv", ANIO, ESCALA))
)

tabla_modelos_tfm <- copy(tabla_modelos)
tabla_modelos_tfm[, DIC := sprintf("%.1f", DIC)]
tabla_modelos_tfm[, WAIC := sprintf("%.1f", WAIC)]
tabla_modelos_tfm[, RMSE := sprintf("%.3f", RMSE)]
tabla_modelos_tfm[, `Cov95 (%)` := sprintf("%.1f", Cov95)]
tabla_modelos_tfm[, `DIC difference` := fifelse(
  Model == resultado_m2$modelo, "Reference", sprintf("%+.1f", `DIC difference`)
)]
tabla_modelos_tfm[, `WAIC difference` := fifelse(
  Model == resultado_m2$modelo, "Reference", sprintf("%+.1f", `WAIC difference`)
)]
tabla_modelos_tfm[, `RMSE gain (%)` := fifelse(
  Model == resultado_m2$modelo, "Reference", sprintf("%+.1f", `RMSE gain (%)`)
)]
tabla_modelos_tfm <- tabla_modelos_tfm[, .(
  Model, Mesh, DIC, `DIC difference`, WAIC, `WAIC difference`,
  RMSE, `RMSE gain (%)`, `Cov95 (%)`
)]

ruta_tabla_modelos_tfm <- file.path(
  DIR_FIG_MALLAS,
  sprintf("tabla_comparacion_modelos_%d_%s.png", ANIO, ESCALA)
)
booktabs_png(
  tabla_modelos_tfm,
  ruta_tabla_modelos_tfm,
  title = "Assessment of the spatial effect",
  subtitle = "M1 includes an SPDE spatial field; M2 contains the same fixed covariates only",
  note = paste0(
    "All models use the same observations and leave-one-station-out groups. ",
    "DIC and WAIC differences are M1 minus M2 (negative favours M1). RMSE gain ",
    "is relative to M2 (positive favours M1); Cov95 should be close to 95%."
  ),
  widths = c(2.15, 0.70, 0.68, 1.10, 0.68, 1.18, 0.70, 1.02, 0.82),
  align = c("left", "left", rep("right", 7)),
  group_starts = 2L,
  font_size = 8.5,
  row_height = 0.28
)
cat("Tabla academica de modelos guardada en:", ruta_tabla_modelos_tfm, "\n")
