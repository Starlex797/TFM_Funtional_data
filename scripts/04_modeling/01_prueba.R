# ==============================================================================
# Estudio 2019-2021
# ==============================================================================

library(INLA)
library(data.table)
library(here)

# ------------------------------------------------------------------------------
# CONFIGURACION
# ------------------------------------------------------------------------------

COVARIABLES <- c("intensidad", "Temperatura", "Velocidad_Viento")

# El analisis de sensibilidad ajusta el mismo modelo espacial en cada malla.
# La comparacion M1-M2 se hace despues contra un unico modelo sin campo espacial.
TIPOS <- c("espacial")

DIR_OUT <- here("outputs", "modelo", "Seleccion_variables" )
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PASO 1: DATOS
# ==============================================================================
ANIOS <- 2019:2021

df <- rbindlist(
    lapply(ANIOS, function(anio) {
        ruta <- here(
            "data", "processed", "Maestro", as.character(anio),
            sprintf("dataset_maestro_inla_%d_DIARIO.rds", anio)
        )
        as.data.table(readRDS(ruta))
    }),
    use.names = TRUE,
    fill = TRUE
)

df[, FECHA := as.Date(FECHA)]

# La respuesta ya viene transformada: LOG_NO2_DIARIO es log1p(DATO_DIARIO).
estaciones <- sort(unique(df$ESTACION))
filas_por_estacion <- split(seq_len(nrow(df)), df$ESTACION)
grupos <- filas_por_estacion[as.character(df$ESTACION)]

df[, y := LOG_NO2_DIARIO]
df[, NO2 := DATO_DIARIO]

# ==============================================================================
# PASO 2-8: UN MODELO POR (MALLA, TIPO)
# ==============================================================================
    
mesh <- readRDS("malla_spde_madrid_media.rds")

spde <- inla.spde2.pcmatern(
      mesh = mesh,
      prior.range = c(9.3, 0.5), # P(rango < 9.3 km) = 0.5
      prior.sigma = c(0.6, 0.2)  # P(sigma > 0.5) = 0.5
    )

# --- PASOS 3, 4 y 6: indice, matriz A y formula -----------------------------
   

# Un unico campo espacial compartido por los 60 dias.

s.index <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)# 
A.estud <- inla.spde.make.A(mesh = mesh, loc =  as.matrix(df[, .(X_km, Y_km)]) )# UTM 30N en km, ya vienen del Paso 2
 
# Una replica del campo espacial por dia (n.group = n_dias), encadenadas
# por un AR1. Incognitas del campo: nodos x dias.


# --- PASO 5: inla.stack ----------------------------------------------------

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
m0 <- inla(y_response ~ -1 + Intercept + f(Temperatura, model = "rw2", scale.model = TRUE) + f(Velocidad_Viento, model = "rw2", scale.model = TRUE) +
    Presion_Barometrica + Humedad_Relativa + Radiacion_Solar + Precipitacion+ f(spatial.field, model = spde),
      data = inla.stack.data(stack.est, spde = spde),
      family = "gaussian", # la respuesta ya esta en log: gaussiana, no gamma
      control.predictor = list(A = inla.stack.A(stack.est), compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, return.marginals.predictor = TRUE, openmp.strategy = "huge"),
      control.inla = list( int.strategy = "eb"), # Empirical Bayes: mas rapido
      num.threads = 5

    )
minutos <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# --- PASO 7: rango y varianza en escala fisica (km) ------------------------
# INLA estima el campo en escala logaritmica interna (Theta1, Theta2).
# inla.spde2.result lo devuelve a rango (km) y varianza marginal.
spde.result <- inla.spde2.result(
inla = m0, name = "spatial.field",
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
gcv <- inla.group.cv(result = m0, groups = grupos)

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
lcpo <- mean(-log(m0$cpo$cpo), na.rm = TRUE)
pit_mean <- mean(m0$cpo$pit, na.rm = TRUE)
pit_sd <- sd(m0$cpo$pit, na.rm = TRUE)

   
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
      DIC        = m0$dic$dic,
      WAIC       = m0$waic$waic,
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