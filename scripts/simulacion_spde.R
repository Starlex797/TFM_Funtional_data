# ==============================================================================
# SIMULACIÓN INLA-SPDE: TAMAÑO DE MALLA Y COMPARACIÓN DE MODELOS
# ==============================================================================
# Objetivo: Seleccionar 19 días, dividir 15 train / 4 test, y comparar:
#   - 3 tamaños de malla (gruesa, media, fina)
#   - 2 modelos: solo espacial vs espacio-temporal (AR1)
# Métrica principal: RMSE, MAE y Coverage 95% sobre los 4 días de test.
# ==============================================================================

library(INLA)
library(data.table)
library(sf)
library(here)
library(ggplot2)

set.seed(7291)

carpeta_sim <- here("outputs", "simulacion")
dir.create(carpeta_sim, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. CARGA DE DATOS Y SELECCIÓN DE 19 DÍAS CONSECUTIVOS DE NOVIEMBRE
# ==============================================================================

dt_maestro <- readRDS(here("data", "processed",
                           "dataset_maestro_inla_2025_DIARIO.rds"))
setDT(dt_maestro)
View()

fechas_disponibles <- sort(unique(dt_maestro$FECHA))

# Filtrar solo las fechas que corresponden a noviembre (mes "11")
fechas_noviembre <- fechas_disponibles[format(fechas_disponibles, "%m") == "11"]

# Verificar que haya al menos 19 días en noviembre
if(length(fechas_noviembre) < 19) {
  stop("No hay suficientes días en noviembre en tu dataset.")
}

# Seleccionar los primeros 19 días consecutivos de noviembre
fechas_sim <- fechas_noviembre[1:19]

cat("Días consecutivos seleccionados para la simulación:\n")
print(data.table(
  idx        = seq_along(fechas_sim),
  fecha      = fechas_sim,
  dia_sem    = weekdays(fechas_sim),
  conjunto   = ifelse(seq_along(fechas_sim) <= 15, "TRAIN", "TEST") # Wong usa 14 y 5, pero 15 y 4 está bien.
))

fechas_train <- fechas_sim[1:15]
fechas_test  <- fechas_sim[16:19]

# Filtrar y reindexar (el ID_TIEMPO_SIM debe ir del 1 al 19 secuencialmente)
dt_sim <- dt_maestro[FECHA %in% fechas_sim]
dt_sim[, ID_TIEMPO_SIM := match(FECHA, fechas_sim)]
setorder(dt_sim, ID_TIEMPO_SIM, ESTACION)

n_train <- 15L
n_test  <- 4L
ndays   <- 19L

cat(sprintf("\nFilas totales: %d | Estaciones: %d | Días: %d (train=%d, test=%d)\n",
            nrow(dt_sim), uniqueN(dt_sim$ESTACION), ndays, n_train, n_test))
# ==============================================================================
# 2. COORDENADAS UTM (km) Y BOUNDARIES
# ==============================================================================

coords_unicas <- unique(dt_sim[, .(ESTACION, LONGITUD, LATITUD)])
coords_sf <- st_as_sf(coords_unicas,
                       coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)

coords_unicas[, X_km := st_coordinates(coords_sf)[, 1] / 1000]
coords_unicas[, Y_km := st_coordinates(coords_sf)[, 2] / 1000]

dt_sim <- merge(dt_sim, coords_unicas[, .(ESTACION, X_km, Y_km)],
                by = "ESTACION", all.x = TRUE)
setorder(dt_sim, ID_TIEMPO_SIM, ESTACION)

coords_puntos <- as.matrix(dt_sim[, .(X_km, Y_km)])
coords_matriz <- as.matrix(unique(dt_sim[, .(X_km, Y_km)]))

# Boundaries para las mallas
bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.05, resolution = 50)
bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.2)

# ==============================================================================
# 3. DEFINIR COVARIABLES Y VECTOR RESPUESTA (NA en test)
# ==============================================================================

# Respuesta: NA para los 4 días de test (INLA predice automáticamente)
y_full <- dt_sim$LOG_NO2_DIARIO
y_train_test <- ifelse(dt_sim$FECHA %in% fechas_train,
                       dt_sim$LOG_NO2_DIARIO, NA)

covariables <- list(
  trafico_intensidad  = dt_sim$intensidad,
  trafico_carga       = dt_sim$carga,
  temperatura         = dt_sim$Temperatura,
  humedad             = dt_sim$Humedad_Relativa,
  precipitacion       = dt_sim$Precipitaciones,
  presion_barometrica = dt_sim$`Presion Barométrica`,
  radiacion_solar     = dt_sim$`Radiación Solar`,
  velocidad_viento    = dt_sim$`Velocidad Viento`

)

# Índices de las filas de test para extraer predicciones
idx_test_filas <- which(dt_sim$FECHA %in% fechas_test)
y_test_real    <- dt_sim$LOG_NO2_DIARIO[idx_test_filas]

# ==============================================================================
# 4. DEFINIR LAS 3 MALLAS A COMPARAR
# ==============================================================================

config_mallas <- list(
  gruesa = list(max.edge = c(8, 12), cutoff = 0.5,  label = "Gruesa (8 km)"),
  media  = list(max.edge = c(4, 8),  cutoff = 0.5,  label = "Media (4 km)"),
  fina   = list(max.edge = c(1, 4),  cutoff = 0.25, label = "Fina (1 km)")
)

# ==============================================================================
# 5. FÓRMULAS (compartidas por todas las iteraciones)
# ==============================================================================

efectos_fijos <- y ~ 0 + intercept +
  trafico_intensidad + temperatura + precipitacion + velocidad_viento

formula_st <- update(efectos_fijos,
  . ~ . + f(campo_espacial,
            model = spde,
            group = campo_espacial.group,
            control.group = list(model = "ar1")))

formula_s <- update(efectos_fijos,
  . ~ . + f(campo_espacial_s, model = spde))

# ==============================================================================
# 6. BUCLE PRINCIPAL: 3 MALLAS × 2 MODELOS
# ==============================================================================

resultados <- list()

for (malla_nombre in names(config_mallas)) {

  cfg <- config_mallas[[malla_nombre]]
  cat("\n", strrep("=", 70), "\n")
  cat(" MALLA:", cfg$label, "\n")
  cat(strrep("=", 70), "\n")

  # --- Crear malla ---
  malla <- inla.mesh.2d(
    loc      = coords_matriz,
    boundary = list(bnd_inner, bnd_outer),
    max.edge = cfg$max.edge,
    cutoff   = cfg$cutoff
  )
  cat(sprintf("  Nodos: %d\n", malla$n))

  # --- SPDE ---
  spde <- inla.spde2.matern(mesh = malla, alpha = 2)

  # =========================================================================
  # 6a. MODELO ESPACIO-TEMPORAL (AR1)
  # =========================================================================
  cat("\n  [ST] Ajustando modelo espacio-temporal...\n")

  indice_st <- inla.spde.make.index(
    name = "campo_espacial", n.spde = spde$n.spde, n.group = ndays
  )
  A_st <- inla.spde.make.A(
    mesh = malla, loc = coords_puntos,
    group = dt_sim$ID_TIEMPO_SIM, n.group = ndays
  )

  stk_st <- inla.stack(
    tag     = "sim_st",
    data    = list(y = y_train_test),
    A       = list(A_st, 1),
    effects = list(
      c(indice_st, list(intercept = 1)),
      covariables
    ),
    compress = FALSE
  )

  t0 <- Sys.time()
  modelo_st <- tryCatch(
    inla(
      formula           = formula_st,
      data              = inla.stack.data(stk_st, spde = spde),
      family            = "gaussian",
      control.predictor = list(A = inla.stack.A(stk_st), compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE),
      control.inla      = list(strategy = "laplace"),
      verbose           = FALSE
    ),
    error = function(e) { message("  ERROR ST: ", e$message); NULL }
  )
  t_st <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (!is.null(modelo_st)) {
    idx_st_data  <- inla.stack.index(stk_st, tag = "sim_st")$data
    pred_st_mean <- modelo_st$summary.fitted.values$mean[idx_st_data][idx_test_filas]
    pred_st_sd   <- modelo_st$summary.fitted.values$sd[idx_st_data][idx_test_filas]

    # Posición relativa de test dentro del stack
    # idx_test_filas son posiciones en dt_sim → posiciones en el stack data
    pred_st_mean <- modelo_st$summary.fitted.values$mean[idx_st_data]
    pred_st_sd   <- modelo_st$summary.fitted.values$sd[idx_st_data]

    pred_st_mean_test <- pred_st_mean[idx_test_filas]
    pred_st_sd_test   <- pred_st_sd[idx_test_filas]

    rmse_st <- sqrt(mean((pred_st_mean_test - y_test_real)^2, na.rm = TRUE))
    mae_st  <- mean(abs(pred_st_mean_test - y_test_real), na.rm = TRUE)

    # Coverage 95%
    lo95 <- pred_st_mean_test - 1.96 * pred_st_sd_test
    hi95 <- pred_st_mean_test + 1.96 * pred_st_sd_test
    cov95_st <- mean(y_test_real >= lo95 & y_test_real <= hi95, na.rm = TRUE)

    cat(sprintf("  [ST] RMSE=%.4f | MAE=%.4f | Cov95=%.1f%% | %.1f min\n",
                rmse_st, mae_st, cov95_st * 100, t_st))

    resultados[[length(resultados) + 1]] <- data.table(
      Malla    = cfg$label,
      Modelo   = "Espacio-temporal (AR1)",
      n_nodos  = malla$n,
      DIC      = modelo_st$dic$dic,
      WAIC     = modelo_st$waic$waic,
      RMSE     = rmse_st,
      MAE      = mae_st,
      Cov95    = cov95_st,
      Tiempo_min = round(t_st, 2)
    )
  }

  # =========================================================================
  # 6b. MODELO SOLO ESPACIAL
  # =========================================================================
  cat("\n  [S]  Ajustando modelo solo espacial...\n")

  indice_s_solo <- inla.spde.make.index(
    name = "campo_espacial_s", n.spde = spde$n.spde
  )
  A_s <- inla.spde.make.A(mesh = malla, loc = coords_puntos)

  stk_s <- inla.stack(
    tag     = "sim_s",
    data    = list(y = y_train_test),
    A       = list(A_s, 1),
    effects = list(
      c(indice_s_solo, list(intercept = 1)),
      covariables
    ),
    compress = FALSE
  )

  t0 <- Sys.time()
  modelo_s <- tryCatch(
    inla(
      formula           = formula_s,
      data              = inla.stack.data(stk_s, spde = spde),
      family            = "gaussian",
      control.predictor = list(A = inla.stack.A(stk_s), compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE),
      control.inla      = list(strategy = "laplace"),
      verbose           = FALSE
    ),
    error = function(e) { message("  ERROR S: ", e$message); NULL }
  )
  t_s <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (!is.null(modelo_s)) {
    idx_s_data <- inla.stack.index(stk_s, tag = "sim_s")$data
    pred_s_mean <- modelo_s$summary.fitted.values$mean[idx_s_data]
    pred_s_sd   <- modelo_s$summary.fitted.values$sd[idx_s_data]

    pred_s_mean_test <- pred_s_mean[idx_test_filas]
    pred_s_sd_test   <- pred_s_sd[idx_test_filas]

    rmse_s <- sqrt(mean((pred_s_mean_test - y_test_real)^2, na.rm = TRUE))
    mae_s  <- mean(abs(pred_s_mean_test - y_test_real), na.rm = TRUE)

    lo95 <- pred_s_mean_test - 1.96 * pred_s_sd_test
    hi95 <- pred_s_mean_test + 1.96 * pred_s_sd_test
    cov95_s <- mean(y_test_real >= lo95 & y_test_real <= hi95, na.rm = TRUE)

    cat(sprintf("  [S]  RMSE=%.4f | MAE=%.4f | Cov95=%.1f%% | %.1f min\n",
                rmse_s, mae_s, cov95_s * 100, t_s))

    resultados[[length(resultados) + 1]] <- data.table(
      Malla    = cfg$label,
      Modelo   = "Solo espacial",
      n_nodos  = malla$n,
      DIC      = modelo_s$dic$dic,
      WAIC     = modelo_s$waic$waic,
      RMSE     = rmse_s,
      MAE      = mae_s,
      Cov95    = cov95_s,
      Tiempo_min = round(t_s, 2)
    )
  }

  rm(malla, spde, modelo_st, modelo_s); gc()
}

# ==============================================================================
# 7. TABLA COMPARATIVA FINAL
# ==============================================================================

tabla_comparativa <- rbindlist(resultados)
setorder(tabla_comparativa, RMSE)

cat("\n", strrep("=", 80), "\n")
cat("   RESULTADOS DE LA SIMULACIÓN: 3 MALLAS × 2 MODELOS\n")
cat(strrep("=", 80), "\n")
cat(sprintf("   Train: %d días | Test: %d días | Estaciones: %d\n",
            n_train, n_test, uniqueN(dt_sim$ESTACION)))
cat(strrep("-", 80), "\n")
print(tabla_comparativa, row.names = FALSE)

write.csv(tabla_comparativa,
          file.path(carpeta_sim, "comparacion_mallas_modelos.csv"),
          row.names = FALSE)

# ==============================================================================
# 8. VISUALIZACIÓN
# ==============================================================================

tabla_comparativa[, Malla := factor(Malla,
  levels = c("Gruesa (8 km)", "Media (4 km)", "Fina (1 km)"))]

ggplot(tabla_comparativa, aes(x = Malla, y = RMSE, fill = Modelo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(RMSE, 4)),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Espacio-temporal (AR1)" = "#2166AC",
                                "Solo espacial"          = "#B2182B")) +
  labs(
    title    = "Comparación de Modelos INLA-SPDE por Tamaño de Malla",
    subtitle = sprintf("RMSE sobre %d días de test | %d días de entrenamiento | Madrid 2025",
                       n_test, n_train),
    x = "Resolución de malla", y = "RMSE (log NO₂)", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "top"
  )

ggsave(file.path(carpeta_sim, "comparacion_rmse_malla_modelo.png"),
       width = 10, height = 6, dpi = 300)

  # Gráfico de tiempo de cómputo
ggplot(tabla_comparativa, aes(x = Malla, y = Tiempo_min, fill = Modelo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = paste0(Tiempo_min, " min")),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Espacio-temporal (AR1)" = "#2166AC",
                                "Solo espacial"          = "#B2182B")) +
  labs(
    title = "Tiempo de Cómputo por Malla y Modelo",
    x = "Resolución de malla", y = "Tiempo (minutos)", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

ggsave(file.path(carpeta_sim, "comparacion_tiempo_malla_modelo.png"),
       width = 10, height = 6, dpi = 300)

cat("\n✅ Simulación completada. Resultados en:", carpeta_sim, "\n")

