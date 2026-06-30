# ==============================================================================
# SIMULACIÓN INLA-SPDE: MÚLTIPLES FRECUENCIAS × MALLAS × MODELOS
# ==============================================================================
# Objetivo:
#   - Comparar 3 frecuencias temporales: HORARIA, DIARIA, MENSUAL
#   - Para cada frecuencia:
#       • Selección de covariables por VIF (sin stepwise)
#       • Horario y Mensual: solo modelo espacial (3 mallas)
#       • Diario: modelo espacial + espacio-temporal AR1 (3 mallas)
#   - Guardar resultados como tablas PNG y gráficos
# ==============================================================================

library(INLA)
library(fmesher)
library(data.table)
library(sf)
library(here)
library(ggplot2)
library(car)
library(gt)

set.seed(4827)

# ==============================================================================
# CONFIGURACIÓN GENERAL
# ==============================================================================

VIF_UMBRAL     <- 5
DIC_MEJORA_MIN <- 2

config_mallas <- list(
  gruesa = list(max.edge = c(8, 12), cutoff = 0.5,  label = "Gruesa (8 km)"),
  media  = list(max.edge = c(4,  8), cutoff = 0.5,  label = "Media (4 km)"),
  fina   = list(max.edge = c(1,  4), cutoff = 0.25, label = "Fina (1 km)")
)

# Covariables candidatas (nombres estandarizados presentes en los 3 datasets)
COVS_NOMBRES <- c("intensidad", "Temperatura",
                   "Precipitaciones", "Presion Barométrica",
                   "Velocidad Viento")
COVS_ALIAS   <- c("trafico_intensidad", "temperatura",
                   "precipitacion", "presion_barometrica",
                   "velocidad_viento")

# Columnas raw para re-estandarizar en mensual
COLS_RAW <- c("intensidad_raw", "Temperatura_raw",
              "Precipitaciones_raw", "Presion Barométrica_raw",
              "Velocidad Viento_raw")

# ==============================================================================
# FUNCIONES AUXILIARES
# ==============================================================================

# --- Guardar tabla gt como PNG ---
guardar_tabla_png <- function(df, titulo, subtitulo = NULL, ruta_png,
                              ancho_px = 800) {
  tbl <- gt(df) |>
    tab_header(title = titulo, subtitle = subtitulo) |>
    tab_options(
      table.font.size      = px(12),
      heading.title.font.size = px(16),
      heading.subtitle.font.size = px(12),
      column_labels.font.weight = "bold",
      table.border.top.color = "black",
      table.border.bottom.color = "black",
      heading.border.bottom.color = "black",
      column_labels.border.bottom.color = "black"
    ) |>
    opt_horizontal_padding(scale = 2)

  gtsave(tbl, filename = ruta_png, vwidth = ancho_px)
  cat("  Tabla guardada:", basename(ruta_png), "\n")
}

# --- Preparar datos según frecuencia y periodo estacional ---
preparar_datos <- function(frecuencia, periodo = "noviembre") {

  if (frecuencia == "horario") {
    dt <- readRDS(here("data", "processed",
                       "dataset_maestro_inla_2025_HORARIO.rds"))
    setDT(dt)
    setnames(dt, "LOG_NO2_HORARIO", "LOG_NO2")

    todas_fechas <- sort(unique(dt$FECHA))

    if (periodo == "noviembre") {
      fechas_mes <- todas_fechas[format(todas_fechas, "%m") == "11"]
      fechas_sel <- fechas_mes[13:17]
    } else if (periodo == "marzo") {
      fechas_mes <- todas_fechas[format(todas_fechas, "%m") == "03"]
      fechas_sel <- fechas_mes[7:11]
    }

    dt <- dt[FECHA %in% fechas_sel]
    dt[, ID_TIEMPO_SIM := as.integer(factor(
      paste(FECHA, sprintf("%02d", HORA))))]
    setorder(dt, ID_TIEMPO_SIM, ESTACION)
    n_periodos <- uniqueN(dt$ID_TIEMPO_SIM)
    n_test_p   <- 2L * 24L
    n_train_p  <- n_periodos - n_test_p
    periodos   <- sort(unique(dt$ID_TIEMPO_SIM))
    periodos_train <- periodos[seq_len(n_train_p)]
    periodos_test  <- periodos[(n_train_p + 1L):n_periodos]
    dt[, es_train := ID_TIEMPO_SIM %in% periodos_train]
    lab_train <- sprintf("%s → %s (%d horas)",
                         min(fechas_sel), fechas_sel[3], n_train_p)
    lab_test  <- sprintf("%s → %s (%d horas)",
                         fechas_sel[4], max(fechas_sel), n_test_p)

  } else if (frecuencia == "diario") {
    dt <- readRDS(here("data", "processed",
                       "dataset_maestro_inla_2025_DIARIO.rds"))
    setDT(dt)
    setnames(dt, "LOG_NO2_DIARIO", "LOG_NO2")

    todas_fechas <- sort(unique(dt$FECHA))

    if (periodo == "noviembre") {
      fechas_mes <- todas_fechas[format(todas_fechas, "%m") == "11"]
      fechas_sel <- fechas_mes[1:19]
    } else if (periodo == "marzo") {
      fechas_mes <- todas_fechas[format(todas_fechas, "%m") == "03"]
      fechas_sel <- fechas_mes[3:21]
    }

    dt <- dt[FECHA %in% fechas_sel]
    dt[, ID_TIEMPO_SIM := match(FECHA, fechas_sel)]
    setorder(dt, ID_TIEMPO_SIM, ESTACION)
    n_periodos <- 19L
    n_test_p   <- 4L
    n_train_p  <- 15L
    periodos_train <- 1L:n_train_p
    periodos_test  <- (n_train_p + 1L):n_periodos
    dt[, es_train := ID_TIEMPO_SIM %in% periodos_train]
    lab_train <- sprintf("%s → %s (%d días)",
                         fechas_sel[1], fechas_sel[15], n_train_p)
    lab_test  <- sprintf("%s → %s (%d días)",
                         fechas_sel[16], fechas_sel[19], n_test_p)

  } else if (frecuencia == "mensual") {
    dt_diario <- readRDS(here("data", "processed",
                              "dataset_maestro_inla_2025_DIARIO.rds"))
    setDT(dt_diario)
    dt_diario[, ANIO_MES := as.Date(format(FECHA, "%Y-%m-01"))]

    dt <- dt_diario[, .(
      DATO_NO2                  = mean(DATO_DIARIO,             na.rm = TRUE),
      intensidad_raw            = mean(intensidad_raw,          na.rm = TRUE),
      Temperatura_raw           = mean(Temperatura_raw,         na.rm = TRUE),
      Precipitaciones_raw       = sum(Precipitaciones_raw,      na.rm = TRUE),
      `Presion Barométrica_raw` = mean(`Presion Barométrica_raw`, na.rm = TRUE),
      `Velocidad Viento_raw`    = mean(`Velocidad Viento_raw`,  na.rm = TRUE),
      n_dias                    = .N
    ), by = .(ESTACION, ANIO_MES, barrio, distrito, LONGITUD, LATITUD,
              ID_DISTRITO)]

    setnames(dt, "ANIO_MES", "FECHA")
    dt[, LOG_NO2 := log(DATO_NO2)]

    # Re-estandarizar a escala mensual
    for (i in seq_along(COLS_RAW)) {
      dt[, (COVS_NOMBRES[i]) := scale(get(COLS_RAW[i]))[, 1]]
    }

    fechas_disp <- sort(unique(dt$FECHA))
    n_periodos  <- length(fechas_disp)
    n_test_p    <- 2L
    n_train_p   <- n_periodos - n_test_p
    dt[, ID_TIEMPO_SIM := match(FECHA, fechas_disp)]
    setorder(dt, ID_TIEMPO_SIM, ESTACION)
    periodos_train <- seq_len(n_train_p)
    periodos_test  <- (n_train_p + 1L):n_periodos
    dt[, es_train := ID_TIEMPO_SIM %in% periodos_train]
    lab_train <- sprintf("%s → %s (%d meses)",
                         fechas_disp[1], fechas_disp[n_train_p], n_train_p)
    lab_test  <- sprintf("%s → %s (%d meses)",
                         fechas_disp[n_train_p + 1], fechas_disp[n_periodos],
                         n_test_p)
    rm(dt_diario)
  }

  # Coordenadas UTM en km
  coords_unicas <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
  coords_sf <- st_as_sf(coords_unicas,
                         coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
    st_transform(25830)
  coords_unicas[, X_km := st_coordinates(coords_sf)[, 1] / 1000]
  coords_unicas[, Y_km := st_coordinates(coords_sf)[, 2] / 1000]
  dt <- merge(dt, coords_unicas[, .(ESTACION, X_km, Y_km)],
              by = "ESTACION", all.x = TRUE)
  setorder(dt, ID_TIEMPO_SIM, ESTACION)

  list(
    dt           = dt,
    n_periodos   = uniqueN(dt$ID_TIEMPO_SIM),
    n_train      = n_train_p,
    n_test       = n_test_p,
    lab_train    = lab_train,
    lab_test     = lab_test
  )
}

# --- Heatmap de correlación (datos no estandarizados) ---
generar_heatmap_correlacion <- function(dt, carpeta) {

  # Usar todas las covariables raw disponibles (incluyendo las excluidas)
  cols_raw_todas <- c("intensidad_raw", "carga_raw", "Temperatura_raw",
                      "Humedad_Relativa_raw", "Precipitaciones_raw",
                      "Presion Barométrica_raw", "Radiación Solar_raw",
                      "Velocidad Viento_raw")
  cols_raw_disp <- intersect(cols_raw_todas, names(dt))

  etiquetas <- c(
    intensidad_raw              = "Intensidad tráfico",
    carga_raw                   = "Carga tráfico",
    Temperatura_raw             = "Temperatura",
    Humedad_Relativa_raw        = "Humedad relativa",
    Precipitaciones_raw         = "Precipitaciones",
    `Presion Barométrica_raw`   = "Presión barométrica",
    `Radiación Solar_raw`       = "Radiación solar",
    `Velocidad Viento_raw`      = "Velocidad viento"
  )

  dt_train <- dt[es_train == TRUE]
  mat <- dt_train[, c(..cols_raw_disp)]
  mat[, LOG_NO2 := dt_train$LOG_NO2]

  nombres_plot <- c(etiquetas[cols_raw_disp], "log NO2")
  setnames(mat, names(mat), nombres_plot)

  cor_mat <- cor(mat, use = "pairwise.complete.obs")

  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Var1", "Var2", "Correlacion")

  niveles <- c(setdiff(nombres_plot, "log NO2"), "log NO2")
  cor_df$Var1 <- factor(cor_df$Var1, levels = niveles)
  cor_df$Var2 <- factor(cor_df$Var2, levels = rev(niveles))

  ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlacion)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", Correlacion)),
              size = 3, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1),
                         name = "r") +
    labs(title = "Correlación entre covariables y log NO2",
         subtitle = "Datos no estandarizados (train)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title   = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(file.path(carpeta, "heatmap_correlacion_raw.png"),
         width = 9, height = 7, dpi = 300)
  cat("  Heatmap de correlación guardado: heatmap_correlacion_raw.png\n")
}

# --- Selección de covariables: solo VIF hacia atrás ---
seleccionar_covariables <- function(datos, carpeta) {

  dt <- datos$dt

  # 0. Heatmap de correlación (con todas las variables, antes de excluir)
  generar_heatmap_correlacion(dt, carpeta)

  # Construir lista de covariables candidatas (ya sin carga, humedad, radiación)
  covs_disponibles <- intersect(COVS_NOMBRES, names(dt))
  alias_map <- setNames(COVS_ALIAS[match(covs_disponibles, COVS_NOMBRES)],
                        covs_disponibles)

  covariables_candidatas <- setNames(
    lapply(covs_disponibles, function(nm) dt[[nm]]),
    alias_map[covs_disponibles]
  )

  # Pre-filtro: descartar covariables con varianza ≈ 0 en train
  vars_pool <- names(covariables_candidatas)
  for (v in vars_pool) {
    sd_v <- sd(covariables_candidatas[[v]][dt$es_train], na.rm = TRUE)
    if (is.na(sd_v) || sd_v < 1e-10) {
      cat(sprintf("    Descartada '%s' (sd ≈ 0 en train)\n", v))
      vars_pool <- setdiff(vars_pool, v)
    }
  }
  covariables_candidatas <- covariables_candidatas[vars_pool]

  # ---------- VIF hacia atrás ----------
  df_vif <- as.data.frame(
    lapply(covariables_candidatas, function(x) x[dt$es_train])
  )
  df_vif$LOG_NO2 <- dt$LOG_NO2[dt$es_train]
  df_vif <- na.omit(df_vif)
  vars_vif <- names(covariables_candidatas)

  cat("\n  --- VIF hacia atrás (umbral =", VIF_UMBRAL, ") ---\n")

  repeat {
    fml <- as.formula(paste("LOG_NO2 ~", paste(vars_vif, collapse = " + ")))
    lm_fit <- lm(fml, data = df_vif)
    vif_vals <- car::vif(lm_fit)

    cat(sprintf("    Variables: %s\n", paste(vars_vif, collapse = ", ")))

    vif_max <- max(vif_vals)
    if (vif_max <= VIF_UMBRAL) {
      cat(sprintf("    → Todas VIF ≤ %.1f. Fin.\n", VIF_UMBRAL))
      break
    }

    var_elim <- names(which.max(vif_vals))
    cat(sprintf("    → Eliminando '%s' (VIF = %.3f)\n", var_elim, vif_max))
    vars_vif <- setdiff(vars_vif, var_elim)

    if (length(vars_vif) < 2) {
      warning("VIF dejó menos de 2 variables.")
      break
    }
  }

  # Tabla VIF final
  tabla_vif <- data.frame(
    Variable = names(vif_vals),
    VIF      = round(as.numeric(vif_vals), 4),
    row.names = NULL
  )
  tabla_vif <- tabla_vif[order(tabla_vif$VIF, decreasing = TRUE), ]

  guardar_tabla_png(
    tabla_vif,
    titulo    = "Selección de Covariables — VIF",
    subtitulo = sprintf("Umbral VIF = %g | %d variables retenidas",
                        VIF_UMBRAL, length(vars_vif)),
    ruta_png  = file.path(carpeta, "tabla_seleccion_vif.png")
  )

  cat(sprintf("  Variables tras VIF: %s\n", paste(vars_vif, collapse = ", ")))

  # Preparar objetos comunes para significancia y stepwise
  coords_puntos <- as.matrix(dt[, .(X_km, Y_km)])
  coords_matriz <- as.matrix(unique(dt[, .(X_km, Y_km)]))
  bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.05, resolution = 50)
  bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.2)

  malla_ref <- inla.mesh.2d(
    loc      = coords_matriz,
    boundary = list(bnd_inner, bnd_outer),
    max.edge = c(4, 8), cutoff = 0.5
  )
  spde_ref <- inla.spde2.matern(mesh = malla_ref, alpha = 2)
  indice_ref <- inla.spde.make.index("campo_espacial",
                                      n.spde = spde_ref$n.spde)
  A_ref <- inla.spde.make.A(mesh = malla_ref, loc = coords_puntos)

  y_train_test   <- ifelse(dt$es_train, dt$LOG_NO2, NA_real_)
  idx_test_filas <- which(!dt$es_train)
  y_test_real    <- dt$LOG_NO2[idx_test_filas]

  # ---------- Significancia bayesiana (malla referencia) ----------
  formula_full <- as.formula(paste(
    "y ~ 0 + intercept +",
    paste(vars_vif, collapse = " + "),
    "+ f(campo_espacial, model = spde_ref)"
  ))

  stk_full <- inla.stack(
    tag = "full", data = list(y = y_train_test),
    A = list(A_ref, 1),
    effects = list(c(indice_ref, list(intercept = 1)),
                   covariables_candidatas[vars_vif]),
    compress = FALSE
  )

  cat("  Ajustando modelo completo (significancia)...\n")
  modelo_full <- inla(
    formula           = formula_full,
    data              = inla.stack.data(stk_full, spde = spde_ref),
    family            = "gaussian",
    control.predictor = list(A = inla.stack.A(stk_full), compute = TRUE),
    control.compute   = list(dic = TRUE, waic = FALSE, cpo = FALSE),
    control.inla      = list(strategy = "laplace"),
    verbose           = FALSE
  )

  sf_full <- modelo_full$summary.fixed
  sf_vars <- sf_full[rownames(sf_full) != "intercept", ]
  sig_mask <- sf_vars[, "0.025quant"] > 0 | sf_vars[, "0.975quant"] < 0
  vars_sig <- rownames(sf_vars)[sig_mask]

  tabla_sig <- data.frame(
    Variable  = rownames(sf_vars),
    Media     = round(sf_vars$mean, 4),
    SD        = round(sf_vars$sd, 4),
    `Q2.5`    = round(sf_vars[, "0.025quant"], 4),
    `Q97.5`   = round(sf_vars[, "0.975quant"], 4),
    Sig_95    = ifelse(sig_mask, "Sí", "No"),
    row.names = NULL, check.names = FALSE
  )

  guardar_tabla_png(
    tabla_sig,
    titulo    = "Significancia Bayesiana (IC 95%)",
    subtitulo = sprintf("%d significativas de %d",
                        sum(sig_mask), nrow(sf_vars)),
    ruta_png  = file.path(carpeta, "tabla_significancia_bayesiana.png")
  )

  cat(sprintf("  Significativas: %s\n", paste(vars_sig, collapse = ", ")))

  # ---------- Stepwise DIC ----------
  ajustar_sw <- function(vars_mod) {
    fml_sw <- as.formula(paste(
      "y ~ 0 + intercept +",
      paste(vars_mod, collapse = " + "),
      "+ f(campo_espacial, model = spde_ref)"
    ))
    stk_sw <- inla.stack(
      tag = "sw", data = list(y = y_train_test),
      A = list(A_ref, 1),
      effects = list(c(indice_ref, list(intercept = 1)),
                     covariables_candidatas[vars_mod]),
      compress = FALSE
    )
    mod <- tryCatch(
      inla(formula = fml_sw,
           data = inla.stack.data(stk_sw, spde = spde_ref),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(stk_sw), compute = TRUE),
           control.compute   = list(dic = TRUE, waic = FALSE, cpo = FALSE),
           control.inla      = list(strategy = "laplace"),
           verbose = FALSE),
      error = function(e) { message("  ERROR: ", e$message); NULL }
    )
    if (is.null(mod)) return(list(DIC = Inf, RMSE = NA_real_, MAE = NA_real_))
    idx_d <- inla.stack.index(stk_sw, tag = "sw")$data
    pred  <- mod$summary.fitted.values$mean[idx_d][idx_test_filas]
    list(DIC  = mod$dic$dic,
         RMSE = sqrt(mean((pred - y_test_real)^2, na.rm = TRUE)),
         MAE  = mean(abs(pred - y_test_real), na.rm = TRUE))
  }

  vars_actuales <- if (length(vars_sig) > 0) vars_sig else vars_vif
  tabla_sw <- list()
  res_actual <- ajustar_sw(vars_actuales)

  cat(sprintf("  Stepwise inicial: %s (DIC=%.2f)\n",
              paste(vars_actuales, collapse = " + "), res_actual$DIC))

  tabla_sw[[1]] <- data.frame(
    Iter = 0L, Accion = "inicial", Variable = "-",
    Variables = paste(vars_actuales, collapse = " + "),
    DIC  = round(res_actual$DIC, 2),
    RMSE = round(res_actual$RMSE, 4),
    MAE  = round(res_actual$MAE, 4),
    stringsAsFactors = FALSE
  )

  for (iter in seq_len(length(vars_vif))) {
    cands <- data.frame(accion = character(), variable = character(),
                        DIC = numeric(), stringsAsFactors = FALSE)

    for (v in setdiff(vars_vif, vars_actuales)) {
      r <- ajustar_sw(c(vars_actuales, v))
      cands <- rbind(cands, data.frame(accion = "añadir", variable = v,
                                       DIC = r$DIC, RMSE = r$RMSE,
                                       MAE = r$MAE, stringsAsFactors = FALSE))
    }
    if (length(vars_actuales) >= 2) {
      for (v in vars_actuales) {
        r <- ajustar_sw(setdiff(vars_actuales, v))
        cands <- rbind(cands, data.frame(accion = "eliminar", variable = v,
                                         DIC = r$DIC, RMSE = r$RMSE,
                                         MAE = r$MAE, stringsAsFactors = FALSE))
      }
    }

    mejor  <- cands[which.min(cands$DIC), ]
    mejora <- res_actual$DIC - mejor$DIC

    if (mejora < DIC_MEJORA_MIN) {
      cat(sprintf("  Stepwise convergido (ΔDIC=%.2f)\n", mejora))
      break
    }

    vars_actuales <- if (mejor$accion == "añadir")
      c(vars_actuales, mejor$variable) else
      setdiff(vars_actuales, mejor$variable)

    res_actual <- list(DIC = mejor$DIC, RMSE = mejor$RMSE, MAE = mejor$MAE)

    cat(sprintf("  Iter %d: %s '%s' (ΔDIC=%.2f)\n",
                iter, mejor$accion, mejor$variable, mejora))

    tabla_sw[[iter + 1]] <- data.frame(
      Iter = iter, Accion = mejor$accion, Variable = mejor$variable,
      Variables = paste(vars_actuales, collapse = " + "),
      DIC  = round(mejor$DIC, 2),
      RMSE = round(mejor$RMSE, 4),
      MAE  = round(mejor$MAE, 4),
      stringsAsFactors = FALSE
    )
  }

  tabla_stepwise <- do.call(rbind, tabla_sw)

  guardar_tabla_png(
    tabla_stepwise,
    titulo    = "Stepwise DIC — Selección de Covariables",
    subtitulo = sprintf("Criterio: ΔDIC mínimo = %g | Variables finales: %s",
                        DIC_MEJORA_MIN,
                        paste(vars_actuales, collapse = ", ")),
    ruta_png  = file.path(carpeta, "tabla_stepwise_dic.png"),
    ancho_px  = 1000
  )

  cat(sprintf("  Variables finales: %s\n",
              paste(vars_actuales, collapse = ", ")))

  list(
    vars_finales           = vars_actuales,
    covariables_candidatas = covariables_candidatas,
    y_train_test           = y_train_test,
    idx_test_filas         = idx_test_filas,
    y_test_real            = y_test_real,
    coords_puntos          = coords_puntos,
    coords_matriz          = coords_matriz,
    bnd_inner              = bnd_inner,
    bnd_outer              = bnd_outer
  )
}

# --- Ajustar y evaluar un modelo (espacial o espacio-temporal) ---
ajustar_evaluar <- function(tipo_modelo, cfg_malla, sel, datos) {

  dt       <- datos$dt
  n_per    <- datos$n_periodos
  covs     <- sel$covariables_candidatas[sel$vars_finales]
  y_tt     <- sel$y_train_test
  idx_test <- sel$idx_test_filas
  y_real   <- sel$y_test_real

  malla <- inla.mesh.2d(
    loc      = sel$coords_matriz,
    boundary = list(sel$bnd_inner, sel$bnd_outer),
    max.edge = cfg_malla$max.edge,
    cutoff   = cfg_malla$cutoff
  )
  spde <- inla.spde2.matern(mesh = malla, alpha = 2)

  ef <- paste("y ~ 0 + intercept +", paste(sel$vars_finales, collapse = " + "))

  if (tipo_modelo == "espacial") {

    indice <- inla.spde.make.index("campo_espacial", n.spde = spde$n.spde)
    A <- inla.spde.make.A(mesh = malla, loc = sel$coords_puntos)
    formula_mod <- as.formula(paste(ef, "+ f(campo_espacial, model = spde)"))
    tag <- "sim_s"

  } else {
    indice <- inla.spde.make.index("campo_espacial",
                                    n.spde = spde$n.spde,
                                    n.group = n_per)
    A <- inla.spde.make.A(mesh = malla, loc = sel$coords_puntos,
                           group = dt$ID_TIEMPO_SIM, n.group = n_per)
    formula_mod <- as.formula(paste(
      ef,
      "+ f(campo_espacial, model = spde,",
      "    group = campo_espacial.group,",
      "    control.group = list(model = 'ar1'))"
    ))
    tag <- "sim_st"
  }

  stk <- inla.stack(
    tag     = tag,
    data    = list(y = y_tt),
    A       = list(A, 1),
    effects = list(c(indice, list(intercept = 1)), covs),
    compress = FALSE
  )

  t0 <- Sys.time()
  modelo <- tryCatch(
    inla(
      formula           = formula_mod,
      data              = inla.stack.data(stk, spde = spde),
      family            = "gaussian",
      control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE),
      control.inla      = list(strategy = "laplace"),
      verbose           = FALSE
    ),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )
  t_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (is.null(modelo)) {
    return(data.table(
      Malla = cfg_malla$label, Modelo = tipo_modelo, n_nodos = malla$n,
      DIC = NA_real_, WAIC = NA_real_, RMSE = NA_real_, MAE = NA_real_,
      Cov95 = NA_real_, Rango_km = NA_real_, Sigma2 = NA_real_,
      Tiempo_min = round(t_min, 2)
    ))
  }

  idx_data <- inla.stack.index(stk, tag = tag)$data
  pred_mean <- modelo$summary.fitted.values$mean[idx_data][idx_test]
  pred_sd   <- modelo$summary.fitted.values$sd[idx_data][idx_test]

  rmse  <- sqrt(mean((pred_mean - y_real)^2, na.rm = TRUE))
  mae   <- mean(abs(pred_mean - y_real), na.rm = TRUE)
  lo95  <- pred_mean - 1.96 * pred_sd
  hi95  <- pred_mean + 1.96 * pred_sd
  cov95 <- mean(y_real >= lo95 & y_real <= hi95, na.rm = TRUE)

  spde_res   <- inla.spde.result(modelo, "campo_espacial", spde)
  rango_post <- exp(spde_res$summary.log.range.nominal$mean)
  sigma_post <- exp(spde_res$summary.log.variance.nominal$mean)

  etiq_modelo <- ifelse(tipo_modelo == "espacial",
                        "Solo espacial", "Espacio-temporal (AR1)")

  cat(sprintf("    [%s | %s] RMSE=%.4f MAE=%.4f Cov95=%.1f%% %.1f min\n",
              cfg_malla$label, etiq_modelo, rmse, mae, cov95 * 100, t_min))

  data.table(
    Malla      = cfg_malla$label,
    Modelo     = etiq_modelo,
    n_nodos    = malla$n,
    DIC        = round(modelo$dic$dic, 2),
    WAIC       = round(modelo$waic$waic, 2),
    RMSE       = round(rmse, 4),
    MAE        = round(mae, 4),
    Cov95      = round(cov95 * 100, 1),
    Rango_km   = round(rango_post, 2),
    Sigma2     = round(sigma_post, 4),
    Tiempo_min = round(t_min, 2)
  )
}

# ==============================================================================
# BUCLE PRINCIPAL: POR FRECUENCIA
# ==============================================================================

ejecuciones <- list(
  list(freq = "horario", periodo = "noviembre"),
  list(freq = "horario", periodo = "marzo"),
  list(freq = "diario",  periodo = "noviembre"),
  list(freq = "diario",  periodo = "marzo"),
  list(freq = "mensual", periodo = "anual")
)

# Modelos por frecuencia
modelos_por_freq <- list(
  horario = "espacial",
  diario  = c("espacial", "espacio-temporal"),
  mensual = "espacial"
)

resumen_global <- list()

for (ejec in ejecuciones) {

  freq    <- ejec$freq
  periodo <- ejec$periodo
  etiqueta <- sprintf("%s — %s", toupper(freq), toupper(periodo))

  cat("\n", strrep("#", 80), "\n")
  cat(sprintf("  %s\n", etiqueta))
  cat(strrep("#", 80), "\n")

  carpeta <- here("outputs", "simulacion", freq, periodo)
  dir.create(carpeta, showWarnings = FALSE, recursive = TRUE)

  # 1. Preparar datos
  cat("  Preparando datos...\n")
  datos <- preparar_datos(freq, periodo)
  cat(sprintf("  %d filas | %d estaciones | %d periodos (train=%d, test=%d)\n",
              nrow(datos$dt), uniqueN(datos$dt$ESTACION),
              datos$n_periodos, datos$n_train, datos$n_test))

  # 2. Selección de covariables (solo VIF)
  sel <- seleccionar_covariables(datos, carpeta)

  # 3. Ajustar modelos: 3 mallas × modelos según frecuencia
  tipos_modelo <- modelos_por_freq[[freq]]
  resultados <- list()

  for (malla_nombre in names(config_mallas)) {
    cfg <- config_mallas[[malla_nombre]]
    cat(sprintf("\n  === Malla: %s ===\n", cfg$label))

    for (tipo in tipos_modelo) {
      res <- ajustar_evaluar(tipo, cfg, sel, datos)
      resultados[[length(resultados) + 1]] <- res
    }

    gc()
  }

  tabla_comp <- rbindlist(resultados)
  setorder(tabla_comp, RMSE)

  # 4. Guardar tabla comparativa como PNG
  guardar_tabla_png(
    tabla_comp,
    titulo    = sprintf("Comparación de Modelos INLA-SPDE — %s", etiqueta),
    subtitulo = sprintf("Train: %s | Test: %s",
                        datos$lab_train, datos$lab_test),
    ruta_png  = file.path(carpeta, "tabla_comparacion_mallas_modelos.png"),
    ancho_px  = 1100
  )

  # 5. Gráfico RMSE
  tabla_comp[, Malla := factor(Malla,
    levels = c("Gruesa (8 km)", "Media (4 km)", "Fina (1 km)"))]

  ggplot(tabla_comp, aes(x = Malla, y = RMSE, fill = Modelo)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = round(RMSE, 4)),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Solo espacial" = "#B2182B",
                                  "Espacio-temporal (AR1)" = "#2166AC")) +
    labs(
      title = sprintf("RMSE por Malla y Modelo — %s", etiqueta),
      subtitle = sprintf("Test: %s", datos$lab_test),
      x = "Resolución de malla", y = "RMSE (log NO₂)", fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")

  ggsave(file.path(carpeta, "grafico_rmse.png"),
         width = 10, height = 6, dpi = 300)

  # 6. Gráfico tiempo de cómputo
  ggplot(tabla_comp, aes(x = Malla, y = Tiempo_min, fill = Modelo)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = paste0(Tiempo_min, " min")),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Solo espacial" = "#B2182B",
                                  "Espacio-temporal (AR1)" = "#2166AC")) +
    labs(
      title = sprintf("Tiempo de Cómputo — %s", etiqueta),
      x = "Resolución de malla", y = "Tiempo (minutos)", fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")

  ggsave(file.path(carpeta, "grafico_tiempo.png"),
         width = 10, height = 6, dpi = 300)

  # Acumular para resumen global
  tabla_comp[, Frecuencia := toupper(freq)]
  tabla_comp[, Periodo := toupper(periodo)]
  resumen_global[[length(resumen_global) + 1]] <- tabla_comp

  cat(sprintf("\n  %s completado. Resultados en: %s\n", etiqueta, carpeta))
}

# ==============================================================================
# TABLA RESUMEN GLOBAL (todas las frecuencias)
# ==============================================================================

tabla_global <- rbindlist(resumen_global)
setorder(tabla_global, Frecuencia, Periodo, RMSE)

mejor_por_config <- tabla_global[, .SD[which.min(RMSE)],
                                  by = .(Frecuencia, Periodo)]

guardar_tabla_png(
  mejor_por_config[, .(Frecuencia, Periodo, Malla, Modelo, RMSE, MAE, Cov95,
                       Rango_km, Sigma2, Tiempo_min)],
  titulo    = "Mejor Configuración por Frecuencia y Periodo",
  subtitulo = "Modelo con menor RMSE en cada combinación",
  ruta_png  = here("outputs", "simulacion", "tabla_resumen_global.png"),
  ancho_px  = 1100
)

guardar_tabla_png(
  tabla_global[, .(Frecuencia, Periodo, Malla, Modelo, n_nodos, DIC, WAIC,
                   RMSE, MAE, Cov95, Tiempo_min)],
  titulo    = "Resultados Completos — Todas las Frecuencias y Periodos",
  subtitulo = "Horario/Mensual: solo espacial | Diario: espacial + espacio-temporal",
  ruta_png  = here("outputs", "simulacion", "tabla_resultados_completos.png"),
  ancho_px  = 1300
)

cat("\n", strrep("=", 80), "\n")
cat("  SIMULACIÓN COMPLETADA\n")
cat(strrep("=", 80), "\n")
cat("  Resultados guardados en:", here("outputs", "simulacion"), "\n")

