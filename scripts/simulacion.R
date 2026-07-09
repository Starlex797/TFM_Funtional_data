# ==============================================================================
# SIMULACIÓN INLA-SPDE: BLOQUES 0–7 POR FRECUENCIA (MENSUAL → DIARIO → HORARIO)
# ==============================================================================
# Objetivo:
#   - Comparar 3 frecuencias temporales: MENSUAL, DIARIA, HORARIA
#   - Para cada frecuencia (Bloques 0–7):
#       • B0: Preparación de datos
#       • B1: Selección de covariables (VIF + significancia + stepwise DIC)
#       • B2: Modelo espacial SPDE (4 mallas)
#       • B3: Diagnóstico de residuos del modelo espacial
#       • B4: Modelo espacio-temporal AR1 (si B3 lo justifica)
#       • B5: Diagnóstico de residuos del AR1
#       • B6: Comparación y selección de modelos
#       • B7: Tiempo de recuperación de zonas contaminadas
#   - Las variables candidatas se definen de forma independiente por frecuencia
#     en la sección CONFIGURACIÓN GENERAL (COVS_MENSUAL/DIARIO/HORARIO_NOMBRES).
#   - Outputs organizados en subcarpetas por bloque:
#       seleccion/ | ACF/ | variograma/ | QQ/ | comparacion/ | recuperacion/
# ==============================================================================

library(INLA)
library(fmesher)
library(data.table)
library(sf)
library(here)
library(ggplot2)
library(car)
library(gt)
library(gstat)

set.seed(4827)

# ==============================================================================
# CONFIGURACIÓN GENERAL
# ==============================================================================

VIF_UMBRAL     <- 5
DIC_MEJORA_MIN <- 2

config_mallas <- list(
  muy_gruesa = list(max.edge = c(12, 18), cutoff = 1.0,  label = "Muy gruesa (12 km)"),
  gruesa     = list(max.edge = c(8, 12),  cutoff = 0.5,  label = "Gruesa (8 km)"),
  media      = list(max.edge = c(4,  8),  cutoff = 0.5,  label = "Media (4 km)"),
  fina       = list(max.edge = c(1,  4),  cutoff = 0.25, label = "Fina (1 km)")
)

# ==============================================================================
# COVARIABLES CANDIDATAS POR FRECUENCIA
# Modificar aquí para cambiar las variables de cada escala temporal.
# Nombres : columnas estandarizadas presentes en el dataset.
# Alias   : nombres internos que usará INLA (sin espacios ni caracteres especiales).
# Mantener el mismo orden en Nombres y Alias.
# ==============================================================================

COVS_MENSUAL_NOMBRES  <- c("intensidad", "Temperatura",
                             "Precipitaciones", "Velocidad Viento")
COVS_MENSUAL_ALIAS    <- c("trafico_intensidad", "temperatura",
                             "precipitacion", "velocidad_viento")

COVS_DIARIO_NOMBRES   <- c("intensidad", "Temperatura",
                             "Precipitaciones", "Velocidad Viento")
COVS_DIARIO_ALIAS     <- c("trafico_intensidad", "temperatura",
                             "precipitacion", "velocidad_viento")

COVS_HORARIO_NOMBRES  <- c("intensidad", "Temperatura",
                             "Precipitaciones", "Velocidad Viento")
COVS_HORARIO_ALIAS    <- c("trafico_intensidad", "temperatura",
                             "precipitacion", "velocidad_viento")

# Columnas raw para re-estandarizar en mensual (no modificar)
COLS_RAW_MENSUAL <- c("intensidad_raw", "Temperatura_raw",
                       "Precipitaciones_raw", "Velocidad Viento_raw")

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
    dt <- readRDS(here("data", "processed", "Maestro", "horario",
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
    dt <- readRDS(here("data", "processed", "Maestro", "diario",
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
    dt_diario <- readRDS(here("data", "processed", "Maestro", "diario",
                              "dataset_maestro_inla_2025_DIARIO.rds"))
    setDT(dt_diario)
    dt_diario[, ANIO_MES := as.Date(format(FECHA, "%Y-%m-01"))]

    dt <- dt_diario[, .(
      DATO_NO2                  = mean(DATO_DIARIO,             na.rm = TRUE),
      intensidad_raw            = mean(intensidad_raw,          na.rm = TRUE),
      Temperatura_raw           = mean(Temperatura_raw,         na.rm = TRUE),
      Precipitaciones_raw       = sum(Precipitaciones_raw,      na.rm = TRUE),
      `Velocidad Viento_raw`    = mean(`Velocidad Viento_raw`,  na.rm = TRUE),
      n_dias                    = .N
    ), by = .(ESTACION, ANIO_MES, barrio, distrito, LONGITUD, LATITUD,
              ID_DISTRITO)]

    setnames(dt, "ANIO_MES", "FECHA")
    dt[, LOG_NO2 := log(DATO_NO2)]

    # Re-estandarizar a escala mensual
    for (i in seq_along(COLS_RAW_MENSUAL)) {
      dt[, (COVS_MENSUAL_NOMBRES[i]) := scale(get(COLS_RAW_MENSUAL[i]))[, 1]]
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

  # NO2 en escala original (exp del log) para análisis exploratorio
  dt[, NO2_raw := exp(LOG_NO2)]

  list(
    dt           = dt,
    n_periodos   = uniqueN(dt$ID_TIEMPO_SIM),
    n_train      = n_train_p,
    n_test       = n_test_p,
    lab_train    = lab_train,
    lab_test     = lab_test
  )
}

# --- Heatmap de correlación (datos no estandarizados) + Análisis temporal ---

# --- Crear subcarpetas de outputs por bloque ---
crear_subcarpetas <- function(carpeta_base) {
  subs <- list(
    seleccion    = file.path(carpeta_base, "seleccion"),
    ACF          = file.path(carpeta_base, "ACF"),
    variograma   = file.path(carpeta_base, "variograma"),
    QQ           = file.path(carpeta_base, "QQ"),
    comparacion  = file.path(carpeta_base, "comparacion"),
    recuperacion = file.path(carpeta_base, "recuperacion"),
    modelos      = carpeta_base
  )
  lapply(subs, dir.create, showWarnings = FALSE, recursive = TRUE)
  subs
}

# ==============================================================================
# BLOQUE 1 — SELECCIÓN DE COVARIABLES (independiente por frecuencia)
# ==============================================================================
# covs_nombres : nombres de columnas estandarizadas en el dataset
# covs_alias   : aliases sin espacios para INLA (mismo orden)
# carpeta_sel  : subcarpeta donde guardar tablas de seleccion
# ==============================================================================

seleccionar_covariables <- function(datos, carpeta_sel,
                                    covs_nombres, covs_alias) {
  dt <- datos$dt

  covs_disponibles <- intersect(covs_nombres, names(dt))
  alias_map <- setNames(covs_alias[match(covs_disponibles, covs_nombres)],
                        covs_disponibles)
  todas_covs <- setNames(
    lapply(covs_disponibles, function(nm) dt[[nm]]),
    alias_map[covs_disponibles]
  )

  vars_pool <- names(todas_covs)

  # Pre-filtro varianza cero
  for (v in vars_pool) {
    sd_v <- sd(todas_covs[[v]][dt$es_train], na.rm = TRUE)
    if (is.na(sd_v) || sd_v < 1e-10) {
      cat(sprintf("    Descartada '%s' (sd ~ 0 en train)\n", v))
      vars_pool <- setdiff(vars_pool, v)
    }
  }
  todas_covs    <- todas_covs[vars_pool]
  covs_candidatas <- todas_covs

  cat(sprintf("  Vars candidatas: %s\n", paste(vars_pool, collapse = ", ")))

  y_train_test    <- ifelse(dt$es_train, dt$LOG_NO2, NA_real_)
  idx_test_filas  <- which(!dt$es_train)
  idx_train_filas <- which( dt$es_train)
  y_test_real     <- dt$LOG_NO2[idx_test_filas]

  # --- VIF backward ---
  vars_vif <- vars_pool
  if (length(vars_pool) >= 2) {
    df_vif <- as.data.frame(lapply(covs_candidatas, function(x) x[dt$es_train]))
    df_vif$LOG_NO2 <- dt$LOG_NO2[dt$es_train]
    df_vif <- na.omit(df_vif)
    vif_vals <- NULL
    cat("\n  --- VIF hacia atras (umbral =", VIF_UMBRAL, ") ---\n")
    repeat {
      if (length(vars_vif) < 2) break
      fml      <- as.formula(paste("LOG_NO2 ~", paste(vars_vif, collapse = " + ")))
      lm_fit   <- lm(fml, data = df_vif)
      vif_vals <- car::vif(lm_fit)
      cat(sprintf("    Vars: %s\n", paste(vars_vif, collapse = ", ")))
      if (max(vif_vals) <= VIF_UMBRAL) { cat("    -> Fin VIF.\n"); break }
      var_elim <- names(which.max(vif_vals))
      cat(sprintf("    -> Eliminando '%s' (VIF=%.3f)\n", var_elim, max(vif_vals)))
      vars_vif <- setdiff(vars_vif, var_elim)
    }
    if (!is.null(vif_vals)) {
      tabla_vif <- data.frame(Variable = names(vif_vals),
                               VIF = round(as.numeric(vif_vals), 4), row.names = NULL)
      guardar_tabla_png(
        tabla_vif[order(tabla_vif$VIF, decreasing = TRUE), ],
        titulo    = "Seleccion de Covariables - VIF",
        subtitulo = sprintf("Umbral VIF = %g | %d variables retenidas",
                            VIF_UMBRAL, length(vars_vif)),
        ruta_png  = file.path(carpeta_sel, "tabla_seleccion_vif.png")
      )
    }
  }
  cat(sprintf("  Variables tras VIF: %s\n", paste(vars_vif, collapse = ", ")))

  # --- SPDE de referencia (malla media) ---
  coords_puntos <- as.matrix(dt[, .(X_km, Y_km)])
  coords_matriz <- as.matrix(unique(dt[, .(X_km, Y_km)]))
  bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.05, resolution = 50)
  bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.2)
  malla_ref  <- inla.mesh.2d(loc = coords_matriz, boundary = list(bnd_inner, bnd_outer),
                               max.edge = c(4, 8), cutoff = 0.5)
  spde_ref   <- inla.spde2.matern(mesh = malla_ref, alpha = 2)
  indice_ref <- inla.spde.make.index("campo_espacial", n.spde = spde_ref$n.spde)
  A_ref      <- inla.spde.make.A(mesh = malla_ref, loc = coords_puntos)

  formula_full <- as.formula(paste(
    "y ~ 0 + intercept +", paste(vars_vif, collapse = " + "),
    "+ f(campo_espacial, model = spde_ref)"
  ))
  stk_full <- inla.stack(
    tag = "full", data = list(y = y_train_test), A = list(A_ref, 1),
    effects = list(c(indice_ref, list(intercept = 1)), todas_covs[vars_vif]),
    compress = FALSE
  )
  cat("  Ajustando modelo completo (significancia bayesiana)...\n")
  modelo_full <- inla(
    formula = formula_full, data = inla.stack.data(stk_full, spde = spde_ref),
    family = "gaussian",
    control.predictor = list(A = inla.stack.A(stk_full), compute = TRUE),
    control.compute   = list(dic = TRUE, waic = FALSE, cpo = FALSE),
    control.inla      = list(strategy = "laplace"), verbose = FALSE
  )

  sf_full  <- modelo_full$summary.fixed
  sf_vars  <- sf_full[rownames(sf_full) != "intercept", , drop = FALSE]
  sig_mask <- sf_vars[, "0.025quant"] > 0 | sf_vars[, "0.975quant"] < 0
  vars_sig <- rownames(sf_vars)[sig_mask]

  tabla_sig <- data.frame(
    Variable = rownames(sf_vars),
    Media    = round(sf_vars$mean,           4),
    SD       = round(sf_vars$sd,             4),
    Q2.5     = round(sf_vars[, "0.025quant"],4),
    Q97.5    = round(sf_vars[, "0.975quant"],4),
    Sig_95   = ifelse(sig_mask, "Si", "No"),
    row.names = NULL, check.names = FALSE
  )
  guardar_tabla_png(
    tabla_sig,
    titulo    = "Significancia Bayesiana (IC 95%)",
    subtitulo = sprintf("%d significativas de %d", sum(sig_mask), nrow(sf_vars)),
    ruta_png  = file.path(carpeta_sel, "tabla_significancia_bayesiana.png")
  )
  cat(sprintf("  Significativas: %s\n", paste(vars_sig, collapse = ", ")))

  # --- Stepwise DIC ---
  ajustar_sw <- function(vars_mod) {
    fml_sw <- as.formula(paste("y ~ 0 + intercept +", paste(vars_mod, collapse = " + "),
                                "+ f(campo_espacial, model = spde_ref)"))
    stk_sw <- inla.stack(
      tag = "sw", data = list(y = y_train_test), A = list(A_ref, 1),
      effects = list(c(indice_ref, list(intercept = 1)), todas_covs[vars_mod]),
      compress = FALSE
    )
    mod <- tryCatch(
      inla(formula = fml_sw, data = inla.stack.data(stk_sw, spde = spde_ref),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(stk_sw), compute = TRUE),
           control.compute   = list(dic = TRUE, waic = FALSE, cpo = FALSE),
           control.inla      = list(strategy = "laplace"), verbose = FALSE),
      error = function(e) { message("  ERROR sw: ", e$message); NULL }
    )
    if (is.null(mod)) return(list(DIC = Inf, RMSE = NA_real_, MAE = NA_real_))
    idx_d <- inla.stack.index(stk_sw, tag = "sw")$data
    pred  <- mod$summary.fitted.values$mean[idx_d][idx_test_filas]
    list(DIC  = mod$dic$dic,
         RMSE = sqrt(mean((pred - y_test_real)^2, na.rm = TRUE)),
         MAE  = mean(abs(pred - y_test_real), na.rm = TRUE))
  }

  vars_sw_actual <- if (length(vars_sig) > 0) vars_sig else vars_vif
  res_actual <- ajustar_sw(vars_sw_actual)
  tabla_sw   <- list()
  tabla_sw[[1]] <- data.frame(
    Iter = 0L, Accion = "inicial", Variable = "-",
    Variables = paste(vars_sw_actual, collapse = " + "),
    DIC = round(res_actual$DIC, 2), RMSE = round(res_actual$RMSE, 4),
    MAE = round(res_actual$MAE, 4), stringsAsFactors = FALSE
  )
  cat(sprintf("  Stepwise inicial: %s (DIC=%.2f)\n",
              paste(vars_sw_actual, collapse = " + "), res_actual$DIC))

  for (iter in seq_len(length(vars_vif) + 1L)) {
    cands <- data.frame(accion = character(), variable = character(),
                        DIC = numeric(), RMSE = numeric(), MAE = numeric(),
                        stringsAsFactors = FALSE)
    for (v in setdiff(vars_vif, vars_sw_actual)) {
      r <- ajustar_sw(c(vars_sw_actual, v))
      cands <- rbind(cands, data.frame(accion = "anadir", variable = v,
                                       DIC = r$DIC, RMSE = r$RMSE, MAE = r$MAE,
                                       stringsAsFactors = FALSE))
    }
    if (length(vars_sw_actual) >= 2) {
      for (v in vars_sw_actual) {
        r <- ajustar_sw(setdiff(vars_sw_actual, v))
        cands <- rbind(cands, data.frame(accion = "eliminar", variable = v,
                                         DIC = r$DIC, RMSE = r$RMSE, MAE = r$MAE,
                                         stringsAsFactors = FALSE))
      }
    }
    if (nrow(cands) == 0) break
    mejor  <- cands[which.min(cands$DIC), ]
    mejora <- res_actual$DIC - mejor$DIC
    if (mejora < DIC_MEJORA_MIN) {
      cat(sprintf("  Stepwise convergido (DDIC=%.2f)\n", mejora)); break
    }
    vars_sw_actual <- if (mejor$accion == "anadir") c(vars_sw_actual, mejor$variable)
                      else setdiff(vars_sw_actual, mejor$variable)
    res_actual <- list(DIC = mejor$DIC, RMSE = mejor$RMSE, MAE = mejor$MAE)
    cat(sprintf("  Iter %d: %s '%s' (DDIC=%.2f)\n",
                iter, mejor$accion, mejor$variable, mejora))
    tabla_sw[[iter + 1]] <- data.frame(
      Iter = iter, Accion = mejor$accion, Variable = mejor$variable,
      Variables = paste(vars_sw_actual, collapse = " + "),
      DIC = round(mejor$DIC, 2), RMSE = round(mejor$RMSE, 4),
      MAE = round(mejor$MAE, 4), stringsAsFactors = FALSE
    )
  }

  guardar_tabla_png(
    do.call(rbind, tabla_sw),
    titulo    = "Stepwise DIC - Seleccion de Covariables",
    subtitulo = sprintf("DDIC min = %g | Vars finales: %s",
                        DIC_MEJORA_MIN, paste(vars_sw_actual, collapse = ", ")),
    ruta_png  = file.path(carpeta_sel, "tabla_stepwise_dic.png"), ancho_px = 1000
  )
  cat(sprintf("  Variables finales: %s\n", paste(vars_sw_actual, collapse = ", ")))

  list(
    vars_finales    = vars_sw_actual,
    todas_covs      = todas_covs,
    y_train_test    = y_train_test,
    idx_test_filas  = idx_test_filas,
    idx_train_filas = idx_train_filas,
    y_test_real     = y_test_real,
    y_train_real    = dt$LOG_NO2[dt$es_train],
    coords_puntos   = coords_puntos,
    coords_matriz   = coords_matriz,
    bnd_inner       = bnd_inner,
    bnd_outer       = bnd_outer
  )
}

# ==============================================================================
# BLOQUES 2 / 4 — AJUSTAR Y EVALUAR MODELO
# ==============================================================================

ajustar_evaluar <- function(tipo_modelo, cfg_malla, sel, datos) {

  dt    <- datos$dt
  n_per <- datos$n_periodos

  covs_modelo <- sel$todas_covs[sel$vars_finales[sel$vars_finales %in% names(sel$todas_covs)]]

  malla <- inla.mesh.2d(
    loc = sel$coords_matriz, boundary = list(sel$bnd_inner, sel$bnd_outer),
    max.edge = cfg_malla$max.edge, cutoff = cfg_malla$cutoff
  )
  spde <- inla.spde2.matern(mesh = malla, alpha = 2)
  ef   <- paste("y ~ 0 + intercept +", paste(sel$vars_finales, collapse = " + "))

  if (tipo_modelo == "espacial") {
    indice      <- inla.spde.make.index("campo_espacial", n.spde = spde$n.spde)
    A           <- inla.spde.make.A(mesh = malla, loc = sel$coords_puntos)
    formula_mod <- as.formula(paste(ef, "+ f(campo_espacial, model = spde)"))
    tag         <- "sim_s"
  } else {
    indice <- inla.spde.make.index("campo_espacial",
                                    n.spde = spde$n.spde, n.group = n_per)
    A <- inla.spde.make.A(mesh = malla, loc = sel$coords_puntos,
                           group = dt$ID_TIEMPO_SIM, n.group = n_per)
    formula_mod <- as.formula(paste(
      ef, "+ f(campo_espacial, model = spde,",
      "    group = campo_espacial.group,",
      "    control.group = list(model = 'ar1'))"
    ))
    tag <- "sim_st"
  }

  stk <- inla.stack(
    tag = tag, data = list(y = sel$y_train_test), A = list(A, 1),
    effects = list(c(indice, list(intercept = 1)), covs_modelo), compress = FALSE
  )

  t0 <- Sys.time()
  modelo <- tryCatch(
    inla(
      formula = formula_mod, data = inla.stack.data(stk, spde = spde),
      family = "gaussian",
      control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
      control.inla      = list(strategy = "laplace"), verbose = FALSE
    ),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )
  t_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (is.null(modelo)) {
    return(list(
      metricas = data.table(Malla = cfg_malla$label, Modelo = tipo_modelo,
                             n_nodos = malla$n, DIC = NA_real_, WAIC = NA_real_,
                             RMSE = NA_real_, MAE = NA_real_, Cov95 = NA_real_,
                             Rango_km = NA_real_, Sigma2 = NA_real_,
                             Tiempo_min = round(t_min, 2)),
      modelo = NULL, residuos = NULL, stk = stk, tag = tag, malla = malla, spde = spde
    ))
  }

  idx_data  <- inla.stack.index(stk, tag = tag)$data
  pred_all  <- modelo$summary.fitted.values$mean[idx_data]
  pred_test <- pred_all[sel$idx_test_filas]
  pred_sd   <- modelo$summary.fitted.values$sd[idx_data][sel$idx_test_filas]

  residuos_train <- sel$y_train_real - pred_all[sel$idx_train_filas]

  rmse  <- sqrt(mean((pred_test - sel$y_test_real)^2, na.rm = TRUE))
  mae   <- mean(abs(pred_test  - sel$y_test_real),    na.rm = TRUE)
  cov95 <- mean(sel$y_test_real >= pred_test - 1.96 * pred_sd &
                sel$y_test_real <= pred_test + 1.96 * pred_sd, na.rm = TRUE)

  spde_res   <- inla.spde.result(modelo, "campo_espacial", spde)
  rango_post <- exp(spde_res$summary.log.range.nominal$mean)
  sigma_post <- exp(spde_res$summary.log.variance.nominal$mean)
  etiq       <- ifelse(tipo_modelo == "espacial", "Solo espacial", "Espacio-temporal (AR1)")

  cat(sprintf("    [%s | %s] RMSE=%.4f MAE=%.4f Cov95=%.1f%% %.1fmin\n",
              cfg_malla$label, etiq, rmse, mae, cov95 * 100, t_min))

  list(
    metricas = data.table(
      Malla = cfg_malla$label, Modelo = etiq, n_nodos = malla$n,
      DIC = round(modelo$dic$dic, 2), WAIC = round(modelo$waic$waic, 2),
      RMSE = round(rmse, 4), MAE = round(mae, 4), Cov95 = round(cov95 * 100, 1),
      Rango_km = round(rango_post, 2), Sigma2 = round(sigma_post, 4),
      Tiempo_min = round(t_min, 2)
    ),
    modelo = modelo, residuos = residuos_train,
    stk = stk, tag = tag, malla = malla, spde = spde
  )
}

# ==============================================================================
# BLOQUES 3 / 5 — DIAGNOSTICO DE RESIDUOS
# ==============================================================================
# B3: verifica si modelo espacial deja autocorrelacion temporal/espacial.
#     Si p < 0.05 en Ljung-Box -> AR1 justificado.
# B5: verifica que AR1 elimino la autocorrelacion.
# ==============================================================================

diagnostico_residuos <- function(residuos, datos, sel,
                                  dir_acf, dir_var, dir_qq,
                                  etiqueta = "", lag_max = 10) {

  dt_train <- datos$dt[es_train == TRUE]

  # --- ACF + Ljung-Box (carpeta ACF/) ---
  res_temporal <- tapply(residuos, dt_train$ID_TIEMPO_SIM, mean, na.rm = TRUE)
  lag_uso      <- min(lag_max, length(res_temporal) - 1L)

  png(file.path(dir_acf, paste0("acf_", etiqueta, ".png")),
      width = 800, height = 500)
  acf(res_temporal, lag.max = lag_uso,
      main = sprintf("ACF residuos - %s", etiqueta),
      col = "#2166AC", lwd = 2)
  dev.off()

  lj      <- Box.test(res_temporal, lag = lag_uso, type = "Ljung-Box")
  acf_sig <- lj$p.value < 0.05
  cat(sprintf("  [%s] Ljung-Box p=%.4f -> %s\n", etiqueta, lj$p.value,
              ifelse(acf_sig, "autocorrelacion DETECTADA (AR1 justificado)",
                              "sin autocorrelacion significativa")))

  # --- Variograma de residuos (carpeta variograma/) ---
  res_por_est <- tapply(residuos, dt_train$ESTACION, mean, na.rm = TRUE)
  coords_est  <- unique(dt_train[, .(ESTACION, X_km, Y_km)])
  setkey(coords_est, ESTACION)
  coords_est  <- coords_est[names(res_por_est)]

  res_sf <- st_as_sf(
    data.frame(x = coords_est$X_km, y = coords_est$Y_km,
               res = as.numeric(res_por_est)),
    coords = c("x", "y")
  )
  v_emp <- tryCatch(variogram(res ~ 1, data = res_sf, cutoff = 20, width = 4),
                    error = function(e) NULL)
  if (!is.null(v_emp)) {
    png(file.path(dir_var, paste0("variograma_", etiqueta, ".png")),
        width = 700, height = 500)
    print(plot(v_emp, main = sprintf("Variograma residuos - %s", etiqueta),
               xlab = "Distancia (km)", ylab = "Semivarianza"))
    dev.off()
  }

  # --- Q-Q plot (carpeta QQ/) ---
  png(file.path(dir_qq, paste0("qq_", etiqueta, ".png")),
      width = 600, height = 600)
  qqnorm(residuos, main = sprintf("Q-Q residuos - %s", etiqueta),
         pch = 16, cex = 0.6, col = "#666666")
  qqline(residuos, col = "#B2182B", lwd = 2)
  dev.off()

  invisible(list(acf_significativa = acf_sig, ljung_box_p = lj$p.value))
}

# ==============================================================================
# BLOQUE 6 — COMPARACION DE SIGNIFICANCIA (espacial vs AR1)
# ==============================================================================

comparar_significancia <- function(modelo_esp, modelo_ar1, vars_finales,
                                    carpeta_comp, etiqueta) {
  sf_esp <- modelo_esp$summary.fixed
  sf_ar1 <- if (!is.null(modelo_ar1)) modelo_ar1$summary.fixed else NULL

  vars_mostrar <- vars_finales[vars_finales != "intercept"]

  tabla <- lapply(vars_mostrar, function(v) {
    row_e <- if (v %in% rownames(sf_esp)) sf_esp[v, ] else NULL
    row_a <- if (!is.null(sf_ar1) && v %in% rownames(sf_ar1)) sf_ar1[v, ] else NULL
    sig_e <- if (!is.null(row_e)) row_e[, "0.025quant"] > 0 | row_e[, "0.975quant"] < 0 else NA
    sig_a <- if (!is.null(row_a)) row_a[, "0.025quant"] > 0 | row_a[, "0.975quant"] < 0 else NA

    interp <- if (is.na(sig_e) || is.na(sig_a)) "Solo en un modelo"
              else if ( sig_e &&  sig_a) "Efecto genuino"
              else if ( sig_e && !sig_a) "Proxy autocorrelacion"
              else if (!sig_e &&  sig_a) "Emerge con AR1"
              else                        "No significativa"

    data.frame(
      Variable       = v,
      Coef_Espacial  = if (!is.null(row_e)) round(row_e$mean, 4) else NA,
      Sig_Espacial   = ifelse(is.na(sig_e), "---", ifelse(sig_e, "Si", "No")),
      Coef_AR1       = if (!is.null(row_a)) round(row_a$mean, 4) else NA,
      Sig_AR1        = if (is.null(sf_ar1)) "---"
                       else ifelse(is.na(sig_a), "---", ifelse(sig_a, "Si", "No")),
      Interpretacion = interp,
      row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE
    )
  })

  guardar_tabla_png(
    do.call(rbind, tabla),
    titulo    = sprintf("Significancia Comparada - %s", etiqueta),
    subtitulo = "Proxy autocorrelacion: sig. en espacial pero no en AR1 | Efecto genuino: sig. en ambos",
    ruta_png  = file.path(carpeta_comp, "tabla_significancia_comparada.png"),
    ancho_px  = 1000
  )
  cat(sprintf("  Tabla de significancia comparada guardada en comparacion/ - %s\n", etiqueta))
}

# ==============================================================================
# BLOQUE 7 — TIEMPO DE RECUPERACION DE ZONAS CONTAMINADAS
# ==============================================================================
# Usa rho (parametro AR1) para calcular tau = -1/log(rho) en unidades temporales.
# ==============================================================================

calcular_recuperacion <- function(modelo_ar1, spde_obj, datos, sel,
                                   stk_gan, tag_gan, carpeta_rec, etiqueta) {

  if (is.null(modelo_ar1)) {
    cat("  [B7] Sin modelo AR1: se omite calculo de recuperacion.\n")
    return(invisible(NULL))
  }

  hp      <- modelo_ar1$summary.hyperpar
  rho_row <- grep("GroupRho|Rho for", rownames(hp), value = TRUE)
  if (length(rho_row) == 0) {
    cat("  [B7] Parametro AR1 no encontrado en hyperpar.\n")
    return(invisible(NULL))
  }

  rho    <- hp[rho_row[1], "mean"]
  rho_lo <- hp[rho_row[1], "0.025quant"]
  rho_hi <- hp[rho_row[1], "0.975quant"]
  tau    <- -1 / log(abs(rho))
  tau_lo <- -1 / log(abs(rho_hi))
  tau_hi <- -1 / log(abs(rho_lo))

  cat(sprintf("  [B7] rho = %.3f [%.3f, %.3f] | tau = %.2f periodos [%.2f, %.2f]\n",
              rho, rho_lo, rho_hi, tau, tau_lo, tau_hi))

  guardar_tabla_png(
    data.frame(
      Parametro = c("rho (AR1)", "tau (periodos hasta 1/e)"),
      Media     = round(c(rho,    tau),    3),
      IC_2.5    = round(c(rho_lo, tau_lo), 3),
      IC_97.5   = round(c(rho_hi, tau_hi), 3),
      row.names = NULL
    ),
    titulo    = sprintf("Tiempo de Recuperacion - %s", etiqueta),
    subtitulo = "tau = -1/log(rho): periodos para decaer al 37% del exceso inicial",
    ruta_png  = file.path(carpeta_rec, "tabla_recuperacion.png")
  )

  # Mapa de hotspots
  idx_data  <- inla.stack.index(stk_gan, tag = tag_gan)$data
  pred_mean <- modelo_ar1$summary.fitted.values$mean[idx_data]
  dt_train  <- datos$dt[es_train == TRUE]

  pred_dt <- data.table(ESTACION = dt_train$ESTACION,
                         pred_log = pred_mean)
  pred_est <- pred_dt[, .(NO2_pred = mean(exp(pred_log), na.rm = TRUE)), by = ESTACION]
  pred_est <- merge(pred_est, unique(datos$dt[, .(ESTACION, X_km, Y_km)]), by = "ESTACION")

  ggplot(pred_est, aes(x = X_km, y = Y_km, color = NO2_pred, size = NO2_pred)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(
      low = "#2166AC", mid = "#FFFFBF", high = "#B2182B",
      midpoint = median(pred_est$NO2_pred, na.rm = TRUE),
      name = "NO2 pred\n(ug/m3)"
    ) +
    scale_size_continuous(guide = "none", range = c(2, 8)) +
    labs(
      title    = sprintf("Hotspots de NO2 - %s", etiqueta),
      subtitle = sprintf("rho = %.3f | tau ~ %.1f periodos de recuperacion", rho, tau),
      x = "X (km UTM)", y = "Y (km UTM)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(carpeta_rec, "mapa_hotspots_recuperacion.png"),
         width = 9, height = 7, dpi = 300)

  invisible(list(rho = rho, tau = tau))
}

# ==============================================================================
# BUCLE PRINCIPAL — MENSUAL -> DIARIO -> HORARIO (Bloques 0-7)
# ==============================================================================

ejecuciones <- list(
  list(freq = "mensual", periodo = "anual"),
  list(freq = "diario",  periodo = "noviembre"),
  list(freq = "diario",  periodo = "marzo"),
  list(freq = "horario", periodo = "noviembre"),
  list(freq = "horario", periodo = "marzo")
)

resumen_global <- list()

for (ejec in ejecuciones) {

  freq     <- ejec$freq
  periodo  <- ejec$periodo
  etiqueta <- sprintf("%s - %s", toupper(freq), toupper(periodo))

  cat("\n", strrep("#", 80), "\n")
  cat(sprintf("  %s\n", etiqueta))
  cat(strrep("#", 80), "\n")

  carpeta <- here("outputs", "simulacion", freq, periodo)
  subs    <- crear_subcarpetas(carpeta)

  # Covariables según frecuencia (modificar COVS_*_NOMBRES en config)
  covs_cfg <- switch(freq,
    mensual = list(nombres = COVS_MENSUAL_NOMBRES, alias = COVS_MENSUAL_ALIAS),
    diario  = list(nombres = COVS_DIARIO_NOMBRES,  alias = COVS_DIARIO_ALIAS),
    horario = list(nombres = COVS_HORARIO_NOMBRES,  alias = COVS_HORARIO_ALIAS)
  )

  # ── B0: Preparar datos ────────────────────────────────────────────────────
  cat("  [B0] Preparando datos...\n")
  datos <- preparar_datos(freq, periodo)
  cat(sprintf("  %d filas | %d estaciones | %d periodos (train=%d, test=%d)\n",
              nrow(datos$dt), uniqueN(datos$dt$ESTACION),
              datos$n_periodos, datos$n_train, datos$n_test))

  # ── B1: Seleccion de covariables ─────────────────────────────────────────
  cat("  [B1] Seleccion de covariables...\n")
  sel <- seleccionar_covariables(
    datos        = datos,
    carpeta_sel  = subs$seleccion,
    covs_nombres = covs_cfg$nombres,
    covs_alias   = covs_cfg$alias
  )

  # ── B2: Modelo espacial (4 mallas) ───────────────────────────────────────
  cat("  [B2] Modelo espacial (4 mallas)...\n")
  res_espacial <- list()
  for (malla_nombre in names(config_mallas)) {
    cfg <- config_mallas[[malla_nombre]]
    cat(sprintf("  === Malla: %s ===\n", cfg$label))
    res_espacial[[malla_nombre]] <- ajustar_evaluar("espacial", cfg, sel, datos)
    gc()
  }

  rmse_esp        <- sapply(res_espacial, function(r) ifelse(is.na(r$metricas$RMSE), Inf, r$metricas$RMSE))
  malla_mejor_esp <- names(which.min(rmse_esp))
  res_esp_gan     <- res_espacial[[malla_mejor_esp]]

  # ── B3: Diagnostico residuos espacial ────────────────────────────────────
  cat("  [B3] Diagnostico de residuos del modelo espacial...\n")
  etiq_esp <- paste0("espacial_", gsub("[^A-Za-z0-9]", "_", malla_mejor_esp))
  diag_esp <- diagnostico_residuos(
    residuos = res_esp_gan$residuos, datos = datos, sel = sel,
    dir_acf  = subs$ACF,
    dir_var  = subs$variograma,
    dir_qq   = subs$QQ,
    etiqueta = etiq_esp
  )

  # ── B4 + B5: AR1 (solo si B3 lo justifica) ───────────────────────────────
  res_ar1         <- list()
  modelo_ar1_gan  <- NULL
  spde_ar1_gan    <- res_esp_gan$spde
  stk_gan         <- res_esp_gan$stk
  tag_gan         <- res_esp_gan$tag
  malla_mejor_ar1 <- malla_mejor_esp

  if (diag_esp$acf_significativa) {
    cat("  [B4] Autocorrelacion detectada -> ajustando AR1 (4 mallas)...\n")
    for (malla_nombre in names(config_mallas)) {
      cfg <- config_mallas[[malla_nombre]]
      cat(sprintf("  === AR1 | Malla: %s ===\n", cfg$label))
      res_ar1[[malla_nombre]] <- ajustar_evaluar("espacio-temporal", cfg, sel, datos)
      gc()
    }
    rmse_ar1        <- sapply(res_ar1, function(r) ifelse(is.na(r$metricas$RMSE), Inf, r$metricas$RMSE))
    malla_mejor_ar1 <- names(which.min(rmse_ar1))
    modelo_ar1_gan  <- res_ar1[[malla_mejor_ar1]]$modelo
    spde_ar1_gan    <- res_ar1[[malla_mejor_ar1]]$spde
    stk_gan         <- res_ar1[[malla_mejor_ar1]]$stk
    tag_gan         <- res_ar1[[malla_mejor_ar1]]$tag

    cat("  [B5] Diagnostico de residuos del modelo AR1...\n")
    etiq_ar1 <- paste0("ar1_", gsub("[^A-Za-z0-9]", "_", malla_mejor_ar1))
    diagnostico_residuos(
      residuos = res_ar1[[malla_mejor_ar1]]$residuos, datos = datos, sel = sel,
      dir_acf  = subs$ACF,
      dir_var  = subs$variograma,
      dir_qq   = subs$QQ,
      etiqueta = etiq_ar1
    )
  } else {
    cat("  [B4] Sin autocorrelacion -> modelo espacial suficiente.\n")
  }

  # ── B6: Comparacion y seleccion ──────────────────────────────────────────
  cat("  [B6] Comparando modelos...\n")
  metricas_esp <- rbindlist(lapply(res_espacial, `[[`, "metricas"))
  metricas_ar1 <- if (length(res_ar1) > 0)
                    rbindlist(lapply(res_ar1, `[[`, "metricas")) else data.table()

  tabla_comp <- rbindlist(list(metricas_esp, metricas_ar1), fill = TRUE)
  setorder(tabla_comp, RMSE)

  guardar_tabla_png(
    tabla_comp,
    titulo    = sprintf("Comparacion de Modelos INLA-SPDE - %s", etiqueta),
    subtitulo = sprintf("Train: %s | Test: %s", datos$lab_train, datos$lab_test),
    ruta_png  = file.path(subs$comparacion, "tabla_comparacion_mallas_modelos.png"),
    ancho_px  = 1100
  )

  niveles_malla <- c("Muy gruesa (12 km)", "Gruesa (8 km)", "Media (4 km)", "Fina (1 km)")
  tabla_comp[, Malla := factor(Malla, levels = intersect(niveles_malla, unique(Malla)))]

  ggplot(tabla_comp, aes(x = Malla, y = RMSE, fill = Modelo)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%.4f", RMSE)),
              position = position_dodge(width = 0.7), vjust = -0.4, size = 3.2) +
    scale_fill_manual(values = c("Solo espacial" = "#B2182B",
                                  "Espacio-temporal (AR1)" = "#2166AC"),
                      na.value = "grey70") +
    labs(title    = sprintf("RMSE - %s", etiqueta),
         subtitle = sprintf("Test: %s", datos$lab_test),
         x = "Resolucion de malla", y = "RMSE (log NO2)", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  ggsave(file.path(subs$comparacion, "grafico_rmse.png"),
         width = 10, height = 6, dpi = 300)

  ggplot(tabla_comp, aes(x = Malla, y = Tiempo_min, fill = Modelo)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = paste0(Tiempo_min, " min")),
              position = position_dodge(width = 0.7), vjust = -0.4, size = 3.2) +
    scale_fill_manual(values = c("Solo espacial" = "#B2182B",
                                  "Espacio-temporal (AR1)" = "#2166AC"),
                      na.value = "grey70") +
    labs(title = sprintf("Tiempo de computo - %s", etiqueta),
         x = "Resolucion de malla", y = "Tiempo (min)", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  ggsave(file.path(subs$comparacion, "grafico_tiempo.png"),
         width = 10, height = 6, dpi = 300)

  if (!is.null(modelo_ar1_gan) && !is.null(res_esp_gan$modelo)) {
    comparar_significancia(
      modelo_esp   = res_esp_gan$modelo,
      modelo_ar1   = modelo_ar1_gan,
      vars_finales = sel$vars_finales,
      carpeta_comp = subs$comparacion,
      etiqueta     = etiqueta
    )
  }

  # ── B7: Tiempo de recuperacion ────────────────────────────────────────────
  cat("  [B7] Calculando tiempo de recuperacion...\n")
  calcular_recuperacion(
    modelo_ar1  = modelo_ar1_gan,
    spde_obj    = spde_ar1_gan,
    datos       = datos,
    sel         = sel,
    stk_gan     = stk_gan,
    tag_gan     = tag_gan,
    carpeta_rec = subs$recuperacion,
    etiqueta    = etiqueta
  )

  tabla_comp[, Frecuencia := toupper(freq)]
  tabla_comp[, Periodo    := toupper(periodo)]
  resumen_global[[length(resumen_global) + 1]] <- tabla_comp

  cat(sprintf("  %s completado. Resultados en: %s\n", etiqueta, carpeta))
}

# ==============================================================================
# TABLA RESUMEN GLOBAL
# ==============================================================================

tabla_global <- rbindlist(resumen_global, fill = TRUE)
setorder(tabla_global, Frecuencia, Periodo, RMSE)

mejor_por_config <- tabla_global[, .SD[which.min(RMSE)], by = .(Frecuencia, Periodo)]

guardar_tabla_png(
  mejor_por_config[, .(Frecuencia, Periodo, Malla, Modelo,
                        RMSE, MAE, Cov95, Rango_km, Sigma2, Tiempo_min)],
  titulo    = "Mejor Configuracion por Frecuencia y Periodo",
  subtitulo = "Modelo con menor RMSE en test por frecuencia y periodo",
  ruta_png  = here("outputs", "simulacion", "tabla_resumen_global.png"),
  ancho_px  = 1100
)

guardar_tabla_png(
  tabla_global[, .(Frecuencia, Periodo, Malla, Modelo, n_nodos,
                   DIC, WAIC, RMSE, MAE, Cov95, Tiempo_min)],
  titulo    = "Resultados Completos - Todas las Frecuencias y Periodos",
  subtitulo = "B3 decide si AR1 | Variables independientes por frecuencia",
  ruta_png  = here("outputs", "simulacion", "tabla_resultados_completos.png"),
  ancho_px  = 1300
)

cat("\n", strrep("=", 80), "\n")
cat("  SIMULACION COMPLETADA\n")
cat(strrep("=", 80), "\n")
cat("  Resultados guardados en:", here("outputs", "simulacion"), "\n")

