# ==============================================================================
# SIMULACIÓN INLA-SPDE: BLOQUES 0–7 POR FRECUENCIA (MENSUAL → DIARIO → HORARIO)
# ==============================================================================
# Objetivo:
#   - Comparar 3 frecuencias temporales: MENSUAL (2019-2025), DIARIA, HORARIA
#   - Mensual: identificar covariables MACRO (efecto estacional/interanual)
#   - Diario/Horario: identificar covariables MICRO (efecto a corto plazo)
#   - Contraste estacional diario/horario: INVIERNO (ene+feb) vs VERANO (jul+ago)
#   - Para cada frecuencia (Bloques 0–7):
#       • B0:   Preparación de datos
#       • B0.5: Selección por residuos parciales (Y - beta*X iterativo)
#       • B1:   Selección de covariables (VIF + significancia + stepwise DIC)
#       • B2:   Modelo espacial SPDE (4 mallas)
#       • B3:   Diagnóstico de residuos del modelo espacial (ACF + QQ)
#       • B4:   Modelo espacio-temporal AR1 (si B3 lo justifica)
#       • B5:   Diagnóstico de residuos del AR1
#       • B6:   Comparación y selección de modelos
#       • B6.5: Tablas de diagnóstico (efectos fijos, hiperparámetros, SPDE)
#       • B7:   Tiempo de recuperación de zonas contaminadas
#       • B8:   (solo horario) Modelo de residuos diario -> horario:
#               (1) modelo DIARIO sobre el año completo -> beta_hat de las
#               covariables significativas; (2) pseudo-respuesta horaria
#               y*_sh = y_sh - sum(beta_hat_k * x_k,sh) regresada sobre las
#               covariables NO significativas del diario + retardo h-1 del
#               tráfico, con filtro VIF previo entre ellas (sin correlación
#               entre covariables), para identificar variables importantes
#               que hayan quedado en los residuos
#   - NOTA: No se calcula variograma empírico (24 estaciones = 276 pares,
#     insuficientes). El SPDE estima rango y varianza vía Matérn.
#   - Las variables candidatas se definen de forma independiente por frecuencia
#     en la sección CONFIGURACIÓN GENERAL (COVS_DIARIO/HORARIO_NOMBRES).
#   - Outputs organizados en subcarpetas por bloque:
#       seleccion/ | ACF/ | QQ/ | comparacion/ | diagnosticos/ | recuperacion/
# ==============================================================================

library(INLA)
library(fmesher)
library(data.table)
library(sf)
library(here)
library(ggplot2)
library(car)
library(gridExtra)
library(grid)
library(gstat)

set.seed(4827)

# ==============================================================================
# CONFIGURACIÓN GENERAL
# ==============================================================================

VIF_UMBRAL <- 5
DIC_MEJORA_MIN <- 2

# Perfil horario mensual: meses consecutivos desde enero y holdout al final.
N_MESES_PERFIL <- 9L
N_MESES_TEST_PERFIL <- 2L

config_mallas <- list(
  muy_gruesa = list(max.edge = c(12, 18), cutoff = 1.0, label = "Muy gruesa (12 km)"),
  gruesa     = list(max.edge = c(8, 12), cutoff = 0.5, label = "Gruesa (8 km)"),
  media      = list(max.edge = c(4, 8), cutoff = 0.5, label = "Media (4 km)"),
  fina       = list(max.edge = c(1, 4), cutoff = 0.25, label = "Fina (1 km)")
)

crear_spde <- function(mesh) {
  inla.spde2.pcmatern(
    mesh = mesh, alpha = 2,
    prior.range = c(5, 0.5),
    prior.sigma = c(1, 0.01)
  )
}

# ==============================================================================
# COVARIABLES CANDIDATAS POR FRECUENCIA
# Modificar aquí para cambiar las variables de cada escala temporal.
# Nombres : columnas estandarizadas presentes en el dataset.
# Alias   : nombres internos que usará INLA (sin espacios ni caracteres especiales).
# Mantener el mismo orden en Nombres y Alias.
# ==============================================================================

# --- MENSUAL (2019-2025, 7 anos): todas las covariables disponibles ---
# A escala mensual podemos incluir mas covariables porque la agregacion
# reduce el ruido y la multicolinealidad. Esto permite identificar
# cuales son "macro" (significativas a escala mensual/estacional)
# vs "micro" (solo significativas a escala diaria/horaria).
COVS_MENSUAL_NOMBRES <- c(
  "intensidad", "carga", "Temperatura",
  "Humedad_Relativa", "Precipitaciones",
  "Presion Barométrica", "Radiación Solar",
  "Velocidad Viento"
)
COVS_MENSUAL_ALIAS <- c(
  "trafico_intensidad", "trafico_carga", "temperatura",
  "humedad_relativa", "precipitacion",
  "presion_barometrica", "radiacion_solar",
  "velocidad_viento"
)

# --- DIARIO y HORARIO (2025, 1 ano) ---
COVS_DIARIO_NOMBRES <- c(
  "intensidad", "Temperatura",
  "Precipitaciones", "Velocidad Viento"
)
COVS_DIARIO_ALIAS <- c(
  "trafico_intensidad", "temperatura",
  "precipitacion", "velocidad_viento"
)

COVS_HORARIO_NOMBRES <- c(
  "intensidad", "Temperatura",
  "Precipitaciones", "Velocidad Viento"
)
COVS_HORARIO_ALIAS <- c(
  "trafico_intensidad", "temperatura",
  "precipitacion", "velocidad_viento"
)

# ==============================================================================
# FUNCIONES AUXILIARES
# ==============================================================================

# --- Guardar tabla como PNG (sin navegador: gridExtra + ggsave) ---
guardar_tabla_png <- function(df, titulo, subtitulo = NULL, ruta_png,
                              ancho_px = 800) {
  # Convertir a data.frame plano si es data.table
  df <- as.data.frame(df)

  # Construir el grob de la tabla con estilo mínimo
  tbl_theme <- gridExtra::ttheme_minimal(
    core = list(
      fg_params = list(fontsize = 9),
      bg_params = list(fill = c("white", "#F5F5F5"), col = NA)
    ),
    colhead = list(
      fg_params = list(fontsize = 9, fontface = "bold"),
      bg_params = list(fill = "#DDEEFF", col = NA)
    )
  )
  tbl_grob <- gridExtra::tableGrob(df, rows = NULL, theme = tbl_theme)

  # Título y subtítulo como grobs de texto
  titulo_grob <- grid::textGrob(
    titulo,
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )

  if (!is.null(subtitulo)) {
    sub_grob <- grid::textGrob(
      subtitulo,
      gp = grid::gpar(fontsize = 9, col = "grey40")
    )
    combinado <- gridExtra::arrangeGrob(
      titulo_grob, sub_grob, tbl_grob,
      nrow = 3,
      heights = grid::unit(c(0.55, 0.35, nrow(df) * 0.35 + 0.5), "inches")
    )
    alto_in <- 0.55 + 0.35 + nrow(df) * 0.35 + 0.5
  } else {
    combinado <- gridExtra::arrangeGrob(
      titulo_grob, tbl_grob,
      nrow = 2,
      heights = grid::unit(c(0.55, nrow(df) * 0.35 + 0.5), "inches")
    )
    alto_in <- 0.55 + nrow(df) * 0.35 + 0.5
  }

  ancho_in <- max(6, ancho_px / 96)
  alto_in <- max(3, alto_in)

  ggplot2::ggsave(
    filename = ruta_png,
    plot     = combinado,
    width    = ancho_in,
    height   = alto_in,
    dpi      = 96,
    bg       = "white"
  )
  cat("  Tabla guardada:", basename(ruta_png), "\n")
}

# --- Preparar datos según frecuencia y periodo estacional ---
# Periodos: "invierno" (ene+feb) vs "verano" (jul+ago)
# Contraste estacional: invierno = alta contaminacion (inversiones termicas,
# calefaccion, menor dispersion) vs verano = baja (mayor fotolisis, dispersion)
preparar_datos <- function(frecuencia, periodo = "invierno") {
  if (frecuencia == "horario") {
    dt <- readRDS(here(
      "data", "processed", "Maestro", "horario",
      "dataset_maestro_inla_2025_HORARIO.rds"
    ))
    setDT(dt)
    setnames(dt, "LOG_NO2_HORARIO", "LOG_NO2")

    todas_fechas <- sort(unique(dt$FECHA))

    # 5 dias totales: 4 de entrenamiento + 1 de test
    N_DIAS_HORARIO <- 5L
    N_DIAS_TEST_H <- 1L
    if (periodo == "invierno") {
      fechas_mes <- todas_fechas[format(todas_fechas, "%m") %in% c("01", "02")]
      n_dias_disp <- length(fechas_mes)
      n_dias_sel <- min(N_DIAS_HORARIO, n_dias_disp)
      fechas_sel <- fechas_mes[seq_len(n_dias_sel)]
    } else if (periodo == "verano") {
      fechas_mes <- todas_fechas[format(todas_fechas, "%m") %in% c("07", "08")]
      n_dias_disp <- length(fechas_mes)
      n_dias_sel <- min(N_DIAS_HORARIO, n_dias_disp)
      fechas_sel <- fechas_mes[seq_len(n_dias_sel)]
    }

    dt <- dt[FECHA %in% fechas_sel]
    dt[, ID_TIEMPO_SIM := as.integer(factor(
      paste(FECHA, sprintf("%02d", HORA))
    ))]
    setorder(dt, ID_TIEMPO_SIM, ESTACION)
    n_periodos <- uniqueN(dt$ID_TIEMPO_SIM)
    # 4 dias de test (96 horas) al final de la ventana
    n_test_p <- min(N_DIAS_TEST_H * 24L, floor(n_periodos * 0.30))
    n_train_p <- n_periodos - n_test_p
    periodos <- sort(unique(dt$ID_TIEMPO_SIM))
    periodos_train <- periodos[seq_len(n_train_p)]
    periodos_test <- periodos[(n_train_p + 1L):n_periodos]
    dt[, es_train := ID_TIEMPO_SIM %in% periodos_train]
    fecha_corte_train <- fechas_sel[ceiling(n_train_p / 24)]
    fecha_corte_test <- fechas_sel[ceiling(n_train_p / 24) + 1]
    lab_train <- sprintf(
      "%s → %s (%d horas)",
      min(fechas_sel), fecha_corte_train, n_train_p
    )
    lab_test <- sprintf(
      "%s → %s (%d horas)",
      fecha_corte_test, max(fechas_sel), n_test_p
    )
  } else if (frecuencia == "horario_perfil") {
    # Perfil medio de 24 horas por estacion, mes y tipo de dia.
    dt <- readRDS(here(
      "data", "processed", "Maestro", "horario",
      "dataset_maestro_inla_2025_HORARIO.rds"
    ))
    setDT(dt)
    setnames(dt, "LOG_NO2_HORARIO", "LOG_NO2")

    if (N_MESES_PERFIL < 3L) {
      stop("N_MESES_PERFIL debe ser al menos 3.")
    }
    if (N_MESES_TEST_PERFIL < 1L || N_MESES_TEST_PERFIL >= N_MESES_PERFIL) {
      stop("N_MESES_TEST_PERFIL debe estar entre 1 y N_MESES_PERFIL - 1.")
    }

    dt[, MES := as.integer(format(FECHA, "%m"))]
    meses_disponibles <- sort(unique(dt$MES))
    meses_seleccionados <- head(meses_disponibles, N_MESES_PERFIL)
    if (length(meses_seleccionados) < N_MESES_PERFIL) {
      stop(sprintf(
        "Se solicitaron %d meses para horario_perfil, pero solo hay %d.",
        N_MESES_PERFIL, length(meses_seleccionados)
      ))
    }
    dt <- dt[MES %in% meses_seleccionados]
    dt[, dow := as.integer(format(FECHA, "%u"))]
    dt[, tipo_dia := fifelse(dow <= 5, "laborable", "finde")]

    covs_cols <- intersect(COVS_HORARIO_NOMBRES, names(dt))
    dt <- dt[!is.na(LOG_NO2), c(
      list(
        LOG_NO2 = mean(LOG_NO2, na.rm = TRUE),
        LONGITUD = LONGITUD[1L],
        LATITUD = LATITUD[1L]
      ),
      lapply(.SD, function(x) mean(x, na.rm = TRUE))
    ), by = .(ESTACION, MES, tipo_dia, HORA), .SDcols = covs_cols]

    dt[, tipo_ord := fifelse(tipo_dia == "laborable", 1L, 2L)]
    clave <- unique(dt[, .(MES, tipo_ord, HORA)])
    setorder(clave, MES, tipo_ord, HORA)
    clave[, ID_TIEMPO_SIM := .I]
    dt <- merge(dt, clave, by = c("MES", "tipo_ord", "HORA"))
    dt <- dt[complete.cases(dt[, c("LOG_NO2", covs_cols), with = FALSE])]
    setorder(dt, ID_TIEMPO_SIM, ESTACION)

    meses_test <- tail(meses_seleccionados, N_MESES_TEST_PERFIL)
    meses_train <- setdiff(meses_seleccionados, meses_test)
    periodos_train <- sort(unique(dt[MES %in% meses_train, ID_TIEMPO_SIM]))
    periodos_test <- sort(unique(dt[MES %in% meses_test, ID_TIEMPO_SIM]))
    n_train_p <- length(periodos_train)
    n_test_p <- length(periodos_test)
    dt[, es_train := ID_TIEMPO_SIM %in% periodos_train]

    # Un unico campo espacio-temporal compartido de 24 horas. Esto evita crear
    # 432 campos y no concatena meses/tipos de dia en una falsa serie temporal.
    horas_perfil <- sort(unique(dt$HORA))
    if (length(horas_perfil) != 24L) {
      stop(sprintf("El perfil debe contener 24 horas; se encontraron %d.", length(horas_perfil)))
    }
    dt[, ID_GRUPO_MODELO := match(HORA, horas_perfil)]
    dt[, tipo_finde := as.numeric(tipo_dia == "finde")]
    dt[, hora_sin24 := sin(2 * pi * as.numeric(HORA) / 24)]
    dt[, hora_cos24 := cos(2 * pi * as.numeric(HORA) / 24)]
    dt[, hora_sin12 := sin(2 * pi * as.numeric(HORA) / 12)]
    dt[, hora_cos12 := cos(2 * pi * as.numeric(HORA) / 12)]
    dt[, finde_sin24 := tipo_finde * hora_sin24]
    dt[, finde_cos24 := tipo_finde * hora_cos24]
    dt[, finde_sin12 := tipo_finde * hora_sin12]
    dt[, finde_cos12 := tipo_finde * hora_cos12]
    media_mes_train <- mean(dt[es_train == TRUE, MES])
    sd_mes_train <- sd(dt[es_train == TRUE, MES])
    dt[, mes_tendencia := (MES - media_mes_train) / sd_mes_train]

    efectos_base <- list(
      tipo_finde = dt$tipo_finde,
      mes_tendencia = dt$mes_tendencia,
      finde_sin24 = dt$finde_sin24,
      finde_cos24 = dt$finde_cos24,
      finde_sin12 = dt$finde_sin12,
      finde_cos12 = dt$finde_cos12
    )

    n_periodos <- uniqueN(dt$ID_TIEMPO_SIM)
    lab_train <- sprintf(
      "Meses %s | laborable+finde x 24h (%d perfiles-hora)",
      paste(range(meses_train), collapse = "-"), n_train_p
    )
    lab_test <- sprintf(
      "Meses %s | laborable+finde x 24h (%d perfiles-hora)",
      paste(range(meses_test), collapse = "-"), n_test_p
    )
  } else if (frecuencia == "diario") {
    dt <- readRDS(here(
      "data", "processed", "Maestro", "diario",
      "dataset_maestro_inla_2025_DIARIO.rds"
    ))
    setDT(dt)
    setnames(dt, "LOG_NO2_DIARIO", "LOG_NO2")

    todas_fechas <- sort(unique(dt$FECHA))

    # 1 mes train + 1 semana test (reducido para memoria)
    if (periodo == "invierno") {
      fechas_mes <- todas_fechas[format(todas_fechas, "%m") %in% c("01", "02")]
      fechas_sel <- fechas_mes[seq_len(min(35L, length(fechas_mes)))]
    } else if (periodo == "verano") {
      fechas_mes <- todas_fechas[format(todas_fechas, "%m") %in% c("07", "08")]
      fechas_sel <- fechas_mes[seq_len(min(35L, length(fechas_mes)))]
    }

    dt <- dt[FECHA %in% fechas_sel]
    dt[, ID_TIEMPO_SIM := match(FECHA, fechas_sel)]
    setorder(dt, ID_TIEMPO_SIM, ESTACION)
    n_periodos <- length(fechas_sel)
    n_test_p <- 7L
    n_train_p <- n_periodos - n_test_p
    periodos_train <- 1L:n_train_p
    periodos_test <- (n_train_p + 1L):n_periodos
    dt[, es_train := ID_TIEMPO_SIM %in% periodos_train]
    lab_train <- sprintf(
      "%s → %s (%d días)",
      fechas_sel[1], fechas_sel[n_train_p], n_train_p
    )
    lab_test <- sprintf(
      "%s → %s (%d días)",
      fechas_sel[n_train_p + 1], fechas_sel[n_periodos],
      n_test_p
    )
  } else if (frecuencia == "mensual") {
    # Dataset mensual multi-anual (2019-2025)
    # Train: 2019-2024 (~72 meses) | Test: 2025 (12 meses)
    dt <- readRDS(here(
      "data", "processed", "Maestro", "mensual",
      "dataset_maestro_inla_2019_2025_MENSUAL.rds"
    ))
    setDT(dt)

    # Eliminar meses con NO2 = NA
    dt <- dt[!is.na(LOG_NO2)]

    fechas_disp <- sort(unique(dt$FECHA))
    dt[, ID_TIEMPO_SIM := match(FECHA, fechas_disp)]
    setorder(dt, ID_TIEMPO_SIM, ESTACION)

    # Train: 2022-2024 | Test: 2025 (reducido de 2019-2024 para memoria)
    dt <- dt[year(FECHA) >= 2022]
    fechas_disp <- sort(unique(dt$FECHA))
    dt[, ID_TIEMPO_SIM := match(FECHA, fechas_disp)]
    setorder(dt, ID_TIEMPO_SIM, ESTACION)

    n_periodos <- length(fechas_disp)
    fechas_test <- fechas_disp[year(fechas_disp) == 2025]
    n_test_p <- length(fechas_test)
    n_train_p <- n_periodos - n_test_p

    periodos_train <- seq_len(n_train_p)
    periodos_test <- (n_train_p + 1L):n_periodos
    dt[, es_train := ID_TIEMPO_SIM %in% periodos_train]

    lab_train <- sprintf(
      "%s → %s (%d meses, 2019-2024)",
      fechas_disp[1], fechas_disp[n_train_p], n_train_p
    )
    lab_test <- sprintf(
      "%s → %s (%d meses, 2025)",
      fechas_disp[n_train_p + 1], fechas_disp[n_periodos],
      n_test_p
    )
  }

  if (!exists("efectos_base", inherits = FALSE)) efectos_base <- list()
  if (!"ID_GRUPO_MODELO" %in% names(dt)) dt[, ID_GRUPO_MODELO := ID_TIEMPO_SIM]

  # Coordenadas UTM en km
  coords_unicas <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
  coords_sf <- st_as_sf(coords_unicas,
    coords = c("LONGITUD", "LATITUD"), crs = 4326
  ) |>
    st_transform(25830)
  coords_unicas[, X_km := st_coordinates(coords_sf)[, 1] / 1000]
  coords_unicas[, Y_km := st_coordinates(coords_sf)[, 2] / 1000]
  dt <- merge(dt, coords_unicas[, .(ESTACION, X_km, Y_km)],
    by = "ESTACION", all.x = TRUE
  )
  setorder(dt, ID_TIEMPO_SIM, ESTACION)

  # NO2 en escala original (exp del log) para análisis exploratorio
  dt[, NO2_raw := exp(LOG_NO2)]

  list(
    dt           = dt,
    n_periodos   = uniqueN(dt$ID_TIEMPO_SIM),
    n_train      = n_train_p,
    n_test       = n_test_p,
    n_grupos_modelo = uniqueN(dt$ID_GRUPO_MODELO),
    lab_train    = lab_train,
    lab_test     = lab_test,
    efectos_base = efectos_base,
    frecuencia   = frecuencia
  )
}

# --- Heatmap de correlación (datos no estandarizados) + Análisis temporal ---

# --- Crear subcarpetas de outputs por bloque ---
crear_subcarpetas <- function(carpeta_base) {
  subs <- list(
    seleccion    = file.path(carpeta_base, "seleccion"),
    ACF          = file.path(carpeta_base, "ACF"),
    QQ           = file.path(carpeta_base, "QQ"),
    comparacion  = file.path(carpeta_base, "comparacion"),
    recuperacion = file.path(carpeta_base, "recuperacion"),
    diagnosticos = file.path(carpeta_base, "diagnosticos"),
    modelos      = carpeta_base
  )
  lapply(subs, dir.create, showWarnings = FALSE, recursive = TRUE)
  subs
}

guardar_perfil_horario <- function(datos, carpeta) {
  if (!identical(datos$frecuencia, "horario_perfil")) return(invisible(NULL))

  perfil <- datos$dt[, .(
    NO2_media = mean(NO2_raw, na.rm = TRUE),
    NO2_sd = sd(NO2_raw, na.rm = TRUE),
    n = .N
  ), by = .(
    Muestra = fifelse(es_train, "Entrenamiento", "Prueba"),
    tipo_dia, HORA
  )]
  perfil[, `:=`(
    IC_2.5 = NO2_media - 1.96 * NO2_sd / sqrt(n),
    IC_97.5 = NO2_media + 1.96 * NO2_sd / sqrt(n)
  )]
  fwrite(perfil, file.path(carpeta, "perfil_observado_laborable_finde.csv"))

  p <- ggplot(perfil, aes(
    x = HORA, y = NO2_media, color = tipo_dia,
    fill = tipo_dia, group = tipo_dia
  )) +
    geom_ribbon(aes(ymin = IC_2.5, ymax = IC_97.5), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(laborable = "#2166AC", finde = "#B2182B")) +
    scale_fill_manual(values = c(laborable = "#2166AC", finde = "#B2182B")) +
    scale_x_continuous(breaks = seq(min(perfil$HORA), max(perfil$HORA), by = 3)) +
    facet_wrap(~ Muestra, ncol = 1, scales = "free_y") +
    labs(
      title = "Perfil horario medio de NO2",
      subtitle = "Comparacion entre dias laborables y fines de semana",
      x = "Hora", y = "NO2 medio (ug/m3)", color = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  ggsave(
    file.path(carpeta, "perfil_observado_laborable_finde.png"),
    p, width = 10, height = 8, dpi = 300
  )
  invisible(perfil)
}

# ==============================================================================
# BLOQUE 0.5 — SELECCIÓN POR RESIDUOS PARCIALES (INLA-SPDE)
# ==============================================================================
# Enfoque iterativo usando modelos INLA-SPDE univariantes:
#   En cada paso se ajusta un modelo INLA con campo espacial Matérn + una
#   covariable candidata. Se extrae el beta posterior medio y se resta su
#   contribución (Y <- Y - beta_SPDE * X). Se elige la covariable que más
#   reduce la SD del residuo. Los betas incorporan la estructura espacial,
#   evitando sesgos por confusión espacial que tendría una regresión OLS.
# Resultado: orden de importancia marginal y tabla de reducciones.
# ==============================================================================

seleccion_por_residuos_parciales <- function(datos, covs_nombres, covs_alias,
                                             carpeta_sel) {
  dt <- datos$dt
  alias_map <- setNames(covs_alias, covs_nombres)

  covs_disp <- intersect(covs_nombres, names(dt))
  todas_covs <- setNames(
    lapply(covs_disp, function(nm) dt[[nm]]),
    alias_map[covs_disp]
  )
  efectos_base <- datos$efectos_base
  rhs_rp <- paste(c("cov_test", names(efectos_base)), collapse = " + ")
  fml_rp <- as.formula(paste(
    "y ~ 0 + intercept +", rhs_rp,
    "+ f(campo_espacial, model = spde_rp)"
  ))

  vars_pool <- names(todas_covs)
  y_actual <- dt$LOG_NO2
  idx_train <- which(dt$es_train)

  # --- Construir malla y SPDE para residuos parciales (malla media) ---
  coords_puntos <- as.matrix(dt[, .(X_km, Y_km)])
  coords_matriz <- as.matrix(unique(dt[, .(X_km, Y_km)]))
  bnd_inner_rp <- inla.nonconvex.hull(coords_matriz, convex = -0.05, resolution = 50)
  bnd_outer_rp <- inla.nonconvex.hull(coords_matriz, convex = -0.2)
  malla_rp <- inla.mesh.2d(
    loc = coords_matriz,
    boundary = list(bnd_inner_rp, bnd_outer_rp),
    max.edge = c(4, 8), cutoff = 0.5
  )
  spde_rp <- crear_spde(malla_rp)
  indice_rp <- inla.spde.make.index("campo_espacial", n.spde = spde_rp$n.spde)
  A_rp <- inla.spde.make.A(mesh = malla_rp, loc = coords_puntos)

  # Respuesta enmascarada (solo train)
  y_train_mask <- ifelse(dt$es_train, y_actual, NA_real_)

  resultados <- list()
  vars_seleccionadas <- character(0)
  sd_original <- sd(y_actual[idx_train], na.rm = TRUE)

  cat(sprintf(
    "  [Residuos parciales INLA-SPDE] SD original de LOG_NO2 = %.4f\n",
    sd_original
  ))

  for (paso in seq_along(vars_pool)) {
    vars_restantes <- setdiff(vars_pool, vars_seleccionadas)
    if (length(vars_restantes) == 0) break

    # Probar cada covariable restante con modelo INLA-SPDE univariante
    errores <- sapply(vars_restantes, function(v) {
      x_v <- todas_covs[[v]]

      stk_rp <- inla.stack(
        tag = "rp", data = list(y = y_train_mask),
        A = list(A_rp, 1),
        effects = list(
          c(indice_rp, list(intercept = 1)),
          c(list(cov_test = x_v), efectos_base)
        ),
        compress = FALSE
      )

      mod_rp <- tryCatch(
        inla(
          formula = fml_rp,
          data = inla.stack.data(stk_rp, spde = spde_rp),
          family = "gaussian",
          control.predictor = list(A = inla.stack.A(stk_rp), compute = FALSE),
          control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
          control.inla = list(strategy = "gaussian"),
          verbose = FALSE
        ),
        error = function(e) NULL
      )

      if (is.null(mod_rp)) {
        return(Inf)
      }

      beta <- mod_rp$summary.fixed["cov_test", "mean"]
      residuo <- y_actual - beta * x_v
      sd(residuo[idx_train], na.rm = TRUE)
    })

    mejor_var <- names(which.min(errores))
    x_mejor <- todas_covs[[mejor_var]]

    # Ajustar INLA-SPDE con la mejor covariable para obtener beta posterior
    stk_mejor <- inla.stack(
      tag = "rp", data = list(y = y_train_mask),
      A = list(A_rp, 1),
      effects = list(
        c(indice_rp, list(intercept = 1)),
        c(list(cov_test = x_mejor), efectos_base)
      ),
      compress = FALSE
    )

    mod_mejor <- inla(
      formula = fml_rp,
      data = inla.stack.data(stk_mejor, spde = spde_rp),
      family = "gaussian",
      control.predictor = list(A = inla.stack.A(stk_mejor), compute = FALSE),
      control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
      control.inla = list(strategy = "gaussian"),
      verbose = FALSE
    )

    beta_mejor <- mod_mejor$summary.fixed["cov_test", "mean"]

    error_antes <- sd(y_actual[idx_train], na.rm = TRUE)
    y_actual <- y_actual - beta_mejor * x_mejor
    # Actualizar respuesta enmascarada para siguiente iteracion
    y_train_mask <- ifelse(dt$es_train, y_actual, NA_real_)
    error_despues <- sd(y_actual[idx_train], na.rm = TRUE)
    reduccion_pct <- (1 - error_despues / error_antes) * 100
    reduccion_acum <- (1 - error_despues / sd_original) * 100

    resultados[[paso]] <- data.frame(
      Paso = paso,
      Variable = mejor_var,
      Beta = round(beta_mejor, 4),
      SD_antes = round(error_antes, 4),
      SD_despues = round(error_despues, 4),
      Reduccion_pct = round(reduccion_pct, 2),
      Reduccion_acumulada_pct = round(reduccion_acum, 2)
    )

    vars_seleccionadas <- c(vars_seleccionadas, mejor_var)

    cat(sprintf(
      "  Paso %d: restar '%s' (beta_SPDE=%.4f) | SD: %.4f -> %.4f (%.1f%% paso | %.1f%% acum)\n",
      paso, mejor_var, beta_mejor, error_antes, error_despues,
      reduccion_pct, reduccion_acum
    ))

    # Liberar memoria entre pasos
    rm(mod_mejor, stk_mejor)
    gc(verbose = FALSE)
  }

  tabla_res <- do.call(rbind, resultados)

  guardar_tabla_png(
    tabla_res,
    titulo = "Seleccion por Residuos Parciales (INLA-SPDE)",
    subtitulo = sprintf(
      "Y_nuevo = Y - beta_SPDE*X en cada paso | SD original = %.4f",
      sd_original
    ),
    ruta_png = file.path(carpeta_sel, "tabla_residuos_parciales.png"),
    ancho_px = 1000
  )

  # Grafico de reduccion acumulada
  ggplot(tabla_res, aes(x = reorder(Variable, Paso), y = Reduccion_acumulada_pct)) +
    geom_col(fill = "#2166AC", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", Reduccion_acumulada_pct)),
      vjust = -0.4, size = 3.5
    ) +
    labs(
      title = "Reduccion acumulada de variabilidad por covariable",
      subtitle = "Orden de importancia marginal (betas de INLA-SPDE)",
      x = NULL, y = "Reduccion acumulada (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
  ggsave(file.path(carpeta_sel, "grafico_residuos_parciales.png"),
    width = 8, height = 5, dpi = 300
  )

  list(
    orden_variables    = vars_seleccionadas,
    tabla              = tabla_res,
    y_residual_final   = y_actual,
    sd_original        = sd_original
  )
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
  alias_map <- setNames(
    covs_alias[match(covs_disponibles, covs_nombres)],
    covs_disponibles
  )
  todas_covs <- setNames(
    lapply(covs_disponibles, function(nm) dt[[nm]]),
    alias_map[covs_disponibles]
  )
  efectos_base <- datos$efectos_base
  vars_base <- names(efectos_base)

  vars_pool <- names(todas_covs)

  # Pre-filtro varianza cero
  for (v in vars_pool) {
    sd_v <- sd(todas_covs[[v]][dt$es_train], na.rm = TRUE)
    if (is.na(sd_v) || sd_v < 1e-10) {
      cat(sprintf("    Descartada '%s' (sd ~ 0 en train)\n", v))
      vars_pool <- setdiff(vars_pool, v)
    }
  }
  todas_covs <- todas_covs[vars_pool]
  covs_candidatas <- todas_covs

  cat(sprintf("  Vars candidatas: %s\n", paste(vars_pool, collapse = ", ")))

  y_train_test <- ifelse(dt$es_train, dt$LOG_NO2, NA_real_)
  idx_test_filas <- which(!dt$es_train)
  idx_train_filas <- which(dt$es_train)
  y_test_real <- dt$LOG_NO2[idx_test_filas]

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
      fml <- as.formula(paste("LOG_NO2 ~", paste(vars_vif, collapse = " + ")))
      lm_fit <- lm(fml, data = df_vif)
      vif_vals <- car::vif(lm_fit)
      cat(sprintf("    Vars: %s\n", paste(vars_vif, collapse = ", ")))
      if (max(vif_vals) <= VIF_UMBRAL) {
        cat("    -> Fin VIF.\n")
        break
      }
      var_elim <- names(which.max(vif_vals))
      cat(sprintf("    -> Eliminando '%s' (VIF=%.3f)\n", var_elim, max(vif_vals)))
      vars_vif <- setdiff(vars_vif, var_elim)
    }
    if (!is.null(vif_vals)) {
      tabla_vif <- data.frame(
        Variable = names(vif_vals),
        VIF = round(as.numeric(vif_vals), 4), row.names = NULL
      )
      guardar_tabla_png(
        tabla_vif[order(tabla_vif$VIF, decreasing = TRUE), ],
        titulo = "Seleccion de Covariables - VIF",
        subtitulo = sprintf(
          "Umbral VIF = %g | %d variables retenidas",
          VIF_UMBRAL, length(vars_vif)
        ),
        ruta_png = file.path(carpeta_sel, "tabla_seleccion_vif.png")
      )
    }
  }
  cat(sprintf("  Variables tras VIF: %s\n", paste(vars_vif, collapse = ", ")))

  # --- SPDE de referencia (malla media) ---
  coords_puntos <- as.matrix(dt[, .(X_km, Y_km)])
  coords_matriz <- as.matrix(unique(dt[, .(X_km, Y_km)]))
  bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.05, resolution = 50)
  bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.2)
  malla_ref <- inla.mesh.2d(
    loc = coords_matriz, boundary = list(bnd_inner, bnd_outer),
    max.edge = c(4, 8), cutoff = 0.5
  )
  spde_ref <- crear_spde(malla_ref)
  indice_ref <- inla.spde.make.index("campo_espacial", n.spde = spde_ref$n.spde)
  A_ref <- inla.spde.make.A(mesh = malla_ref, loc = coords_puntos)

  formula_full <- as.formula(paste(
    "y ~ 0 + intercept +", paste(c(vars_base, vars_vif), collapse = " + "),
    "+ f(campo_espacial, model = spde_ref)"
  ))
  stk_full <- inla.stack(
    tag = "full", data = list(y = y_train_test), A = list(A_ref, 1),
    effects = list(
      c(indice_ref, list(intercept = 1)),
      c(efectos_base, todas_covs[vars_vif])
    ),
    compress = FALSE
  )
  cat("  Ajustando modelo completo (significancia bayesiana)...\n")
  modelo_full <- inla(
    formula = formula_full, data = inla.stack.data(stk_full, spde = spde_ref),
    family = "gaussian",
    control.predictor = list(A = inla.stack.A(stk_full), compute = TRUE),
    control.compute = list(dic = TRUE, waic = FALSE, cpo = FALSE),
    control.inla = list(strategy = "laplace"), verbose = FALSE
  )

  sf_full <- modelo_full$summary.fixed
  sf_vars <- sf_full[intersect(vars_vif, rownames(sf_full)), , drop = FALSE]
  sig_mask <- sf_vars[, "0.025quant"] > 0 | sf_vars[, "0.975quant"] < 0
  vars_sig <- rownames(sf_vars)[sig_mask]

  tabla_sig <- data.frame(
    Variable = rownames(sf_vars),
    Media = round(sf_vars$mean, 4),
    SD = round(sf_vars$sd, 4),
    Q2.5 = round(sf_vars[, "0.025quant"], 4),
    Q97.5 = round(sf_vars[, "0.975quant"], 4),
    Sig_95 = ifelse(sig_mask, "Si", "No"),
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
    fml_sw <- as.formula(paste(
      "y ~ 0 + intercept +", paste(c(vars_base, vars_mod), collapse = " + "),
      "+ f(campo_espacial, model = spde_ref)"
    ))
    stk_sw <- inla.stack(
      tag = "sw", data = list(y = y_train_test), A = list(A_ref, 1),
      effects = list(
        c(indice_ref, list(intercept = 1)),
        c(efectos_base, todas_covs[vars_mod])
      ),
      compress = FALSE
    )
    mod <- tryCatch(
      inla(
        formula = fml_sw, data = inla.stack.data(stk_sw, spde = spde_ref),
        family = "gaussian",
        control.predictor = list(A = inla.stack.A(stk_sw), compute = TRUE),
        control.compute = list(dic = TRUE, waic = FALSE, cpo = FALSE),
        control.inla = list(strategy = "laplace"), verbose = FALSE
      ),
      error = function(e) {
        message("  ERROR sw: ", e$message)
        NULL
      }
    )
    if (is.null(mod)) {
      return(list(DIC = Inf, RMSE = NA_real_, MAE = NA_real_))
    }
    idx_d <- inla.stack.index(stk_sw, tag = "sw")$data
    pred <- mod$summary.fitted.values$mean[idx_d][idx_test_filas]
    list(
      DIC = mod$dic$dic,
      RMSE = sqrt(mean((pred - y_test_real)^2, na.rm = TRUE)),
      MAE = mean(abs(pred - y_test_real), na.rm = TRUE)
    )
  }

  vars_sw_actual <- if (length(vars_sig) > 0) vars_sig else vars_vif
  res_actual <- ajustar_sw(vars_sw_actual)
  tabla_sw <- list()
  tabla_sw[[1]] <- data.frame(
    Iter = 0L, Accion = "inicial", Variable = "-",
    Variables = paste(vars_sw_actual, collapse = " + "),
    DIC = round(res_actual$DIC, 2), RMSE = round(res_actual$RMSE, 4),
    MAE = round(res_actual$MAE, 4), stringsAsFactors = FALSE
  )
  cat(sprintf(
    "  Stepwise inicial: %s (DIC=%.2f)\n",
    paste(vars_sw_actual, collapse = " + "), res_actual$DIC
  ))

  for (iter in seq_len(length(vars_vif) + 1L)) {
    cands <- data.frame(
      accion = character(), variable = character(),
      DIC = numeric(), RMSE = numeric(), MAE = numeric(),
      stringsAsFactors = FALSE
    )
    for (v in setdiff(vars_vif, vars_sw_actual)) {
      r <- ajustar_sw(c(vars_sw_actual, v))
      cands <- rbind(cands, data.frame(
        accion = "anadir", variable = v,
        DIC = r$DIC, RMSE = r$RMSE, MAE = r$MAE,
        stringsAsFactors = FALSE
      ))
    }
    if (length(vars_sw_actual) >= 2) {
      for (v in vars_sw_actual) {
        r <- ajustar_sw(setdiff(vars_sw_actual, v))
        cands <- rbind(cands, data.frame(
          accion = "eliminar", variable = v,
          DIC = r$DIC, RMSE = r$RMSE, MAE = r$MAE,
          stringsAsFactors = FALSE
        ))
      }
    }
    if (nrow(cands) == 0) break
    mejor <- cands[which.min(cands$DIC), ]
    mejora <- res_actual$DIC - mejor$DIC
    if (mejora < DIC_MEJORA_MIN) {
      cat(sprintf("  Stepwise convergido (DDIC=%.2f)\n", mejora))
      break
    }
    vars_sw_actual <- if (mejor$accion == "anadir") {
      c(vars_sw_actual, mejor$variable)
    } else {
      setdiff(vars_sw_actual, mejor$variable)
    }
    res_actual <- list(DIC = mejor$DIC, RMSE = mejor$RMSE, MAE = mejor$MAE)
    cat(sprintf(
      "  Iter %d: %s '%s' (DDIC=%.2f)\n",
      iter, mejor$accion, mejor$variable, mejora
    ))
    tabla_sw[[iter + 1]] <- data.frame(
      Iter = iter, Accion = mejor$accion, Variable = mejor$variable,
      Variables = paste(vars_sw_actual, collapse = " + "),
      DIC = round(mejor$DIC, 2), RMSE = round(mejor$RMSE, 4),
      MAE = round(mejor$MAE, 4), stringsAsFactors = FALSE
    )
  }

  guardar_tabla_png(
    do.call(rbind, tabla_sw),
    titulo = "Stepwise DIC - Seleccion de Covariables",
    subtitulo = sprintf(
      "DDIC min = %g | Vars finales: %s",
      DIC_MEJORA_MIN, paste(vars_sw_actual, collapse = ", ")
    ),
    ruta_png = file.path(carpeta_sel, "tabla_stepwise_dic.png"), ancho_px = 1000
  )
  cat(sprintf("  Variables finales: %s\n", paste(vars_sw_actual, collapse = ", ")))

  list(
    vars_finales    = c(vars_base, vars_sw_actual),
    vars_base       = vars_base,
    vars_vif        = vars_vif,
    vars_primarias  = vars_sw_actual,
    todas_covs      = c(efectos_base, todas_covs),
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
  dt <- datos$dt
  n_grupos <- datos$n_grupos_modelo
  if (any(dt$ID_GRUPO_MODELO < 1L | dt$ID_GRUPO_MODELO > n_grupos)) {
    stop("ID_GRUPO_MODELO contiene indices fuera de 1:n_grupos_modelo.")
  }

  covs_modelo <- sel$todas_covs[sel$vars_finales[sel$vars_finales %in% names(sel$todas_covs)]]

  malla <- inla.mesh.2d(
    loc = sel$coords_matriz, boundary = list(sel$bnd_inner, sel$bnd_outer),
    max.edge = cfg_malla$max.edge, cutoff = cfg_malla$cutoff
  )
  spde <- crear_spde(malla)
  ef <- paste("y ~ 0 + intercept +", paste(sel$vars_finales, collapse = " + "))

  if (tipo_modelo == "espacial") {
    indice <- inla.spde.make.index("campo_espacial", n.spde = spde$n.spde)
    A <- inla.spde.make.A(mesh = malla, loc = sel$coords_puntos)
    formula_mod <- as.formula(paste(
      ef, "+ f(campo_espacial, model = spde)"
    ))
    tag <- "sim_s"
  } else {
    indice <- inla.spde.make.index("campo_espacial",
      n.spde = spde$n.spde, n.group = n_grupos
    )
    A <- inla.spde.make.A(
      mesh = malla, loc = sel$coords_puntos,
      group = dt$ID_GRUPO_MODELO, n.group = n_grupos
    )
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
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
      control.inla = list(strategy = "laplace"), verbose = FALSE
    ),
    error = function(e) {
      message("  ERROR: ", e$message)
      NULL
    }
  )
  t_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (is.null(modelo)) {
    return(list(
      metricas = data.table(
        Malla = cfg_malla$label, Modelo = tipo_modelo,
        n_nodos = malla$n, DIC = NA_real_, WAIC = NA_real_,
        RMSE = NA_real_, MAE = NA_real_, Cov95 = NA_real_,
        Rango_km = NA_real_, Sigma2 = NA_real_,
        Tiempo_min = round(t_min, 2)
      ),
      modelo = NULL, residuos = NULL, stk = stk, tag = tag, malla = malla, spde = spde
    ))
  }

  idx_data <- inla.stack.index(stk, tag = tag)$data
  pred_all <- modelo$summary.fitted.values$mean[idx_data]
  pred_test <- pred_all[sel$idx_test_filas]
  pred_sd <- modelo$summary.fitted.values$sd[idx_data][sel$idx_test_filas]
  fila_precision <- grep(
    "Precision for the Gaussian observations",
    rownames(modelo$summary.hyperpar), value = TRUE
  )
  if (length(fila_precision) == 0L) {
    stop("No se encontro la precision gaussiana para calcular Cov95.")
  }
  var_observacion <- 1 / modelo$summary.hyperpar[fila_precision[1], "mean"]
  pred_sd_total <- sqrt(pred_sd^2 + var_observacion)

  residuos_train <- sel$y_train_real - pred_all[sel$idx_train_filas]

  rmse <- sqrt(mean((pred_test - sel$y_test_real)^2, na.rm = TRUE))
  mae <- mean(abs(pred_test - sel$y_test_real), na.rm = TRUE)
  cov95 <- mean(sel$y_test_real >= pred_test - 1.96 * pred_sd_total &
    sel$y_test_real <= pred_test + 1.96 * pred_sd_total, na.rm = TRUE)

  spde_res <- inla.spde.result(modelo, "campo_espacial", spde)
  rango_post <- exp(spde_res$summary.log.range.nominal$mean)
  sigma_post <- exp(spde_res$summary.log.variance.nominal$mean)
  etiq <- ifelse(tipo_modelo == "espacial", "Solo espacial", "Espacio-temporal (AR1)")

  cat(sprintf(
    "    [%s | %s] RMSE=%.4f MAE=%.4f Cov95=%.1f%% %.1fmin\n",
    cfg_malla$label, etiq, rmse, mae, cov95 * 100, t_min
  ))

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
                                 dir_acf, dir_qq,
                                 etiqueta = "", lag_max = 10) {
  # NOTA: No se calcula variograma empirico porque con 24 estaciones solo hay
  # 276 pares unicos, insuficientes para estimar la semivarianza de forma
  # estable. El modelo SPDE ya estima el rango y la varianza del campo
  # espacial via la covarianza de Matern, que es mas fiable en este caso.

  dt_train <- datos$dt[es_train == TRUE]

  # --- ACF + Ljung-Box (carpeta ACF/) ---
  id_diagnostico <- if (identical(datos$frecuencia, "horario_perfil")) {
    dt_train$ID_GRUPO_MODELO
  } else {
    dt_train$ID_TIEMPO_SIM
  }
  res_temporal <- tapply(residuos, id_diagnostico, mean, na.rm = TRUE)
  lag_uso <- min(lag_max, length(res_temporal) - 1L)

  png(file.path(dir_acf, paste0("acf_", etiqueta, ".png")),
    width = 800, height = 500
  )
  acf(res_temporal,
    lag.max = lag_uso,
    main = sprintf("ACF residuos - %s", etiqueta),
    col = "#2166AC", lwd = 2
  )
  dev.off()

  lj <- Box.test(res_temporal, lag = lag_uso, type = "Ljung-Box")
  acf_sig <- lj$p.value < 0.05
  cat(sprintf(
    "  [%s] Ljung-Box p=%.4f -> %s\n", etiqueta, lj$p.value,
    ifelse(acf_sig, "autocorrelacion DETECTADA (AR1 justificado)",
      "sin autocorrelacion significativa"
    )
  ))

  # --- Q-Q plot (carpeta QQ/) ---
  png(file.path(dir_qq, paste0("qq_", etiqueta, ".png")),
    width = 600, height = 600
  )
  qqnorm(residuos,
    main = sprintf("Q-Q residuos - %s", etiqueta),
    pch = 16, cex = 0.6, col = "#666666"
  )
  qqline(residuos, col = "#B2182B", lwd = 2)
  dev.off()

  invisible(list(acf_significativa = acf_sig, ljung_box_p = lj$p.value))
}

# ==============================================================================
# BLOQUE 8 — MODELO DE RESIDUOS DIARIO -> HORARIO (solo horario)
# ==============================================================================
# Esquema en dos etapas:
#   Etapa 1 (DIARIO, año completo): y_sd = beta_0 + beta_1*x1_sd + ...
#     Regresión sobre todas las covariables candidatas (filtro VIF previo).
#     Se retienen las significativas al 5% con sus beta_hat.
#   Etapa 2 (HORARIO, un mes): pseudo-respuesta
#     y*_sh = y_sh - sum_k beta_hat_k * x_k,sh
#     (se resta la contribución de las significativas del diario, evaluadas
#     en sus valores horarios) y se regresa y* sobre las covariables NO
#     significativas de la etapa 1 más el retardo h-1 de la intensidad de
#     tráfico. Filtro VIF previo entre las candidatas para que la correlación
#     entre ellas no infle ni enmascare ningún efecto.
#   Si alguna resulta significativa sobre y*, es una variable importante que
#   ha quedado en los residuos del modelo diario (dinámica intra-diaria).
# Ambas etapas usan las columnas *_raw (escala original): los beta_hat del
# diario solo pueden restarse a la escala horaria si las unidades coinciden
# (las columnas estandarizadas usan media/SD distintas en cada dataset).
# ==============================================================================

modelo_residuos_diario_horario <- function(mes_horario = "01",
                                           alpha_sig = 0.05) {
  carpeta <- here("outputs", "simulacion", "residuos_diario_horario")
  dir.create(carpeta, showWarnings = FALSE, recursive = TRUE)

  covs_raw <- c(
    intensidad_trafico  = "intensidad_raw",
    carga_trafico       = "carga_raw",
    temperatura         = "Temperatura_raw",
    humedad_relativa    = "Humedad_Relativa_raw",
    precipitacion       = "Precipitaciones_raw",
    presion_barometrica = "Presion Barométrica_raw",
    radiacion_solar     = "Radiación Solar_raw",
    velocidad_viento    = "Velocidad Viento_raw"
  )

  vif_backward_b8 <- function(df, respuesta, vars, titulo, ruta_png) {
    vif_vals <- NULL
    repeat {
      if (length(vars) < 2) break
      vif_vals <- car::vif(lm(reformulate(vars, respuesta), data = df))
      if (max(vif_vals) <= VIF_UMBRAL) break
      elim <- names(which.max(vif_vals))
      cat(sprintf(
        "  [B8] Eliminando '%s' por multicolinealidad (VIF=%.2f)\n",
        elim, max(vif_vals)
      ))
      vars <- setdiff(vars, elim)
    }
    if (!is.null(vif_vals)) {
      tabla_vif <- data.frame(
        Variable = names(vif_vals),
        VIF = round(as.numeric(vif_vals), 3), row.names = NULL
      )
      guardar_tabla_png(
        tabla_vif[order(tabla_vif$VIF, decreasing = TRUE), ],
        titulo = titulo,
        subtitulo = sprintf(
          "Umbral VIF = %g | %d variables retenidas",
          VIF_UMBRAL, length(vars)
        ),
        ruta_png = ruta_png
      )
    }
    vars
  }

  tabla_coefs_b8 <- function(modelo_lm) {
    cf <- summary(modelo_lm)$coefficients
    data.frame(
      Variable = rownames(cf),
      Beta = signif(cf[, "Estimate"], 4),
      SE = signif(cf[, "Std. Error"], 3),
      t = round(cf[, "t value"], 2),
      p_valor = signif(cf[, "Pr(>|t|)"], 3),
      Sig_5pct = ifelse(cf[, "Pr(>|t|)"] < alpha_sig, "Si", "No"),
      row.names = NULL, check.names = FALSE
    )
  }

  # ── Etapa 1: modelo DIARIO sobre el año completo ───────────────────────────
  dt_d <- readRDS(here(
    "data", "processed", "Maestro", "diario",
    "dataset_maestro_inla_2025_DIARIO.rds"
  ))
  setDT(dt_d)
  df1 <- data.table(y = dt_d$LOG_NO2_DIARIO)
  for (a in names(covs_raw)) df1[, (a) := dt_d[[covs_raw[[a]]]]]
  df1 <- na.omit(df1)
  cat(sprintf(
    "  [B8-Etapa1] Modelo diario: %d obs (%d estaciones x año completo)\n",
    nrow(df1), uniqueN(dt_d$ESTACION)
  ))

  vars1 <- vif_backward_b8(
    df1, "y", names(covs_raw),
    titulo = "B8 Etapa 1 - VIF covariables diarias",
    ruta_png = file.path(carpeta, "tabla_vif_etapa1_diario.png")
  )

  lm1 <- lm(reformulate(vars1, "y"), data = df1)
  tabla1 <- tabla_coefs_b8(lm1)
  r2_1 <- summary(lm1)$r.squared
  vars_sig <- setdiff(tabla1$Variable[tabla1$Sig_5pct == "Si"], "(Intercept)")
  betas_sig <- coef(lm1)[vars_sig]

  guardar_tabla_png(
    tabla1,
    titulo = "B8 Etapa 1 - Modelo diario (año completo)",
    subtitulo = sprintf(
      "log(NO2)_sd ~ covariables (escala original) | R2 = %.4f | Significativas: %s",
      r2_1, paste(vars_sig, collapse = ", ")
    ),
    ruta_png = file.path(carpeta, "tabla_etapa1_modelo_diario.png"),
    ancho_px = 1000
  )
  fwrite(tabla1, file.path(carpeta, "etapa1_modelo_diario.csv"))
  cat(sprintf(
    "  [B8-Etapa1] R2 = %.4f | Significativas: %s\n",
    r2_1, paste(vars_sig, collapse = ", ")
  ))

  # ── Etapa 2: pseudo-respuesta horaria sobre un mes ─────────────────────────
  dt_h <- readRDS(here(
    "data", "processed", "Maestro", "horario",
    "dataset_maestro_inla_2025_HORARIO.rds"
  ))
  setDT(dt_h)
  dt_h <- dt_h[format(FECHA, "%m") == mes_horario]

  df2 <- data.table(
    ESTACION = dt_h$ESTACION, FECHA = dt_h$FECHA, HORA = dt_h$HORA,
    y = dt_h$LOG_NO2_HORARIO
  )
  for (a in names(covs_raw)) df2[, (a) := dt_h[[covs_raw[[a]]]]]
  setorder(df2, ESTACION, FECHA, HORA)
  df2[, intensidad_trafico_lag1 := shift(intensidad_trafico, 1L), by = ESTACION]
  df2 <- na.omit(df2)
  cat(sprintf(
    "  [B8-Etapa2] Mes %s: %d obs horarias tras retardo y NA\n",
    mes_horario, nrow(df2)
  ))

  # Restar la contribución de las significativas del diario (beta_hat fijos)
  df2[, y_res := y]
  for (v in vars_sig) df2[, y_res := y_res - betas_sig[[v]] * get(v)]

  # Candidatas: no significativas en la etapa 1 + retardo del tráfico
  vars2 <- c(setdiff(names(covs_raw), vars_sig), "intensidad_trafico_lag1")
  cat(sprintf(
    "  [B8-Etapa2] Candidatas sobre residuos: %s\n",
    paste(vars2, collapse = ", ")
  ))

  vars2 <- vif_backward_b8(
    df2, "y_res", vars2,
    titulo = "B8 Etapa 2 - VIF covariables no significativas + lag tráfico",
    ruta_png = file.path(carpeta, "tabla_vif_etapa2_horario.png")
  )

  # --- Modelo de residuos mejorado -------------------------------------------
  # (1) Efectos fijos por estación: los offsets sistemáticos de nivel entre
  #     estaciones (heterogeneidad local, misalignment) no contaminan los
  #     coeficientes; cada efecto se estima con la variación DENTRO de cada
  #     estación.
  # (2) Errores estándar robustos por clúster (estación): los residuos
  #     horarios están autocorrelacionados; los SE clásicos de OLS serían
  #     optimistas. El clúster por estación admite autocorrelación arbitraria
  #     dentro de cada estación.
  # (3) Delta R2 parcial: pérdida de R2 al eliminar cada covariable
  #     (manteniendo el resto), como medida de importancia práctica.
  df2[, ESTACION := factor(ESTACION)]
  lm2 <- lm(reformulate(c(vars2, "ESTACION"), "y_res"), data = df2)
  r2_2 <- summary(lm2)$r.squared

  vcov_cl <- sandwich::vcovCL(lm2, cluster = df2$ESTACION)
  ct2 <- lmtest::coeftest(lm2, vcov. = vcov_cl)

  delta_r2 <- sapply(vars2, function(v) {
    r2_sin <- summary(lm(
      reformulate(c(setdiff(vars2, v), "ESTACION"), "y_res"), data = df2
    ))$r.squared
    r2_2 - r2_sin
  })

  sd_yres <- sd(df2$y_res)
  tabla2 <- data.frame(
    Variable = vars2,
    Beta = signif(ct2[vars2, "Estimate"], 4),
    Beta_std = round(
      ct2[vars2, "Estimate"] *
        sapply(vars2, function(v) sd(df2[[v]])) / sd_yres, 4
    ),
    SE_cluster = signif(ct2[vars2, "Std. Error"], 3),
    t = round(ct2[vars2, "t value"], 2),
    p_valor = signif(ct2[vars2, "Pr(>|t|)"], 3),
    Sig_5pct = ifelse(ct2[vars2, "Pr(>|t|)"] < alpha_sig, "Si", "No"),
    Delta_R2 = round(delta_r2[vars2], 4),
    row.names = NULL, check.names = FALSE
  )

  # Diagnóstico: ¿queda autocorrelación en los residuos del modelo de errores?
  df2[, res_modelo := residuals(lm2)]
  res_t <- df2[, .(res = mean(res_modelo)), by = .(FECHA, HORA)]
  setorder(res_t, FECHA, HORA)
  lj2 <- Box.test(res_t$res, lag = 24, type = "Ljung-Box")
  png(file.path(carpeta, "acf_etapa2_residuos.png"), width = 800, height = 500)
  acf(res_t$res, lag.max = 48,
      main = sprintf(
        "ACF residuos etapa 2 (media por hora) | Ljung-Box p = %.3g",
        lj2$p.value
      ),
      col = "#2166AC", lwd = 2)
  dev.off()

  guardar_tabla_png(
    tabla2,
    titulo = "B8 Etapa 2 - Residuos horarios vs covariables no significativas",
    subtitulo = sprintf(
      "y* ~ covariables + FE estacion | SE cluster por estacion | Mes %s | R2 = %.4f (FE incluidos) | Beta_std = betas por SD",
      mes_horario, r2_2
    ),
    ruta_png = file.path(carpeta, "tabla_etapa2_residuos_horario.png"),
    ancho_px = 1100
  )
  fwrite(tabla2, file.path(carpeta, "etapa2_residuos_horario.csv"))

  vars_alerta <- setdiff(
    tabla2$Variable[tabla2$Sig_5pct == "Si"], "(Intercept)"
  )
  if (length(vars_alerta) > 0) {
    cat(sprintf(
      "  [B8] Variables importantes que quedan en los residuos: %s (R2=%.4f)\n",
      paste(vars_alerta, collapse = ", "), r2_2
    ))
  } else {
    cat(sprintf(
      "  [B8] Ninguna covariable explica los residuos horarios (R2=%.4f).\n",
      r2_2
    ))
  }

  invisible(list(
    tabla_etapa1 = tabla1, r2_etapa1 = r2_1,
    vars_significativas_diario = vars_sig,
    betas_diario = betas_sig,
    tabla_etapa2 = tabla2, r2_etapa2 = r2_2,
    ljung_box_p_etapa2 = lj2$p.value,
    vars_en_residuos = vars_alerta
  ))
}

# ==============================================================================
# BLOQUE 6a — TABLAS DE DIAGNOSTICO DEL MODELO (efectos fijos, hiper, SPDE)
# ==============================================================================
# Genera y guarda tablas PNG con:
#   1. Efectos fijos (summary.fixed): coeficientes + IC 95%
#   2. Hiperparametros (summary.hyperpar): precisiones, rho AR1, etc.
#   3. Resultados fisicos del SPDE: rango espacial (km) y varianza del campo
# ==============================================================================

guardar_diagnosticos_modelo <- function(modelo, spde_obj, carpeta_diag,
                                        etiqueta, tipo_modelo) {
  if (is.null(modelo)) {
    cat(sprintf("  [Diagnosticos] Sin modelo %s para %s.\n", tipo_modelo, etiqueta))
    return(invisible(NULL))
  }

  prefijo <- gsub("[^A-Za-z0-9]", "_", tolower(tipo_modelo))

  # --- 1. Tabla de efectos fijos ---
  tabla_fijos <- modelo$summary.fixed
  df_fijos <- data.frame(
    Efecto = rownames(tabla_fijos),
    Media = round(tabla_fijos$mean, 4),
    SD = round(tabla_fijos$sd, 4),
    Q2.5 = round(tabla_fijos[, "0.025quant"], 4),
    Mediana = round(tabla_fijos[, "0.5quant"], 4),
    Q97.5 = round(tabla_fijos[, "0.975quant"], 4),
    Sig_95 = ifelse(tabla_fijos[, "0.025quant"] > 0 |
      tabla_fijos[, "0.975quant"] < 0, "Si", "No"),
    row.names = NULL, check.names = FALSE
  )

  guardar_tabla_png(
    df_fijos,
    titulo = sprintf("Efectos Fijos - %s | %s", tipo_modelo, etiqueta),
    subtitulo = "Coeficientes posteriores con IC 95% credible",
    ruta_png = file.path(
      carpeta_diag,
      sprintf("tabla_efectos_fijos_%s.png", prefijo)
    ),
    ancho_px = 900
  )
  fwrite(df_fijos, file.path(
    carpeta_diag, sprintf("efectos_fijos_%s.csv", prefijo)
  ))

  # --- 2. Tabla de hiperparametros ---
  tabla_hiper <- modelo$summary.hyperpar
  df_hiper <- data.frame(
    Hiperparametro = rownames(tabla_hiper),
    Media = round(tabla_hiper$mean, 4),
    SD = round(tabla_hiper$sd, 4),
    Q2.5 = round(tabla_hiper[, "0.025quant"], 4),
    Mediana = round(tabla_hiper[, "0.5quant"], 4),
    Q97.5 = round(tabla_hiper[, "0.975quant"], 4),
    row.names = NULL, check.names = FALSE
  )

  guardar_tabla_png(
    df_hiper,
    titulo = sprintf("Hiperparametros - %s | %s", tipo_modelo, etiqueta),
    subtitulo = "Incluye precision gaussiana, parametros SPDE y rho AR1 (si aplica)",
    ruta_png = file.path(
      carpeta_diag,
      sprintf("tabla_hiperparametros_%s.png", prefijo)
    ),
    ancho_px = 1000
  )
  fwrite(df_hiper, file.path(
    carpeta_diag, sprintf("hiperparametros_%s.csv", prefijo)
  ))

  # --- 3. Resultados fisicos del SPDE (rango y varianza) ---
  # Un fallo diagnostico no debe cancelar el resto de ejecuciones.
  df_spde <- tryCatch({
    spde_res <- inla.spde.result(modelo, "campo_espacial", spde_obj)
    sr <- spde_res$summary.log.range.nominal
    sv <- spde_res$summary.log.variance.nominal
    rango_mean <- exp(sr[1, "mean"])
    var_mean <- exp(sv[1, "mean"])
    rango_q <- exp(sr[1, c("0.025quant", "0.975quant")])
    var_q <- exp(sv[1, c("0.025quant", "0.975quant")])

    tabla <- data.frame(
      Parametro = c("Rango espacial (km)", "Varianza del campo (sill)"),
      Media = round(c(rango_mean, var_mean), 4),
      IC_2.5 = round(c(rango_q[[1]], var_q[[1]]), 4),
      IC_97.5 = round(c(rango_q[[2]], var_q[[2]]), 4),
      row.names = NULL
    )
    guardar_tabla_png(
      tabla,
      titulo = sprintf("Parametros Fisicos SPDE - %s | %s", tipo_modelo, etiqueta),
      subtitulo = "Rango: distancia de correlacion ~0.13 | Varianza: intensidad del campo",
      ruta_png = file.path(
        carpeta_diag,
        sprintf("tabla_spde_fisico_%s.png", prefijo)
      ),
      ancho_px = 900
    )
    fwrite(tabla, file.path(
      carpeta_diag, sprintf("spde_fisico_%s.csv", prefijo)
    ))
    tabla
  }, error = function(e) {
    warning(sprintf(
      "No se pudo extraer el diagnostico SPDE de %s: %s",
      tipo_modelo, conditionMessage(e)
    ))
    NULL
  })

  cat(sprintf("  [Diagnosticos] Tablas de %s guardadas en diagnosticos/\n", tipo_modelo))

  invisible(list(
    efectos_fijos = df_fijos,
    hiperparametros = df_hiper,
    spde_fisico = df_spde,
    rango_km = if (!is.null(df_spde)) df_spde$Media[1] else NA_real_,
    varianza_espacial = if (!is.null(df_spde)) df_spde$Media[2] else NA_real_
  ))
}

# ==============================================================================
# BLOQUE 6b — COMPARACION DE SIGNIFICANCIA (espacial vs AR1)
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
    mismo_signo <- if (!is.null(row_e) && !is.null(row_a)) {
      sign(row_e[, "mean"]) == sign(row_a[, "mean"])
    } else {
      NA
    }

    interp <- if (is.na(sig_e) || is.na(sig_a)) {
      "Solo en un modelo"
    } else if (sig_e && sig_a && !mismo_signo) {
      "Signo inestable"
    } else if (sig_e && sig_a) {
      "Efecto genuino"
    } else if (sig_e && !sig_a) {
      "Proxy autocorrelacion"
    } else if (!sig_e && sig_a) {
      "Emerge con AR1"
    } else {
      "No significativa"
    }

    data.frame(
      Variable = v,
      Coef_Espacial = if (!is.null(row_e)) round(row_e$mean, 4) else NA,
      Sig_Espacial = ifelse(is.na(sig_e), "---", ifelse(sig_e, "Si", "No")),
      Coef_AR1 = if (!is.null(row_a)) round(row_a$mean, 4) else NA,
      Sig_AR1 = if (is.null(sf_ar1)) {
        "---"
      } else {
        ifelse(is.na(sig_a), "---", ifelse(sig_a, "Si", "No"))
      },
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

  hp <- modelo_ar1$summary.hyperpar
  rho_row <- grep("GroupRho|Rho for", rownames(hp), value = TRUE)
  if (length(rho_row) == 0) {
    cat("  [B7] Parametro AR1 no encontrado en hyperpar.\n")
    return(invisible(NULL))
  }

  rho <- hp[rho_row[1], "mean"]
  rho_lo <- hp[rho_row[1], "0.025quant"]
  rho_hi <- hp[rho_row[1], "0.975quant"]
  tau_desde_rho <- function(x) -1 / log(pmax(abs(x), .Machine$double.eps))
  tau <- tau_desde_rho(rho)
  tau_lo <- tau_desde_rho(rho_lo)
  tau_hi <- tau_desde_rho(rho_hi)
  es_perfil <- identical(datos$frecuencia, "horario_perfil")
  unidad_tau <- if (es_perfil) "horas" else "periodos"

  cat(sprintf(
    "  [B7] rho = %.3f [%.3f, %.3f] | tau = %.2f %s [%.2f, %.2f]\n",
    rho, rho_lo, rho_hi, tau, unidad_tau, tau_lo, tau_hi
  ))

  guardar_tabla_png(
    data.frame(
      Parametro = c("rho (AR1)", sprintf("tau (%s hasta 1/e)", unidad_tau)),
      Media     = round(c(rho, tau), 3),
      IC_2.5    = round(c(rho_lo, tau_lo), 3),
      IC_97.5   = round(c(rho_hi, tau_hi), 3),
      row.names = NULL
    ),
    titulo = sprintf(
      "%s - %s",
      if (es_perfil) "Persistencia Horaria" else "Tiempo de Recuperacion",
      etiqueta
    ),
    subtitulo = sprintf(
      "tau = -1/log(rho): %s para decaer al 37%%",
      unidad_tau
    ),
    ruta_png = file.path(carpeta_rec, "tabla_recuperacion.png")
  )

  if (es_perfil) {
    return(invisible(list(rho = rho, tau = tau)))
  }

  # Mapa de hotspots
  idx_data <- inla.stack.index(stk_gan, tag = tag_gan)$data
  pred_mean <- modelo_ar1$summary.fitted.values$mean[idx_data][sel$idx_train_filas]
  dt_train <- datos$dt[es_train == TRUE]

  pred_dt <- data.table(
    ESTACION = dt_train$ESTACION,
    pred_log = pred_mean
  )
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
      title = sprintf("Hotspots de NO2 - %s", etiqueta),
      subtitle = sprintf("rho = %.3f | tau ~ %.1f periodos de recuperacion", rho, tau),
      x = "X (km UTM)", y = "Y (km UTM)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(carpeta_rec, "mapa_hotspots_recuperacion.png"),
    width = 9, height = 7, dpi = 300
  )

  invisible(list(rho = rho, tau = tau))
}

# ==============================================================================
# BUCLE PRINCIPAL — MENSUAL -> DIARIO -> HORARIO (Bloques 0-7)
# ==============================================================================
# Mensual (2019-2025): identifica covariables MACRO (efecto estacional/interanual)
# Diario/Horario (2025): identifica covariables MICRO (efecto a corto plazo)
# La tabla final compara significancia entre escalas para clasificar cada variable.
# ==============================================================================

ejecuciones <- list(
  # --- Estudio completo: descomenta para incluir otras escalas ---
  # list(freq = "mensual", periodo = "multianual"),
  # list(freq = "diario",  periodo = "invierno"),
  # list(freq = "diario",  periodo = "verano"),
  # list(freq = "horario", periodo = "invierno"),
  # Perfil horario agregado: ejecucion aislada para evitar bloqueos previos.
  list(
    freq = "horario_perfil",
    periodo = sprintf("2025_%dmeses", N_MESES_PERFIL)
  )
)

resumen_global <- list()
resumen_parciales <- list() # Para tabla macro/micro

for (ejec in ejecuciones) {
  freq <- ejec$freq
  periodo <- ejec$periodo
  etiqueta <- sprintf("%s - %s", toupper(freq), toupper(periodo))

  cat("\n", strrep("#", 80), "\n")
  cat(sprintf("  %s\n", etiqueta))
  cat(strrep("#", 80), "\n")

  carpeta <- here("outputs", "simulacion", freq, periodo)
  subs <- crear_subcarpetas(carpeta)

  # Covariables según frecuencia (modificar COVS_*_NOMBRES en config)
  covs_cfg <- switch(freq,
    mensual        = list(nombres = COVS_MENSUAL_NOMBRES, alias = COVS_MENSUAL_ALIAS),
    diario         = list(nombres = COVS_DIARIO_NOMBRES, alias = COVS_DIARIO_ALIAS),
    horario        = list(nombres = COVS_HORARIO_NOMBRES, alias = COVS_HORARIO_ALIAS),
    horario_perfil = list(nombres = COVS_HORARIO_NOMBRES, alias = COVS_HORARIO_ALIAS)
  )

  # ── B0: Preparar datos ────────────────────────────────────────────────────
  cat("  [B0] Preparando datos...\n")
  datos <- preparar_datos(freq, periodo)
  cat(sprintf(
    "  %d filas | %d estaciones | %d periodos | %d grupos latentes (train=%d, test=%d)\n",
    nrow(datos$dt), uniqueN(datos$dt$ESTACION),
    datos$n_periodos, datos$n_grupos_modelo, datos$n_train, datos$n_test
  ))
  guardar_perfil_horario(datos, carpeta)

  # ── B0.5: Residuos parciales (orden de importancia marginal) ────────────
  cat("  [B0.5] Seleccion por residuos parciales...\n")
  res_parciales <- seleccion_por_residuos_parciales(
    datos        = datos,
    covs_nombres = covs_cfg$nombres,
    covs_alias   = covs_cfg$alias,
    carpeta_sel  = subs$seleccion
  )
  cat(sprintf(
    "  Orden de importancia: %s\n",
    paste(res_parciales$orden_variables, collapse = " > ")
  ))

  # Guardar para tabla macro/micro final
  if (!is.null(res_parciales$tabla)) {
    tbl_parc <- as.data.table(res_parciales$tabla)
    tbl_parc[, Frecuencia := if (freq == "horario_perfil") "HORARIO" else toupper(freq)]
    tbl_parc[, Periodo := toupper(periodo)]
    resumen_parciales[[length(resumen_parciales) + 1]] <- tbl_parc
  }

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

  rmse_esp <- sapply(res_espacial, function(r) ifelse(is.na(r$metricas$RMSE), Inf, r$metricas$RMSE))
  if (!any(is.finite(rmse_esp))) {
    stop(sprintf("Ningun modelo espacial convergio para %s.", etiqueta))
  }
  malla_mejor_esp <- names(rmse_esp)[which.min(rmse_esp)]
  res_esp_gan <- res_espacial[[malla_mejor_esp]]

  # ── B3: Diagnostico residuos espacial ────────────────────────────────────
  cat("  [B3] Diagnostico de residuos del modelo espacial...\n")
  etiq_esp <- paste0("espacial_", gsub("[^A-Za-z0-9]", "_", malla_mejor_esp))
  diag_esp <- diagnostico_residuos(
    residuos = res_esp_gan$residuos, datos = datos, sel = sel,
    dir_acf = subs$ACF,
    dir_qq = subs$QQ,
    etiqueta = etiq_esp
  )

  # ── B4 + B5: AR1 (solo si B3 lo justifica) ───────────────────────────────
  res_ar1 <- list()
  modelo_ar1_gan <- NULL
  spde_ar1_gan <- res_esp_gan$spde
  stk_gan <- res_esp_gan$stk
  tag_gan <- res_esp_gan$tag
  malla_mejor_ar1 <- malla_mejor_esp

  if (diag_esp$acf_significativa) {
    cat("  [B4] Autocorrelacion detectada -> ajustando AR1 (4 mallas)...\n")
    for (malla_nombre in names(config_mallas)) {
      cfg <- config_mallas[[malla_nombre]]
      cat(sprintf("  === AR1 | Malla: %s ===\n", cfg$label))
      res_ar1[[malla_nombre]] <- ajustar_evaluar("espacio-temporal", cfg, sel, datos)
      gc()
    }
    rmse_ar1 <- sapply(res_ar1, function(r) ifelse(is.na(r$metricas$RMSE), Inf, r$metricas$RMSE))
    if (any(is.finite(rmse_ar1))) {
      malla_mejor_ar1 <- names(rmse_ar1)[which.min(rmse_ar1)]
      modelo_ar1_gan <- res_ar1[[malla_mejor_ar1]]$modelo
      spde_ar1_gan <- res_ar1[[malla_mejor_ar1]]$spde
      stk_gan <- res_ar1[[malla_mejor_ar1]]$stk
      tag_gan <- res_ar1[[malla_mejor_ar1]]$tag

      cat("  [B5] Diagnostico de residuos del modelo AR1...\n")
      etiq_ar1 <- paste0("ar1_", gsub("[^A-Za-z0-9]", "_", malla_mejor_ar1))
      diagnostico_residuos(
        residuos = res_ar1[[malla_mejor_ar1]]$residuos, datos = datos, sel = sel,
        dir_acf = subs$ACF,
        dir_qq = subs$QQ,
        etiqueta = etiq_ar1
      )
    } else {
      warning(sprintf("Ningun modelo AR1 convergio para %s; se conserva el espacial.", etiqueta))
    }
  } else {
    cat("  [B4] Sin autocorrelacion -> modelo espacial suficiente.\n")
  }

  # ── B6: Comparacion y seleccion ──────────────────────────────────────────
  cat("  [B6] Comparando modelos...\n")
  metricas_esp <- rbindlist(lapply(res_espacial, `[[`, "metricas"))
  metricas_ar1 <- if (length(res_ar1) > 0) {
    rbindlist(lapply(res_ar1, `[[`, "metricas"))
  } else {
    data.table()
  }

  tabla_comp <- rbindlist(list(metricas_esp, metricas_ar1), fill = TRUE)
  setorder(tabla_comp, RMSE)
  fwrite(tabla_comp, file.path(subs$comparacion, "resultados_modelos.csv"))

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
      position = position_dodge(width = 0.7), vjust = -0.4, size = 3.2
    ) +
    scale_fill_manual(
      values = c(
        "Solo espacial" = "#B2182B",
        "Espacio-temporal (AR1)" = "#2166AC"
      ),
      na.value = "grey70"
    ) +
    labs(
      title = sprintf("RMSE - %s", etiqueta),
      subtitle = sprintf("Test: %s", datos$lab_test),
      x = "Resolucion de malla", y = "RMSE (log NO2)", fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  ggsave(file.path(subs$comparacion, "grafico_rmse.png"),
    width = 10, height = 6, dpi = 300
  )

  ggplot(tabla_comp, aes(x = Malla, y = Tiempo_min, fill = Modelo)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = paste0(Tiempo_min, " min")),
      position = position_dodge(width = 0.7), vjust = -0.4, size = 3.2
    ) +
    scale_fill_manual(
      values = c(
        "Solo espacial" = "#B2182B",
        "Espacio-temporal (AR1)" = "#2166AC"
      ),
      na.value = "grey70"
    ) +
    labs(
      title = sprintf("Tiempo de computo - %s", etiqueta),
      x = "Resolucion de malla", y = "Tiempo (min)", fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  ggsave(file.path(subs$comparacion, "grafico_tiempo.png"),
    width = 10, height = 6, dpi = 300
  )

  if (!is.null(modelo_ar1_gan) && !is.null(res_esp_gan$modelo)) {
    comparar_significancia(
      modelo_esp   = res_esp_gan$modelo,
      modelo_ar1   = modelo_ar1_gan,
      vars_finales = sel$vars_finales,
      carpeta_comp = subs$comparacion,
      etiqueta     = etiqueta
    )
  }

  # ── B6.5: Tablas de diagnostico (efectos fijos, hiper, SPDE) ─────────────
  cat("  [B6.5] Generando tablas de diagnostico del modelo...\n")
  guardar_diagnosticos_modelo(
    modelo = res_esp_gan$modelo,
    spde_obj = res_esp_gan$spde,
    carpeta_diag = subs$diagnosticos,
    etiqueta = etiqueta,
    tipo_modelo = "Solo espacial"
  )
  if (!is.null(modelo_ar1_gan)) {
    guardar_diagnosticos_modelo(
      modelo = modelo_ar1_gan,
      spde_obj = spde_ar1_gan,
      carpeta_diag = subs$diagnosticos,
      etiqueta = etiqueta,
      tipo_modelo = "Espacio-temporal (AR1)"
    )
  }

  # ── B7: Tiempo de recuperacion ────────────────────────────────────────────
  cat("  [B7] Calculando tiempo de recuperacion...\n")
  tryCatch(
    calcular_recuperacion(
      modelo_ar1  = modelo_ar1_gan,
      spde_obj    = spde_ar1_gan,
      datos       = datos,
      sel         = sel,
      stk_gan     = stk_gan,
      tag_gan     = tag_gan,
      carpeta_rec = subs$recuperacion,
      etiqueta    = etiqueta
    ),
    error = function(e) warning(sprintf(
      "No se pudo completar B7 para %s: %s", etiqueta, conditionMessage(e)
    ))
  )

  tabla_comp[, Frecuencia := toupper(freq)]
  tabla_comp[, Periodo := toupper(periodo)]
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
  mejor_por_config[, .(
    Frecuencia, Periodo, Malla, Modelo,
    RMSE, MAE, Cov95, Rango_km, Sigma2, Tiempo_min
  )],
  titulo = "Mejor Configuracion por Frecuencia y Periodo",
  subtitulo = "Modelo con menor RMSE en test por frecuencia y periodo",
  ruta_png = here("outputs", "simulacion", "tabla_resumen_global.png"),
  ancho_px = 1100
)

guardar_tabla_png(
  tabla_global[, .(
    Frecuencia, Periodo, Malla, Modelo, n_nodos,
    DIC, WAIC, RMSE, MAE, Cov95, Tiempo_min
  )],
  titulo = "Resultados Completos - Todas las Frecuencias y Periodos",
  subtitulo = "B3 decide si AR1 | Variables independientes por frecuencia",
  ruta_png = here("outputs", "simulacion", "tabla_resultados_completos.png"),
  ancho_px = 1300
)

# ==============================================================================
# TABLA CLASIFICACION MACRO vs MICRO
# ==============================================================================
# Compara la reduccion de variabilidad de cada covariable entre escalas:
#   - MACRO: variable significativa a escala MENSUAL (efecto estacional/interanual)
#   - MICRO: variable significativa solo a escala DIARIA/HORARIA (efecto a corto plazo)
#   - AMBOS: significativa en ambas escalas
# ==============================================================================

if (length(resumen_parciales) > 0) {
  tabla_parciales <- rbindlist(resumen_parciales, fill = TRUE)

  # Pivotear: una fila por variable, columnas = reduccion por frecuencia
  tabla_macro_micro <- dcast(
    tabla_parciales,
    Variable ~ Frecuencia,
    value.var = "Reduccion_pct",
    fun.aggregate = mean
  )

  # Clasificar: umbral = 2% de reduccion como "significativa"
  UMBRAL_RELEVANCIA <- 2.0

  # Detectar columnas de frecuencia presentes
  tiene_mensual <- "MENSUAL" %in% names(tabla_macro_micro)
  tiene_diario <- "DIARIO" %in% names(tabla_macro_micro)
  tiene_horario <- "HORARIO" %in% names(tabla_macro_micro)

  tabla_macro_micro[, Clasificacion := {
    sig_mensual <- if (tiene_mensual) !is.na(MENSUAL) & MENSUAL >= UMBRAL_RELEVANCIA else FALSE
    sig_diario <- if (tiene_diario) !is.na(DIARIO) & DIARIO >= UMBRAL_RELEVANCIA else FALSE
    sig_horario <- if (tiene_horario) !is.na(HORARIO) & HORARIO >= UMBRAL_RELEVANCIA else FALSE
    sig_micro <- sig_diario | sig_horario

    fifelse(
      sig_mensual & sig_micro, "AMBOS (macro + micro)",
      fifelse(
        sig_mensual, "MACRO (estacional)",
        fifelse(
          sig_micro, "MICRO (corto plazo)",
          "NO RELEVANTE"
        )
      )
    )
  }, by = Variable]

  # Redondear columnas numéricas
  cols_num <- intersect(c("MENSUAL", "DIARIO", "HORARIO"), names(tabla_macro_micro))
  for (col in cols_num) {
    tabla_macro_micro[, (col) := round(get(col), 2)]
  }

  guardar_tabla_png(
    tabla_macro_micro,
    titulo = "Clasificacion Macro vs Micro de Covariables",
    subtitulo = sprintf(
      "MACRO = sig. a escala mensual (7 anos) | MICRO = sig. a escala diaria/horaria | Umbral = %.0f%% reduccion",
      UMBRAL_RELEVANCIA
    ),
    ruta_png = here("outputs", "simulacion", "tabla_macro_micro_covariables.png"),
    ancho_px = 1000
  )

  cat("\n--- Clasificacion de covariables ---\n")
  print(tabla_macro_micro)
  cat("\nTabla guardada: tabla_macro_micro_covariables.png\n")
}

# ==============================================================================
# BLOQUE 8 — MODELO DE RESIDUOS DIARIO -> HORARIO (solo horario)
# ==============================================================================
# Independiente del bucle principal: carga sus propios datos (diario año
# completo + horario de un mes). Mes de la etapa horaria configurable.

cat("\n[B8] Modelo de residuos diario -> horario...\n")
res_b8 <- tryCatch(
  modelo_residuos_diario_horario(mes_horario = "01"),
  error = function(e) {
    warning(sprintf("No se pudo completar B8: %s", conditionMessage(e)))
    NULL
  }
)

cat("\n", strrep("=", 80), "\n")
cat("  SIMULACION COMPLETADA\n")
cat(strrep("=", 80), "\n")
cat("  Resultados guardados en:", here("outputs", "simulacion"), "\n")
