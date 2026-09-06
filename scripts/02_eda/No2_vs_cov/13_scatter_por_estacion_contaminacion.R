# ==============================================================================
# SCATTER NO₂ vs COVARIABLES — POR ESTACIÓN DE CONTAMINACIÓN
# Una gráfica por estación · 8 paneles (uno por covariable) · Color = estación año
# Escala diaria y horaria · Madrid 2025
# Outputs: outputs/analysis/scatter_por_estacion/diario/ y /horario/
# ==============================================================================
# ==============================================================================
# 13_scatter_por_estacion_contaminacion.R
# ==============================================================================
#
# QUÉ HACE:
#
# - Carga los datasets maestros diario y horario de 2025.
# - Elimina las observaciones sin valores válidos de NO2.
# - Clasifica las observaciones según la estación del año.
# - Localiza dinámicamente las columnas de las ocho covariables.
#
# - Genera una figura independiente para cada estación de contaminación.
# - Cada figura contiene ocho paneles, uno por covariable:
#     · Temperatura.
#     · Humedad relativa.
#     · Precipitaciones.
#     · Presión barométrica.
#     · Radiación solar.
#     · Velocidad del viento.
#     · Intensidad del tráfico.
#     · Carga de tráfico.
#
# - El eje vertical contiene las concentraciones de NO2.
# - El eje horizontal contiene el valor de cada covariable.
# - Los puntos se colorean según la estación del año.
# - Las escalas del eje horizontal son independientes para cada covariable.
# - No se añaden líneas de tendencia.
#
# - El proceso se realiza por separado para:
#     · Datos diarios.
#     · Datos horarios.
#
# FINALIDAD PARA EL TFM:
#
# Este script resume en una única figura el comportamiento de todas las
# covariables para una estación concreta.
#
# Permite realizar una comparación rápida entre predictores y detectar:
#
# - Qué covariables parecen estar más relacionadas con el NO2.
# - Si la relación cambia entre estaciones del año.
# - Si existen agrupaciones estacionales.
# - Si aparecen valores extremos o distribuciones anómalas.
# - Si la relación es diferente en escala diaria y horaria.
# - Si una estación tiene un comportamiento diferente al resto.
#
# Es una herramienta útil para el diagnóstico exploratorio por estación y para
# decidir si los modelos necesitan efectos espaciales o coeficientes diferentes
# según la localización.
#
# SALIDAS:
#
# outputs/analysis/scatter_por_estacion/
#     · diario/
#     · horario/
#
# Se genera una figura diaria y otra horaria para cada estación de contaminación.


library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. FUNCIONES AUXILIARES
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  tolower(gsub("^_|_$", "", x))
}

agregar_estacion_anio <- function(dt) {
  dt[, mes := as.integer(format(FECHA, "%m"))]
  dt[, estacion_anio := fcase(
    mes %in% c(12L, 1L, 2L),  "Invierno",
    mes %in% c(3L, 4L, 5L),   "Primavera",
    mes %in% c(6L, 7L, 8L),   "Verano",
    mes %in% c(9L, 10L, 11L), "Oto\u00f1o"
  )]
  dt[, estacion_anio := factor(
    estacion_anio,
    levels = c("Invierno", "Primavera", "Verano", "Oto\u00f1o")
  )]
  dt
}

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

cat("Cargando maestro DIARIO 2025...\n")
dt_d <- readRDS(here(
  "data", "processed", "Maestro", "diario",
  "dataset_maestro_inla_2025_DIARIO.rds"
))
dt_d <- dt_d[!is.na(DATO_DIARIO)]
dt_d <- agregar_estacion_anio(dt_d)

cat("Cargando maestro HORARIO 2025...\n")
dt_h <- readRDS(here(
  "data", "processed", "Maestro", "horario",
  "dataset_maestro_inla_2025_HORARIO.rds"
))
dt_h <- dt_h[!is.na(DATO)]
dt_h <- agregar_estacion_anio(dt_h)

cat(sprintf(
  "  Diario  : %d filas · %d estaciones\n",
  nrow(dt_d), uniqueN(dt_d$ESTACION)
))
cat(sprintf(
  "  Horario : %d filas · %d estaciones\n\n",
  nrow(dt_h), uniqueN(dt_h$ESTACION)
))

# ==============================================================================
# 3. COLUMNAS DE COVARIABLES (resolución dinámica para nombres con tildes)
# ==============================================================================

resolver_cols_raw <- function(dt) {
  todas <- names(dt)[grepl("_raw$", names(dt))]

  list(
    Temperatura      = intersect("Temperatura_raw", todas),
    Humedad          = intersect("Humedad_Relativa_raw", todas),
    Precipitaciones  = intersect("Precipitaciones_raw", todas),
    Presion          = grep("^Presion", todas, value = TRUE),
    Radiacion        = grep("^Radiaci", todas, value = TRUE),
    Viento           = grep("^Velocidad", todas, value = TRUE),
    Intensidad       = intersect("intensidad_raw", todas),
    Carga            = intersect("carga_raw", todas)
  )
}

cols_d <- resolver_cols_raw(dt_d)
cols_h <- resolver_cols_raw(dt_h)

# Etiquetas legibles para los paneles
etiquetas_cov <- c(
  Temperatura     = "Temperatura (\u00b0C)",
  Humedad         = "Humedad relativa (%)",
  Precipitaciones = "Precipitaciones (mm)",
  Presion         = "Presi\u00f3n barom\u00e9trica (mbar)",
  Radiacion       = "Radiaci\u00f3n solar (W/m\u00b2)",
  Viento          = "Velocidad viento (m/s)",
  Intensidad      = "Intensidad tr\u00e1fico (veh/h)",
  Carga           = "Carga tr\u00e1fico"
)

# ==============================================================================
# 4. PALETA Y TEMA
# ==============================================================================

paleta_estaciones <- c(
  "Invierno" = "#2980b9",
  "Primavera" = "#27ae60",
  "Verano" = "#e67e22",
  "Oto\u00f1o" = "#8e44ad"
)

tema_scatter_est <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "gray40", size = 9),
      plot.caption = element_text(color = "gray55", size = 7.5),
      strip.text = element_text(
        face = "bold", size = 9,
        margin = margin(t = 4, b = 4)
      ),
      strip.background = element_rect(fill = "gray96", color = "gray80"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.6, "cm"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.9, "lines"),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5)
    )
}

# ==============================================================================
# 5. FUNCIÓN PRINCIPAL: SCATTER MULTI-PANEL PARA UNA ESTACIÓN
# ==============================================================================

scatter_estacion <- function(dt, col_no2, cols_cov, etiquetas,
                             nombre_est, escala, punto_size, punto_alpha) {
  # Seleccionar covariables disponibles
  cols_validas <- Filter(function(c) length(c) > 0 && c %in% names(dt), cols_cov)

  if (length(cols_validas) == 0) {
    return(invisible(NULL))
  }

  # Subconjunto de la estación
  dt_est <- dt[ESTACION == nombre_est,
    c("estacion_anio", col_no2, unlist(cols_validas)),
    with = FALSE
  ]

  # Renombrar NO2 a una columna uniforme
  setnames(dt_est, col_no2, "NO2")

  # Pivotar a largo: una fila por (observación, covariable)
  dt_long <- melt(
    dt_est,
    id.vars       = c("NO2", "estacion_anio"),
    measure.vars  = unlist(cols_validas),
    variable.name = "cov_col",
    value.name    = "cov_valor"
  )

  # Asignar etiqueta legible a cada covariable
  # Mapeo: nombre_clave → col_real → etiqueta
  mapa_etiq <- data.table(
    cov_col  = unlist(cols_validas),
    etiqueta = etiquetas[names(cols_validas)]
  )
  dt_long <- merge(dt_long, mapa_etiq, by = "cov_col", all.x = TRUE)
  dt_long[, etiqueta := factor(etiqueta, levels = etiquetas)]

  n_pts <- nrow(dt_long[!is.na(cov_valor) & !is.na(NO2)])

  ggplot(
    dt_long[!is.na(cov_valor) & !is.na(NO2)],
    aes(x = cov_valor, y = NO2, color = estacion_anio)
  ) +
    geom_point(size = punto_size, alpha = punto_alpha, na.rm = TRUE) +
    facet_wrap(~etiqueta, nrow = 4, ncol = 2, scales = "free_x") +
    scale_color_manual(
      values = paleta_estaciones,
      name = "Estaci\u00f3n del a\u00f1o"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    labs(
      title = sprintf("NO\u2082 vs Covariables \u2014 %s", nombre_est),
      subtitle = sprintf(
        "Madrid 2025  \u00b7  Escala %s  \u00b7  %s observaciones v\u00e1lidas",
        escala, format(n_pts, big.mark = ".")
      ),
      x = "Valor de la covariable",
      y = "NO\u2082 (\u00b5g/m\u00b3)",
      caption = "Cada punto: una observaci\u00f3n  \u00b7  Sin l\u00ednea de tendencia  \u00b7  Color: estaci\u00f3n del a\u00f1o"
    ) +
    tema_scatter_est() +
    guides(color = guide_legend(
      nrow = 1,
      override.aes = list(size = 3, alpha = 0.9)
    ))
}

# ==============================================================================
# 6. DIRECTORIOS DE SALIDA
# ==============================================================================

dir_diario <- here("outputs", "figures", "EDA", "NO2", "Scatter", "diario")
dir_horario <- here("outputs", "figures", "EDA", "NO2", "Scatter", "horario")
dir.create(dir_diario, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_horario, recursive = TRUE, showWarnings = FALSE)
cat(sprintf(
  "Carpeta de salida: %s\n\n",
  here("outputs", "figures", "EDA", "NO2", "Scatter")
))

# ==============================================================================
# 7. BUCLE: UNA ESTACIÓN A LA VEZ — DIARIO Y HORARIO
# ==============================================================================

estaciones_d <- sort(unique(dt_d$ESTACION))
estaciones_h <- sort(unique(dt_h$ESTACION))
estaciones_todas <- union(estaciones_d, estaciones_h)

cat(sprintf("Estaciones a procesar: %d\n\n", length(estaciones_todas)))

for (est in estaciones_todas) {
  nombre_arch <- limpiar_nombre(est)

  # ── DIARIO ──────────────────────────────────────────────────────────────────
  if (est %in% estaciones_d) {
    p <- scatter_estacion(
      dt          = dt_d,
      col_no2     = "DATO_DIARIO",
      cols_cov    = cols_d,
      etiquetas   = etiquetas_cov,
      nombre_est  = est,
      escala      = "diaria",
      punto_size  = 1.3,
      punto_alpha = 0.60
    )

    if (!is.null(p)) {
      archivo <- file.path(
        dir_diario,
        sprintf("scatter_diario_%s_2025.png", nombre_arch)
      )
      ggsave(archivo, plot = p, width = 13, height = 14, dpi = 180, bg = "white")
      cat(sprintf("  \u2713 [diario]  %s\n", basename(archivo)))
    }
  }

  # ── HORARIO ─────────────────────────────────────────────────────────────────
  if (est %in% estaciones_h) {
    p <- scatter_estacion(
      dt          = dt_h,
      col_no2     = "DATO",
      cols_cov    = cols_h,
      etiquetas   = etiquetas_cov,
      nombre_est  = est,
      escala      = "horaria",
      punto_size  = 0.2,
      punto_alpha = 0.15
    )

    if (!is.null(p)) {
      archivo <- file.path(
        dir_horario,
        sprintf("scatter_horario_%s_2025.png", nombre_arch)
      )
      ggsave(archivo, plot = p, width = 13, height = 14, dpi = 180, bg = "white")
      cat(sprintf("  \u2713 [horario] %s\n", basename(archivo)))
    }
  }

  cat("\n")
}

# ==============================================================================
# 8. RESUMEN
# ==============================================================================

n_d <- length(list.files(dir_diario, pattern = "\\.png$"))
n_h <- length(list.files(dir_horario, pattern = "\\.png$"))

cat("==============================================================\n")
cat(sprintf("  PNG diarios  : %d  \u2192  %s\n", n_d, dir_diario))
cat(sprintf("  PNG horarios : %d  \u2192  %s\n", n_h, dir_horario))
cat("==============================================================\n")
