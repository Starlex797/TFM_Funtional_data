# ==============================================================================
# SCATTER PLOTS: NO₂ vs COVARIABLES — DIARIO Y HORARIO — MADRID 2025
# Facetado por estación del año · Sin línea de tendencia
# Color por estación de medición
# Outputs: outputs/scatter_no2_covariables/diario/ y /horario/
# ==============================================================================
# QUÉ HACE:
#
# - Carga los datasets maestros de 2025 en dos resoluciones:
#     · Diaria.
#     · Horaria.
#
# - Elimina las observaciones sin una concentración válida de NO2.
# - Añade la estación del año a ambos datasets.
# - Resuelve dinámicamente los nombres de columnas con tildes o caracteres
#   especiales.
#
# - Analiza las siguientes covariables:
#     · Temperatura.
#     · Humedad relativa.
#     · Precipitaciones.
#     · Presión barométrica.
#     · Radiación solar.
#     · Velocidad del viento.
#     · Intensidad del tráfico.
#     · Carga de tráfico.
#
# - Genera para cada covariable:
#     · Un scatter plot con datos diarios.
#     · Un scatter plot con datos horarios.
#
# - Los puntos se colorean según la estación de contaminación.
# - Los gráficos se dividen en cuatro paneles según la estación del año.
# - Informa del número de observaciones válidas utilizadas.
# - No añade rectas de tendencia ni calcula correlaciones.
#
# FINALIDAD PARA EL TFM:
#
# Este script ofrece una visión conjunta de la relación entre el NO2 y las
# covariables, diferenciando:
#
# - Las estaciones de contaminación.
# - Las estaciones del año.
# - La resolución diaria y horaria.
#
# Permite detectar visualmente:
#
# - Relaciones lineales o no lineales.
# - Cambios estacionales.
# - Diferencias entre estaciones.
# - Valores extremos.
# - Acumulaciones o patrones particulares de los datos.
# - Diferencias entre la relación diaria y la relación horaria.
#
# Es especialmente útil para decidir si una covariable debe modelizarse de forma
# lineal, no lineal o mediante interacciones estacionales.
#
# SALIDAS:
#
# outputs/scatter_no2_covariables/
#     · diario/
#     · horario/
#
# Se genera un gráfico diario y otro horario para cada covariable.


library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. FUNCIÓN AUXILIAR
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  tolower(gsub("^_|_$", "", x))
}

estacion_anio_var <- function(dt) {
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
dt_d <- estacion_anio_var(dt_d)
cat(sprintf(
  "  Diario  : %d filas · %d estaciones\n",
  nrow(dt_d), uniqueN(dt_d$ESTACION)
))

cat("Cargando maestro HORARIO 2025...\n")
dt_h <- readRDS(here(
  "data", "processed", "Maestro", "horario",
  "dataset_maestro_inla_2025_HORARIO.rds"
))
dt_h <- dt_h[!is.na(DATO)]
dt_h <- estacion_anio_var(dt_h)
cat(sprintf(
  "  Horario : %d filas · %d estaciones\n",
  nrow(dt_h), uniqueN(dt_h$ESTACION)
))

# ==============================================================================
# 3. RESOLUCIÓN DE NOMBRES DE COLUMNAS CON CARACTERES ESPECIALES
# ==============================================================================

col_presion_d <- grep("^Presion", names(dt_d), value = TRUE)
col_presion_d <- col_presion_d[grepl("_raw", col_presion_d)]

col_presion_h <- grep("^Presion", names(dt_h), value = TRUE)
col_presion_h <- col_presion_h[grepl("_raw", col_presion_h)]

col_radiacion_d <- grep("^Radiaci", names(dt_d), value = TRUE)
col_radiacion_d <- col_radiacion_d[grepl("_raw", col_radiacion_d)]

col_radiacion_h <- grep("^Radiaci", names(dt_h), value = TRUE)
col_radiacion_h <- col_radiacion_h[grepl("_raw", col_radiacion_h)]

col_viento_d <- grep("^Velocidad", names(dt_d), value = TRUE)
col_viento_d <- col_viento_d[grepl("_raw", col_viento_d)]

col_viento_h <- grep("^Velocidad", names(dt_h), value = TRUE)
col_viento_h <- col_viento_h[grepl("_raw", col_viento_h)]

# ==============================================================================
# 4. DEFINICIÓN DE COVARIABLES (diario y horario comparten las mismas)
# ==============================================================================

vars_cov <- list(
  list(
    col_d = "Temperatura_raw", col_h = "Temperatura_raw",
    label = "Temperatura (\u00b0C)",
    file = "temperatura"
  ),
  list(
    col_d = "Humedad_Relativa_raw", col_h = "Humedad_Relativa_raw",
    label = "Humedad relativa (%)",
    file = "humedad"
  ),
  list(
    col_d = "Precipitaciones_raw", col_h = "Precipitaciones_raw",
    label = "Precipitaciones (mm)",
    file = "precipitacion"
  ),
  list(
    col_d = col_presion_d, col_h = col_presion_h,
    label = "Presi\u00f3n barom\u00e9trica (mbar)",
    file = "presion"
  ),
  list(
    col_d = col_radiacion_d, col_h = col_radiacion_h,
    label = "Radiaci\u00f3n solar (W/m\u00b2)",
    file = "radiacion"
  ),
  list(
    col_d = col_viento_d, col_h = col_viento_h,
    label = "Velocidad del viento (m/s)",
    file = "viento"
  ),
  list(
    col_d = "intensidad_raw", col_h = "intensidad_raw",
    label = "Intensidad tr\u00e1fico (veh/h)",
    file = "intensidad"
  ),
  list(
    col_d = "carga_raw", col_h = "carga_raw",
    label = "Carga tr\u00e1fico",
    file = "carga"
  )
)

# ==============================================================================
# 5. PALETAS Y TEMA
# ==============================================================================

estaciones_ord <- sort(union(unique(dt_d$ESTACION), unique(dt_h$ESTACION)))
n_est <- length(estaciones_ord)
paleta_est <- setNames(scales::hue_pal()(n_est), estaciones_ord)

tema_scatter <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "gray40", size = 9.5),
      plot.caption = element_text(color = "gray55", size = 8),
      strip.text = element_text(
        face = "bold", size = 11,
        margin = margin(t = 5, b = 5)
      ),
      strip.background = element_rect(fill = "gray96", color = "gray80"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 7.5),
      legend.key.size = unit(0.55, "cm"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.title = element_text(size = 10)
    )
}

# ==============================================================================
# 6. DIRECTORIOS DE SALIDA
# ==============================================================================

dir_d <- here("outputs", "scatter_no2_covariables", "diario")
dir_h <- here("outputs", "scatter_no2_covariables", "horario")
dir.create(dir_d, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_h, recursive = TRUE, showWarnings = FALSE)
cat(sprintf(
  "\nDirectorio de salida: %s\n\n",
  here("outputs", "scatter_no2_covariables")
))

# ==============================================================================
# 7. FUNCIÓN GENERADORA DE SCATTER
# ==============================================================================

hacer_scatter <- function(datos, col_x, col_y, label_x, label_y,
                          titulo, subtitulo, caption,
                          punto_size, punto_alpha) {
  dt_p <- datos[
    !is.na(get(col_x)) & !is.na(get(col_y)),
    .(ESTACION, estacion_anio,
      X = get(col_x), Y = get(col_y)
    )
  ]

  n_pts <- nrow(dt_p)

  ggplot(dt_p, aes(x = X, y = Y, color = ESTACION)) +
    geom_point(size = punto_size, alpha = punto_alpha, na.rm = TRUE) +
    facet_wrap(~estacion_anio, nrow = 2, ncol = 2) +
    scale_color_manual(values = paleta_est, name = "Estaci\u00f3n") +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    labs(
      title = titulo,
      subtitle = sprintf(
        "%s  \u00b7  %s observaciones v\u00e1lidas",
        subtitulo, format(n_pts, big.mark = ".")
      ),
      x = label_x,
      y = label_y,
      caption = caption
    ) +
    tema_scatter() +
    guides(color = guide_legend(
      ncol = 4,
      override.aes = list(size = 2.5, alpha = 0.9)
    ))
}

# ==============================================================================
# 8. BUCLE PRINCIPAL
# ==============================================================================

cat("--- Generando scatter plots DIARIOS ---\n")

for (vi in vars_cov) {
  # ── DIARIO ──────────────────────────────────────────────────────────────────
  col_x <- vi$col_d
  if (length(col_x) == 0 || !col_x %in% names(dt_d)) {
    cat(sprintf("  [SKIP diario] %s\n", vi$file))
    next
  }

  p_d <- hacer_scatter(
    datos       = dt_d,
    col_x       = col_x,
    col_y       = "DATO_DIARIO",
    label_x     = vi$label,
    label_y     = "NO\u2082 (\u00b5g/m\u00b3)",
    titulo      = sprintf("NO\u2082 vs %s \u2014 Escala diaria, Madrid 2025", vi$label),
    subtitulo   = "Todas las estaciones de medici\u00f3n  \u00b7  Facetado por estaci\u00f3n del a\u00f1o",
    caption     = "Cada punto: observaci\u00f3n diaria de una estaci\u00f3n  \u00b7  Sin l\u00ednea de tendencia",
    punto_size  = 1.2,
    punto_alpha = 0.55
  )

  archivo_d <- file.path(
    dir_d,
    sprintf("scatter_diario_no2_vs_%s_2025.png", vi$file)
  )
  ggsave(archivo_d, plot = p_d, width = 13, height = 9, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 [diario]  %s\n", basename(archivo_d)))

  # ── HORARIO ─────────────────────────────────────────────────────────────────
  col_x_h <- vi$col_h
  if (length(col_x_h) == 0 || !col_x_h %in% names(dt_h)) {
    cat(sprintf("  [SKIP horario] %s\n", vi$file))
    next
  }

  p_h <- hacer_scatter(
    datos       = dt_h,
    col_x       = col_x_h,
    col_y       = "DATO",
    label_x     = vi$label,
    label_y     = "NO\u2082 (\u00b5g/m\u00b3)",
    titulo      = sprintf("NO\u2082 vs %s \u2014 Escala horaria, Madrid 2025", vi$label),
    subtitulo   = "Todas las estaciones de medici\u00f3n  \u00b7  Facetado por estaci\u00f3n del a\u00f1o",
    caption     = "Cada punto: observaci\u00f3n horaria de una estaci\u00f3n  \u00b7  Sin l\u00ednea de tendencia",
    punto_size  = 0.25,
    punto_alpha = 0.18
  )

  archivo_h <- file.path(
    dir_h,
    sprintf("scatter_horario_no2_vs_%s_2025.png", vi$file)
  )
  ggsave(archivo_h, plot = p_h, width = 13, height = 9, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 [horario] %s\n\n", basename(archivo_h)))
}

# ==============================================================================
# 9. RESUMEN
# ==============================================================================

cat("==============================================================\n")
cat(sprintf("  Covariables procesadas : %d\n", length(vars_cov)))
cat(sprintf("  PNG diarios   (%d)  \u2192  %s\n", length(vars_cov), dir_d))
cat(sprintf("  PNG horarios  (%d)  \u2192  %s\n", length(vars_cov), dir_h))
cat("==============================================================\n")
