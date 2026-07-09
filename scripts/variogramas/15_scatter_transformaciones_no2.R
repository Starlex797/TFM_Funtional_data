# ==============================================================================
# SCATTER: COMPARACIÓN NO₂ vs log(NO₂) · EFECTO DE X vs X²
# Verifica visualmente si la transformación log + término cuadrático
# mejoran la linealidad · Temperatura y Velocidad del viento
# Outputs: outputs/scatter_no2_covariables/transformaciones/
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA DE DATOS
# ==============================================================================

cat("Cargando datos...\n")

dt_d <- readRDS(here("data", "processed", "Maestro", "diario",
                     "dataset_maestro_inla_2025_DIARIO.rds"))
dt_d <- dt_d[!is.na(DATO_DIARIO)]

dt_h <- readRDS(here("data", "processed", "Maestro", "horario",
                     "dataset_maestro_inla_2025_HORARIO.rds"))
dt_h <- dt_h[!is.na(DATO)]

# Estación del año
agregar_ea <- function(dt) {
  dt[, mes := as.integer(format(FECHA, "%m"))]
  dt[, estacion_anio := fcase(
    mes %in% c(12L,1L,2L),  "Invierno",
    mes %in% c(3L,4L,5L),   "Primavera",
    mes %in% c(6L,7L,8L),   "Verano",
    mes %in% c(9L,10L,11L), "Oto\u00f1o"
  )]
  dt[, estacion_anio := factor(estacion_anio,
     levels = c("Invierno","Primavera","Verano","Oto\u00f1o"))]
  dt
}

dt_d <- agregar_ea(dt_d)
dt_h <- agregar_ea(dt_h)

cat(sprintf("  Diario  : %d filas\n  Horario : %d filas\n\n",
            nrow(dt_d), nrow(dt_h)))

# ==============================================================================
# 2. RESOLVER COLUMNA VIENTO
# ==============================================================================

col_viento_d <- grep("^Velocidad.*_raw$", names(dt_d), value = TRUE)
col_viento_h <- grep("^Velocidad.*_raw$", names(dt_h), value = TRUE)

# ==============================================================================
# 3. DEFINICIÓN DE VARIABLES A COMPARAR
# ==============================================================================

vars_comparar <- list(
  list(
    col_d     = "Temperatura_raw",
    col_h     = "Temperatura_raw",
    label_x   = "Temperatura (\u00b0C)",
    label_x2  = "Temperatura\u00b2 (\u00b0C\u00b2)",
    file      = "temperatura"
  ),
  list(
    col_d     = col_viento_d,
    col_h     = col_viento_h,
    label_x   = "Velocidad viento (m/s)",
    label_x2  = "Velocidad viento\u00b2 (m/s)\u00b2",
    file      = "viento"
  )
)

# ==============================================================================
# 4. PALETA Y TEMA
# ==============================================================================

paleta_ea <- c("Invierno"="#2980b9", "Primavera"="#27ae60",
               "Verano"="#e67e22",   "Oto\u00f1o"="#8e44ad")

tema_comp <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(color = "gray40", size = 9.5),
      plot.caption     = element_text(color = "gray55", size = 8),
      strip.text       = element_text(face = "bold", size = 10,
                                      margin = margin(t = 4, b = 4)),
      strip.background = element_rect(fill = "gray96", color = "gray80"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1.2, "lines"),
      axis.title       = element_text(size = 10)
    )
}

# ==============================================================================
# 5. DIRECTORIO DE SALIDA
# ==============================================================================

dir_salida <- here("outputs", "scatter_no2_covariables", "transformaciones")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Directorio: %s\n\n", dir_salida))

# ==============================================================================
# 6. FUNCIÓN: GENERAR LA COMPARACIÓN DE 3 PANELES
# ==============================================================================

scatter_comparacion <- function(dt, col_x, col_no2, col_log_no2,
                                 label_x, label_x2,
                                 titulo, escala,
                                 pt_size, pt_alpha) {

  dt_p <- dt[!is.na(get(col_x)) & !is.na(get(col_no2)) & !is.na(get(col_log_no2)),
             .(estacion_anio,
               X      = get(col_x),
               X2     = get(col_x)^2,
               NO2    = get(col_no2),
               logNO2 = get(col_log_no2))]

  n_pts <- nrow(dt_p)

  # Panel 1: NO₂ vs X
  p1 <- data.table(estacion_anio = dt_p$estacion_anio,
                   x_val = dt_p$X, y_val = dt_p$NO2,
                   panel = sprintf("NO\u2082 vs %s", label_x))

  # Panel 2: log(NO₂) vs X
  p2 <- data.table(estacion_anio = dt_p$estacion_anio,
                   x_val = dt_p$X, y_val = dt_p$logNO2,
                   panel = sprintf("log(NO\u2082) vs %s", label_x))

  # Panel 3: log(NO₂) vs X²
  p3 <- data.table(estacion_anio = dt_p$estacion_anio,
                   x_val = dt_p$X2, y_val = dt_p$logNO2,
                   panel = sprintf("log(NO\u2082) vs %s", label_x2))

  dt_all <- rbindlist(list(p1, p2, p3))
  dt_all[, panel := factor(panel, levels = c(
    sprintf("NO\u2082 vs %s", label_x),
    sprintf("log(NO\u2082) vs %s", label_x),
    sprintf("log(NO\u2082) vs %s", label_x2)
  ))]

  ggplot(dt_all, aes(x = x_val, y = y_val, color = estacion_anio)) +

    geom_point(size = pt_size, alpha = pt_alpha, na.rm = TRUE) +

    facet_wrap(~ panel, nrow = 1, scales = "free") +

    scale_color_manual(values = paleta_ea,
                       name   = "Estaci\u00f3n del a\u00f1o") +

    labs(
      title    = titulo,
      subtitle = sprintf(
        "Escala %s  \u00b7  %s observaciones  \u00b7  Todas las estaciones de medici\u00f3n",
        escala, format(n_pts, big.mark = ".")
      ),
      x        = NULL,
      y        = NULL,
      caption  = paste0(
        "Izquierda: relaci\u00f3n original (NO\u2082 bruto)  \u00b7  ",
        "Centro: tras log-transformar la respuesta  \u00b7  ",
        "Derecha: log(NO\u2082) vs t\u00e9rmino cuadr\u00e1tico"
      )
    ) +
    tema_comp() +
    guides(color = guide_legend(
      nrow = 1,
      override.aes = list(size = 3, alpha = 0.9)
    ))
}

# ==============================================================================
# 7. BUCLE PRINCIPAL
# ==============================================================================

for (vi in vars_comparar) {

  nombre <- vi$file

  # ── DIARIO ──────────────────────────────────────────────────────────────────
  if (length(vi$col_d) > 0 && vi$col_d %in% names(dt_d)) {

    p <- scatter_comparacion(
      dt           = dt_d,
      col_x        = vi$col_d,
      col_no2      = "DATO_DIARIO",
      col_log_no2  = "LOG_NO2_DIARIO",
      label_x      = vi$label_x,
      label_x2     = vi$label_x2,
      titulo       = sprintf(
        "Efecto de la transformaci\u00f3n log y cuadr\u00e1tica \u2014 %s",
        vi$label_x),
      escala       = "diaria",
      pt_size      = 1.1,
      pt_alpha     = 0.50
    )

    archivo <- file.path(dir_salida,
                          sprintf("comparacion_diario_%s_2025.png", nombre))
    ggsave(archivo, plot = p, width = 16, height = 6, dpi = 200, bg = "white")
    cat(sprintf("  \u2713 [diario]  %s\n", basename(archivo)))
  }

  # ── HORARIO ─────────────────────────────────────────────────────────────────
  if (length(vi$col_h) > 0 && vi$col_h %in% names(dt_h)) {

    p <- scatter_comparacion(
      dt           = dt_h,
      col_x        = vi$col_h,
      col_no2      = "DATO",
      col_log_no2  = "LOG_NO2_HORARIO",
      label_x      = vi$label_x,
      label_x2     = vi$label_x2,
      titulo       = sprintf(
        "Efecto de la transformaci\u00f3n log y cuadr\u00e1tica \u2014 %s",
        vi$label_x),
      escala       = "horaria",
      pt_size      = 0.15,
      pt_alpha     = 0.12
    )

    archivo <- file.path(dir_salida,
                          sprintf("comparacion_horario_%s_2025.png", nombre))
    ggsave(archivo, plot = p, width = 16, height = 6, dpi = 200, bg = "white")
    cat(sprintf("  \u2713 [horario] %s\n\n", basename(archivo)))
  }
}

cat(sprintf("\u2713 Archivos guardados en:\n  %s\n", dir_salida))
