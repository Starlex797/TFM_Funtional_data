# ==============================================================================
# SCATTER NO₂ vs COVARIABLES — POR ESTACIÓN DE MEDICIÓN · MADRID 2025
# Relación mensual y estacional entre cada covariable y el NO₂ bruto
# Una carpeta por estación · Un PNG por covariable y escala temporal
# Outputs: outputs/covariables_vs_no2_por_estacion/{mensual|estacional}/{estacion}/
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA Y PREPARACIÓN DEL DATASET MAESTRO
# ==============================================================================

cat("Cargando dataset maestro diario 2025...\n")

dt_m <- readRDS(here("data", "processed", "Maestro", "diario",
                     "dataset_maestro_inla_2025_DIARIO.rds"))

# Eliminar misalignment: solo filas con NO₂ válido
dt_m <- dt_m[!is.na(DATO_DIARIO)]

# Mes y nombre de mes
meses_es <- c("Enero","Febrero","Marzo","Abril","Mayo","Junio",
               "Julio","Agosto","Septiembre","Octubre","Noviembre","Diciembre")
dt_m[, mes_num := as.integer(format(FECHA, "%m"))]
dt_m[, mes_nom := factor(meses_es[mes_num], levels = meses_es)]

# Estación del año
dt_m[, estacion_anio := factor(fcase(
  mes_num %in% c(12L, 1L, 2L),  "Invierno",
  mes_num %in% c(3L, 4L, 5L),   "Primavera",
  mes_num %in% c(6L, 7L, 8L),   "Verano",
  mes_num %in% c(9L, 10L, 11L), "Oto\u00f1o"
), levels = c("Invierno", "Primavera", "Verano", "Oto\u00f1o"))]

cat(sprintf("  Filas v\u00e1lidas : %d\n", nrow(dt_m)))
cat(sprintf("  Estaciones   : %d\n", uniqueN(dt_m$ESTACION)))

# ==============================================================================
# 2. DEFINICIÓN DE COVARIABLES (nombre exacto resuelto dinámicamente)
# ==============================================================================

col_presion_raw   <- grep("^Presion",      names(dt_m), value = TRUE)
col_presion_raw   <- col_presion_raw[grepl("_raw", col_presion_raw)]
col_radiacion_raw <- grep("^Radiaci",      names(dt_m), value = TRUE)
col_radiacion_raw <- col_radiacion_raw[grepl("_raw", col_radiacion_raw)]
col_viento_raw    <- grep("^Velocidad",    names(dt_m), value = TRUE)
col_viento_raw    <- col_viento_raw[grepl("_raw", col_viento_raw)]

vars_cov <- list(
  list(col = "Temperatura_raw",       label = "Temperatura (\u00b0C)",              file = "temperatura"),
  list(col = "Humedad_Relativa_raw",  label = "Humedad relativa (%)",               file = "humedad"),
  list(col = "Precipitaciones_raw",   label = "Precipitaciones (mm)",               file = "precipitacion"),
  list(col = col_presion_raw,         label = "Presi\u00f3n barom\u00e9trica (mbar)", file = "presion"),
  list(col = col_radiacion_raw,       label = "Radiaci\u00f3n solar (W/m\u00b2)",    file = "radiacion"),
  list(col = col_viento_raw,          label = "Velocidad viento (m/s)",              file = "viento"),
  list(col = "intensidad_raw",        label = "Intensidad tr\u00e1fico (veh/h)",    file = "intensidad"),
  list(col = "carga_raw",             label = "Carga tr\u00e1fico",                 file = "carga")
)

# ==============================================================================
# 3. FUNCIÓN AUXILIAR: nombre de carpeta limpio
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  tolower(gsub("^_|_$", "", x))
}

# ==============================================================================
# 4. DIRECTORIOS BASE DE SALIDA
# ==============================================================================

dir_base       <- here("outputs", "covariables_vs_no2_por_estacion")
dir_mensual    <- file.path(dir_base, "mensual")
dir_estacional <- file.path(dir_base, "estacional")
dir.create(dir_mensual,    recursive = TRUE, showWarnings = FALSE)
dir.create(dir_estacional, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nDirectorio de salida: %s\n\n", dir_base))

# ==============================================================================
# 5. PARÁMETROS ESTÉTICOS COMPARTIDOS
# ==============================================================================

tema_scatter <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(color = "gray40", size = 9),
      plot.caption      = element_text(color = "gray55", size = 7.5),
      strip.text        = element_text(face = "bold", size = 9.5,
                                       margin = margin(t = 4, b = 4)),
      strip.background  = element_rect(fill = "gray96", color = "gray80"),
      panel.grid.minor  = element_blank(),
      panel.spacing     = unit(0.8, "lines"),
      axis.title        = element_text(size = 9)
    )
}

# ==============================================================================
# 6. FUNCIÓN GENERADORA DE SCATTER
# ==============================================================================

hacer_scatter <- function(dt_est, col_cov, label_cov, titulo, facet_var,
                          nrow_f, ncol_f) {

  # Seleccionar y limpiar
  dt_p <- dt_est[!is.na(get(col_cov)),
                 .(FECHA,
                   NO2 = DATO_DIARIO,
                   COV = get(col_cov),
                   FACET = get(facet_var))]

  if (nrow(dt_p) < 5L) return(NULL)

  # Correlación global (para subtítulo)
  r_global <- round(cor(dt_p$COV, dt_p$NO2, use = "complete.obs"), 2)
  signo    <- ifelse(r_global >= 0, "positiva", "negativa")

  ggplot(dt_p, aes(x = COV, y = NO2)) +

    # Nube de puntos
    geom_point(size = 0.9, alpha = 0.45, color = "#4a90d9") +

    # Línea de tendencia lineal con banda de confianza
    geom_smooth(method = "lm", formula = y ~ x,
                se = TRUE, color = "#e74c3c",
                fill = "#e74c3c", alpha = 0.12, linewidth = 0.9) +

    # Facetas por mes o estación
    facet_wrap(~ FACET, nrow = nrow_f, ncol = ncol_f, scales = "free") +

    scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +

    labs(
      title    = titulo,
      subtitle = sprintf(
        "Relaci\u00f3n global: r = %.2f (%s)  \u00b7  L\u00ednea roja: regresi\u00f3n lineal  \u00b7  Banda: IC 95%%",
        r_global, signo
      ),
      x       = label_cov,
      y       = "NO\u2082 (\u00b5g/m\u00b3)",
      caption = "Cada punto: observaci\u00f3n diaria  \u00b7  Pendiente positiva (+): NO\u2082 sube con la covariable  \u00b7  Negativa (\u2212): baja"
    ) +
    tema_scatter()
}

# ==============================================================================
# 7. BUCLE PRINCIPAL — ESTACIONES × COVARIABLES × ESCALA
# ==============================================================================

estaciones_unicas <- sort(unique(dt_m$ESTACION))
n_est  <- length(estaciones_unicas)
n_cov  <- length(vars_cov)
total  <- n_est * n_cov * 2L

cat(sprintf("Generando %d gr\u00e1ficas (%d estaciones \u00d7 %d covariables \u00d7 2 escalas)...\n\n",
            total, n_est, n_cov))

contador <- 0L

for (est in estaciones_unicas) {

  est_limpio <- limpiar_nombre(est)

  # Crear subcarpetas para esta estación
  dir_est_m <- file.path(dir_mensual,    est_limpio)
  dir_est_e <- file.path(dir_estacional, est_limpio)
  dir.create(dir_est_m, showWarnings = FALSE)
  dir.create(dir_est_e, showWarnings = FALSE)

  dt_est <- dt_m[ESTACION == est]

  for (vi in vars_cov) {

    col_cov <- vi$col

    if (!col_cov %in% names(dt_m)) next

    tit_base <- sprintf("NO\u2082 vs %s\n%s \u2014 Madrid 2025", vi$label, est)

    # ── Mensual (12 paneles) ──────────────────────────────────────────────────
    p_mes <- hacer_scatter(
      dt_est   = dt_est,
      col_cov  = col_cov,
      label_cov = vi$label,
      titulo   = tit_base,
      facet_var = "mes_nom",
      nrow_f   = 3L, ncol_f = 4L
    )
    if (!is.null(p_mes)) {
      archivo_m <- file.path(dir_est_m,
                             sprintf("no2_vs_%s_mensual.png", vi$file))
      ggsave(archivo_m, plot = p_mes,
             width = 13, height = 9, dpi = 180, bg = "white")
    }

    # ── Estacional (4 paneles) ────────────────────────────────────────────────
    p_est <- hacer_scatter(
      dt_est   = dt_est,
      col_cov  = col_cov,
      label_cov = vi$label,
      titulo   = tit_base,
      facet_var = "estacion_anio",
      nrow_f   = 2L, ncol_f = 2L
    )
    if (!is.null(p_est)) {
      archivo_e <- file.path(dir_est_e,
                             sprintf("no2_vs_%s_estacional.png", vi$file))
      ggsave(archivo_e, plot = p_est,
             width = 10, height = 8, dpi = 180, bg = "white")
    }

    contador <- contador + 2L
  }

  cat(sprintf("  \u2713 %-35s  [%d gr\u00e1ficas]\n", est, n_cov * 2L))
}

# ==============================================================================
# 8. RESUMEN
# ==============================================================================

cat(sprintf("\n\u2713 %d gr\u00e1ficas guardadas en:\n  %s\n", contador, dir_base))
cat(sprintf("  Mensual    : %s\n", dir_mensual))
cat(sprintf("  Estacional : %s\n", dir_estacional))
