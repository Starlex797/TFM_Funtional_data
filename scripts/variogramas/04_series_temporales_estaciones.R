# ==============================================================================
# SERIES TEMPORALES — NO₂ DIARIO Y HORARIO POR ESTACIÓN
# Madrid 2025 · Delimitadores mensuales · Un PNG por estación y escala
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA DE DATOS
# ==============================================================================

dt_diario <- readRDS(here("data", "processed", "contaminacion", "diario",
                           "aire_madrid_2025_No2_trans_diarios.rds"))

dt_horario <- readRDS(here("data", "processed", "contaminacion", "horario",
                            "aire_madrid_2025_No2_horarios.rds"))

# Construir DATETIME continuo para la serie horaria
# H01 → 00:00, H02 → 01:00, ..., H24 → 23:00
dt_horario[, Hora_num := as.integer(sub("H0?", "", as.character(HORA)))]
dt_horario[, DATETIME := as.POSIXct(FECHA) + (Hora_num - 1L) * 3600L]

# Ordenar para que la línea sea continua
setorder(dt_diario, ESTACION, FECHA)
setorder(dt_horario, ESTACION, DATETIME)

# ==============================================================================
# 2. DELIMITADORES MENSUALES
# ==============================================================================

# Primer día de cada mes (para vlines y etiquetas)
primeros_meses      <- seq(as.Date("2025-01-01"), as.Date("2025-12-01"), by = "1 month")
primeros_meses_ct   <- as.POSIXct(primeros_meses)   # versión POSIXct para gráficos horarios
etiquetas_meses     <- format(primeros_meses, "%b")  # "Ene", "Feb", ...

# Posición central de cada mes (para leyenda de fondo alternado)
midpoint_meses      <- primeros_meses + 15L
midpoint_meses_ct   <- as.POSIXct(midpoint_meses)

# Franjas alternadas de fondo (rectángulos mes par / mes impar)
# --- versión Date (diario) ---
franjas_diario <- data.frame(
  xmin  = primeros_meses,
  xmax  = c(primeros_meses[-1], as.Date("2026-01-01")),
  shade = (seq_along(primeros_meses) %% 2L == 0L)
) |> subset(shade)

# --- versión POSIXct (horario) ---
meses_limites_ct <- c(primeros_meses_ct, as.POSIXct("2026-01-01"))
franjas_horario <- data.frame(
  xmin  = meses_limites_ct[-length(meses_limites_ct)],
  xmax  = meses_limites_ct[-1],
  shade = (seq_along(primeros_meses) %% 2L == 0L)
) |> subset(shade)

# ==============================================================================
# 3. DIRECTORIOS DE SALIDA
# ==============================================================================

dir_diario  <- here("outputs", "series_temporales", "diario")
dir_horario <- here("outputs", "series_temporales", "horario")
dir.create(dir_diario,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_horario, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 4. FUNCIÓN AUXILIAR: nombre de archivo limpio
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", tolower(x))
}

# ==============================================================================
# 5. LISTA DE ESTACIONES
# ==============================================================================

estaciones_unicas <- sort(unique(dt_diario$ESTACION))
cat(sprintf("Estaciones detectadas: %d\n", length(estaciones_unicas)))
cat(paste0("  · ", estaciones_unicas, collapse = "\n"), "\n\n")

# ==============================================================================
# 6. BUCLE PRINCIPAL — DIARIO
# ==============================================================================

cat("--- Generando series DIARIAS ---\n")

for (est in estaciones_unicas) {

  dt_est <- dt_diario[ESTACION == est & !is.na(DATO_DIARIO)]

  if (nrow(dt_est) < 2L) {
    cat(sprintf("  [SKIP] %s — sin datos diarios\n", est))
    next
  }

  n_dias <- nrow(dt_est)
  pct_ok <- round(100 * n_dias / 365, 1)

  p <- ggplot(dt_est, aes(x = FECHA, y = DATO_DIARIO)) +

    # Franjas alternadas de fondo por mes
    geom_rect(
      data        = franjas_diario,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill        = "gray90", alpha = 0.5
    ) +

    # Líneas divisorias mensuales
    geom_vline(
      xintercept = as.numeric(primeros_meses),
      color      = "gray60", linewidth = 0.4, linetype = "dashed"
    ) +

    # Serie temporal
    geom_line(color = "#2c3e50", linewidth = 0.55, alpha = 0.85) +
    geom_point(size = 0.9, color = "#2980b9", alpha = 0.6) +

    # Umbral legal OMS (40 µg/m³ media anual)
    geom_hline(yintercept = 40, color = "#e74c3c",
               linewidth = 0.7, linetype = "longdash", alpha = 0.8) +
    annotate("text", x = as.Date("2025-01-05"), y = 41.5,
             label = "Límite OMS (40 µg/m³)", color = "#e74c3c",
             size = 2.8, hjust = 0, fontface = "italic") +

    # Eje X: tick en inicio de cada mes con nombre abreviado
    scale_x_date(
      breaks       = primeros_meses,
      labels       = etiquetas_meses,
      minor_breaks = NULL,
      expand       = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +

    labs(
      title    = paste0("Serie Temporal Diaria de NO\u2082 \u2014 ", est),
      subtitle = sprintf("Madrid 2025  |  %d d\u00edas v\u00e1lidos (%.1f%%)", n_dias, pct_ok),
      x        = NULL,
      y        = "NO\u2082 (\u00b5g/m\u00b3)",
      caption  = "Puntos azules: observaci\u00f3n diaria  \u00b7  L\u00ednea roja: l\u00edmite anual OMS"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "gray40", size = 10),
      plot.caption     = element_text(color = "gray55", size = 8),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )

  nombre_file <- file.path(dir_diario,
                            paste0("diario_", limpiar_nombre(est), ".png"))
  ggsave(nombre_file, plot = p, width = 13, height = 5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 7. BUCLE PRINCIPAL — HORARIO
# ==============================================================================

cat("\n--- Generando series HORARIAS ---\n")

for (est in estaciones_unicas) {

  dt_est <- dt_horario[ESTACION == est & !is.na(DATO)]

  if (nrow(dt_est) < 2L) {
    cat(sprintf("  [SKIP] %s — sin datos horarios\n", est))
    next
  }

  n_horas <- nrow(dt_est)
  pct_ok  <- round(100 * n_horas / 8760, 1)

  p <- ggplot(dt_est, aes(x = DATETIME, y = DATO)) +

    # Franjas alternadas de fondo por mes
    geom_rect(
      data        = franjas_horario,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill        = "gray90", alpha = 0.5
    ) +

    # Líneas divisorias mensuales
    geom_vline(
      xintercept  = primeros_meses_ct,
      color       = "gray60", linewidth = 0.4, linetype = "dashed"
    ) +

    # Serie horaria
    geom_line(color = "#2c3e50", linewidth = 0.3, alpha = 0.7) +

    # Umbral legal OMS
    geom_hline(yintercept = 40, color = "#e74c3c",
               linewidth = 0.6, linetype = "longdash", alpha = 0.75) +

    # Eje X: tick en inicio de cada mes con nombre abreviado
    scale_x_datetime(
      breaks       = primeros_meses_ct,
      labels       = etiquetas_meses,
      minor_breaks = NULL,
      expand       = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +

    labs(
      title    = paste0("Serie Temporal Horaria de NO\u2082 \u2014 ", est),
      subtitle = sprintf("Madrid 2025  |  %d horas v\u00e1lidas (%.1f%% del a\u00f1o)", n_horas, pct_ok),
      x        = NULL,
      y        = "NO\u2082 (\u00b5g/m\u00b3)",
      caption  = "Cada punto = medici\u00f3n horaria  \u00b7  L\u00ednea roja discontinua: l\u00edmite anual OMS (40 \u00b5g/m\u00b3)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "gray40", size = 10),
      plot.caption       = element_text(color = "gray55", size = 8),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )

  nombre_file <- file.path(dir_horario,
                            paste0("horario_", limpiar_nombre(est), ".png"))
  ggsave(nombre_file, plot = p, width = 14, height = 5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 8. PANEL RESUMEN — TODAS LAS ESTACIONES EN UN SOLO GRÁFICO (FACET)
# ==============================================================================

cat("\n--- Generando paneles resumen (todas las estaciones) ---\n")

# Panel diario
p_panel_diario <- ggplot(dt_diario[!is.na(DATO_DIARIO)],
                          aes(x = FECHA, y = DATO_DIARIO)) +

  geom_rect(
    data        = franjas_diario,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill        = "gray90", alpha = 0.5
  ) +
  geom_vline(xintercept  = as.numeric(primeros_meses),
             color = "gray65", linewidth = 0.3, linetype = "dashed") +
  geom_line(color = "#2c3e50", linewidth = 0.4, alpha = 0.8) +
  geom_hline(yintercept = 40, color = "#e74c3c",
             linewidth = 0.5, linetype = "longdash", alpha = 0.7) +
  scale_x_date(breaks = primeros_meses, labels = etiquetas_meses,
               minor_breaks = NULL, expand = expansion(mult = c(0.01, 0.01))) +
  facet_wrap(~ ESTACION, ncol = 4, scales = "free_y") +
  labs(
    title    = "Series Temporales Diarias de NO\u2082 \u2014 Todas las Estaciones",
    subtitle = "Madrid 2025  \u00b7  L\u00ednea roja: l\u00edmite anual OMS (40 \u00b5g/m\u00b3)",
    x        = NULL,
    y        = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "gray40", size = 10),
    strip.text         = element_text(face = "bold", size = 7.5),
    axis.text.x        = element_text(angle = 60, hjust = 1, size = 6),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(0.8, "lines")
  )

ggsave(
  file.path(here("outputs", "series_temporales"), "panel_diario_todas_estaciones.png"),
  plot = p_panel_diario, width = 22, height = 18, dpi = 200, bg = "white"
)
cat("  \u2713 panel_diario_todas_estaciones.png\n")

# Panel horario
p_panel_horario <- ggplot(dt_horario[!is.na(DATO)],
                           aes(x = DATETIME, y = DATO)) +

  geom_rect(
    data        = franjas_horario,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill        = "gray90", alpha = 0.5
  ) +
  geom_vline(xintercept = primeros_meses_ct,
             color = "gray65", linewidth = 0.3, linetype = "dashed") +
  geom_line(color = "#2c3e50", linewidth = 0.25, alpha = 0.65) +
  geom_hline(yintercept = 40, color = "#e74c3c",
             linewidth = 0.45, linetype = "longdash", alpha = 0.7) +
  scale_x_datetime(breaks = primeros_meses_ct, labels = etiquetas_meses,
                   minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))) +
  facet_wrap(~ ESTACION, ncol = 4, scales = "free_y") +
  labs(
    title    = "Series Temporales Horarias de NO\u2082 \u2014 Todas las Estaciones",
    subtitle = "Madrid 2025  \u00b7  L\u00ednea roja: l\u00edmite anual OMS (40 \u00b5g/m\u00b3)",
    x        = NULL,
    y        = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "gray40", size = 10),
    strip.text         = element_text(face = "bold", size = 7.5),
    axis.text.x        = element_text(angle = 60, hjust = 1, size = 6),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(0.8, "lines")
  )

ggsave(
  file.path(here("outputs", "series_temporales"), "panel_horario_todas_estaciones.png"),
  plot = p_panel_horario, width = 22, height = 18, dpi = 200, bg = "white"
)
cat("  \u2713 panel_horario_todas_estaciones.png\n")

# ==============================================================================
# RESUMEN FINAL
# ==============================================================================

n_d <- length(list.files(dir_diario,  pattern = "\\.png$"))
n_h <- length(list.files(dir_horario, pattern = "\\.png$"))

cat(sprintf("\n\u2705 Script completado.\n"))
cat(sprintf("   Diario  : %d PNGs en outputs/series_temporales/diario/\n",  n_d))
cat(sprintf("   Horario : %d PNGs en outputs/series_temporales/horario/\n", n_h))
cat(sprintf("   Paneles : 2 PNGs en outputs/series_temporales/\n"))

