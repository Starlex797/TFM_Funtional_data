# ==============================================================================
# SERIES TEMPORALES CLIMÁTICAS — DIARIO Y HORARIO POR ESTACIÓN
# Madrid 2025 · 6 variables · Delimitadores mensuales
# Outputs: outputs/series_temporales_clima/
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CONFIGURACIÓN DE VARIABLES
# ==============================================================================

vars_clima <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  "Presion Barométrica", "Radiación Solar", "Velocidad Viento"
)

labels_clima <- c(
  "Temperatura (°C)",
  "Humedad Relativa (%)",
  "Precipitaciones (mm)",
  "Presión Barométrica (hPa)",
  "Radiación Solar (W/m²)",
  "Velocidad Viento (m/s)"
)

# Mapa nombre variable → etiqueta legible
dict_labels <- setNames(labels_clima, vars_clima)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

dt_diario <- readRDS(here("data", "processed", "Clima", "diario",
                           "meteo_madrid_2025_diario.rds"))

dt_horario <- readRDS(here("data", "processed", "Clima", "horario",
                            "meteo_madrid_2025_horario.rds"))

# Variables disponibles en los datos (intersección con vars definidas)
vars_diario  <- intersect(vars_clima, names(dt_diario))
vars_horario <- intersect(vars_clima, names(dt_horario))

cat(sprintf("Variables diario  : %s\n", paste(vars_diario,  collapse = ", ")))
cat(sprintf("Variables horario : %s\n", paste(vars_horario, collapse = ", ")))

# ==============================================================================
# 3. TRANSFORMAR A FORMATO LARGO
# ==============================================================================

# --- Diario ---
dt_diario_largo <- melt(
  dt_diario,
  id.vars      = c("ESTACION", "FECHA"),
  measure.vars = vars_diario,
  variable.name = "Variable",
  value.name    = "Valor"
)
dt_diario_largo[, Variable   := as.character(Variable)]
dt_diario_largo[, Label      := dict_labels[Variable]]
dt_diario_largo[, Label      := factor(Label, levels = labels_clima)]
setorder(dt_diario_largo, ESTACION, Variable, FECHA)

# --- Horario ---
dt_horario[, Hora_num := as.integer(sub("H0?", "", as.character(HORA)))]
dt_horario[, DATETIME := as.POSIXct(FECHA) + (Hora_num - 1L) * 3600L]

dt_horario_largo <- melt(
  dt_horario,
  id.vars      = c("ESTACION", "FECHA", "HORA", "DATETIME"),
  measure.vars = vars_horario,
  variable.name = "Variable",
  value.name    = "Valor"
)
dt_horario_largo[, Variable := as.character(Variable)]
dt_horario_largo[, Label    := dict_labels[Variable]]
dt_horario_largo[, Label    := factor(Label, levels = labels_clima)]
setorder(dt_horario_largo, ESTACION, Variable, DATETIME)

# ==============================================================================
# 4. DELIMITADORES MENSUALES
# ==============================================================================

primeros_meses    <- seq(as.Date("2025-01-01"), as.Date("2025-12-01"), by = "1 month")
primeros_meses_ct <- as.POSIXct(primeros_meses)
etiquetas_meses   <- format(primeros_meses, "%b")

# Franjas alternadas (mes par con fondo gris)
franjas_diario <- data.frame(
  xmin  = primeros_meses,
  xmax  = c(primeros_meses[-1], as.Date("2026-01-01")),
  shade = (seq_along(primeros_meses) %% 2L == 0L)
) |> subset(shade)

meses_limites_ct <- c(primeros_meses_ct, as.POSIXct("2026-01-01"))
franjas_horario <- data.frame(
  xmin  = meses_limites_ct[-length(meses_limites_ct)],
  xmax  = meses_limites_ct[-1],
  shade = (seq_along(primeros_meses) %% 2L == 0L)
) |> subset(shade)

# ==============================================================================
# 5. DIRECTORIOS DE SALIDA
# ==============================================================================

dir_base    <- here("outputs", "series_temporales_clima")
dir_diario  <- file.path(dir_base, "diario")
dir_horario <- file.path(dir_base, "horario")
dir_paneles <- file.path(dir_base, "paneles")

for (d in c(dir_diario, dir_horario, dir_paneles))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 6. FUNCIÓN AUXILIAR
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", tolower(x))
}

# Paleta fija por variable (orden según labels_clima)
paleta_vars <- c(
  "Temperatura (°C)"          = "#e74c3c",
  "Humedad Relativa (%)"      = "#2980b9",
  "Precipitaciones (mm)"      = "#27ae60",
  "Presión Barométrica (hPa)" = "#8e44ad",
  "Radiación Solar (W/m²)"    = "#f39c12",
  "Velocidad Viento (m/s)"    = "#16a085"
)

# ==============================================================================
# 7. BUCLE — SERIES DIARIAS POR ESTACIÓN (1 PNG por estación, 6 facetas)
# ==============================================================================

estaciones_unicas <- sort(unique(dt_diario_largo$ESTACION))
cat(sprintf("\nEstaciones: %d\n", length(estaciones_unicas)))

cat("\n--- Generando series DIARIAS por estación ---\n")

for (est in estaciones_unicas) {

  dt_est <- dt_diario_largo[ESTACION == est & !is.na(Valor)]

  if (nrow(dt_est) == 0L) {
    cat(sprintf("  [SKIP] %s\n", est))
    next
  }

  n_vars_ok <- uniqueN(dt_est$Variable)

  p <- ggplot(dt_est, aes(x = FECHA, y = Valor, color = Label)) +

    # Franjas mensuales de fondo
    geom_rect(
      data        = franjas_diario,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.55
    ) +

    # Líneas divisorias mensuales
    geom_vline(
      xintercept  = as.numeric(primeros_meses),
      color       = "gray60", linewidth = 0.35, linetype = "dashed"
    ) +

    # Serie temporal
    geom_line(linewidth = 0.55, alpha = 0.85) +
    geom_point(size = 0.8, alpha = 0.5) +

    scale_color_manual(values = paleta_vars, guide = "none") +

    scale_x_date(
      breaks       = primeros_meses,
      labels       = etiquetas_meses,
      minor_breaks = NULL,
      expand       = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +

    facet_wrap(~ Label, ncol = 2, scales = "free_y") +

    labs(
      title    = paste0("Series Temporales Diarias \u2014 ", est),
      subtitle = "Madrid 2025  \u00b7  Variables clim\u00e1ticas",
      x        = NULL,
      y        = "Valor",
      caption  = "Franjas grises: meses pares  \u00b7  L\u00edneas discontinuas: inicio de mes"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "gray40", size = 9),
      plot.caption       = element_text(color = "gray55", size = 7.5),
      strip.text         = element_text(face = "bold", size = 9),
      strip.background   = element_rect(fill = "gray96", color = "gray80"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(1, "lines")
    )

  nombre_file <- file.path(dir_diario,
                            paste0("diario_", limpiar_nombre(est), ".png"))
  ggsave(nombre_file, plot = p, width = 14, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 8. BUCLE — SERIES HORARIAS POR ESTACIÓN (1 PNG por estación, 6 facetas)
# ==============================================================================

cat("\n--- Generando series HORARIAS por estación ---\n")

for (est in estaciones_unicas) {

  dt_est <- dt_horario_largo[ESTACION == est & !is.na(Valor)]

  if (nrow(dt_est) == 0L) {
    cat(sprintf("  [SKIP] %s\n", est))
    next
  }

  n_horas <- nrow(dt_est)

  p <- ggplot(dt_est, aes(x = DATETIME, y = Valor, color = Label)) +

    # Franjas mensuales de fondo
    geom_rect(
      data        = franjas_horario,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.55
    ) +

    # Líneas divisorias mensuales
    geom_vline(
      xintercept  = primeros_meses_ct,
      color       = "gray60", linewidth = 0.3, linetype = "dashed"
    ) +

    # Serie horaria
    geom_line(linewidth = 0.3, alpha = 0.75) +

    scale_color_manual(values = paleta_vars, guide = "none") +

    scale_x_datetime(
      breaks       = primeros_meses_ct,
      labels       = etiquetas_meses,
      minor_breaks = NULL,
      expand       = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +

    facet_wrap(~ Label, ncol = 2, scales = "free_y") +

    labs(
      title    = paste0("Series Temporales Horarias \u2014 ", est),
      subtitle = "Madrid 2025  \u00b7  Variables clim\u00e1ticas",
      x        = NULL,
      y        = "Valor",
      caption  = "Franjas grises: meses pares  \u00b7  L\u00edneas discontinuas: inicio de mes"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "gray40", size = 9),
      plot.caption       = element_text(color = "gray55", size = 7.5),
      strip.text         = element_text(face = "bold", size = 9),
      strip.background   = element_rect(fill = "gray96", color = "gray80"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(1, "lines")
    )

  nombre_file <- file.path(dir_horario,
                            paste0("horario_", limpiar_nombre(est), ".png"))
  ggsave(nombre_file, plot = p, width = 14, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 9. PANELES RESUMEN — TODAS LAS ESTACIONES POR VARIABLE (diario y horario)
# ==============================================================================

cat("\n--- Generando paneles resumen por variable ---\n")

for (v in vars_diario) {

  lbl <- dict_labels[[v]]

  # ── Panel diario ────────────────────────────────────────────────────────────
  dt_v <- dt_diario_largo[Variable == v & !is.na(Valor)]

  if (nrow(dt_v) > 0L) {
    p_d <- ggplot(dt_v, aes(x = FECHA, y = Valor)) +

      geom_rect(
        data        = franjas_diario,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "gray88", alpha = 0.55
      ) +
      geom_vline(xintercept = as.numeric(primeros_meses),
                 color = "gray60", linewidth = 0.3, linetype = "dashed") +
      geom_line(color = paleta_vars[[lbl]], linewidth = 0.45, alpha = 0.85) +
      scale_x_date(breaks = primeros_meses, labels = etiquetas_meses,
                   minor_breaks = NULL, expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      facet_wrap(~ ESTACION, ncol = 4, scales = "free_y") +
      labs(
        title    = paste0("Serie Diaria: ", lbl, " \u2014 Todas las Estaciones"),
        subtitle = "Madrid 2025",
        x        = NULL,
        y        = lbl
      ) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(color = "gray40", size = 9),
        strip.text         = element_text(face = "bold", size = 7),
        axis.text.x        = element_text(angle = 60, hjust = 1, size = 6),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.spacing      = unit(0.7, "lines")
      )

    nombre_file <- file.path(dir_paneles,
                              paste0("panel_diario_", limpiar_nombre(v), ".png"))
    ggsave(nombre_file, plot = p_d, width = 22, height = 18, dpi = 200, bg = "white")
    cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
  }

  # ── Panel horario ───────────────────────────────────────────────────────────
  dt_vh <- dt_horario_largo[Variable == v & !is.na(Valor)]

  if (nrow(dt_vh) > 0L) {
    p_h <- ggplot(dt_vh, aes(x = DATETIME, y = Valor)) +

      geom_rect(
        data        = franjas_horario,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "gray88", alpha = 0.55
      ) +
      geom_vline(xintercept = primeros_meses_ct,
                 color = "gray60", linewidth = 0.25, linetype = "dashed") +
      geom_line(color = paleta_vars[[lbl]], linewidth = 0.25, alpha = 0.75) +
      scale_x_datetime(breaks = primeros_meses_ct, labels = etiquetas_meses,
                       minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      facet_wrap(~ ESTACION, ncol = 4, scales = "free_y") +
      labs(
        title    = paste0("Serie Horaria: ", lbl, " \u2014 Todas las Estaciones"),
        subtitle = "Madrid 2025",
        x        = NULL,
        y        = lbl
      ) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(color = "gray40", size = 9),
        strip.text         = element_text(face = "bold", size = 7),
        axis.text.x        = element_text(angle = 60, hjust = 1, size = 6),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.spacing      = unit(0.7, "lines")
      )

    nombre_file <- file.path(dir_paneles,
                              paste0("panel_horario_", limpiar_nombre(v), ".png"))
    ggsave(nombre_file, plot = p_h, width = 22, height = 18, dpi = 200, bg = "white")
    cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
  }
}

# ==============================================================================
# 10. RESUMEN FINAL
# ==============================================================================

n_d <- length(list.files(dir_diario,  pattern = "\\.png$"))
n_h <- length(list.files(dir_horario, pattern = "\\.png$"))
n_p <- length(list.files(dir_paneles, pattern = "\\.png$"))

cat(sprintf("\n\u2705 Script completado.\n"))
cat(sprintf("   Diario   : %d PNGs en outputs/series_temporales_clima/diario/\n",  n_d))
cat(sprintf("   Horario  : %d PNGs en outputs/series_temporales_clima/horario/\n", n_h))
cat(sprintf("   Paneles  : %d PNGs en outputs/series_temporales_clima/paneles/\n", n_p))
cat(sprintf("   TOTAL    : %d archivos PNG\n", n_d + n_h + n_p))

