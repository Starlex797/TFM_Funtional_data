# ==============================================================================
# NO2 y covariables climáticas 2019-2025 (semanal y mensual)
# ==============================================================================
# Igual que el análisis de tráfico pero con el clima. El clima es casi uniforme
# en toda la ciudad, así que se usa la media de Madrid (no se separa por
# estación). Se distinguen dos grupos:
#   - MACRO (ciclo estacional, se ven en medias mensuales):
#       temperatura, humedad, radiación, presión.
#   - MICRO (efecto a escala diaria/de evento, apenas visible en medias):
#       precipitación, velocidad del viento.
# Gráfica: paneles apilados (NO2 arriba + covariables), compartiendo el eje
# temporal, con el confinamiento sombreado.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

DIR_SALIDA <- here("outputs", "figures", "no2_clima_2019_2025")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
ANIOS <- 2019:2025
COVID_INI <- as.Date("2020-03-14")
COVID_FIN <- as.Date("2020-06-21")

CLIMVARS <- c(
  "Temperatura_raw", "Humedad_Relativa_raw", "Radiación Solar_raw",
  "Presion Barométrica_raw", "Precipitaciones_raw", "Velocidad Viento_raw"
)
VARS_MACRO <- c("Temperatura_raw", "Humedad_Relativa_raw",
                "Radiación Solar_raw", "Presion Barométrica_raw")
VARS_MICRO <- c("Precipitaciones_raw", "Velocidad Viento_raw")

ETQ <- c(
  "NO2" = "NO₂ (µg/m³)",
  "Temperatura_raw" = "Temperatura (°C)",
  "Humedad_Relativa_raw" = "Humedad rel. (%)",
  "Radiación Solar_raw" = "Radiación (W/m²)",
  "Presion Barométrica_raw" = "Presión (hPa)",
  "Precipitaciones_raw" = "Precipitación (mm/día)",
  "Velocidad Viento_raw" = "Viento (m/s)"
)
COLORES <- c(
  "NO₂ (µg/m³)" = "#b2182b", "Temperatura (°C)" = "#ef8a62",
  "Humedad rel. (%)" = "#2166ac", "Radiación (W/m²)" = "#e08214",
  "Presión (hPa)" = "#7570b3", "Precipitación (mm/día)" = "#4575b4",
  "Viento (m/s)" = "#1b7837"
)

# ------------------------------------------------------------------------------
# 1. Cargar y apilar los maestros diarios (NO2 diario + clima)
# ------------------------------------------------------------------------------

lista <- lapply(ANIOS, function(a) {
  f <- here("data", "processed", "Maestro", "diario",
            sprintf("dataset_maestro_inla_%d_DIARIO.rds", a))
  if (!file.exists(f)) return(NULL)
  d <- as.data.table(readRDS(f))
  d[, c("FECHA", "DATO_DIARIO", CLIMVARS), with = FALSE]
})
dt <- rbindlist(lista, use.names = TRUE)

dt[, semana := as.Date(cut(FECHA, breaks = "week"))]
dt[, mes := as.Date(format(FECHA, "%Y-%m-01"))]

# ------------------------------------------------------------------------------
# 2. Agregación: media de Madrid por periodo (NO2 + clima)
# ------------------------------------------------------------------------------

agregar_long <- function(dt, col_periodo) {
  a <- dt[, c(
    list(NO2 = mean(DATO_DIARIO, na.rm = TRUE)),
    lapply(.SD, mean, na.rm = TRUE)
  ), by = col_periodo, .SDcols = CLIMVARS]
  setnames(a, col_periodo, "periodo")
  melt(a, id.vars = "periodo", variable.name = "Var", value.name = "valor")
}

long_semana <- agregar_long(dt, "semana")
long_mes <- agregar_long(dt, "mes")

# ------------------------------------------------------------------------------
# 3. Gráfica de paneles apilados
# ------------------------------------------------------------------------------

grafica <- function(long, vars_orden, titulo, subtitulo) {
  d <- long[Var %in% c("NO2", vars_orden)]
  niveles <- c("NO2", vars_orden)
  d[, Var := factor(Var, levels = niveles, labels = ETQ[niveles])]

  ggplot(d, aes(periodo, valor, color = Var)) +
    annotate("rect", xmin = COVID_INI, xmax = COVID_FIN,
             ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.15) +
    geom_line(linewidth = 0.55) +
    facet_grid(Var ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = COLORES, guide = "none") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title = titulo, subtitle = subtitulo,
      x = NULL, y = NULL,
      caption = "Franja gris = confinamiento (mar-jun 2020). Media de Madrid (todas las estaciones)."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

guardar <- function(long, vars, sufijo, titulo, subtitulo, alto) {
  p <- grafica(long, vars, titulo, subtitulo)
  ggsave(file.path(DIR_SALIDA, sprintf("no2_clima_%s.png", sufijo)),
         plot = p, width = 12, height = alto, dpi = 200, bg = "white")
}

# MACRO (temperatura, humedad, radiación, presión)
guardar(long_mes, VARS_MACRO, "macro_mensual",
        "NO₂ y clima macro (2019-2025)",
        "Media mensual — variables de ciclo estacional", 10)
guardar(long_semana, VARS_MACRO, "macro_semanal",
        "NO₂ y clima macro (2019-2025)",
        "Media semanal — variables de ciclo estacional", 10)

# MICRO (precipitación, viento)
guardar(long_mes, VARS_MICRO, "micro_mensual",
        "NO₂ y clima micro (2019-2025)",
        "Media mensual — su efecto es a escala diaria, apenas se aprecia aquí", 6.5)
guardar(long_semana, VARS_MICRO, "micro_semanal",
        "NO₂ y clima micro (2019-2025)",
        "Media semanal — su efecto es a escala diaria/de evento", 6.5)

fwrite(long_mes, file.path(DIR_SALIDA, "no2_clima_mensual.csv"))
fwrite(long_semana, file.path(DIR_SALIDA, "no2_clima_semanal.csv"))

cat("Gráficas guardadas en:\n", DIR_SALIDA, "\n")
