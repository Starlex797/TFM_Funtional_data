# ==============================================================================
# Scatter log(NO2) vs radiación y precipitación, descartando los ceros triviales
# ==============================================================================
# Sin transformaciones. Se filtran:
#   - Radiación solar: solo valores > 0 (se descartan las horas nocturnas).
#   - Precipitación:   solo valores > 0 (se descartan las horas sin lluvia).
# Así el scatter muestra la relación en las horas que realmente informan.
# Escala horaria, todas las estaciones, Madrid 2025.
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

DIR_SALIDA <- here("outputs", "analysis", "scatter_no2_covariables", "sin_ceros")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

d <- as.data.table(readRDS(here(
  "data", "processed", "Maestro", "horario",
  "dataset_maestro_inla_2025_HORARIO.rds"
)))
d[, logNO2 := LOG_NO2_HORARIO]

rad <- d[!is.na(logNO2) & !is.na(`Radiación Solar_raw`) & `Radiación Solar_raw` > 0,
         .(x = `Radiación Solar_raw`, logNO2, panel = "Radiación solar (solo horas de día, >0 W/m²)")]
prec <- d[!is.na(logNO2) & !is.na(Precipitaciones_raw) & Precipitaciones_raw > 0,
          .(x = Precipitaciones_raw, logNO2, panel = "Precipitación (solo horas con lluvia, >0 mm)")]

cat(sprintf("Radiación: %d horas de día (de %d totales)\n",
            nrow(rad), sum(!is.na(d$`Radiación Solar_raw`))))
cat(sprintf("Precipitación: %d horas con lluvia (de %d totales)\n",
            nrow(prec), sum(!is.na(d$Precipitaciones_raw))))

# Correlaciones sobre el subconjunto filtrado
cat("\nCorrelaciones (subconjunto sin ceros):\n")
cat(sprintf("  Radiación : r = %.3f | rho = %.3f\n",
            cor(rad$x, rad$logNO2), cor(rad$x, rad$logNO2, method = "spearman")))
cat(sprintf("  Precipit. : r = %.3f | rho = %.3f\n",
            cor(prec$x, prec$logNO2), cor(prec$x, prec$logNO2, method = "spearman")))

largo <- rbindlist(list(rad, prec))
largo[, panel := factor(panel, levels = c(
  "Radiación solar (solo horas de día, >0 W/m²)",
  "Precipitación (solo horas con lluvia, >0 mm)"
))]

p <- ggplot(largo, aes(x, logNO2)) +
  geom_point(color = "#0057FF", size = 0.5, alpha = 0.05) +
  facet_wrap(~ panel, scales = "free_x") +
  labs(
    title = "log(NO₂) vs radiación y precipitación, sin los ceros triviales",
    subtitle = "Sin transformar. Radiación: solo horas de día. Precipitación: solo horas con lluvia.",
    x = NULL, y = "log(NO₂ + 1)",
    caption = "Cada punto es una hora. Todas las estaciones. Madrid 2025."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ruta <- file.path(DIR_SALIDA, "scatter_radiacion_precip_sin_ceros_2025.png")
ggsave(ruta, p, width = 12, height = 6, dpi = 200, bg = "white")
cat("\nGuardado:", ruta, "\n")
