# ==============================================================================
# Transformaciones de covariables vs log(NO2): crudo vs transformado
# ==============================================================================
# Radiación solar y Precipitación -> log(x+1) (muchos ceros, muy asimétricas).
# Velocidad del viento -> se prueban inversa 1/(x+1) y log(x+1) (relación
# decreciente y convexa; se busca cuál linealiza mejor).
# Se comparan lado a lado el scatter CRUDO y el TRANSFORMADO (log NO2 en Y).
# Escala horaria, todas las estaciones, Madrid 2025.
# Output: outputs/analysis/scatter_no2_covariables/transformaciones/
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

DIR_SALIDA <- here("outputs", "analysis", "scatter_no2_covariables", "transformaciones")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

d <- as.data.table(readRDS(here(
  "data", "processed", "Maestro", "horario",
  "dataset_maestro_inla_2025_HORARIO.rds"
)))
d[, logNO2 := LOG_NO2_HORARIO]

base <- d[!is.na(logNO2), .(
  logNO2,
  rad  = `Radiación Solar_raw`,
  prec = Precipitaciones_raw,
  vie  = `Velocidad Viento_raw`
)]

# Formato largo: cada covariable en su versión cruda y transformada
largo <- rbindlist(list(
  base[!is.na(rad),  .(x = rad,        logNO2, panel = "Radiación solar (cruda)")],
  base[!is.na(rad),  .(x = log(rad+1), logNO2, panel = "Radiación solar: log(x+1)")],
  base[!is.na(prec), .(x = prec,        logNO2, panel = "Precipitación (cruda)")],
  base[!is.na(prec), .(x = log(prec+1), logNO2, panel = "Precipitación: log(x+1)")],
  base[!is.na(vie),  .(x = vie,          logNO2, panel = "Viento (crudo)")],
  base[!is.na(vie),  .(x = 1 / (vie+1),  logNO2, panel = "Viento: inversa 1/(x+1)")],
  base[!is.na(vie),  .(x = log(vie+1),   logNO2, panel = "Viento: log(x+1)")]
))
largo[, panel := factor(panel, levels = c(
  "Radiación solar (cruda)",  "Radiación solar: log(x+1)",
  "Precipitación (cruda)",    "Precipitación: log(x+1)",
  "Viento (crudo)",           "Viento: inversa 1/(x+1)",
  "Viento: log(x+1)"
))]

p <- ggplot(largo, aes(x, logNO2)) +
  geom_point(color = "#0057FF", size = 0.4, alpha = 0.03) +
  facet_wrap(~ panel, scales = "free_x", ncol = 2) +
  labs(
    title = "Transformaciones de covariables frente a log(NO₂) — horario 2025",
    subtitle = paste(
      "Radiación y precipitación: cruda vs log(x+1).",
      "Viento: crudo, inversa 1/(x+1) y log(x+1)."
    ),
    x = NULL, y = "log(NO₂ + 1)",
    caption = "Cada punto es una hora. Todas las estaciones. Eje X libre por panel."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ruta <- file.path(DIR_SALIDA, "transformaciones_covariables_logNO2_2025.png")
ggsave(ruta, p, width = 12, height = 13, dpi = 200, bg = "white")
cat("Guardado:", ruta, "\n")
