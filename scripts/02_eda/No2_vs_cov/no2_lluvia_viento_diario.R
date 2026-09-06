# ==============================================================================
# Efecto de la lluvia y el viento sobre el NO2 (escala DIARIA, 2019-2025)
# ==============================================================================
# La precipitación y el viento actúan a escala de evento (día), no estacional,
# por lo que una serie temporal agregada no muestra su patrón. Aquí se usan:
#   - Escala DIARIA (media de Madrid por día).
#   - Gráficas de RELACIÓN (boxplots), no series temporales:
#       * NO2 en días de lluvia vs secos.
#       * NO2 por cuantiles de viento.
#   - Versión con ANOMALÍA de NO2 (día - media del mes) para eliminar el
#     confusor estacional y aislar el efecto del evento.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

DIR_SALIDA <- here("outputs", "figures", "no2_clima_2019_2025")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
ANIOS <- 2019:2025
UMBRAL_LLUVIA <- 1 # mm/día (media de Madrid) para considerar día de lluvia

# ------------------------------------------------------------------------------
# 1. Cargar maestros diarios y agregar a media de Madrid por día
# ------------------------------------------------------------------------------

lista <- lapply(ANIOS, function(a) {
  f <- here("data", "processed", "Maestro", "diario",
            sprintf("dataset_maestro_inla_%d_DIARIO.rds", a))
  if (!file.exists(f)) return(NULL)
  d <- as.data.table(readRDS(f))
  d[, .(FECHA, NO2 = DATO_DIARIO,
        precip = Precipitaciones_raw, viento = `Velocidad Viento_raw`)]
})
dt <- rbindlist(lista, use.names = TRUE)

dia <- dt[, .(
  NO2 = mean(NO2, na.rm = TRUE),
  precip = mean(precip, na.rm = TRUE),
  viento = mean(viento, na.rm = TRUE)
), by = FECHA]
dia <- dia[is.finite(NO2) & is.finite(precip) & is.finite(viento)]

# ------------------------------------------------------------------------------
# 2. Categorías y anomalía de NO2 (sin ciclo estacional)
# ------------------------------------------------------------------------------

dia[, anio_mes := format(FECHA, "%Y-%m")]
dia[, NO2_anom := NO2 - mean(NO2), by = anio_mes] # día - media de su mes

dia[, Lluvia := factor(
  fifelse(precip >= UMBRAL_LLUVIA, "Lluvia (≥1 mm)", "Seco (<1 mm)"),
  levels = c("Seco (<1 mm)", "Lluvia (≥1 mm)")
)]
dia[, Viento_q := cut(
  viento, quantile(viento, 0:4 / 4, na.rm = TRUE),
  labels = c("Q1 flojo", "Q2", "Q3", "Q4 fuerte"), include.lowest = TRUE
)]

cat(sprintf("Días: %d | de lluvia: %d | secos: %d\n",
            nrow(dia), sum(dia$Lluvia == "Lluvia (≥1 mm)"),
            sum(dia$Lluvia == "Seco (<1 mm)")))
cat("\nMediana de NO2 (crudo / anomalía) por categoría:\n")
print(dia[, .(NO2 = round(median(NO2), 1), Anom = round(median(NO2_anom), 1)), by = Lluvia])
print(dia[, .(NO2 = round(median(NO2), 1), Anom = round(median(NO2_anom), 1)), by = Viento_q][order(Viento_q)])

# ------------------------------------------------------------------------------
# 3. Gráfica de 4 paneles
# ------------------------------------------------------------------------------

col_lluvia <- c("Seco (<1 mm)" = "#fdae61", "Lluvia (≥1 mm)" = "#4575b4")
box_theme <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "none",
        panel.grid.major.x = element_blank())

p1 <- ggplot(dia, aes(Lluvia, NO2, fill = Lluvia)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6) +
  scale_fill_manual(values = col_lluvia) +
  labs(title = "NO₂ según lluvia", x = NULL, y = "NO₂ (µg/m³)") + box_theme

p2 <- ggplot(dia, aes(Viento_q, NO2, fill = Viento_q)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6) +
  scale_fill_brewer(palette = "Greens") +
  labs(title = "NO₂ según viento", x = NULL, y = NULL) + box_theme

p3 <- ggplot(dia, aes(Lluvia, NO2_anom, fill = Lluvia)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_boxplot(outlier.size = 0.4, width = 0.6) +
  scale_fill_manual(values = col_lluvia) +
  labs(title = "Anomalía de NO₂ según lluvia", x = NULL,
       y = "Anomalía NO₂ (µg/m³)") + box_theme

p4 <- ggplot(dia, aes(Viento_q, NO2_anom, fill = Viento_q)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_boxplot(outlier.size = 0.4, width = 0.6) +
  scale_fill_brewer(palette = "Greens") +
  labs(title = "Anomalía de NO₂ según viento", x = NULL, y = NULL) + box_theme

fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Efecto de la lluvia y el viento sobre el NO₂ (diario, 2019-2025)",
    subtitle = paste(
      "Arriba: NO₂ crudo. Abajo: anomalía (día − media de su mes),",
      "que elimina el ciclo estacional y aísla el efecto del evento."
    ),
    caption = "Media de Madrid por día. Viento en cuartiles (Q1 flojo → Q4 fuerte).",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey30"))
  )

ggsave(file.path(DIR_SALIDA, "no2_lluvia_viento_diario.png"),
       plot = fig, width = 11, height = 8, dpi = 200, bg = "white")
fwrite(dia, file.path(DIR_SALIDA, "no2_lluvia_viento_diario.csv"))
cat("\nGráfica guardada en:\n", DIR_SALIDA, "\n")
