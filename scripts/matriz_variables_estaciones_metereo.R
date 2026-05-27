library(data.table)
library(ggplot2)
library(here)
library(data.table)
vars_clima <- setdiff(names(datos_metereo_horarios), c("ESTACION", "FECHA", "HORA"))

# % de datos disponibles por estación y variable
mat <- datos_metereo_horarios[, lapply(.SD, function(x) round(mean(!is.na(x)) * 100, 1)),
                               by = ESTACION, .SDcols = vars_clima]

# Formato largo para el heatmap
mat_long <- melt(mat, id.vars = "ESTACION", variable.name = "Variable", value.name = "pct_disponible")

grafico_estaciones_metereo_variables<-ggplot(mat_long, aes(x = Variable, y = ESTACION, fill = pct_disponible)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(pct_disponible, "%")), size = 2.8) +
  scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                       midpoint = 50, limits = c(0, 100),
                       name = "% disponible") +
  labs(title = "Disponibilidad de variables meteorológicas por estación",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        axis.text.y = element_text(size = 8))

# Guardar el gráfico

# ==============================================================================
# GUARDADO AUTOMÁTICO DEL GRÁFICO EN LA CARPETA OUTPUT
# ==============================================================================

# 1. Aseguramos que la carpeta de destino existe (si no existe, R la crea)
dir.create(here("output", "figures"), showWarnings = FALSE, recursive = TRUE)

# 2. Guardar en formato PNG (Ideal para presentaciones de PowerPoint o GitHub)
ggsave(
  filename = here("output", "figures", "disponibilidad_variables_meteo.png"),
  plot = grafico_estaciones_metereo_variables,
  width = 9,       # Ancho en pulgadas (ajusta según veas el solapamiento del texto)
  height = 7,      # Alto en pulgadas
  dpi = 300,       # Alta resolución (calidad de imprenta requerida para el TFM)
  bg = "white"     # Fondo blanco para evitar transparencias raras
)


cat("💾 ¡Gráficos guardados con éxito en la carpeta 'output/figures/'!\n")
