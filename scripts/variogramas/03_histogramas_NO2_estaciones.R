# ==============================================================================
# HISTOGRAMAS DE NO₂ POR ESTACIÓN — ESCALA CRUDA VS LOG
# Madrid 2025 · Escala diaria
# Salida: un PNG por estación con dos paneles (crudo | log)
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA Y PREPARACIÓN
# ==============================================================================

dt_no2 <- readRDS(here("data", "processed", "contaminacion", "diario",
                        "aire_madrid_2025_No2_trans_diarios.rds"))

# Formato largo: una fila por (estación, día, escala)
dt_largo <- melt(
  dt_no2[!is.na(DATO_DIARIO) & !is.na(LOG_NO2_DIARIO)],
  id.vars      = "ESTACION",
  measure.vars = c("DATO_DIARIO", "LOG_NO2_DIARIO"),
  variable.name = "Escala",
  value.name    = "Valor"
)

dt_largo[, Escala := factor(
  Escala,
  levels = c("DATO_DIARIO", "LOG_NO2_DIARIO"),
  labels = c("NO\u2082 crudo (\u00b5g/m\u00b3)", "log(NO\u2082)")
)]

# Estadísticas por estación y escala
dt_stats <- dt_largo[, .(
  Media   = mean(Valor),
  Mediana = median(Valor),
  N       = .N,
  SD      = round(sd(Valor), 2)
), by = .(ESTACION, Escala)]

colores_escala <- c(
  "NO\u2082 crudo (\u00b5g/m\u00b3)" = "#2980b9",
  "log(NO\u2082)"                    = "#27ae60"
)

# ==============================================================================
# 2. CARPETA DE SALIDA
# ==============================================================================

carpeta_out <- here("outputs", "histograma")
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 3. LOOP: UN PNG POR ESTACIÓN
# ==============================================================================

estaciones_lista <- unique(dt_largo$ESTACION)
n_estaciones     <- length(estaciones_lista)

cat(sprintf("Generando %d histogramas...\n", n_estaciones))

for (est in estaciones_lista) {

  dt_est <- dt_largo[ESTACION == est]
  dt_st  <- dt_stats[ESTACION == est]

  p <- ggplot(dt_est, aes(x = Valor, fill = Escala)) +
    geom_histogram(
      bins      = 25,
      color     = "white",
      alpha     = 0.85,
      linewidth = 0.3
    ) +
    # Media
    geom_vline(
      data      = dt_st,
      aes(xintercept = Media),
      color     = "#c0392b",
      linewidth = 0.9,
      linetype  = "dashed"
    ) +
    # Mediana
    geom_vline(
      data      = dt_st,
      aes(xintercept = Mediana),
      color     = "#e67e22",
      linewidth = 0.9,
      linetype  = "dotted"
    ) +
    # Anotación estadísticas en cada panel
    geom_label(
      data = dt_st,
      aes(
        x     = Inf,
        y     = Inf,
        label = paste0(
          "Media:   ", round(Media,   2), "\n",
          "Mediana: ", round(Mediana, 2), "\n",
          "DE:      ", SD,              "\n",
          "n = ",      N
        )
      ),
      hjust        = 1.05,
      vjust        = 1.05,
      size         = 3.5,
      label.size   = 0.3,
      fill         = "white",
      color        = "gray30",
      inherit.aes  = FALSE
    ) +
    scale_fill_manual(values = colores_escala, guide = "none") +
    facet_wrap(~ Escala, scales = "free", ncol = 2) +
    labs(
      title    = paste("Estaci\u00f3n:", est),
      subtitle = "Distribuci\u00f3n de NO\u2082 \u00b7 Madrid 2025 \u00b7 Escala cruda vs transformaci\u00f3n log\nL\u00ednea roja discontinua: media  \u00b7  L\u00ednea naranja punteada: mediana",
      x        = NULL,
      y        = "Frecuencia (d\u00edas)",
      caption  = "Fuente: Red de Monitoreo de Calidad del Aire de Madrid \u00b7 Media diaria de registros horarios v\u00e1lidos"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title        = element_text(face = "bold", size = 16),
      plot.subtitle     = element_text(color = "gray40", size = 10, lineheight = 1.4),
      plot.caption      = element_text(color = "gray55", size = 8),
      strip.text        = element_text(face = "bold", size = 13),
      strip.background  = element_rect(fill = "gray93", color = "gray75"),
      panel.grid.minor  = element_blank(),
      panel.spacing     = unit(1.5, "lines"),
      axis.text         = element_text(size = 10)
    )

  nombre_archivo <- paste0("hist_NO2_", gsub("[^a-zA-Z0-9]", "_", est), ".png")

  ggsave(
    filename = file.path(carpeta_out, nombre_archivo),
    plot     = p,
    width    = 12,
    height   = 6,
    dpi      = 200
  )
}

cat(sprintf("\u2705 %d histogramas guardados en:\n   %s\n", n_estaciones, carpeta_out))

