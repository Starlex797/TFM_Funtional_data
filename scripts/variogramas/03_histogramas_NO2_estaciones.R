# ==============================================================================
# HISTOGRAMAS DE NO₂ POR ESTACIÓN — ESCALA CRUDA VS LOG
# Madrid 2025 · Escala diaria
# Salida: un PNG por estación con dos paneles (crudo | log)
# ==============================================================================
# Resumen

# Carga los datos diarios de NO₂ de 2025.
# Conserva observaciones que tengan disponibles tanto:el NO₂ crudo;
# el logaritmo del NO₂.

# Transforma los datos a formato largo para comparar ambas escalas.
# Calcula media, mediana, desviación típica y número de observaciones por estación y escala.
# Genera dos clases de gráficos:
# Histograma globalAgrupa las observaciones de todas las estaciones.
# Compara la distribución cruda con la distribución logarítmica.
# Calcula asimetría y curtosis.
# Superpone una distribución normal con la misma media y desviación típica.
# Sirve como diagnóstico visual de normalidad y asimetría.

# Histogramas por estaciónCrea un PNG para cada estación.
# Presenta dos paneles: NO₂ crudo y log(NO₂).
# Señala la media y la mediana.
# Añade las estadísticas principales dentro del gráfico.
# ================================================================================
# Finalidad: # Justificación si es conveniente aplicar una transformación logarítmica antes de modelizar.



library(data.table)
library(ggplot2)
library(here)
library(moments) # skewness(), kurtosis()

# ==============================================================================
# 1. CARGA Y PREPARACIÓN
# ==============================================================================

dt_no2 <- readRDS(here(
  "data", "processed", "contaminacion", "diario",
  "aire_madrid_2025_No2_trans_diarios.rds"
))

# Formato largo: una fila por (estación, día, escala)
dt_largo <- melt(
  dt_no2[!is.na(DATO_DIARIO) & !is.na(LOG_NO2_DIARIO)],
  id.vars = "ESTACION",
  measure.vars = c("DATO_DIARIO", "LOG_NO2_DIARIO"),
  variable.name = "Escala",
  value.name = "Valor"
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
# 2b. HISTOGRAMA GLOBAL — ¿HACE FALTA TRANSFORMAR?
# Todas las estaciones agrupadas · crudo vs log · asimetría + curva normal
# ==============================================================================

# Estadísticas globales (todas las estaciones, ambas escalas)
dt_stats_global <- dt_largo[, .(
  Media     = mean(Valor),
  Mediana   = median(Valor),
  SD        = sd(Valor),
  Asimetria = skewness(Valor),
  Curtosis  = kurtosis(Valor),
  N         = .N
), by = Escala]

cat("\n=============================================================\n")
cat(" DIAGNÓSTICO DE ASIMETRÍA — NO₂ (todas las estaciones)\n")
cat("=============================================================\n")
print(dt_stats_global, digits = 4)
cat("-------------------------------------------------------------\n")
cat(" |Asimetría| < 0.5 aprox. simétrica · > 1 fuertemente asimétrica\n")
cat("-------------------------------------------------------------\n\n")

# Curva normal de referencia (misma media/SD que los datos), precomputada
# por panel para que respete el rango libre de cada facet ("scales = free")
curva_normal <- dt_largo[,
  {
    rango <- range(Valor)
    xs <- seq(rango[1], rango[2], length.out = 200)
    media <- mean(Valor)
    de <- sd(Valor)
    .(x = xs, y = dnorm(xs, mean = media, sd = de))
  },
  by = Escala
]

p_global <- ggplot(dt_largo, aes(x = Valor)) +
  geom_histogram(
    aes(y = after_stat(density), fill = Escala),
    bins = 40,
    color = "white",
    alpha = 0.85,
    linewidth = 0.25
  ) +
  # Curva normal de referencia (misma media/SD que los datos observados)
  geom_line(
    data = curva_normal,
    aes(x = x, y = y),
    color = "#c0392b",
    linewidth = 0.9,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_vline(
    data = dt_stats_global,
    aes(xintercept = Media),
    color = "#c0392b", linewidth = 0.7, linetype = "solid"
  ) +
  geom_label(
    data = dt_stats_global,
    aes(
      x = Inf,
      y = Inf,
      label = paste0(
        "Asimetría: ", round(Asimetria, 2), "\n",
        "Curtosis:  ", round(Curtosis, 2),  "\n",
        "Media:     ", round(Media, 1),     "\n",
        "n = ",        N
      )
    ),
    hjust = 1.05, vjust = 1.05, size = 3.6,
    label.size = 0.3, fill = "white", color = "gray30",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colores_escala, guide = "none") +
  facet_wrap(~Escala, scales = "free", ncol = 2) +
  labs(
    title    = "Distribución global de NO₂ — ¿hace falta transformar?",
    subtitle = "Todas las estaciones agrupadas · Madrid 2025 · Curva roja: normal con la misma media/DE",
    x        = NULL,
    y        = "Densidad",
    caption  = "Fuente: Red de Monitoreo de Calidad del Aire de Madrid · Media diaria de registros horarios válidos"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(color = "gray40", size = 10, lineheight = 1.4),
    plot.caption     = element_text(color = "gray55", size = 8),
    strip.text       = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray93", color = "gray75"),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1.5, "lines"),
    axis.text        = element_text(size = 10)
  )

ggsave(
  filename = file.path(carpeta_out, "00_histograma_global_NO2_crudo_vs_log.png"),
  plot     = p_global,
  width    = 12,
  height   = 6,
  dpi      = 200
)

cat(
  "✅ Histograma global guardado en:\n   ",
  file.path(carpeta_out, "00_histograma_global_NO2_crudo_vs_log.png"), "\n\n"
)


# ==============================================================================
# 3. LOOP: UN PNG POR ESTACIÓN
# ==============================================================================

estaciones_lista <- unique(dt_largo$ESTACION)
n_estaciones <- length(estaciones_lista)

cat(sprintf("Generando %d histogramas...\n", n_estaciones))

for (est in estaciones_lista) {
  dt_est <- dt_largo[ESTACION == est]
  dt_st <- dt_stats[ESTACION == est]

  p <- ggplot(dt_est, aes(x = Valor, fill = Escala)) +
    geom_histogram(
      bins      = 25,
      color     = "white",
      alpha     = 0.85,
      linewidth = 0.3
    ) +
    # Media
    geom_vline(
      data = dt_st,
      aes(xintercept = Media),
      color = "#c0392b",
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    # Mediana
    geom_vline(
      data = dt_st,
      aes(xintercept = Mediana),
      color = "#e67e22",
      linewidth = 0.9,
      linetype = "dotted"
    ) +
    # Anotación estadísticas en cada panel
    geom_label(
      data = dt_st,
      aes(
        x = Inf,
        y = Inf,
        label = paste0(
          "Media:   ", round(Media, 2), "\n",
          "Mediana: ", round(Mediana, 2), "\n",
          "DE:      ", SD, "\n",
          "n = ", N
        )
      ),
      hjust = 1.05,
      vjust = 1.05,
      size = 3.5,
      label.size = 0.3,
      fill = "white",
      color = "gray30",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = colores_escala, guide = "none") +
    facet_wrap(~Escala, scales = "free", ncol = 2) +
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
