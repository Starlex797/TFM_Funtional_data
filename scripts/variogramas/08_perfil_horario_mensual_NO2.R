# ==============================================================================
# PERFIL HORARIO DE NO₂ — COMPARACIÓN MENSUAL (LABORABLES vs FIN DE SEMANA)
# Madrid 2025 · Todas las estaciones de medición · Escala: hora del día (0–23)
# Outputs: outputs/estudio_a_nivel_horario/
# ==============================================================================
# QUÉ HACE:
#
# - Carga los registros horarios de NO2 de 2025.
# - Convierte las variables H01–H24 en horas comprendidas entre 0 y 23.
# - Clasifica cada observación como:
#     · Laborable.
#     · Fin de semana.
# - Identifica el mes correspondiente a cada observación.
# - Calcula el perfil horario medio para cada estación de medición.
# - Agrega posteriormente los perfiles de las estaciones utilizando:
#     · Mediana.
#     · Primer cuartil.
#     · Tercer cuartil.
#
# El script genera tres gráficos:
#
# 1. Perfiles horarios mensuales para los días laborables.
# 2. Perfiles horarios mensuales para los fines de semana.
# 3. Panel de doce meses comparando laborables y fines de semana.
#
# Las bandas de los gráficos representan el rango intercuartílico entre
# estaciones y, por tanto, ofrecen una medida de la variabilidad espacial.
#
# FINALIDAD PARA EL TFM:
#
# Este script estudia con mayor resolución la evolución del ciclo horario del
# NO2 a lo largo del año.
#
# Permite comprobar:
#
# - Si las horas de los picos de contaminación cambian según el mes.
# - Si la magnitud del NO2 presenta un patrón estacional mensual.
# - Si la diferencia entre laborables y fines de semana es estable.
# - Qué meses presentan mayor variabilidad entre estaciones.
# - Si existen meses con perfiles horarios especialmente atípicos.
#
# Este análisis complementa al script 07:
#
# - El script 07 proporciona una visión estacional resumida.
# - El script 08 proporciona una visión mensual más detallada.
#
# OBSERVACIÓN METODOLÓGICA:
#
# Los festivos que caen entre lunes y viernes quedan clasificados como
# laborables.
#
# SALIDAS:
#
# outputs/estudio_a_nivel_horario/
#     · perfil_mensual_laborables_NO2_2025.png
#     · perfil_mensual_finde_NO2_2025.png
#     · perfil_mensual_comparacion_NO2_2025.png


library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA DE DATOS
# ==============================================================================

cat("Cargando datos horarios 2025...\n")

dt_h <- readRDS(here(
  "data", "processed", "contaminacion", "horario",
  "aire_madrid_2025_No2_horarios.rds"
))

# Hora del día (0–23): H01 → 0, H02 → 1, …, H24 → 23
dt_h[, Hora_num := as.integer(sub("H0?", "", as.character(HORA))) - 1L]

# Tipo de día
dt_h[, dia_semana := as.integer(format(FECHA, "%u"))]
dt_h[, tipo_dia := fifelse(dia_semana <= 5L, "Laborable", "Fin de semana")]
dt_h[, tipo_dia := factor(tipo_dia, levels = c("Laborable", "Fin de semana"))]

# Mes numérico y nombre en español
dt_h[, mes := as.integer(format(FECHA, "%m"))]
dt_h[, mes_nom := factor(
  format(FECHA, "%B"),
  levels = format(
    seq(as.Date("2025-01-01"), by = "1 month", length.out = 12),
    "%B"
  )
)]

cat(sprintf("  Registros cargados : %d\n", nrow(dt_h)))
cat(sprintf("  Estaciones         : %d\n", uniqueN(dt_h$ESTACION)))
cat(sprintf("  Meses              : %d\n", uniqueN(dt_h$mes)))

# ==============================================================================
# 2. PERFILES HORARIOS
# ==============================================================================

# 2a. Por estación de medición (para los plots individuales con ribbons)
perfil_est <- dt_h[
  !is.na(DATO),
  .(NO2_medio = mean(DATO, na.rm = TRUE)),
  by = .(ESTACION, tipo_dia, mes, mes_nom, Hora_num)
]

# 2b. Agregado entre estaciones: mediana + IQR (para líneas centrales)
perfil_mes <- perfil_est[, .(
  NO2_mediana = median(NO2_medio, na.rm = TRUE),
  NO2_q25     = quantile(NO2_medio, 0.25, na.rm = TRUE),
  NO2_q75     = quantile(NO2_medio, 0.75, na.rm = TRUE)
), by = .(tipo_dia, mes, mes_nom, Hora_num)]

setorder(perfil_mes, tipo_dia, mes, Hora_num)

# ==============================================================================
# 3. PARÁMETROS ESTÉTICOS COMPARTIDOS
# ==============================================================================

breaks_hora <- seq(0L, 23L, by = 3L)
etiquetas_hora <- sprintf("%02d:00", breaks_hora)

# Paleta de 12 meses (espectro continuo azul → verde → naranja → rojo)
paleta_meses <- setNames(
  colorRampPalette(c("#1a6faf", "#2ecc71", "#e67e22", "#c0392b"))(12L),
  levels(perfil_mes$mes_nom)
)

paleta_tipo <- c(
  "Laborable" = "#2980b9",
  "Fin de semana" = "#27ae60"
)

tema_mensual <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40", size = 10),
      plot.caption = element_text(color = "gray55", size = 8),
      strip.text = element_text(
        face = "bold", size = 10,
        margin = margin(b = 5, t = 5)
      ),
      strip.background = element_rect(fill = "gray96", color = "gray80"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8.5),
      legend.key.width = unit(1.4, "cm"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.9, "lines"),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
}

# ==============================================================================
# 4. DIRECTORIO DE SALIDA
# ==============================================================================

dir_salida <- here("outputs", "estudio_a_nivel_horario")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 5. PLOT 1 — LABORABLES: 12 líneas coloreadas por mes (panel único)
# ==============================================================================

cat("\nGenerando gráfico 1: LABORABLES — 12 meses superpuestos...\n")

datos_lab_mes <- perfil_mes[tipo_dia == "Laborable"]

p_lab_mes <- ggplot(
  datos_lab_mes,
  aes(
    x = Hora_num, y = NO2_mediana,
    color = mes_nom, group = mes_nom
  )
) +

  # Banda IQR por mes
  geom_ribbon(
    aes(ymin = NO2_q25, ymax = NO2_q75, fill = mes_nom),
    alpha = 0.10, color = NA
  ) +

  # Línea central (mediana entre estaciones)
  geom_line(linewidth = 0.9, alpha = 0.9) +

  # Referencia OMS 24h
  geom_hline(
    yintercept = 25, color = "#e74c3c",
    linewidth = 0.6, linetype = "longdash", alpha = 0.7
  ) +
  annotate("text",
    x = 0.2, y = 26.8,
    label = "Gu\u00eda OMS 24h (25 \u00b5g/m\u00b3)",
    color = "#e74c3c", size = 2.7, hjust = 0, fontface = "italic"
  ) +
  scale_x_continuous(
    breaks = breaks_hora, labels = etiquetas_hora,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  scale_color_manual(values = paleta_meses, name = "Mes") +
  scale_fill_manual(values = paleta_meses, guide = "none") +
  labs(
    title    = "Perfil Horario de NO\u2082 en D\u00edas Laborables por Mes \u2014 Madrid 2025",
    subtitle = "Todas las estaciones de medici\u00f3n  \u00b7  Mediana entre estaciones  \u00b7  Banda: rango IQR (Q25\u2013Q75)",
    x        = "Hora del d\u00eda",
    y        = "NO\u2082 (\u00b5g/m\u00b3)",
    caption  = "Cada l\u00ednea: mediana horaria de NO\u2082 entre estaciones  \u00b7  Lunes\u2013viernes"
  ) +
  tema_mensual() +
  guides(color = guide_legend(nrow = 2, override.aes = list(linewidth = 1.8)))

archivo_lab_mes <- file.path(dir_salida, "perfil_mensual_laborables_NO2_2025.png")
ggsave(archivo_lab_mes,
  plot = p_lab_mes,
  width = 13, height = 7, dpi = 200, bg = "white"
)
cat(sprintf("  \u2713 %s\n", basename(archivo_lab_mes)))

# ==============================================================================
# 6. PLOT 2 — FIN DE SEMANA: 12 líneas coloreadas por mes (panel único)
# ==============================================================================

cat("Generando gráfico 2: FIN DE SEMANA — 12 meses superpuestos...\n")

datos_fds_mes <- perfil_mes[tipo_dia == "Fin de semana"]

p_fds_mes <- ggplot(
  datos_fds_mes,
  aes(
    x = Hora_num, y = NO2_mediana,
    color = mes_nom, group = mes_nom
  )
) +
  geom_ribbon(
    aes(ymin = NO2_q25, ymax = NO2_q75, fill = mes_nom),
    alpha = 0.10, color = NA
  ) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_hline(
    yintercept = 25, color = "#e74c3c",
    linewidth = 0.6, linetype = "longdash", alpha = 0.7
  ) +
  annotate("text",
    x = 0.2, y = 26.8,
    label = "Gu\u00eda OMS 24h (25 \u00b5g/m\u00b3)",
    color = "#e74c3c", size = 2.7, hjust = 0, fontface = "italic"
  ) +
  scale_x_continuous(
    breaks = breaks_hora, labels = etiquetas_hora,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  scale_color_manual(values = paleta_meses, name = "Mes") +
  scale_fill_manual(values = paleta_meses, guide = "none") +
  labs(
    title    = "Perfil Horario de NO\u2082 en Fin de Semana por Mes \u2014 Madrid 2025",
    subtitle = "Todas las estaciones de medici\u00f3n  \u00b7  Mediana entre estaciones  \u00b7  Banda: rango IQR (Q25\u2013Q75)",
    x        = "Hora del d\u00eda",
    y        = "NO\u2082 (\u00b5g/m\u00b3)",
    caption  = "Cada l\u00ednea: mediana horaria de NO\u2082 entre estaciones  \u00b7  S\u00e1bado y domingo"
  ) +
  tema_mensual() +
  guides(color = guide_legend(nrow = 2, override.aes = list(linewidth = 1.8)))

archivo_fds_mes <- file.path(dir_salida, "perfil_mensual_finde_NO2_2025.png")
ggsave(archivo_fds_mes,
  plot = p_fds_mes,
  width = 13, height = 7, dpi = 200, bg = "white"
)
cat(sprintf("  \u2713 %s\n", basename(archivo_fds_mes)))

# ==============================================================================
# 7. PLOT 3 — COMPARACIÓN: 12 facetas (mes), laborable vs fin de semana
# ==============================================================================

cat("Generando gráfico 3: COMPARACIÓN mensual laborable vs fin de semana...\n")

p_comp_mes <- ggplot(
  perfil_mes,
  aes(x = Hora_num, color = tipo_dia, fill = tipo_dia)
) +
  geom_ribbon(
    aes(ymin = NO2_q25, ymax = NO2_q75),
    alpha = 0.15, color = NA
  ) +
  geom_line(aes(y = NO2_mediana), linewidth = 1.0) +
  geom_hline(
    yintercept = 25, color = "#e74c3c",
    linewidth = 0.5, linetype = "longdash", alpha = 0.65
  ) +
  facet_wrap(~mes_nom, nrow = 3, ncol = 4) +
  scale_x_continuous(
    breaks = seq(0L, 21L, by = 6L),
    labels = sprintf("%02d:00", seq(0L, 21L, by = 6L)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_color_manual(values = paleta_tipo, name = NULL) +
  scale_fill_manual(values = paleta_tipo, name = NULL) +
  labs(
    title    = "Perfil Horario de NO\u2082 por Mes: Laborables vs Fin de Semana \u2014 Madrid 2025",
    subtitle = "Todas las estaciones de medici\u00f3n  \u00b7  Mediana entre estaciones  \u00b7  Banda: rango IQR (Q25\u2013Q75)  \u00b7  L\u00ednea roja: gu\u00eda OMS 24h",
    x        = "Hora del d\u00eda",
    y        = "NO\u2082 (\u00b5g/m\u00b3)",
    caption  = "Azul: lunes\u2013viernes  \u00b7  Verde: s\u00e1bado y domingo"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 9.5),
    plot.caption = element_text(color = "gray55", size = 8.5),
    strip.text = element_text(
      face = "bold", size = 10.5,
      margin = margin(b = 5, t = 5)
    ),
    strip.background = element_rect(fill = "gray96", color = "gray80"),
    legend.position = "top",
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.85, "lines"),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

archivo_comp_mes <- file.path(dir_salida, "perfil_mensual_comparacion_NO2_2025.png")
ggsave(archivo_comp_mes,
  plot = p_comp_mes,
  width = 15, height = 11, dpi = 200, bg = "white"
)
cat(sprintf("  \u2713 %s\n", basename(archivo_comp_mes)))

# ==============================================================================
# 8. RESUMEN
# ==============================================================================

cat("\n==============================================================\n")
cat("RESUMEN\n")
cat("==============================================================\n")
cat(sprintf("  Estaciones de medici\u00f3n : %d\n", uniqueN(dt_h$ESTACION)))
cat(sprintf("  Meses                  : 12 (enero\u2013diciembre 2025)\n"))
cat(sprintf("  Archivos guardados en:\n  %s\n", dir_salida))
cat(sprintf("  \u00b7 %s\n", basename(archivo_lab_mes)))
cat(sprintf("  \u00b7 %s\n", basename(archivo_fds_mes)))
cat(sprintf("  \u00b7 %s\n", basename(archivo_comp_mes)))
cat("==============================================================\n")
