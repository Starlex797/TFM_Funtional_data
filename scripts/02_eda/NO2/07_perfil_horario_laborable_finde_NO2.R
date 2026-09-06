# ==============================================================================
# PERFIL HORARIO DE NO₂ — LABORABLES vs FIN DE SEMANA POR ESTACIÓN DEL AÑO
# Madrid 2025 · Todas las estaciones de medición · Escala: hora del día (0–23)
# Outputs: outputs/analysis/estudio_a_nivel_horario/
# ==============================================================================
# QUÉ HACE:
#
# - Carga los registros horarios de NO2 de 2025.
# - Convierte las variables H01–H24 en horas comprendidas entre 0 y 23.
# - Clasifica cada observación según:
#     · Día laborable.
#     · Fin de semana.
# - Clasifica también cada observación según la estación del año:
#     · Invierno.
#     · Primavera.
#     · Verano.
#     · Otoño.
# - Calcula el NO2 horario medio para cada combinación de:
#     · Estación de medición.
#     · Tipo de día.
#     · Estación del año.
#     · Hora del día.
#
# El script genera tres gráficos:
#
# 1. Perfil horario de NO2 durante los días laborables.
# 2. Perfil horario de NO2 durante los fines de semana.
# 3. Comparación directa entre laborables y fines de semana.
#
# En los gráficos:
#
# - Se muestran las curvas de las diferentes estaciones de medición.
# - Se calcula la mediana del NO2 entre estaciones.
# - Se representa el rango intercuartílico como medida de variabilidad espacial.
# - Se separan las cuatro estaciones del año mediante facetas.
#
# FINALIDAD PARA EL TFM:
#
# Este script permite identificar la estructura cíclica diaria del NO2 y su
# relación con los patrones de actividad humana.
#
# Ayuda a detectar:
#
# - Picos de NO2 relacionados con las horas punta de tráfico.
# - Diferencias entre jornadas laborales y fines de semana.
# - Cambios estacionales en la forma del ciclo diario.
# - Cambios estacionales en la intensidad de las concentraciones.
# - Heterogeneidad entre estaciones de medición.
#
# Los resultados pueden justificar la inclusión en los modelos de:
#
# - La hora del día.
# - El tipo de día.
# - La estación del año.
# - Interacciones entre hora, tipo de día y estación del año.
#
# OBSERVACIÓN METODOLÓGICA:
#
# El script considera laborables todos los días de lunes a viernes. Los festivos
# que caen entre semana no se excluyen y quedan clasificados como laborables.
#
# SALIDAS:
#
# outputs/analysis/estudio_a_nivel_horario/
#     · perfil_horario_laborables_NO2_2025.png
#     · perfil_horario_finde_NO2_2025.png
#     · perfil_horario_comparacion_NO2_2025.png

library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA DE DATOS
# ==============================================================================

cat("Cargando datos horarios 2025...\n")

dt_h <- readRDS(here(
  "data", "processed", "contaminacion", "horario",
  "aire_madrid_2025_No2_horarios1.rds"
))

# Hora del día (0–23): H01 → 0, H02 → 1, …, H24 → 23
dt_h[, Hora_num := as.integer(sub("H0?", "", as.character(HORA))) - 1L]

# ==============================================================================
# 2. CLASIFICACIÓN: DÍA DE LA SEMANA Y ESTACIÓN DEL AÑO
# ==============================================================================

# Día de la semana (1 = lunes, 7 = domingo)
dt_h[, dia_semana := as.integer(format(FECHA, "%u"))]

# Tipo de día
dt_h[, tipo_dia := fifelse(dia_semana <= 5L, "Weekday", "Weekend")]

# Estación del año (calendario meteorológico)
dt_h[, mes := as.integer(format(FECHA, "%m"))]
dt_h[, estacion_anio := fcase(
  mes %in% c(12L, 1L, 2L),  "Winter",
  mes %in% c(3L, 4L, 5L),   "Spring",
  mes %in% c(6L, 7L, 8L),   "Summer",
  mes %in% c(9L, 10L, 11L), "Autumn"
)]

# Orden de las estaciones del año (circular, empezando en invierno)
dt_h[, estacion_anio := factor(
  estacion_anio,
  levels = c("Winter", "Spring", "Summer", "Autumn")
)]

cat(sprintf("  Total registros: %d\n", nrow(dt_h)))
cat(sprintf("  Estaciones de medición: %d\n", uniqueN(dt_h$ESTACION)))

# ==============================================================================
# 3. PERFIL HORARIO MEDIO — por tipo de día, estación del año y estación
# ==============================================================================

perfil <- dt_h[
  !is.na(DATO),
  .(
    NO2_medio = mean(DATO, na.rm = TRUE),
    n_obs = .N
  ),
  by = .(ESTACION, tipo_dia, estacion_anio, Hora_num)
]

cat(sprintf("  Perfiles calculados: %d filas\n", nrow(perfil)))

# ==============================================================================
# 4. DIRECTORIO DE SALIDA
# ==============================================================================

dir_salida <- here("outputs", "figures", "EDA", "NO2", "Perfil_horario")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Directorio de salida: %s\n\n", dir_salida))

# ==============================================================================
# 5. PARÁMETROS DE ESTILO COMPARTIDOS
# ==============================================================================

# Paleta de estaciones de medición (se asigna automáticamente)
estaciones_unicas <- sort(unique(perfil$ESTACION))
n_est <- length(estaciones_unicas)

# Paleta de colores para las estaciones de medición
set.seed(2408)
paleta_est <- setNames(
  scales::hue_pal()(n_est),
  estaciones_unicas
)

# Etiquetas del eje X cada 3 horas
breaks_hora <- seq(0L, 23L, by = 3L)
etiquetas_hora <- sprintf("%02d:00", breaks_hora)

# Tema base
tema_perfil <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40", size = 10),
      plot.caption = element_text(color = "gray55", size = 8),
      strip.text = element_text(
        face = "bold", size = 11,
        margin = margin(b = 6, t = 6)
      ),
      strip.background = element_rect(fill = "gray96", color = "gray80"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.key.width = unit(1.5, "cm"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.2, "lines"),
      axis.title = element_text(size = 10)
    )
}

# ==============================================================================
# 6. FUNCIÓN GENERADORA DE GRÁFICO
# ==============================================================================

grafico_perfil <- function(datos, tipo, titulo, subtitulo, color_ribbon) {
  # Banda de confianza: media ± 1 sd entre estaciones (variabilidad espacial)
  banda <- datos[, .(
    NO2_q25 = quantile(NO2_medio, 0.25, na.rm = TRUE),
    NO2_q75 = quantile(NO2_medio, 0.75, na.rm = TRUE),
    NO2_med = median(NO2_medio, na.rm = TRUE)
  ), by = .(estacion_anio, Hora_num)]

  ggplot() +

    # Banda intercuartílica (variabilidad entre estaciones)
    geom_ribbon(
      data = banda,
      aes(x = Hora_num, ymin = NO2_q25, ymax = NO2_q75),
      fill = color_ribbon, alpha = 0.18, inherit.aes = FALSE
    ) +

    # Mediana global entre estaciones
    geom_line(
      data = banda,
      aes(x = Hora_num, y = NO2_med),
      color = color_ribbon, linewidth = 1.1, linetype = "dashed",
      inherit.aes = FALSE
    ) +

    # Perfil por estación de medición
    geom_line(
      data = datos,
      aes(x = Hora_num, y = NO2_medio, color = ESTACION),
      linewidth = 0.55, alpha = 0.75
    ) +

    # Facetas por estación del año
    facet_wrap(~estacion_anio, nrow = 2, ncol = 2) +
    scale_x_continuous(
      breaks = breaks_hora,
      labels = etiquetas_hora,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    scale_color_manual(values = paleta_est, name = "Monitoring station") +
    labs(
      title = titulo,
      subtitle = subtitulo,
      x = "Hour of day",
      y = "NO\u2082 (\u00b5g/m\u00b3)",
      caption = paste0(
        "Coloured lines: mean hourly profile by monitoring station  \u00b7  ",
        "Dashed line: network median  \u00b7  ",
        "Band: interquartile range (Q25\u2013Q75)"
      )
    ) +
    tema_perfil() +
    guides(color = guide_legend(ncol = 4, override.aes = list(linewidth = 1.5)))
}

# ==============================================================================
# 7. GRÁFICO 1 — LABORABLES
# ==============================================================================

cat("Generando gráfico: LABORABLES...\n")

datos_lab <- perfil[tipo_dia == "Weekday"]

p_lab <- grafico_perfil(
  datos = datos_lab,
  tipo = "Weekday",
  titulo = "Hourly NO\u2082 Profile on Weekdays \u2014 Madrid, 2025",
  subtitulo = paste0(
    "All monitoring stations  \u00b7  ",
    "Mean hourly concentration by season  \u00b7  ",
    "Monday to Friday"
  ),
  color_ribbon = "#2980b9"
)

archivo_lab <- file.path(dir_salida, "perfil_horario_laborables_NO2_2025.png")
ggsave(archivo_lab,
  plot = p_lab,
  width = 14, height = 10, dpi = 200, bg = "white"
)
cat(sprintf("  \u2713 Guardado: %s\n", basename(archivo_lab)))

# ==============================================================================
# 8. GRÁFICO 2 — FIN DE SEMANA
# ==============================================================================

cat("Generando gráfico: FIN DE SEMANA...\n")

datos_fds <- perfil[tipo_dia == "Weekend"]

p_fds <- grafico_perfil(
  datos = datos_fds,
  tipo = "Weekend",
  titulo = "Hourly NO\u2082 Profile on Weekends \u2014 Madrid, 2025",
  subtitulo = paste0(
    "All monitoring stations  \u00b7  ",
    "Mean hourly concentration by season  \u00b7  ",
    "Saturday and Sunday"
  ),
  color_ribbon = "#27ae60"
)

archivo_fds <- file.path(dir_salida, "perfil_horario_finde_NO2_2025.png")
ggsave(archivo_fds,
  plot = p_fds,
  width = 14, height = 10, dpi = 200, bg = "white"
)
cat(sprintf("  \u2713 Guardado: %s\n", basename(archivo_fds)))

# ==============================================================================
# 9. GRÁFICO 3 — COMPARACIÓN SUPERPUESTA (LABORABLE vs FIN DE SEMANA)
# ==============================================================================

cat("Generando gráfico: COMPARACIÓN SUPERPUESTA...\n")

# Agregado por {tipo_dia, estacion_anio, Hora_num}: mediana y rango IQR entre estaciones
perfil_comp <- perfil[, .(
  NO2_mediana = median(NO2_medio, na.rm = TRUE),
  NO2_q25     = quantile(NO2_medio, 0.25, na.rm = TRUE),
  NO2_q75     = quantile(NO2_medio, 0.75, na.rm = TRUE)
), by = .(tipo_dia, estacion_anio, Hora_num)]

perfil_comp[, estacion_anio := factor(
  estacion_anio,
  levels = c("Winter", "Spring", "Summer", "Autumn")
)]
perfil_comp[, tipo_dia := factor(
  tipo_dia,
  levels = c("Weekday", "Weekend")
)]

paleta_tipo <- c(
  "Weekday" = "#2980b9",
  "Weekend" = "#27ae60"
)

p_comp <- ggplot(
  perfil_comp,
  aes(x = Hora_num, color = tipo_dia, fill = tipo_dia)
) +

  # Banda IQR por tipo de día
  geom_ribbon(
    aes(ymin = NO2_q25, ymax = NO2_q75),
    alpha = 0.15, color = NA
  ) +

  # Línea de mediana
  geom_line(aes(y = NO2_mediana), linewidth = 1.1) +
  facet_wrap(~estacion_anio, nrow = 2, ncol = 2) +
  scale_x_continuous(
    breaks = breaks_hora,
    labels = etiquetas_hora,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  scale_color_manual(values = paleta_tipo, name = NULL) +
  scale_fill_manual(values = paleta_tipo, name = NULL) +
  labs(
    title = "Hourly NO\u2082 Profile: Weekdays vs Weekends \u2014 Madrid, 2025",
    subtitle = paste0(
      "All monitoring stations  \u00b7  ",
      "Median across stations  \u00b7  ",
      "Band: interquartile range (Q25\u2013Q75)"
    ),
    x = "Hour of day",
    y = "NO\u2082 (\u00b5g/m\u00b3)",
    caption = "Blue: Monday\u2013Friday  \u00b7  Green: Saturday and Sunday"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 10),
    plot.caption = element_text(color = "gray55", size = 8.5),
    strip.text = element_text(
      face = "bold", size = 12,
      margin = margin(b = 6, t = 6)
    ),
    strip.background = element_rect(fill = "gray96", color = "gray80"),
    legend.position = "top",
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    axis.title = element_text(size = 10)
  )

archivo_comp <- file.path(dir_salida, "perfil_horario_comparacion_NO2_2025.png")
ggsave(archivo_comp,
  plot = p_comp,
  width = 13, height = 9, dpi = 200, bg = "white"
)
cat(sprintf("  \u2713 Guardado: %s\n", basename(archivo_comp)))

# ==============================================================================
# 10. GRAFICO 4 - PERFIL GENERAL SIN DIVISION POR ESTACIONES
# ==============================================================================

cat("Generando grafico: PERFIL GENERAL...\n")

# Primero se estima el perfil anual de cada estacion de medicion. Despues se
# resume su distribucion mediante la mediana y el rango intercuartilico.
perfil_estacion_general <- dt_h[
  !is.na(DATO),
  .(NO2_medio = mean(DATO)),
  by = .(ESTACION, tipo_dia, Hora_num)
]

perfil_general <- perfil_estacion_general[, .(
  NO2_mediana = median(NO2_medio),
  NO2_q25 = quantile(NO2_medio, 0.25),
  NO2_q75 = quantile(NO2_medio, 0.75)
), by = .(tipo_dia, Hora_num)]

perfil_general[, tipo_dia := factor(
  tipo_dia,
  levels = c("Weekday", "Weekend")
)]

p_general <- ggplot(
  perfil_general,
  aes(x = Hora_num, color = tipo_dia, fill = tipo_dia)
) +
  geom_ribbon(
    aes(ymin = NO2_q25, ymax = NO2_q75),
    alpha = 0.15, color = NA
  ) +
  geom_line(aes(y = NO2_mediana), linewidth = 1.15) +
  scale_x_continuous(
    breaks = breaks_hora,
    labels = etiquetas_hora,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  scale_color_manual(values = paleta_tipo, name = NULL) +
  scale_fill_manual(values = paleta_tipo, name = NULL) +
  labs(
    title = "Overall Hourly NO\u2082 Profile: Weekdays vs Weekends \u2014 Madrid, 2025",
    subtitle = paste0(
      "All seasons and monitoring stations  \u00b7  ",
      "Median across stations  \u00b7  ",
      "Band: interquartile range (Q25\u2013Q75)"
    ),
    x = "Hour of day",
    y = "NO\u2082 (\u00b5g/m\u00b3)",
    caption = "Blue: Monday\u2013Friday  \u00b7  Green: Saturday and Sunday"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 10),
    plot.caption = element_text(color = "gray55", size = 8.5),
    legend.position = "top",
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10)
  )

archivo_general <- file.path(
  dir_salida,
  "perfil_horario_general_laborable_finde_NO2_2025.png"
)
ggsave(
  archivo_general,
  plot = p_general,
  width = 13, height = 6.5, dpi = 200, bg = "white"
)
cat(sprintf("  \u2713 Guardado: %s\n", basename(archivo_general)))

# ==============================================================================
# 11. RESUMEN
# ==============================================================================

cat("\n==============================================================\n")
cat("RESUMEN\n")
cat("==============================================================\n")
cat(sprintf("  Estaciones de medición incluidas : %d\n", n_est))
cat(sprintf("  Estaciones del año               : Winter, Spring, Summer, Autumn\n"))
cat(sprintf(
  "  Registros válidos (laborables)   : %d\n",
  dt_h[tipo_dia == "Weekday" & !is.na(DATO), .N]
))
cat(sprintf(
  "  Registros válidos (finde)        : %d\n",
  dt_h[tipo_dia == "Weekend" & !is.na(DATO), .N]
))
cat(sprintf("\nArchivos guardados en:\n  %s\n", dir_salida))
cat(sprintf("  · %s\n", basename(archivo_lab)))
cat(sprintf("  · %s\n", basename(archivo_fds)))
cat(sprintf("  · %s\n", basename(archivo_comp)))
cat(sprintf("  · %s\n", basename(archivo_general)))
cat("==============================================================\n")
