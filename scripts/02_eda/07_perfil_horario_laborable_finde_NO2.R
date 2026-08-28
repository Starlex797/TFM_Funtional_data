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
# - Se incluye una línea de referencia de 25 µg/m³.
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
  "aire_madrid_2025_No2_horarios.rds"
))

# Hora del día (0–23): H01 → 0, H02 → 1, …, H24 → 23
dt_h[, Hora_num := as.integer(sub("H0?", "", as.character(HORA))) - 1L]

# ==============================================================================
# 2. CLASIFICACIÓN: DÍA DE LA SEMANA Y ESTACIÓN DEL AÑO
# ==============================================================================

# Día de la semana (1 = lunes, 7 = domingo)
dt_h[, dia_semana := as.integer(format(FECHA, "%u"))]

# Tipo de día
dt_h[, tipo_dia := fifelse(dia_semana <= 5L, "Laborable", "Fin de semana")]

# Estación del año (calendario meteorológico)
dt_h[, mes := as.integer(format(FECHA, "%m"))]
dt_h[, estacion_anio := fcase(
  mes %in% c(12L, 1L, 2L),  "Invierno",
  mes %in% c(3L, 4L, 5L),   "Primavera",
  mes %in% c(6L, 7L, 8L),   "Verano",
  mes %in% c(9L, 10L, 11L), "Otoño"
)]

# Orden de las estaciones del año (circular, empezando en invierno)
dt_h[, estacion_anio := factor(
  estacion_anio,
  levels = c("Invierno", "Primavera", "Verano", "Otoño")
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

dir_salida <- here("outputs", "analysis", "estudio_a_nivel_horario")
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

    # Umbral OMS referencia (25 µg/m³ media 24h — guía 2021)
    geom_hline(
      yintercept = 25, color = "#e74c3c",
      linewidth = 0.6, linetype = "longdash", alpha = 0.7
    ) +
    annotate("text",
      x = 0.2, y = 26.5,
      label = "Guía OMS 24h (25 µg/m³)",
      color = "#e74c3c", size = 2.6, hjust = 0, fontface = "italic"
    ) +

    # Facetas por estación del año
    facet_wrap(~estacion_anio, nrow = 2, ncol = 2) +
    scale_x_continuous(
      breaks = breaks_hora,
      labels = etiquetas_hora,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    scale_color_manual(values = paleta_est, name = "Estación de medición") +
    labs(
      title = titulo,
      subtitle = subtitulo,
      x = "Hora del día",
      y = "NO\u2082 (\u00b5g/m\u00b3)",
      caption = paste0(
        "L\u00edneas de color: perfil horario medio por estaci\u00f3n de medici\u00f3n  \u00b7  ",
        "L\u00ednea discontinua: mediana global  \u00b7  ",
        "Banda: rango intercuart\u00edlico (Q25\u2013Q75)"
      )
    ) +
    tema_perfil() +
    guides(color = guide_legend(ncol = 4, override.aes = list(linewidth = 1.5)))
}

# ==============================================================================
# 7. GRÁFICO 1 — LABORABLES
# ==============================================================================

cat("Generando gráfico: LABORABLES...\n")

datos_lab <- perfil[tipo_dia == "Laborable"]

p_lab <- grafico_perfil(
  datos = datos_lab,
  tipo = "Laborable",
  titulo = "Perfil Horario de NO\u2082 en D\u00edas Laborables \u2014 Madrid 2025",
  subtitulo = paste0(
    "Todas las estaciones de medici\u00f3n  \u00b7  ",
    "Media horaria por estaci\u00f3n del a\u00f1o  \u00b7  ",
    "Lunes a viernes"
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

datos_fds <- perfil[tipo_dia == "Fin de semana"]

p_fds <- grafico_perfil(
  datos = datos_fds,
  tipo = "Fin de semana",
  titulo = "Perfil Horario de NO\u2082 en Fin de Semana \u2014 Madrid 2025",
  subtitulo = paste0(
    "Todas las estaciones de medici\u00f3n  \u00b7  ",
    "Media horaria por estaci\u00f3n del a\u00f1o  \u00b7  ",
    "S\u00e1bado y domingo"
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
  levels = c("Invierno", "Primavera", "Verano", "Otoño")
)]
perfil_comp[, tipo_dia := factor(
  tipo_dia,
  levels = c("Laborable", "Fin de semana")
)]

paleta_tipo <- c(
  "Laborable" = "#2980b9",
  "Fin de semana" = "#27ae60"
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

  # Referencia OMS 24h
  geom_hline(
    yintercept = 25, color = "#e74c3c",
    linewidth = 0.6, linetype = "longdash", alpha = 0.7
  ) +
  annotate("text",
    x = 0.2, y = 26.5,
    label = "Gu\u00eda OMS 24h (25 \u00b5g/m\u00b3)",
    color = "#e74c3c", size = 2.7, hjust = 0, fontface = "italic"
  ) +
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
    title = "Perfil Horario de NO\u2082: Laborables vs Fin de Semana \u2014 Madrid 2025",
    subtitle = paste0(
      "Todas las estaciones de medici\u00f3n  \u00b7  ",
      "Mediana entre estaciones  \u00b7  ",
      "Banda: rango intercuart\u00edlico (Q25\u2013Q75)"
    ),
    x = "Hora del d\u00eda",
    y = "NO\u2082 (\u00b5g/m\u00b3)",
    caption = "Azul: lunes\u2013viernes  \u00b7  Verde: s\u00e1bado y domingo  \u00b7  L\u00ednea roja: gu\u00eda OMS 24h"
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
# 10. RESUMEN
# ==============================================================================

cat("\n==============================================================\n")
cat("RESUMEN\n")
cat("==============================================================\n")
cat(sprintf("  Estaciones de medición incluidas : %d\n", n_est))
cat(sprintf("  Estaciones del año               : Invierno, Primavera, Verano, Otoño\n"))
cat(sprintf(
  "  Registros válidos (laborables)   : %d\n",
  dt_h[tipo_dia == "Laborable" & !is.na(DATO), .N]
))
cat(sprintf(
  "  Registros válidos (finde)        : %d\n",
  dt_h[tipo_dia == "Fin de semana" & !is.na(DATO), .N]
))
cat(sprintf("\nArchivos guardados en:\n  %s\n", dir_salida))
cat(sprintf("  · %s\n", basename(archivo_lab)))
cat(sprintf("  · %s\n", basename(archivo_fds)))
cat(sprintf("  · %s\n", basename(archivo_comp)))
cat("==============================================================\n")
