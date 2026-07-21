# ==============================================================================
# SERIES TEMPORALES HISTÓRICAS DE NO₂ — 2019 A 2025
# Madrid · Escala diaria, mensual y horaria · Por estación
# Outputs: outputs/variogramas/series_temporales/
# ==============================================================================
# QUÉ HACE:
#
# - Carga y combina los datos de NO2 correspondientes al periodo 2019–2025.
# - Trabaja con tres resoluciones temporales:
#     · Horaria.
#     · Diaria.
#     · Mensual.
# - Añade el año de procedencia de cada observación.
# - Construye fechas continuas para las series mensuales y horarias.
# - Genera para cada estación:
#     · Una serie histórica diaria.
#     · Una serie histórica mensual.
#     · Una serie histórica horaria.
# - Señala visualmente los comienzos de año y de mes.
# - Incluye una línea horizontal de referencia en 40 µg/m³.
# - Calcula la cobertura temporal disponible para cada estación.
# - Genera paneles conjuntos con todas las estaciones para las tres escalas.
#
# FINALIDAD PARA EL TFM:
#
# Este script proporciona el contexto histórico del NO2 en Madrid. Permite
# estudiar tendencias a largo plazo, estacionalidad, periodos anómalos y
# diferencias persistentes entre estaciones.
#
# También permite comprobar si el comportamiento observado durante 2025 es
# representativo de los años anteriores o presenta características particulares.
#
# Puede utilizarse para analizar:
#
# - Tendencias temporales entre 2019 y 2025.
# - Posibles efectos asociados al periodo de la pandemia.
# - Estacionalidad recurrente del NO2.
# - Años o episodios especialmente contaminados.
# - Diferencias estructurales entre estaciones.
# - Cobertura y continuidad de las series históricas.
#
# SALIDAS:
#
# outputs/historicos_No2/
#     · diario/
#     · mensual/
#     · horario/
#     · panel_diario_todas_estaciones_2019_2025.png
#     · panel_mensual_todas_estaciones_2019_2025.png
#     · panel_horario_todas_estaciones_2019_2025.png


library(data.table)
library(ggplot2)
library(here)

ANIOS <- 2019:2025

# ==============================================================================
# 1. CARGA Y APILADO DE DATOS (7 años)
# ==============================================================================

cargar_anios <- function(escala, patron_archivo, patron_cols) {
  rbindlist(
    lapply(ANIOS, function(a) {
      ruta <- here(
        "data", "processed", "contaminacion", escala,
        sprintf(patron_archivo, a)
      )
      if (!file.exists(ruta)) {
        message("  No encontrado: ", basename(ruta))
        return(NULL)
      }
      dt <- readRDS(ruta)
      dt[, ANO := as.integer(a)]
      dt
    }),
    fill = TRUE
  )
}

cat("Cargando datos diarios...\n")
dt_diario <- cargar_anios("diario", "aire_madrid_%d_No2_trans_diarios.rds")

cat("Cargando datos mensuales...\n")
dt_mensual <- cargar_anios("mensual", "aire_madrid_%d_log_No2_mensuales.rds")

cat("Cargando datos horarios...\n")
dt_horario <- cargar_anios("horario", "aire_madrid_%d_No2_horarios.rds")

# Convertir MES "YYYY-MM" → Date (día 1 de cada mes) para el eje temporal
dt_mensual[, FECHA_MES := as.Date(paste0(MES, "-01"))]

# Construir DATETIME continuo para la serie horaria
dt_horario[, Hora_num := as.integer(sub("H0?", "", as.character(HORA)))]
dt_horario[, DATETIME := as.POSIXct(FECHA) + (Hora_num - 1L) * 3600L]

# Ordenar
setorder(dt_diario, ESTACION, FECHA)
setorder(dt_mensual, ESTACION, FECHA_MES)
setorder(dt_horario, ESTACION, DATETIME)

cat(sprintf(
  "  Diario  : %d obs · %d estaciones\n",
  nrow(dt_diario), uniqueN(dt_diario$ESTACION)
))
cat(sprintf(
  "  Mensual : %d obs · %d estaciones\n",
  nrow(dt_mensual), uniqueN(dt_mensual$ESTACION)
))
cat(sprintf(
  "  Horario : %d obs · %d estaciones\n",
  nrow(dt_horario), uniqueN(dt_horario$ESTACION)
))

# ==============================================================================
# 2. DELIMITADORES TEMPORALES — AÑOS Y MESES
# ==============================================================================

# Primero de cada año (delimitador principal)
primeros_anios <- seq(as.Date("2019-01-01"), as.Date("2026-01-01"), by = "1 year")
primeros_anios_ct <- as.POSIXct(primeros_anios)
etiquetas_anios <- format(primeros_anios, "%Y")

# Primero de cada mes (delimitador secundario) — para diario y mensual
primeros_meses <- seq(as.Date("2019-01-01"), as.Date("2025-12-01"), by = "1 month")
primeros_meses_ct <- as.POSIXct(primeros_meses)

# Franjas alternadas por AÑO (par = gris claro)
anios_limites <- primeros_anios
anios_limites_ct <- primeros_anios_ct

franjas_anios_d <- data.frame(
  xmin  = anios_limites[-length(anios_limites)],
  xmax  = anios_limites[-1],
  shade = (seq_along(ANIOS) %% 2L == 0L)
) |> subset(shade)

franjas_anios_ct <- data.frame(
  xmin  = anios_limites_ct[-length(anios_limites_ct)],
  xmax  = anios_limites_ct[-1],
  shade = (seq_along(ANIOS) %% 2L == 0L)
) |> subset(shade)

# ==============================================================================
# 3. DIRECTORIOS DE SALIDA
# ==============================================================================

dir_base <- here("outputs", "historicos_No2")
dir_diario <- file.path(dir_base, "diario")
dir_mensual <- file.path(dir_base, "mensual")
dir_horario <- file.path(dir_base, "horario")

for (d in c(dir_diario, dir_mensual, dir_horario)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 4. FUNCIÓN AUXILIAR
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", tolower(x))
}

# Tema base compartido
tema_ts <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "gray40", size = 9),
      plot.caption       = element_text(color = "gray55", size = 7.5),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

tema_panel <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = 14),
      plot.subtitle      = element_text(color = "gray40", size = 9),
      strip.text         = element_text(face = "bold", size = 7),
      strip.background   = element_rect(fill = "gray96", color = "gray80"),
      axis.text.x        = element_text(angle = 60, hjust = 1, size = 6),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(0.8, "lines")
    )
}

estaciones_unicas <- sort(unique(dt_diario$ESTACION))
cat(sprintf("\nEstaciones: %d\n\n", length(estaciones_unicas)))

# ==============================================================================
# 5. SERIES DIARIAS — POR ESTACIÓN
# ==============================================================================

cat("--- Generando series DIARIAS por estación (2019-2025) ---\n")

for (est in estaciones_unicas) {
  dt_est <- dt_diario[ESTACION == est & !is.na(DATO_DIARIO)]
  if (nrow(dt_est) < 2L) {
    cat(sprintf("  [SKIP diario] %s\n", est))
    next
  }

  n_dias <- nrow(dt_est)
  pct_ok <- round(100 * n_dias / (365.25 * 7), 1)

  p <- ggplot(dt_est, aes(x = FECHA, y = DATO_DIARIO)) +

    # Franjas alternadas por año
    geom_rect(
      data = franjas_anios_d,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.5
    ) +

    # Líneas divisorias de año (gruesas)
    geom_vline(
      xintercept = as.numeric(primeros_anios),
      color = "gray40", linewidth = 0.6, linetype = "solid"
    ) +

    # Líneas divisorias de mes (finas, secundarias)
    geom_vline(
      xintercept = as.numeric(primeros_meses),
      color = "gray75", linewidth = 0.2, linetype = "dashed"
    ) +

    # Serie temporal
    geom_line(color = "#2c3e50", linewidth = 0.4, alpha = 0.8) +

    # Umbral OMS (40 µg/m³)
    geom_hline(
      yintercept = 40, color = "#e74c3c",
      linewidth = 0.65, linetype = "longdash", alpha = 0.8
    ) +
    annotate("text",
      x = as.Date("2019-02-01"), y = 42,
      label = "Límite OMS 40 µg/m³", color = "#e74c3c",
      size = 2.6, hjust = 0, fontface = "italic"
    ) +

    # Etiquetas de año sobre las líneas divisorias
    annotate("text",
      x = as.Date(paste0(ANIOS, "-07-01")),
      y = Inf,
      label = as.character(ANIOS),
      vjust = 1.5, size = 3.2, color = "gray35", fontface = "bold"
    ) +
    scale_x_date(
      breaks = primeros_meses,
      labels = function(d) {
        ifelse(format(d, "%m") == "01",
          format(d, "%b\n%Y"),
          format(d, "%b")
        )
      },
      minor_breaks = NULL,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    labs(
      title    = paste0("Serie Diaria de NO\u2082 \u2014 ", est),
      subtitle = sprintf("Madrid 2019\u20132025  |  %d d\u00edas v\u00e1lidos (%.1f%%)", n_dias, pct_ok),
      x        = NULL,
      y        = "NO\u2082 (\u00b5g/m\u00b3)",
      caption  = "Franjas grises: a\u00f1os pares  \u00b7  L\u00edneas verticales finas: inicio de mes  \u00b7  L\u00ednea roja: l\u00edmite OMS"
    ) +
    tema_ts()

  nombre_file <- file.path(
    dir_diario,
    paste0("diario_", limpiar_nombre(est), "_2019_2025.png")
  )
  ggsave(nombre_file, plot = p, width = 16, height = 5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 6. SERIES MENSUALES — POR ESTACIÓN
# ==============================================================================

cat("\n--- Generando series MENSUALES por estación (2019-2025) ---\n")

for (est in estaciones_unicas) {
  dt_est <- dt_mensual[ESTACION == est & !is.na(DATO_MENSUAL)]
  if (nrow(dt_est) < 2L) {
    cat(sprintf("  [SKIP mensual] %s\n", est))
    next
  }

  p <- ggplot(dt_est, aes(x = FECHA_MES, y = DATO_MENSUAL)) +

    # Franjas alternadas por año
    geom_rect(
      data = franjas_anios_d,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.5
    ) +

    # Líneas divisorias de año
    geom_vline(
      xintercept = as.numeric(primeros_anios),
      color = "gray40", linewidth = 0.6
    ) +

    # Líneas de mes (secundarias)
    geom_vline(
      xintercept = as.numeric(primeros_meses),
      color = "gray75", linewidth = 0.25, linetype = "dashed"
    ) +

    # Área bajo la curva
    geom_area(fill = "#2980b9", alpha = 0.15) +
    geom_line(color = "#2980b9", linewidth = 0.75) +
    geom_point(color = "#2980b9", size = 1.8, shape = 21, fill = "white", stroke = 1) +

    # Umbral OMS
    geom_hline(
      yintercept = 40, color = "#e74c3c",
      linewidth = 0.65, linetype = "longdash", alpha = 0.8
    ) +
    annotate("text",
      x = as.Date("2019-02-01"), y = 42,
      label = "Límite OMS 40 µg/m³", color = "#e74c3c",
      size = 2.6, hjust = 0, fontface = "italic"
    ) +

    # Etiquetas de año
    annotate("text",
      x = as.Date(paste0(ANIOS, "-07-01")),
      y = Inf,
      label = as.character(ANIOS),
      vjust = 1.5, size = 3.2, color = "gray35", fontface = "bold"
    ) +
    scale_x_date(
      breaks = primeros_meses,
      labels = function(d) {
        ifelse(format(d, "%m") == "01",
          format(d, "%b\n%Y"),
          format(d, "%b")
        )
      },
      minor_breaks = NULL,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      title    = paste0("Serie Mensual de NO\u2082 \u2014 ", est),
      subtitle = "Madrid 2019\u20132025  |  Media mensual de medias diarias v\u00e1lidas",
      x        = NULL,
      y        = "NO\u2082 (\u00b5g/m\u00b3)",
      caption  = "Franjas grises: a\u00f1os pares  \u00b7  L\u00ednea roja: l\u00edmite anual OMS (40 \u00b5g/m\u00b3)"
    ) +
    tema_ts()

  nombre_file <- file.path(
    dir_mensual,
    paste0("mensual_", limpiar_nombre(est), "_2019_2025.png")
  )
  ggsave(nombre_file, plot = p, width = 16, height = 5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 7. SERIES HORARIAS — POR ESTACIÓN
# ==============================================================================

cat("\n--- Generando series HORARIAS por estación (2019-2025) ---\n")

for (est in estaciones_unicas) {
  dt_est <- dt_horario[ESTACION == est & !is.na(DATO)]
  if (nrow(dt_est) < 2L) {
    cat(sprintf("  [SKIP horario] %s\n", est))
    next
  }

  n_horas <- nrow(dt_est)
  pct_ok <- round(100 * n_horas / (8760 * 7), 1)

  p <- ggplot(dt_est, aes(x = DATETIME, y = DATO)) +

    # Franjas alternadas por año
    geom_rect(
      data = franjas_anios_ct,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.5
    ) +

    # Líneas divisorias de año
    geom_vline(
      xintercept = primeros_anios_ct,
      color = "gray40", linewidth = 0.6
    ) +

    # Líneas de mes (secundarias)
    geom_vline(
      xintercept = primeros_meses_ct,
      color = "gray75", linewidth = 0.15, linetype = "dashed"
    ) +

    # Serie horaria
    geom_line(color = "#2c3e50", linewidth = 0.25, alpha = 0.65) +

    # Umbral OMS
    geom_hline(
      yintercept = 40, color = "#e74c3c",
      linewidth = 0.55, linetype = "longdash", alpha = 0.75
    ) +

    # Etiquetas de año
    annotate("text",
      x = as.POSIXct(paste0(ANIOS, "-07-01")),
      y = Inf,
      label = as.character(ANIOS),
      vjust = 1.5, size = 3.2, color = "gray35", fontface = "bold"
    ) +
    scale_x_datetime(
      breaks       = primeros_anios_ct,
      labels       = etiquetas_anios,
      minor_breaks = primeros_meses_ct,
      expand       = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    labs(
      title    = paste0("Serie Horaria de NO\u2082 \u2014 ", est),
      subtitle = sprintf("Madrid 2019\u20132025  |  %d horas v\u00e1lidas (%.1f%% del per\u00edodo)", n_horas, pct_ok),
      x        = NULL,
      y        = "NO\u2082 (\u00b5g/m\u00b3)",
      caption  = "Franjas grises: a\u00f1os pares  \u00b7  Ticks menores: inicio de mes  \u00b7  L\u00ednea roja: l\u00edmite OMS"
    ) +
    tema_ts()

  nombre_file <- file.path(
    dir_horario,
    paste0("horario_", limpiar_nombre(est), "_2019_2025.png")
  )
  ggsave(nombre_file, plot = p, width = 16, height = 5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 8. PANELES RESUMEN — TODAS LAS ESTACIONES (facet)
# ==============================================================================

cat("\n--- Generando paneles resumen (todas las estaciones) ---\n")

# ── Panel diario ──────────────────────────────────────────────────────────────
p_panel_d <- ggplot(
  dt_diario[!is.na(DATO_DIARIO)],
  aes(x = FECHA, y = DATO_DIARIO)
) +
  geom_rect(
    data = franjas_anios_d,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_vline(
    xintercept = as.numeric(primeros_anios),
    color = "gray45", linewidth = 0.4
  ) +
  geom_line(color = "#2c3e50", linewidth = 0.3, alpha = 0.75) +
  geom_hline(
    yintercept = 40, color = "#e74c3c",
    linewidth = 0.4, linetype = "longdash", alpha = 0.7
  ) +
  scale_x_date(
    breaks = primeros_anios, labels = etiquetas_anios,
    minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))
  ) +
  facet_wrap(~ESTACION, ncol = 4, scales = "free_y") +
  labs(
    title = "Series Diarias de NO\u2082 \u2014 Todas las Estaciones",
    subtitle = "Madrid 2019\u20132025  \u00b7  L\u00ednea roja: l\u00edmite OMS (40 \u00b5g/m\u00b3)",
    x = NULL, y = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  tema_panel()

ggsave(file.path(dir_base, "panel_diario_todas_estaciones_2019_2025.png"),
  plot = p_panel_d, width = 22, height = 18, dpi = 200, bg = "white"
)
cat("  \u2713 panel_diario_todas_estaciones_2019_2025.png\n")

# ── Panel mensual ─────────────────────────────────────────────────────────────
p_panel_m <- ggplot(
  dt_mensual[!is.na(DATO_MENSUAL)],
  aes(x = FECHA_MES, y = DATO_MENSUAL)
) +
  geom_rect(
    data = franjas_anios_d,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_vline(
    xintercept = as.numeric(primeros_anios),
    color = "gray45", linewidth = 0.4
  ) +
  geom_area(fill = "#2980b9", alpha = 0.15) +
  geom_line(color = "#2980b9", linewidth = 0.5) +
  geom_point(
    color = "#2980b9", size = 0.9, shape = 21,
    fill = "white", stroke = 0.7
  ) +
  geom_hline(
    yintercept = 40, color = "#e74c3c",
    linewidth = 0.4, linetype = "longdash", alpha = 0.7
  ) +
  scale_x_date(
    breaks = primeros_anios, labels = etiquetas_anios,
    minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))
  ) +
  facet_wrap(~ESTACION, ncol = 4, scales = "free_y") +
  labs(
    title = "Series Mensuales de NO\u2082 \u2014 Todas las Estaciones",
    subtitle = "Madrid 2019\u20132025  \u00b7  L\u00ednea roja: l\u00edmite OMS (40 \u00b5g/m\u00b3)",
    x = NULL, y = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  tema_panel()

ggsave(file.path(dir_base, "panel_mensual_todas_estaciones_2019_2025.png"),
  plot = p_panel_m, width = 22, height = 18, dpi = 200, bg = "white"
)
cat("  \u2713 panel_mensual_todas_estaciones_2019_2025.png\n")

# ── Panel horario ─────────────────────────────────────────────────────────────
p_panel_h <- ggplot(
  dt_horario[!is.na(DATO)],
  aes(x = DATETIME, y = DATO)
) +
  geom_rect(
    data = franjas_anios_ct,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_vline(
    xintercept = primeros_anios_ct,
    color = "gray45", linewidth = 0.35
  ) +
  geom_line(color = "#2c3e50", linewidth = 0.2, alpha = 0.6) +
  geom_hline(
    yintercept = 40, color = "#e74c3c",
    linewidth = 0.35, linetype = "longdash", alpha = 0.7
  ) +
  scale_x_datetime(
    breaks = primeros_anios_ct, labels = etiquetas_anios,
    minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))
  ) +
  facet_wrap(~ESTACION, ncol = 4, scales = "free_y") +
  labs(
    title = "Series Horarias de NO\u2082 \u2014 Todas las Estaciones",
    subtitle = "Madrid 2019\u20132025  \u00b7  L\u00ednea roja: l\u00edmite OMS (40 \u00b5g/m\u00b3)",
    x = NULL, y = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  tema_panel()

ggsave(file.path(dir_base, "panel_horario_todas_estaciones_2019_2025.png"),
  plot = p_panel_h, width = 22, height = 18, dpi = 200, bg = "white"
)
cat("  \u2713 panel_horario_todas_estaciones_2019_2025.png\n")

# ==============================================================================
# 9. RESUMEN FINAL
# ==============================================================================

n_d <- length(list.files(dir_diario, pattern = "\\.png$"))
n_m <- length(list.files(dir_mensual, pattern = "\\.png$"))
n_h <- length(list.files(dir_horario, pattern = "\\.png$"))

cat(sprintf("\n\u2705 Script completado.\n"))
cat(sprintf("   Diario   : %d PNGs en outputs/variogramas/series_temporales/diario/\n", n_d))
cat(sprintf("   Mensual  : %d PNGs en outputs/variogramas/series_temporales/mensual/\n", n_m))
cat(sprintf("   Horario  : %d PNGs en outputs/variogramas/series_temporales/horario/\n", n_h))
cat(sprintf("   Paneles  : 3 PNGs en outputs/variogramas/series_temporales/\n"))
cat(sprintf("   TOTAL    : %d archivos PNG\n", n_d + n_m + n_h + 3L))
