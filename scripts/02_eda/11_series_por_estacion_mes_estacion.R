# ==============================================================================
# SERIES TEMPORALES DE NO₂ POR ESTACIÓN SELECCIONADA — MADRID 2025
# Vista mensual (12 paneles) y por estación del año (4 paneles)
# Outputs: outputs/analysis/series_por_estacion/mensual/ y /estacion_anio/
# ==============================================================================
# QUÉ HACE:
#
# - Carga el dataset maestro diario de 2025.
# - Selecciona siete estaciones de interés:
#     · El Pardo.
#     · Casa de Campo.
#     · Paseo de la Castellana.
#     · Plaza Elíptica.
#     · Moratalaz.
#     · Arturo Soria.
#     · Vallecas.
#
# - Elimina las filas sin una concentración diaria válida de NO2.
# - Construye variables para:
#     · Mes.
#     · Día del mes.
#     · Estación del año.
#
# - Genera dos gráficos para cada estación seleccionada.
#
# GRÁFICO MENSUAL:
#
# - Divide la serie diaria en doce paneles, uno por mes.
# - Representa la evolución del NO2 a lo largo de los días de cada mes.
# - Añade una línea horizontal con la mediana mensual.
# - Añade una referencia horizontal en 40 µg/m³.
#
# GRÁFICO ESTACIONAL:
#
# - Divide la serie en cuatro paneles:
#     · Invierno.
#     · Primavera.
#     · Verano.
#     · Otoño.
# - Representa la evolución diaria del NO2 dentro de cada estación del año.
# - Añade una referencia horizontal en 40 µg/m³.
#
# FINALIDAD PARA EL TFM:
#
# Este script permite examinar detalladamente la evolución diaria del NO2 en un
# conjunto representativo de estaciones.
#
# Las estaciones elegidas incluyen entornos de características diferentes, lo
# que permite estudiar:
#
# - Diferencias entre estaciones urbanas y suburbanas.
# - Episodios concretos de contaminación.
# - Persistencia de concentraciones elevadas.
# - Variabilidad dentro de cada mes.
# - Cambios entre estaciones del año.
# - Diferencias estacionales según la localización.
#
# Complementa los paneles generales porque permite observar con mayor claridad
# el comportamiento individual de algunas estaciones relevantes para el TFM.
#
# OBSERVACIÓN SOBRE EL SCRIPT:
#
# En el gráfico mensual se dibuja una línea azul correspondiente a la mediana
# mensual. Sin embargo, el texto del gráfico menciona una banda de rango
# intercuartílico. Actualmente no se genera esa banda, por lo que sería
# conveniente cambiar el texto o añadir realmente el rango intercuartílico.
#
# SALIDAS:
#
# outputs/analysis/series_por_estacion/
#     · mensual/
#     · estacion_anio/
#
# Se generan dos PNG por estación seleccionada.


library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA Y FILTRADO
# ==============================================================================

cat("Cargando dataset maestro diario 2025...\n")

dt_m <- readRDS(here("data", "processed", "Maestro", "diario",
                     "dataset_maestro_inla_2025_DIARIO.rds"))

# Estaciones de interés
estaciones_sel <- c(
  "El Pardo",
  "Casa de Campo",
  "P\u00ba Castellana",
  "Plaza El\u00edptica",
  "Moratalaz",
  "Arturo Soria",
  "Vallecas"
)

# Filtrar y eliminar filas sin NO₂ (misalignment)
dt_sel <- dt_m[ESTACION %in% estaciones_sel & !is.na(DATO_DIARIO)]

# Respetar el orden de presentación elegido
dt_sel[, ESTACION := factor(ESTACION, levels = estaciones_sel)]

cat(sprintf("  Estaciones cargadas : %d\n", uniqueN(dt_sel$ESTACION)))
cat(sprintf("  Filas válidas       : %d\n", nrow(dt_sel)))

# ==============================================================================
# 2. VARIABLES TEMPORALES AUXILIARES
# ==============================================================================

# Mes
nombres_mes <- c("Enero","Febrero","Marzo","Abril","Mayo","Junio",
                 "Julio","Agosto","Septiembre","Octubre","Noviembre","Diciembre")
dt_sel[, mes_num := as.integer(format(FECHA, "%m"))]
dt_sel[, mes_nom := factor(nombres_mes[mes_num], levels = nombres_mes)]
dt_sel[, dia_mes := as.integer(format(FECHA, "%d"))]

# Estación del año
dt_sel[, estacion_anio := fcase(
  mes_num %in% c(12L, 1L, 2L),  "Invierno",
  mes_num %in% c(3L, 4L, 5L),   "Primavera",
  mes_num %in% c(6L, 7L, 8L),   "Verano",
  mes_num %in% c(9L, 10L, 11L), "Oto\u00f1o"
)]
dt_sel[, estacion_anio := factor(
  estacion_anio,
  levels = c("Invierno", "Primavera", "Verano", "Oto\u00f1o")
)]

# ==============================================================================
# 3. PALETAS Y TEMAS
# ==============================================================================

paleta_meses <- setNames(
  colorRampPalette(c("#1a6faf","#2ecc71","#e67e22","#c0392b"))(12L),
  nombres_mes
)

paleta_estaciones <- c(
  "Invierno"  = "#2980b9",
  "Primavera" = "#27ae60",
  "Verano"    = "#e67e22",
  "Oto\u00f1o" = "#8e44ad"
)

tema_serie <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold", size = 14),
      plot.subtitle     = element_text(color = "gray40", size = 9.5),
      plot.caption      = element_text(color = "gray55", size = 8),
      strip.text        = element_text(face = "bold", size = 10,
                                       margin = margin(t = 5, b = 5)),
      strip.background  = element_rect(fill = "gray96", color = "gray80"),
      legend.position   = "none",
      panel.grid.minor  = element_blank(),
      panel.spacing     = unit(0.8, "lines"),
      axis.title        = element_text(size = 10),
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 8)
    )
}

# ==============================================================================
# 4. DIRECTORIOS DE SALIDA
# ==============================================================================

dir_mensual    <- here("outputs", "analysis", "series_por_estacion", "mensual")
dir_estac_anio <- here("outputs", "analysis", "series_por_estacion", "estacion_anio")
dir.create(dir_mensual,    recursive = TRUE, showWarnings = FALSE)
dir.create(dir_estac_anio, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nCarpeta de salida: %s\n\n",
            here("outputs", "analysis", "series_por_estacion")))

# ==============================================================================
# 5. FUNCIÓN AUXILIAR: nombre de archivo limpio
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  tolower(gsub("^_|_$", "", x))
}

# ==============================================================================
# 6. BUCLE — UNA ESTACIÓN A LA VEZ
# ==============================================================================

for (est in estaciones_sel) {

  dt_est <- dt_sel[ESTACION == est]

  if (nrow(dt_est) < 5L) {
    cat(sprintf("  [SKIP] %s — sin datos\n", est)); next
  }

  n_dias <- nrow(dt_est)
  nombre_archivo <- limpiar_nombre(est)

  # ------------------------------------------------------------------
  # 6a. GRÁFICO MENSUAL — 12 paneles (uno por mes)
  # ------------------------------------------------------------------

  # Mediana mensual (línea horizontal de referencia dentro de cada panel)
  mediana_mes <- dt_est[, .(mediana = median(DATO_DIARIO, na.rm = TRUE)), by = mes_nom]

  p_mes <- ggplot(dt_est, aes(x = dia_mes, y = DATO_DIARIO)) +

    # Línea de mediana mensual (referencia horizontal por panel)
    geom_hline(
      data        = mediana_mes,
      aes(yintercept = mediana),
      color       = "steelblue", linewidth = 0.5,
      linetype    = "dotted", alpha = 0.8,
      inherit.aes = FALSE
    ) +

    # Serie diaria coloreada por mes
    geom_line(aes(color = mes_nom), linewidth = 0.75, na.rm = TRUE) +
    geom_point(aes(color = mes_nom), size = 1.4, alpha = 0.75, na.rm = TRUE) +

    # Límite OMS 40 µg/m³ (media anual de referencia visual)
    geom_hline(yintercept = 40, color = "#e74c3c",
               linewidth = 0.55, linetype = "longdash", alpha = 0.75) +

    facet_wrap(~ mes_nom, nrow = 3, ncol = 4) +

    scale_x_continuous(
      breaks = c(1, 10, 20, 31),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    scale_color_manual(values = paleta_meses) +

    labs(
      title    = sprintf("NO\u2082 diario por mes \u2014 %s", est),
      subtitle = sprintf(
        "Madrid 2025  \u00b7  %d d\u00edas v\u00e1lidos  \u00b7  Cada panel: un mes del a\u00f1o",
        n_dias
      ),
      x       = "D\u00eda del mes",
      y       = "NO\u2082 (\u00b5g/m\u00b3)",
      caption = "L\u00ednea roja discontinua: l\u00edmite OMS 40 \u00b5g/m\u00b3  \u00b7  Banda azul: rango IQR del mes"
    ) +
    tema_serie()

  archivo_mes <- file.path(dir_mensual,
                           sprintf("mensual_no2_%s_2025.png", nombre_archivo))
  ggsave(archivo_mes, plot = p_mes,
         width = 14, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  [mensual]      \u2713 %s\n", basename(archivo_mes)))

  # ------------------------------------------------------------------
  # 6b. GRÁFICO POR ESTACIÓN DEL AÑO — 4 paneles
  # ------------------------------------------------------------------

  # Numeración de día dentro de cada estación (evita el gap en Invierno)
  dt_est[, dia_estacion := rowid(estacion_anio)]

  p_est <- ggplot(dt_est, aes(x = dia_estacion, y = DATO_DIARIO,
                               color = estacion_anio)) +

    # Serie diaria
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_point(size = 0.9, alpha = 0.65, na.rm = TRUE) +

    # Límite OMS
    geom_hline(yintercept = 40, color = "#e74c3c",
               linewidth = 0.55, linetype = "longdash", alpha = 0.75) +
    annotate("text", x = 2, y = 42,
             label = "OMS 40 \u00b5g/m\u00b3", color = "#e74c3c",
             size = 2.8, hjust = 0, fontface = "italic") +

    facet_wrap(~ estacion_anio, nrow = 2, ncol = 2) +

    scale_x_continuous(
      breaks = c(1, 30, 60, 90),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    scale_color_manual(values = paleta_estaciones) +

    labs(
      title    = sprintf("NO\u2082 diario por estaci\u00f3n del a\u00f1o \u2014 %s", est),
      subtitle = sprintf(
        "Madrid 2025  \u00b7  %d d\u00edas v\u00e1lidos  \u00b7  Cada panel: una estaci\u00f3n del a\u00f1o",
        n_dias
      ),
      x       = "D\u00eda dentro de la estaci\u00f3n",
      y       = "NO\u2082 (\u00b5g/m\u00b3)",
      caption = "L\u00ednea roja discontinua: l\u00edmite OMS 40 \u00b5g/m\u00b3  \u00b7  Invierno: ene\u2013feb + dic"
    ) +
    tema_serie()

  archivo_est <- file.path(dir_estac_anio,
                            sprintf("estacion_anio_no2_%s_2025.png", nombre_archivo))
  ggsave(archivo_est, plot = p_est,
         width = 13, height = 9, dpi = 200, bg = "white")
  cat(sprintf("  [estacion_año] \u2713 %s\n\n", basename(archivo_est)))
}

# ==============================================================================
# 7. RESUMEN
# ==============================================================================

cat("==============================================================\n")
cat(sprintf("  Estaciones procesadas : %d\n", length(estaciones_sel)))
cat(sprintf("  PNG mensuales        : %d  →  %s\n",
            length(estaciones_sel), dir_mensual))
cat(sprintf("  PNG estación del año : %d  →  %s\n",
            length(estaciones_sel), dir_estac_anio))
cat("==============================================================\n")
