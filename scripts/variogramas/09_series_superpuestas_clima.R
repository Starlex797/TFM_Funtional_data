# ==============================================================================
# SERIES TEMPORALES SUPERPUESTAS — ESTACIONES CLIMATOLÓGICAS MADRID 2025
# Compara visualmente todas las estaciones en cada variable meteorológica
# Outputs: outputs/comparacion_estaciones_clima/
# ==============================================================================
QUÉ HACE:
#
# - Carga las observaciones meteorológicas diarias de 2025.
# - Analiza las siguientes variables:
#     · Temperatura.
#     · Humedad relativa.
#     · Presión barométrica.
#     · Radiación solar.
#     · Velocidad del viento.
#     · Precipitaciones.
#
# Para cada variable meteorológica:
#
# - Elimina las observaciones ausentes.
# - Identifica las estaciones con información disponible.
# - Superpone las series temporales de todas las estaciones.
# - Utiliza un color diferente para cada estación.
# - Añade franjas y separadores mensuales.
# - Genera un archivo PNG independiente.
#
# FINALIDAD PARA EL TFM:
#
# Este script permite evaluar visualmente la similitud de las mediciones entre
# las diferentes estaciones meteorológicas.
#
# Sirve para comprobar:
#
# - Si las estaciones siguen un patrón temporal común.
# - Qué variables presentan mayor variabilidad espacial.
# - Si existen estaciones meteorológicas atípicas.
# - Si aparecen cambios bruscos o posibles errores de medición.
# - Si existen periodos con información incompleta.
# - Si es razonable asignar los datos meteorológicos de una estación climática
#   a una estación de calidad del aire próxima.
#
# Esta comprobación es importante para evaluar la calidad del emparejamiento
# espacial entre los datos de contaminación y las variables meteorológicas.
#
# SALIDAS:
#
# outputs/comparacion_estaciones_clima/
#
# Se genera un PNG para cada variable meteorológica.


library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA DE DATOS
# ==============================================================================

cat("Cargando datos meteorológicos diarios 2025...\n")

dt_meteo_d <- readRDS(here("data", "processed", "clima", "diario",
                            "meteo_madrid_2025_diario.rds"))

cat(sprintf("  Filas: %d  ·  Estaciones: %d\n",
            nrow(dt_meteo_d), uniqueN(dt_meteo_d$ESTACION)))

# ==============================================================================
# 2. DEFINICIÓN DE VARIABLES Y ETIQUETAS
# ==============================================================================

vars_info <- list(
  list(col  = "Temperatura",
       ylab = "Temperatura (°C)",
       tit  = "Temperatura diaria",
       file = "temp"),

  list(col  = "Humedad_Relativa",
       ylab = "Humedad relativa (%)",
       tit  = "Humedad relativa diaria",
       file = "humedad"),

  list(col  = "Presion Barométrica",
       ylab = "Presión barométrica (mbar)",
       tit  = "Presión barométrica diaria",
       file = "presion"),

  list(col  = "Radiación Solar",
       ylab = expression("Radiación solar (W/m"^2*")"),
       tit  = "Radiación solar diaria",
       file = "radiacion"),

  list(col  = "Velocidad Viento",
       ylab = "Velocidad del viento (m/s)",
       tit  = "Velocidad del viento diaria",
       file = "viento"),

  list(col  = "Precipitaciones",
       ylab = "Precipitaciones (mm)",
       tit  = "Precipitaciones diarias",
       file = "precipitacion")
)

# ==============================================================================
# 3. DIRECTORIO DE SALIDA
# ==============================================================================

dir_salida <- here("outputs", "comparacion_estaciones_clima")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Directorio: %s\n\n", dir_salida))

# ==============================================================================
# 4. PARÁMETROS ESTÉTICOS COMPARTIDOS
# ==============================================================================

# Delimitadores mensuales
primeros_meses    <- seq(as.Date("2025-01-01"), as.Date("2025-12-01"), by = "1 month")
etiquetas_meses   <- format(primeros_meses, "%b")

# Franjas alternadas de fondo
franjas <- data.frame(
  xmin  = primeros_meses,
  xmax  = c(primeros_meses[-1], as.Date("2026-01-01")),
  shade = seq_along(primeros_meses) %% 2L == 0L
) |> subset(shade)

tema_sup <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "gray40", size = 9.5),
      plot.caption     = element_text(color = "gray55", size = 8),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 8),
      legend.key.width = unit(1.2, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 9)
    )
}

# ==============================================================================
# 5. BUCLE: UN PNG POR VARIABLE
# ==============================================================================

for (vi in vars_info) {

  col  <- vi$col
  ylab <- vi$ylab
  tit  <- vi$tit
  file <- vi$file

  # Filtrar estaciones con datos en esta variable
  dt_var <- dt_meteo_d[!is.na(get(col)), .(ESTACION, FECHA, VALOR = get(col))]

  n_est <- uniqueN(dt_var$ESTACION)

  if (n_est == 0L) {
    cat(sprintf("  [SKIP] %s — sin datos\n", col))
    next
  }

  # Paleta: tantos colores como estaciones
  paleta <- setNames(
    scales::hue_pal()(n_est),
    sort(unique(dt_var$ESTACION))
  )

  n_col_leyenda <- ifelse(n_est <= 6L, n_est,
                   ifelse(n_est <= 12L, 3L, 4L))

  p <- ggplot(dt_var, aes(x = FECHA, y = VALOR,
                           color = ESTACION, group = ESTACION)) +

    # Franjas de fondo
    geom_rect(
      data        = franjas,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray92", alpha = 0.6
    ) +

    # Líneas divisorias mensuales
    geom_vline(
      xintercept = as.numeric(primeros_meses),
      color = "gray70", linewidth = 0.3, linetype = "dashed"
    ) +

    # Series superpuestas
    geom_line(linewidth = 0.55, alpha = 0.80, na.rm = TRUE) +

    scale_x_date(
      breaks       = primeros_meses,
      labels       = etiquetas_meses,
      minor_breaks = NULL,
      expand       = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    scale_color_manual(values = paleta, name = "Estación") +

    labs(
      title    = paste0(tit, " \u2014 Madrid 2025"),
      subtitle = sprintf("Todas las estaciones superpuestas  \u00b7  %d estaciones con datos",
                         n_est),
      x        = NULL,
      y        = ylab,
      caption  = "Cada l\u00ednea representa una estaci\u00f3n de medici\u00f3n  \u00b7  Datos diarios"
    ) +
    tema_sup() +
    guides(color = guide_legend(ncol = n_col_leyenda,
                                override.aes = list(linewidth = 1.5)))

  # Altura proporcional al nº de filas de leyenda
  n_filas_leyenda <- ceiling(n_est / n_col_leyenda)
  alto <- 5.5 + n_filas_leyenda * 0.35

  archivo <- file.path(dir_salida,
                       sprintf("series_superpuestas_%s_2025.png", file))
  ggsave(archivo, plot = p, width = 13, height = alto, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %-14s  [%2d estaciones]  %s\n",
              col, n_est, basename(archivo)))
}

cat("\n\u2713 Todos los archivos guardados en:\n  ", dir_salida, "\n")
