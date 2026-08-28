# ==============================================================================
# SERIES TEMPORALES: COVARIABLES vs NO₂ — MADRID 2025 (ESCALA DIARIA)
# Compara visualmente cada covariable con el NO₂ bruto por estación
# Todas las estaciones superpuestas · Sin misalignment (filas completas)
# Outputs: outputs/analysis/covariables_vs_no2/
# ==============================================================================
# QUÉ HACE:
#
# - Carga el dataset maestro diario utilizado posteriormente en los modelos INLA.
# - Elimina las filas que no tienen una medición válida de NO2.
# - Esta eliminación evita asociaciones incorrectas entre contaminación,
#   meteorología, tráfico y localización espacial.
#
# El script compara el NO2 con ocho covariables:
#
# Variables meteorológicas:
#
# - Temperatura.
# - Humedad relativa.
# - Precipitaciones.
# - Presión barométrica.
# - Radiación solar.
# - Velocidad del viento.
#
# Variables de tráfico:
#
# - Intensidad del tráfico.
# - Carga de tráfico.
#
# Para cada covariable:
#
# - Conserva únicamente las filas donde existen tanto el NO2 como la covariable.
# - Calcula el porcentaje de cobertura conjunta.
# - Transforma los datos a formato largo.
# - Genera un gráfico con dos paneles:
#     · Panel superior: serie diaria de NO2.
#     · Panel inferior: serie diaria de la covariable.
# - Superpone las series de todas las estaciones.
# - Utiliza escalas verticales independientes, ya que las variables tienen
#   unidades y rangos diferentes.
#
# FINALIDAD PARA EL TFM:
#
# Este script realiza una primera exploración temporal de la relación entre el
# NO2 y sus posibles variables explicativas.
#
# Permite detectar:
#
# - Ciclos estacionales comunes.
# - Relaciones temporales directas o inversas.
# - Episodios meteorológicos asociados con reducciones del NO2.
# - Aumentos de contaminación coincidentes con una mayor intensidad de tráfico.
# - Posibles retrasos entre una covariable y la respuesta del NO2.
# - Periodos con baja cobertura conjunta.
# - Diferencias en estos patrones entre estaciones.
#
# Este análisis es un paso previo a:
#
# - Los análisis de correlación.
# - La selección de covariables.
# - La formulación de los modelos estadísticos.
# - La interpretación de los coeficientes estimados posteriormente.
#
# SALIDAS:
#
# outputs/analysis/covariables_vs_no2/
library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA Y LIMPIEZA DEL DATASET MAESTRO
# ==============================================================================

cat("Cargando dataset maestro diario 2025...\n")

dt_m <- readRDS(here("data", "processed", "Maestro", "diario",
                     "dataset_maestro_inla_2025_DIARIO.rds"))

# Eliminación del misalignment: conservar solo filas con NO₂ válido
# (filas sin DATO_DIARIO implican que la estación no tiene medición ese día,
#  lo que produce una asociación incorrecta con el barrio de tráfico)
dt_m <- dt_m[!is.na(DATO_DIARIO)]

cat(sprintf("  Filas tras eliminar misalignment: %d\n", nrow(dt_m)))
cat(sprintf("  Estaciones                       : %d\n", uniqueN(dt_m$ESTACION)))

# ==============================================================================
# 2. DEFINICIÓN DE COVARIABLES
# ==============================================================================

# Resolver nombre exacto de columnas con caracteres especiales
col_presion_raw <- grep("^Presion", names(dt_m), value = TRUE)
col_presion_raw <- col_presion_raw[grepl("_raw", col_presion_raw)]

vars_cov <- list(

  # Climatológicas
  list(col   = "Temperatura_raw",
       label = "Temperatura (\u00b0C)",
       file  = "temperatura"),

  list(col   = "Humedad_Relativa_raw",
       label = "Humedad relativa (%)",
       file  = "humedad"),

  list(col   = "Precipitaciones_raw",
       label = "Precipitaciones (mm)",
       file  = "precipitacion"),

  list(col   = col_presion_raw,
       label = "Presi\u00f3n barom\u00e9trica (mbar)",
       file  = "presion"),

  list(col   = "Radiaci\u00f3n Solar_raw",
       label = "Radiaci\u00f3n solar (W/m\u00b2)",
       file  = "radiacion"),

  list(col   = "Velocidad Viento_raw",
       label = "Velocidad del viento (m/s)",
       file  = "viento"),

  # Tráfico
  list(col   = "intensidad_raw",
       label = "Intensidad tr\u00e1fico (veh/h)",
       file  = "intensidad"),

  list(col   = "carga_raw",
       label = "Carga tr\u00e1fico",
       file  = "carga")
)

# ==============================================================================
# 3. PARÁMETROS ESTÉTICOS COMPARTIDOS
# ==============================================================================

primeros_meses  <- seq(as.Date("2025-01-01"), as.Date("2025-12-01"), by = "1 month")
etiquetas_meses <- format(primeros_meses, "%b")

franjas <- data.frame(
  xmin  = primeros_meses,
  xmax  = c(primeros_meses[-1], as.Date("2026-01-01")),
  shade = seq_along(primeros_meses) %% 2L == 0L
) |> subset(shade)

n_est    <- uniqueN(dt_m$ESTACION)
estaciones_ord <- sort(unique(dt_m$ESTACION))
paleta   <- setNames(scales::hue_pal()(n_est), estaciones_ord)
n_col_lg <- ifelse(n_est <= 8L, 2L, 4L)

tema_cov <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title        = element_text(face = "bold", size = 13),
      plot.subtitle     = element_text(color = "gray40", size = 9.5),
      plot.caption      = element_text(color = "gray55", size = 8),
      strip.text        = element_text(face = "bold", size = 10.5,
                                       margin = margin(t = 5, b = 5)),
      strip.background  = element_rect(fill = "gray96", color = "gray80"),
      legend.position   = "bottom",
      legend.title      = element_text(face = "bold", size = 9),
      legend.text       = element_text(size = 7.5),
      legend.key.width  = unit(1.1, "cm"),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank(),
      axis.text.x         = element_text(angle = 45, hjust = 1, size = 8.5),
      panel.spacing       = unit(1, "lines")
    )
}

# ==============================================================================
# 4. DIRECTORIO DE SALIDA
# ==============================================================================

dir_salida <- here("outputs", "analysis", "covariables_vs_no2")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nDirectorio de salida: %s\n\n", dir_salida))

# ==============================================================================
# 5. BUCLE PRINCIPAL — UN PNG POR COVARIABLE
# ==============================================================================

for (vi in vars_cov) {

  col_cov <- vi$col

  # Comprobar que la columna existe
  if (!col_cov %in% names(dt_m)) {
    cat(sprintf("  [SKIP] %-35s — columna no encontrada\n", col_cov))
    next
  }

  # Seleccionar y filtrar: solo filas con ambas variables presentes
  dt_plot <- dt_m[
    !is.na(get(col_cov)),
    .(ESTACION, FECHA,
      NO2 = DATO_DIARIO,
      COV = get(col_cov))
  ]

  n_validas <- nrow(dt_plot)
  if (n_validas < 2L) {
    cat(sprintf("  [SKIP] %-35s — sin datos suficientes\n", col_cov))
    next
  }

  # Pivotar a formato largo para facet_grid con escalas libres
  dt_long <- melt(
    dt_plot,
    id.vars      = c("ESTACION", "FECHA"),
    measure.vars = c("NO2", "COV"),
    variable.name = "panel_id",
    value.name    = "valor"
  )

  # Etiquetas de panel (NO₂ siempre arriba)
  label_no2 <- "NO\u2082 (\u00b5g/m\u00b3)"
  dt_long[, panel := fifelse(panel_id == "NO2", label_no2, vi$label)]
  dt_long[, panel := factor(panel, levels = c(label_no2, vi$label))]

  # Porcentaje de cobertura de la covariable
  pct_cov <- round(100 * n_validas / nrow(dt_m), 1)

  p <- ggplot(dt_long, aes(x = FECHA, y = valor,
                            color = ESTACION, group = ESTACION)) +

    # Franjas de fondo
    geom_rect(
      data        = franjas,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray92", alpha = 0.55
    ) +

    # Líneas divisorias mensuales
    geom_vline(
      xintercept  = as.numeric(primeros_meses),
      color       = "gray68", linewidth = 0.3, linetype = "dashed"
    ) +

    # Series superpuestas por estación
    geom_line(linewidth = 0.5, alpha = 0.80, na.rm = TRUE) +

    # Paneles independientes con escala Y libre
    facet_grid(panel ~ ., scales = "free_y", switch = "y") +

    scale_x_date(
      breaks       = primeros_meses,
      labels       = etiquetas_meses,
      minor_breaks = NULL,
      expand       = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.06))) +
    scale_color_manual(values = paleta, name = "Estaci\u00f3n") +

    labs(
      title    = sprintf("NO\u2082 vs %s \u2014 Madrid 2025", vi$label),
      subtitle = sprintf(
        "Series diarias superpuestas  \u00b7  %d estaciones  \u00b7  %.1f%% filas v\u00e1lidas en ambas variables",
        n_est, pct_cov
      ),
      x       = NULL,
      y       = NULL,
      caption = "Panel superior: NO\u2082 bruto (\u00b5g/m\u00b3)  \u00b7  Panel inferior: covariable en unidades originales  \u00b7  Cada l\u00ednea: una estaci\u00f3n de medici\u00f3n"
    ) +
    tema_cov() +
    guides(color = guide_legend(ncol = n_col_lg,
                                override.aes = list(linewidth = 1.5)))

  archivo <- file.path(dir_salida,
                       sprintf("series_no2_vs_%s_2025.png", vi$file))
  ggsave(archivo, plot = p, width = 13, height = 8, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %-35s  [%d filas]  %s\n",
              col_cov, n_validas, basename(archivo)))
}

cat(sprintf("\n\u2713 Todos los archivos guardados en:\n  %s\n", dir_salida))
