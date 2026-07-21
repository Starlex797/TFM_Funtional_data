# ==============================================================================
# SCATTERPLOTS: NO2 vs COVARIABLES  ·  UNA SOLA ESTACION
# Compara la respuesta cruda y log-transformada frente a cada covariable.
# Al fijar una estacion se elimina la variacion ENTRE estaciones y se ve la
# relacion DENTRO de la estacion (menos difusa).
# Escalas diaria y horaria, Madrid 2025.
# Outputs: outputs/scatter_no2_covariables/no2_vs_log_no2/<estacion>/
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

set.seed(4827)

# ------------------------------------------------------------------------------
# PARAMETRO: estacion a analizar (debe coincidir con la columna ESTACION)
# ------------------------------------------------------------------------------
ESTACION_FILTRO <- "Plaza Elíptica"

# ==============================================================================
# 1. DATOS
# ==============================================================================

dt_d <- as.data.table(readRDS(here(
  "data", "processed", "Maestro", "diario",
  "dataset_maestro_inla_2025_DIARIO.rds"
)))
dt_h <- as.data.table(readRDS(here(
  "data", "processed", "Maestro", "horario",
  "dataset_maestro_inla_2025_HORARIO.rds"
)))

# --- Filtrar a una sola estacion ---
if (!ESTACION_FILTRO %in% dt_h$ESTACION)
  stop(sprintf("La estacion '%s' no existe en los datos.", ESTACION_FILTRO))
dt_d <- dt_d[ESTACION == ESTACION_FILTRO]
dt_h <- dt_h[ESTACION == ESTACION_FILTRO]
cat(sprintf("Estacion: %s | filas diario: %d | filas horario: %d\n",
            ESTACION_FILTRO, nrow(dt_d), nrow(dt_h)))

# Slug seguro para nombres de archivo/carpeta (sin tildes ni caracteres raros)
estacion_slug <- gsub("[^A-Za-z0-9]+", "_",
                      iconv(ESTACION_FILTRO, to = "ASCII//TRANSLIT"))
estacion_slug <- gsub("^_|_$", "", estacion_slug)

resolver_raw <- function(dt, patron) {
  candidatos <- grep(patron, names(dt), value = TRUE)
  candidatos <- candidatos[grepl("_raw$", candidatos)]
  if (length(candidatos) == 0L) return(NA_character_)
  candidatos[1L]
}

vars_cov <- list(
  list(
    col_d = "Temperatura_raw", col_h = "Temperatura_raw",
    label = "Temperatura (C)", file = "temperatura"
  ),
  list(
    col_d = "Humedad_Relativa_raw", col_h = "Humedad_Relativa_raw",
    label = "Humedad relativa (%)", file = "humedad"
  ),
  list(
    col_d = "Precipitaciones_raw", col_h = "Precipitaciones_raw",
    label = "Precipitaciones (mm)", file = "precipitacion"
  ),
  list(
    col_d = resolver_raw(dt_d, "^Presion"),
    col_h = resolver_raw(dt_h, "^Presion"),
    label = "Presion barometrica (mbar)", file = "presion"
  ),
  list(
    col_d = resolver_raw(dt_d, "^Radiaci.*Solar"),
    col_h = resolver_raw(dt_h, "^Radiaci.*Solar"),
    label = "Radiacion solar (W/m2)", file = "radiacion"
  ),
  list(
    col_d = resolver_raw(dt_d, "^Velocidad"),
    col_h = resolver_raw(dt_h, "^Velocidad"),
    label = "Velocidad del viento (m/s)", file = "viento"
  ),
  list(
    col_d = "intensidad_raw", col_h = "intensidad_raw",
    label = "Intensidad de trafico (veh/h)", file = "intensidad"
  ),
  list(
    col_d = "carga_raw", col_h = "carga_raw",
    label = "Carga de trafico", file = "carga"
  )
)

dir_salida <- here(
  "outputs", "scatter_no2_covariables", "no2_vs_log_no2", estacion_slug
)
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. PANELES DE RESPUESTA
# ==============================================================================

crear_paneles_y <- function(dt_plot) {
  rbindlist(
    list(
      data.table(
        X = dt_plot$X,
        Y = dt_plot$NO2,
        Panel = "NO2 sin transformar"
      ),
      data.table(
        X = dt_plot$X,
        Y = dt_plot$logNO2,
        Panel = "log(NO2 + 1)"
      )
    )
  )
}

# ==============================================================================
# 3. GRAFICO
# ==============================================================================

crear_scatter <- function(dt, col_x, col_y_raw, col_y_log, label_x,
                          escala, max_puntos = 40000L) {
  dt_p <- dt[
    !is.na(get(col_x)) & !is.na(get(col_y_raw)) & !is.na(get(col_y_log)),
    .(X = get(col_x), NO2 = get(col_y_raw), logNO2 = get(col_y_log))
  ]
  dt_plot <- if (nrow(dt_p) > max_puntos) {
    dt_p[sample.int(nrow(dt_p), max_puntos)]
  } else {
    copy(dt_p)
  }

  paneles <- crear_paneles_y(dt_plot)

  ggplot(paneles, aes(X, Y)) +
    geom_point(color = "#377EB8", size = 0.35, alpha = 0.16) +
    facet_wrap(~ Panel, nrow = 1, scales = "free_y") +
    labs(
      title = paste("Relacion entre NO2 y", label_x),
      subtitle = sprintf(
        "Estacion %s | escala %s | %s observaciones validas",
        ESTACION_FILTRO, escala, format(nrow(dt_p), big.mark = ".")
      ),
      x = label_x, y = NULL,
      caption = "Comparacion directa de la variable objetivo sin transformar y log-transformada."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      plot.caption = element_text(color = "gray35", size = 9),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.2, "lines")
    )
}

# ==============================================================================
# 4. EJECUCION
# ==============================================================================

archivos <- character()

for (vi in vars_cov) {
  configuraciones <- list(
    list(
      dt = dt_d, col_x = vi$col_d,
      col_y_raw = "DATO_DIARIO", col_y_log = "LOG_NO2_DIARIO",
      escala = "diaria", prefijo = "diario"
    ),
    list(
      dt = dt_h, col_x = vi$col_h,
      col_y_raw = "DATO", col_y_log = "LOG_NO2_HORARIO",
      escala = "horaria", prefijo = "horario"
    )
  )

  for (cfg in configuraciones) {
    if (is.na(cfg$col_x) || !cfg$col_x %in% names(cfg$dt)) {
      cat(sprintf("[SKIP] %s - %s\n", cfg$escala, vi$label))
      next
    }

    p <- crear_scatter(
      dt = cfg$dt, col_x = cfg$col_x,
      col_y_raw = cfg$col_y_raw, col_y_log = cfg$col_y_log,
      label_x = vi$label,
      escala = cfg$escala
    )

    archivo <- file.path(
      dir_salida,
      sprintf("scatter_%s_no2_log_no2_vs_%s_2025.png", cfg$prefijo, vi$file)
    )
    ggsave(archivo, p, width = 13, height = 6, dpi = 220, bg = "white")
    archivos <- c(archivos, archivo)
    cat(sprintf("[OK] %s\n", basename(archivo)))
  }
}

cat(sprintf(
  "\n%d scatterplots guardados en:\n%s\n",
  length(archivos), dir_salida
))
