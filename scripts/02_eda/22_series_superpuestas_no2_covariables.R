# ==============================================================================
# SERIES TEMPORALES SUPERPUESTAS: NO2 vs CADA COVARIABLE  ·  POR ESTACION
# Para 3 estaciones (Plaza Eliptica, El Pardo, Retiro) y 3 escalas temporales
# (horaria, diaria, mensual) se superpone la serie de la variable objetivo
# (log NO2) con la de cada covariable. Ambas se ESTANDARIZAN (z-score) para que
# sean comparables y se aprecie el co-movimiento (relacion temporal).
#
# Una faceta por covariable; dos lineas: NO2 (objetivo) y la covariable.
# Escala horaria: por legibilidad se muestra una VENTANA representativa
#   (por defecto 2 semanas de invierno); el ano completo a nivel horario es
#   ilegible para ver el ciclo diario.
#
# Outputs: outputs/analysis/series_superpuestas/<estacion>/<escala>.png
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

set.seed(4827)

# ------------------------------------------------------------------------------
# PARAMETROS
# ------------------------------------------------------------------------------
ESTACIONES <- c("Plaza Elíptica", "El Pardo", "Retiro")

COL_NO2 <- "#111111"   # variable objetivo (negro)
COL_COV <- "#0057FF"   # covariable (azul vivo)

# Ventana para la escala horaria (para ver el ciclo diario sin saturar)
VENTANA_INI  <- as.Date("2025-01-13")
VENTANA_DIAS <- 14L

# Covariables (columna cruda -> etiqueta)
COVS <- list(
  list(col = "intensidad_raw",          lab = "Intensidad de trafico"),
  list(col = "carga_raw",               lab = "Carga de trafico"),
  list(col = "Temperatura_raw",         lab = "Temperatura (C)"),
  list(col = "Humedad_Relativa_raw",    lab = "Humedad relativa (%)"),
  list(col = "Precipitaciones_raw",     lab = "Precipitaciones (mm)"),
  list(col = "Presion Barométrica_raw", lab = "Presion barometrica (mbar)"),
  list(col = "Radiación Solar_raw",     lab = "Radiacion solar (W/m2)"),
  list(col = "Velocidad Viento_raw",    lab = "Velocidad del viento (m/s)")
)

# Configuracion por escala: dataset + columnas de NO2 + como se arma el tiempo
ESCALAS <- list(
  horaria = list(
    rds     = here("data", "processed", "Maestro", "horario",
                   "dataset_maestro_inla_2025_HORARIO.rds"),
    no2_log = "LOG_NO2_HORARIO",
    tiempo  = function(d) as.POSIXct(d$FECHA, tz = "UTC") + d$HORA * 3600,
    x_lab   = "Fecha y hora"
  ),
  diaria = list(
    rds     = here("data", "processed", "Maestro", "diario",
                   "dataset_maestro_inla_2025_DIARIO.rds"),
    no2_log = "LOG_NO2_DIARIO",
    tiempo  = function(d) as.POSIXct(d$FECHA, tz = "UTC"),
    x_lab   = "Fecha"
  ),
  mensual = list(
    rds     = here("data", "processed", "Maestro", "mensual",
                   "dataset_maestro_inla_2025_MENSUAL.rds"),
    no2_log = "LOG_NO2_MENSUAL",
    tiempo  = function(d) as.POSIXct(d$FECHA, tz = "UTC"),
    x_lab   = "Mes"
  )
)

# ------------------------------------------------------------------------------
# Utilidades
# ------------------------------------------------------------------------------
zscore <- function(x) as.numeric(scale(x))

slug <- function(s) {
  s <- gsub("[^A-Za-z0-9]+", "_", iconv(s, to = "ASCII//TRANSLIT"))
  gsub("^_|_$", "", s)
}

# Construye el data.table largo (una fila por tiempo x covariable x serie)
construir_largo <- function(dt, no2_log, covs_disp) {
  dt <- dt[order(Tiempo)]
  rbindlist(lapply(covs_disp, function(v) {
    ok <- !is.na(dt[[v$col]]) & !is.na(dt[[no2_log]])
    d  <- dt[ok]
    if (nrow(d) < 3) return(NULL)
    rbind(
      data.table(Tiempo = d$Tiempo, Valor = zscore(d[[no2_log]]),
                 Serie = "NO2 (objetivo)", Covariable = v$lab),
      data.table(Tiempo = d$Tiempo, Valor = zscore(d[[v$col]]),
                 Serie = v$lab,             Covariable = v$lab)
    )
  }))
}

crear_grafico <- function(largo, estacion, escala, x_lab, n_obs, sub_extra = "") {
  niveles <- vapply(COVS, function(v) v$lab, character(1))
  largo[, Covariable := factor(Covariable, levels = niveles)]
  # La serie objetivo se dibuja en todas las facetas; su etiqueta de leyenda es fija
  largo[, TipoSerie := ifelse(Serie == "NO2 (objetivo)",
                              "NO2 (objetivo)", "Covariable")]
  ggplot(largo, aes(Tiempo, Valor, color = TipoSerie, linewidth = TipoSerie)) +
    geom_line(alpha = 0.85) +
    facet_wrap(~ Covariable, ncol = 2, scales = "free_y") +
    scale_color_manual(values = c("NO2 (objetivo)" = COL_NO2,
                                  "Covariable"     = COL_COV), name = NULL) +
    scale_linewidth_manual(values = c("NO2 (objetivo)" = 0.55,
                                      "Covariable" = 0.45), guide = "none") +
    labs(
      title = sprintf("Series superpuestas NO2 vs covariables  -  %s", estacion),
      subtitle = sprintf(
        "Escala %s | %s obs | series estandarizadas (z-score)%s",
        escala, format(n_obs, big.mark = "."), sub_extra),
      x = x_lab, y = "Valor estandarizado (z)",
      caption = "Ambas series en z-score para comparar su co-movimiento."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      plot.caption  = element_text(color = "gray35", size = 8),
      strip.text    = element_text(face = "bold", size = 9),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

# ==============================================================================
# BUCLE PRINCIPAL
# ==============================================================================
for (esc_nombre in names(ESCALAS)) {
  cfg <- ESCALAS[[esc_nombre]]
  if (!file.exists(cfg$rds)) { cat(sprintf("[SKIP escala] %s (falta %s)\n",
                                           esc_nombre, cfg$rds)); next }
  dt_all <- as.data.table(readRDS(cfg$rds))
  covs_disp <- Filter(function(v) v$col %in% names(dt_all), COVS)

  for (estacion in ESTACIONES) {
    if (!estacion %in% dt_all$ESTACION) {
      cat(sprintf("[SKIP] %s no existe (%s)\n", estacion, esc_nombre)); next
    }
    dt <- dt_all[ESTACION == estacion]
    dt[, Tiempo := cfg$tiempo(dt)]

    sub_extra <- ""
    if (esc_nombre == "horaria") {
      ini <- as.POSIXct(VENTANA_INI, tz = "UTC")
      fin <- ini + VENTANA_DIAS * 86400
      dt  <- dt[Tiempo >= ini & Tiempo < fin]
      sub_extra <- sprintf(" | ventana %s (+%d dias)",
                           format(VENTANA_INI), VENTANA_DIAS)
    }
    if (nrow(dt) < 3) { cat(sprintf("[SKIP] %s/%s sin datos\n",
                                    estacion, esc_nombre)); next }

    largo <- construir_largo(dt, cfg$no2_log, covs_disp)
    if (is.null(largo) || nrow(largo) == 0) {
      cat(sprintf("[SKIP] %s/%s sin series validas\n", estacion, esc_nombre))
      next
    }
    p <- crear_grafico(largo, estacion, esc_nombre, cfg$x_lab,
                       nrow(dt), sub_extra)

    dir_out <- here("outputs", "analysis", "series_superpuestas", slug(estacion))
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    archivo <- file.path(dir_out, sprintf("series_%s_%s.png",
                                          slug(estacion), esc_nombre))
    ggsave(archivo, p, width = 12, height = 9, dpi = 200, bg = "white")
    cat(sprintf("[OK] %s/%s -> %s\n", estacion, esc_nombre, basename(archivo)))
  }
}

cat("\nListo. Series en outputs/analysis/series_superpuestas/<estacion>/\n")
