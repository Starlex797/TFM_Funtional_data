# ==============================================================================
# PANEL DE SCATTERS: NO2 vs TODAS LAS COVARIABLES  ·  UNA ESTACION POR TIPOLOGIA
# Para cada tipologia de estacion (Suburbana, Urbana fondo, Urbana trafico) se
# toma UNA estacion representativa y se genera un panel con una faceta por
# covariable:
#   - PNG 1: NO2 sin transformar
#   - PNG 2: log(NO2 + 1)
# Puntos en color vivo. Escala horaria, Madrid 2025.
# Outputs: outputs/analysis/scatter_no2_covariables/panel_estacion/<estacion>/
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

set.seed(4827)

# ------------------------------------------------------------------------------
# PARAMETROS
# ------------------------------------------------------------------------------
# Una estacion representativa por tipologia (NOM_TIPO):
ESTACIONES_TIPO <- c(
  "Suburbana"      = "Casa de Campo",
  "Urbana fondo"   = "Retiro",
  "Urbana tráfico" = "Plaza Elíptica"
)
COLOR_VIVO      <- "#0057FF"   # azul vivo
ESCALA          <- "horaria"   # dataset horario

# Covariables (columna cruda -> etiqueta con unidades)
COVS <- list(
  list(col = "intensidad_raw",          lab = "Intensidad de trafico (veh/h)"),
  list(col = "carga_raw",               lab = "Carga de trafico"),
  list(col = "Temperatura_raw",         lab = "Temperatura (C)"),
  list(col = "Humedad_Relativa_raw",    lab = "Humedad relativa (%)"),
  list(col = "Precipitaciones_raw",     lab = "Precipitaciones (mm)"),
  list(col = "Presion Barométrica_raw", lab = "Presion barometrica (mbar)"),
  list(col = "Radiación Solar_raw",     lab = "Radiacion solar (W/m2)"),
  list(col = "Velocidad Viento_raw",    lab = "Velocidad del viento (m/s)")
)

# ==============================================================================
# 1. DATOS
# ==============================================================================
dt_all <- as.data.table(readRDS(here(
  "data", "processed", "Maestro", "horario",
  "dataset_maestro_inla_2025_HORARIO.rds"
)))

# Respuesta cruda y log-transformada
dt_all[, NO2    := DATO]
dt_all[, logNO2 := LOG_NO2_HORARIO]

# --- Eliminar outliers de tráfico (por barrio) --------------------------------
# Son errores del dato de tráfico (valores de miles de veh/h) que comprimen el
# eje y distorsionan la relacion. Se detectan POR BARRIO (el nivel de trafico
# depende del barrio): valor > Tukey (Q3+3*IQR) Y > 2x el percentil 99 del
# barrio. Se anulan intensidad y carga en esas horas.
marcar_outlier_trafico <- function(x) {
  q <- quantile(x, c(.25, .75, .99), na.rm = TRUE)
  x > (q[2] + 3 * (q[2] - q[1])) & x > 2 * q[3]
}
if ("barrio" %in% names(dt_all)) {
  dt_all[, es_out_trafico := FALSE]
  dt_all[!is.na(intensidad_raw),
         es_out_trafico := marcar_outlier_trafico(intensidad_raw), by = barrio]
  cat(sprintf("Outliers de trafico eliminados: %d\n", sum(dt_all$es_out_trafico)))
  dt_all[es_out_trafico == TRUE,
         `:=`(intensidad_raw = NA_real_, carga_raw = NA_real_)]
  dt_all[, es_out_trafico := NULL]
}

covs_disp <- Filter(function(v) v$col %in% names(dt_all), COVS)
niveles   <- vapply(covs_disp, function(v) v$lab, character(1))

# ==============================================================================
# 2. FORMATO LARGO: una fila por (observacion, covariable)
# ==============================================================================
construir_largo <- function(dt, col_y) {
  rbindlist(lapply(covs_disp, function(v) {
    d <- dt[!is.na(get(v$col)) & !is.na(get(col_y)),
            .(X = get(v$col), Y = get(col_y))]
    d[, Covariable := v$lab]
    d
  }))
}

# ==============================================================================
# 3. FUNCION DE GRAFICO (panel facetado)
# ==============================================================================
crear_panel <- function(dt, col_y, y_lab, titulo_y, estacion, tipo) {
  largo <- construir_largo(dt, col_y)
  largo[, Covariable := factor(Covariable, levels = niveles)]

  ggplot(largo, aes(X, Y)) +
    geom_point(color = COLOR_VIVO, size = 0.6, alpha = 0.28) +
    facet_wrap(~ Covariable, scales = "free_x", ncol = 4) +
    labs(
      title = sprintf("NO2 (%s) frente a las covariables  -  %s [%s]",
                      titulo_y, estacion, tipo),
      subtitle = sprintf(
        "Escala %s | %s observaciones | una faceta por covariable",
        ESCALA, format(nrow(dt), big.mark = ".")),
      x = NULL, y = y_lab,
      caption = "Cada punto es una hora. Eje X libre por covariable."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title     = element_text(face = "bold", size = 15),
      plot.subtitle  = element_text(color = "gray40"),
      plot.caption   = element_text(color = "gray35", size = 9),
      strip.text     = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.spacing  = unit(1, "lines")
    )
}

# ==============================================================================
# 4. BUCLE SOBRE UNA ESTACION POR TIPOLOGIA -> 2 PNG cada una
# ==============================================================================
for (tipo in names(ESTACIONES_TIPO)) {
  estacion <- ESTACIONES_TIPO[[tipo]]
  if (!estacion %in% dt_all$ESTACION) {
    cat(sprintf("[SKIP] '%s' no existe en los datos.\n", estacion)); next
  }
  dt <- dt_all[ESTACION == estacion]

  estacion_slug <- gsub("[^A-Za-z0-9]+", "_",
                        iconv(estacion, to = "ASCII//TRANSLIT"))
  estacion_slug <- gsub("^_|_$", "", estacion_slug)

  dir_salida <- here("outputs", "analysis", "scatter_no2_covariables",
                     "panel_estacion", estacion_slug)
  dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n[%s] %s | filas: %d\n", tipo, estacion, nrow(dt)))

  p_raw <- crear_panel(dt, "NO2", "NO2 (ug/m3)", "sin transformar",
                       estacion, tipo)
  f_raw <- file.path(dir_salida, sprintf(
    "panel_%s_NO2_crudo_vs_covariables_2025.png", estacion_slug))
  ggsave(f_raw, p_raw, width = 14, height = 7, dpi = 220, bg = "white")
  cat(sprintf("  [OK] %s\n", basename(f_raw)))

  p_log <- crear_panel(dt, "logNO2", "log(NO2 + 1)", "log-transformado",
                       estacion, tipo)
  f_log <- file.path(dir_salida, sprintf(
    "panel_%s_logNO2_vs_covariables_2025.png", estacion_slug))
  ggsave(f_log, p_log, width = 14, height = 7, dpi = 220, bg = "white")
  cat(sprintf("  [OK] %s\n", basename(f_log)))
}

cat("\nListo. Paneles en outputs/analysis/scatter_no2_covariables/panel_estacion/<estacion>/\n")
