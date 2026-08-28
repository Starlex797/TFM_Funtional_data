# ==============================================================================
# DIAS-EVENTO A NIVEL HORARIO: ¿COMO AFECTAN LLUVIA Y VIENTO AL NO2?
# Se toman dos dias concretos (uno lluvioso, uno ventoso) y se dibujan sus 24 h
# superponiendo NO2 con cada covariable. Graficos APARTE por evento y estacion.
# Mismo formato que los perfiles ciclicos: version estandarizada (z-score) y
# version raw (doble eje, unidades reales).
#
# Estaciones: Plaza Eliptica (trafico), El Pardo (suburbana), Retiro (fondo).
# Outputs: outputs/analysis/dias_evento/<estacion>/<evento>_<fecha>_{estandarizado,raw}.png
# ==============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(here)

ESTACIONES <- c("Plaza Elíptica", "El Pardo", "Retiro")
COL_NO2 <- "#111111"
COL_COV <- "#0057FF"

# Dias-evento seleccionados (24 h completas en las 3 estaciones):
#   - lluvia: 2025-05-02 (max lluvia media, viento bajo -> aisla la lluvia)
#   - viento: 2025-01-27 (max viento, casi sin lluvia -> aisla el viento)
DIAS_EVENTO <- list(
  lluvia = as.Date("2025-05-02"),
  viento = as.Date("2025-01-27")
)

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

ORDEN <- 1:24
ETIQ  <- as.character(1:24)
NO2_COL <- "DATO"

slug <- function(s) {
  s <- gsub("[^A-Za-z0-9]+", "_", iconv(s, to = "ASCII//TRANSLIT"))
  gsub("^_|_$", "", s)
}
zscore <- function(x) { s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - mean(x, na.rm = TRUE))
  (x - mean(x, na.rm = TRUE)) / s }

# ------------------------------------------------------------------------------
# Perfil (24 h de un dia): NO2 + cada covariable por hora
# ------------------------------------------------------------------------------
perfil_dia <- function(dt, covs_disp) {
  prof <- data.table(pos = ORDEN)
  prof[, NO2 := sapply(ORDEN, function(g) {
    v <- dt[HORA == g][[NO2_COL]]; if (length(v)) v[1] else NA_real_ })]
  for (v in covs_disp)
    prof[, (v$lab) := sapply(ORDEN, function(g) {
      x <- dt[HORA == g][[v$col]]; if (length(x)) x[1] else NA_real_ })]
  prof
}

# (a) Version estandarizada
fig_std <- function(prof, covs_disp, estacion, titulo, agg) {
  niveles <- vapply(covs_disp, function(v) v$lab, character(1))
  largo <- rbindlist(lapply(covs_disp, function(v) rbind(
    data.table(pos = prof$pos, Valor = zscore(prof$NO2),
               Serie = "NO2 (objetivo)", Covariable = v$lab),
    data.table(pos = prof$pos, Valor = zscore(prof[[v$lab]]),
               Serie = "Covariable",     Covariable = v$lab))))
  largo[, Covariable := factor(Covariable, levels = niveles)]
  ggplot(largo, aes(pos, Valor, color = Serie)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
    facet_wrap(~ Covariable, ncol = 2) +
    scale_x_continuous(breaks = ORDEN, labels = ETIQ) +
    scale_color_manual(values = c("NO2 (objetivo)" = COL_NO2,
                                  "Covariable" = COL_COV), name = NULL) +
    labs(title = sprintf("NO2 vs covariables  -  %s  -  %s", titulo, estacion),
         subtitle = sprintf("Estandarizado (z-score) | %s", agg),
         x = "Hora del dia", y = "Valor estandarizado (z)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "gray40"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "top", panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 7))
}

# (b) Version raw (doble eje por covariable)
mini_dual <- function(prof, v) {
  df <- data.table(pos = prof$pos, no2 = prof$NO2, cov = prof[[v$lab]])
  r_no2 <- range(df$no2, na.rm = TRUE); r_cov <- range(df$cov, na.rm = TRUE)
  d_no2 <- diff(r_no2); d_cov <- diff(r_cov)
  if (!is.finite(d_cov) || d_cov == 0) d_cov <- 1
  if (!is.finite(d_no2) || d_no2 == 0) d_no2 <- 1
  df[, cov_s := (cov - r_cov[1]) / d_cov * d_no2 + r_no2[1]]
  ggplot(df, aes(pos)) +
    geom_line(aes(y = no2), color = COL_NO2, linewidth = 0.6) +
    geom_point(aes(y = no2), color = COL_NO2, size = 0.9) +
    geom_line(aes(y = cov_s), color = COL_COV, linewidth = 0.6) +
    geom_point(aes(y = cov_s), color = COL_COV, size = 0.9) +
    scale_x_continuous(breaks = ORDEN, labels = ETIQ) +
    scale_y_continuous(name = "NO2 (ug/m3)",
      sec.axis = sec_axis(~ (. - r_no2[1]) / d_no2 * d_cov + r_cov[1],
                          name = v$lab)) +
    labs(title = v$lab, x = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9),
          axis.title.y.left  = element_text(color = COL_NO2, size = 8),
          axis.title.y.right = element_text(color = COL_COV, size = 8),
          axis.text = element_text(size = 6.5),
          panel.grid.minor = element_blank())
}

fig_raw <- function(prof, covs_disp, estacion, titulo, agg) {
  minis <- lapply(covs_disp, function(v) mini_dual(prof, v))
  wrap_plots(minis, ncol = 2) +
    plot_annotation(
      title = sprintf("NO2 vs covariables (valores reales)  -  %s  -  %s",
                      titulo, estacion),
      subtitle = sprintf("Doble eje: NO2 negro (izq.), covariable azul (der.) | %s", agg),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(color = "gray40")))
}

# ==============================================================================
# BUCLE PRINCIPAL
# ==============================================================================
dt_all <- as.data.table(readRDS(here("data","processed","Maestro","horario",
                                     "dataset_maestro_inla_2025_HORARIO.rds")))
dt_all[, DIA := as.Date(FECHA)]
covs_disp <- Filter(function(v) v$col %in% names(dt_all), COVS)

for (evento in names(DIAS_EVENTO)) {
  fecha <- DIAS_EVENTO[[evento]]
  titulo <- sprintf("dia %s %s (%s)", evento, format(fecha),
                    format(fecha, "%A"))
  agg <- sprintf("dia %s: %s", evento, format(fecha))

  for (estacion in ESTACIONES) {
    dt <- dt_all[ESTACION == estacion & DIA == fecha]
    if (nrow(dt) < 2) { cat(sprintf("[SKIP] %s/%s sin datos\n",
                                    estacion, evento)); next }
    prof <- perfil_dia(dt, covs_disp)

    dir_out <- here("outputs", "analysis", "dias_evento", slug(estacion))
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    base <- sprintf("%s_%s", evento, format(fecha))

    ggsave(file.path(dir_out, sprintf("%s_estandarizado.png", base)),
           fig_std(prof, covs_disp, estacion, titulo, agg),
           width = 11, height = 9, dpi = 200, bg = "white")
    ggsave(file.path(dir_out, sprintf("%s_raw.png", base)),
           fig_raw(prof, covs_disp, estacion, titulo, agg),
           width = 12, height = 9, dpi = 200, bg = "white")
    cat(sprintf("[OK] %s / %s -> estandarizado + raw\n", estacion, evento))
  }
}

cat("\nListo. Dias-evento en outputs/analysis/dias_evento/<estacion>/\n")
