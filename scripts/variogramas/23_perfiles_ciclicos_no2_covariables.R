# ==============================================================================
# PERFILES CICLICOS: NO2 vs COVARIABLES  ·  POR ESTACION
# No es la serie cronologica, sino el PATRON MEDIO (perfil) en tres ciclos:
#   - Horario : hora del dia (1..24)     -> ciclo diurno
#   - Diario  : dia de la semana (Lun..Dom) -> ciclo semanal
#   - Mensual : mes (Ene..Dic)           -> ciclo estacional/anual
#
# Para cada estacion y ciclo se promedia NO2 y cada covariable por posicion del
# ciclo y se superponen. Dos versiones:
#   (a) ESTANDARIZADA (z-score del perfil): compara la FORMA / co-movimiento.
#   (b) RAW (unidades reales, doble eje): lee valores reales de cada variable.
#
# Estaciones: Plaza Eliptica (trafico), El Pardo (suburbana), Retiro (fondo).
# Outputs: outputs/perfiles_ciclicos/<estacion>/perfil_<ciclo>_<version>.png
# ==============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(here)

ESTACIONES <- c("Plaza Elíptica", "El Pardo", "Retiro")
COL_NO2 <- "#111111"   # variable objetivo (negro)
COL_COV <- "#0057FF"   # covariable (azul vivo)

# A nivel horario NO se promedia: se toma un DIA CONCRETO y se dibujan sus 24 h.
DIA_CONCRETO <- as.Date("2025-01-15")   # miercoles de invierno con 24 h completas

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

# Configuracion por ciclo: dataset, columna NO2 cruda, como se define la
# posicion del ciclo, etiquetas del eje y titulo.
CICLOS <- list(
  horario = list(
    rds   = here("data","processed","Maestro","horario",
                 "dataset_maestro_inla_2025_HORARIO.rds"),
    no2   = "DATO",
    pos   = function(d) as.integer(d$HORA),
    orden = 1:24, etiquetas = as.character(1:24),
    x_lab = "Hora del dia",
    titulo = sprintf("dia %s (%s)", format(DIA_CONCRETO),
                     format(DIA_CONCRETO, "%A")),
    filtro = function(d) d[as.Date(FECHA) == DIA_CONCRETO],  # un dia concreto
    agg = sprintf("valores de un unico dia (%s)", format(DIA_CONCRETO))
  ),
  diario = list(
    rds   = here("data","processed","Maestro","diario",
                 "dataset_maestro_inla_2025_DIARIO.rds"),
    no2   = "DATO_DIARIO",
    pos   = function(d) as.integer(format(d$FECHA, "%u")),   # 1=Lun..7=Dom
    orden = 1:7,
    etiquetas = c("Lun","Mar","Mie","Jue","Vie","Sab","Dom"),
    x_lab = "Dia de la semana", titulo = "ciclo semanal (Lun-Dom)",
    filtro = NULL, agg = "patron medio por posicion del ciclo"
  ),
  mensual = list(
    rds   = here("data","processed","Maestro","mensual",
                 "dataset_maestro_inla_2025_MENSUAL.rds"),
    no2   = "DATO_MENSUAL",
    pos   = function(d) as.integer(format(d$FECHA, "%m")),   # 1..12
    orden = 1:12,
    etiquetas = c("Ene","Feb","Mar","Abr","May","Jun",
                  "Jul","Ago","Sep","Oct","Nov","Dic"),
    x_lab = "Mes", titulo = "ciclo estacional (Ene-Dic)",
    filtro = NULL, agg = "patron medio por posicion del ciclo"
  )
)

slug <- function(s) {
  s <- gsub("[^A-Za-z0-9]+", "_", iconv(s, to = "ASCII//TRANSLIT"))
  gsub("^_|_$", "", s)
}
zscore <- function(x) { s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - mean(x, na.rm = TRUE)); (x - mean(x, na.rm = TRUE)) / s }

# ------------------------------------------------------------------------------
# Perfil medio por posicion del ciclo (NO2 + cada covariable)
# ------------------------------------------------------------------------------
calcular_perfil <- function(dt, no2_col, covs_disp, orden) {
  prof <- data.table(pos = orden)
  prof[, NO2 := sapply(orden, function(g)
    mean(dt[POS == g][[no2_col]], na.rm = TRUE))]
  for (v in covs_disp)
    prof[, (v$lab) := sapply(orden, function(g)
      mean(dt[POS == g][[v$col]], na.rm = TRUE))]
  prof
}

# ------------------------------------------------------------------------------
# (a) Figura ESTANDARIZADA: faceta por covariable, dos lineas z-score
# ------------------------------------------------------------------------------
fig_estandarizada <- function(prof, covs_disp, cfg, estacion, n_obs) {
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
    scale_x_continuous(breaks = cfg$orden, labels = cfg$etiquetas) +
    scale_color_manual(values = c("NO2 (objetivo)" = COL_NO2,
                                  "Covariable" = COL_COV), name = NULL) +
    labs(title = sprintf("NO2 vs covariables  -  %s  -  %s", cfg$titulo, estacion),
         subtitle = sprintf("Estandarizado (z-score) | %s obs | %s",
                            format(n_obs, big.mark = "."), cfg$agg),
         x = cfg$x_lab, y = "Valor estandarizado (z)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "gray40"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "top", panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 7))
}

# ------------------------------------------------------------------------------
# (b) Figura RAW: un mini-grafico de doble eje por covariable (patchwork)
# ------------------------------------------------------------------------------
mini_dual <- function(prof, v, cfg) {
  df <- data.table(pos = prof$pos, no2 = prof$NO2, cov = prof[[v$lab]])
  r_no2 <- range(df$no2, na.rm = TRUE); r_cov <- range(df$cov, na.rm = TRUE)
  d_no2 <- diff(r_no2); d_cov <- diff(r_cov)
  if (!is.finite(d_cov) || d_cov == 0) d_cov <- 1
  if (!is.finite(d_no2) || d_no2 == 0) d_no2 <- 1
  df[, cov_s := (cov - r_cov[1]) / d_cov * d_no2 + r_no2[1]]  # cov -> escala NO2

  ggplot(df, aes(pos)) +
    geom_line(aes(y = no2), color = COL_NO2, linewidth = 0.6) +
    geom_point(aes(y = no2), color = COL_NO2, size = 0.9) +
    geom_line(aes(y = cov_s), color = COL_COV, linewidth = 0.6) +
    geom_point(aes(y = cov_s), color = COL_COV, size = 0.9) +
    scale_x_continuous(breaks = cfg$orden, labels = cfg$etiquetas) +
    scale_y_continuous(
      name = "NO2 (ug/m3)",
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

fig_raw <- function(prof, covs_disp, cfg, estacion, n_obs) {
  minis <- lapply(covs_disp, function(v) mini_dual(prof, v, cfg))
  wrap_plots(minis, ncol = 2) +
    plot_annotation(
      title = sprintf("NO2 vs covariables (valores reales)  -  %s  -  %s",
                      cfg$titulo, estacion),
      subtitle = sprintf("Doble eje: NO2 en negro (izq.), covariable en azul (der.) | %s obs | %s",
                         format(n_obs, big.mark = "."), cfg$agg),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(color = "gray40")))
}

# ==============================================================================
# BUCLE PRINCIPAL
# ==============================================================================
for (ciclo in names(CICLOS)) {
  cfg <- CICLOS[[ciclo]]
  if (!file.exists(cfg$rds)) { cat(sprintf("[SKIP ciclo] %s\n", ciclo)); next }
  dt_all <- as.data.table(readRDS(cfg$rds))
  covs_disp <- Filter(function(v) v$col %in% names(dt_all), COVS)

  for (estacion in ESTACIONES) {
    if (!estacion %in% dt_all$ESTACION) {
      cat(sprintf("[SKIP] %s (%s)\n", estacion, ciclo)); next }
    dt <- dt_all[ESTACION == estacion]
    if (!is.null(cfg$filtro)) dt <- cfg$filtro(dt)
    if (nrow(dt) < 2) { cat(sprintf("[SKIP] %s/%s sin datos\n",
                                    estacion, ciclo)); next }
    dt[, POS := cfg$pos(dt)]
    n_obs <- nrow(dt)

    prof <- calcular_perfil(dt, cfg$no2, covs_disp, cfg$orden)

    dir_out <- here("outputs", "perfiles_ciclicos", slug(estacion))
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

    p_std <- fig_estandarizada(prof, covs_disp, cfg, estacion, n_obs)
    f_std <- file.path(dir_out, sprintf("perfil_%s_estandarizado.png", ciclo))
    ggsave(f_std, p_std, width = 11, height = 9, dpi = 200, bg = "white")

    p_raw <- fig_raw(prof, covs_disp, cfg, estacion, n_obs)
    f_raw <- file.path(dir_out, sprintf("perfil_%s_raw.png", ciclo))
    ggsave(f_raw, p_raw, width = 12, height = 9, dpi = 200, bg = "white")

    cat(sprintf("[OK] %s / %s -> estandarizado + raw\n", estacion, ciclo))
  }
}

cat("\nListo. Perfiles en outputs/perfiles_ciclicos/<estacion>/\n")
