# ==============================================================================
# SCRIPT: EXPLORACIÓN VISUAL COMPLETA ANTES DE MODELADO
# ==============================================================================
# Objetivo:
#   - Cargar datasets (horario, diario, mensual)
#   - Visualizar: distribuciones, correlaciones (NO2 raw), linealidad
#   - Crear un reporte PNG consolidado por frecuencia
#   - Detectar anomalías, valores faltantes, outliers antes de INLA
# ==============================================================================

library(INLA)
library(fmesher)
library(data.table)
library(sf)
library(here)
library(ggplot2)
library(car)
library(mgcv)
library(patchwork)
library(gt)

set.seed(4827)

# ==============================================================================
# CONFIG
# ==============================================================================

VIF_UMBRAL     <- 5
DIC_MEJORA_MIN <- 2

# Covariables estandarizadas y raw
COVS_NOMBRES <- c("intensidad", "Temperatura",
                  "Precipitaciones", "Presion Barométrica",
                  "Velocidad Viento")
COVS_ALIAS   <- c("trafico_intensidad", "temperatura",
                  "precipitacion", "presion_barometrica",
                  "velocidad_viento")

# ==============================================================================
# CARGAR DATOS
# ==============================================================================

cat("\n", strrep("=", 80), "\n")
cat("  CARGA DE DATOS\n")
cat(strrep("=", 80), "\n")

dt_h <- readRDS(here("data", "processed", "Maestro", "horario", "dataset_maestro_inla_2025_HORARIO.rds"))
dt_d <- readRDS(here("data", "processed", "Maestro", "diario", "dataset_maestro_inla_2025_DIARIO.rds"))

setDT(dt_h)
setDT(dt_d)

setnames(dt_h, "LOG_NO2_HORARIO", "LOG_NO2")
setnames(dt_d, "LOG_NO2_DIARIO", "LOG_NO2")

dt_h[, NO2_raw := exp(LOG_NO2)]
dt_d[, NO2_raw := exp(LOG_NO2)]

# Datos mensuales
dt_d[, ANIO_MES := as.Date(format(FECHA, "%Y-%m-01"))]
dt_m <- dt_d[, .(
  DATO_NO2                    = mean(DATO_DIARIO,             na.rm = TRUE),
  intensidad_raw              = mean(intensidad_raw,          na.rm = TRUE),
  Temperatura_raw             = mean(Temperatura_raw,         na.rm = TRUE),
  Precipitaciones_raw         = sum(Precipitaciones_raw,      na.rm = TRUE),
  `Presion Barométrica_raw`   = mean(`Presion Barométrica_raw`, na.rm = TRUE),
  `Velocidad Viento_raw`      = mean(`Velocidad Viento_raw`,  na.rm = TRUE),
  n_dias                      = .N
), by = .(ESTACION, ANIO_MES, barrio, distrito, LONGITUD, LATITUD, ID_DISTRITO)]

setnames(dt_m, "ANIO_MES", "FECHA")
dt_m[, LOG_NO2 := log(DATO_NO2)]
dt_m[, NO2_raw := exp(LOG_NO2)]

# Re-estandarizar covariables a escala mensual
cols_raw_m <- c("intensidad_raw", "Temperatura_raw",
                "Precipitaciones_raw", "Presion Barométrica_raw",
                "Velocidad Viento_raw")
for (i in seq_along(COVS_NOMBRES)) {
  dt_m[, (COVS_NOMBRES[i]) := scale(get(cols_raw_m[i]))[, 1]]
}

cat("  Horario: ", nrow(dt_h), "observaciones |", uniqueN(dt_h$ESTACION),
    "estaciones\n")
cat("  Diario:  ", nrow(dt_d), "observaciones |", uniqueN(dt_d$ESTACION),
    "estaciones\n")
cat("  Mensual: ", nrow(dt_m), "observaciones |", uniqueN(dt_m$ESTACION),
    "estaciones\n")

# ==============================================================================
# FUNCIÓN: CREAR RESUMEN EXPLORATORIO POR FRECUENCIA
# ==============================================================================

explorar_frecuencia <- function(dt, frecuencia, periodo_sel = NULL) {

  if (!is.null(periodo_sel)) {
    if (frecuencia %in% c("horario", "diario")) {
      mes <- as.numeric(format(dt$FECHA, "%m"))[1]
      if (periodo_sel == "noviembre") {
        dt <- dt[format(FECHA, "%m") == "11"]
      } else if (periodo_sel == "marzo") {
        dt <- dt[format(FECHA, "%m") == "03"]
      }
    }
  }

  dt <- dt[complete.cases(dt[, .(NO2_raw)])]

  carpeta <- here("outputs", "exploracion_visual", frecuencia)
  if (!is.null(periodo_sel)) {
    carpeta <- file.path(carpeta, tolower(periodo_sel))
  }
  dir.create(carpeta, showWarnings = FALSE, recursive = TRUE)

  cat("\n  ├─", toupper(frecuencia))
  if (!is.null(periodo_sel)) cat(" — ", toupper(periodo_sel))
  cat("\n")

  # ---------- 1. Distribuciones univariantes ----------
  cat("  │  ├─ Distribuciones univariantes...\n")

  cols_raw_disp <- intersect(c("intensidad_raw", "Temperatura_raw",
                               "Precipitaciones_raw", "Presion Barométrica_raw",
                               "Velocidad Viento_raw"),
                             names(dt))
  etiquetas <- c(
    intensidad_raw            = "Intensidad tráfico",
    Temperatura_raw           = "Temperatura (°C)",
    Precipitaciones_raw       = "Precipitaciones (mm)",
    `Presion Barométrica_raw` = "Presión (hPa)",
    `Velocidad Viento_raw`    = "Velocidad viento (m/s)"
  )

  plots_dist <- lapply(c("NO2_raw", cols_raw_disp), function(v) {
    lab <- ifelse(v == "NO2_raw", "NO2 (µg/m³)", etiquetas[v])
    ggplot(dt, aes(x = !!rlang::sym(v))) +
      geom_histogram(bins = 30, fill = "#2166AC", alpha = 0.7, color = "white") +
      labs(title = lab, x = NULL, y = "Frecuencia") +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 9))
  })

  n_plots <- length(plots_dist)
  n_cols  <- min(3, n_plots)
  n_rows  <- ceiling(n_plots / n_cols)

  p_dist <- Reduce(`+`, plots_dist) +
    plot_layout(ncol = n_cols) +
    plot_annotation(title = "Distribuciones univariantes",
                    theme = theme(plot.title = element_text(face = "bold")))

  ggsave(file.path(carpeta, "01_distribuciones.png"),
         plot = p_dist, width = n_cols * 4, height = n_rows * 3.5, dpi = 300)

  # ---------- 2. Heatmap correlación (NO2 raw vs covariables) ----------
  cat("  │  ├─ Correlación (NO2 raw vs covariables)...\n")

  mat <- dt[, c(..cols_raw_disp)]
  mat[, NO2_raw := dt$NO2_raw]

  nombres_plot <- c(etiquetas[cols_raw_disp], "NO2 (µg/m³)")
  setnames(mat, names(mat), nombres_plot)

  cor_mat <- cor(mat, use = "pairwise.complete.obs")
  cor_df  <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Var1", "Var2", "Correlacion")

  niveles <- c(setdiff(nombres_plot, "NO2 (µg/m³)"), "NO2 (µg/m³)")
  cor_df$Var1 <- factor(cor_df$Var1, levels = niveles)
  cor_df$Var2 <- factor(cor_df$Var2, levels = rev(niveles))

  p_cor <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlacion)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", Correlacion)),
              size = 3, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = "Correlación: NO2 (escala original) vs covariables",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title  = element_text(face = "bold", size = 10))

  ggsave(file.path(carpeta, "02_heatmap_correlacion.png"),
         plot = p_cor, width = 8, height = 6, dpi = 300)

  # ---------- 3. Linealidad: scatter + suavizadores ----------
  cat("  │  ├─ Análisis de linealidad...\n")

  plots_lin <- lapply(cols_raw_disp, function(v) {
    df_p <- na.omit(data.frame(cov = dt[[v]], no2 = dt$NO2_raw))
    if (nrow(df_p) < 10) return(NULL)

    ggplot(df_p, aes(x = cov, y = no2)) +
      geom_point(alpha = 0.2, size = 0.8, color = "grey40") +
      geom_smooth(aes(color = "Lineal (OLS)"),
                  method = "lm",    se = TRUE, linewidth = 0.9) +
      geom_smooth(aes(color = "Smooth (LOESS)"),
                  method = "loess", se = TRUE, linewidth = 0.9,
                  linetype = "dashed") +
      scale_color_manual(
        values = c("Lineal (OLS)" = "#B2182B", "Smooth (LOESS)" = "#2166AC")
      ) +
      labs(title = etiquetas[v],
           x = etiquetas[v], y = "NO2 (µg/m³)", color = NULL) +
      theme_minimal(base_size = 9) +
      theme(plot.title    = element_text(face = "bold", size = 9),
            legend.position = "bottom",
            legend.text     = element_text(size = 7))
  })

  plots_lin <- Filter(Negate(is.null), plots_lin)
  if (length(plots_lin) > 0) {
    n_plots_lin <- length(plots_lin)
    n_cols_lin  <- min(3, n_plots_lin)
    n_rows_lin  <- ceiling(n_plots_lin / n_cols_lin)

    p_lin <- Reduce(`+`, plots_lin) +
      plot_layout(ncol = n_cols_lin, guides = "collect") +
      plot_annotation(
        title = "Linealidad: NO2 vs covariables (rojo=OLS, azul=LOESS)",
        theme = theme(plot.title = element_text(face = "bold", size = 11))
      ) &
      theme(legend.position = "bottom")

    ggsave(file.path(carpeta, "03_linealidad_scatter.png"),
           plot = p_lin, width = n_cols_lin * 4.5, height = n_rows_lin * 4,
           dpi = 300)
  }

  # ---------- 4. Tabla de linealidad (Pearson/Spearman/AIC) ----------
  cat("  │  ├─ Tabla de linealidad...\n")

  tabla_lin <- lapply(cols_raw_disp, function(v) {
    df_v <- na.omit(data.frame(cov = dt[[v]], no2 = dt$NO2_raw))
    if (nrow(df_v) < 10) return(NULL)

    r_p <- cor(df_v$cov, df_v$no2, method = "pearson")
    r_s <- cor(df_v$cov, df_v$no2, method = "spearman")

    tryCatch({
      lm_fit  <- lm(no2 ~ cov, data = df_v)
      gam_fit <- mgcv::gam(no2 ~ s(cov, k = 4), data = df_v, method = "REML")

      aic_lm  <- AIC(lm_fit)
      aic_gam <- AIC(gam_fit)
      delta   <- round(aic_lm - aic_gam, 1)

      no_lineal <- delta > 2 | abs(r_p - r_s) > 0.10

      data.frame(
        Variable   = etiquetas[v],
        r_Pearson  = round(r_p, 3),
        r_Spearman = round(r_s, 3),
        Dif_abs    = round(abs(r_p - r_s), 3),
        AIC_LM     = round(aic_lm, 1),
        AIC_GAM    = round(aic_gam, 1),
        ΔAIC       = delta,
        Tipo       = ifelse(no_lineal, "No lineal", "Lineal"),
        row.names  = NULL,
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)
  })

  tabla_lin <- do.call(rbind, Filter(Negate(is.null), tabla_lin))

  if (!is.null(tabla_lin) && nrow(tabla_lin) > 0) {
    tbl <- gt::gt(tabla_lin) |>
      gt::tab_header(
        title    = "Análisis de Linealidad",
        subtitle = "Criterios: |Pearson - Spearman| > 0.10 o ΔAIC > 2 → No lineal"
      ) |>
      gt::tab_options(
        table.font.size = px(10),
        heading.title.font.size = px(13),
        heading.subtitle.font.size = px(10),
        column_labels.font.weight = "bold"
      ) |>
      gt::opt_horizontal_padding(scale = 1.5)

    gt::gtsave(tbl, file.path(carpeta, "04_tabla_linealidad.png"))
  }

  # ---------- 5. Resumen de valores faltantes y outliers ----------
  cat("  │  ├─ Resumen de datos faltantes y outliers...\n")

  resumen_datos <- data.frame(
    Variable = c("NO2_raw", cols_raw_disp),
    N_total = sapply(c("NO2_raw", cols_raw_disp), function(v) nrow(dt)),
    N_validos = sapply(c("NO2_raw", cols_raw_disp),
                       function(v) sum(!is.na(dt[[v]]))),
    Pct_faltante = round(sapply(c("NO2_raw", cols_raw_disp),
                                 function(v) 100 * mean(is.na(dt[[v]]))), 2),
    Media = round(sapply(c("NO2_raw", cols_raw_disp),
                         function(v) mean(dt[[v]], na.rm = TRUE)), 2),
    Mediana = round(sapply(c("NO2_raw", cols_raw_disp),
                           function(v) median(dt[[v]], na.rm = TRUE)), 2),
    SD = round(sapply(c("NO2_raw", cols_raw_disp),
                      function(v) sd(dt[[v]], na.rm = TRUE)), 2),
    Min = round(sapply(c("NO2_raw", cols_raw_disp),
                       function(v) min(dt[[v]], na.rm = TRUE)), 2),
    Max = round(sapply(c("NO2_raw", cols_raw_disp),
                       function(v) max(dt[[v]], na.rm = TRUE)), 2),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  tbl_datos <- gt::gt(resumen_datos) |>
    gt::tab_header(
      title = "Resumen descriptivo de variables"
    ) |>
    gt::tab_options(
      table.font.size = px(9),
      heading.title.font.size = px(12),
      column_labels.font.weight = "bold"
    ) |>
    gt::opt_horizontal_padding(scale = 1.2)

  gt::gtsave(tbl_datos, file.path(carpeta, "05_resumen_datos.png"))

  cat("  │  └─ ✓ Completado\n")

  invisible(list(
    n_obs = nrow(dt),
    n_est = uniqueN(dt$ESTACION),
    tabla_lin = tabla_lin,
    resumen = resumen_datos
  ))
}

# ==============================================================================
# EJECUCIÓN: EXPLORACIÓN POR FRECUENCIA Y PERIODO
# ==============================================================================

cat("\n", strrep("=", 80), "\n")
cat("  EXPLORACIÓN VISUAL POR FRECUENCIA\n")
cat(strrep("=", 80), "\n")

# Horario: noviembre y marzo
cat("\n  HORARIO\n")
res_h_nov <- explorar_frecuencia(dt_h, "horario", "noviembre")
res_h_mar <- explorar_frecuencia(dt_h, "horario", "marzo")

# Diario: noviembre y marzo
cat("\n  DIARIO\n")
res_d_nov <- explorar_frecuencia(dt_d, "diario", "noviembre")
res_d_mar <- explorar_frecuencia(dt_d, "diario", "marzo")

# Mensual: año completo
cat("\n  MENSUAL\n")
res_m <- explorar_frecuencia(dt_m, "mensual")

cat("\n", strrep("=", 80), "\n")
cat("  EXPLORACIÓN COMPLETADA\n")
cat(strrep("=", 80), "\n")
cat("  Outputs en: outputs/exploracion_visual/\n")
cat("    └─ horario/  (noviembre, marzo)\n")
cat("    └─ diario/   (noviembre, marzo)\n")
cat("    └─ mensual/\n")
