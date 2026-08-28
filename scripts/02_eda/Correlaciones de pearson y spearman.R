
# =============================================================================
# SIMULACION_INTENTO_2.R
# Análisis de correlación entre covariables y NO2
# Dataset Maestro Horario 2025
# =============================================================================

library(tidyverse)
library(ggcorrplot)

# -----------------------------------------------------------------------------
# Configuración general
# -----------------------------------------------------------------------------

ESTACIONES_INTERES <- c(
  "Pza. del Carmen",
  "Casa de Campo",
  "Pº Castellana",
  "Plaza Elíptica",
  "Arturo Soria",
  "Vallecas",
  "Barrio del Pilar"   # Estación más próxima a Peñagrande en datos 2025
)

ETIQUETAS_ESTACIONES <- c(
  "Pza. del Carmen"  = "Plaza del Carmen",
  "Casa de Campo"    = "Casa de Campo",
  "Pº Castellana"    = "Castellana",
  "Plaza Elíptica"   = "Plaza Elíptica",
  "Arturo Soria"     = "Arturo Soria",
  "Vallecas"         = "Vallecas",
  "Barrio del Pilar" = "Peñagrande (B.Pilar)"
)

# Carpetas de resultados
DIR_BASE       <- file.path("outputs", "analysis", "correlacion_scatter", "horario")
DIR_PEARSON    <- file.path(DIR_BASE, "correlacion de pearson")
DIR_SPEARMAN   <- file.path(DIR_BASE, "spearman")
DIR_SCATTER    <- file.path(DIR_BASE, "scatter plots")
DIR_COMPARATIVA <- file.path(DIR_BASE, "comparativa estaciones")

for (d in c(DIR_PEARSON, DIR_SPEARMAN, DIR_SCATTER, DIR_COMPARATIVA)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# 1. CARGA Y PREPARACIÓN DE DATOS
# =============================================================================

datos_h <- readRDS("data/processed/Maestro/horario/dataset_maestro_inla_2025_HORARIO.rds")

covariables_raw <- c(
  "Temperatura_raw",
  "Humedad_Relativa_raw",
  "Presion Barométrica_raw",
  "Radiación Solar_raw",
  "Velocidad Viento_raw",
  "intensidad_raw",
  "carga_raw",
  "Precipitaciones_raw"
)

df_sel <- datos_h |>
  filter(ESTACION %in% ESTACIONES_INTERES) |>
  mutate(
    ESTACION_LABEL = ETIQUETAS_ESTACIONES[ESTACION],
    NO2            = DATO,
    log_NO2        = LOG_NO2_HORARIO
  ) |>
  select(ESTACION, ESTACION_LABEL, FECHA, HORA, NO2, log_NO2,
         all_of(covariables_raw)) |>
  rename(
    "Temperatura (°C)"          = Temperatura_raw,
    "Humedad Relativa (%)"      = Humedad_Relativa_raw,
    "Presión Barométrica (hPa)" = `Presion Barométrica_raw`,
    "Radiación Solar (W/m²)"    = `Radiación Solar_raw`,
    "Velocidad Viento (m/s)"    = `Velocidad Viento_raw`,
    "Intensidad Tráfico"        = intensidad_raw,
    "Carga Tráfico"             = carga_raw,
    "Precipitaciones (mm)"      = Precipitaciones_raw
  ) |>
  drop_na()

covs_label <- c(
  "Temperatura (°C)",
  "Humedad Relativa (%)",
  "Presión Barométrica (hPa)",
  "Radiación Solar (W/m²)",
  "Velocidad Viento (m/s)",
  "Intensidad Tráfico",
  "Carga Tráfico",
  "Precipitaciones (mm)"
)

cat("Observaciones tras filtrado:", nrow(df_sel), "\n")
cat("Estaciones:", paste(unique(df_sel$ESTACION_LABEL), collapse = ", "), "\n")

# =============================================================================
# FUNCIÓN AUXILIAR: genera y guarda un heatmap de correlación
# =============================================================================

guardar_heatmap <- function(mat, metodo, titulo, dir_out, nombre_archivo) {

  color_titulo <- if (metodo == "pearson") "#1A5276" else "#6E2F91"

  p <- ggcorrplot(
    mat,
    method   = "square",
    type     = "full",
    lab      = TRUE,
    lab_size = 3.2,
    tl.cex   = 9,
    colors   = c("#4575B4", "white", "#D73027"),
    title    = titulo,
    ggtheme  = theme_minimal(base_size = 10)
  ) +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 11, color = color_titulo
      )
    )

  ggsave(file.path(dir_out, nombre_archivo), p, width = 10, height = 9, dpi = 180)
  cat("  Guardado:", file.path(dir_out, nombre_archivo), "\n")
}

# Helper: matriz de correlación incluyendo NO2 y log(NO2)
build_mat <- function(df, metodo) {
  df |>
    select(NO2, log_NO2, all_of(covs_label)) |>
    rename("NO2" = NO2, "log(NO2)" = log_NO2) |>
    cor(method = metodo, use = "complete.obs")
}

# =============================================================================
# 2. NIVEL HORARIO — HEATMAPS DE CORRELACIÓN
# =============================================================================

## ── 2a. Pearson — global ─────────────────────────────────────────────────────
cat("\n--- Pearson: heatmap global ---\n")
guardar_heatmap(
  mat            = build_mat(df_sel, "pearson"),
  metodo         = "pearson",
  titulo         = "Correlación de Pearson — Nivel Horario 2025\nTodas las estaciones seleccionadas",
  dir_out        = DIR_PEARSON,
  nombre_archivo = "heatmap_pearson_global_horario.png"
)

## ── 2b. Pearson — por estación ───────────────────────────────────────────────
cat("--- Pearson: heatmaps por estación ---\n")
for (est in unique(df_sel$ESTACION_LABEL)) {
  df_est <- df_sel |> filter(ESTACION_LABEL == est)
  if (nrow(df_est) < 30) next
  guardar_heatmap(
    mat            = build_mat(df_est, "pearson"),
    metodo         = "pearson",
    titulo         = paste0("Correlación de Pearson — Nivel Horario 2025\n", est),
    dir_out        = DIR_PEARSON,
    nombre_archivo = paste0("heatmap_pearson_", gsub("[^[:alnum:]]", "_", est), "_horario.png")
  )
}

## ── 2c. Spearman — global ────────────────────────────────────────────────────
cat("\n--- Spearman: heatmap global ---\n")
guardar_heatmap(
  mat            = build_mat(df_sel, "spearman"),
  metodo         = "spearman",
  titulo         = "Correlación de Spearman — Nivel Horario 2025\nTodas las estaciones seleccionadas",
  dir_out        = DIR_SPEARMAN,
  nombre_archivo = "heatmap_spearman_global_horario.png"
)

## ── 2d. Spearman — por estación ──────────────────────────────────────────────
cat("--- Spearman: heatmaps por estación ---\n")
for (est in unique(df_sel$ESTACION_LABEL)) {
  df_est <- df_sel |> filter(ESTACION_LABEL == est)
  if (nrow(df_est) < 30) next
  guardar_heatmap(
    mat            = build_mat(df_est, "spearman"),
    metodo         = "spearman",
    titulo         = paste0("Correlación de Spearman — Nivel Horario 2025\n", est),
    dir_out        = DIR_SPEARMAN,
    nombre_archivo = paste0("heatmap_spearman_", gsub("[^[:alnum:]]", "_", est), "_horario.png")
  )
}

# =============================================================================
# 3. NIVEL HORARIO — SCATTER PLOTS (NO2 y log(NO2) vs covariables)
# =============================================================================
# Un archivo PNG por covariable; dos filas (NO2 / log NO2) × 7 columnas (estaciones)

cat("\n--- Scatter plots ---\n")

for (cov in covs_label) {

  df_long <- df_sel |>
    select(ESTACION_LABEL, all_of(cov), NO2, log_NO2) |>
    pivot_longer(
      cols      = c(NO2, log_NO2),
      names_to  = "tipo_NO2",
      values_to = "valor_NO2"
    ) |>
    mutate(
      tipo_NO2 = if_else(
        tipo_NO2 == "NO2",
        "NO\u2082 (\u03bcg/m\u00b3)",
        "log(NO\u2082)"
      )
    )

  p_scatter <- ggplot(df_long, aes(x = .data[[cov]], y = valor_NO2)) +
    geom_point(alpha = 0.07, size = 0.6) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.9,
                color = "#D73027") +
    facet_grid(tipo_NO2 ~ ESTACION_LABEL, scales = "free_y") +
    labs(
      title   = paste0("Relación entre ", cov, " y NO\u2082 — Nivel Horario 2025"),
      x       = cov,
      y       = NULL,
      caption = "Línea roja: suavizado LOESS con banda de confianza al 95%"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 10),
      strip.text   = element_text(size = 7.5),
      axis.text.x  = element_text(angle = 30, hjust = 1, size = 7),
      plot.caption = element_text(size = 7, hjust = 1)
    )

  nombre_archivo <- paste0(
    "scatter_", gsub("[^[:alnum:]]", "_", cov), "_horario.png"
  )
  ggsave(file.path(DIR_SCATTER, nombre_archivo), p_scatter,
         width = 18, height = 6, dpi = 160)
  cat("  Guardado:", nombre_archivo, "\n")
}

# =============================================================================
# 4. NIVEL HORARIO — COMPARATIVA ENTRE ESTACIONES (gráfico de líneas)
# =============================================================================
# Para cada método (Pearson / Spearman) y cada variable respuesta (NO2 / log NO2):
#   X = covariable  |  Y = coeficiente de correlación  |  líneas = estación
# Se generan 4 gráficos:
#   • pearson_NO2        • pearson_logNO2
#   • spearman_NO2       • spearman_logNO2
# Y uno combinado con los 4 paneles juntos.

cat("\n--- Comparativa entre estaciones ---\n")

# Calcular correlaciones por estación, método y variable respuesta
cor_df <- expand.grid(
  ESTACION_LABEL = unique(df_sel$ESTACION_LABEL),
  Covariable     = covs_label,
  Metodo         = c("pearson", "spearman"),
  Respuesta      = c("NO2", "log_NO2"),
  stringsAsFactors = FALSE
) |>
  rowwise() |>
  mutate(
    cor_val = {
      # Capturar el valor escalar de la fila actual con nombre distinto
      # para evitar que filter() lo confunda con la columna de df_sel
      est_actual <- ESTACION_LABEL
      cov_actual <- Covariable
      met_actual <- Metodo
      res_actual <- Respuesta

      df_est <- df_sel[df_sel$ESTACION_LABEL == est_actual, ]
      x <- df_est[[cov_actual]]
      y <- if (res_actual == "NO2") df_est$NO2 else df_est$log_NO2
      suppressWarnings(cor(x, y, method = met_actual, use = "complete.obs"))
    }
  ) |>
  ungroup() |>
  mutate(
    Metodo_label   = if_else(Metodo == "pearson", "Pearson", "Spearman"),
    Respuesta_label = if_else(Respuesta == "NO2",
                              "NO\u2082 (\u03bcg/m\u00b3)",
                              "log(NO\u2082)")
  )

# Orden de covariables en el eje X (de mayor a menor correlación absoluta media)
orden_covs <- cor_df |>
  group_by(Covariable) |>
  summarise(abs_mean = mean(abs(cor_val), na.rm = TRUE)) |>
  arrange(desc(abs_mean)) |>
  pull(Covariable)

cor_df <- cor_df |>
  mutate(Covariable = factor(Covariable, levels = orden_covs))

# Paleta de colores para las 7 estaciones
paleta_estaciones <- c(
  "Plaza del Carmen"      = "#E41A1C",
  "Casa de Campo"         = "#377EB8",
  "Castellana"            = "#4DAF4A",
  "Plaza Elíptica"        = "#984EA3",
  "Arturo Soria"          = "#FF7F00",
  "Vallecas"              = "#A65628",
  "Peñagrande (B.Pilar)"  = "#F781BF"
)

# ── Función para generar un gráfico de líneas de comparativa ─────────────────
plot_comparativa <- function(datos, titulo, color_titulo = "black") {
  ggplot(datos, aes(x = Covariable, y = cor_val,
                    color = ESTACION_LABEL, group = ESTACION_LABEL)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    scale_color_manual(values = paleta_estaciones, name = "Estación") +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    labs(
      title = titulo,
      x     = NULL,
      y     = "Coeficiente de correlación"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 11,
                                  color = color_titulo),
      axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
}

# ── 4a. Gráficos individuales (un PNG por combinación método × respuesta) ─────
combinaciones <- list(
  list(metodo = "pearson",  respuesta = "NO2",
       label_r = "NO\u2082 (\u03bcg/m\u00b3)",
       titulo = "Comparativa de correlación de Pearson vs NO\u2082 — Nivel Horario 2025",
       archivo = "comparativa_pearson_NO2_horario.png",
       col_tit = "#1A5276"),
  list(metodo = "pearson",  respuesta = "log_NO2",
       label_r = "log(NO\u2082)",
       titulo = "Comparativa de correlación de Pearson vs log(NO\u2082) — Nivel Horario 2025",
       archivo = "comparativa_pearson_logNO2_horario.png",
       col_tit = "#1A5276"),
  list(metodo = "spearman", respuesta = "NO2",
       label_r = "NO\u2082 (\u03bcg/m\u00b3)",
       titulo = "Comparativa de correlación de Spearman vs NO\u2082 — Nivel Horario 2025",
       archivo = "comparativa_spearman_NO2_horario.png",
       col_tit = "#6E2F91"),
  list(metodo = "spearman", respuesta = "log_NO2",
       label_r = "log(NO\u2082)",
       titulo = "Comparativa de correlación de Spearman vs log(NO\u2082) — Nivel Horario 2025",
       archivo = "comparativa_spearman_logNO2_horario.png",
       col_tit = "#6E2F91")
)

for (cfg in combinaciones) {
  datos_sub <- cor_df |>
    filter(Metodo == cfg$metodo, Respuesta == cfg$respuesta)

  p <- plot_comparativa(datos_sub, cfg$titulo, cfg$col_tit)
  ggsave(file.path(DIR_COMPARATIVA, cfg$archivo), p,
         width = 11, height = 6, dpi = 180)
  cat("  Guardado:", cfg$archivo, "\n")
}

# ── 4b. Panel combinado (2×2: método × variable respuesta) ───────────────────
p_panel <- cor_df |>
  mutate(
    Metodo_label    = factor(Metodo_label,    levels = c("Pearson", "Spearman")),
    Respuesta_label = factor(Respuesta_label,
                             levels = c("NO\u2082 (\u03bcg/m\u00b3)", "log(NO\u2082)"))
  ) |>
  ggplot(aes(x = Covariable, y = cor_val,
             color = ESTACION_LABEL, group = ESTACION_LABEL)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_grid(Metodo_label ~ Respuesta_label) +
  scale_color_manual(values = paleta_estaciones, name = "Estación") +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(
    title   = "Comparativa de correlación entre estaciones — Nivel Horario 2025",
    x       = NULL,
    y       = "Coeficiente de correlación",
    caption = "Pearson (lineal) vs Spearman (rango); NO\u2082 en escala original y logarítmica"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 7.5),
    legend.position = "right",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8),
    strip.text      = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption    = element_text(size = 7, hjust = 1)
  )

ggsave(file.path(DIR_COMPARATIVA, "comparativa_panel_2x2_horario.png"),
       p_panel, width = 14, height = 9, dpi = 180)
cat("  Guardado: comparativa_panel_2x2_horario.png\n")

# =============================================================================
# 5. ESTRATIFICACIÓN POR ESTACIÓN DEL AÑO — HEATMAPS ESTACIONALES
# =============================================================================
# Objetivo: comparar cómo cambian las correlaciones entre covariables y NO2
# a lo largo del año (Invierno / Primavera / Verano / Otoño).
# Se generan heatmaps globales (las 7 estaciones conjuntas) y por estación
# meteorológica. Carpeta separada: estratificacion estacional/

cat("\n--- Estratificación estacional ---\n")

DIR_ESTACIONAL         <- file.path(DIR_BASE, "estratificacion estacional")
DIR_ESTACIONAL_PEARSON <- file.path(DIR_ESTACIONAL, "pearson")
DIR_ESTACIONAL_SPEARMAN<- file.path(DIR_ESTACIONAL, "spearman")

for (d in c(DIR_ESTACIONAL_PEARSON, DIR_ESTACIONAL_SPEARMAN)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Definición de estaciones del año (hemisferio norte, Madrid)
ESTACIONES_ANNO <- list(
  Invierno  = c(12L, 1L, 2L),
  Primavera = c(3L, 4L, 5L),
  Verano    = c(6L, 7L, 8L),
  "Oto\u00f1o"  = c(9L, 10L, 11L)
)

# Añadir mes y estación del año al data frame filtrado
df_estacional <- df_sel |>
  mutate(
    Mes = as.integer(format(as.Date(FECHA), "%m")),
    Estacion_Anno = case_when(
      Mes %in% ESTACIONES_ANNO$Invierno  ~ "Invierno",
      Mes %in% ESTACIONES_ANNO$Primavera ~ "Primavera",
      Mes %in% ESTACIONES_ANNO$Verano    ~ "Verano",
      Mes %in% ESTACIONES_ANNO[["Oto\u00f1o"]] ~ "Oto\u00f1o"
    ),
    Estacion_Anno = factor(
      Estacion_Anno,
      levels = c("Invierno", "Primavera", "Verano", "Oto\u00f1o")
    )
  )

# Colores de título por estación del año
colores_estacion <- c(
  Invierno  = "#1F618D",   # azul frío
  Primavera = "#1E8449",   # verde
  Verano    = "#B7770D",   # naranja cálido
  "Oto\u00f1o"    = "#784212"    # marrón
)

# ── 5a. Heatmaps globales por estación del año (Pearson + Spearman) ──────────

for (estacion_anno in levels(df_estacional$Estacion_Anno)) {

  df_tmp <- df_estacional |> filter(Estacion_Anno == estacion_anno)
  n_obs  <- nrow(df_tmp)
  cat("  ", estacion_anno, "— n =", n_obs, "\n")
  if (n_obs < 30) next

  for (metodo in c("pearson", "spearman")) {

    mat_tmp <- df_tmp |>
      select(NO2, log_NO2, all_of(covs_label)) |>
      rename("NO2" = NO2, "log(NO2)" = log_NO2) |>
      cor(method = metodo, use = "complete.obs")

    dir_tmp <- if (metodo == "pearson") DIR_ESTACIONAL_PEARSON else DIR_ESTACIONAL_SPEARMAN
    metodo_label <- if (metodo == "pearson") "Pearson" else "Spearman"
    col_tit <- colores_estacion[[estacion_anno]]

    guardar_heatmap(
      mat            = mat_tmp,
      metodo         = metodo,
      titulo         = paste0(
        "Correlaci\u00f3n de ", metodo_label,
        " — ", estacion_anno,
        " 2025\nTodas las estaciones seleccionadas"
      ),
      dir_out        = dir_tmp,
      nombre_archivo = paste0(
        "heatmap_", metodo, "_",
        tolower(gsub("\u00f1", "n", estacion_anno)),
        "_global_horario.png"
      )
    )
  }
}

# ── 5b. Heatmaps por estación meteorológica × estación del año ───────────────

for (estacion_anno in levels(df_estacional$Estacion_Anno)) {
  df_tmp <- df_estacional |> filter(Estacion_Anno == estacion_anno)

  for (est in unique(df_tmp$ESTACION_LABEL)) {
    df_est_tmp <- df_tmp |> filter(ESTACION_LABEL == est)
    if (nrow(df_est_tmp) < 30) next

    for (metodo in c("pearson", "spearman")) {

      mat_tmp <- df_est_tmp |>
        select(NO2, log_NO2, all_of(covs_label)) |>
        rename("NO2" = NO2, "log(NO2)" = log_NO2) |>
        cor(method = metodo, use = "complete.obs")

      dir_tmp      <- if (metodo == "pearson") DIR_ESTACIONAL_PEARSON else DIR_ESTACIONAL_SPEARMAN
      metodo_label <- if (metodo == "pearson") "Pearson" else "Spearman"

      guardar_heatmap(
        mat            = mat_tmp,
        metodo         = metodo,
        titulo         = paste0(
          "Correlaci\u00f3n de ", metodo_label,
          " — ", estacion_anno, " 2025\n", est
        ),
        dir_out        = dir_tmp,
        nombre_archivo = paste0(
          "heatmap_", metodo, "_",
          tolower(gsub("\u00f1", "n", estacion_anno)), "_",
          gsub("[^[:alnum:]]", "_", est),
          "_horario.png"
        )
      )
    }
  }
}

# ── 5c. Panel comparativo 1×4 (las 4 estaciones del año en una sola imagen) ──
# Muestra los coeficientes de correlación de cada covariable con NO2 y log(NO2)
# por estación del año → visualmente muy potente para la comparativa en el TFM

cat("\n  Generando panel comparativo estacional...\n")

# Calcular correlaciones por (estación del año × covariable × método × respuesta)
# para las 7 estaciones de medición juntas
cor_estacional <- expand.grid(
  Estacion_Anno = levels(df_estacional$Estacion_Anno),
  Covariable    = covs_label,
  Metodo        = c("pearson", "spearman"),
  Respuesta     = c("NO2", "log_NO2"),
  stringsAsFactors = FALSE
) |>
  rowwise() |>
  mutate(
    cor_val = {
      ea  <- Estacion_Anno
      cov <- Covariable
      met <- Metodo
      res <- Respuesta
      df_tmp <- df_estacional[df_estacional$Estacion_Anno == ea, ]
      x <- df_tmp[[cov]]
      y <- if (res == "NO2") df_tmp$NO2 else df_tmp$log_NO2
      suppressWarnings(cor(x, y, method = met, use = "complete.obs"))
    }
  ) |>
  ungroup() |>
  mutate(
    Metodo_label    = if_else(Metodo == "pearson", "Pearson", "Spearman"),
    Respuesta_label = if_else(Respuesta == "NO2",
                              "NO\u2082 (\u03bcg/m\u00b3)", "log(NO\u2082)"),
    Estacion_Anno   = factor(
      Estacion_Anno,
      levels = c("Invierno", "Primavera", "Verano", "Oto\u00f1o")
    ),
    Covariable = factor(Covariable, levels = covs_label)
  )

# Paleta y forma para los 4 periodos estacionales
paleta_anno <- c(
  Invierno  = "#1F618D",
  Primavera = "#1E8449",
  Verano    = "#B7770D",
  "Oto\u00f1o"  = "#784212"
)

# Panel 2×2: filas = método, columnas = variable respuesta
# Dentro de cada panel: línea por estación del año, X = covariable
p_estacional_panel <- cor_estacional |>
  mutate(
    Metodo_label    = factor(Metodo_label,    levels = c("Pearson", "Spearman")),
    Respuesta_label = factor(Respuesta_label,
                             levels = c("NO\u2082 (\u03bcg/m\u00b3)", "log(NO\u2082)"))
  ) |>
  ggplot(aes(x = Covariable, y = cor_val,
             color = Estacion_Anno, group = Estacion_Anno)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_grid(Metodo_label ~ Respuesta_label) +
  scale_color_manual(values = paleta_anno, name = "Estaci\u00f3n del a\u00f1o") +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(
    title   = "Comparativa estacional de correlaci\u00f3n con NO\u2082 — Nivel Horario 2025",
    x       = NULL,
    y       = "Coeficiente de correlaci\u00f3n",
    caption = paste0(
      "Invierno: dic-feb \u2022 Primavera: mar-may \u2022 ",
      "Verano: jun-ago \u2022 Oto\u00f1o: sep-nov\n",
      "Todas las estaciones de medici\u00f3n seleccionadas conjuntas"
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
    legend.position = "right",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    strip.text      = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption    = element_text(size = 7.5, hjust = 1, color = "grey40")
  )

ggsave(
  file.path(DIR_ESTACIONAL, "panel_comparativo_estacional_horario.png"),
  p_estacional_panel, width = 14, height = 9, dpi = 180
)
cat("  Guardado: panel_comparativo_estacional_horario.png\n")

# ── 5d. Panel visual 1×4: cuatro heatmaps lado a lado (cowplot) ──────────────
# Un heatmap por estación del año (Pearson, NO2 + log_NO2 juntos)
# Muy impactante para presentaciones y TFM

if (requireNamespace("cowplot", quietly = TRUE)) {

  library(cowplot)

  recoger_heatmap_estacional <- function(estacion_anno, metodo) {

    df_tmp   <- df_estacional[df_estacional$Estacion_Anno == estacion_anno, ]
    mat_tmp  <- df_tmp |>
      select(NO2, log_NO2, all_of(covs_label)) |>
      rename("NO2" = NO2, "log(NO2)" = log_NO2) |>
      cor(method = metodo, use = "complete.obs")

    metodo_label <- if (metodo == "pearson") "Pearson" else "Spearman"
    col_tit <- colores_estacion[[estacion_anno]]

    ggcorrplot(
      mat_tmp,
      method   = "square",
      type     = "full",
      lab      = TRUE,
      lab_size = 2.4,
      tl.cex   = 7.5,
      colors   = c("#4575B4", "white", "#D73027"),
      title    = paste0(estacion_anno),
      ggtheme  = theme_minimal(base_size = 8)
    ) +
      theme(
        plot.title = element_text(
          hjust = 0.5, face = "bold", size = 10, color = col_tit
        ),
        plot.margin = margin(4, 4, 4, 4)
      )
  }

  for (metodo in c("pearson", "spearman")) {

    lista_plots <- lapply(
      levels(df_estacional$Estacion_Anno),
      recoger_heatmap_estacional,
      metodo = metodo
    )

    metodo_label <- if (metodo == "pearson") "Pearson" else "Spearman"
    titulo_panel <- paste0(
      "Heatmaps estacionales de correlaci\u00f3n de ", metodo_label,
      " — Nivel Horario 2025"
    )

    p_grid <- plot_grid(
      plotlist = lista_plots,
      nrow     = 1,
      ncol     = 4,
      labels   = NULL
    )
    p_titulo <- ggdraw() +
      draw_label(titulo_panel, fontface = "bold", size = 11, hjust = 0.5)

    p_final <- plot_grid(p_titulo, p_grid, ncol = 1, rel_heights = c(0.06, 1))

    archivo_grid <- paste0("heatmaps_1x4_", metodo, "_estacional_horario.png")
    ggsave(
      file.path(DIR_ESTACIONAL, archivo_grid),
      p_final, width = 20, height = 6, dpi = 180
    )
    cat("  Guardado:", archivo_grid, "\n")
  }

} else {
  message("Instala 'cowplot' para generar el panel 1x4 de heatmaps: install.packages('cowplot')")
}

cat("\n=== An\u00e1lisis de correlaci\u00f3n (nivel horario) completado. ===\n")
cat("Resultados en:", DIR_BASE, "\n")

# =============================================================================
# FUNCIÓN GENÉRICA — reutilizada para nivel DIARIO y MENSUAL
# =============================================================================
# Parámetros:
#   datos       → data frame ya cargado
#   col_no2     → nombre de columna con el valor de NO2
#   col_log_no2 → nombre de columna con log(NO2)
#   col_fecha   → nombre de columna de fecha (Date o character YYYY-MM)
#   nivel       → "diario" | "mensual"  (se usa en títulos y nombres de archivo)
#   dir_nivel   → carpeta raíz de salida para este nivel

run_analisis_correlacion <- function(datos, col_no2, col_log_no2,
                                     col_fecha, nivel, dir_nivel) {

  # ── Subcarpetas ─────────────────────────────────────────────────────────────
  DIR_P   <- file.path(dir_nivel, "correlacion de pearson")
  DIR_S   <- file.path(dir_nivel, "spearman")
  DIR_SC  <- file.path(dir_nivel, "scatter plots")
  DIR_CMP <- file.path(dir_nivel, "comparativa estaciones")
  DIR_EST <- file.path(dir_nivel, "estratificacion estacional")
  DIR_EP  <- file.path(DIR_EST, "pearson")
  DIR_ES  <- file.path(DIR_EST, "spearman")
  for (d in c(DIR_P, DIR_S, DIR_SC, DIR_CMP, DIR_EP, DIR_ES)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # ── Preparación del data frame ───────────────────────────────────────────────
  # Detectar nombres de columnas con encoding variable mediante grep
  col_presion_raw <- grep("^Presion Bar.*_raw$", names(datos), value = TRUE)[1]
  col_solar_raw   <- grep("Radi.*Solar_raw$",    names(datos), value = TRUE)[1]
  col_viento_raw  <- grep("Velocidad Viento_raw", names(datos), value = TRUE)[1]

  df <- datos |>
    filter(ESTACION %in% ESTACIONES_INTERES) |>
    mutate(
      ESTACION_LABEL = ETIQUETAS_ESTACIONES[ESTACION],
      NO2     = .data[[col_no2]],
      log_NO2 = .data[[col_log_no2]],
      FECHA_  = as.Date(if (col_fecha == "MES") paste0(.data[[col_fecha]], "-01")
                        else .data[[col_fecha]])
    ) |>
    select(ESTACION, ESTACION_LABEL, FECHA_, NO2, log_NO2,
           all_of(covariables_raw)) |>
    rename(
      "Temperatura (\u00b0C)"          = Temperatura_raw,
      "Humedad Relativa (%)"           = Humedad_Relativa_raw,
      "Presi\u00f3n Barom\u00e9trica (hPa)" = !!sym(col_presion_raw),
      "Radiaci\u00f3n Solar (W/m\u00b2)"    = !!sym(col_solar_raw),
      "Velocidad Viento (m/s)"         = !!sym(col_viento_raw),
      "Intensidad Tr\u00e1fico"        = intensidad_raw,
      "Carga Tr\u00e1fico"             = carga_raw,
      "Precipitaciones (mm)"           = Precipitaciones_raw
    ) |>
    drop_na()

  n <- nrow(df)
  cat("\n[", toupper(nivel), "] Observaciones:", n,
      "| Estaciones:", n_distinct(df$ESTACION_LABEL), "\n")
  if (n < 10) { warning("Muy pocos datos para el nivel ", nivel); return(invisible(NULL)) }

  # Función auxiliar local: construir matriz de correlación
  bmat <- function(d, met) {
    d |> select(NO2, log_NO2, all_of(covs_label)) |>
      rename("NO2" = NO2, "log(NO2)" = log_NO2) |>
      cor(method = met, use = "complete.obs")
  }

  nivel_lbl <- switch(nivel,
    diario  = "Nivel Diario 2025",
    mensual = "Nivel Mensual 2025"
  )

  # ── Heatmaps Pearson + Spearman (global y por estación) ─────────────────────
  for (met in c("pearson", "spearman")) {
    dir_m  <- if (met == "pearson") DIR_P else DIR_S
    met_lb <- if (met == "pearson") "Pearson" else "Spearman"

    # Global
    guardar_heatmap(
      mat = bmat(df, met), metodo = met,
      titulo = paste0("Correlaci\u00f3n de ", met_lb, " \u2014 ", nivel_lbl,
                      "\nTodas las estaciones seleccionadas"),
      dir_out = dir_m,
      nombre_archivo = paste0("heatmap_", met, "_global_", nivel, ".png")
    )

    # Por estación
    for (est in unique(df$ESTACION_LABEL)) {
      df_e <- df[df$ESTACION_LABEL == est, ]
      if (nrow(df_e) < 10) next
      guardar_heatmap(
        mat = bmat(df_e, met), metodo = met,
        titulo = paste0("Correlaci\u00f3n de ", met_lb, " \u2014 ", nivel_lbl, "\n", est),
        dir_out = dir_m,
        nombre_archivo = paste0("heatmap_", met, "_",
                                gsub("[^[:alnum:]]", "_", est), "_", nivel, ".png")
      )
    }
    cat("  Heatmaps", met_lb, "guardados.\n")
  }

  # ── Scatter plots (NO2 y log NO2 vs cada covariable) ────────────────────────
  for (cov in covs_label) {
    df_long <- df |>
      select(ESTACION_LABEL, all_of(cov), NO2, log_NO2) |>
      pivot_longer(c(NO2, log_NO2), names_to = "tipo_NO2", values_to = "valor_NO2") |>
      mutate(tipo_NO2 = if_else(tipo_NO2 == "NO2",
                                "NO\u2082 (\u03bcg/m\u00b3)", "log(NO\u2082)"))

    p_sc <- ggplot(df_long, aes(x = .data[[cov]], y = valor_NO2)) +
      geom_point(alpha = 0.15, size = 0.8) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 0.9, color = "#D73027") +
      facet_grid(tipo_NO2 ~ ESTACION_LABEL, scales = "free_y") +
      labs(
        title   = paste0("Relaci\u00f3n entre ", cov, " y NO\u2082 \u2014 ", nivel_lbl),
        x = cov, y = NULL,
        caption = "L\u00ednea roja: suavizado LOESS con banda de confianza al 95%"
      ) +
      theme_bw(base_size = 9) +
      theme(
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 10),
        strip.text   = element_text(size = 7.5),
        axis.text.x  = element_text(angle = 30, hjust = 1, size = 7),
        plot.caption = element_text(size = 7, hjust = 1)
      )
    ggsave(file.path(DIR_SC,
                     paste0("scatter_", gsub("[^[:alnum:]]", "_", cov), "_", nivel, ".png")),
           p_sc, width = 18, height = 6, dpi = 160)
  }
  cat("  Scatter plots guardados.\n")

  # ── Comparativa entre estaciones (gráfico de líneas) ────────────────────────
  cor_cmp <- expand.grid(
    ESTACION_LABEL = unique(df$ESTACION_LABEL),
    Covariable     = covs_label,
    Metodo         = c("pearson", "spearman"),
    Respuesta      = c("NO2", "log_NO2"),
    stringsAsFactors = FALSE
  ) |>
    rowwise() |>
    mutate(cor_val = {
      ea  <- ESTACION_LABEL; cov_ <- Covariable
      met <- Metodo;         res  <- Respuesta
      d_  <- df[df$ESTACION_LABEL == ea, ]
      x   <- d_[[cov_]]
      y   <- if (res == "NO2") d_$NO2 else d_$log_NO2
      suppressWarnings(cor(x, y, method = met, use = "complete.obs"))
    }) |>
    ungroup() |>
    mutate(
      Metodo_label    = if_else(Metodo == "pearson", "Pearson", "Spearman"),
      Respuesta_label = if_else(Respuesta == "NO2", "NO\u2082 (\u03bcg/m\u00b3)", "log(NO\u2082)"),
      Covariable      = factor(Covariable, levels = covs_label)
    )

  p_cmp <- cor_cmp |>
    mutate(
      Metodo_label    = factor(Metodo_label,    levels = c("Pearson", "Spearman")),
      Respuesta_label = factor(Respuesta_label, levels = c("NO\u2082 (\u03bcg/m\u00b3)", "log(NO\u2082)"))
    ) |>
    ggplot(aes(x = Covariable, y = cor_val,
               color = ESTACION_LABEL, group = ESTACION_LABEL)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    facet_grid(Metodo_label ~ Respuesta_label) +
    scale_color_manual(values = paleta_estaciones, name = "Estaci\u00f3n") +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    labs(
      title   = paste0("Comparativa de correlaci\u00f3n entre estaciones \u2014 ", nivel_lbl),
      x = NULL, y = "Coeficiente de correlaci\u00f3n",
      caption = "Pearson (lineal) vs Spearman (rango); NO\u2082 original y logar\u00edtmica"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 7.5),
      legend.position  = "right",
      strip.text       = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 7, hjust = 1)
    )
  ggsave(file.path(DIR_CMP, paste0("comparativa_panel_2x2_", nivel, ".png")),
         p_cmp, width = 14, height = 9, dpi = 180)

  # Gráficos individuales por método × respuesta
  for (met in c("pearson", "spearman")) {
    for (res in c("NO2", "log_NO2")) {
      met_lb  <- if (met == "pearson") "Pearson" else "Spearman"
      res_lb  <- if (res == "NO2") "NO2" else "logNO2"
      res_tit <- if (res == "NO2") "NO\u2082" else "log(NO\u2082)"
      p_ind <- cor_cmp |>
        filter(Metodo == met, Respuesta == res) |>
        ggplot(aes(x = Covariable, y = cor_val,
                   color = ESTACION_LABEL, group = ESTACION_LABEL)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2.2) +
        scale_color_manual(values = paleta_estaciones, name = "Estaci\u00f3n") +
        scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
        labs(
          title = paste0("Comparativa ", met_lb, " vs ", res_tit, " \u2014 ", nivel_lbl),
          x = NULL, y = "Coeficiente de correlaci\u00f3n"
        ) +
        theme_bw(base_size = 10) +
        theme(
          plot.title       = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
          legend.position  = "right",
          panel.grid.minor = element_blank()
        )
      ggsave(file.path(DIR_CMP,
                       paste0("comparativa_", met, "_", res_lb, "_", nivel, ".png")),
             p_ind, width = 11, height = 6, dpi = 180)
    }
  }
  cat("  Comparativa entre estaciones guardada.\n")

  # ── Estratificación estacional ───────────────────────────────────────────────
  df_seas <- df |>
    mutate(
      Mes_ = as.integer(format(FECHA_, "%m")),
      Estacion_Anno = case_when(
        Mes_ %in% c(12L, 1L, 2L) ~ "Invierno",
        Mes_ %in% c(3L, 4L, 5L)  ~ "Primavera",
        Mes_ %in% c(6L, 7L, 8L)  ~ "Verano",
        Mes_ %in% c(9L, 10L, 11L)~ "Oto\u00f1o"
      ),
      Estacion_Anno = factor(Estacion_Anno,
                             levels = c("Invierno", "Primavera", "Verano", "Oto\u00f1o"))
    )

  colores_seas <- c(Invierno = "#1F618D", Primavera = "#1E8449",
                    Verano = "#B7770D", "Oto\u00f1o" = "#784212")

  # Heatmaps globales y por estación
  for (estacion_anno in levels(df_seas$Estacion_Anno)) {
    df_t <- df_seas[df_seas$Estacion_Anno == estacion_anno, ]
    if (nrow(df_t) < 10) next

    for (met in c("pearson", "spearman")) {
      dir_m <- if (met == "pearson") DIR_EP else DIR_ES
      met_lb <- if (met == "pearson") "Pearson" else "Spearman"

      # Global estacional
      guardar_heatmap(
        mat = bmat(df_t, met), metodo = met,
        titulo = paste0("Correlaci\u00f3n de ", met_lb, " \u2014 ",
                        estacion_anno, " | ", nivel_lbl,
                        "\nTodas las estaciones seleccionadas"),
        dir_out = dir_m,
        nombre_archivo = paste0("heatmap_", met, "_",
                                tolower(gsub("\u00f1", "n", estacion_anno)),
                                "_global_", nivel, ".png")
      )

      # Por estación de medición
      for (est in unique(df_t$ESTACION_LABEL)) {
        df_te <- df_t[df_t$ESTACION_LABEL == est, ]
        if (nrow(df_te) < 5) next
        guardar_heatmap(
          mat = bmat(df_te, met), metodo = met,
          titulo = paste0("Correlaci\u00f3n de ", met_lb, " \u2014 ",
                          estacion_anno, " | ", nivel_lbl, "\n", est),
          dir_out = dir_m,
          nombre_archivo = paste0("heatmap_", met, "_",
                                  tolower(gsub("\u00f1", "n", estacion_anno)), "_",
                                  gsub("[^[:alnum:]]", "_", est), "_", nivel, ".png")
        )
      }
    }
  }

  # Panel comparativo estacional (líneas)
  cor_seas <- expand.grid(
    Estacion_Anno = levels(df_seas$Estacion_Anno),
    Covariable    = covs_label,
    Metodo        = c("pearson", "spearman"),
    Respuesta     = c("NO2", "log_NO2"),
    stringsAsFactors = FALSE
  ) |>
    rowwise() |>
    mutate(cor_val = {
      ea  <- Estacion_Anno; cov_ <- Covariable
      met <- Metodo;        res  <- Respuesta
      d_  <- df_seas[df_seas$Estacion_Anno == ea, ]
      x   <- d_[[cov_]]
      y   <- if (res == "NO2") d_$NO2 else d_$log_NO2
      suppressWarnings(cor(x, y, method = met, use = "complete.obs"))
    }) |>
    ungroup() |>
    mutate(
      Metodo_label    = factor(if_else(Metodo == "pearson", "Pearson", "Spearman"),
                               levels = c("Pearson", "Spearman")),
      Respuesta_label = factor(if_else(Respuesta == "NO2",
                                       "NO\u2082 (\u03bcg/m\u00b3)", "log(NO\u2082)"),
                               levels = c("NO\u2082 (\u03bcg/m\u00b3)", "log(NO\u2082)")),
      Estacion_Anno   = factor(Estacion_Anno,
                               levels = c("Invierno", "Primavera", "Verano", "Oto\u00f1o")),
      Covariable      = factor(Covariable, levels = covs_label)
    )

  p_seas <- ggplot(cor_seas,
                   aes(x = Covariable, y = cor_val,
                       color = Estacion_Anno, group = Estacion_Anno)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_grid(Metodo_label ~ Respuesta_label) +
    scale_color_manual(values = colores_seas, name = "Estaci\u00f3n del a\u00f1o") +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    labs(
      title   = paste0("Comparativa estacional de correlaci\u00f3n \u2014 ", nivel_lbl),
      x = NULL, y = "Coeficiente de correlaci\u00f3n",
      caption = "Invierno: dic-feb  \u2022  Primavera: mar-may  \u2022  Verano: jun-ago  \u2022  Oto\u00f1o: sep-nov"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
      legend.position  = "right",
      legend.title     = element_text(size = 10, face = "bold"),
      strip.text       = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 7.5, hjust = 1, color = "grey40")
    )
  ggsave(file.path(DIR_EST, paste0("panel_comparativo_estacional_", nivel, ".png")),
         p_seas, width = 14, height = 9, dpi = 180)

  # Grid 1×4 cowplot
  if (requireNamespace("cowplot", quietly = TRUE)) {
    library(cowplot)
    for (met in c("pearson", "spearman")) {
      met_lb <- if (met == "pearson") "Pearson" else "Spearman"
      lista <- lapply(levels(df_seas$Estacion_Anno), function(ea) {
        df_t <- df_seas[df_seas$Estacion_Anno == ea, ]
        mat_ <- bmat(df_t, met)
        col_ <- colores_seas[[ea]]
        ggcorrplot(mat_, method = "square", type = "full",
                   lab = TRUE, lab_size = 2.4, tl.cex = 7.5,
                   colors = c("#4575B4", "white", "#D73027"),
                   title = ea, ggtheme = theme_minimal(base_size = 8)) +
          theme(plot.title  = element_text(hjust = 0.5, face = "bold",
                                           size = 10, color = col_),
                plot.margin = margin(4, 4, 4, 4))
      })
      p_tit <- ggdraw() +
        draw_label(paste0("Heatmaps estacionales \u2014 Correlaci\u00f3n de ", met_lb,
                          " \u2014 ", nivel_lbl),
                   fontface = "bold", size = 11, hjust = 0.5)
      p_fin <- plot_grid(p_tit,
                         plot_grid(plotlist = lista, nrow = 1, ncol = 4),
                         ncol = 1, rel_heights = c(0.06, 1))
      ggsave(file.path(DIR_EST,
                       paste0("heatmaps_1x4_", met, "_estacional_", nivel, ".png")),
             p_fin, width = 20, height = 6, dpi = 180)
    }
  }
  cat("  Estratificaci\u00f3n estacional guardada.\n")
  cat("=== NIVEL", toupper(nivel), "completado. Resultados en:", dir_nivel, "===\n")
}

# =============================================================================
# NIVEL DIARIO
# =============================================================================

cat("\n\n====================================================\n")
cat("  EJECUTANDO ANÁLISIS A NIVEL DIARIO\n")
cat("====================================================\n")

datos_d <- readRDS("data/processed/Maestro/diario/dataset_maestro_inla_2025_DIARIO.rds")

DIR_BASE_DIARIO <- file.path("outputs", "analysis", "correlacion_scatter", "diario")

run_analisis_correlacion(
  datos       = datos_d,
  col_no2     = "DATO_DIARIO",
  col_log_no2 = "LOG_NO2_DIARIO",
  col_fecha   = "FECHA",
  nivel       = "diario",
  dir_nivel   = DIR_BASE_DIARIO
)

# =============================================================================
# NIVEL MENSUAL
# =============================================================================

cat("\n\n====================================================\n")
cat("  EJECUTANDO ANÁLISIS A NIVEL MENSUAL\n")
cat("====================================================\n")

# Crear dataset maestro mensual si no existe
ruta_mensual <- "data/processed/Maestro/mensual/dataset_maestro_inla_2025_MENSUAL.rds"

if (!file.exists(ruta_mensual)) {
  cat("  Creando dataset maestro mensual...\n")
  datos_mensual_tmp <- readRDS("data/processed/Maestro/diario/dataset_maestro_inla_2025_DIARIO.rds")
  datos_m <- datos_mensual_tmp |>
    mutate(MES = format(as.Date(FECHA), "%Y-%m")) |>
    group_by(ESTACION, MES, barrio, LONGITUD, LATITUD, distrito, ID_DISTRITO) |>
    summarise(
      DATO_MENSUAL              = mean(DATO_DIARIO,              na.rm = TRUE),
      LOG_NO2_MENSUAL           = mean(LOG_NO2_DIARIO,           na.rm = TRUE),
      intensidad_raw            = mean(intensidad_raw,            na.rm = TRUE),
      carga_raw                 = mean(carga_raw,                 na.rm = TRUE),
      Temperatura_raw           = mean(Temperatura_raw,           na.rm = TRUE),
      Humedad_Relativa_raw      = mean(Humedad_Relativa_raw,      na.rm = TRUE),
      Precipitaciones_raw       = sum(Precipitaciones_raw,        na.rm = TRUE),
      `Presion Barométrica_raw` = mean(`Presion Barométrica_raw`, na.rm = TRUE),
      `Radiación Solar_raw`     = mean(`Radiación Solar_raw`,     na.rm = TRUE),
      `Velocidad Viento_raw`    = mean(`Velocidad Viento_raw`,    na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(FECHA = as.Date(paste0(MES, "-01")))
  saveRDS(datos_m, ruta_mensual)
  cat("  Dataset mensual guardado en:", ruta_mensual, "\n")
} else {
  datos_m <- readRDS(ruta_mensual)
  cat("  Dataset mensual cargado desde:", ruta_mensual, "\n")
}

DIR_BASE_MENSUAL <- file.path("outputs", "analysis", "correlacion_scatter", "mensual")

run_analisis_correlacion(
  datos       = datos_m,
  col_no2     = "DATO_MENSUAL",
  col_log_no2 = "LOG_NO2_MENSUAL",
  col_fecha   = "FECHA",
  nivel       = "mensual",
  dir_nivel   = DIR_BASE_MENSUAL
)

cat("\n=== ANÁLISIS COMPLETO (horario + diario + mensual) FINALIZADO ===\n")

