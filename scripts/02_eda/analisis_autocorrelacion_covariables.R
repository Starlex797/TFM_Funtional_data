# ==============================================================================
# ANÁLISIS DE AUTOCORRELACIÓN: COVARIABLES vs NO2 (DATOS RAW, SIN LOG)
# ==============================================================================
# Objetivo: Entender la estructura temporal de las correlaciones entre
# covariables clima/tráfico y NO2 raw (sin transformación logarítmica)
# ==============================================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

set.seed(4827)

# ==============================================================================
# 1. CARGAR DATOS
# ==============================================================================

cat("\n1. Cargando datos...\n")

# Datos horarios
dt_horario <- readRDS(here("data", "processed", "Maestro", "horario", "dataset_maestro_inla_2025_HORARIO.rds"))
setDT(dt_horario)

# Datos diarios
dt_diario <- readRDS(here("data", "processed", "Maestro", "diario", "dataset_maestro_inla_2025_DIARIO.rds"))
setDT(dt_diario)

# Datos mensuales (construir a partir de diarios)
dt_diario[, ANIO_MES := as.Date(format(FECHA, "%Y-%m-01"))]
dt_mensual <- dt_diario[, .(
  DATO_NO2                  = mean(DATO_DIARIO, na.rm = TRUE),
  intensidad_raw            = mean(intensidad_raw, na.rm = TRUE),
  Temperatura_raw           = mean(Temperatura_raw, na.rm = TRUE),
  Precipitaciones_raw       = sum(Precipitaciones_raw, na.rm = TRUE),
  `Presion Barométrica_raw` = mean(`Presion Barométrica_raw`, na.rm = TRUE),
  `Velocidad Viento_raw`    = mean(`Velocidad Viento_raw`, na.rm = TRUE),
  n_dias                    = .N
), by = .(ESTACION, ANIO_MES, barrio, distrito, LONGITUD, LATITUD, ID_DISTRITO)]
setnames(dt_mensual, "ANIO_MES", "FECHA")

cat("   Horario:  ", nrow(dt_horario), "registros\n")
cat("   Diario:   ", nrow(dt_diario), "registros\n")
cat("   Mensual:  ", nrow(dt_mensual), "registros\n")

# ==============================================================================
# 2. ANÁLISIS POR FRECUENCIA: HORARIA
# ==============================================================================

cat("\n2. FRECUENCIA HORARIA\n")
cat("   Período: Noviembre (5 días)\n")

todas_fechas_h <- sort(unique(dt_horario$FECHA))
fechas_mes_h <- todas_fechas_h[format(todas_fechas_h, "%m") == "11"]
fechas_sel_h <- fechas_mes_h[13:17]
dt_h_subset <- dt_horario[FECHA %in% fechas_sel_h]
dt_h_subset[, ID_TIEMPO := as.integer(factor(paste(FECHA, sprintf("%02d", HORA))))]

# Crear carpeta de salida
carpeta_h <- here("outputs", "analysis", "analisis_autocorrelacion", "horario")
dir.create(carpeta_h, showWarnings = FALSE, recursive = TRUE)

cols_raw_h <- c("intensidad_raw", "Temperatura_raw", "Precipitaciones_raw",
                "Presion Barométrica_raw", "Velocidad Viento_raw")
cols_raw_h_disp <- intersect(cols_raw_h, names(dt_h_subset))

etiquetas_h <- c(
  intensidad_raw = "Intensidad tráfico",
  Temperatura_raw = "Temperatura",
  Precipitaciones_raw = "Precipitaciones",
  `Presion Barométrica_raw` = "Presión barométrica",
  `Velocidad Viento_raw` = "Velocidad viento"
)

# Correlaciones globales
mat_h <- dt_h_subset[, c(..cols_raw_h_disp)]
mat_h[, NO2_raw := dt_h_subset$DATO_HORARIO]
cor_h <- cor(mat_h, use = "pairwise.complete.obs")

cat("\n   Matriz de correlaciones (Pearson):\n")
print(round(cor_h[, "NO2_raw"], 4))

# Correlaciones por hora del día (estacionalidad intradía)
cor_por_hora <- data.frame()
for (h in sort(unique(dt_h_subset$HORA))) {
  dt_h_hora <- dt_h_subset[HORA == h]
  if (nrow(dt_h_hora) >= 2) {
    for (cov in cols_raw_h_disp) {
      cor_val <- cor(dt_h_hora[[cov]], dt_h_hora$DATO_HORARIO, use = "complete.obs")
      if (!is.na(cor_val)) {
        cor_por_hora <- rbind(cor_por_hora, data.frame(
          Hora = h,
          Covariable = etiquetas_h[cov],
          Correlacion = cor_val,
          N = nrow(dt_h_hora)
        ))
      }
    }
  }
}

# Gráfico: correlaciones por hora del día
p_h <- ggplot(cor_por_hora, aes(x = Hora, y = Correlacion, color = Covariable, group = Covariable)) +
  geom_line(alpha = 0.7, linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Covariable, scales = "free_y", nrow = 2) +
  labs(
    title = "Autocorrelación Horaria: Covariables vs NO2 (raw)",
    subtitle = "Correlación de Pearson por hora del día (Noviembre 5 días)",
    x = "Hora del día", y = "Correlación con NO2 raw"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(carpeta_h, "correlaciones_por_hora.png"),
       width = 12, height = 8, dpi = 300)
cat("   Gráfico guardado: correlaciones_por_hora.png\n")

# ==============================================================================
# 3. ANÁLISIS POR FRECUENCIA: DIARIA
# ==============================================================================

cat("\n3. FRECUENCIA DIARIA\n")
cat("   Período: Noviembre (19 días)\n")

todas_fechas_d <- sort(unique(dt_diario$FECHA))
fechas_mes_d <- todas_fechas_d[format(todas_fechas_d, "%m") == "11"]
fechas_sel_d <- fechas_mes_d[1:19]
dt_d_subset <- dt_diario[FECHA %in% fechas_sel_d]
dt_d_subset[, ID_TIEMPO := match(FECHA, fechas_sel_d)]

carpeta_d <- here("outputs", "analysis", "analisis_autocorrelacion", "diario")
dir.create(carpeta_d, showWarnings = FALSE, recursive = TRUE)

cols_raw_d <- c("intensidad_raw", "Temperatura_raw", "Precipitaciones_raw",
                "Presion Barométrica_raw", "Velocidad Viento_raw")
cols_raw_d_disp <- intersect(cols_raw_d, names(dt_d_subset))

etiquetas_d <- c(
  intensidad_raw = "Intensidad tráfico",
  Temperatura_raw = "Temperatura",
  Precipitaciones_raw = "Precipitaciones",
  `Presion Barométrica_raw` = "Presión barométrica",
  `Velocidad Viento_raw` = "Velocidad viento"
)

# Correlaciones globales
mat_d <- dt_d_subset[, c(..cols_raw_d_disp)]
mat_d[, NO2_raw := dt_d_subset$DATO_DIARIO]
cor_d <- cor(mat_d, use = "pairwise.complete.obs")

cat("\n   Matriz de correlaciones (Pearson):\n")
print(round(cor_d[, "NO2_raw"], 4))

# Correlaciones por día
cor_por_dia <- data.frame()
for (dia in sort(unique(dt_d_subset$ID_TIEMPO))) {
  dt_d_dia <- dt_d_subset[ID_TIEMPO == dia]
  if (nrow(dt_d_dia) >= 2) {
    for (cov in cols_raw_d_disp) {
      cor_val <- cor(dt_d_dia[[cov]], dt_d_dia$DATO_DIARIO, use = "complete.obs")
      if (!is.na(cor_val)) {
        cor_por_dia <- rbind(cor_por_dia, data.frame(
          Dia = dia,
          Covariable = etiquetas_d[cov],
          Correlacion = cor_val,
          N = nrow(dt_d_dia)
        ))
      }
    }
  }
}

# Gráfico: correlaciones por día
p_d <- ggplot(cor_por_dia, aes(x = Dia, y = Correlacion, color = Covariable, group = Covariable)) +
  geom_line(alpha = 0.7, linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Covariable, scales = "free_y", nrow = 2) +
  labs(
    title = "Autocorrelación Diaria: Covariables vs NO2 (raw)",
    subtitle = "Correlación de Pearson por día (Noviembre 19 días)",
    x = "Día dentro del período", y = "Correlación con NO2 raw"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(carpeta_d, "correlaciones_por_dia.png"),
       width = 12, height = 8, dpi = 300)
cat("   Gráfico guardado: correlaciones_por_dia.png\n")

# ==============================================================================
# 4. ANÁLISIS POR FRECUENCIA: MENSUAL
# ==============================================================================

cat("\n4. FRECUENCIA MENSUAL\n")
cat("   Período: Año completo\n")

carpeta_m <- here("outputs", "analysis", "analisis_autocorrelacion", "mensual")
dir.create(carpeta_m, showWarnings = FALSE, recursive = TRUE)

cols_raw_m <- c("intensidad_raw", "Temperatura_raw", "Precipitaciones_raw",
                "Presion Barométrica_raw", "Velocidad Viento_raw")
cols_raw_m_disp <- intersect(cols_raw_m, names(dt_mensual))

etiquetas_m <- c(
  intensidad_raw = "Intensidad tráfico",
  Temperatura_raw = "Temperatura",
  Precipitaciones_raw = "Precipitaciones",
  `Presion Barométrica_raw` = "Presión barométrica",
  `Velocidad Viento_raw` = "Velocidad viento"
)

# Correlaciones globales
mat_m <- dt_mensual[, c(..cols_raw_m_disp)]
mat_m[, NO2_raw := dt_mensual$DATO_NO2]
cor_m <- cor(mat_m, use = "pairwise.complete.obs")

cat("\n   Matriz de correlaciones (Pearson):\n")
print(round(cor_m[, "NO2_raw"], 4))

# Agregar mes para análisis de estacionalidad
dt_mensual[, Mes := month(FECHA)]
dt_mensual[, Mes_nombre := c("Ene", "Feb", "Mar", "Abr", "May", "Jun",
                             "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")[Mes]]

cor_por_mes <- data.frame()
for (mes in sort(unique(dt_mensual$Mes))) {
  dt_m_mes <- dt_mensual[Mes == mes]
  if (nrow(dt_m_mes) >= 2) {
    for (cov in cols_raw_m_disp) {
      cor_val <- cor(dt_m_mes[[cov]], dt_m_mes$DATO_NO2, use = "complete.obs")
      if (!is.na(cor_val)) {
        mes_nom <- c("Ene", "Feb", "Mar", "Abr", "May", "Jun",
                     "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")[mes]
        cor_por_mes <- rbind(cor_por_mes, data.frame(
          Mes = mes,
          Mes_nombre = mes_nom,
          Covariable = etiquetas_m[cov],
          Correlacion = cor_val,
          N = nrow(dt_m_mes)
        ))
      }
    }
  }
}

# Gráfico: correlaciones por mes (estacionalidad)
if (nrow(cor_por_mes) > 0) {
  p_m <- ggplot(cor_por_mes, aes(x = factor(Mes_nombre, levels = c("Ene", "Feb", "Mar", "Abr", "May", "Jun",
                                                                       "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")),
                                  y = Correlacion, color = Covariable, group = Covariable)) +
    geom_line(alpha = 0.7, linewidth = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~Covariable, scales = "free_y", nrow = 2) +
    labs(
      title = "Autocorrelación Mensual: Covariables vs NO2 (raw) — Ciclo Estacional",
      subtitle = "Correlación de Pearson por mes del año (datos promediados mensualmente)",
      x = "Mes", y = "Correlación con NO2 raw"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(carpeta_m, "correlaciones_por_mes_estacionalidad.png"),
         width = 12, height = 8, dpi = 300)
  cat("   Gráfico guardado: correlaciones_por_mes_estacionalidad.png\n")
}

# ==============================================================================
# 5. RESUMEN COMPARATIVO
# ==============================================================================

cat("\n5. RESUMEN COMPARATIVO\n")
cat("\n   HORARIO (global):\n")
print(round(cor_h[, "NO2_raw"], 4))

cat("\n   DIARIO (global):\n")
print(round(cor_d[, "NO2_raw"], 4))

cat("\n   MENSUAL (global):\n")
print(round(cor_m[, "NO2_raw"], 4))

cat("\n\n   Interpretación:\n")
cat("   - Variabilidad en correlaciones por período sugiere AUTOCORRELACIÓN TEMPORAL\n")
cat("   - Correlaciones que cambian de signo indican CONFUSIÓN ESTACIONAL\n")
cat("   - Patrón estacional en correlaciones mensuales es clave para validar AR1\n")

cat("\n   Análisis completado. Resultados en: ", here("outputs", "analysis", "analisis_autocorrelacion"), "\n\n")
