# ==============================================================================
# EDA 2025
#===============================================================================
# Structure of the analysis:
#   1. Variogramas experimentales (variables con n ≥ 20 estaciones)
#   2. Índice de Moran Global: MC (999 perm.) + analítico (todas las variables)
#   3. Análisis de robustez: k-vecinos vs distancia euclidiana (5 km)
#   4. Guardar todos los resultados
# ==============================================================================


# ==============================================================================
# Library imports
# ==============================================================================
library(data.table)
library(sf)
library(gstat)
library(ggplot2)
library(ggrepel)
library(spdep)
library(dplyr)
library(gridExtra)
library(grid)
library(here)

# Hyperparameters for Moran's I and variogram calculations
set.seed(4712)
NSIM      <- 999   # Montecarlo Permutations for Moran's I
CUTOFF_KM <- 20    # distance max for  variograma (km)
WIDTH_KM  <- 4     # width for variograma (km)
N_MIN_VAR <- 20    # min stations to calculate variogram

carpeta_correlacion <- here("outputs", "correlacion_espacial")

dir.create(carpeta_correlacion, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. CARGA Y PREPARACIÓN DE DATOS
# ==============================================================================

# NO2 diario → media anual por estación
dt_aire_raw  <- readRDS(here("data", "processed", "contaminacion", "mensual",
                             "aire_madrid_2025_log_No2_mensuales.rds"))
dt_aire_diario<- readRDS(here("data", "processed", "Contaminacion", "diario",
                             "aire_madrid_2025_No2_trans_diarios.rds"))
prueba<- View(dt_aire_raw)
dt_no2_anual <- dt_aire_raw[!is.na(DATO_MENSUAL), .(
  valor    = mean(DATO_MENSUAL, na.rm = TRUE),
  LONGITUD = unique(LONGITUD),
  LATITUD  = unique(LATITUD)
), by = ESTACION]

# Meteorología mensual → media anual por estación
dt_meteo <- readRDS(here("data", "processed", "Clima", "mensual",
                         "meteo_madrid_2025_mensual.rds"))

# Tráfico mensual por distrito (para scatter NO₂ vs tráfico)
dt_traf_mensual_dist <- readRDS(here("data", "processed",
                                     "trafico_madrid_2025_mensual_distrito.rds"))

vars_clima   <- c("Temperatura", "Humedad_Relativa","Precipitaciones", "Presion Barométrica","Radiación Solar","Velocidad Viento")
labels_clima <- c("Temperatura (°C)", "Humedad Relativa (%)", "Precipitaciones (mm)", "Presión Barométrica (hPa)","Radiación Solar (W/h)","Velocidad Viento (m/s)")

# We are calculating the mean annual value for each climate variable by station, excluding NA values. 
# The result is stored in a list where each element corresponds to a climate variable and contains a data.table 
# with the mean value, longitude, and latitude for each station.
lista_clima <- setNames(
  lapply(vars_clima, function(v) {
    dt_meteo[!is.na(get(v)), .(
      valor    = mean(get(v), na.rm = TRUE),
      LONGITUD = unique(LONGITUD),
      LATITUD  = unique(LATITUD)
    ), by = ESTACION]
  }),
  labels_clima
)

# Dataset unificado: NO2 + todas las variables climáticas
todos_datos <- c(list("NO2 (µg/m³)" = dt_no2_anual), lista_clima)
print(head(todos_datos))
cat("Número de estaciones por variable:\n")
# Number of stations per climate variable 
for (nm in names(todos_datos)) {
  elegible_var <- if (nrow(todos_datos[[nm]]) >= N_MIN_VAR) "  ← variograma" else ""
  cat(sprintf("  %-30s n = %2d%s\n", nm, nrow(todos_datos[[nm]]), elegible_var))
}


#===============================================================================
# Block 1: Studying the temporal evolution of climate variables and the spatial variability between stations
# and contamination levels, using line plots and coefficient of variation (CV) between stations for each month.
# ==============================================================================
# Montly assessment of climate variables (line plots + CV between stations)
# ==============================================================================
# Line plots are faceted by variable, and a separate plot shows the coefficient of variation (CV)
# between stations for each month.
# Goal: to assess the temporal evolution of climate variables and the spatial variability between stations.
# ==============================================================================

dt_largo <- melt(
  dt_meteo,
  id.vars    = c("ESTACION", "MES"),
  measure.vars = vars_clima,
  variable.name = "Variable",
  value.name    = "Valor"
)

setDT(dt_largo)
etiquetas_lineas <- setNames(labels_clima, vars_clima)
dt_largo[, Variable_label := etiquetas_lineas[as.character(Variable)]]
dt_largo[, Mes_num := as.integer(sub("2025-", "", MES))]

meses_lab <- c("Ene","Feb","Mar","Abr","May","Jun",
               "Jul","Ago","Sep","Oct","Nov","Dic")

plots_lineas <- list()

for (var_label in unique(dt_largo$Variable_label)) {

  dt_sub    <- dt_largo[Variable_label == var_label & !is.na(Valor)]
  n_est_var <- uniqueN(dt_sub$ESTACION)

  p <- ggplot(dt_sub, aes(x = Mes_num, y = Valor,
                           color = ESTACION, group = ESTACION)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_x_continuous(breaks = 1:12, labels = meses_lab) +
    labs(
      title    = paste("Evolución Mensual —", var_label),
      subtitle = paste0(n_est_var, " estaciones · Madrid 2025"),
      x        = "Mes",
      y        = var_label,
      color    = "Estación"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "gray40", size = 10),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.4, "cm")
    )

  plots_lineas[[var_label]] <- p
}

# Coeficiente de variación entre estaciones por mes
dt_cv <- dt_largo[!is.na(Valor), .(
  media = mean(Valor),
  sd    = sd(Valor),
  n_est = uniqueN(ESTACION)
), by = .(Variable_label, Mes_num)]
dt_cv[, CV := (sd / abs(media)) * 100]

plot_cv_facet <- ggplot(dt_cv, aes(x = Mes_num, y = CV)) +
  geom_line(linewidth = 0.7, color = "#2c3e50") +
  geom_point(size = 2, color = "#2c3e50") +
  scale_x_continuous(breaks = 1:12, labels = meses_lab) +
  facet_wrap(~ Variable_label, scales = "free_y", ncol = 2) +
  labs(
    title    = "Coeficiente de Variación entre Estaciones por Mes",
    subtitle = "CV (%) = sd / |media| × 100 · Madrid 2025",
    x        = "Mes",
    y        = "CV (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "gray40", size = 10),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = "gray85"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(plot_cv_facet)

# ==============================================================================
# GUARDAR GRÁFICOS — Evolución mensual por variable + CV entre estaciones
# ==============================================================================

# Función auxiliar: convierte el label de una variable en nombre de fichero válido
sanitizar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")   # quita tildes / caracteres especiales
  x <- gsub("[^[:alnum:]_]", "_", x)      # reemplaza todo lo que no sea alfanum.
  x <- gsub("_+", "_", x)                 # colapsa guiones bajos múltiples
  x <- gsub("^_|_$", "", x)              # elimina guiones al inicio/final
  tolower(x)
}

# 1. Una PNG por variable climática (evolución mensual por estación)
for (nm in names(plots_lineas)) {
  nombre_archivo <- paste0("evolucion_mensual_", sanitizar_nombre(nm), ".png")
  ggsave(
    filename = file.path(carpeta_correlacion, nombre_archivo),
    plot     = plots_lineas[[nm]],
    width    = 10, height = 5.5, dpi = 150
  )
}

# 2. Gráfico del Coeficiente de Variación (facet con todas las variables)
ggsave(
  filename = file.path(carpeta_correlacion, "cv_entre_estaciones_por_mes.png"),
  plot     = plot_cv_facet,
  width    = 10, height = 8, dpi = 150
)

cat("✅ Gráficos de líneas y CV guardados en:\n  ", carpeta_correlacion, "\n")

# ==============================================================================
# 2c. PERFILES TEMPORALES DE CONTAMINACIÓN — NO₂ A TRES ESCALAS
# ==============================================================================
# Objetivo: caracterizar el comportamiento del NO₂ en Madrid 2025 a escala
#   mensual  → ciclo estacional (calefacción, inversión térmica)
#   semanal  → laborable vs fin de semana (peso del tráfico)
#   horaria  → patrón de "doble joroba" en horas punta
# ==============================================================================

# ── 2c.1 Escala mensual ────────────────────────────────────────────────────────
dt_diario_no2 <- dt_aire_diario[!is.na(DATO_DIARIO)]
dt_diario_no2[, Mes_num := as.integer(format(FECHA, "%m"))]

plot_perfil_mensual <- ggplot(dt_diario_no2,
                              aes(x = factor(Mes_num), y = DATO_DIARIO)) +
  geom_boxplot(fill = "#3498db", color = "#2c3e50", alpha = 0.7,
               outlier.size = 0.8, outlier.alpha = 0.4) +
  scale_x_discrete(labels = meses_lab) +
  labs(
    title    = "Perfil Mensual de NO\u2082 \u2014 Ciclo Estacional",
    subtitle = "Distribuci\u00f3n de medias diarias (todas las estaciones) \u00b7 Madrid 2025",
    x        = "Mes",
    y        = "NO\u2082 (\u00b5g/m\u00b3)",
    caption  = "Cada observaci\u00f3n es la media diaria de una estaci\u00f3n"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "gray40", size = 10),
    panel.grid.minor = element_blank()
  )

print(plot_perfil_mensual)

# ── 2c.2 Escala semanal ────────────────────────────────────────────────────────
# format "%u" → 1 = lunes, …, 7 = domingo  (independiente del locale)
dias_lab <- c("Lun", "Mar", "Mi\u00e9", "Jue", "Vie", "S\u00e1b", "Dom")

dt_diario_no2[, Dia_num  := as.integer(format(FECHA, "%u"))]
dt_diario_no2[, Tipo_dia := fifelse(Dia_num >= 6L, "Fin de semana", "Laborable")]
dt_diario_no2[, Dia_fac  := factor(Dia_num, levels = 1:7, labels = dias_lab)]

plot_perfil_semanal <- ggplot(dt_diario_no2,
                              aes(x = Dia_fac, y = DATO_DIARIO, fill = Tipo_dia)) +
  geom_boxplot(color = "#2c3e50", alpha = 0.75,
               outlier.size = 0.8, outlier.alpha = 0.4) +
  scale_fill_manual(
    values = c("Laborable" = "#2980b9", "Fin de semana" = "#e67e22"),
    name   = NULL
  ) +
  labs(
    title    = "Perfil Semanal de NO\u2082 \u2014 Laborable vs Fin de Semana",
    subtitle = "Peso del tr\u00e1fico laboral en la contaminaci\u00f3n \u00b7 Madrid 2025",
    x        = "D\u00eda de la semana",
    y        = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "gray40", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "top"
  )

print(plot_perfil_semanal)

# ── 2c.3 Escala horaria ────────────────────────────────────────────────────────
dt_horario_no2 <- readRDS(here("data", "processed", "contaminacion", "horario",
                               "aire_madrid_2025_No2_horarios.rds"))

# HORA es factor "H01"…"H24" → entero 1…24
dt_horario_no2[, Hora_num := as.integer(sub("H0?", "", as.character(HORA)))]
dt_horario_no2[, Dia_num  := as.integer(format(FECHA, "%u"))]
dt_horario_no2[, Tipo_dia := fifelse(Dia_num >= 6L, "Fin de semana", "Laborable")]

# Mediana y rango intercuartílico por hora y tipo de día
dt_hora_resumen <- dt_horario_no2[!is.na(DATO), .(
  mediana = median(DATO),
  p25     = quantile(DATO, 0.25),
  p75     = quantile(DATO, 0.75)
), by = .(Hora_num, Tipo_dia)]

plot_perfil_horario <- ggplot(dt_hora_resumen,
                              aes(x = Hora_num, y = mediana,
                                  color = Tipo_dia, fill = Tipo_dia)) +
  geom_ribbon(aes(ymin = p25, ymax = p75), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(2, 24, by = 2),
                     labels = function(x) sprintf("%02d:00", x)) +
  scale_color_manual(
    values = c("Laborable" = "#2980b9", "Fin de semana" = "#e67e22"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = c("Laborable" = "#2980b9", "Fin de semana" = "#e67e22"),
    name   = NULL
  ) +
  labs(
    title    = "Perfil Horario de NO\u2082 \u2014 Patr\u00f3n de Doble Joroba",
    subtitle = "Mediana \u00b1 IQR por hora del d\u00eda \u00b7 Laborable vs Fin de Semana \u00b7 Madrid 2025",
    x        = "Hora del d\u00eda",
    y        = "NO\u2082 (\u00b5g/m\u00b3)",
    caption  = "Banda: rango intercuart\u00edlico (P25\u2013P75) sobre todas las estaciones"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "gray40", size = 10),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position  = "top"
  )

print(plot_perfil_horario)


# ==============================================================================
# Correlation heatmap: NO2, climate variables and traffic variables
# ==============================================================================
# Uses the daily master dataset (dataset_maestro_inla_2025_DIARIO.rds) built in
# Paso 2. Pearson correlation on raw-scale variables (result is identical to
# standardised variables because correlation is invariant to linear transforms).
# ==============================================================================



library(reshape2)

dt_maestro_eda <- readRDS(here("data", "processed",
                               "dataset_maestro_inla_2025_DIARIO.rds"))

cols_no2_hm     <- "DATO_DIARIO"
cols_trafico_hm <- c("intensidad_raw", "carga_raw")
cols_clima_hm   <- c("Temperatura_raw", "Humedad_Relativa_raw",
                      "Precipitaciones_raw", "Presion Barométrica_raw",
                      "Radiación Solar_raw", "Velocidad Viento_raw")

all_cols_hm <- c(cols_no2_hm, cols_trafico_hm, cols_clima_hm)

cor_mat_hm <- cor(dt_maestro_eda[, ..all_cols_hm], use = "pairwise.complete.obs")

labels_hm <- c("NO2", "Intensidad", "Carga",
               "Temperatura", "Humedad Rel.", "Precipitaciones",
               "Presión Barom.", "Radiación Solar", "Vel. Viento")
rownames(cor_mat_hm) <- labels_hm
colnames(cor_mat_hm) <- labels_hm

cor_long_hm <- as.data.table(melt(cor_mat_hm,
                                   varnames = c("Var1", "Var2"),
                                   value.name = "cor"))
cor_long_hm[, Var1 := factor(Var1, levels = rev(labels_hm))]
cor_long_hm[, Var2 := factor(Var2, levels = labels_hm)]

plot_heatmap_cor <- ggplot(cor_long_hm, aes(Var2, Var1, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(cor, 2)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Correlación") +
  labs(title = "Correlación: NO2, Variables Climáticas y Tráfico",
       subtitle = "Datos diarios · Madrid 2025",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_heatmap_cor)

ggsave(file.path(carpeta_correlacion, "00j_Heatmap_Correlacion_NO2_Clima_Trafico.png"),
       plot = plot_heatmap_cor, width = 10, height = 8, dpi = 300)

rm(dt_maestro_eda)





















# ==============================================================================
# 3. VARIOGRAMAS EXPERIMENTALES (n ≥ 20 estaciones)
# ==============================================================================
# Se calcula variograma experimental + ajuste de modelo teórico esférico para
# las variables con suficiente soporte espacial.
# Variables elegibles automáticamente según N_MIN_VAR.
# ==============================================================================

# ==============================================================================
# VARIOGRAMA DE RESIDUOS DEL NO2 (DÍA TÍPICO)
# ==============================================================================

library(data.table)
library(sf)
library(gstat)
library(ggplot2)
library(here)

# 1. Cargar la base de datos maestra
dt_maestro <- readRDS(here("data","processed","dataset_maestro_inla_2025_DIARIO.rds"))
setDT(dt_maestro)

# 2. Seleccionar un día específico (Corte Transversal Espacial)
# Para la estadística espacial clásica, no puedes meter todos los días a la vez.
# Elegimos una fecha representativa (ej. un día de entre semana en un mes estándar).
# Cambia esta fecha por el "Día Típico" que calculaste en el paso anterior.
dia_objetivo <- as.Date("2025-05-14") 
dt_dia <- dt_maestro[FECHA == dia_objetivo]

# Eliminamos filas con NAs en las covariables para que la regresión no falle
covs <- c("intensidad", "Temperatura", "Humedad_Relativa",
          "Precipitaciones", "Presion Barométrica", "Radiación Solar")
dt_dia <- dt_dia[complete.cases(dt_dia[, ..covs])]

# Renombramos columnas con espacios/tildes para que gstat las acepte en fórmulas
setnames(dt_dia,
         old = c("Presion Barométrica", "Radiación Solar"),
         new = c("Presion_Barometrica", "Radiacion_Solar"),
         skip_absent = TRUE)

# 3. Convertir a objeto espacial proyectado (UTM 30N, metros)
sf_dia <- st_as_sf(dt_dia,
                   coords = c("LONGITUD", "LATITUD"),
                   crs = 4326) |>
  st_transform(25830)

# 4. Calcular el Variograma Empírico de los RESIDUOS
# gstat ajusta internamente una regresión OLS sobre las covariables y calcula
# el variograma sobre los residuos (autocorrelación espacial neta).
fml_residuos <- LOG_NO2_DIARIO ~ intensidad + Temperatura + Humedad_Relativa +
  Precipitaciones + Presion_Barometrica + Radiacion_Solar

# cutoff: Distancia máxima a mirar (ej. 20 km = 20000 metros)
# width: Intervalo de búsqueda (ej. 4 km = 4000 metros para garantizar +30 parejas por punto)
var_empirico <- variogram(fml_residuos, 
                          data   = sf_dia, 
                          cutoff = 20000, 
                          width  = 4000)

# 5. Ajustar un Modelo Teórico (Curva)
# Inicializamos el modelo con unos valores a ojo:
# psill = varianza esperada, model = "Sph" (Esférico), range = rango inicial estimado, nugget = error base
# Usamos la varianza de los residuos OLS como psill inicial
resid_var <- var(residuals(lm(fml_residuos, data = sf_dia)))
modelo_inicial <- vgm(psill = resid_var,
                      model = "Sph",
                      range = 1000,
                      nugget = 0)

# El algoritmo ajusta la curva matemáticamente a los puntos experimentales
var_teorico <- fit.variogram(var_empirico, modelo_inicial)

# 6. Mostrar el Rango Espacial Real (En Kilómetros)
rango_metros <- var_teorico[var_teorico$model == "Sph", "range"]
cat("\n======================================================\n")
cat(" RANGO ESPACIAL DE LOS RESIDUOS DE NO2:", round(rango_metros / 1000, 2), "Kilómetros\n")
cat("======================================================\n")

# 7. Graficar el resultado final
plot_variograma <- plot(var_empirico, var_teorico, 
                        main = paste("Variograma de Residuos Espaciales -", dia_objetivo),
                        xlab = "Distancia (metros)", 
                        ylab = "Semivarianza")
print(plot_variograma)

# ==============================================================================
# 4. ÍNDICE DE MORAN GLOBAL — TODAS LAS VARIABLES
# ==============================================================================
# Para cada variable:
#   a) Test Monte Carlo (999 permutaciones) — robusto para n pequeño
#   b) Test analítico (asunción de normalidad) — comparación
#   c) Diagrama de dispersión de Moran
#
# Selección de k-vecinos automática:
#   n ≥ 20 → k = 4  |  n < 20 → k = 3
# ==============================================================================

moran_global_plot <- function(dt, variable_label, k = NULL, nsim = NSIM) {
  n <- nrow(dt)
  if (n < 4) {
    message("Saltando '", variable_label, "': n = ", n, " (insuficiente)")
    return(NULL)
  }
  if (is.null(k)) k <- if (n >= 20) 4L else 3L

  sf_obj <- st_as_sf(dt, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
    st_transform(crs = 25830)
  coords <- st_coordinates(sf_obj)
  nb     <- knn2nb(knearneigh(coords, k = k))
  lw     <- nb2listw(nb, style = "W")

  # Test Monte Carlo
  mc    <- moran.mc(dt$valor, lw, nsim = nsim)
  I     <- round(mc$statistic[[1]], 4)
  pv_mc <- round(mc$p.value, 4)

  # Test analítico
  ta    <- moran.test(dt$valor, lw, alternative = "greater")
  pv_an <- round(ta$p.value, 4)

  sig <- case_when(
    pv_mc < 0.01 ~ "***",
    pv_mc < 0.05 ~ "**",
    pv_mc < 0.10 ~ "*",
    TRUE         ~ "n.s."
  )

  # Scatter plot de Moran
  z       <- as.numeric(scale(dt$valor))
  Wz      <- lag.listw(lw, z)
  df_plot <- data.frame(z = z, Wz = Wz, estacion = dt$ESTACION)

  x_r <- quantile(z,  0.88); x_l <- quantile(z,  0.12)
  y_t <- quantile(Wz, 0.88); y_b <- quantile(Wz, 0.12)

  label_ann <- paste0(
    "I = ", I, "  ", sig,
    "\np(MC) = ", pv_mc,
    "\np(an) = ", pv_an,
    "\nn = ", n, "  k = ", k
  )

  p <- ggplot(df_plot, aes(x = z, y = Wz)) +
    # Fondo por cuadrantes
    annotate("rect", xmin =    0, xmax =  Inf, ymin =    0, ymax =  Inf, fill = "#c0392b", alpha = 0.07) +
    annotate("rect", xmin = -Inf, xmax =    0, ymin = -Inf, ymax =    0, fill = "#2980b9", alpha = 0.07) +
    annotate("rect", xmin =    0, xmax =  Inf, ymin = -Inf, ymax =    0, fill = "#f39c12", alpha = 0.04) +
    annotate("rect", xmin = -Inf, xmax =    0, ymin =    0, ymax =  Inf, fill = "#f39c12", alpha = 0.04) +
    # Ejes de referencia
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
    # Línea de regresión (pendiente = I de Moran)
    geom_abline(slope = I, intercept = 0, color = "#e74c3c", linewidth = 0.9) +
    # Puntos y etiquetas
    geom_point(size = 3.5, color = "#2c3e50", alpha = 0.9) +
    geom_text_repel(aes(label = estacion), size = 2.5, color = "gray30",
                    max.overlaps = 25, box.padding = 0.3) +
    # Etiquetas de cuadrante
    annotate("text", x = x_r, y = y_t, label = "HH", color = "#c0392b", size = 4, fontface = "bold", alpha = 0.8) +
    annotate("text", x = x_l, y = y_b, label = "LL", color = "#2980b9", size = 4, fontface = "bold", alpha = 0.8) +
    annotate("text", x = x_r, y = y_b, label = "HL", color = "gray45",  size = 4, fontface = "bold", alpha = 0.8) +
    annotate("text", x = x_l, y = y_t, label = "LH", color = "gray45",  size = 4, fontface = "bold", alpha = 0.8) +
    # Anotación con estadísticos
    annotate("label", x = Inf, y = Inf, hjust = 1.08, vjust = 1.3,
             label = label_ann, size = 3.8, fontface = "bold",
             fill = "white", color = "#2c3e50", label.size = 0.4) +
    labs(
      title    = paste("Diagrama de Dispersión de Moran —", variable_label),
      subtitle = paste0(n, " estaciones · k = ", k, " vecinos más cercanos · Madrid 2025"),
      x        = paste0(variable_label, "  (z-score)"),
      y        = "Rezago espacial Wz  (media ponderada de vecinos)",
      caption  = "Línea roja: pendiente = I de Moran · p(MC): 999 perm. · p(an): normalidad · UTM 30N"
    ) +
    theme_minimal() +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "gray40", size = 10),
      panel.grid.minor = element_blank()
    )

  list(plot = p, I = I, p_mc = pv_mc, p_analitico = pv_an, signif = sig, n = n, k = k)
}

# Calcular Moran para todas las variables
cat("\n--- Calculando Moran's I Global ---\n")
resultados_moran <- setNames(
  lapply(names(todos_datos), function(nm) {
    r <- moran_global_plot(todos_datos[[nm]], variable_label = nm)
    if (!is.null(r))
      cat(sprintf("  %-30s I = %7.4f  p(MC) = %.4f  p(an) = %.4f  %s\n",
                  nm, r$I, r$p_mc, r$p_analitico, r$signif))
    r
  }),
  names(todos_datos)
)

# Tabla resumen: MC + analítico + diferencia
tabla_moran <- do.call(rbind, lapply(names(resultados_moran), function(nm) {
  r <- resultados_moran[[nm]]
  if (is.null(r)) return(NULL)
  data.frame(
    Variable    = nm,
    N           = r$n,
    k           = r$k,
    I_Moran     = r$I,
    p_MC        = r$p_mc,
    p_Analitico = r$p_analitico,
    Dif_p       = round(abs(r$p_mc - r$p_analitico), 4),
    Signif      = r$signif,
    stringsAsFactors = FALSE
  )
}))

cat("\n================================================================\n")
cat("   ÍNDICE DE MORAN GLOBAL — Monte Carlo vs Analítico\n")
cat("================================================================\n")
print(tabla_moran, row.names = FALSE)
cat("\nSignif (basado en MC): *** p<0.01  ** p<0.05  * p<0.10  n.s. = no significativo\n")
cat("AVISO: Variables con n < 12 tienen baja potencia estadística.\n")

# Imprimir scatter plots de Moran
for (nm in names(resultados_moran)) {
  r <- resultados_moran[[nm]]
  if (!is.null(r)) print(r$plot)
}

# Gráfico de barras resumen del I de Moran
p_barras_resumen <- ggplot(tabla_moran,
                           aes(x = I_Moran,
                               y = reorder(Variable, I_Moran),
                               fill = Signif)) +
  geom_col(width = 0.75, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(aes(label = Signif,
                x = I_Moran + sign(I_Moran) * 0.005),
            size = 4.5, fontface = "bold", hjust = -0.2) +
  scale_fill_manual(
    values = c("***" = "#c0392b", "**" = "#e74c3c",
               "*"   = "#f39c12", "n.s." = "#bdc3c7"),
    name = "Significación"
  ) +
  labs(
    title    = "Índice de Moran Global por Variable",
    subtitle = "Autocorrelación espacial — Medias anuales 2025, Madrid",
    x        = "I de Moran",
    y        = NULL,
    caption  = "Test Monte Carlo 999 permutaciones · UTM 30N"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 15, margin = margin(b = 8)),
    plot.subtitle    = element_text(color = "gray40", size = 11),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

print(p_barras_resumen)


# ==============================================================================
# 5. ANÁLISIS DE ROBUSTEZ: k-VECINOS vs DISTANCIA EUCLIDIANA (5 km)
# ==============================================================================
# Comprueba si las conclusiones de Moran se mantienen cuando se cambia la
# definición de vecindad de k-vecinos a distancia realista (≤5 km).
# ==============================================================================

cat("\n\n")
cat("==============================================================================\n")
cat("VALIDACIÓN DE ROBUSTEZ: k-VECINOS vs DISTANCIA EUCLIDIANA (5 km)\n")
cat("==============================================================================\n")

moran_por_distancia <- function(dt, dist_km = 5, nsim = 999) {
  n <- nrow(dt)
  if (n < 4) return(NULL)
  sf_obj <- st_as_sf(dt, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
    st_transform(crs = 25830)
  coords        <- st_coordinates(sf_obj)
  nb            <- dnearneigh(coords, d1 = 0, d2 = dist_km * 1000)
  n_sin_vecinos <- sum(card(nb) == 0)
  # Si más del 30% de estaciones queda sin vecinos, la distancia es demasiado restrictiva
  if ((n - n_sin_vecinos) < n * 0.7)
    return(list(I = NA, p = NA, signif = "FALLO", n = n, dist_km = dist_km))
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  mc  <- moran.mc(dt$valor, lw, nsim = nsim, zero.policy = TRUE)
  I   <- round(mc$statistic[[1]], 4)
  p   <- round(mc$p.value, 4)
  sig <- case_when(p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.10 ~ "*", TRUE ~ "n.s.")
  list(I = I, p = p, signif = sig, n = n, dist_km = dist_km,
       n_sin_vecinos = n_sin_vecinos)
}

moran_por_kvecinos <- function(dt, k = NULL, nsim = 999) {
  n <- nrow(dt)
  if (n < 4) return(NULL)
  if (is.null(k)) k <- if (n >= 20) 4L else 3L
  sf_obj <- st_as_sf(dt, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
    st_transform(crs = 25830)
  coords   <- st_coordinates(sf_obj)
  nb       <- knn2nb(knearneigh(coords, k = k))
  lw       <- nb2listw(nb, style = "W")
  mc       <- moran.mc(dt$valor, lw, nsim = nsim)
  I        <- round(mc$statistic[[1]], 4)
  p        <- round(mc$p.value, 4)
  sig      <- case_when(p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.10 ~ "*", TRUE ~ "n.s.")
  dist_mat <- as.matrix(dist(as.data.frame(coords)))
  dist_med <- mean(apply(dist_mat, 1, function(d) sort(d)[k + 1])) / 1000
  list(I = I, p = p, signif = sig, n = n, k = k,
       dist_media_km = round(dist_med, 1))
}

tabla_robustez <- NULL
for (nm in names(todos_datos)) {
  dt   <- todos_datos[[nm]]
  r_k  <- moran_por_kvecinos(dt)
  r_d5 <- moran_por_distancia(dt, dist_km = 5)
  if (!is.null(r_k) && !is.null(r_d5)) {
    tabla_robustez <- rbind(tabla_robustez, data.frame(
      Variable      = nm,
      n             = r_k$n,
      I_k_vecinos   = r_k$I,
      p_k_vecinos   = r_k$p,
      Sig_k         = r_k$signif,
      k             = r_k$k,
      dist_media_km = r_k$dist_media_km,
      I_dist5km     = r_d5$I,
      p_dist5km     = r_d5$p,
      Sig_d5        = r_d5$signif,
      stringsAsFactors = FALSE
    ))
  }
}

cat("\nCOMPARACIÓN LADO A LADO:\n\n")
for (i in seq_len(nrow(tabla_robustez))) {
  row <- tabla_robustez[i, ]
  cat(sprintf("%-30s (n=%2d)\n", row$Variable, row$n))
  cat(sprintf("  k-vecinos (k=%d, dist media=%.1f km): I = %7.4f  p = %.4f  %s\n",
              row$k, row$dist_media_km, row$I_k_vecinos, row$p_k_vecinos, row$Sig_k))
  cat(sprintf("  Distancia 5 km:                       I = %7.4f  p = %.4f  %s",
              row$I_dist5km, row$p_dist5km, row$Sig_d5))
  acuerdo <- (row$Sig_k == "n.s." && row$Sig_d5 == "n.s.") ||
             (row$Sig_k != "n.s." && row$Sig_d5 != "n.s.")
  cat(if (acuerdo) "  \u2713 ACUERDO" else "  \u26a0 DISCORDANCIA")
  cat("\n\n")
}


# ==============================================================================
# 6. GUARDAR TODOS LOS RESULTADOS
# ==============================================================================
cat("--- Guardando resultados en", carpeta_correlacion, "---\n")

# Tablas CSV
write.csv(tabla_moran,
          file.path(carpeta_correlacion, "00_Resumen_Moran_Global.csv"),
          row.names = FALSE)

if (!is.null(tabla_variogramas) && nrow(tabla_variogramas) > 0) {
  write.csv(tabla_variogramas,
            file.path(carpeta_correlacion, "03_Parametros_Variogramas.csv"),
            row.names = FALSE)
}

write.csv(tabla_robustez,
          file.path(carpeta_correlacion, "04_Robustez_kVecinos_vs_Distancia5km.csv"),
          row.names = FALSE)

# Gráfico de barras Moran
ggsave(file.path(carpeta_correlacion, "01_Resumen_Barras_Moran_Global.png"),
       plot = p_barras_resumen, width = 11, height = 8, dpi = 300)

# Líneas temporales por estación (datos mensuales)
for (var_label in names(plots_lineas)) {
  nombre_safe <- gsub("[^a-zA-Z0-9]", "", var_label)
  ggsave(
    file.path(carpeta_correlacion, paste0("00b_Lineas_Mensuales_", nombre_safe, ".png")),
    plot  = plots_lineas[[var_label]],
    width = 12, height = 7, dpi = 300
  )
}

# CV entre estaciones por mes
ggsave(file.path(carpeta_correlacion, "00c_CV_Estaciones_Mes.png"),
       plot = plot_cv_facet, width = 12, height = 9, dpi = 300)

# Perfiles temporales de NO2 (mensual, semanal, horario)
ggsave(file.path(carpeta_correlacion, "00d_Perfil_Mensual_NO2.png"),
       plot = plot_perfil_mensual, width = 11, height = 7, dpi = 300)

ggsave(file.path(carpeta_correlacion, "00e_Perfil_Semanal_NO2.png"),
       plot = plot_perfil_semanal, width = 11, height = 7, dpi = 300)

ggsave(file.path(carpeta_correlacion, "00f_Perfil_Horario_NO2.png"),
       plot = plot_perfil_horario, width = 12, height = 7, dpi = 300)

# Scatter plots de correlación inicial
ggsave(file.path(carpeta_correlacion, "00g_Scatter_NO2_vs_Viento.png"),
       plot = plot_scatter_viento, width = 11, height = 7, dpi = 300)

ggsave(file.path(carpeta_correlacion, "00h_Scatter_NO2_vs_Trafico.png"),
       plot = plot_scatter_traf, width = 11, height = 7, dpi = 300)

ggsave(file.path(carpeta_correlacion, "00i_Panel_Correlaciones_NO2.png"),
       plot = plot_panel_correlaciones, width = 14, height = 7, dpi = 300)

# Variograma comparativo
ggsave(file.path(carpeta_correlacion, "02_Variogramas_Comparativo.png"),
       plot = plot_variogramas_comparativo, width = 14, height = 5, dpi = 300)

# Variogramas individuales
for (nm in names(resultados_variogramas)) {
  nombre_safe <- gsub("[^a-zA-Z0-9]", "", nm)
  ggsave(
    file.path(carpeta_correlacion, paste0("02_Variograma_", nombre_safe, ".png")),
    plot  = crear_plot_variograma(resultados_variogramas[[nm]], nm),
    width = 10, height = 7, dpi = 300
  )
}

# Scatter plots de Moran (uno por variable)
for (i in seq_along(names(resultados_moran))) {
  nm <- names(resultados_moran)[i]
  r  <- resultados_moran[[i]]
  if (!is.null(r)) {
    nombre_archivo <- paste0(
      sprintf("%02d", i + 1), "_Moran_Scatter_",
      gsub("[^a-zA-Z0-9]", "", nm), ".png"
    )
    ggsave(file.path(carpeta_correlacion, nombre_archivo),
           plot = r$plot, width = 12, height = 9, dpi = 300, scale = 1.1)
  }
}

# Homogeneidad espacial — gráfico resumen
ggsave(file.path(carpeta_correlacion, "05_Homogeneidad_Espacial_Resumen.png"),
       plot = plot_homogeneidad_resumen, width = 11, height = 7, dpi = 300)

# Homogeneidad espacial — heatmaps por variable
for (v in names(plots_heatmaps)) {
  nombre_safe <- gsub("[^a-zA-Z0-9]", "", v)
  ggsave(
    file.path(carpeta_correlacion, paste0("05_Heatmap_Cor_", nombre_safe, ".png")),
    plot  = plots_heatmaps[[v]],
    width = 10, height = 8, dpi = 300
  )
}

# Homogeneidad espacial — tabla CSV
write.csv(tabla_homogeneidad,
          file.path(carpeta_correlacion, "05_Homogeneidad_Espacial_Comparativa.csv"),
          row.names = FALSE)

# README
vars_con_variograma <- paste(names(resultados_variogramas), collapse = ", ")
readme_lineas <- c(
  "# ANÁLISIS DE CORRELACIÓN ESPACIAL — MADRID 2025",
  "",
  "## Metodología",
  "",
  "### Variograma experimental",
  paste0("- Variables con n >= ", N_MIN_VAR, " estaciones: ", vars_con_variograma),
  paste0("- Cutoff: ", CUTOFF_KM, " km  |  Anchura de bin: ", WIDTH_KM, " km"),
  "- Modelo teórico ajustado: esférico (fit.variogram)",
  "",
  "### Índice de Moran Global",
  paste0("- Test Monte Carlo (", NSIM, " permutaciones) + test analítico (comparación)"),
  "- Vecindario: k-vecinos más cercanos (k=4 si n>=20; k=3 si n<20)",
  "- Proyección: UTM 30N (EPSG:25830)",
  "- Hipótesis nula: valores distribuidos aleatoriamente en el espacio",
  "",
  "### Análisis de robustez",
  "- Comparación k-vecinos vs vecindad por distancia euclídea (5 km)",
  "",
  "### Homogeneidad espacial (correlaciones temporales entre estaciones)",
  "- Correlación de Pearson entre pares de estaciones usando serie horaria completa",
  "- Calculado para 2022 y 2025 para verificar estabilidad interanual",
  "- Decisión: r >= 0.9 → media ciudad | r 0.7-0.9 → zona gris | r < 0.7 → interpolar",
  "",
  "## Archivos generados",
  "- 00_Resumen_Moran_Global.csv  : I, p(MC), p(analítico), Dif_p por variable",
  "- 01_Resumen_Barras_Moran_Global.png",
  "- 02_Variogramas_Comparativo.png",
  "- 02_Variograma_*.png          : variograma individual por variable",
  "- 03_Parametros_Variogramas.csv: nugget, sill, range por variable",
  "- 0X_Moran_Scatter_*.png       : scatter plot de Moran por variable",
  "- 04_Robustez_kVecinos_vs_Distancia5km.csv",
  "- 05_Homogeneidad_Espacial_Resumen.png : gráfico barras resumen",
  "- 05_Heatmap_Cor_*.png         : heatmap correlación por variable",
  "- 05_Homogeneidad_Espacial_Comparativa.csv : tabla 2022 vs 2025",
  "",
  "## Resultados principales",
  "- Significativo   : Humedad Relativa (p~0.04, **)",
  "- Marginal        : Presión Barométrica (p~0.09, *)",
  "- No significativo: NO2, Temperatura y resto",
  "",
  "### Homogeneidad espacial",
  "- Media ciudad: Temperatura (r~0.99), Presión Barométrica (r~0.99), Humedad Relativa (r~0.97)",
  "- Zona gris: Velocidad Viento (r~0.72-0.79)",
  "- Interpolar: Dir.Viento (r~0.51), Precipitaciones (r~0.59)",
  "- Patrones estables entre 2022 y 2025",
  "",
  paste0("Generado: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(readme_lineas, file.path(carpeta_correlacion, "README.txt"))

cat("\n\u2713 RESULTADOS GUARDADOS EN:\n")
cat(carpeta_correlacion, "\n\n")
cat("Archivos generados:\n  \u2022 ")
list.files(carpeta_correlacion) |> sort() |> cat(sep = "\n  \u2022 ")
cat("\n")

# ==============================================================================
# 7. GRÁFICA ESTRELLA — PERFIL HORARIO NO₂ + TRÁFICO POR TEMPORADA Y TIPO DE DÍA
# ==============================================================================
# Inspirada en Fig. 11 de Laña (2024). Cuatro paneles (temporadas) × dos versiones
# (martes laborable / domingo festivo).
# Eje izquierdo (línea roja): concentración horaria de NO₂ (µg/m³)
# Eje derecho  (barras azules): intensidad de tráfico (veh/hora)
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

# ── 7a. Preparar NO₂ horario ──────────────────────────────────────────────────
dt_no2_h <- copy(dt_horario_no2)

if (!"Hora_num" %in% names(dt_no2_h))
  dt_no2_h[, Hora_num := as.integer(gsub("H0?", "", as.character(HORA)))]

# Temporada según mes (hemisferio norte)
dt_no2_h[, mes_n := as.integer(format(FECHA, "%m"))]
dt_no2_h[, Temporada := fcase(
  mes_n %in% c(12L, 1L, 2L),  "Invierno",
  mes_n %in% c(3L,  4L, 5L),  "Primavera",
  mes_n %in% c(6L,  7L, 8L),  "Verano",
  mes_n %in% c(9L, 10L, 11L), "Oto\u00f1o"
)]

# Tipo de día: solo martes (2) y domingo (7)
dt_no2_h[, dia_sem := as.integer(format(FECHA, "%u"))]
dt_no2_h[, Tipo_dia2 := fcase(
  dia_sem == 2L, "Martes (laborable)",
  dia_sem == 7L, "Domingo (festivo)",
  default = NA_character_
)]

# Media horaria de NO₂ (todas las estaciones)
no2_perfil <- dt_no2_h[!is.na(DATO) & !is.na(Tipo_dia2), .(
  NO2_media = mean(DATO, na.rm = TRUE)
), by = .(Hora_num, Temporada, Tipo_dia2)]


dt_th <- readRDS(here("data", "processed", "trafico_madrid_2025_horario.rds"))

# Hora_num: convierte timestamp 00:00→1, 01:00→2 … 23:00→24 (== convenci\u00f3n H01-H24)
dt_th[, Hora_num := as.integer(format(fecha_hora, "%H")) + 1L]
dt_th[, mes_n   := as.integer(format(FECHA, "%m"))]
dt_th[, Temporada := fcase(
  mes_n %in% c(12L, 1L, 2L),  "Invierno",
  mes_n %in% c(3L,  4L, 5L),  "Primavera",
  mes_n %in% c(6L,  7L, 8L),  "Verano",
  mes_n %in% c(9L, 10L, 11L), "Oto\u00f1o"
)]
dt_th[, dia_sem := as.integer(format(FECHA, "%u"))]
dt_th[, Tipo_dia2 := fcase(
  dia_sem == 2L, "Martes (laborable)",
  dia_sem == 7L, "Domingo (festivo)",
  default = NA_character_
)]

traf_perfil <- dt_th[!is.na(intensidad) & !is.na(Tipo_dia2), .(
  Traf_media = mean(intensidad, na.rm = TRUE)
), by = .(Hora_num, Temporada, Tipo_dia2)]

rm(dt_th); gc()
cat("   \u2705 Tr\u00e1fico agregado y memoria liberada.\n")

# ── 7c. Unir y calcular factor de escala para el eje doble ───────────────────
dt_estrella <- merge(no2_perfil, traf_perfil,
                     by = c("Hora_num", "Temporada", "Tipo_dia2"))

dt_estrella[, Temporada := factor(Temporada,
  levels = c("Invierno", "Primavera", "Verano", "Oto\u00f1o"))]
dt_estrella[, Tipo_dia2 := factor(Tipo_dia2,
  levels = c("Martes (laborable)", "Domingo (festivo)"))]

# Factor de escala global: mapea el rango de tráfico al rango de NO₂
escala <- max(dt_estrella$NO2_media,  na.rm = TRUE) /
          max(dt_estrella$Traf_media, na.rm = TRUE)

horas_labels <- setNames(
  sprintf("%02d:00", 0:23),   # 00:00 … 23:00
  1:24                         # Hora_num 1 … 24
)

# ── 7d. Función de construcción ──────────────────────────────────────────────
hacer_estrella <- function(datos, subtitulo) {

  ggplot(datos, aes(x = Hora_num)) +

    # Barras de tráfico (escaladas al eje izquierdo)
    geom_col(aes(y = Traf_media * escala),
             fill = "#aed6f1", color = NA, alpha = 0.8, width = 0.85) +

    # Línea y puntos de NO₂
    geom_line(aes(y = NO2_media),
              color = "#c0392b", linewidth = 1.1) +
    geom_point(aes(y = NO2_media),
               color = "#c0392b", size = 1.8, shape = 21,
               fill = "white", stroke = 1) +

    # Ejes dobles
    scale_y_continuous(
      name   = expression(NO[2] ~ "(" * mu * "g/m"^3 * ")"),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = sec_axis(
        transform = ~ . / escala,
        name      = "Intensidad Tr\u00e1fico (veh\u00edculos/hora)"
      )
    ) +
    scale_x_continuous(
      breaks = c(1, 5, 9, 13, 17, 21),
      labels = c("01:00", "05:00", "09:00", "13:00", "17:00", "21:00"),
      minor_breaks = seq(1, 24, 1)
    ) +
    facet_wrap(~ Temporada, ncol = 2) +
    labs(
      title    = subtitulo,
      subtitle = paste0("Media horaria sobre todas las estaciones \u00b7 Madrid 2025\n",
                        "L\u00ednea roja: NO\u2082 (\u00b5g/m\u00b3)  |  Barras azules: Intensidad tr\u00e1fico (veh/hora)"),
      x        = "Hora del d\u00eda",
      caption  = paste0("NO\u2082: media de todas las estaciones de medida  \u00b7  ",
                        "Tr\u00e1fico: media de todos los distritos de Madrid")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "gray35", size = 9, lineheight = 1.4),
      plot.caption       = element_text(color = "gray55", size = 8),
      strip.text         = element_text(face = "bold", size = 11),
      strip.background   = element_rect(fill = "gray95", color = "gray80"),
      axis.title.y.left  = element_text(color = "#c0392b", face = "bold", size = 10),
      axis.text.y.left   = element_text(color = "#c0392b"),
      axis.title.y.right = element_text(color = "#2980b9", face = "bold", size = 10),
      axis.text.y.right  = element_text(color = "#2980b9"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor.y = element_blank(),
      panel.spacing      = unit(1.4, "lines")
    )
}

# ── 7e. Generar las dos versiones ─────────────────────────────────────────────
plot_estrella_martes <- hacer_estrella(
  dt_estrella[Tipo_dia2 == "Martes (laborable)"],
  "Perfil Horario de NO\u2082 y Tr\u00e1fico por Temporada \u2014 Martes (laborable)"
)

plot_estrella_domingo <- hacer_estrella(
  dt_estrella[Tipo_dia2 == "Domingo (festivo)"],
  "Perfil Horario de NO\u2082 y Tr\u00e1fico por Temporada \u2014 Domingo (festivo)"
)

print(plot_estrella_martes)
print(plot_estrella_domingo)

ggsave(file.path(carpeta_correlacion, "06_Estrella_Perfil_NO2_Trafico_Martes.png"),
       plot = plot_estrella_martes, width = 14, height = 8, dpi = 300)
ggsave(file.path(carpeta_correlacion, "06_Estrella_Perfil_NO2_Trafico_Domingo.png"),
       plot = plot_estrella_domingo, width = 14, height = 8, dpi = 300)

cat("\u2705 Gr\u00e1ficas estrella guardadas en:\n  ", carpeta_correlacion, "\n")

