# ==============================================================================
# ÍNDICE DE MORAN GLOBAL Y LOCAL (LISA) — ¿EXISTE CORRELACIÓN ESPACIAL?
# Madrid 2025 · Contaminación (NO₂) · Clima · Intensidad de tráfico
# Escala: media anual por estación (24 estaciones, dataset maestro mensual)
# ==============================================================================
# Estructura:
#   1. Carga y agregación anual por estación (NO₂, clima, tráfico)
#   2. Índice de Moran GLOBAL: ¿hay autocorrelación espacial en conjunto?
#        -> tabla resumen + gráfico de barras + scatter plots de Moran
#   3. Índice de Moran LOCAL (LISA): ¿DÓNDE están los clusters/outliers?
#        -> mapas de cluster (Alto-Alto, Bajo-Bajo, Alto-Bajo, Bajo-Alto)
#   4. Guardar todo en outputs/analysis/moran_espacial/
# ==============================================================================
# ==============================================================================
# 16_indice_moran_global_local.R
# ==============================================================================
#
# QUÉ HACE:
#
# - Carga el dataset maestro mensual correspondiente al periodo 2019–2025.
# - Selecciona exclusivamente las observaciones de 2025.
# - Calcula una media anual para cada estación.
#
# - Analiza nueve variables:
#     · NO2.
#     · Intensidad del tráfico.
#     · Carga de tráfico.
#     · Temperatura.
#     · Humedad relativa.
#     · Precipitaciones.
#     · Presión barométrica.
#     · Radiación solar.
#     · Velocidad del viento.
#
# - Construye una estructura de vecindad espacial mediante los vecinos más
#   próximos:
#     · k = 4 cuando hay al menos 20 estaciones.
#     · k = 3 cuando hay menos de 20 estaciones.
#
# - Calcula el índice de Moran global.
# - Utiliza 999 permutaciones Monte Carlo para evaluar su significación.
# - Calcula también el contraste analítico del índice de Moran.
#
# MORAN GLOBAL:
#
# - Evalúa si los valores similares se encuentran espacialmente próximos.
# - Genera una tabla resumen con:
#     · Índice de Moran.
#     · p-valor de Monte Carlo.
#     · p-valor analítico.
#     · Número de estaciones.
#     · Número de vecinos.
#     · Nivel de significación.
#
# - Genera gráficos de Moran donde se compara el valor estandarizado de cada
#   estación con el rezago espacial de sus vecinas.
#
# MORAN LOCAL O LISA:
#
# - Identifica dónde aparecen agrupaciones espaciales y valores atípicos.
# - Clasifica las estaciones como:
#     · Alto-Alto: valor alto rodeado de valores altos.
#     · Bajo-Bajo: valor bajo rodeado de valores bajos.
#     · Alto-Bajo: valor alto rodeado de valores bajos.
#     · Bajo-Alto: valor bajo rodeado de valores altos.
#     · No significativo.
#
# - Representa los resultados LISA sobre el mapa de Madrid.
#
# FINALIDAD PARA EL TFM:
#
# Este script comprueba formalmente si existe dependencia espacial en el NO2,
# en el tráfico y en las variables meteorológicas.
#
# El índice de Moran global responde a la pregunta:
#
#     ¿Existe autocorrelación espacial en el conjunto de Madrid?
#
# El índice de Moran local responde a la pregunta:
#
#     ¿En qué estaciones se localizan los clusters y los outliers espaciales?
#
# Una autocorrelación espacial significativa indica que las observaciones no son
# espacialmente independientes. Esto ayuda a justificar la utilización de un
# modelo espacial como INLA-SPDE.
#
# Los mapas LISA permiten detectar:
#
# - Zonas con concentraciones elevadas rodeadas de valores elevados.
# - Zonas con concentraciones bajas espacialmente agrupadas.
# - Estaciones atípicas respecto a su entorno.
# - Patrones espaciales comunes entre NO2, tráfico y meteorología.
#
# OBSERVACIÓN METODOLÓGICA:
#
# Los resultados dependen de la definición de vecindad elegida. En este caso se
# emplean vecinos más próximos y no contigüidad entre polígonos, ya que las
# unidades analizadas son estaciones puntuales.
#
# SALIDAS:
#
# outputs/analysis/moran_espacial/
#
# Entre las principales salidas se encuentran:
#
# - Resumen del índice de Moran global en CSV.
# - Tabla formateada del índice de Moran global.
# - Gráfico de barras con los índices de Moran.
# - Panel de scatter plots de Moran.
# - Un mapa LISA para cada variable.
# - Tabla resumen con los conteos de clusters LISA.


library(data.table)
library(sf)
library(spdep)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(gt)
library(here)

set.seed(4712)
NSIM <- 999 # permutaciones Monte Carlo (Moran global)
ALPHA <- 0.05 # significación para clasificar clusters LISA

carpeta_out <- here("outputs", "analysis", "moran_espacial")
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. CARGA Y AGREGACIÓN ANUAL POR ESTACIÓN
# ==============================================================================
# El dataset maestro ya asocia, por estación y mes, el NO₂, las variables
# climáticas y el tráfico (intensidad/carga) de su área de influencia.
# Se agrega a media anual 2025 para tener un único valor por estación.
# ==============================================================================

dt_maestro <- readRDS(here(
  "data", "processed", "Maestro", "mensual",
  "dataset_maestro_inla_2019_2025_MENSUAL.rds"
))
setDT(dt_maestro)

dt_2025 <- dt_maestro[ANIO == 2025]

vars_moran <- c(
  "DATO_NO2", "intensidad_raw", "carga_raw",
  "Temperatura_raw", "Humedad_Relativa_raw", "Precipitaciones_raw",
  "Presion Barométrica_raw", "Radiación Solar_raw", "Velocidad Viento_raw"
)
labels_moran <- c(
  "NO₂ (µg/m³)", "Intensidad de tráfico", "Carga de tráfico",
  "Temperatura (°C)", "Humedad Relativa (%)", "Precipitaciones (mm)",
  "Presión Barométrica (hPa)", "Radiación Solar (W/h)", "Velocidad Viento (m/s)"
)
grupo_moran <- c(
  "Contaminación", "Tráfico", "Tráfico",
  "Clima", "Clima", "Clima", "Clima", "Clima", "Clima"
)

lista_moran <- setNames(
  lapply(vars_moran, function(v) {
    dt_2025[!is.na(get(v)), .(
      valor    = mean(get(v), na.rm = TRUE),
      LONGITUD = unique(LONGITUD),
      LATITUD  = unique(LATITUD)
    ), by = ESTACION]
  }),
  labels_moran
)

cat("Estaciones por variable (media anual 2025):\n")
for (nm in names(lista_moran)) {
  cat(sprintf("  %-30s n = %d\n", nm, nrow(lista_moran[[nm]])))
}

# Distritos (fondo cartográfico para los mapas LISA)
distritos <- st_read(here("data", "raw", "geometrias", "DISTRITOS.shp"), quiet = TRUE)


# ==============================================================================
# 2. FUNCIONES: PREPARAR VECINDAD ESPACIAL + MORAN GLOBAL + MORAN LOCAL
# ==============================================================================

# Matriz de vecindad k-NN (misma regla que el resto del proyecto: n≥20 -> k=4)
preparar_vecindad <- function(dt, k = NULL) {
  n <- nrow(dt)
  if (is.null(k)) k <- if (n >= 20) 4L else 3L

  sf_obj <- st_as_sf(dt, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
    st_transform(crs = 25830)
  coords <- st_coordinates(sf_obj)
  nb <- knn2nb(knearneigh(coords, k = k))
  lw <- nb2listw(nb, style = "W")

  list(sf_obj = sf_obj, coords = coords, lw = lw, n = n, k = k)
}

# --- Moran GLOBAL: ¿hay autocorrelación espacial en conjunto? --------------
moran_global <- function(dt, prep, variable_label, nsim = NSIM) {
  mc <- moran.mc(dt$valor, prep$lw, nsim = nsim)
  I <- round(mc$statistic[[1]], 4)
  pv_mc <- round(mc$p.value, 4)

  ta <- moran.test(dt$valor, prep$lw, alternative = "greater")
  pv_an <- round(ta$p.value, 4)

  sig <- if (pv_mc < 0.01) "***" else if (pv_mc < 0.05) "**" else if (pv_mc < 0.10) "*" else "n.s."

  # Scatter plot de Moran (z vs. rezago espacial Wz)
  z <- as.numeric(scale(dt$valor))
  Wz <- lag.listw(prep$lw, z)
  df_plot <- data.frame(z = z, Wz = Wz, estacion = dt$ESTACION)

  label_ann <- paste0("I = ", I, "  ", sig, "\np(MC) = ", pv_mc)

  p <- ggplot(df_plot, aes(x = z, y = Wz)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#c0392b", alpha = 0.07) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#1a5276", alpha = 0.07) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.3) +
    geom_abline(slope = I, intercept = 0, color = "#e74c3c", linewidth = 0.8) +
    geom_point(size = 2.3, color = "#2c3e50", alpha = 0.9) +
    annotate("label",
      x = Inf, y = Inf, hjust = 1.08, vjust = 1.2,
      label = label_ann, size = 3.2, fontface = "bold",
      fill = "white", color = "#2c3e50", label.size = 0.3
    ) +
    labs(title = variable_label, x = "z-score", y = "Wz (rezago espacial)") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )

  list(
    I = I, p_mc = pv_mc, p_analitico = pv_an, signif = sig,
    n = prep$n, k = prep$k, plot = p
  )
}

# --- Moran LOCAL (LISA): ¿DÓNDE están los clusters? -------------------------
moran_local <- function(dt, prep, variable_label, alpha = ALPHA) {
  loc <- localmoran(dt$valor, prep$lw, alternative = "two.sided")
  col_p <- grep("^Pr", colnames(loc), value = TRUE)[1]
  p_val <- loc[, col_p]

  z <- as.numeric(scale(dt$valor))
  Wz <- lag.listw(prep$lw, z)

  cuadrante <- ifelse(z >= 0 & Wz >= 0, "Alto-Alto (HH)",
    ifelse(z < 0 & Wz < 0, "Bajo-Bajo (LL)",
      ifelse(z >= 0 & Wz < 0, "Alto-Bajo (HL)",
        "Bajo-Alto (LH)"
      )
    )
  )
  cluster <- ifelse(p_val < alpha, cuadrante, "No significativo")

  sf_res <- prep$sf_obj
  sf_res$ESTACION <- dt$ESTACION
  sf_res$valor <- dt$valor
  sf_res$Ii <- loc[, "Ii"]
  sf_res$p_valor <- p_val
  sf_res$Cluster <- factor(cluster,
    levels = c(
      "Alto-Alto (HH)", "Bajo-Bajo (LL)",
      "Alto-Bajo (HL)", "Bajo-Alto (LH)",
      "No significativo"
    )
  )
  sf_res$X <- st_coordinates(sf_res)[, 1]
  sf_res$Y <- st_coordinates(sf_res)[, 2]

  n_clusters <- sum(cluster != "No significativo")
  cat(sprintf(
    "  %-30s %d/%d estaciones en cluster/outlier significativo (α=%.2f)\n",
    variable_label, n_clusters, prep$n, alpha
  ))

  sf_res
}

colores_lisa <- c(
  "Alto-Alto (HH)" = "#c0392b",
  "Bajo-Bajo (LL)" = "#1a5276",
  "Alto-Bajo (HL)" = "#f39c12",
  "Bajo-Alto (LH)" = "#8e44ad",
  "No significativo" = "grey75"
)

mapa_lisa <- function(sf_res, variable_label) {
  ggplot() +
    geom_sf(data = distritos, fill = "grey97", colour = "grey85", linewidth = 0.35) +
    geom_point(
      data = sf_res,
      aes(x = X, y = Y, fill = Cluster),
      shape = 21, size = 5, colour = "grey20", stroke = 0.5
    ) +
    geom_text_repel(
      data = sf_res,
      aes(x = X, y = Y, label = ESTACION),
      size = 2.4, colour = "grey20", bg.color = "white", bg.r = 0.12,
      box.padding = 0.4, point.padding = 0.3, segment.size = 0.3,
      max.overlaps = 30, seed = 7231
    ) +
    scale_fill_manual(values = colores_lisa, name = "Cluster LISA", drop = FALSE) +
    labs(
      title    = paste("Moran Local (LISA) —", variable_label),
      subtitle = sprintf("Madrid 2025 · media anual por estación · significación α = %.2f (999 perm.)", ALPHA),
      caption  = "HH: foco de valores altos · LL: foco de valores bajos · HL/LH: outliers espaciales"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(colour = "grey40", size = 9),
      plot.caption     = element_text(colour = "grey55", size = 8),
      axis.title       = element_blank(),
      axis.text        = element_text(size = 7, colour = "grey60"),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    coord_sf(crs = st_crs(distritos))
}


# ==============================================================================
# 3. CALCULAR MORAN GLOBAL Y LOCAL PARA TODAS LAS VARIABLES
# ==============================================================================

cat("\n--- Calculando Moran's I Global ---\n")
resultados_global <- list()
resultados_local <- list()

for (nm in names(lista_moran)) {
  dt <- lista_moran[[nm]]
  prep <- preparar_vecindad(dt)

  rg <- moran_global(dt, prep, nm)
  resultados_global[[nm]] <- rg
  cat(sprintf("  %-30s I = %7.4f  p(MC) = %.4f  %s\n", nm, rg$I, rg$p_mc, rg$signif))

  resultados_local[[nm]] <- moran_local(dt, prep, nm)
}


# ==============================================================================
# 4. MORAN GLOBAL — TABLA RESUMEN + BARRAS + SCATTER (grid)
# ==============================================================================

tabla_global <- rbindlist(lapply(names(resultados_global), function(nm) {
  r <- resultados_global[[nm]]
  data.table(
    Grupo = grupo_moran[match(nm, labels_moran)],
    Variable = nm,
    N = r$n,
    k = r$k,
    I_Moran = r$I,
    p_MC = r$p_mc,
    p_Analitico = r$p_analitico,
    Significacion = r$signif
  )
}))
setorder(tabla_global, Grupo, -I_Moran)

cat("\n================================================================\n")
cat("   ÍNDICE DE MORAN GLOBAL — NO₂, Clima y Tráfico (Madrid 2025)\n")
cat("================================================================\n")
print(tabla_global)
cat("\nSignif: *** p<0.01  ** p<0.05  * p<0.10  n.s. = no significativo\n")

fwrite(tabla_global, file.path(carpeta_out, "00_Resumen_Moran_Global.csv"))

tbl_gt <- tabla_global |>
  gt(groupname_col = "Grupo") |>
  tab_header(
    title    = md("**Índice de Moran Global**"),
    subtitle = md("NO₂, clima y tráfico · media anual por estación · Madrid 2025")
  ) |>
  cols_label(
    Variable = md("**Variable**"), N = md("**N**"), k = md("**k**"),
    I_Moran = md("**I de Moran**"), p_MC = md("**p (MC)**"),
    p_Analitico = md("**p (anal.)**"), Significacion = md("**Sig.**")
  ) |>
  fmt_number(columns = c(I_Moran, p_MC, p_Analitico), decimals = 4) |>
  data_color(columns = I_Moran, method = "numeric", palette = c("#d4e9f7", "#08519c")) |>
  tab_source_note(source_note = md("Test Monte Carlo (999 perm.) · k-NN, k=4 · UTM 30N (EPSG:25830)")) |>
  tab_style(
    style = list(cell_fill(color = "#1a3a5c"), cell_text(color = "white", weight = "bold", size = px(15))),
    locations = cells_title(groups = "title")
  ) |>
  tab_style(
    style = list(cell_fill(color = "#1a3a5c"), cell_text(color = "#d0e4f5", size = px(11))),
    locations = cells_title(groups = "subtitle")
  ) |>
  tab_style(
    style = list(cell_fill(color = "#2c5f8a"), cell_text(color = "white", weight = "bold")),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = list(cell_fill(color = "#e8eef4"), cell_text(weight = "bold")),
    locations = cells_row_groups()
  )

gtsave(tbl_gt, filename = file.path(carpeta_out, "00_Tabla_Moran_Global.png"), zoom = 2, expand = 20)

# Barras resumen
p_barras <- ggplot(
  tabla_global,
  aes(x = I_Moran, y = reorder(Variable, I_Moran), fill = Significacion)
) +
  geom_col(width = 0.72, color = "white", linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_text(aes(label = Significacion, x = I_Moran + sign(I_Moran) * 0.01),
    size = 4, fontface = "bold", hjust = -0.15
  ) +
  scale_fill_manual(
    values = c(
      "***" = "#c0392b", "**" = "#e74c3c",
      "*" = "#f39c12", "n.s." = "#bdc3c7"
    ),
    name = "Significación"
  ) +
  labs(
    title = "Índice de Moran Global por Variable",
    subtitle = "NO₂, clima y tráfico · medias anuales 2025, Madrid",
    x = "I de Moran", y = NULL,
    caption = "Test Monte Carlo 999 permutaciones · k-NN k=4 · UTM 30N"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    panel.grid.minor = element_blank(), legend.position = "bottom"
  )

ggsave(file.path(carpeta_out, "01_Barras_Moran_Global.png"),
  plot = p_barras, width = 9, height = 6, dpi = 300, bg = "white"
)

# Grid de scatter plots de Moran (una figura con las 9 variables)
grid_scatter <- arrangeGrob(grobs = lapply(resultados_global, `[[`, "plot"), ncol = 3)
ggsave(file.path(carpeta_out, "02_Scatter_Moran_Global_Todas_Variables.png"),
  plot = grid_scatter, width = 13, height = 11, dpi = 250, bg = "white"
)

cat("\n✅ Moran global guardado (tabla, barras, scatter grid) en:\n   ", carpeta_out, "\n")


# ==============================================================================
# 5. MORAN LOCAL (LISA) — MAPAS DE CLUSTER POR VARIABLE
# ==============================================================================

nombres_archivo <- gsub("[^a-zA-Z0-9]+", "_", names(resultados_local))

for (i in seq_along(resultados_local)) {
  nm <- names(resultados_local)[i]
  p <- mapa_lisa(resultados_local[[nm]], nm)
  arc <- sprintf("03_LISA_%02d_%s.png", i, nombres_archivo[i])
  ggsave(file.path(carpeta_out, arc), plot = p, width = 11, height = 9, dpi = 250, bg = "white")
}

cat(sprintf("✅ %d mapas LISA guardados en:\n   %s\n", length(resultados_local), carpeta_out))

# Resumen de conteos de cluster por variable
tabla_lisa <- rbindlist(lapply(names(resultados_local), function(nm) {
  tb <- table(resultados_local[[nm]]$Cluster)
  data.table(Variable = nm, Cluster = names(tb), N_estaciones = as.integer(tb))
}))
tabla_lisa_ancha <- dcast(tabla_lisa, Variable ~ Cluster, value.var = "N_estaciones", fill = 0)

fwrite(tabla_lisa_ancha, file.path(carpeta_out, "04_Resumen_LISA_Conteos.csv"))

tbl_lisa_gt <- tabla_lisa_ancha |>
  gt() |>
  tab_header(
    title    = md("**Resumen de Clusters LISA por Variable**"),
    subtitle = md("Número de estaciones en cada tipo de cluster/outlier · Madrid 2025")
  ) |>
  tab_source_note(source_note = md("HH: foco alto · LL: foco bajo · HL/LH: outliers espaciales · α = 0.05")) |>
  tab_style(
    style = list(cell_fill(color = "#1a3a5c"), cell_text(color = "white", weight = "bold", size = px(15))),
    locations = cells_title(groups = "title")
  ) |>
  tab_style(
    style = list(cell_fill(color = "#1a3a5c"), cell_text(color = "#d0e4f5", size = px(11))),
    locations = cells_title(groups = "subtitle")
  ) |>
  tab_style(
    style = list(cell_fill(color = "#2c5f8a"), cell_text(color = "white", weight = "bold")),
    locations = cells_column_labels()
  )

gtsave(tbl_lisa_gt, filename = file.path(carpeta_out, "04_Tabla_LISA_Conteos.png"), zoom = 2, expand = 20)

cat("\n✅ Resumen de clusters LISA guardado (csv + png) en:\n   ", carpeta_out, "\n")

cat("\n================================================================\n")
cat(" PROCESO COMPLETADO — Moran Global y Local (LISA)\n")
cat(" Salidas en:", carpeta_out, "\n")
cat("================================================================\n")
