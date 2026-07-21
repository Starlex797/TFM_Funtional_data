#===============================================================================
# Interpolation climate variables 
#===============================================================================

# ==============================================================================
# EJECUCIÓN
# ==============================================================================
library(here)
source(here("R", "FUNCIONES_INTERPOLACION.R"))
library(ggplot2)
library(viridis)


mis_variables_clima <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  "Presion_Barometrica", "Radiacion_Solar", "Velocidad Viento"
)
etiquetas_variables <- c(
  Temperatura = "Temperatura",
  Humedad_Relativa = "Humedad relativa",
  Precipitaciones = "Precipitaciones",
  Presion_Barometrica = "Presion barometrica",
  Radiacion_Solar = "Radiacion solar",
  `Velocidad Viento` = "Velocidad del viento"
)
dt_meteo <- readRDS(here("data", "processed", "Clima", "diario", "meteo_madrid_2025_diario.rds"))
setDT(dt_meteo)

# Normalizar columnas con acentos para no depender del locale de Windows/R.
col_presion <- grep("^Presion Barom", names(dt_meteo), value = TRUE)
col_radiacion <- grep("Solar$", names(dt_meteo), value = TRUE)
if (length(col_presion) == 1L) setnames(dt_meteo, col_presion, "Presion_Barometrica")
if (length(col_radiacion) == 1L) setnames(dt_meteo, col_radiacion, "Radiacion_Solar")

# Día con lluvia y viento seleccionado para el análisis
# (máxima precipitación media con 9 estaciones completas: 27.9 mm, viento 1.16 m/s)
FECHA_ANALISIS <- as.Date("2025-04-03")

dt_meteo_dia <- dt_meteo[FECHA == FECHA_ANALISIS]
cat(sprintf("Día de análisis: %s | Estaciones: %d\n",
            FECHA_ANALISIS, nrow(dt_meteo_dia)))
cat(sprintf("  Precipitación media: %.2f mm | Viento medio: %.2f m/s\n",
            mean(dt_meteo_dia$Precipitaciones, na.rm = TRUE),
            mean(dt_meteo_dia$`Velocidad Viento`, na.rm = TRUE)))

tabla_rmse_loocv <- comparar_interpolaciones_loocv(
  dt_meteo   = dt_meteo_dia,
  variables  = mis_variables_clima,
  k_vecinos  = 5
)

# ==============================================================================
# TABLA Y GRÁFICO COMPARATIVO DE MÉTODOS (PNG)
# ==============================================================================
library(gt)

dir.create(here("output", "figures", "interpolacion_clima"),
           recursive = TRUE, showWarnings = FALSE)

# --- Tabla gt como PNG ---
tabla_display <- copy(tabla_rmse_loocv)
tabla_display[, Mejor_Metodo := NULL]

tbl <- gt(tabla_display) |>
  tab_header(
    title    = "Comparación de Métodos de Interpolación — LOOCV",
    subtitle = sprintf("Fecha: %s | k vecinos = 5", format(FECHA_ANALISIS, "%d %b %Y"))
  ) |>
  tab_spanner(label = "RMSE", columns = starts_with("RMSE_")) |>
  tab_spanner(label = "MAE",  columns = starts_with("MAE_")) |>
  cols_label(
    RMSE_Media = "Media", RMSE_NN = "1-NN", RMSE_kNN = "kNN", RMSE_IDW = "IDW",
    MAE_Media  = "Media", MAE_NN  = "1-NN", MAE_kNN  = "kNN", MAE_IDW  = "IDW"
  ) |>
  tab_options(
    table.font.size      = px(12),
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(12),
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.border.bottom.color = "black",
    column_labels.border.bottom.color = "black"
  ) |>
  opt_horizontal_padding(scale = 2)

ruta_tabla <- here(
  "output", "figures", "interpolacion_clima",
  "tabla_comparacion_metodos.png"
)
tryCatch(
  gtsave(tbl, filename = ruta_tabla, vwidth = 1000),
  error = function(e) {
    warning("gtsave no disponible; se usa gridExtra: ", conditionMessage(e))
    tabla_grob <- gridExtra::tableGrob(
      as.data.frame(tabla_display), rows = NULL,
      theme = gridExtra::ttheme_minimal(base_size = 8)
    )
    titulo_grob <- grid::textGrob(
      "Comparacion de metodos de interpolacion - LOOCV",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    combinado <- gridExtra::arrangeGrob(
      titulo_grob, tabla_grob, ncol = 1,
      heights = grid::unit(c(0.6, 4.8), "in")
    )
    ggplot2::ggsave(
      ruta_tabla, combinado,
      width = 12, height = 5.4, dpi = 150, bg = "white"
    )
  }
)
cat("Tabla comparativa guardada: tabla_comparacion_metodos.png\n")

# --- Gráfico de barras RMSE por variable y método ---
library(tidyr)

df_rmse <- tabla_rmse_loocv[, .(Variable, Media = RMSE_Media, `1-NN` = RMSE_NN,
                                 `kNN (k=5)` = RMSE_kNN, `IDW (p=2)` = RMSE_IDW)]
df_long <- pivot_longer(df_rmse, cols = -Variable,
                        names_to = "Metodo", values_to = "RMSE")
df_long$Metodo <- factor(df_long$Metodo,
                         levels = c("Media", "1-NN", "kNN (k=5)", "IDW (p=2)"))

ggplot(df_long, aes(x = Variable, y = RMSE, fill = Metodo)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", RMSE)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 2.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "RMSE por Variable y Método de Interpolación (LOOCV)",
    subtitle = sprintf("Madrid — %s", format(FECHA_ANALISIS, "%d %b %Y")),
    x = NULL, y = "RMSE", fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1),
    plot.title    = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(here("output", "figures", "interpolacion_clima",
            "grafico_rmse_metodos.png"),
       width = 10, height = 6, dpi = 300, bg = "white")
cat("Gráfico RMSE comparativo guardado: grafico_rmse_metodos.png\n")

# ==============================================================================
# MAPAS DE INTERPOLACIÓN
# ==============================================================================
library(ggplot2)
library(viridis)

# --- 1. Corregir coordenadas corruptas ---
dt_meteo[LATITUD  < 35, LATITUD  := LATITUD  * 10^floor(log10(40 / LATITUD) + 1)]
dt_meteo[LONGITUD > -1, LONGITUD := LONGITUD * 10]

# --- 2. Rejilla espacial sobre Madrid (CRS unificado UTM 25830) ---
cat("Creando la rejilla espacial sobre Madrid...\n")

mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE) |>
  st_make_valid()

# Transformar a UTM ANTES de crear rejilla e intersección (evita CRS mismatch)
mapa_distritos_utm <- st_transform(mapa_distritos, 25830)

rejilla <- st_make_grid(mapa_distritos_utm, n = c(100, 100), what = "centers") |>
  st_as_sf()
# Ambos ya están en 25830 → no hay error de CRS
rejilla_madrid <- st_intersection(rejilla, st_union(mapa_distritos_utm))

coords_rejilla <- st_coordinates(rejilla_madrid)
dx <- diff(sort(unique(coords_rejilla[, 1])))[1]
dy <- diff(sort(unique(coords_rejilla[, 2])))[1]

# Bordes de distritos para superponer en los mapas
bordes_utm    <- st_geometry(mapa_distritos_utm) |> st_cast("MULTILINESTRING") |> st_cast("LINESTRING")
bordes_coords <- as.data.frame(st_coordinates(bordes_utm))

# Usar el mismo día con lluvia y viento
fecha_mapa <- FECHA_ANALISIS
cat("Fecha seleccionada para los mapas:", format(fecha_mapa), "\n")

# --- 3. Bucle: interpolar + mapa comparativo por variable ---
f_mapa <- Z ~ 1

dir.create(here("output", "figures", "interpolacion_clima"), recursive = TRUE, showWarnings = FALSE)

# Incertidumbre empirica del IDW:
# 1) prediccion leave-one-out en cada estacion;
# 2) error cuadratico de cada estacion;
# 3) RMSE local sobre la rejilla usando los mismos pesos IDW (p = 2).
# No es una varianza posterior, sino un diagnostico espacial de error LOOCV.
incertidumbre_idw_loocv <- function(sf_estaciones, coords_rejilla, idp = 2) {
  coords_est <- st_coordinates(sf_estaciones)
  valores <- sf_estaciones$Z
  n_est <- nrow(coords_est)

  residuos_loo <- vapply(seq_len(n_est), function(i) {
    idx_train <- setdiff(seq_len(n_est), i)
    distancias <- sqrt(
      (coords_est[i, 1] - coords_est[idx_train, 1])^2 +
        (coords_est[i, 2] - coords_est[idx_train, 2])^2
    )
    pesos <- 1 / pmax(distancias, 1)^idp
    pred_loo <- sum(pesos * valores[idx_train]) / sum(pesos)
    valores[i] - pred_loo
  }, numeric(1))

  dx_grid <- outer(coords_rejilla[, 1], coords_est[, 1], "-")
  dy_grid <- outer(coords_rejilla[, 2], coords_est[, 2], "-")
  dist_grid <- sqrt(dx_grid^2 + dy_grid^2)
  pesos_grid <- 1 / pmax(dist_grid, 1)^idp
  pesos_grid <- pesos_grid / rowSums(pesos_grid)

  rmse_local <- sqrt(as.vector(pesos_grid %*% (residuos_loo^2)))
  distancia_min_km <- apply(dist_grid, 1, min) / 1000

  list(
    rmse_local = rmse_local,
    rmse_global = sqrt(mean(residuos_loo^2)),
    residuos_loo = residuos_loo,
    distancia_min_km = distancia_min_km
  )
}

for (variable_mapa in mis_variables_clima) {
  variable_label <- unname(etiquetas_variables[variable_mapa])
  
  cat(sprintf("\n--- Generando mapa: %s ---\n", variable_mapa))
  
  dt_instante <- dt_meteo[FECHA == fecha_mapa & !is.na(get(variable_mapa))]
  
  if (nrow(dt_instante) < 3) {
    cat("  Menos de 3 estaciones disponibles. Saltando.\n")
    next
  }
  
  sf_estaciones <- st_as_sf(dt_instante, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
    st_transform(25830)
  sf_estaciones$Z <- sf_estaciones[[variable_mapa]]
  coords_est <- as.data.frame(st_coordinates(sf_estaciones))
  
  # 4 superficies de interpolación
  media_val <- mean(sf_estaciones$Z, na.rm = TRUE)
  mod_nn    <- idw(f_mapa, sf_estaciones, rejilla_madrid, nmax = 1, debug.level = 0)
  k_real    <- min(5, nrow(sf_estaciones))
  mod_knn   <- idw(f_mapa, sf_estaciones, rejilla_madrid, nmax = k_real, idp = 0, debug.level = 0)
  mod_idw   <- idw(f_mapa, sf_estaciones, rejilla_madrid, idp  = 2, debug.level = 0)
  
  df_plot <- rbind(
    data.frame(X = coords_rejilla[, 1], Y = coords_rejilla[, 2],
               Valor = media_val,            Metodo = "1. Media (Baseline)"),
    data.frame(X = coords_rejilla[, 1], Y = coords_rejilla[, 2],
               Valor = mod_nn$var1.pred,     Metodo = "2. Vecino Más Cercano (1-NN)"),
    data.frame(X = coords_rejilla[, 1], Y = coords_rejilla[, 2],
               Valor = mod_knn$var1.pred,    Metodo = "3. kNN (k=5)"),
    data.frame(X = coords_rejilla[, 1], Y = coords_rejilla[, 2],
               Valor = mod_idw$var1.pred,    Metodo = "4. IDW (p=2)")
  )
  df_plot$Metodo <- factor(df_plot$Metodo,
                           levels = c("1. Media (Baseline)", "2. Vecino Más Cercano (1-NN)",
                                      "3. kNN (k=5)", "4. IDW (p=2)"))
  
  p <- ggplot(df_plot, aes(X, Y, fill = Valor)) +
    geom_tile(width = dx, height = dy) +
    geom_path(data = bordes_coords, aes(x = X, y = Y, group = L1),
              color = "white", linewidth = 0.3, inherit.aes = FALSE) +
    geom_point(data = coords_est, aes(X, Y),
               color = "black", size = 1.5, shape = 16, inherit.aes = FALSE) +
    facet_wrap(~ Metodo, ncol = 2) +
    scale_fill_viridis_c(option = "turbo", name = variable_label) +
    labs(
      title    = paste("Interpolación espacial:", variable_label),
      subtitle = paste("Madrid —", format(fecha_mapa, "%d %b %Y"),
                       "| Rejilla 100x100 | Puntos negros = estaciones reales"),
      caption  = "Métodos: Media global, 1-NN, kNN (k=5), IDW (p=2)"
    ) +
    theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      axis.text       = element_blank(),
      axis.title      = element_blank(),
      strip.text      = element_text(face = "bold", size = 11),
      plot.title      = element_text(face = "bold", size = 13),
      legend.position = "right"
    )
  
  print(p)
  
  nombre_archivo <- gsub("[^a-zA-Z0-9]", "_", variable_mapa)
  ggsave(
    here("output", "figures", "interpolacion_clima",
         paste0("interpolacion_", nombre_archivo, ".png")),
    plot = p, width = 10, height = 8, dpi = 200, bg = "white"
  )
  cat(sprintf("  Guardado: interpolacion_%s.png\n", nombre_archivo))

  # Mapa independiente de incertidumbre empirica del IDW.
  incertidumbre <- incertidumbre_idw_loocv(
    sf_estaciones = sf_estaciones,
    coords_rejilla = coords_rejilla,
    idp = 2
  )
  df_incertidumbre <- data.frame(
    X = coords_rejilla[, 1],
    Y = coords_rejilla[, 2],
    Incertidumbre = incertidumbre$rmse_local
  )

  p_incertidumbre <- ggplot(
    df_incertidumbre,
    aes(X, Y, fill = Incertidumbre)
  ) +
    geom_tile(width = dx, height = dy) +
    geom_path(
      data = bordes_coords,
      aes(x = X, y = Y, group = L1),
      color = "white", linewidth = 0.35, inherit.aes = FALSE
    ) +
    geom_point(
      data = coords_est, aes(X, Y),
      color = "black", fill = "white", size = 2,
      shape = 21, stroke = 0.7, inherit.aes = FALSE
    ) +
    scale_fill_viridis_c(
      option = "magma", direction = -1,
      name = paste0("RMSE local\n", variable_label)
    ) +
    coord_equal() +
    labs(
      title = paste("Incertidumbre empirica IDW:", variable_label),
      subtitle = sprintf(
        "Madrid - %s | LOOCV global RMSE = %.3f | %d estaciones",
        format(fecha_mapa, "%d %b %Y"),
        incertidumbre$rmse_global,
        nrow(sf_estaciones)
      ),
      caption = paste(
        "RMSE LOOCV local ponderado con IDW (p=2).",
        "No representa una varianza posterior."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "right"
    )

  ggsave(
    here(
      "output", "figures", "interpolacion_clima",
      paste0("incertidumbre_IDW_", nombre_archivo, ".png")
    ),
    plot = p_incertidumbre,
    width = 9, height = 7, dpi = 300, bg = "white"
  )
  cat(sprintf("  Guardado: incertidumbre_IDW_%s.png\n", nombre_archivo))
}

cat("\nMapas guardados en output/figures/interpolacion_clima/\n")

