#===============================================================================
# Interpolation climate variables 
#===============================================================================

# ==============================================================================
# EJECUCIÓN
# ==============================================================================
source(here("R","FUNCIONES_INTERPOLACION.R"))


mis_variables_clima <- c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                         "Presion Barométrica", "Radiación Solar", "Velocidad Viento")
dt_meteo <- readRDS(here("data", "processed", "Clima", "diario", "meteo_madrid_2025_diario.rds"))

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

for (variable_mapa in mis_variables_clima) {
  
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
    scale_fill_viridis_c(option = "turbo", name = variable_mapa) +
    labs(
      title    = paste("Interpolación espacial:", variable_mapa),
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
}

cat("\nMapas guardados en output/figures/interpolacion_clima/\n")