# ==============================================================================
# Kriging ordinario del tráfico vs agregación a áreas (barrio)
# ==============================================================================
# Los medidores de tráfico son datos PUNTUALES (id + utm_x/utm_y). En el
# preprocesamiento se transformaron a datos de ÁREA (media por barrio). Aquí:
#   1. Se visualizan las localizaciones de los medidores.
#   2. Se compara, por LOOCV, si predice mejor el KRIGING ORDINARIO (puntual) o
#      la media por barrio (área) la intensidad de un medidor oculto.
# La intensidad se agrega a media anual 2025 por medidor (patrón espacial estable).
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(gstat)
  library(ggplot2)
  library(viridis)
})
sf_use_s2(FALSE)

ANIO <- 2025L
DIR_SALIDA <- here("outputs", "figures", "kriging_trafico")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
RUTA_CACHE <- here("data", "processed", "Trafico",
                   sprintf("sensores_%d_media.rds", ANIO))

# ------------------------------------------------------------------------------
# 1. Agregar intensidad media anual por medidor (con cache)
# ------------------------------------------------------------------------------

if (file.exists(RUTA_CACHE)) {
  sensores <- readRDS(RUTA_CACHE)
  cat("Sensores leídos de cache:", nrow(sensores), "\n")
} else {
  dir_horario <- here("data", "raw", "Datos_trafico",
                      sprintf("Datos_limpios_%d", ANIO), "Horario")
  fs <- list.files(dir_horario, pattern = "\\.rds$", full.names = TRUE)
  cat("Agregando", length(fs), "ficheros mensuales...\n")

  acc <- vector("list", length(fs))
  for (i in seq_along(fs)) {
    d <- as.data.table(readRDS(fs[i]))
    acc[[i]] <- d[!is.na(intensidad),
      .(s = sum(as.numeric(intensidad)), n = .N),
      by = .(id, utm_x, utm_y, barrio, distrito)
    ]
    rm(d)
    gc()
    cat("  ", i, "/", length(fs), basename(fs[i]), "\n")
  }
  agg <- rbindlist(acc)
  sensores <- agg[
    , .(
      intensidad = sum(s) / sum(n), n_obs = sum(n),
      utm_x = first(utm_x), utm_y = first(utm_y),
      barrio = first(barrio), distrito = first(distrito)
    ),
    by = id
  ]
  saveRDS(sensores, RUTA_CACHE)
  cat("Sensores agregados y cacheados:", nrow(sensores), "\n")
}

# Limpieza: coordenadas válidas dentro de la zona UTM de Madrid.
sensores <- sensores[
  !is.na(utm_x) & !is.na(utm_y) & !is.na(intensidad) &
    utm_x > 400000 & utm_x < 470000 & utm_y > 4450000 & utm_y < 4500000
]
cat(sprintf(
  "Medidores válidos: %d | intensidad media: %.0f | mediana: %.0f\n",
  nrow(sensores), mean(sensores$intensidad), median(sensores$intensidad)
))

sf_sensores <- st_as_sf(sensores, coords = c("utm_x", "utm_y"), crs = 25830)

# ------------------------------------------------------------------------------
# 2. Mapa de localizaciones de los medidores
# ------------------------------------------------------------------------------

mapa_utm <- st_transform(st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE
)), 25830)

p_loc <- ggplot() +
  geom_sf(data = mapa_utm, fill = "grey96", color = "grey70", linewidth = 0.3) +
  geom_sf(data = sf_sensores, aes(color = intensidad), size = 0.7, alpha = 0.8) +
  scale_color_viridis_c(
    option = "inferno", trans = "sqrt", direction = -1,
    name = "Intensidad\nmedia (veh/h)"
  ) +
  labs(
    title = sprintf("Medidores de tráfico de Madrid (%d) — %d medidores",
                    ANIO, nrow(sensores)),
    subtitle = "Intensidad media anual por medidor (escala sqrt)",
    caption = "Puntos = detectores de tráfico | polígonos = distritos"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey35"),
    legend.position = "right"
  )

ggsave(file.path(DIR_SALIDA, "localizacion_medidores_trafico.png"),
       plot = p_loc, width = 9, height = 8, dpi = 200, bg = "white")
cat("Mapa de localizaciones guardado.\n")
