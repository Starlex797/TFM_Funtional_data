# ==============================================================================
# Superficie krigeada de la intensidad de tráfico
# ==============================================================================
# Reutiliza el cache de medidores (media anual por medidor) y produce:
#   - Mapa de puntos (medidores)
#   - Superficie de intensidad por kriging ordinario
#   - Mapa de incertidumbre (desviación típica de kriging)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(gstat)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})
sf_use_s2(FALSE)

ANIO <- 2025L
DIR_SALIDA <- here("outputs", "figures", "kriging_trafico")

sensores <- readRDS(here("data", "processed", "Trafico",
                         sprintf("sensores_%d_media.rds", ANIO)))
sensores <- sensores[
  !is.na(utm_x) & !is.na(utm_y) & !is.na(intensidad) &
    utm_x > 400000 & utm_x < 470000 & utm_y > 4450000 & utm_y < 4500000
]
sf_sensores <- st_as_sf(sensores, coords = c("utm_x", "utm_y"), crs = 25830)

mapa_utm <- st_transform(st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE
)), 25830)

# Rejilla de kriging.
rejilla <- st_as_sf(st_make_grid(mapa_utm, n = c(150, 150), what = "centers"))
rejilla_madrid <- st_intersection(rejilla, st_union(mapa_utm))

# Variograma + kriging a la rejilla (vecindario local).
vario_emp <- variogram(intensidad ~ 1, sf_sensores)
vario_mod <- fit.variogram(vario_emp, vgm(c("Exp", "Sph", "Gau")))
krig <- krige(intensidad ~ 1, sf_sensores, rejilla_madrid,
              model = vario_mod, nmax = 40, debug.level = 0)

coords_g <- st_coordinates(rejilla_madrid)
df_k <- data.frame(
  X = coords_g[, 1], Y = coords_g[, 2],
  pred = krig$var1.pred, sd = sqrt(krig$var1.var)
)
df_pts <- as.data.frame(st_coordinates(sf_sensores))
df_pts$intensidad <- sensores$intensidad
bordes <- as.data.frame(st_coordinates(
  st_cast(st_cast(st_geometry(mapa_utm), "MULTILINESTRING"), "LINESTRING")
))

tema <- theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        legend.position = "right", legend.key.size = unit(0.4, "cm"))

p_pts <- ggplot() +
  geom_path(data = bordes, aes(X, Y, group = L1), color = "grey70", linewidth = 0.3) +
  geom_point(data = df_pts, aes(X, Y, color = intensidad), size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "inferno", trans = "sqrt", direction = -1,
                        name = "veh/h", limits = c(0, 1500), oob = scales::squish) +
  coord_equal() + labs(title = "Medidores (puntos)") + tema

p_surf <- ggplot() +
  geom_tile(data = df_k, aes(X, Y, fill = pred)) +
  geom_path(data = bordes, aes(X, Y, group = L1), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", trans = "sqrt", direction = -1,
                       name = "veh/h", limits = c(0, 1500), oob = scales::squish) +
  coord_equal() + labs(title = "Superficie krigeada") + tema

p_sd <- ggplot() +
  geom_tile(data = df_k, aes(X, Y, fill = sd)) +
  geom_path(data = bordes, aes(X, Y, group = L1), color = "white", linewidth = 0.3) +
  geom_point(data = df_pts, aes(X, Y), color = "black", size = 0.15, alpha = 0.5) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "SD\n(veh/h)") +
  coord_equal() + labs(title = "Incertidumbre (desv. típica)") + tema

fig <- (p_pts | p_surf | p_sd) +
  plot_annotation(
    title = "Kriging ordinario de la intensidad de tráfico — media anual 2025",
    subtitle = paste(
      "Superficie continua desde los medidores puntuales.",
      "La incertidumbre crece donde no hay medidores (NW natural)."
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, color = "grey30"))
  )

ggsave(file.path(DIR_SALIDA, "superficie_krigeada_trafico.png"),
       plot = fig, width = 16, height = 6.5, dpi = 200, bg = "white")
cat("Superficie krigeada guardada.\n")
