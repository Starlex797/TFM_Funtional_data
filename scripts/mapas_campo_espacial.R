# ==============================================================================
# MAPAS DE LA SUPERFICIE CONTINUA DE NO2 EN MADRID (Modelo Espacio-Temporal)
# ==============================================================================
# Genera tres tipos de mapas sobre una rejilla regular (300 × 300 celdas):
#
#   Mapa 1 — Media Posterior del NO2 (µg/m³)
#             Superficie predicha anual: exp(intercepto + campo_espacial_medio).
#             Muestra el patrón de contaminación a lo largo de Madrid.
#
#   Mapa 2 — Incertidumbre (Desviación Estándar Posterior)
#             SD del campo latente promediado sobre el año.
#             Donde es alta, el modelo necesita más datos para ser preciso.
#
#   Mapa 3 — Probabilidad de Excedencia del umbral regulatorio (40 µg/m³)
#             P(NO2 > 40) en cada celda, usando la distribución normal posterior.
#             Umbral: valor límite anual de la UE para NO2 (Directiva 2008/50/CE).
#
# CRS de trabajo: UTM Zona 30N en km (EPSG:25830 / 1000)
# ==============================================================================

library(INLA)
library(data.table)
library(sf)
library(ggplot2)
library(dplyr)
library(here)

# ------------------------------------------------------------------------------
# 0. CARGAR OBJETOS NECESARIOS
# ------------------------------------------------------------------------------
modelo_st    <- readRDS(here("data", "processed", "modelo_final_no2_madrid.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid.rds"))
dt_maestro   <- readRDS(here("data", "processed", "dataset_maestro_inla_2025.rds"))
setDT(dt_maestro)

# Mapa de distritos en km (mismo CRS que la malla)
mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"),
                           quiet = TRUE)
mapa_distritos_km <- st_transform(mapa_distritos, 25830)
st_geometry(mapa_distritos_km) <- st_geometry(mapa_distritos_km) / 1000
madrid_union <- st_union(mapa_distritos_km)

# Extraer bordes como data.frame XY para geom_path (evita conflicto geom_sf / coord_equal)
bordes_distritos <- st_cast(mapa_distritos_km, "POLYGON")
bordes_coords    <- as.data.frame(st_coordinates(bordes_distritos))
# Columnas: X, Y, L1 (anillo), L2 (polígono) → group = interaction(L1, L2)

# Coordenadas de las estaciones de medición (para superponer en los mapas)
estaciones_km <- unique(dt_maestro[, .(ESTACION, LONGITUD, LATITUD)]) |>
  st_as_sf(coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)
st_geometry(estaciones_km) <- st_geometry(estaciones_km) / 1000

# ------------------------------------------------------------------------------
# 1. EXTRAER Y PROMEDIAR EL CAMPO LATENTE SOBRE EL TIEMPO
# ------------------------------------------------------------------------------
n_nodos <- malla_madrid$n
ndays   <- length(unique(dt_maestro$ID_TIEMPO))

cat("Nodos de la malla:", n_nodos, "| Días:", ndays, "\n")
cat("Longitud de summary.random$campo_espacial:",
    nrow(modelo_st$summary.random$campo_espacial), "\n")

# Reshape a matriz (nodos × días)
campo_mean_mat <- matrix(
  modelo_st$summary.random$campo_espacial$mean,
  nrow = n_nodos, ncol = ndays
)
campo_sd_mat <- matrix(
  modelo_st$summary.random$campo_espacial$sd,
  nrow = n_nodos, ncol = ndays
)

# Media y SD anuales en cada nodo de la malla
campo_mean_anual <- rowMeans(campo_mean_mat)   # Promedio temporal
campo_sd_anual   <- rowMeans(campo_sd_mat)     # Incertidumbre media

# ------------------------------------------------------------------------------
# 2. PROYECTAR SOBRE REJILLA REGULAR (300 × 300)
# ------------------------------------------------------------------------------
proyector <- inla.mesh.projector(malla_madrid, dims = c(300, 300))

grid_campo_medio <- inla.mesh.project(proyector, campo_mean_anual)
grid_campo_sd    <- inla.mesh.project(proyector, campo_sd_anual)

# Intercepto posterior para obtener predicción en escala log(NO2)
beta0 <- modelo_st$summary.fixed["intercept", "mean"]

# Predicción en escala original: exp(intercepto + campo_espacial)
grid_no2_log <- beta0 + grid_campo_medio         # Log(NO2) predicho
grid_no2_ugm3 <- exp(grid_no2_log)               # NO2 en µg/m³

# Probabilidad de excedencia: P(log(NO2) > log(40)) en cada celda
# Usando aprox. normal: N(media_campo + beta0, sd_campo^2)
umbral_log <- log(40)   # 40 µg/m³ = límite anual UE (Directiva 2008/50/CE)
grid_exceedance <- pnorm(umbral_log,
                          mean = grid_no2_log,
                          sd   = grid_campo_sd,
                          lower.tail = FALSE)   # P(X > umbral)

# ------------------------------------------------------------------------------
# 3. ENSAMBLAR DATA FRAME Y ENMASCARAR CON EL POLÍGONO DE MADRID
# ------------------------------------------------------------------------------
grid_df <- expand.grid(
  X_km = proyector$x,
  Y_km = proyector$y
) |>
  mutate(
    no2_ugm3   = as.vector(grid_no2_ugm3),
    campo_sd   = as.vector(grid_campo_sd),
    exceedance = as.vector(grid_exceedance)
  )

# Convertir a sf para comprobar intersección con Madrid
grid_sf <- st_as_sf(
  grid_df |> filter(!is.na(no2_ugm3)),
  coords = c("X_km", "Y_km"),
  crs = st_crs(mapa_distritos_km)
)
dentro <- st_within(grid_sf, madrid_union, sparse = FALSE)[, 1]

# Poner NA fuera de Madrid
idx_validos <- which(!is.na(grid_df$no2_ugm3))
grid_df$no2_ugm3[idx_validos[!dentro]]   <- NA
grid_df$campo_sd[idx_validos[!dentro]]   <- NA
grid_df$exceedance[idx_validos[!dentro]] <- NA

# Eliminamos celdas vacías para acelerar el renderizado
grid_df <- grid_df |> filter(!is.na(no2_ugm3))

cat("Celdas dentro de Madrid:", nrow(grid_df), "\n")

# ------------------------------------------------------------------------------
# 4. TEMA COMÚN PARA LOS TRES MAPAS
# ------------------------------------------------------------------------------
tema_mapa <- theme_minimal() +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_text(size = 7, color = "grey50"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    legend.position  = "right",
    legend.key.width = unit(0.4, "cm"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

# Coordenadas de estaciones como data.frame (para geom_point)
est_coords <- as.data.frame(st_coordinates(estaciones_km))
colnames(est_coords) <- c("X_km", "Y_km")

# ------------------------------------------------------------------------------
# MAPA 1: MEDIA POSTERIOR DEL NO2 (µg/m³)
# ------------------------------------------------------------------------------
cat("Generando Mapa 1: Media Posterior del NO2...\n")

mapa1 <- ggplot() +
  geom_raster(data = grid_df, aes(x = X_km, y = Y_km, fill = no2_ugm3)) +
  scale_fill_gradientn(
    colours  = c("#2166ac", "#74add1", "#e0f3f8", "#fee090", "#f46d43", "#a50026"),
    name     = "NO₂ (µg/m³)",
    na.value = "transparent",
    guide = guide_colorbar(barheight = 10)
  ) +
  geom_path(data = bordes_coords,
            aes(x = X, y = Y, group = interaction(L1, L2)),
            colour = "white", linewidth = 0.4) +
  geom_point(data = est_coords, aes(x = X_km, y = Y_km),
             shape = 21, fill = "white", colour = "black", size = 1.8) +
  labs(
    title    = "Mapa 1: Superficie Predicha de NO₂ en Madrid (2025)",
    subtitle = "Media posterior anual — exp(intercepto + campo espacial latente)"
  ) +
  coord_equal() +
  tema_mapa

# ------------------------------------------------------------------------------
# MAPA 2: INCERTIDUMBRE (DESVIACIÓN ESTÁNDAR POSTERIOR)
# ------------------------------------------------------------------------------
cat("Generando Mapa 2: Incertidumbre del campo latente...\n")

mapa2 <- ggplot() +
  geom_raster(data = grid_df, aes(x = X_km, y = Y_km, fill = campo_sd)) +
  scale_fill_gradientn(
    colours  = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    name     = "SD Posterior\n(escala log)",
    na.value = "transparent",
    guide = guide_colorbar(barheight = 10)
  ) +
  geom_path(data = bordes_coords,
            aes(x = X, y = Y, group = interaction(L1, L2)),
            colour = "white", linewidth = 0.4) +
  geom_point(data = est_coords, aes(x = X_km, y = Y_km),
             shape = 21, fill = "white", colour = "black", size = 1.8) +
  labs(
    title    = "Mapa 2: Incertidumbre de la Estimación (SD Posterior)",
    subtitle = "Alta SD → el modelo tiene menos información en esa zona"
  ) +
  coord_equal() +
  tema_mapa

# ------------------------------------------------------------------------------
# MAPA 3: PROBABILIDAD DE EXCEDENCIA DEL LÍMITE REGULATORIO (40 µg/m³)
# ------------------------------------------------------------------------------
cat("Generando Mapa 3: Probabilidad de excedencia (> 40 µg/m³)...\n")

mapa3 <- ggplot() +
  geom_raster(data = grid_df, aes(x = X_km, y = Y_km, fill = exceedance)) +
  scale_fill_gradientn(
    colours  = c("#1a9850", "#a6d96a", "#ffffbf", "#fdae61", "#d73027"),
    limits   = c(0, 1),
    labels   = scales::percent_format(accuracy = 1),
    name     = "P(NO₂ > 40)",
    na.value = "transparent",
    guide = guide_colorbar(barheight = 10)
  ) +
  geom_sf(data = mapa_distritos_km, fill = NA, colour = "white",
          linewidth = 0.4, alpha = 0.8) +
  geom_point(data = est_coords, aes(x = X_km, y = Y_km),
             shape = 21, fill = "white", colour = "black", size = 1.8) +
  labs(
    title    = "Mapa 3: Probabilidad de Superar el Límite Anual de la UE",
    subtitle = "Umbral: 40 µg/m³ (Directiva 2008/50/CE) — Rojo = mayor riesgo"
  ) +
  coord_equal() +
  tema_mapa

# ------------------------------------------------------------------------------
# 5. GUARDAR LOS TRES MAPAS
# ------------------------------------------------------------------------------
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(here("output", "figures", "mapa1_media_posterior_no2.png"),
       plot = mapa1, width = 7, height = 6, dpi = 300, bg = "white")

ggsave(here("output", "figures", "mapa2_incertidumbre_sd.png"),
       plot = mapa2, width = 7, height = 6, dpi = 300, bg = "white")

ggsave(here("output", "figures", "mapa3_exceedance_40ugm3.png"),
       plot = mapa3, width = 7, height = 6, dpi = 300, bg = "white")

cat("\n✅ Los tres mapas se han guardado en output/figures/\n")
cat("   - mapa1_media_posterior_no2.png\n")
cat("   - mapa2_incertidumbre_sd.png\n")
cat("   - mapa3_exceedance_40ugm3.png\n")

# ------------------------------------------------------------------------------
# (OPCIONAL) Vista interactiva con mapview
# ------------------------------------------------------------------------------
# Descomenta este bloque para explorar los mapas de forma interactiva en RStudio.
# Requiere: install.packages("mapview")
#
library(mapview)
library(stars)
#
grid_stars <- grid_df |>
   select(X_km, Y_km, no2_ugm3) |>
   filter(!is.na(no2_ugm3)) |>
   st_as_sf(coords = c("X_km", "Y_km"), crs = st_crs(mapa_distritos_km)) |>
   st_rasterize()

mapview(grid_stars, layer.name = "NO2 medio (µg/m³)", col.regions = hcl.colors(100, "RdYlBu", rev = TRUE))

