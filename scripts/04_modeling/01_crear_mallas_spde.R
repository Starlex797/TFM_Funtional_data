# ==============================================================================
# STEP 1: CREATION OF THE SPATIAL MESH FOR NO2  (MESH)
# ==============================================================================

library(INLA)
library(sf)
library(data.table)
library(here)
library(ggplot2) # se necesita ya en el Bloque 3 para el diagnóstico de boundaries

# 1. We load the NO2 data for 2025 (daily) and convert it to a data.table for efficient processing.
dt_no2_2025 <- readRDS(here("data", "processed", "Contaminacion", "diario", "aire_madrid_2025_No2_trans_diarios1.rds"))
dt_clima_2025 <- readRDS(here("data", "processed", "Clima", "diario", "meteo_madrid_2025_diario5.rds"))

Variable_estudio <- "NO2" # Variable de estudio: NO2, Radiacion_Solar, etc.

setDT(dt_no2_2025)
setDT(dt_clima_2025)

# 2. extract unique coordinates and project them to UTM 30N (EPSG:25830) in kilometers.
# Se extraen las coordenadas únicas y se proyectan a UTM 30N (EPSG:25830) en kilómetros.
coords_estaciones <- unique(dt_no2_2025[, .(ESTACION, LONGITUD, LATITUD)])

# !is.na(Variable_estudio) comprueba el TEXTO "Radiacion_Solar" (nunca es NA):
# no filtraba nada. get(Variable_estudio) recupera la columna que ese texto
# nombra, para quedarnos solo con las estaciones que sí miden esa variable.
# coords_estaciones_clima <- unique(
# dt_clima_2025[!is.na(get(Variable_estudio)), .(ESTACION, LONGITUD, LATITUD)]
# )

if (Variable_estudio %in% c("NO2", "DATO_DIARIO", "LOG_NO2_DIARIO")) {
  datos_a_proyectar <- coords_estaciones
} else {
  if (!Variable_estudio %in% names(dt_clima_2025)) {
    stop("La variable '", Variable_estudio, "' no existe en dt_clima_2025.")
  }
  datos_a_proyectar <- unique(
    dt_clima_2025[!is.na(get(Variable_estudio)), .(ESTACION, LONGITUD, LATITUD)]
  )
}

coords_sf <- st_as_sf(datos_a_proyectar, coords = c("LONGITUD", "LATITUD"), crs = 4326) # CRS original en WGS84 (grados)
coords_utm <- st_transform(coords_sf, 25830) # proyection in utm 30N (EPSG:25830)

# Matrix of coordinates in kilometers (for INLA)
coords_matriz <- st_coordinates(coords_utm) / 1000
carpeta_figuras_mallas <- here("outputs", "figures", "modelo", "mallas", Variable_estudio)
if (!dir.exists(carpeta_figuras_mallas)) dir.create(carpeta_figuras_mallas, recursive = TRUE)

# Map of Madrid districts in UTM 30N (EPSG:25830) in kilometers
mapa_distritos_km <- st_transform(
  st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE),
  25830
)
# 3.Define the boundaries (Boundaries)

bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.3, resolution = 50) # Contorno interno ajustado, con resolución fija para que no se autointerseque.
bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.6, resolution = 50) # Contorno externo de amortiguación, para evitar efectos de borde.

# 4. Explore different mesh resolutions (coarse, medium, fine) and save them to disk.

# Cutoff por variable climática = distancia entre las DOS estaciones más
distancia_minima_variable <- function(variable) {
  coords_var <- unique(
    dt_clima_2025[!is.na(get(variable)), .(ESTACION, LONGITUD, LATITUD)]
  )
  if (nrow(coords_var) < 2) {
    stop("Menos de 2 estaciones miden '", variable, "': no se puede calcular un cutoff.")
  }
  xy_km <- st_coordinates(st_transform(
    st_as_sf(coords_var, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830
  )) / 1000
  min(dist(xy_km))
}

# NO2
edge_gruesa_NO2 <- c(8, 10)
edge_media_NO2 <- c(3, 8)
edge_fina_NO2 <- c(1.86, 10)
cutoff_NO2 <- 0.79 # Cut off to avoid creating very small triangles between nearby stations. This merges close points into one, preventing tiny triangles.

# Radiacion_Solar
edge_gruesa_RS <- c(5, 6)
edge_media_RS <- c(3, 6)
edge_fina_RS <- c(1.86, 6)
cutoff_RS <- distancia_minima_variable("Radiacion_Solar") # distancia mínima entre las 2 estaciones más cercanas que miden radiación solar
# Precipitacion
edge_gruesa_Prec <- c(5, 6)
edge_media_Prec <- c(3, 6)
edge_fina_Prec <- c(1.86, 6)
cutoff_Prec <- distancia_minima_variable("Precipitaciones") # ídem, precipitaciones
# Velocidad_Viento
edge_gruesa_VV <- c(5, 6)
edge_media_VV <- c(3, 6)
edge_fina_VV <- c(1.86, 86)
cutoff_VV <- distancia_minima_variable("Velocidad_Viento") # ídem, velocidad del viento
# Presion_Barometrica
edge_gruesa_PB <- c(5, 6)
edge_media_PB <- c(3, 6)
edge_fina_PB <- c(1.86, 6)
cutoff_PB <- distancia_minima_variable("Presion_Barometrica") # ídem, presión barométrica

cat(sprintf(
  "\nCutoff calculado por variable (km): RS=%.3f | Prec=%.3f | VV=%.3f | PB=%.3f\n",
  cutoff_RS, cutoff_Prec, cutoff_VV, cutoff_PB
))

# Selector: coge los parámetros del bloque que corresponde a Variable_estudio
# (arriba) y los deja con nombre genérico. Antes las mallas usaban siempre
# edge_*_NO2/cutoff_NO2 sin mirar Variable_estudio -> los bloques _RS/_Prec/
# _VV/_PB quedaban definidos pero nunca se usaban.
sufijo_variable <- c(
  NO2 = "NO2",
  DATO_DIARIO = "NO2",
  LOG_NO2_DIARIO = "NO2",
  Radiacion_Solar = "RS",
  Precipitaciones = "Prec",
  Velocidad_Viento = "VV",
  Presion_Barometrica = "PB"
)
sufijo <- sufijo_variable[[Variable_estudio]]
if (is.null(sufijo)) {
  stop("No hay parámetros de malla (edge/cutoff) definidos para '", Variable_estudio, "'.")
}
edge_gruesa <- get(paste0("edge_gruesa_", sufijo))
edge_media <- get(paste0("edge_media_", sufijo))
edge_fina <- get(paste0("edge_fina_", sufijo))
cutoff <- get(paste0("cutoff_", sufijo))

# Opcion A: coarse mesh (low resolution, faster)
malla_gruesa <- inla.mesh.2d(
  loc = coords_matriz, # The coordinates of the stations are used to define the mesh nodes. The mesh will be constructed around these points, ensuring that the spatial model can capture the variability in NO2 levels across the study area.
  boundary = list(bnd_inner, bnd_outer), # we define the inner and outer boundaries of the mesh. The inner boundary is a tighter fit around the stations, while the outer boundary provides a buffer zone to avoid edge effects in the spatial model.
  max.edge = edge_gruesa, # define the spatial range as the maximum distance between mesh nodes. Smaller values create a finer mesh with more nodes, while larger values create a coarser mesh with fewer nodes.
  cutoff = cutoff # Cut off to avoid creating very small triangles between nearby stations. This merges close points into one, preventing tiny triangles.
)

# Option b: medium mesh (medium resolution, slower)
malla_media <- inla.mesh.2d(
  loc = coords_matriz,
  boundary = list(bnd_inner, bnd_outer),
  max.edge = edge_media,
  cutoff = cutoff
)

# Option c: fine mesh (high resolution, slowest)
malla_fina <- inla.mesh.2d(
  loc = coords_matriz,
  boundary = list(bnd_inner, bnd_outer),
  max.edge = edge_fina,
  cutoff = cutoff
)


# 5. Save the meshes to disk for later use in modeling. The meshes are saved as RDS files in the "data/processed/Malla" directory. The directory is created if it does not exist. Each mesh is saved with a descriptive filename indicating its resolution (coarse, medium, fine).
# saveRDS() no crea directorios: hay que asegurar la carpeta antes de guardar.
carpeta_mallas <- here("data", "processed", "Malla", Variable_estudio)
dir.create(carpeta_mallas, recursive = TRUE, showWarnings = FALSE)

# carpeta_mallas ya incluye Variable_estudio (línea de arriba): no hace falta
# repetirlo en el nombre de fichero. La versión anterior metía Variable_estudio
# como un segmento de RUTA en vez de parte del nombre, así que intentaba
# guardar dentro de una subcarpeta que no existe (saveRDS no la crea) y habría
# fallado con "cannot open the connection".
saveRDS(malla_gruesa, file.path(carpeta_mallas, "malla_spde_madrid_gruesa.rds"))
saveRDS(malla_media, file.path(carpeta_mallas, "malla_spde_madrid_media.rds"))
saveRDS(malla_fina, file.path(carpeta_mallas, "malla_spde_madrid_fina.rds"))


# Visualization of the number of vertices in each mesh. This provides a quick overview of the complexity of each mesh, which can impact computational time and model performance.
cat("Vértices Malla Gruesa:", malla_gruesa$n, "\n")
cat("Vértices Malla Media:", malla_media$n, "\n")
cat("Vértices Malla Fina:", malla_fina$n, "\n")

# ==============================================================================
# VISUALIZACIÓN: ESTACIONES SOBRE LAS MALLAS SPDE
# ==============================================================================
library(ggplot2)
library(fmesher)
library(gridExtra)

carpeta_figuras_mallas <- here("outputs", "figures", "modelo", "mallas", Variable_estudio)
if (!dir.exists(carpeta_figuras_mallas)) dir.create(carpeta_figuras_mallas, recursive = TRUE)

# Mapa de distritos en km (mismo CRS que coords_matriz)
st_geometry(mapa_distritos_km) <- st_geometry(mapa_distritos_km) / 1000

# Estaciones como data.frame en km. coords_matriz sale de datos_a_proyectar
# (las estaciones de Variable_estudio, filtradas m\u00e1s arriba), as\u00ed que las
# etiquetas tienen que salir de esas mismas estaciones y no de coords_estaciones
# (que son las de NO2: otro conjunto, con otro recuento -> longitudes distintas).
coords_df <- as.data.frame(coords_matriz)
colnames(coords_df) <- c("X", "Y")
coords_df$ESTACION <- datos_a_proyectar$ESTACION

# Paleta y etiquetas compartidas
etiqueta_estacion <- paste(gsub("_", " ", Variable_estudio), "station")
pal_leyenda <- setNames(
  c("#d73027", "#2166ac", "#8c510a"), # brown (BrBG), contrasts with the mesh blue
  c(etiqueta_estacion, "SPDE mesh", "Madrid districts")
)

# Función reutilizable: una malla + estaciones superpuestas + leyenda
crear_mapa_malla <- function(malla, titulo, subtitulo = NULL) {
  # Datos vacíos para registrar "SPDE mesh" en la leyenda (geom_fm no lo hace)
  df_dummy <- data.frame(x = NA_real_, y = NA_real_)

  ggplot() +
    geom_sf(
      data = mapa_distritos_km, fill = NA,
      aes(color = "Madrid districts"), linewidth = 0.4
    ) +
    # geom_fm no registra aes(color) en la escala -> color fijo
    geom_fm(data = malla, color = "#2166ac", linewidth = 0.2, alpha = 0.45) +
    # Capa invisible que crea la entrada de la malla en la leyenda
    geom_line(data = df_dummy, aes(x = x, y = y, color = "SPDE mesh")) +
    geom_point(
      data = coords_df, aes(
        x = X, y = Y,
        color = etiqueta_estacion
      ),
      shape = 17, size = 2.5
    ) +
    scale_color_manual(name = NULL, values = pal_leyenda) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape     = c(NA, 17, NA),
          linetype  = c(1, 0, 1),
          linewidth = c(0.6, 0, 0.6),
          size      = c(NA, 2.5, NA)
        )
      )
    ) +
    labs(
      title = titulo,
      subtitle = subtitulo %||% sprintf("%d vertices | UTM 30N (km)", malla$n),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "white", color = "gray85"),
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(color = "gray40", size = 9),
      axis.text        = element_text(size = 7, color = "gray50"),
      legend.position  = "bottom",
      legend.text      = element_text(size = 9)
    )
}

# Operador nulo-coalescencia (disponible en R ≥ 4.4, pero lo definimos por si acaso)
`%||%` <- function(a, b) if (!is.null(a)) a else b

mapa_gruesa <- crear_mapa_malla(malla_gruesa, sprintf("Coarse mesh  (max.edge = %g km)", edge_gruesa[1]))
mapa_media <- crear_mapa_malla(malla_media, sprintf("Medium mesh  (max.edge = %g km)", edge_media[1]))
mapa_fina <- crear_mapa_malla(malla_fina, sprintf("Fine mesh    (max.edge = %g km)", edge_fina[1]))

# Guardado con verificación.
#
guardar_png <- function(nombre, grafico, width, height) {
  ruta <- file.path(carpeta_figuras_mallas, nombre)
  mtime_previo <- if (file.exists(ruta)) file.mtime(ruta) else as.POSIXct(NA)

  ggsave(ruta, grafico, width = width, height = height, dpi = 200, bg = "white")

  escrito <- file.exists(ruta) &&
    (is.na(mtime_previo) || file.mtime(ruta) > mtime_previo)
  if (!escrito) {
    stop(
      "No se pudo reescribir '", nombre, "'. Suele ser porque el PNG está ",
      "abierto en un visor o en la vista previa del IDE: ciérralo y vuelve a ",
      "ejecutar. Ruta: ", ruta
    )
  }
  cat("  guardado:", nombre, "\n")
  invisible(ruta)
}

# Guardar cada malla individualmente
guardar_png("malla_gruesa.png", mapa_gruesa, width = 7, height = 6)
guardar_png("malla_media.png", mapa_media, width = 7, height = 6)
guardar_png("malla_fina.png", mapa_fina, width = 7, height = 6)

# Guardar mapa comparativo (3 paneles)
comparacion <- arrangeGrob(mapa_gruesa, mapa_media, mapa_fina, ncol = 3)
guardar_png("comparacion_mallas.png", comparacion, width = 18, height = 6)

cat("Mapas guardados en:", carpeta_figuras_mallas, "\n")
grid.arrange(mapa_gruesa, mapa_media, mapa_fina, ncol = 3)
