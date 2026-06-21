# ==============================================================================
# PASO 1: CREACIÓN DE LA MALLA ESPACIAL (MESH)
# ==============================================================================

library(INLA)
library(sf)
library(data.table)
library(here)

# 1. Cargar datos de NO2 (variable respuesta)
dt_no2_2025 <- readRDS(here("data", "processed", "aire_madrid_2025_No2_trans_diarios.rds"))
setDT(dt_no2_2025)

# 2. Extraer y proyectar coordenadas espaciales
# Se extraen las coordenadas únicas y se proyectan a UTM 30N (EPSG:25830) en kilómetros.
coords_estaciones <- unique(dt_no2_2025[, .(ESTACION, LONGITUD, LATITUD)])
coords_sf <- st_as_sf(coords_estaciones, coords = c("LONGITUD", "LATITUD"), crs = 4326) # CRS original en WGS84 (grados)
coords_utm <- st_transform(coords_sf, 25830) # Proyección a UTM zona 30N (metros)

# Matriz de coordenadas en kilómetros (requerido por INLA para estabilidad numérica)
coords_matriz <- st_coordinates(coords_utm) / 1000 

# 3. Definir los límites geométricos (Boundaries)
# Se genera un contorno interno ajustado a las estaciones y un contorno externo de amortiguación.
# La función convex define cómo se dibuja la línea de la frontera alrededor de las estaciones.
bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.05, resolution = 50) # Se deja un 5% de margen adicional para incluir todas las estaciones dentro del contorno. Esto ayuda a asegurar que el modelo espacial tenga en cuenta todas las ubicaciones de las estaciones sin ser demasiado ajustado.
bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.2) # Se deja un 20% de margen adicional para amortiguación fuera de las estaciones. Esto ayuda a evitar problemas de borde en el modelado espacial.

# 4. Construir las Mallas (Meshes) de diferentes resoluciones

# OPCIÓN A: Malla Gruesa (Baja resolución, muy rápida)
# Triángulos internos de hasta 4 km, externos de 10 km
malla_gruesa <- inla.mesh.2d(
  loc = coords_matriz, 
  boundary = list(bnd_inner, bnd_outer), 
  max.edge = c(8, 12),  
  cutoff = 0.5         
)

# OPCIÓN B: Malla Media (Tu configuración original)
# Triángulos internos de hasta 2 km, externos de 5 km
malla_media <- inla.mesh.2d(
  loc = coords_matriz, 
  boundary = list(bnd_inner, bnd_outer), 
  max.edge = c(4, 8),  
  cutoff = 0.5         
)

# OPCIÓN C: Malla Fina (Alta resolución, más lenta)
# Triángulos internos de hasta 1 km, externos de 3 km
# Bajamos el cutoff a 0.25 para permitir que los puntos estén más cerca sin agruparse
malla_fina <- inla.mesh.2d(
  loc = coords_matriz, 
  boundary = list(bnd_inner, bnd_outer), 
  max.edge = c(1, 4),  
  cutoff = 0.25         
)

# 5. Exportar los tres archivos al disco duro
saveRDS(malla_gruesa, here("data", "processed", "malla_spde_madrid_gruesa.rds"))
saveRDS(malla_media, here("data", "processed", "malla_spde_madrid_media.rds"))
saveRDS(malla_fina, here("data", "processed", "malla_spde_madrid_fina.rds"))

# Visualizar cuántos vértices (nodos) tiene cada una para compararlas en el TFM
cat("Vértices Malla Gruesa:", malla_gruesa$n, "\n")
cat("Vértices Malla Media:", malla_media$n, "\n")
cat("Vértices Malla Fina:", malla_fina$n, "\n")

# ==============================================================================
# VISUALIZACIÓN: ESTACIONES DE NO2 SOBRE LAS MALLAS SPDE
# ==============================================================================
library(ggplot2)
library(fmesher)
library(gridExtra)

carpeta_figuras_mallas <- here("output", "figures", "mallas")
if (!dir.exists(carpeta_figuras_mallas)) dir.create(carpeta_figuras_mallas, recursive = TRUE)

# Mapa de distritos en km (mismo CRS que coords_matriz)
mapa_distritos_km <- st_transform(
  st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE),
  25830
)
st_geometry(mapa_distritos_km) <- st_geometry(mapa_distritos_km) / 1000

# Estaciones como data.frame en km
coords_df <- as.data.frame(coords_matriz)
colnames(coords_df) <- c("X", "Y")
coords_df$ESTACION <- coords_estaciones$ESTACION

# Paleta y etiquetas compartidas
pal_leyenda <- c(
  "Estaci\u00f3n NO\u2082" = "#d73027",
  "Malla SPDE"             = "#2166ac",
  "Distritos de Madrid"    = "gray50"
)

# Función reutilizable: una malla + estaciones superpuestas + leyenda
crear_mapa_malla <- function(malla, titulo, subtitulo = NULL) {

  # Datos vacíos para registrar "Malla SPDE" en la leyenda (geom_fm no lo hace)
  df_dummy <- data.frame(x = NA_real_, y = NA_real_)

  ggplot() +
    geom_sf(data = mapa_distritos_km, fill = NA,
            aes(color = "Distritos de Madrid"), linewidth = 0.4) +
    # geom_fm no registra aes(color) en la escala -> color fijo
    geom_fm(data = malla, color = "#2166ac", linewidth = 0.2, alpha = 0.45) +
    # Capa invisible que crea la entrada de la malla en la leyenda
    geom_line(data = df_dummy, aes(x = x, y = y, color = "Malla SPDE")) +
    geom_point(data = coords_df, aes(x = X, y = Y,
               color = "Estaci\u00f3n NO\u2082"),
               shape = 17, size = 2.5) +
    scale_color_manual(name = NULL, values = pal_leyenda) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape     = c(NA,   17,   NA),
          linetype  = c(1,    0,    1),
          linewidth = c(0.6,  0,    0.6),
          size      = c(NA,   2.5,  NA)
        )
      )
    ) +
    labs(
      title    = titulo,
      subtitle = subtitulo %||% sprintf("%d v\u00e9rtices | UTM 30N (km)", malla$n),
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

mapa_gruesa <- crear_mapa_malla(malla_gruesa, "Malla Gruesa  (max.edge = 8 km)")
mapa_media  <- crear_mapa_malla(malla_media,  "Malla Media   (max.edge = 4 km)")
mapa_fina   <- crear_mapa_malla(malla_fina,   "Malla Fina    (max.edge = 1 km)")

# Guardar cada malla individualmente
ggsave(file.path(carpeta_figuras_mallas, "malla_gruesa.png"),
       mapa_gruesa, width = 7, height = 6, dpi = 200, bg = "white")
ggsave(file.path(carpeta_figuras_mallas, "malla_media.png"),
       mapa_media,  width = 7, height = 6, dpi = 200, bg = "white")
ggsave(file.path(carpeta_figuras_mallas, "malla_fina.png"),
       mapa_fina,   width = 7, height = 6, dpi = 200, bg = "white")

# Guardar mapa comparativo (3 paneles)
comparacion <- arrangeGrob(mapa_gruesa, mapa_media, mapa_fina, ncol = 3)
ggsave(file.path(carpeta_figuras_mallas, "comparacion_mallas.png"),
       comparacion, width = 18, height = 6, dpi = 200, bg = "white")

cat("Mapas guardados en:", carpeta_figuras_mallas, "\n")
grid.arrange(mapa_gruesa, mapa_media, mapa_fina, ncol = 3)

