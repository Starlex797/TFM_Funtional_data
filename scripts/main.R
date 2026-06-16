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

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# VISUALIZACIÓN COMPARATIVA DE MALLAS
#-----------------------------------------------------------------------------------
library(ggplot2)
library(sf)
library(fmesher)
library(here)
library(gridExtra) # Para poner los gráficos uno al lado del otro

# 1. Cargar el mapa de distritos original
mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"))

# 2. Asegurarnos de que está en el CRS oficial (25830 - metros) y pasarlo a KILÓMETROS
mapa_distritos_utm <- st_transform(mapa_distritos, 25830)
st_geometry(mapa_distritos_utm) <- st_geometry(mapa_distritos_utm) / 1000

# 3. Preparar los puntos de las estaciones en formato data.frame (en km)
coords_df <- as.data.frame(coords_matriz)
colnames(coords_df) <- c("X", "Y")

# 4. Crear una función base de ggplot para no repetir código
crear_mapa_malla <- function(malla_input, titulo_malla) {
  ggplot() +
    # Capa 1: Los distritos de Madrid de fondo (sin relleno, solo el borde)
    geom_sf(data = mapa_distritos_utm, fill = NA, color = "gray40", linewidth = 0.5) +
    
    # Capa 2: La malla de triángulos de INLA que le pasemos
    geom_fm(data = malla_input, color = "blue", linewidth = 0.2, alpha = 0.4) +
    
    # Capa 3: Tus estaciones de contaminación en rojo
    geom_point(data = coords_df, aes(x = X, y = Y), color = "red", size = 2) +
    
    labs(
      title = titulo_malla,
      subtitle = "UTM Zona 30N (km)",
      x = "", y = "" # Quitamos texto de ejes para que quede más limpio al juntarlos
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = "lightgrey"),
      plot.title = element_text(face = "bold", size = 12)
    )
}

# 5. Generar los tres gráficos usando la función
grafico_gruesa <- crear_mapa_malla(malla_gruesa, "Malla Gruesa (4 km)")
grafico_media  <- crear_mapa_malla(malla_media, "Malla Media (2 km) - ACTUAL")
grafico_fina   <- crear_mapa_malla(malla_fina, "Malla Fina (1 km)")

# 6. Guardar el súper mapa comparativo combinando los tres
pdf(here("output", "figures", "comparacion_mallas_madrid.pdf"), width = 15, height = 5)
grid.arrange(grafico_gruesa, grafico_media, grafico_fina, ncol = 3)
dev.off()

# También lo puedes imprimir en la pestaña de Plots de RStudio para verlo ahora mismo:
grid.arrange(grafico_gruesa, grafico_media, grafico_fina, ncol = 3)

