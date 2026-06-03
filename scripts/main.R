# ==============================================================================
# PASO 1: CREACIÓN DE LA MALLA ESPACIAL (MESH)
# ==============================================================================

library(INLA)
library(sf)
library(data.table)
library(here)

# 1. Cargar datos de NO2 (variable respuesta)
dt_no2 <- readRDS(here("data", "processed", "aire_madrid_2025_No2_trans_diarios.rds"))
setDT(dt_no2)

View(dt_no2)

# 2. Extraer y proyectar coordenadas espaciales
# Se extraen las coordenadas únicas y se proyectan a UTM 30N (EPSG:25830) en kilómetros.
coords_estaciones <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])
coords_sf <- st_as_sf(coords_estaciones, coords = c("LONGITUD", "LATITUD"), crs = 4326) # CRS original en WGS84 (grados)
coords_utm <- st_transform(coords_sf, 25830) # Proyección a UTM zona 30N (metros)

# Matriz de coordenadas en kilómetros (requerido por INLA para estabilidad numérica)
coords_matriz <- st_coordinates(coords_utm) / 1000 

# 3. Definir los límites geométricos (Boundaries)
# Se genera un contorno interno ajustado a las estaciones y un contorno externo de amortiguación.
# La función convex define cómo se dibuja la línea de la frontera alrededor de las estaciones.
bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.05, resolution = 50) # Se deja un 5% de margen adicional para incluir todas las estaciones dentro del contorno. Esto ayuda a asegurar que el modelo espacial tenga en cuenta todas las ubicaciones de las estaciones sin ser demasiado ajustado.
bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.2) # Se deja un 20% de margen adicional para amortiguación fuera de las estaciones. Esto ayuda a evitar problemas de borde en el modelado espacial.

# 4. Construir la Malla (Mesh). También llamado discretizar ( es el proceso de dividir un espacio continuo espacial en regiones discretas)
malla_madrid <- inla.mesh.2d(
  loc = coords_matriz, # Ubicaciones de las estaciones (en km) #ES la matriz de puntos de las localizaciones usado para determinar el dominio. 
  boundary = list(bnd_inner, bnd_outer), # Contornos interno y externo
  max.edge = c(2, 5),  # 2 km de resolución interna, 5 km externa. El triángulo más pequeño tendrá lados de 2 km, y el más grande de 5 km. Esto controla la densidad de la malla.
  cutoff = 0.5         # Nos indica que si hay dos estaciones a menos de 500 metros de distancia, se trata como si fuera el mismo punto. 
)

# 5. Exportar y visualizar
saveRDS(malla_madrid, here("data", "processed", "malla_spde_madrid.rds"))

pdf(here("output", "figures", "malla_espacial_madrid.pdf"))
plot(malla_madrid, main = "Malla SPDE - Madrid 2025 (Km)")
points(coords_matriz, col = "red", pch = 16, cex = 1.2)
dev.off()


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
library(ggplot2)
library(sf)
library(fmesher)
library(here)

# 1. Cargar el mapa de distritos original
mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"))

# 2. Asegurarnos de que está en el CRS oficial (25830 - metros) y pasarlo a KILÓMETROS
mapa_distritos_utm <- st_transform(mapa_distritos, 25830)
st_geometry(mapa_distritos_utm) <- st_geometry(mapa_distritos_utm) / 1000

# 3. Preparar los puntos de las estaciones en formato data.frame (en km)
coords_df <- as.data.frame(coords_matriz)
colnames(coords_df) <- c("X", "Y")

# 4. El súper mapa combinando todo en Kilómetros
ggplot() +
  # Capa 1: Los distritos de Madrid de fondo (sin relleno, solo el borde)
  geom_sf(data = mapa_distritos_utm, fill = NA, color = "gray40", linewidth = 0.5) +
  
  # Capa 2: La malla de triángulos de INLA
  geom_fm(data = malla_madrid, color = "blue", linewidth = 0.2, alpha = 0.4) +
  
  # Capa 3: Tus estaciones de contaminación en rojo
  geom_point(data = coords_df, aes(x = X, y = Y), color = "red", size = 2) +
  
  labs(
    title = "Malla SPDE e Infraestructura de Distritos en Madrid",
    subtitle = "Escala unificada en Kilómetros (UTM Zona 30N)",
    x = "Eje X (km)", y = "Eje Y (km)"
  ) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "lightgrey"))

