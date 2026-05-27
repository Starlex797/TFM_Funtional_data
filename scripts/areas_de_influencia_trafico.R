library(sf)
library(data.table)

#--------------------------Lectura de los datos -----------------------------------------
dt_ubicacion_trafico <- fread("C:\\Users\\HP\\Desktop\\TFM\\Base_de_datos\\Datos_trafico\\Detectores\\2025_Enero_medidores.csv")
df_estaciones_contami <- readRDS("estaciones_contami.rds")
names(df_estaciones_contami)
#----------------------------------------------------------------------------------------
# 2. Convertir a objeto espacial (sf)

geo_trafico <- st_as_sf(dt_ubicacion_trafico, coords = c("longitud", "latitud"), crs = 4326)

# 3. Cargar estaciones de contaminación 
geo_contaminacion <- st_as_sf(df_estaciones_contami, coords = c("longitud", "latitud"), crs = 4326)

# 4. Proyectar a un sistema métrico para Madrid 

geo_trafico <- st_transform(geo_trafico, 25830)
geo_contami <- st_transform(geo_contaminacion, 25830)
#-----------------------------------------------------------
# Creamos el área de influencia (500 metros)
buffer_contami <- st_buffer(geo_contami, dist = 500)

# Esto unirá a cada medidor de tráfico el nombre de la estación de contaminación más cercana
puente_trafico_contami <- st_join(geo_trafico, buffer_contami, join = st_intersects)

# Filtramos para quedarnos solo con los que cayeron dentro de algún área
puente_final <- puente_trafico_contami[!is.na(puente_trafico_contami$estacion), ]

# Convertimos a data.table para trabajar más rápido después
setDT(puente_final)
#---------------------Número de medidores bajo el área de influencia---------
conteo_medidores <- puente_final[, .(
  num_medidores = uniqueN(id),               # Cuenta IDs únicos de medidores 
  medidores_m30 = sum(tipo_elem == "M30"),     # Cuántos son de la M30 
  medidores_urb = sum(tipo_elem == "URB")   # Cuántos son urbanos 
), by = .(estacion, nom_tipo)]          # Agrupado por estación y su tipo

# Ordenar de mayor a menor para ver cuáles tienen más cobertura
conteo_medidores <- conteo_medidores[order(-num_medidores)]

# Visualizar el resultado
print(conteo_medidores)
#------------------------------Mapa------------------------------------------
library(ggplot2)
ggplot() +
  # 1. Medidores de tráfico totales (fondo)
  # Diferenciamos M30 de URB. Los puntos M30 dibujarán el trazado de la autopista.
  geom_sf(data = geo_trafico, aes(color = tipo_elem), size = 0.5, alpha = 0.4) +
  
  # 2. Áreas de influencia de 500m
  geom_sf(data = buffer_contami, fill = "blue", alpha = 0.1, color = "blue", linetype = "dashed") +
  
  # 3. Medidores SELECCIONADOS (dentro de los buffers)
  # Conservan su identificación única e invariable.
  geom_sf(data = st_as_sf(puente_final), color = "purple", size = 1, alpha = 0.8) +
  
  # 4. Estaciones de contaminación diferenciadas por tipo (nom_tipo)
  # Usamos shape 23 (diamante con relleno) para separar la leyenda de colores.
  geom_sf(data = geo_contami, aes(fill = nom_tipo), shape = 23, size = 4, color = "black") +
  
  # Configuración de escalas: M30 en naranja para que resalte
  scale_color_manual(values = c("M30" = "orange", "URB" = "darkred"), name = "Red de Tráfico") +
  scale_fill_brewer(palette = "Set1", name = "Tipo Estación Contam.") +
  
  labs(
    title = "Análisis Espacial: Tráfico (M30 vs Urbano) y Contaminación",
    subtitle = "La M30 (naranja) es el único entorno con datos de velocidad media (vmed) ",
    caption = "Fuente: Ayuntamiento de Madrid. La Carga en M30 suele ser cero.",
    x = "Longitud (UTM)",
    y = "Latitud (UTM)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")

#-----------------------Mapa Interactivo-----------------------------------
library(leaflet)
library(sf)

library(leaflet)
library(sf)

# 1. Transformamos todo a coordenadas geográficas (WGS84 - EPSG:4326)
geo_trafico_4326 <- st_transform(geo_trafico, 4326)
geo_contami_4326 <- st_transform(geo_contami, 4326)
puente_final_4326 <- st_transform(st_as_sf(puente_final), 4326)
buffer_contami_4326 <- st_transform(buffer_contami, 4326)

# 2. Definimos las paletas de colores
# Paleta para los medidores de tráfico
pal_trafico <- colorFactor(palette = c("orange", "darkred"), domain = geo_trafico_4326$tipo_elem)

# Paleta NUEVA para las estaciones de contaminación según su tipo
# Usamos colores que contrasten bien: Azul, Verde y Negro
pal_contami <- colorFactor(palette = c("blue", "green", "black"), domain = geo_contami_4326$nom_tipo)

# 3. Creamos el mapa
leaflet() %>%
  addTiles() %>%  # Fondo de mapa estándar
  
  # Capa 1: Áreas de influencia (Buffers)
  addPolygons(data = buffer_contami_4326, color = "blue", weight = 1, 
              fillOpacity = 0.1, group = "Áreas de Influencia (500m)") %>%
  
  # Capa 2: Todos los medidores de tráfico de Madrid
  addCircleMarkers(data = geo_trafico_4326, radius = 2, 
                   color = ~pal_trafico(tipo_elem), stroke = FALSE, fillOpacity = 0.4,
                   popup = ~paste("ID:", id, "<br>Tipo:", tipo_elem),
                   group = "Todos los Medidores") %>%
  
  # Capa 3: Medidores SELECCIONADOS (Resaltados en morado)
  addCircleMarkers(data = puente_final_4326, radius = 4, 
                   popup = ~paste("<b>Medidor Seleccionado</b><br>ID:", id, 
                                  "<br>Tipo:", tipo_elem,
                                  "<br>Estación asociada:", estacion),
                   group = "Medidores Seleccionados") %>%
  
  # Capa 4: Estaciones de Contaminación (AHORA CON COLORES POR TIPO)
  addCircleMarkers(data = geo_contami_4326, 
                   radius = 8, # Las hacemos más grandes para que destaquen
                   color = ~pal_contami(nom_tipo), 
                   weight = 3, fillOpacity = 1, opacity = 1,
                   popup = ~paste("<b>Estación de Contaminación</b><br>", 
                                  estacion, "<br>Tipo:", nom_tipo),
                   group = "Estaciones Contaminación") %>%
  
  # Control de capas
  addLayersControl(
    overlayGroups = c("Estaciones Contaminación", "Medidores Seleccionados", 
                      "Áreas de Influencia (500m)", "Todos los Medidores"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  
  # Leyenda 1: Red de Tráfico (Abajo a la derecha)
  addLegend(pal = pal_trafico, values = geo_trafico_4326$tipo_elem, 
            title = "Red de Tráfico", position = "bottomright") %>%
  
  # Leyenda 2: Tipo de Estación de Contaminación (Arriba a la derecha)
  addLegend(pal = pal_contami, values = geo_contami_4326$nom_tipo, 
            title = "Tipo Estación Contam.", position = "topright")

