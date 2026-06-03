# ==============================================================================
# PASO 2: ASOCIACIÓN DE POLÍGONOS, TRÁFICO Y METEOROLOGÍA (EL CRUCE MAESTRO)
# ==============================================================================

library(data.table)
library(sf)
library(here)
source(here("R", "cleaning_functions.R"))

# 0. EL TRUCO PARA EL ERROR DE "DUPLICATE VERTEX"
# Apagamos el motor de geometría esférica (S2) temporalmente
sf_use_s2(FALSE) 

# 1. Cargar los tres datasets diarios limpios desde 'data/processed'
dt_no2     <- readRDS(here("data", "processed", "aire_madrid_2025_No2_trans_diarios.rds"))
dt_trafico <- readRDS(here("data", "processed", "trafico_madrid_2025_anual.rds"))
dt_meteo   <- readRDS(here("data", "processed", "meteo_madrid_2025_diario.rds"))


setDT(dt_no2)
setDT(dt_trafico)
setDT(dt_meteo)

# Normalizamos los nombres de distrito en tráfico
dt_trafico[, distrito := limpiar_nombres(distrito)]

# 2. Cargar el mapa oficial de distritos, REPARARLO y pasarlo a UTM (metros)
mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE)

# Reparamos geometrías corruptas o con vértices duplicados
mapa_distritos <- st_make_valid(mapa_distritos) 

# Lo proyectamos a UTM 30N (EPSG:25830) para evitar problemas esféricos
mapa_distritos <- st_transform(mapa_distritos, 25830)

# Normalizamos los nombres de distrito en el mapa
mapa_distritos$distrito <- limpiar_nombres(mapa_distritos$name)

# 3. CRUCE ESPACIAL: ¿En qué distrito se ubica físicamente cada estación de NO2?
estaciones_coords <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])

# Le decimos a R que las coordenadas originales vienen en WGS84 (EPSG:4326)
estaciones_sf <- st_as_sf(estaciones_coords, coords = c("LONGITUD", "LATITUD"), crs = 4326)

# Las transformamos a UTM 30N (EPSG:25830) para que encajen perfecto con los distritos
estaciones_sf <- st_transform(estaciones_sf, 25830)

# Spatial Join: Intersecamos las estaciones con los polígonos de los distritos
estaciones_con_distrito <- st_join(estaciones_sf, mapa_distritos[, "distrito"], join = st_intersects)
dt_estaciones_distrito  <- as.data.table(estaciones_con_distrito)
dt_estaciones_distrito[, geometry := NULL]

# Sincronizamos el nombre del distrito de vuelta en el dataset de contaminación
dt_no2 <- merge(dt_no2, dt_estaciones_distrito[, .(ESTACION, distrito)], by = "ESTACION", all.x = TRUE)

# 4. INTEGRACIÓN DE LA METEOROLOGÍA 
# Calculamos la media diaria de la ciudad para las variables climáticas
cols_clima <- setdiff(names(dt_meteo), c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "FECHA", "ID_TIEMPO")) 
#Set differ se utiliza para 

dt_meteo_ciudad <- dt_meteo[, lapply(.SD, mean, na.rm = TRUE), by = .(FECHA), .SDcols = cols_clima]

# 5. EL CRUCE MAESTRO (Unificación Espacio-Temporal)
# Primero: Cruzamos el Tráfico usando como llaves la FECHA y el DISTRITO
dt_maestro <- merge(dt_no2, dt_trafico, by = c("FECHA", "distrito"), all.x = TRUE)

# Segundo: Cruzamos la Meteorología usando únicamente la FECHA como llave
dt_maestro <- merge(dt_maestro, dt_meteo_ciudad, by = "FECHA", all.x = TRUE)

# Creamos un ID numérico para los distritos
dt_maestro[, ID_DISTRITO := .GRP, by = distrito]

# 5b. ESTANDARIZACIÓN DE COVARIABLES DE TRÁFICO
# Escalamos a media 0 y desviación típica 1 para estabilidad numérica en INLA
# y facilitar la interpretación de los coeficientes (en unidades de sd).
dt_maestro[, intensidad_raw := intensidad]  # Guardamos los valores originales
dt_maestro[, carga_raw      := carga]

dt_maestro[, intensidad := scale(intensidad)[, 1]]
dt_maestro[, carga      := scale(carga)[, 1]]

cat("--- Estandarización de tráfico ---\n")
cat("intensidad: media =", round(mean(dt_maestro$intensidad, na.rm = TRUE), 4),
    "| sd =", round(sd(dt_maestro$intensidad, na.rm = TRUE), 4), "\n")
cat("carga:      media =", round(mean(dt_maestro$carga, na.rm = TRUE), 4),
    "| sd =", round(sd(dt_maestro$carga, na.rm = TRUE), 4), "\n")

# Ordenación cronológica y espacial estricta (Indispensable para el AR1 de INLA) 
setorder(dt_maestro, ID_TIEMPO, ESTACION)
View(dt_maestro)
# 6. Guardar el archivo final unificado
saveRDS(dt_maestro, here("data", "processed", "dataset_maestro_inla_2025.rds"))

# Restauramos el motor S2 para no afectar a otros scripts que puedas tener
sf_use_s2(TRUE) 

# ------------------------------------------------------------------------------
# CONTROL DE CALIDAD
# ------------------------------------------------------------------------------
cat("✅ ¡Paso 2 completado con éxito!\n")
cat("Filas totales del Dataset Maestro (Días × Estaciones):", nrow(dt_maestro), "\n")
cat("¿Existen NAs en el ID_TIEMPO?:", any(is.na(dt_maestro$ID_TIEMPO)), "\n")
print(head(dt_maestro[, .(FECHA, ESTACION, distrito, LOG_NO2_DIARIO, ID_TIEMPO, ID_DISTRITO)]))

