# ==============================================================================
# PASO 2: ASOCIACIÓN DE POLÍGONOS, TRÁFICO Y METEOROLOGÍA (EL CRUCE MAESTRO)
# ==============================================================================

library(data.table)
library(sf)
library(gstat)
library(here)
source(here("R", "cleaning_functions.R"))
source(here("R", "FUNCIONES_INTERPOLACION.R"))

# 0. EL TRUCO PARA EL ERROR DE "DUPLICATE VERTEX"
# Apagamos el motor de geometría esférica (S2) temporalmente
sf_use_s2(FALSE) 

# 1. Cargar los tres datasets diarios limpios desde 'data/processed'
dt_no2     <- readRDS(here("data", "processed","Contaminacion","diario","aire_madrid_2025_No2_trans_diarios.rds"))
dt_trafico <- readRDS(here("data", "processed","trafico_madrid_2025_diario_barrio.rds"))
dt_meteo   <- readRDS(here("data", "processed","Clima","diario", "meteo_madrid_2025_diario.rds"))
setDT(dt_no2)
setDT(dt_trafico)
setDT(dt_meteo)

# Filtramos solo el año 2025
dt_no2     <- dt_no2[year(FECHA) == 2025]
dt_trafico <- dt_trafico[year(FECHA) == 2025]
dt_meteo   <- dt_meteo[year(FECHA) == 2025]

# Normalizamos los nombres de distrito en tráfico
dt_trafico[, distrito := limpiar_nombres(distrito)]

# 2. Cargar geometrías de distritos y barrios, reparar y proyectar a UTM
mapa_distritos <- st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE) |>
  st_make_valid() |> st_transform(25830)
mapa_distritos$distrito <- limpiar_nombres(mapa_distritos$name)

mapa_barrios <- st_read(here("data", "raw", "Geometrias", "BARRIOS.shp"), quiet = TRUE) |>
  st_make_valid() |> st_transform(25830)
mapa_barrios$barrio <- tolower(trimws(mapa_barrios$NOMBRE))

# Normalizamos el nombre de barrio en tráfico al mismo formato
dt_trafico[, barrio := tolower(trimws(barrio))]

# 3. CRUCE ESPACIAL: ¿En qué distrito y barrio se ubica cada estación de NO2?
estaciones_coords <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])
estaciones_sf <- st_as_sf(estaciones_coords, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)

# Asignamos distrito
est_distrito <- st_join(estaciones_sf, mapa_distritos[, "distrito"], join = st_intersects)
# Asignamos barrio
est_barrio <- st_join(estaciones_sf, mapa_barrios[, "barrio"], join = st_intersects)

dt_est_geo <- merge(
  as.data.table(est_distrito)[, .(ESTACION, distrito)],
  as.data.table(est_barrio)[, .(ESTACION, barrio)],
  by = "ESTACION"
)

# Sincronizamos distrito y barrio en el dataset de contaminación
dt_no2 <- merge(dt_no2, dt_est_geo, by = "ESTACION", all.x = TRUE)

# 4. INTEGRACIÓN DE LA METEOROLOGÍA (IDW diario por estación)
# Interpolamos las variables climáticas desde las estaciones meteo a las
# ubicaciones de las estaciones de NO2, día a día, usando IDW.
coords_no2 <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])

dt_clima_interp <- interpolar_idw_clima(
  dt_meteo    = dt_meteo,
  dt_objetivo = coords_no2,
  variables   = c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                  "Presion Barométrica", "Radiación Solar")
)

# Guardamos una copia para no recalcular
saveRDS(dt_clima_interp,
        here("data", "processed", "Clima", "diario", "clima_interpolado_diario_2025.rds"))

# 5. EL CRUCE MAESTRO (Unificación Espacio-Temporal)
# Primero: Cruzamos el Tráfico usando como llaves la FECHA y el BARRIO
dt_maestro <- merge(dt_no2, dt_trafico[, .(barrio, FECHA, intensidad, ocupacion, carga)],
                    by = c("FECHA", "barrio"), all.x = TRUE)

# Segundo: Cruzamos la Meteorología interpolada por ESTACION y FECHA
dt_maestro <- merge(dt_maestro, dt_clima_interp, by = c("ESTACION", "FECHA"), all.x = TRUE)

# Creamos un ID numérico para los distritos
dt_maestro[, ID_DISTRITO := .GRP, by = distrito]

# 5b. ESTANDARIZACIÓN DE COVARIABLES (TRÁFICO + CLIMA)
# Escalamos a media 0 y desviación típica 1 para estabilidad numérica en INLA
# y facilitar la interpretación de los coeficientes (en unidades de sd).

# Tráfico
dt_maestro[, intensidad_raw := intensidad]
dt_maestro[, carga_raw      := carga]
dt_maestro[, intensidad := scale(intensidad)[, 1]]
dt_maestro[, carga      := scale(carga)[, 1]]

# Clima
cols_clima_std <- c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                    "Presion Barométrica", "Radiación Solar")
cols_clima_std <- intersect(cols_clima_std, names(dt_maestro))

for (v in cols_clima_std) {
  raw_name <- paste0(v, "_raw")
  dt_maestro[, (raw_name) := get(v)]
  dt_maestro[, (v) := scale(get(v))[, 1]]
}

cat("--- Estandarización de covariables ---\n")
cols_std <- c("intensidad", "carga", cols_clima_std)
for (v in cols_std) {
  cat(sprintf("%-25s media = %7.4f | sd = %6.4f\n",
              v, mean(dt_maestro[[v]], na.rm = TRUE), sd(dt_maestro[[v]], na.rm = TRUE)))
}

# Ordenación cronológica y espacial estricta (Indispensable para el AR1 de INLA) 
setorder(dt_maestro, ID_TIEMPO, ESTACION)
# 6. Guardar el archivo final unificado
saveRDS(dt_maestro, here("data", "processed", "dataset_maestro_inla_2025_DIARIO.rds"))

# Restauramos el motor S2 para no afectar a otros scripts que puedas tener
sf_use_s2(TRUE) 

revision_dataset<-View(dt_maestro)

# ------------------------------------------------------------------------------
# CONTROL DE CALIDAD
# ------------------------------------------------------------------------------
cat("✅ ¡Paso 2 completado con éxito!\n")
cat("Filas totales del Dataset Maestro (Días × Estaciones):", nrow(dt_maestro), "\n")
cat("¿Existen NAs en el ID_TIEMPO?:", any(is.na(dt_maestro$ID_TIEMPO)), "\n")
print(head(dt_maestro[, .(FECHA, ESTACION, distrito, LOG_NO2_DIARIO, ID_TIEMPO, ID_DISTRITO)]))


