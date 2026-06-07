#-------------------------------------------------------------------------------
#----------------------Spatial Correlation--------------------------------------
#-------------------------------------------------------------------------------


#----------------------Air quality monitoring stations--------------------------
# ==============================================================================
# ANÁLISIS DE CORRELACIÓN ESPACIAL PARA DATOS PUNTUALES (VARIOGRAMA)
# ==============================================================================
library(data.table)
library(sf)
library(gstat)
library(ggplot2)
library(here)

# 1. CARGAR DATOS
# Cargamos tu dataset diario transformado (usando la ruta correcta a tu archivo RDS)
ruta_datos <- here("data", "processed", "contaminacion", "diario", "aire_madrid_2025_No2_trans_diarios.rds")
dt_aire <- readRDS(ruta_datos)

# 2. AGREGACIÓN ANUAL
# Calculamos el nivel medio anual de NO2 (en su escala original) para cada estación
dt_estaciones <- dt_aire[!is.na(DATO_DIARIO), .(
  media_no2 = mean(DATO_DIARIO, na.rm = TRUE),
  X = unique(LONGITUD),
  Y = unique(LATITUD)
), by = ESTACION]

# 3. CONVERSIÓN A OBJETO ESPACIAL
# R necesita saber que la X y la Y son coordenadas reales de la Tierra.
# Asumimos que tus coordenadas (LONGITUD/LATITUD) están en WGS84 (EPSG:4326)
sf_estaciones <- st_as_sf(dt_estaciones, 
                          coords = c("X", "Y"), 
                          crs = 4326)

# Las proyectamos a UTM 30N (metros) para que el variograma calcule la distancia real
sf_estaciones <- st_transform(sf_estaciones, crs = 25830)

# 4. CÁLCULO DEL VARIOGRAMA EXPERIMENTAL
# La fórmula (media_no2 ~ 1) significa que evaluamos la media espacial constante.
# 'cutoff' es la distancia máxima a evaluar (e.g., 20.000 metros = 20 km)
variograma_exp <- variogram(media_no2 ~ 1, 
                            data = sf_estaciones, 
                            cutoff = 20000, 
                            width = 2000) # Agrupamos las distancias de 2 en 2 km

# 5. VISUALIZACIÓN (ggplot2)
# Nota: NO se usa geom_smooth sobre el variograma — ya es una estimación suavizada
# por bins de distancia. Añadir loess encima puede crear patrones artificiales.
plot_variograma <- ggplot(variograma_exp, aes(x = dist / 1000, y = gamma)) +
  geom_line(color = "#34495e", linetype = "dashed") +
  geom_point(size = 3, color = "#2c3e50") +
  labs(
    title = "Variograma Experimental: Correlación Espacial del NO2",
    subtitle = "Ciudad de Madrid (Media Anual)",
    x = "Distancia entre estaciones (Kilómetros)",
    y = "Semivarianza",
    caption = "Nota: Una curva ascendente indica presencia de autocorrelación espacial positiva."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11),
    panel.grid.minor = element_blank()
  )

print(plot_variograma)


# ==============================================================================
# ANÁLISIS DE CORRELACIÓN ESPACIAL PARA VARIABLES CLIMÁTICAS (VARIOGRAMA)
# ==============================================================================
# ANÁLISIS DE CORRELACIÓN ESPACIAL PARA VARIABLES CLIMÁTICAS (VARIOGRAMA)
# ==============================================================================
library(data.table)
library(sf)
library(gstat)
library(ggplot2)
library(here)

# 1. CARGAR DATOS METEOROLÓGICOS DIARIOS
ruta_meteo <- here("data", "processed","Clima", "mensual", "meteo_madrid_2025_mensual.rds")
dt_meteo <- readRDS(ruta_meteo)

# 2. AGREGACIÓN ANUAL POR ESTACIÓN
# Calculamos la media anual usando los nombres exactos de tu diccionario
dt_estaciones_meteo <- dt_meteo[, .(
  Temperatura      = mean(Temperatura, na.rm = TRUE),
  Humedad_Relativa = mean(Humedad_Relativa, na.rm = TRUE),
  LONGITUD         = unique(LONGITUD),
  LATITUD          = unique(LATITUD)
), by = ESTACION]

# 3. CONVERSIÓN A OBJETO ESPACIAL PROYECTADO (UTM 30N - metros)
sf_meteo <- st_as_sf(dt_estaciones_meteo, 
                     coords = c("LONGITUD", "LATITUD"), 
                     crs = 4326)
sf_meteo <- st_transform(sf_meteo, crs = 25830)
# ==============================================================================
# 4. CÁLCULO DE VARIOGRAMAS EXPERIMENTALES (FILTRADO POR VARIABLE)
# ==============================================================================

# --- Temperatura ---
# Seleccionamos solo las estaciones que SÍ tienen datos de Temperatura
sf_temp <- sf_meteo[!is.na(sf_meteo$Temperatura), ]

cat("Estaciones válidas para Temperatura:", nrow(sf_temp), "\n")
v_temp <- variogram(Temperatura ~ 1, data = sf_temp, cutoff = 20000, width = 2000)
v_temp$Variable <- "Temperatura (°C)"


# --- Humedad Relativa ---
# Seleccionamos solo las estaciones que SÍ tienen datos de Humedad Relativa
sf_hum <- sf_meteo[!is.na(sf_meteo$Humedad_Relativa), ]

cat("Estaciones válidas para Humedad Relativa:", nrow(sf_hum), "\n")
v_hum <- variogram(Humedad_Relativa ~ 1, data = sf_hum, cutoff = 20000, width = 2000)
v_hum$Variable <- "Humedad Relativa (%)"


# Combinamos ambos resultados en una única tabla para el gráfico
v_clima_total <- rbind(v_temp, v_hum)

# ==============================================================================
# 5. VISUALIZACIÓN ACADÉMICA COMPARATIVA (ggplot2)
# ==============================================================================
plot_clima <- ggplot(v_clima_total, aes(x = dist / 1000, y = gamma)) +
  geom_line(color = "#34495e", linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.5, color = "#2c3e50") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  labs(
    title = "Variograma Experimental de Covariables Meteorológicas",
    subtitle = "Ciudad de Madrid (Medias Anuales)",
    x = "Distancia entre estaciones meteorológicas (Kilómetros)",
    y = "Semivarianza",
    caption = "Nota: El comportamiento ascendente demuestra la continuidad espacial de la variable."
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title        = element_text(face = "bold", size = 14, margin = margin(b = 5)),
    plot.subtitle     = element_text(color = "gray40", size = 11, margin = margin(b = 12)),
    strip.text        = element_text(face = "bold", size = 11, color = "#2c3e50"),
    strip.background  = element_rect(fill = "gray95", color = "gray85"),
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(1.5, "lines")
  )

print(plot_clima)


# ==============================================================================
# ANÁLISIS DE CORRELACIÓN ESPACIAL DIARIA CON ENVOLVENTES DE MONTE CARLO
# ==============================================================================
library(data.table)
library(sf)
library(gstat)
library(ggplot2)
library(here)

# Parámetro global para Monte Carlo
num_simulaciones <- 99

cat("\n==================================================\n")
cat("1. VARIOGRAMA Y MONTE CARLO PARA NO2 (DÍA ESPECÍFICO)\n")
cat("==================================================\n")

# 1.1 CARGA Y PREPARACIÓN NO2
ruta_datos <- here("data", "processed", "contaminacion", "diario", "aire_madrid_2025_No2_trans_diarios.rds")
dt_aire <- readRDS(ruta_datos)

# 🔥 CAMBIO CLAVE: Filtramos un día concreto de alta contaminación (Invierno)
dia_no2 <- as.Date("2025-01-15")

dt_estaciones <- dt_aire[as.Date(FECHA) == dia_no2 & !is.na(DATO_DIARIO), .(
  valor_no2 = DATO_DIARIO, # Tomamos el valor exacto de ese día, no la media
  X = unique(LONGITUD),
  Y = unique(LATITUD)
), by = ESTACION]

sf_estaciones <- st_as_sf(dt_estaciones, coords = c("X", "Y"), crs = 4326)
sf_estaciones <- st_transform(sf_estaciones, crs = 25830)

# 1.2 VARIOGRAMA EMPÍRICO NO2
variograma_exp <- variogram(valor_no2 ~ 1, data = sf_estaciones, cutoff = 20000, width = 2000)

# 1.3 ENVOLVENTE MONTE CARLO NO2
cat("   -> Ejecutando", num_simulaciones, "permutaciones espaciales para NO2...\n")
set.seed(8271)
sf_simulacion  <- sf_estaciones
sim_resultados <- matrix(NA, nrow = nrow(variograma_exp), ncol = num_simulaciones)

for (i in seq_len(num_simulaciones)) {
  sf_simulacion$valor_no2 <- sample(sf_estaciones$valor_no2)
  v_sim <- variogram(valor_no2 ~ 1, data = sf_simulacion, cutoff = 20000, width = 2000)
  if (nrow(v_sim) == nrow(variograma_exp))
    sim_resultados[, i] <- v_sim$gamma
}

variograma_exp$limite_inferior <- apply(sim_resultados, 1, quantile, probs = 0.025, na.rm = TRUE)
variograma_exp$limite_superior <- apply(sim_resultados, 1, quantile, probs = 0.975, na.rm = TRUE)

# 1.4 VISUALIZACIÓN NO2
plot_variograma_no2 <- ggplot(variograma_exp, aes(x = dist / 1000)) +
  geom_ribbon(aes(ymin = limite_inferior, ymax = limite_superior), fill = "gray80", alpha = 0.5) +
  geom_line(aes(y = gamma), color = "#34495e", linetype = "dashed") +
  geom_point(aes(y = gamma), size = 3, color = "#2c3e50") +
  labs(
    title = "Variograma Experimental con Envolvente de Monte Carlo (95%)",
    subtitle = paste("Autocorrelación Espacial del NO2 -", format(dia_no2, "%d de %B de %Y")),
    x = "Distancia entre estaciones (Kilómetros)",
    y = "Semivarianza",
    caption = "Nota: Curva empírica fuera de la franja gris confirma correlación espacial no aleatoria."
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40", size = 11),
        panel.grid.minor = element_blank())

print(plot_variograma_no2)


cat("\n==================================================\n")
cat("2. VARIOGRAMAS Y MONTE CARLO PARA CLIMATOLOGÍA (DÍA ESPECÍFICO)\n")
cat("==================================================\n")

# 2.1 CARGA Y PREPARACIÓN CLIMA
ruta_meteo <- here("data", "processed","Clima", "diario", "meteo_madrid_2025_diario.rds")
dt_meteo <- readRDS(ruta_meteo)

# 🔥 CAMBIO CLAVE: Filtramos un día concreto de verano (Isla de Calor)
dia_clima <- as.Date("2025-07-15")

dt_estaciones_meteo <- dt_meteo[as.Date(FECHA) == dia_clima, .(
  Temperatura      = Temperatura,
  Humedad_Relativa = Humedad_Relativa,
  LONGITUD         = unique(LONGITUD),
  LATITUD          = unique(LATITUD)
), by = ESTACION]

sf_meteo <- st_as_sf(dt_estaciones_meteo, coords = c("LONGITUD", "LATITUD"), crs = 4326)
sf_meteo <- st_transform(sf_meteo, crs = 25830)

# 2.2 TEMPERATURA: Variograma y Monte Carlo
cat("   -> Ejecutando", num_simulaciones, "permutaciones espaciales para Temperatura...\n")
sf_temp <- sf_meteo[!is.na(sf_meteo$Temperatura), ]
v_temp  <- variogram(Temperatura ~ 1, data = sf_temp, cutoff = 20000, width = 2000)

set.seed(8271)
sf_sim_temp  <- sf_temp
sim_temp     <- matrix(NA, nrow = nrow(v_temp), ncol = num_simulaciones)
for (i in seq_len(num_simulaciones)) {
  sf_sim_temp$Temperatura <- sample(sf_temp$Temperatura)
  v_sim <- variogram(Temperatura ~ 1, data = sf_sim_temp, cutoff = 20000, width = 2000)
  if (nrow(v_sim) == nrow(v_temp)) sim_temp[, i] <- v_sim$gamma
}
v_temp$Variable        <- "Temperatura (°C)"
v_temp$limite_inferior <- apply(sim_temp, 1, quantile, probs = 0.025, na.rm = TRUE)
v_temp$limite_superior <- apply(sim_temp, 1, quantile, probs = 0.975, na.rm = TRUE)

# 2.3 HUMEDAD RELATIVA: Variograma y Monte Carlo
cat("   -> Ejecutando", num_simulaciones, "permutaciones espaciales para Humedad...\n")
sf_hum <- sf_meteo[!is.na(sf_meteo$Humedad_Relativa), ]
v_hum  <- variogram(Humedad_Relativa ~ 1, data = sf_hum, cutoff = 20000, width = 2000)

set.seed(8271)
sf_sim_hum  <- sf_hum
sim_hum     <- matrix(NA, nrow = nrow(v_hum), ncol = num_simulaciones)
for (i in seq_len(num_simulaciones)) {
  sf_sim_hum$Humedad_Relativa <- sample(sf_hum$Humedad_Relativa)
  v_sim <- variogram(Humedad_Relativa ~ 1, data = sf_sim_hum, cutoff = 20000, width = 2000)
  if (nrow(v_sim) == nrow(v_hum)) sim_hum[, i] <- v_sim$gamma
}
v_hum$Variable        <- "Humedad Relativa (%)"
v_hum$limite_inferior <- apply(sim_hum, 1, quantile, probs = 0.025, na.rm = TRUE)
v_hum$limite_superior <- apply(sim_hum, 1, quantile, probs = 0.975, na.rm = TRUE)

# 2.4 VISUALIZACIÓN COMBINADA CLIMA
v_clima_total <- rbind(v_temp, v_hum)

plot_clima_mc <- ggplot(v_clima_total, aes(x = dist / 1000)) +
  geom_ribbon(aes(ymin = limite_inferior, ymax = limite_superior), fill = "gray80", alpha = 0.5) +
  geom_line(aes(y = gamma), color = "#34495e", linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(y = gamma), size = 2.5, color = "#2c3e50") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  labs(
    title = "Variograma Experimental de Covariables (con Monte Carlo 95%)",
    subtitle = paste("Ciudad de Madrid -", format(dia_clima, "%d de %B de %Y")),
    x = "Distancia entre estaciones meteorológicas (Kilómetros)",
    y = "Semivarianza",
    caption = "Nota: La franja gris representa la hipótesis nula de aleatoriedad espacial."
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title        = element_text(face = "bold", size = 14, margin = margin(b = 5)),
    plot.subtitle     = element_text(color = "gray40", size = 11, margin = margin(b = 12)),
    strip.text        = element_text(face = "bold", size = 11, color = "#2c3e50"),
    strip.background  = element_rect(fill = "gray95", color = "gray85"),
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(1.5, "lines")
  )

print(plot_clima_mc)
