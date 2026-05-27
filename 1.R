library(dplyr)
library(lubridate)

# 1. URL para el año 2016 completo (Barajas, Cuatro Vientos, Getafe, Torrejón)
# Pedimos todas las estaciones y todas las variables posibles
url_todo_2016 <- "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=LEMD&station=LECU&station=LEGT&station=LETO&data=all&year1=2016&month1=1&day1=1&year2=2016&month2=12&day2=31&tz=Etc%2FUTC&format=comma&latlon=no&direct=yes"

# 2. Descargar y leer (na.strings = "M" porque IEM usa M para datos faltantes)
df_2016_completo <- read.csv(url_todo_2016, skip = 5, na.strings = "M")

# 3. Limpiar y convertir unidades (IEM viene en sistema americano)
df_final_2016 <- df_2016_completo %>%
  mutate(
    fecha_hora = as.POSIXct(valid, format="%Y-%m-%d %H:%M", tz="UTC"),
    temp_c = round((tmpf - 32) * 5/9, 1),      # Temperatura en Celsius
    rocio_c = round((dwpf - 32) * 5/9, 1),     # Punto de rocío
    humedad_rel = relh,                        # Humedad relativa %
    viento_dir = drct,                         # Dirección viento (0-360)
    viento_kmh = round(sknt * 1.852, 1),       # Velocidad viento en km/h
    presion_mb = mslp,                         # Presión al nivel del mar
    lluvia_mm = p01i * 25.4                    # Pulgadas a Milímetros
  ) %>%
  select(estacion = station, fecha_hora, temp_c, rocio_c, humedad_rel, 
         viento_dir, viento_kmh, presion_mb, lluvia_mm)

# 4. Verificar que tenemos el año entero
print(summary(df_final_2016))
print(table(df_final_2016$estacion)) # Ver cuántas horas hay por estación

