
#----------------------------Librerías-----------------------------
library(tidyverse)
library(readxl)
library(data.table)
library(dplyr)
library(lubridate)
#----------------------------Lectura de los datos-------------------
ruta_trafico<-"C:\\Users\\HP\\Desktop\\TFM\\Base_de_datos\\Datos_trafico\\Datos_trafico_2025\\Enero_2025\\2025_Enero.csv"
#-------------------------------------------------------------------
limpiar_trafico_dt <- function(dt) {
  
  # Asegurarnos de que R trate el objeto estrictamente como un data.table
  # Esto elimina el aviso de la copia superficial
  setDT(dt)
  
  # 1. Transformación de Fecha
  dt[, fecha := as.POSIXct(fecha, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
  dt[, fecha_hora := floor_date(fecha, "hour")]
  
  # 2. Limpieza de Errores y Negativos 
  dt[error %in% c("E", "S") | intensidad < 0, intensidad := NA]
  dt[error %in% c("E", "S") | ocupacion < 0, ocupacion := NA]
  dt[error %in% c("E", "S") | carga < 0, carga := NA]
  dt[error %in% c("E", "S") | vmed < 0, vmed := NA]
  
  # 3. Agregación Horaria
  dt_horario <- dt[, .(
    intensidad = mean(intensidad, na.rm = TRUE),
    ocupacion  = mean(ocupacion, na.rm = TRUE),
    carga      = mean(carga, na.rm = TRUE),
    vmed       = mean(vmed, na.rm = TRUE),
    registros_validos = sum(!is.na(intensidad))
  ), by = .(id, tipo_elem, fecha_hora)]
  
  # 4. Limpiar los NaN 
  for(col in c("intensidad", "ocupacion", "carga", "vmed")) {
    set(dt_horario, which(is.nan(dt_horario[[col]])), col, NA)
  }
  
  return(dt_horario)
}
#-----------------------------Visualización de los datos---------------------
data_raw_trafico<- fread(ruta_trafico, sep = ";",drop="registros_validos")
names(data_raw_trafico)
head(data_raw_trafico)
#--------------Cogemos un subconjunto para ver si la limpieza se ha hecho de manera correcta---
# Filtra solo unos pocos IDs para probar tu código sin que colapse R---------------------------
ids_prueba <- unique(data_raw_trafico$id)[1:4]
ids_prueba
df_mini <- data_raw_trafico %>% filter(id %in% ids_prueba)
df_mini_limpio <- limpiar_trafico_dt(df_mini)
df_mini_limpio
#------------------------Diagnóstico-----------------------------------------
head(df_mini_limpio$fecha) # Está bien hecho la media 
