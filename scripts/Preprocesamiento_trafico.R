#------------Preprocesamiento de los datos de tráfico--------------------------------------------
library(data.table)

# ----------------- 1. DEFINICIÓN DE LA FUNCIÓN MEJORADA -----------------
limpiar_trafico_dt <- function(dt, dt_ubicaciones = NULL,frecuencia="months") {
  
  # Usamos copy() para evitar modificar el dataframe original
  dt_work <- copy(dt)
  setDT(dt_work)
  
  # 1. Transformación de Fecha
  dt_work[, fecha := as.POSIXct(fecha, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Madrid")]
  
  if (frecuencia == "months") {
    #Eliminamos horas/minutos  y fijamos el día al 1 de cada mes en formato Date
    dt_work[,month_year := as.Date(trunc(fecha,"months"))]
    time_column <- "month_year" # Name of the column to group by 
  } 
  else if (frecuencia == "hours"){
  dt_work[, fecha_hora := as.POSIXct(trunc(fecha, "hours"))]
  time_column<-"date_hour" 
    }
    else {
      stop(" The frequence is nor valid. Choose months or hours")
    }
  
  # 2. Limpieza de Errores y Negativos 
  cols_medidas <- c("intensidad", "ocupacion", "carga", "vmed")
  
  for (col in cols_medidas) {
    dt_work[error %in% c("E", "S") | get(col) < 0, (col) := NA]
  }
  
  # Función auxiliar para calcular medias evitando los NaN
  mean_na <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }
  
  # 3. Agregación Horaria
  dt_horario <- dt_work[, .(
    intensidad        = mean_na(intensidad),
    ocupacion         = mean_na(ocupacion),
    carga             = mean_na(carga),
    vmed              = mean_na(vmed),
    registros_validos = sum(!is.na(intensidad))
  ), by = .(id, tipo_elem)]
  
  # 4. CRUCE ESPACIAL: Añadir las ubicaciones si se proporciona el dataset
  if (!is.null(dt_ubicaciones)) {
    
    # Aseguramos que dt_ubicaciones sea un data.table
    setDT(dt_ubicaciones)
    
    # Definimos las columnas que queremos traernos del archivo de ubicaciones
    # (Excluimos 'tipo_elem' porque ya lo tenemos en dt_horario)
    cols_ubicacion <- c("id", "distrito", "nombre", "utm_x", "utm_y", "longitud", "latitud")
    
    # Verificamos qué columnas existen realmente en el archivo cargado para evitar errores
    cols_existentes <- intersect(cols_ubicacion, names(dt_ubicaciones))
    
    # Hacemos un Left Join basado en el 'id'
    dt_horario <- merge(dt_horario, 
                        dt_ubicaciones[, ..cols_existentes], 
                        by = "id", 
                        all.x = TRUE)
  }
  
  return(dt_horario)
}

# ----------------- 2. CARGA DE DATOS Y DIAGNÓSTICO -----------------
ruta_trafico_2025<-"C:\\Users\\HP\\Desktop\\TFM\\Base_de_datos\\Datos_trafico\\Datos_trafico_2025\\Enero_2025\\Enero_2025.csv"
# 1. Cargar el histórico de tráfico
data_raw_trafico <- fread(ruta_trafico_2025, sep = ";")


# 2. Cargar el fichero mensual de ubicaciones (asegúrate de poner tu ruta correcta)
# Ojo: comprueba si el separador es ";" o "," en este archivo
ruta_ubicaciones_Enero_2025 <- "C:\\Users\\HP\\Desktop\\TFM\\Base_de_datos\\Datos_trafico\\Detectores_2025\\Enero.csv"
data_ubicaciones <- fread(ruta_ubicaciones_Enero_2025, sep = ";")

# 3. Seleccionar unos IDs de prueba
ids_prueba <- unique(data_raw_trafico$id)[1:4]
df_mini <- data_raw_trafico[id %in% ids_prueba]

# 4. Aplicar la función de limpieza Y AÑADIR UBICACIONES
df_mini_limpio <- limpiar_trafico_dt(df_mini, dt_ubicaciones = data_ubicaciones, frecuencia="months")

# ------------------------ Ver Resultados -----------------------------------------
print("Columnas tras el cruce espacial:")
print(names(df_mini_limpio))

print("Primeras filas con datos georreferenciados:")
print(head(df_mini_limpio))

print("Número de valores faltantes (NAs) por columna:")
print(data_raw_trafico[, lapply(.SD, function(x) sum(is.na(x)))])


#----------------------- Tablas limpias mes y horas ------------------------------------------


# 1. Definimos los meses que quieres procesar
meses <- c("Enero","Febrero", "Marzo", "Abril", "Mayo", "Junio", 
           "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre")

# 2. Creamos una carpeta nueva para guardar los archivos limpios (cámbiala si quieres)
carpeta_salida <- "C:\\Users\\HP\\Desktop\\TFM\\Base_de_datos\\Datos_trafico\\Datos_limpios_2025\\"

# Si la carpeta no existe, R la creará por ti automáticamente
if(!dir.exists(carpeta_salida)) dir.create(carpeta_salida)

# 3. Iniciamos el Bucle para procesar mes a mes
for (mes in meses) {
  
  cat("\n--------------------------------------------------\n")
  print(paste("Iniciando procesamiento del mes:", mes))
  
  # A) Construimos las rutas dinámicamente usando 'paste0'
  # Ajusta las rutas si tuvieran alguna diferencia en mayúsculas/minúsculas
  ruta_trafico_mes <- paste0("C:\\Users\\HP\\Desktop\\TFM\\Base_de_datos\\Datos_trafico\\Datos_trafico_2025\\", mes, "_2025\\", mes, "_2025.csv")
  ruta_ubica_mes   <- paste0("C:\\Users\\HP\\Desktop\\TFM\\Base_de_datos\\Datos_trafico\\Detectores_2025\\", mes, ".csv")
  
  # B) Comprobamos si los archivos existen en tu ordenador para que el bucle no dé error
  if (file.exists(ruta_trafico_mes) & file.exists(ruta_ubica_mes)) {
    
    # 1. Leer los datos brutos
    dt_raw <- fread(ruta_trafico_mes, sep = ";")
    dt_ubi <- fread(ruta_ubica_mes, sep = ";")
    
    # 2. Aplicar tu función de limpieza y cruce espacial
    dt_limpio <- limpiar_trafico_dt(dt_raw, dt_ubicaciones = dt_ubi)
    
    # 3. Guardar el archivo limpio. 
    # Usamos saveRDS() porque conserva el formato exacto de fecha (POSIXct) y comprime muy bien.
    ruta_guardado <- paste0(carpeta_salida, "Trafico_Limpio_", mes, "_2025.rds")
    saveRDS(dt_limpio, ruta_guardado)
    
    # Si prefieres guardarlo en CSV tradicional, descomenta la siguiente línea:
    # fwrite(dt_limpio, paste0(carpeta_salida, "Trafico_Limpio_", mes, "_2025.csv"), sep = ";")
    
    print(paste("ÉXITO: Mes de", mes, "limpiado y guardado en", ruta_guardado))
    
    # 4. PASO CRÍTICO: Liberar memoria RAM
    # Borramos los objetos pesados del entorno
    rm(dt_raw, dt_ubi, dt_limpio)
    # Forzamos a R a que pase el "camión de la basura" (Garbage Collector) para vaciar la RAM
    gc()
    
  } else {
    print(paste("⚠️ AVISO: Faltan archivos (tráfico o medidores) para el mes de", mes, ". Saltando al siguiente..."))
  }
}

print("PROCESAMIENTO TOTAL COMPLETADO.")


