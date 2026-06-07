#-------------------------------------------------------------------------------
#----------------------1) Funtion to clean the air data-------------------------
#------------------------------------------------------------------------------
library(data.table)

limpiar_aire_madrid <- function(dt_bruto, dt_ubica) {
  
  # I create a copy of the original data to avoid modifying it by reference (data.table modifica por referencia, así que es importante no tocar el original)
  dt <- copy(dt_bruto)
  setDT(dt) # We ensure it's a data.table, por si acaso viene como data.frame u otro formato.
  
  # 1) It is necessary to merge the coordinates of the stations with the raw data, in order to have the longitude and latitude of each station in the same table.
  cols_utiles <- dt_ubica[, c("CODIGO_CORTO", "LONGITUD", "LATITUD", "NOM_TIPO")]
  dt <- merge(cols_utiles, dt, by.x = "CODIGO_CORTO", by.y = "ESTACION") 
  setDT(dt)
  
  # Rename the station code column to "ESTACION" for consistency with the rest of the code and dictionaries
  setnames(dt, "CODIGO_CORTO", "ESTACION")
  
  # 2). Drop redundant columns that we won't use in the analysis. 
  cols_borrar <- intersect(names(dt), c("PROVINCIA", "MUNICIPIO", "PUNTO_MUESTREO","Radicación Solar"))
  if(length(cols_borrar) > 0) dt[, (cols_borrar) := NULL]
  
  # 3). Select wich variables we want to keep, based on the dictionaries we have in the global environment. This way we ensure that we are working only with the stations and magnitudes that we are interested in.
  dt <- dt[MAGNITUD %in% as.numeric(names(nombres_magnitudes)) & 
             ESTACION %in% as.numeric(names(nombres_estaciones_aire))]
  
  # 4.) The data base have two columns for each hour: one with the value of the measurement (H01, H02, ..., H24) and another with the validation code (V01, V02, ..., V24).
  #We only want to keep the values that have a validation code of "V" (valid), and set the rest to NA. We do this in a loop to avoid writing 24 lines de código manualmente.
  for (i in 1:24) {
    h_col <- paste0("H", sprintf("%02d", i)) #paste0("H", sprintf("%02d", i)) genera "H01", "H02", ..., "H24" de forma automática
    v_col <- paste0("V", sprintf("%02d", i))
    
    if (v_col %in% names(dt) && h_col %in% names(dt)) {
      filas_invalidas <- which(dt[[v_col]] != "V")
      if (length(filas_invalidas) > 0)
        set(dt, i = filas_invalidas, j = h_col, value = NA_real_)
    }
  }
  # 5.) Change the codes of magnitudes and stations to their names, and create a date column from the year, 
  #month and day columns. We use set() in data.table to avoid the warning of recycling with named vectors.
  
  # unname() evita que data.table dispare el warning de reciclado con vectores nombrados
  
  set(dt, j = "MAGNITUD", value = unname(nombres_magnitudes[as.character(dt$MAGNITUD)]))
  set(dt, j = "ESTACION", value = unname(nombres_estaciones_aire[as.character(dt$ESTACION)]))
  # Our data had separate columns for year, month and day (ANO, MES, DIA). We combine them into a single date column called "FECHA" 
  #using sprintf to format the date as "YYYY-MM-DD" and then converting it to Date class.
  set(dt, j = "FECHA",    value = as.Date(sprintf("%04d-%02d-%02d", dt$ANO, dt$MES, dt$DIA)))
  
  # 6.) We have already created the date column, so we can remove the original year, month and day columns, 
  #as well as the validation code columns (V01 to V24) that we no longer need. We use set() in a loop to avoid writing 24 lines of code manually, 
  #and also to avoid the warning of .mapply when deleting multiple columns with NULL.
  
  cols_v <- names(dt)[grep("^V[0-9]{2}$", names(dt))]
  for (col in c("ANO", "MES", "DIA", cols_v)) set(dt, j = col, value = NULL)
  
  # 7.) Finally, we reshape the data from wide format (one column per hour) to long format (one row per station-date-hour), using melt() from data.table.
  
  h_cols <- paste0("H", sprintf("%02d", 1:24)) # I am generating the name 
  
  dt_long <- melt(dt, 
                  id.vars = c("ESTACION", "LONGITUD", "LATITUD", "NOM_TIPO", "MAGNITUD", "FECHA"), 
                  measure.vars = h_cols,
                  variable.name = "HORA", 
                  value.name = "DATO")
  
  # (Opcional) Si prefieres que la columna HORA sea numérica (1, 2, ..., 24) en vez de texto ("H01"),
  # puedes deserializarla descomentando la siguiente línea:
  # dt_long[, HORA := as.integer(gsub("H", "", HORA))]
  
  return(dt_long) 
}


#-------------------------------------------------------------------------------
#--------------- FUNCTION TO COMBINE MONTHLY FILES INTO A SINGLE ANNUAL FILE ---
#-------------------------------------------------------------------------------
# Reads all CSV or Excel files in a folder (one per month), combines them and
# returns a single data.table sorted by ESTACION.
#
# Arguments:
#   carpeta      – Path to the folder containing the monthly files.
#   sep          – Field separator for CSV files (default ";").
#   patron       – Regex to filter files (default: all .csv and .xlsx/.xls).
#   ordenar_por  – Column(s) to sort the result (default: ESTACION).
#
# Returns a data.table with all months combined, ordered by station.

combinar_meses_anual <- function(carpeta,
                                  sep       = ";",
                                  patron    = "\\.(csv|xlsx|xls)$",
                                  ordenar_por = "ESTACION") {
  
  # -- Dependencies -----------------------------------------------------------
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("El paquete 'data.table' es necesario.")
  
  # -- Discover files ---------------------------------------------------------
  archivos <- list.files(carpeta, pattern = patron, full.names = TRUE,
                         ignore.case = TRUE)
  
  if (length(archivos) == 0)
    stop("No se encontraron archivos CSV ni Excel en: ", carpeta)
  
  message("  -> Encontrados ", length(archivos), " archivo(s) en: ", carpeta)
  
  # -- Read each file ---------------------------------------------------------
  lista_meses <- lapply(archivos, function(ruta) {
    
    ext <- tolower(tools::file_ext(ruta))
    
    dt <- tryCatch({
      if (ext == "csv") {
        data.table::fread(ruta, sep = sep, encoding = "Latin-1")
      } else {
        # Excel: requires readxl
        if (!requireNamespace("readxl", quietly = TRUE))
          stop("El paquete 'readxl' es necesario para leer archivos Excel.")
        data.table::as.data.table(readxl::read_excel(ruta))
      }
    }, error = function(e) {
      warning("Error leyendo ", basename(ruta), ": ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(dt))
      message("    OK  ", basename(ruta), " (", nrow(dt), " filas)")
    
    dt
  })
  
  # -- Drop files that failed to load -----------------------------------------
  lista_meses <- Filter(Negate(is.null), lista_meses)
  
  if (length(lista_meses) == 0)
    stop("No se pudo leer ningún archivo correctamente.")
  
  # -- Combine all months -----------------------------------------------------
  dt_anual <- data.table::rbindlist(lista_meses, use.names = TRUE, fill = TRUE)
  
  # -- Sort by station (and optionally other columns) -------------------------
  cols_orden <- intersect(ordenar_por, names(dt_anual))
  if (length(cols_orden) > 0)
    data.table::setorderv(dt_anual, cols_orden)
  
  message("  -> Total filas combinadas: ", nrow(dt_anual),
          " | Estaciones únicas: ",
          if ("ESTACION" %in% names(dt_anual)) data.table::uniqueN(dt_anual$ESTACION) else "N/A")
  
  return(dt_anual)
}




#--------------------------FUNCIÓN PARA PASAR DE HORAS A DIAS----------------------
#----------------------------------------------------------------------------------
library(data.table)

# Función para agregar datos horarios a diarios con control de NAs
agregar_a_diario <- function(dt, col_grupo = "ESTACION", col_fecha = "FECHA", col_valor = "DATO",
                             col_longitud = "LONGITUD", col_latitud = "LATITUD", umbral_na = 0.2) {
  
  dt <- as.data.table(dt)
  
  # Columnas de agrupación: estación + coordenadas (constantes por estación) + fecha
  cols_by <- c(col_grupo, col_longitud, col_latitud, col_fecha)
  # Solo incluimos las que realmente existen en la tabla
  cols_by <- intersect(cols_by, names(dt))
  
  dt_diario <- dt[, .(
    DATO_DIARIO = if (sum(is.na(.SD[[col_valor]])) / .N >= umbral_na) {
      NA_real_
    } else {
      mean(.SD[[col_valor]], na.rm = TRUE)
    }
  ), by = cols_by]
  
  return(dt_diario)
}

#-------------------------------------------------------------------------------
#------------------- Function to daily to monthly with NA control --------------
#-------------------------------------------------------------------------------
agregar_a_mensual<- function(dt, col_grupo = "ESTACION", col_fecha = "FECHA", col_valor = "DATO",
                             col_longitud = "LONGITUD", col_latitud = "LATITUD", umbral_na = 0.2) {
  
  dt <- as.data.table(dt)
  
  # Aseguramos que la columna de fecha es de tipo Date
  dt[, (col_fecha) := as.Date(get(col_fecha))]
  
  # Extraemos el año y mes para agrupar por mes
  dt[, MES := format(get(col_fecha), "%Y-%m")]
  
  # Columnas de agrupación: estación + coordenadas (constantes por estación) + mes
  cols_by <- c(col_grupo, col_longitud, col_latitud, "MES")
  cols_by <- intersect(cols_by, names(dt))
  
  dt_mensual <- dt[, .(
    DATO_MENSUAL = if (sum(is.na(.SD[[col_valor]])) / .N >= umbral_na) {
      NA_real_
    } else {
      mean(.SD[[col_valor]], na.rm = TRUE)
    }
  ), by = cols_by]
  
  return(dt_mensual)
}






# ==============================================================================
# FUNCIÓN DE LIMPIEZA DE METEOROLOGÍA (Formato Ancho por Covariable)
# ==============================================================================

# ==============================================================================
# FUNCIÓN DE LIMPIEZA DE METEOROLOGÍA (Adaptada a tu CSV real)
# ==============================================================================

limpiar_datos_metereo <- function(dt_bruto, dt_ubica) {
  
  dt <- copy(dt_bruto)
  setDT(dt)
  
  # 1. CRUCE ESPACIAL
  # LONGITUD y LATITUD del CSV de estaciones vienen codificadas como enteros con
  # puntos de miles europeos (ej: -37.122.567 = -3.712257 grados).
  # Transformación: quitar puntos → entero → dividir entre 10^7
  cols_utiles <- dt_ubica[, .(
    CODIGO_CORTO,
    LONGITUD  = as.numeric(gsub("\\.", "", LONGITUD))  / 1e7,
    LATITUD   = as.numeric(gsub("\\.", "", LATITUD))   / 1e7,
    X_km      = as.numeric(gsub(",", ".", COORDENADA_X_ETRS89)) / 1000,
    Y_km      = as.numeric(gsub(",", ".", COORDENADA_Y_ETRS89)) / 1000
  )]
  
  dt <- merge(cols_utiles, dt, by.x = "CODIGO_CORTO", by.y = "ESTACION")
  setDT(dt)
  setnames(dt, "CODIGO_CORTO", "ESTACION")
  
  # 2. Eliminamos columnas redundantes
  cols_borrar <- intersect(names(dt), c("PROVINCIA", "MUNICIPIO", "PUNTO_MUESTREO"))
  if (length(cols_borrar) > 0) for (col in cols_borrar) set(dt, j = col, value = NULL)
  
  # 3. Filtrado con los diccionarios del entorno global
  dt <- dt[MAGNITUD %in% as.numeric(names(nombres_magnitudes_clima))]
  dt <- dt[ESTACION %in% as.numeric(names(nombres_estaciones_clima))]
  
  # 4. Validación (Solo aceptamos código "V")
  for (i in 1:24) {
    h_col <- paste0("H", sprintf("%02d", i))
    v_col <- paste0("V", sprintf("%02d", i))
    if (v_col %in% names(dt) && h_col %in% names(dt)) {
      filas_invalidas <- which(dt[[v_col]] != "V")
      if (length(filas_invalidas) > 0)
        set(dt, i = filas_invalidas, j = h_col, value = NA_real_)
    }
  }
  
  # 5. Traducción y fecha
  set(dt, j = "MAGNITUD", value = unname(nombres_magnitudes_clima[as.character(dt$MAGNITUD)]))
  set(dt, j = "ESTACION", value = unname(nombres_estaciones_clima[as.character(dt$ESTACION)]))
  set(dt, j = "FECHA",    value = as.Date(sprintf("%04d-%02d-%02d", dt$ANO, dt$MES, dt$DIA)))
  
  # =========================================================================
  # 5b. NUEVO: FILTRADO DE ESTACIONES Y ELIMINACIÓN DE RADIACIÓN SOLAR
  # =========================================================================
  estaciones_objetivo <- c("Plaza España", "Ensanche de Vallecas", "Escuelas Aguirre", "Urb. Embajada (Barajas)", 
                           "Arturo Soria", "Plaza Elíptica", "Farolillo", "Sanchinarro", 
                           "Casa de Campo", "El Pardo", "Barajas", "Juan Carlos I", 
                           "Plaza del Carmen", "Tres Olivos", "Moratalaz", "J.M.D Moratalaz", 
                           "Cuatro Caminos", "J.M.D Villaverde", "Barrio del Pilar", "E.D.A.R La China", 
                           "Méndez Álvaro", "Centro Mpal. De Acústica", "Pº Castellana", "J.M.D Hortaleza", 
                           "Retiro", "Peñagrande", "Plaza Castilla", "J.M.D Chamberí", 
                           "J.M.D Centro", "J.M.D Chamartín", "J.M.D Vallecas 1", "J.M.D Vallecas 2", 
                           "Matadero 01", "Matadero 02")
   
  
  # Nos quedamos solo con las estaciones de la lista
  dt <- dt[ESTACION %in% estaciones_objetivo]
  
  # Eliminamos la radiación solar (pasando a minúsculas por si lleva tilde en el diccionario)
  dt <- dt[!tolower(MAGNITUD) %in% c("radiacion solar", "radiación solar")]
  # =========================================================================
  
  # 6. Borrado de columnas sobrantes
  cols_v <- names(dt)[grep("^V[0-9]{2}$", names(dt))]
  for (col in c("ANO", "MES", "DIA", cols_v)) set(dt, j = col, value = NULL)
  
  # 7. Formato largo (melt por hora)
  h_cols <- paste0("H", sprintf("%02d", 1:24))
  dt_long <- melt(dt,
                  id.vars      = c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "MAGNITUD", "FECHA"),
                  measure.vars = h_cols,
                  variable.name = "HORA",
                  value.name   = "VALOR_METEO")
  
  # 8. Formato ancho (una columna por variable climática)
  dt_wide <- dcast(dt_long,
                   ESTACION + LONGITUD + LATITUD + X_km + Y_km + FECHA + HORA ~ MAGNITUD,
                   value.var = "VALOR_METEO")
  
  return(dt_wide)
}

#-----------------------------------------------------------------------------
# ----------------- Cleaning Function for Hourly Traffic Data with Spatial Imputation ----------------
#-----------------------------------------------------------------------------
limpiar_trafico_espacial_horario <- function(dt, dt_ubicaciones, mapa_distritos, mapa_barrios) {
  
  # Create a working copy of the data to avoid modifying the original by reference (data.table modifies by reference, so it's important not to touch the original)
  dt_work <- copy(dt)
  setDT(dt_work)
  
  # 1. Temporal Formatting: Combine date and hour into a single POSIXct column for easier handling of hourly data.
  dt_work[, fecha_hora := as.POSIXct(fecha, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Madrid")] # It passes the date and hour to a single column of type POSIXct, which is more convenient for handling hourly data. We specify the format of the input date and time, and set the timezone to Madrid.
  dt_work[, FECHA := as.Date(fecha_hora)] 
  
  # 2. We only keep these variables 
  cols_medidas <- c("intensidad", "ocupacion", "carga")
  
  # 2b. We set to NA the measurements that have error codes "E" (error) or "S" (suspect), as well as any negative values which are not physically meaningful for these variables.
  for (col in cols_medidas) {
    dt_work[error %in% c("E", "S") | get(col) < 0, (col) := NA]
  }
  
  # 3. Spatial Coordinates: We need to assign UTM coordinates to each sensor based on
  #the provided location data. This is crucial for the spatial imputation step later on.
  
  setDT(dt_ubicaciones)
  coor_medidores <- unique(dt_ubicaciones[, .(id, utm_x, utm_y)], by = "id")
  # We convert the coordinates to an sf object for spatial operations. We only keep sensors that have valid UTM coordinates (non-NA).
  
# We are ignoring sensors without valid UTM coordinates because they cannot be 
  #used for spatial imputation. If we kept them, they would not have a defined
  #location in space, which would cause problems when trying to assign them to districts
  #or neighborhoods, and they would also be excluded from the nearest neighbor imputation step 
  #since they don't have coordinates to calculate distances. By filtering out sensors with NA UTM
  #coordinates at this stage, we ensure that all remaining sensors have valid spatial information 
  #that can be used for the subsequent spatial operations.
  
  sf_sensores <- st_as_sf(coor_medidores[!is.na(utm_x)], 
                          coords = c("utm_x", "utm_y"), 
                          crs = 25830, 
                          remove = FALSE)
  
  cat("Sensores con coordenadas UTM válidas: ", nrow(sf_sensores), "\n")
  cat("Sensores sin coordenadas UTM válidas (excluidos): ", nrow(coor_medidores) - nrow(sf_sensores), "\n")
  
  # 4a. Assignation of districts based on spatial intersection. We transform the district map to 
  #the same CRS as the sensor coordinates (UTM 30N, EPSG:25830) to ensure that the spatial join works 
  #correctly. Then we use st_join with st_intersects to assign each sensor to the district it falls into. 
  #This will add a new column "distrito" to sf_sensores_dt with the name of the district for each sensor.
  mapa_distritos <- st_transform(mapa_distritos, crs = 25830)
  sf_sensores_dt <- st_join(sf_sensores, mapa_distritos[, "distrito"], join = st_intersects)
  
  # 4b. The same procedure but with neighborhood map.
  mapa_barrios <- st_transform(mapa_barrios, crs = 25830)
  sf_sensores_ba <- st_join(sf_sensores, mapa_barrios[, "barrio"], join = st_intersects)
  
  # Combine the relevant columns from both spatial joins into a single data.table for easier merging 
  #with the hourly data. We keep the sensor ID, UTM coordinates, district and neighborhood information.
  dt_sensores_fijos <- as.data.table(sf_sensores_dt)[, .(id, utm_x, utm_y, distrito)]
  dt_barrios_fijos  <- as.data.table(sf_sensores_ba)[, .(id, barrio)]
  dt_sensores_fijos <- merge(dt_sensores_fijos, dt_barrios_fijos, by = "id", all.x = TRUE)
  
  # Cross-setting the spatial information with the hourly data. We merge the spatial information
  #(UTM coordinates, district and neighborhood) into the main hourly data table based on the sensor ID. 
  #This will allow us to perform spatial imputation later on.
  dt_work <- merge(dt_work, dt_sensores_fijos[, .(id, utm_x, utm_y, distrito, barrio)],
                   by = "id", all.x = TRUE)
  
  # Discard rows with missing spatial information (sensors without valid UTM coordinates or without district assignment) 
  #since they cannot be used for spatial imputation or aggregation by district/neighborhood.
  #This ensures that the subsequent analysis is based on sensors with complete spatial information.
  dt_work <- dt_work[!is.na(utm_x) & !is.na(distrito) & distrito != ""]
  
  # 5. Imputation with nearest neighbor spatial interpolation. 
  #For each hour and each variable, we identify the sensors that have missing values.
  
  # Adding a summary of the number of NAs before imputation for each variable, to have a reference of how many values are being imputed. This can help us understand the extent of missing data and the impact of the imputation process.
  cat("Conteo de valores NA",sum(is.na(dt_work$intensidad)), "en intensidad,",
      sum(is.na(dt_work$ocupacion)), "en ocupación y",
      sum(is.na(dt_work$carga)), "en carga.\n")
  
  cat("   -> Iniciando interpolación espacial a nivel horario...\n")
  
  # Tally the total number of NAs across all three variables to determine if we need to perform 
  #imputation. If there are no NAs, we can skip the imputation step entirely, which saves time and 
  #computational resources.
  total_nas <- sum(is.na(dt_work$intensidad)) + sum(is.na(dt_work$ocupacion)) + sum(is.na(dt_work$carga))
  
  if (total_nas > 0) {
    
    cat("   -> Iniciando interpolación espacial a nivel horario...\n")
    
    for (col in cols_medidas) {
      
      horas_con_na <- dt_work[is.na(get(col)), unique(fecha_hora)]
      
      if (length(horas_con_na) > 0) {
        cat("      Imputando", length(horas_con_na), "horas defectuosas en la variable:", col, "\n")
      }
      
      for (h in horas_con_na) {
        
        dt_h <- dt_work[fecha_hora == h, .(id, utm_x, utm_y, valor = get(col))]
        
        faltantes <- dt_h[is.na(valor)]
        validos   <- dt_h[!is.na(valor)]
        
        if (nrow(validos) > 0 && nrow(faltantes) > 0) {
          
          sf_falt <- st_as_sf(faltantes, coords = c("utm_x", "utm_y"), crs = 25830)
          sf_val  <- st_as_sf(validos,   coords = c("utm_x", "utm_y"), crs = 25830)
          
          nn_indices        <- st_nearest_feature(sf_falt, sf_val)
          valores_imputados <- validos$valor[nn_indices]
          
          ids_faltantes <- faltantes$id
          dt_work[fecha_hora == h & id %in% ids_faltantes,
                  (col) := valores_imputados[match(id, ids_faltantes)]]
        }
      }
    }
    cat("   -> ✅ Interpolación completada con éxito.\n")
    
  } else {
    # Si total_nas es 0, salta directamente aquí
    cat("   -> ✅ ¡Datos perfectos! No hay valores nulos, se omite la interpolación.\n")
  }
  # 6. Selection of relevant columns for the output. We keep the original date and time column 
  #(fecha_hora), the date (FECHA), sensor ID, UTM coordinates, district and neighborhood information, 
  #and the traffic variables (intensidad, ocupacion, carga).
  cols_output <- c("fecha_hora", "FECHA", "id", "utm_x", "utm_y",
                   "distrito", "barrio", "intensidad", "ocupacion", "carga")
  dt_horario <- dt_work[, ..cols_output]
  
  return(dt_horario)
}


#------------------------------------------------------------------------------
# FUNCIÓN: CONVERSIÓN FLEXIBLE DE RESOLUCIÓN TEMPORAL PARA CUALQUIER DATASET
#------------------------------------------------------------------------------
# Convierte un data.table de resolución fina (horaria/diaria) a diaria o mensual.
#
# Argumentos:
#   dt           – data.table de entrada
#   a            – destino: "diario" o "mensual"
#   col_fecha    – nombre de la columna de fecha (Date o POSIXct)
#   cols_grupo   – vector de columnas de agrupación (ej. c("ESTACION","LONGITUD","LATITUD"))
#   cols_valores – vector de columnas numéricas a agregar
#   umbral_na    – proporción máxima de NAs permitida antes de devolver NA (default 0.2)
#   fun_agr      – función de agregación aplicada a cada columna (default mean)
#
# Devuelve un data.table con una fila por grupo × período.

convertir_resolucion <- function(dt,
                                 a            = c("diario", "mensual"),
                                 col_fecha    = "FECHA",
                                 cols_grupo   = "ESTACION",
                                 cols_valores = NULL,
                                 umbral_na    = 0.2,
                                 fun_agr      = mean) {
  
  a  <- match.arg(a)
  dt <- as.data.table(dt)
  
  # Aseguramos que la columna de fecha sea Date
  dt[, (col_fecha) := as.Date(get(col_fecha))]
  
  # Creamos la columna de período según la resolución pedida
  col_periodo <- if (a == "diario") "FECHA_AGR" else "MES"
  dt[, (col_periodo) := if (a == "diario") {
    get(col_fecha)
  } else {
    format(get(col_fecha), "%Y-%m")
  }]
  
  # Columnas de agrupación finales
  cols_by <- intersect(c(cols_grupo, col_periodo), names(dt))
  
  # Si no se especifican columnas de valores, tomamos todas las numéricas restantes
  if (is.null(cols_valores)) {
    cols_valores <- setdiff(names(dt), c(cols_grupo, col_fecha, col_periodo))
    cols_valores <- cols_valores[sapply(dt[, ..cols_valores], is.numeric)]
  }
  
  dt_agr <- dt[, lapply(.SD, function(x) {
    if (!is.numeric(x)) return(NA_real_)
    if (sum(is.na(x)) / .N >= umbral_na) NA_real_ else fun_agr(x, na.rm = TRUE)
  }), by = cols_by, .SDcols = cols_valores]
  
  # Renombramos FECHA_AGR -> FECHA para mantener consistencia si es diario
  if (a == "diario" && "FECHA_AGR" %in% names(dt_agr))
    setnames(dt_agr, "FECHA_AGR", "FECHA")
  
  setorderv(dt_agr, cols_by)
  return(dt_agr)
}


#------------------------------------------------------------------------------
# FUNCIÓN: AGREGACIÓN DE DATOS HORARIOS DE TRÁFICO A NIVEL DIARIO POR DISTRITO
#------------------------------------------------------------------------------
# Toma el output de limpiar_trafico_espacial_horario() y agrega por distrito y
# fecha, calculando la media de cada variable e indicando cuántos sensores
# contribuyen a cada combinación.
#
# Argumentos:
#   dt_horario  – data.table con datos horarios limpios (output de la función anterior)
#
# Devuelve un data.table con una fila por (distrito, FECHA).

agregar_trafico_diario <- function(dt_horario) {
  
  dt <- as.data.table(dt_horario)
  
  # Agregación a nivel de DISTRITO
  dt_diario_distrito <- dt[, .(
    intensidad    = mean(intensidad, na.rm = TRUE),
    ocupacion     = mean(ocupacion,  na.rm = TRUE),
    carga         = mean(carga,      na.rm = TRUE),
    num_medidores = uniqueN(id)
  ), by = .(distrito, FECHA)]
  setorder(dt_diario_distrito, FECHA, distrito)
  
  # Agregación a nivel de BARRIO (excluimos sensores sin barrio asignado)
  # Se incluye distrito en el by para conservar la jerarquía barrio → distrito
  dt_diario_barrio <- dt[!is.na(barrio) & barrio != "", .(
    intensidad    = mean(intensidad, na.rm = TRUE),
    ocupacion     = mean(ocupacion,  na.rm = TRUE),
    carga         = mean(carga,      na.rm = TRUE),
    num_medidores = uniqueN(id)
  ), by = .(distrito, barrio, FECHA)]
  setorder(dt_diario_barrio, FECHA, distrito, barrio)
  
  return(list(distrito = dt_diario_distrito, barrio = dt_diario_barrio))
}


#------------------------------------------------------------------------------
# FUNCIÓN: AGREGACIÓN DE DATOS DE TRÁFICO A NIVEL MENSUAL POR DISTRITO Y BARRIO
#------------------------------------------------------------------------------
# Toma el output de agregar_trafico_diario() (lista con $distrito y $barrio) o
# directamente el data.table diario por distrito/barrio, y calcula la media
# mensual de intensidad, ocupación y carga.
#
# Argumentos:
#   dt_diario_distrito – data.table diario por distrito
#   dt_diario_barrio   – data.table diario por barrio
#   umbral_na          – proporción máxima de NAs para devolver NA (default 0.2)
#
# Devuelve una lista con $distrito y $barrio, igual que agregar_trafico_diario().

agregar_trafico_mensual <- function(dt_diario_distrito,
                                    dt_diario_barrio  = NULL,
                                    umbral_na         = 0.2) {
  
  vars_trafico <- c("intensidad", "ocupacion", "carga")
  
  agregar_un_nivel <- function(dt, cols_grupo) {
    dt <- as.data.table(dt)
    dt[, FECHA := as.Date(FECHA)]
    dt[, MES   := format(FECHA, "%Y-%m")]
    
    cols_by <- intersect(c(cols_grupo, "MES"), names(dt))
    
    dt_mens <- dt[, c(
      lapply(.SD, function(x) {
        if (sum(is.na(x)) / .N >= umbral_na) NA_real_
        else mean(x, na.rm = TRUE)
      }),
      list(num_dias = .N)
    ), by = cols_by, .SDcols = vars_trafico]
    
    setorderv(dt_mens, cols_by)
    return(dt_mens)
  }
  
  dt_mens_distrito <- agregar_un_nivel(dt_diario_distrito, cols_grupo = "distrito")
  
  resultado <- list(distrito = dt_mens_distrito)
  
  if (!is.null(dt_diario_barrio)) {
    dt_mens_barrio   <- agregar_un_nivel(dt_diario_barrio,   cols_grupo = c("distrito", "barrio"))
    resultado$barrio <- dt_mens_barrio
  }
  
  return(resultado)
}

#----------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
library(stringi) # ¡NUEVA LIBRERÍA! Para quitar tildes fácilmente

# Función auxiliar para limpiar textos (quita tildes, minúsculas, espacios extra)
limpiar_nombres <- function(texto) {
  texto <- tolower(trimws(texto))
  texto <- stri_trans_general(texto, "Latin-ASCII") # Quita tildes
  # Reemplaza múltiples espacios o guiones raros por un solo espacio
  texto <- gsub("[[:punct:]]+", " ", texto) 
  texto <- gsub("\\s+", " ", texto)
  return(trimws(texto))
}


# ==============================================================================
# FUNCIÓN: Cargar datos crudos meteorológicos para un año dado
#
# Busca primero un CSV anual único con patrón "{anio}_datos_metereo.csv".
# Si no existe, combina todos los CSVs mensuales encontrados recursivamente
# dentro de la carpeta del año.
# ==============================================================================
cargar_datos_metereo_anual <- function(anio, carpeta_base) {
  
  carpeta_anio <- file.path(carpeta_base, as.character(anio))
  
  if (!dir.exists(carpeta_anio))
    stop("No existe la carpeta para el año ", anio, ": ", carpeta_anio)
  
  # Intento 1: archivo CSV anual único
  ruta_anual <- file.path(carpeta_anio, paste0(anio, "_datos_metereo.csv"))
  
  if (file.exists(ruta_anual)) {
    message("  -> Cargando archivo anual único: ", basename(ruta_anual))
    return(fread(ruta_anual, sep = ";"))
  }
  
  # Intento 2: archivos mensuales en subcarpetas
  archivos_csv <- list.files(carpeta_anio, pattern = "\\.csv$",
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = TRUE)
  
  if (length(archivos_csv) == 0)
    stop("No se encontraron archivos CSV en: ", carpeta_anio)
  
  message("  -> Encontrados ", length(archivos_csv),
          " archivo(s) mensual(es). Combinando...")
  
  lista_meses <- lapply(archivos_csv, function(ruta) {
    dt <- tryCatch(
      fread(ruta, sep = ";", encoding = "Latin-1"),
      error = function(e) {
        warning("Error leyendo ", basename(ruta), ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(dt))
      message("    OK  ", basename(ruta), " (", nrow(dt), " filas)")
    dt
  })
  
  lista_meses <- Filter(Negate(is.null), lista_meses)
  
  if (length(lista_meses) == 0)
    stop("No se pudo leer ningún archivo correctamente.")
  
  dt_anual <- rbindlist(lista_meses, use.names = TRUE, fill = TRUE)
  message("  -> Total filas combinadas: ", nrow(dt_anual))
  return(dt_anual)
}

# ==============================================================================
# FUNCIÓN: Procesar un año completo → devuelve lista con $diario y $mensual
# ==============================================================================
procesar_anio_meteo <- function(anio, carpeta_base, ruta_estaciones) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("  Procesando año:", anio, "\n")
  cat(strrep("=", 60), "\n")
  
  # 1. Carga de datos brutos
  cat("⏳ Cargando datos...\n")
  dt_ubicaciones <- as.data.table(
    read.csv(ruta_estaciones, sep = ";", fileEncoding = "latin1")
  )
  data_raw <- cargar_datos_metereo_anual(anio, carpeta_base)
  
  # 2. Limpieza y reestructuración
  cat("⏳ Limpiando y reestructurando...\n")
  datos_horarios <- limpiar_datos_metereo(data_raw, dt_ubicaciones)
  
  cols_clima <- setdiff(
    names(datos_horarios),
    c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "FECHA", "HORA")
  )
  
  # 3. Agregación diaria (umbral 20% NA)
  cat("⏳ Agregando a escala diaria...\n")
  datos_diarios <- datos_horarios[, lapply(.SD, function(x) {
    if (sum(is.na(x)) / .N >= 0.2) NA_real_ else mean(x, na.rm = TRUE)
  }), by = .(ESTACION, LONGITUD, LATITUD, X_km, Y_km, FECHA),
  .SDcols = cols_clima]
  
  datos_diarios <- datos_diarios[year(FECHA) == anio]
  datos_diarios[, ANO := anio]
  setorder(datos_diarios, FECHA)
  
  cat("  ✅ Diario:", nrow(datos_diarios), "filas |",
      uniqueN(datos_diarios$ESTACION), "estaciones |",
      uniqueN(datos_diarios$FECHA), "días\n")
  
  # 4. Agregación mensual (umbral 20% NA)
  cat("⏳ Agregando a escala mensual...\n")
  datos_mensuales <- convertir_resolucion(
    dt           = datos_horarios,
    a            = "mensual",
    col_fecha    = "FECHA",
    cols_grupo   = c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km"),
    cols_valores = cols_clima,
    umbral_na    = 0.2
  )
  
  datos_mensuales <- datos_mensuales[substr(MES, 1, 4) == as.character(anio)]
  datos_mensuales[, ANO := anio]
  setorder(datos_mensuales, MES, ESTACION)
  
  cat("  ✅ Mensual:", nrow(datos_mensuales), "filas |",
      uniqueN(datos_mensuales$ESTACION), "estaciones |",
      uniqueN(datos_mensuales$MES), "meses\n")
  
  return(list(horario=datos_horarios,diario = datos_diarios, mensual = datos_mensuales))
}

