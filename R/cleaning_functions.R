#-------------------------Función de limpieza de la contaminacion--------------------
library(data.table)

limpiar_aire_madrid <- function(dt_bruto, dt_ubica) {
  
  # Usamos copy() por seguridad
  dt <- copy(dt_bruto)
  setDT(dt)
  
  # 1. CRUCE ESPACIAL
  cols_utiles <- dt_ubica[, c("CODIGO_CORTO", "LONGITUD", "LATITUD", "NOM_TIPO")]
  dt <- merge(cols_utiles, dt, by.x = "CODIGO_CORTO", by.y = "ESTACION") 
  setDT(dt)
  
  # Renombramos solo CODIGO_CORTO → ESTACION (el resto ya está en mayúsculas)
  setnames(dt, "CODIGO_CORTO", "ESTACION")
  
  # 2. LIMPIEZA DE COLUMNAS REDUNDANTES
  cols_borrar <- intersect(names(dt), c("PROVINCIA", "MUNICIPIO", "PUNTO_MUESTREO","Radicación Solar"))
  if(length(cols_borrar) > 0) dt[, (cols_borrar) := NULL]
  
  # 3. FILTRADO DE ESTACIONES Y CONTAMINANTES
  dt <- dt[MAGNITUD %in% as.numeric(names(nombres_magnitudes)) & 
             ESTACION %in% as.numeric(names(nombres_estaciones_aire))]
  
  # 4. CRITERIO DE VALIDACIÓN (Código "V")
  for (i in 1:24) {
    h_col <- paste0("H", sprintf("%02d", i))
    v_col <- paste0("V", sprintf("%02d", i))
    
    if (v_col %in% names(dt) && h_col %in% names(dt)) {
      filas_invalidas <- which(dt[[v_col]] != "V")
      if (length(filas_invalidas) > 0)
        set(dt, i = filas_invalidas, j = h_col, value = NA_real_)
    }
  }
  # 5. TRADUCCIÓN Y REESTRUCTURACIÓN TEMPORAL
  # unname() evita que data.table dispare el warning de reciclado con vectores nombrados
  set(dt, j = "MAGNITUD", value = unname(nombres_magnitudes[as.character(dt$MAGNITUD)]))
  set(dt, j = "ESTACION", value = unname(nombres_estaciones_aire[as.character(dt$ESTACION)]))
  set(dt, j = "FECHA",    value = as.Date(sprintf("%04d-%02d-%02d", dt$ANO, dt$MES, dt$DIA)))
  
  # 6. LIMPIEZA FINAL (Basura)
  # set() en bucle evita el warning de .mapply al borrar múltiples columnas con NULL
  cols_v <- names(dt)[grep("^V[0-9]{2}$", names(dt))]
  for (col in c("ANO", "MES", "DIA", cols_v)) set(dt, j = col, value = NULL)
  
  # 7. PASO A FORMATO LARGO (Melt)
  # Generamos el vector con los nombres de las columnas de horas (H01 a H24)
  h_cols <- paste0("H", sprintf("%02d", 1:24))
  
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

# --- EJEMPLO DE USO ---
# Asumiendo que tu tabla larga de NO2 se llama 'dt_no2'
# dt_diario_no2 <- agregar_a_diario(dt_no2, col_grupo = "ESTACION", col_fecha = "FECHA", col_valor = "DATO")



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
  estaciones_objetivo <- c(
    "Plaza Elíptica", "Peñagrande", "Juan Carlos I", "J.M.D Villaverde", 
    "J.M.D Moratalaz", "J.M.D Hortaleza", "Centro Mpal. De Acústica", "Casa de Campo"
  )
  
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
#----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

# ----------------- 1. DEFINICIÓN DE LA FUNCIÓN ESPACIAL DIARIA -----------------
limpiar_trafico_espacial_horario <- function(dt, dt_ubicaciones, mapa_poligonos) {
  
  dt_work <- copy(dt)
  setDT(dt_work)
  
  # 1. Transformación Temporal
  dt_work[, fecha_hora := as.POSIXct(fecha, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Madrid")]
  dt_work[, FECHA := as.Date(fecha_hora)] 
  
  # 2. Limpieza de Errores y Negativos (ELIMINAMOS 'vmed' DE LA LISTA)
  cols_medidas <- c("intensidad", "ocupacion", "carga")
  
  for (col in cols_medidas) {
    dt_work[error %in% c("E", "S") | get(col) < 0, (col) := NA]
  }
  
  # 3. PREPARACIÓN GEOMÉTRICA DE LOS SENSORES
  setDT(dt_ubicaciones)
  coor_medidores <- unique(dt_ubicaciones[, .(id, utm_x, utm_y)], by = "id")
  
  sf_sensores <- st_as_sf(coor_medidores[!is.na(utm_x)], 
                          coords = c("utm_x", "utm_y"), 
                          crs = 25830, 
                          remove = FALSE)
  
  # 4. ASIGNACIÓN GEOMÉTRICA DE DISTRITOS
  mapa_poligonos <- st_transform(mapa_poligonos, crs = 25830)
  sf_sensores_distrito <- st_join(sf_sensores, mapa_poligonos[, "distrito"], join = st_intersects)
  dt_sensores_fijos <- as.data.table(sf_sensores_distrito)
  
  # Cruzamos coordenadas y distrito con los datos horarios
  dt_work <- merge(dt_work, dt_sensores_fijos[, .(id, utm_x, utm_y, distrito)], by = "id", all.x = TRUE)
  
  # Limpiamos datos que caen fuera del mapa o sin coordenadas
  dt_work <- dt_work[!is.na(utm_x) & !is.na(distrito) & distrito != ""]
  
  # 5. IMPUTACIÓN ESPACIAL HORARIA: VECINO MÁS CERCANO
  cat("   -> Iniciando interpolación espacial a nivel horario...\n")
  
  for (col in cols_medidas) {
    
    # Identificamos qué horas exactas tienen algún NA en esta variable
    horas_con_na <- dt_work[is.na(get(col)), unique(fecha_hora)]
    
    if (length(horas_con_na) > 0) {
      cat("      Imputando", length(horas_con_na), "horas defectuosas en la variable:", col, "\n")
    }
    
    for (h in horas_con_na) {
      
      # Filtramos los datos solo de esa hora exacta
      dt_h <- dt_work[fecha_hora == h, .(id, utm_x, utm_y, valor = get(col))]
      
      faltantes <- dt_h[is.na(valor)]
      validos   <- dt_h[!is.na(valor)]
      
      # Si hay al menos un sensor válido para copiar en esa hora
      if (nrow(validos) > 0 && nrow(faltantes) > 0) {
        
        sf_falt <- st_as_sf(faltantes, coords = c("utm_x", "utm_y"), crs = 25830)
        sf_val  <- st_as_sf(validos, coords = c("utm_x", "utm_y"), crs = 25830)
        
        # Matemáticas: ¿Quién está más cerca en este instante?
        nn_indices <- st_nearest_feature(sf_falt, sf_val)
        valores_imputados <- validos$valor[nn_indices]
        
        # Inyección a prueba de fallos
        ids_faltantes <- faltantes$id
        dt_work[fecha_hora == h & id %in% ids_faltantes, (col) := valores_imputados[match(id, ids_faltantes)]]
      }
    }
  }
  
  # 6. AGREGACIÓN DIARIA POR DISTRITO (Sin 'vmed')
  dt_areas_diario <- dt_work[, .(
    intensidad        = mean(intensidad, na.rm = TRUE),
    ocupacion         = mean(ocupacion, na.rm = TRUE),
    carga             = mean(carga, na.rm = TRUE),
    num_medidores     = uniqueN(id)
  ), by = .(distrito, FECHA)]
  
  return(dt_areas_diario)
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

