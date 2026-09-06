#===============================================================================
# Functions for cleaning and processing raw data from Madrid's air quality and meteorological stations.
#===============================================================================
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
                             col_longitud = "LONGITUD", col_latitud = "LATITUD", umbral_na = 0.3) {
  
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
                             col_longitud = "LONGITUD", col_latitud = "LATITUD", umbral_na = 0.3) {
  
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
# Function to clean hourly meteorological data with spatial imputation and filtering
# ==============================================================================

limpiar_datos_metereo <- function(dt_bruto, dt_ubica) {
  
  dt <- copy(dt_bruto)
  setDT(dt)
  
  # 1. Merge with station coordinates and convert to numeric
  # We extract the useful columns from the station locations table and convert them 
  # to numeric, handling any formatting issues (e.g., commas, dots).
  
  cols_utiles <- dt_ubica[, .(
    CODIGO_CORTO,
    LONGITUD  = as.numeric(gsub("\\.", "", LONGITUD))  / 1e7,
    LATITUD   = as.numeric(gsub("\\.", "", LATITUD))   / 1e7,
    X_km      = as.numeric(gsub(",", ".", COORDENADA_X_ETRS89)) / 1000,
    Y_km      = as.numeric(gsub(",", ".", COORDENADA_Y_ETRS89)) / 1000
  )]
  
  # We merge the useful columns with the raw data, using the station code as the key.
  
  dt <- merge(cols_utiles, dt, by.x = "CODIGO_CORTO", by.y = "ESTACION")
  setDT(dt)
  setnames(dt, "CODIGO_CORTO", "ESTACION")
  
  # 2. Remove redundant columns that we won't use in the analysis. We check if they 
  # exist before trying to delete them.
  cols_borrar <- intersect(names(dt), c("PROVINCIA", "MUNICIPIO", "PUNTO_MUESTREO"))
  if (length(cols_borrar) > 0) for (col in cols_borrar) set(dt, j = col, value = NULL)
  
  # 3. Filter to keep only the magnitudes and stations we are interested in,
  # based on the dictionaries defined in the global environment. 
  # This ensures that we are working only with relevant data.
  dt <- dt[MAGNITUD %in% as.numeric(names(nombres_magnitudes_clima))]
  dt <- dt[ESTACION %in% as.numeric(names(nombres_estaciones_clima))]
  
  # 4. Set to NA the hourly values that do not have a valid validation code ("V").
  for (i in 1:24) {
    h_col <- paste0("H", sprintf("%02d", i))
    v_col <- paste0("V", sprintf("%02d", i))
    # We check if both the value column and the validation column exist before trying to access them.
    if (v_col %in% names(dt) && h_col %in% names(dt)) {
      filas_invalidas <- which(dt[[v_col]] != "V")
      if (length(filas_invalidas) > 0)
        set(dt, i = filas_invalidas, j = h_col, value = NA_real_)
    }
  }
  
  # 5. Change the codes of magnitudes and stations to their names,
  # and create a date column from the year, month and day columns. 
  # We use set() in data.table to avoid the warning of recycling with named vectors.
  set(dt, j = "MAGNITUD", value = unname(nombres_magnitudes_clima[as.character(dt$MAGNITUD)]))
  set(dt, j = "ESTACION", value = unname(nombres_estaciones_clima[as.character(dt$ESTACION)]))
  set(dt, j = "FECHA",    value = as.Date(sprintf("%04d-%02d-%02d", dt$ANO, dt$MES, dt$DIA)))
  
 
  # 5b. Filter to keep only the stations of interest (those that are relevant for our analysis).
  estaciones_objetivo <- c("Plaza España", "Ensanche de Vallecas", "Escuelas Aguirre", "Urb. Embajada (Barajas)", 
                           "Arturo Soria", "Plaza Elíptica", "Farolillo", "Sanchinarro", 
                           "Casa de Campo", "El Pardo", "Barajas", "Juan Carlos I", 
                           "Plaza del Carmen", "Tres Olivos", "Moratalaz", "J.M.D Moratalaz", 
                           "Cuatro Caminos", "J.M.D Villaverde", "Barrio del Pilar", "E.D.A.R La China", 
                           "Méndez Álvaro", "Centro Mpal. De Acústica", "Pº Castellana", "J.M.D Hortaleza", 
                           "Retiro", "Peñagrande", "Plaza Castilla", "J.M.D Chamberí", 
                           "J.M.D Centro", "J.M.D Chamartín", "J.M.D Vallecas 1", "J.M.D Vallecas 2", 
                           "Matadero 01", "Matadero 02")
   
  
  # Filtering the data to keep only the stations of interest. 
  # This is important to avoid processing unnecessary data and 
  # to focus on the relevant stations for our analysis.
  dt <- dt[ESTACION %in% estaciones_objetivo]
  
  # Exclude variables that are not relevant for our analysis (e.g., wind direction).
  # Estos nombres deben coincidir con los de nombres_magnitudes_clima en
  # minusculas. Antes la lista decia "radiacion.uv", que no coincidia con
  # "Radiacion Ultravioleta", asi que la UV nunca llegaba a excluirse.
  vars_excluir <- c("direccion_viento", "radiacion_ultravioleta")
  dt <- dt[!tolower(MAGNITUD) %in% vars_excluir]
 
  
  # 6. Remove the original year, month and day columns, as well as the validation 
  # code columns (V01 to V24) that we no longer need.
  cols_v <- names(dt)[grep("^V[0-9]{2}$", names(dt))]
  for (col in c("ANO", "MES", "DIA", cols_v)) set(dt, j = col, value = NULL)
  
  # 7. Reshape the data from wide format (one column per hour) to long format 
  # (one row per station-date-hour), using melt() from data.table.
  h_cols <- paste0("H", sprintf("%02d", 1:24))
  dt_long <- melt(dt,
                  id.vars      = c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "MAGNITUD", "FECHA"),
                  measure.vars = h_cols,
                  variable.name = "HORA",
                  value.name   = "VALOR_METEO")
  
  # 8. Reshape back to wide format, with one column per magnitude, using dcast() from data.table.
  dt_wide <- dcast(dt_long,
                   ESTACION + LONGITUD + LATITUD + X_km + Y_km + FECHA + HORA ~ MAGNITUD,
                   value.var = "VALOR_METEO")
  
  return(dt_wide)
}

# ==============================================================================
#  Function to impute missing values (NAs) in hourly meteorological data, 
# using a two-step strategy:
# ==============================================================================
# Strategy for imputing missing values (NAs) in hourly meteorological data:
#
#   1. Linear interpolation: for each station, the missing values are filled 
#      using linear interpolation
#
#   2. Fallback to nearest station: for any remaining missing values, the value 
#      from the geographically closest station (by Euclidean distance in km) that has 
#      a valid reading at the same timestamp is used.
#
# Note: A station is considered NOT to be measuring a variable when 100% of its values ~
# are NA (e.g., a station without a rain gauge). These stations are completely ignored, 
# and their NAs are left intact.

# Arguments:
#   dt_horario – data.table from limpieza de datos 
#
#   cols_clima – Vector of climatic column names to be imputed. If NULL, 
#                columns are automatically detected as all numeric columns not included 
#                in the identifier set.
#
#   maxgap     – Maximum number of consecutive NA hours that can be filled via linear interpolation (default is 3).
#
# ==============================================================================

# ==============================================================================
# Inventario oficial de sensores por estación
# ==============================================================================
# El fichero estaciones.csv del Ayuntamiento incluye una columna por sensor
# ("VV (81)", "T (83)", "HR (86)"...) con una "X" cuando la estación dispone de
# ese aparato. Esa es la fuente de verdad para distinguir un hueco ESTRUCTURAL
# (la estación no mide esa variable) de un hueco TEMPORAL (mide, pero falta el
# dato). Sin este inventario ambos casos son un NA indistinguible.
#
# Devuelve un data.table ESTACION x Variable con la marca tiene_sensor.
#
# Nota: read.csv convierte "VV (81)" en "VV..81.", así que las columnas se
# localizan por la abreviatura MÁS el código de magnitud, y no por el nombre
# literal. Se evita buscar solo el código porque "COORDENADA_Y_ETRS89" también
# contiene "89".
# ==============================================================================

# Estados posibles de una celda climática, en orden de calidad decreciente.
ESTADOS_CLIMA <- c("OK", "IMPUTADO", "FALLO", "AUSENTE", "SIN_SENSOR")

# Proporción de horas imputadas a partir de la cual un periodo agregado (día o
# mes) se considera IMPUTADO.
#
# NO se usa "alguna hora imputada": ese criterio se acumula al subir de escala.
# Con un 0,1 % de horas imputadas en 2019, encadenar "algún" hora -> día y
# "algún" día -> mes producía un 2,5 % de días y un 36,5 % de meses marcados
# como imputados. El dato es el mismo; lo que se degrada es el indicador. Por
# eso cada escala cuenta SUS PROPIAS horas y ninguna hereda de la anterior.
UMBRAL_IMPUTADO_CLIMA <- 0.3

inventario_sensores_clima <- function(dt_ubica) {

  dt_ubica <- as.data.table(dt_ubica)

  abrev_a_codigo <- c(VV = "81", DV = "82", T = "83", HR = "86",
                      PB = "87", RS = "88", P = "89")

  filas <- list()

  for (ab in names(abrev_a_codigo)) {

    codigo <- abrev_a_codigo[[ab]]

    # La variable debe estar en el diccionario; si no, no se modela.
    if (!codigo %in% names(nombres_magnitudes_clima)) next
    variable <- nombres_magnitudes_clima[[codigo]]

    patron <- paste0("^", ab, "[^A-Za-z0-9]*", codigo)
    col <- grep(patron, names(dt_ubica), value = TRUE)

    if (length(col) != 1L) {
      warning("No se localiza una única columna de sensor para ", variable,
              " (patrón '", patron, "'). Se omite del inventario.")
      next
    }

    filas[[length(filas) + 1]] <- data.table(
      CODIGO_CORTO = dt_ubica$CODIGO_CORTO,
      Variable     = variable,
      tiene_sensor = toupper(trimws(as.character(dt_ubica[[col]]))) == "X"
    )
  }

  if (length(filas) == 0L)
    stop("No se ha podido construir el inventario de sensores desde estaciones.csv")

  inv <- rbindlist(filas)

  # Se traduce el código de estación al nombre, que es la clave que usan los
  # datos horarios una vez limpiados.
  inv[, ESTACION := unname(nombres_estaciones_clima[as.character(CODIGO_CORTO)])]
  inv <- inv[!is.na(ESTACION)]

  inv[, .(ESTACION, Variable, tiene_sensor)]
}


imputar_na_horario <- function(dt_horario,
                               cols_clima = NULL, # It is optional
                               maxgap     = 3L,   # maxgap is the maximum number of consecutive NA
                               usar_vecino = FALSE,
                               marcar_imputados = FALSE,
                               marcar_estado = FALSE,
                               inventario = NULL) {
  # hours that can be filled via linear interpolation (default is 3).
  #
  # marcar_imputados: si es TRUE, añade una columna lógica "imp_<variable>" por
  # cada variable climática, con TRUE en las celdas que estaban vacías (NA) en el
  # dato original y que la imputación ha rellenado. Sirve para poder distinguir
  # después las medidas REALES de las imputadas (p. ej. para validar la
  # interpolación solo sobre medidas reales).
  
  # If the 'zoo' package is not installed, we stop the function and inform the user to install it.
  if (!requireNamespace("zoo", quietly = TRUE))
    stop("El paquete 'zoo' es necesario para la interpolación lineal. ",
         "Instálalo con install.packages('zoo').")

  dt <- copy(dt_horario)
  setDT(dt)
  
  # Define the identifier columns that are not to be imputed. 
  # These columns are used to identify each unique observation (station, coordinates, date, hour).
  id_cols <- c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km",
               "FECHA", "HORA")
  
  # If cols_clima is NULL, we automatically detect the climate variable columns 
  # as all numeric columns not included in the identifier set.
  if (is.null(cols_clima)) {
    cols_clima <- setdiff(names(dt), id_cols)
    cols_clima <- cols_clima[sapply(dt[, ..cols_clima], is.numeric)]
  }

  if (!all(c("X_km", "Y_km") %in% names(dt)))
    stop("Las columnas X_km e Y_km son necesarias para la imputación espacial.")

  # Numeric hour for correct ordering (H01 → 1, …, H24 → 24)
  dt[, HORA_NUM := as.integer(gsub("H", "", as.character(HORA)))]
  setorder(dt, ESTACION, FECHA, HORA_NUM)

  # --- Etiqueta de procedencia de cada celda ---------------------------------
  # El valor numérico NO se toca: sigue siendo el número o NA. La razón por la
  # que una celda está vacía viaja en una columna paralela "<variable>_estado".
  # Se evita a propósito un valor centinela (-999) dentro de la columna
  # numérica, porque is.na() dejaría de verlo y contaminaría medias y umbrales.
  #
  #   OK         medida real del sensor
  #   IMPUTADO   estaba vacía y la imputación le puso valor
  #   FALLO      estaba vacía y sigue vacía
  #   SIN_SENSOR la estación no dispone de ese aparato (según estaciones.csv)
  #
  # na_inicial es la foto de qué celdas estaban vacías ANTES de imputar esa
  # variable; dt mantiene un orden de filas fijo, así que sigue alineada.
  # Las filas anadidas para completar la rejilla vienen marcadas; si el llamante
  # no aporta la columna, se asume que todas las filas son reales.
  if (!"fila_ausente" %in% names(dt)) dt[, fila_ausente := FALSE]

  asignar_estado <- function(v, na_inicial) {

    if (!marcar_estado) return(invisible(NULL))

    # REGLA DE AUSENCIA DE SENSOR (empírica, no de metadatos):
    # una estación no dispone del aparato si NO aporta ni una sola lectura real
    # de esa variable en todo el año. Se evalúa sobre na_inicial, es decir ANTES
    # de imputar, para que el relleno no borre la evidencia.
    #
    # Se descarta usar estaciones.csv como fuente: ese fichero describe la red
    # ACTUAL y contradice a los años antiguos (en 2019 da E.D.A.R La China por
    # estación de solo viento cuando reporta temperatura, humedad, presión y
    # radiación). El dato es la única fuente coherente año a año.
    #
    # LIMITACIÓN CONOCIDA: la regla es anual, así que un sensor retirado a
    # mitad de año (los 5 barómetros de septiembre de 2022) no cuenta como
    # SIN_SENSOR sino como FALLO durante los meses restantes.
    est_sin_sensor <- dt[, .(vacia = all(na_inicial[.I])),
                         by = ESTACION][vacia == TRUE, ESTACION]

    tiene <- !dt$ESTACION %in% est_sin_sensor

    # SIN_SENSOR manda sobre AUSENTE: si la estacion no tiene el aparato, lo
    # informativo es eso, no que falte la fila de ese dia.
    dt[, (paste0(v, "_estado")) := factor(
      fifelse(!tiene,                 "SIN_SENSOR",
      fifelse(fila_ausente %in% TRUE, "AUSENTE",
      fifelse(!na_inicial,            "OK",
      fifelse(!is.na(get(v)),         "IMPUTADO", "FALLO")))),
      levels = ESTADOS_CLIMA
    )]

    # Contraste con el inventario oficial, si se aporta. Es solo diagnóstico:
    # NO decide la etiqueta, pero avisa de discrepancias entre la red declarada
    # y la observada, que son en sí un resultado del capítulo de calidad.
    if (!is.null(inventario) && v %in% inventario$Variable) {
      con_sensor <- inventario[Variable == v & tiene_sensor == TRUE, ESTACION]
      solo_dato <- setdiff(setdiff(dt[, unique(ESTACION)], est_sin_sensor), con_sensor)
      solo_meta <- setdiff(con_sensor, setdiff(dt[, unique(ESTACION)], est_sin_sensor))
      if (length(solo_dato) || length(solo_meta))
        message("   · discrepancia inventario/dato en ", v,
                if (length(solo_dato)) paste0(" | miden sin figurar: ", paste(solo_dato, collapse = ", ")) else "",
                if (length(solo_meta)) paste0(" | figuran sin medir: ", paste(solo_meta, collapse = ", ")) else "")
    }

    invisible(NULL)
  }

  if (!is.null(inventario)) inventario <- as.data.table(inventario)

  cat(if (isTRUE(usar_vecino))
        "⏳ Imputación de NAs horarios (lineal + estación más cercana)...
"
      else
        "⏳ Imputación de NAs horarios (solo interpolación lineal intra-estación)...
")
  
  # Loop over each climate variable to impute missing values
  for (v in cols_clima) {

    # Snapshot of which cells were originally missing (before imputing this var).
    # dt keeps a fixed row order during the whole loop (setorder above), so this
    # logical vector stays aligned with the column when we build the flag below.
    na_inicial <- is.na(dt[[v]])
    na_antes <- sum(na_inicial)
    if (na_antes == 0L) {
      if (marcar_imputados) dt[, (paste0("imp_", v)) := FALSE]
      asignar_estado(v, na_inicial)
      next
    }

    # Identify stations that truly measure this variable (< 100% NA)
    est_mide <- dt[, .(pct_na = sum(is.na(get(v))) / .N), by = ESTACION][
      pct_na < 1, ESTACION]

    if (length(est_mide) == 0L) {
      if (marcar_imputados) dt[, (paste0("imp_", v)) := FALSE]
      asignar_estado(v, na_inicial)
      next
    }

    # --- Step 1: linear interpolation (maxgap) within each station ----------
    dt[ESTACION %in% est_mide,
       (v) := zoo::na.approx(get(v), maxgap = maxgap, na.rm = FALSE),
       by = ESTACION]

    na_post_interp <- sum(is.na(dt[[v]]))
    n_interp <- na_antes - na_post_interp

    # --- Paso 2 (DESACTIVADO POR DEFECTO): estacion mas cercana -------------
    # Copia el valor de la estacion valida mas proxima en el mismo instante.
    #
    # POR QUE ESTA DESACTIVADO. Esto NO es imputacion temporal: es una
    # interpolacion espacial 1-NN sin ponderar por distancia, hecha ANTES de la
    # interpolacion espacial buena (IDW/kriging) que resuelve el
    # desalineamiento. Interpolar dos veces en el espacio tiene tres efectos
    # indeseables para un modelo SPDE:
    #   1. Doble conteo: si la estacion donante tambien esta entre los vecinos
    #      del IDW, su valor entra dos veces en la media ponderada.
    #   2. Validacion circular: un LOOCV sobre celdas rellenadas por 1-NN
    #      favorece artificialmente a 1-NN e IDW frente al kriging.
    #   3. Correlacion espacial fabricada: copiar entre vecinos sesga al alza el
    #      rango espacial que el SPDE trata precisamente de estimar.
    #
    # Y lo que se gana es poco: a escala DIARIA el umbral del 30 % de NA ya
    # absorbe los huecos (se pierden 1-3 dias al ano). El coste es mayor solo a
    # escala horaria en presion y radiacion, que tienen 8 estaciones frente a un
    # umbral IDW de 7.
    #
    # Se conserva tras el argumento usar_vecino para poder reactivarlo y
    # comparar (p. ej. el rango espacial posterior con y sin este paso).
    if (isTRUE(usar_vecino)) {
    
    # We create a table with the coordinates of the stations that measure this variable, 
    # to avoid recalculating distances for each timestamp.
    coords_est <- unique(dt[ESTACION %in% est_mide, .(ESTACION, X_km, Y_km)])
    
    dt_na_rows <- dt[ESTACION %in% est_mide & is.na(get(v)), .(ESTACION, FECHA, HORA)]
    
    # If there are any timestamps with NA values, we proceed to find the nearest station with a valid value.
    if (nrow(dt_na_rows) > 0) {
      timestamps_na <- unique(dt_na_rows[, .(FECHA, HORA)])
      
      # For each timestamp with NA values, we find the nearest station with a valid value and impute it.
      for (i in seq_len(nrow(timestamps_na))) {
        f <- timestamps_na$FECHA[i]
        h <- timestamps_na$HORA[i]
        
        # We create a temporary table with the values of the variable v for all stations 
        # at the current timestamp (FECHA × HORA).
        dt_ts <- dt[ESTACION %in% est_mide & FECHA == f & HORA == h,
                    .(ESTACION, val = get(v))]

        missing_ests <- dt_ts[is.na(val), ESTACION]
        valid_dt     <- merge(dt_ts[!is.na(val)], coords_est, by = "ESTACION")
        
        # If there are no missing stations or no valid stations, we skip to the next timestamp.
        if (length(missing_ests) == 0L || nrow(valid_dt) == 0L) next
        
        # For each station with a missing value, we calculate the Euclidean distance to all valid 
        # stations and take the value of the closest one.
        for (est in missing_ests) {
          cx <- coords_est[ESTACION == est, X_km]
          cy <- coords_est[ESTACION == est, Y_km]
          
          # If the coordinates of the station are missing or NA, we skip to the next station.
          if (length(cx) == 0L || is.na(cx)) next
          
          # Calculate distances to all valid stations and find the closest one
          dists    <- sqrt((valid_dt$X_km - cx)^2 + (valid_dt$Y_km - cy)^2)
          
          # We take the value of the closest station (the one with the minimum distance) 
          # and assign it to the missing station.
          best_val <- valid_dt$val[which.min(dists)]

          dt[ESTACION == est & FECHA == f & HORA == h, (v) := best_val]
        }
      }
    }
    
    }  # fin de if (isTRUE(usar_vecino))

    # After both steps, we count the remaining NAs for this variable and print 
    # a summary of the imputation process.
    na_post_spatial <- sum(is.na(dt[[v]]))
    n_spatial <- na_post_interp - na_post_spatial

    # Flag: celdas que estaban vacías y ahora tienen valor -> fueron imputadas.
    if (marcar_imputados)
      dt[, (paste0("imp_", v)) := na_inicial & !is.na(get(v))]

    asignar_estado(v, na_inicial)

    if (isTRUE(usar_vecino)) {
      cat(sprintf("   %-22s | antes: %6d NA → interp: -%d, est. cercana: -%d → quedan: %d NA
",
                  v, na_antes, n_interp, n_spatial, na_post_spatial))
    } else {
      cat(sprintf("   %-22s | antes: %6d NA → interp lineal: -%d → quedan: %d NA
",
                  v, na_antes, n_interp, na_post_spatial))
    }
  }

  dt[, HORA_NUM := NULL]
  cat("   ✅ Imputación completada.\n")

  return(dt)
}


# ==============================================================================
# Function to load annual meteorological data, either from a single annual CSV or 
# by combining monthly CSVs
#
# First, it checks if a single annual CSV exists. If not, it looks for monthly 
# CSVs in subfolders and combines them.
# ==============================================================================

cargar_datos_metereo_anual <- function(anio, carpeta_base) {
  
  # Check if the base folder exists
  carpeta_anio <- file.path(carpeta_base, as.character(anio))
  
  if (!dir.exists(carpeta_anio))
    stop("No existe la carpeta para el año ", anio, ": ", carpeta_anio)
  
  # First attempt: single annual CSV
  ruta_anual <- file.path(carpeta_anio, paste0(anio, "_datos_metereo.csv"))
  
  # If the annual CSV exists, load it directly
  if (file.exists(ruta_anual)) {
    message("  -> Cargando archivo anual único: ", basename(ruta_anual))
    return(fread(ruta_anual, sep = ";"))
  }
  
  # Second attempt: combine monthly CSVs
  archivos_csv <- list.files(carpeta_anio, pattern = "\\.csv$",
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = TRUE)
  
  # If no monthly CSVs are found, throw an error
  if (length(archivos_csv) == 0)
    stop("No se encontraron archivos CSV en: ", carpeta_anio)
  
  message("  -> Encontrados ", length(archivos_csv),
          " archivo(s) mensual(es). Combinando...")
  
  # Read each monthly CSV and combine them into a single data.table
  # Fread is used with tryCatch to handle any errors in reading individual files, 
  # and a warning is issued for any file that fails to load.
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
  
  # Remove any NULL entries from the list of monthly data.tables (files that failed to load)
  lista_meses <- Filter(Negate(is.null), lista_meses)
  
  if (length(lista_meses) == 0)
    stop("No se pudo leer ningún archivo correctamente.")
  
  dt_anual <- rbindlist(lista_meses, use.names = TRUE, fill = TRUE)
  message("  -> Total filas combinadas: ", nrow(dt_anual))
  return(dt_anual)
}

# ==============================================================================
# Function to process meteorological data for a given year
# ==============================================================================

# ==============================================================================
# Estado de un periodo agregado (día o mes) a partir de las horas que lo componen
# ==============================================================================
# Recibe el horario ya etiquetado y devuelve una tabla con una columna
# <variable>_estado por cada variable climática, resumida al periodo indicado.
#
#   SIN_SENSOR  todas las horas del periodo carecen de aparato
#   FALLO       faltan >= umbral_na de las horas -> el agregado sale NA
#   IMPUTADO    >= UMBRAL_IMPUTADO_CLIMA de las horas CON dato son imputadas
#   OK          resto
#
# Las dos escalas agregadas cuentan horas directamente, sin encadenarse entre
# sí, para que el criterio no se acumule.
estado_periodo_clima <- function(dt_horario, cols_clima, cols_grupo,
                                 umbral_na = 0.3) {

  cols_estado <- intersect(paste0(cols_clima, "_estado"), names(dt_horario))
  if (length(cols_estado) == 0) return(NULL)

  dt_horario[, lapply(.SD, function(s) {
    s       <- as.character(s)
    n       <- length(s)
    n_sin   <- sum(s == "SIN_SENSOR")
    n_aus   <- sum(s == "AUSENTE")
    n_fallo <- sum(s == "FALLO")
    n_imp   <- sum(s == "IMPUTADO")
    n_ok    <- sum(s == "OK")

    if (n_sin == n) {
      "SIN_SENSOR"
    } else if (n_aus == n) {
      "AUSENTE"
    } else if ((n_fallo + n_aus + n_sin) / n >= umbral_na) {
      "FALLO"
    } else if (n_imp / max(n_ok + n_imp, 1L) >= UMBRAL_IMPUTADO_CLIMA) {
      "IMPUTADO"
    } else {
      "OK"
    }
  }), by = cols_grupo, .SDcols = cols_estado]
}


procesar_anio_meteo <- function(anio, carpeta_base, ruta_estaciones) {
  
  # 1. Load raw data
  cat("⏳ Cargando datos...\n")
  dt_ubicaciones <- as.data.table(
    read.csv(ruta_estaciones, sep = ";")
  )

  # En Windows, forzar fileEncoding = "latin1" puede truncar silenciosamente
  # este CSV en la primera estacion cuando la locale nativa ya es CP1252. Una
  # lectura parcial haria que todo el ano se redujera a una sola estacion.
  if (nrow(dt_ubicaciones) < 20L) {
    stop(
      "El inventario de estaciones solo contiene ", nrow(dt_ubicaciones),
      " filas; revisa la codificacion de ", basename(ruta_estaciones), "."
    )
  }

  # Inventario oficial de sensores: permite marcar el hueco ESTRUCTURAL antes de
  # imputar, en lugar de deducirlo después a partir de los NA.
  inventario <- inventario_sensores_clima(dt_ubicaciones)

  # Load the annual meteorological data, either from a single CSV or by combining monthly CSVs
  data_raw <- cargar_datos_metereo_anual(anio, carpeta_base)
  
  # 2. Cleaning and restructuring
  cat("Limpiando y reestructurando...\n")
  datos_horarios <- limpiar_datos_metereo(data_raw, dt_ubicaciones)
  
  # 2a. Detect columns to impute (all numeric except identifiers)
  # Setdiff is used to exclude the identifier columns from the list of climate variable columns.
  cols_clima <- setdiff(
    names(datos_horarios),
    c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "FECHA", "HORA",
      "fila_ausente")
  )
  
  # 2a-bis. Rejilla completa estacion x dia x hora ---------------------------
  # El proveedor no publica fila cuando la estacion esta caida el dia entero:
  # no escribe 24 registros vacios, sencillamente no escribe nada. Esos dias no
  # son NA, no existen, y desaparecen del denominador sin dejar rastro. Aqui se
  # materializan como filas con las variables a NA.
  #
  # El right join contra la rejilla conserva las filas presentes y crea las que
  # faltan: no hace falta anti-join ni rbind.
  datos_horarios[, fila_ausente := FALSE]

  rejilla <- CJ(
    ESTACION = unique(datos_horarios$ESTACION),
    FECHA    = seq(as.Date(paste0(anio, "-01-01")),
                   as.Date(paste0(anio, "-12-31")), by = "day"),
    HORA     = unique(datos_horarios$HORA)
  )

  datos_horarios <- datos_horarios[rejilla, on = .(ESTACION, FECHA, HORA)]
  datos_horarios[is.na(fila_ausente), fila_ausente := TRUE]

  # Las coordenadas son constantes por estacion y se recuperan de sus otras
  # filas. Sin esto, la agregacion diaria -que agrupa por estacion Y
  # coordenadas- crearia grupos aparte para las filas nuevas, en silencio.
  cols_coord <- c("LONGITUD", "LATITUD", "X_km", "Y_km")
  datos_horarios[, (cols_coord) := lapply(.SD, function(x) x[!is.na(x)][1L]),
                 by = ESTACION, .SDcols = cols_coord]

  cat(sprintf("   Rejilla completada: +%d filas (%d dias-estacion ausentes)\n",
              sum(datos_horarios$fila_ausente),
              sum(datos_horarios$fila_ausente) / uniqueN(datos_horarios$HORA)))

  # 2b. Imputation of NAs in hourly data
  cat("⏳ Imputando NAs horarios...\n")
  datos_horarios <- imputar_na_horario(datos_horarios,
                                       cols_clima = cols_clima,
                                       maxgap     = 3L,
                                       marcar_estado = TRUE,
                                       inventario = inventario)

  cols_estado <- paste0(cols_clima, "_estado")
  cols_estado <- intersect(cols_estado, names(datos_horarios))

  # 2c. Keep only the target year at hourly scale too.
  # The raw folder of a year may include spill-over days of the next year
  # (e.g. 2025 dragged Jan-Mar 2026). The daily/monthly outputs are filtered
  # below, but the hourly output must be filtered here so all three levels are
  # consistent and no foreign year leaks into meteo_madrid_<anio>_horario.rds.
  datos_horarios <- datos_horarios[year(FECHA) == anio]

  # 3. Aggregation to daily scale (threshold 30% NA)
  # Source: datos_horarios (already filtered and with NAs resolved at hourly scale)
  cat("⏳ Agregando a escala diaria...\n")
  by_diario  <- c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "FECHA")
  cols_suma  <- intersect("Precipitaciones", cols_clima)
  cols_media <- setdiff(cols_clima, cols_suma)
  
  # We calculate the daily mean for variables that are averaged (e.g., temperature)
  # and the daily sum for precipitation,
  dt_media <- if (length(cols_media) > 0) {
    datos_horarios[, lapply(.SD, function(x) {
      # For each variable, if the proportion of NAs is greater than or equal to 30%,
      # we return NA for the daily mean. Otherwise, we calculate the mean of the available values, 
      # ignoring NAs.
      if (sum(is.na(x)) / .N >= 0.3) NA_real_ else mean(x, na.rm = TRUE)
    }), by = by_diario, .SDcols = cols_media]
  } else {
    unique(datos_horarios[, ..by_diario])
  }

  # We calculate the daily sum for precipitation, following the same logic as for the mean.
  datos_diarios <- if (length(cols_suma) > 0) {
    dt_suma <- datos_horarios[, lapply(.SD, function(x) {
      # For each variable, if the proportion of NAs is greater than or equal to 30%,
      if (sum(is.na(x)) / .N >= 0.3) NA_real_ else sum(x, na.rm = TRUE)
    }), by = by_diario, .SDcols = cols_suma]
    # We merge the daily mean and daily sum tables by the grouping columns (station, coordinates, date).
    merge(dt_media, dt_suma, by = by_diario)
  } else {
    dt_media
  }
  
  # Estado del dato diario. Se hereda del horario para que el fichero diario
  # —que es el que consumen los análisis de calidad— también distinga medida,
  # imputación, fallo y ausencia de sensor. Un día que mezcla horas medidas e
  # imputadas se etiqueta IMPUTADO: basta una hora inventada para que la media
  # deje de ser una medida pura.
  if (length(cols_estado) > 0) {

    estado_diario <- estado_periodo_clima(datos_horarios, cols_clima,
                                          cols_grupo = by_diario,
                                          umbral_na  = 0.3)

    datos_diarios <- merge(datos_diarios, estado_diario,
                           by = by_diario, all.x = TRUE)

    # Red de seguridad: si la agregación descartó el día por umbral de NA, el
    # valor sale NA y el estado tiene que ser FALLO aunque el recuento de horas
    # no lo hubiera detectado (p. ej. por redondeo en el umbral).
    #
    # SIN_SENSOR y AUSENTE quedan fuera: en ambos el valor es NA por definición
    # y ya están correctamente etiquetados. Sin esta excepción, un día entero
    # ausente se reclasificaría como avería y se perderia la distincion.
    for (v in cols_clima) {
      ce <- paste0(v, "_estado")
      if (!ce %in% names(datos_diarios)) next
      datos_diarios[is.na(get(v)) & !get(ce) %in% c("SIN_SENSOR", "AUSENTE"),
              (ce) := "FALLO"]
      set(datos_diarios, j = ce,
          value = factor(datos_diarios[[ce]], levels = ESTADOS_CLIMA))
    }
  }

  # We filter the daily data to keep only the rows corresponding to the specified year.
  datos_diarios <- datos_diarios[year(FECHA) == anio]
  datos_diarios[, ANO := anio]
  setorder(datos_diarios, FECHA)
  
  cat("  ✅ Diario:", nrow(datos_diarios), "filas |",
      uniqueN(datos_diarios$ESTACION), "estaciones |",
      uniqueN(datos_diarios$FECHA), "días\n")
  
  # 4. Aggregation to monthly scale (threshold 20% NA)
  
  cat("⏳ Agregando a escala mensual...\n")
  # We define the grouping columns for monthly aggregation, which include station and coordinates.
  cols_grupo_m <- c("ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km")
  
  # We calculate the monthly mean for variables that are averaged (e.g., temperature) 
  # and the monthly sum for precipitation.
  mens_media <- if (length(cols_media) > 0) {
    convertir_resolucion(
      dt           = datos_diarios,
      a            = "mensual",
      col_fecha    = "FECHA",
      cols_grupo   = cols_grupo_m,
      cols_valores = cols_media,
      umbral_na    = 0.3,
      fun_agr      = mean
    )
  } else NULL
  # We calculate the monthly sum for precipitation, following the same logic as for the mean.
  datos_mensuales <- if (length(cols_suma) > 0) {
    mens_suma <- convertir_resolucion(
      dt           = datos_diarios,
      a            = "mensual",
      col_fecha    = "FECHA",
      cols_grupo   = cols_grupo_m,
      cols_valores = cols_suma,
      umbral_na    = 0.3,
      fun_agr      = sum
    )
    # We merge the monthly mean and monthly sum tables by the grouping columns (station, coordinates, month).
    if (!is.null(mens_media)) merge(mens_media, mens_suma, by = c(cols_grupo_m, "MES"))
    else mens_suma
  } else {
    mens_media
  }
  
  # Estado del dato mensual. Se calcula desde las HORAS del mes, no desde los
  # estados diarios: encadenar día -> mes reintroduciría el efecto acumulativo.
  if (length(cols_estado) > 0) {

    horario_mes <- copy(datos_horarios)
    horario_mes[, MES := format(FECHA, "%Y-%m")]

    estado_mensual <- estado_periodo_clima(horario_mes, cols_clima,
                                           cols_grupo = c(cols_grupo_m, "MES"),
                                           umbral_na  = 0.3)

    datos_mensuales <- merge(datos_mensuales, estado_mensual,
                             by = c(cols_grupo_m, "MES"), all.x = TRUE)

    for (v in cols_clima) {
      ce <- paste0(v, "_estado")
      if (!ce %in% names(datos_mensuales)) next
      datos_mensuales[is.na(get(v)) & !get(ce) %in% c("SIN_SENSOR", "AUSENTE"),
                (ce) := "FALLO"]
      set(datos_mensuales, j = ce,
          value = factor(datos_mensuales[[ce]], levels = ESTADOS_CLIMA))
    }
  }

  # We filter the monthly data to keep only the rows corresponding to the specified year.
  datos_mensuales <- datos_mensuales[substr(MES, 1, 4) == as.character(anio)]
  datos_mensuales[, ANO := anio]
  setorder(datos_mensuales, MES, ESTACION)
  
  cat("  ✅ Mensual:", nrow(datos_mensuales), "filas |",
      uniqueN(datos_mensuales$ESTACION), "estaciones |",
      uniqueN(datos_mensuales$MES), "meses\n")
  
  # We return a list containing the hourly, daily, and monthly data.tables.
  return(list(horario=datos_horarios,diario = datos_diarios, mensual = datos_mensuales))
}

#===============================================================================
# Cleaning Function for Hourly Traffic Data with Spatial Imputation 
#===============================================================================

limpiar_trafico_espacial_horario <- function(dt, dt_ubicaciones, mapa_distritos, mapa_barrios) {
  
  # Create a working copy of the data to avoid modifying the original by reference 
  # (data.table modifies by reference, so it's important not to touch the original)
  
  dt_work <- copy(dt)
  setDT(dt_work)
  
  # 1. Temporal Formatting: Combine date and hour into a single POSIXct column 
  # for easier handling of hourly data.
  
  dt_work[, fecha_hora := as.POSIXct(fecha, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Madrid")] 
  
  # It passes the date and hour to a single column of type POSIXct, which is more convenient for 
  # handling hourly data. We specify the format of the input date and time, and set the timezone to Madrid.
  
  dt_work[, FECHA := as.Date(fecha_hora)] 
  
  # 2. We only keep these variables 
  cols_medidas <- c("intensidad", "carga")
  
  # 2b. We set to NA the measurements that have error codes "E" (error) or "S" (suspect),
  # as well as any negative values which are not physically meaningful for these variables.
  
  for (col in cols_medidas) {
    dt_work[error %in% c("E", "S") | get(col) < 0, (col) := NA]
  }
  
  # 3. Spatial Coordinates: We need to assign UTM coordinates to each sensor based on
  #the provided location data. This is crucial for the spatial imputation step later on.
  # dt_ubicaciones is the file where we have the coordinates of the sensors
  
  setDT(dt_ubicaciones)
  coor_medidores <- unique(dt_ubicaciones[, .(id, utm_x, utm_y)], by = "id")
  
  # We convert the coordinates to an sf object for spatial operations. We only keep 
  # sensors that have valid UTM coordinates (non-NA).
  
  # We are ignoring sensors without valid UTM coordinates because they cannot be 
  # used for spatial imputation. If we kept them, they would not have a defined
  # location in space, which would cause problems when trying to assign them to districts
  # or neighborhoods, and they would also be excluded from the nearest neighbor imputation step 
  # since they don't have coordinates to calculate distances. By filtering out sensors with NA UTM
  # coordinates at this stage, we ensure that all remaining sensors have valid spatial information 
  # that can be used for the subsequent spatial operations.
  
  sf_sensores <- st_as_sf(coor_medidores[!is.na(utm_x)], 
                          coords = c("utm_x", "utm_y"), 
                          crs = 25830, 
                          remove = FALSE)
  
  cat("Sensores con coordenadas UTM válidas: ", nrow(sf_sensores), "\n")
  cat("Sensores sin coordenadas UTM válidas (excluidos): ", nrow(coor_medidores) - nrow(sf_sensores), "\n")
  
  # 4a. Assignation of districts based on spatial intersection. We transform the district map to 
  # the same CRS as the sensor coordinates to ensure that the spatial join works 
  # correctly. Then we use st_join with st_intersects to assign each sensor to the district it falls into. 
  # This will add a new column "distrito" to sf_sensores_dt with the name of the district for each sensor.
  
  mapa_distritos <- st_transform(mapa_distritos, crs = 25830)
  sf_sensores_dt <- st_join(sf_sensores, mapa_distritos[, "distrito"], join = st_intersects)
  
  # 4b. The same procedure but with neighborhood map.
  mapa_barrios <- st_transform(mapa_barrios, crs = 25830)
  sf_sensores_ba <- st_join(sf_sensores, mapa_barrios[, "barrio"], join = st_intersects)
  
  # Combine the relevant columns from both spatial joins into a single data.table for easier merging 
  # with the hourly data. We keep the sensor ID, UTM coordinates, district and neighborhood information.
  
  dt_sensores_fijos <- as.data.table(sf_sensores_dt)[, .(id, utm_x, utm_y, distrito)]
  dt_barrios_fijos  <- as.data.table(sf_sensores_ba)[, .(id, barrio)]
  dt_sensores_fijos <- merge(dt_sensores_fijos, dt_barrios_fijos, by = "id", all.x = TRUE)
  
  # Cross-setting the spatial information with the hourly data. We merge the spatial information
  # (UTM coordinates, district and neighborhood) into the main hourly data table based on the sensor ID. 
  # This will allow us to perform spatial imputation later on.
  
  dt_work <- merge(dt_work, dt_sensores_fijos[, .(id, utm_x, utm_y, distrito, barrio)],
                   by = "id", all.x = TRUE)
  
  # Discard rows with missing spatial information (sensors without valid UTM coordinates or without district assignment) 
  # since they cannot be used for spatial imputation or aggregation by district/neighborhood.
  # This ensures that the subsequent analysis is based on sensors with complete spatial information.
  
  dt_work <- dt_work[!is.na(utm_x) & !is.na(distrito) & distrito != ""]
  
  # 5. Imputation with nearest neighbor spatial interpolation. 
  # For each hour and each variable, we identify the sensors that have missing values.
  
  # Adding a summary of the number of NAs before imputation for each variable, to have
  # a reference of how many values are being imputed. This can help us understand the extent 
  # of missing data and the impact of the imputation process.
  
  cat("Conteo de valores NA",sum(is.na(dt_work$intensidad)), "en intensidad,",
      sum(is.na(dt_work$carga)), "en carga.\n")
  
  # Tally the total number of NAs across all three variables to determine if we need to perform 
  # imputation. If there are no NAs, we can skip the imputation step entirely.
  
  total_nas <- sum(is.na(dt_work$intensidad)) + sum(is.na(dt_work$carga))
  
  if (total_nas > 0) {
    
    cat("   -> Iniciando imputación por vecinos más cercanos (st_nearest_feature) a nivel horario...\n")
    
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
          
          # Convert the data.tables to sf objects for spatial operations.
          
          sf_falt <- st_as_sf(faltantes, coords = c("utm_x", "utm_y"), crs = 25830)
          sf_val  <- st_as_sf(validos,   coords = c("utm_x", "utm_y"), crs = 25830)
          
          # We use st_nearest_feature to find the index of the nearest valid sensor for each missing sensor.
          nn_indices        <- st_nearest_feature(sf_falt, sf_val)
          valores_imputados <- validos$valor[nn_indices]
          
          ids_faltantes <- faltantes$id
          dt_work[fecha_hora == h & id %in% ids_faltantes,
                  (col) := valores_imputados[match(id, ids_faltantes)]]
        }
      }
    }
    
    cat("   -> ✅ Imputación por vecinos más cercanos completada con éxito.\n")
    
  } else {
    
    cat("   -> ✅ ¡Datos perfectos! No hay valores nulos, se omite la imputación.\n")
  }
  
  # 6. Selection of relevant columns for the output.
  # We keep the original date and time column 
  # (fecha_hora), the date (FECHA), sensor ID, UTM coordinates, district and neighborhood information, 
  # and the traffic variables (intensidad, ocupacion, carga).
  
  cols_output <- c("fecha_hora", "FECHA", "id", "utm_x", "utm_y",
                   "distrito", "barrio", "intensidad", "carga")
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
#   umbral_na    – proporción máxima de NAs permitida antes de devolver NA (default 0.3)
#   fun_agr      – función de agregación aplicada a cada columna (default mean)
#
# Devuelve un data.table con una fila por grupo × período.

convertir_resolucion <- function(dt,
                                 a            = c("diario", "mensual"),
                                 col_fecha    = "FECHA",
                                 cols_grupo   = "ESTACION",
                                 cols_valores = NULL,
                                 umbral_na    = 0.3,
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


#===============================================================================
# Function: Aggregation to Daily Traffic Data by District and Neighborhood
#===============================================================================

# The input is the hourly traffic data (output of limpiar_trafico_espacial_horario), 
# and it calculates the daily average of intensity, occupancy, and load for each district and neighborhood.
# It also counts the number of unique sensors contributing to each daily average.

agregar_trafico_diario <- function(dt_horario) {
  
  dt <- as.data.table(dt_horario)
  
  # Traffic data for each district
  # (we include the district in the by to maintain the hierarchy neighborhood → district)
  
  dt_diario_distrito <- dt[, .(
    intensidad    = mean(intensidad, na.rm = TRUE),
    carga         = mean(carga,      na.rm = TRUE),
    num_medidores = uniqueN(id)
  ), by = .(distrito, FECHA)]
  setorder(dt_diario_distrito, FECHA, distrito)
  
  # Aggregation by neighborhood (barrio) within each district
  
  dt_diario_barrio <- dt[!is.na(barrio) & barrio != "", .(
    intensidad    = mean(intensidad, na.rm = TRUE),
    carga         = mean(carga,      na.rm = TRUE),
    num_medidores = uniqueN(id)
  ), by = .(distrito, barrio, FECHA)]
  setorder(dt_diario_barrio, FECHA, distrito, barrio)
  
  # Return a list with both daily data.tables for district and neighborhood
  
  return(list(distrito = dt_diario_distrito, barrio = dt_diario_barrio))
}


#===============================================================================
# Function to aggregate daily traffic data to monthly averages by district and neighborhood
#===============================================================================

agregar_trafico_mensual <- function(dt_diario_distrito,
                                    dt_diario_barrio  = NULL, # This argument is optional. 
                                    umbral_na         = 0.3) {
  
  vars_trafico <- c("intensidad", "carga")
  
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


#===============================================================================
# Function: Aggregation of Hourly Traffic Data by District and Neighborhood
#===============================================================================

# The input is the cleaned hourly traffic data (output of limpiar_trafico_espacial_horario).
# It extracts the hour from fecha_hora and calculates the mean of intensity, occupancy,
# and load for each district and neighborhood × date × hour combination.
# It also counts the number of unique sensors contributing to each average.

#-------------------------------------------------------------------------------
# Rejilla horaria completa para los agregados de trafico por zona
#-------------------------------------------------------------------------------
# agregar_trafico_horario_zona() agrupa con by = .(zona, FECHA, HORA), y un "by"
# SOLO crea grupos para las combinaciones presentes en los datos. Si ningun
# sensor de un barrio reporta a las 03:00, esa hora no genera fila: no queda
# como NA, sencillamente no existe, y desaparece de cualquier denominador.
#
# En 2019 esto afectaba a 852 horas del agregado por barrio y 64 del de
# distrito, concentradas en la madrugada, cuando el trafico es minimo y los
# sensores pueden no transmitir.
#
# Las filas anadidas llevan num_medidores = 0, que es el dato honesto (cero
# sensores reportaron) y ademas sirve de indicador de fiabilidad por celda: hoy
# una media apoyada en 28 sensores y otra en 1 entran al modelo con el mismo
# peso, y con esta columna se pueden distinguir.
#
# Argumentos:
#   dt        - agregado por zona (salida de agregar_trafico_horario_zona)
#   col_zona  - "distrito" o "barrio"
#   dias      - vector de fechas que debe cubrir la rejilla
#
# Devuelve el mismo data.table con las filas que faltaban y una columna ESTADO
# de dos niveles: OK / AUSENTE. No hay FALLO porque la agregacion usa
# mean(na.rm = TRUE) y nunca produce NA en las filas que si existen.
completar_rejilla_trafico <- function(dt, col_zona, dias) {

  dt <- as.data.table(dt)

  rejilla <- CJ(zona = unique(dt[[col_zona]]), FECHA = dias, HORA = 1:24)
  setnames(rejilla, "zona", col_zona)

  dt[, fila_ausente := FALSE]
  dt <- dt[rejilla, on = c(col_zona, "FECHA", "HORA")]
  dt[is.na(fila_ausente), `:=`(fila_ausente = TRUE, num_medidores = 0L)]

  # El barrio arrastra su distrito: en las filas nuevas hay que heredarlo para
  # no romper la jerarquia barrio -> distrito.
  if (col_zona == "barrio" && "distrito" %in% names(dt)) {
    jer <- unique(dt[!is.na(distrito), .(barrio, distrito)])
    dt[, distrito := jer$distrito[match(barrio, jer$barrio)]]
  }

  dt[, ESTADO := factor(fifelse(fila_ausente, "AUSENTE", "OK"),
                        levels = c("OK", "AUSENTE"))]
  dt[]
}


agregar_trafico_horario_zona <- function(dt_horario) {

  dt <- as.data.table(dt_horario)

  # Extract numeric hour (1–24) from the datetime column
  dt[, HORA := hour(fecha_hora) + 1L]

  # Aggregation by district × date × hour
  dt_hora_distrito <- dt[, .(
    intensidad    = mean(intensidad, na.rm = TRUE),
    carga         = mean(carga,      na.rm = TRUE),
    num_medidores = uniqueN(id)
  ), by = .(distrito, FECHA, HORA)]
  setorder(dt_hora_distrito, FECHA, HORA, distrito)

  # Aggregation by neighborhood (barrio) within each district × date × hour
  dt_hora_barrio <- dt[!is.na(barrio) & barrio != "", .(
    intensidad    = mean(intensidad, na.rm = TRUE),
    carga         = mean(carga,      na.rm = TRUE),
    num_medidores = uniqueN(id)
  ), by = .(distrito, barrio, FECHA, HORA)]
  setorder(dt_hora_barrio, FECHA, HORA, distrito, barrio)

  return(list(distrito = dt_hora_distrito, barrio = dt_hora_barrio))
}


#----------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
library(stringi) # Para quitar tildes fácilmente

# Función auxiliar para limpiar textos (quita tildes, minúsculas, espacios extra)
limpiar_nombres <- function(texto) {
  texto <- tolower(trimws(texto))
  texto <- stri_trans_general(texto, "Latin-ASCII") # Quita tildes
  # Reemplaza múltiples espacios o guiones raros por un solo espacio
  texto <- gsub("[[:punct:]]+", " ", texto) 
  texto <- gsub("\\s+", " ", texto)
  return(trimws(texto))
}




