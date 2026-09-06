# ==============================================================================
# Cleaning functions for the official daily meteorological files
# ==============================================================================
#
# These functions are deliberately independent from the hourly climate pipeline.
# The D01...D31 values supplied by the data provider are treated as final daily
# observations: no hourly-to-daily or daily-to-monthly aggregation is performed.

library(data.table)

ESTADOS_CLIMA_DIARIO <- c("OK", "IMPUTADO", "FALLO", "AUSENTE", "SIN_SENSOR")

VARIABLES_CLIMA_DIARIAS <- unname(nombres_magnitudes_clima[
  c("81", "83", "86", "87", "88", "89")
])


leer_meteo_diario_anual <- function(anio, carpeta_diarios) {

  patron <- paste0("^", anio, "_meteorologicos-diarios-csv[.]csv$")
  archivos <- list.files(
    carpeta_diarios,
    pattern = patron,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(archivos) != 1L) {
    stop(
      "Se esperaba un unico CSV diario para ", anio, " y se encontraron ",
      length(archivos), ": ", paste(basename(archivos), collapse = ", ")
    )
  }

  cat("   Leyendo ", basename(archivos), "...\n", sep = "")
  fread(archivos, sep = ";", encoding = "Latin-1")
}


coordenadas_estaciones_diarias <- function(dt_ubicaciones) {

  requeridas <- c(
    "CODIGO_CORTO", "LONGITUD", "LATITUD",
    "COORDENADA_X_ETRS89", "COORDENADA_Y_ETRS89"
  )
  ausentes <- setdiff(requeridas, names(dt_ubicaciones))
  if (length(ausentes)) {
    stop(
      "Faltan columnas en el catalogo de estaciones: ",
      paste(ausentes, collapse = ", ")
    )
  }

  ubicaciones <- as.data.table(copy(dt_ubicaciones))[, .(
    CODIGO_CORTO = as.integer(CODIGO_CORTO),
    LONGITUD = suppressWarnings(
      as.numeric(gsub("\\.", "", as.character(LONGITUD))) / 1e7
    ),
    LATITUD = suppressWarnings(
      as.numeric(gsub("\\.", "", as.character(LATITUD))) / 1e7
    ),
    X_km = suppressWarnings(
      as.numeric(gsub(",", ".", as.character(COORDENADA_X_ETRS89))) / 1000
    ),
    Y_km = suppressWarnings(
      as.numeric(gsub(",", ".", as.character(COORDENADA_Y_ETRS89))) / 1000
    )
  )]

  ubicaciones[, ESTACION := unname(
    nombres_estaciones_clima[as.character(CODIGO_CORTO)]
  )]
  ubicaciones <- ubicaciones[!is.na(ESTACION)]

  duplicadas <- ubicaciones[duplicated(ESTACION), unique(ESTACION)]
  if (length(duplicadas)) {
    stop(
      "Hay estaciones duplicadas en el catalogo: ",
      paste(duplicadas, collapse = ", ")
    )
  }

  ubicaciones[, .(ESTACION, LONGITUD, LATITUD, X_km, Y_km)]
}


reestructurar_meteo_diario <- function(dt_bruto, anio) {

  dt <- as.data.table(copy(dt_bruto))
  requeridas <- c(
    "ESTACION", "MAGNITUD", "ANO", "MES",
    paste0("D", sprintf("%02d", 1:31)),
    paste0("V", sprintf("%02d", 1:31))
  )
  ausentes <- setdiff(requeridas, names(dt))
  if (length(ausentes)) {
    stop(
      "Faltan columnas en el CSV diario: ",
      paste(ausentes, collapse = ", ")
    )
  }

  anios_archivo <- sort(unique(dt$ANO))
  if (!anio %in% anios_archivo) {
    stop("El CSV no contiene observaciones del anio ", anio, ".")
  }
  if (length(anios_archivo) > 1L) {
    message(
      "   El CSV tambien contiene los anios ",
      paste(setdiff(anios_archivo, anio), collapse = ", "),
      "; se conservara solo ", anio, "."
    )
  }

  dt <- dt[ANO == anio]
  dt[, `:=`(
    ESTACION = as.integer(ESTACION),
    MAGNITUD = as.integer(MAGNITUD),
    MES = as.integer(MES)
  )]

  codigos_variables <- c("81", "83", "86", "87", "88", "89")
  codigos_estaciones <- as.integer(names(nombres_estaciones_clima))
  dt <- dt[
    MAGNITUD %in% as.integer(codigos_variables) &
      ESTACION %in% codigos_estaciones
  ]

  if (!nrow(dt)) {
    stop("No quedan observaciones tras filtrar estaciones y magnitudes.")
  }

  claves <- c("ESTACION", "MAGNITUD", "ANO", "MES")
  duplicadas <- dt[duplicated(dt, by = claves)]
  if (nrow(duplicadas)) {
    stop(
      "Hay ", nrow(duplicadas),
      " filas duplicadas por estacion, magnitud, anio y mes."
    )
  }

  dt[, id_fila := .I]
  columnas_dato <- paste0("D", sprintf("%02d", 1:31))
  columnas_validez <- paste0("V", sprintf("%02d", 1:31))

  valores <- melt(
    dt,
    id.vars = c("id_fila", "ESTACION", "MAGNITUD", "ANO", "MES"),
    measure.vars = columnas_dato,
    variable.name = "DIA_COLUMNA",
    value.name = "VALOR_FUENTE"
  )
  validez <- melt(
    dt,
    id.vars = "id_fila",
    measure.vars = columnas_validez,
    variable.name = "DIA_COLUMNA",
    value.name = "VALIDACION"
  )

  valores[, DIA := as.integer(sub("D", "", DIA_COLUMNA))]
  validez[, DIA := as.integer(sub("V", "", DIA_COLUMNA))]
  valores[, VALIDACION := as.character(validez$VALIDACION)]
  valores[, DIA_COLUMNA := NULL]

  codigos_desconocidos <- setdiff(
    unique(na.omit(valores$VALIDACION)), c("V", "N", "")
  )
  if (length(codigos_desconocidos)) {
    stop(
      "Codigos de validacion no reconocidos: ",
      paste(codigos_desconocidos, collapse = ", ")
    )
  }

  valores[, FECHA := suppressWarnings(as.IDate(sprintf(
    "%04d-%02d-%02d", ANO, MES, DIA
  )))]

  # D29-D31 exist in every row, including months where those dates do not
  # exist. They are placeholders, not missing observations.
  valores <- valores[!is.na(FECHA)]
  valores[, VALOR_FUENTE := suppressWarnings(as.numeric(VALOR_FUENTE))]
  valores[, VALOR := fifelse(VALIDACION == "V", VALOR_FUENTE, NA_real_)]

  valores[, `:=`(
    ESTACION = unname(nombres_estaciones_clima[as.character(ESTACION)]),
    MAGNITUD = unname(nombres_magnitudes_clima[as.character(MAGNITUD)])
  )]

  if (anyNA(valores$ESTACION) || anyNA(valores$MAGNITUD)) {
    stop("No se pudieron traducir todos los codigos de estacion o magnitud.")
  }

  dias_presentes <- unique(valores[, .(ESTACION, FECHA)])
  dias_presentes[, presente_fuente := TRUE]

  datos <- dcast(
    valores,
    ESTACION + FECHA ~ MAGNITUD,
    value.var = "VALOR"
  )
  datos <- merge(
    datos,
    dias_presentes,
    by = c("ESTACION", "FECHA"),
    all.x = TRUE,
    sort = FALSE
  )

  for (variable in VARIABLES_CLIMA_DIARIAS) {
    if (!variable %in% names(datos)) datos[, (variable) := NA_real_]
  }

  list(
    datos = datos,
    estaciones = sort(unique(valores$ESTACION))
  )
}


completar_rejilla_meteo_diaria <- function(datos, estaciones, anio,
                                            coordenadas) {

  calendario <- seq(
    as.IDate(paste0(anio, "-01-01")),
    as.IDate(paste0(anio, "-12-31")),
    by = "day"
  )
  rejilla <- CJ(ESTACION = estaciones, FECHA = calendario)

  dt <- datos[rejilla, on = .(ESTACION, FECHA)]
  dt[, fila_ausente := is.na(presente_fuente)]
  dt[, presente_fuente := NULL]

  dt <- merge(
    dt,
    coordenadas,
    by = "ESTACION",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(dt, ESTACION, FECHA)

  if (dt[, anyNA(LONGITUD) || anyNA(LATITUD) || anyNA(X_km) || anyNA(Y_km)]) {
    sin_coordenadas <- dt[
      is.na(LONGITUD) | is.na(LATITUD) | is.na(X_km) | is.na(Y_km),
      unique(ESTACION)
    ]
    stop(
      "Faltan coordenadas para: ",
      paste(sin_coordenadas, collapse = ", ")
    )
  }

  dt[]
}


etiquetar_e_imputar_meteo_diario <- function(dt, maxgap_dias = 3L,
                                              imputar_fallos = TRUE) {

  datos <- as.data.table(copy(dt))
  setorder(datos, ESTACION, FECHA)

  if (isTRUE(imputar_fallos) && !requireNamespace("zoo", quietly = TRUE)) {
    stop(
      "El paquete 'zoo' es necesario para la interpolacion lineal diaria."
    )
  }
  if (length(maxgap_dias) != 1L || is.na(maxgap_dias) || maxgap_dias < 0) {
    stop("maxgap_dias debe ser un entero no negativo.")
  }
  maxgap_dias <- as.integer(maxgap_dias)

  for (variable in VARIABLES_CLIMA_DIARIAS) {
    estado <- paste0(variable, "_estado")

    estaciones_sin_sensor <- datos[, .(
      sin_sensor = all(is.na(get(variable)))
    ), by = ESTACION][sin_sensor == TRUE, ESTACION]

    datos[, (estado) := fcase(
      ESTACION %in% estaciones_sin_sensor, "SIN_SENSOR",
      fila_ausente, "AUSENTE",
      !is.na(get(variable)), "OK",
      default = "FALLO"
    )]

    if (isTRUE(imputar_fallos) && maxgap_dias > 0L) {
      datos[, .valor_interpolado := {
        x <- get(variable)
        if (sum(!is.na(x)) < 2L) x else zoo::na.approx(
          x, na.rm = FALSE, maxgap = maxgap_dias
        )
      }, by = ESTACION]

      imputables <- which(
        datos[[estado]] == "FALLO" & !is.na(datos$.valor_interpolado)
      )
      if (length(imputables)) {
        set(
          datos,
          i = imputables,
          j = variable,
          value = datos$.valor_interpolado[imputables]
        )
        set(datos, i = imputables, j = estado, value = "IMPUTADO")
      }
      datos[, .valor_interpolado := NULL]
    }

    set(
      datos,
      j = estado,
      value = factor(datos[[estado]], levels = ESTADOS_CLIMA_DIARIO)
    )
  }

  # A day absent from the source must remain absent even when it is surrounded
  # by observations that would permit a numerical interpolation.
  datos[
    fila_ausente == TRUE,
    (VARIABLES_CLIMA_DIARIAS) := lapply(.SD, function(x) NA_real_),
    .SDcols = VARIABLES_CLIMA_DIARIAS
  ]

  datos[, ANO := as.integer(format(FECHA, "%Y"))]
  columnas_estado <- paste0(VARIABLES_CLIMA_DIARIAS, "_estado")
  setcolorder(datos, c(
    "ESTACION", "LONGITUD", "LATITUD", "X_km", "Y_km", "FECHA",
    VARIABLES_CLIMA_DIARIAS, columnas_estado, "fila_ausente", "ANO"
  ))
  setorder(datos, FECHA, ESTACION)
  datos[]
}


validar_meteo_diario_procesado <- function(dt, anio) {

  datos <- as.data.table(dt)
  calendario <- seq(
    as.IDate(paste0(anio, "-01-01")),
    as.IDate(paste0(anio, "-12-31")),
    by = "day"
  )
  n_estaciones <- uniqueN(datos$ESTACION)
  esperado <- n_estaciones * length(calendario)

  if (nrow(datos) != esperado) {
    stop("Filas obtenidas: ", nrow(datos), "; esperadas: ", esperado, ".")
  }
  if (datos[, anyDuplicated(paste(ESTACION, FECHA)) > 0L]) {
    stop("Hay claves estacion-fecha duplicadas.")
  }
  if (!setequal(unique(datos$FECHA), calendario)) {
    stop("El calendario diario no esta completo.")
  }
  if (any(datos$ANO != anio)) {
    stop("La salida contiene fechas ajenas al anio ", anio, ".")
  }

  if (datos[
    fila_ausente == TRUE,
    any(!is.na(unlist(.SD))),
    .SDcols = VARIABLES_CLIMA_DIARIAS
  ]) {
    stop("Alguna fila ausente contiene valores climaticos no NA.")
  }

  for (variable in VARIABLES_CLIMA_DIARIAS) {
    estado <- paste0(variable, "_estado")
    estados <- as.character(datos[[estado]])
    valor <- datos[[variable]]

    if (anyNA(estados) || !all(estados %in% ESTADOS_CLIMA_DIARIO)) {
      stop("Estados desconocidos o NA en ", estado, ".")
    }
    if (any(estados %in% c("OK", "IMPUTADO") & is.na(valor))) {
      stop("Estado con valor NA incoherente en ", variable, ".")
    }
    if (any(estados %in% c("FALLO", "AUSENTE", "SIN_SENSOR") & !is.na(valor))) {
      stop("Estado sin dato con valor numerico incoherente en ", variable, ".")
    }
  }

  invisible(TRUE)
}


resumen_calidad_meteo_diario <- function(dt, anio) {

  datos <- as.data.table(dt)
  estados <- paste0(VARIABLES_CLIMA_DIARIAS, "_estado")
  largo <- melt(
    datos,
    id.vars = c("ESTACION", "FECHA"),
    measure.vars = estados,
    variable.name = "VARIABLE",
    value.name = "ESTADO"
  )
  largo[, VARIABLE := sub("_estado$", "", VARIABLE)]

  resumen <- dcast(
    largo[, .N, by = .(VARIABLE, ESTADO)],
    VARIABLE ~ ESTADO,
    value.var = "N",
    fill = 0
  )
  for (estado in ESTADOS_CLIMA_DIARIO) {
    if (!estado %in% names(resumen)) resumen[, (estado) := 0L]
  }
  setcolorder(resumen, c("VARIABLE", ESTADOS_CLIMA_DIARIO))

  cat(
    "   ", anio, ": ", nrow(datos), " filas | ",
    uniqueN(datos$ESTACION), " estaciones | ",
    uniqueN(datos$FECHA), " dias | ",
    datos[, sum(fila_ausente)], " dias-estacion ausentes\n",
    sep = ""
  )
  print(resumen)
  invisible(resumen)
}


procesar_anio_meteo_diario_fuente <- function(anio, carpeta_diarios,
                                               ruta_estaciones,
                                               maxgap_dias = 3L,
                                               imputar_fallos = TRUE) {

  ubicaciones_raw <- as.data.table(read.csv(ruta_estaciones, sep = ";"))
  if (nrow(ubicaciones_raw) < 20L) {
    stop(
      "El catalogo de estaciones solo contiene ", nrow(ubicaciones_raw),
      " filas; revisa la codificacion de ", basename(ruta_estaciones), "."
    )
  }
  coordenadas <- coordenadas_estaciones_diarias(ubicaciones_raw)

  bruto <- leer_meteo_diario_anual(anio, carpeta_diarios)
  reestructurado <- reestructurar_meteo_diario(bruto, anio)
  completo <- completar_rejilla_meteo_diaria(
    reestructurado$datos,
    reestructurado$estaciones,
    anio,
    coordenadas
  )
  resultado <- etiquetar_e_imputar_meteo_diario(
    completo,
    maxgap_dias = maxgap_dias,
    imputar_fallos = imputar_fallos
  )

  validar_meteo_diario_procesado(resultado, anio)
  resumen_calidad_meteo_diario(resultado, anio)
  resultado
}
