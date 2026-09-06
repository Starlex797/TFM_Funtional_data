# Preparacion de datos para la comparacion manual de interpoladores.
# No ejecuta analisis al hacer source(). Dependencias: data.table, sf, ggplot2.

variables_comparacion_clima <- c(
  'Temperatura', 'Humedad_Relativa', 'Precipitaciones',
  'Presion_Barometrica', 'Radiacion_Solar', 'Velocidad_Viento'
)

archivo_clima_5 <- function(raiz, escala, anio = 2025L) {
  if (!escala %in% c('horario', 'diario', 'mensual')) stop('Escala desconocida.')
  carpeta <- file.path(raiz, 'data', 'processed', 'Clima', escala)
  # Los archivos actuales terminan en diario5, no diario_5. Ambos son aceptados,
  # pero nunca se escoge silenciosamente entre dos versiones existentes.
  candidatos <- file.path(carpeta, sprintf('meteo_madrid_%d_%s%s.rds',
                                           anio, escala, c('5', '_5')))
  encontrados <- candidatos[file.exists(candidatos)]
  if (length(encontrados) != 1L) stop('Se esperaba un unico archivo: ',
                                     paste(candidatos, collapse = ' / '))
  encontrados
}

normalizar_entrada_comparacion <- function(dt, escala, anio) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  requeridas <- c('ESTACION', 'LONGITUD', 'LATITUD', variables_comparacion_clima,
                  paste0(variables_comparacion_clima, '_estado'))
  if (!all(requeridas %in% names(dt))) stop('Faltan columnas: ',
      paste(setdiff(requeridas, names(dt)), collapse = ', '))
  dt[, ESTACION := enc2utf8(as.character(ESTACION))]
  if (escala == 'mensual') {
    if (!'MES' %in% names(dt)) stop('El archivo mensual debe contener MES.')
    dt[, FECHA := as.Date(paste0(as.character(MES), '-01'))]
  } else {
    dt[, FECHA := as.Date(FECHA)]
  }
  if (anyNA(dt$FECHA)) stop('Fechas invalidas.')
  dt <- dt[format(FECHA, '%Y') == as.character(anio)]
  if (!nrow(dt)) stop('No hay datos del anio solicitado.')
  if (escala == 'horario') {
    dt[, HORA := as.integer(sub('^H', '', as.character(HORA)))]
    if (anyNA(dt$HORA) || any(!dt$HORA %in% 1:24)) stop('Horas fuera de 1:24.')
  } else dt[, HORA := 0L]
  dt[, MES := format(FECHA, '%Y-%m')]
  if (anyDuplicated(dt[, .(ESTACION, FECHA, HORA)])) stop('Estacion-tiempo duplicados.')
  dt[]
}

coordenadas_comparacion <- function(dt) {
  coords <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
  if (anyDuplicated(coords$ESTACION)) stop('Coordenadas no constantes por estacion.')
  coords[, `:=`(LONGITUD_ORIGINAL = LONGITUD, LATITUD_ORIGINAL = LATITUD)]
  # Reparacion decimal ya utilizada por los scripts del proyecto. Se conserva
  # cada valor original y se comprueba el resultado contra X_km/Y_km si existen.
  i <- which(is.finite(coords$LATITUD) & coords$LATITUD > 0 & coords$LATITUD < 35)
  coords[i, LATITUD := LATITUD * 10^round(log10(40.4 / LATITUD))]
  i <- which(is.finite(coords$LONGITUD) & coords$LONGITUD < 0 & coords$LONGITUD > -1)
  coords[i, LONGITUD := LONGITUD * 10]
  if (any(!is.finite(coords$LATITUD) | !is.finite(coords$LONGITUD) |
          coords$LATITUD < 40.2 | coords$LATITUD > 40.7 |
          coords$LONGITUD < -4 | coords$LONGITUD > -3.3)) stop('Coordenadas fuera de Madrid.')
  xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(coords,
          coords = c('LONGITUD', 'LATITUD'), crs = 4326), 25830))
  coords[, `:=`(X = xy[, 1], Y = xy[, 2])]
  if (all(c('X_km', 'Y_km') %in% names(dt))) {
    originales <- unique(dt[, .(ESTACION, X_km, Y_km)])
    if (anyDuplicated(originales$ESTACION)) stop('Coordenadas proyectadas inconsistentes.')
    coords <- merge(coords, originales, by = 'ESTACION', sort = FALSE)
    coords[, Diferencia_UTM_m := sqrt((X - 1000 * X_km)^2 + (Y - 1000 * Y_km)^2)]
    if (any(!is.finite(coords$Diferencia_UTM_m) | coords$Diferencia_UTM_m > 20)) {
      stop('Las coordenadas corregidas no coinciden con X_km/Y_km (20 m).')
    }
  }
  # Un pliegue por ubicacion: no se valida usando otro sensor colocalizado.
  coords[, Ubicacion := .GRP, by = .(LONGITUD, LATITUD)]
  coords[]
}

preparar_variable_comparacion <- function(dt, horario, coords, variable, escala) {
  d <- dt[, .(ESTACION, FECHA, HORA, MES, Valor = get(variable),
              Estado = as.character(get(paste0(variable, '_estado'))))]
  estados <- c('OK', 'IMPUTADO', 'FALLO', 'AUSENTE', 'SIN_SENSOR')
  if (anyNA(d$Estado) || any(!d$Estado %in% estados)) stop('Estado desconocido en ', variable)
  d[, Hora_imputada := FALSE]
  if (escala != 'horario') {
    llave <- if (escala == 'diario') 'FECHA' else 'MES'
    flags <- horario[, .(Hora_imputada = any(
      as.character(get(paste0(variable, '_estado'))) == 'IMPUTADO')),
      by = c('ESTACION', llave)]
    d[, Hora_imputada := NULL]
    d <- merge(d, flags, by = c('ESTACION', llave), all.x = TRUE, sort = FALSE)
    if (anyNA(d$Hora_imputada)) stop('Falta trazabilidad horaria para ', escala)
  }
  flag <- paste0('imp_', variable)
  if (flag %in% names(dt)) {
    # Compatibilidad adicional con marcas antiguas; no sustituye a *_estado.
    ff <- dt[, .(ESTACION, FECHA, HORA, Imputado_legacy = get(flag))]
    d <- merge(d, ff, by = c('ESTACION', 'FECHA', 'HORA'), sort = FALSE)
  } else d[, Imputado_legacy := FALSE]
  d[, Valido := is.finite(Valor) & Estado == 'OK' & !Hora_imputada &
      !is.na(Imputado_legacy) & !Imputado_legacy]
  audit <- d[, .(N_filas = .N, N_finitos_archivo = sum(is.finite(Valor)),
      N_OK_archivo = sum(Estado == 'OK' & is.finite(Valor)),
      N_OK_con_hora_imputada = sum(Estado == 'OK' & is.finite(Valor) & Hora_imputada),
      N_validos = sum(Valido)), by = ESTACION]
  d <- merge(d[Valido == TRUE], coords[, .(ESTACION, Ubicacion)], by = 'ESTACION')
  if (!nrow(d)) stop('Sin observaciones reales para ', variable, '/', escala)
  # Si hay dos sensores validos colocalizados se usa su media como observacion
  # de esa ubicacion. Ambos quedan excluidos juntos durante la validacion.
  d <- d[, .(Valor = mean(Valor)), by = .(Ubicacion, FECHA, HORA, MES)]
  sitios <- unique(coords[Ubicacion %in% d$Ubicacion, .(Ubicacion, X, Y)])
  data.table::setorder(sitios, Ubicacion)
  tiempos <- unique(dt[, .(FECHA, HORA, MES)])
  data.table::setorder(tiempos, FECHA, HORA)
  tiempos[, Tiempo := .I]
  d <- merge(d, tiempos, by = c('FECHA', 'HORA', 'MES'))
  Y <- matrix(NA_real_, nrow(tiempos), nrow(sitios))
  Y[cbind(d$Tiempo, match(d$Ubicacion, sitios$Ubicacion))] <- d$Valor
  D <- as.matrix(dist(as.matrix(sitios[, .(X, Y)])))
  diag(D) <- Inf
  list(Y = Y, D = D, sitios = sitios, tiempos = tiempos, auditoria = audit)
}

