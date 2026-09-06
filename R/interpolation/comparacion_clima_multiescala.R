# Sensibilidad de RMSE a k. No selecciona ningun k ni metodo automaticamente.
# Cada ubicacion se predice usando las restantes en el mismo instante.

# 1. Definir una muestra comun para que todos los k comparen los mismos casos.
muestra_comun_rmse <- function(Y, variable, cobertura_minima = 0.90,
                              umbral_lluvia = 0.1) {
  donantes <- rowSums(is.finite(Y)) - 1L
  objetivos <- is.finite(Y) & (donantes >= 1L)
  # Evaluacion CONDICIONADA a lluvia observada en el objetivo.
  # Los vecinos secos conservan su valor cero.
  if (variable == 'Precipitaciones') objetivos <- objetivos & (Y > umbral_lluvia)
  objetivos[is.na(objetivos)] <- FALSE
  if (!sum(objetivos)) stop('No hay objetivos evaluables para ', variable)

  cobertura <- rbindlist(lapply(seq_len(ncol(Y) - 1L), function(k) {
    incluidos <- objetivos & (donantes >= k)
    base_estacion <- colSums(objetivos)
    fraccion <- colSums(incluidos)[base_estacion > 0] / base_estacion[base_estacion > 0]
    data.table(k = k, N = sum(incluidos), N_base = sum(objetivos),
               Cobertura = sum(incluidos) / sum(objetivos),
               Cobertura_min_ubicacion = min(fraccion))
  }))
  # Esto limita el rango REPRESENTABLE con los datos; no elige el mejor k.
  k_max <- max(cobertura[Cobertura >= cobertura_minima &
                         Cobertura_min_ubicacion >= cobertura_minima, k])
  mascara <- objetivos & (donantes >= k_max)
  list(mascara = mascara, k_max = k_max, cobertura = cobertura)
}

# 2. Predecir la ubicacion retirada para k=1,...,k_max.
# Se ordenan las estaciones por distancia. En cada instante se cuentan
# solamente vecinos con observacion valida.
predicciones_rmse <- function(Y, D, objetivo, k_max) {
  candidatos <- setdiff(order(D[objetivo, ]), objetivo)
  if (any(!is.finite(D[objetivo, candidatos]) | D[objetivo, candidatos] <= 0)) {
    stop('Agrupar primero los sensores colocalizados.')
  }
  z <- Y[, candidatos, drop = FALSE]
  n <- nrow(Y)
  pred <- list(
    Media = rowMeans(z, na.rm = TRUE),
    NN = rep(NA_real_, n),
    KNN = matrix(NA_real_, n, k_max),
    IDW1 = matrix(NA_real_, n, k_max),
    IDW2 = matrix(NA_real_, n, k_max)
  )
  cuenta <- suma <- suma1 <- suma2 <- peso1 <- peso2 <- numeric(n)
  for (j in seq_along(candidatos)) {
    ok <- is.finite(z[, j])
    distancia <- D[objetivo, candidatos[j]]
    cuenta <- cuenta + ok
    suma[ok] <- suma[ok] + z[ok, j]
    suma1[ok] <- suma1[ok] + z[ok, j] / distancia
    suma2[ok] <- suma2[ok] + z[ok, j] / distancia^2
    peso1[ok] <- peso1[ok] + 1 / distancia
    peso2[ok] <- peso2[ok] + 1 / distancia^2
    for (k in seq_len(min(j, k_max))) {
      filas <- which(ok & cuenta == k)
      pred$KNN[filas, k] <- suma[filas] / k
      pred$IDW1[filas, k] <- suma1[filas] / peso1[filas]
      pred$IDW2[filas, k] <- suma2[filas] / peso2[filas]
      if (k == 1L) pred$NN[filas] <- z[filas, j]
    }
  }
  pred
}

# 3. Calcular RMSE para cada metodo y k.
# Primero MSE por ubicacion; despues raiz de su media.
# Asi cada ubicacion tiene el mismo peso aunque le falten algunos periodos.
calcular_curva_rmse <- function(red, variable, cobertura_minima = 0.90,
                                umbral_lluvia = 0.1) {
  Y <- red$Y
  muestra <- muestra_comun_rmse(Y, variable, cobertura_minima, umbral_lluvia)
  errores <- list()
  i <- 0L
  for (s in seq_len(ncol(Y))) {
    filas <- which(muestra$mascara[, s])
    if (!length(filas)) next
    observado <- Y[filas, s]
    pred <- predicciones_rmse(Y, red$D, s, muestra$k_max)
    for (k in seq_len(muestra$k_max)) {
      estimado <- cbind(
        Media = pred$Media[filas], NN = pred$NN[filas],
        KNN = pred$KNN[filas, k], IDW1 = pred$IDW1[filas, k],
        IDW2 = pred$IDW2[filas, k]
      )
      stopifnot(all(is.finite(estimado)))
      i <- i + 1L
      errores[[i]] <- data.table(
        Ubicacion = red$sitios$Ubicacion[s], k = k,
        Metodo = c('Media', '1-NN', 'KNN', 'IDW p=1', 'IDW p=2'),
        N = length(filas), MSE = colMeans((estimado - observado)^2)
      )
    }
    # Referencia cero: permite ver si, en lluvia, se mejora predecir siempre 0.
    if (variable == 'Precipitaciones') {
      for (k in seq_len(muestra$k_max)) {
        i <- i + 1L
        errores[[i]] <- data.table(Ubicacion = red$sitios$Ubicacion[s],
          k = k, Metodo = 'Siempre cero', N = length(filas), MSE = mean(observado^2))
      }
    }
  }
  por_ubicacion <- rbindlist(errores)
  curva <- por_ubicacion[, .(RMSE = sqrt(mean(MSE)), N = sum(N),
                             N_ubicaciones = .N), by = .(Metodo, k)]
  curva[, Variable := variable]
  curva[, Muestra := if (variable == 'Precipitaciones')
    paste0('Objetivo con lluvia > ', umbral_lluvia, ' mm') else 'Todos los periodos validos']
  list(curva = curva, errores = por_ubicacion, cobertura = muestra$cobertura)
}

# 4. Una figura por variable: RMSE en Y, numero entero de vecinos en X.
# Las lineas horizontales son referencias, no recomendaciones.
grafico_rmse_k <- function(curva, escala, unidad, k_elegido = NA_integer_) {
  variable <- unique(curva$Variable)
  locales <- curva[Metodo %in% c('KNN', 'IDW p=1', 'IDW p=2')]
  referencias <- unique(curva[Metodo %in% c('Media', '1-NN', 'Siempre cero'),
                              .(Metodo, RMSE)])
  k_max <- max(curva$k)
  subtitulo <- paste0(unique(curva$Muestra), ' | ', unique(curva$N),
                     ' observaciones | ', unique(curva$N_ubicaciones), ' ubicaciones')
  nota <- if (k_max == 1L)
    'Solo hay soporte comun para k=1: esta muestra no permite comparar otros k.' else
    'Todos los puntos usan la misma muestra. El valor de k lo elige el usuario.'
  g <- ggplot(locales, aes(k, RMSE, color = Metodo)) +
    geom_point(size = 2.8) +
    geom_hline(data = referencias, aes(yintercept = RMSE, linetype = Metodo),
               color = 'grey35', linewidth = 0.6) +
    scale_color_manual(values = c('KNN' = '#1673B1', 'IDW p=1' = '#D55E00',
                                   'IDW p=2' = '#009E73')) +
    scale_linetype_manual(values = c('Media' = 'dashed', '1-NN' = 'dotted',
                                      'Siempre cero' = 'dotdash')) +
    scale_x_continuous(breaks = seq_len(k_max)) +
    labs(title = paste(gsub('_', ' ', variable), '-', escala, '2025'),
         subtitle = subtitulo, x = 'Numero de vecinos (k)',
         y = paste0('RMSE (', unidad, ')'), color = 'Metodo',
         linetype = 'Referencia', caption = nota) +
    theme_minimal(base_size = 13) +
    theme(legend.position = 'bottom', panel.grid.minor = element_blank(),
          plot.title = element_text(face = 'bold'),
          plot.caption = element_text(hjust = 0))
  if (k_max > 1L) g <- g + geom_line(linewidth = 0.9)
  if (k_max == 1L) g <- g + expand_limits(y = 0) +
    labs(caption = paste(nota, 'Con k=1, KNN e IDW coinciden con 1-NN.'))
  if (!is.na(k_elegido)) g <- g +
    geom_vline(xintercept = k_elegido, color = 'black', linewidth = 0.6) +
    labs(caption = paste0('k elegido por el usuario: ', k_elegido, '. ', nota))
  g
}

# 5. Ejecutar sensibilidad. Solo guarda curvas, errores y auditoria.
ejecutar_sensibilidad_manual <- function(raiz, escalas, anio, cobertura_minima,
                                         umbral_lluvia, k_elegido) {
  directorio <- file.path(raiz, 'outputs', 'interpolacion',
    paste0('rmse_k_manual_', anio, '_5_', format(Sys.time(), '%Y%m%d_%H%M%S')))
  if (dir.exists(directorio)) stop('La carpeta de salida ya existe.')
  dir.create(directorio, recursive = TRUE)
  horario_path <- archivo_clima_5(raiz, 'horario', anio)
  horario <- normalizar_entrada_comparacion(readRDS(horario_path), 'horario', anio)
  unidades <- c(Temperatura = 'grados C', Humedad_Relativa = '%',
    Precipitaciones = 'mm', Presion_Barometrica = 'hPa',
    Radiacion_Solar = 'W/m2', Velocidad_Viento = 'm/s')
  entradas <- auditorias <- list()
  for (escala in escalas) {
    ruta <- archivo_clima_5(raiz, escala, anio)
    dt <- if (escala == 'horario') horario else
      normalizar_entrada_comparacion(readRDS(ruta), escala, anio)
    coords <- coordenadas_comparacion(dt)
    salida <- file.path(directorio, escala)
    dir.create(salida)
    fwrite(coords, file.path(salida, 'coordenadas.csv'))
    curvas <- list()
    for (v in variables_comparacion_clima) {
      message(escala, ' / ', v)
      red <- preparar_variable_comparacion(dt, horario, coords, v, escala)
      resultado <- calcular_curva_rmse(red, v, cobertura_minima, umbral_lluvia)
      curvas[[v]] <- resultado$curva
      elegido <- k_elegido[k_elegido$Variable == v, escala]
      validar_k_manual(elegido, resultado$curva$k, v, escala)
      g <- grafico_rmse_k(resultado$curva, escala, unidades[[v]], elegido)
      for (extension in c('png', 'pdf')) {
        ggsave(file.path(salida, paste0('RMSE_k_', v, '.', extension)),
               g, width = 10, height = 6, dpi = 180, bg = 'white')
      }
      fwrite(resultado$cobertura, file.path(salida, paste0('cobertura_', v, '.csv')))
      fwrite(resultado$errores, file.path(salida, paste0('errores_', v, '.csv')))
      red$auditoria[, `:=`(Escala = escala, Variable = v)]
      auditorias[[paste(escala, v)]] <- red$auditoria
    }
    tabla <- rbindlist(curvas)
    fwrite(tabla, file.path(salida, 'curvas_rmse.csv'))
    # Tabla completa para consultar RMSE de cualquier k antes de decidir.
    sensibilidad <- dcast(tabla, Variable + Muestra + k + N + N_ubicaciones ~ Metodo,
                           value.var = 'RMSE')
    fwrite(sensibilidad, file.path(salida, paste0('sensibilidad_rmse_', escala, '.csv')))
    entradas[[escala]] <- data.table(Escala = escala, Archivo = normalizePath(ruta, winslash = '/'),
                                     MD5 = unname(tools::md5sum(ruta)))
  }
  fwrite(rbindlist(entradas), file.path(directorio, 'archivos_entrada.csv'))
  fwrite(rbindlist(auditorias), file.path(directorio, 'auditoria_datos.csv'))
  writeLines(c(
    paste('Anio:', anio), paste('Escalas:', paste(escalas, collapse = ', ')),
    'Metodo: retirar una ubicacion completa y predecir con las demas en el mismo instante.',
    'RMSE = raiz de la media de MSE por ubicacion; todas las ubicaciones pesan igual.',
    paste('Precipitacion: objetivo observado >', umbral_lluvia, 'mm. Vecinos secos conservados.'),
    'La evaluacion de lluvia es condicional: no evalua la deteccion de lluvia en periodos secos.',
    'Siempre cero es una referencia de error, no un interpolador candidato.',
    paste('Cobertura minima para el rango k (global y por ubicacion):', cobertura_minima),
    'El rango se limita por disponibilidad; no hay seleccion automatica de k ni de metodo.',
    'Todos los k dibujados se evaluan sobre exactamente la misma muestra.',
    'Mensual: doce periodos y exclusiones por imputacion pueden limitar k a 1.',
    'Solo estado OK, sin imputacion. En agregados se excluye cualquier periodo con horas imputadas.',
    paste('Trazabilidad horaria:', normalizePath(horario_path, winslash = '/')),
    paste('MD5 horario:', unname(tools::md5sum(horario_path))),
    'Los agregados no se recalculan. Pueden seguir siendo incompletos si el preprocesamiento los acepto.',
    'Colocalizados: media observada por ubicacion e instante; se retiran juntos.',
    'La seleccion posterior es exploratoria: el RMSE usado para elegir k no es una evaluacion independiente.',
    'No se ajustan ensambles, no se generan rankings ni se cruzan escalas en una tabla.',
    'Entrada sin modificar. Consultar README_comparacion_2025.md.'
  ), file.path(directorio, 'METODOLOGIA.txt'))
  writeLines(capture.output(sessionInfo()), file.path(directorio, 'sessionInfo.txt'))
  message('Resultados: ', normalizePath(directorio, winslash = '/'))
  directorio
}
