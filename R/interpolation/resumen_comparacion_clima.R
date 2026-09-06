# Tablas por escala usando exclusivamente el k escrito por el usuario.
validar_k_manual <- function(k, disponibles, variable, escala) {
  if (length(k) != 1L) stop('Falta una eleccion unica para ', variable, '/', escala)
  if (is.na(k)) return(invisible(NULL))
  if (!is.numeric(k) || !is.finite(k) || k != as.integer(k) || !k %in% disponibles) {
    stop('k no representado para ', variable, '/', escala,
         '. Valores disponibles: ', paste(sort(unique(disponibles)), collapse = ', '))
  }
}

tabla_escala_manual <- function(curvas, k_elegido, escala) {
  rbindlist(lapply(unique(curvas$Variable), function(v) {
    d <- curvas[Variable == v]
    seleccionado <- k_elegido[k_elegido$Variable == v, escala]
    validar_k_manual(seleccionado, d$k, v, escala)
    referencia <- d[k == min(d$k)]
    rmse_local <- function(metodo) {
      if (is.na(seleccionado)) return(NA_real_)
      d[Metodo == metodo & k == seleccionado, RMSE]
    }
    data.table(
      Variable = v, k_elegido = as.integer(seleccionado),
      Estado = if (is.na(seleccionado)) 'Pendiente' else 'Elegido por el usuario',
      N = unique(d$N), N_ubicaciones = unique(d$N_ubicaciones),
      RMSE_Media = referencia[Metodo == 'Media', RMSE],
      RMSE_1NN = referencia[Metodo == '1-NN', RMSE],
      RMSE_KNN = rmse_local('KNN'), RMSE_IDW_p1 = rmse_local('IDW p=1'),
      RMSE_IDW_p2 = rmse_local('IDW p=2'),
      RMSE_siempre_cero = if (v == 'Precipitaciones')
        referencia[Metodo == 'Siempre cero', RMSE] else NA_real_
    )
  }))
}

exportar_tablas_manuales <- function(raiz, k_elegido, directorio = NULL) {
  if (is.null(directorio)) {
    carpetas <- list.dirs(file.path(raiz, 'outputs', 'interpolacion'), recursive = FALSE)
    carpetas <- carpetas[grepl('^rmse_k_manual_', basename(carpetas))]
    carpetas <- carpetas[file.exists(file.path(carpetas, 'METODOLOGIA.txt'))]
    completas <- vapply(carpetas, function(p) all(file.exists(
      file.path(p, c('horario', 'diario', 'mensual'), 'curvas_rmse.csv'))), logical(1))
    carpetas <- carpetas[completas]
    if (!length(carpetas)) stop('Ejecuta primero la sensibilidad RMSE-k de las tres escalas.')
    directorio <- carpetas[which.max(file.info(file.path(carpetas, 'METODOLOGIA.txt'))$mtime)]
  }
  unidades <- c(Temperatura = 'grados C', Humedad_Relativa = '%',
    Precipitaciones = 'mm', Presion_Barometrica = 'hPa',
    Radiacion_Solar = 'W/m2', Velocidad_Viento = 'm/s')
  # Cada escala se lee y exporta por separado.
  for (escala in c('horario', 'diario', 'mensual')) {
    salida <- file.path(directorio, escala)
    ruta <- file.path(salida, 'curvas_rmse.csv')
    if (!file.exists(ruta)) next
    curvas <- fread(ruta)
    tabla <- tabla_escala_manual(curvas, k_elegido, escala)
    fwrite(tabla, file.path(salida, paste0('tabla_rmse_', escala, '.csv')))
    visual <- copy(tabla)[, c('N_ubicaciones', 'Estado') := NULL]
    visual[, Variable := gsub('_', ' ', Variable)]
    for (col in grep('^RMSE', names(visual), value = TRUE)) {
      set(visual, j = col, value = ifelse(is.na(visual[[col]]), '-', sprintf('%.3f', visual[[col]])))
    }
    visual[, k_elegido := ifelse(is.na(k_elegido), '-', as.character(k_elegido))]
    setnames(visual,
      c('k_elegido', 'RMSE_Media', 'RMSE_1NN', 'RMSE_KNN', 'RMSE_IDW_p1',
        'RMSE_IDW_p2', 'RMSE_siempre_cero'),
      c('k', 'Media', '1-NN', 'KNN', 'IDW p=1', 'IDW p=2', 'Siempre 0'))
    # El valor cero solo es una referencia pertinente para la fila de lluvia.
    grob <- gridExtra::tableGrob(as.data.frame(visual), rows = NULL,
                              theme = gridExtra::ttheme_minimal(base_size = 12))
    figura <- gridExtra::arrangeGrob(
      grid::textGrob(paste('Comparacion de RMSE - escala', escala, '- 2025'),
                     gp = grid::gpar(fontsize = 15, fontface = 'bold')),
      grob,
      grid::textGrob(paste('Precipitacion: solo objetivos con lluvia.',
                          'Guion: k pendiente o referencia no aplicable.'),
                     gp = grid::gpar(fontsize = 10)),
      ncol = 1, heights = c(0.6, 3.8, 0.6))
    for (extension in c('png', 'pdf')) {
      ggsave(file.path(salida, paste0('tabla_rmse_', escala, '.', extension)), figura,
             width = 12.5, height = 4.5, dpi = 180, bg = 'white')
    }
    # Si se cambia el k elegido, actualizar la marca vertical sin recalcular CV.
    for (v in unique(curvas$Variable)) {
      k <- k_elegido[k_elegido$Variable == v, escala]
      g <- grafico_rmse_k(curvas[Variable == v], escala, unidades[[v]], k)
      for (extension in c('png', 'pdf')) {
        ggsave(file.path(salida, paste0('RMSE_k_', v, '.', extension)),
               g, width = 10, height = 6, dpi = 180, bg = 'white')
      }
    }
    message('Tabla ', escala, ': ', file.path(salida, paste0('tabla_rmse_', escala, '.csv')))
  }
  fwrite(as.data.table(k_elegido), file.path(directorio, 'k_elegido_usuario.csv'))
  invisible(directorio)
}
