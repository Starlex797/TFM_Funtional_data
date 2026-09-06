# ==============================================================================
# RELACIÓN DIARIA/MENSUAL ENTRE log(NO2) Y COVARIABLES POR TIPOLOGÍA
# Madrid, estaciones representativas: Plaza Elíptica, Retiro y Casa de Campo
#
# Diseño:
#   - una covariable por fila y una estación por columna;
#   - puntos diarios semitransparentes para covariables continuas;
#   - boxplots para comparar días con lluvia frente a días sin lluvia;
#   - escala X común dentro de cada fila y libre entre covariables;
#   - tamaños de los grupos de lluvia anotados en los boxplots.
#
# Uso:
#   Rscript <este_script> monthly
#   Rscript <este_script> daily
#
# Salidas:
#   outputs/figures/eda/no2_covariables_tipologia_diario/
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(here)
  library(patchwork)
})

# ==============================================================================
# 1. CONFIGURACIÓN
# ==============================================================================

# La escala se recibe como primer argumento; si se omite, se conserva la versión
# mensual como comportamiento por defecto.
argumentos <- commandArgs(trailingOnly = TRUE)
escala_solicitada <- if (length(argumentos) >= 1L) {
  tolower(argumentos[1L])
} else {
  "monthly"
}
alias_escala <- c(
  daily = "diaria", monthly = "mensual",
  diaria = "diaria", mensual = "mensual"
)
if (!escala_solicitada %in% names(alias_escala)) {
  stop("Invalid scale: ", escala_solicitada, ". Use daily or monthly.")
}
ESCALA <- unname(alias_escala[escala_solicitada])

CONFIGURACIONES <- list(
  diaria = list(
    archivo = c(
      "data", "processed", "Maestro", "diario",
      "dataset_maestro_inla_2025_DIARIO.rds"
    ),
    respuesta = "dato_diario",
    tiempo = "fecha",
    escala_adjetivo = "daily",
    unidad_observacion = "day",
    alpha = 0.30,
    tamano = 0.55,
    lluvia_no = "No rain",
    lluvia_si = "Rain",
    x_lluvia = "Daily precipitation",
    unidad_precipitacion = "mm/day"
  ),
  mensual = list(
    archivo = c(
      "data", "processed", "Maestro", "mensual",
      "dataset_maestro_inla_2019_2025_MENSUAL.rds"
    ),
    respuesta = "dato_no2",
    tiempo = "fecha",
    escala_adjetivo = "monthly",
    unidad_observacion = "month",
    alpha = 0.48,
    tamano = 0.85,
    lluvia_no = "Dry month",
    lluvia_si = "Month with rain",
    x_lluvia = "Monthly precipitation",
    unidad_precipitacion = "mm/month"
  )
)

CFG <- CONFIGURACIONES[[ESCALA]]

# Regla de Tukey para detectar valores atípicos de tráfico dentro de cada
# estación. Se usa 1,5 * IQR, la convención habitual para boxplots.
FACTOR_IQR_TRAFICO <- 1.5

# Estaciones y orden de las columnas. La etiqueta incluye la tipología para que
# la comparación sustantiva sea inmediata.
ESTACIONES <- c(
  "plaza_eliptica" = "Plaza El\u00edptica\nUrban traffic",
  "retiro" = "Retiro\nUrban background",
  "casa_de_campo" = "Casa de Campo\nSuburban"
)

# Columnas crudas del maestro: conservan las unidades originales.
COVARIABLES <- data.table(
  id = c(
    "temperatura", "humedad", "precipitacion", "viento",
    "radiacion", "presion", "trafico"
  ),
  columna = c(
    "temperatura_raw",
    "humedad_relativa_raw",
    "precipitaciones_raw",
    "velocidad_viento_raw",
    "radiacion_solar_raw",
    "presion_barometrica_raw",
    "intensidad_raw"
  ),
  etiqueta = c(
    "Temperature (\u00b0C)",
    "Relative humidity (%)",
    "Precipitation (mm/day)",
    "Wind speed (m/s)",
    "Solar radiation (W/m\u00b2)",
    "Barometric pressure (mbar)",
    "Traffic intensity (veh/h)"
  )
)
COVARIABLES[id == "precipitacion", etiqueta := paste0(
  "Precipitation (", CFG$unidad_precipitacion, ")"
)]

DIR_SALIDA <- here(
  "outputs", "figures", "eda", "no2_covariables_tipologia_multiescala"
)
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. CARGA Y VALIDACIÓN DEL MAESTRO DIARIO
# ==============================================================================

archivos <- do.call(here, as.list(CFG$archivo))

faltan_archivos <- archivos[!file.exists(archivos)]
if (length(faltan_archivos) > 0L) {
  stop(
    "The following master data files were not found:\n",
    paste(faltan_archivos, collapse = "\n")
  )
}

# Convierte texto con tildes (o con marcadores <U+00XX> producidos por algunas
# configuraciones de Windows) en identificadores ASCII comparables.
normalizar_id <- function(x) {
  y <- as.character(x)
  reemplazos <- c(
    "<U+00C1>" = "A", "<U+00C9>" = "E", "<U+00CD>" = "I",
    "<U+00D3>" = "O", "<U+00DA>" = "U", "<U+00DC>" = "U",
    "<U+00D1>" = "N", "<U+00E1>" = "a", "<U+00E9>" = "e",
    "<U+00ED>" = "i", "<U+00F3>" = "o", "<U+00FA>" = "u",
    "<U+00FC>" = "u", "<U+00F1>" = "n"
  )
  for (patron in names(reemplazos)) {
    y <- gsub(patron, reemplazos[[patron]], y, fixed = TRUE)
  }

  y_ascii <- suppressWarnings(iconv(y, from = "UTF-8", to = "ASCII//TRANSLIT"))
  y[!is.na(y_ascii)] <- y_ascii[!is.na(y_ascii)]
  y <- tolower(gsub("[^A-Za-z0-9]+", "_", y))
  gsub("^_+|_+$", "", y)
}

# Se usa copy() para no modificar por referencia el objeto leído.
datos <- rbindlist(
  lapply(archivos, function(f) as.data.table(copy(readRDS(f)))),
  use.names = TRUE,
  fill = TRUE
)

# Los identificadores internos no dependen de las tildes ni de la locale.
nombres_normalizados <- normalizar_id(names(datos))
if (anyDuplicated(nombres_normalizados)) {
  stop("Column-name normalization produced duplicate names.")
}
setnames(datos, nombres_normalizados)

columnas_requeridas <- c(
  "estacion", CFG$tiempo, CFG$respuesta, COVARIABLES$columna
)
faltan_columnas <- setdiff(columnas_requeridas, names(datos))
if (length(faltan_columnas) > 0L) {
  stop(
    "The data frame is missing these required columns: ",
    paste(faltan_columnas, collapse = ", ")
  )
}

# Normalización mínima de tipos. suppressWarnings permite detectar valores no
# convertibles y tratarlos como NA durante el filtrado de pares completos.
datos[, fecha := as.Date(fecha)]
for (v in c(CFG$respuesta, COVARIABLES$columna)) {
  datos[, (v) := suppressWarnings(as.numeric(get(v)))]
}
datos[, no2 := get(CFG$respuesta)]

datos[, estacion_id := normalizar_id(estacion)]

# Comprobar la unidad observacional antes de representar "un punto por día".
duplicados <- datos[
  estacion_id %chin% names(ESTACIONES),
  .N,
  by = c("estacion_id", CFG$tiempo)
][N > 1L]
if (nrow(duplicados) > 0L) {
  stop(
    "Duplicate observations exist for the selected time unit. Clean the data ",
    "before creating the figure."
  )
}

estaciones_ausentes <- setdiff(names(ESTACIONES), unique(datos$estacion_id))
if (length(estaciones_ausentes) > 0L) {
  stop(
    "These monitoring stations were not found: ",
    paste(estaciones_ausentes, collapse = ", ")
  )
}

datos <- datos[estacion_id %chin% names(ESTACIONES)]

# El tráfico tiene niveles estructuralmente distintos según la estación. Por
# ello, sus atípicos se detectan DENTRO de cada estación y no sobre el conjunto
# mezclado. Se sustituyen por NA solo para esta figura; el maestro no se modifica.
limites_trafico <- datos[is.finite(intensidad_raw),
  {
    q1 <- quantile(intensidad_raw, 0.25, na.rm = TRUE, names = FALSE)
    q3 <- quantile(intensidad_raw, 0.75, na.rm = TRUE, names = FALSE)
    iqr <- q3 - q1
    .(
      limite_inferior = q1 - FACTOR_IQR_TRAFICO * iqr,
      limite_superior = q3 + FACTOR_IQR_TRAFICO * iqr
    )
  },
  by = estacion_id
]

datos[
  limites_trafico,
  on = "estacion_id",
  `:=`(
    limite_inferior = i.limite_inferior,
    limite_superior = i.limite_superior
  )
]
datos[, outlier_trafico :=
  is.finite(intensidad_raw) &
    (intensidad_raw < limite_inferior | intensidad_raw > limite_superior)]
n_out_trafico <- datos[outlier_trafico == TRUE, .N]
datos[outlier_trafico == TRUE, intensidad_raw := NA_real_]
datos[, c("limite_inferior", "limite_superior", "outlier_trafico") := NULL]

# log(NO2) solo está definido para concentraciones positivas. Los NA, infinitos,
# ceros y eventuales valores negativos no se dibujan.
n_no2_no_positivo <- datos[
  !is.na(no2) & is.finite(no2) & no2 <= 0,
  .N
]
if (n_no2_no_positivo > 0L) {
  warning(
    sprintf(
      "%d observations with NO2 <= 0 were excluded because log(NO2) is undefined.",
      n_no2_no_positivo
    )
  )
}

datos[
  is.finite(no2) & no2 > 0,
  log_no2 := log(no2)
]

# ==============================================================================
# 3. FORMATO LARGO Y GRUPOS DE LLUVIA
# ==============================================================================

largo <- melt(
  datos[, c("estacion_id", "fecha", "log_no2", COVARIABLES$columna), with = FALSE],
  id.vars = c("estacion_id", "fecha", "log_no2"),
  measure.vars = COVARIABLES$columna,
  variable.name = "columna_covariable",
  value.name = "x",
  variable.factor = FALSE,
  na.rm = FALSE
)

# Pares completos y finitos: así los valores faltantes no detienen el gráfico.
largo <- largo[is.finite(x) & is.finite(log_no2)]
largo <- merge(
  largo,
  COVARIABLES,
  by.x = "columna_covariable",
  by.y = "columna",
  all.x = TRUE,
  sort = FALSE
)

largo[, estacion_panel := factor(
  unname(ESTACIONES[estacion_id]),
  levels = unname(ESTACIONES)
)]
largo[, id := factor(id, levels = COVARIABLES$id)]

# La precipitación se trata como una condición binaria, coherente con la
# pregunta "¿cambia el NO2 los días de lluvia?". Se conservan ambos grupos.
# Los valores negativos, si existieran por un error de origen, se excluyen.
n_precipitacion_negativa <- largo[id == "precipitacion" & x < 0, .N]
if (n_precipitacion_negativa > 0L) {
  warning(sprintf(
    "%d negative precipitation values were excluded.",
    n_precipitacion_negativa
  ))
  largo <- largo[!(id == "precipitacion" & x < 0)]
}
largo[, grupo_lluvia := factor(
  fifelse(id == "precipitacion" & x > 0, CFG$lluvia_si, CFG$lluvia_no),
  levels = c(CFG$lluvia_no, CFG$lluvia_si)
)]

if (nrow(largo) == 0L) {
  stop("No complete and finite pairs remain for plotting.")
}

# Only the two rainfall-group sample sizes are annotated.
anotaciones_lluvia <- largo[id == "precipitacion", .(
  n_no_llueve = sum(grupo_lluvia == CFG$lluvia_no),
  n_llueve = sum(grupo_lluvia == CFG$lluvia_si)
), by = .(id, estacion_panel)]
anotaciones_lluvia[, etiqueta_n := sprintf(
  "n(dry) = %d\nn(rain) = %d", n_no_llueve, n_llueve
)]

anotaciones <- anotaciones_lluvia[, .(id, estacion_panel, etiqueta_n)]

# ==============================================================================
# 4. UNA FILA DE PANELES POR COVARIABLE
# ==============================================================================

crear_fila <- function(id_covariable, mostrar_cabeceras = FALSE) {
  datos_fila <- largo[id == id_covariable]
  anotacion_fila <- anotaciones[id == id_covariable]
  etiqueta_x <- COVARIABLES[id == id_covariable, etiqueta]

  # Al crear cada covariable como un ggplot independiente, facet_wrap() conserva
  # una única escala X para sus tres estaciones. Al apilar después los siete
  # ggplots, cada fila puede tener un rango X distinto.
  if (id_covariable == "precipitacion") {
    p <- ggplot(
      datos_fila,
      aes(x = grupo_lluvia, y = log_no2, fill = grupo_lluvia)
    ) +
      geom_boxplot(
        width = 0.56, linewidth = 0.45, alpha = 0.88,
        outlier.alpha = 0.35, outlier.size = 0.75, na.rm = TRUE
      ) +
      scale_fill_manual(
        values = setNames(
          c("#D9D9D9", "#6BAED6"),
          c(CFG$lluvia_no, CFG$lluvia_si)
        ),
        guide = "none"
      ) +
      scale_x_discrete(drop = FALSE) +
      labs(x = CFG$x_lluvia)
  } else {
    p <- ggplot(datos_fila, aes(x = x, y = log_no2)) +
      geom_point(
        color = "#2C5C85", size = CFG$tamano, alpha = CFG$alpha,
        shape = 16, na.rm = TRUE
      ) +
      scale_x_continuous(
        labels = scales::label_number(decimal.mark = ",", big.mark = "."),
        expand = expansion(mult = c(0.035, 0.055))
      ) +
      labs(x = etiqueta_x)
  }

  p <- p +
    geom_label(
      data = anotacion_fila,
      aes(x = Inf, y = Inf, label = etiqueta_n),
      inherit.aes = FALSE,
      hjust = 1.08, vjust = 1.10,
      size = 2.25, lineheight = 0.90,
      label.padding = unit(0.10, "lines"),
      label.r = unit(0.08, "lines"),
      linewidth = 0,
      fill = "white", color = "#333333"
    ) +
    facet_wrap(
      ~estacion_panel,
      nrow = 1, scales = "fixed", drop = FALSE
    ) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 0.5, decimal.mark = ","),
      expand = expansion(mult = c(0.035, 0.13))
    ) +
    labs(y = expression(log(NO[2]))) +
    theme_minimal(base_size = 8.5) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.border = element_rect(color = "grey72", fill = NA, linewidth = 0.35),
      axis.title.x = element_text(size = 7.8, margin = margin(t = 2)),
      axis.title.y = element_text(size = 7.8, margin = margin(r = 3)),
      axis.text = element_text(size = 6.8, color = "grey25"),
      axis.ticks = element_line(color = "grey55", linewidth = 0.25),
      axis.ticks.length = unit(1.2, "pt"),
      strip.background = element_rect(
        fill = "#E8EEF3", color = "grey72", linewidth = 0.35
      ),
      strip.text = element_text(
        face = "bold", size = 7.7, lineheight = 0.95,
        margin = margin(t = 3, b = 3)
      ),
      panel.spacing.x = unit(10, "pt"),
      plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 3)
    )

  # Las estaciones funcionan como cabeceras de columna y solo se muestran en
  # la primera fila para ahorrar espacio vertical.
  if (!mostrar_cabeceras) {
    p <- p + theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )
  }

  p
}

filas <- lapply(seq_len(nrow(COVARIABLES)), function(i) {
  crear_fila(
    id_covariable = COVARIABLES$id[i],
    mostrar_cabeceras = i == 1L
  )
})

fecha_min <- min(datos$fecha, na.rm = TRUE)
fecha_max <- max(datos$fecha, na.rm = TRUE)
periodo <- if (format(fecha_min, "%Y") == format(fecha_max, "%Y")) {
  format(fecha_min, "%Y")
} else {
  sprintf(
    "%s to %s", format(fecha_min, "%Y-%m"), format(fecha_max, "%Y-%m")
  )
}
periodo_archivo <- if (format(fecha_min, "%Y") == format(fecha_max, "%Y")) {
  format(fecha_min, "%Y")
} else {
  paste(format(fecha_min, "%Y"), format(fecha_max, "%Y"), sep = "_")
}

figura <- wrap_plots(filas, ncol = 1, heights = rep(1, length(filas))) +
  plot_annotation(
    title = paste0(
      if (ESCALA == "mensual") "Monthly" else "Daily",
      " relationship between NO\u2082 and its covariates in Madrid"
    ),
    subtitle = paste0(
      "Period: ", periodo,
      "  \u00b7  Source: Madrid City Council (air quality, weather and traffic)"
    ),
    caption = paste0(
      "Each point represents one ", CFG$unidad_observacion,
      "; no smoothing curves are shown. Binary precipitation: the observation ",
      "is classified as rainy when precipitation is > 0 ", CFG$unidad_precipitacion, ".\n",
      "Boxplots show the median, interquartile range and 1.5\u00d7IQR whiskers."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, color = "#1E2D3D"),
      plot.subtitle = element_text(size = 8.4, color = "grey35", margin = margin(b = 5)),
      plot.caption = element_text(
        size = 7.4, color = "grey35", hjust = 0,
        margin = margin(t = 5)
      ),
      plot.margin = margin(t = 8, r = 8, b = 7, l = 6)
    )
  )

# ==============================================================================
# 5. EXPORTACIÓN PARA LA MEMORIA
# ==============================================================================

archivo_png <- file.path(
  DIR_SALIDA,
  sprintf(
    "panel_%s_log_no2_covariates_station_type_%s_no_smoother_binary_rainfall.png",
    ESCALA, periodo_archivo
  )
)
archivo_pdf <- sub("\\.png$", ".pdf", archivo_png)

# A4 vertical: adecuado para una página completa de la memoria. El PNG se
# exporta a 320 ppp y el PDF conserva líneas y texto vectoriales.
ggsave(
  filename = archivo_png, plot = figura,
  width = 210, height = 297, units = "mm", dpi = 320, bg = "white"
)

dispositivo_pdf <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
ggsave(
  filename = archivo_pdf, plot = figura, device = dispositivo_pdf,
  width = 210, height = 297, units = "mm", bg = "white"
)

cat("\nFigure generated successfully:\n")
cat("  PNG: ", archivo_png, "\n", sep = "")
cat("  PDF: ", archivo_pdf, "\n", sep = "")
cat("  Traffic outliers excluded (1.5*IQR rule within station): ",
  n_out_trafico, "\n",
  sep = ""
)
