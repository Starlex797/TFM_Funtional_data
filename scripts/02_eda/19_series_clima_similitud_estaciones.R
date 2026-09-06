# ==============================================================================
# SIMILITUD ENTRE ESTACIONES CLIMATOLÓGICAS — HORARIO / DIARIO / MENSUAL
# Madrid 2025 · 6 variables meteorológicas · Series superpuestas por estación
#
# Objetivo: estudiar cómo de parecidos son los valores climatológicos entre
# estaciones en tres escalas temporales, como SUSTITUTO de la correlación
# espacial (con 8-25 estaciones el variograma/Moran es inestable).
#
# ¿Por qué esto y no el variograma? El variograma empírico se estima solo
# con la dimensión ESPACIAL: con n estaciones hay n(n-1)/2 pares únicos
# (8 est -> 28 pares; 25 est -> 300), que repartidos en intervalos de
# distancia dejan un puñado de pares por bin: la semivarianza resultante es
# muy inestable (la literatura recomienda ~100+ localizaciones; Webster y
# Oliver, 2007). Además exige supuestos de estacionariedad e isotropía que
# no se pueden comprobar con tan pocos puntos, y los efectos locales
# sistemáticos (exposición del anemómetro, isla de calor) contaminan la
# estimación. Este enfoque, en cambio, explota la dimensión TEMPORAL como
# replicación: cada par de estaciones comparte hasta 365 días u 8.760 horas
# de observación, de modo que su correlación se estima con miles de réplicas
# en lugar de con una sola "distancia". El resultado (r̄ por pares, series
# superpuestas y mapas por día) responde a la misma pregunta que el
# variograma —¿cuánto se parecen dos localizaciones?— pero con precisión,
# por escala temporal y sin supuestos estructurales. La estructura espacial
# formal (rango, varianza) la estima después el propio SPDE vía Matérn
# dentro del modelo, donde sí se aprovecha toda la información conjunta.
#   - Panel horario : semana de invierno + semana de verano (legible)
#   - Panel diario  : año completo
#   - Panel mensual : agregado del diario (media; suma para precipitación)
# Cuantificación: correlación de Pearson media (y mínima) entre todos los
# pares de estaciones, calculada sobre el AÑO COMPLETO en cada escala.
#   r̄ alto  -> campo espacialmente homogéneo (una estación representa bien
#              a las demás; el misalignment es poco grave para esa variable)
#   r̄ bajo  -> heterogeneidad espacial (la interpolación aporta más error)
# Mapas por variable (2 paneles): visualizan la CORRELACIÓN ESPACIAL, es
# decir, si las estaciones registran o no los mismos valores:
#   (a) media anual por estación  -> ¿coinciden los NIVELES?
#   (b) r̄ diaria de cada estación con el resto -> ¿covarían en el tiempo?
# Outputs: outputs/analysis/similitud_estaciones_clima/
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)
library(gridExtra)
library(grid)
library(sf)

# ==============================================================================
# 1. CONFIGURACIÓN
# ==============================================================================

SEMANA_INVIERNO <- seq(as.Date("2025-01-13"), as.Date("2025-01-19"), by = "day")
SEMANA_VERANO <- seq(as.Date("2025-07-14"), as.Date("2025-07-20"), by = "day")

vars_info <- list(
  list(
    col = "Temperatura", ylab = "Temperatura (°C)",
    ylab_mes = "Temperatura media (°C)", file = "temperatura", agg = "media"
  ),
  list(
    col = "Humedad_Relativa", ylab = "Humedad relativa (%)",
    ylab_mes = "Humedad relativa media (%)", file = "humedad", agg = "media"
  ),
  list(
    col = "Precipitaciones", ylab = "Precipitaciones (mm)",
    ylab_mes = "Precipitación acumulada (mm/mes)", file = "precipitacion", agg = "suma"
  ),
  list(
    col = "Presion Barométrica", ylab = "Presión barométrica (mbar)",
    ylab_mes = "Presión barométrica media (mbar)", file = "presion", agg = "media"
  ),
  list(
    col = "Radiación Solar", ylab = "Radiación solar (W/m²)",
    ylab_mes = "Radiación solar media (W/m²)", file = "radiacion", agg = "media"
  ),
  list(
    col = "Velocidad Viento", ylab = "Velocidad del viento (m/s)",
    ylab_mes = "Velocidad media del viento (m/s)", file = "viento", agg = "media"
  )
)

dir_salida <- here("outputs", "EDA", "Clima", "similitud_estaciones_clima")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. CARGA DE DATOS (red meteorológica original, SIN interpolar)
# ==============================================================================

dt_d <- readRDS(here(
  "data", "processed", "Clima", "diario",
  "meteo_madrid_2025_diario.rds"
))
setDT(dt_d)

dt_h <- readRDS(here(
  "data", "processed", "Clima", "horario",
  "meteo_madrid_2025_horario.rds"
))
setDT(dt_h)
dt_h[, Hora_num := as.integer(sub("H0?", "", as.character(HORA)))]
dt_h[, DATETIME := as.POSIXct(FECHA, tz = "UTC") + (Hora_num - 1L) * 3600L]

cat(sprintf("Diario : %d filas · %d estaciones\n", nrow(dt_d), uniqueN(dt_d$ESTACION)))
cat(sprintf("Horario: %d filas · %d estaciones\n\n", nrow(dt_h), uniqueN(dt_h$ESTACION)))

# Cartografía de distritos (misma capa que el resto de mapas del proyecto)
distritos <- st_read(here("data", "raw", "geometrias", "DISTRITOS.shp"),
  quiet = TRUE
)
# OJO: LONGITUD/LATITUD vienen corruptas para algunas estaciones en el RDS
# (p.ej. Centro Mpal. De Acústica, E.D.A.R La China); X_km/Y_km (UTM 30N)
# son correctas para las 26, así que la geometría se construye desde ellas.
coords_est <- unique(dt_d[, .(ESTACION, X_m = X_km * 1000, Y_m = Y_km * 1000)])

# ==============================================================================
# 3. AUXILIARES
# ==============================================================================

# Correlación media/mínima entre pares de estaciones (año completo).
# dt_largo: columnas <clave temporal>, ESTACION, valor
correlacion_pares <- function(dt_largo, clave_tiempo) {
  ancho <- dcast(dt_largo, paste(clave_tiempo, "~ ESTACION"), value.var = "valor")
  m <- as.matrix(ancho[, -1])
  if (ncol(m) < 2) {
    return(list(
      media = NA_real_, minima = NA_real_, n_est = ncol(m),
      r_por_estacion = NULL
    ))
  }
  cm <- suppressWarnings(cor(m, use = "pairwise.complete.obs"))
  pares <- cm[lower.tri(cm)]
  pares <- pares[!is.na(pares)]
  diag(cm) <- NA
  list(
    media = mean(pares),
    minima = min(pares),
    n_est = ncol(m),
    # r media de CADA estación con todas las demás (para el mapa)
    r_por_estacion = rowMeans(cm, na.rm = TRUE)
  )
}

tema_similitud <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = 11),
      plot.subtitle      = element_text(color = "gray40", size = 8.5),
      strip.text         = element_text(face = "bold", size = 9.5),
      strip.background   = element_rect(fill = "gray96", color = "gray80"),
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8)
    )
}

guardar_tabla_png_local <- function(df, titulo, subtitulo, ruta_png) {
  tbl_grob <- gridExtra::tableGrob(
    as.data.frame(df),
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core = list(
        fg_params = list(fontsize = 9),
        bg_params = list(fill = c("white", "#F5F5F5"), col = NA)
      ),
      colhead = list(
        fg_params = list(fontsize = 9, fontface = "bold"),
        bg_params = list(fill = "#DDEEFF", col = NA)
      )
    )
  )
  combinado <- gridExtra::arrangeGrob(
    grid::textGrob(titulo, gp = grid::gpar(fontsize = 13, fontface = "bold")),
    grid::textGrob(subtitulo, gp = grid::gpar(fontsize = 9, col = "grey40")),
    tbl_grob,
    nrow = 3,
    heights = grid::unit(c(0.5, 0.35, nrow(df) * 0.35 + 0.5), "inches")
  )
  ggsave(ruta_png, combinado,
    width = 10,
    height = 0.85 + nrow(df) * 0.35 + 0.5, dpi = 150, bg = "white"
  )
}

primeros_meses <- seq(as.Date("2025-01-01"), as.Date("2025-12-01"), by = "1 month")
etiquetas_meses <- format(primeros_meses, "%b")

# ==============================================================================
# 4. BUCLE PRINCIPAL — UNA FIGURA (3 ESCALAS) POR VARIABLE + CORRELACIONES
# ==============================================================================

resumen <- list()

for (vi in vars_info) {
  v <- vi$col
  if (!v %in% names(dt_d)) {
    cat(sprintf("  [SKIP] %s — columna no encontrada\n", v))
    next
  }

  # ---- Datos por escala (solo estaciones que miden la variable) -------------
  d_diario <- dt_d[!is.na(get(v)), .(ESTACION, FECHA, valor = get(v))]
  d_horario <- dt_h[
    !is.na(get(v)),
    .(ESTACION, FECHA, DATETIME, valor = get(v))
  ]
  d_mensual <- if (vi$agg == "suma") {
    d_diario[, .(valor = sum(valor)), by = .(ESTACION, MES = month(FECHA))]
  } else {
    d_diario[, .(valor = mean(valor)), by = .(ESTACION, MES = month(FECHA))]
  }

  estaciones <- sort(unique(d_diario$ESTACION))
  paleta <- setNames(scales::hue_pal()(length(estaciones)), estaciones)

  # ---- Correlaciones entre pares (año completo, por escala) -----------------
  cor_h <- correlacion_pares(
    d_horario[, .(FECHA, DATETIME, ESTACION, valor)],
    "DATETIME"
  )
  cor_d <- correlacion_pares(d_diario, "FECHA")
  cor_m <- correlacion_pares(d_mensual, "MES")

  resumen[[v]] <- data.table(
    Variable = vi$ylab,
    N_estaciones = cor_d$n_est,
    r_medio_horario = round(cor_h$media, 3),
    r_medio_diario = round(cor_d$media, 3),
    r_medio_mensual = round(cor_m$media, 3),
    r_min_horario = round(cor_h$minima, 3),
    r_min_diario = round(cor_d$minima, 3),
    r_min_mensual = round(cor_m$minima, 3)
  )

  # ---- Panel horario: semana invierno + semana verano ------------------------
  d_sem <- rbind(
    d_horario[FECHA %in% SEMANA_INVIERNO][, periodo := "Invierno (13-19 ene)"],
    d_horario[FECHA %in% SEMANA_VERANO][, periodo := "Verano (14-20 jul)"]
  )
  p_hor <- ggplot(d_sem, aes(DATETIME, valor, color = ESTACION)) +
    geom_line(linewidth = 0.35, alpha = 0.75, na.rm = TRUE) +
    facet_wrap(~periodo, nrow = 1, scales = "free_x") +
    scale_x_datetime(date_breaks = "1 day", date_labels = "%d %b") +
    scale_color_manual(values = paleta) +
    labs(
      title = sprintf("Escala HORARIA — r̄ entre estaciones (año completo) = %.2f", cor_h$media),
      x = NULL, y = vi$ylab
    ) +
    tema_similitud()

  # ---- Panel diario: año completo (leyenda común aquí) ------------------------
  n_col_lg <- ifelse(length(estaciones) <= 10L, 3L, 5L)
  p_dia <- ggplot(d_diario, aes(FECHA, valor, color = ESTACION)) +
    geom_line(linewidth = 0.4, alpha = 0.75, na.rm = TRUE) +
    scale_x_date(breaks = primeros_meses, labels = etiquetas_meses) +
    scale_color_manual(values = paleta, name = "Estación") +
    labs(
      title = sprintf("Escala DIARIA — r̄ = %.2f", cor_d$media),
      x = NULL, y = vi$ylab
    ) +
    tema_similitud() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 6.5),
      legend.key.width = unit(0.8, "cm")
    ) +
    guides(color = guide_legend(
      ncol = n_col_lg,
      override.aes = list(linewidth = 1.4)
    ))

  # ---- Panel mensual ----------------------------------------------------------
  p_mes <- ggplot(d_mensual, aes(MES, valor, color = ESTACION)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_x_continuous(breaks = 1:12, labels = etiquetas_meses) +
    scale_color_manual(values = paleta) +
    labs(
      title = sprintf("Escala MENSUAL — r̄ = %.2f", cor_m$media),
      x = NULL, y = vi$ylab_mes
    ) +
    tema_similitud()

  # ---- Componer figura ---------------------------------------------------------
  titulo <- grid::textGrob(
    sprintf("%s — Similitud entre estaciones · Madrid 2025", vi$ylab),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  subtitulo <- grid::textGrob(
    sprintf(
      "%d estaciones · r̄ = correlación de Pearson media entre pares (sustituto de la correlación espacial) · r mín: horario %.2f | diario %.2f | mensual %.2f",
      cor_d$n_est, cor_h$minima, cor_d$minima, cor_m$minima
    ),
    gp = grid::gpar(fontsize = 9, col = "grey40")
  )

  fig <- gridExtra::arrangeGrob(
    titulo, subtitulo, p_hor, p_dia, p_mes,
    ncol = 1, heights = c(0.45, 0.3, 3.2, 4.6, 3.0)
  )

  archivo <- file.path(
    dir_salida,
    sprintf("similitud_%s_horario_diario_mensual.png", vi$file)
  )
  ggsave(archivo, fig, width = 13, height = 12.5, dpi = 200, bg = "white")
  cat(sprintf(
    "  ✓ %-22s [%d est] r̄ horario=%.2f diario=%.2f mensual=%.2f  %s\n",
    v, cor_d$n_est, cor_h$media, cor_d$media, cor_m$media, basename(archivo)
  ))

  # ---- Mapas: ¿registran las estaciones los mismos valores? -------------------
  # Panel (a): media anual -> compara NIVELES entre estaciones.
  # Panel (b): r̄ diaria de cada estación con el resto -> compara COVARIACIÓN
  #            (lectura directa de la correlación espacial de la variable).
  est_stats <- d_diario[, .(media_anual = mean(valor)), by = ESTACION]
  est_stats[, r_media := cor_d$r_por_estacion[ESTACION]]

  sf_mapa <- st_as_sf(
    merge(est_stats, coords_est, by = "ESTACION"),
    coords = c("X_m", "Y_m"), crs = st_crs(distritos)
  )

  tema_mapa <- theme_void(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 8.5),
      legend.text     = element_text(size = 8)
    )

  p_map_val <- ggplot() +
    geom_sf(
      data = distritos, fill = "grey97", colour = "grey80",
      linewidth = 0.3
    ) +
    geom_sf(
      data = sf_mapa, aes(color = media_anual), size = 4.2,
      alpha = 0.95
    ) +
    geom_sf_text(
      data = sf_mapa,
      aes(label = sprintf("%.1f", media_anual)),
      size = 2.7, nudge_y = 1300, color = "grey15"
    ) +
    scale_color_viridis_c(name = vi$ylab) +
    labs(title = "(a) Media anual por estación — ¿mismos niveles?") +
    tema_mapa

  p_map_cor <- ggplot() +
    geom_sf(
      data = distritos, fill = "grey97", colour = "grey80",
      linewidth = 0.3
    ) +
    geom_sf(data = sf_mapa, aes(color = r_media), size = 4.2, alpha = 0.95) +
    geom_sf_text(
      data = sf_mapa,
      aes(label = sprintf("%.2f", r_media)),
      size = 2.7, nudge_y = 1300, color = "grey15"
    ) +
    scale_color_gradient(
      low = "#B2182B", high = "#2166AC",
      name = "r̄ con el\nresto (diaria)"
    ) +
    labs(title = "(b) Correlación media con el resto — ¿covarían?") +
    tema_mapa

  fig_mapa <- gridExtra::arrangeGrob(
    grid::textGrob(
      sprintf("%s — ¿Miden las estaciones lo mismo? · Madrid 2025", vi$ylab),
      gp = grid::gpar(fontsize = 13, fontface = "bold")
    ),
    grid::textGrob(
      paste(
        "(a) niveles medios anuales · (b) correlación diaria de cada",
        "estación con las demás: sustituto de la correlación espacial"
      ),
      gp = grid::gpar(fontsize = 9, col = "grey40")
    ),
    gridExtra::arrangeGrob(p_map_val, p_map_cor, ncol = 2),
    ncol = 1, heights = c(0.4, 0.3, 7)
  )

  archivo_mapa <- file.path(
    dir_salida,
    sprintf("mapa_similitud_%s.png", vi$file)
  )
  ggsave(archivo_mapa, fig_mapa,
    width = 13, height = 7, dpi = 200,
    bg = "white"
  )
  cat(sprintf("    ✓ mapa: %s\n", basename(archivo_mapa)))

  # ---- Mapas por día estacional (invierno/primavera/verano/otoño) -------------
  # Un día concreto de cada estación del año con el valor DIARIO observado en
  # cada estación de medida. Es la lectura más directa de la CORRELACIÓN
  # ESPACIAL: si en un mismo día todas las estaciones registran valores
  # similares (colores homogéneos), el campo es espacialmente coherente y una
  # estación representa bien a las demás; si hay gradientes fuertes, la
  # variable es heterogénea y la interpolación introduce más error.
  # Selección del día:
  #   - Precipitaciones : día de la temporada con lluvia en MÁS estaciones
  #                       (desempate: mayor acumulado) -> garantiza valores > 0
  #   - Velocidad Viento: día más ventoso (mayor media entre estaciones)
  #   - Resto           : día fijo a mitad de temporada (15 ene/abr/jul/oct)
  temporadas <- list(
    Invierno  = list(meses = c(12L, 1L, 2L), objetivo = as.Date("2025-01-15")),
    Primavera = list(meses = 3L:5L, objetivo = as.Date("2025-04-15")),
    Verano    = list(meses = 6L:8L, objetivo = as.Date("2025-07-15")),
    Otono     = list(meses = 9L:11L, objetivo = as.Date("2025-10-15"))
  )
  etiquetas_temporada <- c(
    Invierno = "Invierno", Primavera = "Primavera",
    Verano = "Verano", Otono = "Otoño"
  )

  elegir_dia <- function(d, meses, objetivo, modo) {
    dd <- d[month(FECHA) %in% meses]
    if (nrow(dd) == 0L) {
      return(as.Date(NA))
    }
    if (modo == "lluvia") {
      por_dia <- dd[, .(n_con = sum(valor > 0), total = sum(valor)), by = FECHA]
      por_dia <- por_dia[n_con > 0]
      if (nrow(por_dia) == 0L) {
        return(as.Date(NA))
      }
      setorder(por_dia, -n_con, -total)
      return(por_dia$FECHA[1])
    }
    if (modo == "max_media") {
      por_dia <- dd[, .(m = mean(valor)), by = FECHA]
      return(por_dia$FECHA[which.max(por_dia$m)])
    }
    fechas <- unique(dd$FECHA)
    fechas[which.min(abs(as.numeric(fechas - objetivo)))]
  }

  modo_dia <- if (v == "Precipitaciones") {
    "lluvia"
  } else if (v == "Velocidad Viento") {
    "max_media"
  } else {
    "fijo"
  }

  paneles_dia <- list()
  for (nom_temp in names(temporadas)) {
    tp <- temporadas[[nom_temp]]
    dia <- elegir_dia(d_diario, tp$meses, tp$objetivo, modo_dia)
    if (is.na(dia)) {
      cat(sprintf("    [AVISO] %s: sin día válido en %s\n", v, nom_temp))
      next
    }
    d_dia <- d_diario[FECHA == dia]
    sf_dia <- st_as_sf(
      merge(d_dia, coords_est, by = "ESTACION"),
      coords = c("X_m", "Y_m"), crs = st_crs(distritos)
    )
    rng <- range(d_dia$valor)

    paneles_dia[[nom_temp]] <- ggplot() +
      geom_sf(
        data = distritos, fill = "grey97", colour = "grey80",
        linewidth = 0.3
      ) +
      geom_sf(data = sf_dia, aes(color = valor), size = 4.4, alpha = 0.95) +
      geom_sf_text(
        data = sf_dia, aes(label = sprintf("%.1f", valor)),
        size = 2.6, nudge_y = 1300, color = "grey15"
      ) +
      scale_color_viridis_c(name = NULL) +
      labs(
        title = sprintf(
          "%s · %s", etiquetas_temporada[[nom_temp]],
          format(dia, "%d %b %Y")
        ),
        subtitle = sprintf(
          "rango: %.1f – %.1f  |  sd entre estaciones = %.2f",
          rng[1], rng[2], sd(d_dia$valor)
        )
      ) +
      tema_mapa +
      theme(plot.subtitle = element_text(
        size = 8, hjust = 0.5,
        color = "grey40"
      ))
  }

  if (length(paneles_dia) > 0) {
    fig_dias <- gridExtra::arrangeGrob(
      grid::textGrob(
        sprintf(
          "%s — Valores diarios por estación en un día de cada temporada · Madrid 2025",
          vi$ylab
        ),
        gp = grid::gpar(fontsize = 13, fontface = "bold")
      ),
      grid::textGrob(
        paste(
          "Lectura directa de la correlación espacial: colores/valores",
          "homogéneos = las estaciones registran lo mismo ese día;",
          "gradientes = heterogeneidad espacial"
        ),
        gp = grid::gpar(fontsize = 9, col = "grey40")
      ),
      gridExtra::arrangeGrob(grobs = paneles_dia, ncol = 2),
      ncol = 1, heights = c(0.4, 0.3, 11.5)
    )
    archivo_dias <- file.path(
      dir_salida, sprintf("mapa_dias_estacionales_%s.png", vi$file)
    )
    ggsave(archivo_dias, fig_dias,
      width = 12, height = 12.5, dpi = 200,
      bg = "white"
    )
    cat(sprintf("    ✓ mapa días estacionales: %s\n", basename(archivo_dias)))
  }
}

# ==============================================================================
# 5. TABLA RESUMEN DE CORRELACIONES POR VARIABLE Y ESCALA
# ==============================================================================

tabla_resumen <- rbindlist(resumen)
setorder(tabla_resumen, -r_medio_diario)
fwrite(tabla_resumen, file.path(dir_salida, "tabla_similitud_correlaciones.csv"))

guardar_tabla_png_local(
  tabla_resumen,
  titulo = "Similitud entre estaciones climáticas — correlación media entre pares",
  subtitulo = "Pearson sobre el año completo 2025 en cada escala · sustituto de la correlación espacial (pocas estaciones para variograma)",
  ruta_png = file.path(dir_salida, "tabla_similitud_correlaciones.png")
)

cat(sprintf("\n✓ Resultados en: %s\n", dir_salida))
