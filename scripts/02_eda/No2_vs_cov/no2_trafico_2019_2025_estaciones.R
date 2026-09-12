# ==============================================================================
# NO2 y tráfico 2019-2025: una estación de cada tipología (semanal y mensual)
# ==============================================================================
# Igual que el análisis por tipología, pero mostrando UNA estación representativa
# de cada tipo (en vez de la media del grupo):
#   - Tráfico:   Barrio del Pilar
#   - Fondo:     Retiro
#   - Suburbana: El Pardo
# Semanal: SOLO días laborables (L-V), para que el fin de semana no diluya la
# señal (el tráfico cae mucho en sábado/domingo). Mensual: todos los días.
# Gráfica: panel doble (NO2 arriba, tráfico abajo) por estación, con el
# confinamiento sombreado.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(sf)
  library(viridis)
})

DIR_SALIDA <- here("outputs", "figures", "no2_trafico_2019_2025")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
ANIOS <- 2019:2025
COVID_INI <- as.Date("2020-03-14")
COVID_FIN <- as.Date("2020-06-21")

# Estación representativa por tipología (fácil de cambiar).
SELECCION <- c(
  "Barrio del Pilar" = "Tráfico",
  "Retiro" = "Fondo",
  "El Pardo" = "Suburbana"
)

# ------------------------------------------------------------------------------
# 1. Cargar y apilar los maestros diarios (NO2 diario + tráfico por barrio)
# ------------------------------------------------------------------------------

lista <- lapply(ANIOS, function(a) {
  f <- here(
    "data", "processed", "Maestro", "diario",
    sprintf("dataset_maestro_inla_%d_DIARIO.rds", a)
  )
  if (!file.exists(f)) {
    return(NULL)
  }
  d <- as.data.table(readRDS(f))
  d[, .(ESTACION, FECHA, NO2 = DATO_DIARIO, Trafico = intensidad_raw)]
})
dt <- rbindlist(lista, use.names = TRUE)

# ------------------------------------------------------------------------------
# 2. Seleccionar las estaciones y etiquetarlas con su tipología
# ------------------------------------------------------------------------------

dt <- dt[ESTACION %in% names(SELECCION)]
dt[, Tipo := SELECCION[ESTACION]]
dt[, Etiqueta := paste0(ESTACION, " (", Tipo, ")")]
faltan <- setdiff(names(SELECCION), unique(dt$ESTACION))
if (length(faltan) > 0) warning("No encontradas: ", paste(faltan, collapse = ", "))

# ------------------------------------------------------------------------------
# 3. Agregación semanal (solo laborables) y mensual (todos los días)
# ------------------------------------------------------------------------------

dt[, laborable := as.integer(format(FECHA, "%u")) <= 5L] # 1=lun ... 7=dom
dt[, semana := as.Date(cut(FECHA, breaks = "week"))]
dt[, mes := as.Date(format(FECHA, "%Y-%m-01"))]

agregar <- function(d, col_periodo) {
  out <- d[, .(
    NO2 = mean(NO2, na.rm = TRUE),
    Trafico = mean(Trafico, na.rm = TRUE)
  ), by = c("Etiqueta", col_periodo)]
  setnames(out, col_periodo, "periodo")
  out[]
}

agg_semana <- agregar(dt[laborable == TRUE], "semana")
agg_mes <- agregar(dt, "mes")

# ------------------------------------------------------------------------------
# 4. Gráfica: panel doble NO2 / tráfico por estación
# ------------------------------------------------------------------------------

PALETA <- c(
  "Barrio del Pilar (Tráfico)" = "#d73027",
  "Retiro (Fondo)" = "#fc8d59",
  "El Pardo (Suburbana)" = "#4575b4"
)
NIVELES <- names(PALETA)

grafica_doble <- function(d, subtitulo) {
  dl <- melt(d,
    id.vars = c("Etiqueta", "periodo"),
    measure.vars = c("NO2", "Trafico"),
    variable.name = "Var", value.name = "valor"
  )
  dl[, Var := factor(Var,
    levels = c("NO2", "Trafico"),
    labels = c("NO₂ (µg/m³)", "Tráfico (veh/h)")
  )]
  dl[, Etiqueta := factor(Etiqueta, levels = NIVELES)]

  ggplot(dl, aes(periodo, valor, color = Etiqueta)) +
    annotate("rect",
      xmin = COVID_INI, xmax = COVID_FIN,
      ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.15
    ) +
    geom_line(linewidth = 0.55) +
    facet_grid(Var ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = PALETA) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title = "NO₂ y tráfico por estación y tipología (2019-2025)",
      subtitle = subtitulo,
      x = NULL, y = NULL, color = "Estación (tipología)",
      caption = "Franja gris = confinamiento y desescalada (mar-jun 2020). Tráfico = intensidad media por barrio."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.placement = "outside",
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

p_semana <- grafica_doble(
  agg_semana, "Media semanal — SOLO días laborables (L-V)"
)
p_mes <- grafica_doble(agg_mes, "Media mensual — todos los días")

ggsave(file.path(DIR_SALIDA, "no2_trafico_semanal_estaciones.png"),
  plot = p_semana, width = 12, height = 7, dpi = 200, bg = "white"
)
ggsave(file.path(DIR_SALIDA, "no2_trafico_mensual_estaciones.png"),
  plot = p_mes, width = 12, height = 7, dpi = 200, bg = "white"
)
fwrite(agg_semana, file.path(DIR_SALIDA, "no2_trafico_semanal_estaciones.csv"))
fwrite(agg_mes, file.path(DIR_SALIDA, "no2_trafico_mensual_estaciones.csv"))

cat("Estaciones usadas:\n")
print(SELECCION)
cat("\nGráficas guardadas en:\n", DIR_SALIDA, "\n")

# ------------------------------------------------------------------------------
# 5. NO2 Y TRAFICO HORARIOS DURANTE VARIOS DIAS CONSECUTIVOS
# ------------------------------------------------------------------------------

# Parametros editables. FECHA_INICIO_HORARIA debe pertenecer al ANIO_HORARIO.
ANIO_HORARIO <- 2019L
ESTACION_HORARIA <- "Barrio del Pilar"
FECHA_INICIO_HORARIA <- as.Date("2019-10-14")
NUM_DIAS_HORARIO <- 7L
MAX_DESFASE_HORAS <- 12L
FECHA_DIA_HORARIO <- as.Date("2019-10-14")

archivo_horario <- here(
  "data", "processed", "Maestro", as.character(ANIO_HORARIO),
  sprintf("dataset_maestro_inla_%d_HORARIO.rds", ANIO_HORARIO)
)
if (!file.exists(archivo_horario)) {
  archivo_horario <- here(
    "data", "processed", "Maestro", "horario",
    sprintf("dataset_maestro_inla_%d_HORARIO.rds", ANIO_HORARIO)
  )
}
if (!file.exists(archivo_horario)) {
  stop("No se encuentra el maestro horario para ", ANIO_HORARIO, ".")
}

maestro_horario <- as.data.table(readRDS(archivo_horario))
columnas_horarias <- c("ESTACION", "FECHA", "HORA", "DATO", "intensidad_raw")
faltan_horarias <- setdiff(columnas_horarias, names(maestro_horario))
if (length(faltan_horarias) > 0L) {
  stop(
    "Faltan columnas en el maestro horario: ",
    paste(faltan_horarias, collapse = ", ")
  )
}
if (!ESTACION_HORARIA %in% maestro_horario$ESTACION) {
  stop("La estacion no existe en el maestro horario: ", ESTACION_HORARIA)
}

fecha_fin_horaria <- FECHA_INICIO_HORARIA + NUM_DIAS_HORARIO - 1L
horario <- maestro_horario[
  ESTACION == ESTACION_HORARIA &
    FECHA >= FECHA_INICIO_HORARIA & FECHA <= fecha_fin_horaria,
  .(
    FECHA = as.Date(FECHA),
    HORA = as.integer(HORA),
    NO2 = as.numeric(DATO),
    Trafico = as.numeric(intensidad_raw)
  )
]
setorder(horario, FECHA, HORA)

if (nrow(horario) == 0L) {
  stop("No hay observaciones en el periodo horario seleccionado.")
}
if (uniqueN(horario$FECHA) != NUM_DIAS_HORARIO) {
  warning(
    "El periodo contiene ", uniqueN(horario$FECHA), " dias con datos de los ",
    NUM_DIAS_HORARIO, " solicitados."
  )
}

# HORA=1 representa el intervalo que comienza a las 00:00.
horario[, fecha_hora :=
  as.POSIXct(FECHA, tz = "UTC") + (HORA - 1L) * 3600]
horario[, `:=`(
  cambio_NO2 = NO2 - shift(NO2),
  cambio_Trafico = Trafico - shift(Trafico)
)]

# Para k > 0 se correlaciona el trafico en t con el NO2 en t+k. Por tanto, un
# maximo en k horas indica una asociacion descriptiva con respuesta posterior.
correlaciones_desfase <- rbindlist(lapply(0:MAX_DESFASE_HORAS, function(k) {
  no2_futuro <- shift(horario$NO2, n = k, type = "lead")
  cambio_no2_futuro <- shift(horario$cambio_NO2, n = k, type = "lead")
  completos_nivel <- complete.cases(horario$Trafico, no2_futuro)
  completos_cambio <- complete.cases(
    horario$cambio_Trafico, cambio_no2_futuro
  )

  data.table(
    desfase_horas = k,
    correlacion_niveles = if (sum(completos_nivel) >= 3L) {
      cor(horario$Trafico[completos_nivel], no2_futuro[completos_nivel])
    } else {
      NA_real_
    },
    correlacion_incrementos = if (sum(completos_cambio) >= 3L) {
      cor(
        horario$cambio_Trafico[completos_cambio],
        cambio_no2_futuro[completos_cambio]
      )
    } else {
      NA_real_
    },
    n_niveles = sum(completos_nivel),
    n_incrementos = sum(completos_cambio)
  )
}))

mejor_nivel <- correlaciones_desfase[which.max(correlacion_niveles)]
mejor_incremento <- correlaciones_desfase[which.max(correlacion_incrementos)]

# Dos lineas separadas, con el mismo eje temporal, evitan comparar erroneamente
# unidades distintas mediante un eje Y doble.
lineas_horarias <- melt(
  horario,
  id.vars = "fecha_hora",
  measure.vars = c("NO2", "Trafico"),
  variable.name = "Variable",
  value.name = "valor"
)
lineas_horarias[, Variable := factor(
  Variable,
  levels = c("NO2", "Trafico"),
  labels = c("NO2 (ug/m3)", "Traffic intensity (veh/h)")
)]

colores_horarios <- c(
  "NO2 (ug/m3)" = "#B13A3A",
  "Traffic intensity (veh/h)" = "#246B8E"
)

p_horario <- ggplot(
  lineas_horarias,
  aes(x = fecha_hora, y = valor, color = Variable)
) +
  geom_vline(
    xintercept = as.POSIXct(
      seq(FECHA_INICIO_HORARIA, fecha_fin_horaria + 1L, by = "day"),
      tz = "UTC"
    ),
    color = "grey82", linewidth = 0.3
  ) +
  geom_line(linewidth = 0.65, na.rm = TRUE) +
  geom_point(size = 0.8, alpha = 0.75, na.rm = TRUE) +
  facet_grid(Variable ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = colores_horarios, guide = "none") +
  scale_x_datetime(
    date_breaks = "12 hours",
    date_labels = "%d %b\n%H:%M",
    expand = expansion(mult = c(0.005, 0.01)),
    timezone = "UTC"
  ) +
  labs(
    title = "Hourly evolution of NO2 and traffic intensity",
    subtitle = sprintf(
      "%s | %s to %s | Highest level correlation: lag %d h (r = %.2f)",
      ESTACION_HORARIA,
      format(FECHA_INICIO_HORARIA, "%d %b %Y"),
      format(fecha_fin_horaria, "%d %b %Y"),
      mejor_nivel$desfase_horas,
      mejor_nivel$correlacion_niveles
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Both panels share the time axis but use their own y-axis. ",
      "A positive lag means that NO2 is observed after traffic."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey35", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = "#EEF2F4", color = "grey75"),
    strip.text.y.left = element_text(face = "bold", angle = 90),
    axis.text.x = element_text(size = 7.5),
    plot.margin = margin(8, 10, 8, 8)
  )

correlaciones_largo <- melt(
  correlaciones_desfase,
  id.vars = "desfase_horas",
  measure.vars = c("correlacion_niveles", "correlacion_incrementos"),
  variable.name = "metrica",
  value.name = "correlacion"
)
correlaciones_largo[, metrica := factor(
  metrica,
  levels = c("correlacion_niveles", "correlacion_incrementos"),
  labels = c("Hourly levels", "Hourly increases")
)]

p_desfase <- ggplot(
  correlaciones_largo,
  aes(x = desfase_horas, y = correlacion, color = metrica)
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.8) +
  scale_color_manual(
    values = c("Hourly levels" = "#6A3D9A", "Hourly increases" = "#D07A25")
  ) +
  scale_x_continuous(breaks = 0:MAX_DESFASE_HORAS) +
  labs(
    title = "Lagged correlation between traffic and subsequent NO2",
    subtitle = sprintf(
      "Maximum for increases: lag %d h (r = %.2f)",
      mejor_incremento$desfase_horas,
      mejor_incremento$correlacion_incrementos
    ),
    x = "NO2 lag after traffic (hours)",
    y = "Pearson correlation",
    color = NULL,
    caption = paste0(
      "Descriptive correlations for the selected seven-day window; ",
      "they do not establish a causal traffic effect."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey35", hjust = 0),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.margin = margin(8, 10, 8, 8)
  )

sufijo_horario <- sprintf(
  "%s_%s_%ddias",
  gsub("[^a-z0-9]+", "_", tolower(ESTACION_HORARIA)),
  format(FECHA_INICIO_HORARIA, "%Y%m%d"),
  NUM_DIAS_HORARIO
)
ruta_lineas_horarias <- file.path(
  DIR_SALIDA,
  paste0("no2_trafico_lineas_horarias_", sufijo_horario, ".png")
)
ruta_desfases <- file.path(
  DIR_SALIDA,
  paste0("no2_trafico_correlacion_desfases_", sufijo_horario, ".png")
)
ruta_correlaciones <- file.path(
  DIR_SALIDA,
  paste0("no2_trafico_correlacion_desfases_", sufijo_horario, ".csv")
)

ggsave(
  ruta_lineas_horarias,
  plot = p_horario, width = 12, height = 7.2, dpi = 300, bg = "white"
)
ggsave(
  ruta_desfases,
  plot = p_desfase, width = 8.5, height = 5.2, dpi = 300, bg = "white"
)
fwrite(correlaciones_desfase, ruta_correlaciones)

cat("\nAnalisis horario:\n")
cat("  Estacion: ", ESTACION_HORARIA, "\n", sep = "")
cat(
  "  Periodo: ", format(FECHA_INICIO_HORARIA), " a ",
  format(fecha_fin_horaria), "\n",
  sep = ""
)
cat(sprintf(
  "  Niveles: correlacion maxima con NO2 %d h despues (r = %.3f).\n",
  mejor_nivel$desfase_horas, mejor_nivel$correlacion_niveles
))
cat(sprintf(
  "  Incrementos: correlacion maxima con NO2 %d h despues (r = %.3f).\n",
  mejor_incremento$desfase_horas, mejor_incremento$correlacion_incrementos
))
cat("  Lineas: ", ruta_lineas_horarias, "\n", sep = "")
cat("  Desfases: ", ruta_desfases, "\n", sep = "")
cat("  Tabla: ", ruta_correlaciones, "\n", sep = "")

# ------------------------------------------------------------------------------
# 6. COMPARACION DE LAS 24 HORAS DE UN UNICO DIA LABORABLE
# ------------------------------------------------------------------------------

if (as.integer(format(FECHA_DIA_HORARIO, "%u")) > 5L) {
  warning("FECHA_DIA_HORARIO no corresponde a un dia laborable.")
}

horario_dia <- maestro_horario[
  ESTACION == ESTACION_HORARIA & FECHA == FECHA_DIA_HORARIO,
  .(
    hora = as.integer(HORA) - 1L,
    NO2 = as.numeric(DATO),
    Trafico = as.numeric(intensidad_raw)
  )
]
setorder(horario_dia, hora)

if (nrow(horario_dia) != 24L || !identical(horario_dia$hora, 0:23)) {
  warning(
    "El dia seleccionado no contiene exactamente las 24 horas esperadas. ",
    "Se representaran las horas disponibles."
  )
}

dia_largo <- melt(
  horario_dia,
  id.vars = "hora",
  measure.vars = c("NO2", "Trafico"),
  variable.name = "Variable",
  value.name = "valor"
)
dia_largo[, valor_estandarizado := {
  desviacion <- sd(valor, na.rm = TRUE)
  if (is.finite(desviacion) && desviacion > 0) {
    (valor - mean(valor, na.rm = TRUE)) / desviacion
  } else {
    NA_real_
  }
}, by = Variable]
dia_largo[, Variable := factor(
  Variable,
  levels = c("NO2", "Trafico"),
  labels = c("NO2", "Traffic intensity")
)]

p_dia_horario <- ggplot(
  dia_largo,
  aes(x = hora, y = valor_estandarizado, color = Variable)
) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(
    values = c("NO2" = "#B13A3A", "Traffic intensity" = "#246B8E")
  ) +
  scale_x_continuous(
    breaks = seq(0, 23, by = 2),
    minor_breaks = 0:23,
    labels = sprintf("%02d:00", seq(0, 23, by = 2)),
    limits = c(0, 23),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Hourly NO2 and traffic intensity during one working day",
    subtitle = sprintf(
      "%s | %s | All 24 hourly observations",
      ESTACION_HORARIA,
      format(FECHA_DIA_HORARIO, "%d %b %Y")
    ),
    x = "Hour of day",
    y = "Standardized value (z-score)",
    color = NULL,
    caption = paste0(
      "Each series is standardized separately. Values above zero are above ",
      "that variable's daily mean; the comparison concerns timing, not magnitude."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9.5, color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey35", hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    legend.position = "top",
    axis.text.x = element_text(size = 8),
    plot.margin = margin(8, 10, 8, 8)
  )

ruta_dia_horario <- file.path(
  DIR_SALIDA,
  sprintf(
    "no2_trafico_un_dia_horario_%s_%s.png",
    gsub("[^a-z0-9]+", "_", tolower(ESTACION_HORARIA)),
    format(FECHA_DIA_HORARIO, "%Y%m%d")
  )
)
ggsave(
  ruta_dia_horario,
  plot = p_dia_horario, width = 11, height = 5.8, dpi = 300, bg = "white"
)

cat("\nGrafica de un dia laborable:\n")
cat("  Fecha: ", format(FECHA_DIA_HORARIO), "\n", sep = "")
cat("  Horas representadas: ", nrow(horario_dia), "\n", sep = "")
cat("  PNG: ", ruta_dia_horario, "\n", sep = "")

# ------------------------------------------------------------------------------
# 6. MAPA HORARIO DE UN DIA LABORABLE
# ------------------------------------------------------------------------------

# Dia editable. Se muestran sus 24 horas en una cuadricula de mapas.
FECHA_MAPA <- as.Date("2019-10-16")

if (as.integer(format(FECHA_MAPA, "%u")) > 5L) {
  warning("FECHA_MAPA corresponde a un fin de semana.")
}
if (as.integer(format(FECHA_MAPA, "%Y")) != ANIO_HORARIO) {
  stop("FECHA_MAPA debe pertenecer a ANIO_HORARIO.")
}

archivo_trafico_horario <- here(
  "data", "processed", "Trafico", "Horario_Barrio",
  as.character(ANIO_HORARIO),
  sprintf("trafico_madrid_%d_horario_barrio1.rds", ANIO_HORARIO)
)
archivo_barrios <- here("data", "raw", "geometrias", "BARRIOS.shp")
if (!file.exists(archivo_trafico_horario)) {
  stop("No se encuentra el trafico horario por barrio: ", archivo_trafico_horario)
}
if (!file.exists(archivo_barrios)) {
  stop("No se encuentra la geometria de barrios: ", archivo_barrios)
}

normalizar_barrio <- function(x) {
  y <- iconv(tolower(trimws(as.character(x))), to = "ASCII//TRANSLIT")
  gsub("[^a-z0-9]", "", y)
}

trafico_horario <- as.data.table(readRDS(archivo_trafico_horario))
trafico_dia <- trafico_horario[
  FECHA == FECHA_MAPA,
  .(Trafico = mean(intensidad, na.rm = TRUE)),
  by = .(barrio, HORA)
]
trafico_dia[, clave_barrio := normalizar_barrio(barrio)]
trafico_dia[, Hora_panel := factor(
  HORA,
  levels = 1:24,
  labels = sprintf("H%02d", 1:24)
)]

if (uniqueN(trafico_dia$HORA) != 24L) {
  warning("El trafico no contiene las 24 horas para FECHA_MAPA.")
}

barrios_sf <- st_make_valid(st_read(archivo_barrios, quiet = TRUE))
barrios_sf$clave_barrio <- normalizar_barrio(barrios_sf$NOMBRE)
barrios_dia_sf <- merge(
  barrios_sf,
  trafico_dia[, .(clave_barrio, Hora_panel, Trafico)],
  by = "clave_barrio",
  all.x = FALSE,
  sort = FALSE
)

no2_dia <- maestro_horario[
  FECHA == FECHA_MAPA,
  .(
    FECHA = as.Date(FECHA),
    ESTACION,
    LONGITUD = as.numeric(LONGITUD),
    LATITUD = as.numeric(LATITUD),
    HORA = as.integer(HORA),
    NO2 = as.numeric(DATO)
  )
][is.finite(LONGITUD) & is.finite(LATITUD) & is.finite(NO2)]
no2_dia[, Hora_panel := factor(
  HORA,
  levels = 1:24,
  labels = sprintf("H%02d", 1:24)
)]

if (nrow(no2_dia) == 0L) {
  stop("No hay mediciones de NO2 para FECHA_MAPA.")
}
if (uniqueN(no2_dia$HORA) != 24L) {
  warning("El NO2 no contiene las 24 horas para FECHA_MAPA.")
}

no2_dia_sf <- st_transform(
  st_as_sf(
    no2_dia,
    coords = c("LONGITUD", "LATITUD"),
    crs = 4326,
    remove = FALSE
  ),
  st_crs(barrios_dia_sf)
)

p_mapa_horario <- ggplot() +
  geom_sf(
    data = barrios_dia_sf,
    aes(fill = Trafico),
    color = "white",
    linewidth = 0.05
  ) +
  geom_sf(
    data = no2_dia_sf,
    aes(size = NO2),
    shape = 21,
    fill = "#D1495B",
    color = "white",
    stroke = 0.25,
    alpha = 0.88
  ) +
  facet_wrap(~Hora_panel, ncol = 6, drop = FALSE) +
  scale_fill_viridis_c(
    option = "C",
    direction = -1,
    trans = "sqrt",
    name = "Traffic\n(veh/h)",
    na.value = "grey90"
  ) +
  scale_size_continuous(
    range = c(0.35, 3.2),
    limits = range(no2_dia$NO2, na.rm = TRUE),
    name = "NO2\n(ug/m3)"
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Hourly spatial evolution of NO2 and traffic intensity in Madrid",
    subtitle = paste0(
      format(FECHA_MAPA, "%A, %d %B %Y"),
      " | Neighbourhood colour: traffic intensity | Circles: monitoring-station NO2"
    ),
    caption = paste0(
      "H01 represents 00:00-01:00 and H24 represents 23:00-24:00. ",
      "Traffic is aggregated by neighbourhood."
    )
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#1E2D3D"),
    plot.subtitle = element_text(size = 9.2, color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey35", hjust = 0),
    strip.background = element_rect(
      fill = "#EEF2F4", color = "grey72", linewidth = 0.3
    ),
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = unit(2.5, "mm"),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

ruta_mapa_horario <- file.path(
  DIR_SALIDA,
  paste0(
    "mapa_no2_trafico_24_horas_",
    format(FECHA_MAPA, "%Y%m%d"),
    ".png"
  )
)
ruta_datos_mapa <- sub("\\.png$", ".csv", ruta_mapa_horario)

ggsave(
  ruta_mapa_horario,
  plot = p_mapa_horario,
  width = 16, height = 12, dpi = 300, bg = "white"
)
fwrite(
  no2_dia[, .(FECHA, HORA, ESTACION, LONGITUD, LATITUD, NO2)],
  ruta_datos_mapa
)

cat("\nMapa horario de un dia laborable:\n")
cat("  Fecha: ", format(FECHA_MAPA), "\n", sep = "")
cat("  Horas de trafico: ", uniqueN(trafico_dia$HORA), "\n", sep = "")
cat("  Horas de NO2: ", uniqueN(no2_dia$HORA), "\n", sep = "")
cat("  Mapa: ", ruta_mapa_horario, "\n", sep = "")
