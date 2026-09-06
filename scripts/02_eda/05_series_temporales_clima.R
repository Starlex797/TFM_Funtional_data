# ==============================================================================
# SERIES TEMPORALES CLIMÁTICAS — DIARIO Y HORARIO POR ESTACIÓN
# Madrid 2025 · 6 variables · Delimitadores mensuales
# Outputs: outputs/analysis/series_temporales_clima/
# ==============================================================================
# Describe la evolución temporal de las covariables climáticas que pueden explicar
# el comportamiento del No2. Ayuda a detectar estacionalidad, episodios de lluvia, ciclos de temperatura. Etc.


library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CONFIGURACIÓN DE VARIABLES
# ==============================================================================

vars_clima <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  "Presion Barométrica", "Radiación Solar", "Velocidad Viento"
)

labels_clima <- c(
  "Temperatura (°C)",
  "Humedad Relativa (%)",
  "Precipitaciones (mm)",
  "Presión Barométrica (hPa)",
  "Radiación Solar (W/m²)",
  "Velocidad Viento (m/s)"
)

# Mapa nombre variable → etiqueta legible
dict_labels <- setNames(labels_clima, vars_clima)

# --- Emparejado de nombres entre convenciones ---------------------------------
# En el proyecto conviven dos formas de nombrar las mismas variables:
#   "Presion_Barometrica"  -> R/utilities/dictionaries.R, que es lo que acaba en
#                             meteo_madrid_<anio>_diario<N>.rds
#   "Presion Barométrica"  -> dataset maestro, y lo que declara vars_clima arriba
# Un intersect() literal entre ambas devuelve solo 3 de las 6 variables y las
# otras 3 desaparecen sin aviso. Se emparejan ignorando acentos, espacios y
# guiones bajos para que el script funcione con cualquiera de las dos.
normalizar_var <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  tolower(gsub("[^[:alnum:]]", "", x))
}

# Devuelve un vector con nombres = nombre canonico y valores = nombre real en
# los datos, para poder traducir en los dos sentidos.
resolver_vars <- function(canonicas, disponibles) {
  idx <- match(normalizar_var(canonicas), normalizar_var(disponibles))
  ok <- !is.na(idx)
  setNames(disponibles[idx[ok]], canonicas[ok])
}

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

dt_diario <- readRDS(here(
  "data", "processed", "Clima", "diario",
  "meteo_madrid_2025_diario5.rds"
))

dt_horario <- readRDS(here(
  "data", "processed", "Clima", "horario",
  "meteo_madrid_2025_horario5.rds"
))

# Variables disponibles en los datos, emparejadas por nombre normalizado
map_diario <- resolver_vars(vars_clima, names(dt_diario))
map_horario <- resolver_vars(vars_clima, names(dt_horario))

vars_diario <- unname(map_diario)
vars_horario <- unname(map_horario)

# Etiqueta legible indexada por el nombre REAL de la columna, que es el que
# viajara en la columna Variable despues del melt().
labels_diario <- setNames(dict_labels[names(map_diario)], vars_diario)
labels_horario <- setNames(dict_labels[names(map_horario)], vars_horario)

cat(sprintf("Variables diario  : %s\n", paste(vars_diario, collapse = ", ")))
cat(sprintf("Variables horario : %s\n", paste(vars_horario, collapse = ", ")))

if (length(vars_diario) < length(vars_clima)) {
  cat(sprintf(
    "  Aviso: no se encontraron en el diario: %s\n",
    paste(setdiff(vars_clima, names(map_diario)), collapse = ", ")
  ))
}

# ==============================================================================
# 3. TRANSFORMAR A FORMATO LARGO
# ==============================================================================

# --- Diario ---
dt_diario_largo <- melt(
  dt_diario,
  id.vars = c("ESTACION", "FECHA"),
  measure.vars = vars_diario,
  variable.name = "Variable",
  value.name = "Valor"
)
dt_diario_largo[, Variable := as.character(Variable)]
dt_diario_largo[, Label := labels_diario[Variable]]
dt_diario_largo[, Label := factor(Label, levels = labels_clima)]
setorder(dt_diario_largo, ESTACION, Variable, FECHA)

# --- Horario ---
dt_horario[, Hora_num := as.integer(sub("H0?", "", as.character(HORA)))]
dt_horario[, DATETIME := as.POSIXct(FECHA) + (Hora_num - 1L) * 3600L]

dt_horario_largo <- melt(
  dt_horario,
  id.vars = c("ESTACION", "FECHA", "HORA", "DATETIME"),
  measure.vars = vars_horario,
  variable.name = "Variable",
  value.name = "Valor"
)
dt_horario_largo[, Variable := as.character(Variable)]
dt_horario_largo[, Label := labels_horario[Variable]]
dt_horario_largo[, Label := factor(Label, levels = labels_clima)]
setorder(dt_horario_largo, ESTACION, Variable, DATETIME)

# ==============================================================================
# 4. DELIMITADORES MENSUALES
# ==============================================================================

primeros_meses <- seq(as.Date("2025-01-01"), as.Date("2025-12-01"), by = "1 month")
primeros_meses_ct <- as.POSIXct(primeros_meses)
etiquetas_meses <- format(primeros_meses, "%b")

# Franjas alternadas (mes par con fondo gris)
franjas_diario <- data.frame(
  xmin  = primeros_meses,
  xmax  = c(primeros_meses[-1], as.Date("2026-01-01")),
  shade = (seq_along(primeros_meses) %% 2L == 0L)
) |> subset(shade)

meses_limites_ct <- c(primeros_meses_ct, as.POSIXct("2026-01-01"))
franjas_horario <- data.frame(
  xmin  = meses_limites_ct[-length(meses_limites_ct)],
  xmax  = meses_limites_ct[-1],
  shade = (seq_along(primeros_meses) %% 2L == 0L)
) |> subset(shade)

# ==============================================================================
# 5. DIRECTORIOS DE SALIDA
# ==============================================================================

dir_base <- here("outputs", "EDA", "Clima", "series_temporales_clima")
dir_diario <- file.path(dir_base, "diario")
dir_horario <- file.path(dir_base, "horario")
dir_paneles <- file.path(dir_base, "paneles")

for (d in c(dir_diario, dir_horario, dir_paneles)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 6. FUNCIÓN AUXILIAR
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", tolower(x))
}

# Paleta fija por variable (orden según labels_clima)
paleta_vars <- c(
  "Temperatura (°C)"          = "#e74c3c",
  "Humedad Relativa (%)"      = "#2980b9",
  "Precipitaciones (mm)"      = "#27ae60",
  "Presión Barométrica (hPa)" = "#8e44ad",
  "Radiación Solar (W/m²)"    = "#f39c12",
  "Velocidad Viento (m/s)"    = "#16a085"
)

# ==============================================================================
# 7. BUCLE — SERIES DIARIAS POR ESTACIÓN (1 PNG por estación, 6 facetas)
# ==============================================================================

estaciones_unicas <- sort(unique(dt_diario_largo$ESTACION))
cat(sprintf("\nEstaciones: %d\n", length(estaciones_unicas)))

cat("\n--- Generando series DIARIAS por estación ---\n")

for (est in estaciones_unicas) {
  dt_est <- dt_diario_largo[ESTACION == est & !is.na(Valor)]

  if (nrow(dt_est) == 0L) {
    cat(sprintf("  [SKIP] %s\n", est))
    next
  }

  n_vars_ok <- uniqueN(dt_est$Variable)

  p <- ggplot(dt_est, aes(x = FECHA, y = Valor, color = Label)) +

    # Franjas mensuales de fondo
    geom_rect(
      data = franjas_diario,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.55
    ) +

    # Líneas divisorias mensuales
    geom_vline(
      xintercept = as.numeric(primeros_meses),
      color = "gray60", linewidth = 0.35, linetype = "dashed"
    ) +

    # Serie temporal
    geom_line(linewidth = 0.55, alpha = 0.85) +
    geom_point(size = 0.8, alpha = 0.5) +
    scale_color_manual(values = paleta_vars, guide = "none") +
    scale_x_date(
      breaks       = primeros_meses,
      labels       = etiquetas_meses,
      minor_breaks = NULL,
      expand       = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    facet_wrap(~Label, ncol = 2, scales = "free_y") +
    labs(
      title    = paste0("Series Temporales Diarias \u2014 ", est),
      subtitle = "Madrid 2025  \u00b7  Variables clim\u00e1ticas",
      x        = NULL,
      y        = "Valor",
      caption  = "Franjas grises: meses pares  \u00b7  L\u00edneas discontinuas: inicio de mes"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "gray40", size = 9),
      plot.caption       = element_text(color = "gray55", size = 7.5),
      strip.text         = element_text(face = "bold", size = 9),
      strip.background   = element_rect(fill = "gray96", color = "gray80"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(1, "lines")
    )

  nombre_file <- file.path(
    dir_diario,
    paste0("diario_", limpiar_nombre(est), ".png")
  )
  ggsave(nombre_file, plot = p, width = 14, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 8. BUCLE — SERIES HORARIAS POR ESTACIÓN (1 PNG por estación, 6 facetas)
# ==============================================================================

cat("\n--- Generando series HORARIAS por estación ---\n")

for (est in estaciones_unicas) {
  dt_est <- dt_horario_largo[ESTACION == est & !is.na(Valor)]

  if (nrow(dt_est) == 0L) {
    cat(sprintf("  [SKIP] %s\n", est))
    next
  }

  n_horas <- nrow(dt_est)

  p <- ggplot(dt_est, aes(x = DATETIME, y = Valor, color = Label)) +

    # Franjas mensuales de fondo
    geom_rect(
      data = franjas_horario,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.55
    ) +

    # Líneas divisorias mensuales
    geom_vline(
      xintercept = primeros_meses_ct,
      color = "gray60", linewidth = 0.3, linetype = "dashed"
    ) +

    # Serie horaria
    geom_line(linewidth = 0.3, alpha = 0.75) +
    scale_color_manual(values = paleta_vars, guide = "none") +
    scale_x_datetime(
      breaks       = primeros_meses_ct,
      labels       = etiquetas_meses,
      minor_breaks = NULL,
      expand       = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    facet_wrap(~Label, ncol = 2, scales = "free_y") +
    labs(
      title    = paste0("Series Temporales Horarias \u2014 ", est),
      subtitle = "Madrid 2025  \u00b7  Variables clim\u00e1ticas",
      x        = NULL,
      y        = "Valor",
      caption  = "Franjas grises: meses pares  \u00b7  L\u00edneas discontinuas: inicio de mes"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "gray40", size = 9),
      plot.caption       = element_text(color = "gray55", size = 7.5),
      strip.text         = element_text(face = "bold", size = 9),
      strip.background   = element_rect(fill = "gray96", color = "gray80"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(1, "lines")
    )

  nombre_file <- file.path(
    dir_horario,
    paste0("horario_", limpiar_nombre(est), ".png")
  )
  ggsave(nombre_file, plot = p, width = 14, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 9. PANELES RESUMEN — TODAS LAS ESTACIONES POR VARIABLE (diario y horario)
# ==============================================================================

cat("\n--- Generando paneles resumen por variable ---\n")

for (v in vars_diario) {
  lbl <- labels_diario[[v]]

  # ── Panel diario ────────────────────────────────────────────────────────────
  dt_v <- dt_diario_largo[Variable == v & !is.na(Valor)]

  if (nrow(dt_v) > 0L) {
    p_d <- ggplot(dt_v, aes(x = FECHA, y = Valor)) +
      geom_rect(
        data = franjas_diario,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "gray88", alpha = 0.55
      ) +
      geom_vline(
        xintercept = as.numeric(primeros_meses),
        color = "gray60", linewidth = 0.3, linetype = "dashed"
      ) +
      geom_line(color = paleta_vars[[lbl]], linewidth = 0.45, alpha = 0.85) +
      scale_x_date(
        breaks = primeros_meses, labels = etiquetas_meses,
        minor_breaks = NULL, expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      facet_wrap(~ESTACION, ncol = 4, scales = "free_y") +
      labs(
        title    = paste0("Serie Diaria: ", lbl, " \u2014 Todas las Estaciones"),
        subtitle = "Madrid 2025",
        x        = NULL,
        y        = lbl
      ) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(color = "gray40", size = 9),
        strip.text         = element_text(face = "bold", size = 7),
        axis.text.x        = element_text(angle = 60, hjust = 1, size = 6),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.spacing      = unit(0.7, "lines")
      )

    nombre_file <- file.path(
      dir_paneles,
      paste0("panel_diario_", limpiar_nombre(v), ".png")
    )
    ggsave(nombre_file, plot = p_d, width = 22, height = 18, dpi = 200, bg = "white")
    cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
  }

  # ── Panel horario ───────────────────────────────────────────────────────────
  dt_vh <- dt_horario_largo[Variable == v & !is.na(Valor)]

  if (nrow(dt_vh) > 0L) {
    p_h <- ggplot(dt_vh, aes(x = DATETIME, y = Valor)) +
      geom_rect(
        data = franjas_horario,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "gray88", alpha = 0.55
      ) +
      geom_vline(
        xintercept = primeros_meses_ct,
        color = "gray60", linewidth = 0.25, linetype = "dashed"
      ) +
      geom_line(color = paleta_vars[[lbl]], linewidth = 0.25, alpha = 0.75) +
      scale_x_datetime(
        breaks = primeros_meses_ct, labels = etiquetas_meses,
        minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      facet_wrap(~ESTACION, ncol = 4, scales = "free_y") +
      labs(
        title    = paste0("Serie Horaria: ", lbl, " \u2014 Todas las Estaciones"),
        subtitle = "Madrid 2025",
        x        = NULL,
        y        = lbl
      ) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(color = "gray40", size = 9),
        strip.text         = element_text(face = "bold", size = 7),
        axis.text.x        = element_text(angle = 60, hjust = 1, size = 6),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.spacing      = unit(0.7, "lines")
      )

    nombre_file <- file.path(
      dir_paneles,
      paste0("panel_horario_", limpiar_nombre(v), ".png")
    )
    ggsave(nombre_file, plot = p_h, width = 22, height = 18, dpi = 200, bg = "white")
    cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
  }
}

# ==============================================================================
# 10. ESTRUCTURA ESPACIAL: OFFSET FIJO DE ESTACIÓN vs. RUIDO METEOROLÓGICO
# ==============================================================================
# Las series de las secciones anteriores muestran que las estaciones no se
# solapan: hay uno o dos grados de separacion permanente entre unas y otras.
# Aqui se cuantifica que parte de esa separacion es una CONSTANTE de la estacion
# (cota, entorno urbano) y que parte es variacion meteorologica del dia.
#
# Descomposicion, por variable:
#
#   Valor(estacion, dia) = media de Madrid(dia) + offset(estacion) + residuo
#
#   offset(estacion) = media anual de su desviacion respecto a la media diaria
#                      -> componente ESTRUCTURAL, fija todo el ano
#   residuo          = lo que queda dia a dia
#                      -> componente METEOROLOGICA, lo unico realmente aleatorio
#
# Para que sirve: fija la escala con la que juzgar cualquier diferencia. Si el
# 96 % de la separacion entre dos estaciones es un offset constante, una brecha
# de 4 hPa no es un error del sensor ni de la interpolacion, es geografia. Lo
# anomalo es apartarse del offset habitual de ESA estacion.

cat("\n--- Descomposicion offset / ruido por variable ---\n")

# Estaciones minimas por dia para que la media diaria de Madrid sea comparable,
# y dias minimos por estacion para que su offset sea estable.
MIN_ESTACIONES_DIA <- 3L
MIN_DIAS_OFFSET <- 100L

dir_estructura <- file.path(dir_base, "estructura_espacial")
dir.create(dir_estructura, recursive = TRUE, showWarnings = FALSE)

dt_desc <- dt_diario_largo[!is.na(Valor)]

# Media de Madrid del dia, solo con los dias que tienen suficientes estaciones.
dt_desc[, N_EST_DIA := .N, by = .(Variable, FECHA)]
dt_desc <- dt_desc[N_EST_DIA >= MIN_ESTACIONES_DIA]
dt_desc[, MEDIA_DIA := mean(Valor), by = .(Variable, FECHA)]
dt_desc[, DESVIACION := Valor - MEDIA_DIA]

# Offset de cada estacion: su desviacion media a lo largo del ano, y cuanto
# oscila esa desviacion (= el ruido que no explica el offset).
offsets <- dt_desc[, .(
  n_dias = .N,
  offset = mean(DESVIACION),
  sd_offset = sd(DESVIACION)
), by = .(Variable, Label, ESTACION)][n_dias >= MIN_DIAS_OFFSET]

# Residuo de cada dia una vez descontado el offset de su estacion.
dt_desc <- merge(dt_desc, offsets[, .(Variable, ESTACION, offset)],
  by = c("Variable", "ESTACION")
)
dt_desc[, RESIDUO := DESVIACION - offset]

# Dispersion espacial observada: media cuadratica de la sd entre estaciones de
# cada dia. Se usa la media CUADRATICA, no la aritmetica, porque es la que
# resulta comparable con un RMSE de interpolacion.
sd_espacial <- dt_desc[, .(sd_dia = sd(Valor)), by = .(Variable, FECHA)][
  !is.na(sd_dia), .(sd_espacial = sqrt(mean(sd_dia^2))),
  by = Variable
]

descomposicion <- offsets[, .(
  n_estaciones = .N,
  # Componente estructural: cuanto se separan entre si los offsets.
  offset_fijo = sd(offset),
  # Componente meteorologica: oscilacion tipica alrededor del propio offset.
  ruido = sqrt(mean(sd_offset^2)),
  offset_min = min(offset),
  offset_max = max(offset),
  est_min = ESTACION[which.min(offset)],
  est_max = ESTACION[which.max(offset)]
), by = .(Variable, Label)]

descomposicion <- merge(descomposicion, sd_espacial, by = "Variable")

# Reparto de varianza entre las dos componentes. Se calcula sobre la suma de
# ambas (y no sobre sd_espacial) para que sea un porcentaje propio y no pueda
# pasarse de 100 por diferencias de denominador.
descomposicion[, pct_estructural := round(
  100 * offset_fijo^2 / (offset_fijo^2 + ruido^2), 1
)]
# Umbral de alarma realista: 3 veces el ruido residual. Por debajo de esto una
# diferencia es meteorologia normal; por encima, merece mirarse.
descomposicion[, umbral_alarma := round(3 * ruido, 2)]

cols_num <- c("sd_espacial", "offset_fijo", "ruido", "offset_min", "offset_max")
descomposicion[, (cols_num) := lapply(.SD, round, 2), .SDcols = cols_num]

setorder(descomposicion, -pct_estructural)
setcolorder(descomposicion, c(
  "Label", "n_estaciones", "sd_espacial", "offset_fijo", "ruido",
  "pct_estructural", "umbral_alarma", "offset_min", "est_min",
  "offset_max", "est_max"
))

print(descomposicion[, .(
  Variable = Label, n_est = n_estaciones, sd_espacial,
  offset_fijo, ruido, pct_estructural, umbral_alarma
)])

cat("\n  sd_espacial    : separacion tipica entre estaciones el mismo dia\n")
cat("  offset_fijo    : parte de esa separacion que es constante todo el ano\n")
cat("  ruido          : oscilacion real dia a dia alrededor del propio offset\n")
cat("  pct_estructural: % de la varianza espacial que explica el offset fijo\n")
cat("  umbral_alarma  : 3 x ruido; por debajo, diferencia normal\n")

fwrite(descomposicion, file.path(dir_estructura, "descomposicion_offset_ruido.csv"))
fwrite(offsets[order(Variable, -offset)], file.path(dir_estructura, "offsets_por_estacion.csv"))

# ==============================================================================
# 11. FIGURA — OFFSET DE CADA ESTACIÓN CON SU OSCILACIÓN
# ==============================================================================
# Un punto por estacion (su offset) y una barra (+-1 sd de ese offset). La
# lectura es inmediata: si los puntos estan muy separados y las barras son
# cortas, la variable esta dominada por la geografia; si los puntos se apinan
# y las barras son largas, manda el tiempo atmosferico.

cat("\n--- Generando figura de offsets por estacion ---\n")

# Truco estandar para ordenar dentro de cada faceta: se antepone la variable a
# la etiqueta, se ordena globalmente y luego se recorta el prefijo en el eje.
offsets_plot <- copy(offsets)
offsets_plot[, ETIQUETA := paste(Variable, ESTACION, sep = "___")]
setorder(offsets_plot, Variable, offset)
offsets_plot[, ETIQUETA := factor(ETIQUETA, levels = ETIQUETA)]

p_off <- ggplot(offsets_plot, aes(x = offset, y = ETIQUETA, color = Label)) +
  geom_vline(xintercept = 0, color = "gray45", linewidth = 0.4) +
  geom_linerange(
    aes(xmin = offset - sd_offset, xmax = offset + sd_offset),
    linewidth = 0.6, alpha = 0.55
  ) +
  geom_point(size = 1.9) +
  scale_color_manual(values = paleta_vars, guide = "none") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  facet_wrap(~Label, ncol = 3, scales = "free") +
  labs(
    title = "Desnivel permanente de cada estación respecto a la media de Madrid",
    subtitle = paste(
      "Punto: desviación media del año (offset)  · ",
      "Barra: ±1 sd de esa desviación (ruido diario)"
    ),
    x = "Desviación respecto a la media diaria de Madrid",
    y = NULL,
    caption = paste(
      "Barras cortas y puntos separados = la diferencia entre estaciones es",
      "estructural (cota, entorno), no meteorológica"
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 9),
    plot.caption = element_text(color = "gray55", size = 7.5),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray96", color = "gray80"),
    axis.text.y = element_text(size = 6.5),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.9, "lines")
  )

ggsave(file.path(dir_estructura, "offsets_por_estacion.png"),
  plot = p_off, width = 15, height = 11, dpi = 200, bg = "white"
)
cat("  ✓ offsets_por_estacion.png\n")

# ==============================================================================
# 12. FIGURA — LA DISPERSIÓN ANTES Y DESPUÉS DE QUITAR EL OFFSET
# ==============================================================================
# Misma magnitud en las dos columnas y misma escala por fila, asi que el ancho
# de la banda es directamente comparable. Si al pasar a la derecha la banda se
# estrecha mucho, la separacion entre estaciones era el offset; si apenas
# cambia, era meteorologia.

cat("\n--- Generando figura antes/despues del offset ---\n")

etq_antes <- "Desviación respecto a la media diaria"
etq_despues <- "Tras restar el offset de la estación"

dt_comp <- melt(
  dt_desc,
  id.vars = c("ESTACION", "FECHA", "Label"),
  measure.vars = c("DESVIACION", "RESIDUO"),
  variable.name = "Tipo", value.name = "V"
)
dt_comp[, Tipo := factor(
  fifelse(Tipo == "DESVIACION", etq_antes, etq_despues),
  levels = c(etq_antes, etq_despues)
)]

p_anom <- ggplot(dt_comp, aes(x = FECHA, y = V, group = ESTACION, color = Label)) +
  geom_hline(yintercept = 0, color = "gray45", linewidth = 0.35) +
  geom_line(linewidth = 0.22, alpha = 0.35) +
  scale_color_manual(values = paleta_vars, guide = "none") +
  scale_x_date(
    breaks = primeros_meses, labels = etiquetas_meses,
    minor_breaks = NULL, expand = expansion(mult = c(0.01, 0.01))
  ) +
  facet_grid(Label ~ Tipo, scales = "free_y", switch = "y") +
  labs(
    title = "Cuánto de la separación entre estaciones es un desnivel fijo",
    subtitle = "Una línea por estación  ·  Madrid 2025  ·  misma escala en las dos columnas",
    x = NULL, y = NULL,
    caption = paste(
      "El estrechamiento de la banda al pasar a la derecha es la parte de la",
      "diferencia que explica la posición de la estación"
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 9),
    plot.caption = element_text(color = "gray55", size = 7.5),
    strip.text.x = element_text(face = "bold", size = 9),
    strip.text.y.left = element_text(face = "bold", size = 7.5, angle = 90),
    strip.background = element_rect(fill = "gray96", color = "gray80"),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.7, "lines")
  )

ggsave(file.path(dir_estructura, "dispersion_antes_despues_offset.png"),
  plot = p_anom, width = 14, height = 13, dpi = 200, bg = "white"
)
cat("  ✓ dispersion_antes_despues_offset.png\n")

# ==============================================================================
# 13. RESUMEN FINAL
# ==============================================================================

n_d <- length(list.files(dir_diario, pattern = "\\.png$"))
n_h <- length(list.files(dir_horario, pattern = "\\.png$"))
n_p <- length(list.files(dir_paneles, pattern = "\\.png$"))

n_e <- length(list.files(dir_estructura, pattern = "\\.png$"))

cat(sprintf("\n\u2705 Script completado.\n"))
cat(sprintf("   Salida en: %s\n", dir_base))
cat(sprintf("   Diario     : %d PNGs en diario/\n", n_d))
cat(sprintf("   Horario    : %d PNGs en horario/\n", n_h))
cat(sprintf("   Paneles    : %d PNGs en paneles/\n", n_p))
cat(sprintf("   Estructura : %d PNGs + 2 CSV en estructura_espacial/\n", n_e))
cat(sprintf("   TOTAL      : %d archivos PNG\n", n_d + n_h + n_p + n_e))
