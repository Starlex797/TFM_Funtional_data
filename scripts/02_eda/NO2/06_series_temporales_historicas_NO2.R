# ==============================================================================
# SERIES TEMPORALES HISTÓRICAS DE NO₂ — 2019 A 2025
# Madrid · Escala diaria, mensual y horaria · Por estación
# Outputs: outputs/analysis/variogramas/series_temporales/
# ==============================================================================
# QUÉ HACE:
#
# - Carga y combina los datos de NO2 correspondientes al periodo 2019–2025.
# - Trabaja con tres resoluciones temporales:
#     · Horaria.
#     · Diaria.
#     · Mensual.
# - Añade el año de procedencia de cada observación.
# - Construye fechas continuas para las series mensuales y horarias.
# - Genera para cada estación:
#     · Una serie histórica diaria.
#     · Una serie histórica mensual.
#     · Una serie histórica horaria.
# - Señala visualmente los comienzos de año y de mes.
# - Incluye una línea horizontal de referencia en 40 µg/m³.
# - Calcula la cobertura temporal disponible para cada estación.
# - Genera paneles conjuntos con todas las estaciones para las tres escalas.
#
# FINALIDAD PARA EL TFM:
#
# Este script proporciona el contexto histórico del NO2 en Madrid. Permite
# estudiar tendencias a largo plazo, estacionalidad, periodos anómalos y
# diferencias persistentes entre estaciones.
#
# También permite comprobar si el comportamiento observado durante 2025 es
# representativo de los años anteriores o presenta características particulares.
#
# Puede utilizarse para analizar:
#
# - Tendencias temporales entre 2019 y 2025.
# - Posibles efectos asociados al periodo de la pandemia.
# - Estacionalidad recurrente del NO2.
# - Años o episodios especialmente contaminados.
# - Diferencias estructurales entre estaciones.
# - Cobertura y continuidad de las series históricas.
#
# SALIDAS:
#
# outputs/analysis/historicos_No2/
#     · diario/
#     · mensual/
#     · horario/
#     · panel_diario_todas_estaciones_2019_2025.png
#     · panel_mensual_todas_estaciones_2019_2025.png
#     · panel_horario_todas_estaciones_2019_2025.png


library(data.table)
library(ggplot2)
library(gt)
library(here)

ANIOS <- 2019:2025

# ==============================================================================
# 1. CARGA Y APILADO DE DATOS (7 años)
# ==============================================================================

cargar_anios <- function(escala, patron_archivo, patron_cols) {
  rbindlist(
    lapply(ANIOS, function(a) {
      ruta <- here(
        "data", "processed", "contaminacion", escala,
        sprintf(patron_archivo, a)
      )
      if (!file.exists(ruta)) {
        message("  No encontrado: ", basename(ruta))
        return(NULL)
      }
      dt <- readRDS(ruta)
      dt[, ANO := as.integer(a)]
      dt
    }),
    fill = TRUE
  )
}

cat("Cargando datos diarios...\n")
dt_diario <- cargar_anios("diario", "aire_madrid_%d_No2_trans_diarios1.rds")

cat("Cargando datos mensuales...\n")
dt_mensual <- cargar_anios("mensual", "aire_madrid_%d_log_No2_mensuales1.rds")

cat("Cargando datos horarios...\n")
dt_horario <- cargar_anios("horario", "aire_madrid_%d_No2_horarios.rds")

# Convertir MES "YYYY-MM" → Date (día 1 de cada mes) para el eje temporal
dt_mensual[, FECHA_MES := as.Date(paste0(MES, "-01"))]

# Construir DATETIME continuo para la serie horaria
dt_horario[, Hora_num := as.integer(sub("H0?", "", as.character(HORA)))]
dt_horario[, DATETIME := as.POSIXct(FECHA) + (Hora_num - 1L) * 3600L]

# Ordenar
setorder(dt_diario, ESTACION, FECHA)
setorder(dt_mensual, ESTACION, FECHA_MES)
setorder(dt_horario, ESTACION, DATETIME)

cat(sprintf(
  "  Diario  : %d obs · %d estaciones\n",
  nrow(dt_diario), uniqueN(dt_diario$ESTACION)
))
cat(sprintf(
  "  Mensual : %d obs · %d estaciones\n",
  nrow(dt_mensual), uniqueN(dt_mensual$ESTACION)
))
cat(sprintf(
  "  Horario : %d obs · %d estaciones\n",
  nrow(dt_horario), uniqueN(dt_horario$ESTACION)
))

# ==============================================================================
# 2. DELIMITADORES TEMPORALES — AÑOS Y MESES
# ==============================================================================

# Primero de cada año (delimitador principal)
primeros_anios <- seq(as.Date("2019-01-01"), as.Date("2026-01-01"), by = "1 year")
primeros_anios_ct <- as.POSIXct(primeros_anios)
etiquetas_anios <- format(primeros_anios, "%Y")

# Primero de cada mes (delimitador secundario) — para diario y mensual
primeros_meses <- seq(as.Date("2019-01-01"), as.Date("2025-12-01"), by = "1 month")
primeros_meses_ct <- as.POSIXct(primeros_meses)

# Franjas alternadas por AÑO (par = gris claro)
anios_limites <- primeros_anios
anios_limites_ct <- primeros_anios_ct

franjas_anios_d <- data.frame(
  xmin  = anios_limites[-length(anios_limites)],
  xmax  = anios_limites[-1],
  shade = (seq_along(ANIOS) %% 2L == 0L)
) |> subset(shade)

franjas_anios_ct <- data.frame(
  xmin  = anios_limites_ct[-length(anios_limites_ct)],
  xmax  = anios_limites_ct[-1],
  shade = (seq_along(ANIOS) %% 2L == 0L)
) |> subset(shade)

# ==============================================================================
# 3. DIRECTORIOS DE SALIDA
# ==============================================================================

dir_base <- here("outputs", "figures", "EDA", "NO2", "Series_historicas")
dir_diario <- file.path(dir_base, "diario")
dir_mensual <- file.path(dir_base, "mensual")
dir_horario <- file.path(dir_base, "horario")
dir_general <- file.path(dir_base, "general")

for (d in c(dir_diario, dir_mensual, dir_horario, dir_general)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 4. FUNCIÓN AUXILIAR
# ==============================================================================

limpiar_nombre <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", tolower(x))
}

# Tema base compartido
tema_ts <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "gray40", size = 9),
      plot.caption       = element_text(color = "gray55", size = 7.5),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

tema_panel <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = 14),
      plot.subtitle      = element_text(color = "gray40", size = 9),
      strip.text         = element_text(face = "bold", size = 7),
      strip.background   = element_rect(fill = "gray96", color = "gray80"),
      axis.text.x        = element_text(angle = 60, hjust = 1, size = 6),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(0.8, "lines")
    )
}

estaciones_unicas <- sort(unique(dt_diario$ESTACION))
cat(sprintf("\nEstaciones: %d\n\n", length(estaciones_unicas)))

# ==============================================================================
# 5. SERIES GENERALES DE LA RED - HORARIA, DIARIA Y MENSUAL
# ==============================================================================

# Cada observacion general es la media aritmetica de todas las estaciones con
# un valor valido en ese instante. N_ESTACIONES documenta su disponibilidad.
general_horario <- dt_horario[
  !is.na(DATO),
  .(
    NO2_GENERAL = mean(DATO),
    N_ESTACIONES = uniqueN(ESTACION)
  ),
  by = DATETIME
][order(DATETIME)]

general_diario <- dt_diario[
  !is.na(DATO_DIARIO),
  .(
    NO2_GENERAL = mean(DATO_DIARIO),
    N_ESTACIONES = uniqueN(ESTACION)
  ),
  by = FECHA
][order(FECHA)]

general_mensual <- dt_mensual[
  !is.na(DATO_MENSUAL),
  .(
    NO2_GENERAL = mean(DATO_MENSUAL),
    N_ESTACIONES = uniqueN(ESTACION)
  ),
  by = FECHA_MES
][order(FECHA_MES)]

if (!nrow(general_horario) || !nrow(general_diario) || !nrow(general_mensual)) {
  stop("No se pudieron calcular las series generales de NO2.")
}

rango_estaciones <- function(x) {
  paste(range(x$N_ESTACIONES), collapse = "-")
}

# Centros aproximados de las estaciones meteorologicas. Se muestran cuatro
# etiquetas por anio para evitar saturar el eje de una serie de siete anios.
centros_temporadas <- seq(
  as.Date("2019-01-15"), as.Date("2025-10-15"),
  by = "3 months"
)
centros_temporadas_ct <- as.POSIXct(centros_temporadas, tz = "UTC")
limites_temporadas <- seq(
  as.Date("2019-03-01"), as.Date("2025-12-01"),
  by = "3 months"
)
limites_temporadas_ct <- as.POSIXct(limites_temporadas, tz = "UTC")
etiquetas_temporadas <- rep(
  c("Winter", "Spring", "Summer", "Autumn"),
  length.out = length(centros_temporadas)
)
es_invierno <- etiquetas_temporadas == "Winter"
etiquetas_temporadas[es_invierno] <- paste0(
  etiquetas_temporadas[es_invierno], "\n",
  format(centros_temporadas[es_invierno], "%Y")
)

limites_anuales <- data.frame(
  Referencia = factor(
    c("current_40", "future_20"),
    levels = c("current_40", "future_20")
  ),
  Valor = c(40, 20)
)

confinamiento_d <- data.frame(
  xmin = as.Date("2020-03-14"),
  xmax = as.Date("2020-06-21")
)
confinamiento_ct <- data.frame(
  xmin = as.POSIXct("2020-03-14", tz = "UTC"),
  xmax = as.POSIXct("2020-06-21", tz = "UTC")
)

escalas_limites <- list(
  scale_color_manual(
    name = "Annual-mean references",
    values = c(
      "current_40" = "#b2182b",
      "future_20" = "#6a3d9a"
    ),
    labels = c(
      "current_40" = "Current annual limit: 40 ug/m3",
      "future_20" = "2030 annual limit: 20 ug/m3"
    )
  ),
  scale_linetype_manual(
    name = "Annual-mean references",
    values = c(
      "current_40" = "longdash",
      "future_20" = "dotdash"
    ),
    labels = c(
      "current_40" = "Current annual limit: 40 ug/m3",
      "future_20" = "2030 annual limit: 20 ug/m3"
    )
  )
)

cat("--- Generando series GENERALES de Madrid (2019-2025) ---\n")

p_general_h <- ggplot(general_horario, aes(DATETIME, NO2_GENERAL)) +
  geom_rect(
    data = franjas_anios_ct,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_rect(
    data = confinamiento_ct,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#fdae61", alpha = 0.28
  ) +
  geom_vline(
    xintercept = primeros_anios_ct,
    color = "gray45", linewidth = 0.45
  ) +
  geom_vline(
    xintercept = limites_temporadas_ct,
    color = "gray62", linewidth = 0.3, linetype = "dashed"
  ) +
  geom_hline(
    data = limites_anuales,
    aes(yintercept = Valor, color = Referencia, linetype = Referencia),
    linewidth = 0.65
  ) +
  geom_line(color = "#2c3e50", linewidth = 0.2, alpha = 0.65) +
  annotate(
    "text",
    x = as.POSIXct("2020-05-02", tz = "UTC"), y = Inf,
    label = "First state\nof alarm", vjust = 1.15,
    size = 2.4, fontface = "bold", color = "#8c510a"
  ) +
  scale_x_datetime(
    breaks = centros_temporadas_ct,
    labels = etiquetas_temporadas,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  escalas_limites +
  labs(
    title = "Network-wide hourly NO\u2082 concentration - Madrid, 2019-2025",
    subtitle = "Arithmetic mean across monitoring stations with valid data at each hour",
    x = NULL,
    y = "NO\u2082 (\u00b5g/m\u00b3)",
    caption = sprintf(
      paste(
        "Each observation is an hourly network mean (%s stations available).",
        "Dashed vertical lines delimit meteorological seasons; orange band: first state of alarm (14 Mar-21 Jun 2020).",
        "Horizontal lines are annual-mean references, not hourly exceedance thresholds.",
        sep = "\n"
      ),
      rango_estaciones(general_horario)
    )
  ) +
  tema_ts() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(dir_general, "serie_general_horaria_NO2_2019_2025.png"),
  plot = p_general_h, width = 16, height = 5.8, dpi = 200, bg = "white"
)
cat("  [OK] serie_general_horaria_NO2_2019_2025.png\n")

p_general_d <- ggplot(general_diario, aes(FECHA, NO2_GENERAL)) +
  geom_rect(
    data = franjas_anios_d,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_rect(
    data = confinamiento_d,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#fdae61", alpha = 0.28
  ) +
  geom_vline(
    xintercept = primeros_anios,
    color = "gray45", linewidth = 0.45
  ) +
  geom_vline(
    xintercept = limites_temporadas,
    color = "gray62", linewidth = 0.3, linetype = "dashed"
  ) +
  geom_hline(
    data = limites_anuales,
    aes(yintercept = Valor, color = Referencia, linetype = Referencia),
    linewidth = 0.65
  ) +
  geom_line(color = "#2166ac", linewidth = 0.4, alpha = 0.8) +
  annotate(
    "text",
    x = as.Date("2020-05-02"), y = Inf,
    label = "First state\nof alarm", vjust = 1.15,
    size = 2.4, fontface = "bold", color = "#8c510a"
  ) +
  scale_x_date(
    breaks = centros_temporadas,
    labels = etiquetas_temporadas,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  escalas_limites +
  labs(
    title = "Network-wide daily NO\u2082 concentration - Madrid, 2019-2025",
    subtitle = "Arithmetic mean across monitoring stations with valid data on each day",
    x = NULL,
    y = "NO\u2082 (\u00b5g/m\u00b3)",
    caption = sprintf(
      paste(
        "Each observation is a daily network mean (%s stations available).",
        "Dashed vertical lines delimit meteorological seasons; orange band: first state of alarm (14 Mar-21 Jun 2020).",
        "Horizontal lines are annual-mean references, not daily exceedance thresholds.",
        sep = "\n"
      ),
      rango_estaciones(general_diario)
    )
  ) +
  tema_ts() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(dir_general, "serie_general_diaria_NO2_2019_2025.png"),
  plot = p_general_d, width = 16, height = 5.8, dpi = 200, bg = "white"
)
cat("  [OK] serie_general_diaria_NO2_2019_2025.png\n")

p_general_m <- ggplot(general_mensual, aes(FECHA_MES, NO2_GENERAL)) +
  geom_rect(
    data = franjas_anios_d,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_rect(
    data = confinamiento_d,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#fdae61", alpha = 0.28
  ) +
  geom_vline(
    xintercept = primeros_anios,
    color = "gray45", linewidth = 0.45
  ) +
  geom_vline(
    xintercept = limites_temporadas,
    color = "gray62", linewidth = 0.3, linetype = "dashed"
  ) +
  geom_hline(
    data = limites_anuales,
    aes(yintercept = Valor, color = Referencia, linetype = Referencia),
    linewidth = 0.65
  ) +
  geom_area(fill = "#2ca25f", alpha = 0.16) +
  geom_line(color = "#238b45", linewidth = 0.75) +
  geom_point(
    color = "#238b45", fill = "white", shape = 21,
    size = 1.8, stroke = 0.8
  ) +
  annotate(
    "text",
    x = as.Date("2020-05-02"), y = Inf,
    label = "First state\nof alarm", vjust = 1.15,
    size = 2.4, fontface = "bold", color = "#8c510a"
  ) +
  scale_x_date(
    breaks = centros_temporadas,
    labels = etiquetas_temporadas,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  escalas_limites +
  labs(
    title = "Network-wide monthly NO\u2082 concentration - Madrid, 2019-2025",
    subtitle = "Arithmetic mean across monitoring stations with a valid monthly mean",
    x = NULL,
    y = "NO\u2082 (\u00b5g/m\u00b3)",
    caption = sprintf(
      paste(
        "Each observation is a monthly network mean (%s stations available).",
        "Dashed vertical lines delimit meteorological seasons; orange band: first state of alarm (14 Mar-21 Jun 2020).",
        "Horizontal lines are annual-mean references, not monthly exceedance thresholds.",
        sep = "\n"
      ),
      rango_estaciones(general_mensual)
    )
  ) +
  tema_ts() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(dir_general, "serie_general_mensual_NO2_2019_2025.png"),
  plot = p_general_m, width = 16, height = 5.8, dpi = 200, bg = "white"
)
cat("  [OK] serie_general_mensual_NO2_2019_2025.png\n")

# Tabla academica: primer estado de alarma frente al mismo periodo de 2019 ----
resumen_alarma <- rbindlist(lapply(c(2019L, 2020L), function(anio) {
  fecha_inicio <- as.Date(sprintf("%d-03-14", anio))
  fecha_fin <- as.Date(sprintf("%d-06-20", anio))
  x <- general_diario[FECHA >= fecha_inicio & FECHA <= fecha_fin]

  data.table(
    Year = anio,
    Days = nrow(x),
    Stations_min = min(x$N_ESTACIONES),
    Stations_max = max(x$N_ESTACIONES),
    Mean_NO2 = mean(x$NO2_GENERAL)
  )
}))

if (any(resumen_alarma$Days != 99L)) {
  stop("La comparacion del estado de alarma no contiene 99 dias en ambos anios.")
}

media_2019 <- resumen_alarma[Year == 2019L, Mean_NO2]
resumen_alarma[, Absolute_change := Mean_NO2 - media_2019]
resumen_alarma[, Relative_change := 100 * Absolute_change / media_2019]
resumen_alarma[Year == 2019L, `:=`(
  Absolute_change = NA_real_,
  Relative_change = NA_real_
)]
resumen_alarma[, Stations_per_day := paste0(Stations_min, "-", Stations_max)]

tabla_alarma_gt <- resumen_alarma[, .(
  Year, Days, Stations_per_day, Mean_NO2, Absolute_change, Relative_change
)] |>
  as.data.frame() |>
  gt() |>
  tab_header(
    title = html(
      "<strong>Table 3.2. Network-wide NO<sub>2</sub> during the first COVID-19 state of alarm</strong>"
    ),
    subtitle = html(
      "Madrid &middot; Equivalent period: 14 March&ndash;20 June, 2019 and 2020"
    )
  ) |>
  cols_label(
    Year = html("<strong>Year</strong>"),
    Days = html("<strong>Days</strong>"),
    Stations_per_day = html("<strong>Stations<br>per day</strong>"),
    Mean_NO2 = html("<strong>Mean NO<sub>2</sub><br>(&micro;g/m<sup>3</sup>)</strong>"),
    Absolute_change = html(
      "<strong>Absolute difference<br>(&micro;g/m<sup>3</sup>)</strong>"
    ),
    Relative_change = html(
      "<strong>Relative difference<br>versus 2019 (%)</strong>"
    )
  ) |>
  tab_spanner(
    label = html("<strong>Data availability</strong>"),
    columns = c(Days, Stations_per_day)
  ) |>
  tab_spanner(
    label = html("<strong>Network-wide concentration</strong>"),
    columns = c(Mean_NO2, Absolute_change, Relative_change)
  ) |>
  fmt_integer(columns = c(Year, Days), use_seps = FALSE) |>
  fmt_number(
    columns = c(Mean_NO2, Absolute_change),
    decimals = 2, dec_mark = ".", sep_mark = ","
  ) |>
  fmt_number(
    columns = Relative_change,
    decimals = 1, dec_mark = ".", sep_mark = ",",
    pattern = "{x}%"
  ) |>
  sub_missing(
    columns = c(Absolute_change, Relative_change),
    missing_text = "Reference"
  ) |>
  cols_align(align = "center", columns = c(Year, Days, Stations_per_day)) |>
  cols_align(
    align = "right",
    columns = c(Mean_NO2, Absolute_change, Relative_change)
  ) |>
  tab_source_note(
    source_note = html(
      paste0(
        "<strong>Source.</strong> Madrid City Council Air Quality Monitoring ",
        "Network; period defined from Spanish Royal Decree 463/2020."
      )
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold", size = px(16)),
    locations = cells_title(groups = "title")
  ) |>
  tab_style(
    style = cell_text(size = px(11), color = "#333333"),
    locations = cells_title(groups = "subtitle")
  ) |>
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_borders(sides = "bottom", color = "#333333", weight = px(1.5))
    ),
    locations = list(cells_column_spanners(), cells_column_labels())
  ) |>
  tab_style(
    style = cell_borders(sides = "bottom", color = "#D9D9D9", weight = px(0.5)),
    locations = cells_body(rows = 1)
  ) |>
  cols_width(
    Year ~ px(70),
    Days ~ px(70),
    Stations_per_day ~ px(105),
    Mean_NO2 ~ px(120),
    Absolute_change ~ px(145),
    Relative_change ~ px(155)
  ) |>
  tab_options(
    table.font.names = "Arial",
    table.font.size = 12,
    table.width = px(790),
    heading.align = "left",
    heading.border.bottom.color = "#111111",
    heading.border.bottom.width = px(2),
    table.border.top.color = "#111111",
    table.border.top.width = px(2),
    table.border.bottom.color = "#111111",
    table.border.bottom.width = px(2),
    data_row.padding = px(8),
    source_notes.font.size = 9,
    source_notes.padding = px(7)
  )

# Version PNG reproducible sin navegador, con el mismo lenguaje visual que
# Table 3.1: fondo blanco, reglas horizontales y ausencia de celdas coloreadas.
x_columnas <- c(0.065, 0.205, 0.335, 0.49, 0.68, 0.875)
fila_2019 <- c(
  "2019", "99", "22-24", "26.69", "Reference", "Reference"
)
fila_2020 <- c(
  "2020", "99", "22-24", "13.45", "-13.24", "-49.6%"
)

tabla_alarma_png <- ggplot() +
  xlim(0, 1) +
  ylim(0.35, 1) +
  annotate(
    "segment",
    x = 0.02, xend = 0.98, y = 0.975, yend = 0.975,
    linewidth = 0.8, color = "#111111"
  ) +
  annotate(
    "text",
    x = 0.025, y = 0.93,
    label = "bold('Table 3.2. Network-wide NO'[2]~'during the first COVID-19 state of alarm')",
    parse = TRUE,
    hjust = 0, fontface = "bold", family = "Arial", size = 5.6,
    color = "#1F1F1F"
  ) +
  annotate(
    "text",
    x = 0.025, y = 0.865,
    label = "Madrid, equivalent period: 14 March-20 June, 2019 and 2020",
    hjust = 0, family = "Arial", size = 3.8, color = "#1F1F1F"
  ) +
  annotate(
    "segment",
    x = 0.02, xend = 0.98, y = 0.825, yend = 0.825,
    linewidth = 0.8, color = "#111111"
  ) +
  annotate(
    "text",
    x = 0.27, y = 0.785, label = "Data availability",
    fontface = "bold", family = "Arial", size = 3.8, color = "#1F1F1F"
  ) +
  annotate(
    "text",
    x = 0.695, y = 0.785, label = "Network-wide concentration",
    fontface = "bold", family = "Arial", size = 3.8, color = "#1F1F1F"
  ) +
  annotate(
    "segment",
    x = 0.145, xend = 0.405, y = 0.758, yend = 0.758,
    linewidth = 0.55, color = "#C9C9C9"
  ) +
  annotate(
    "segment",
    x = 0.425, xend = 0.975, y = 0.758, yend = 0.758,
    linewidth = 0.55, color = "#C9C9C9"
  ) +
  annotate(
    "text",
    x = x_columnas[c(1, 2, 3, 6)], y = 0.715,
    label = c("Year", "Days", "Stations\nper day", "Relative difference\nversus 2019 (%)"),
    fontface = "bold", family = "Arial", size = 3.65, color = "#1F1F1F"
  ) +
  annotate(
    "text",
    x = x_columnas[4], y = 0.728,
    label = "bold('Mean NO'[2])", parse = TRUE,
    fontface = "bold", family = "Arial", size = 3.65, color = "#1F1F1F"
  ) +
  annotate(
    "text",
    x = x_columnas[4], y = 0.695,
    label = "bold(mu*'g/m'^3)", parse = TRUE,
    family = "Arial", size = 3.65, color = "#1F1F1F"
  ) +
  annotate(
    "text",
    x = x_columnas[5], y = 0.728,
    label = "Absolute difference",
    fontface = "bold", family = "Arial", size = 3.65, color = "#1F1F1F"
  ) +
  annotate(
    "text",
    x = x_columnas[5], y = 0.695,
    label = "bold(mu*'g/m'^3)", parse = TRUE,
    family = "Arial", size = 3.65, color = "#1F1F1F"
  ) +
  annotate(
    "segment",
    x = 0.02, xend = 0.98, y = 0.665, yend = 0.665,
    linewidth = 0.6, color = "#C9C9C9"
  ) +
  annotate(
    "text",
    x = rep(x_columnas, 2),
    y = rep(c(0.615, 0.555), each = length(x_columnas)),
    label = c(fila_2019, fila_2020),
    family = "Arial", size = 3.7, color = "#1F1F1F"
  ) +
  annotate(
    "segment",
    x = 0.02, xend = 0.98, y = 0.585, yend = 0.585,
    linewidth = 0.3, color = "#D9D9D9"
  ) +
  annotate(
    "segment",
    x = 0.02, xend = 0.98, y = 0.515, yend = 0.515,
    linewidth = 0.6, color = "#C9C9C9"
  ) +
  annotate(
    "text",
    x = 0.025, y = 0.455, label = "Source.",
    hjust = 0, fontface = "bold", family = "Arial", size = 2.8
  ) +
  annotate(
    "text",
    x = 0.078, y = 0.455,
    label = paste0(
      "Madrid City Council Air Quality Monitoring Network; period defined from Spanish ",
      "Royal Decree 463/2020."
    ),
    hjust = 0, family = "Arial", size = 2.8
  ) +
  annotate(
    "segment",
    x = 0.02, xend = 0.98, y = 0.405, yend = 0.405,
    linewidth = 0.8, color = "#111111"
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 8, 5, 8)
  )

ggsave(
  filename = file.path(dir_general, "table_NO2_first_state_alarm_2019_2020.png"),
  plot = tabla_alarma_png,
  width = 9.7, height = 3.25, dpi = 200, bg = "white"
)
cat("  [OK] table_NO2_first_state_alarm_2019_2020.png\n")

# ==============================================================================
# 6. SERIES DIARIAS — POR ESTACIÓN
# ==============================================================================

cat("--- Generando series DIARIAS por estación (2019-2025) ---\n")

for (est in estaciones_unicas) {
  dt_est <- dt_diario[ESTACION == est & !is.na(DATO_DIARIO)]
  if (nrow(dt_est) < 2L) {
    cat(sprintf("  [SKIP diario] %s\n", est))
    next
  }

  n_dias <- nrow(dt_est)
  pct_ok <- round(100 * n_dias / (365.25 * 7), 1)

  p <- ggplot(dt_est, aes(x = FECHA, y = DATO_DIARIO)) +

    # Franjas alternadas por año
    geom_rect(
      data = franjas_anios_d,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.5
    ) +

    # Líneas divisorias de año (gruesas)
    geom_vline(
      xintercept = as.numeric(primeros_anios),
      color = "gray40", linewidth = 0.6, linetype = "solid"
    ) +

    # Líneas divisorias de mes (finas, secundarias)
    geom_vline(
      xintercept = as.numeric(primeros_meses),
      color = "gray75", linewidth = 0.2, linetype = "dashed"
    ) +

    # Serie temporal
    geom_line(color = "#2c3e50", linewidth = 0.4, alpha = 0.8) +

    # Umbral OMS (40 µg/m³)
    geom_hline(
      yintercept = 40, color = "#e74c3c",
      linewidth = 0.65, linetype = "longdash", alpha = 0.8
    ) +
    annotate("text",
      x = as.Date("2019-02-01"), y = 42,
      label = "Límite OMS 40 µg/m³", color = "#e74c3c",
      size = 2.6, hjust = 0, fontface = "italic"
    ) +

    # Etiquetas de año sobre las líneas divisorias
    annotate("text",
      x = as.Date(paste0(ANIOS, "-07-01")),
      y = Inf,
      label = as.character(ANIOS),
      vjust = 1.5, size = 3.2, color = "gray35", fontface = "bold"
    ) +
    scale_x_date(
      breaks = primeros_meses,
      labels = function(d) {
        ifelse(format(d, "%m") == "01",
          format(d, "%b\n%Y"),
          format(d, "%b")
        )
      },
      minor_breaks = NULL,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    labs(
      title    = paste0("Serie Diaria de NO\u2082 \u2014 ", est),
      subtitle = sprintf("Madrid 2019\u20132025  |  %d d\u00edas v\u00e1lidos (%.1f%%)", n_dias, pct_ok),
      x        = NULL,
      y        = "NO\u2082 (\u00b5g/m\u00b3)",
      caption  = "Franjas grises: a\u00f1os pares  \u00b7  L\u00edneas verticales finas: inicio de mes  \u00b7  L\u00ednea roja: l\u00edmite OMS"
    ) +
    tema_ts()

  nombre_file <- file.path(
    dir_diario,
    paste0("diario_", limpiar_nombre(est), "_2019_2025.png")
  )
  ggsave(nombre_file, plot = p, width = 16, height = 5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 7. SERIES MENSUALES — POR ESTACIÓN
# ==============================================================================

cat("\n--- Generando series MENSUALES por estación (2019-2025) ---\n")

for (est in estaciones_unicas) {
  dt_est <- dt_mensual[ESTACION == est & !is.na(DATO_MENSUAL)]
  if (nrow(dt_est) < 2L) {
    cat(sprintf("  [SKIP mensual] %s\n", est))
    next
  }

  p <- ggplot(dt_est, aes(x = FECHA_MES, y = DATO_MENSUAL)) +

    # Franjas alternadas por año
    geom_rect(
      data = franjas_anios_d,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.5
    ) +

    # Líneas divisorias de año
    geom_vline(
      xintercept = as.numeric(primeros_anios),
      color = "gray40", linewidth = 0.6
    ) +

    # Líneas de mes (secundarias)
    geom_vline(
      xintercept = as.numeric(primeros_meses),
      color = "gray75", linewidth = 0.25, linetype = "dashed"
    ) +

    # Área bajo la curva
    geom_area(fill = "#2980b9", alpha = 0.15) +
    geom_line(color = "#2980b9", linewidth = 0.75) +
    geom_point(color = "#2980b9", size = 1.8, shape = 21, fill = "white", stroke = 1) +

    # Umbral OMS
    geom_hline(
      yintercept = 40, color = "#e74c3c",
      linewidth = 0.65, linetype = "longdash", alpha = 0.8
    ) +
    annotate("text",
      x = as.Date("2019-02-01"), y = 42,
      label = "Límite OMS 40 µg/m³", color = "#e74c3c",
      size = 2.6, hjust = 0, fontface = "italic"
    ) +

    # Etiquetas de año
    annotate("text",
      x = as.Date(paste0(ANIOS, "-07-01")),
      y = Inf,
      label = as.character(ANIOS),
      vjust = 1.5, size = 3.2, color = "gray35", fontface = "bold"
    ) +
    scale_x_date(
      breaks = primeros_meses,
      labels = function(d) {
        ifelse(format(d, "%m") == "01",
          format(d, "%b\n%Y"),
          format(d, "%b")
        )
      },
      minor_breaks = NULL,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      title    = paste0("Serie Mensual de NO\u2082 \u2014 ", est),
      subtitle = "Madrid 2019\u20132025  |  Media mensual de medias diarias v\u00e1lidas",
      x        = NULL,
      y        = "NO\u2082 (\u00b5g/m\u00b3)",
      caption  = "Franjas grises: a\u00f1os pares  \u00b7  L\u00ednea roja: l\u00edmite anual OMS (40 \u00b5g/m\u00b3)"
    ) +
    tema_ts()

  nombre_file <- file.path(
    dir_mensual,
    paste0("mensual_", limpiar_nombre(est), "_2019_2025.png")
  )
  ggsave(nombre_file, plot = p, width = 16, height = 5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 8. SERIES HORARIAS — POR ESTACIÓN
# ==============================================================================

cat("\n--- Generando series HORARIAS por estación (2019-2025) ---\n")

for (est in estaciones_unicas) {
  dt_est <- dt_horario[ESTACION == est & !is.na(DATO)]
  if (nrow(dt_est) < 2L) {
    cat(sprintf("  [SKIP horario] %s\n", est))
    next
  }

  n_horas <- nrow(dt_est)
  pct_ok <- round(100 * n_horas / (8760 * 7), 1)

  p <- ggplot(dt_est, aes(x = DATETIME, y = DATO)) +

    # Franjas alternadas por año
    geom_rect(
      data = franjas_anios_ct,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gray88", alpha = 0.5
    ) +

    # Líneas divisorias de año
    geom_vline(
      xintercept = primeros_anios_ct,
      color = "gray40", linewidth = 0.6
    ) +

    # Líneas de mes (secundarias)
    geom_vline(
      xintercept = primeros_meses_ct,
      color = "gray75", linewidth = 0.15, linetype = "dashed"
    ) +

    # Serie horaria
    geom_line(color = "#2c3e50", linewidth = 0.25, alpha = 0.65) +

    # Umbral OMS
    geom_hline(
      yintercept = 40, color = "#e74c3c",
      linewidth = 0.55, linetype = "longdash", alpha = 0.75
    ) +

    # Etiquetas de año
    annotate("text",
      x = as.POSIXct(paste0(ANIOS, "-07-01")),
      y = Inf,
      label = as.character(ANIOS),
      vjust = 1.5, size = 3.2, color = "gray35", fontface = "bold"
    ) +
    scale_x_datetime(
      breaks       = primeros_anios_ct,
      labels       = etiquetas_anios,
      minor_breaks = primeros_meses_ct,
      expand       = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    labs(
      title    = paste0("Serie Horaria de NO\u2082 \u2014 ", est),
      subtitle = sprintf("Madrid 2019\u20132025  |  %d horas v\u00e1lidas (%.1f%% del per\u00edodo)", n_horas, pct_ok),
      x        = NULL,
      y        = "NO\u2082 (\u00b5g/m\u00b3)",
      caption  = "Franjas grises: a\u00f1os pares  \u00b7  Ticks menores: inicio de mes  \u00b7  L\u00ednea roja: l\u00edmite OMS"
    ) +
    tema_ts()

  nombre_file <- file.path(
    dir_horario,
    paste0("horario_", limpiar_nombre(est), "_2019_2025.png")
  )
  ggsave(nombre_file, plot = p, width = 16, height = 5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 %s\n", basename(nombre_file)))
}

# ==============================================================================
# 9. PANELES RESUMEN — TODAS LAS ESTACIONES (facet)
# ==============================================================================

cat("\n--- Generando paneles resumen (todas las estaciones) ---\n")

# ── Panel diario ──────────────────────────────────────────────────────────────
p_panel_d <- ggplot(
  dt_diario[!is.na(DATO_DIARIO)],
  aes(x = FECHA, y = DATO_DIARIO)
) +
  geom_rect(
    data = franjas_anios_d,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_vline(
    xintercept = as.numeric(primeros_anios),
    color = "gray45", linewidth = 0.4
  ) +
  geom_line(color = "#2c3e50", linewidth = 0.3, alpha = 0.75) +
  geom_hline(
    yintercept = 40, color = "#e74c3c",
    linewidth = 0.4, linetype = "longdash", alpha = 0.7
  ) +
  scale_x_date(
    breaks = primeros_anios, labels = etiquetas_anios,
    minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))
  ) +
  facet_wrap(~ESTACION, ncol = 4, scales = "free_y") +
  labs(
    title = "Series Diarias de NO\u2082 \u2014 Todas las Estaciones",
    subtitle = "Madrid 2019\u20132025  \u00b7  L\u00ednea roja: l\u00edmite OMS (40 \u00b5g/m\u00b3)",
    x = NULL, y = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  tema_panel()

ggsave(file.path(dir_base, "panel_diario_todas_estaciones_2019_2025.png"),
  plot = p_panel_d, width = 22, height = 18, dpi = 200, bg = "white"
)
cat("  \u2713 panel_diario_todas_estaciones_2019_2025.png\n")

# ── Panel mensual ─────────────────────────────────────────────────────────────
p_panel_m <- ggplot(
  dt_mensual[!is.na(DATO_MENSUAL)],
  aes(x = FECHA_MES, y = DATO_MENSUAL)
) +
  geom_rect(
    data = franjas_anios_d,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_vline(
    xintercept = as.numeric(primeros_anios),
    color = "gray45", linewidth = 0.4
  ) +
  geom_area(fill = "#2980b9", alpha = 0.15) +
  geom_line(color = "#2980b9", linewidth = 0.5) +
  geom_point(
    color = "#2980b9", size = 0.9, shape = 21,
    fill = "white", stroke = 0.7
  ) +
  geom_hline(
    yintercept = 40, color = "#e74c3c",
    linewidth = 0.4, linetype = "longdash", alpha = 0.7
  ) +
  scale_x_date(
    breaks = primeros_anios, labels = etiquetas_anios,
    minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))
  ) +
  facet_wrap(~ESTACION, ncol = 4, scales = "free_y") +
  labs(
    title = "Series Mensuales de NO\u2082 \u2014 Todas las Estaciones",
    subtitle = "Madrid 2019\u20132025  \u00b7  L\u00ednea roja: l\u00edmite OMS (40 \u00b5g/m\u00b3)",
    x = NULL, y = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  tema_panel()

ggsave(file.path(dir_base, "panel_mensual_todas_estaciones_2019_2025.png"),
  plot = p_panel_m, width = 22, height = 18, dpi = 200, bg = "white"
)
cat("  \u2713 panel_mensual_todas_estaciones_2019_2025.png\n")

# ── Panel horario ─────────────────────────────────────────────────────────────
p_panel_h <- ggplot(
  dt_horario[!is.na(DATO)],
  aes(x = DATETIME, y = DATO)
) +
  geom_rect(
    data = franjas_anios_ct,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "gray88", alpha = 0.5
  ) +
  geom_vline(
    xintercept = primeros_anios_ct,
    color = "gray45", linewidth = 0.35
  ) +
  geom_line(color = "#2c3e50", linewidth = 0.2, alpha = 0.6) +
  geom_hline(
    yintercept = 40, color = "#e74c3c",
    linewidth = 0.35, linetype = "longdash", alpha = 0.7
  ) +
  scale_x_datetime(
    breaks = primeros_anios_ct, labels = etiquetas_anios,
    minor_breaks = NULL, expand = expansion(mult = c(0.005, 0.005))
  ) +
  facet_wrap(~ESTACION, ncol = 4, scales = "free_y") +
  labs(
    title = "Series Horarias de NO\u2082 \u2014 Todas las Estaciones",
    subtitle = "Madrid 2019\u20132025  \u00b7  L\u00ednea roja: l\u00edmite OMS (40 \u00b5g/m\u00b3)",
    x = NULL, y = "NO\u2082 (\u00b5g/m\u00b3)"
  ) +
  tema_panel()

ggsave(file.path(dir_base, "panel_horario_todas_estaciones_2019_2025.png"),
  plot = p_panel_h, width = 22, height = 18, dpi = 200, bg = "white"
)
cat("  \u2713 panel_horario_todas_estaciones_2019_2025.png\n")

# ==============================================================================
# 10. RESUMEN FINAL
# ==============================================================================

n_d <- length(list.files(dir_diario, pattern = "\\.png$"))
n_m <- length(list.files(dir_mensual, pattern = "\\.png$"))
n_h <- length(list.files(dir_horario, pattern = "\\.png$"))

cat(sprintf("\n\u2705 Script completado.\n"))
cat(sprintf("   Diario   : %d PNGs en outputs/analysis/variogramas/series_temporales/diario/\n", n_d))
cat(sprintf("   Mensual  : %d PNGs en outputs/analysis/variogramas/series_temporales/mensual/\n", n_m))
cat(sprintf("   Horario  : %d PNGs en outputs/analysis/variogramas/series_temporales/horario/\n", n_h))
cat(sprintf("   Paneles  : 3 PNGs en outputs/analysis/variogramas/series_temporales/\n"))
cat(sprintf("   TOTAL    : %d archivos PNG\n", n_d + n_m + n_h + 3L))
