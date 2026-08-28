# ==============================================================================
# Eventos aislados de viento y de lluvia: NO2, precipitación y viento HORARIOS
# ==============================================================================
# Para separar los efectos, se buscan automáticamente:
#   - Un día VENTOSO SIN lluvia  -> aísla el efecto del viento (dispersión).
#   - Un día de LLUVIA SIN viento -> aísla el efecto de la lluvia (lavado).
# Cada evento se muestra en una ventana de 5 días (día ± 2), a escala horaria,
# con el día del evento marcado.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

ANIO <- 2025L
LADO <- 2L # días a cada lado -> ventana de 5 días
DIR_SALIDA <- here("outputs", "figures", "no2_clima_2019_2025")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

a_hora <- function(x) as.integer(gsub("[^0-9]", "", as.character(x)))

# ------------------------------------------------------------------------------
# 1. NO2 + meteo horarios (media de Madrid por FECHA+HORA)
# ------------------------------------------------------------------------------

no2 <- as.data.table(readRDS(
  here("data", "processed", "Contaminacion", "horario",
       sprintf("aire_madrid_%d_No2_horarios.rds", ANIO))
))
no2[, HORA := a_hora(HORA)]
no2_h <- no2[, .(NO2 = mean(DATO, na.rm = TRUE)), by = .(FECHA, HORA)]

met <- as.data.table(readRDS(
  here("data", "processed", "Clima", "horario",
       sprintf("meteo_madrid_%d_horario.rds", ANIO))
))
met[, HORA := a_hora(HORA)]
met_h <- met[, .(
  precip = mean(Precipitaciones, na.rm = TRUE),
  viento = mean(`Velocidad Viento`, na.rm = TRUE)
), by = .(FECHA, HORA)]

dth <- merge(no2_h, met_h, by = c("FECHA", "HORA"))
dth <- dth[FECHA >= as.Date(sprintf("%d-01-01", ANIO)) &
             FECHA <= as.Date(sprintf("%d-12-31", ANIO))]
dth[, dt := as.POSIXct(FECHA, tz = "UTC") + (HORA - 1L) * 3600]

# ------------------------------------------------------------------------------
# 2. Selección del tramo de 5 días con más lluvia y con más viento
# ------------------------------------------------------------------------------

dias <- dth[, .(
  precip = sum(precip, na.rm = TRUE), viento = mean(viento, na.rm = TRUE)
), by = FECHA][order(FECHA)]

VIENTO_ALTO <- quantile(dias$viento, 0.90, na.rm = TRUE)
LLUVIA_MIN <- 3 # mm totales/día para marcar un día como lluvioso

# Tramos de 5 días consecutivos, AISLADOS (para separar los efectos):
#   - Lluvioso pero con poco viento (viento del tramo bajo la mediana).
#   - Ventoso pero seco (lluvia del tramo casi nula).
dias[, precip5 := frollsum(precip, 5, align = "center")]
dias[, viento5 := frollmean(viento, 5, align = "center")]
VIENTO5_MED <- median(dias$viento5, na.rm = TRUE)
SECO5 <- 2 # mm en 5 días ~ tramo prácticamente seco

# Tramo lluvioso y calmado: más lluvia entre los tramos de viento flojo.
p5_calmo <- dias$precip5
p5_calmo[is.na(dias$viento5) | dias$viento5 > VIENTO5_MED] <- NA
i_lluvia <- which.max(p5_calmo)
# Tramo ventoso y seco: más viento entre los tramos secos.
v5_seco <- dias$viento5
v5_seco[is.na(dias$precip5) | dias$precip5 > SECO5] <- NA
i_viento <- which.max(v5_seco)

win_lluvia <- dias[(i_lluvia - 2L):(i_lluvia + 2L)]
win_viento <- dias[(i_viento - 2L):(i_viento + 2L)]
dias_lluviosos <- win_lluvia[precip >= LLUVIA_MIN, FECHA]
dias_ventosos <- win_viento[viento >= VIENTO_ALTO, FECHA]

cat(sprintf(
  "Tramo LLUVIOSO y calmado: %s a %s | lluvia total %.1f mm, viento medio %.2f m/s\n",
  min(win_lluvia$FECHA), max(win_lluvia$FECHA),
  sum(win_lluvia$precip), mean(win_lluvia$viento)
))
cat(sprintf(
  "Tramo VENTOSO y seco: %s a %s | viento medio %.2f m/s, lluvia total %.1f mm\n",
  min(win_viento$FECHA), max(win_viento$FECHA),
  mean(win_viento$viento), sum(win_viento$precip)
))

# ------------------------------------------------------------------------------
# 3. Función que dibuja un tramo de 5 días marcando los días de evento
# ------------------------------------------------------------------------------

grafica_evento <- function(win_dias, dias_ev, color_ev, etiqueta_ev, sufijo) {
  ini <- min(win_dias$FECHA)
  fin <- max(win_dias$FECHA)
  dwin <- dth[FECHA >= ini & FECHA <= fin]
  ev <- data.table(
    xmin = as.POSIXct(dias_ev, tz = "UTC"),
    xmax = as.POSIXct(dias_ev, tz = "UTC") + 86400
  )
  capa <- geom_rect(data = ev, aes(xmin = xmin, xmax = xmax,
                                   ymin = -Inf, ymax = Inf),
                    fill = color_ev, alpha = 0.18, inherit.aes = FALSE)
  eje <- scale_x_datetime(date_breaks = "1 day", date_labels = "%d-%b")
  base <- theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())

  p1 <- ggplot(dwin, aes(dt, NO2)) + capa +
    geom_line(color = "#b2182b", linewidth = 0.6) + eje +
    labs(title = "NO₂ horario", x = NULL, y = "µg/m³") + base
  p2 <- ggplot(dwin, aes(dt, precip)) + capa +
    geom_col(fill = "#2166ac", width = 3000) + eje +
    labs(title = "Precipitación horaria", x = NULL, y = "mm") + base
  p3 <- ggplot(dwin, aes(dt, viento)) + capa +
    geom_line(color = "#1b7837", linewidth = 0.6) + eje +
    labs(title = "Viento horario", x = NULL, y = "m/s") + base

  fig <- (p1 / p2 / p3) + plot_layout(heights = c(1.2, 1, 1)) +
    plot_annotation(
      title = sprintf("NO₂, lluvia y viento a escala horaria — %s a %s",
                      format(ini, "%d %b"), format(fin, "%d %b %Y")),
      subtitle = sprintf("Franja sombreada = %s. Media de Madrid.", etiqueta_ev),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 10, color = "grey30"))
    )
  ruta <- file.path(DIR_SALIDA, sprintf("no2_evento_%s_horario.png", sufijo))
  ggsave(ruta, plot = fig, width = 11, height = 7.5, dpi = 200, bg = "white")
  cat("Guardado:", ruta, "\n")
}

grafica_evento(win_lluvia, dias_lluviosos, "#4575b4",
               "días de lluvia (tramo con poco viento)", "lluvia")
grafica_evento(win_viento, dias_ventosos, "#1b7837",
               "días ventosos (tramo seco)", "viento")
