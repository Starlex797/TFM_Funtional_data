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
  f <- here("data", "processed", "Maestro", "diario",
            sprintf("dataset_maestro_inla_%d_DIARIO.rds", a))
  if (!file.exists(f)) return(NULL)
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
  dl <- melt(d, id.vars = c("Etiqueta", "periodo"),
             measure.vars = c("NO2", "Trafico"),
             variable.name = "Var", value.name = "valor")
  dl[, Var := factor(Var, levels = c("NO2", "Trafico"),
                     labels = c("NO₂ (µg/m³)", "Tráfico (veh/h)"))]
  dl[, Etiqueta := factor(Etiqueta, levels = NIVELES)]

  ggplot(dl, aes(periodo, valor, color = Etiqueta)) +
    annotate("rect", xmin = COVID_INI, xmax = COVID_FIN,
             ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.15) +
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
       plot = p_semana, width = 12, height = 7, dpi = 200, bg = "white")
ggsave(file.path(DIR_SALIDA, "no2_trafico_mensual_estaciones.png"),
       plot = p_mes, width = 12, height = 7, dpi = 200, bg = "white")
fwrite(agg_semana, file.path(DIR_SALIDA, "no2_trafico_semanal_estaciones.csv"))
fwrite(agg_mes, file.path(DIR_SALIDA, "no2_trafico_mensual_estaciones.csv"))

cat("Estaciones usadas:\n"); print(SELECCION)
cat("\nGráficas guardadas en:\n", DIR_SALIDA, "\n")
