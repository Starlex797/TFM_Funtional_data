# ==============================================================================
# NO2 y tráfico 2019-2025 por tipología de estación (semanal y mensual)
# ==============================================================================
# Relaciona la evolución del NO2 con la intensidad de tráfico (por barrio, la del
# preprocesamiento inicial) a lo largo de 2019-2025, distinguiendo la tipología
# de estación: tráfico, fondo y suburbana.
#   - Semanal: para ver con nitidez el confinamiento de 2020.
#   - Mensual: para ver la tendencia de fondo y la estacionalidad.
# Gráfica: panel doble (NO2 arriba, tráfico abajo) compartiendo el eje temporal,
# una línea por tipología, con el confinamiento sombreado.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

DIR_SALIDA <- here("outputs", "figures", "no2_trafico_2019_2025")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
ANIOS <- 2019:2025

# Estado de alarma / confinamiento estricto y desescalada (marzo-junio 2020).
COVID_INI <- as.Date("2020-03-14")
COVID_FIN <- as.Date("2020-06-21")

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
# 2. Tipología de estación (desde el NO2 crudo) y etiquetas
# ------------------------------------------------------------------------------

no2_raw <- as.data.table(readRDS(
  here("data", "processed", "Contaminacion", "horario",
       "aire_madrid_2025_No2_horarios.rds")
))
tipo <- unique(no2_raw[, .(ESTACION, NOM_TIPO)])
tipo[, Tipo := fifelse(
  grepl("tráfico|trafico", NOM_TIPO, ignore.case = TRUE), "Tráfico",
  fifelse(grepl("fondo", NOM_TIPO, ignore.case = TRUE), "Fondo", "Suburbana")
)]

dt <- merge(dt, tipo[, .(ESTACION, Tipo)], by = "ESTACION", all.x = TRUE)
dt <- dt[!is.na(Tipo)]
cat("Estaciones por tipología:\n")
print(unique(dt[, .(ESTACION, Tipo)])[, .N, by = Tipo])

# ------------------------------------------------------------------------------
# 3. Agregación (media por estación y luego media por tipología)
# ------------------------------------------------------------------------------

agregar <- function(dt, col_periodo) {
  por_est <- dt[, .(
    NO2 = mean(NO2, na.rm = TRUE),
    Trafico = mean(Trafico, na.rm = TRUE)
  ), by = c("ESTACION", "Tipo", col_periodo)]
  out <- por_est[, .(
    NO2 = mean(NO2, na.rm = TRUE),
    Trafico = mean(Trafico, na.rm = TRUE)
  ), by = c("Tipo", col_periodo)]
  setnames(out, col_periodo, "periodo")
  out[]
}

dt[, semana := as.Date(cut(FECHA, breaks = "week"))]
dt[, mes := as.Date(format(FECHA, "%Y-%m-01"))]

agg_semana <- agregar(dt, "semana")
agg_mes <- agregar(dt, "mes")

# ------------------------------------------------------------------------------
# 4. Gráfica: panel doble NO2 / tráfico por tipología
# ------------------------------------------------------------------------------

PALETA <- c("Tráfico" = "#d73027", "Fondo" = "#fc8d59", "Suburbana" = "#4575b4")

grafica_doble <- function(d, subtitulo) {
  dl <- melt(d, id.vars = c("Tipo", "periodo"),
             measure.vars = c("NO2", "Trafico"),
             variable.name = "Var", value.name = "valor")
  dl[, Var := factor(Var, levels = c("NO2", "Trafico"),
                     labels = c("NO₂ (µg/m³)", "Tráfico (veh/h)"))]

  ggplot(dl, aes(periodo, valor, color = Tipo)) +
    annotate("rect", xmin = COVID_INI, xmax = COVID_FIN,
             ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.15) +
    geom_line(linewidth = 0.55) +
    facet_grid(Var ~ ., scales = "free_y", switch = "y") +
    scale_color_manual(values = PALETA) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title = "NO₂ y tráfico por tipología de estación (2019-2025)",
      subtitle = subtitulo,
      x = NULL, y = NULL, color = "Tipología",
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

p_semana <- grafica_doble(agg_semana, "Media semanal")
p_mes <- grafica_doble(agg_mes, "Media mensual")

ggsave(file.path(DIR_SALIDA, "no2_trafico_semanal_tipologia.png"),
       plot = p_semana, width = 12, height = 7, dpi = 200, bg = "white")
ggsave(file.path(DIR_SALIDA, "no2_trafico_mensual_tipologia.png"),
       plot = p_mes, width = 12, height = 7, dpi = 200, bg = "white")

fwrite(agg_semana, file.path(DIR_SALIDA, "no2_trafico_semanal_tipologia.csv"))
fwrite(agg_mes, file.path(DIR_SALIDA, "no2_trafico_mensual_tipologia.csv"))

cat("\nGráficas guardadas en:\n", DIR_SALIDA, "\n")
