# ==============================================================================
# SCATTER PLOTS: NO₂ vs COVARIABLES CLIMÁTICAS E INTENSIDAD DE TRÁFICO
# Madrid 2019-2025 · Escala diaria · Un único PNG con todos los paneles
# Output: outputs/figures/eda/no2_cov/
# ==============================================================================
# QUÉ HACE:
#
# - Carga el maestro diario conjunto 2019-2025
#   (data/processed/Maestro/diario/dataset_maestro_inla_2019_2025_DIARIO.rds).
# - Elimina las observaciones sin una concentración válida de NO2.
# - Representa el NO2 diario frente a cada covariable en su escala original
#   (columnas *_raw):
#     · Temperatura.
#     · Humedad relativa.
#     · Precipitaciones.
#     · Presión barométrica.
#     · Radiación solar.
#     · Velocidad del viento.
#     · Intensidad del tráfico.
# - Todos los paneles comparten el eje Y (NO2) y cada uno tiene su propio eje X.
# - Los puntos se colorean por tipología de estación.
# - Etiquetas en inglés y pie de figura explicativo para el TFM.
# - Sin línea de suavizado (LOESS) ni coeficientes de correlación.
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CONFIGURACIÓN
# ==============================================================================

NOMBRE_ARCHIVO_DATOS <- "dataset_maestro_inla_2019_2025_DIARIO.rds"
RUTA_DATOS <- here("data", "processed", "Maestro", "diario", NOMBRE_ARCHIVO_DATOS)

DIR_SALIDA <- here("outputs", "figures", "eda", "no2_cov")
ARCHIVO_PNG <- file.path(DIR_SALIDA, "scatter_no2_covariates_daily_2019_2025.png")

PUNTO_SIZE <- 0.35
PUNTO_ALPHA <- 0.15

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

if (!file.exists(RUTA_DATOS)) {
  stop("No se encuentra el maestro: ", RUTA_DATOS)
}

cat("Cargando maestro DIARIO 2019-2025...\n")
dt <- as.data.table(readRDS(RUTA_DATOS))

# Unifica nombres con tildes o espacios (p. ej. "Presión Barométrica_raw" ->
# "Presion_Barometrica_raw") para que coincidan con los usados más abajo.
normalizar <- function(x) iconv(x, to = "ASCII//TRANSLIT")
nombres_limpios <- gsub("\\s+", "_", normalizar(names(dt)))
setnames(dt, names(dt), make.unique(nombres_limpios))

dt[, FECHA := as.Date(FECHA)]
dt <- dt[!is.na(DATO_DIARIO)]
cat(sprintf(
  "  %d filas · %d estaciones · %s a %s\n",
  nrow(dt), uniqueN(dt$ESTACION), format(min(dt$FECHA)), format(max(dt$FECHA))
))

# ==============================================================================
# 3. COVARIABLES Y ETIQUETAS EN INGLÉS
# ==============================================================================

# Nombre de columna -> etiqueta del panel (orden de aparición en la figura)
vars_cov <- c(
  Temperatura_raw         = "Temperature (°C)",
  Humedad_Relativa_raw    = "Relative humidity (%)",
  Precipitaciones_raw     = "Precipitation (mm)",
  Presion_Barometrica_raw = "Barometric pressure (mbar)",
  Radiacion_Solar_raw     = "Solar radiation (W/m²)",
  Velocidad_Viento_raw    = "Wind speed (m/s)",
  intensidad_raw          = "Traffic intensity (veh/h)"
)

faltan <- setdiff(names(vars_cov), names(dt))
if (length(faltan)) {
  cat("  [AVISO] Covariables ausentes, se omiten:", paste(faltan, collapse = ", "), "\n")
  vars_cov <- vars_cov[setdiff(names(vars_cov), faltan)]
}
if (!length(vars_cov)) stop("Ninguna covariable encontrada en el maestro.")
if (!"NOM_TIPO" %in% names(dt)) {
  stop("El maestro no tiene la columna NOM_TIPO (tipología de estación).")
}

# Tipología de estación en inglés (se normalizan tildes para el emparejamiento)
tipo_en <- c(
  "urbana trafico" = "Urban traffic",
  "urbana fondo"   = "Urban background",
  "suburbana"      = "Suburban"
)
dt[, station_type := tipo_en[tolower(normalizar(NOM_TIPO))]]
dt[is.na(station_type), station_type := as.character(NOM_TIPO)]
dt[, station_type := factor(station_type, levels = unique(c(tipo_en, station_type)))]
dt[, station_type := droplevels(station_type)]

paleta_tipo <- c(
  "Urban traffic"    = "#d73027",
  "Urban background" = "#4575b4",
  "Suburban"         = "#1a9850"
)

# ==============================================================================
# 4. FORMATO LARGO
# ==============================================================================

dt_long <- melt(
  dt[, c("ESTACION", "FECHA", "station_type", "DATO_DIARIO", names(vars_cov)), with = FALSE],
  id.vars = c("ESTACION", "FECHA", "station_type", "DATO_DIARIO"),
  measure.vars = names(vars_cov),
  variable.name = "covariable",
  value.name = "valor"
)
dt_long <- dt_long[!is.na(valor)]
dt_long[, panel := factor(vars_cov[as.character(covariable)], levels = vars_cov)]

# Orden aleatorio de filas: evita que una tipología tape sistemáticamente a otra
set.seed(2019)
dt_long <- dt_long[sample(.N)]

n_por_panel <- dt_long[, .N, by = panel][order(panel)]
cat("\nObservaciones válidas por panel:\n")
print(n_por_panel)

# ==============================================================================
# 5. FIGURA
# ==============================================================================

pie_figura <- paste(
  "Figure: Daily mean NO₂ concentration (µg/m³) plotted against each",
  "meteorological covariate and traffic intensity at the Madrid air-quality",
  sprintf(
    "monitoring network, %s–%s (%d stations).",
    format(min(dt$FECHA), "%Y"), format(max(dt$FECHA), "%Y"), uniqueN(dt$ESTACION)
  ),
  "Each point is one station-day; colours indicate the station type (urban",
  "traffic, urban background, suburban). All panels share the vertical axis,",
  "while each horizontal axis is shown in the original units of its covariate.",
  "No smoothing curve or correlation coefficient is added, so the panels show",
  "the raw marginal relationship, its dispersion, possible non-linearities and",
  "extreme values before any modelling."
)

p <- ggplot(dt_long, aes(x = valor, y = DATO_DIARIO, color = station_type)) +
  geom_point(size = PUNTO_SIZE, alpha = PUNTO_ALPHA, stroke = 0) +
  facet_wrap(~panel, ncol = 2, scales = "free_x", strip.position = "bottom") +
  scale_color_manual(values = paleta_tipo, name = "Station type", drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.03)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  labs(
    x = NULL,
    y = "Daily NO₂ (µg/m³)",
    caption = paste(strwrap(pie_figura, width = 115), collapse = "\n")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_minimal(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 10, margin = margin(t = 2, b = 6)),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80"),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8.5),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9.5),
    plot.caption = element_text(
      hjust = 0, size = 9, color = "gray20", lineheight = 1.15,
      margin = margin(t = 12)
    ),
    plot.caption.position = "plot",
    plot.margin = margin(10, 14, 10, 10)
  )

# ==============================================================================
# 6. GUARDADO
# ==============================================================================

dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
ggsave(ARCHIVO_PNG, plot = p, width = 8.5, height = 11, dpi = 300, bg = "white")

cat("==============================================================\n")
cat(sprintf("  Covariables representadas : %d\n", length(vars_cov)))
cat(sprintf("  PNG -> %s\n", ARCHIVO_PNG))
cat("==============================================================\n")
