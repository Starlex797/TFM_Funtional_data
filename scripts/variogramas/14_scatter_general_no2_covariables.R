# ==============================================================================
# SCATTER GENERAL: NO₂ vs COVARIABLES — SIN AGRUPACIÓN
# Todos los datos juntos · Diario y horario · Madrid 2025
# Outputs: outputs/scatter_no2_covariables/general/
# ==============================================================================
 QUÉ HACE:
#
# - Carga los datasets maestros diario y horario de 2025.
# - Elimina las filas sin una concentración válida de NO2.
# - Analiza las ocho covariables meteorológicas y de tráfico.
#
# - Agrupa todas las estaciones y todas las observaciones en una única nube de
#   puntos para cada covariable.
# - No diferencia por:
#     · Estación de contaminación.
#     · Mes.
#     · Estación del año.
#     · Tipo de día.
#
# - Genera dos gráficos para cada covariable:
#     · Scatter plot general diario.
#     · Scatter plot general horario.
#
# - Los datos diarios se representan en azul.
# - Los datos horarios se representan en rojo y con mayor transparencia debido
#   al elevado número de observaciones.
# - No añade líneas de tendencia.
#
# FINALIDAD PARA EL TFM:
#
# Este script proporciona la representación más general y directa de la relación
# entre el NO2 y cada covariable.
#
# Sirve como diagnóstico inicial para comprobar:
#
# - La forma global de la relación.
# - La presencia de valores extremos.
# - Posibles relaciones lineales o no lineales.
# - Heterocedasticidad.
# - Zonas con alta concentración de observaciones.
# - Diferencias visuales entre escala diaria y horaria.
#
# No obstante, al mezclar estaciones y periodos temporales, una relación global
# puede estar condicionada por diferencias espaciales o estacionales. Por ello,
# estos gráficos deben interpretarse junto con los scripts 11, 12 y 13.
#
# SALIDAS:
#
# outputs/scatter_no2_covariables/general/
#
# Se genera un scatter plot diario y otro horario para cada covariable.


# ==============================================================================
# SECUENCIA GENERAL DE ESTE BLOQUE
# ==============================================================================
#
# Relación NO2-covariables por estación, mes y estación del año
#     ->
# Evolución temporal detallada de estaciones seleccionadas
#     ->
# Relación conjunta diferenciada por estación y época del año
#     ->
# Resumen de todas las covariables para cada estación de contaminación
#     ->
# Relación general sin ningún tipo de agrupación
#
# En conjunto, estos scripts estudian la relación entre el NO2 y las covariables
# desde diferentes niveles de agregación.
#
# El script 14 ofrece la visión más general.
# El script 12 diferencia por estación del año y estación de medición.
# Los scripts 11 y 13 permiten estudiar cada estación individualmente.
#
# Esta comparación ayuda a evitar que una asociación global se interprete como
# una relación homogénea cuando en realidad puede depender de la localización,
# de la época del año o de la escala temporal utilizada.
library(data.table)
library(ggplot2)
library(here)

# ==============================================================================
# 1. CARGA DE DATOS
# ==============================================================================

cat("Cargando maestro DIARIO 2025...\n")
dt_d <- readRDS(here("data", "processed", "Maestro", "diario",
                     "dataset_maestro_inla_2025_DIARIO.rds"))
dt_d <- dt_d[!is.na(DATO_DIARIO)]

cat("Cargando maestro HORARIO 2025...\n")
dt_h <- readRDS(here("data", "processed", "Maestro", "horario",
                     "dataset_maestro_inla_2025_HORARIO.rds"))
dt_h <- dt_h[!is.na(DATO)]

cat(sprintf("  Diario  : %d filas\n", nrow(dt_d)))
cat(sprintf("  Horario : %d filas\n\n", nrow(dt_h)))

# ==============================================================================
# 2. COVARIABLES
# ==============================================================================

resolver_col <- function(dt, patron) {
  c <- grep(patron, names(dt), value = TRUE)
  c[grepl("_raw$", c)][1]
}

vars_cov <- list(
  list(col   = "Temperatura_raw",
       label = "Temperatura (\u00b0C)",
       file  = "temperatura"),

  list(col   = "Humedad_Relativa_raw",
       label = "Humedad relativa (%)",
       file  = "humedad"),

  list(col   = "Precipitaciones_raw",
       label = "Precipitaciones (mm)",
       file  = "precipitacion"),

  list(col   = resolver_col(dt_d, "^Presion"),
       label = "Presi\u00f3n barom\u00e9trica (mbar)",
       file  = "presion"),

  list(col   = resolver_col(dt_d, "^Radiaci"),
       label = "Radiaci\u00f3n solar (W/m\u00b2)",
       file  = "radiacion"),

  list(col   = resolver_col(dt_d, "^Velocidad"),
       label = "Velocidad del viento (m/s)",
       file  = "viento"),

  list(col   = "intensidad_raw",
       label = "Intensidad tr\u00e1fico (veh/h)",
       file  = "intensidad"),

  list(col   = "carga_raw",
       label = "Carga tr\u00e1fico",
       file  = "carga")
)

# ==============================================================================
# 3. DIRECTORIO DE SALIDA
# ==============================================================================

dir_salida <- here("outputs", "scatter_no2_covariables", "general")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Directorio: %s\n\n", dir_salida))

# ==============================================================================
# 4. TEMA
# ==============================================================================

tema_gen <- theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 9.5),
    plot.caption  = element_text(color = "gray55", size = 8),
    panel.grid.minor = element_blank(),
    axis.title    = element_text(size = 11)
  )

# ==============================================================================
# 5. BUCLE
# ==============================================================================

for (vi in vars_cov) {

  col <- vi$col
  if (is.null(col) || is.na(col) || col == "") {
    cat(sprintf("  [SKIP] %s\n", vi$file)); next
  }

  # ── DIARIO ──────────────────────────────────────────────────────────────────
  if (col %in% names(dt_d)) {

    dt_p <- dt_d[!is.na(get(col)), .(X = get(col), Y = DATO_DIARIO)]
    n    <- nrow(dt_p)

    p_d <- ggplot(dt_p, aes(X, Y)) +
      geom_point(size = 1.0, alpha = 0.45, color = "#2c7bb6", na.rm = TRUE) +
      labs(
        title    = sprintf("NO\u2082 vs %s \u2014 Madrid 2025", vi$label),
        subtitle = sprintf("Escala diaria  \u00b7  %s observaciones  \u00b7  Todas las estaciones",
                           format(n, big.mark = ".")),
        x        = vi$label,
        y        = "NO\u2082 (\u00b5g/m\u00b3)",
        caption  = "Sin l\u00ednea de tendencia  \u00b7  Datos diarios 2025"
      ) +
      tema_gen

    archivo_d <- file.path(dir_salida,
                            sprintf("scatter_general_diario_no2_vs_%s_2025.png", vi$file))
    ggsave(archivo_d, plot = p_d, width = 8, height = 6, dpi = 200, bg = "white")
    cat(sprintf("  \u2713 [diario]  %s\n", basename(archivo_d)))
  }

  # ── HORARIO ─────────────────────────────────────────────────────────────────
  col_h <- if (col %in% names(dt_h)) col else resolver_col(dt_h, sub("_raw", "", col))

  if (!is.null(col_h) && !is.na(col_h) && col_h %in% names(dt_h)) {

    dt_p <- dt_h[!is.na(get(col_h)), .(X = get(col_h), Y = DATO)]
    n    <- nrow(dt_p)

    p_h <- ggplot(dt_p, aes(X, Y)) +
      geom_point(size = 0.15, alpha = 0.12, color = "#d7191c", na.rm = TRUE) +
      labs(
        title    = sprintf("NO\u2082 vs %s \u2014 Madrid 2025", vi$label),
        subtitle = sprintf("Escala horaria  \u00b7  %s observaciones  \u00b7  Todas las estaciones",
                           format(n, big.mark = ".")),
        x        = vi$label,
        y        = "NO\u2082 (\u00b5g/m\u00b3)",
        caption  = "Sin l\u00ednea de tendencia  \u00b7  Datos horarios 2025"
      ) +
      tema_gen

    archivo_h <- file.path(dir_salida,
                            sprintf("scatter_general_horario_no2_vs_%s_2025.png", vi$file))
    ggsave(archivo_h, plot = p_h, width = 8, height = 6, dpi = 200, bg = "white")
    cat(sprintf("  \u2713 [horario] %s\n\n", basename(archivo_h)))
  }
}

# ==============================================================================
# 6. RESUMEN
# ==============================================================================

n_png <- length(list.files(dir_salida, pattern = "\\.png$"))
cat(sprintf("\n\u2713 %d PNG guardados en:\n  %s\n", n_png, dir_salida))
