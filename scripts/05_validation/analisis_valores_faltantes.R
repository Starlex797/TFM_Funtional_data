# ==============================================================================
# ESTUDIO DE VALORES FALTANTES DEL CLIMA (2019-2025)
# ==============================================================================
# Sustituye a la sección 1 de analisis_calidad_datos.R. La diferencia es que
# aquí NO se deduce el tipo de hueco a posteriori: se lee directamente de las
# columnas <variable>_estado que genera el preprocesamiento.
#
# Los cuatro estados:
#   OK          medida real del sensor
#   IMPUTADO    estaba vacía y el relleno le puso valor
#   FALLO       estaba vacía y sigue vacía  -> avería
#   SIN_SENSOR  la estación no dispone del aparato ese año
#
# DENOMINADOR. Todas las tasas se calculan sobre las celdas MEDIBLES, es decir
# excluyendo SIN_SENSOR. No se le puede reprochar a una estación sin barómetro
# que no mida la presión. Esto es lo que separa "la presión tiene un 30 % de
# completitud" (falso, mezcla estaciones sin aparato) de "los 8 barómetros de
# la red miden al 99 %" (cierto).
#
#   completitud_medida = OK                / medibles
#   tasa_imputacion    = IMPUTADO          / medibles
#   tasa_averia        = FALLO             / medibles
#   completitud_final  = (OK + IMPUTADO)   / medibles   <- lo que ve el modelo
#
# Entradas : data/processed/Clima/{horario,diario}/meteo_madrid_<anio>_<escala><SUFIJO>.rds
# Salidas  : outputs/tables/valores_faltantes/   y   outputs/figures/valores_faltantes/
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0. Parámetros
# ------------------------------------------------------------------------------
ANIOS  <- 2019:2025
SUFIJO <- "2"          # versión de los ficheros con columnas _estado
ESCALA <- "horario"    # "horario" = donde ocurren las averías; "diario" = escala del modelo

VARS_CLIMA <- c("Temperatura", "Humedad_Relativa", "Precipitaciones",
                "Presion Barométrica", "Radiación Solar", "Velocidad Viento")

# Rótulos en inglés, coherentes con los mapas ya generados para la memoria.
ETIQUETAS <- c(
  "Temperatura"          = "Temperature",
  "Humedad_Relativa"     = "Relative humidity",
  "Precipitaciones"      = "Precipitation",
  "Presion Barométrica"  = "Barometric pressure",
  "Radiación Solar"      = "Solar radiation",
  "Velocidad Viento"     = "Wind speed"
)

ESTADOS <- c("OK", "IMPUTADO", "FALLO", "SIN_SENSOR")

COLOR_ESTADO <- c(
  "Measured"    = "#2C6E4E",
  "Imputed"     = "#C08A1E",
  "Failure"     = "#A8383A",
  "No sensor"   = "#B9C2C9"
)

ESTADO_EN <- c("OK" = "Measured", "IMPUTADO" = "Imputed",
               "FALLO" = "Failure", "SIN_SENSOR" = "No sensor")

DIR_TAB <- here("outputs", "tables",  "valores_faltantes")
DIR_FIG <- here("outputs", "figures", "valores_faltantes")
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

tema <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey35", size = 9.5),
    plot.caption     = element_text(color = "grey45", size = 7.5, hjust = 0),
    strip.text       = element_text(face = "bold", size = 9)
  )

# ------------------------------------------------------------------------------
# 1. Carga: una fila por estación x tiempo x variable, con su estado
# ------------------------------------------------------------------------------
ruta_anio <- function(anio) {
  here("data", "processed", "Clima", ESCALA,
       sprintf("meteo_madrid_%d_%s%s.rds", anio, ESCALA, SUFIJO))
}

largo <- list()
anios_ausentes <- integer(0)

for (anio in ANIOS) {

  ruta <- ruta_anio(anio)
  if (!file.exists(ruta)) { anios_ausentes <- c(anios_ausentes, anio); next }

  dt <- as.data.table(readRDS(ruta))

  cols_estado <- paste0(VARS_CLIMA, "_estado")
  faltan <- setdiff(cols_estado, names(dt))
  if (length(faltan) > 0)
    stop("El fichero de ", anio, " no tiene columnas de estado: ",
         paste(faltan, collapse = ", "),
         "\nRegenera el preprocesamiento con marcar_estado = TRUE.")

  # Formato largo: una fila por celda, que es la unidad natural del estudio.
  m <- melt(dt[, c("ESTACION", "FECHA", cols_estado), with = FALSE],
            id.vars       = c("ESTACION", "FECHA"),
            measure.vars  = cols_estado,
            variable.name = "Variable",
            value.name    = "Estado")

  m[, Variable := sub("_estado$", "", as.character(Variable))]
  m[, `:=`(Anio = anio, Estado = factor(as.character(Estado), levels = ESTADOS))]

  largo[[length(largo) + 1]] <- m[, .(Anio, ESTACION, FECHA, Variable, Estado)]
}

if (length(anios_ausentes) > 0)
  warning("Años sin fichero (se omiten): ", paste(anios_ausentes, collapse = ", "))
if (length(largo) == 0)
  stop("No se ha cargado ningún año. Revisa SUFIJO y ESCALA.")

celdas <- rbindlist(largo)
celdas[, Variable_lab := factor(ETIQUETAS[Variable], levels = unname(ETIQUETAS))]

cat("Celdas analizadas:", format(nrow(celdas), big.mark = ","),
    "| escala:", ESCALA, "| años:", uniqueN(celdas$Anio), "\n\n")

# ------------------------------------------------------------------------------
# 2. Tabla A: reparto por variable y año
# ------------------------------------------------------------------------------
resumir <- function(dt, ...) {
  dt[, .(
    n_celdas   = .N,
    n_ok       = sum(Estado == "OK"),
    n_imputado = sum(Estado == "IMPUTADO"),
    n_fallo    = sum(Estado == "FALLO"),
    n_sin      = sum(Estado == "SIN_SENSOR")
  ), by = ...][
    , medibles := n_ok + n_imputado + n_fallo][
    , `:=`(
      completitud_medida = round(100 * n_ok / medibles, 2),
      tasa_imputacion    = round(100 * n_imputado / medibles, 2),
      tasa_averia        = round(100 * n_fallo / medibles, 2),
      completitud_final  = round(100 * (n_ok + n_imputado) / medibles, 2)
    )][]
}

tabla_var <- resumir(celdas, .(Anio, Variable = Variable_lab))
setorder(tabla_var, Variable, Anio)

# Nº de estaciones que realmente disponen del sensor cada año.
estaciones_activas <- celdas[Estado != "SIN_SENSOR",
                             .(n_estaciones = uniqueN(ESTACION)),
                             by = .(Anio, Variable = Variable_lab)]
tabla_var <- merge(tabla_var, estaciones_activas,
                   by = c("Anio", "Variable"), all.x = TRUE)
tabla_var[is.na(n_estaciones), n_estaciones := 0L]

cat("--- A) Reparto por variable y año ---\n")
print(tabla_var[, .(Anio, Variable, n_estaciones, completitud_medida,
                    tasa_imputacion, tasa_averia)])
fwrite(tabla_var, file.path(DIR_TAB, sprintf("faltantes_por_variable_%s.csv", ESCALA)))

# ------------------------------------------------------------------------------
# 3. Tabla B: por estación (solo donde la estación tiene el sensor)
# ------------------------------------------------------------------------------
tabla_est <- resumir(celdas[Estado != "SIN_SENSOR"],
                     .(Anio, ESTACION, Variable = Variable_lab))
setorder(tabla_est, Variable, Anio, completitud_medida)

fwrite(tabla_est, file.path(DIR_TAB, sprintf("faltantes_por_estacion_%s.csv", ESCALA)))

cat("\n--- B) Las 12 peores combinaciones estación x variable x año ---\n")
print(head(tabla_est[, .(Anio, ESTACION, Variable, completitud_medida,
                         tasa_imputacion, tasa_averia)], 12))

# ------------------------------------------------------------------------------
# 4. Tabla C: inventario OBSERVADO de sensores
# ------------------------------------------------------------------------------
# Qué mide cada estación cada año, deducido del dato y no de estaciones.csv.
# Es la tabla que documenta la reducción de la red.
inventario_obs <- celdas[, .(tiene_sensor = any(Estado != "SIN_SENSOR")),
                         by = .(Anio, ESTACION, Variable = Variable_lab)]

fwrite(dcast(inventario_obs, ESTACION + Variable ~ Anio, value.var = "tiene_sensor"),
       file.path(DIR_TAB, "inventario_sensores_observado.csv"))

cat("\n--- C) Estaciones con sensor por variable y año ---\n")
print(dcast(inventario_obs[tiene_sensor == TRUE, .N, by = .(Variable, Anio)],
            Variable ~ Anio, value.var = "N", fill = 0L))

# ------------------------------------------------------------------------------
# 5. Figuras
# ------------------------------------------------------------------------------

# 5.1 Reparto de estados (barras apiladas) -------------------------------------
prop <- celdas[, .N, by = .(Anio, Variable_lab, Estado)]
prop[, pct := 100 * N / sum(N), by = .(Anio, Variable_lab)]
prop[, Estado_lab := factor(ESTADO_EN[as.character(Estado)],
                            levels = names(COLOR_ESTADO))]

g1 <- ggplot(prop, aes(x = factor(Anio), y = pct, fill = Estado_lab)) +
  geom_col(width = 0.78) +
  facet_wrap(~ Variable_lab, ncol = 3) +
  scale_fill_manual(values = COLOR_ESTADO, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "Where every climate cell comes from",
    subtitle = sprintf("Share of %s cells by origin, Madrid 2019-2025", ESCALA),
    x = NULL, y = "% of cells",
    caption  = "'No sensor' is a structural gap: the station never carried that instrument that year."
  ) + tema + theme(legend.position = "bottom")

ggsave(file.path(DIR_FIG, sprintf("reparto_estados_%s.png", ESCALA)),
       g1, width = 10, height = 6, dpi = 300, bg = "white")

# 5.2 Completitud medida, con denominador correcto ------------------------------
g2 <- ggplot(tabla_var, aes(x = factor(Anio), y = Variable, fill = completitud_medida)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.1f", completitud_medida)), size = 3.1, color = "grey12") +
  scale_fill_gradient(low = "#F5D9D9", high = "#2C6E4E",
                      limits = c(80, 100), oob = scales::squish,
                      name = "% measured") +
  labs(
    title    = "Completeness among stations that carry the sensor",
    subtitle = "Real measurements over measurable cells; stations without the instrument are excluded",
    x = NULL, y = NULL,
    caption  = "Excluding 'no sensor' from the denominator is what separates a sparse network from a faulty one."
  ) + tema

ggsave(file.path(DIR_FIG, sprintf("completitud_medida_%s.png", ESCALA)),
       g2, width = 9, height = 4.2, dpi = 300, bg = "white")

# 5.3 Tamaño de la red: estaciones con sensor ----------------------------------
g3 <- ggplot(tabla_var, aes(x = Anio, y = n_estaciones, color = Variable)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "grey45") +
  annotate("text", x = min(ANIOS), y = 7.5, hjust = 0, size = 2.9, color = "grey45",
           label = "IDW minimum = 7") +
  facet_wrap(~ Variable, ncol = 3) +
  scale_x_continuous(breaks = ANIOS) +
  labs(
    title    = "Size of the measuring network by variable",
    subtitle = "Stations actually carrying each instrument, derived from the data itself",
    x = NULL, y = "Stations with sensor"
  ) + tema + theme(legend.position = "none",
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))

ggsave(file.path(DIR_FIG, "tamano_red_por_variable.png"),
       g3, width = 10, height = 5.5, dpi = 300, bg = "white")

# 5.4 Averías por estación -----------------------------------------------------
# Solo las estaciones que tienen el sensor: aquí sí se puede hablar de fallo.
g4 <- ggplot(tabla_est, aes(x = factor(Anio), y = ESTACION, fill = tasa_averia)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_wrap(~ Variable, ncol = 3) +
  scale_fill_gradient(low = "#EDF1F4", high = "#A8383A", name = "% failure") +
  labs(
    title    = "Outages by station",
    subtitle = "Share of measurable cells left empty after imputation",
    x = NULL, y = NULL,
    caption  = "Blank cells: the station did not carry that instrument that year."
  ) + tema +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 7))

ggsave(file.path(DIR_FIG, sprintf("averias_por_estacion_%s.png", ESCALA)),
       g4, width = 12, height = 8, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# 6. Comprobaciones de coherencia
# ------------------------------------------------------------------------------
# El estudio no vale nada si las categorías no cuadran con el total de celdas.
cat("\n--- COHERENCIA ---\n")

err <- tabla_var[n_ok + n_imputado + n_fallo + n_sin != n_celdas]
cat("Filas donde los 4 estados no suman el total:", nrow(err), "\n")

sin_estado <- celdas[is.na(Estado), .N]
cat("Celdas sin estado asignado:", sin_estado, "\n")

sin_medibles <- tabla_var[medibles == 0, .N]
cat("Combinaciones variable-año sin ninguna estación con sensor:", sin_medibles, "\n")

if (nrow(err) == 0 && sin_estado == 0)
  cat("\nOK: el reparto cuadra en todas las combinaciones.\n")
else
  warning("Hay incoherencias: revisa el preprocesamiento antes de usar estas cifras.")

cat("\nTablas  ->", DIR_TAB, "\n")
cat("Figuras ->", DIR_FIG, "\n")
