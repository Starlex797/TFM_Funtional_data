# ==============================================================================
# ESTUDIO DE CALIDAD DE LOS DATOS (2019 - 2025)
# ==============================================================================
# Objetivo:
#   Cuantificar la calidad de los datos PROCESADOS que realmente entran en el
#   modelo, para el periodo 2019-2025, distinguiendo:
#
#     (A) COVARIABLES CLIMATICAS seleccionadas (6):
#         Temperatura, Humedad_Relativa, Precipitaciones,
#         Presion Barometrica, Radiacion Solar, Velocidad Viento.
#         (NO se estudian las descartadas: Radiacion Ultravioleta -solo 2019-
#          ni Direccion del Viento -excluida en la limpieza-).
#
#     (B) CONTAMINACION seleccionada (1):
#         NO2 (media diaria, log-transformada). No se estudian el resto de
#         magnitudes medidas (SO2, CO, NO, PM2.5, PM10, NOx, O3) porque no
#         se usan como variable respuesta del modelo.
#
# Se cuantifica la calidad en DOS niveles del pipeline:
#   1. DATO MEDIDO (estaciones): completitud real de la observacion, separando
#      huecos ESTRUCTURALES (la estacion no mide esa variable) de huecos
#      TEMPORALES (la estacion mide pero le falta un dia).
#   2. DATO FINAL DEL MODELO (tras relleno/interpolacion): completitud del
#      dato que efectivamente alimenta el modelo.
#
# METODO DE RELLENO DE HUECOS (ver detalle en la seccion final y en el .md):
#   - Clima (horario): imputacion en 2 pasos -> (1) interpolacion lineal
#     intra-estacion (zoo::na.approx, maxgap = 3 h) y (2) estacion mas cercana
#     (distancia euclidea) en el mismo instante. Estaciones 100% NA se respetan.
#   - Agregacion temporal con umbral de NA (dia NA si >=30% horas faltan;
#     mes NA si >=30% dias faltan). NO2: umbral del 20%.
#   - Clima -> ubicaciones de NO2: interpolacion espacial IDW (beta=2, min 7
#     estaciones), validada por LOOCV.
#   - NO2: los huecos NO se rellenan; los dias bajo umbral se marcan NA y el
#     modelo INLA los trata como respuesta faltante.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0. Parametros globales
# ------------------------------------------------------------------------------
ANIOS <- 2019

# Covariables climaticas SELECCIONADAS (nombre real en los .rds)
VARS_CLIMA <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  "Presion Barométrica", "Radiación Solar", "Velocidad Viento"
)

# Etiquetas legibles para tablas y figuras
ETIQUETAS <- c(
  "Temperatura"          = "Temperatura",
  "Humedad_Relativa"     = "Humedad relativa",
  "Precipitaciones"      = "Precipitaciones",
  "Presion Barométrica"  = "Presion barometrica",
  "Radiación Solar"      = "Radiacion solar",
  "Velocidad Viento"     = "Velocidad viento"
)

# Rutas de salida
DIR_TABLAS <- here("outputs", "tables", "calidad_datos")
DIR_FIGURAS <- here("outputs", "figures", "calidad_datos")
dir.create(DIR_TABLAS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIGURAS, recursive = TRUE, showWarnings = FALSE)

# Tema comun para las figuras
tema_calidad <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey30"),
    axis.text.x = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

cat("\n", strrep("=", 78), "\n", sep = "")
cat("ESTUDIO DE CALIDAD DE DATOS 2019-2025 | variables SELECCIONADAS\n")
cat(strrep("=", 78), "\n", sep = "")

# ==============================================================================
# 1. CALIDAD DE LAS COVARIABLES CLIMATICAS (DATO MEDIDO EN ESTACIONES)
# ==============================================================================
# Fuente: data/processed/Clima/diario/meteo_madrid_<anio>_diario.rds
# Este dato ya incorpora la imputacion horaria en 2 pasos y la agregacion a
# escala diaria (umbral 30% NA).
#
# MATIZ IMPORTANTE: cuando una estacion deja de disponer de sensor para una
# variable, NO aparece como fila con NA, sino que directamente NO tiene fila
# ese dia. Por eso la completitud se mide contra la REJILLA COMPLETA
# estacion x dia (todas las estaciones que reportan la variable en algun
# momento del anio x todos los dias del anio). Asi se capturan tambien las
# ausencias por estaciones que se apagan a mitad de anio (p.ej. la red de
# presion/radiacion que cae a 3 estaciones desde sep-2022).
#
# Se distinguen:
#   - COBERTURA ESPACIAL: nº de estaciones que reportan la variable (diaria).
#   - COMPLETITUD DE REJILLA: obs. validas / (estaciones_anio x dias).
# ------------------------------------------------------------------------------

resumen_clima <- list()
cobertura_estaciones <- list() # detalle estacion x variable x anio
cobertura_diaria <- list() # nº estaciones por dia (para figura temporal)

for (anio in ANIOS) {
  ruta <- here(
    "data", "processed", "Clima", "diario",
    sprintf("meteo_madrid_%d_diario.rds", anio)
  )
  if (!file.exists(ruta)) {
    warning("No existe: ", ruta)
    next
  }

  dt <- readRDS(ruta)
  setDT(dt)

  n_est_total <- uniqueN(dt$ESTACION) # estaciones meteo con alguna fila el anio
  dias_anio <- sort(unique(dt$FECHA))
  n_dias <- length(dias_anio)

  for (v in VARS_CLIMA) {
    if (!v %in% names(dt)) next

    # Estaciones que reportan la variable en algun momento del anio
    est_ano <- dt[!is.na(get(v)), unique(ESTACION)]
    n_est_ano <- length(est_ano)

    # Nº de estaciones con dato valido por dia (cobertura espacial diaria)
    por_dia <- dt[ESTACION %in% est_ano,
      .(n_validas = sum(!is.na(get(v)))),
      by = FECHA
    ]
    # Rellenar dias sin ninguna estacion (no aparecen) con 0
    por_dia <- por_dia[data.table(FECHA = dias_anio), on = "FECHA"]
    por_dia[is.na(n_validas), n_validas := 0L]

    # Completitud de REJILLA: obs validas / (estaciones_anio x dias)
    obs_validas <- sum(por_dia$n_validas)
    celdas_grid <- n_est_ano * n_dias
    compl_rejilla <- if (celdas_grid > 0) 100 * obs_validas / celdas_grid else NA_real_

    resumen_clima[[length(resumen_clima) + 1]] <- data.table(
      Anio                    = anio,
      Variable                = ETIQUETAS[[v]],
      N_estaciones_reportan   = n_est_ano,
      Estaciones_dia_min      = min(por_dia$n_validas),
      Estaciones_dia_mediana  = as.numeric(median(por_dia$n_validas)),
      Estaciones_dia_max      = max(por_dia$n_validas),
      Cobertura_espacial_pct  = round(100 * n_est_ano / n_est_total, 1),
      Completitud_rejilla_pct = round(compl_rejilla, 1),
      N_dias                  = n_dias
    )

    # Detalle estacion x variable (completitud de cada estacion sobre el anio)
    det_est <- dt[ESTACION %in% est_ano,
      .(dias_reportados = sum(!is.na(get(v)))),
      by = ESTACION
    ]
    det_est[, `:=`(
      Anio = anio, Variable = ETIQUETAS[[v]],
      completitud_pct = round(100 * dias_reportados / n_dias, 1)
    )]
    cobertura_estaciones[[length(cobertura_estaciones) + 1]] <-
      det_est[, .(Anio, Variable, ESTACION, dias_reportados, completitud_pct)]

    # Serie diaria de cobertura (para figura)
    por_dia[, `:=`(Anio = anio, Variable = ETIQUETAS[[v]])]
    cobertura_diaria[[length(cobertura_diaria) + 1]] <- por_dia
  }
}

tabla_clima_medido <- rbindlist(resumen_clima)
setorder(tabla_clima_medido, Variable, Anio)
tabla_cobertura_diaria <- rbindlist(cobertura_diaria)

cat("\n--- (A) COVARIABLES CLIMATICAS: dato MEDIDO en estaciones ---\n\n")
print(tabla_clima_medido)

fwrite(
  tabla_clima_medido,
  file.path(DIR_TABLAS, "calidad_clima_medido.csv")
)
fwrite(
  rbindlist(cobertura_estaciones),
  file.path(DIR_TABLAS, "cobertura_estaciones_clima.csv")
)

# ==============================================================================
# 2. CALIDAD DE LAS COVARIABLES CLIMATICAS (DATO FINAL: TRAS IDW)
# ==============================================================================
# Fuente: data/processed/Clima/diario/clima_interpolado_diario_<anio>.rds
# Clima interpolado por IDW a las ubicaciones de las estaciones de NO2. Aqui
# el hueco estructural desaparece: TODAS las ubicaciones reciben un valor
# estimado para las 6 covariables cada dia.
# ------------------------------------------------------------------------------

resumen_interp <- list()

for (anio in ANIOS) {
  ruta <- here(
    "data", "processed", "Clima", "diario",
    sprintf("clima_interpolado_diario_%d.rds", anio)
  )
  if (!file.exists(ruta)) {
    warning("No existe interpolado: ", ruta)
    next
  }

  dt <- readRDS(ruta)
  setDT(dt)

  for (v in VARS_CLIMA) {
    if (!v %in% names(dt)) next
    compl <- 100 * (nrow(dt) - sum(is.na(dt[[v]]))) / nrow(dt)
    resumen_interp[[length(resumen_interp) + 1]] <- data.table(
      Anio = anio,
      Variable = ETIQUETAS[[v]],
      N_ubicaciones = uniqueN(dt$ESTACION),
      N_dias = uniqueN(dt$FECHA),
      Completitud_pct = round(compl, 1)
    )
  }
}

tabla_clima_interp <- rbindlist(resumen_interp)
setorder(tabla_clima_interp, Variable, Anio)

cat("\n--- (A') COVARIABLES CLIMATICAS: dato FINAL tras interpolacion IDW ---\n\n")
print(tabla_clima_interp)
fwrite(
  tabla_clima_interp,
  file.path(DIR_TABLAS, "calidad_clima_interpolado.csv")
)

# ==============================================================================
# 3. CALIDAD DE LA CONTAMINACION (NO2 DIARIO)
# ==============================================================================
# Fuente: data/processed/Contaminacion/diario/aire_madrid_<anio>_No2_trans_diarios.rds
# El dato diario ya aplica: filtro de validacion horaria ("V"), agregacion a
# media diaria (dia NA si >=20% horas faltan) y log(x+1). Los huecos NO se
# rellenan: los dias bajo umbral quedan NA y el modelo los trata como faltantes.
# ------------------------------------------------------------------------------

resumen_no2 <- list()
compl_est_no2 <- list()

for (anio in ANIOS) {
  ruta <- here(
    "data", "processed", "Contaminacion", "diario",
    sprintf("aire_madrid_%d_No2_trans_diarios.rds", anio)
  )
  if (!file.exists(ruta)) {
    warning("No existe NO2: ", ruta)
    next
  }

  dt <- readRDS(ruta)
  setDT(dt)

  n_reg <- nrow(dt)
  n_na <- sum(is.na(dt$DATO_DIARIO))
  n_validos <- n_reg - n_na

  resumen_no2[[length(resumen_no2) + 1]] <- data.table(
    Anio               = anio,
    N_estaciones       = uniqueN(dt$ESTACION),
    N_dias             = uniqueN(dt$FECHA),
    N_registros        = n_reg,
    N_dias_validos     = n_validos,
    N_dias_NA          = n_na,
    Completitud_pct    = round(100 * n_validos / n_reg, 1)
  )

  # Completitud por estacion (para ver estaciones problematicas)
  por_est <- dt[, .(completitud_pct = round(100 * mean(!is.na(DATO_DIARIO)), 1)),
    by = ESTACION
  ][order(completitud_pct)]
  por_est[, Anio := anio]
  compl_est_no2[[length(compl_est_no2) + 1]] <- por_est
}

tabla_no2 <- rbindlist(resumen_no2)
setorder(tabla_no2, Anio)

cat("\n--- (B) CONTAMINACION NO2: dato diario ---\n\n")
print(tabla_no2)
fwrite(tabla_no2, file.path(DIR_TABLAS, "calidad_no2.csv"))
fwrite(
  rbindlist(compl_est_no2),
  file.path(DIR_TABLAS, "calidad_no2_por_estacion.csv")
)

# ==============================================================================
# 4. FIGURAS
# ==============================================================================

# --- 4.1 Heatmap: completitud de rejilla del clima medido (var x anio) --------
hm <- copy(tabla_clima_medido)
hm[, Anio := factor(Anio)]

g1 <- ggplot(hm, aes(x = Anio, y = Variable, fill = Completitud_rejilla_pct)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f", Completitud_rejilla_pct)),
    size = 3.4, color = "grey15"
  ) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "#F7F7A0", high = "#1A9850",
    midpoint = 55, limits = c(20, 100),
    name = "% rejilla\ncompleta"
  ) +
  labs(
    title = "Completitud del dato climatico medido (rejilla estacion x dia)",
    subtitle = "% de celdas con observacion valida sobre la red que reporta cada variable (2019-2025)",
    x = NULL, y = NULL
  ) +
  tema_calidad
ggsave(file.path(DIR_FIGURAS, "heatmap_completitud_rejilla_clima.png"),
  g1,
  width = 9, height = 4.5, dpi = 150
)

# --- 4.2 Heatmap: cobertura espacial (nº estaciones que reportan) -------------
g2 <- ggplot(hm, aes(x = Anio, y = Variable, fill = N_estaciones_reportan)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = N_estaciones_reportan),
    size = 3.6, color = "grey15"
  ) +
  scale_fill_gradient(
    low = "#FEE0D2", high = "#3182BD",
    name = "Estaciones\nque reportan"
  ) +
  labs(
    title = "Cobertura espacial de las covariables climaticas",
    subtitle = "Numero de estaciones meteorologicas que reportan cada variable (de 26 posibles)",
    x = NULL, y = NULL
  ) +
  tema_calidad
ggsave(file.path(DIR_FIGURAS, "heatmap_cobertura_espacial_clima.png"),
  g2,
  width = 9, height = 4.5, dpi = 150
)

# --- 4.2b Serie temporal: nº de estaciones por dia (revela la caida de la red) -
cd <- copy(tabla_cobertura_diaria)
cd[, FECHA := as.Date(FECHA)]
g2b <- ggplot(cd, aes(x = FECHA, y = n_validas, color = Variable)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "grey40") +
  annotate("text",
    x = min(cd$FECHA), y = 7.6, hjust = 0,
    label = "min. IDW = 7", size = 3, color = "grey40"
  ) +
  facet_wrap(~Variable, ncol = 2) +
  scale_x_date(date_labels = "%Y") +
  labs(
    title = "Cobertura diaria de estaciones meteorologicas por variable (2019-2025)",
    subtitle = "Nº de estaciones con dato valido cada dia. La caida bajo 7 impide la interpolacion IDW",
    x = NULL, y = "Estaciones con dato"
  ) +
  tema_calidad +
  theme(legend.position = "none")
ggsave(file.path(DIR_FIGURAS, "serie_cobertura_diaria_clima.png"),
  g2b,
  width = 11, height = 7, dpi = 150
)

# --- 4.3 NO2: completitud diaria por anio -------------------------------------
g3 <- ggplot(tabla_no2, aes(x = factor(Anio), y = Completitud_pct)) +
  geom_col(fill = "#3182BD", width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", Completitud_pct)),
    vjust = -0.5, size = 3.6
  ) +
  coord_cartesian(ylim = c(90, 100)) +
  labs(
    title = "Completitud del NO2 diario",
    subtitle = "% de registros estacion-dia validos tras el control de calidad (umbral 20% horas faltantes)",
    x = NULL, y = "% dias validos"
  ) +
  tema_calidad
ggsave(file.path(DIR_FIGURAS, "barras_completitud_no2.png"),
  g3,
  width = 8, height = 4.5, dpi = 150
)

# --- 4.4 Completitud del producto final (IDW) por variable y anio -------------
hi <- copy(tabla_clima_interp)
hi[, Anio := factor(Anio)]
g4 <- ggplot(hi, aes(x = Anio, y = Variable, fill = Completitud_pct)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f", Completitud_pct)),
    size = 3.4, color = "grey15"
  ) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "#F7F7A0", high = "#1A9850",
    midpoint = 80, limits = c(60, 100),
    name = "% dias\nestimados"
  ) +
  labs(
    title = "Completitud del producto FINAL de covariables (tras IDW)",
    subtitle = "% de dias con valor estimado en las 24 ubicaciones de NO2. Los huecos = dias con < 7 estaciones",
    x = NULL, y = NULL
  ) +
  tema_calidad
ggsave(file.path(DIR_FIGURAS, "heatmap_completitud_clima_interpolado.png"),
  g4,
  width = 9, height = 4.5, dpi = 150
)

cat("\nFiguras guardadas en:", DIR_FIGURAS, "\n")
cat("Tablas guardadas en:", DIR_TABLAS, "\n")

# ==============================================================================
# 5. RESUMEN GLOBAL EN CONSOLA
# ==============================================================================
cat("\n", strrep("=", 78), "\n", sep = "")
cat("RESUMEN GLOBAL 2019-2025\n")
cat(strrep("=", 78), "\n", sep = "")

cat("\nCLIMA (medido) - por variable (media 2019-2025):\n")
print(tabla_clima_medido[, .(
  Completitud_rejilla_media = round(mean(Completitud_rejilla_pct), 1),
  Estaciones_reportan_min   = min(N_estaciones_reportan),
  Estaciones_reportan_max   = max(N_estaciones_reportan)
),
by = Variable
])

cat(
  "\nCLIMA (final tras IDW) - completitud media:",
  round(mean(tabla_clima_interp$Completitud_pct), 1), "%\n"
)
cat(
  "  (ojo: 2022 sep-dic la red de presion/radiacion cae a 3 estaciones",
  "-> IDW imposible esos dias)\n"
)

cat(
  "\nNO2 - completitud diaria media:",
  round(mean(tabla_no2$Completitud_pct), 1), "% |",
  "rango:", round(min(tabla_no2$Completitud_pct), 1), "-",
  round(max(tabla_no2$Completitud_pct), 1), "%\n"
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("Estudio de calidad completado.\n")
cat(strrep("=", 78), "\n", sep = "")
