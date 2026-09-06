# ==============================================================================
# 29 — DISTANCIAS ENTRE LAS 24 ESTACIONES DE NO2 -> PRIOR PC DEL RANGO (SPDE)
# ==============================================================================
# Objetivo: para fijar el prior PC del RANGO del campo Matérn (SPDE) hay que
# conocer la ESCALA ESPACIAL real de la red. Este script:
#   1) Toma las 24 estaciones de NO2 y sus coordenadas oficiales en UTM 30N
#      (EPSG:25830) desde el fichero de estaciones (COORDENADA_X/Y_ETRS89),
#      pasadas a km; no se reproyecta desde LONGITUD/LATITUD.
#   2) Calcula la matriz de distancias reales (276 pares) y sus resúmenes:
#      mínima, vecino más cercano, mediana, máxima.
#   3) Fija el diámetro del dominio en D = 28 km (extensión del área de estudio)
#      y de ahí el rango a priori rho = D/3.
#   4) Sugiere valores para prior.range del PC prior a partir de esas distancias.
#
# Recordatorio del PC prior: prior.range = c(r0, p) codifica P(rango < r0) = p.
#   - p pequeño con r0 pequeño  -> penaliza rangos CORTOS (campo suave).
#   - p = 0.5 con r0            -> mediana del rango = r0.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
})
source(here("R", "utilities", "academic_quality_tables.R"))
source(here("R", "utilities", "dictionaries.R"))

# Diámetro del dominio de estudio (km). Se fija a priori: es la extensión del
# área cubierta por la red, no la distancia máxima observada entre dos estaciones.
DIAMETRO_KM <- 28

DIR_FIG <- here("outputs", "figures", "distancias_estaciones")
DIR_TAB <- here("outputs", "tables", "distancias_estaciones")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("DISTANCIAS ENTRE ESTACIONES DE NO2 -> PRIOR PC DEL RANGO\n")
cat(strrep("=", 70), "\n", sep = "")

# ==============================================================================
# 1. Coordenadas UTM de las 24 estaciones de NO2
# ==============================================================================
# El .rds diario solo guarda LONGITUD/LATITUD, pero el fichero de estaciones ya
# trae las coordenadas oficiales en ETRS89 / UTM 30N (EPSG:25830), en metros y
# con coma decimal: COORDENADA_X_ETRS89 / COORDENADA_Y_ETRS89. Se usan esas —
# las mismas del padrón municipal — en lugar de reproyectar desde lon/lat.
f_no2 <- here(
  "data", "processed", "Contaminacion", "diario",
  "aire_madrid_2025_No2_trans_diarios1.rds"
)
d <- readRDS(f_no2)
setDT(d)
est_no2 <- sort(unique(d$ESTACION))

f_ubic <- here("data", "raw", "Datos_contaminacion", "Estaciones", "datos.csv")
ubic <- fread(f_ubic)
setDT(ubic)
cols_req <- c("CODIGO_CORTO", "COORDENADA_X_ETRS89", "COORDENADA_Y_ETRS89")
if (!all(cols_req %in% names(ubic))) {
  stop(
    "Faltan columnas UTM en el fichero de estaciones: ",
    paste(setdiff(cols_req, names(ubic)), collapse = ", ")
  )
}

# Metros -> km. El nombre de estación se recupera del diccionario, que es el que
# se aplicó al generar los datos procesados.
num_utm <- function(x) as.numeric(gsub(",", ".", as.character(x)))
est <- ubic[, .(
  ESTACION = unname(nombres_estaciones_aire[as.character(CODIGO_CORTO)]),
  X_km = num_utm(COORDENADA_X_ETRS89) / 1000,
  Y_km = num_utm(COORDENADA_Y_ETRS89) / 1000
)][ESTACION %in% est_no2]
est <- unique(est, by = "ESTACION")
setorder(est, ESTACION)

faltan <- setdiff(est_no2, est$ESTACION)
if (length(faltan) > 0) {
  stop("Estaciones de NO2 sin coordenadas UTM: ", paste(faltan, collapse = ", "))
}
cat(sprintf("\nEstaciones de NO2: %d\n", nrow(est)))

# Validación: coordenadas UTM 30N plausibles para Madrid (km)
fuera <- est[!is.finite(X_km) | !is.finite(Y_km) |
  X_km < 400 | X_km > 480 | Y_km < 4440 | Y_km > 4520]
if (nrow(fuera) > 0) {
  cat("  AVISO: coordenadas UTM fuera del rango esperado de Madrid:\n")
  print(fuera)
}

# ==============================================================================
# 2. Extensión del dominio (caja envolvente, km)
# ==============================================================================
ancho_km <- diff(range(est$X_km))
alto_km <- diff(range(est$Y_km))
diag_km <- sqrt(ancho_km^2 + alto_km^2)

# ==============================================================================
# 3. Matriz de distancias reales (km) y resúmenes
# ==============================================================================
D <- as.matrix(dist(est[, .(X_km, Y_km)])) # 24 x 24
rownames(D) <- colnames(D) <- est$ESTACION
pares <- D[upper.tri(D)] # 276 distancias únicas

# Distancia al vecino más cercano de cada estación
diag(D) <- Inf
nn <- apply(D, 1, min)
diag(D) <- 0

resumen <- data.table(
  Medida = c(
    "Nº estaciones", "Nº pares",
    "Dist. mínima (par)", "Vecino más cercano (mediana)",
    "Vecino más cercano (máx)", "Distancia mediana (pares)",
    "Distancia media (pares)", "Distancia máxima (par)",
    "Ancho dominio (X)", "Alto dominio (Y)", "Diagonal dominio"
  ),
  Valor_km = round(c(
    nrow(est), length(pares),
    min(pares), median(nn), max(nn), median(pares),
    mean(pares), max(pares), ancho_km, alto_km, diag_km
  ), 2)
)
cat("\n--- Resumen de distancias (km) ---\n")
print(resumen)
fwrite(resumen, file.path(DIR_TAB, "resumen_distancias.csv"))

# Pares extremos: qué estaciones marcan la distancia mínima y la máxima
idx <- which(upper.tri(D), arr.ind = TRUE)
par_min <- idx[which.min(pares), ]
par_max <- idx[which.max(pares), ]
nom_min <- paste(est$ESTACION[par_min["row"]], "-", est$ESTACION[par_min["col"]])
nom_max <- paste(est$ESTACION[par_max["row"]], "-", est$ESTACION[par_max["col"]])

d_min <- min(pares) # distancia mínima observada entre dos estaciones
d_max <- max(pares) # distancia máxima observada entre dos estaciones

# Diámetro del dominio: valor fijado (extensión del área de estudio), no la
# distancia máxima observada. Al no conocerse el rango real del campo Matérn,
# se toma rho = D/3.
diametro <- DIAMETRO_KM
rho_prior <- diametro / 3

# Cuantiles de las distancias entre pares
q <- quantile(pares, probs = c(0, .1, .25, .5, .75, .9, 1))
cat("\n--- Cuantiles de las distancias entre pares (km) ---\n")
print(round(q, 2))
fwrite(
  data.table(Cuantil = names(q), km = round(as.numeric(q), 2)),
  file.path(DIR_TAB, "cuantiles_distancias.csv")
)

# ==============================================================================
# 4. TABLA ACADÉMICA: escala espacial de la red y rango a priori
# ==============================================================================
# Al desconocerse el rango real del campo Matérn, se fija rho = D/3, donde D es
# el diámetro del dominio de estudio (28 km). Es la regla habitual: el campo
# decorrelaciona a un tercio de la extensión del dominio, lo que evita tanto un
# campo plano (rho ~ D) como uno puramente local (rho ~ d_min).
tabla_acad <- data.table(
  Medida = c(
    "Distancia mínima entre estaciones",
    "Distancia máxima entre estaciones",
    "Diámetro del dominio de estudio",
    "Rango a priori del campo Matérn"
  ),
  Símbolo = c("d_min", "d_max", "D", "rho = D/3"),
  Estaciones = c(nom_min, nom_max, "--", "--"),
  `Valor (km)` = sprintf("%.2f", c(d_min, d_max, diametro, rho_prior))
)
cat("\n--- Tabla académica: escala espacial y rango a priori ---\n")
print(tabla_acad)
fwrite(tabla_acad, file.path(DIR_TAB, "tabla_escala_espacial_rho.csv"))

# Versión LaTeX (booktabs) lista para el documento
tex <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Escala espacial de la red de estaciones de NO$_2$ y rango a priori del campo espacial.}",
  "\\label{tab:escala-espacial-rho}",
  "\\begin{tabular}{lllr}", "\\toprule",
  "Medida & S\\'imbolo & Estaciones & Valor (km) \\\\", "\\midrule",
  sprintf("Distancia m\\'inima entre estaciones & $d_{\\min}$ & %s & %.2f \\\\", nom_min, d_min),
  sprintf("Distancia m\\'axima entre estaciones & $d_{\\max}$ & %s & %.2f \\\\", nom_max, d_max),
  sprintf("Di\\'ametro del dominio de estudio & $D$ & --- & %.2f \\\\", diametro),
  sprintf("Rango a priori & $\\rho = D/3$ & --- & %.2f \\\\", rho_prior),
  "\\bottomrule", "\\end{tabular}",
  sprintf(paste0(
    "\\begin{tablenotes}\\small\\item Nota: %d estaciones (%d pares). ",
    "Coordenadas oficiales ETRS89 / UTM 30N (EPSG:25830) del fichero de ",
    "estaciones (COORDENADA\\_X/Y\\_ETRS89), en km. Al desconocerse el rango ",
    "real, se fija $\\rho$ en un tercio del di\\'ametro del dominio.",
    "\\end{tablenotes}"
  ), nrow(est), length(pares)),
  "\\end{table}"
)
writeLines(tex, file.path(DIR_TAB, "tabla_escala_espacial_rho.tex"))

# Versión PNG con estilo booktabs (misma utilidad que las tablas de calidad)
booktabs_png(
  tabla_acad, file.path(DIR_FIG, "tabla_escala_espacial_rho.png"),
  title = "Escala espacial de la red de estaciones de NO2",
  subtitle = "Distancias entre estaciones y rango a priori del campo espacial (SPDE)",
  note = sprintf(
    paste0(
      "Nota: %d estaciones (%d pares). Coordenadas oficiales ETRS89 / UTM 30N ",
      "(EPSG:25830) del fichero de estaciones, en km, y distancias euclídeas. ",
      "El diámetro D = %.0f km es la extensión del dominio de estudio. Al ",
      "desconocerse el rango real del campo Matérn, se fija rho = D/3."
    ),
    nrow(est), length(pares), diametro
  ),
  widths = c(3.60, 0.95, 2.30, 1.05),
  align = c("left", "center", "left", "right")
)

# ==============================================================================
# 5. SUGERENCIA de prior PC para el rango
# ==============================================================================
# Criterios habituales:
#   - No tiene sentido un rango menor que la separación típica entre estaciones
#     (no se puede resolver más fino) -> usar el vecino más cercano como cota.
#   - Ni mayor que el dominio -> la mediana razonable está entre la mediana de
#     los pares y ~1/3 del diámetro.
r0_pen <- round(median(nn), 1) # penalizar rangos < separación típica
r0_med <- round(median(pares) / 2, 1) # mediana del rango ~ mitad de la mediana de pares
r0_dom <- round(rho_prior, 1) # 1/3 del diámetro del dominio (D = 28 km)

sugerencias <- data.table(
  Opcion = c(
    "A) Penalizar rangos cortos", "B) Mediana en ~1/2 mediana de pares",
    "C) Mediana en ~1/3 del diámetro"
  ),
  prior.range = c(
    sprintf("c(%.1f, 0.1)", r0_pen),
    sprintf("c(%.1f, 0.5)", r0_med),
    sprintf("c(%.1f, 0.5)", r0_dom)
  ),
  Interpretacion = c(
    sprintf("P(rango < %.1f km) = 0.1: rango casi seguro > separación típica", r0_pen),
    sprintf("mediana del rango = %.1f km", r0_med),
    sprintf("mediana del rango = %.1f km", r0_dom)
  )
)
cat("\n--- Sugerencias de prior.range (PC prior) ---\n")
print(sugerencias)
fwrite(sugerencias, file.path(DIR_TAB, "sugerencias_prior_rango.csv"))

cat(sprintf(
  paste0(
    "\n>> Escala de la red: d_min ~%.2f km | d_max ~%.2f km | ",
    "diámetro D = %.0f km | rho = D/3 ~%.2f km\n"
  ),
  d_min, d_max, diametro, rho_prior
))
cat(sprintf(
  ">> Regla de malla: max.edge del SPDE debe ser < rango/3 (~%.1f km)\n",
  rho_prior / 3
))

# ==============================================================================
# 6. Figuras: histograma de distancias + mapa de estaciones
# ==============================================================================
g_hist <- ggplot(data.table(d = pares), aes(d)) +
  geom_histogram(bins = 30, fill = "#2166AC", color = "white") +
  geom_vline(xintercept = median(nn), color = "#1A9850", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = median(pares), color = "#B2182B", linewidth = 1) +
  geom_vline(xintercept = rho_prior, color = "#762A83", linewidth = 1, linetype = "dotted") +
  annotate("text",
    x = median(nn), y = Inf, label = "vecino cercano (mediana)",
    vjust = 2, hjust = -0.05, color = "#1A9850", size = 3
  ) +
  annotate("text",
    x = median(pares), y = Inf, label = "mediana de pares",
    vjust = 3.5, hjust = -0.05, color = "#B2182B", size = 3
  ) +
  annotate("text",
    x = rho_prior, y = Inf, label = "rho = D/3",
    vjust = 5, hjust = -0.05, color = "#762A83", size = 3
  ) +
  labs(
    title = "Distancias entre las 24 estaciones de NO2",
    subtitle = sprintf(
      "%d pares | líneas: vecino más cercano (mediana), mediana de pares y rango a priori",
      length(pares)
    ),
    x = "Distancia (km)", y = "Nº de pares"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_FIG, "histograma_distancias.png"), g_hist, width = 9, height = 5, dpi = 150)

g_map <- ggplot(est, aes(X_km, Y_km)) +
  geom_point(color = "#2166AC", size = 3) +
  geom_text(aes(label = ESTACION), size = 2.4, vjust = -1, check_overlap = TRUE) +
  coord_equal() +
  labs(
    title = "Estaciones de NO2 (UTM 30N, km)",
    subtitle = sprintf("Dominio ~ %.0f × %.0f km", ancho_km, alto_km),
    x = "X (km)", y = "Y (km)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_FIG, "mapa_estaciones.png"), g_map, width = 8, height = 8, dpi = 150)

fwrite(
  est[, .(ESTACION, X_km = round(X_km, 3), Y_km = round(Y_km, 3))],
  file.path(DIR_TAB, "estaciones_coordenadas.csv")
)

cat("\nFiguras en:", DIR_FIG, "\nTablas en:", DIR_TAB, "\n")
cat("--- Completado ---\n")
