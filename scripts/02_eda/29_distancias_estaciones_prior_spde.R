# ==============================================================================
# 29 — DISTANCIAS ENTRE LAS 24 ESTACIONES DE NO2 -> PRIOR PC DEL RANGO (SPDE)
# ==============================================================================
# Objetivo: para fijar el prior PC del RANGO del campo Matérn (SPDE) hay que
# conocer la ESCALA ESPACIAL real de la red. Este script:
#   1) Toma las 24 estaciones de NO2 (datos ya procesados).
#   2) Proyecta sus coordenadas a UTM 30N (EPSG:25830) en km — IGUAL que el
#      modelo (X_km/Y_km).
#   3) Calcula la matriz de distancias reales (276 pares) y sus resúmenes:
#      mínima, vecino más cercano, mediana, máxima (diámetro del dominio).
#   4) Sugiere valores para prior.range del PC prior a partir de esas distancias.
#
# Recordatorio del PC prior: prior.range = c(r0, p) codifica P(rango < r0) = p.
#   - p pequeño con r0 pequeño  -> penaliza rangos CORTOS (campo suave).
#   - p = 0.5 con r0            -> mediana del rango = r0.
# ==============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(sf); library(ggplot2)
})

DIR_FIG <- here("outputs", "figures", "distancias_estaciones")
DIR_TAB <- here("outputs", "tables", "distancias_estaciones")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("DISTANCIAS ENTRE ESTACIONES DE NO2 -> PRIOR PC DEL RANGO\n")
cat(strrep("=", 70), "\n", sep = "")

# ==============================================================================
# 1. Coordenadas de las 24 estaciones de NO2 (datos procesados)
# ==============================================================================
f_no2 <- here("data", "processed", "Contaminacion", "diario",
              "aire_madrid_2025_No2_trans_diarios.rds")
d <- readRDS(f_no2); setDT(d)
est <- unique(d[, .(ESTACION, LONGITUD, LATITUD)])
setorder(est, ESTACION)
cat(sprintf("\nEstaciones de NO2: %d\n", nrow(est)))

# Validación: coordenadas dentro del rango de Madrid
fuera <- est[LATITUD < 40.2 | LATITUD > 40.7 | LONGITUD < -4.0 | LONGITUD > -3.3]
if (nrow(fuera) > 0) {
  cat("  AVISO: coordenadas fuera del rango de Madrid:\n"); print(fuera)
}

# ==============================================================================
# 2. Proyección a UTM 30N (EPSG:25830) en km — igual que el modelo
# ==============================================================================
sf_est <- st_transform(st_as_sf(est, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830)
xy_km <- st_coordinates(sf_est) / 1000
est[, `:=`(X_km = xy_km[, 1], Y_km = xy_km[, 2])]

# Extensión del dominio (caja envolvente)
ancho_km <- diff(range(est$X_km))
alto_km  <- diff(range(est$Y_km))
diag_km  <- sqrt(ancho_km^2 + alto_km^2)

# ==============================================================================
# 3. Matriz de distancias reales (km) y resúmenes
# ==============================================================================
D <- as.matrix(dist(est[, .(X_km, Y_km)]))         # 24 x 24
rownames(D) <- colnames(D) <- est$ESTACION
pares <- D[upper.tri(D)]                             # 276 distancias únicas

# Distancia al vecino más cercano de cada estación
diag(D) <- Inf
nn <- apply(D, 1, min)
diag(D) <- 0

resumen <- data.table(
  Medida = c("Nº estaciones", "Nº pares",
             "Dist. mínima (par)", "Vecino más cercano (mediana)",
             "Vecino más cercano (máx)", "Distancia mediana (pares)",
             "Distancia media (pares)", "Distancia máxima (diámetro)",
             "Ancho dominio (X)", "Alto dominio (Y)", "Diagonal dominio"),
  Valor_km = round(c(nrow(est), length(pares),
             min(pares), median(nn), max(nn), median(pares),
             mean(pares), max(pares), ancho_km, alto_km, diag_km), 2))
cat("\n--- Resumen de distancias (km) ---\n"); print(resumen)
fwrite(resumen, file.path(DIR_TAB, "resumen_distancias.csv"))

# Cuantiles de las distancias entre pares
q <- quantile(pares, probs = c(0, .1, .25, .5, .75, .9, 1))
cat("\n--- Cuantiles de las distancias entre pares (km) ---\n")
print(round(q, 2))
fwrite(data.table(Cuantil = names(q), km = round(as.numeric(q), 2)),
       file.path(DIR_TAB, "cuantiles_distancias.csv"))

# ==============================================================================
# 4. SUGERENCIA de prior PC para el rango
# ==============================================================================
# Criterios habituales:
#   - No tiene sentido un rango menor que la separación típica entre estaciones
#     (no se puede resolver más fino) -> usar el vecino más cercano como cota.
#   - Ni mayor que el dominio -> la mediana razonable está entre la mediana de
#     los pares y ~1/2 del diámetro.
r0_pen  <- round(median(nn), 1)          # penalizar rangos < separación típica
r0_med  <- round(median(pares) / 2, 1)   # mediana del rango ~ mitad de la mediana de pares
r0_dom  <- round(diag_km / 3, 1)         # ~1/3 del diámetro del dominio

sugerencias <- data.table(
  Opcion = c("A) Penalizar rangos cortos", "B) Mediana en ~1/2 mediana de pares",
             "C) Mediana en ~1/3 del diámetro"),
  prior.range = c(sprintf("c(%.1f, 0.1)", r0_pen),
                  sprintf("c(%.1f, 0.5)", r0_med),
                  sprintf("c(%.1f, 0.5)", r0_dom)),
  Interpretacion = c(
    sprintf("P(rango < %.1f km) = 0.1: rango casi seguro > separación típica", r0_pen),
    sprintf("mediana del rango = %.1f km", r0_med),
    sprintf("mediana del rango = %.1f km", r0_dom)))
cat("\n--- Sugerencias de prior.range (PC prior) ---\n"); print(sugerencias)
fwrite(sugerencias, file.path(DIR_TAB, "sugerencias_prior_rango.csv"))

cat(sprintf("\n>> Escala de la red: separación típica ~%.1f km | diámetro ~%.1f km\n",
            median(nn), max(pares)))
cat(sprintf(">> Regla de malla: max.edge del SPDE debe ser < rango/3 (~%.1f km)\n",
            r0_med / 3))

# ==============================================================================
# 5. Figuras: histograma de distancias + mapa de estaciones
# ==============================================================================
g_hist <- ggplot(data.table(d = pares), aes(d)) +
  geom_histogram(bins = 30, fill = "#2166AC", color = "white") +
  geom_vline(xintercept = median(nn), color = "#1A9850", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = median(pares), color = "#B2182B", linewidth = 1) +
  annotate("text", x = median(nn), y = Inf, label = "vecino cercano (mediana)",
           vjust = 2, hjust = -0.05, color = "#1A9850", size = 3) +
  annotate("text", x = median(pares), y = Inf, label = "mediana de pares",
           vjust = 2, hjust = -0.05, color = "#B2182B", size = 3) +
  labs(title = "Distancias entre las 24 estaciones de NO2",
       subtitle = "276 pares | líneas: vecino más cercano (mediana) y mediana de pares",
       x = "Distancia (km)", y = "Nº de pares") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_FIG, "histograma_distancias.png"), g_hist, width = 9, height = 5, dpi = 150)

g_map <- ggplot(est, aes(X_km, Y_km)) +
  geom_point(color = "#2166AC", size = 3) +
  geom_text(aes(label = ESTACION), size = 2.4, vjust = -1, check_overlap = TRUE) +
  coord_equal() +
  labs(title = "Estaciones de NO2 (UTM 30N, km)",
       subtitle = sprintf("Dominio ~ %.0f × %.0f km", ancho_km, alto_km),
       x = "X (km)", y = "Y (km)") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_FIG, "mapa_estaciones.png"), g_map, width = 8, height = 8, dpi = 150)

fwrite(est[, .(ESTACION, LONGITUD, LATITUD, X_km = round(X_km, 3), Y_km = round(Y_km, 3))],
       file.path(DIR_TAB, "estaciones_coordenadas.csv"))

cat("\nFiguras en:", DIR_FIG, "\nTablas en:", DIR_TAB, "\n")
cat("--- Completado ---\n")
