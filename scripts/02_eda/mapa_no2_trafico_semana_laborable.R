# ==============================================================================
# Mapa: NO2 vs intensidad de trafico — una semana de dias laborables
# ==============================================================================
# Relaciona la contaminacion (NO2) con el trafico en una semana concreta,
# considerando solo los dias laborables (Lun-Vie, sin festivos):
#   - Fondo coropletico: intensidad media de trafico por DISTRITO.
#   - Burbujas: NO2 medio por estacion (tamaño y color ~ concentracion).
#   - Panel de dispersion: NO2 de cada estacion vs intensidad de su distrito,
#     con correlaciones de Pearson y Spearman.
# Trafico agregado por distrito (veh/h); NO2 de las estaciones de calidad del aire.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})
sf_use_s2(FALSE)

# --- Parametros: semana a analizar (lunes de la semana) -----------------------
SEMANA_INI <- as.Date("2025-03-10")             # lunes; semana tipica sin festivos
DIR_SALIDA <- here("outputs", "figures", "no2_trafico_semana")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

FESTIVOS_2025 <- as.Date(c(
  "2025-01-01", "2025-01-06", "2025-04-17", "2025-04-18", "2025-05-01",
  "2025-05-02", "2025-05-15", "2025-07-25", "2025-08-15", "2025-10-12",
  "2025-11-01", "2025-11-09", "2025-12-06", "2025-12-08", "2025-12-25"
))

# Dias laborables de la semana (Lun-Vie menos festivos).
dias <- SEMANA_INI + 0:6
dias_lab <- dias[as.integer(format(dias, "%u")) <= 5 & !(dias %in% FESTIVOS_2025)]
cat("Dias laborables considerados:\n"); print(dias_lab)

# Normaliza nombres de distrito (minusculas, sin acentos ni signos) para el join.
norm_dist <- function(x) {
  x <- tolower(trimws(x))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^a-z]", "", x)
  x <- gsub("sanblascanillejas", "sanblas", x)   # trafico trae el nombre largo
  x
}

# ------------------------------------------------------------------------------
# 1. Trafico: intensidad media por distrito en los dias laborables de la semana
# ------------------------------------------------------------------------------

f_traf <- list.files(
  here("data", "raw", "Datos_trafico", "Datos_limpios_2025", "Diario"),
  pattern = "\\.rds$", full.names = TRUE)
traf <- rbindlist(lapply(f_traf, function(f) as.data.table(readRDS(f))),
                  use.names = TRUE, fill = TRUE)
traf_sem <- traf[FECHA %in% dias_lab,
                 .(intensidad = mean(intensidad, na.rm = TRUE)), by = distrito]
traf_sem[, key := norm_dist(distrito)]
cat(sprintf("\nDistritos con trafico: %d | intensidad media: %.0f veh/h\n",
            nrow(traf_sem), mean(traf_sem$intensidad)))

# ------------------------------------------------------------------------------
# 2. NO2: media por estacion en los mismos dias laborables
# ------------------------------------------------------------------------------

no2 <- as.data.table(readRDS(here(
  "data", "processed", "Contaminacion", "diario",
  "aire_madrid_2025_No2_trans_diarios.rds")))
no2_sem <- no2[FECHA %in% dias_lab,
               .(NO2 = mean(DATO_DIARIO, na.rm = TRUE)),
               by = .(ESTACION, LONGITUD, LATITUD)][!is.na(NO2)]
sf_no2 <- st_as_sf(no2_sem, coords = c("LONGITUD", "LATITUD"), crs = 4326)
cat(sprintf("Estaciones NO2: %d | NO2 medio: %.1f ug/m3\n",
            nrow(no2_sem), mean(no2_sem$NO2)))

# ------------------------------------------------------------------------------
# 3. Geometria de distritos + union con el trafico
# ------------------------------------------------------------------------------

dist_sf <- st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE))
dist_sf$key <- norm_dist(dist_sf$name)
dist_sf <- merge(dist_sf, traf_sem[, .(key, distrito, intensidad)],
                 by = "key", all.x = TRUE)
sin_match <- dist_sf$name[is.na(dist_sf$intensidad)]
if (length(sin_match)) cat("Distritos sin trafico casado:", paste(sin_match, collapse=", "), "\n")

# Cada estacion NO2 hereda la intensidad de su distrito (join espacial).
sf_no2 <- st_join(sf_no2, dist_sf[, c("name", "intensidad")], join = st_within)
setDT(no2_sem)
no2_sem[, intensidad := sf_no2$intensidad]
no2_sem[, distrito := sf_no2$name]
xy <- st_coordinates(st_centroid(sf_no2))   # solo para etiquetas si hiciera falta
no2_rel <- no2_sem[!is.na(intensidad)]

r_p <- cor(no2_rel$NO2, no2_rel$intensidad, use = "complete.obs")
r_s <- cor(no2_rel$NO2, no2_rel$intensidad, method = "spearman", use = "complete.obs")
cat(sprintf("\nRelacion NO2 ~ intensidad (n=%d): Pearson=%.2f | Spearman=%.2f\n",
            nrow(no2_rel), r_p, r_s))

# ------------------------------------------------------------------------------
# 4a. Mapa: coropletico de trafico + burbujas de NO2
# ------------------------------------------------------------------------------

etq_sem <- sprintf("%s a %s", format(min(dias_lab), "%d/%m"),
                   format(max(dias_lab), "%d/%m/%Y"))

p_mapa <- ggplot() +
  geom_sf(data = dist_sf, aes(fill = intensidad), color = "white", linewidth = 0.3) +
  geom_sf(data = sf_no2, aes(size = NO2), shape = 21,
          fill = "#d73027", color = "white", stroke = 0.6, alpha = 0.9) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Tráfico\n(veh/h)",
                       na.value = "grey85") +
  scale_size_continuous(range = c(2, 9), name = "NO₂\n(µg/m³)") +
  labs(title = "NO₂ e intensidad de tráfico en Madrid",
       subtitle = sprintf("Días laborables %s. Fondo = tráfico por distrito; burbujas = NO₂ por estación.", etq_sem)) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey30"),
        legend.position = "right")

# ------------------------------------------------------------------------------
# 4b. Dispersion: NO2 vs intensidad del distrito
# ------------------------------------------------------------------------------

p_disp <- ggplot(no2_rel, aes(intensidad, NO2)) +
  geom_smooth(method = "lm", se = FALSE, color = "#d73027", linewidth = 0.8) +
  geom_point(color = "#2166ac", size = 2.6, alpha = 0.85) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Pearson r = %.2f\nSpearman ρ = %.2f", r_p, r_s),
           hjust = -0.1, vjust = 1.3, fontface = "bold", size = 3.6) +
  labs(title = "Relación NO₂ — tráfico",
       subtitle = "Cada punto es una estación (tráfico de su distrito)",
       x = "Intensidad de tráfico del distrito (veh/h)", y = "NO₂ medio (µg/m³)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey30"))

fig <- p_mapa + p_disp + plot_layout(widths = c(1.6, 1))
ggsave(file.path(DIR_SALIDA, "mapa_no2_trafico_semana.png"), fig,
       width = 15, height = 7.5, dpi = 200, bg = "white")

fwrite(no2_rel[order(-NO2)], file.path(DIR_SALIDA, "no2_trafico_estaciones.csv"))
cat("\nSalidas guardadas en:\n", DIR_SALIDA, "\n")
