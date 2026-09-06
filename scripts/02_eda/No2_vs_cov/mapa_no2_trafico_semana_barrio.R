# ==============================================================================
# Mapa: NO2 vs intensidad de trafico POR BARRIO — semana de dias laborables
# ==============================================================================
# Version a nivel de BARRIO (mas fina que distrito) para que la relacion
# trafico-NO2 se aprecie mejor. Para los dias laborables de una semana:
#   - Fondo coropletico: intensidad media de trafico por BARRIO.
#   - Burbujas: NO2 medio por estacion (tamaño y color ~ concentracion).
#   - Panel de dispersion: NO2 de cada estacion vs intensidad de SU barrio.
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

SEMANA_INI <- as.Date("2025-03-10")             # lunes; semana tipica sin festivos
DIR_SALIDA <- here("outputs", "figures", "no2_trafico_semana")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

FESTIVOS_2025 <- as.Date(c(
  "2025-01-01", "2025-01-06", "2025-04-17", "2025-04-18", "2025-05-01",
  "2025-05-02", "2025-05-15", "2025-07-25", "2025-08-15", "2025-10-12",
  "2025-11-01", "2025-11-09", "2025-12-06", "2025-12-08", "2025-12-25"
))

dias <- SEMANA_INI + 0:6
dias_lab <- dias[as.integer(format(dias, "%u")) <= 5 & !(dias %in% FESTIVOS_2025)]
cat("Dias laborables considerados:\n"); print(dias_lab)

# Normaliza nombres de barrio (minusculas, sin acentos ni signos) para el join.
norm_bar <- function(x) {
  x <- tolower(trimws(x)); x <- iconv(x, to = "ASCII//TRANSLIT"); gsub("[^a-z]", "", x)
}

# ------------------------------------------------------------------------------
# 1. Trafico: intensidad media por barrio en los dias laborables de la semana
# ------------------------------------------------------------------------------

f_traf <- list.files(
  here("data", "raw", "Datos_trafico", "Datos_limpios_2025", "Diario_Barrio"),
  pattern = "rds$", full.names = TRUE)
traf <- rbindlist(lapply(f_traf, function(f) as.data.table(readRDS(f))),
                  use.names = TRUE, fill = TRUE)
traf_sem <- traf[FECHA %in% dias_lab, .(
  intensidad = mean(intensidad, na.rm = TRUE),
  ocupacion  = mean(ocupacion,  na.rm = TRUE),
  carga      = mean(carga,      na.rm = TRUE)
), by = barrio]
traf_sem[, key := norm_bar(barrio)]
cat(sprintf("\nBarrios con trafico: %d | intensidad media: %.0f veh/h\n",
            nrow(traf_sem), mean(traf_sem$intensidad)))

# ------------------------------------------------------------------------------
# 2. Geometria de barrios (UTM 25830) + union con el trafico
# ------------------------------------------------------------------------------

bar_sf <- st_make_valid(st_read(
  here("data", "raw", "geometrias", "BARRIOS.shp"), quiet = TRUE))
bar_sf$key <- norm_bar(bar_sf$NOMBRE)
bar_sf <- merge(bar_sf, traf_sem[, .(key, intensidad, ocupacion, carga)],
                by = "key", all.x = TRUE)
sin_traf <- sum(is.na(bar_sf$intensidad))
if (sin_traf) cat("Barrios del shp sin trafico (grises):", sin_traf, "\n")

# ------------------------------------------------------------------------------
# 3. NO2: media por estacion en los mismos dias laborables + barrio de cada una
# ------------------------------------------------------------------------------

no2 <- as.data.table(readRDS(here(
  "data", "processed", "Contaminacion", "diario",
  "aire_madrid_2025_No2_trans_diarios.rds")))
no2_sem <- no2[FECHA %in% dias_lab,
               .(NO2 = mean(DATO_DIARIO, na.rm = TRUE)),
               by = .(ESTACION, LONGITUD, LATITUD)][!is.na(NO2)]
sf_no2 <- st_transform(
  st_as_sf(no2_sem, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830)

# Cada estacion hereda la intensidad de su barrio (join espacial).
sf_no2 <- st_join(sf_no2, bar_sf[, c("NOMBRE", "intensidad", "ocupacion", "carga")],
                  join = st_within)
setDT(no2_sem)
no2_sem[, `:=`(barrio = sf_no2$NOMBRE, intensidad = sf_no2$intensidad,
               ocupacion = sf_no2$ocupacion, carga = sf_no2$carga)]
no2_rel <- no2_sem[!is.na(intensidad)]

# ------------------------------------------------------------------------------
# 4. Mapa + dispersion para cada metrica de trafico
#    - intensidad: lo pedido (veh/h). NO se relaciona con el NO2.
#    - carga: congestion. SI se relaciona (la intensidad no distingue autovia
#      fluida de calle atascada; la congestion si).
# ------------------------------------------------------------------------------

etq_sem <- sprintf("%s a %s", format(min(dias_lab), "%d/%m"),
                   format(max(dias_lab), "%d/%m/%Y"))

METRICAS <- list(
  list(col = "intensidad", suf = "barrio",
       ejeuds = "Intensidad de tráfico del barrio (veh/h)",
       leyenda = "Tráfico\n(veh/h)", trans = "sqrt",
       titulo = "NO₂ e intensidad de tráfico en Madrid, por barrio"),
  list(col = "carga", suf = "barrio_carga",
       ejeuds = "Carga / congestión del barrio (%)",
       leyenda = "Carga\n(%)", trans = "identity",
       titulo = "NO₂ y congestión de tráfico (carga) en Madrid, por barrio")
)

for (mt in METRICAS) {
  v <- mt$col
  r_p <- cor(no2_rel$NO2, no2_rel[[v]], use = "complete.obs")
  r_s <- cor(no2_rel$NO2, no2_rel[[v]], method = "spearman", use = "complete.obs")
  cat(sprintf("\n[%s] NO2 ~ %s (n=%d): Pearson=%.2f | Spearman=%.2f\n",
              mt$suf, v, nrow(no2_rel), r_p, r_s))

  p_mapa <- ggplot() +
    geom_sf(data = bar_sf, aes(fill = .data[[v]]), color = "white", linewidth = 0.15) +
    geom_sf(data = sf_no2, aes(size = NO2), shape = 21,
            fill = "#d73027", color = "white", stroke = 0.6, alpha = 0.9) +
    scale_fill_viridis_c(option = "mako", direction = -1, name = mt$leyenda,
                         na.value = "grey88", trans = mt$trans) +
    scale_size_continuous(range = c(2, 9), name = "NO₂\n(µg/m³)") +
    labs(title = mt$titulo,
         subtitle = sprintf("Días laborables %s. Fondo = tráfico por barrio; burbujas = NO₂ por estación.", etq_sem)) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey30"),
          legend.position = "right")

  p_disp <- ggplot(no2_rel, aes(.data[[v]], NO2)) +
    geom_smooth(method = "lm", se = FALSE, color = "#d73027", linewidth = 0.8) +
    geom_point(color = "#2166ac", size = 2.6, alpha = 0.85) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("Pearson r = %.2f\nSpearman ρ = %.2f", r_p, r_s),
             hjust = -0.1, vjust = 1.3, fontface = "bold", size = 3.6) +
    labs(title = "Relación NO₂ — tráfico (barrio)",
         subtitle = "Cada punto es una estación (tráfico de su barrio)",
         x = mt$ejeuds, y = "NO₂ medio (µg/m³)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey30"))

  fig <- p_mapa + p_disp + plot_layout(widths = c(1.6, 1))
  ggsave(file.path(DIR_SALIDA, sprintf("mapa_no2_trafico_semana_%s.png", mt$suf)),
         fig, width = 15, height = 7.5, dpi = 200, bg = "white")
}

fwrite(no2_rel[order(-NO2)], file.path(DIR_SALIDA, "no2_trafico_estaciones_barrio.csv"))
cat("\nSalidas guardadas en:\n", DIR_SALIDA, "\n")
