# ==============================================================================
# MAPAS DE TRÁFICO POR BARRIO — ¿DÓNDE ESTÁ LA MAYOR INTENSIDAD Y CARGA?
# Madrid · Intensidad (veh/h) y Carga (%) · Nivel barrio
# ==============================================================================
# Cortes generados:
#   1. Hora punta mañana  (7-9h)   · 2025 · laborables
#   2. Hora punta tarde   (18-20h) · 2025 · laborables   (siguiente hora punta)
#   3. Hora valle/no punta(4-6h)   · 2025 · laborables   (mínimo de tráfico)
#   4. Anual por barrio — 2019     · media de todos los días
#   5. Anual por barrio — 2025     · media de todos los días
#
# Las ventanas horarias 7-9 / 18-20 / 4-6 se justifican con el perfil horario
# medio de intensidad entre semana (ver bloque 1): mínimo en 5-6h, meseta alta
# sostenida entre 10h y 22h con picos relativos en la mañana y en la tarde.
# ==============================================================================
# QUÉ HACE:
#
# - Carga la cartografía de los barrios de Madrid.
# - Normaliza sus nombres para poder cruzarlos con los datos de tráfico.
# - Calcula puntos interiores de los barrios para colocar etiquetas.
#
# - Carga los datos horarios de tráfico de 2025.
# - Calcula el perfil horario medio de intensidad durante los días laborables.
# - Utiliza ese perfil para justificar tres ventanas horarias:
#     · Punta de mañana: 7:00–9:00.
#     · Punta de tarde: 18:00–20:00.
#     · Periodo valle: 4:00–6:00.
#
# - Calcula para cada barrio:
#     · Intensidad media de tráfico.
#     · Carga media de tráfico.
#     · Número de observaciones disponibles.
#
# - Realiza las agregaciones correspondientes a:
#     · Punta de mañana de 2025.
#     · Punta de tarde de 2025.
#     · Periodo valle de 2025.
#     · Media anual de 2019.
#     · Media anual de 2025.
#
# - Utiliza escalas de color comunes para que las franjas horarias sean
#   comparables entre sí.
# - Utiliza también escalas comunes para comparar 2019 con 2025.
# - Destaca y etiqueta los seis barrios con valores más elevados.
#
# - Genera:
#     · Un gráfico del perfil horario medio.
#     · Un mapa doble para la punta de mañana.
#     · Un mapa doble para la punta de tarde.
#     · Un mapa doble para el periodo valle.
#     · Una figura comparativa anual 2019–2025.
#
# - Cada mapa doble muestra:
#     · Intensidad de tráfico.
#     · Carga de tráfico.
#
# - Genera además una tabla con los diez barrios de mayor intensidad y carga
#   para cada periodo.
#
# FINALIDAD PARA EL TFM:
#
# Este script describe la distribución espacial del tráfico, una de las fuentes
# principales de emisiones de NO2 en Madrid.
#
# Permite identificar:
#
# - Barrios con mayor intensidad de circulación.
# - Barrios con mayor grado de saturación o carga.
# - Diferencias entre horas punta y horas valle.
# - Cambios espaciales del tráfico entre 2019 y 2025.
# - Posibles coincidencias entre zonas de tráfico elevado y concentraciones
#   elevadas de NO2.
#
# Los resultados ayudan a justificar la inclusión de intensidad y carga como
# covariables del modelo.
#
# También permiten interpretar espacialmente los efectos del tráfico y estudiar
# si las estaciones situadas cerca de barrios con tráfico intenso presentan
# mayores concentraciones de contaminación.
#
# OBSERVACIÓN METODOLÓGICA:
#
# La selección de las ventanas horarias no se realiza de forma arbitraria. El
# propio script genera previamente el perfil horario medio para respaldar la
# elección de los periodos de punta y valle.
#
# SALIDAS:
#
# outputs/analysis/mapas/trafico/
#
# Salidas principales:
#
# - 00_perfil_horario_ventanas.png
# - 00_Resumen_Top10_Barrios_Trafico.csv
# - 01_trafico_punta_manana_7_9h.png
# - 02_trafico_punta_tarde_18_20h.png
# - 03_trafico_valle_4_6h.png
# - 04_trafico_anual_2019_vs_2025.png


library(data.table)
library(sf)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(scales)
library(here)

source(here("R", "cleaning", "cleaning_functions.R"))  # limpiar_nombres()

TOP_N <- 6  # nº de barrios etiquetados como "mayor tráfico" en cada mapa

carpeta_out <- here("outputs", "analysis", "mapas", "trafico")
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 0. GEOMETRÍA DE BARRIOS
# ==============================================================================

mapa_barrios <- st_read(here("data", "raw", "geometrias", "BARRIOS.shp"), quiet = TRUE) |>
  st_transform(25830)
mapa_barrios$barrio_norm <- limpiar_nombres(mapa_barrios$NOMBRE)

# Centroides (dentro del polígono) para colocar las etiquetas de texto
pts_centroide <- st_point_on_surface(mapa_barrios)
coords_centroide <- st_coordinates(pts_centroide)
mapa_barrios$X <- coords_centroide[, 1]
mapa_barrios$Y <- coords_centroide[, 2]


# ==============================================================================
# 1. PERFIL HORARIO MEDIO (LABORABLES) — JUSTIFICA LAS VENTANAS ELEGIDAS
# ==============================================================================

dt_horario_2025 <- readRDS(here("data", "processed", "Trafico", "horario",
                                "trafico_madrid_2025_horario_barrio.rds"))
setDT(dt_horario_2025)
dt_horario_2025[, dow := as.integer(format(FECHA, "%u"))]

perfil_horario <- dt_horario_2025[dow <= 5, .(
  intensidad_media = mean(intensidad, na.rm = TRUE)
), by = HORA][order(HORA)]

cat("Perfil horario medio de intensidad (laborables, 2025):\n")
print(perfil_horario)

p_perfil <- ggplot(perfil_horario, aes(x = HORA, y = intensidad_media)) +
  annotate("rect", xmin = 7, xmax = 9,   ymin = -Inf, ymax = Inf, fill = "#c0392b", alpha = 0.12) +
  annotate("rect", xmin = 18, xmax = 20, ymin = -Inf, ymax = Inf, fill = "#e67e22", alpha = 0.12) +
  annotate("rect", xmin = 4, xmax = 6,   ymin = -Inf, ymax = Inf, fill = "#2980b9", alpha = 0.12) +
  geom_line(color = "#2c3e50", linewidth = 1) +
  geom_point(color = "#2c3e50", size = 1.8) +
  scale_x_continuous(breaks = seq(0, 24, 2)) +
  labs(
    title    = "Perfil horario medio de intensidad de tráfico — Madrid 2025",
    subtitle = "Laborables · todos los barrios · rojo: punta mañana (7-9h) · naranja: punta tarde (18-20h) · azul: valle (4-6h)",
    x = "Hora del día", y = "Intensidad media (veh/h)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey40", size = 9),
        panel.grid.minor = element_blank())

ggsave(file.path(carpeta_out, "00_perfil_horario_ventanas.png"),
      plot = p_perfil, width = 10, height = 5.5, dpi = 250, bg = "white")


# ==============================================================================
# 2. AGREGACIÓN DE TRÁFICO POR BARRIO
# ==============================================================================

# --- Media por barrio en una franja horaria concreta (laborables) -----------
agregar_trafico_horario <- function(anio, horas, solo_laborables = TRUE) {
  dt <- readRDS(here("data", "processed", "Trafico", "horario",
                     sprintf("trafico_madrid_%d_horario_barrio.rds", anio)))
  setDT(dt)
  dt <- dt[HORA %in% horas]
  if (solo_laborables) {
    dt[, dow := as.integer(format(FECHA, "%u"))]
    dt <- dt[dow <= 5]
  }
  dt[, barrio_norm := limpiar_nombres(barrio)]
  dt[!is.na(intensidad) & !is.na(carga), .(
    intensidad = mean(intensidad, na.rm = TRUE),
    carga      = mean(carga,      na.rm = TRUE),
    n_obs      = .N
  ), by = barrio_norm]
}

# --- Media anual por barrio (todos los días, escala diaria) -----------------
agregar_trafico_anual <- function(anio) {
  dt <- readRDS(here("data", "processed", "Trafico", "diario",
                     sprintf("trafico_madrid_%d_diario_barrio.rds", anio)))
  setDT(dt)
  dt[, barrio_norm := limpiar_nombres(barrio)]
  dt[!is.na(intensidad) & !is.na(carga), .(
    intensidad = mean(intensidad, na.rm = TRUE),
    carga      = mean(carga,      na.rm = TRUE),
    n_obs      = .N
  ), by = barrio_norm]
}

dt_punta_manana <- agregar_trafico_horario(2025, 7:9)
dt_punta_tarde  <- agregar_trafico_horario(2025, 18:20)
dt_valle        <- agregar_trafico_horario(2025, 4:6)
dt_anual_2019   <- agregar_trafico_anual(2019)
dt_anual_2025   <- agregar_trafico_anual(2025)

# Límites de color comunes: comparables entre las 3 franjas horarias,
# y por separado, comparables entre los 2 años (antes/después)
horarias <- rbindlist(list(dt_punta_manana, dt_punta_tarde, dt_valle))
LIM_INT_HORARIA   <- range(horarias$intensidad, na.rm = TRUE)
LIM_CARGA_HORARIA <- range(horarias$carga,      na.rm = TRUE)

anuales <- rbindlist(list(dt_anual_2019, dt_anual_2025))
LIM_INT_ANUAL   <- range(anuales$intensidad, na.rm = TRUE)
LIM_CARGA_ANUAL <- range(anuales$carga,      na.rm = TRUE)


# ==============================================================================
# 3. FUNCIÓN DE MAPEO: INTENSIDAD + CARGA EN UNA SOLA FIGURA
# ==============================================================================

mapa_metrica <- function(sf_datos, top_datos, columna, etiqueta, unidad, limites, paleta) {
  ggplot(sf_datos) +
    geom_sf(aes(fill = .data[[columna]]), colour = "white", linewidth = 0.08) +
    geom_sf(data = top_datos, fill = NA, colour = "#c0392b", linewidth = 0.8) +
    geom_text_repel(data = top_datos,
                    aes(x = X, y = Y, label = NOMBRE),
                    size = 2.5, fontface = "bold", colour = "#7b241c",
                    bg.color = "white", bg.r = 0.12,
                    box.padding = 0.4, min.segment.length = 0,
                    max.overlaps = 30, seed = 7231) +
    scale_fill_viridis_c(option = paleta, limits = limites, oob = scales::squish,
                         na.value = "grey85",
                         name = paste0(etiqueta, "\n(", unidad, ")")) +
    labs(title = etiqueta) +
    theme_void(base_size = 10) +
    theme(legend.position  = "right",
          plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
          plot.background  = element_rect(fill = "white", colour = NA))
}

generar_mapa_periodo <- function(dt_periodo, titulo, subtitulo, archivo,
                                 lim_intensidad, lim_carga, top_n = TOP_N) {

  sf_p <- merge(mapa_barrios, dt_periodo, by = "barrio_norm", all.x = TRUE)

  top_intensidad <- sf_p[order(-sf_p$intensidad), ][seq_len(top_n), ]
  top_carga      <- sf_p[order(-sf_p$carga), ][seq_len(top_n), ]

  p_int <- mapa_metrica(sf_p, top_intensidad, "intensidad", "Intensidad",
                        "veh/h", lim_intensidad, "plasma")
  p_car <- mapa_metrica(sf_p, top_carga, "carga", "Carga",
                        "%", lim_carga, "magma")

  titulo_grob <- textGrob(
    label = paste0(titulo, "\n", subtitulo),
    gp = gpar(fontface = "bold", fontsize = 13, col = "grey15"),
    just = "left", x = 0.02
  )

  p_final <- arrangeGrob(p_int, p_car, ncol = 2, top = titulo_grob)
  ggsave(file.path(carpeta_out, archivo), plot = p_final,
        width = 13, height = 7.5, dpi = 250, bg = "white")

  cat(sprintf("✅ %s guardado\n", archivo))

  # Ranking impreso en consola (top 5 por métrica)
  cat(sprintf("   Top 5 intensidad: %s\n",
             paste(head(as.data.frame(top_intensidad)$NOMBRE, 5), collapse = ", ")))
  cat(sprintf("   Top 5 carga:      %s\n",
             paste(head(as.data.frame(top_carga)$NOMBRE, 5), collapse = ", ")))

  list(top_intensidad = top_intensidad, top_carga = top_carga)
}


# ==============================================================================
# 4. GENERAR LOS 5 MAPAS
# ==============================================================================

r1 <- generar_mapa_periodo(dt_punta_manana,
     "Hora punta mañana (7h-9h)", "Madrid 2025 · laborables · media por barrio",
     "01_trafico_punta_manana_7_9h.png", LIM_INT_HORARIA, LIM_CARGA_HORARIA)

r2 <- generar_mapa_periodo(dt_punta_tarde,
     "Hora punta tarde (18h-20h)", "Madrid 2025 · laborables · media por barrio",
     "02_trafico_punta_tarde_18_20h.png", LIM_INT_HORARIA, LIM_CARGA_HORARIA)

r3 <- generar_mapa_periodo(dt_valle,
     "Hora valle / no punta (4h-6h)", "Madrid 2025 · laborables · media por barrio",
     "03_trafico_valle_4_6h.png", LIM_INT_HORARIA, LIM_CARGA_HORARIA)

# --- Comparativa anual 2019 vs 2025: una sola figura, 2 años x 2 métricas ---
generar_mapa_anual_comparativo <- function(dt_2019, dt_2025, archivo,
                                           lim_intensidad, lim_carga, top_n = TOP_N) {

  sf_2019 <- merge(mapa_barrios, dt_2019, by = "barrio_norm", all.x = TRUE)
  sf_2025 <- merge(mapa_barrios, dt_2025, by = "barrio_norm", all.x = TRUE)

  top_int_2019 <- sf_2019[order(-sf_2019$intensidad), ][seq_len(top_n), ]
  top_car_2019 <- sf_2019[order(-sf_2019$carga), ][seq_len(top_n), ]
  top_int_2025 <- sf_2025[order(-sf_2025$intensidad), ][seq_len(top_n), ]
  top_car_2025 <- sf_2025[order(-sf_2025$carga), ][seq_len(top_n), ]

  p1 <- mapa_metrica(sf_2019, top_int_2019, "intensidad", "2019 · Intensidad", "veh/h", lim_intensidad, "plasma")
  p2 <- mapa_metrica(sf_2019, top_car_2019, "carga", "2019 · Carga", "%", lim_carga, "magma")
  p3 <- mapa_metrica(sf_2025, top_int_2025, "intensidad", "2025 · Intensidad", "veh/h", lim_intensidad, "plasma")
  p4 <- mapa_metrica(sf_2025, top_car_2025, "carga", "2025 · Carga", "%", lim_carga, "magma")

  titulo_grob <- textGrob(
    label = "Tráfico anual por barrio — 2019 vs 2025\nMedia de todos los días del año · misma escala de color entre años",
    gp = gpar(fontface = "bold", fontsize = 13, col = "grey15"),
    just = "left", x = 0.02
  )

  p_final <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2, top = titulo_grob)
  ggsave(file.path(carpeta_out, archivo), plot = p_final,
        width = 13, height = 14, dpi = 250, bg = "white")

  cat(sprintf("✅ %s guardado\n", archivo))
  cat(sprintf("   Top 5 intensidad 2019: %s\n", paste(head(as.data.frame(top_int_2019)$NOMBRE, 5), collapse = ", ")))
  cat(sprintf("   Top 5 intensidad 2025: %s\n", paste(head(as.data.frame(top_int_2025)$NOMBRE, 5), collapse = ", ")))
  cat(sprintf("   Top 5 carga 2019:      %s\n", paste(head(as.data.frame(top_car_2019)$NOMBRE, 5), collapse = ", ")))
  cat(sprintf("   Top 5 carga 2025:      %s\n", paste(head(as.data.frame(top_car_2025)$NOMBRE, 5), collapse = ", ")))

  list(top_int_2019 = top_int_2019, top_car_2019 = top_car_2019,
       top_int_2025 = top_int_2025, top_car_2025 = top_car_2025)
}

r4 <- generar_mapa_anual_comparativo(dt_anual_2019, dt_anual_2025,
     "04_trafico_anual_2019_vs_2025.png", LIM_INT_ANUAL, LIM_CARGA_ANUAL)


# ==============================================================================
# 5. TABLA RESUMEN — TOP 10 BARRIOS POR PERIODO Y MÉTRICA
# ==============================================================================

periodos <- list(
  "Punta mañana (7-9h) 2025" = dt_punta_manana,
  "Punta tarde (18-20h) 2025" = dt_punta_tarde,
  "Valle (4-6h) 2025"         = dt_valle,
  "Anual 2019"                = dt_anual_2019,
  "Anual 2025"                = dt_anual_2025
)

tabla_resumen <- rbindlist(lapply(names(periodos), function(nm) {
  dt <- merge(as.data.table(mapa_barrios)[, .(barrio_norm, NOMBRE, NOMDIS)],
             periodos[[nm]], by = "barrio_norm")
  top_i <- dt[order(-intensidad)][1:10, .(Periodo = nm, Metrica = "Intensidad",
                                          Barrio = NOMBRE, Distrito = NOMDIS,
                                          Valor = round(intensidad, 1))]
  top_c <- dt[order(-carga)][1:10, .(Periodo = nm, Metrica = "Carga",
                                     Barrio = NOMBRE, Distrito = NOMDIS,
                                     Valor = round(carga, 1))]
  rbind(top_i, top_c)
}))

fwrite(tabla_resumen, file.path(carpeta_out, "00_Resumen_Top10_Barrios_Trafico.csv"))

cat("\n================================================================\n")
cat(" PROCESO COMPLETADO — Mapas de tráfico por barrio\n")
cat(" Salidas en:", carpeta_out, "\n")
cat("================================================================\n")
