# ==============================================================================
# Tráfico en las estaciones de NO2: tres métodos y su variación por estación
# ==============================================================================
# Se asigna la intensidad de tráfico a cada estación de contaminación con tres
# métodos y se compara el valor que da cada uno, estación por estación:
#   1. Media de barrio (areal)
#   2. Buffer de 300 m (media de los medidores del entorno)
#   3. Kriging ordinario
# Salidas: mapa de los tres métodos y comparación por estación (variación).
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(gstat)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(scales)
})
sf_use_s2(FALSE)

ANIO <- 2025L
DIA <- as.Date("2025-04-03") # día concreto (jueves laborable)
RADIO_BUFFER <- 300
LIM_INT <- c(0, 1500) # escala común de intensidad para los mapas
DIR_SALIDA <- here("outputs", "figures", "kriging_trafico")
MESES_ES <- c("Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio",
              "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre")
fecha_str <- format(DIA, "%Y-%m-%d")
fecha_txt <- format(DIA, "%d %b %Y")

# ------------------------------------------------------------------------------
# 1. Datos: intensidad media del DÍA concreto por medidor
# ------------------------------------------------------------------------------

mes_nombre <- MESES_ES[as.integer(format(DIA, "%m"))]
ruta_dia <- here("data", "raw", "Datos_trafico",
                 sprintf("Datos_limpios_%d", ANIO), "Horario",
                 sprintf("Trafico_Horario_%s_%d.rds", mes_nombre, ANIO))
dt_raw <- as.data.table(readRDS(ruta_dia))
sensores <- dt_raw[FECHA == DIA & !is.na(intensidad),
  .(intensidad = mean(intensidad)),
  by = .(id, utm_x, utm_y, barrio, distrito)
]
sensores <- sensores[
  !is.na(utm_x) & !is.na(utm_y) &
    utm_x > 400000 & utm_x < 470000 & utm_y > 4450000 & utm_y < 4500000
]
cat(sprintf("Día %s: %d medidores | intensidad media %.0f veh/h\n",
            fecha_str, nrow(sensores), mean(sensores$intensidad)))
sf_sensores <- st_as_sf(sensores, coords = c("utm_x", "utm_y"), crs = 25830)

distritos <- st_transform(st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE
)), 25830)
barrios <- st_transform(st_make_valid(st_read(
  here("data", "raw", "geometrias", "BARRIOS.shp"), quiet = TRUE
)), 25830)

dt_no2 <- as.data.table(readRDS(
  here("data", "processed", "Contaminacion", "horario",
       sprintf("aire_madrid_%d_No2_horarios.rds", ANIO))
))
no2_est <- unique(dt_no2[, .(ESTACION, LONGITUD, LATITUD)])
sf_no2 <- st_transform(
  st_as_sf(no2_est, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830
)

# ------------------------------------------------------------------------------
# 2. Método 1: media de barrio
# ------------------------------------------------------------------------------

idx_sb <- st_intersects(sf_sensores, barrios)
sens_barrio <- vapply(idx_sb, function(i) {
  if (length(i) == 0L) NA_integer_ else i[1]
}, integer(1))
dt_sb <- data.table(intensidad = sensores$intensidad, b = sens_barrio)
media_barrio <- dt_sb[!is.na(b), .(int_barrio = mean(intensidad)), by = b]
barrios$int_barrio <- media_barrio$int_barrio[match(seq_len(nrow(barrios)), media_barrio$b)]

idx_nb <- st_intersects(sf_no2, barrios)
val_barrio <- vapply(idx_nb, function(i) {
  if (length(i) == 0L) NA_real_ else barrios$int_barrio[i[1]]
}, numeric(1))

# ------------------------------------------------------------------------------
# 3. Método 2: buffer de 300 m
# ------------------------------------------------------------------------------

nb_est <- st_is_within_distance(sf_no2, sf_sensores, dist = RADIO_BUFFER)
val_buffer <- vapply(nb_est, function(idx) {
  if (length(idx) == 0L) NA_real_ else mean(sensores$intensidad[idx])
}, numeric(1))

# ------------------------------------------------------------------------------
# 4. Método 3: kriging ordinario
# ------------------------------------------------------------------------------

vario_emp <- variogram(intensidad ~ 1, sf_sensores)
vario_mod <- fit.variogram(vario_emp, vgm(c("Exp", "Sph", "Gau")))
val_krig <- krige(intensidad ~ 1, sf_sensores, sf_no2,
                  model = vario_mod, nmax = 40, debug.level = 0)$var1.pred

# ------------------------------------------------------------------------------
# 5. Superficies para los mapas
# ------------------------------------------------------------------------------

rejilla <- st_as_sf(st_make_grid(distritos, n = c(110, 110), what = "centers"))
rejilla_madrid <- st_intersection(rejilla, st_union(distritos))
coords_g <- st_coordinates(rejilla_madrid)
dx <- diff(sort(unique(coords_g[, 1])))[1]
dy <- diff(sort(unique(coords_g[, 2])))[1]

# Buffer sobre la rejilla.
nb_g <- st_is_within_distance(rejilla_madrid, sf_sensores, dist = RADIO_BUFFER)
surf_buffer <- vapply(nb_g, function(idx) {
  if (length(idx) == 0L) NA_real_ else mean(sensores$intensidad[idx])
}, numeric(1))

# Kriging sobre la rejilla.
surf_krig <- krige(intensidad ~ 1, sf_sensores, rejilla_madrid,
                   model = vario_mod, nmax = 40, debug.level = 0)$var1.pred

df_g <- data.table(X = coords_g[, 1], Y = coords_g[, 2],
                   Buffer = surf_buffer, Kriging = surf_krig)
df_no2 <- data.table(X = st_coordinates(sf_no2)[, 1],
                     Y = st_coordinates(sf_no2)[, 2])
bordes <- as.data.frame(st_coordinates(
  st_cast(st_cast(st_geometry(distritos), "MULTILINESTRING"), "LINESTRING")
))

escala_int <- function(nombre) {
  scale_fill_viridis_c(option = "inferno", trans = "sqrt", direction = -1,
                       limits = LIM_INT, oob = squish, name = nombre)
}
capa_no2 <- list(
  geom_point(data = df_no2, aes(X, Y), color = "black", size = 1.6),
  geom_point(data = df_no2, aes(X, Y), color = "cyan", size = 0.9)
)
tema <- theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        legend.position = "right", legend.key.size = unit(0.4, "cm"))

p_barrio <- ggplot() +
  geom_sf(data = barrios, aes(fill = int_barrio), color = "white", linewidth = 0.1) +
  geom_sf(data = sf_no2, color = "black", size = 1.6) +
  geom_sf(data = sf_no2, color = "cyan", size = 0.9) +
  escala_int("veh/h") +
  labs(title = "1. Media de barrio") + tema

p_buffer <- ggplot() +
  geom_tile(data = df_g, aes(X, Y, fill = Buffer), width = dx, height = dy) +
  geom_path(data = bordes, aes(X, Y, group = L1), color = "white", linewidth = 0.2) +
  capa_no2 + escala_int("veh/h") + coord_equal() +
  labs(title = "2. Buffer 300 m") + tema

p_krig <- ggplot() +
  geom_tile(data = df_g, aes(X, Y, fill = Kriging), width = dx, height = dy) +
  geom_path(data = bordes, aes(X, Y, group = L1), color = "white", linewidth = 0.2) +
  capa_no2 + escala_int("veh/h") + coord_equal() +
  labs(title = "3. Kriging ordinario") + tema

fig_mapas <- (p_barrio | p_buffer | p_krig) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Intensidad de tráfico por los tres métodos — ", fecha_txt),
    subtitle = "Puntos cian = estaciones de NO₂. Escala de intensidad común.",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, color = "grey30")))

ggsave(file.path(DIR_SALIDA, sprintf("mapa_tres_metodos_trafico_%s.png", fecha_str)),
       plot = fig_mapas, width = 16, height = 6.5, dpi = 200, bg = "white")
cat("Mapa de los tres métodos guardado.\n")

# ------------------------------------------------------------------------------
# 6. Comparación por estación: valores y variación
# ------------------------------------------------------------------------------

dt_est <- data.table(
  ESTACION = no2_est$ESTACION,
  Barrio = round(val_barrio, 0),
  Buffer = round(val_buffer, 0),
  Kriging = round(val_krig, 0)
)
m <- as.matrix(dt_est[, .(Barrio, Buffer, Kriging)])
dt_est[, `:=`(
  Media_3 = round(rowMeans(m, na.rm = TRUE), 0),
  Rango = round(apply(m, 1, function(x) diff(range(x, na.rm = TRUE))), 0),
  CV_pct = round(100 * apply(m, 1, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  }), 1)
)]
setorder(dt_est, -Rango)
cat("\nTráfico asignado a cada estación (ordenado por variación entre métodos):\n")
print(dt_est)
fwrite(dt_est, file.path(DIR_SALIDA,
       sprintf("trafico_por_estacion_tres_metodos_%s.csv", fecha_str)))

# Figura: los tres valores por estación (longitud de la línea = desacuerdo).
dl <- melt(dt_est, id.vars = "ESTACION",
           measure.vars = c("Barrio", "Buffer", "Kriging"),
           variable.name = "Metodo", value.name = "Intensidad")
orden_est <- dt_est[order(Rango), ESTACION]
dl[, ESTACION := factor(ESTACION, levels = orden_est)]

p_est <- ggplot(dl, aes(Intensidad, ESTACION)) +
  geom_line(aes(group = ESTACION), color = "grey75", linewidth = 0.8) +
  geom_point(aes(color = Metodo), size = 2.6) +
  scale_color_manual(values = c(
    Barrio = "#1b9e77", Buffer = "#d95f02", Kriging = "#7570b3"
  )) +
  labs(
    title = paste0("Tráfico asignado a cada estación de NO₂ según el método — ", fecha_txt),
    subtitle = paste(
      "Cada línea une los tres valores de una estación.",
      "Cuanto más larga, mayor el desacuerdo entre métodos."
    ),
    x = "Intensidad (veh/h)", y = NULL, color = "Método"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "grey93"))

ggsave(file.path(DIR_SALIDA, sprintf("trafico_por_estacion_comparacion_%s.png", fecha_str)),
       plot = p_est, width = 10, height = 8, dpi = 200, bg = "white")

cat(sprintf(
  "\nVariación entre métodos: rango medio %.0f veh/h | máx %.0f (%s)\n",
  mean(dt_est$Rango, na.rm = TRUE),
  max(dt_est$Rango, na.rm = TRUE),
  dt_est[which.max(Rango), ESTACION]
))
cat("Salidas en:\n", DIR_SALIDA, "\n")
