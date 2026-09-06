# ==============================================================================
# PASO 1: CREACIÓN DE LA MALLA ESPACIAL (MESH)
# ==============================================================================

library(INLA)
library(sf)
library(data.table)
library(here)

# 1. Cargar datos de NO2 (variable respuesta)
dt_no2_2025 <- readRDS(here("data", "processed", "Contaminacion", "diario", "aire_madrid_2025_No2_trans_diarios1.rds"))
setDT(dt_no2_2025)

# 2. Extraer y proyectar coordenadas espaciales
# Se extraen las coordenadas únicas y se proyectan a UTM 30N (EPSG:25830) en kilómetros.
coords_estaciones <- unique(dt_no2_2025[, .(ESTACION, LONGITUD, LATITUD)])
coords_sf <- st_as_sf(coords_estaciones, coords = c("LONGITUD", "LATITUD"), crs = 4326) # CRS original en WGS84 (grados)
coords_utm <- st_transform(coords_sf, 25830) # Proyección a UTM zona 30N (metros)

# Matriz de coordenadas en kilómetros (requerido por INLA para estabilidad numérica)
coords_matriz <- st_coordinates(coords_utm) / 1000
carpeta_figuras_mallas <- here("outputs", "figures", "modelo", "mallas")
if (!dir.exists(carpeta_figuras_mallas)) dir.create(carpeta_figuras_mallas, recursive = TRUE)

# Mapa de distritos en km (mismo CRS que coords_matriz)
mapa_distritos_km <- st_transform(
  st_read(here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE),
  25830
)
# 3. Definir los límites geométricos (Boundaries)
# Se genera un contorno interno ajustado a las estaciones y un contorno externo de amortiguación.
# La función convex define cómo se dibuja la línea de la frontera alrededor de las estaciones.
# bnd_inner <- as.inla.mesh.segment(mapa_distritos_km) # Contorno interno ajustado a las estaciones
bnd_inner <- inla.nonconvex.hull(coords_matriz, convex = -0.3) # Se deja un 5% de margen adicional para incluir todas las estaciones dentro del contorno. Esto ayuda a asegurar que el modelo espacial tenga en cuenta todas las ubicaciones de las estaciones sin ser demasiado ajustado.
bnd_outer <- inla.nonconvex.hull(coords_matriz, convex = -0.6) # Se deja un 20% de margen adicional para amortiguación fuera de las estaciones. Esto ayuda a evitar problemas de borde en el modelado espacial.

# 4. Construir las Mallas (Meshes) de diferentes resoluciones
# max.edge en km, c(interior, exterior). Se declaran aquí una sola vez para que
# los títulos de las figuras no se desincronicen de los parámetros reales.
edge_gruesa <- c(5, 6)
edge_media <- c(3, 6)
edge_fina <- c(1.86, 6)

# OPCIÓN A: Malla Gruesa (Baja resolución, muy rápida)
malla_gruesa <- inla.mesh.2d(
  loc = coords_matriz,
  boundary = list(bnd_inner, bnd_outer), # Se definen los contornos
  max.edge = edge_gruesa,
  cutoff = 0.79
)

# OPCIÓN B: Malla Media (Tu configuración original)

malla_media <- inla.mesh.2d(
  loc = coords_matriz,
  boundary = list(bnd_inner, bnd_outer),
  max.edge = edge_media,
  cutoff = 0.79
)

# OPCIÓN C: Malla Fina (Alta resolución, más lenta)

malla_fina <- inla.mesh.2d(
  loc = coords_matriz,
  boundary = list(bnd_inner, bnd_outer),
  max.edge = edge_fina,
  cutoff = 0.79
) # distancia minima entre nodo. En km
# Por ejemplo, esto se usa cuando si dos estaciones están muy pegadas ,
# la triangulación crearía triangulos diminutos y alargados entre ellas. Entonces lo que
# hace es fusionar esos dos puntos en uno solo, para que no haya triangulos diminutos.


# 5. Exportar los tres archivos al disco duro
saveRDS(malla_gruesa, here("data", "processed", "malla_spde_madrid_gruesa.rds"))
saveRDS(malla_media, here("data", "processed", "malla_spde_madrid_media.rds"))
saveRDS(malla_fina, here("data", "processed", "malla_spde_madrid_fina.rds"))


# Visualizar cuántos vértices (nodos) tiene cada una para compararlas en el TFM
cat("Vértices Malla Gruesa:", malla_gruesa$n, "\n")
cat("Vértices Malla Media:", malla_media$n, "\n")
cat("Vértices Malla Fina:", malla_fina$n, "\n")

# ==============================================================================
# VISUALIZACIÓN: ESTACIONES DE NO2 SOBRE LAS MALLAS SPDE
# ==============================================================================
library(ggplot2)
library(fmesher)
library(gridExtra)

carpeta_figuras_mallas <- here("outputs", "figures", "modelo", "mallas")
if (!dir.exists(carpeta_figuras_mallas)) dir.create(carpeta_figuras_mallas, recursive = TRUE)

# Mapa de distritos en km (mismo CRS que coords_matriz)
st_geometry(mapa_distritos_km) <- st_geometry(mapa_distritos_km) / 1000

# Estaciones como data.frame en km
coords_df <- as.data.frame(coords_matriz)
colnames(coords_df) <- c("X", "Y")
coords_df$ESTACION <- coords_estaciones$ESTACION

# Paleta y etiquetas compartidas
pal_leyenda <- c(
  "NO\u2082 station" = "#d73027",
  "SPDE mesh"        = "#2166ac",
  "Madrid districts" = "#8c510a" # brown (BrBG), contrasts with the mesh blue
)

# Función reutilizable: una malla + estaciones superpuestas + leyenda
crear_mapa_malla <- function(malla, titulo, subtitulo = NULL) {
  # Datos vacíos para registrar "SPDE mesh" en la leyenda (geom_fm no lo hace)
  df_dummy <- data.frame(x = NA_real_, y = NA_real_)

  ggplot() +
    geom_sf(
      data = mapa_distritos_km, fill = NA,
      aes(color = "Madrid districts"), linewidth = 0.4
    ) +
    # geom_fm no registra aes(color) en la escala -> color fijo
    geom_fm(data = malla, color = "#2166ac", linewidth = 0.2, alpha = 0.45) +
    # Capa invisible que crea la entrada de la malla en la leyenda
    geom_line(data = df_dummy, aes(x = x, y = y, color = "SPDE mesh")) +
    geom_point(
      data = coords_df, aes(
        x = X, y = Y,
        color = "NO\u2082 station"
      ),
      shape = 17, size = 2.5
    ) +
    scale_color_manual(name = NULL, values = pal_leyenda) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape     = c(NA, 17, NA),
          linetype  = c(1, 0, 1),
          linewidth = c(0.6, 0, 0.6),
          size      = c(NA, 2.5, NA)
        )
      )
    ) +
    labs(
      title = titulo,
      subtitle = subtitulo %||% sprintf("%d vertices | UTM 30N (km)", malla$n),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "white", color = "gray85"),
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(color = "gray40", size = 9),
      axis.text        = element_text(size = 7, color = "gray50"),
      legend.position  = "bottom",
      legend.text      = element_text(size = 9)
    )
}

# Operador nulo-coalescencia (disponible en R ≥ 4.4, pero lo definimos por si acaso)
`%||%` <- function(a, b) if (!is.null(a)) a else b

mapa_gruesa <- crear_mapa_malla(malla_gruesa, sprintf("Coarse mesh  (max.edge = %g km)", edge_gruesa[1]))
mapa_media <- crear_mapa_malla(malla_media, sprintf("Medium mesh  (max.edge = %g km)", edge_media[1]))
mapa_fina <- crear_mapa_malla(malla_fina, sprintf("Fine mesh    (max.edge = %g km)", edge_fina[1]))

# Guardado con verificación.
# El device PNG por defecto (ragg) NO da error si el fichero está bloqueado
# —p. ej. abierto en un visor o en la vista previa del IDE—: solo emite el aviso
# "agg could not write to the given file" y sigue. El resultado es una figura
# vieja en disco que no refleja los parámetros actuales. Aquí se comprueba que
# el fichero se ha reescrito de verdad y, si no, se para con un mensaje claro.
guardar_png <- function(nombre, grafico, width, height) {
  ruta <- file.path(carpeta_figuras_mallas, nombre)
  mtime_previo <- if (file.exists(ruta)) file.mtime(ruta) else as.POSIXct(NA)

  ggsave(ruta, grafico, width = width, height = height, dpi = 200, bg = "white")

  escrito <- file.exists(ruta) &&
    (is.na(mtime_previo) || file.mtime(ruta) > mtime_previo)
  if (!escrito) {
    stop(
      "No se pudo reescribir '", nombre, "'. Suele ser porque el PNG está ",
      "abierto en un visor o en la vista previa del IDE: ciérralo y vuelve a ",
      "ejecutar. Ruta: ", ruta
    )
  }
  cat("  guardado:", nombre, "\n")
  invisible(ruta)
}

# Guardar cada malla individualmente
guardar_png("malla_gruesa.png", mapa_gruesa, width = 7, height = 6)
guardar_png("malla_media.png", mapa_media, width = 7, height = 6)
guardar_png("malla_fina.png", mapa_fina, width = 7, height = 6)

# Guardar mapa comparativo (3 paneles)
comparacion <- arrangeGrob(mapa_gruesa, mapa_media, mapa_fina, ncol = 3)
guardar_png("comparacion_mallas.png", comparacion, width = 18, height = 6)

cat("Mapas guardados en:", carpeta_figuras_mallas, "\n")
grid.arrange(mapa_gruesa, mapa_media, mapa_fina, ncol = 3)
