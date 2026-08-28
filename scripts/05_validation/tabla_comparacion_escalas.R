# ==============================================================================
# Tabla comparativa de métodos de interpolación: escala x muestra
# ==============================================================================
# Reproduce el formato de la "Tabla 4" del TFM pero cruzando:
#   - Escala:  diaria vs horaria
#   - Muestra: toda la muestra vs un único día
# El ganador (menor RMSE) de cada fila aparece en negrita.
#
# Lee los CSV generados por los scripts de simulación diaria y horaria en sus
# carpetas versionadas más recientes. Si el CSV del horario completo aún no
# existe (solo se exportó como PNG en ejecuciones antiguas), usa como respaldo
# los valores de la última ejecución conocida.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(gridExtra)
  library(grid)
  library(ggplot2)
})

VARIABLES <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  "Presion_Barometrica", "Radiacion_Solar", "Velocidad Viento"
)

ETIQUETA_VAR <- c(
  Temperatura = "Temperatura",
  Humedad_Relativa = "Humedad relativa",
  Precipitaciones = "Precipitaciones",
  Presion_Barometrica = "Presión barom.",
  Radiacion_Solar = "Radiación solar",
  `Velocidad Viento` = "Velocidad viento"
)

# ------------------------------------------------------------------------------
# 1. Localizar las carpetas versionadas más recientes
# ------------------------------------------------------------------------------

ultimo_dir <- function(patron) {
  dirs <- list.dirs(here("outputs", "figures"), recursive = FALSE)
  dirs <- dirs[grepl(patron, basename(dirs))]
  if (length(dirs) == 0L) {
    stop("No se encontró ninguna carpeta que cumpla: ", patron)
  }
  dirs[which.max(file.info(dirs)$mtime)]
}

dir_diario <- ultimo_dir("^interpolacion_clima_[0-9]{8}")
dir_horario <- ultimo_dir("^interpolacion_clima_horaria_[0-9]{8}")

cat("Diario :", basename(dir_diario), "\n")
cat("Horario:", basename(dir_horario), "\n")

# ------------------------------------------------------------------------------
# 2. Cargar las cuatro tablas
# ------------------------------------------------------------------------------

diario_full <- fread(file.path(dir_diario, "tabla_comparacion_metodos.csv"))
diario_dia <- fread(file.path(dir_diario, "tabla_comparacion_metodos_un_dia.csv"))
horario_dia <- fread(file.path(dir_horario, "tabla_comparacion_metodos_un_dia.csv"))

f_horario_full <- file.path(dir_horario, "tabla_comparacion_metodos.csv")
if (file.exists(f_horario_full)) {
  horario_full <- fread(f_horario_full)
} else {
  message("CSV del horario completo no encontrado: se usan valores de respaldo.")
  horario_full <- data.table(
    Variable = VARIABLES,
    Muestra = c(
      "Todas las horas 2025", "Todas las horas 2025", "Horas de lluvia",
      "Todas las horas 2025", "Horas radiación >10 W/m²", "Horas viento >=P75"
    ),
    N_Predicciones = c(195044, 177222, 5021, 68368, 33914, 21690),
    RMSE_Media = c(1.954, 5.762, 0.963, 4.552, 54.799, 0.830),
    RMSE_NN = c(2.505, 5.846, 1.159, 5.693, 64.784, 1.207),
    RMSE_kNN = c(2.008, 5.301, 0.965, 3.909, 54.004, 1.003),
    RMSE_IDW_b1 = c(2.081, 5.161, 0.999, 4.085, 54.328, 1.024),
    RMSE_IDW_b2 = c(2.145, 5.174, 1.041, 4.292, 56.308, 1.058),
    RMSE_Ensemble = c(2.100, 5.112, 0.998, 4.263, 54.963, 1.045),
    Mejor_Metodo = c("Media", "Ensemble", "Media", "kNN", "kNN", "Media")
  )
}

# ------------------------------------------------------------------------------
# 3. Estandarizar a un formato común
# ------------------------------------------------------------------------------

mapear_ganador <- function(x) {
  fifelse(
    x %in% c("Vecino Cercano", "Closest observation"), "1-NN",
    fifelse(
      x == "IDW beta=1", "IDW b1",
      fifelse(x == "IDW beta=2", "IDW b2", x)
    )
  )
}

estandarizar <- function(dt, escala, es_un_dia) {
  dt <- copy(dt)
  if (es_un_dia) {
    muestra <- paste0("Un día ", format(as.Date(dt$Dia), "%d-%b"))
  } else {
    muestra <- dt$Muestra
  }
  out <- data.table(
    Escala = escala,
    Muestra = muestra,
    Variable = unname(ETIQUETA_VAR[dt$Variable]),
    Pred = dt$N_Predicciones,
    Media = dt$RMSE_Media,
    `1-NN` = dt$RMSE_NN,
    kNN = dt$RMSE_kNN,
    `IDW b1` = dt$RMSE_IDW_b1,
    `IDW b2` = dt$RMSE_IDW_b2,
    Ensemble = dt$RMSE_Ensemble,
    Ganador = mapear_ganador(dt$Mejor_Metodo),
    orden = match(dt$Variable, VARIABLES)
  )
  setorder(out, orden)
  out[, orden := NULL]
  out
}

tabla <- rbindlist(list(
  estandarizar(diario_full, "Diaria", FALSE),
  estandarizar(diario_dia, "Diaria", TRUE),
  estandarizar(horario_full, "Horaria", FALSE),
  estandarizar(horario_dia, "Horaria", TRUE)
), use.names = TRUE)

# ------------------------------------------------------------------------------
# 4. Formato de presentación
# ------------------------------------------------------------------------------

cols_metodo <- c("Media", "1-NN", "kNN", "IDW b1", "IDW b2", "Ensemble")

tabla_fmt <- copy(tabla)
tabla_fmt[, Pred := format(Pred, big.mark = ".", trim = TRUE)]
for (col in cols_metodo) {
  tabla_fmt[[col]] <- sprintf("%.2f", tabla[[col]])
}

df_render <- as.data.frame(tabla_fmt)

# ------------------------------------------------------------------------------
# 5. tableGrob con el ganador en negrita
# ------------------------------------------------------------------------------

tg <- tableGrob(
  df_render,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 8,
    core = list(fg_params = list(hjust = 1, x = 0.92)),
    colhead = list(fg_params = list(fontface = "bold"))
  )
)

# Negrita en la celda del método ganador de cada fila.
for (i in seq_len(nrow(df_render))) {
  jcol <- match(df_render$Ganador[i], names(df_render))
  if (!is.na(jcol)) {
    ind <- which(
      tg$layout$name == "core-fg" &
        tg$layout$t == i + 1L &
        tg$layout$l == jcol
    )
    if (length(ind) == 1L) {
      tg$grobs[[ind]] <- editGrob(
        tg$grobs[[ind]],
        gp = gpar(fontface = "bold")
      )
    }
  }
}

# Líneas horizontales separando los cuatro bloques (cada 6 filas).
for (fila_bloque in c(0L, 6L, 12L, 18L, 24L)) {
  tg <- gtable::gtable_add_grob(
    tg,
    grobs = segmentsGrob(
      x0 = unit(0, "npc"), x1 = unit(1, "npc"),
      y0 = unit(0, "npc"), y1 = unit(0, "npc"),
      gp = gpar(lwd = 1.4, col = "grey40")
    ),
    t = fila_bloque + 1L, b = fila_bloque + 1L,
    l = 1, r = ncol(tg)
  )
}

titulo <- grid::textGrob(
  "RMSE por método, variable, escala y muestra (validación 2025)",
  gp = grid::gpar(fontsize = 14, fontface = "bold")
)
subtitulo <- grid::textGrob(
  paste(
    "En negrita, el método ganador (menor RMSE).",
    "Ensemble = pesos 1/RMSE de 1-NN, IDW b1 y kNN.",
    "k=3 | IDW beta=1 y 2."
  ),
  gp = grid::gpar(fontsize = 9)
)

composicion <- arrangeGrob(
  titulo, subtitulo, tg,
  ncol = 1,
  heights = grid::unit(c(0.5, 0.4, 7.2), "in")
)

DIR_SALIDA <- here("outputs", "figures", "comparacion_escalas")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
ruta_png <- file.path(DIR_SALIDA, "tabla_comparacion_escalas.png")
ruta_csv <- file.path(DIR_SALIDA, "tabla_comparacion_escalas.csv")

ggsave(
  ruta_png,
  plot = composicion,
  width = 12,
  height = 9,
  dpi = 200,
  bg = "white"
)
fwrite(tabla, ruta_csv)

cat("\nTabla combinada:\n")
print(tabla)
cat("\nPNG guardada en:\n", ruta_png, "\n")
