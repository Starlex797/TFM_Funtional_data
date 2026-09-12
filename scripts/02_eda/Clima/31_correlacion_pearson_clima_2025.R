# ==============================================================================
# CORRELACION DE PEARSON ENTRE VARIABLES CLIMATOLOGICAS - ANO 2025
# ==============================================================================
#
# Se calcula una tabla para cada escala: horaria, diaria y mensual.
# Como no todas las estaciones miden las mismas variables, cada correlacion se
# calcula con los pares disponibles para esas dos variables.
#
# Salidas:
#   outputs/figures/EDA/correlacion_pearson_clima_2025/
#
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)


# 1. Configuracion -------------------------------------------------------------

ANIO <- 2025
ESCALAS <- c("horario", "diario", "mensual")

VARIABLES <- c(
  "Temperatura",
  "Humedad_Relativa",
  "Precipitaciones",
  "Presion_Barometrica",
  "Radiacion_Solar",
  "Velocidad_Viento"
)

ETIQUETAS <- c(
  Temperatura = "Temperatura",
  Humedad_Relativa = "Humedad relativa",
  Precipitaciones = "Precipitaciones",
  Presion_Barometrica = "Presion barometrica",
  Radiacion_Solar = "Radiacion solar",
  Velocidad_Viento = "Velocidad del viento"
)

CARPETA_SALIDA <- here(
  "outputs", "figures", "EDA", "correlacion_pearson_clima_2025"
)
dir.create(CARPETA_SALIDA, recursive = TRUE, showWarnings = FALSE)


# 2. Funcion para calcular una escala -----------------------------------------

calcular_correlacion <- function(escala) {

  # Archivo de 2025 con sufijo 5.
  archivo <- here(
    "data", "processed", "Clima", escala,
    paste0("meteo_madrid_", ANIO, "_", escala, "5.rds")
  )

  if (!file.exists(archivo)) {
    stop("No se encuentra el archivo: ", archivo)
  }

  datos <- as.data.table(readRDS(archivo))

  if (!all(VARIABLES %in% names(datos))) {
    stop("Faltan variables climatologicas en el archivo ", basename(archivo))
  }

  # Se excluyen los valores imputados o con problemas de calidad.
  # Las versiones con sufijo 5 incluyen una columna de estado por variable.
  for (variable in VARIABLES) {
    columna_estado <- paste0(variable, "_estado")

    if (columna_estado %in% names(datos)) {
      datos[get(columna_estado) != "OK", (variable) := NA_real_]
    }
  }

  clima <- datos[, ..VARIABLES]

  # Correlacion de Pearson usando todos los pares disponibles.
  matriz_cor <- cor(clima, use = "pairwise.complete.obs", method = "pearson")

  # Numero de observaciones utilizadas en cada correlacion.
  matriz_n <- matrix(
    0L,
    nrow = length(VARIABLES),
    ncol = length(VARIABLES),
    dimnames = list(VARIABLES, VARIABLES)
  )

  for (i in seq_along(VARIABLES)) {
    for (j in seq_along(VARIABLES)) {
      matriz_n[i, j] <- sum(
        complete.cases(clima[[VARIABLES[i]]], clima[[VARIABLES[j]]])
      )
    }
  }

  # Guardar las matrices en CSV.
  tabla_cor <- data.table(Variable = ETIQUETAS[rownames(matriz_cor)])
  tabla_cor <- cbind(tabla_cor, as.data.table(round(matriz_cor, 3)))
  setnames(tabla_cor, VARIABLES, ETIQUETAS[VARIABLES])

  tabla_n <- data.table(Variable = ETIQUETAS[rownames(matriz_n)])
  tabla_n <- cbind(tabla_n, as.data.table(matriz_n))
  setnames(tabla_n, VARIABLES, ETIQUETAS[VARIABLES])

  fwrite(
    tabla_cor,
    file.path(CARPETA_SALIDA, paste0("correlacion_pearson_", escala, ".csv"))
  )
  fwrite(
    tabla_n,
    file.path(CARPETA_SALIDA, paste0("numero_pares_", escala, ".csv"))
  )

  # Convertir la matriz a formato largo para dibujarla.
  tabla_grafico <- as.data.table(as.table(matriz_cor))
  setnames(tabla_grafico, c("Variable_1", "Variable_2", "Correlacion"))

  tabla_grafico[, Variable_1 := factor(
    ETIQUETAS[as.character(Variable_1)],
    levels = rev(ETIQUETAS[VARIABLES])
  )]
  tabla_grafico[, Variable_2 := factor(
    ETIQUETAS[as.character(Variable_2)],
    levels = ETIQUETAS[VARIABLES]
  )]

  grafico <- ggplot(
    tabla_grafico,
    aes(x = Variable_2, y = Variable_1, fill = Correlacion)
  ) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", Correlacion)), size = 4) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson r"
    ) +
    coord_equal() +
    labs(
      title = paste("Correlacion de Pearson - escala", escala),
      subtitle = "Variables climatologicas de Madrid, 2025",
      x = NULL,
      y = NULL,
      caption = "Cada correlacion usa los pares disponibles de observaciones con estado OK."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(CARPETA_SALIDA, paste0("correlacion_pearson_", escala, ".png")),
    grafico,
    width = 10,
    height = 8,
    dpi = 200,
    bg = "white"
  )

  # Tabla con las correlaciones altas, sin repetir cada pareja dos veces.
  tabla_pares <- rbindlist(lapply(seq_len(length(VARIABLES) - 1L), function(i) {
    rbindlist(lapply((i + 1L):length(VARIABLES), function(j) {
      data.table(
        Escala = escala,
        Variable_1 = ETIQUETAS[VARIABLES[i]],
        Variable_2 = ETIQUETAS[VARIABLES[j]],
        Pearson_r = matriz_cor[i, j],
        N = matriz_n[i, j]
      )
    }))
  }))

  tabla_pares[, Posible_colinealidad := abs(Pearson_r) >= 0.7]
  tabla_pares[, Pearson_r := round(Pearson_r, 3)]

  fwrite(
    tabla_pares,
    file.path(CARPETA_SALIDA, paste0("pares_correlacion_", escala, ".csv"))
  )

  cat("Escala", escala, "terminada.\n")
  tabla_pares
}


# 3. Ejecutar las tres escalas -------------------------------------------------

resultados <- rbindlist(lapply(ESCALAS, calcular_correlacion))

# Resumen con las parejas que alcanzan |r| >= 0.70.
correlaciones_altas <- resultados[Posible_colinealidad == TRUE]
setorder(correlaciones_altas, Escala, -abs(Pearson_r))

fwrite(
  correlaciones_altas,
  file.path(CARPETA_SALIDA, "resumen_posible_colinealidad.csv")
)

cat("\nResultados guardados en:\n", CARPETA_SALIDA, "\n")
print(correlaciones_altas)
