# ==============================================================================
# CORRELACION DE PEARSON ENTRE VARIABLES CLIMATOLOGICAS, POR ESTACION
# Escala diaria - Ano 2025
# ==============================================================================
#
# El script 31 calcula UNA matriz de correlacion juntando todas las estaciones.
# Este la calcula POR ESTACION, para cinco de ellas, y las dibuja una al lado de
# otra. La pregunta que responde es distinta: no "como se relacionan las
# variables en Madrid", sino "esa relacion es la misma en todos los puntos de la
# ciudad o cambia segun donde midas".
#
# Importa para el modelo INLA-SPDE: si la estructura de correlacion entre
# covariables fuera muy distinta de una estacion a otra, un unico coeficiente
# fijo por covariable para toda la ciudad seria una simplificacion discutible.
#
# ESTACIONES. Se eligen las cinco pedidas. Ojo con Moratalaz: en los datos hay
# dos estaciones con ese nombre y solo "J.M.D Moratalaz" tiene los seis
# sensores; la llamada "Moratalaz" a secas solo mide temperatura y humedad, asi
# que no se puede correlacionar nada con ella. Se usa "J.M.D Moratalaz".
#
# Salidas: outputs/figures/EDA/correlacion_clima_estaciones_2025/
# ==============================================================================

library(data.table)
library(ggplot2)
library(here)


# 1. Configuracion -------------------------------------------------------------

ANIO <- 2025
ESCALA <- "diario"

ESTACIONES <- c(
  "J.M.D Moratalaz",
  "J.M.D Hortaleza",
  "J.M.D Villaverde",
  "Casa de Campo",
  "Peñagrande"
)

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
  "outputs", "figures", "EDA", "correlacion_clima_estaciones_2025"
)
dir.create(CARPETA_SALIDA, recursive = TRUE, showWarnings = FALSE)


# 2. Datos ---------------------------------------------------------------------

archivo <- here(
  "data", "processed", "Clima", ESCALA,
  paste0("meteo_madrid_", ANIO, "_", ESCALA, "5.rds")
)
if (!file.exists(archivo)) stop("No se encuentra el archivo: ", archivo)

datos <- as.data.table(readRDS(archivo))

faltan <- setdiff(ESTACIONES, unique(datos$ESTACION))
if (length(faltan) > 0) {
  stop("Estaciones no encontradas en el fichero: ", paste(faltan, collapse = ", "))
}

datos <- datos[ESTACION %in% ESTACIONES]

# Se anulan los valores que no estan marcados como OK. Las versiones con sufijo
# 5 traen una columna de estado por variable (OK / IMPUTADO / FALLO / AUSENTE /
# SIN_SENSOR). Correlacionar valores imputados inflaria la correlacion, porque
# la imputacion se construye precisamente a partir de las otras variables.
for (variable in VARIABLES) {
  columna_estado <- paste0(variable, "_estado")
  if (columna_estado %in% names(datos)) {
    datos[get(columna_estado) != "OK", (variable) := NA_real_]
  }
}


# 3. Una matriz de correlacion por estacion ------------------------------------

# pairwise.complete.obs: cada pareja de variables usa los dias en que las dos
# estan disponibles, en lugar de tirar el dia entero si falla un sensor.
correlaciones <- rbindlist(lapply(ESTACIONES, function(estacion) {

  clima <- datos[ESTACION == estacion, VARIABLES, with = FALSE]
  matriz_cor <- cor(clima, use = "pairwise.complete.obs", method = "pearson")

  tabla <- as.data.table(as.table(matriz_cor))
  setnames(tabla, c("Variable_1", "Variable_2", "Correlacion"))
  tabla[, ESTACION := estacion]

  # Numero de dias que sostiene cada correlacion: sin esto, un 0.9 calculado
  # con 20 dias pareceria tan solido como uno calculado con 360.
  tabla[, N := mapply(
    function(v1, v2) sum(complete.cases(clima[[v1]], clima[[v2]])),
    as.character(Variable_1), as.character(Variable_2)
  )]
  tabla[]
}))

fwrite(
  correlaciones[, .(ESTACION, Variable_1, Variable_2,
                    Pearson_r = round(Correlacion, 3), N)],
  file.path(CARPETA_SALIDA, "correlaciones_por_estacion.csv")
)


# 4. Heatmap -------------------------------------------------------------------

correlaciones[, Variable_1 := factor(
  ETIQUETAS[as.character(Variable_1)], levels = rev(ETIQUETAS[VARIABLES])
)]
correlaciones[, Variable_2 := factor(
  ETIQUETAS[as.character(Variable_2)], levels = ETIQUETAS[VARIABLES]
)]
correlaciones[, ESTACION := factor(ESTACION, levels = ESTACIONES)]

grafico <- ggplot(correlaciones, aes(Variable_2, Variable_1, fill = Correlacion)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Correlacion)), size = 2.9) +
  # Escala fija en [-1, 1] y centrada en 0: sin esto, cada panel se colorearia
  # con su propio rango y los paneles no serian comparables entre si, que es
  # justo lo que se quiere mirar.
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), name = "Pearson r"
  ) +
  facet_wrap(~ESTACION, ncol = 3) +
  coord_equal() +
  labs(
    title = "Correlacion de Pearson entre variables climatologicas, por estacion",
    subtitle = paste0("Madrid, ", ANIO, " - escala ", ESCALA),
    x = NULL, y = NULL,
    caption = paste(
      "Solo observaciones con estado OK.",
      "Cada correlacion usa los dias con las dos variables disponibles."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(CARPETA_SALIDA, "heatmap_correlacion_clima_estaciones.png"),
  grafico, width = 14, height = 9.5, dpi = 200, bg = "white"
)


# 5. Cuanto varia cada correlacion entre estaciones ----------------------------

# La figura se lee de un vistazo, pero para citar un numero en la memoria hace
# falta el rango: si el recorrido de r entre las cinco estaciones es pequeno,
# la estructura de correlacion es homogenea en la ciudad.
pares <- correlaciones[as.integer(Variable_2) < (length(VARIABLES) + 1L - as.integer(Variable_1))]

dispersion <- pares[, .(
  r_min = round(min(Correlacion), 3),
  r_max = round(max(Correlacion), 3),
  recorrido = round(max(Correlacion) - min(Correlacion), 3),
  N_min = min(N)
), by = .(Variable_1, Variable_2)]
setorder(dispersion, -recorrido)

fwrite(dispersion, file.path(CARPETA_SALIDA, "dispersion_entre_estaciones.csv"))

cat("\n--- Parejas cuya correlacion mas cambia entre estaciones ---\n")
print(head(dispersion, 10))

cat("\nGuardado en:", CARPETA_SALIDA, "\n")
