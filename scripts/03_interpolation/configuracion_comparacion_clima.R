# 1. Datos: anio 2025, archivos con sufijo 5.
ANIO <- 2025L

# 2. Muestra de lluvia: solo se puntuan objetivos observados > 0.1 mm.
# Los ceros de las estaciones VECINAS se mantienen para hacer la prediccion.
UMBRAL_LLUVIA <- 0.1

# 3. Rango representable de k: conservar al menos el 90% de los objetivos
# evaluables de cada ubicacion. Todos los k usan una muestra comun.
# Este parametro controla cobertura, NO selecciona el k con menor error.
COBERTURA_MINIMA <- 0.90

# 4. TU ELECCION, despues de mirar los graficos RMSE-k.
# Sustituye cada NA por el k que quieras usar en esa variable y escala.
# Ese k se aplica a KNN, IDW p=1 e IDW p=2; Media y 1-NN son referencias.
# Mientras quede NA, la tabla indica "Pendiente" y no inventa una eleccion.
K_ELEGIDO <- data.frame(
  Variable = c(
    "Temperatura", "Humedad_Relativa", "Precipitaciones",
    "Presion_Barometrica", "Radiacion_Solar", "Velocidad_Viento"
  ),
  horario = c(3, 3, 4, 3, 3, 4),
  diario = c(3, 3, 3, 3, 3, 2),
  mensual = c(3, 3, 2, 1, 1, 2)
)
