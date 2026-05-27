# scripts/02_eda_and_plots.R
# Objetivo: Transformación a formato largo de NO2, filtrado de días laborables y gráficos.

# 1. Cargar librerías -----------------------------------------------------
library(here)
library(data.table)

# 2. Cargar datos procesados del Script 01 --------------------------------
ruta_entrada <- here("data", "processed", "aire_madrid_2025_limpio.rds")

if (!file.exists(ruta_entrada)) {
  stop("❌ El archivo de datos procesados de 2025 no existe. Ejecuta primero el script 01.")
}

dt_2025 <- readRDS(ruta_entrada)

# 3. Reestructuración de Datos (De Ancho a Largo) -------------------------
cols_horas <- paste0("H", sprintf("%02d", 1:24))

dt_largo <- melt(dt_2025, 
                 id.vars = c("ESTACION", "LONGITUD", "LATITUD", "NOM_TIPO", "MAGNITUD", "FECHA"),
                 measure.vars = cols_horas,
                 variable.name = "HORA",
                 value.name = "NIVEL_CONTAMINACION")

# Filtramos solo NO2 y eliminamos registros sin datos (las horas sin código "V")
dt_objetivo <- dt_largo[MAGNITUD == "NO2" & !is.na(NIVEL_CONTAMINACION)]

# wday(FECHA) devuelve el día de la semana (1 = Domingo, 7 = Sábado)
dt_objetivo[, DIA_SEMANA := wday(FECHA)]

# Filtramos: nos quedamos solo con días laborables (excluimos 1 y 7)
dt_objetivo <- dt_objetivo[!(DIA_SEMANA %in% c(1, 7))]
dt_objetivo[, DIA_SEMANA := NULL] # Borramos columna auxiliar

# Aplicamos la transformación logarítmica para estabilizar la varianza
dt_objetivo[, LOG_NIVEL := log(NIVEL_CONTAMINACION + 1)]
# Visualizacion de la tabla 
View(dt_objetivo)
# Guardamos este dataset en formato largo; lo usaremos en el paso de geoestadística
saveRDS(dt_objetivo, here("data", "processed", "no2_2025_largo.rds"))


# 4. Generación y Guardado de Gráficos ------------------------------------

# --- 4.1 HISTOGRAMAS ---
# Abrimos un dispositivo PNG para guardar la imagen con buena resolución
png(here("output", "figures", "01_no2_histogramas.png"), width = 1000, height = 500, res = 120)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# Histograma Original
hist(dt_objetivo$NIVEL_CONTAMINACION,
     col = "tomato",
     main = "Histograma de NO2\n(Datos Originales)",
     xlab = expression("NO"[2]*" ("*mu*"g/m"^3*")"),
     ylab = "Frecuencia",
     breaks = 50)

# Histograma Transformado
hist(dt_objetivo$LOG_NIVEL,
     col = "lightblue",
     main = "Histograma de Ln(NO2)\n(Datos Transformados)",
     xlab = "Ln(NO2 + 1)",
     ylab = "Frecuencia",
     breaks = 50)

dev.off() # Cerramos el dispositivo para escribir el archivo en el disco


# --- 4.2 BOXPLOTS ---
png(here("output", "figures", "02_no2_boxplots.png"), width = 1200, height = 600, res = 120)

par(mfrow = c(1, 2), mar = c(8, 4, 4, 2) + 0.1) # Más margen inferior para los nombres de las estaciones

# Boxplot Datos Originales
boxplot(NIVEL_CONTAMINACION ~ ESTACION, data = dt_objetivo,
        col = "tomato",
        main = "Niveles de NO2 (Originales)",
        ylab = expression("NO"[2]*" ("*mu*"g/m"^3*")"), 
        las = 2,           
        cex.axis = 0.6,    
        outline = TRUE)    

# Boxplot Datos Transformados
boxplot(LOG_NIVEL ~ ESTACION, data = dt_objetivo,
        col = "lightblue",
        main = "Niveles de Ln(NO2) (Transformados)",
        ylab = "Ln(NO2 + 1)", 
        las = 2, 
        cex.axis = 0.6,
        outline = TRUE)

par(mfrow = c(1, 1)) # Resetear los márgenes gráficos por defecto
dev.off() 

cat("✅ Script 02 completado con éxito. Gráficos guardados en output/figures/\n")

#----------------------------------------------------------------------


# scripts/03_spatial_analysis.R
# Objetivo: Agregación mensual y modelización de variogramas espaciales (Enero 2025).

# 1. Cargar librerías -----------------------------------------------------
library(here)
library(data.table)
library(sf)
library(gstat)

# 2. Cargar datos en formato largo del paso anterior ----------------------
ruta_entrada <- here("data", "processed", "no2_2025_largo.rds")

if (!file.exists(ruta_entrada)) {
  stop("❌ El archivo 'no2_2025_largo.rds' no existe. Ejecuta primero el script 02.")
}

dt_objetivo <- readRDS(ruta_entrada)

# 3. Mensualización (El "Día Laborable Típico") ---------------------------
dt_objetivo[, MES_ANIO := as.Date(trunc(FECHA, "months"))]

datos_mensuales <- dt_objetivo[, .(
  NO2_MEDIA = mean(NIVEL_CONTAMINACION, na.rm = TRUE),
  NO2_LOG   = log(mean(NIVEL_CONTAMINACION, na.rm = TRUE) + 1)
), by = .(ESTACION, LONGITUD, LATITUD, NOM_TIPO, MES_ANIO)]


# 4. Análisis Geoestadístico Espacial (Enero 2025) -----------------------
dt_enero <- datos_mensuales[MES_ANIO == as.Date("2025-01-01")]

# Conversión a objeto espacial sf (Coordenadas geográficas WGS84)
sf_enero <- st_as_sf(dt_enero, coords = c("LONGITUD", "LATITUD"), crs = 4326)

# Proyección a coordenadas métricas (UTM Zona 30N para Madrid)
sf_enero_utm <- st_transform(sf_enero, crs = 25830)


# --- 4.1 Nube de Variograma ---
nube_variograma <- gstat::variogram(NO2_LOG ~ 1, data = sf_enero_utm, cloud = TRUE)

# Guardamos el gráfico de la nube directamente en formato PNG
png(here("output", "figures", "03_no2_nube_variograma.png"), width = 700, height = 500, res = 120)

# Al usar la función plot nativa sobre un objeto de gstat, se genera la gráfica correctamente
print(plot(nube_variograma,
           main = "Nube de Variograma - Ln(NO2) Enero 2025",
           col = "darkorange", pch = 16,
           xlab = "Distancia de separación (metros)",
           ylab = "Semivarianza (Diferencias al cuadrado / 2)"))

dev.off()


# --- 4.2 Variograma Empírico y Ajuste Teórico ---
var_empirico <- gstat::variogram(NO2_LOG ~ 1, data = sf_enero_utm)

# Ajuste del modelo teórico Exponencial a partir del empírico
var_ajustado_v2 <- gstat::fit.variogram(
  var_empirico,
  model = gstat::vgm(psill = 0.04, model = "Exp", range = 8000, nugget = 0.01)
)

cat("\n📊 Parámetros estimados del Variograma Ajustado:\n")
print(var_ajustado_v2)

# Guardamos la comparación del variograma empírico vs ajustado en PNG
png(here("output", "figures", "04_no2_variograma_ajustado.png"), width = 700, height = 500, res = 120)

print(plot(var_empirico, var_ajustado_v2,
           main = "Variograma Espacial de Ln(NO2) — Enero 2025 (Laborables)",
           col  = "darkred",
           xlab = "Distancia de separación (metros)",
           ylab = "Semivarianza"))

dev.off()

cat("✅ Script 03 completado con éxito. Variogramas calculados y guardados.\n")

