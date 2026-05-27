# 1. Instalar y cargar
if(!require(meteostat)) install.packages("meteostat")
library(meteostat)
library(dplyr)

# 2. Obtener datos de Barajas (08221) - Es la más fiable para datos horarios antiguos
# Usamos 'hourly' para alta frecuencia
madrid_barajas <- hourly("08221", start = "2014-01-01", end = "2018-12-31")

# 3. Obtener datos de Cuatro Vientos (08223)
madrid_4vientos <- hourly("08223", start = "2014-01-01", end = "2018-12-31")

# Ver los datos (vienen limpios y con nombres de columna claros)
head(madrid_barajas)
