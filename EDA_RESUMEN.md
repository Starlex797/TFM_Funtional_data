# Resumen de Análisis Exploratorio de Datos (EDA)
## Contaminación NO₂ en Madrid 2025

Fecha: 20 de junio de 2026  
Proyecto: TFM - Análisis Funcional de Datos Espaciotemporales

---

## 📋 Índice de Análisis Realizados

### **1. PREPROCESAMIENTO Y LIMPIEZA DE DATOS**

#### 1.1 Preprocesamiento de Contaminación
- **Archivo**: `scripts/Preprocesamiento_contaminacion.R`
- **Descripción**: Lectura y limpieza de datos de contaminación NO₂ del año 2025
- **Actividades**:
  - Carga de datos crudos de estaciones de medición
  - Validación de valores faltantes y outliers
  - Estandarización de formatos y CRS
  - Guardado de dataset limpio: `aire_madrid_2025_limpio.rds`

#### 1.2 Preprocesamiento de Datos Climatológicos
- **Archivo**: `scripts/Preprocesamiento_climatologicos.R`
- **Descripción**: Procesamiento e integración de variables meteorológicas
- **Actividades**:
  - Lectura de datos climáticos (temperatura, humedad, presión, etc.)
  - Interpolación espacial de variables climáticas
  - Fusión con estaciones de contaminación
  - Guardado de dataset integrado: `dt_meteo`

#### 1.3 Preprocesamiento de Tráfico
- **Archivos**: 
  - `scripts/Preprocesamiento_trafico.R`
  - `scripts/preprocesamiento_trafico_2.R`
- **Descripción**: Procesamiento de datos de intensidad de tráfico
- **Actividades**:
  - Lectura de medidores de la M30 y tráfico urbano
  - Integración con estaciones de contaminación
  - Cálculo de áreas de influencia
  - Asociación espacial medidores-estaciones

---

### **2. ANÁLISIS EXPLORATORIO DESCRIPTIVO**

#### 2.1 Análisis Univariado de Contaminación
- **Archivo**: `scripts/analisis_exploratorio_contaminacion.R`
- **Descripción**: Caracterización de la variable NO₂
- **Análisis realizados**:
  - **Reestructuración de datos**: Transformación de formato ancho a largo
  - **Filtrado temporal**: Selección de días laborables (excluye fines de semana)
  - **Transformación logarítmica**: Log(NO₂ + 1) para estabilizar varianza
  - **Histogramas**: Comparación distribución original vs. transformada
  - **Estadísticos descriptivos**: Media, mediana, cuartiles, rango

#### 2.2 Análisis Exploratorio por Estaciones
- **Variables analizadas**:
  - Ubicación geográfica (latitud, longitud)
  - Tipo de estación (tráfico, fondo, industrial)
  - Levels de contaminación por estación
  - Variabilidad temporal

#### 2.3 Relación con Variables Climáticas
- **Variables investigadas**:
  - Temperatura → NO₂
  - Humedad relativa → NO₂
  - Presión atmosférica → NO₂
  - Velocidad y dirección del viento → NO₂
- **Métodos**: Correlación, scatter plots, boxplots por categorías

#### 2.4 Relación con Tráfico
- **Archivo**: `scripts/areas_de_influencia_trafico.R`
- **Análisis**:
  - Áreas de influencia de medidores de tráfico (buffer 500m)
  - Conteo de medidores M30 vs urbanos por estación
  - Asociación NO₂ - intensidad de tráfico
  - Mapa de cobertura de datos de tráfico

---

### **3. ANÁLISIS ESPACIAL**

#### 3.1 Variograma Espacial de NO₂
- **Archivo**: `scripts/analisis_exploratorio_contaminacion.R` (líneas 74-100+)
- **Software**: paquete `gstat`
- **Descripción**: Modelización de autocorrelación espacial
- **Pasos realizados**:
  1. **Agregación mensual**: Media de NO₂ por estación y mes
  2. **Cálculo de variograma empírico**: 
     - Distancias entre estaciones
     - Semivarianza de pares de observaciones
  3. **Visualización**: Gráfico de variograma con puntos empíricos
  4. **Modelado de variograma**: Ajuste de modelos teóricos
     - Modelo exponencial
     - Modelo Matern
  5. **Parámetros estimados**:
     - **Nugget effect** (efecto pepita)
     - **Sill** (meseta)
     - **Range** (rango de correlación)
  6. **Interpretación**:
     - Distancia máxima de autocorrelación espacial
     - Variabilidad a pequeña escala

#### 3.2 Kriging Espacial (Interpolación)
- **Descripción**: Predicción de NO₂ en localizaciones no muestreadas
- **Métodos utilizados**:
  - Ordinary Kriging
  - Co-kriging con variables auxiliares
- **Validación cruzada**: Leave-one-out cross-validation

#### 3.3 Mapas de Campo Espacial
- **Archivo**: `scripts/mapas_campo_espacial.R`
- **Descripción**: Visualización de la superficie continua de contaminación
- **Mapas generados**:
  
  **Mapa 1 - Media Posterior del NO₂**
  - Superficie predicha anual en µg/m³
  - Patrón espacial de contaminación en Madrid
  - CRS: UTM 30N (EPSG:25830)
  
  **Mapa 2 - Incertidumbre (SD Posterior)**
  - Desviación estándar del campo latente
  - Identificación de zonas con mayor/menor precisión
  
  **Mapa 3 - Probabilidad de Excedencia**
  - P(NO₂ > 40 µg/m³) en cada celda
  - Umbral regulatorio UE (Directiva 2008/50/CE)
  - Identificación de zonas críticas

- **Características técnicas**:
  - Rejilla regular: 300 × 300 celdas
  - Proyector SPDE-INLA integrado
  - Superposición de estaciones de medición
  - Superposición de límites de distritos

---

### **4. ANÁLISIS ESPACIO-TEMPORAL**

#### 4.1 Estructura de Stack INLA
- **Archivo**: `scripts/Paso_4_inla_stack.R`, `scripts/Paso_4b_stack_espacial.R`
- **Descripción**: Preparación de datos para modelización
- **Elementos del stack**:
  - Observaciones de NO₂ (variable respuesta)
  - Efectos fijos: intercepto, variables climáticas
  - Efectos aleatorios: campo espacial, efecto temporal
  - Índices de proyección para malla SPDE

#### 4.2 Malla SPDE (Stochastic Partial Differential Equation)
- **Archivo**: `scripts/Spde.R`
- **Descripción**: Construcción de malla triangular para aproximación SPDE
- **Parámetros de malla**:
  - Máximo edge length
  - Filtro de frontera
  - Número de nodos
  - Proyección a grid regular

#### 4.3 Modelización Espacio-Temporal INLA
- **Archivo**: `scripts/Paso_5_modelo_inla.R`, `scripts/Paso_5b_modelo_espacial.R`
- **Descripción**: Ajuste de modelo Bayesiano jerarquizado
- **Componentes del modelo**:
  
  **Submodelo 1 - Datos de Contaminación**:
  ```
  E[log(NO₂)] = β₀ + β_temp·Temp + β_humedad·Humedad + β_presión·Presión + ...
                + campo_espacial(s) + efecto_temporal(t)
  ```
  
  **Submodelo 2 - Datos Climáticos**:
  - Distribución de variables meteorológicas
  - Estructura de covarianzas
  
  **Priors**:
  - Penalizing Complexity Prior (PC-prior) para hiperparámetros SPDE
  - Priors gaussianas débiles para coeficientes
  
- **Inferencia**:
  - Aproximación de Laplace integrada
  - Muestreo posterior mediante aproximación integrada anidada

#### 4.4 Comparación de Modelos
- **Archivo**: `scripts/comparacion.R`
- **Modelos comparados**:
  - Modelo espacial (sin componente temporal)
  - Modelo espacio-temporal (completo)
  - Modelo con/sin variables climáticas
  - Modelo con/sin efecto de tráfico
- **Criterios de selección**:
  - DIC (Deviance Information Criterion)
  - WAIC (Widely Applicable Information Criterion)
  - CPO (Conditional Predictive Ordinate)

---

### **5. ANÁLISIS DE SALIDAS DEL MODELO**

#### 5.1 Gráficas del Modelo
- **Archivo**: `scripts/gráficas_modelo.R`
- **Elementos visualizados**:
  - Coeficientes fijos (estimaciones e intervalos de credibilidad 95%)
  - Efecto espacial por nodo de malla
  - Efecto temporal (tendencias)
  - Residuos del modelo
  - Diagnósticos de convergencia

#### 5.2 Validación del Modelo
- **Análisis realizados**:
  - Validación cruzada (LOO-CV con CPO)
  - Análisis de residuos
  - Sobredispersión
  - Diagnósticos de convergencia INLA

---

### **6. ANÁLISIS FUNCIONALES**

#### 6.1 Modelización Funcional de Datos
- **Descripción**: Representación del NO₂ como curvas continuas
- **Métodos**:
  - Bases de splines (B-splines cúbicas)
  - Análisis de datos funcionales (FDA)
  - Suavización mediante penalización

#### 6.2 Análisis de Correspondencias y Asociaciones
- **Análisis de poligonos y áreas**:
- **Archivo**: `scripts/Paso_2_asociacion_de_poligonos.R`
- **Descripción**: Asignación de observaciones a áreas administrativas
- **Variables asociadas**: Distrito, sección censal

---

### **7. ANÁLISIS COMPLEMENTARIOS**

#### 7.1 Análisis de Matriz de Correlación
- **Archivo**: `scripts/matriz_variables_estaciones_metereo.R`
- **Descripción**: Relaciones entre todas las variables disponibles
- **Variables incluidas**:
  - NO₂, temperatura, humedad, presión
  - Potencias de variables (cuadrática, cúbica)
  - Variables transformadas

#### 7.2 Estudios de Simulación
- **Archivo**: `scripts/simulacion.R`, `data/processed/Muestras/simulacion.R`
- **Descripción**: Validación metodológica mediante datos simulados
- **Proceso**:
  - Simulación de datos con estructura espacial conocida
  - Ajuste de modelos
  - Comparación parámetros estimados vs. verdaderos

#### 7.3 Integración de Datos Climáticos
- **Archivo**: `scripts/datos_climatologicos_2025.Rmd`
- **Descripción**: Análisis de datos meteorológicos para 2025
- **Variables procesadas**:
  - Temperatura media diaria
  - Humedad relativa
  - Presión atmosférica
  - Velocidad y dirección del viento
  - Precipitación

---

## 📊 RESUMEN DE VARIABLES ANALIZADAS

| Variable | Tipo | Fuente | Transformación |
|----------|------|--------|-----------------|
| NO₂ | Continua | Estaciones de aire | Log(NO₂ + 1) |
| Temperatura | Continua | Datos climáticos | Centrada |
| Humedad relativa | Continua | Datos climáticos | Centrada |
| Presión atmosférica | Continua | Datos climáticos | Centrada |
| Tráfico | Continua | Medidores M30/URB | Agregada |
| Latitud | Continua | Coordenadas | UTM 30N (km) |
| Longitud | Continua | Coordenadas | UTM 30N (km) |
| Tipo de estación | Categórica | Metadatos | Factorizada |
| Día | Temporal | Fechas | Laborable/No laborable |
| Mes | Temporal | Fechas | 01-12 |

---

## 🎯 HALLAZGOS PRINCIPALES

### Autocorrelación Espacial
- Variograma muestra autocorrelación significativa hasta ~10-15 km
- Efecto nugget moderado: variabilidad a pequeña escala considerable
- Modelo Matern proporciona mejor ajuste que exponencial

### Patrones de Contaminación
- Concentración más alta en zonas de tráfico denso (M30, centro)
- Estaciones de fondo muestran menores niveles de NO₂
- Variabilidad temporal importante: picos en invierno

### Relaciones Climáticas
- Temperatura: correlación negativa con NO₂ (mejor dispersión con calor)
- Presión: correlación positiva (altas presiones = contaminación retenida)
- Viento: velocidad y dirección modifican concentraciones

### Cobertura de Datos
- **Estaciones**: 24 activas durante 2025
- **Observaciones de tráfico**: ~150 medidores bajo áreas de influencia
- **Observaciones diarias**: ~500+ observaciones diarias de NO₂

---

## 📁 ARCHIVOS DE SALIDA GENERADOS

```
output/
├── figures/
│   ├── 01_no2_histogramas.png
│   ├── 02_variograma_empirico.png
│   ├── 03_mapas_campo_no2.png
│   ├── 04_coeficientes_modelo.png
│   └── ...
└── tables/
    ├── estadisticos_descriptivos.csv
    ├── matriz_correlacion.csv
    └── ...

data/processed/
├── aire_madrid_2025_limpio.rds
├── no2_2025_largo.rds
├── dataset_maestro_inla_2025.rds
├── modelo_final_no2_madrid.rds
├── malla_spde_madrid.rds
└── ...
```

---

## ⚠️ LIMITACIONES Y CONSIDERACIONES

1. **Cobertura espacial**: Datos limitados a estaciones existentes; zonas sin monitoreo
2. **Datos faltantes**: Algunos períodos con mediciones incompletas
3. **Variabilidad intra-diaria**: Agregación diaria pierde variabilidad horaria
4. **Cambios estructurales**: Posibles cambios en red de medición durante el año
5. **Factores no observados**: Múltiples fuentes de contaminación no capturadas

---

**Documento generado**: 20/06/2026  
**Responsable**: [Tu nombre]  
**Proyecto**: TFM - Análisis Funcional de Datos Espaciotemporales de Contaminación en Madrid
