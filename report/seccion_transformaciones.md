# Análisis exploratorio de las relaciones covariable–respuesta y justificación de las transformaciones

## 1. Objetivo

Antes de especificar el modelo espacio-temporal es necesario decidir **la forma en que
entran las variables**: si la variable respuesta (la concentración de NO₂) debe modelarse
en su escala original o transformada, y si las covariables mantienen una relación
aproximadamente lineal con la respuesta o requieren alguna transformación o un tratamiento
no lineal. Esta sección documenta el análisis gráfico que fundamenta esas decisiones.

El análisis se apoya en diagramas de dispersión (*scatterplots*) del NO₂ frente a cada una
de las covariables candidatas —tráfico (intensidad y carga), temperatura, humedad relativa,
precipitación, presión barométrica, radiación solar y velocidad del viento—, a resolución
**horaria** durante el año 2025.

## 2. Estrategia: estratificación por tipología de estación

Un diagrama de dispersión que agregue todas las estaciones de la ciudad mezcla dos fuentes
de variación muy distintas: la variación **entre estaciones** (cada emplazamiento tiene un
nivel base de NO₂ diferente según su entorno) y la variación **dentro de cada estación** (la
respuesta de la contaminación a las condiciones de tráfico y meteorología a lo largo del
tiempo). Esa mezcla difumina las relaciones y puede inducir a error (falacia ecológica): una
nube difusa en el conjunto agregado no implica que la relación sea débil dentro de cada
estación.

Para aislar la relación *dentro de la estación* se estratifica por la **tipología** de la
estación de medida, que la red de vigilancia de la calidad del aire de Madrid clasifica en
tres clases:

- **Suburbana** (representada por *Casa de Campo*),
- **Urbana de fondo** (representada por *Retiro*),
- **Urbana de tráfico** (representada por *Plaza Elíptica*).

Se selecciona una estación representativa de cada tipología y se comparan, para cada
covariable, dos versiones de la respuesta: el **NO₂ sin transformar** y su
**transformación logarítmica** `log(NO₂ + 1)`. Al fijar la estación desaparece la variación
entre emplazamientos y la relación con cada covariable se lee con mucha más nitidez
(Figuras 1–6).

## 3. La variable respuesta: necesidad de la transformación logarítmica

En las tres tipologías, la distribución del **NO₂ en su escala original es marcadamente
asimétrica a la derecha**: la masa de observaciones se concentra en valores bajos
(aproximadamente 0–50 µg/m³) y se prolonga en una cola de episodios de contaminación alta
(hasta ~150 µg/m³ en la estación de tráfico, ~100 µg/m³ en fondo y ~120 µg/m³ en suburbana).
Esta asimetría es la esperable en una variable de concentración: es no negativa y se genera
por procesos **multiplicativos** (acumulación y dispersión que actúan en proporción, no en
términos aditivos).

La transformación `log(NO₂ + 1)` **simetriza** la nube de puntos y estabiliza la dispersión
a lo largo del rango de las covariables. Esto es precisamente lo que exige la hipótesis de
**errores gaussianos homocedásticos** sobre la que se construye el modelo (la verosimilitud
gaussiana de INLA). El término `+1` evita el problema de `log(0)` en las horas con NO₂ nulo.

> **Decisión.** Se modela la respuesta en escala logarítmica, `log(NO₂ + 1)`. Es la elección
> estándar para concentraciones de contaminantes y queda respaldada visualmente por la
> simetrización observada en todas las tipologías. Esta decisión se confirma con los
> diagnósticos de residuos del modelo ajustado (histograma y gráfico cuantil-cuantil de los
> residuos, que deben aproximarse a la normalidad) y, formalmente, con un análisis de
> Box–Cox, cuyo parámetro óptimo se espera próximo a λ ≈ 0.

**Un matiz importante.** El diagrama de dispersión marginal informa sobre todo de la
**linealidad de cada covariable**; la decisión sobre la transformación de la *respuesta*
depende de la distribución de los **residuos** del modelo, no de la nube cruda. Ambas piezas
apuntan en la misma dirección (a favor del logaritmo), pero conviene distinguir los dos
objetivos: transformar la respuesta corrige la asimetría y la heterocedasticidad de los
errores; transformar un predictor corrige la no-linealidad de su relación con la respuesta.

## 4. Las covariables: forma de la relación

Del examen de las facetas por covariable, comunes a las tres tipologías, se extraen los
siguientes patrones.

### 4.1 Relaciones fuertes y no lineales

- **Velocidad del viento.** Es la covariable meteorológica con la señal más nítida: relación
  claramente **decreciente y curva** (forma hiperbólica), coherente con el mecanismo físico de
  **dispersión** —la concentración disminuye aproximadamente como el inverso de la velocidad
  del viento—. La caída es abrupta a vientos bajos y se aplana a vientos altos.
- **Radiación solar.** Relación **decreciente y curva**, con una acumulación de puntos en el
  valor cero (horas nocturnas). La radiación actúa además como *proxy* del ciclo diurno: de
  noche (radiación nula) el NO₂ se acumula, y de día la mezcla atmosférica y la fotoquímica lo
  reducen.
- **Tráfico (intensidad y carga).** Relación **creciente** pero con forma de cuña: el tráfico
  fija el *techo* del NO₂ (a mayor tráfico, mayores máximos posibles) más que su valor exacto,
  de modo que la dispersión vertical es amplia. La intensidad presenta además una distribución
  sesgada a la derecha.

En estos tres casos la relación **no es una recta**. La transformación logarítmica de la
respuesta atenúa la curvatura pero **no la elimina**: en los paneles de `log(NO₂ + 1)` el
viento y la radiación siguen mostrando una relación curva. Por tanto, la no-linealidad debe
tratarse en el lado de los **predictores**.

### 4.2 Relaciones moderadas o aproximadamente lineales

- **Temperatura.** Relación débil, con ligera curvatura (mayores concentraciones en el rango
  frío, compatible con la mayor estabilidad atmosférica invernal). Admite un tratamiento
  lineal, a lo sumo con un término cuadrático.
- **Presión barométrica.** Relación **creciente**: la presión alta se asocia a situaciones
  anticiclónicas de estancamiento que favorecen la acumulación de contaminantes. Es
  aproximadamente monótona y no requiere transformación.

### 4.3 Relaciones débiles o con estructura especial

- **Humedad relativa.** Nube **difusa** en las tres tipologías, incluso tras aislar la
  estación. Que la dispersión persista *dentro* de un mismo emplazamiento confirma que se
  trata de una covariable genuinamente débil, y no de un artefacto de agregación.
- **Precipitación.** Distribución fuertemente **inflada en cero** (la inmensa mayoría de las
  horas no registra lluvia). Este problema **no se resuelve con una transformación de escala**:
  la variable es casi degenerada. Su tratamiento natural sería, en su caso, un **indicador
  binario** (hora con o sin precipitación) más que una covariable continua, o bien su descarte
  como predictor de efecto apreciable.

### 4.4 Un artefacto a tener en cuenta

En los paneles de la respuesta logarítmica aparecen **bandas horizontales** discretas. No
constituyen señal: son consecuencia de que el NO₂ se registra en valores enteros de µg/m³, de
modo que `log(entero + 1)` genera un conjunto discreto de niveles. Es un artefacto de la
resolución de medida y no afecta a la interpretación de las relaciones.

## 5. Diferencias entre tipologías

La estratificación revela que **el peso del tráfico local depende de la tipología**, lo que
constituye el hallazgo con mayor relevancia para la especificación del modelo:

- En la estación de **tráfico** (*Plaza Elíptica*) los niveles de NO₂ son los más elevados y
  el tráfico del propio emplazamiento gobierna claramente el techo de la concentración: el
  sensor "respira" su vía.
- En la estación de **fondo** (*Retiro*) y en la **suburbana** (*Casa de Campo*) el tráfico
  inmediato es bajo y de rango estrecho, y explica poco de la variación; la concentración
  responde sobre todo a la **meteorología** (viento, radiación, presión) y refleja la
  contaminación de fondo urbano y el transporte regional.

En cambio, el efecto de las variables **meteorológicas** es **consistente entre tipologías**
(viento y radiación siempre reducen el NO₂; la presión lo aumenta), lo que sugiere que la
componente atmosférica es relativamente homogénea a escala de ciudad, mientras que la
componente de emisión es marcadamente **local**.

Esta heterogeneidad espacial del efecto del tráfico justifica que el modelo no se limite a
efectos fijos globales, sino que incorpore un **campo espacial** (que absorbe el nivel base
distinto de cada punto) y la **tipología/distrito** de la estación: el mismo valor de tráfico
no significa lo mismo en un cañón de tráfico que en un parque.

## 6. Síntesis y decisiones adoptadas

| Elemento | Diagnóstico gráfico | Decisión |
|---|---|---|
| **Respuesta NO₂** | Asimetría derecha; el log simetriza (todas las tipologías) | Modelar `log(NO₂ + 1)`; confirmar con residuos y Box–Cox |
| **Viento** | Decreciente, curva (≈ 1/viento) | Efecto **no lineal** (p. ej. *spline* `rw2`) o transformación |
| **Radiación solar** | Decreciente, curva, con exceso de ceros (noche) | Efecto **no lineal**; controlar el ciclo horario |
| **Tráfico (intensidad/carga)** | Creciente, sesgado, en cuña | Lineal; valorar `log` para interpretar el coeficiente |
| **Temperatura** | Débil, leve curvatura | Lineal (o término cuadrático) |
| **Presión** | Creciente, monótona | Lineal, sin transformación |
| **Humedad** | Difusa dentro de la estación | Covariable **débil**; efecto lineal pequeño |
| **Precipitación** | Inflada en cero | No transformar; indicador binario o descarte |

**Conclusión.** El análisis exploratorio estratificado por tipología justifica **una única
transformación de la variable respuesta** —el logaritmo— y aconseja tratar la no-linealidad
de las covariables meteorológicas dominantes (viento y radiación) mediante **términos no
lineales** en lugar de transformaciones *ad hoc* de escala, opción que en el marco INLA-SPDE
resulta más flexible y evita fijar a priori una forma funcional. Las covariables débiles
(humedad, precipitación) no mejoran con ninguna transformación, lo que anticipa su papel
secundario en el modelo. Estas decisiones, tomadas sobre la evidencia visual, se ratifican a
posteriori con los diagnósticos de residuos y los criterios de ajuste del modelo estimado.

---

### Figuras

- **Figura 1.** NO₂ sin transformar frente a las covariables — *Casa de Campo* [Suburbana].
- **Figura 2.** `log(NO₂ + 1)` frente a las covariables — *Casa de Campo* [Suburbana].
- **Figura 3.** NO₂ sin transformar frente a las covariables — *Retiro* [Urbana de fondo].
- **Figura 4.** `log(NO₂ + 1)` frente a las covariables — *Retiro* [Urbana de fondo].
- **Figura 5.** NO₂ sin transformar frente a las covariables — *Plaza Elíptica* [Urbana de tráfico].
- **Figura 6.** `log(NO₂ + 1)` frente a las covariables — *Plaza Elíptica* [Urbana de tráfico].

*Fuente: elaboración propia. Escala horaria, Madrid 2025. Cada punto es una hora; el eje
horizontal es libre en cada faceta. Archivos en*
`outputs/scatter_no2_covariables/panel_estacion/<estacion>/`.
