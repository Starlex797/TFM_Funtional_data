# Limitaciones del análisis de correlación y del VIF: no estacionariedad y dependencia de la escala

## Planteamiento

Un aspecto metodológico que conviene explicitar es que **ni el coeficiente de
correlación ni el Factor de Inflación de la Varianza (VIF) son medidas absolutas
y estables**: ambos son *resúmenes estadísticos calculados sobre un conjunto de
datos concreto*, y su valor cambia con la escala temporal de agregación, con la
época del año y con la estación de medida considerada. Ignorar esto lleva a
conclusiones frágiles sobre qué covariables importan. A continuación se explican
las carencias de cada herramienta y por qué las correlaciones de las variables
climáticas con el NO₂ resultan tan inestables.

## Carencias del VIF como criterio de selección

El VIF cuantifica la redundancia de un predictor a partir del coeficiente de
determinación de su regresión **lineal** sobre el resto de predictores
(`VIF = 1/(1 − R²)`). De esta definición se derivan sus principales limitaciones:

1. **Detecta únicamente colinealidad lineal.** Dos variables ligadas de forma no
   lineal —como la temperatura y la radiación solar a través del ciclo diurno—
   pueden presentar un VIF bajo y ser, sin embargo, redundantes. El VIF no
   captura ese tipo de dependencia.
2. **No considera la variable respuesta.** El VIF mide la redundancia *entre
   predictores*, no su utilidad para explicar el NO₂. Una covariable puede tener
   VIF bajo y ser irrelevante, o VIF alto y ser físicamente esencial; eliminar
   por VIF puede descartar una variable importante.
3. **El umbral es una regla empírica** (habitualmente 5 ó 10), sin una base
   teórica que lo fundamente.
4. **No captura la confusión estructural del ciclo compartido.** Temperatura,
   radiación y humedad comparten un mismo motor —el ciclo diurno y estacional—;
   esa redundancia mediada por el ciclo no queda bien resumida en un único valor
   de VIF.
5. **No es invariante: depende de la muestra.** El VIF calculado a resolución
   horaria difiere del calculado a resolución diaria, y ambos difieren del de un
   único mes o estación del año. No es una propiedad intrínseca de las variables,
   sino un retrato de ese subconjunto concreto.
6. **Asume observaciones independientes.** Se calcula sobre una regresión
   ordinaria que ignora la autocorrelación espacial y temporal de los datos, por
   lo que su información efectiva está sobreestimada.

En consecuencia, el VIF es un **filtro lineal, estático y dependiente de la
muestra**, útil como primer cribado de redundancia, pero insuficiente para
decidir la relevancia de una covariable.

## Por qué las correlaciones climáticas cambian con la escala y la estación

El hecho de que el coeficiente de Pearson entre el NO₂ y una variable climática
varíe —en magnitud e incluso en signo— según se calcule a escala horaria, diaria
o mensual, o según la época del año, **no es un error ni una inconsistencia de
los datos**: es la consecuencia esperable de que una correlación sea un resumen
condicionado al diseño de muestreo. Intervienen cuatro fenómenos.

### 1. La escala determina qué mecanismo físico se mide (cambio de soporte)

Cada nivel de agregación temporal expone una fuente de variabilidad distinta:

- A escala **horaria** domina el **ciclo diurno**: la radiación solar actúa casi
  como un reloj (nula de noche, alta de día) y el NO₂ se acumula en las horas sin
  mezcla atmosférica; la correlación refleja, sobre todo, la alternancia
  día/noche.
- A escala **diaria** el ciclo diurno se promedia y aflora la **variabilidad
  sinóptica** (paso de frentes, episodios de viento).
- A escala **mensual** domina el **ciclo estacional** (acumulación invernal frente
  a dispersión estival).

Como cada escala mide un mecanismo diferente, el mismo par de variables produce
correlaciones distintas. Además, agregar promedia el ruido y realza la señal
lenta, por lo que los coeficientes mensuales tienden a ser considerablemente
mayores que los horarios. Este fenómeno es un caso del **problema del cambio de
soporte** (o de la unidad de análisis modificable): correlaciones calculadas a
distintos niveles de agregación no son comparables porque no miden lo mismo.

### 2. La relación es no estacionaria (depende de la época del año)

El vínculo entre el NO₂ y la meteorología no es constante a lo largo del año. En
**invierno**, con atmósfera estable, emisiones de calefacción y escasa mezcla
vertical, la temperatura y la radiación aparecen fuertemente anticorreladas con
el NO₂. En **verano**, la mayor mezcla y la actividad fotoquímica modifican esa
relación. Un coeficiente de correlación resume una relación que además es curva
sobre el rango de valores observado; al restringir ese rango a una única estación
del año se observa un tramo distinto de la curva y el ajuste lineal resultante
cambia. No es que la física varíe, sino que **cambia qué mecanismo domina la
varianza observada** en cada régimen.

### 3. El coeficiente de Pearson solo capta la parte lineal de relaciones curvas

Las relaciones del NO₂ con el viento (aproximadamente inversa), la radiación
(decreciente y curva) o la temperatura (en forma de U) no son lineales. Al
muestrear esas curvas sobre rangos distintos (estaciones) o agregaciones
distintas (escalas), la pendiente lineal que devuelve Pearson varía
apreciablemente. De ahí su inestabilidad.

### 4. El ciclo compartido infla y desestabiliza la correlación

Buena parte de la correlación climática no procede de un vínculo directo, sino de
que ambas series **siguen el mismo ciclo** (diurno o estacional). El peso de esa
componente compartida difiere según la escala y la estación del año. Su efecto es
tan grande que, al eliminarlo, las correlaciones cambian drásticamente: en este
trabajo, la correlación mensual del NO₂ con la velocidad del viento pasó de
−0.76 (sobre datos crudos) a ≈ +0.32 (sobre anomalías, tras retirar el ciclo
estacional), con un intervalo de confianza que roza el cero.

## Consecuencia metodológica y solución adoptada

El VIF y la correlación marginal fallan por la misma razón de fondo: proporcionan
**un único número sobre un conjunto de datos concreto**, ocultando que la
estructura de dependencia es **no estacionaria** (varía con la estación del año),
**dependiente de la escala** (varía con la agregación) y **no lineal**. Por ello
sus valores "bailan" y no deben tomarse como propiedades absolutas de las
variables.

Para obtener un diagnóstico robusto del tipo de relación se adoptó un
procedimiento que neutraliza cada uno de estos sesgos:

- **Cálculo sobre anomalías** (retirando el ciclo diurno/estacional de ambas
  series) para eliminar la correlación espuria inducida por el reloj común.
- **Estratificación por estación de medida y por escala temporal**, respetando la
  no estacionariedad y el cambio de soporte en lugar de promediarlos.
- **Uso conjunto de Pearson (asociación lineal), Spearman (monótona) y la
  distancia de correlación (cualquier dependencia)**, que permite separar la
  parte lineal de la no lineal.
- **Intervalos de confianza por *block bootstrap***, que preservan la
  autocorrelación temporal y distinguen la señal del ruido.

Y para la **selección de covariables**, en lugar de apoyarse únicamente en el
VIF, se complementó con criterios de **importancia condicional** (residuos
parciales en dos etapas y selección por DIC), que evalúan la contribución real de
cada variable al NO₂ controlando el efecto del resto. Este enfoque es el que
permite, por ejemplo, identificar que la velocidad del viento es un efecto
robusto mientras que parte del efecto aparente del tráfico funcionaba como
*proxy* de la autocorrelación temporal.
