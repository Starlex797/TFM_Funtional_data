# Modelo de los residuos horarios: metodología paso a paso

Este apartado describe cómo se construye y evalúa el **modelo de los residuos a
nivel horario** sobre una muestra de 5 días, y por qué se toma cada decisión. El
objetivo no es solo ajustar covariables, sino **determinar qué estructura
adicional necesita el residuo** una vez retirado el efecto de las covariables
principales: ¿basta con las covariables?, ¿hace falta un término **temporal**
(y de qué orden, AR(1) o AR(2))?, ¿y un término **espacial** (campo SPDE)?

## Idea general y justificación

Tras ajustar las covariables dominantes, lo que queda —el **residuo**— no suele
ser ruido blanco: conserva la **autocorrelación temporal** (una hora contaminada
va seguida de otra) y la **autocorrelación espacial** (estaciones próximas se
parecen). Modelar explícitamente esa estructura es lo que distingue un modelo
geoestadístico espacio-temporal de una simple regresión. El experimento compara,
sobre el residuo, cinco especificaciones que **añaden estructura de forma
incremental**, para decidir con criterios objetivos cuál es necesaria.

## Paso 1 — Selección de la muestra (5 días, invierno)

Se toman los primeros **5 días de enero** (2025-01-01 a 2025-01-05), las 24
estaciones, a resolución horaria (2.880 observaciones, 120 instantes de tiempo).

*Por qué:* una ventana corta permite **comparar varias especificaciones de
modelo** (incluidas las espacio-temporales, costosas) en un tiempo de cómputo
razonable, manteniendo el ciclo diurno completo (5 × 24 h). Se elige **invierno**
porque es cuando el NO₂ y su dinámica temporal son más marcados, de modo que la
estructura a detectar es más clara.

## Paso 2 — Estandarización de covariables

Cada covariable candidata (tráfico: intensidad y carga; clima: temperatura,
humedad, precipitación, presión, radiación, viento) se **estandariza** (media 0,
desviación 1) dentro de la muestra.

*Por qué:* deja los coeficientes en una **escala común** y permite usar un
**umbral de tamaño de efecto** (|β estandarizado| ≥ 0.10) para separar
covariables "fuertes" de "débiles". Esto es necesario porque, con N grande, la
significancia estadística no discrimina (casi todo sale significativo); el
criterio operativo es la **magnitud** del efecto.

## Paso 3 — Coordenadas en km y retardo temporal

Las coordenadas se proyectan a **UTM en kilómetros**, y se construye el **retardo
h-1 del tráfico** (`intensidad` de la hora anterior, por estación).

*Por qué:* el campo espacial SPDE (Matérn) se define sobre **distancias reales**,
por lo que se necesitan coordenadas métricas, no grados. El retardo h-1 se añade
como candidata porque el efecto del tráfico sobre el NO₂ **no es instantáneo**:
la contaminación de una hora arrastra la emisión de la anterior.

## Paso 4 — Índices de tiempo y de estación

Se crea un índice temporal `t` (1…120, el instante horario) y un índice de
estación.

*Por qué:* son la base de los términos estructurales. El **AR** se define sobre
`t` (replicado por estación → cada estación tiene su propia serie temporal), y el
componente **espacio-temporal** agrupa el campo espacial por instante `t`.

## Paso 5 — Partición train / test (temporal)

**Entrenamiento = días 1-4 (96 h)**, **test = día 5 (24 h)**. En el ajuste, la
respuesta de las horas de test se pone a `NA`.

*Por qué:* evalúa la capacidad **predictiva fuera de muestra** (predecir un día
futuro), no solo el ajuste. Poner `NA` en el test permite que INLA **prediga**
esas horas usando el modelo entrenado con el resto, de forma limpia y sin fuga de
información.

## Paso 6 — Construcción del residuo y\* (dos etapas)

Se ajusta una regresión de `log(NO₂)` sobre todas las covariables, se identifican
las **fuertes** (|β estandarizado| ≥ 0.10) y se forma la **pseudo-respuesta**:

**y\*_sh = log(NO₂)_sh − Σ_j β̂_j · X_j,sh   (j = covariables fuertes)**

*Por qué:* al restar la contribución de las covariables dominantes, `y*` es
**"lo que queda por explicar"**. Modelar `y*` aísla la pregunta de interés: de esa
variabilidad restante, ¿cuánta es estructura temporal, cuánta espacial, y cuánta
explican todavía las covariables débiles? Es un **residuo de dos etapas**; nota
metodológica: los β̂ de la primera etapa se tratan como conocidos (problema del
"regresor generado"), por lo que la inferencia de la segunda etapa es orientativa.

## Paso 7 — Malla y campo SPDE (Matérn)

Se construye una **malla triangular** sobre Madrid (resolución media, `max.edge`
≈ 4-8 km) y un **campo Matérn** vía SPDE (`alpha = 2`).

*Por qué:* el SPDE aproxima un campo gaussiano continuo como solución de una
ecuación diferencial estocástica, lo que permite a INLA estimarlo con rapidez.
Representa la **dependencia espacial** del residuo (parches de contaminación de
unos km). La malla **media** es el compromiso entre precisión y coste (en el
análisis previo daba el mismo RMSE que la fina en una fracción del tiempo).

## Paso 8 — Los cinco modelos comparados

Sobre `y*`, todos comparten los mismos **efectos fijos** (covariables débiles +
retardo del tráfico) y solo cambian en la **estructura añadida**:

| Modelo | Estructura añadida | Qué comprueba |
|---|---|---|
| **A** | ninguna (solo covariables) | ¿bastan las covariables débiles? |
| **B** | + **AR(1)** temporal (por estación) | ¿hay persistencia temporal de orden 1? |
| **C** | + **AR(2)** temporal (por estación) | ¿mejora un orden temporal 2? |
| **D** | + **campo espacial SPDE** | ¿hay estructura espacial? |
| **E** | **SPDE ⊗ AR(1)** (espacio-temporal) | ¿conviene combinar ambas? |

*Por qué esta batería:* aislando cada término se puede atribuir la mejora a su
causa. Comparar **B vs C** responde directamente a la pregunta de si un **AR(2)**
aporta sobre el **AR(1)** (relevante si en la ACF del AR(1) quedaban picos en los
retardos 2-3). Comparar **A/B/C vs D** contrasta lo temporal frente a lo
espacial, y **E** evalúa la combinación.

## Paso 9 — Métricas de comparación

Para cada modelo se calculan:
- **DIC y WAIC**: ajuste penalizando la complejidad (menor = mejor). Miden lo bien
  que el modelo describe los datos de entrenamiento sin premiar el sobreajuste.
- **RMSE y MAE** sobre el test: error de **predicción puntual** del día 5.
- **Cobertura 95%**: porcentaje de observaciones de test dentro del intervalo
  predictivo del 95% → mide la **calibración de la incertidumbre**.

*Por qué las tres a la vez:* son complementarias y evitan conclusiones erróneas.
El RMSE solo mira el punto y **puede favorecer a un modelo mal calibrado**; el
DIC/WAIC mira el ajuste global; la cobertura comprueba si los intervalos son
honestos. Un buen modelo debe equilibrar las tres (RMSE bajo **y** cobertura
cercana al 95%). El intervalo predictivo se construye con la desviación
**predictiva completa**, sumando la varianza del ajuste y la del ruido de
observación (1/precisión gaussiana).

## Paso 10 — Diagnóstico de residuos del mejor modelo

Del modelo con menor WAIC se examinan:
- **ACF** de los residuos (media por instante): comprueba si **queda
  autocorrelación temporal** sin modelar (retardos dentro de las bandas = ruido
  blanco → estructura temporal bien capturada).
- **QQ-plot**: comprueba la **normalidad** de los residuos (supuesto del modelo
  gaussiano); desviaciones en las colas indican que los eventos extremos se
  ajustan peor.

*Por qué:* cierran el círculo. La selección por WAIC dice **qué modelo ajusta
mejor**; los diagnósticos confirman que ese modelo **deja un residuo limpio**
(sin estructura temporal ni desvíos graves de la normalidad), que es la condición
para que su inferencia sea válida.

## Resumen del razonamiento

> Se aísla el residuo del NO₂ horario (`y*`), y se le añade estructura de forma
> incremental —temporal (AR1/AR2) y espacial (SPDE)— comparando con criterios que
> miden a la vez **ajuste** (DIC/WAIC), **precisión** (RMSE/MAE) y **calibración**
> (cobertura), y validando con **diagnósticos de residuos** (ACF, QQ). Así se
> decide, con evidencia y no por defecto, **qué complejidad merece el modelo** a
> escala horaria: si el residuo exige un término temporal, de qué orden, y si
> además necesita el campo espacial.
