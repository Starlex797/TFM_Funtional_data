# Modelo INLA-SPDE de los residuos del modelo horario

## Objetivo y justificación

Una vez ajustado el modelo horario espacio-temporal, la pregunta que se plantea
es si **queda en sus residuos alguna covariable relevante que el modelo no haya
incorporado**. Si la hubiera, sería una señal de que el modelo horario está
**mal especificado por omisión** y de que esa covariable debería añadirse. Si no
la hay, se confirma que el modelo ya recoge toda la señal que las covariables
disponibles pueden aportar y que lo que queda es ruido sin estructura explicable.

Para responderlo de forma rigurosa no basta con mirar correlaciones simples: hay
que **ajustar un modelo geoestadístico sobre los residuos** que separe lo que una
covariable explica de lo que es mera proximidad espacial. Por eso el modelo de
los residuos es a su vez un **INLA-SPDE** (campo gaussiano Matérn), y no una
regresión ordinaria.

## Paso 1 — Punto de partida: el modelo horario

El análisis parte del **modelo horario de invierno** ya seleccionado:

- **Muestra:** primeros 5 días de enero (2025-01-01 a 2025-01-05), las 24
  estaciones, a resolución horaria; log(NO₂) como respuesta. Entrenamiento =
  días 1-4 (96 h); test = día 5 (24 h).
- **Covariables** (seleccionadas mediante filtro VIF, significancia bayesiana y
  *stepwise* por DIC): **intensidad de tráfico, temperatura y velocidad del
  viento**. La precipitación se descartó por no ser significativa (su intervalo
  de credibilidad al 95 % contenía el 0), de modo que **no forma parte** del
  modelo horario.
- **Estructura:** campo espacial **Matérn (SPDE)** combinado con un término
  **AR(1) temporal** por instante horario. Para este diagnóstico el campo se
  ajusta sobre la malla triangular **gruesa (~8 km)**, la misma que se usa en el
  modelo de los residuos, de modo que ambos comparten resolución espacial. (La
  resolución de la malla no altera la conclusión: el ajuste queda dominado por el
  término AR(1), que explica casi toda la variabilidad con independencia de la
  malla.)

## Paso 2 — Construcción de los residuos

Se obtiene el residuo de cada observación de entrenamiento como

**r_sh = log(NO₂)_sh − ŷ_sh**,

donde ŷ_sh es la predicción del modelo horario (efecto de las covariables usadas
más el campo espacio-temporal). El día de test se **enmascara** (respuesta a
`NA`) en el ajuste, de forma que los residuos analizados son estrictamente los
del **conjunto de entrenamiento**, coherentes con el diagnóstico de residuos del
propio modelo horario.

El residuo `r` representa **la parte de la contaminación horaria que el modelo no
explica**: ni por sus covariables ni por su estructura espacio-temporal.

## Paso 3 — Qué covariables se introducen (y cuáles no)

El principio metodológico central es que **el modelo de los residuos no puede
contener las covariables que el modelo horario ya usa**. Reintroducir
`intensidad`, `temperatura` o `viento` sería redundante: su efecto ya está
retirado del residuo y su coeficiente sería, por construcción, prácticamente
nulo. Por tanto, en el modelo de residuos se introducen **únicamente las
covariables NO usadas**:

- `carga` de tráfico,
- `humedad` relativa,
- `presión` barométrica,
- `radiación` solar,
- `precipitación` (que quedó fuera del modelo horario),
- el **retardo h-1** de la intensidad de tráfico (información intra-horaria que un
  modelo no puede capturar sin un término explícito de rezago).

Antes de introducirlas se aplica un **filtro VIF** (umbral 5) entre ellas, para
garantizar que **no exista colinealidad** entre las covariables candidatas: si dos
estuvieran muy correlacionadas, una podría inflar o enmascarar el efecto de la
otra y falsear la conclusión.

## Paso 4 — El modelo INLA-SPDE de los residuos

Sobre el residuo `r` se ajusta el modelo

**r_s = β₀ + Σ_k β_k · x_{k,s} + ω(s) + ε_s**,

donde:

- **β_k** son los efectos de las covariables no usadas,
- **ω(s)** es un **campo gaussiano Matérn** (SPDE) sobre la malla **gruesa
  (~8 km)**, que recoge la estructura espacial que aún compartan los residuos
  entre estaciones próximas,
- **ε_s** es el ruido de observación.

El campo espacial ω es esencial: **separa** la variabilidad que una covariable
explica de la que es simple vecindad espacial, evitando atribuir a una covariable
un efecto que en realidad es autocorrelación espacial del residuo.

## Paso 5 — Criterios de decisión

- **Significancia:** una covariable se considera importante si su **intervalo de
  credibilidad bayesiano al 95 %** (cuantiles 2.5 % y 97.5 % de la posterior de
  β_k) **no contiene el 0**. Es el mismo criterio de significancia usado en la
  selección del modelo horario. No se emplean *p*-valores porque el ajuste es
  bayesiano (INLA).
- **Tamaño del efecto:** se acompaña del **coeficiente estandarizado**
  β_std = β_k / SD(r), que expresa el efecto en desviaciones típicas del residuo
  y permite comparar covariables en una escala común.
- **Estructura espacial remanente:** el **rango** y la **varianza** del campo ω
  indican si en los residuos queda dependencia espacial no captada.
- **Varianza explicada:** el porcentaje del residuo que explica el modelo
  completo (campo + covariables) y el que explican las covariables por sí solas.
- **Autocorrelación remanente:** **ACF** y test de **Ljung-Box** sobre el residuo
  final, para comprobar que no queda estructura temporal sin modelar.

## Paso 6 — Interpretación

- Si **alguna** covariable no usada resulta significativa (IC 95 % sin el 0) y con
  tamaño de efecto apreciable → es **información que el modelo horario omitió** y
  que convendría incorporar.
- Si **ninguna** lo es → el modelo horario **ya captura toda la señal explicable**
  con las covariables disponibles; el residuo es ruido sin contenido atribuible a
  covariables, y no procede ampliar el modelo.

En ambos casos, el rango/varianza del campo ω y el diagnóstico de autocorrelación
informan sobre si el residuo conserva estructura (espacial o temporal) que ninguna
covariable explica.

## Resultados

Efectos de las covariables **no usadas** en el modelo SPDE de los residuos
(ordenadas por tamaño de efecto):

| Covariable no usada | β | β_std | IC 95 % | ¿Importante? |
|---|---|---|---|---|
| presión | 0.0021 | 0.043 | [−0.004, 0.008] | no |
| humedad | 0.0016 | 0.032 | [−0.001, 0.004] | no |
| precipitación | −0.0015 | −0.030 | [−0.011, 0.009] | no |
| carga de tráfico | 0.0012 | 0.025 | [−0.002, 0.004] | no |
| radiación | 0.0005 | 0.011 | [−0.002, 0.003] | no |
| retardo h-1 tráfico | −0.0002 | −0.005 | [−0.003, 0.002] | no |

- **R² del modelo horario (train):** 99.2 % (SD del residuo = 0.049 frente a
  SD de log NO₂ = 0.544).
- **Campo espacial de los residuos:** rango ≈ 101 km, σ² ≈ 1·10⁻⁵. Un rango
  mayor que la propia extensión del área de estudio con varianza prácticamente
  nula indica que **no queda estructura espacial** en los residuos.
- **% del residuo explicado:** modelo completo (campo + covariables) 0.07 %;
  covariables solas 0.07 %. Es decir, las covariables no usadas **no explican
  prácticamente nada** del residuo.
- **Ljung-Box (autocorrelación remanente):** p ≈ 2·10⁻¹¹.

> **Conclusión.** Ninguna de las covariables no incorporadas al modelo horario
> (carga de tráfico, humedad, presión, radiación, precipitación ni el retardo
> h-1 del tráfico) resulta significativa dentro de sus residuos: todos los
> intervalos de credibilidad al 95 % contienen el 0 y los tamaños de efecto son
> despreciables (|β_std| < 0.05). El modelo horario espacio-temporal **ya recoge
> toda la señal explicable** con las covariables disponibles; el residuo carece
> de estructura espacial (campo Matérn de varianza nula) y las covariables
> omitidas no aportan información. **No procede, por tanto, ampliar el modelo
> horario con ninguna de estas covariables.**
>
> *Matiz metodológico.* El modelo horario explica el 99.2 % de la varianza en
> entrenamiento gracias a la gran flexibilidad del campo espacio-temporal
> (SPDE ⊗ AR1), de modo que el residuo analizado es muy pequeño. Esto implica que
> el contraste tiene, por construcción, **poca potencia** para detectar efectos
> débiles: el campo espacio-temporal absorbe buena parte de la variabilidad que
> de otro modo podría atribuirse a covariables. La conclusión "ninguna covariable
> queda en los residuos" debe leerse como *"el modelo espacio-temporal ya captura
> esa estructura"*, no como que dichas covariables carezcan de relación física con
> el NO₂. El leve residuo de autocorrelación temporal (Ljung-Box) confirma que
> queda una pequeña estructura temporal, pero de magnitud despreciable
> (0.07 % de varianza).
