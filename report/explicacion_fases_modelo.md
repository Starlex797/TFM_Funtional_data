# Modelización espacio-temporal del NO₂ en Madrid: explicación de las fases

Este documento recorre, fase a fase, la lógica del análisis: para cada etapa se
explica la **teoría** (por qué se hace y qué garantiza) y se **interpretan los
resultados** obtenidos en el proyecto. El hilo es: *qué mido → qué lo explica y
con qué forma → cómo lo encajo en un modelo espacio-temporal → qué estima el
modelo → me fío de él → está limpio de estructura no modelada*.

---

## FASE 0 — Contexto y objetivo

**Teoría.** El NO₂ urbano se genera por **emisión** (sobre todo tráfico) y se
elimina por **dispersión y transformación atmosférica** (viento, mezcla vertical,
fotoquímica). Además presenta fuerte **dependencia espacial** (estaciones
próximas se parecen) y **temporal** (persistencia de hora a hora y de día a día).
Un modelo adecuado debe combinar covariables físicas con un término que capture
esa estructura espacio-temporal residual. Elegimos el enfoque **INLA-SPDE**:
inferencia bayesiana aproximada (rápida) con un campo espacial continuo de tipo
**Matérn** representado como solución de una EDP estocástica (SPDE), y un término
temporal **AR(1)**.

**Datos.** NO₂ de 24 estaciones de Madrid (2025), en tres resoluciones
(horaria, diaria, mensual), cruzado con covariables de tráfico (intensidad,
carga) y clima interpolado (temperatura, humedad, precipitación, presión,
radiación, viento).

---

## FASE 1 — La variable respuesta

**Teoría.** El modelo gaussiano de INLA asume **errores simétricos y
homocedásticos**. Una concentración de contaminante es positiva, asimétrica a la
derecha y de naturaleza **multiplicativa** (acumulación/dilución proporcionales),
por lo que rara vez cumple esos supuestos en su escala original. La
transformación logarítmica los restablece.

**Resultados e interpretación.** En las tres tipologías de estación el NO₂ crudo
está claramente sesgado a la derecha (masa en 0–50 µg/m³ y cola hasta ~130), y
`log(NO₂+1)` **simetriza** la distribución y estabiliza la dispersión. Por tanto
se modela `log(NO₂+1)`. Es coherente con que el NO₂ sea positivo y multiplicativo,
y se confirma a posteriori con el gráfico cuantil-cuantil de los residuos (deben
aproximarse a la normal) y con Box–Cox (λ óptimo ≈ 0). Matiz metodológico
importante: el diagrama de dispersión sirve para decidir la forma de los
*predictores*; la transformación de la *respuesta* se valida con los **residuos**,
no con la nube cruda.

---

## FASE 2 — Selección de variables

Es el corazón del trabajo: decidir **qué covariables entran, con qué forma
(lineal / no lineal) y cuáles se descartan**. Se hace en capas, de lo más simple
a lo más riguroso.

### 2.1 Multicolinealidad (VIF)
**Teoría.** Si dos covariables están muy correlacionadas entre sí, sus
coeficientes se vuelven inestables y no interpretables (se "reparten" el efecto).
El **Factor de Inflación de la Varianza (VIF)** mide esa redundancia; se eliminan
hacia atrás las variables con VIF por encima de un umbral (5).

**Resultado.** A nivel horario, tras el filtro VIF **ninguna covariable se
elimina** (todas por debajo de 5): las ocho candidatas son suficientemente
independientes entre sí para coexistir en el modelo.

### 2.2 Forma de la relación con la respuesta (exploratorio)
**Teoría.** Antes de fijar términos lineales conviene ver la forma de cada
relación y a qué **escala temporal** actúa cada variable, porque una nube
marginal difusa puede esconder una relación real enmascarada por la mezcla de
estaciones o por el ciclo diario/estacional.

**Resultados e interpretación.**
- Aislando **una estación** (para quitar la variación entre emplazamientos), las
  relaciones se leen mejor: **viento** y **radiación** decrecientes y curvas,
  **presión** creciente, **tráfico** creciente en forma de cuña, y
  **humedad/precipitación** difusas incluso dentro de la estación.
- Los **perfiles cíclicos** revelan una separación **macro/micro**:
  - *Macro (estacional):* temperatura, radiación, humedad y presión tienen ciclo
    anual limpio que gobierna la envolvente del NO₂ (alto en invierno, bajo en
    verano).
  - *Micro (horaria/episódica):* el **tráfico** es casi plano a lo largo del año
    (solo cae en agosto) pero domina el **ciclo horario** (doble pico mañana y
    noche); **viento** y **precipitación** son **episódicos** (actúan día a día).
- Los **días-evento** confirman el mecanismo físico: en un día **lluvioso** o
  **ventoso** el NO₂ se desploma (de ~130 a ~24–45 µg/m³) pese a ser laborable de
  invierno → la meteorología "apaga" el NO₂ por dispersión/lavado.

### 2.3 Cuantificación no sesgada de la dependencia
**Teoría.** Una correlación cruda entre NO₂ y una covariable está contaminada
por: (i) el **ciclo compartido** (ambos siguen el reloj → correlación espuria),
(ii) la **mezcla de estaciones** (falacia ecológica), (iii) la
**autocorrelación** (infla la significancia) y (iv) la **no-linealidad** (Pearson
la infravalora). La solución: medir sobre **anomalías** (quitando el ciclo),
**por estación**, con **Pearson (lineal) + Spearman (monótona) + distancia de
correlación (cualquier dependencia)** e **intervalos por *block bootstrap***.

**Resultados e interpretación (Plaza Elíptica).** La corrección cambia el
diagnóstico de forma decisiva:
- **El tráfico mensual pasó de −0.38 (crudo, negativo y absurdo) a +0.65 / +0.50
  (positivo y coherente)**: el negativo era confusión estacional. Además
  intensidad y carga dejan de contradecirse.
- **El viento mensual cayó de −0.76 a +0.32 con IC que roza el 0**: casi toda
  aquella correlación gigante era el ciclo estacional compartido.
- La **radiación horaria** se desploma a ≈−0.08: su "efecto" era el reloj diurno
  que hemos retirado.
- Tras desestacionalizar, **casi todas las relaciones quedan "lineales"** y más
  débiles: **buena parte de la no-linealidad de la tabla cruda era el ciclo**, no
  una curvatura real. Las señales que sobreviven robustas: **viento** (fuerte,
  negativo, ~lineal a nivel horario/diario) y **presión** (positiva). Humedad y
  precipitación quedan débiles.

### 2.4 Diagnóstico en dos etapas (residuo parcial)
**Teoría.** Con N enorme (horario ≈ 17.500 obs/mes) **todo** sale significativo
por intervalo de confianza; el criterio operativo es el **tamaño del efecto**. Se
ajusta un modelo con las candidatas, se separan **fuertes** (|β estandarizado| ≥
0.10) y **débiles**, se resta la contribución de las fuertes para formar el
residuo `y*` y se comprueba si las **descartadas** aún explican `y*`.

**Resultados e interpretación (enero, horario).**
- **Fuertes:** viento (β = −0.45), carga (0.36), temperatura (−0.19), intensidad
  (0.13). Al restarlas, la desviación típica baja de `SD(y)=0.83` a
  `SD(y*)=0.60` → explican ~**48 %** de la varianza.
- **Débiles:** presión, humedad, precipitación, radiación. En la segunda etapa
  siguen siendo "significativas" (por el N), pero **reducen la SD solo un
  0.68 %** → aporte práctico **despreciable**.
- Con el modelo combinado se añade el **retardo de tráfico (h-1)**: la
  **intensidad(h-1)** resulta significativa (β = −0.022) pero también con un
  efecto minúsculo.
- **Conclusión de selección:** entran como efectos fijos las **fuertes** (tráfico
  + viento + temperatura); presión se mantiene por su interés físico/estadístico
  moderado; humedad y precipitación se descartan como *drivers*; lo que queda sin
  explicar es **estructura espacial + ruido**, no covariables olvidadas —
  justamente lo que capturará el campo SPDE.

---

## FASE 3 — Construcción del modelo INLA-SPDE

**Teoría.**
- **Malla (*mesh*).** El campo espacial continuo se aproxima sobre una malla
  triangular. `max.edge` controla la resolución (fina dentro del dominio, gruesa
  en el borde) y `cutoff` evita triángulos degenerados junto a estaciones
  próximas.
- **SPDE Matérn.** Un campo gaussiano Matérn es solución de una SPDE; INLA
  resuelve esa versión para ganar velocidad. Con `alpha = 2` la suavidad es la
  habitual (ν = 1). Sus dos hiperparámetros interpretables son el **rango**
  (distancia a la que la correlación cae a ~0.13) y la **varianza marginal**.
- **Proyección A + `inla.stack`.** La matriz `A` conecta cada observación con los
  vértices de la malla; el *stack* organiza respuesta, efectos fijos y campo
  latente. Para el modelo espacio-temporal, `A` y el índice llevan un **grupo
  temporal** y el campo se replica en el tiempo con estructura **AR(1)**:
  `Corr(w_t, w_{t-1}) = ρ`.
- **Dos fórmulas comparables.** Se mantiene idéntica la parte fija y solo cambia
  el término latente: `w(s)` (solo espacial) frente a `w(s,t)` con AR(1)
  (espacio-temporal). Así la comparación aísla el efecto de añadir la dinámica
  temporal.

**En el proyecto.** Se construyó un dataset de "día tipo" (24 estaciones × 7 días
de la semana) para medir el coste de cómputo y comparar ambos modelos antes de
lanzar el año completo.

---

## FASE 4 — Resultados del INLA

**Teoría de la lectura.** Un ajuste INLA se interpreta en tres bloques:
1. **Efectos fijos** (`summary.fixed`): media posterior e IC 95 % de cada β, **en
   escala log** (un β positivo → la covariable aumenta el NO₂; su magnitud es el
   cambio en log-NO₂ por unidad de covariable estandarizada).
2. **Hiperparámetros**: rango espacial (km), varianza del campo, **ρ del AR(1)**
   (persistencia temporal) y precisión del ruido gaussiano (**nugget**).
3. **Descomposición de la varianza**: qué parte explica cada componente
   (covariables / campo espacial / ruido).

**Resultados e interpretación.**
- **Signos de los efectos fijos** (coherentes con el análisis de correlaciones
  corregidas): **tráfico +** (más tráfico, más NO₂), **viento −** y
  **temperatura −** (dispersión y mezcla), **presión +** (estancamiento). Son los
  esperados físicamente.
- **Estructura latente.** El campo Matérn recoge la heterogeneidad espacial que
  las covariables no explican (nivel base distinto por emplazamiento: una
  estación de tráfico "respira" su calle, una de fondo refleja la contaminación
  urbana general). El **AR(1)** captura la persistencia temporal (una hora/día
  contaminado tiende a seguir contaminado), coherente con el ciclo diario y con
  la memoria que veíamos en las series.
- **Reparto de varianza** (diagnóstico de residuos, enero): las covariables
  fuertes explican ~48 % de la variabilidad de `log(NO₂)`; el resto se reparte
  entre **campo espacial** (autocorrelación real, confirmada por el test de
  Moran) y **ruido**. Es decir, el modelo espacio-temporal es necesario porque
  queda estructura espacial y temporal que las covariables no capturan.

---

## FASE 5 — Comparación y validación

**Teoría.** Hay que distinguir dos preguntas:
- **Ajuste dentro de muestra**: **DIC** y **WAIC** penalizan la complejidad;
  menor es mejor. Miden lo bien que el modelo describe los datos con los que se
  entrenó.
- **Capacidad predictiva fuera de muestra**: **validación cruzada espacial por
  bloques** (se ocultan estaciones **completas** y se predicen desde las vecinas)
  con **RMSE, MAE, cobertura del IC 95 % y CRPS**. Mide la interpolación a
  lugares **no observados**.

**Resultados e interpretación (día tipo).**

| Modelo | DIC | WAIC | p.eff | RMSE (in-sample) |
|---|---|---|---|---|
| Solo espacial | −401.7 | −399.0 | 27.9 | 0.062 |
| Espacio-temporal AR(1) | −750.7 | −749.6 | 77.3 | 0.015 |

- **Dentro de muestra el espacio-temporal gana con claridad** (DIC/WAIC mucho
  menores, RMSE 4× menor): el AR(1) añade flexibilidad que ajusta la dinámica
  temporal de las estaciones observadas. Su coste es más **parámetros efectivos**
  (p.eff 77 vs 28).
- **Matiz crucial (validación espacial).** Si en cambio se evalúa la
  **interpolación a estaciones no observadas** (CV por bloques), el
  espacio-temporal puede **empeorar** el RMSE/MAE respecto al solo espacial: el
  AR(1) no ayuda a un punto **sin ningún dato** (no hay serie temporal local que
  propagar), pero sí añade parámetros que no generalizan → paga complejidad sin
  beneficio. Conclusión honesta: **si el objetivo es rellenar/prever en
  estaciones observadas, el espacio-temporal es superior; si el objetivo es
  interpolar a ubicaciones nuevas, el AR(1) no aporta e incluso resta**.

---

## FASE 6 — Diagnóstico del modelo final

**Teoría.** Un buen modelo deja **residuos sin estructura**: (i) **QQ-plot** ≈
recta (normalidad), (ii) **ACF** de los residuos intra-estación ≈ 0 (sin memoria
temporal sobrante → el AR(1) ya la capturó), (iii) **variograma / índice de
Moran** sin autocorrelación espacial residual (→ el campo SPDE ya la capturó).

**Resultados e interpretación.** El diagnóstico en dos etapas mostró que, tras
retirar las covariables fuertes, el residuo **no esconde covariables
descartadas** (su aporte adicional es < 1 %) y sí presenta **autocorrelación
espacial significativa** (Moran) y **memoria temporal** (ACF de lag corto) —
exactamente lo que **campo SPDE + AR(1)** están diseñados para modelar. Es la
justificación empírica de la arquitectura del modelo.

---

## FASE 7 — Conclusiones

- El `log(NO₂)` en Madrid queda descrito por **cuatro componentes
  complementarias**: **tráfico local** (emisión, escala horaria), **meteorología**
  (dispersión: viento −, temperatura −, presión +; escala horaria y estacional),
  **campo espacial** (nivel base heterogéneo por emplazamiento) y **persistencia
  temporal AR(1)**.
- El fenómeno es **multi-escala**: a nivel horario mandan tráfico y dispersión; a
  nivel estacional manda la meteorología; el tráfico apenas tiene ciclo anual.
- Las covariables **higrométricas** (humedad, precipitación) son **débiles** como
  *drivers* continuos y se descartan; la precipitación, además, está inflada en
  cero.
- **Limitaciones**: el residuo de dos etapas trata los β de la primera etapa como
  conocidos (regresor generado); la validación cruzada espacial revela los
  límites del AR(1) para interpolar a puntos nuevos; el clima es interpolado.
- **Líneas futuras**: términos no lineales (splines) para viento/radiación,
  interacción tipología × tráfico, y validación **temporal** (ocultar horas en
  estaciones observadas) para cuantificar dónde el espacio-temporal sí gana.

---

### Trazabilidad (scripts y salidas)
- Transformación de la respuesta y forma de covariables: `15_`, `21_`, `23_`,
  `24_scatter/perfiles`, `report/seccion_transformaciones.md`.
- Correlaciones corregidas (anomalías + IC): `25_correlaciones_corregidas_anomalias.R`.
- Diagnóstico en dos etapas: `18_`, `20_residuos`.
- Malla/SPDE/stack/modelos: `Spde.R`, `Paso_4_inla_stack.R`, `Paso_5_modelo_inla.R`,
  `Paso_5b_modelo_espacial.R`.
- Comparación y validación: `comparacion.R`,
  `comparacion_modelos_dia_tipo.csv`.
