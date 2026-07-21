# Metodologia del modelo INLA-SPDE

Este apartado describe el procedimiento seguido para ajustar los modelos
INLA-SPDE utilizados en el analisis espacio-temporal de las concentraciones de
NO2 en Madrid. Aunque el archivo principal se denomina `simulacion.R`, el
procedimiento no corresponde a una simulacion Monte Carlo, sino a una
validacion temporal sobre datos reales mediante holdout. Es decir, el modelo se
ajusta con un tramo temporal de entrenamiento y se evalua sobre un tramo final
no utilizado durante el ajuste.

El objetivo del procedimiento es doble. Por un lado, se busca identificar las
covariables meteorologicas y de trafico asociadas a la variabilidad del NO2 en
distintas escalas temporales. Por otro lado, se pretende construir un modelo
capaz de capturar la estructura espacial residual mediante un campo gaussiano
de Matern aproximado por SPDE, y, cuando es necesario, incorporar dependencia
temporal mediante una estructura autorregresiva AR(1).

## 1. Variable respuesta y estructura general del modelo

La variable de estudio es la concentracion de dioxido de nitrogeno, NO2. Dado
que la distribucion de NO2 presenta asimetria positiva y valores extremos, se
trabaja con la transformacion logaritmica:

$$
y_i = \log(NO_{2,i}+1),
$$

donde \(y_i\) representa la observacion transformada en la localizacion
\(s_i\) y en el instante temporal \(t_i\).

El modelo utilizado es un modelo bayesiano jerarquico gaussiano. La primera
capa del modelo, correspondiente a la verosimilitud, se define como:

$$
y_i \mid \eta_i, \sigma_\varepsilon^2
\sim
\mathcal{N}(\eta_i,\sigma_\varepsilon^2),
\qquad i=1,\dots,n.
$$

El predictor lineal se escribe como:

$$
\eta_i =
\beta_0
+
\sum_{k=1}^{K}\beta_k X_k(s_i,t_i)
+
u(s_i,t_i),
$$

donde:

- \(\beta_0\) es el intercepto.
- \(X_k(s_i,t_i)\) son las covariables meteorologicas o de trafico.
- \(\beta_k\) mide el efecto fijo de la covariable \(k\).
- \(u(s_i,t_i)\) es el efecto espacial o espacio-temporal no observado.
- \(\varepsilon_i \sim \mathcal{N}(0,\sigma_\varepsilon^2)\) es el error de observacion o nugget.

En forma compacta, el modelo puede expresarse como:

$$
\mathbf{y}
=
\mathbf{X}\boldsymbol{\beta}
+
\mathbf{A}\mathbf{w}
+
\boldsymbol{\varepsilon}.
$$

Aqui, \(\mathbf{X}\) es la matriz de covariables, \(\boldsymbol{\beta}\) el vector
de efectos fijos, \(\mathbf{w}\) el vector de valores latentes del campo espacial
en los nodos de la malla, y \(\mathbf{A}\) la matriz de proyeccion que lleva el
campo desde la malla hasta las estaciones de medicion.

## 2. Preparacion de datos y particion train-test

El primer paso del procedimiento consiste en preparar los datos de acuerdo con
la escala temporal considerada. El script trabaja con varias escalas:

- Escala mensual.
- Escala diaria.
- Escala horaria.
- Perfil horario laborable/finde.

En todos los casos se construye una variable respuesta comun, denominada en el
script `LOG_NO2`, que corresponde a \(\log(NO2+1)\). Tambien se preparan las
coordenadas proyectadas de las estaciones, \(X_{km}\) e \(Y_{km}\), utilizadas
para construir la malla espacial.

La validacion del modelo se realiza mediante holdout temporal. Esto significa
que se deja fuera el tramo final de la serie como conjunto de test. Formalmente,
si \(T_{train}\) es el conjunto de tiempos de entrenamiento y \(T_{test}\) el
conjunto de tiempos de evaluacion, entonces:

$$
T = T_{train} \cup T_{test},
\qquad
T_{train} \cap T_{test} = \emptyset.
$$

Durante el ajuste, las observaciones del conjunto test se introducen como
valores faltantes:

$$
y_i^{obs} =
\begin{cases}
y_i, & t_i \in T_{train},\\
NA, & t_i \in T_{test}.
\end{cases}
$$

INLA estima el modelo con las observaciones de entrenamiento y obtiene la
distribucion predictiva posterior para las observaciones del test.

Esta particion temporal es preferible a una particion aleatoria porque las
observaciones de contaminacion presentan autocorrelacion temporal. Una division
aleatoria podria mezclar informacion del futuro en el entrenamiento y producir
una evaluacion demasiado optimista.

## 3. Estandarizacion de covariables

Las covariables candidatas se estandarizan antes de entrar en el modelo:

$$
\widetilde{X}_{k,i}
=
\frac{X_{k,i}-\bar{X}_k}{s_k},
$$

donde \(\bar{X}_k\) y \(s_k\) son la media y desviacion tipica de la covariable
\(k\) en el conjunto de entrenamiento.

Esta transformacion tiene dos ventajas. En primer lugar, mejora la estabilidad
numerica del ajuste. En segundo lugar, permite comparar los coeficientes
\(\beta_k\), ya que cada uno representa el cambio esperado en \(\log(NO2+1)\)
ante un incremento de una desviacion tipica en la covariable correspondiente.

## 4. Seleccion inicial por residuos parciales

Antes del ajuste final se realiza una ordenacion de las covariables mediante un
procedimiento de residuos parciales. La idea es identificar que covariables
reducen mas la variabilidad de la respuesta cuando se estiman dentro de un
modelo INLA-SPDE.

Se parte del residuo inicial:

$$
r^{(0)} = y.
$$

En cada paso \(m\), para cada covariable candidata \(X_j\), se ajusta un modelo
univariante de la forma:

$$
r^{(m-1)}_i
=
\alpha
+
\beta_j X_{j,i}
+
u_j(s_i)
+
\varepsilon_i.
$$

Despues se calcula el residuo que quedaria al retirar la contribucion estimada
de esa covariable:

$$
r^{(m)}_j =
r^{(m-1)}
-
\widehat{\beta}_j X_j.
$$

La covariable seleccionada en el paso \(m\) es aquella que mas reduce la
desviacion tipica del residuo:

$$
j^\ast
=
\arg\min_j
\operatorname{sd}
\left(
r^{(m-1)}
-
\widehat{\beta}_j X_j
\right).
$$

Entonces se actualiza:

$$
r^{(m)}
=
r^{(m-1)}
-
\widehat{\beta}_{j^\ast}X_{j^\ast}.
$$

Este paso no se realiza con una regresion lineal ordinaria, sino con un modelo
INLA-SPDE univariante. La razon es que el coeficiente de cada covariable se
estima teniendo en cuenta la estructura espacial residual. De esta forma se
reduce el riesgo de atribuir a una covariable un efecto que en realidad
corresponde a un patron espacial no observado.

## 5. Control de multicolinealidad mediante VIF

Tras la ordenacion inicial, se aplica un filtro de multicolinealidad. Para cada
covariable \(X_k\), se calcula el factor de inflacion de la varianza:

$$
VIF_k = \frac{1}{1-R_k^2},
$$

donde \(R_k^2\) es el coeficiente de determinacion obtenido al regresar
\(X_k\) sobre el resto de covariables.

Un valor alto de \(VIF_k\) indica que \(X_k\) esta fuertemente explicada por
otras covariables, lo que puede inflar la varianza de \(\widehat{\beta}_k\) y
hacer inestable su interpretacion. En el procedimiento se utiliza el umbral:

$$
VIF_k \leq 5.
$$

Cuando alguna covariable supera este umbral, se elimina iterativamente la
covariable con mayor VIF hasta que todas las covariables retenidas cumplen el
criterio.

## 6. Seleccion bayesiana de covariables

Una vez eliminada la multicolinealidad excesiva, se ajusta un modelo INLA-SPDE
con las covariables restantes. Para cada efecto fijo \(\beta_k\), INLA
proporciona su distribucion posterior marginal:

$$
\pi(\beta_k \mid \mathbf{y}).
$$

Una covariable se considera significativa si su intervalo creible posterior al
95 % no contiene el cero:

$$
0 \notin
\left[
q_{0.025}(\beta_k),
q_{0.975}(\beta_k)
\right].
$$

Este criterio permite decidir si existe evidencia posterior suficiente de que
el efecto de una covariable es positivo o negativo.

Adicionalmente, se realiza una seleccion stepwise basada en DIC. Se aceptan
cambios en el conjunto de covariables si la mejora en DIC es al menos:

$$
\Delta DIC \geq 2.
$$

El objetivo de este paso es evitar modelos innecesariamente complejos cuando la
mejora estadistica es pequena.

## 7. Construccion de la malla espacial

El enfoque SPDE requiere discretizar el dominio espacial mediante una malla
triangular. En este trabajo, la malla se construye a partir de las coordenadas
de las estaciones de medicion en kilometros.

El script compara varias resoluciones de malla:

- Malla muy gruesa.
- Malla gruesa.
- Malla media.
- Malla fina.

La eleccion de la malla afecta al equilibrio entre precision y coste
computacional. Una malla mas fina aproxima mejor el campo espacial continuo,
pero incrementa el numero de nodos \(m\) y, por tanto, el numero de variables
latentes que debe estimar el modelo.

## 8. Campo espacial Matern y formulacion SPDE

El componente espacial no observado se modeliza como un campo gaussiano con
covarianza Matern:

$$
\operatorname{Cov}
\left(
u(s),u(s')
\right)
=
\sigma^2
\frac{2^{1-\nu}}{\Gamma(\nu)}
\left(
\kappa ||s-s'||
\right)^\nu
K_\nu
\left(
\kappa ||s-s'||
\right),
$$

donde:

- \(\sigma^2\) es la varianza marginal del campo espacial.
- \(\kappa\) controla la escala espacial.
- \(\nu\) es el parametro de suavidad.
- \(K_\nu\) es la funcion de Bessel modificada de segunda especie.
- \(||s-s'||\) es la distancia euclidea entre dos localizaciones.

El rango espacial se define como:

$$
r = \frac{\sqrt{8\nu}}{\kappa}.
$$

Este rango se interpreta como la distancia a partir de la cual la correlacion
espacial se vuelve pequena. En terminos practicos, indica hasta que distancia
dos estaciones tienden a compartir informacion espacial residual.

El metodo SPDE utiliza el resultado de Lindgren, Rue y Lindstrom, segun el cual
un campo de Matern puede representarse como solucion de la ecuacion diferencial
parcial estocastica:

$$
\left(
\kappa^2 - \Delta
\right)^{\alpha/2}
\left(
\tau u(s)
\right)
=
\mathcal{W}(s),
$$

donde \(\Delta\) es el operador laplaciano, \(\mathcal{W}(s)\) es ruido blanco
gaussiano y:

$$
\alpha = \nu + \frac{d}{2}.
$$

Como el dominio espacial es bidimensional, \(d=2\). En el script se utiliza
`alpha = 2`, por lo que:

$$
\nu = \alpha - \frac{d}{2}
=
2 - 1
=
1.
$$

Por tanto, el modelo ajusta un campo Matern con suavidad \(\nu=1\).

La aproximacion por elementos finitos representa el campo continuo como:

$$
u(s)
\approx
\sum_{j=1}^{m}
\psi_j(s)w_j,
$$

donde \(\psi_j(s)\) son funciones base asociadas a los nodos de la malla y
\(w_j\) son pesos latentes. En forma matricial:

$$
\mathbf{u}
=
\mathbf{A}\mathbf{w}.
$$

El vector de pesos latentes sigue una distribucion gaussiana:

$$
\mathbf{w}
\sim
\mathcal{N}
\left(
\mathbf{0},
\mathbf{Q}(\kappa,\tau)^{-1}
\right).
$$

La matriz \(\mathbf{Q}(\kappa,\tau)\) es una matriz de precision dispersa. Esta
propiedad es esencial, ya que permite realizar inferencia bayesiana de forma
eficiente con INLA.

En el codigo se utilizan priors de complejidad penalizada para el SPDE:

$$
P(r < 5) = 0.5,
\qquad
P(\sigma > 1) = 0.01.
$$

Estos priors regularizan el rango espacial y la varianza marginal del campo,
evitando soluciones excesivamente complejas o inestables.

## 9. Matriz de proyeccion e INLA stack

Las observaciones no se localizan en los nodos de la malla, sino en las
estaciones de medicion. Por ello se construye una matriz de proyeccion
\(\mathbf{A}\) que interpola el campo desde los nodos hasta las localizaciones
observadas.

Cada elemento de \(\mathbf{A}\) puede escribirse como:

$$
A_{ij} = \psi_j(s_i),
$$

de forma que:

$$
u(s_i)
=
\sum_{j=1}^{m}
A_{ij}w_j.
$$

En INLA, la respuesta, la matriz de proyeccion, los indices del campo latente y
las covariables se ensamblan mediante `inla.stack`. Conceptualmente, el stack
construye el modelo:

$$
\boldsymbol{\eta}
=
\mathbf{X}\boldsymbol{\beta}
+
\mathbf{A}\mathbf{w}.
$$

Esto permite que el modelo estime simultaneamente los efectos fijos de las
covariables y el campo espacial residual.

## 10. Modelo solo espacial

El primer modelo ajustado en cada escala es un modelo solo espacial:

$$
y(s_i,t_i)
=
\beta_0
+
\sum_{k=1}^{K}
\beta_k X_k(s_i,t_i)
+
u(s_i)
+
\varepsilon_i.
$$

En este modelo, el campo espacial \(u(s_i)\) no cambia con el tiempo. Por tanto,
captura diferencias espaciales persistentes entre zonas o estaciones que no son
explicadas por las covariables.

Este modelo es util como primera aproximacion porque permite responder a la
pregunta:

> Existe estructura espacial residual en el NO2 una vez controladas las
> covariables meteorologicas y de trafico?

Si el campo espacial tiene varianza relevante o mejora claramente la prediccion,
significa que existe informacion espacial no capturada por las covariables
incluidas.

## 11. Diagnostico de residuos del modelo espacial

Tras ajustar el modelo solo espacial, se analizan los residuos:

$$
\widehat{\varepsilon}_i
=
y_i - \widehat{y}_i.
$$

El objetivo es comprobar si el modelo ha capturado adecuadamente la estructura
de los datos. Para ello se utilizan dos diagnosticos:

1. La funcion de autocorrelacion de los residuos.
2. El grafico Q-Q de normalidad.

La autocorrelacion muestral en el retardo \(h\) se define como:

$$
\widehat{\rho}(h)
=
\frac{
\sum_t
\left(
r_t-\bar{r}
\right)
\left(
r_{t-h}-\bar{r}
\right)
}{
\sum_t
\left(
r_t-\bar{r}
\right)^2
}.
$$

Ademas, se aplica el test de Ljung-Box:

$$
Q
=
n(n+2)
\sum_{h=1}^{H}
\frac{\widehat{\rho}(h)^2}{n-h}.
$$

Si el test detecta autocorrelacion significativa, se considera que el modelo
solo espacial no es suficiente, porque queda estructura temporal en los
residuos. En ese caso se ajusta el modelo espacio-temporal con AR(1).

## 12. Modelo espacio-temporal con AR(1)

Cuando los residuos del modelo espacial presentan autocorrelacion temporal, se
introduce una estructura AR(1) en el campo latente. El modelo pasa a ser:

$$
y(s_i,t_i)
=
\beta_0
+
\sum_{k=1}^{K}
\beta_k X_k(s_i,t_i)
+
u(s_i,t_i)
+
\varepsilon_i.
$$

El campo espacio-temporal evoluciona segun:

$$
u(s,t)
=
\rho u(s,t-1)
+
\xi(s,t),
$$

donde \(\rho\) es el parametro de autocorrelacion temporal y \(\xi(s,t)\) es
una innovacion espacial con estructura Matern.

En terminos del vector latente de la malla:

$$
\mathbf{w}_t
=
\rho \mathbf{w}_{t-1}
+
\boldsymbol{\xi}_t.
$$

La matriz de precision conjunta puede escribirse como un producto de Kronecker:

$$
\mathbf{Q}_{ST}
=
\mathbf{Q}_{AR1}(\rho)
\otimes
\mathbf{Q}_{SPDE}(\kappa,\tau).
$$

Este modelo permite capturar simultaneamente:

- Correlacion espacial entre estaciones cercanas.
- Persistencia temporal de los episodios de contaminacion.
- Efectos explicativos de las covariables.

## 13. Perfil horario laborable/finde

Ademas del modelo horario directo, se construye un modelo de perfil horario. En
este caso, los datos horarios se agregan por:

$$
\text{estacion} \times \text{mes} \times \text{tipo de dia} \times \text{hora}.
$$

El tipo de dia distingue entre dias laborables y fines de semana:

$$
D =
\begin{cases}
0, & \text{laborable},\\
1, & \text{fin de semana}.
\end{cases}
$$

Para capturar el ciclo diario se introducen terminos armonicos:

$$
\sin\left(\frac{2\pi h}{24}\right),
\qquad
\cos\left(\frac{2\pi h}{24}\right),
$$

y tambien terminos de ciclo de 12 horas:

$$
\sin\left(\frac{2\pi h}{12}\right),
\qquad
\cos\left(\frac{2\pi h}{12}\right).
$$

Ademas, se incluyen interacciones entre el indicador de fin de semana y los
terminos armonicos, por ejemplo:

$$
D \cdot \sin\left(\frac{2\pi h}{24}\right),
\qquad
D \cdot \cos\left(\frac{2\pi h}{24}\right).
$$

El predictor lineal del perfil horario puede escribirse como:

$$
\eta(s,m,d,h)
=
\beta_0
+
\sum_{k=1}^{K}
\beta_k X_k(s,m,d,h)
+
\gamma_1 D
+
\gamma_2 M
+
f_{24}(h)
+
f_{12}(h)
+
D f_{24}(h)
+
D f_{12}(h)
+
u(s,h),
$$

donde \(M\) representa la tendencia mensual estandarizada, \(h\) la hora del dia
y \(u(s,h)\) el campo espacio-temporal asociado al perfil de 24 horas.

En este caso, la estructura AR(1) se interpreta sobre las 24 horas del perfil:

$$
u(s,h)
=
\rho u(s,h-1)
+
\xi(s,h).
$$

Este planteamiento reduce la dimension del problema horario y permite estudiar
la diferencia entre patrones laborables y de fin de semana sin ajustar un campo
independiente para cada hora real del ano.

## 14. Inferencia mediante INLA

El modelo pertenece a la familia de modelos latentes gaussianos. Sea:

$$
\mathbf{x}
=
\left(
\boldsymbol{\beta},
\mathbf{w}
\right)
$$

el campo latente, y sea \(\boldsymbol{\theta}\) el vector de hiperparametros:

$$
\boldsymbol{\theta}
=
\left(
\kappa,
\tau,
\sigma_\varepsilon^2,
\rho
\right).
$$

El objetivo bayesiano es obtener las distribuciones posteriores marginales:

$$
\pi(x_j \mid \mathbf{y})
=
\int
\pi(x_j \mid \boldsymbol{\theta},\mathbf{y})
\pi(\boldsymbol{\theta}\mid \mathbf{y})
d\boldsymbol{\theta}.
$$

INLA aproxima estas distribuciones sin utilizar MCMC. En lugar de simular
cadenas, utiliza aproximaciones de Laplace anidadas para calcular de forma
determinista las marginales posteriores de los efectos fijos, del campo latente
y de los hiperparametros.

Esto resulta especialmente adecuado para modelos SPDE, ya que la matriz de
precision del campo latente es dispersa.

## 15. Evaluacion predictiva

Una vez ajustado el modelo, se obtienen las predicciones posteriores medias
\(\widehat{y}_i\) para las observaciones del conjunto test. Se calculan tres
medidas principales.

El error cuadratico medio:

$$
RMSE
=
\sqrt{
\frac{1}{n_{test}}
\sum_{i \in test}
\left(
\widehat{y}_i-y_i
\right)^2
}.
$$

El error absoluto medio:

$$
MAE
=
\frac{1}{n_{test}}
\sum_{i \in test}
\left|
\widehat{y}_i-y_i
\right|.
$$

Y la cobertura posterior al 95 %:

$$
Cob_{95}
=
\frac{1}{n_{test}}
\sum_{i \in test}
\mathbf{1}
\left\{
y_i
\in
\left[
\widehat{y}_i - 1.96\widehat{\sigma}_i,
\widehat{y}_i + 1.96\widehat{\sigma}_i
\right]
\right\}.
$$

El RMSE y el MAE evaluan la precision puntual de la prediccion. La cobertura
evalua si la incertidumbre posterior esta bien calibrada.

Tambien se calculan criterios de informacion como DIC y WAIC. El DIC se define
como:

$$
DIC = \bar{D} + p_D,
$$

donde \(\bar{D}\) es la devianza posterior media y \(p_D\) el numero efectivo
de parametros.

El WAIC se define como:

$$
WAIC
=
-2
\left(
lppd - p_{WAIC}
\right).
$$

Estos criterios penalizan la complejidad del modelo y permiten comparar
alternativas, aunque la seleccion final prioriza la capacidad predictiva en el
conjunto test.

## 16. Interpretacion de los parametros del SPDE

El modelo proporciona parametros fisicos interpretables del campo espacial.

El rango espacial \(r\) indica la distancia aproximada hasta la que existe
correlacion espacial residual:

$$
r = \frac{\sqrt{8\nu}}{\kappa}.
$$

Un rango grande indica que las estaciones separadas por varios kilometros
comparten informacion residual. Un rango pequeno indica que la dependencia
espacial es mas local.

La varianza marginal \(\sigma^2\) del campo indica la intensidad de la
heterogeneidad espacial no explicada por las covariables:

$$
\sigma^2 = \operatorname{Var}\{u(s)\}.
$$

Si \(\sigma^2\) es elevada, significa que las covariables incluidas no explican
toda la estructura espacial del NO2.

En el modelo AR(1), el parametro \(\rho\) mide la persistencia temporal:

$$
\operatorname{Corr}\{u(s,t),u(s,t+h)\}
=
\rho^h.
$$

A partir de \(\rho\), se calcula un tiempo de recuperacion o decorrelacion:

$$
\tau
=
-\frac{1}{\log(\rho)}.
$$

Este valor representa el numero de periodos necesarios para que la correlacion
de un episodio de contaminacion caiga aproximadamente a \(1/e\), es decir, al
37 % de su valor inicial.

## 17. Modelo de residuos diario-horario

El esquema tambien incluye una etapa adicional para estudiar si las variables
que no son significativas a escala diaria pueden explicar variabilidad residual
a escala horaria.

En una primera etapa se ajusta un modelo diario:

$$
y_{sd}
=
\beta_0
+
\sum_{k \in S}
\beta_k X_{k,sd}
+
\varepsilon_{sd},
$$

donde \(s\) representa la estacion y \(d\) el dia. Se identifican las
covariables significativas a escala diaria y se guardan sus coeficientes
\(\widehat{\beta}_k\).

En la segunda etapa, se construye una pseudo-respuesta horaria:

$$
y^\ast_{sh}
=
y_{sh}
-
\sum_{k \in S}
\widehat{\beta}_k X_{k,sh}.
$$

Esta pseudo-respuesta representa la parte del NO2 horario que no queda explicada
por los efectos diarios significativos. Posteriormente se regresa \(y^\ast\)
sobre las covariables no significativas en el modelo diario y sobre posibles
retardos del trafico:

$$
y^\ast_{sh}
=
\alpha
+
\sum_{\ell \in R}
\delta_\ell Z_{\ell,sh}
+
\varepsilon_{sh}.
$$

Si alguna variable resulta significativa en esta segunda etapa, se interpreta
como evidencia de que dicha variable contiene informacion horaria que no se
detecta en la escala diaria.

## 18. Resumen del procedimiento

El flujo completo del modelo puede resumirse en los siguientes pasos:

1. Cargar el dataset correspondiente a la escala temporal de analisis.
2. Definir la respuesta \(y=\log(NO2+1)\).
3. Seleccionar la ventana temporal y separar entrenamiento y test mediante
   holdout temporal.
4. Estandarizar las covariables usando el conjunto de entrenamiento.
5. Ordenar las covariables mediante residuos parciales con modelos INLA-SPDE
   univariantes.
6. Eliminar multicolinealidad mediante VIF.
7. Retener covariables con significancia bayesiana y mejorar el conjunto final
   mediante stepwise DIC.
8. Construir la malla espacial y el objeto SPDE.
9. Ensamblar el modelo mediante `inla.stack`.
10. Ajustar el modelo solo espacial para varias mallas.
11. Diagnosticar residuos mediante ACF, Ljung-Box y Q-Q plot.
12. Ajustar el modelo espacio-temporal AR(1) si queda autocorrelacion temporal.
13. Comparar modelos mediante RMSE, MAE, cobertura, DIC, WAIC y tiempo de
    ejecucion.
14. Extraer e interpretar efectos fijos, hiperparametros, rango espacial,
    varianza espacial y persistencia temporal.

## 19. Formulacion final

La formulacion general del modelo utilizado en la tesis puede sintetizarse como:

$$
\boxed{
\log(NO2(s,t)+1)
=
\beta_0
+
\sum_{k=1}^{K}
\beta_k X_k(s,t)
+
u(s,t)
+
\varepsilon(s,t)
}
$$

con:

$$
\varepsilon(s,t)
\sim
\mathcal{N}(0,\sigma_\varepsilon^2),
$$

$$
u(s,t)
\sim
\text{campo gaussiano Matern aproximado mediante SPDE},
$$

y, cuando se detecta dependencia temporal:

$$
u(s,t)
=
\rho u(s,t-1)
+
\xi(s,t).
$$

Este enfoque permite separar la variabilidad del NO2 en tres componentes:

- La parte explicada por covariables observadas.
- La variabilidad espacial residual.
- La dependencia temporal no explicada por las covariables.

Por ello, el modelo INLA-SPDE resulta adecuado para estudiar la contaminacion
por NO2 en presencia de misalignment espacial, estructura temporal y
heterogeneidad no observada entre estaciones.

