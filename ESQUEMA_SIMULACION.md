# Esquema metodológico de `simulacion.R`

**TFM — Interpolación espacio-temporal de NO₂ en Madrid con modelos INLA-SPDE**
Documento explicativo del flujo, el fundamento matemático, la justificación de cada
decisión y posibles mejoras. Nivel: Máster en Estadística.

---

## 0. Idea general y una precisión de nomenclatura

El script recorre **tres escalas temporales** del NO₂ (mensual, diaria, horaria) y, para cada
una, ajusta y valida un modelo bayesiano jerárquico espacio-temporal **INLA-SPDE** con
covariables meteorológicas y de tráfico. Objetivo doble: **seleccionar** covariables/estructura
por escala y **validar** la capacidad predictiva dejando fuera un tramo temporal.

> **Precisión importante.** El fichero se llama `simulacion.R` pero **no es una simulación Monte
> Carlo** (no genera datos con parámetros *verdaderos* conocidos para medir su recuperación). Es
> **validación cruzada temporal por holdout** sobre datos reales. Conviene nombrarlo así en la
> memoria. En §9 se propone añadir una simulación con verdad conocida.

---

## 1. Fundamento matemático del modelo

### 1.1 El modelo jerárquico bayesiano en tres niveles

Sea $y_i = \log(\text{NO}_2)$ observado en la estación $s_i$ en el instante $t_i$, $i=1,\dots,n$.

**Nivel 1 — Verosimilitud (datos):**

$$y_i \mid \eta_i, \theta \;\sim\; \mathcal{N}\!\left(\eta_i,\; \sigma_\varepsilon^2\right),
\qquad \sigma_\varepsilon^2 = 1/\tau_\varepsilon$$

donde $\tau_\varepsilon$ es la **precisión del ruido de observación (nugget)**. El predictor
lineal es

$$\eta_i \;=\; \beta_0 \;+\; \sum_{k=1}^{K}\beta_k\,X_{k}(s_i,t_i)\;+\; u(s_i,t_i).$$

**Nivel 2 — Campo latente gaussiano** $\mathbf{x} = (\beta_0,\boldsymbol\beta,\mathbf{w})$, con
$\mathbf{w}$ los valores del campo espacial en los nodos de la malla:

$$\mathbf{x}\mid\boldsymbol\theta \;\sim\; \mathcal{N}\!\left(\mathbf{0},\,
\mathbf{Q}(\boldsymbol\theta)^{-1}\right),$$

un **campo aleatorio markoviano gaussiano (GMRF)** con **matriz de precisión dispersa**
$\mathbf{Q}$ (la inversa de la covarianza). La dispersión de $\mathbf{Q}$ es lo que hace viable
la inferencia.

**Nivel 3 — Hiperparámetros** $\boldsymbol\theta = (\kappa,\tau,\tau_\varepsilon,\rho)$ con
priors $\pi(\boldsymbol\theta)$.

### 1.2 El campo espacial: covarianza de Matérn

El campo continuo $u(s)$ es un **campo gaussiano estacionario e isótropo** con función de
covarianza de **Matérn**:

$$\operatorname{Cov}\!\big(u(s),u(s')\big) \;=\;
\sigma^2\,\frac{2^{1-\nu}}{\Gamma(\nu)}\,(\kappa\,\lVert s-s'\rVert)^{\nu}\,
K_\nu\!\big(\kappa\,\lVert s-s'\rVert\big),$$

donde:
- $\lVert s-s'\rVert$ = distancia euclídea (en km),
- $\kappa>0$ = **parámetro de escala** (inverso del rango),
- $\nu>0$ = **suavidad** del campo,
- $\sigma^2$ = **varianza marginal**,
- $K_\nu$ = función de Bessel modificada de segunda especie.

El **rango empírico** (distancia a la que la correlación cae a $\approx 0{,}13$) es

$$r \;=\; \frac{\sqrt{8\nu}}{\kappa}.$$

En el código se estima $r$ en km (bloque B6.5) para interpretar "hasta dónde llega" la
correlación espacial.

### 1.3 El truco SPDE: de campo continuo a GMRF disperso

Lindgren, Rue & Lindström (2011) demostraron que un campo de Matérn es la **solución
estacionaria** de la ecuación diferencial parcial estocástica (SPDE):

$$\big(\kappa^2 - \Delta\big)^{\alpha/2}\,\big(\tau\,u(s)\big) \;=\; \mathcal{W}(s),
\qquad \alpha = \nu + d/2,$$

donde $\Delta$ es el **Laplaciano**, $\mathcal{W}$ un **ruido blanco gaussiano** y $d$ la
dimensión. En 2D ($d=2$) con `alpha = 2` del código:

$$\alpha = 2 \;\Longrightarrow\; \nu = \alpha - d/2 = 2 - 1 = 1.$$

Es decir, se ajusta una **Matérn de suavidad $\nu=1$**.

**Aproximación por elementos finitos.** El campo se representa en una **malla triangular** con
funciones base $\psi_j$ lineales a trozos (valen 1 en su nodo, 0 en los demás):

$$u(s) \;\approx\; \sum_{j=1}^{m} \psi_j(s)\,w_j,
\qquad \mathbf{w}=(w_1,\dots,w_m)^\top \sim \mathcal{N}(\mathbf{0},\mathbf{Q}^{-1}).$$

La precisión $\mathbf{Q}$ se construye a partir de dos matrices dispersas de la malla:
- la **matriz de masa** $\mathbf{C}$ con $C_{ij}=\int \psi_i(s)\psi_j(s)\,ds$,
- la **matriz de rigidez** $\mathbf{G}$ con $G_{ij}=\int \nabla\psi_i(s)\cdot\nabla\psi_j(s)\,ds$.

Para $\alpha=2$:

$$\mathbf{Q}(\kappa,\tau) \;=\; \tau^2\,\big(\kappa^2\mathbf{C}+\mathbf{G}\big)\,
\mathbf{C}^{-1}\,\big(\kappa^2\mathbf{C}+\mathbf{G}\big).$$

$\mathbf{Q}$ es **dispersa** (cada nodo solo se relaciona con sus vecinos en la malla), lo que
permite factorizaciones de Cholesky rápidas — la clave computacional del método.

### 1.4 Matriz de proyección $\mathbf{A}$

Las observaciones no están en los nodos, sino en las 24 estaciones. La **matriz de proyección**
$\mathbf{A}$ ($n \times m$) evalúa las funciones base en las localizaciones:

$$A_{ij} = \psi_j(s_i),\qquad u(s_i) = (\mathbf{A}\mathbf{w})_i.$$

El predictor lineal en forma matricial es $\boldsymbol\eta = \mathbf{X}\boldsymbol\beta +
\mathbf{A}\mathbf{w}$. En INLA esto se ensambla con `inla.stack`.

### 1.5 Extensión espacio-temporal: AR1 separable (Bloque B4)

El campo se replica en cada instante $t$ y los campos consecutivos se enlazan con un proceso
**autorregresivo de orden 1**:

$$w(s,t) \;=\; \rho\,w(s,t-1) \;+\; \xi(s,t),\qquad |\rho|<1,$$

siendo $\xi(\cdot,t)$ campos espaciales de Matérn independientes. Esto genera una **covarianza
separable (producto de Kronecker)** espacio × tiempo. La precisión conjunta es

$$\mathbf{Q}_{st} \;=\; \mathbf{Q}_{\text{AR1}}(\rho)\;\otimes\;\mathbf{Q}_{\text{esp}}(\kappa,\tau),$$

donde la precisión del AR1 para $T$ instantes es la tridiagonal

$$\mathbf{Q}_{\text{AR1}}(\rho)=\frac{1}{1-\rho^2}
\begin{pmatrix}
1 & -\rho & & &\\
-\rho & 1+\rho^2 & -\rho & &\\
 & \ddots & \ddots & \ddots &\\
 & & -\rho & 1+\rho^2 & -\rho\\
 & & & -\rho & 1
\end{pmatrix}.$$

**Consecuencia práctica (memoria):** el campo latente tiene dimensión $m \times T$ (nodos ×
instantes). A escala horaria esto crece muy rápido (p. ej. $304 \times 120 \approx 36.500$
parámetros para 5 días), y factorizar $\mathbf{Q}_{st}$ es lo que consume RAM.

### 1.6 Inferencia: INLA

Se busca la marginal posterior de cada componente latente,

$$\pi(x_i\mid\mathbf{y}) \;=\; \int \pi(x_i\mid\boldsymbol\theta,\mathbf{y})\,
\pi(\boldsymbol\theta\mid\mathbf{y})\,d\boldsymbol\theta.$$

INLA (Rue, Martino & Chopin, 2009) aproxima esto **sin MCMC** mediante **aproximaciones de
Laplace anidadas**:

1. Aproxima la posterior de los hiperparámetros:
$$\tilde\pi(\boldsymbol\theta\mid\mathbf{y}) \;\propto\;
\left.\frac{\pi(\mathbf{x},\boldsymbol\theta,\mathbf{y})}
{\tilde\pi_G(\mathbf{x}\mid\boldsymbol\theta,\mathbf{y})}\right|_{\mathbf{x}=\mathbf{x}^*(\boldsymbol\theta)},$$
con $\tilde\pi_G$ la aproximación gaussiana en la moda $\mathbf{x}^*$.
2. Aproxima $\pi(x_i\mid\boldsymbol\theta,\mathbf{y})$ (estrategia *gaussian*, *simplified
   laplace* o *laplace*).
3. **Integra numéricamente** sobre una rejilla de $\boldsymbol\theta$.

Es determinista, rápido y sin diagnósticos de convergencia de cadenas: el estándar para GMRF de
gran dimensión.

---

## 2. Configuración — los "mandos"

| Elemento | Valor | Fundamento |
|---|---|---|
| `VIF_UMBRAL` | 5 | Factor de inflación de varianza. Se descarta $X_k$ si $\text{VIF}_k>5$ (criterio conservador). |
| `DIC_MEJORA_MIN` | 2 | $\Delta\text{DIC}<2$ se considera no informativo (regla estándar de comparación de modelos). |
| `config_mallas` | 4 mallas (12/8/4/1 km) | El resultado del SPDE depende de la discretización; comparar resoluciones elige el equilibrio sesgo–coste con criterio empírico. |
| `COVS_*` | 8 (mensual) / 4 (diario, horario) | Más covariables donde la agregación reduce ruido y colinealidad (macro); menos donde interesa el corto plazo (micro). |

**Factor de inflación de la varianza:**

$$\text{VIF}_k \;=\; \frac{1}{1-R_k^2},$$

donde $R_k^2$ es el coeficiente de determinación de regresar $X_k$ sobre las demás covariables.
Mide cuánto se infla $\operatorname{Var}(\hat\beta_k)$ por colinealidad.

---

## 3. Preparación de datos (Bloque B0)

Cada escala carga su dataset y produce respuesta, covariables **estandarizadas**
$\tilde X_k = (X_k-\bar X_k)/s_{X_k}$, coordenadas UTM en km, índice temporal `ID_TIEMPO_SIM` e
indicador `es_train`.

| Escala | Ventana | Train / Test |
|---|---|---|
| Mensual | 2019–2025 (usa 2022+) | Train 2022–24, Test 2025 |
| Diario | 35 días (invierno/verano) | Test = últimos 7 días |
| Horario | 5 días (invierno) | Train 4 días, Test 1 día |
| Horario-perfil | perfil laborable/finde por mes (ene–sep) | Train meses 1–7, Test 8–9 |

**Por qué holdout temporal (no aleatorio ni k-fold):** con autocorrelación temporal, un *split*
aleatorio filtra información del futuro (fuga de datos). Dejar fuera el **tramo final** imita
predecir un periodo no observado y estima honestamente el error de generalización.

**Por qué log(NO₂):** la distribución cruda es asimétrica (asimetría $\approx 0{,}88$); el
logaritmo la simetriza ($\approx 0$) y estabiliza la varianza, requisito de la verosimilitud
gaussiana. Formalmente, si $\text{NO}_2\sim\text{lognormal}$, entonces $\log(\text{NO}_2)$ es
gaussiana.

**Por qué estandarizar:** deja $\hat\beta_k$ en "efecto por desviación típica", comparables entre
sí, y mejora el condicionamiento numérico.

---

## 4. Selección de covariables (Bloques B0.5 y B1)

### 4.1 B0.5 — Residuos parciales (ordenación por importancia)

Iterativo. Partiendo de $r^{(0)}=y$, en el paso $m$:

1. Para cada covariable candidata $j$ se ajusta un SPDE univariante y se obtiene $\hat\beta_j$.
2. Se elige
$$j^* = \arg\min_{j}\; \operatorname{sd}\!\big(r^{(m-1)} - \hat\beta_j\,X_j\big).$$
3. Se actualiza el residuo:
$$r^{(m)} = r^{(m-1)} - \hat\beta_{j^*}\,X_{j^*}.$$

**Por qué dentro de un SPDE y no OLS:** al estimar $\hat\beta_j$ **con el campo espacial
presente**, los coeficientes están corregidos de **confusión espacial** (dos covariables con
patrón espacial similar no se roban el efecto mutuamente).

### 4.2 B1 — Selección formal (tres filtros)

1. **VIF hacia atrás** (elimina colinealidad; $\text{VIF}>5$).
2. **Significancia bayesiana**: se retiene $X_k$ si su IC 95 % excluye 0,
$$0 \notin \big[\,q_{0{,}025}(\beta_k),\; q_{0{,}975}(\beta_k)\,\big].$$
3. **Stepwise por DIC**: se añade/quita la variable que más reduce el DIC mientras
   $\Delta\text{DIC}\ge 2$.

---

## 5. Ajuste y evaluación (Bloques B2 y B4)

**B2** ajusta el modelo solo-espacial en las 4 mallas; **B4** el AR1, solo si B3 lo justifica.
Sobre el **test** no visto se calculan, con $\hat y_i$ la predicción posterior media y
$\hat\sigma_i$ su desviación:

$$\text{RMSE}=\sqrt{\frac{1}{n_{\text{test}}}\sum_{i}(\hat y_i-y_i)^2},\qquad
\text{MAE}=\frac{1}{n_{\text{test}}}\sum_i |\hat y_i-y_i|,$$

$$\text{Cobertura}_{95}=\frac{1}{n_{\text{test}}}\sum_i
\mathbf{1}\!\big\{\,y_i\in[\hat y_i \pm 1{,}96\,\hat\sigma_i]\,\big\}.$$

**Por qué la cobertura:** RMSE mide el error del punto, pero un modelo bayesiano debe calibrar la
**incertidumbre**; la cobertura debería aproximarse a 0,95. Es una comprobación que los modelos
puramente predictivos suelen omitir.

**Criterios de información (dentro de muestra).**

$$\text{DIC} = \bar D + p_D,\qquad \bar D = \mathbb{E}_{\mathbf{x}\mid\mathbf y}[D(\mathbf x)],
\quad p_D=\bar D - D(\bar{\mathbf x}),\quad D(\mathbf x)=-2\log p(\mathbf y\mid\mathbf x),$$

$$\text{WAIC} = -2\big(\text{lppd} - p_{\text{WAIC}}\big),\quad
\text{lppd}=\sum_i\log \mathbb{E}[p(y_i\mid\mathbf x)],\quad
p_{\text{WAIC}}=\sum_i\operatorname{Var}\big[\log p(y_i\mid\mathbf x)\big].$$

$p_D$ y $p_{\text{WAIC}}$ son el **número efectivo de parámetros** (penalización por
complejidad). El WAIC es preferible por usar toda la posterior predictiva.

**CPO (validación *leave-one-out* interna):**

$$\text{CPO}_i = \pi(y_i\mid \mathbf{y}_{-i}),\qquad
\text{LPML}=\sum_i \log \text{CPO}_i.$$

Valores bajos de $\text{CPO}_i$ señalan observaciones influyentes o mal predichas.

**Parámetros físicos del campo** (bloque B6.5): rango $r=\sqrt{8\nu}/\kappa$ y varianza marginal
$\sigma^2$, obtenidos de la posterior de los hiperparámetros vía `inla.spde.result`.

---

## 6. Diagnóstico de residuos y "gate" del AR1 (Bloques B3, B5)

Sobre los residuos se calcula la **función de autocorrelación** y el contraste de
**Ljung–Box** a $L$ retardos:

$$\hat\rho_k = \frac{\sum_{t}(r_t-\bar r)(r_{t-k}-\bar r)}{\sum_t (r_t-\bar r)^2},\qquad
Q_{LB}=n(n+2)\sum_{k=1}^{L}\frac{\hat\rho_k^2}{n-k}\;\sim\;\chi^2_L .$$

Si $Q_{LB}$ es significativo (queda autocorrelación temporal) → **se activa el AR1 (B4)**. El
Q-Q plot comprueba la normalidad de los residuos.

**Por qué este *gate*:** parsimonia — no se añade la complejidad temporal salvo que los datos la
exijan. **Matiz crítico:** el test se aplica sobre la *media espacial* de los residuos por
instante, lo que **pierde potencia**; un diagnóstico intra-estación (ver §9) da
$\hat\rho_1\approx 0{,}81$, es decir, el AR1 sí está justificado a nivel horario.

---

## 7. Comparación, interpretación y recuperación (Bloques B6–B7)

- **B6:** tabla y gráficos de RMSE/DIC/WAIC/tiempo → se elige el modelo de menor RMSE.
- **B6.5:** tablas de efectos fijos (con IC), hiperparámetros y parámetros físicos del SPDE.
- **B6b:** compara efectos entre el modelo espacial y el AR1 (¿alguna covariable era un
  artefacto de no modelar el tiempo?).
- **B7 — Tiempo de recuperación.** En un AR1, $\operatorname{corr}(t,t+k)=\rho^{k}$. El tiempo de
  *decorrelación* $\tau$ (caída a $1/e$) cumple $\rho^\tau = e^{-1}$, de donde

$$\boxed{\;\tau = -\dfrac{1}{\ln \rho}\;}$$

Interpretación: nº de periodos que tarda un episodio de contaminación en "olvidarse". Útil para
gestión de calidad del aire.

**Clasificación MACRO vs MICRO:** cruzando la importancia de cada covariable entre escalas se
etiqueta como **macro** (significativa a escala mensual → estacional/interanual) o **micro** (solo
diaria/horaria → corto plazo). Es la síntesis que responde a la pregunta científica del TFM.

---

## 8. Limitaciones conocidas (para anticipar al tribunal)

1. **"Simulación" = validación por holdout**, no Monte Carlo con verdad conocida (§0, §9.1).
2. **Inestabilidad numérica horaria.** El modelo solo-espacial (B0.5/B2) amontona miles de
   observaciones horarias sobre 24 puntos en un único campo estático; el sistema queda **mal
   condicionado** (no separa campo de nugget) y el solver **se cae** por encima de ~4.600
   observaciones.
3. **Potencia del *gate* AR1** (test sobre la media espacial): el AR1 casi nunca se dispara pese a
   que sí es necesario.
4. **Holdout único** → RMSE sin intervalo de confianza.
5. **Priors por defecto** en `inla.spde2.matern`.

---

## 9. Posibles mejoras (ordenadas por impacto)

### 9.1 Añadir un estudio de simulación real (recuperación de parámetros)
Generar un campo con **verdad conocida** $(\kappa_0,\tau_0,\rho_0,\boldsymbol\beta_0)$:
$$\mathbf{w}^{\text{true}}\sim\mathcal{N}(\mathbf 0, \mathbf Q(\kappa_0,\tau_0)^{-1})
\;\;(\texttt{inla.qsample}),\qquad
y_i = \mathbf{x}_i^\top\boldsymbol\beta_0 + (\mathbf A\mathbf w^{\text{true}})_i + \varepsilon_i,$$
ajustar el modelo y medir **sesgo** $\mathbb E[\hat\theta]-\theta_0$ y **cobertura** de los IC,
repetido en $R$ réplicas. Responde: *¿con 24 estaciones el SPDE recupera bien el rango?*

### 9.2 Resolver el crash horario con AR1 desde el inicio
Sustituir el solo-espacial por la estructura espacio-temporal desde B0.5/B2 (evita el mal
condicionamiento), o trabajar sobre el **perfil laborable/finde** (reduce ~10.000 a ~430 filas).

### 9.3 Diagnóstico intra-estación + descomposición de varianza
Sustituir el Ljung-Box sobre la media espacial por **ACF por estación**. Prueba realizada:
$\hat\rho_1=0{,}81$, las 24 estaciones rechazan ruido blanco. Descomponer además
$$\operatorname{Var}(y)\approx
\underbrace{\operatorname{Var}(\mathbf X\boldsymbol\beta)}_{\text{covariables}}
+\underbrace{\sigma^2_{\text{esp}}}_{\text{campo}}
+\underbrace{\sigma^2_\varepsilon}_{\text{nugget}}.$$
(En enero salió ≈ 47 % / 11 % / 42 %.)

### 9.4 Selección por tamaño de efecto a nivel horario
Con $N\approx 17.500$ todo es significativo; usar un umbral de $|\hat\beta_{\text{estand.}}|$
para distinguir lo relevante de lo despreciable.

### 9.5 Priors de Complejidad Penalizada (PC-priors)
Con `inla.spde2.pcmatern`, priors de la forma
$$P(r < r_0)=\alpha_r,\qquad P(\sigma>\sigma_0)=\alpha_\sigma,$$
que penalizan la desviación respecto a un modelo base, regularizan y evitan que
$\tau_\varepsilon$ diverja (una causa del crash).

### 9.6 Predecir fuera del *stack* de estimación
Ajustar solo con las observaciones y proyectar después con `inla.mesh.projector`; reduce mucho la
memoria del AR1.

### 9.7 Validación con réplicas (*rolling-origin*)
Repetir el holdout con varios orígenes temporales para dar **IC del RMSE**.

### 9.8 Modo compacto / PARDISO
`inla.setOption(inla.mode="compact")` y el solver PARDISO aceleran y ahorran memoria.

### 9.9 Covariables retardadas (resultado negativo, probado)
Añadir tráfico en $h-1$ **no** resultó significativo una vez incluido el contemporáneo (colineales
entre sí). La persistencia del NO₂ no viene de la memoria del tráfico, sino del proceso
atmosférico → se modela con **AR1**, no con retardos.

---

## 10. Referencias

- Lindgren, F., Rue, H., Lindström, J. (2011). *An explicit link between Gaussian fields and
  Gaussian Markov random fields: the SPDE approach.* JRSS-B, 73(4), 423–498.
- Rue, H., Martino, S., Chopin, N. (2009). *Approximate Bayesian inference for latent Gaussian
  models using integrated nested Laplace approximations.* JRSS-B, 71(2), 319–392.
- Krainski, E. T. et al. (2019). *Advanced Spatial Modeling with Stochastic Partial Differential
  Equations using R and INLA.* Chapman & Hall/CRC.
- Blangiardo, M., Cameletti, M. (2015). *Spatial and Spatio-temporal Bayesian Models with R-INLA.*
  Wiley.
- Simpson, D. et al. (2017). *Penalising model component complexity: A principled, practical
  approach to constructing priors (PC-priors).* Statistical Science, 32(1), 1–28.
- Watanabe, S. (2010). *Asymptotic equivalence of Bayes cross validation and WAIC.* JMLR, 11.
