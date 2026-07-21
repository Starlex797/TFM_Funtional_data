# Esquema del Marco Teórico

> **Modelización espacio-temporal bayesiana del NO₂ en Madrid mediante INLA-SPDE**

Este esquema ordena el marco teórico del TFM siguiendo el hilo lógico del propio
análisis: de lo **determinista y local** (interpolar clima) a lo **estocástico y
jerárquico** (campo espacial + inferencia bayesiana). Cada apartado indica **qué
teoría desarrollar**, **por qué es necesaria en el proyecto** y las **referencias
canónicas**. La numeración es propuesta de índice.

---

## 1. Fundamentos de la estadística espacial

Establece el vocabulario y la taxonomía sobre la que se apoya todo lo demás.

- **1.1. Tipos de datos espaciales.** Datos geoestadísticos (medidos en puntos
  fijos: las 24 estaciones), datos de área/lattice (barrios, tráfico agregado) y
  patrones de puntos. Justificar que el NO₂ es un **proceso geoestadístico**.
- **1.2. El proceso estocástico espacial** $\{Z(s): s \in D \subset \mathbb{R}^2\}$.
  Realización única, necesidad de supuestos que sustituyan a la repetición.
- **1.3. Estacionariedad e isotropía.** Estacionariedad estricta, débil y
  **intrínseca**; isotropía. Por qué se necesitan para poder estimar con una sola
  realización.
- **1.4. Dependencia espacial y el "primer axioma de la geografía".** Autocorrelación
  positiva a corta distancia → motiva tanto la interpolación como el campo latente.
- **1.5. Problemas específicos del dato:** *misalignment* espacial (clima y NO₂ no se
  miden en los mismos puntos → hace falta interpolar), *change of support* y el
  sesgo ecológico (tráfico por barrio vs. NO₂ puntual).

*Referencias:* Cressie (1993); Banerjee, Carlin & Gelfand (2015); Diggle & Ribeiro (2007).

---

## 2. Interpolación espacial determinista (covariables climáticas)

Base teórica de la fase de **construcción de covariables**: el clima se mide en pocas
estaciones meteorológicas y hay que llevarlo a cada estación de NO₂.

- **2.1. El problema de la predicción espacial.** Estimar $\hat Z(s_0)$ en un punto no
  observado a partir de $Z(s_1),\dots,Z(s_n)$.
- **2.2. Métodos de referencia (benchmark).** Media global, vecino más próximo (1-NN)
  y k-vecinos (kNN) como interpoladores triviales contra los que comparar.
- **2.3. Ponderación por distancia inversa (IDW).** Formulación
  $\hat Z(s_0)=\sum_i w_i Z(s_i)$ con $w_i \propto \lVert s_0-s_i\rVert^{-p}$;
  el papel del **exponente $p$** (se usa $p=2$), la interpolación exacta y sus
  limitaciones (efecto "ojo de buey", sin modelo de incertidumbre nativo).
- **2.4. Cuantificación empírica de la incertidumbre del IDW.** LOOCV local ponderado
  como sustituto del error de predicción del que IDW carece por construcción.
- **2.5. Criterio de selección del interpolador.** Validación **LOOCV** con
  **RMSE/MAE**; justificación de por qué IDW ($p=2$) se elige como método operativo
  para generar las covariables climáticas del modelo.
- **2.6. (Opcional) Puente hacia lo estocástico.** IDW como caso particular
  determinista frente al **kriging**, que sí modela la covarianza → enlaza con §4.

*Referencias:* Shepard (1968); Li & Heap (2014); Cressie (1993, cap. kriging).

---

## 3. Modelos jerárquicos bayesianos

Núcleo conceptual del trabajo: por qué un modelo en **capas** es la forma natural de
separar señal física, estructura espacial y ruido.

- **3.1. Inferencia bayesiana: fundamentos.** Verosimilitud, previa, posterior,
  teorema de Bayes; interpretación de la incertidumbre como distribución.
- **3.2. La estructura jerárquica en tres niveles** (Berliner):
  - **Nivel de datos (verosimilitud):** $y_i \mid \eta_i,\theta \sim \mathcal{N}(\eta_i,\sigma_\varepsilon^2)$.
  - **Nivel de proceso latente:** predictor lineal
    $\eta_i = \beta_0 + \sum_k \beta_k X_k(s_i,t_i) + u(s_i,t_i)$.
  - **Nivel de hiperparámetros:** priors sobre $\theta=(\kappa,\tau,\sigma_\varepsilon^2,\rho)$.
- **3.3. Efectos fijos vs. efectos latentes.** Covariables observadas (β) frente al
  campo no observado $u(\cdot)$; qué representa cada componente.
- **3.4. Distribuciones previas.** Priors no informativas vs. informativas; introducir
  los **Penalised Complexity priors (PC-priors)** como priors que penalizan la
  complejidad y regularizan (se usan sobre rango y varianza del SPDE).
- **3.5. Descomposición de la varianza.** Partición: covariables / campo espacial /
  ruido (*nugget*); base para interpretar cuánto explica cada capa.

*Referencias:* Gelman et al. (2013), *Bayesian Data Analysis*; Banerjee et al. (2015);
Blangiardo & Cameletti (2015).

---

## 4. Campos aleatorios gaussianos y geoestadística

Especifica **cómo** se modela la capa espacial latente $u(s)$.

- **4.1. Campos aleatorios gaussianos (GRF).** Definición, caracterización por media y
  función de covarianza; por qué la gaussianidad hace tratable la inferencia.
- **4.2. Funciones de covarianza y la familia Matérn.**
  $$C(h)=\sigma^2\frac{2^{1-\nu}}{\Gamma(\nu)}(\kappa h)^\nu K_\nu(\kappa h)$$
  Parámetros: **varianza marginal** $\sigma^2$, **escala** $\kappa$, **suavidad** $\nu$.
- **4.3. Rango espacial.** $r=\sqrt{8\nu}/\kappa$: distancia a la que la correlación se
  vuelve despreciable; su interpretación física en el NO₂ urbano.
- **4.4. El variograma.** Semivariograma empírico y teórico; *nugget*, *sill*, rango;
  su uso como **diagnóstico** de estructura espacial residual (scripts `variogramas/`).
- **4.5. Kriging (predicción óptima).** Kriging ordinario/universal como el mejor
  predictor lineal insesgado (BLUP); relación con el modelo latente bayesiano.
- **4.6. Autocorrelación espacial en datos de área.** **Índice de Moran** global y
  local (LISA) como test de existencia de estructura espacial que justifica el campo.

*Referencias:* Cressie (1993); Rasmussen & Williams (2006) para Matérn; Stein (1999).

---

## 5. Modelos Latentes Gaussianos (LGM)

Marco unificador que engloba el modelo del TFM y habilita INLA.

- **5.1. Definición de la clase LGM.** Verosimilitud de la familia exponencial +
  campo latente gaussiano $x$ + hiperparámetros $\theta$ (de dimensión baja).
- **5.2. Campos aleatorios de Markov gaussianos (GMRF).** Independencia condicional →
  **matriz de precisión dispersa** $Q$; por qué la dispersión es la clave computacional.
- **5.3. El modelo del TFM como LGM.** Encajar
  $y = X\beta + Aw + \varepsilon$ en la clase: $x=(\beta,w)$, $\theta=(\kappa,\tau,\sigma_\varepsilon^2,\rho)$.

*Referencias:* Rue & Held (2005), *Gaussian Markov Random Fields*; Rue, Martino & Chopin (2009).

---

## 6. El enfoque SPDE (Lindgren–Rue–Lindström)

El puente teórico que convierte un GRF Matérn continuo en un GMRF disperso computable.

- **6.1. La idea central.** Un campo Matérn es solución de la **ecuación diferencial
  parcial estocástica**
  $(\kappa^2-\Delta)^{\alpha/2}(\tau u(s))=\mathcal{W}(s)$,
  con $\alpha=\nu+d/2$; en 2D con $\alpha=2 \Rightarrow \nu=1$.
- **6.2. Discretización por elementos finitos.** Malla triangular (*mesh*),
  funciones base $\psi_j$ y representación $u(s)\approx\sum_j \psi_j(s)\,w_j$.
- **6.3. Construcción de la malla.** `max.edge`, `cutoff`, borde exterior; compromiso
  **resolución ↔ coste computacional**; comparación de mallas (gruesa→fina).
- **6.4. La matriz de proyección $A$ y el `inla.stack`.** $A_{ij}=\psi_j(s_i)$ lleva el
  campo de los nodos a las estaciones; ensamblado de respuesta, β y campo latente.
- **6.5. PC-priors sobre el SPDE.** $P(r<r_0)=p$, $P(\sigma>\sigma_0)=p$; su papel
  regularizador sobre rango y varianza del campo.

*Referencias:* Lindgren, Rue & Lindström (2011, JRSS-B); Bakka et al. (2018);
Krainski et al. (2019), *Advanced Spatial Modeling with SPDEs and R-INLA*.

---

## 7. Inferencia mediante INLA

Cómo se estima el modelo sin MCMC.

- **7.1. El objetivo inferencial.** Marginales posteriores
  $\pi(x_j\mid y)=\int \pi(x_j\mid\theta,y)\,\pi(\theta\mid y)\,d\theta$.
- **7.2. Aproximación de Laplace anidada (INLA).** Aproximación de $\pi(\theta\mid y)$ y
  de las condicionales latentes; exploración de la rejilla de hiperparámetros.
- **7.3. INLA vs. MCMC.** Ventajas (velocidad, determinismo, ideal para GMRF dispersos)
  y límites (clase LGM, θ de baja dimensión); por qué es el método adecuado aquí.
- **7.4. Salidas del ajuste.** `summary.fixed` (β en escala log), hiperparámetros
  (rango, varianza, ρ, precisión del *nugget*) y su lectura.

*Referencias:* Rue, Martino & Chopin (2009, JRSS-B); Blangiardo & Cameletti (2015);
Gómez-Rubio (2020), *Bayesian Inference with INLA*.

---

## 8. Extensión espacio-temporal

Añade la dimensión temporal que exige el dato horario/diario.

- **8.1. Motivación empírica.** Autocorrelación temporal residual (ACF, Ljung-Box):
  el modelo solo espacial deja memoria temporal sin capturar.
- **8.2. Estructura autorregresiva AR(1) en el campo latente.**
  $u(s,t)=\rho\,u(s,t-1)+\xi(s,t)$, con $\xi$ de estructura Matérn.
- **8.3. Separabilidad y producto de Kronecker.**
  $Q_{ST}=Q_{AR1}(\rho)\otimes Q_{SPDE}(\kappa,\tau)$; supuesto de separabilidad
  espacio-temporal y sus implicaciones.
- **8.4. Interpretación de $\rho$.** Persistencia $\text{Corr}=\rho^h$ y tiempo de
  decorrelación $\tau=-1/\log\rho$ (caída a $1/e$).
- **8.5. Modelización del ciclo diario.** Términos **armónicos** (Fourier) de 24 h y
  12 h e interacciones laborable/fin de semana para el perfil horario.

*Referencias:* Cameletti et al. (2013); Blangiardo & Cameletti (2015, cap. espacio-temporal).

---

## 9. Selección de variables y validación

Metodología estadística que rodea al modelo.

- **9.1. Multicolinealidad.** Factor de Inflación de la Varianza (VIF) y eliminación
  hacia atrás (umbral 5).
- **9.2. Selección de covariables.** Residuos parciales y **tamaño del efecto** (β
  estandarizado) frente a la significancia por N grande; significancia bayesiana
  (IC creíble 95 % que no contiene 0) y **stepwise por DIC** (ΔDIC ≥ 2).
- **9.3. Correlaciones no sesgadas.** Anomalías (desestacionalizar), análisis por
  estación, Pearson + Spearman + **distancia de correlación**, IC por *block bootstrap*.
- **9.4. Criterios de información.** **DIC** y **WAIC**: ajuste dentro de muestra
  penalizando complejidad ($\bar D + p_D$; lppd $- p_{WAIC}$).
- **9.5. Validación predictiva fuera de muestra.**
  - *Holdout temporal* (predecir el tramo final como NA).
  - *Validación cruzada espacial por bloques* (ocultar estaciones completas).
  - Métricas: **RMSE, MAE, cobertura del IC 95 %, CRPS**.
- **9.6. Diagnóstico de residuos.** QQ-plot (normalidad), ACF (memoria temporal),
  variograma/Moran (estructura espacial residual): un buen modelo deja residuos
  **sin estructura**.

*Referencias:* Spiegelhalter et al. (2002) DIC; Watanabe (2010) WAIC; Gneiting & Raftery (2007) CRPS.

---

## 10. Contexto de aplicación: contaminación por NO₂

Cierra el marco conectando la estadística con el fenómeno físico.

- **10.1. Química y dinámica del NO₂ urbano.** Emisión (tráfico) vs. eliminación
  (dispersión por viento, mezcla vertical, fotoquímica).
- **10.2. Mecanismos que justifican las covariables.** Viento (−), temperatura (−),
  presión (+), radiación, tráfico (+); naturaleza multi-escala (micro horaria vs.
  macro estacional).
- **10.3. Justificación de la transformación** $\log(NO_2+1)$: positividad, asimetría,
  naturaleza multiplicativa (Box-Cox con λ ≈ 0).
- **10.4. Tipologías de estación** (tráfico, fondo urbano, suburbana) y su papel en la
  heterogeneidad espacial que captura el campo SPDE.
- **10.5. Marco normativo** (opcional): valores límite de NO₂ (UE / OMS) como
  motivación aplicada.

*Referencias:* Seinfeld & Pandis (2016), *Atmospheric Chemistry and Physics*;
normativa EU 2008/50/CE; directrices OMS.

---

## Mapa lógico del marco (resumen)

```
§1 Fundamentos espaciales
        │
        ├─► §2 Interpolación determinista (IDW)  ──►  covariables climáticas
        │
        ├─► §3 Modelos jerárquicos bayesianos  ─┐
        │                                        │
        ├─► §4 Campos gaussianos / Matérn  ──────┤
        │                                        ├─►  §5 LGM  ─►  §6 SPDE  ─►  §7 INLA
        │                                        │                                │
        └─► §8 Extensión espacio-temporal AR(1) ─┘                                │
                                                                                  ▼
                              §9 Selección + validación  ◄──────  modelo ajustado
                                                                                  │
                                                                                  ▼
                                                        §10 Interpretación física (NO₂)
```

**Hilo argumental en una frase:** de interpolar el clima con un método local
determinista (**IDW**, §2), a modelar el NO₂ como un **modelo jerárquico bayesiano**
(§3) cuya capa latente es un **campo gaussiano Matérn** (§4) que, reformulado como
**LGM** (§5) vía **SPDE** (§6), se estima eficientemente con **INLA** (§7), se extiende
en el tiempo con **AR(1)** (§8), y se selecciona/valida (§9) para interpretar
físicamente la contaminación (§10).
