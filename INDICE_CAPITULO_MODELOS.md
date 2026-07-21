# Índice del capítulo de modelización (modelo por modelo)

**TFM — Modelización espacio-temporal del NO₂ en Madrid mediante INLA-SPDE**

Este capítulo se organiza **modelo por modelo**: cada modelo ajustado se presenta como
una sección propia y completa. Para que las secciones sean paralelas y fáciles de
redactar, cada modelo sigue la misma **ficha común**:

> **Ficha común de cada modelo:**
> (a) objetivo y pregunta que responde ·
> (b) datos, escala y partición train/test ·
> (c) covariables y estructura del modelo ·
> (d) malla e hiperparámetros ·
> (e) resultados y diagnóstico ·
> (f) conclusión parcial.

**Leyenda de estado:**
- ✅ material ya generado en el proyecto
- 🟡 parcial
- ❌ pendiente de implementar

---

## Capítulo X. Modelización del NO₂ mediante INLA-SPDE

### X.1. Marco común de modelización
- X.1.1. Modelo jerárquico bayesiano y campo Matérn-SPDE
- X.1.2. Elementos compartidos: log(NO₂+1), estandarización, holdout temporal, PC-priors, inferencia INLA
- X.1.3. Criterios de evaluación comunes (RMSE, MAE, cobertura 95 %, DIC, WAIC)

### X.2. Modelo 0 — Estudio de simulación (recuperación de parámetros) ❌
- Objetivo: validar que el método recupera (κ, σ², ρ, β) con verdad conocida
- Diseño con `inla.qsample`, escenarios (rango, nº de estaciones, nugget)
- Métricas: sesgo y RMSE de los estimadores, cobertura de los IC
- Conclusión: fiabilidad del método con 24 estaciones

### X.3. Modelo 1 — Escala mensual (macro / interanual) ✅
- (a) Pregunta: ¿qué factores estructurales explican el NO₂ a largo plazo?
- (b) Ventana 2022–2024 train / 2025 test
- (c) Covariables seleccionadas + campo solo-espacial / AR(1)
- (d) Malla elegida e hiperparámetros
- (e) Efectos fijos, rango espacial, RMSE/cobertura, ACF/QQ
- (f) Conclusión: componente macro

### X.4. Modelo 2 — Escala diaria (estacional: invierno vs verano) ✅
- (a) Pregunta: ¿cambia la estructura entre estaciones del año?
- (b) Ventanas invierno y verano (test = últimos 7 días); marzo y noviembre como réplica
- (c) Covariables por ventana + estructura
- (d) Comparación de resoluciones de malla
- (e) Resultados por estación del año, diagnóstico de residuos
- (f) Conclusión: estabilidad estacional del patrón

### X.5. Modelo 3 — Escala horaria (corto plazo, 5 días de invierno) ✅
- (a) Pregunta: ¿qué gobierna la variabilidad de alta frecuencia?
- (b) Train días 1–4 / test día 5; nota sobre inestabilidad numérica horaria
- (c) Covariables (intensidad, temperatura, viento) + SPDE⊗AR(1)
- (d) Malla e hiperparámetros; justificación del AR(1) (ρ₁≈0,81)
- (e) Resultados, cobertura, diagnóstico
- (f) Conclusión: componente micro y persistencia temporal

### X.6. Modelo 4 — Perfil horario laborable/finde ✅
- (a) Pregunta: ¿difiere el ciclo diario entre laborables y fin de semana?
- (b) Agregación estación × mes × tipo de día × hora (ene–sep)
- (c) Armónicos sin/cos de 24 h y 12 h + interacciones con indicador de finde + campo SPDE
- (d) Malla e hiperparámetros
- (e) Perfiles estimados, diferencia laborable/finde
- (f) Conclusión: reducción de dimensión y patrón cíclico

### X.7. Modelo 5 — Residuos diario→horario (dos escalas) ✅
- (a) Pregunta: ¿hay información horaria no capturada por los efectos diarios?
- (b) Pseudo-respuesta y* = log(NO₂) − efectos diarios significativos
- (c) Regresión de y* sobre covariables no significativas + retardos de tráfico
- (d)/(e) Resultados de la segunda etapa
- (f) Conclusión + nota del regresor generado

### X.8. Modelo 6 — Diagnóstico en dos etapas del modelo horario ✅ *(script 18)*
- (a) Pregunta: ¿queda efecto de covariables descartadas en el residuo?
- (b) Un mes horario, sin retardos, sin colinealidad
- (c) Etapa 1 (SPDE con significativas) → residuo parcial → Etapa 2 (débiles)
- (e) Resultado: residuo "limpio" o no
- (f) Conclusión parcial

### X.9. Modelo 7 — Modelo combinado de residuos + retardo de tráfico ✅ *(script 20)*
- (a) Pregunta: ¿aportan las débiles + tráfico(h−1) todas juntas?
- (c) Una sola regresión SPDE de y* sobre todas las débiles y el retardo
- (e) Coeficientes γ e IC, % de residuo explicado
- (f) Conclusión: la persistencia viene del proceso atmosférico (AR1), no del retardo

### X.10. Modelo 8 — Comparación de estructuras sobre los residuos ✅ *(script 27)*
- (a) Pregunta: ¿qué estructura conviene añadir a las covariables?
- (c) Cinco modelos: A covariables · B AR(1) · C AR(2) · D SPDE espacial · E SPDE⊗AR(1)
- (e) Tabla comparativa en test (RMSE/MAE/cobertura/DIC/WAIC) + ACF/QQ del mejor
- (f) Conclusión: estructura ganadora

### X.11. Modelo 9 — SPDE sobre los residuos del modelo horario ✅ *(script 28)*
- (a) Pregunta: ¿queda estructura espacial residual tras el modelo bueno?
- (c) SPDE (malla gruesa) sobre los residuos + covariables no usadas
- (e) Rango/σ² del campo de residuos (resultado: ínfimos → sin estructura)
- (f) Conclusión: el modelo horario está bien especificado espacialmente

### X.12. Modelo 10 — Predicción espacial (modelo final aplicado) ❌
- (a) Objetivo: mapear el NO₂ en zonas no monitorizadas
- (c) Proyección de la posterior a rejilla continua de Madrid (`inla.mesh.projector`)
- (e) Mapa de concentración, mapa de incertidumbre, mapa de P(NO₂ > 40)
- (f) Zonas críticas identificadas

### X.13. Síntesis comparativa de los modelos
- X.13.1. Tabla global (escala, covariables, estructura, RMSE, cobertura, rango, ρ)
- X.13.2. Clasificación macro / micro de las covariables
- X.13.3. Modelo(s) recomendado(s) y limitaciones

---

## Notas sobre la organización

- **Modelos X.7–X.11**: son todos de **diagnóstico de la escala horaria**. Si quedan
  repetitivos, pueden agruparse bajo un único epígrafe "Modelos de diagnóstico de
  residuos (escala horaria)" con subsecciones. El desglose anterior los mantiene uno a uno.
- **Síntesis final (X.13)**: imprescindible en el formato modelo-por-modelo. Sin ella el
  capítulo queda como una lista y no como un argumento; aquí es donde se responde la
  pregunta científica del TFM.
- **Alternativa de estructura**: si el capítulo resulta muy largo, separar en
  *Metodología* (X.1–X.2) y *Resultados* (X.3–X.13).

## Correspondencia modelo ↔ código

| Sección | Modelo | Script / salidas |
|---|---|---|
| X.2  | Simulación (recuperación) | *(pendiente)* |
| X.3  | Mensual | `scripts/simulacion.R` · `outputs/simulacion/mensual/` |
| X.4  | Diario | `scripts/simulacion.R` · `outputs/simulacion/diario/` |
| X.5  | Horario | `scripts/simulacion.R` · `outputs/simulacion/horario/` |
| X.6  | Perfil laborable/finde | `scripts/simulacion.R` · `outputs/simulacion/horario_perfil/` |
| X.7  | Residuos diario→horario | `outputs/simulacion/residuos_diario_horario/` |
| X.8  | Diagnóstico dos etapas | `scripts/variogramas/18_residuos_dos_etapas_horario.R` |
| X.9  | Residuo combinado | `scripts/variogramas/20_modelo_residuo_combinado_horario.R` |
| X.10 | Comparación de estructuras | `scripts/variogramas/27_modelo_residuos_horario_5dias.R` |
| X.11 | SPDE sobre residuos | `scripts/variogramas/28_residuos_spde_modelo_horario.R` |
| X.12 | Predicción espacial | *(pendiente)* |
