# ==============================================================================
# DIAGNÓSTICO EN DOS ETAPAS — ¿QUÉ NOS ESTAMOS DEJANDO EN LOS RESIDUOS?
# Nivel horario · un mes · sin retardos · sin multicolinealidad
# ==============================================================================
# Idea:
#   ETAPA 1 — Modelo INLA-SPDE con las covariables candidatas (previamente
#             depuradas de multicolinealidad por VIF). Se identifican las
#             SIGNIFICATIVAS (IC 95% del efecto fijo excluye 0) y se construye
#             el residuo parcial restando SOLO su contribución:
#                 y*_sh = y_sh − Σ_j  β̂_j · X_j,sh     (j = covariables sig.)
#
#   ETAPA 2 — Se regresa y* sobre las covariables que NO salieron significativas
#             (sin retardos). Si alguna de esas variables "descartadas" todavía
#             explica el residuo, es señal que nos estábamos dejando dentro.
#
# Objetivo: comprobar si el residuo del modelo horario está "limpio" de
#   covariables (→ lo que quede es estructura espacial/temporal, que capturan
#   el campo SPDE y el AR1) o si aún esconde efecto de variables descartadas.
#
# Nota metodológica: es un residuo de dos etapas (no un Frisch-Waugh-Lovell
#   completo). Los β̂ de la Etapa 1 se tratan como conocidos, por lo que la
#   significancia de la Etapa 2 es orientativa (problema del regresor generado).
# ==============================================================================
# QUÉ HACE:
#
# Este script realiza un diagnóstico avanzado de los residuos de un modelo
# horario INLA-SPDE para un mes concreto de 2025.
#
# El mes se selecciona mediante:
#
#     MES_ANALISIS <- 1L
#
# Por defecto se analiza enero.
#
# PREPARACIÓN:
#
# - Carga el dataset maestro horario de 2025.
# - Selecciona el mes indicado.
# - Trabaja con log(NO2) como variable respuesta.
# - Estandariza dentro del mes las covariables:
#     · Intensidad del tráfico.
#     · Carga de tráfico.
#     · Temperatura.
#     · Humedad relativa.
#     · Precipitaciones.
#     · Presión barométrica.
#     · Radiación solar.
#     · Velocidad del viento.
#
# CONTROL DE MULTICOLINEALIDAD:
#
# - Calcula el VIF de las covariables candidatas.
# - Aplica una eliminación hacia atrás.
# - Excluye variables hasta que el VIF máximo sea inferior a 5.
#
# CONSTRUCCIÓN ESPACIAL:
#
# - Construye una malla espacial.
# - Define un campo gaussiano Matérn mediante SPDE.
# - Proyecta las observaciones sobre la malla.
#
# ETAPA 1:
#
# - Ajusta un modelo INLA-SPDE con las covariables que superan el filtro VIF.
# - Debido al gran número de observaciones horarias, clasifica las covariables
#   según el tamaño de su coeficiente estandarizado:
#
#     · Fuertes: |beta| >= 0,10.
#     · Débiles:  |beta| < 0,10.
#
# - Construye un residuo parcial eliminando la contribución de las variables
#   consideradas fuertes.
# - Compara la desviación típica original con la del residuo parcial.
#
# ETAPA 2:
#
# - Ajusta otro modelo INLA-SPDE sobre el residuo parcial.
# - Utiliza como predictores las covariables clasificadas como débiles.
# - Comprueba si alguna variable débil todavía explica parte del residuo.
# - Calcula cuánto disminuye la desviación típica después de esta segunda etapa.
#
# DESCOMPOSICIÓN ESPACIAL DEL RESIDUO:
#
# - Estima aproximadamente la contribución de:
#     · Covariables fuertes.
#     · Campo espacial SPDE.
#     · Ruido no estructurado o nugget.
#
# - Extrae el rango espacial estimado por el modelo SPDE.
# - Compara el WAIC de:
#     · Un modelo con campo espacial.
#     · Un modelo sin campo espacial.
#
# - Calcula el índice de Moran de los residuos medios por estación.
# - Comprueba si todavía queda autocorrelación espacial.
#
# ETAPA 3: RETARDO DEL TRÁFICO
#
# - Construye las variables:
#     · Intensidad del tráfico en la hora anterior.
#     · Carga de tráfico en la hora anterior.
#
# - Mantiene únicamente observaciones consecutivas separadas por una hora.
# - Ajusta un nuevo modelo INLA-SPDE sobre el residuo parcial.
# - Evalúa si el tráfico de la hora anterior explica parte del residuo.
# - Compara la reducción de variabilidad conseguida por:
#     · Las covariables débiles.
#     · El tráfico retardado una hora.
#
# DIAGNÓSTICO TEMPORAL:
#
# - Calcula la función de autocorrelación de los residuos hasta 24 horas.
# - Calcula la ACF media entre estaciones.
# - Aplica la prueba de Ljung-Box a cada estación.
# - Determina si el residuo es ruido blanco o conserva estructura temporal.
#
# FINALIDAD PARA EL TFM:
#
# Este script intenta responder a la pregunta:
#
#     ¿Qué estructura queda sin explicar después de incluir las covariables?
#
# Permite determinar si los residuos conservan:
#
# - Efectos de covariables inicialmente consideradas débiles.
# - Dependencia espacial.
# - Efectos retardados del tráfico.
# - Autocorrelación temporal.
# - Ruido no estructurado.
#
# Los resultados permiten justificar componentes adicionales del modelo:
#
# - El campo espacial SPDE, si mejora el WAIC y existe dependencia espacial.
# - Un término AR(1), si los residuos tienen autocorrelación temporal.
# - Retardos del tráfico, si las variables de la hora anterior reducen el
#   residuo.
#
# OBSERVACIÓN METODOLÓGICA:
#
# El análisis en dos etapas es un diagnóstico exploratorio. Los coeficientes de
# la primera etapa se tratan como conocidos cuando se construye el residuo, por
# lo que la significación obtenida en la segunda etapa debe interpretarse de
# forma orientativa.
#
# Aunque el nombre del script hace referencia a dos etapas, el archivo incluye
# también una tercera etapa dedicada al tráfico retardado y un diagnóstico de
# autocorrelación temporal.
#
# SALIDAS:
#
# outputs/analysis/residuos_dos_etapas/horario_{mes}/
#
# Salidas principales:
#
# - Estado de las covariables y resultado del filtro VIF.
# - Resultados de la segunda etapa.
# - Gráfico de reducción de la desviación típica.
# - Descomposición aproximada de la varianza.
# - Resultados del tráfico retardado.
# - Comparación entre covariables débiles y retardo.
# - Gráfico de la ACF temporal del residuo.


library(INLA)
library(data.table)
library(sf)
library(car)
library(spdep)
library(here)
library(ggplot2)
library(gridExtra)
library(grid)

set.seed(4827)

# ------------------------------------------------------------------------------
# PARÁMETROS
# ------------------------------------------------------------------------------
MES_ANALISIS <- 1L # <-- mes a analizar (1 = enero ... 9 = septiembre) · 2025
VIF_UMBRAL <- 5 # umbral de multicolinealidad
# A nivel horario N es enorme (~17.500 obs), así que TODO sale estadísticamente
# significativo (IC95% excluye 0). El filtro operativo es el TAMAÑO DEL EFECTO:
# se consideran "fuertes" (Etapa 1) las covariables con |β estandarizado| >= UMBRAL,
# y el resto ("débiles") pasa a la Etapa 2 para ver si aún explican el residuo.
UMBRAL_BETA <- 0.10 # |beta estandarizado| mínimo para ser "fuerte"

# Covariables candidatas (columnas crudas → se estandarizan dentro del mes)
COVS_RAW <- c(
  "intensidad_raw", "carga_raw", "Temperatura_raw",
  "Humedad_Relativa_raw", "Precipitaciones_raw",
  "Presion Barométrica_raw", "Radiación Solar_raw",
  "Velocidad Viento_raw"
)
COVS_ALIAS <- c(
  "intensidad", "carga", "temperatura", "humedad",
  "precipitacion", "presion", "radiacion", "viento"
)

meses_nombre <- c(
  "enero", "febrero", "marzo", "abril", "mayo", "junio",
  "julio", "agosto", "septiembre", "octubre", "noviembre", "diciembre"
)
mes_lab <- meses_nombre[MES_ANALISIS]

carpeta_out <- here(
  "outputs", "residuos_dos_etapas",
  sprintf("horario_%02d_%s", MES_ANALISIS, mes_lab)
)
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Guardar tabla como PNG (sin navegador)
# ------------------------------------------------------------------------------
guardar_tabla_png <- function(df, titulo, subtitulo = NULL, ruta_png, ancho_px = 900) {
  df <- as.data.frame(df)
  th <- gridExtra::ttheme_minimal(
    core = list(
      fg_params = list(fontsize = 9),
      bg_params = list(fill = c("white", "#F5F5F5"), col = NA)
    ),
    colhead = list(
      fg_params = list(fontsize = 9, fontface = "bold"),
      bg_params = list(fill = "#DDEEFF", col = NA)
    )
  )
  tg <- gridExtra::tableGrob(df, rows = NULL, theme = th)
  tit <- grid::textGrob(titulo, gp = grid::gpar(fontsize = 13, fontface = "bold"))
  if (!is.null(subtitulo)) {
    sub <- grid::textGrob(subtitulo, gp = grid::gpar(fontsize = 9, col = "grey40"))
    comb <- gridExtra::arrangeGrob(tit, sub, tg,
      nrow = 3,
      heights = grid::unit(c(0.55, 0.35, nrow(df) * 0.32 + 0.5), "inches")
    )
    alto <- 0.55 + 0.35 + nrow(df) * 0.32 + 0.5
  } else {
    comb <- gridExtra::arrangeGrob(tit, tg,
      nrow = 2,
      heights = grid::unit(c(0.55, nrow(df) * 0.32 + 0.5), "inches")
    )
    alto <- 0.55 + nrow(df) * 0.32 + 0.5
  }
  ggplot2::ggsave(ruta_png, comb,
    width = max(6, ancho_px / 96),
    height = max(3, alto), dpi = 96, bg = "white"
  )
  cat("  Tabla guardada:", basename(ruta_png), "\n")
}

# ==============================================================================
# 1. CARGA Y PREPARACIÓN (mes horario)
# ==============================================================================
cat(sprintf(
  "\n== Diagnóstico de residuos en dos etapas — %s 2025 (horario) ==\n\n",
  toupper(mes_lab)
))

dt <- readRDS(here(
  "data", "processed", "Maestro", "horario",
  "dataset_maestro_inla_2025_HORARIO.rds"
))
setDT(dt)
setnames(dt, "LOG_NO2_HORARIO", "LOG_NO2")
dt[, MES := as.integer(format(FECHA, "%m"))]
dt <- dt[MES == MES_ANALISIS & !is.na(LOG_NO2)]

# Coordenadas UTM en km
coords_u <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
csf <- st_as_sf(coords_u, coords = c("LONGITUD", "LATITUD"), crs = 4326) |>
  st_transform(25830)
coords_u[, X_km := st_coordinates(csf)[, 1] / 1000]
coords_u[, Y_km := st_coordinates(csf)[, 2] / 1000]
dt <- merge(dt, coords_u[, .(ESTACION, X_km, Y_km)], by = "ESTACION")

# Estandarizar covariables DENTRO del mes (alias sin caracteres especiales)
covs_disp <- character(0)
for (i in seq_along(COVS_RAW)) {
  v <- COVS_RAW[i]
  a <- COVS_ALIAS[i]
  if (v %in% names(dt)) {
    dt[, (a) := as.numeric(scale(get(v)))]
    covs_disp <- c(covs_disp, a)
  }
}

# Filas completas en respuesta + covariables
dt <- dt[complete.cases(dt[, c("LOG_NO2", covs_disp), with = FALSE])]
cat(sprintf(
  "Observaciones: %d | estaciones: %d | horas: %d\n",
  nrow(dt), uniqueN(dt$ESTACION), uniqueN(dt$HORA)
))

# ==============================================================================
# 2. CONTROL DE MULTICOLINEALIDAD (VIF hacia atrás sobre TODAS las candidatas)
# ==============================================================================
cat("\n--- VIF hacia atrás (umbral =", VIF_UMBRAL, ") ---\n")
vars_vif <- covs_disp
vif_final <- NULL
repeat {
  if (length(vars_vif) < 2) break
  fml <- as.formula(paste("LOG_NO2 ~", paste(vars_vif, collapse = " + ")))
  vifs <- car::vif(lm(fml, data = dt))
  vif_final <- vifs
  if (max(vifs) <= VIF_UMBRAL) {
    cat("  -> Sin multicolinealidad.\n")
    break
  }
  drop <- names(which.max(vifs))
  cat(sprintf("  -> Elimino '%s' (VIF = %.2f)\n", drop, max(vifs)))
  vars_vif <- setdiff(vars_vif, drop)
}
vars_excluidas_vif <- setdiff(covs_disp, vars_vif)
cat(sprintf("  Covariables tras VIF: %s\n", paste(vars_vif, collapse = ", ")))
if (length(vars_excluidas_vif) > 0) {
  cat(sprintf(
    "  Excluidas por colinealidad: %s\n",
    paste(vars_excluidas_vif, collapse = ", ")
  ))
}

# ==============================================================================
# 3. MALLA + SPDE (compartida por las dos etapas)
# ==============================================================================
cmat <- as.matrix(unique(dt[, .(X_km, Y_km)]))
bnd_i <- inla.nonconvex.hull(cmat, convex = -0.05, resolution = 50)
bnd_o <- inla.nonconvex.hull(cmat, convex = -0.2)
malla <- inla.mesh.2d(
  loc = cmat, boundary = list(bnd_i, bnd_o),
  max.edge = c(4, 8), cutoff = 0.5
)
spde <- inla.spde2.matern(malla, alpha = 2)
A <- inla.spde.make.A(malla, loc = as.matrix(dt[, .(X_km, Y_km)]))
idx <- inla.spde.make.index("campo", n.spde = spde$n.spde)

ajustar_spde <- function(y, vars) {
  covs_list <- setNames(lapply(vars, function(v) dt[[v]]), vars)
  stk <- inla.stack(
    tag = "d", data = list(y = y), A = list(A, 1),
    effects = list(c(idx, list(intercept = 1)), covs_list)
  )
  fml <- as.formula(paste(
    "y ~ 0 + intercept +", paste(vars, collapse = " + "),
    "+ f(campo, model = spde)"
  ))
  inla(fml,
    data = inla.stack.data(stk, spde = spde), family = "gaussian",
    control.predictor = list(A = inla.stack.A(stk), compute = FALSE),
    control.compute = list(dic = TRUE, waic = TRUE),
    control.inla = list(int.strategy = "eb")
  )
}

es_signif <- function(fixed_row) {
  fixed_row[["0.025quant"]] > 0 | fixed_row[["0.975quant"]] < 0
}

# ==============================================================================
# 4. ETAPA 1 — modelo completo, identificar covariables significativas
# ==============================================================================
cat("\n[ETAPA 1] Ajustando modelo SPDE con las covariables candidatas...\n")
m1 <- ajustar_spde(dt$LOG_NO2, vars_vif)
fixed1 <- m1$summary.fixed
fixed1$sig <- es_signif(fixed1)

# Partición por TAMAÑO DE EFECTO (no por p-valor: con N enorme todo es sig.)
betas_cov <- setNames(fixed1[vars_vif, "mean"], vars_vif)
vars_sig <- vars_vif[abs(betas_cov) >= UMBRAL_BETA] # fuertes -> Etapa 1
vars_nosig <- vars_vif[abs(betas_cov) < UMBRAL_BETA] # débiles -> Etapa 2

cat(sprintf(
  "  Todas significativas por IC95%%: %s\n",
  paste(rownames(fixed1)[fixed1$sig & rownames(fixed1) != "intercept"],
    collapse = ", "
  )
))
cat(sprintf(
  "  FUERTES |β|>=%.2f (Etapa 1): %s\n", UMBRAL_BETA,
  ifelse(length(vars_sig) > 0, paste(vars_sig, collapse = ", "), "ninguna")
))
cat(sprintf(
  "  DÉBILES |β|<%.2f  (Etapa 2): %s\n", UMBRAL_BETA,
  ifelse(length(vars_nosig) > 0, paste(vars_nosig, collapse = ", "), "ninguna")
))

# --- Residuo parcial: y* = y − Σ β̂_sig · X_sig ---
dt[, y_estrella := LOG_NO2]
for (v in vars_sig) dt[, y_estrella := y_estrella - fixed1[v, "mean"] * get(v)]

sd_y <- sd(dt$LOG_NO2)
sd_ye <- sd(dt$y_estrella)
cat(sprintf(
  "  SD(y) = %.4f  ->  SD(y*) tras quitar sig. = %.4f (%.1f%% de la var. original)\n",
  sd_y, sd_ye, 100 * (sd_ye / sd_y)^2
))

# ==============================================================================
# 5. ETAPA 2 — regresar el residuo y* sobre las NO significativas (sin retardos)
# ==============================================================================
tabla_e2 <- NULL
sd_resid2 <- sd_ye
if (length(vars_nosig) >= 1) {
  cat("\n[ETAPA 2] Regresando el residuo y* sobre las no significativas...\n")

  # VIF de control sobre el subconjunto (debe seguir limpio)
  if (length(vars_nosig) >= 2) {
    vif2 <- car::vif(lm(as.formula(paste(
      "y_estrella ~",
      paste(vars_nosig, collapse = " + ")
    )), data = dt))
    cat(sprintf(
      "  VIF Etapa 2 (máx = %.2f): %s\n", max(vif2),
      paste(sprintf("%s=%.2f", names(vif2), vif2), collapse = ", ")
    ))
  }

  m2 <- ajustar_spde(dt$y_estrella, vars_nosig)
  fixed2 <- m2$summary.fixed
  fixed2$sig <- es_signif(fixed2)

  # Residuo tras quitar la contribución de las no-sig
  dt[, resid2 := y_estrella]
  for (v in vars_nosig) dt[, resid2 := resid2 - fixed2[v, "mean"] * get(v)]
  sd_resid2 <- sd(dt$resid2)

  vars_recuperadas <- setdiff(rownames(fixed2)[fixed2$sig], "intercept")
  cat(sprintf(
    "  Variables que TODAVÍA explican el residuo: %s\n",
    ifelse(length(vars_recuperadas) > 0,
      paste(vars_recuperadas, collapse = ", "),
      "ninguna (residuo limpio de covariables)"
    )
  ))
  cat(sprintf("  SD(y*) = %.4f  ->  SD(resid) tras Etapa 2 = %.4f\n", sd_ye, sd_resid2))

  tabla_e2 <- data.frame(
    Variable = rownames(fixed2),
    Media = round(fixed2[["mean"]], 4),
    IC025 = round(fixed2[["0.025quant"]], 4),
    IC975 = round(fixed2[["0.975quant"]], 4),
    Significativa = ifelse(fixed2$sig, "Sí", "no"),
    row.names = NULL
  )
  tabla_e2 <- tabla_e2[tabla_e2$Variable != "intercept", ]
} else {
  cat("\n[ETAPA 2] No hay covariables no significativas que probar.\n")
}

# ==============================================================================
# 5.5 ¿CUÁNTO DEL RESIDUO ES CAMPO ESPACIAL (SPDE) Y CUÁNTO RUIDO?
# ==============================================================================
# Responde a "¿cómo sé que el resto es el campo SPDE?": no se asume, se MIDE.
#   (a) Descomposición de varianza del modelo de la Etapa 1:
#         Var(y) ≈ Var(covariables) + σ²_espacial(SPDE) + σ²_ruido(nugget)
#   (b) Comparación WAIC: modelo CON campo SPDE vs SIN campo (solo covariables).
#   (c) Test de Moran sobre el residuo por estación: ¿hay autocorrelación
#       espacial real que el campo esté capturando (y no ruido)?
# ==============================================================================
cat("\n[5.5] ¿De qué está hecho el residuo? (campo espacial vs ruido)...\n")

spde_res1 <- inla.spde.result(m1, "campo", spde)
sigma2_campo <- exp(spde_res1$summary.log.variance.nominal$mean)
rango_campo <- exp(spde_res1$summary.log.range.nominal$mean)
fila_prec <- grep("^Precision for the Gaussian", rownames(m1$summary.hyperpar))
sigma2_nugget <- 1 / m1$summary.hyperpar$mean[fila_prec]
var_covs <- sd_y^2 - sd_ye^2 # varianza que quitan las covs fuertes
var_total <- var_covs + sigma2_campo + sigma2_nugget

tabla_var <- data.frame(
  Componente = c("Covariables fuertes (Etapa 1)", "Campo espacial SPDE", "Ruido (nugget)"),
  Varianza   = round(c(var_covs, sigma2_campo, sigma2_nugget), 4),
  Pct        = round(100 * c(var_covs, sigma2_campo, sigma2_nugget) / var_total, 1)
)
cat("  --- Descomposición de varianza (aprox.) ---\n")
print(tabla_var)
cat(sprintf("  Rango espacial estimado del campo SPDE: %.1f km\n", rango_campo))

# (b) WAIC: con campo vs sin campo
df0 <- as.data.frame(dt[, c("LOG_NO2", vars_sig), with = FALSE])
m0 <- inla(as.formula(paste("LOG_NO2 ~", paste(vars_sig, collapse = " + "))),
  data = df0, family = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE),
  control.inla = list(int.strategy = "eb")
)
dWAIC <- m0$waic$waic - m1$waic$waic
cat(sprintf(
  "  WAIC sin campo = %.0f | con campo SPDE = %.0f | mejora = %.0f\n",
  m0$waic$waic, m1$waic$waic, dWAIC
))

# (c) Moran sobre el residuo por estación (autocorrelación espacial real)
res_st <- dt[, .(res = mean(y_estrella), X_km = X_km[1], Y_km = Y_km[1]), by = ESTACION]
coords_st <- as.matrix(res_st[, .(X_km, Y_km)])
lw_st <- nb2listw(knn2nb(knearneigh(coords_st, k = 4)), style = "W")
moran <- moran.mc(res_st$res, lw_st, nsim = 999)
cat(sprintf(
  "  Moran I del residuo = %.3f (p = %.3f) -> %s\n",
  moran$statistic, moran$p.value,
  ifelse(moran$p.value < 0.05, "autocorrelacion espacial SIGNIFICATIVA",
    "sin autocorrelacion significativa"
  )
))

guardar_tabla_png(tabla_var,
  titulo = sprintf("¿De qué está hecho el residuo? — horario %s 2025", mes_lab),
  subtitulo = sprintf(
    "Rango SPDE ~%.0f km · Moran I=%.3f (p=%.3f) · WAIC mejora %.0f al añadir el campo",
    rango_campo, moran$statistic, moran$p.value, dWAIC
  ),
  ruta_png = file.path(carpeta_out, "04_descomposicion_varianza.png")
)
fwrite(tabla_var, file.path(carpeta_out, "04_descomposicion_varianza.csv"))

# ==============================================================================
# 5.6 ETAPA 3 — RETARDO DE TRÁFICO (h-1): ¿EXPLICA PARTE DEL RESIDUO?
# ==============================================================================
# Se añade el tráfico de la hora ANTERIOR (intensidad y carga en h-1) sobre el
# residuo y*. Como y* ya no contiene el tráfico de la MISMA hora, cualquier
# reducción aquí es efecto RETARDADO puro. Se compara, sobre el mismo
# subconjunto, con lo que aportaban las covariables débiles (Etapa 2).
# ==============================================================================
cat("\n[5.6] Etapa 3 — retardo de tráfico (h-1)...\n")

dt_lag <- copy(dt)
dt_lag[, datetime := as.POSIXct(FECHA, tz = "UTC") + (HORA - 1) * 3600]
setorder(dt_lag, ESTACION, datetime)
dt_lag[, intensidad_lag1 := shift(intensidad, 1), by = ESTACION]
dt_lag[, carga_lag1 := shift(carga, 1), by = ESTACION]
dt_lag[, dtprev := shift(datetime, 1), by = ESTACION]
dt_lag[, gap_h := as.numeric(difftime(datetime, dtprev, units = "hours"))]
dt_lag <- dt_lag[gap_h == 1 & !is.na(intensidad_lag1) & !is.na(carga_lag1)]

sd_ye_sub <- sd(dt_lag$y_estrella) # baseline sobre el subconjunto con retardo válido

A_lag <- inla.spde.make.A(malla, loc = as.matrix(dt_lag[, .(X_km, Y_km)]))
covs_lag <- list(
  intensidad_lag1 = dt_lag$intensidad_lag1,
  carga_lag1 = dt_lag$carga_lag1
)
stkL <- inla.stack(
  tag = "L", data = list(y = dt_lag$y_estrella), A = list(A_lag, 1),
  effects = list(c(idx, list(intercept = 1)), covs_lag)
)
mL <- inla(y ~ 0 + intercept + intensidad_lag1 + carga_lag1 + f(campo, model = spde),
  data = inla.stack.data(stkL, spde = spde), family = "gaussian",
  control.predictor = list(A = inla.stack.A(stkL), compute = FALSE),
  control.compute = list(waic = TRUE),
  control.inla = list(int.strategy = "eb")
)
fixedL <- mL$summary.fixed
fixedL$sig <- es_signif(fixedL)

dt_lag[, res_lag := y_estrella
  - fixedL["intensidad_lag1", "mean"] * intensidad_lag1
  - fixedL["carga_lag1", "mean"] * carga_lag1]
sd_res_lag <- sd(dt_lag$res_lag)

# Residuo de la Etapa 2 (covariables débiles) sobre el MISMO subconjunto
if (!is.null(tabla_e2)) {
  dt_lag[, res_e2 := y_estrella]
  for (v in vars_nosig) dt_lag[, res_e2 := res_e2 - fixed2[v, "mean"] * get(v)]
  sd_e2_sub <- sd(dt_lag$res_e2)
} else {
  sd_e2_sub <- sd_ye_sub
}

cat(sprintf("  Observaciones con retardo válido: %d\n", nrow(dt_lag)))
cat(sprintf(
  "  β intensidad(h-1) = %.4f [%.4f, %.4f] %s\n",
  fixedL["intensidad_lag1", "mean"], fixedL["intensidad_lag1", "0.025quant"],
  fixedL["intensidad_lag1", "0.975quant"],
  ifelse(fixedL["intensidad_lag1", "sig"], "sig", "no")
))
cat(sprintf(
  "  β carga(h-1)      = %.4f [%.4f, %.4f] %s\n",
  fixedL["carga_lag1", "mean"], fixedL["carga_lag1", "0.025quant"],
  fixedL["carga_lag1", "0.975quant"],
  ifelse(fixedL["carga_lag1", "sig"], "sig", "no")
))
cat(sprintf(
  "  SD(y*)=%.4f | tras covs débiles=%.4f | tras retardo h-1=%.4f\n",
  sd_ye_sub, sd_e2_sub, sd_res_lag
))

tabla_lag <- data.frame(
  Variable      = c("intensidad (h-1)", "carga (h-1)"),
  Beta          = round(c(fixedL["intensidad_lag1", "mean"], fixedL["carga_lag1", "mean"]), 4),
  IC025         = round(c(fixedL["intensidad_lag1", "0.025quant"], fixedL["carga_lag1", "0.025quant"]), 4),
  IC975         = round(c(fixedL["intensidad_lag1", "0.975quant"], fixedL["carga_lag1", "0.975quant"]), 4),
  Significativa = ifelse(c(fixedL["intensidad_lag1", "sig"], fixedL["carga_lag1", "sig"]), "Sí", "no")
)
guardar_tabla_png(tabla_lag,
  titulo = "Etapa 3 — efecto retardado del tráfico (h-1) sobre el residuo",
  subtitulo = sprintf(
    "SD(y*)=%.4f -> tras retardo=%.4f (reduce %.2f%%) · N=%d",
    sd_ye_sub, sd_res_lag, 100 * (1 - sd_res_lag / sd_ye_sub), nrow(dt_lag)
  ),
  ruta_png = file.path(carpeta_out, "05_retardo_trafico.png")
)
fwrite(tabla_lag, file.path(carpeta_out, "05_retardo_trafico.csv"))

# Comparación justa (mismo subconjunto): débiles (Et.2) vs retardo (Et.3)
df_cmp <- data.frame(
  Paso = factor(c("y* (baseline)", "+ covs débiles (Et.2)", "+ tráfico h-1 (Et.3)"),
    levels = c("y* (baseline)", "+ covs débiles (Et.2)", "+ tráfico h-1 (Et.3)")
  ),
  SD = c(sd_ye_sub, sd_e2_sub, sd_res_lag)
)
ggplot(df_cmp, aes(x = Paso, y = SD, fill = Paso)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("%.4f", SD)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("#7fb3d5", "#f39c12", "#c0392b"), guide = "none") +
  coord_cartesian(ylim = c(0, max(df_cmp$SD) * 1.08)) +
  labs(
    title = sprintf("¿Qué explica más el residuo? débiles vs retardo — horario %s 2025", mes_lab),
    subtitle = "SD del residuo (mismo subconjunto): covariables descartadas vs tráfico de la hora anterior",
    x = NULL, y = "Desviación típica del residuo y*"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 9)
  )
ggsave(file.path(carpeta_out, "06_comparacion_retardo_vs_debiles.png"),
  width = 9, height = 6, dpi = 300, bg = "white"
)

# ==============================================================================
# 5.7 ¿EL RUIDO RESTANTE ES BLANCO O TIENE ESTRUCTURA TEMPORAL (AR1)?
# ==============================================================================
# ACF + Ljung-Box del residuo y* DENTRO de cada estación (su serie horaria).
# El campo espacial de la Etapa 1 es constante por estación, así que no afecta
# a la ACF intra-estación: lo que se ve aquí es la memoria temporal del residuo.
#   - ACF de lag-1 alta / Ljung-Box significativo -> queda estructura temporal
#     aprovechable (justifica un término AR1).
#   - ACF ~ 0 -> ruido blanco: el nugget es irreducible, el AR1 no ayudaría.
# ==============================================================================
cat("\n[5.7] ¿El residuo tiene estructura temporal (AR1) o es ruido blanco?...\n")

LAG_MAX_ACF <- 24L
dt[, datetime := as.POSIXct(FECHA, tz = "UTC") + (HORA - 1) * 3600]
setorder(dt, ESTACION, datetime)

series_est <- split(dt$y_estrella, dt$ESTACION)
acf_list <- lapply(series_est, function(s) {
  if (length(s) > LAG_MAX_ACF + 5) {
    as.numeric(acf(s, lag.max = LAG_MAX_ACF, plot = FALSE)$acf)[-1]
  } else {
    NULL
  }
})
acf_mat <- do.call(rbind, acf_list[!vapply(acf_list, is.null, logical(1))])
acf_mean <- colMeans(acf_mat)

lb_p <- vapply(
  series_est, function(s) {
    if (length(s) > LAG_MAX_ACF + 5) {
      Box.test(s, lag = LAG_MAX_ACF, type = "Ljung-Box")$p.value
    } else {
      NA_real_
    }
  },
  numeric(1)
)
n_est_sig <- sum(lb_p < 0.05, na.rm = TRUE)
n_est_tot <- sum(!is.na(lb_p))
n_med <- median(lengths(series_est))
banda <- 1.96 / sqrt(n_med)

cat(sprintf(
  "  ACF lag-1 (media entre estaciones) = %.3f | banda ±%.3f\n",
  acf_mean[1], banda
))
cat(sprintf(
  "  Ljung-Box (lag %d): %d/%d estaciones rechazan ruido blanco (p<0.05)\n",
  LAG_MAX_ACF, n_est_sig, n_est_tot
))

acf_df <- data.frame(lag = seq_len(LAG_MAX_ACF), acf = acf_mean)
ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_hline(yintercept = c(-banda, banda), linetype = "dashed", color = "#c0392b") +
  geom_segment(aes(xend = lag, yend = 0), color = "#1a5276", linewidth = 1) +
  geom_point(color = "#1a5276", size = 2) +
  scale_x_continuous(breaks = seq(0, LAG_MAX_ACF, 2)) +
  labs(
    title = sprintf("ACF temporal del residuo (intra-estación) — horario %s 2025", mes_lab),
    subtitle = sprintf(
      "Media entre estaciones · banda ±%.3f · %d/%d estaciones rechazan ruido blanco (Ljung-Box lag %d)",
      banda, n_est_sig, n_est_tot, LAG_MAX_ACF
    ),
    x = "Retardo (horas)", y = "Autocorrelación"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 9)
  )
ggsave(file.path(carpeta_out, "07_acf_temporal_residuo.png"),
  width = 10, height = 6, dpi = 300, bg = "white"
)

# ==============================================================================
# 6. SALIDAS
# ==============================================================================

# --- Tabla 1: estado de cada covariable candidata ---
estado <- data.frame(
  Variable = covs_disp,
  Estado = ifelse(covs_disp %in% vars_excluidas_vif, "Excluida (colinealidad)",
    ifelse(covs_disp %in% vars_sig, "Fuerte (Etapa 1)",
      "Débil (Etapa 2)"
    )
  ),
  Beta_Etapa1 = sapply(covs_disp, function(v) {
    if (v %in% rownames(fixed1)) round(fixed1[v, "mean"], 4) else NA
  }),
  row.names = NULL
)
estado <- estado[order(-abs(estado$Beta_Etapa1)), ]
guardar_tabla_png(estado,
  titulo = sprintf("Estado de las covariables — horario %s 2025", mes_lab),
  subtitulo = sprintf(
    "Partición por tamaño de efecto |β| >= %.2f · N grande => todo es sig. por IC95%%",
    UMBRAL_BETA
  ),
  ruta_png = file.path(carpeta_out, "01_estado_covariables.png")
)
fwrite(estado, file.path(carpeta_out, "01_estado_covariables.csv"))

# --- Tabla 2: resultados de la Etapa 2 ---
if (!is.null(tabla_e2)) {
  guardar_tabla_png(tabla_e2,
    titulo    = "Etapa 2 — ¿qué queda en el residuo?",
    subtitulo = "Efecto de las covariables NO significativas sobre y* (residuo parcial)",
    ruta_png  = file.path(carpeta_out, "02_etapa2_residuo.png")
  )
  fwrite(tabla_e2, file.path(carpeta_out, "02_etapa2_residuo.csv"))
}

# --- Gráfico: reducción de la desviación típica por etapa ---
df_sd <- data.frame(
  Fase = factor(c("y (original)", "y* (tras Etapa 1)", "resid (tras Etapa 2)"),
    levels = c("y (original)", "y* (tras Etapa 1)", "resid (tras Etapa 2)")
  ),
  SD = c(sd_y, sd_ye, sd_resid2)
)
ggplot(df_sd, aes(x = Fase, y = SD, fill = Fase)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("%.4f", SD)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("#1a5276", "#7fb3d5", "#c0392b"), guide = "none") +
  labs(
    title = sprintf("Reducción de variabilidad por etapa — horario %s 2025", mes_lab),
    subtitle = "Cuánto baja la SD del residuo al quitar covariables significativas y, después, las descartadas",
    x = NULL, y = "Desviación típica de log(NO₂)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 9)
  )
ggsave(file.path(carpeta_out, "03_reduccion_sd_por_etapa.png"),
  width = 9, height = 6, dpi = 300, bg = "white"
)

cat(sprintf("\n✅ Diagnóstico completado. Salidas en:\n   %s\n", carpeta_out))
