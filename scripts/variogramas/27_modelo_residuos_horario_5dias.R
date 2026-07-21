# ==============================================================================
# MODELO DE LOS RESIDUOS HORARIOS (muestra de 5 dias) + ESTRUCTURA
# ------------------------------------------------------------------------------
# Se forma el residuo y* = log(NO2) - contribucion de las covariables FUERTES,
# y se modela y* probando que conviene anadir ADEMAS de las covariables debiles:
#   A) Solo covariables (debiles + retardo h-1 del trafico)   [sin estructura]
#   B) A + AR(1) temporal (por estacion)
#   C) A + AR(2) temporal (por estacion)
#   D) A + campo espacial SPDE
#   E) A + espacio-temporal (SPDE (x) AR1)
# Comparacion en test (dia 5) con DIC/WAIC/RMSE/MAE/Cobertura95 + ACF/QQ.
# Train = dias 1-4 (96 h) | Test = dia 5 (24 h). Malla media (~4 km).
# Salidas: outputs/simulacion/residuos_horario_5dias/
# ==============================================================================

library(INLA)
library(data.table)
library(sf)
library(car)
library(ggplot2)
library(gridExtra); library(grid)
library(here)

set.seed(4827)

FECHA_INI <- as.Date("2025-01-01"); FECHA_FIN <- as.Date("2025-01-05")
DIAS_TRAIN <- 4L                     # dias 1-4 train, dia 5 test
UMBRAL_BETA <- 0.10                  # |beta std| para covariable "fuerte"
carpeta <- here("outputs","simulacion","residuos_horario_5dias")
dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)

COVS_RAW <- c("intensidad_raw","carga_raw","Temperatura_raw","Humedad_Relativa_raw",
              "Precipitaciones_raw","Presion Barométrica_raw","Radiación Solar_raw",
              "Velocidad Viento_raw")
COVS_AL  <- c("intensidad","carga","temperatura","humedad",
              "precipitacion","presion","radiacion","viento")

# ------------------------------------------------------------------------------
# 1. DATOS (5 dias) + estandarizar + coordenadas km + indices
# ------------------------------------------------------------------------------
dt <- as.data.table(readRDS(here("data","processed","Maestro","horario",
        "dataset_maestro_inla_2025_HORARIO.rds")))
setnames(dt, "LOG_NO2_HORARIO", "LOG_NO2")
dt <- dt[as.Date(FECHA) >= FECHA_INI & as.Date(FECHA) <= FECHA_FIN & !is.na(LOG_NO2)]

covs <- character(0)
for (i in seq_along(COVS_RAW)) if (COVS_RAW[i] %in% names(dt)) {
  dt[, (COVS_AL[i]) := as.numeric(scale(get(COVS_RAW[i]))) ]; covs <- c(covs, COVS_AL[i]) }

coords_u <- unique(dt[, .(ESTACION, LONGITUD, LATITUD)])
csf <- st_as_sf(coords_u, coords = c("LONGITUD","LATITUD"), crs = 4326) |> st_transform(25830)
coords_u[, X_km := st_coordinates(csf)[,1]/1000][, Y_km := st_coordinates(csf)[,2]/1000]
dt <- merge(dt, coords_u[, .(ESTACION, X_km, Y_km)], by = "ESTACION")

dt[, datetime := as.POSIXct(FECHA, tz="UTC") + (HORA-1)*3600]
setorder(dt, ESTACION, datetime)
dt[, intensidad_lag1 := shift(intensidad, 1L), by = ESTACION]
dt <- dt[complete.cases(dt[, c("LOG_NO2", covs, "intensidad_lag1"), with=FALSE])]

# indices temporales (1..T) y de estacion
inst <- sort(unique(dt$datetime)); dt[, t_idx := match(datetime, inst)]
dt[, est_idx := as.integer(factor(ESTACION))]
dt[, es_test := t_idx > DIAS_TRAIN*24]        # dia 5 = test

# ------------------------------------------------------------------------------
# 2. RESIDUO y* = log(NO2) - contribucion de covariables FUERTES
# ------------------------------------------------------------------------------
m0 <- lm(reformulate(covs, "LOG_NO2"), data = dt)
b0 <- coef(m0)[covs]
fuertes <- covs[abs(b0) >= UMBRAL_BETA]
debiles <- setdiff(covs, fuertes)
dt[, ystar := LOG_NO2]
for (v in fuertes) dt[, ystar := ystar - b0[[v]]*get(v)]
cat(sprintf("Fuertes (restadas): %s\nDebiles (van al modelo): %s\n",
            paste(fuertes, collapse=", "), paste(debiles, collapse=", ")))
cat(sprintf("SD(log NO2)=%.4f -> SD(y*)=%.4f | obs=%d (train=%d, test=%d)\n",
            sd(dt$LOG_NO2), sd(dt$ystar), nrow(dt), sum(!dt$es_test), sum(dt$es_test)))

vars_fix <- c(debiles, "intensidad_lag1")     # covariables del modelo de residuos

# ------------------------------------------------------------------------------
# 3. MALLA + SPDE (media ~4 km)
# ------------------------------------------------------------------------------
cmat <- as.matrix(unique(dt[, .(X_km, Y_km)]))
malla <- inla.mesh.2d(loc = cmat,
  boundary = list(inla.nonconvex.hull(cmat, convex=-0.05),
                  inla.nonconvex.hull(cmat, convex=-0.2)),
  max.edge = c(4,8), cutoff = 0.5)
spde <- inla.spde2.matern(malla, alpha = 2)

# respuesta con NA en test (para predecir)
y_cv <- ifelse(dt$es_test, NA, dt$ystar)
idx_test <- which(dt$es_test); y_true <- dt$ystar[idx_test]

metricas <- function(pred_m, pred_s, prec) {
  sdp <- sqrt(pred_s^2 + 1/prec); r <- y_true - pred_m
  c(RMSE=sqrt(mean(r^2)), MAE=mean(abs(r)),
    Cov95=100*mean(y_true>=pred_m-1.96*sdp & y_true<=pred_m+1.96*sdp))
}
fila <- function(nombre, mod, pred_m, pred_s) {
  prec <- mod$summary.hyperpar["Precision for the Gaussian observations","mean"]
  mm <- metricas(pred_m, pred_s, prec)
  data.table(Modelo=nombre, DIC=round(mod$dic$dic,1), WAIC=round(mod$waic$waic,1),
             RMSE=round(mm["RMSE"],4), MAE=round(mm["MAE"],4), Cov95=round(mm["Cov95"],1))
}
CC <- list(dic=TRUE, waic=TRUE)

# ------------------------------------------------------------------------------
# 4. MODELOS NO ESPACIALES (data.frame directo)
# ------------------------------------------------------------------------------
df <- copy(dt); df[, y := y_cv]
res <- list()

cat("\n[A] Solo covariables...\n")
fA <- as.formula(paste("y ~", paste(vars_fix, collapse=" + ")))
mA <- inla(fA, family="gaussian", data=df,
           control.predictor=list(compute=TRUE), control.compute=CC)
res[["A"]] <- fila("A. Covariables", mA, mA$summary.fitted.values$mean[idx_test],
                   mA$summary.fitted.values$sd[idx_test])

cat("[B] + AR(1) temporal...\n")
fB <- as.formula(paste("y ~", paste(vars_fix, collapse=" + "),
                       "+ f(t_idx, model='ar1', replicate=est_idx)"))
mB <- inla(fB, family="gaussian", data=df,
           control.predictor=list(compute=TRUE), control.compute=CC)
res[["B"]] <- fila("B. + AR(1) temporal", mB, mB$summary.fitted.values$mean[idx_test],
                   mB$summary.fitted.values$sd[idx_test])

cat("[C] + AR(2) temporal...\n")
fC <- as.formula(paste("y ~", paste(vars_fix, collapse=" + "),
                       "+ f(t_idx, model='ar', order=2, replicate=est_idx)"))
mC <- inla(fC, family="gaussian", data=df,
           control.predictor=list(compute=TRUE), control.compute=CC)
res[["C"]] <- fila("C. + AR(2) temporal", mC, mC$summary.fitted.values$mean[idx_test],
                   mC$summary.fitted.values$sd[idx_test])

# ------------------------------------------------------------------------------
# 5. MODELOS ESPACIALES (stack)
# ------------------------------------------------------------------------------
A_sp <- inla.spde.make.A(malla, loc = as.matrix(dt[, .(X_km, Y_km)]))
idx_s <- inla.spde.make.index("campo", n.spde = spde$n.spde)
covs_list <- c(setNames(lapply(vars_fix, function(v) dt[[v]]), vars_fix),
               list(intercept = rep(1, nrow(dt))))

cat("[D] + espacial SPDE...\n")
stkD <- inla.stack(tag="e", data=list(y=y_cv), A=list(A_sp,1),
                   effects=list(idx_s, covs_list))
fD <- as.formula(paste("y ~ 0 + intercept +", paste(vars_fix, collapse=" + "),
                       "+ f(campo, model=spde)"))
mD <- inla(fD, family="gaussian", data=inla.stack.data(stkD, spde=spde),
           control.predictor=list(A=inla.stack.A(stkD), compute=TRUE), control.compute=CC)
iD <- inla.stack.index(stkD,"e")$data[idx_test]
res[["D"]] <- fila("D. + espacial SPDE", mD, mD$summary.fitted.values$mean[iD],
                   mD$summary.fitted.values$sd[iD])

cat("[E] espacio-temporal SPDE (x) AR1...\n")
idx_st <- inla.spde.make.index("campo", n.spde=spde$n.spde, n.group=max(dt$t_idx))
A_st <- inla.spde.make.A(malla, loc=as.matrix(dt[,.(X_km,Y_km)]),
                         group=dt$t_idx, n.group=max(dt$t_idx))
stkE <- inla.stack(tag="e", data=list(y=y_cv), A=list(A_st,1),
                   effects=list(idx_st, covs_list))
fE <- as.formula(paste("y ~ 0 + intercept +", paste(vars_fix, collapse=" + "),
     "+ f(campo, model=spde, group=campo.group, control.group=list(model='ar1'))"))
mE <- inla(fE, family="gaussian", data=inla.stack.data(stkE, spde=spde),
           control.predictor=list(A=inla.stack.A(stkE), compute=TRUE), control.compute=CC)
iE <- inla.stack.index(stkE,"e")$data[idx_test]
res[["E"]] <- fila("E. SPDE (x) AR1", mE, mE$summary.fitted.values$mean[iE],
                   mE$summary.fitted.values$sd[iE])

# ------------------------------------------------------------------------------
# 6. TABLA COMPARATIVA
# ------------------------------------------------------------------------------
tabla <- rbindlist(res)
print(tabla)
fwrite(tabla, file.path(carpeta, "comparacion_modelos_residuos.csv"))

th <- ttheme_minimal(core=list(fg_params=list(fontsize=10),
        bg_params=list(fill=c("white","#F5F5F5"))),
        colhead=list(fg_params=list(fontsize=10,fontface="bold"),
        bg_params=list(fill="#DDEEFF")))
g <- arrangeGrob(
  textGrob("Modelo de los residuos horarios (5 dias) - ¿que anadir a las covariables?",
           gp=gpar(fontsize=13,fontface="bold")),
  textGrob("Respuesta y* = log(NO2) - covariables fuertes | Train dias 1-4, Test dia 5 | menor DIC/WAIC/RMSE y Cov95~95 mejor",
           gp=gpar(fontsize=9,col="grey40")),
  tableGrob(tabla, rows=NULL, theme=th),
  heights=unit(c(0.5,0.35,nrow(tabla)*0.33+0.5),"inches"))
ggsave(file.path(carpeta,"tabla_comparacion_residuos.png"), g,
       width=10, height=nrow(tabla)*0.4+1.6, dpi=150, bg="white")

# ------------------------------------------------------------------------------
# 7. ACF + QQ del mejor por WAIC
# ------------------------------------------------------------------------------
mejor <- tabla$Modelo[which.min(tabla$WAIC)]
cat(sprintf("\nMejor por WAIC: %s\n", mejor))
mods <- list("A. Covariables"=mA,"B. + AR(1) temporal"=mB,"C. + AR(2) temporal"=mC,
             "D. + espacial SPDE"=list(m=mD,i=inla.stack.index(stkD,"e")$data),
             "E. SPDE (x) AR1"=list(m=mE,i=inla.stack.index(stkE,"e")$data))
obj <- mods[[mejor]]
if (is.list(obj) && !is.null(obj$m)) { fit <- obj$m$summary.fitted.values$mean[obj$i]
} else fit <- obj$summary.fitted.values$mean[seq_len(nrow(dt))]
dt[, res_mejor := ystar - fit]
res_t <- dt[!(es_test), .(r=mean(res_mejor)), by=t_idx][order(t_idx)]
png(file.path(carpeta,"acf_residuos_mejor.png"), width=800, height=500)
acf(res_t$r, lag.max=24, main=sprintf("ACF residuos del mejor modelo: %s", mejor),
    col="#2166AC", lwd=2); dev.off()
png(file.path(carpeta,"qq_residuos_mejor.png"), width=600, height=600)
qqnorm(dt$res_mejor, main=sprintf("Q-Q residuos: %s", mejor)); qqline(dt$res_mejor, col="red", lwd=2)
dev.off()

cat(sprintf("\nListo. Salidas en %s\n", carpeta))
