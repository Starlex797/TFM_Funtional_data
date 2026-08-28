# ==============================================================================
# Kriging Funcional Ordinario (OFK) del NO2 en Madrid
# ==============================================================================
# Replica la estrategia de Montero-Lorenzo, Fernandez-Aviles, Mondejar-Jimenez y
# Vargas-Vargas (2013), "A spatio-temporal geostatistical approach to predicting
# pollution levels: The case of mono-nitrogen oxides in Madrid", Sec. 3.2 y 4.3.2.
#
# Idea: cada estacion se convierte en una CURVA (su perfil temporal de log-NO2) y
# se predice la curva completa en cualquier punto sin estacion como una MEDIA
# PONDERADA de las curvas vecinas. Los pesos (uno por estacion, constante en todo
# t) salen de un kriging ordinario sobre el TRACE-VARIOGRAMA, que integra en el
# tiempo la distancia L2 entre curvas (Ec. 16). No usa covariables ni INLA-SPDE.
#
# Pasos (con las ecuaciones del paper):
#   1. Dato funcional: curva "dia tipo laborable" = 24h x 52 semanas concatenadas.
#   2. Transformacion Ln(NO2) (Sec. 4.2).
#   3. Suavizado B-spline, lambda por GCV (Sec. 4.3.2).
#   4. Trace-variograma empirico (Ec. 16).
#   5. Ajuste de variograma teorico (esferico, como el paper).
#   6. Sistema OFK + prediccion de la curva en la rejilla (Ec. 10, 13-14).
#   7. Trace-variance = incertidumbre global (Ec. 15).
#   8. Validacion LOOCV espacial (RMSE integrado), desglosada por NOM_TIPO.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(sf)
  library(fda)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})
sf_use_s2(FALSE)

ANIO         <- 2025L
NBASIS       <- 252L                      # ~ nudos del paper (L*=248 -> ~252 bases)
LAMBDAS_GCV  <- 10^seq(-2, 4)             # rejilla para elegir la penalizacion
HORA_PUNTA   <- 8L                        # H08, punta de manana
DIR_SALIDA   <- here("outputs", "figures", "kriging_funcional_no2")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

# Festivos 2025 en la ciudad de Madrid (nacionales + CCAA + locales). Los que
# caen en fin de semana son inocuos (ya se descartan los findes).
FESTIVOS_2025 <- as.Date(c(
  "2025-01-01", "2025-01-06", "2025-04-17", "2025-04-18", "2025-05-01",
  "2025-05-02", "2025-05-15", "2025-07-25", "2025-08-15", "2025-10-12",
  "2025-11-01", "2025-11-09", "2025-12-06", "2025-12-08", "2025-12-25"
))

# ------------------------------------------------------------------------------
# 1. Dato funcional: curva "dia tipo laborable" (Sec. 4.2)
# ------------------------------------------------------------------------------

dt <- as.data.table(readRDS(here(
  "data", "processed", "Contaminacion", "horario",
  sprintf("aire_madrid_%d_No2_horarios.rds", ANIO)
)))

dt[, hora := as.integer(sub("H", "", HORA))]
dt[, wday := as.integer(format(FECHA, "%u"))]                  # 1=Lun ... 7=Dom
dt[, laborable := wday <= 5 & !(FECHA %in% FESTIVOS_2025)]

# Media horaria de los dias laborables de cada semana -> "dia tipo" semanal.
lab <- dt[laborable == TRUE & !is.na(DATO)]
lab[, sem_date := as.Date(cut(FECHA, breaks = "week"))]        # semanas Lun-Dom
agg <- lab[, .(m = mean(DATO)),
           by = .(ESTACION, NOM_TIPO, LONGITUD, LATITUD, sem_date, hora)]

# Indice temporal t = 1..N concatenando las 24 h de cada semana en orden.
semanas <- sort(unique(agg$sem_date))
W <- length(semanas)
N <- W * 24L
agg[, sem := match(sem_date, semanas)]
agg[, tidx := (sem - 1L) * 24L + hora]
agg[, vlog := log(pmax(m, 0.1))]                              # Ln(NO2), Sec. 4.2

# Metadatos de estacion y matriz de curvas (N x n) con NA en huecos.
sts <- unique(agg[, .(ESTACION, NOM_TIPO, LONGITUD, LATITUD)])
setorder(sts, ESTACION)
n <- nrow(sts)
Y <- matrix(NA_real_, nrow = N, ncol = n,
            dimnames = list(NULL, sts$ESTACION))
for (k in seq_len(n)) {
  d <- agg[ESTACION == sts$ESTACION[k]]
  Y[d$tidx, k] <- d$vlog
}
cat(sprintf(
  "Estaciones: %d | Semanas: %d | Longitud curva N: %d | %% celdas con dato: %.1f%%\n",
  n, W, N, 100 * mean(!is.na(Y))
))

# ------------------------------------------------------------------------------
# 2-3. Suavizado B-spline con lambda elegido por GCV (Sec. 4.3.2)
# ------------------------------------------------------------------------------

basis <- create.bspline.basis(rangeval = c(1, N), nbasis = NBASIS, norder = 4)

# GCV medio sobre todas las estaciones para cada lambda.
suavizar_estacion <- function(k, fdpar) {
  idx <- which(!is.na(Y[, k]))
  smooth.basis(idx, Y[idx, k], fdpar)
}
gcv_por_lambda <- sapply(LAMBDAS_GCV, function(lmb) {
  fdpar <- fdPar(basis, Lfdobj = 2, lambda = lmb)
  g <- sapply(seq_len(n), function(k) suavizar_estacion(k, fdpar)$gcv)
  mean(g, na.rm = TRUE)
})
lambda_opt <- LAMBDAS_GCV[which.min(gcv_por_lambda)]
cat(sprintf("Lambda de suavizado optimo (GCV): %g\n", lambda_opt))

# Suavizado final: coeficientes en la base comun -> objeto fd -> curvas en malla.
fdpar_opt <- fdPar(basis, Lfdobj = 2, lambda = lambda_opt)
coefs <- matrix(0, nrow = NBASIS, ncol = n)
for (k in seq_len(n)) coefs[, k] <- suavizar_estacion(k, fdpar_opt)$fd$coefs
Xfd <- fd(coefs, basis)
tt <- 1:N
Xmat <- eval.fd(tt, Xfd)                 # N x n : curvas suavizadas de log-NO2
colnames(Xmat) <- sts$ESTACION
hora_de_t <- ((tt - 1L) %% 24L) + 1L     # hora del dia para cada t

# ------------------------------------------------------------------------------
# 4. Trace-variograma empirico (Ec. 16)
# ------------------------------------------------------------------------------

# Coordenadas metricas (UTM 30N ETRS89) para distancias en metros.
sf_sts <- st_transform(
  st_as_sf(sts, coords = c("LONGITUD", "LATITUD"), crs = 4326), 25830
)
XY <- st_coordinates(sf_sts)
D  <- as.matrix(dist(XY))                # n x n : distancias entre estaciones

# Distancia L2 al cuadrado entre cada par de curvas (la integral en t).
ISD <- matrix(0, n, n)
for (i in seq_len(n - 1)) for (j in (i + 1):n) {
  ISD[i, j] <- ISD[j, i] <- sum((Xmat[, i] - Xmat[, j])^2)
}

# Nube por clases de distancia -> gamma(h) = 0.5 * media de la ISD en la clase.
pares <- which(upper.tri(D), arr.ind = TRUE)
dpar  <- D[upper.tri(D)]
ipar  <- ISD[upper.tri(ISD)]
cutoff <- 0.6 * max(dpar)
brk <- seq(0, cutoff, length.out = 13)
cls <- cut(dpar, brk, include.lowest = TRUE)
vgram <- data.table(dist = dpar, isd = ipar, cls = cls)[!is.na(cls),
  .(h = mean(dist), gamma = 0.5 * mean(isd), np = .N), by = cls][order(h)]
cat(sprintf("Trace-variograma: %d clases hasta %.0f m\n", nrow(vgram), cutoff))

# ------------------------------------------------------------------------------
# 5. Ajuste del variograma teorico (esferico, como el paper; se comparan 3)
# ------------------------------------------------------------------------------

m_sph <- function(h, p) p[1] + p[2] * ifelse(
  h < p[3], 1.5 * (h / p[3]) - 0.5 * (h / p[3])^3, 1)
m_exp <- function(h, p) p[1] + p[2] * (1 - exp(-h / p[3]))
m_gau <- function(h, p) p[1] + p[2] * (1 - exp(-(h / p[3])^2))

ajustar <- function(modelo) {
  wss <- function(p) sum(vgram$np * (vgram$gamma - modelo(vgram$h, p))^2)
  ini <- c(min(vgram$gamma), diff(range(vgram$gamma)), 0.4 * cutoff)
  optim(ini, wss, method = "L-BFGS-B",
        lower = c(0, 1e-9, 1e-6),
        upper = c(max(vgram$gamma), 5 * max(vgram$gamma), 3 * cutoff))
}
fits <- list(Esferico = ajustar(m_sph),
             Exponencial = ajustar(m_exp),
             Gaussiano = ajustar(m_gau))
for (nm in names(fits)) cat(sprintf(
  "  %-12s nugget=%.3f  sill=%.3f  rango=%.0f m  WSS=%.3g\n",
  nm, fits[[nm]]$par[1], fits[[nm]]$par[2], fits[[nm]]$par[3], fits[[nm]]$value))

# Se usa el modelo esferico para predecir (fiel al paper).
p_sph <- fits$Esferico$par
gamma_mod <- function(h) m_sph(h, p_sph)

# ------------------------------------------------------------------------------
# 6-7. Sistema OFK: prediccion en la rejilla + trace-variance (Ec. 10, 13-15)
# ------------------------------------------------------------------------------

# Matriz del kriging ordinario (aumentada con la restriccion de insesgadez).
G <- gamma_mod(D); diag(G) <- 0
A <- rbind(cbind(G, 1), c(rep(1, n), 0))
Ainv <- solve(A)

# Rejilla sobre Madrid.
mapa_utm <- st_transform(st_make_valid(st_read(
  here("data", "raw", "geometrias", "madrid_distritos.geojson"), quiet = TRUE
)), 25830)
rej <- st_as_sf(st_make_grid(mapa_utm, n = c(90, 90), what = "centers"))
rej <- st_intersection(rej, st_union(mapa_utm))
GR  <- st_coordinates(rej)
M   <- nrow(GR)
dx  <- diff(sort(unique(GR[, 1])))[1]
dy  <- diff(sort(unique(GR[, 2])))[1]

# Distancias rejilla-estacion y RHS del sistema para todos los nodos a la vez.
D0 <- matrix(0, n, M)
for (k in seq_len(n)) D0[k, ] <- sqrt((GR[, 1] - XY[k, 1])^2 + (GR[, 2] - XY[k, 2])^2)
G0  <- gamma_mod(D0)                                   # n x M
RHS <- rbind(G0, 1)                                    # (n+1) x M
SOL <- Ainv %*% RHS                                    # (n+1) x M
Lam <- SOL[1:n, , drop = FALSE]                        # pesos por nodo
mu  <- SOL[n + 1, ]                                    # multiplicador de Lagrange

# Superficies: cualquier resumen de la curva = combinacion de los pesos.
resumen_estacion <- function(rows) colMeans(Xmat[rows, , drop = FALSE])
m_anual <- resumen_estacion(seq_len(N))                # media anual de log-NO2
m_punta <- resumen_estacion(which(hora_de_t == HORA_PUNTA))
surf_anual <- exp(as.vector(m_anual %*% Lam))          # back-transf. (media geom.)
surf_punta <- exp(as.vector(m_punta %*% Lam))
trace_var  <- colSums(G0 * Lam) + mu                   # Ec. 15 (unidades log^2)

# ------------------------------------------------------------------------------
# 8. Validacion LOOCV espacial (RMSE integrado), desglose por NOM_TIPO
# ------------------------------------------------------------------------------

loocv_rmse <- numeric(n)
for (i in seq_len(n)) {
  j <- setdiff(seq_len(n), i)
  Gi <- gamma_mod(D[j, j]); diag(Gi) <- 0
  Ai <- rbind(cbind(Gi, 1), c(rep(1, n - 1), 0))
  g0 <- c(gamma_mod(D[j, i]), 1)
  sol <- solve(Ai, g0)
  lam <- sol[1:(n - 1)]
  pred <- Xmat[, j, drop = FALSE] %*% lam                # curva predicha
  loocv_rmse[i] <- sqrt(mean((pred - Xmat[, i])^2))      # RMSE integrado (log)
}
res_loocv <- data.table(
  ESTACION = sts$ESTACION, NOM_TIPO = sts$NOM_TIPO, RMSE_log = loocv_rmse
)[order(RMSE_log)]

cat("\nLOOCV RMSE (log-NO2) global:", sprintf("%.4f", mean(loocv_rmse)), "\n")
cat("LOOCV RMSE por tipologia:\n")
print(res_loocv[, .(RMSE_medio = round(mean(RMSE_log), 4), n = .N), by = NOM_TIPO])
fwrite(res_loocv, file.path(DIR_SALIDA, "loocv_rmse.csv"))

# ------------------------------------------------------------------------------
# 9. Graficas
# ------------------------------------------------------------------------------

col_tipo <- c("Urbana tráfico" = "#d73027", "Urbana fondo" = "#4575b4",
              "Suburbana" = "#1a9850")

# (a) Perfil diario medio por estacion (analogo Fig. 3 del paper).
ciclo <- rbindlist(lapply(seq_len(n), function(k) {
  data.table(ESTACION = sts$ESTACION[k], NOM_TIPO = sts$NOM_TIPO[k],
             hora = 1:24,
             NO2 = exp(tapply(Xmat[, k], hora_de_t, mean)))
}))
p_ciclo <- ggplot(ciclo, aes(hora, NO2, group = ESTACION, color = NOM_TIPO)) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  scale_color_manual(values = col_tipo, name = "Tipología") +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(title = "Perfil diario medio de NO₂ por estación (día tipo laborable, 2025)",
       subtitle = "Curvas suavizadas (B-spline) promediadas por hora sobre las 52 semanas",
       x = "Hora del día", y = "NO₂ (µg/m³, back-transf. de log)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_SALIDA, "perfil_diario_medio.png"), p_ciclo,
       width = 10, height = 6, dpi = 200, bg = "white")

# (b) Trace-variograma empirico + modelo esferico ajustado.
hseq <- seq(0, cutoff, length.out = 200)
p_vgm <- ggplot() +
  geom_point(data = vgram, aes(h, gamma, size = np), color = "#2166ac") +
  geom_line(data = data.frame(h = hseq, g = gamma_mod(hseq)),
            aes(h, g), color = "#d73027", linewidth = 1) +
  scale_size_continuous(name = "nº pares") +
  labs(title = "Trace-variograma del NO₂ funcional",
       subtitle = sprintf("Modelo esférico: nugget=%.2f, sill=%.2f, rango=%.0f m",
                          p_sph[1], p_sph[2], p_sph[3]),
       x = "Distancia (m)", y = expression(hat(gamma)(h))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_SALIDA, "trace_variograma.png"), p_vgm,
       width = 9, height = 6, dpi = 200, bg = "white")

# (c) LOOCV: perfil diario observado vs predicho para estaciones representativas.
sel <- res_loocv[, .SD[c(1, .N)], by = NOM_TIPO]$ESTACION
loocv_curvas <- rbindlist(lapply(sel, function(est) {
  i <- match(est, sts$ESTACION)
  j <- setdiff(seq_len(n), i)
  Gi <- gamma_mod(D[j, j]); diag(Gi) <- 0
  Ai <- rbind(cbind(Gi, 1), c(rep(1, n - 1), 0))
  lam <- solve(Ai, c(gamma_mod(D[j, i]), 1))[1:(n - 1)]
  pred <- as.vector(Xmat[, j, drop = FALSE] %*% lam)
  rbind(
    data.table(ESTACION = est, hora = 1:24, tipo = "Observado",
               NO2 = exp(tapply(Xmat[, i], hora_de_t, mean))),
    data.table(ESTACION = est, hora = 1:24, tipo = "OFK (predicho)",
               NO2 = exp(tapply(pred, hora_de_t, mean)))
  )
}))
loocv_curvas[, ESTACION := factor(ESTACION, levels = sel)]
p_loocv <- ggplot(loocv_curvas, aes(hora, NO2, color = tipo)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ESTACION, scales = "free_y") +
  scale_color_manual(values = c("Observado" = "#333333",
                                "OFK (predicho)" = "#d73027"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "Validación LOOCV: perfil diario observado vs. predicho por OFK",
       subtitle = "Estación con mejor y peor RMSE de cada tipología",
       x = "Hora del día", y = "NO₂ (µg/m³)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(DIR_SALIDA, "loocv_curvas.png"), p_loocv,
       width = 11, height = 7, dpi = 200, bg = "white")

# (d) RMSE LOOCV por estacion, coloreado por tipologia.
p_rmse <- ggplot(res_loocv, aes(reorder(ESTACION, RMSE_log), RMSE_log,
                                fill = NOM_TIPO)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = col_tipo, name = "Tipología") +
  labs(title = "RMSE de la validación LOOCV (log-NO₂) por estación",
       subtitle = "El OFK tiende a fallar más en estaciones periféricas/aisladas (paper, Sec. 4.3.2)",
       x = NULL, y = "RMSE integrado (log)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_SALIDA, "loocv_rmse_por_estacion.png"), p_rmse,
       width = 9, height = 7, dpi = 200, bg = "white")

# (e) Superficies OFK: media anual, punta de manana e incertidumbre.
bordes <- st_coordinates(st_cast(st_cast(st_geometry(mapa_utm),
                                         "MULTILINESTRING"), "LINESTRING"))
bordes <- as.data.frame(bordes)
df_est <- data.frame(X = XY[, 1], Y = XY[, 2], NO2 = exp(m_anual))

capa_borde <- geom_path(data = bordes, aes(X, Y, group = L1),
                        color = "white", linewidth = 0.2, inherit.aes = FALSE)
capa_est <- geom_point(data = df_est, aes(X, Y), shape = 21, size = 1.6,
                       fill = "white", color = "black", stroke = 0.4)

mapa_surf <- function(valor, titulo, leyenda, opcion = "plasma") {
  ggplot() +
    geom_tile(data = data.frame(X = GR[, 1], Y = GR[, 2], v = valor),
              aes(X, Y, fill = v), width = dx, height = dy) +
    capa_borde + capa_est +
    scale_fill_viridis_c(option = opcion, name = leyenda) +
    coord_equal() + labs(title = titulo) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.key.size = unit(0.4, "cm"))
}

fig <- (mapa_surf(surf_anual, "Media anual", "NO₂\n(µg/m³)") |
        mapa_surf(surf_punta, sprintf("Punta mañana (H%02d)", HORA_PUNTA), "NO₂\n(µg/m³)") |
        mapa_surf(trace_var, "Incertidumbre (trace-variance)", "σ²\n(log²)", "viridis")) +
  plot_annotation(
    title = "Kriging Funcional Ordinario (OFK) del NO₂ — Madrid 2025",
    subtitle = "Superficies derivadas de la MISMA predicción funcional. Puntos = estaciones.",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(color = "grey30")))
ggsave(file.path(DIR_SALIDA, "mapas_superficie_ofk.png"), fig,
       width = 16, height = 6, dpi = 200, bg = "white")

cat("\nSalidas guardadas en:\n", DIR_SALIDA, "\n")
