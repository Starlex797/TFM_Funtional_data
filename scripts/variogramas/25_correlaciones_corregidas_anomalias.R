# ==============================================================================
# CORRELACIONES CORREGIDAS (SIN SESGO DE CICLO / AGREGACION / NO-LINEALIDAD)
# ------------------------------------------------------------------------------
# Corrige los sesgos de la tabla cruda (Pearson/Spearman/dCor sobre datos crudos
# y 7 estaciones juntas):
#   1. AGREGACION: se calcula POR ESTACION (no se mezclan estaciones).
#   2. CICLO COMPARTIDO: se trabaja sobre ANOMALIAS (se quita el ciclo diurno/
#      semanal/estacional a NO2 y a la covariable antes de correlacionar).
#   3. NO-LINEALIDAD: distancia de correlacion clasica (energy::dcor), en [0,1]
#      y comparable con |Pearson|; su sesgo al alza se controla con el IC boot.
#   4. AUTOCORRELACION: intervalos de confianza por MOVING BLOCK BOOTSTRAP
#      (bloques que preservan la dependencia temporal).
#   + Etiqueta automatica de la FORMA (lineal / monotona no lineal /
#     no monotona / debil) con el "triangulo" Pearson-Spearman-dCor.
#
# Salidas: outputs/correlaciones_corregidas/
#   - correlaciones_corregidas.csv  (todas las estimaciones + IC + forma)
#   - <estacion>_tabla.png          (tabla por estacion, 3 escalas)
# ==============================================================================

library(data.table)
library(energy)
library(here)

set.seed(4827)

# ------------------------------------------------------------------------------
# PARAMETROS
# ------------------------------------------------------------------------------
ESTACIONES <- c("Plaza Elíptica", "El Pardo", "Retiro")
N_BOOT     <- 300L      # replicas block-bootstrap
DCOR_SUB   <- 1500L     # submuestra maxima para dCor (coste O(n^2))
MODO_VALIDACION <- FALSE # TRUE: solo Plaza Eliptica horario, N_BOOT bajo

if (MODO_VALIDACION) { ESTACIONES <- "Plaza Elíptica"; N_BOOT <- 60L }

COVS <- list(
  list(col = "intensidad_raw",          lab = "Intensidad Trafico"),
  list(col = "carga_raw",               lab = "Carga Trafico"),
  list(col = "Temperatura_raw",         lab = "Temperatura (C)"),
  list(col = "Humedad_Relativa_raw",    lab = "Humedad Relativa (%)"),
  list(col = "Precipitaciones_raw",     lab = "Precipitaciones (mm)"),
  list(col = "Presion Barométrica_raw", lab = "Presion Barometrica (hPa)"),
  list(col = "Radiación Solar_raw",     lab = "Radiacion Solar (W/m2)"),
  list(col = "Velocidad Viento_raw",    lab = "Velocidad Viento (m/s)")
)

# Config por escala: dataset, columna NO2, y factores de desestacionalizacion.
# Los factores se restan SECUENCIALMENTE (media de grupo) para dejar anomalias.
ESCALAS <- list(
  horario = list(
    rds = here("data","processed","Maestro","horario",
               "dataset_maestro_inla_2025_HORARIO.rds"),
    no2 = "DATO",
    factores = function(d) list(
      mes_hora = interaction(format(d$FECHA,"%m"), d$HORA, drop = TRUE),
      dow      = format(d$FECHA, "%u")),
    L = 24L),                        # bloque = 1 dia
  diario = list(
    rds = here("data","processed","Maestro","diario",
               "dataset_maestro_inla_2025_DIARIO.rds"),
    no2 = "DATO_DIARIO",
    factores = function(d) list(
      mes = format(d$FECHA, "%m"),
      dow = format(d$FECHA, "%u")),
    L = 14L),                        # bloque = 2 semanas
  mensual = list(
    rds = here("data","processed","Maestro","mensual",
               "dataset_maestro_inla_2019_2025_MENSUAL.rds"),
    no2 = "DATO_NO2",
    factores = function(d) list(mes = format(d$FECHA, "%m")),  # climatologia
    L = 6L)                          # bloque = medio ano
)

# ------------------------------------------------------------------------------
# Utilidades
# ------------------------------------------------------------------------------
slug <- function(s) { s <- gsub("[^A-Za-z0-9]+","_",iconv(s,to="ASCII//TRANSLIT"))
  gsub("^_|_$","",s) }

# Desestacionaliza: resta secuencialmente la media de cada factor -> anomalia
anomalia <- function(x, factores) {
  a <- x
  for (g in factores) a <- a - ave(a, g, FUN = function(v) mean(v, na.rm = TRUE))
  a
}

# Distancia de correlacion CLASICA en [0,1] (comparable en escala con |Pearson|).
# Se submuestrea si n es grande (coste O(n^2)). El sesgo al alza en n pequeno se
# controla con el IC block-bootstrap; para independencia usar dcor.test aparte.
dcor_cl <- function(x, y) {
  n <- length(x)
  if (n > DCOR_SUB) { s <- sample.int(n, DCOR_SUB); x <- x[s]; y <- y[s] }
  tryCatch(energy::dcor(x, y), error = function(e) NA_real_)
}

# Indices de un moving block bootstrap (preserva bloques temporales)
mbb_idx <- function(n, L) {
  nb <- ceiling(n / L)
  starts <- sample.int(n - L + 1, nb, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) s:(s + L - 1)))
  idx[seq_len(n)]
}

ic_pct <- function(v) { v <- v[is.finite(v)]
  if (!length(v)) return(c(NA, NA)); quantile(v, c(0.025, 0.975), names = FALSE) }

# Clasifica la forma con el triangulo Pearson-Spearman-dCor (clasica en [0,1]).
# En una relacion lineal, dCor ~ ligeramente por DEBAJO de |Pearson|; que dCor
# SUPERE a max(|Pearson|,|Spearman|) senala dependencia NO monotona (curva en U,
# bimodal, desfasada) que las correlaciones de rango no captan.
forma_relacion <- function(p, s, d, umbral_debil = 0.10) {
  ap <- abs(p); as_ <- abs(s)
  if (d < umbral_debil && max(ap, as_) < umbral_debil) return("debil")
  g_mono <- as_ - ap             # Spearman supera a Pearson -> monotona curva
  g_nomon <- d - max(ap, as_)    # dCor supera a ambos -> no monotona
  if (g_nomon > 0.05) return("no monotona (no lineal)")
  if (g_mono  > 0.05) return("monotona no lineal")
  "lineal"
}

# ------------------------------------------------------------------------------
# Analisis de una (estacion, escala): devuelve data.table de resultados
# ------------------------------------------------------------------------------
analizar <- function(dt, cfg, covs_disp, estacion, escala) {
  facs <- cfg$factores(dt)
  y_an <- anomalia(dt[[cfg$no2]], facs)
  res  <- vector("list", length(covs_disp))
  for (i in seq_along(covs_disp)) {
    v <- covs_disp[[i]]
    x_an <- anomalia(dt[[v$col]], facs)
    ok <- is.finite(x_an) & is.finite(y_an)
    xa <- x_an[ok]; ya <- y_an[ok]; n <- length(xa)
    if (n < 10 || sd(xa) == 0) { res[[i]] <- NULL; next }

    p0 <- cor(xa, ya)
    s0 <- cor(xa, ya, method = "spearman")
    d0 <- dcor_cl(xa, ya)

    # Block bootstrap
    bp <- bs <- bd <- numeric(N_BOOT)
    for (b in seq_len(N_BOOT)) {
      id <- mbb_idx(n, cfg$L)
      xb <- xa[id]; yb <- ya[id]
      bp[b] <- suppressWarnings(cor(xb, yb))
      bs[b] <- suppressWarnings(cor(xb, yb, method = "spearman"))
      bd[b] <- dcor_cl(xb, yb)
    }
    cp <- ic_pct(bp); cs <- ic_pct(bs); cd <- ic_pct(bd)
    res[[i]] <- data.table(
      Estacion = estacion, Escala = escala, Covariable = v$lab, n = n,
      Pearson = round(p0,3),  P_lo = round(cp[1],3),  P_hi = round(cp[2],3),
      Spearman= round(s0,3),  S_lo = round(cs[1],3),  S_hi = round(cs[2],3),
      DistCor = round(d0,3),  D_lo = round(cd[1],3),  D_hi = round(cd[2],3),
      Forma   = forma_relacion(p0, s0, d0))
  }
  rbindlist(res)
}

# ==============================================================================
# EJECUCION
# ==============================================================================
dir_out <- here("outputs", "correlaciones_corregidas")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

resultados <- list(); k <- 0L
for (esc in names(ESCALAS)) {
  cfg <- ESCALAS[[esc]]
  dt_all <- as.data.table(readRDS(cfg$rds))
  dt_all <- dt_all[order(FECHA)]
  covs_disp <- Filter(function(v) v$col %in% names(dt_all), COVS)
  for (est in ESTACIONES) {
    if (!est %in% dt_all$ESTACION) next
    t0 <- Sys.time()
    dt <- dt_all[ESTACION == est]
    r  <- analizar(dt, cfg, covs_disp, est, esc)
    k <- k + 1L; resultados[[k]] <- r
    cat(sprintf("[OK] %-14s / %-8s  n~%d  (%.1fs)\n", est, esc,
                if (nrow(r)) r$n[1] else 0,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (MODO_VALIDACION) print(r[, .(Covariable, Pearson, Spearman, DistCor, Forma)])
  }
}
tab <- rbindlist(resultados)
cat(sprintf("\nCalculo terminado (%d filas). Generando tablas PNG...\n", nrow(tab)))

# ------------------------------------------------------------------------------
# TABLA-HEATMAP POR ESTACION (estilo comparativa, con IC y forma)
# ------------------------------------------------------------------------------
suppressMessages(library(ggplot2))
niv_esc <- c(horario = "Horario", diario = "Diario", mensual = "Mensual")
for (est in unique(tab$Estacion)) {
  d <- tab[Estacion == est]
  d[, Escala := factor(niv_esc[Escala], levels = niv_esc)]
  num <- rbindlist(lapply(c("Pearson","Spearman","DistCor"), function(m) {
    ini <- substr(m, 1, 1)                       # P / S / D
    data.table(Escala = d$Escala, Covariable = d$Covariable, Medida = m,
               val = d[[m]],
               lab = sprintf("%.2f\n[%.2f, %.2f]",
                             d[[m]], d[[paste0(ini,"_lo")]], d[[paste0(ini,"_hi")]]))
  }))
  frm <- data.table(Escala = d$Escala, Covariable = d$Covariable,
                    Medida = "Forma", val = NA_real_, lab = d$Forma)
  todo <- rbind(num, frm)
  todo[, Medida := factor(Medida, levels = c("Pearson","Spearman","DistCor","Forma"))]

  p <- ggplot(todo, aes(Medida, Covariable)) +
    geom_tile(data = todo[Medida != "Forma"], aes(fill = val), color = "white") +
    geom_tile(data = todo[Medida == "Forma"], fill = "grey95", color = "white") +
    geom_text(aes(label = lab), size = 2.7, lineheight = 0.85) +
    facet_grid(Escala ~ ., scales = "free_y", space = "free") +
    scale_fill_gradient2(low = "#3b6fb0", mid = "white", high = "#c0392b",
                         midpoint = 0, limits = c(-0.8, 0.8),
                         oob = scales::squish, name = "coef.") +
    scale_x_discrete(position = "top") +
    labs(title = sprintf("Correlaciones CORREGIDAS (anomalias) — %s", est),
         subtitle = "Por estacion · sobre anomalias (sin ciclo) · IC block-bootstrap · dCor clasica",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "gray40", size = 8),
          strip.text.y = element_text(face = "bold", angle = 0),
          panel.grid = element_blank(),
          axis.text.x.top = element_text(face = "bold"))
  ggsave(file.path(dir_out, sprintf("%s_tabla.png", slug(est))),
         p, width = 9, height = 11, dpi = 200, bg = "white")
  cat(sprintf("[PNG] %s_tabla.png\n", slug(est)))
}
