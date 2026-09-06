# Pruebas pequenas: pesos, k exacto, lluvia, muestra comun y eleccion manual.
suppressPackageStartupMessages(library(data.table))
source("R/interpolation/preparacion_comparacion_clima.R")
source("R/interpolation/comparacion_clima_multiescala.R")
source("R/interpolation/resumen_comparacion_clima.R")

Y <- rbind(c(100, 10, 20, 40), c(200, NA, 20, 40), c(300, NA, NA, 40))
D <- as.matrix(dist(c(0, 1, 2, 4)))
diag(D) <- Inf
p <- predicciones_rmse(Y, D, 1, 3)
stopifnot(
  identical(p$NN, c(10, 20, 40)), p$KNN[1, 2] == 15,
  isTRUE(all.equal(p$IDW1[1, 2], (10 + 20 / 2) / 1.5)),
  isTRUE(all.equal(p$IDW2[1, 2], (10 + 20 / 4) / 1.25)),
  p$KNN[2, 2] == 30, is.na(p$KNN[3, 2])
)
Y2 <- Y
Y2[, 1] <- -999
stopifnot(identical(p, predicciones_rmse(Y2, D, 1, 3)))

# Durante lluvia en el objetivo, los vecinos secos siguen siendo donantes.
lluvia <- rbind(c(2, 0, 0), c(0, 3, 0), c(0, 0, 0))
D3 <- D[1:3, 1:3]
red <- list(Y = lluvia, D = D3, sitios = data.table(Ubicacion = 1:3))
r <- calcular_curva_rmse(red, "Precipitaciones")
stopifnot(
  r$curva[Metodo == "Siempre cero" & k == 1, RMSE] == sqrt((4 + 9) / 2),
  all(r$curva$N == 2), all(r$curva$N_ubicaciones == 2),
  all(r$curva[Metodo == "KNN", RMSE] == sqrt((4 + 9) / 2)),
  !anyNA(r$curva$RMSE)
)
m <- muestra_comun_rmse(lluvia, "Precipitaciones")
stopifnot(sum(m$mascara) == 2, !any(m$mascara[3, ]), m$k_max == 2)

# Una tabla pendiente nunca escoge por su cuenta el menor RMSE.
k_usuario <- data.frame(Variable = "Precipitaciones", horario = NA)
tabla <- tabla_escala_manual(r$curva, k_usuario, "horario")
stopifnot(
  nrow(tabla) == 1, is.na(tabla$k_elegido),
  is.na(tabla$RMSE_KNN), tabla$Estado == "Pendiente"
)
k_usuario$horario <- 2
tabla <- tabla_escala_manual(r$curva, k_usuario, "horario")
stopifnot(
  nrow(tabla) == 1, tabla$k_elegido == 2,
  tabla$RMSE_KNN == r$curva[Metodo == "KNN" & k == 2, RMSE]
)
stopifnot(inherits(try(validar_k_manual(3, 1:2, "lluvia", "horario"),
  silent = TRUE
), "try-error"))

# La marca agregada OK no basta si existen horas imputadas.
dt <- data.table(
  ESTACION = c("a", "b", "c"), FECHA = as.Date("2025-01-01"),
  HORA = 0L, MES = "2025-01", Temperatura = c(2, 3, 4), Temperatura_estado = "OK"
)
hh <- data.table(
  ESTACION = c("a", "b", "c"), FECHA = as.Date("2025-01-01"),
  MES = "2025-01", Temperatura_estado = c("IMPUTADO", "OK", "OK")
)
coords <- data.table(
  ESTACION = c("a", "b", "c"), Ubicacion = 1:3,
  X = c(0, 1, 2), Y = c(0, 0, 0)
)
red <- preparar_variable_comparacion(dt, hh, coords, "Temperatura", "diario")
stopifnot(ncol(red$Y) == 2, !1 %in% red$sitios$Ubicacion)
hh[, Temperatura_estado := "OK"]
coords[ESTACION == "b", `:=`(Ubicacion = 1L, X = 0)]
red <- preparar_variable_comparacion(dt, hh, coords, "Temperatura", "diario")
stopifnot(ncol(red$Y) == 2, red$Y[1, 1] == 2.5)
cat("OK: RMSE, vecinos, lluvia, muestra comun, k manual e imputaciones.\n")

# Comprobar las salidas reales de la ultima ejecucion, si ya existe.
carpetas <- list.dirs("outputs/interpolacion", recursive = FALSE)
carpetas <- carpetas[grepl("rmse_k_manual_", basename(carpetas))]
carpetas <- carpetas[file.exists(file.path(carpetas, "METODOLOGIA.txt"))]
if (length(carpetas)) {
  carpeta <- carpetas[which.max(file.info(file.path(carpetas, "METODOLOGIA.txt"))$mtime)]
  entrada <- fread(file.path(carpeta, "archivos_entrada.csv"))
  stopifnot(nrow(entrada) == 3, all(unname(tools::md5sum(entrada$Archivo)) == entrada$MD5))
  for (escala in c("horario", "diario", "mensual")) {
    d <- fread(file.path(carpeta, escala, "curvas_rmse.csv"))
    stopifnot(
      uniqueN(d$Variable) == 6, all(is.finite(d$RMSE)),
      all(d[, .(n = uniqueN(N)), by = Variable]$n == 1)
    )
    for (v in unique(d$Variable)) {
      stopifnot(all(file.exists(file.path(
        carpeta, escala,
        paste0("RMSE_k_", v, c(".png", ".pdf"))
      ))))
    }
    elecciones <- data.frame(Variable = unique(d$Variable), valor = 1L)
    names(elecciones)[2] <- escala
    elegida <- tabla_escala_manual(d, elecciones, escala)
    stopifnot(
      nrow(elegida) == 6, all(elegida$k_elegido == 1),
      isTRUE(all.equal(elegida$RMSE_KNN, elegida$RMSE_1NN)),
      isTRUE(all.equal(elegida$RMSE_IDW_p1, elegida$RMSE_1NN)),
      isTRUE(all.equal(elegida$RMSE_IDW_p2, elegida$RMSE_1NN))
    )
  }
  cat("OK: 18 curvas, tres escalas, tablas manuales y archivos de entrada intactos.\n")
}
