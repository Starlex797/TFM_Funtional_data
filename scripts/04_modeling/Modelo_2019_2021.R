# ==============================================================================
# COMPARACION DE DOS MODELOS ESPACIALES INLA-SPDE (2019-2021)
# ==============================================================================
#
# Objetivo metodologico
# ---------------------
# Comparar dos procedimientos sobre el mismo conjunto de datos:
#
#   M1. Modelo espacial base
#       Incluye temperatura, velocidad del viento, lluvia, presion,
#       intensidad y tipologia, ademas de un campo espacial para el NO2.
#
#   M2. Modelo con coeficiente espacialmente variable
#       Conserva exactamente la estructura de M1 y añade un segundo campo SPDE
#       que permite que el efecto de la velocidad del viento cambie
#       espacialmente.
#
# La ecuacion conceptual del segundo modelo es:
#
#   y(s,t) = beta_0 + beta' x(s,t) + w_NO2(s)
#            + Velocidad_Viento(s,t) * w_Viento(s) + error(s,t)
#
# De este modo, la RW2 de velocidad representa su efecto medio no lineal y el
# campo adicional representa su desviacion espacial.
# Todos los demas elementos se mantienen iguales para que la comparacion sea
# atribuible al campo adicional.
#
# Salidas principales
# -------------------
#   1. tabla_comparacion_modelos_2019_2021.csv
#   2. tabla_coeficientes_modelos_2019_2021.csv
#   3. tabla_holdout_espacial_2019_2021.csv
#   4. tabla_holdout_por_estacion_2019_2021.csv
#   5. predicciones_holdout_2019_2021.csv
#
# Nota: RMSE_ajuste se calcula con los mismos datos empleados para ajustar el
# modelo. Sirve como diagnostico descriptivo, no sustituye una
# validacion cruzada o un conjunto de prueba independiente.
# ==============================================================================


# ==============================================================================
# BLOQUE 1. CONFIGURACION DEL ANALISIS
# ==============================================================================

library(INLA)
library(data.table)
library(here)
library(Matrix)


source(here("R", "modeling", "spde_config.R"))
source(here("R", "modeling", "inla_modeling.R"))


# Periodo y escala de los datos.
ANIOS <- 2019:2021
ESCALA <- "DIARIO"
RESPUESTA <- "LOG_NO2_DIARIO"

# Nombre completo del maestro conjunto 2019-2021.
NOMBRE_ARCHIVO_DATOS <-
    "dataset_maestro_inla_20190101_20211231_DIARIO2.rds"
RUTA_DATOS <- here(
    "data", "processed", "Maestro", "diario", NOMBRE_ARCHIVO_DATOS
)

# Las cuatro covariables solicitadas. Llueve es binaria:
#   0 = no llueve; 1 = precipitacion diaria igual o superior a 1 mm.
COVARIABLES <- c(
    "Temperatura",
    "Velocidad_Viento",
    "Llueve",
    "Presion_Barometrica",
    "intensidad"
)

# Covariable cuyo coeficiente se deja variar espacialmente en M2.
# Para estudiar otra covariable basta con cambiar este unico nombre.
COVARIABLE_CAMPO <- "Velocidad_Viento"

# Numero de grupos de los efectos no lineales RW2.
N_GRUPOS_RW2 <- c(
    Temperatura = 40L,
    Velocidad_Viento = 40L
)
METODO_GRUPOS_RW2 <- "quantile"

# Malla y opciones comunes de INLA.
NOMBRE_MALLA <- "fina"
RUTA_MALLA <- here(
    "data", "processed", "Malla", "NO2",
    sprintf("malla_spde_madrid_%s.rds", NOMBRE_MALLA)
)

# PC priors introducidos manualmente para este analisis:
#   prior.range = c(r0, p) -> P(rango < r0) = p
#   prior.sigma = c(s0, p) -> P(sigma > s0) = p
# Cambia aqui los valores si quieres realizar un analisis de sensibilidad.
prior_spde <- list(
    prior.range = c(9.3, 0.5),
    prior.sigma = c(0.6, 0.2)
)

# Validacion espacial comun para M1 y M2. Las respuestas de estas estaciones
# se sustituyen por NA durante el segundo ajuste y se usan solo para evaluar.
CALCULAR_HOLDOUT <- TRUE
ESTACIONES_HOLDOUT <- c(
    "Plaza Castilla",
    "Casa de Campo",
    "Ensanche Vallecas"
)

FAMILIA <- "gaussian"
NUM_THREADS <- 5L
ESTRATEGIA_INTEGRACION <- "eb"
VERBOSE_INLA <- TRUE
GUARDAR_MODELOS <- TRUE

DIR_SALIDA <- here(
    "outputs", "modelo", "Modelo_2019_2021", "comparacion_campos_latentes"
)
DIR_MODELOS <- here(
    "data", "processed", "Modelos", "Modelo_2019_2021"
)
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
if (GUARDAR_MODELOS) {
    dir.create(DIR_MODELOS, recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================
# BLOQUE 2. CARGA Y PREPARACION DE LOS DATOS
# ==============================================================================

if (!file.exists(RUTA_DATOS)) {
    stop("No se encuentra el dataset maestro: ", RUTA_DATOS)
}
if (!file.exists(RUTA_MALLA)) {
    stop("No se encuentra la malla: ", RUTA_MALLA)
}
if (!COVARIABLE_CAMPO %in% COVARIABLES) {
    stop("COVARIABLE_CAMPO debe pertenecer al vector COVARIABLES.")
}

datos <- as.data.table(readRDS(RUTA_DATOS))

columnas_necesarias <- c(
    "ESTACION", "NOM_TIPO", "FECHA", "X_km", "Y_km",
    RESPUESTA, COVARIABLES
)
columnas_ausentes <- setdiff(columnas_necesarias, names(datos))
if (length(columnas_ausentes) > 0L) {
    stop(
        "Faltan columnas necesarias en el maestro: ",
        paste(columnas_ausentes, collapse = ", ")
    )
}

datos[, FECHA := as.Date(FECHA)]
datos <- datos[as.integer(format(FECHA, "%Y")) %in% ANIOS]

# Se emplean casos completos para que M1 y M2 utilicen exactamente las mismas
# observaciones. Esto es imprescindible para comparar WAIC y DIC.
n_antes <- nrow(datos)
datos <- datos[complete.cases(datos[, ..columnas_necesarias])]
n_eliminadas <- n_antes - nrow(datos)

if (nrow(datos) == 0L) {
    stop("No quedan observaciones completas para ajustar los modelos.")
}
if (!all(datos$Llueve %in% c(0, 1))) {
    stop("La variable Llueve debe estar codificada exclusivamente como 0 y 1.")
}

datos[, y := get(RESPUESTA)]

# Suburbana es la referencia. Se utilizan dos indicadores explicitos porque
# las formulas llevan -1 + Intercept; incluir el factor completo junto con el
# intercepto manual generaria una columna redundante.
datos[, NOM_TIPO := relevel(factor(NOM_TIPO), ref = "Suburbana")]
datos[, Tipo_Urbana_fondo := as.integer(NOM_TIPO == "Urbana fondo")]
datos[, Tipo_Urbana_trafico := as.integer(
    NOM_TIPO != "Suburbana" & NOM_TIPO != "Urbana fondo"
)]

# Las formulas utilizan estos dos indices discretos. Deben crearse antes del
# stack y viajar dentro de efectos_fijos para que inla() pueda encontrarlos.
datos[, Temperatura_rw2 := inla.group(
    Temperatura,
    n = N_GRUPOS_RW2[["Temperatura"]],
    method = METODO_GRUPOS_RW2
)]
datos[, Velocidad_Viento_rw2 := inla.group(
    Velocidad_Viento,
    n = N_GRUPOS_RW2[["Velocidad_Viento"]],
    method = METODO_GRUPOS_RW2
)]
setorder(datos, FECHA, ESTACION)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("DATOS DEL ANALISIS\n")
cat(strrep("=", 72), "\n", sep = "")
cat("Archivo:       ", NOMBRE_ARCHIVO_DATOS, "\n", sep = "")
cat(
    "Periodo:       ", as.character(min(datos$FECHA)), " a ",
    as.character(max(datos$FECHA)), "\n",
    sep = ""
)
cat("Observaciones: ", nrow(datos), "\n", sep = "")
cat("Estaciones:    ", uniqueN(datos$ESTACION), "\n", sep = "")
cat("Filas omitidas por datos incompletos: ", n_eliminadas, "\n", sep = "")


# ==============================================================================
# BLOQUE 3. MALLA Y CAMPOS LATENTES SPDE
# ==============================================================================

malla <- readRDS(RUTA_MALLA)

# Campo espacial comun que representa la variacion espacial residual del NO2.
spde_no2 <- crear_spde(
    malla,
    escala = ESCALA,
    prior.range = prior_spde$prior.range,
    prior.sigma = prior_spde$prior.sigma
)
indice_no2 <- inla.spde.make.index(
    name = "campo_no2",
    n.spde = spde_no2$n.spde
)

# Campo adicional e independiente para el coeficiente espacial de la
# covariable elegida. Utiliza los mismos PC priors que el campo del NO2.
# constr = TRUE centra el campo y separa su desviacion espacial del efecto fijo
# medio de la covariable.

spde_covariable <- inla.spde2.pcmatern(
    mesh = malla,
    alpha = 2,
    prior.range = prior_spde$prior.range,
    prior.sigma = prior_spde$prior.sigma,
    constr = TRUE
)
indice_covariable <- inla.spde.make.index(
    name = "campo_covariable",
    n.spde = spde_covariable$n.spde
)

coordenadas <- as.matrix(datos[, .(X_km, Y_km)])
A_no2 <- inla.spde.make.A(mesh = malla, loc = coordenadas)

# Matriz del coeficiente espacialmente variable:
# cada fila de A se multiplica por el valor observado de la covariable.
A_covariable <- Matrix::Diagonal(x = datos[[COVARIABLE_CAMPO]]) %*% A_no2

stopifnot(
    "La matriz A del NO2 tiene un numero de filas incorrecto" =
        nrow(A_no2) == nrow(datos),
    "La matriz A de la covariable tiene un numero de filas incorrecto" =
        nrow(A_covariable) == nrow(datos),
    "La matriz A y la malla tienen dimensiones incompatibles" =
        ncol(A_no2) == malla$n
)

# Efectos fijos compartidos por ambos modelos.
efectos_fijos <- data.frame(
    Intercept = rep(1, nrow(datos)),
    Temperatura_rw2 = datos$Temperatura_rw2,
    Velocidad_Viento_rw2 = datos$Velocidad_Viento_rw2,
    Llueve = datos$Llueve,
    Presion_Barometrica = datos$Presion_Barometrica,
    intensidad = datos$intensidad,
    Tipo_Urbana_fondo = datos$Tipo_Urbana_fondo,
    Tipo_Urbana_trafico = datos$Tipo_Urbana_trafico
)


# ==============================================================================
# BLOQUE 4. FORMULACION Y AJUSTE DE LOS MODELOS
# ==============================================================================

# M1: misma especificacion seleccionada en M5 y campo espacial del NO2.
formula_M1 <- y_response ~ -1 + Intercept +
    f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
    Llueve + Presion_Barometrica + intensidad +
    Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_no2, model = spde_no2)

# M2: M1 mas un campo que permite que el coeficiente del viento cambie
# entre localizaciones. A_covariable introduce el producto
# Velocidad_Viento(s,t) * w_Viento(s) en el predictor.
formula_M2 <- y_response ~ -1 + Intercept +
    f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    Llueve + Presion_Barometrica + intensidad +
    Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_no2, model = spde_no2) +
    f(campo_covariable, model = spde_covariable)

catalogo_modelos <- data.table(
    Modelo = c("M1", "M2"),
    Procedimiento = c(
        "Campo espacial del NO2",
        paste0(
            "Campo espacial del NO2 + coeficiente espacial de ",
            COVARIABLE_CAMPO
        )
    ),
    Formula = c(
        paste(deparse(formula_M1), collapse = " "),
        paste(deparse(formula_M2), collapse = " ")
    )
)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("MODELOS COMPARADOS\n")
cat(strrep("=", 72), "\n", sep = "")
print(catalogo_modelos[, .(Modelo, Procedimiento)], nrows = Inf)

stack_M1 <- crear_stack_inla(
    respuesta = datos$y,
    A = list(A_no2, 1),
    effects = list(indice_no2, efectos_fijos),
    tag = "estimacion"
)
stack_M2 <- crear_stack_inla(
    respuesta = datos$y,
    A = list(A_no2, A_covariable, 1),
    effects = list(indice_no2, indice_covariable, efectos_fijos),
    tag = "estimacion"
)

ajuste_M1 <- ajustar_modelo_cronometrado_inla(
    formula_modelo = formula_M1,
    stack = stack_M1,
    id_modelo = "M1",
    familia = FAMILIA,
    verbose = VERBOSE_INLA,
    num_threads = NUM_THREADS,
    calcular_cpo = TRUE,
    int_strategy = ESTRATEGIA_INTEGRACION
)
ajuste_M2 <- ajustar_modelo_cronometrado_inla(
    formula_modelo = formula_M2,
    stack = stack_M2,
    id_modelo = "M2",
    familia = FAMILIA,
    verbose = VERBOSE_INLA,
    num_threads = NUM_THREADS,
    calcular_cpo = TRUE,
    int_strategy = ESTRATEGIA_INTEGRACION
)


# ==============================================================================
# BLOQUE 5. HOLD-OUT ESPACIAL COMUN PARA M1 Y M2
# ==============================================================================

tabla_holdout <- NULL
tabla_holdout_estacion <- NULL
predicciones_holdout <- NULL
ajustes_holdout <- NULL

if (CALCULAR_HOLDOUT) {
    estaciones_ausentes <- setdiff(
        ESTACIONES_HOLDOUT,
        unique(datos$ESTACION)
    )
    if (length(estaciones_ausentes) > 0L) {
        stop(
            "No existen estas estaciones del hold-out: ",
            paste(estaciones_ausentes, collapse = ", ")
        )
    }

    es_holdout <- datos$ESTACION %in% ESTACIONES_HOLDOUT
    if (!any(es_holdout) || !any(!es_holdout)) {
        stop(
            "El hold-out debe contener observaciones de prueba y de ",
            "entrenamiento."
        )
    }

    y_holdout <- copy(datos$y)
    y_holdout[es_holdout] <- NA_real_

    stack_holdout_M1 <- crear_stack_inla(
        respuesta = y_holdout,
        A = list(A_no2, 1),
        effects = list(indice_no2, efectos_fijos),
        tag = "holdout"
    )
    stack_holdout_M2 <- crear_stack_inla(
        respuesta = y_holdout,
        A = list(A_no2, A_covariable, 1),
        effects = list(indice_no2, indice_covariable, efectos_fijos),
        tag = "holdout"
    )

    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("HOLD-OUT ESPACIAL\n")
    cat(strrep("=", 72), "\n", sep = "")
    print(unique(datos[
        es_holdout,
        .(ESTACION, NOM_TIPO = as.character(NOM_TIPO))
    ])[order(ESTACION)])
    cat("Entrenamiento: ", sum(!es_holdout), " observaciones\n", sep = "")
    cat("Prueba:        ", sum(es_holdout), " observaciones\n", sep = "")

    ajustes_holdout <- list(
        M1 = ajustar_holdout_espacial_inla(
            formula_modelo = formula_M1,
            id_modelo = "M1",
            descripcion = "Campo espacial del NO2",
            stack = stack_holdout_M1,
            indices_observaciones = indices_stack_inla(
                stack_holdout_M1,
                "holdout"
            ),
            datos = datos,
            es_holdout = es_holdout,
            familia = FAMILIA,
            verbose = VERBOSE_INLA,
            num_threads = NUM_THREADS
        ),
        M2 = ajustar_holdout_espacial_inla(
            formula_modelo = formula_M2,
            id_modelo = "M2",
            descripcion = paste0(
                "Campo espacial del NO2 + coeficiente espacial de ",
                COVARIABLE_CAMPO
            ),
            stack = stack_holdout_M2,
            indices_observaciones = indices_stack_inla(
                stack_holdout_M2,
                "holdout"
            ),
            datos = datos,
            es_holdout = es_holdout,
            familia = FAMILIA,
            verbose = VERBOSE_INLA,
            num_threads = NUM_THREADS
        )
    )

    tabla_holdout <- rbindlist(lapply(ajustes_holdout, `[[`, "global"))
    tabla_holdout_estacion <- rbindlist(lapply(
        ajustes_holdout,
        `[[`,
        "por_estacion"
    ))
    predicciones_holdout <- rbindlist(lapply(
        ajustes_holdout,
        `[[`,
        "predicciones"
    ))
    setorder(tabla_holdout, Modelo)
    setorder(tabla_holdout_estacion, Modelo, ESTACION)
    setorder(predicciones_holdout, Modelo, ESTACION, FECHA)
}


# ==============================================================================
# BLOQUE 6. COMPARACION DE LOS DOS PROCEDIMIENTOS
# ==============================================================================

tabla_comparacion <- rbindlist(list(
    resumir_ajuste_inla(
        ajuste = ajuste_M1,
        id_modelo = "M1",
        procedimiento = "Campo espacial del NO2",
        y_observada = datos$y
    ),
    resumir_ajuste_inla(
        ajuste = ajuste_M2,
        id_modelo = "M2",
        procedimiento = paste0(
            "Campo espacial del NO2 + coeficiente espacial de ",
            COVARIABLE_CAMPO
        ),
        y_observada = datos$y
    )
))

if (CALCULAR_HOLDOUT) {
    tabla_comparacion[
        tabla_holdout,
        on = .(Modelo),
        `:=`(
            RMSE_HOLDOUT = i.RMSE_HOLDOUT,
            MAE_HOLDOUT = i.MAE_HOLDOUT,
            Sesgo_HOLDOUT = i.Sesgo_HOLDOUT,
            COV95_HOLDOUT = i.COV95_HOLDOUT,
            Anchura_media_IC95_HOLDOUT = i.Anchura_media_IC95_HOLDOUT,
            Observaciones_HOLDOUT = i.Observaciones_HOLDOUT
        )
    ]
}

# Los deltas se calculan respecto a M1. Un valor negativo indica que M2 reduce
# la metrica correspondiente. Para WAIC, DIC y RMSE, menor es mejor.
referencia_M1 <- tabla_comparacion[Modelo == "M1"]
tabla_comparacion[, `:=`(
    Delta_WAIC_vs_M1 = WAIC - referencia_M1$WAIC,
    Delta_DIC_vs_M1 = DIC - referencia_M1$DIC,
    Delta_RMSE_vs_M1 = RMSE_ajuste - referencia_M1$RMSE_ajuste,
    Orden_WAIC = frank(WAIC, ties.method = "min")
)]
if (CALCULAR_HOLDOUT) {
    tabla_comparacion[, `:=`(
        Delta_RMSE_HOLDOUT_vs_M1 =
            RMSE_HOLDOUT - referencia_M1$RMSE_HOLDOUT,
        Delta_Anchura_HOLDOUT_vs_M1 =
            Anchura_media_IC95_HOLDOUT -
                referencia_M1$Anchura_media_IC95_HOLDOUT
    )]
}

# Tabla conjunta de coeficientes fijos. Las covariables continuas ya proceden
# estandarizadas del maestro; el coeficiente de Llueve compara dias con lluvia
# frente a dias sin lluvia. Un efecto es significativo al 95 % cuando su
# intervalo posterior no contiene cero.
tabla_coeficientes <- rbindlist(list(
    extraer_coeficientes_fijos_inla(
        modelo = ajuste_M1$modelo,
        id_modelo = "M1"
    ),
    extraer_coeficientes_fijos_inla(
        modelo = ajuste_M2$modelo,
        id_modelo = "M2"
    )
))


# ==============================================================================
# BLOQUE 7. EXPORTACION Y PRESENTACION DE RESULTADOS
# ==============================================================================

fwrite(
    tabla_comparacion,
    file.path(DIR_SALIDA, "tabla_comparacion_modelos_2019_2021.csv")
)
fwrite(
    tabla_coeficientes,
    file.path(DIR_SALIDA, "tabla_coeficientes_modelos_2019_2021.csv")
)
if (CALCULAR_HOLDOUT) {
    fwrite(
        tabla_holdout,
        file.path(DIR_SALIDA, "tabla_holdout_espacial_2019_2021.csv")
    )
    fwrite(
        tabla_holdout_estacion,
        file.path(DIR_SALIDA, "tabla_holdout_por_estacion_2019_2021.csv")
    )
    fwrite(
        predicciones_holdout,
        file.path(DIR_SALIDA, "predicciones_holdout_2019_2021.csv")
    )
}
if (GUARDAR_MODELOS) {
    saveRDS(
        ajuste_M1$modelo,
        file.path(DIR_MODELOS, "modelo_M1_espacial_NO2.rds")
    )
    saveRDS(
        ajuste_M2$modelo,
        file.path(DIR_MODELOS, "modelo_M2_campo_velocidad_viento.rds")
    )
}

cat("\n", strrep("=", 72), "\n", sep = "")
cat("TABLA DE COMPARACION\n")
cat(strrep("=", 72), "\n", sep = "")
print(tabla_comparacion, nrows = Inf)
cat("\nInterpretacion de los deltas:\n")
cat("  Delta < 0: M2 mejora la metrica respecto a M1.\n")
cat("  Delta > 0: M2 empeora la metrica respecto a M1.\n")
cat("  COV95 debe aproximarse a 95 %, pero aqui sigue siendo una medida de ajuste.\n")

if (CALCULAR_HOLDOUT) {
    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("RESULTADOS DEL HOLD-OUT ESPACIAL\n")
    cat(strrep("=", 72), "\n", sep = "")
    print(tabla_holdout, nrows = Inf)
    cat("\nResultados por estacion:\n")
    print(tabla_holdout_estacion, nrows = Inf)
    cat(
        "\nEl COV95 del hold-out debe interpretarse junto con RMSE, sesgo ",
        "y anchura.\n",
        sep = ""
    )
}

cat("\n", strrep("=", 72), "\n", sep = "")
cat("TABLA DE COEFICIENTES FIJOS\n")
cat(strrep("=", 72), "\n", sep = "")
print(tabla_coeficientes, nrows = Inf)
cat(
    "\nEl campo de ", COVARIABLE_CAMPO,
    " no tiene un unico coeficiente: es una superficie espacial.\n",
    sep = ""
)

if (interactive()) {
    View(tabla_comparacion, title = "Comparacion de modelos 2019-2021")
    View(tabla_coeficientes, title = "Coeficientes de los modelos 2019-2021")
    if (CALCULAR_HOLDOUT) {
        View(tabla_holdout, title = "Hold-out espacial M1 y M2")
        View(
            tabla_holdout_estacion,
            title = "Hold-out por modelo y estacion"
        )
    }
}

cat("\nResultados guardados en: ", DIR_SALIDA, "\n", sep = "")
