# ==============================================================================
# DESALINEAMIENTO ESPACIAL: COVARIABLE ASIGNADA FRENTE A CAMPO LATENTE ESPACIAL
# ==============================================================================
# M1: covariable del maestro tratada como conocida.
# M2: mediciones meteorologicas z(s,t) = alpha_X + X(s) + error_meteo;
#     NO2: y(s,t) = efectos_comunes + beta_X X(s) + campo_NO2(s) + error_NO2.
# X(s) es un unico campo SPDE Matern espacial, sin AR1 ni grupos temporales.
# Las verosimilitudes de NO2 y meteorologia se ajustan conjuntamente. La copia
# directa del SPDE introduce beta_X X(s) y propaga la incertidumbre del campo.
# La covariable elegida tiene efecto lineal en ambos procedimientos para aislar
# el tratamiento del desalineamiento. El resto de terminos es comun.
#
# ==============================================================================
# BLOQUE 1. CONFIGURACION DEL ANALISIS
# ==============================================================================

library(INLA)
library(data.table)
library(here)
library(Matrix)


source(here("R", "modeling", "spde_config.R"))
source(here("R", "modeling", "inla_modeling.R"))
source(here("R", "modeling", "campo_meteorologico_inla.R"))
source(here("R", "utilities", "academic_quality_tables.R"))


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

# Covariable reconstruida desde sus estaciones meteorologicas.
# Opciones: "Velocidad_Viento", "Presion_Barometrica", "Radiacion_Solar" o
# "Temperatura".
COVARIABLE_CAMPO <- "Temperatura"
OPCIONES_CAMPO <- c(
    "Velocidad_Viento", "Presion_Barometrica", "Radiacion_Solar", "Temperatura"
)
# La covariable del campo se exige completa en el maestro solo cuando se usa:
# asi la radiacion no reduce la muestra de los ajustes de viento o presion.
COVARIABLES <- unique(c(COVARIABLES, COVARIABLE_CAMPO))
ARCHIVOS_CLIMA <- here(
    "data", "processed", "Clima", "diario",
    sprintf("meteo_madrid_%d_diario5.rds", ANIOS)
)
PRIOR_CAMPO_METEO <- list(prior.range = c(9.3, 0.5), prior.sigma = c(3, 0.2))

# Opciones de ejecucion. DIAS_PRUEBA solo acorta los datos para comprobar que
# el codigo y los stacks funcionan; no crea grupos temporales en el SPDE.
SOLO_PREPARAR <- identical(Sys.getenv("INLA_SOLO_PREPARAR"), "1")
DIAS_PRUEBA <- suppressWarnings(as.integer(Sys.getenv("INLA_DIAS_PRUEBA", "0")))
if (is.na(DIAS_PRUEBA) || DIAS_PRUEBA < 0L) DIAS_PRUEBA <- 0L

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
# Malla del campo meteorologico latente. Se toma de la carpeta NO2 porque
# cubre tanto las estaciones meteorologicas como las de contaminacion. Puede
# tener distinta resolucion, pero el campo contiene solo n.spde nodos.
NOMBRE_MALLA_METEO <- "media"
RUTA_MALLA_METEO <- here(
    "data", "processed", "Malla", "NO2",
    sprintf("malla_spde_madrid_%s.rds", NOMBRE_MALLA_METEO)
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
    "Ensanche Vallecas",
    "Villaverde Alto"
)

FAMILIA <- "gaussian"
NUM_THREADS <- 5L
ESTRATEGIA_INTEGRACION <- "eb" # integra la incertidumbre de hiperparametros
VERBOSE_INLA <- TRUE
GUARDAR_MODELOS <- TRUE

DIR_SALIDA <- here(
    "outputs", "modelo", "Modelo_2019_2021", "desalineamiento_campo_espacial",
    COVARIABLE_CAMPO, if (DIAS_PRUEBA > 0) "prueba_corta" else "2019_2021"
)
DIR_MODELOS <- here(
    "data", "processed", "Modelos", "Modelo_2019_2021", "desalineamiento_campo_espacial",
    COVARIABLE_CAMPO, if (DIAS_PRUEBA > 0) "prueba_corta" else "2019_2021"
)
DIR_TABLAS_ACADEMICAS <- file.path(DIR_SALIDA, "academic_tables")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLAS_ACADEMICAS, recursive = TRUE, showWarnings = FALSE)
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
if (!file.exists(RUTA_MALLA_METEO)) {
    stop("No se encuentra la malla meteorologica: ", RUTA_MALLA_METEO)
}
if (!COVARIABLE_CAMPO %in% COVARIABLES) {
    stop("COVARIABLE_CAMPO debe pertenecer al vector COVARIABLES.")
}

datos <- as.data.table(readRDS(RUTA_DATOS))
if (!COVARIABLE_CAMPO %in% OPCIONES_CAMPO) {
    stop("El campo debe ser uno de: ", paste(OPCIONES_CAMPO, collapse = ", "))
}
climatologia <- cargar_meteorologia_espacial(
    ARCHIVOS_CLIMA, ANIOS, COVARIABLE_CAMPO
)
meteo <- climatologia$datos
# Misma escala de la covariable en M1 y M2, sin mezclar estandarizaciones.
columna_raw <- paste0(COVARIABLE_CAMPO, "_raw")
if (!columna_raw %in% names(datos)) stop("Falta la columna fisica ", columna_raw)
datos[, (COVARIABLE_CAMPO) :=
    (get(columna_raw) - climatologia$centro) / climatologia$escala]

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
if (DIAS_PRUEBA > 0) {
    fin_prueba <- as.Date(sprintf("%d-01-01", min(ANIOS))) + DIAS_PRUEBA - 1L
    datos <- datos[FECHA <= fin_prueba]
    meteo <- meteo[FECHA <= fin_prueba]
}

# El filtro distingue entre covariables y respuesta, y el motivo es que INLA
# las trata de forma distinta:
#   - Las covariables NO pueden tener NA: entran en el predictor lineal y una
#     fila sin covariable no puede contribuir al ajuste. Es un filtro duro.
#   - La respuesta SI puede ser NA: INLA la trata como objetivo de prediccion,
#     igual que hace el bloque de hold-out al poner NA en las estaciones
#     reservadas. Esas filas no aportan verosimilitud, asi que no alteran WAIC
#     ni DIC y la comparacion entre M1 y M2 sigue siendo valida; a cambio, se
#     obtiene la prediccion del NO2 en esas estacion-dia y los diagnosticos de
#     residuos trabajan sobre una rejilla mas completa.
# Ambos modelos comparten exactamente las mismas filas, que es lo que exige la
# comparacion de WAIC y DIC.
columnas_obligatorias <- c(
    "ESTACION", "NOM_TIPO", "FECHA", "X_km", "Y_km", COVARIABLES
)

n_antes <- nrow(datos)
datos <- datos[complete.cases(datos[, ..columnas_obligatorias])]
n_eliminadas <- n_antes - nrow(datos)

if (nrow(datos) == 0L) {
    stop("No quedan observaciones completas para ajustar los modelos.")
}

n_respuesta_na <- sum(is.na(datos[[RESPUESTA]]))
if (n_respuesta_na == nrow(datos)) {
    stop("Ninguna fila tiene respuesta observada.")
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
cat("Filas omitidas por covariables incompletas: ", n_eliminadas, "\n", sep = "")
cat(
    "Filas con respuesta NA (se predicen, no se ajustan): ",
    n_respuesta_na, "\n",
    sep = ""
)
cat("Observaciones con respuesta: ", nrow(datos) - n_respuesta_na, "\n", sep = "")


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

# Campo meteorologico espacial: una malla y un unico vector latente X(s).
malla_meteo <- readRDS(RUTA_MALLA_METEO)
campo_meteo <- preparar_campo_meteorologico(
    meteo = meteo,
    datos_no2 = datos,
    prior = PRIOR_CAMPO_METEO,
    malla = malla_meteo
)
spde_meteo <- campo_meteo$spde
coordenadas <- as.matrix(datos[, .(X_km, Y_km)])
A_no2 <- inla.spde.make.A(mesh = malla, loc = coordenadas)
cat("Mediciones meteorologicas validas: ", nrow(meteo), "\n")
cat(
    "Malla meteorologica: ", campo_meteo$malla$n, " nodos; ",
    campo_meteo$n_estaciones_meteo, " estaciones meteorologicas; ",
    campo_meteo$n_estaciones_no2, " estaciones de NO2.\n"
)
fwrite(
    meteo[, .(Observaciones = .N, Estaciones = uniqueN(ESTACION)),
        by = .(Anio = as.integer(format(FECHA, "%Y")))
    ],
    file.path(DIR_SALIDA, "meteorological_data_audit.csv")
)

# Efectos fijos compartidos por ambos modelos.
efectos_fijos <- data.frame(
    Intercept = rep(1, nrow(datos)),
    Temperatura_rw2 = datos$Temperatura_rw2,
    Velocidad_Viento_rw2 = datos$Velocidad_Viento_rw2,
    Velocidad_Viento = datos$Velocidad_Viento,
    Llueve = datos$Llueve,
    Presion_Barometrica = datos$Presion_Barometrica,
    intensidad = datos$intensidad,
    Tipo_Urbana_fondo = datos$Tipo_Urbana_fondo,
    Tipo_Urbana_trafico = datos$Tipo_Urbana_trafico
)
# Covariable del campo en su version estandarizada comun (necesaria en M1
# cuando es Radiacion_Solar; para viento y presion ya esta arriba).
efectos_fijos[[COVARIABLE_CAMPO]] <- datos[[COVARIABLE_CAMPO]]


# ==============================================================================
# BLOQUE 4. FORMULACION Y AJUSTE DE LOS MODELOS
# ==============================================================================

# Formulas editables. La variable elegida entra linealmente en M1 y como
# beta * X(s) en M2; no se introduce ninguna estructura temporal.
if (COVARIABLE_CAMPO == "Velocidad_Viento") {
    formula_M1 <- y_response ~ -1 + Intercept +
        f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
        f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + Presion_Barometrica + Llueve + intensidad +
        Tipo_Urbana_fondo + Tipo_Urbana_trafico +
        f(campo_no2, model = spde_no2)

    formula_M2 <- y_response ~ -1 + Intercept + Intercept_meteo +
        f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
        Presion_Barometrica + Llueve + intensidad +
        Tipo_Urbana_fondo + Tipo_Urbana_trafico +
        f(campo_no2, model = spde_no2) +
        f(campo_meteo, model = spde_meteo) +
        f(campo_meteo_copia,
            copy = "campo_meteo", fixed = FALSE,
            hyper = list(beta = list(prior = "normal", param = c(0, 10)))
        )
} else if (COVARIABLE_CAMPO == "Radiacion_Solar") {
    # La radiacion entra linealmente en M1 y como beta * X(s) en M2; el resto de
    # terminos es identico en ambos modelos.
    formula_M1 <- y_response ~ -1 + Intercept +
        f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
        Radiacion_Solar + Presion_Barometrica + Llueve + intensidad +
        Tipo_Urbana_fondo + Tipo_Urbana_trafico +
        f(campo_no2, model = spde_no2)

    formula_M2 <- y_response ~ -1 + Intercept + Intercept_meteo +
        f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
        Presion_Barometrica + Llueve + intensidad +
        Tipo_Urbana_fondo + Tipo_Urbana_trafico +
        f(campo_no2, model = spde_no2) +
        f(campo_meteo, model = spde_meteo) +
        f(campo_meteo_copia,
            copy = "campo_meteo", fixed = FALSE,
            hyper = list(beta = list(prior = "normal", param = c(0, 10)))
        )
} else if (COVARIABLE_CAMPO == "Temperatura") {
    # La temperatura entra linealmente en M1 y como beta * X(s) en M2 (no como
    # RW2, para que ambos procedimientos estimen un unico coeficiente); el resto
    # de terminos es identico en ambos modelos.
    formula_M1 <- y_response ~ -1 + Intercept +
        f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
        f(Temperatura_rw2, model = "rw2", scale.model = TRUE) + Presion_Barometrica + Llueve + intensidad +
        Tipo_Urbana_fondo + Tipo_Urbana_trafico +
        f(campo_no2, model = spde_no2)

    formula_M2 <- y_response ~ -1 + Intercept + Intercept_meteo +
        f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
        Presion_Barometrica + Llueve + intensidad +
        Tipo_Urbana_fondo + Tipo_Urbana_trafico +
        f(campo_no2, model = spde_no2) +
        f(campo_meteo, model = spde_meteo) +
        f(campo_meteo_copia,
            copy = "campo_meteo", fixed = FALSE,
            hyper = list(beta = list(prior = "normal", param = c(0, 10)))
        )
} else {
    formula_M1 <- y_response ~ -1 + Intercept +
        f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
        Presion_Barometrica + Llueve + intensidad +
        Tipo_Urbana_fondo + Tipo_Urbana_trafico +
        f(campo_no2, model = spde_no2)

    formula_M2 <- y_response ~ -1 + Intercept + Intercept_meteo +
        f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
        f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
        Llueve + intensidad + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
        f(campo_no2, model = spde_no2) +
        f(campo_meteo, model = spde_meteo) +
        f(campo_meteo_copia,
            copy = "campo_meteo", fixed = FALSE,
        )
}
# El parametro del prior normal es la precision; ambos coeficientes usan
# N(0, 1/10), para que la comparacion no cambie el grado de regularizacion.
PRIOR_FIJOS_M1 <- list(
    mean = setNames(list(0), COVARIABLE_CAMPO),
    prec = c(
        setNames(list(10), COVARIABLE_CAMPO),
        list(default = 0.001)
    )
)

catalogo_modelos <- data.table(
    Modelo = c("M1", "M2"),
    Procedimiento = c(
        "Campo espacial del NO2",
        paste0(
            "Campo latente meteorologico espacial de ",
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
attr(stack_M1, "control.fixed") <- PRIOR_FIJOS_M1
stack_M2 <- crear_stack_conjunto_meteo(
    datos$y, A_no2, indice_no2, efectos_fijos, campo_meteo
)
if (SOLO_PREPARAR) {
    saveRDS(
        list(
            n_meteo = nrow(meteo),
            estaciones_meteo = campo_meteo$n_estaciones_meteo,
            estaciones_no2 = campo_meteo$n_estaciones_no2,
            nodos = campo_meteo$malla$n, formulas = catalogo_modelos
        ),
        file.path(DIR_SALIDA, "preparacion_verificada.rds")
    )
    cat("Preparacion completada; no se han ejecutado ajustes.\n")
} else {
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
        familia = rep("gaussian", 2L),
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
        attr(stack_holdout_M1, "control.fixed") <- PRIOR_FIJOS_M1
        # Se oculta NO2 en las cuatro estaciones, pero se conservan todas las
        # mediciones meteorologicas. El campo espacial se estima conjuntamente.
        stack_holdout_M2 <- crear_stack_conjunto_meteo(
            y_holdout, A_no2, indice_no2, efectos_fijos, campo_meteo,
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
                    "Campo latente meteorologico espacial de ",
                    COVARIABLE_CAMPO
                ),
                stack = stack_holdout_M2,
                indices_observaciones = indices_stack_inla(
                    stack_holdout_M2,
                    "holdout"
                ),
                datos = datos,
                es_holdout = es_holdout,
                familia = rep("gaussian", 2L),
                verbose = VERBOSE_INLA,
                num_threads = NUM_THREADS
            )
        )

        saveRDS(
            list(
                modelo = ajustes_holdout$M2$modelo,
                campo = extraer_campo_en_no2(
                    list(modelo = ajustes_holdout$M2$modelo, stack = stack_holdout_M2),
                    datos, climatologia, campo_meteo
                )
            ),
            file.path(DIR_MODELOS, "M2_holdout_campo_espacial.rds")
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
                "Campo latente meteorologico espacial de ",
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
        Delta_CPO_vs_M1 = CPO_mean_neg_log - referencia_M1$CPO_mean_neg_log,
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
            excluir = c("Intercept", "Intercept_meteo"),
            id_modelo = "M2"
        ),
        extraer_beta_campo(ajuste_M2$modelo, COVARIABLE_CAMPO)
    ))


    # ==============================================================================
    # BLOQUE 7. ACADEMIC TABLES IN ENGLISH
    # ==============================================================================

    METHOD_LABELS <- c(
        M1 = "Assigned meteorological covariate",
        M2 = paste(
            "Joint spatial latent",
            c(
                Velocidad_Viento = "wind",
                Presion_Barometrica = "pressure",
                Radiacion_Solar = "solar radiation",
                Temperatura = "temperature"
            )[[COVARIABLE_CAMPO]],
            "field"
        )
    )
    VARIABLE_LABELS <- c(
        Velocidad_Viento = "Wind speed",
        Llueve = "Rain indicator",
        Presion_Barometrica = "Barometric pressure",
        Radiacion_Solar = "Solar radiation",
        Temperatura = "Temperature",
        intensidad = "Traffic intensity",
        Tipo_Urbana_fondo = "Station type: urban background",
        Tipo_Urbana_trafico = "Station type: urban traffic",
        Intercept = "Intercept"
    )

    format_number <- function(x, digits = 2L) {
        ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "--")
    }

    translate_variable <- function(x) {
        translated <- unname(VARIABLE_LABELS[x])
        translated[is.na(translated)] <- gsub("_", " ", x[is.na(translated)])
        translated
    }

    translate_station_type <- function(x) {
        translated <- as.character(x)
        translated[grepl("suburb", translated, ignore.case = TRUE)] <- "Suburban"
        translated[grepl("fondo", translated, ignore.case = TRUE)] <-
            "Urban background"
        translated[grepl("tr.fico|trafico", translated, ignore.case = TRUE)] <-
            "Urban traffic"
        translated
    }

    # Table 1: comparison of the two spatial misalignment procedures.
    comparison_table_english <- tabla_comparacion[, .(
        Model = Modelo,
        Method = unname(METHOD_LABELS[Modelo]),
        WAIC = format_number(WAIC, 1L),
        `Delta WAIC` = format_number(Delta_WAIC_vs_M1, 1L),
        DIC = format_number(DIC, 1L),
        `Delta DIC` = format_number(Delta_DIC_vs_M1, 1L),
        `Mean -log(CPO)` = format_number(CPO_mean_neg_log, 4L),
        `Delta CPO` = format_number(Delta_CPO_vs_M1, 4L),
        `WAIC p_eff` = format_number(WAIC_p_eff, 1L),
        `Fitted RMSE` = format_number(RMSE_ajuste, 3L),
        `Delta fitted RMSE` = format_number(Delta_RMSE_vs_M1, 3L)
    )]

    fwrite(
        comparison_table_english,
        file.path(DIR_TABLAS_ACADEMICAS, "misalignment_method_comparison.csv")
    )
    path_comparison_table <- booktabs_png(
        comparison_table_english,
        file.path(
            DIR_TABLAS_ACADEMICAS,
            "table_misalignment_method_comparison.png"
        ),
        title = "Comparison of spatial misalignment methods",
        subtitle = "Daily NO2 models for Madrid, 2019-2021",
        note = paste0(
            "Delta values are relative to M1; negative values favour M2. Lower WAIC, ",
            "DIC, mean -log(CPO) and RMSE values are preferred. Mean -log(CPO) is the ",
            "average leave-one-out logarithmic score based on valid CPO values. CPO ",
            "failures among observed responses: ",
            paste(
                sprintf("%s=%d", tabla_comparacion$Modelo, tabla_comparacion$CPO_failures),
                collapse = "; "
            ),
            ". All criteria use NO2 observations only; meteorological likelihood terms are excluded. Fitted RMSE is an in-sample diagnostic."
        ),
        widths = c(
            0.50, 3.10, 0.75, 0.82, 0.75, 0.82, 0.98, 0.82, 0.78, 0.88,
            1.05
        ),
        align = c("center", "left", rep("right", 9)),
        font_size = 7.3,
        row_height = 0.34
    )

    # Table 2: common spatial hold-out comparison.
    path_holdout_table <- NULL
    path_holdout_station_table <- NULL
    if (CALCULAR_HOLDOUT) {
        holdout_table_english <- tabla_comparacion[, .(
            Model = Modelo,
            Method = unname(METHOD_LABELS[Modelo]),
            RMSE = format_number(RMSE_HOLDOUT, 3L),
            `Delta RMSE` = format_number(Delta_RMSE_HOLDOUT_vs_M1, 3L),
            MAE = format_number(MAE_HOLDOUT, 3L),
            Bias = format_number(Sesgo_HOLDOUT, 3L),
            `COV95 (%)` = format_number(COV95_HOLDOUT, 1L),
            `Mean 95% width` = format_number(Anchura_media_IC95_HOLDOUT, 3L),
            `Delta width` = format_number(Delta_Anchura_HOLDOUT_vs_M1, 3L),
            Observations = as.character(Observaciones_HOLDOUT)
        )]
        holdout_station_table_english <- tabla_holdout_estacion[, .(
            Model = Modelo,
            Station = ESTACION,
            `Station type` = translate_station_type(NOM_TIPO),
            Observations = as.character(Observaciones),
            RMSE = format_number(RMSE, 3L),
            MAE = format_number(MAE, 3L),
            Bias = format_number(Sesgo, 3L),
            `COV95 (%)` = format_number(COV95, 1L),
            `Mean 95% width` = format_number(Anchura_media_IC95, 3L)
        )]

        fwrite(
            holdout_table_english,
            file.path(DIR_TABLAS_ACADEMICAS, "spatial_holdout_comparison.csv")
        )
        fwrite(
            holdout_station_table_english,
            file.path(DIR_TABLAS_ACADEMICAS, "spatial_holdout_by_station.csv")
        )
        path_holdout_table <- booktabs_png(
            holdout_table_english,
            file.path(
                DIR_TABLAS_ACADEMICAS,
                "table_spatial_holdout_comparison.png"
            ),
            title = "Spatial hold-out comparison of misalignment methods",
            subtitle = paste(
                "Stations excluded jointly:",
                paste(ESTACIONES_HOLDOUT, collapse = ", ")
            ),
            note = paste0(
                "Metrics are calculated only from excluded-station observations. Negative ",
                "deltas favour M2. Predictive COV95 includes observation variance and must ",
                "be interpreted together with RMSE and mean interval width."
            ),
            widths = c(0.50, 3.10, 0.72, 0.82, 0.70, 0.70, 0.78, 0.98, 0.82, 0.88),
            align = c("center", "left", rep("right", 8)),
            font_size = 7.4,
            row_height = 0.34
        )
        path_holdout_station_table <- booktabs_png(
            holdout_station_table_english,
            file.path(
                DIR_TABLAS_ACADEMICAS,
                "table_spatial_holdout_by_station.png"
            ),
            title = "Spatial hold-out performance by station",
            subtitle = "Two spatial misalignment methods; Madrid, 2019-2021",
            note = paste0(
                "The listed stations are excluded simultaneously during fitting. Coverage ",
                "must be assessed together with prediction error and interval width."
            ),
            widths = c(0.50, 1.48, 1.32, 0.82, 0.72, 0.68, 0.68, 0.78, 0.96),
            align = c("center", "left", "left", rep("right", 6)),
            font_size = 7.6,
            row_height = 0.32
        )
    }

    # Tables 3-4: fixed coefficients, separately for M1 and M2.
    coefficient_tables_english <- setNames(lapply(c("M1", "M2"), function(id) {
        table_model <- copy(tabla_coeficientes[Modelo == id])
        table_model[, Variable := translate_variable(Variable)]
        table_model[, Significativa_95 := fifelse(
            Significativa_95 == "Si",
            "Yes",
            "No"
        )]
        table_model[, Modelo := NULL]
        setnames(
            table_model,
            c("Variable", "Coeficiente", "IC95", "Significativa_95"),
            c("Variable", "Posterior mean", "95% credible interval", "Significant (95%)")
        )
        table_model[, `Posterior mean` := format_number(`Posterior mean`, 4L)]
        table_model
    }), c("M1", "M2"))

    coefficient_table_paths <- setNames(vapply(c("M1", "M2"), function(id) {
        table_model <- coefficient_tables_english[[id]]
        fwrite(
            table_model,
            file.path(DIR_TABLAS_ACADEMICAS, sprintf("fixed_effects_%s.csv", id))
        )
        booktabs_png(
            table_model,
            file.path(
                DIR_TABLAS_ACADEMICAS,
                sprintf("table_fixed_effects_%s.png", id)
            ),
            title = sprintf("Fixed effects for model %s", id),
            subtitle = METHOD_LABELS[[id]],
            note = paste0(
                "A fixed effect is labelled significant when its 95% posterior credible ",
                "interval excludes zero. Station-type effects are differences from the ",
                "suburban reference category. In M2, the meteorological coefficient is the ",
                "estimated copy scaling parameter. RW2 curves are reported separately."
            ),
            widths = c(2.55, 1.15, 1.85, 1.25),
            align = c("left", "right", "center", "center"),
            font_size = 8.1,
            row_height = 0.31
        )
    }, character(1)), c("M1", "M2"))

    expected_academic_png <- c(
        path_comparison_table,
        path_holdout_table,
        path_holdout_station_table,
        unname(coefficient_table_paths)
    )
    expected_academic_png <- expected_academic_png[
        !is.na(expected_academic_png) & nzchar(expected_academic_png)
    ]
    if (!all(file.exists(expected_academic_png))) {
        stop("One or more academic PNG tables could not be generated.")
    }


    # ==============================================================================
    # BLOQUE 8. EXPORTACION Y PRESENTACION DE RESULTADOS
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
            file.path(DIR_MODELOS, "modelo_M2_campo_meteorologico_espacial.rds")
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
        " es espacial y su coeficiente beta se incluye en la tabla.\n",
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


    # Proyeccion posterior en unidades fisicas (sin ruido de medicion meteorologica).
    campo_en_no2 <- extraer_campo_en_no2(
        ajuste_M2, datos, climatologia, campo_meteo
    )
    fwrite(campo_en_no2, file.path(DIR_SALIDA, "campo_meteorologico_en_NO2.csv"))
    saveRDS(
        list(
            malla = campo_meteo$malla,
            climatologia = climatologia,
            prior = PRIOR_CAMPO_METEO,
            formulas = catalogo_modelos, estaciones_holdout = ESTACIONES_HOLDOUT
        ),
        file.path(DIR_MODELOS, "configuracion_campo_espacial.rds")
    )
    writeLines(capture.output(sessionInfo()), file.path(DIR_SALIDA, "sessionInfo.txt"))
} # fin de ajustes; SOLO_PREPARAR permite verificar el stack completo
