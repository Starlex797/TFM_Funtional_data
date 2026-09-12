# ==============================================================================
# SELECCION DE VARIABLES PARA EL MODELO INLA-SPDE (2019-2021)
# ==============================================================================
# Objetivo:
#   Comparar las ocho formulas M0-M7 escritas manualmente en el bloque 6.
#   Todos los modelos usan los mismos datos, la misma malla y el mismo campo
#   espacial; solamente cambia la combinacion de covariables.
#
# Metricas principales (calculadas con R/utilities/metricas_predictivas.R):
#   - WAIC y DIC: cuanto menor, mejor.
#   - RMSE_INLA_orientativo: calculado con la media posterior ajustada por INLA.
#     Es rapido, pero no es un RMSE de validacion cruzada.
#   - COV95_INLA_orientativo: porcentaje observado dentro del intervalo formado
#     con mean +/- 1.96 * sd. La sd de summary.fitted.values solo recoge la
#     incertidumbre de la media ajustada, asi que se le suma la varianza
#     residual gaussiana (varianza_residual_gaussiana()) antes de construir el
#     intervalo; sin ese termino el intervalo queda demasiado estrecho y COV95
#     sale muy por debajo de 95 aunque el modelo este bien especificado. Sigue
#     siendo una medida de ajuste orientativa, no una cobertura predictiva LOSO.
#
# Validacion espacial:
#   Se reserva un conjunto fijo de estaciones completas. M5 se reajusta una
#   sola vez sin utilizar sus respuestas y las predicciones se evaluan solo en
#   esas estaciones. Es un hold-out espacial, no una validacion LOSO.
# ==============================================================================

library(INLA)
library(data.table)
library(here)

source(here("R", "modeling", "spde_config.R"))
source(here("R", "utilities", "metricas_predictivas.R"))
source(here("R", "modeling", "inla_modeling.R"))


# ==============================================================================
# 1. CONFIGURACION
# ==============================================================================
ANIOS <- 2019:2021
RESPUESTA <- "LOG_NO2_DIARIO"

# Escriba aqui el nombre completo del unico maestro conjunto que quiere usar,
# incluida la extension .rds. No se construye ni se modifica automaticamente.
NOMBRE_ARCHIVO_DATOS <- "dataset_maestro_inla_20190101_20211231_DIARIO2.rds"
DIRECTORIO_DATOS <- here("data", "processed", "Maestro", "diario")
RUTA_DATOS <- file.path(DIRECTORIO_DATOS, NOMBRE_ARCHIVO_DATOS)

COVARIABLES <- c(
    "intensidad",
    "Temperatura",
    "Velocidad_Viento",
    "Velocidad_Viento_sqrt",
    "Presion_Barometrica",
    "Humedad_Relativa",
    "Precipitaciones",
    "Llueve",
    "Radiacion_Solar"
)

# Numero de grupos para las variables que entran como RW2. Estos valores se
# pueden sustituir despues por los elegidos en el analisis de sensibilidad.
N_GRUPOS_RW2 <- c(
    Temperatura = 40L,
    Velocidad_Viento = 40L
)
METODO_GRUPOS_RW2 <- "quantile"

MALLA <- "fina"
FAMILIA <- "gaussian"
NUM_THREADS <- 5L
VERBOSE_INLA <- TRUE

# Hold-out espacial exclusivo para M5. Los nombres deben coincidir exactamente
# con la columna ESTACION del dataset maestro.
CALCULAR_HOLDOUT <- TRUE
MODELOS_HOLDOUT <- "M5"
ESTACIONES_HOLDOUT <- c(
    "Plaza Castilla", # Urbana trafico
    "Casa de Campo", # Suburbana
    "Ensanche Vallecas" # Urbana fondo
)

DIR_OUT <- here(
    "outputs", "modelo", "Modelo_1", "seleccion_variables", "2019_2021"
)
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)



# ==============================================================================
# 2. CARGA Y PREPARACION DE LOS DATOS
# ==============================================================================
if (!file.exists(RUTA_DATOS)) {
    stop("No se encuentra el maestro conjunto: ", RUTA_DATOS)
}

cat("Maestro utilizado: ", RUTA_DATOS, "\n", sep = "")
df <- as.data.table(readRDS(RUTA_DATOS))

columnas_necesarias <- c(
    "ESTACION", "FECHA", "X_km", "Y_km", RESPUESTA,
    COVARIABLES, "NOM_TIPO"
)
columnas_ausentes <- setdiff(columnas_necesarias, names(df))
if (length(columnas_ausentes) > 0L) {
    stop("Faltan columnas en el maestro: ", paste(columnas_ausentes, collapse = ", "))
}

df[, FECHA := as.Date(FECHA)]
df <- df[as.integer(format(FECHA, "%Y")) %in% ANIOS]
df <- df[complete.cases(df[, ..columnas_necesarias])]
if (!nrow(df)) stop("No quedan observaciones completas para ajustar los modelos.")

# NOM_TIPO ya viene del maestro. Se fija "Suburbana" como categoria de
# referencia. Como la formula utiliza -1 + Intercept, se crean expresamente
# dos indicadores para evitar introducir las tres categorias junto al
# intercepto manual (una combinacion linealmente redundante).
df[, NOM_TIPO := relevel(factor(NOM_TIPO), ref = "Suburbana")]
df[, Tipo_Urbana_fondo := as.integer(NOM_TIPO == "Urbana fondo")]
df[, Tipo_Urbana_trafico := as.integer(
    NOM_TIPO != "Suburbana" & NOM_TIPO != "Urbana fondo"
)]
df[, y := get(RESPUESTA)]
setorder(df, FECHA, ESTACION)

if (interactive()) View(df, title = paste("Maestro:", NOMBRE_ARCHIVO_DATOS))


# =======================================================
# INLA necesita indices discretos para RW2. Las variables originales se
# conservan para poder utilizarlas tambien como efectos lineales.
df[, Temperatura_rw2 := inla.group(
    Temperatura,
    n = N_GRUPOS_RW2[["Temperatura"]],
    method = METODO_GRUPOS_RW2
)]
df[, Velocidad_Viento_rw2 := inla.group(
    Velocidad_Viento,
    n = N_GRUPOS_RW2[["Velocidad_Viento"]],
    method = METODO_GRUPOS_RW2
)]
cat("\nDatos utilizados\n")
cat("Anios:", paste(ANIOS, collapse = ", "), "\n")
cat("Observaciones:", nrow(df), "\n")
cat("Estaciones:", uniqueN(df$ESTACION), "\n")


# ==============================================================================
# 3. MALLA Y CAMPO ESPACIAL COMUN
# ==============================================================================
ruta_malla <- here(
    "data", "processed", "Malla", "NO2",
    sprintf("malla_spde_madrid_%s.rds", MALLA)
)
if (!file.exists(ruta_malla)) stop("No existe la malla: ", ruta_malla)

mesh <- readRDS(ruta_malla)
spde <- crear_spde(mesh, escala = "DIARIO")

indice_espacial <- inla.spde.make.index(
    name = "campo_espacial",
    n.spde = spde$n.spde
)

coordenadas <- as.matrix(df[, .(X_km, Y_km)])
A_espacial <- inla.spde.make.A(mesh = mesh, loc = coordenadas)

stopifnot(
    nrow(A_espacial) == nrow(df),
    ncol(A_espacial) == spde$n.spde
)


# ==============================================================================
# 4. STACK COMUN
# ==============================================================================
efectos_fijos <- data.frame(Intercept = rep(1, nrow(df)))
for (variable in COVARIABLES) {
    efectos_fijos[[variable]] <- df[[variable]]
}
efectos_fijos$Temperatura_rw2 <- df$Temperatura_rw2
efectos_fijos$Velocidad_Viento_rw2 <- df$Velocidad_Viento_rw2
efectos_fijos$Tipo_Urbana_fondo <- df$Tipo_Urbana_fondo
efectos_fijos$Tipo_Urbana_trafico <- df$Tipo_Urbana_trafico

# Un unico stack para todos los modelos, con el mismo campo espacial y
stack_modelo <- crear_stack_inla_spde(
    respuesta = df$y,
    A_espacial = A_espacial,
    indice_espacial = indice_espacial,
    efectos_fijos = efectos_fijos,
    tag = "estimacion"
)

indices_observaciones <- indices_stack_inla(stack_modelo, "estimacion")


# ==============================================================================
# 5. FUNCIONES COMUNES DE MODELIZACION
# ==============================================================================
# crear_stack_inla_spde(), ajustar_modelo_inla() y las funciones de extraccion
# se cargan desde R/modeling/inla_modeling.R. Aqui quedan solamente las
# decisiones particulares de este estudio: datos, formulas y configuracion.

# ==============================================================================
# 6. FORMULAS
# ==============================================================================
# Todas las formulas se escriben completas. Para cambiar un modelo, edita solo
# su formula. Las variables dentro de f(..., model = "rw2") son no lineales.
# Si anades una variable nueva, tambien debe existir en efectos_fijos (bloque 4).

formula_M0 <- y_response ~ -1 + Intercept +
    intensidad + Humedad_Relativa + Radiacion_Solar +
    Temperatura + Velocidad_Viento_sqrt +
    Presion_Barometrica +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M1 <- y_response ~ -1 + Intercept + Humedad_Relativa +
    Velocidad_Viento_sqrt +
    Presion_Barometrica + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    Llueve + intensidad +
    f(campo_espacial, model = spde)

formula_M2 <- y_response ~ -1 + Intercept +
    intensidad + Humedad_Relativa + f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
    Presion_Barometrica +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M3 <- y_response ~ -1 + Intercept +
    Temperatura + Velocidad_Viento_sqrt +
    Presion_Barometrica + intensidad +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M4 <- y_response ~ -1 + Intercept +
    f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    Velocidad_Viento_sqrt + intensidad +
    Presion_Barometrica +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M5 <- y_response ~ -1 + Intercept +
    f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
    Presion_Barometrica + intensidad +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M6 <- y_response ~ -1 + Intercept + intensidad +
    Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M7 <- y_response ~ -1 + Intercept + Humedad_Relativa + Radiacion_Solar +
    Temperatura + Velocidad_Viento_sqrt +
    Presion_Barometrica +
    Llueve +
    f(campo_espacial, model = spde)

# Lista de formulas. Solo sirve para recorrerlas y evitar repetir cuatro veces la
# llamada a inla(); las formulas que debes editar son las escritas arriba.
formulas_modelos <- list(
    M0 = formula_M0,
    M1 = formula_M1,
    M2 = formula_M2,
    M3 = formula_M3,
    M4 = formula_M4,
    M5 = formula_M5,
    M6 = formula_M6,
    M7 = formula_M7
)

descripcion_modelos <- c(
    M0 = "todo",
    M1 = "Modelo con HR",
    M2 = "Modelo HR pero con VV rw2",
    M3 = "Modelo con Temperatura",
    M4 = " Modelo con Temp y vv RW2",
    M5 = " Modelo con temp rw2 y vv rw2",
    M6 = "Modelo sin clima ",
    M7 = "Modelo sin intensidad"
)

catalogo_modelos <- data.table(
    Modelo = names(formulas_modelos),
    Descripcion = unname(descripcion_modelos[names(formulas_modelos)]),
    Formula = vapply(
        formulas_modelos,
        function(x) paste(deparse(x), collapse = " "),
        character(1)
    )
)
cat("\n--- MODELOS QUE SE VAN A AJUSTAR ---\n")
print(catalogo_modelos[, .(Modelo, Descripcion)], nrows = Inf)
if (interactive()) View(catalogo_modelos, title = "Formulas de los modelos")

ajustes <- setNames(
    vector("list", length(formulas_modelos)),
    names(formulas_modelos)
)
for (id_modelo in names(formulas_modelos)) {
    ajustes[[id_modelo]] <- ajustar_modelo_inla(
        formula_modelo = formulas_modelos[[id_modelo]],
        id_modelo = id_modelo,
        descripcion = descripcion_modelos[[id_modelo]],
        stack = stack_modelo,
        spde = spde,
        indices_observaciones = indices_observaciones,
        y_observada = df$y,
        nombres_rw2 = paste0(names(N_GRUPOS_RW2), "_rw2"),
        familia = FAMILIA,
        verbose = VERBOSE_INLA,
        num_threads = NUM_THREADS
    )
}


# ==============================================================================
# 8. TABLAS DE RESULTADOS Y APORTACION DE CADA VARIABLE
# ==============================================================================
modelos <- lapply(ajustes, `[[`, "modelo")
tabla_seleccion <- rbindlist(lapply(ajustes, `[[`, "metricas"))
tablas_coeficientes <- lapply(ajustes, `[[`, "coeficientes")
curvas_rw2_por_modelo <- lapply(ajustes, `[[`, "curvas_rw2")

# Diferencia de cada candidato respecto a M0. Un delta positivo significa que
# el candidato tiene un valor mayor (peor); un delta negativo favorece al
# candidato. El RMSE sigue siendo orientativo porque se calcula en muestra.
metricas_m0 <- tabla_seleccion[Modelo == "M0"]
tabla_seleccion[, `:=`(
    Delta_WAIC_vs_M0 = WAIC - metricas_m0$WAIC,
    Delta_DIC_vs_M0 = DIC - metricas_m0$DIC,
    Delta_RMSE_vs_M0 = RMSE_INLA_orientativo - metricas_m0$RMSE_INLA_orientativo
)]

tabla_comparacion_vs_M0 <- tabla_seleccion[
    Modelo != "M0",
    .(Modelo, Descripcion, Delta_WAIC_vs_M0, Delta_DIC_vs_M0, Delta_RMSE_vs_M0)
][order(-Delta_WAIC_vs_M0)]

fwrite(
    tabla_comparacion_vs_M0,
    file.path(DIR_OUT, "tabla_candidatos_vs_M0_2019_2021.csv")
)


# ==============================================================================
# 9. HOLD-OUT ESPACIAL: TRES ESTACIONES COMPLETAS
# ==============================================================================
tabla_holdout <- NULL
tabla_holdout_estacion <- NULL
predicciones_holdout <- NULL
modelos_holdout <- NULL

if (CALCULAR_HOLDOUT) {
    modelos_ausentes <- setdiff(MODELOS_HOLDOUT, names(formulas_modelos))
    if (length(modelos_ausentes) > 0L) {
        stop(
            "No existen estos modelos solicitados para el hold-out: ",
            paste(modelos_ausentes, collapse = ", ")
        )
    }

    estaciones_ausentes <- setdiff(ESTACIONES_HOLDOUT, unique(df$ESTACION))
    if (length(estaciones_ausentes) > 0L) {
        stop(
            "No existen estas estaciones del hold-out en los datos completos: ",
            paste(estaciones_ausentes, collapse = ", ")
        )
    }

    es_holdout <- df$ESTACION %in% ESTACIONES_HOLDOUT
    if (!any(es_holdout) || !any(!es_holdout)) {
        stop("El hold-out debe dejar observaciones tanto en prueba como en entrenamiento.")
    }

    # Se mantienen las filas de prueba en el stack para obtener su prediccion,
    # pero su respuesta se sustituye por NA y no interviene en el ajuste.
    y_holdout <- copy(df$y)
    y_holdout[es_holdout] <- NA_real_
    stack_holdout <- crear_stack_inla_spde(
        respuesta = y_holdout,
        A_espacial = A_espacial,
        indice_espacial = indice_espacial,
        efectos_fijos = efectos_fijos,
        tag = "holdout"
    )
    indices_stack_holdout <- indices_stack_inla(stack_holdout, "holdout")

    cat("\n--- HOLD-OUT ESPACIAL ---\n")
    cat("Estaciones reservadas:\n")
    print(unique(df[
        es_holdout,
        .(ESTACION, NOM_TIPO = as.character(NOM_TIPO))
    ])[order(ESTACION)])
    cat("Observaciones de entrenamiento: ", sum(!es_holdout), "\n", sep = "")
    cat("Observaciones de prueba: ", sum(es_holdout), "\n", sep = "")

    resultados_holdout <- setNames(
        lapply(MODELOS_HOLDOUT, function(id_modelo) {
            ajustar_holdout_espacial_inla(
                formula_modelo = formulas_modelos[[id_modelo]],
                id_modelo = id_modelo,
                descripcion = descripcion_modelos[[id_modelo]],
                stack = stack_holdout,
                spde = spde,
                indices_observaciones = indices_stack_holdout,
                datos = df,
                es_holdout = es_holdout,
                familia = FAMILIA,
                verbose = VERBOSE_INLA,
                num_threads = NUM_THREADS
            )
        }),
        MODELOS_HOLDOUT
    )

    modelos_holdout <- lapply(resultados_holdout, `[[`, "modelo")
    tabla_holdout <- rbindlist(lapply(resultados_holdout, `[[`, "global"))
    tabla_holdout_estacion <- rbindlist(lapply(
        resultados_holdout,
        `[[`,
        "por_estacion"
    ))
    predicciones_holdout <- rbindlist(lapply(
        resultados_holdout,
        `[[`,
        "predicciones"
    ))
    setorder(tabla_holdout_estacion, Modelo, ESTACION)
    setorder(predicciones_holdout, Modelo, ESTACION, FECHA)

    tabla_seleccion[
        tabla_holdout,
        on = .(Modelo, Descripcion),
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


# Una tabla sencilla por modelo: variable, coeficiente, IC95 y significacion.
for (id_modelo in names(tablas_coeficientes)) {
    fwrite(
        tablas_coeficientes[[id_modelo]],
        file.path(
            DIR_OUT,
            sprintf("tabla_coeficientes_%s_2019_2021.csv", id_modelo)
        )
    )
}
fwrite(
    tabla_seleccion,
    file.path(DIR_OUT, "tabla_comparacion_modelos_2019_2021.csv")
)
if (CALCULAR_HOLDOUT) {
    fwrite(
        tabla_holdout,
        file.path(DIR_OUT, "tabla_holdout_espacial_M5_2019_2021.csv")
    )
    fwrite(
        tabla_holdout_estacion,
        file.path(DIR_OUT, "tabla_holdout_por_estacion_M5_2019_2021.csv")
    )
    fwrite(
        predicciones_holdout,
        file.path(DIR_OUT, "predicciones_holdout_M5_2019_2021.csv")
    )
}


# ==============================================================================
# 10. RESULTADOS
# ==============================================================================
for (id_modelo in names(tablas_coeficientes)) {
    cat("\n--- Tabla de coeficientes ", id_modelo, " ---\n", sep = "")
    print(tablas_coeficientes[[id_modelo]])
}

cat("\n--- Comparacion de modelos ---\n")
print(tabla_seleccion)

cat("\n--- Modelos candidatos frente a M0 ---\n")
print(tabla_comparacion_vs_M0)
cat("Delta > 0: el candidato es peor que M0 en esa metrica.\n")
cat("Delta < 0: el candidato mejora a M0 en esa metrica.\n")

cat("El RMSE_INLA_orientativo es de ajuste y no de validacion cruzada.\n")
cat("El COV95_INLA_orientativo tampoco es validacion cruzada.\n")
cat("Los efectos RW2 se conservan en curvas_rw2_por_modelo y no se mezclan con los coeficientes lineales.\n")

if (CALCULAR_HOLDOUT) {
    cat("\n--- Validacion hold-out espacial de M5 ---\n")
    print(tabla_holdout)
    cat("\n--- Resultados hold-out por estacion y tipologia ---\n")
    print(tabla_holdout_estacion)
    cat("COV95 debe interpretarse junto con RMSE y Anchura_media_IC95.\n")

    if (interactive()) {
        View(tabla_holdout, title = "Hold-out espacial M5")
        View(
            tabla_holdout_estacion,
            title = "Hold-out M5 por estacion y tipologia"
        )
    }
}

cat("\nResultados guardados en: ", DIR_OUT, "\n", sep = "")
modelo_h <- modelos_holdout[["M5"]]

# Coeficientes fijos
View(modelo_h$summary.fixed)

# Efectos RW2
View(modelo_h$summary.random$Temperatura_rw2)
View(modelo_h$summary.random$Velocidad_Viento_rw2)

# Campo espacial
View(modelo_h$summary.random$campo_espacial)

# Rango, desviación espacial y precisión residual
View(modelo_h$summary.hyperpar)
View(modelos_holdout[["M5"]]$summary.hyperpar)
