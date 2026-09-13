# ==============================================================================
# SELECCION DE VARIABLES PARA EL MODELO INLA-SPDE (2019-2021)
# ==============================================================================
# Objetivo:
#   Comparar las catorce formulas M0-M13 escritas manualmente en el bloque 6.
#   M1-M6 eliminan terminos del modelo lineal completo M0; M9-M13 eliminan
#   terminos del modelo con efectos RW2 M8.
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
#   Se reserva un conjunto fijo de estaciones completas. Cada modelo se reajusta
#   una vez sin utilizar sus respuestas y las predicciones se evaluan solo en
#   esas estaciones (RMSE, COV95 y anchura del IC95). Es un hold-out espacial,
#   no una validacion LOSO.
# ==============================================================================

library(INLA)
library(data.table)
library(here)

source(here("R", "modeling", "spde_config.R"))
source(here("R", "utilities", "metricas_predictivas.R"))
source(here("R", "modeling", "inla_modeling.R"))
source(here("R", "utilities", "academic_quality_tables.R"))


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
    Velocidad_Viento = 40L,
    Humedad_Relativa = 40L
)
METODO_GRUPOS_RW2 <- "quantile"

MALLA <- "fina"
FAMILIA <- "gaussian"
NUM_THREADS <- 5L
VERBOSE_INLA <- TRUE

# Hold-out espacial. MODELOS_HOLDOUT = NULL lo aplica a todos los modelos; se
# puede restringir con, por ejemplo, c("M0", "M8"). Los nombres de estaciones
# deben coincidir exactamente con la columna ESTACION del dataset maestro.
CALCULAR_HOLDOUT <- TRUE
MODELOS_HOLDOUT <- NULL
ESTACIONES_HOLDOUT <- c(
    "Plaza Castilla", # Urbana trafico
    "Casa de Campo", # Suburbana
    "Ensanche Vallecas" # Urbana fondo
)

DIR_OUT <- here(
    "outputs", "modelo", "Modelo_1", "seleccion_variables", "2019_2021"
)
DIR_TABLAS_ACADEMICAS <- file.path(DIR_OUT, "academic_tables")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLAS_ACADEMICAS, recursive = TRUE, showWarnings = FALSE)



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
# Se crea <variable>_rw2 para cada variable de N_GRUPOS_RW2.
for (variable in names(N_GRUPOS_RW2)) {
    df[, (paste0(variable, "_rw2")) := inla.group(
        get(variable),
        n = N_GRUPOS_RW2[[variable]],
        method = METODO_GRUPOS_RW2
    )]
}
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
for (variable in names(N_GRUPOS_RW2)) {
    nombre_rw2 <- paste0(variable, "_rw2")
    efectos_fijos[[nombre_rw2]] <- df[[nombre_rw2]]
}
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
    intensidad + Humedad_Relativa + Velocidad_Viento_sqrt +
    Presion_Barometrica +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M1 <- y_response ~ -1 + Intercept + Velocidad_Viento_sqrt +
    Presion_Barometrica + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    Llueve + intensidad +
    f(campo_espacial, model = spde)

formula_M2 <- y_response ~ -1 + Intercept +
    intensidad + Humedad_Relativa +
    Presion_Barometrica +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M3 <- y_response ~ -1 + Intercept +
    Humedad_Relativa + Velocidad_Viento_sqrt + intensidad +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M4 <- y_response ~ -1 + Intercept +
    Humedad_Relativa +
    Velocidad_Viento_sqrt + intensidad +
    Presion_Barometrica +
    Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M5 <- y_response ~ -1 + Intercept +
    Humedad_Relativa +
    Velocidad_Viento_sqrt +
    Presion_Barometrica + intensidad +
    Llueve +
    f(campo_espacial, model = spde)

formula_M6 <- y_response ~ -1 + Intercept + intensidad +
    Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M7 <- y_response ~ -1 + Intercept + Humedad_Relativa + Radiacion_Solar +
    Temperatura + Velocidad_Viento_sqrt +
    Presion_Barometrica +
    Llueve +
    f(campo_espacial, model = spde)

# M8: modelo con temperatura y viento RW2. M9-M13 eliminan terminos de M8.
formula_M8_M5 <- y_response ~ -1 + Intercept + f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
    Presion_Barometrica + intensidad +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M9_M5 <- y_response ~ -1 + Intercept +
    f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    Presion_Barometrica + intensidad +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M10_M5 <- y_response ~ -1 + Intercept +
    f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + intensidad +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M11_M5 <- y_response ~ -1 + Intercept +
    f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + Presion_Barometrica +
    Llueve + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M12_M5 <- y_response ~ -1 + Intercept +
    f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + Presion_Barometrica +
    intensidad + Tipo_Urbana_fondo + Tipo_Urbana_trafico +
    f(campo_espacial, model = spde)

formula_M13_M5 <- y_response ~ -1 + Intercept +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + Llueve +
    intensidad +
    f(campo_espacial, model = spde)

# Lista de formulas. Solo sirve para recorrerlas y evitar repetir la
# llamada a inla(); las formulas que debes editar son las escritas arriba.
formulas_modelos <- list(
    M0 = formula_M0,
    M1 = formula_M1,
    M2 = formula_M2,
    M3 = formula_M3,
    M4 = formula_M4,
    M5 = formula_M5,
    M6 = formula_M6,
    M7 = formula_M7,
    M8 = formula_M8_M5,
    M9 = formula_M9_M5,
    M10 = formula_M10_M5,
    M11 = formula_M11_M5,
    M12 = formula_M12_M5,
    M13 = formula_M13_M5
)

descripcion_modelos <- c(
    M0 = "Full linear model: humidity, sqrt wind, pressure, rain, traffic, station type",
    M1 = "M0 without relative humidity",
    M2 = "M0 without square-root wind speed",
    M3 = "M0 without barometric pressure",
    M4 = "M0 without rain indicator",
    M5 = "M0 without station type",
    M6 = "Traffic intensity and station type only",
    M7 = "Meteorological covariates only (with temperature and solar radiation)",
    M8 = "RW2 temperature and wind speed, pressure, rain, traffic, station type",
    M9 = "M8 without wind speed (RW2)",
    M10 = "M8 without barometric pressure",
    M11 = "M8 without traffic intensity",
    M12 = "M8 without rain indicator",
    M13 = "M8 without temperature (RW2), barometric pressure and station type"
)

# Comparaciones de aportacion: cada modelo reducido frente a su referencia.
# M7 no es una eliminacion de M0 (anade temperatura y radiacion), no se incluye.
referencias_aportacion <- data.table(
    Modelo = c("M1", "M2", "M3", "M4", "M5", "M6", "M9", "M10", "M11", "M12", "M13"),
    Referencia = c(rep("M0", 6L), rep("M8", 5L)),
    Termino_eliminado = c(
        "Relative humidity", "Square-root wind speed", "Barometric pressure",
        "Rain indicator", "Station type", "All meteorological covariates",
        "Wind speed (RW2)", "Barometric pressure", "Traffic intensity",
        "Rain indicator", "Temperature (RW2), barometric pressure, station type"
    )
)
stopifnot(all(c(referencias_aportacion$Modelo, referencias_aportacion$Referencia) %in%
    names(formulas_modelos)))

stopifnot(
    identical(names(formulas_modelos), names(descripcion_modelos)),
    !anyNA(descripcion_modelos),
    all(nzchar(descripcion_modelos))
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

# Catalogo en ingles: puede generarse antes de ejecutar los ajustes.
fwrite(catalogo_modelos, file.path(DIR_OUT, "catalogo_modelos_2019_2021.csv"))
model_catalog_english <- catalogo_modelos[, .(Model = Modelo, Description = Descripcion)]
fwrite(model_catalog_english, file.path(DIR_TABLAS_ACADEMICAS, "model_catalog_2019_2021.csv"))
path_model_catalog <- booktabs_png(
    model_catalog_english,
    file.path(DIR_TABLAS_ACADEMICAS, "table_model_catalog_2019_2021.png"),
    title = "Candidate models for variable selection",
    subtitle = "M0-M13; daily NO2 observations, Madrid, 2019-2021",
    note = paste0(
        "All models include an intercept and the same spatial field. M8-M13 remove ",
        "one term from M5; station type is removed as a two-indicator block. ",
        "RW2 denotes a smooth covariate effect."
    ),
    widths = c(0.6, 4.4), align = c("center", "left"),
    font_size = 8.2, row_height = 0.32
)

# ==============================================================================
# 7. AJUSTE DE LOS CATORCE MODELOS
# ==============================================================================
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
tabla_seleccion[, `:=`(
    p_eff_WAIC = vapply(
        modelos[Modelo],
        function(modelo) modelo$waic$p.eff,
        numeric(1)
    ),
    p_eff_DIC = vapply(
        modelos[Modelo],
        function(modelo) modelo$dic$p.eff,
        numeric(1)
    )
)]

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

# Aportacion de cada termino: modelo reducido menos M5, no menos M0.
# Solo M5 tiene hold-out; estos deltas de RMSE son de ajuste en muestra.
terminos_eliminados <- c(
    M8 = "Temperature (RW2)", M9 = "Wind speed (RW2)",
    M10 = "Barometric pressure", M11 = "Traffic intensity",
    M12 = "Rain indicator", M13 = "Station type"
)
metricas_m5 <- tabla_seleccion[Modelo == "M5"]
tabla_aportacion_vs_M5 <- tabla_seleccion[
    match(names(terminos_eliminados), Modelo),
    .(
        Modelo,
        Termino_eliminado = unname(terminos_eliminados[Modelo]),
        Delta_WAIC_vs_M5 = WAIC - metricas_m5$WAIC,
        Delta_DIC_vs_M5 = DIC - metricas_m5$DIC,
        Delta_RMSE_vs_M5 = RMSE_INLA_orientativo - metricas_m5$RMSE_INLA_orientativo
    )
]
fwrite(
    tabla_aportacion_vs_M5,
    file.path(DIR_OUT, "tabla_aportacion_variables_vs_M5_2019_2021.csv")
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
# 10. ACADEMIC TABLES IN ENGLISH (CSV + PNG)
# ==============================================================================
VARIABLE_LABELS <- c(
    intensidad = "Traffic intensity",
    Temperatura = "Temperature",
    Velocidad_Viento = "Wind speed",
    Velocidad_Viento_sqrt = "Square-root wind speed",
    Presion_Barometrica = "Barometric pressure",
    Humedad_Relativa = "Relative humidity",
    Precipitaciones = "Precipitation",
    Llueve = "Rain indicator",
    Radiacion_Solar = "Solar radiation",
    Tipo_Urbana_fondo = "Station type: urban background",
    Tipo_Urbana_trafico = "Station type: urban traffic",
    Intercept = "Intercept"
)

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

format_number <- function(x, digits = 2L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "--")
}

# Overall model comparison.
model_comparison_english <- tabla_seleccion[, .(
    Model = Modelo,
    Description = Descripcion,
    WAIC = format_number(WAIC, 1L),
    `Delta WAIC` = format_number(Delta_WAIC_vs_M0, 1L),
    DIC = format_number(DIC, 1L),
    `Delta DIC` = format_number(Delta_DIC_vs_M0, 1L),
    `WAIC p_eff` = format_number(p_eff_WAIC, 1L),
    `Fitted RMSE` = format_number(RMSE_INLA_orientativo, 3L)
)]
fwrite(
    model_comparison_english,
    file.path(DIR_TABLAS_ACADEMICAS, "model_comparison_2019_2021.csv")
)
path_model_comparison <- booktabs_png(
    model_comparison_english,
    file.path(
        DIR_TABLAS_ACADEMICAS,
        "table_candidate_comparison_2019_2021.png"
    ),
    title = "Comparison of candidate INLA-SPDE models",
    subtitle = "M0-M13; daily NO2 models for Madrid, 2019-2021",
    note = paste0(
        "Delta values are calculated relative to M0; negative values favour the ",
        "candidate model. Lower WAIC, DIC and RMSE values are preferred. Fitted RMSE ",
        "is an in-sample diagnostic, not a spatial validation metric."
    ),
    widths = c(0.55, 3.15, 0.78, 0.88, 0.78, 0.88, 0.85, 0.90),
    align = c("center", "left", rep("right", 6)),
    font_size = 7.6,
    row_height = 0.32
)

# Removal of individual terms from M5 (M8-M13).
variable_contribution_english <- tabla_aportacion_vs_M5[, .(
    Model = Modelo,
    `Removed term` = Termino_eliminado,
    `Delta WAIC` = format_number(Delta_WAIC_vs_M5, 1L),
    `Delta DIC` = format_number(Delta_DIC_vs_M5, 1L),
    `Delta fitted RMSE` = format_number(Delta_RMSE_vs_M5, 3L)
)]
fwrite(
    variable_contribution_english,
    file.path(DIR_TABLAS_ACADEMICAS, "variable_contribution_vs_M5_2019_2021.csv")
)
path_variable_contribution <- booktabs_png(
    variable_contribution_english,
    file.path(DIR_TABLAS_ACADEMICAS, "table_variable_contribution_vs_M5_2019_2021.png"),
    title = "Variable contribution relative to model M5",
    subtitle = "M8-M13: one term removed at a time",
    note = paste0(
        "Delta = reduced model minus M5. Positive values indicate deterioration ",
        "after removal; negative values favour removal. Fitted RMSE is in-sample. ",
        "Spatial hold-out is performed for M5 only. These comparisons are not ",
        "tests of statistical significance."
    ),
    widths = c(0.55, 2.05, 0.90, 0.90, 1.30),
    align = c("center", "left", rep("right", 3)),
    font_size = 8.0, row_height = 0.34
)

# One fixed-effect table per model. RW2 effects are not assigned a coefficient.
coefficient_tables_english <- setNames(lapply(names(tablas_coeficientes), function(id) {
    table_model <- copy(tablas_coeficientes[[id]])
    table_model[, Variable := translate_variable(Variable)]
    table_model[, Significativa_95 := fifelse(
        Significativa_95 == "Si",
        "Yes",
        "No"
    )]
    setnames(
        table_model,
        c("Variable", "Coeficiente", "IC95", "Significativa_95"),
        c("Variable", "Posterior mean", "95% credible interval", "Significant (95%)")
    )
    table_model[, `Posterior mean` := format_number(`Posterior mean`, 4L)]
    table_model
}), names(tablas_coeficientes))

coefficient_table_paths <- setNames(vapply(
    names(coefficient_tables_english),
    function(id) {
        table_model <- coefficient_tables_english[[id]]
        fwrite(
            table_model,
            file.path(
                DIR_TABLAS_ACADEMICAS,
                sprintf("fixed_effects_%s_2019_2021.csv", id)
            )
        )
        booktabs_png(
            table_model,
            file.path(
                DIR_TABLAS_ACADEMICAS,
                sprintf("table_fixed_effects_%s_2019_2021.png", id)
            ),
            title = sprintf("Fixed effects for model %s", id),
            subtitle = descripcion_modelos[[id]],
            note = paste0(
                "A fixed effect is labelled significant when its 95% posterior credible ",
                "interval excludes zero. ",
                if (any(c("Tipo_Urbana_fondo", "Tipo_Urbana_trafico") %in%
                    all.vars(formulas_modelos[[id]]))) {
                    paste0(
                        "Station-type effects are differences from the ",
                        "suburban reference category. "
                    )
                } else {
                    ""
                },
                "RW2 effects have no single coefficient."
            ),
            widths = c(2.55, 1.15, 1.85, 1.25),
            align = c("left", "right", "center", "center"),
            font_size = 8.2,
            row_height = 0.31
        )
    },
    character(1)
), names(coefficient_tables_english))

# Spatial hold-out tables for M5.
path_holdout_global <- NULL
path_holdout_station <- NULL
if (CALCULAR_HOLDOUT) {
    holdout_global_english <- tabla_holdout[, .(
        Model = Modelo,
        Description = Descripcion,
        RMSE = format_number(RMSE_HOLDOUT, 3L),
        MAE = format_number(MAE_HOLDOUT, 3L),
        Bias = format_number(Sesgo_HOLDOUT, 3L),
        `COV95 (%)` = format_number(COV95_HOLDOUT, 1L),
        `Mean 95% width` = format_number(Anchura_media_IC95_HOLDOUT, 3L),
        Observations = as.character(Observaciones_HOLDOUT)
    )]
    holdout_station_english <- tabla_holdout_estacion[, .(
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
        holdout_global_english,
        file.path(DIR_TABLAS_ACADEMICAS, "spatial_holdout_M5_2019_2021.csv")
    )
    fwrite(
        holdout_station_english,
        file.path(
            DIR_TABLAS_ACADEMICAS,
            "spatial_holdout_by_station_M5_2019_2021.csv"
        )
    )
    path_holdout_global <- booktabs_png(
        holdout_global_english,
        file.path(
            DIR_TABLAS_ACADEMICAS,
            "table_spatial_holdout_M5_2019_2021.png"
        ),
        title = "Spatial hold-out performance of model M5",
        subtitle = paste(
            "Stations excluded jointly:",
            paste(ESTACIONES_HOLDOUT, collapse = ", ")
        ),
        note = paste0(
            "Metrics are calculated only from excluded-station observations. Predictive ",
            "COV95 includes observation variance and must be interpreted together with ",
            "RMSE and mean 95% predictive-interval width."
        ),
        widths = c(0.55, 2.75, 0.72, 0.70, 0.70, 0.78, 1.00, 0.90),
        align = c("center", "left", rep("right", 6)),
        font_size = 8.0,
        row_height = 0.32
    )
    path_holdout_station <- booktabs_png(
        holdout_station_english,
        file.path(
            DIR_TABLAS_ACADEMICAS,
            "table_spatial_holdout_by_station_M5_2019_2021.png"
        ),
        title = "Spatial hold-out performance by station",
        subtitle = "Model M5; daily NO2 observations, 2019-2021",
        note = paste0(
            "The three stations are excluded simultaneously during fitting. Coverage ",
            "must be assessed together with prediction error and interval width."
        ),
        widths = c(0.55, 1.45, 1.30, 0.82, 0.72, 0.68, 0.68, 0.78, 0.96),
        align = c("center", "left", "left", rep("right", 6)),
        font_size = 7.7,
        row_height = 0.34
    )
}

expected_png <- c(
    path_model_catalog,
    path_model_comparison,
    path_variable_contribution,
    unname(coefficient_table_paths),
    path_holdout_global,
    path_holdout_station
)
expected_png <- expected_png[!is.na(expected_png) & nzchar(expected_png)]
if (!all(file.exists(expected_png))) {
    stop("One or more academic PNG tables could not be generated.")
}
cat(
    "\nAcademic PNG tables generated: ",
    length(expected_png),
    "\nDirectory: ",
    DIR_TABLAS_ACADEMICAS,
    "\n",
    sep = ""
)


# ==============================================================================
# 11. RESULTADOS
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
