# ==============================================================================
# SELECCION DE VARIABLES PARA EL MODELO INLA-SPDE (2019-2021)
# ==============================================================================
# Objetivo:
#   Comparar un modelo completo con modelos en los que se quita una covariable
#   cada vez. Todos los modelos usan los mismos datos y el mismo campo espacial.
#
# Metricas principales:
#   - WAIC y DIC: cuanto menor, mejor.
#   - RMSE_INLA_orientativo: calculado con la media posterior ajustada por INLA.
#     Es rapido, pero no es un RMSE de validacion cruzada.
#   - COV95_INLA_orientativo: porcentaje observado dentro del intervalo formado
#     con mean +/- 1.96 * sd de summary.fitted.values, como en simulacion.R.
#     Es una medida de ajuste orientativa, no una cobertura predictiva LOSO.
#
# Validacion opcional:
#   Si CALCULAR_LOSO = TRUE se usa inla.group.cv(), la funcion incluida en INLA,
#   para omitir una estacion completa. Se calculan RMSE_LOSO y COV95_LOSO.
#   No se programan manualmente los folds.
# ==============================================================================

library(INLA)
library(data.table)
library(here)

source(here("R", "modeling", "spde_config.R"))


# ==============================================================================
# 1. CONFIGURACION
# ==============================================================================
FECHA_INICIO <- as.Date("2019-01-01")
FECHA_FIN <- as.Date("2021-12-31")
RESPUESTA <- "LOG_NO2_DIARIO"

# Un unico maestro diario que ya contiene todo el periodo 2019-2021.
# Paso_3_union_maestros_diarios.R guarda los maestros conjuntos en esta carpeta.
RUTA_DATOS <- here(
    "data", "processed", "Maestro", "diario",
    "dataset_maestro_inla_20190101_20211231_DIARIO2.rds"
)

COVARIABLES <- c(
    "intensidad",
    "Temperatura",
    "Velocidad_Viento",
    "Presion_Barometrica",
    "Humedad_Relativa",
    "Radiacion_Solar",
    "Precipitaciones",
    "Llueve"
)

# Numero de grupos para las variables que entran como RW2. Estos valores se
# pueden sustituir despues por los elegidos en el analisis de sensibilidad.
N_GRUPOS_RW2 <- c(
    Temperatura = 40L,
    Velocidad_Viento = 40L,
    Radiacion_Solar = 40L
)
METODO_GRUPOS_RW2 <- "quantile"

MALLA <- "media"
FAMILIA <- "gaussian"
NUM_THREADS <- 5L
VERBOSE_INLA <- FALSE

# Esta es la unica opcion de validacion cruzada del archivo.
CALCULAR_LOSO <- FALSE

DIR_OUT <- here(
    "outputs", "modelo", "Modelo_1", "seleccion_variables", "2019_2021"
)
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. CARGA Y PREPARACION DE LOS DATOS
# ==============================================================================
if (!file.exists(RUTA_DATOS)) {
    stop("No se encuentra el maestro conjunto 2019-2021: ", RUTA_DATOS)
}

df <- as.data.table(readRDS(RUTA_DATOS))

# Indicador binario calculado con la precipitacion original en milimetros.
# No se usa Precipitaciones porque esa columna ya contiene el z-score.
if (!"Precipitaciones_raw" %in% names(df)) {
    stop("Falta Precipitaciones_raw; no se puede aplicar el umbral de 1 mm.")
}
df[, Llueve := fifelse(
    is.na(Precipitaciones_raw),
    NA_integer_,
    as.integer(Precipitaciones_raw >= 1)
)]

columnas_necesarias <- c(
    "ESTACION", "FECHA", "X_km", "Y_km", RESPUESTA, COVARIABLES
)
# ==============================================================================
df[, FECHA := as.Date(FECHA)] # Asegurarse de que FECHA es de clase Date
df <- df[FECHA >= FECHA_INICIO & FECHA <= FECHA_FIN]
if (!nrow(df)) {
    stop("El maestro no contiene observaciones dentro del periodo configurado.")
}
setorder(df, FECHA, ESTACION)

# La respuesta ya viene transformada: LOG_NO2_DIARIO es log1p(DATO_DIARIO).
df[, y := get(RESPUESTA)]

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
df[, Radiacion_Solar_rw2 := inla.group(
    Radiacion_Solar,
    n = N_GRUPOS_RW2[["Radiacion_Solar"]],
    method = METODO_GRUPOS_RW2
)]

cat("\nDatos utilizados\n")
cat("Periodo:", format(FECHA_INICIO), "a", format(FECHA_FIN), "\n")
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
efectos_fijos$Radiacion_Solar_rw2 <- df$Radiacion_Solar_rw2

# Un unico stack para todos los modelos, con el mismo campo espacial y
stack_modelo <- inla.stack(
    data = list(y_response = df$y),
    A = list(A_espacial, 1),
    effects = list(indice_espacial, efectos_fijos), # Efecto fijo es el intercepto
    tag = "estimacion"
)

indices_observaciones <- inla.stack.index(
    stack_modelo,
    tag = "estimacion"
)$data


# ==============================================================================
# 5. MODELOS PARA LA SELECCION DE VARIABLES
# ==============================================================================
# Esta funcion contiene las opciones comunes de INLA. La formula, el nombre y
# la descripcion se indican expresamente al ajustar cada modelo.
ajustar_modelo <- function(formula_modelo, id_modelo, descripcion) {
    cat("\nAjustando ", id_modelo, ": ", descripcion, "...\n", sep = "")

    modelo <- inla(
        formula = formula_modelo,
        data = inla.stack.data(stack_modelo, spde = spde),
        family = FAMILIA,
        verbose = VERBOSE_INLA,
        num.threads = NUM_THREADS,
        control.predictor = list(
            A = inla.stack.A(stack_modelo),
            compute = TRUE
        ),
        control.inla = list(
            strategy = "gaussian",
            int.strategy = "eb"
        ),
        control.compute = list(
            dic = TRUE,
            waic = TRUE,
            cpo = FALSE,
            openmp.strategy = "huge"
        )
    )

    resumen_ajuste <- modelo$summary.fitted.values[
        indices_observaciones, ,
        drop = FALSE
    ]
    prediccion <- resumen_ajuste[, "mean"]
    rmse_inla <- sqrt(mean((prediccion - df$y)^2, na.rm = TRUE))

    # Mismo calculo orientativo que en simulacion.R: intervalo centrado en la
    # media ajustada y construido con la desviacion posterior devuelta por
    # summary.fitted.values. Aqui no se suma la varianza residual gaussiana.
    sd_ajuste <- resumen_ajuste[, "sd"]
    observaciones_validas <- is.finite(df$y) &
        is.finite(prediccion) &
        is.finite(sd_ajuste)
    dentro_intervalo <- df$y[observaciones_validas] >=
        prediccion[observaciones_validas] -
            1.96 * sd_ajuste[observaciones_validas] &
        df$y[observaciones_validas] <=
            prediccion[observaciones_validas] +
                1.96 * sd_ajuste[observaciones_validas]
    cov95_inla <- 100 * mean(dentro_intervalo)

    # Efectos lineales: un coeficiente por variable.
    resumen_fijos <- modelo$summary.fixed
    tabla_fijos <- data.table(
        Variable = rownames(resumen_fijos),
        Coeficiente = resumen_fijos[["mean"]],
        IC95 = sprintf(
            "[%.4f, %.4f]",
            resumen_fijos[["0.025quant"]],
            resumen_fijos[["0.975quant"]]
        ),
        Significativa_95 = fifelse(
            resumen_fijos[["0.025quant"]] > 0 |
                resumen_fijos[["0.975quant"]] < 0,
            "Si",
            "No"
        )
    )
    tabla_fijos <- tabla_fijos[
        Variable != "Intercept",
        .(Variable, Coeficiente, IC95, Significativa_95)
    ]

    # Efectos RW2: no existe un coeficiente unico. Se guarda una fila por nivel
    # de la curva y el criterio de exclusion de cero es solamente puntual.
    posibles_rw2 <- paste0(names(N_GRUPOS_RW2), "_rw2")
    nombres_rw2_modelo <- intersect(
        posibles_rw2,
        names(modelo$summary.random)
    )
    tabla_rw2 <- rbindlist(lapply(nombres_rw2_modelo, function(nombre_rw2) {
        resumen_rw2 <- modelo$summary.random[[nombre_rw2]]
        resultado <- data.table(
            Modelo = id_modelo,
            Descripcion = descripcion,
            Variable = sub("_rw2$", "", nombre_rw2),
            Tipo_efecto = "RW2",
            Nivel_RW2 = resumen_rw2[["ID"]],
            Coeficiente = resumen_rw2[["mean"]],
            Q025 = resumen_rw2[["0.025quant"]],
            Q975 = resumen_rw2[["0.975quant"]]
        )
        resultado[, IC95 := sprintf("[%.4f, %.4f]", Q025, Q975)]
        resultado[, Significativa_95 := fifelse(
            Q025 > 0 | Q975 < 0,
            "Si (puntual)",
            "No (puntual)"
        )]
        resultado
    }), use.names = TRUE, fill = TRUE)

    metricas <- data.table(
        Modelo = id_modelo,
        Descripcion = descripcion,
        WAIC = modelo$waic$waic,
        DIC = modelo$dic$dic,
        RMSE_INLA_orientativo = rmse_inla,
        COV95_INLA_orientativo = cov95_inla
    )

    list(
        modelo = modelo,
        metricas = metricas,
        coeficientes = tabla_fijos,
        curvas_rw2 = tabla_rw2
    )
}

# ------------------------------------------------------------------------------
# MODELOS E
# ------------------------------------------------------------------------------
# M0: referencia con todas las covariables.
formula_M0 <- y_response ~ -1 + Intercept + intensidad + f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
    Presion_Barometrica +
    Precipitaciones + Humedad_Relativa + f(Radiacion_Solar_rw2, model = "rw2", scale.model = "True") +
    f(campo_espacial, model = spde)

ajuste_M0 <- ajustar_modelo(
    formula_modelo = formula_M0,
    id_modelo = "M0",
    descripcion = "Con covariables: todas"
)

# M1: modelo completo inicial.
formula_M1 <- y_response ~ -1 + Intercept + intensidad + f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + Presion_Barometrica + Humedad_Relativa + Precipitaciones + f(campo_espacial, model = spde)

ajuste_M1 <- ajustar_modelo(
    formula_modelo = formula_M1,
    id_modelo = "M1",
    descripcion = "Modelo sin colinealidad"
)

# Escribe los modelos siguientes de la misma forma. Por ejemplo:
## M2: modelo sin Radiacion_Solar
formula_M2 <- y_response ~ -1 + Intercept +
    intensidad +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
    Presion_Barometrica + Precipitaciones +
    f(campo_espacial, model = spde)

ajuste_M2 <- ajustar_modelo(
    formula_modelo = formula_M2,
    id_modelo = "M2",
    descripcion = " sin Humedad_Relativa"
)

formula_M3 <- y_response ~ -1 + Intercept +
    intensidad + Humedad_Relativa +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + Precipitaciones +
    f(campo_espacial, model = spde)

ajuste_M3 <- ajustar_modelo(
    formula_modelo = formula_M3,
    id_modelo = "M3",
    descripcion = "Sin PB"
)
formula_M4 <- y_response ~ -1 + Intercept +
    Presion_Barometrica + Humedad_Relativa +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + Precipitaciones +
    f(campo_espacial, model = spde)
ajuste_M4 <- ajustar_modelo(
    formula_modelo = formula_M4,
    id_modelo = "M4",
    descripcion = "Sin intensidad"
)

formula_M5 <- y_response ~ -1 + Intercept +
    Presion_Barometrica + Humedad_Relativa +
    f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) + intensidad +
    f(campo_espacial, model = spde)

ajuste_M5 <- ajustar_modelo(
    formula_modelo = formula_M5,
    id_modelo = "M5",
    descripcion = "Sin Precipitaciones"
)

formula_M6 <- y_response ~ -1 + Intercept +
    Presion_Barometrica + Humedad_Relativa + intensidad + Precipitaciones +
    f(campo_espacial, model = spde)
ajuste_M6 <- ajustar_modelo(
    formula_modelo = formula_M6,
    id_modelo = "M5",
    descripcion = "Sin VV"
)
formula_M7 <- y_response ~ -1 + Intercept + f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
    Presion_Barometrica + f(Temperatura_rw2, model = "rw2", scale.model = TRUE) + intensidad + Precipitaciones +
    f(campo_espacial, model = spde)
ajuste_M7 <- ajustar_modelo(
    formula_modelo = formula_M7,
    id_modelo = "M5",
    descripcion = "CON TEMP"
)


# Se presentan exactamente cuatro tablas de coeficientes, una por modelo.
ajustes <- list(
    M0 = ajuste_M0,
    M1 = ajuste_M1,
    M2 = ajuste_M2,
    M3 = ajuste_M3,
    M4 = ajuste_M4,
    M5 = ajuste_M5,
    M6 = ajuste_M6,
    M7 = ajuste_M7
)

modelos <- lapply(ajustes, `[[`, "modelo")
tabla_seleccion <- rbindlist(lapply(ajustes, `[[`, "metricas"))
tablas_coeficientes <- lapply(ajustes, `[[`, "coeficientes")
curvas_rw2_por_modelo <- lapply(ajustes, `[[`, "curvas_rw2")


# ==============================================================================
# 7. VALIDACION LOSO OPCIONAL CON LA FUNCION DE INLA
# ==============================================================================
tabla_loso <- NULL

if (CALCULAR_LOSO) {
    filas_por_estacion <- split(seq_len(nrow(df)), df$ESTACION)
    grupos_loso <- unname(filas_por_estacion[as.character(df$ESTACION)])
    maximo_grupo <- max(lengths(filas_por_estacion))

    tabla_loso <- rbindlist(lapply(names(modelos), function(id_modelo) {
        cat("Validacion LOSO de ", id_modelo, " con inla.group.cv()...\n", sep = "")

        resultado_cv <- inla.group.cv(
            result = modelos[[id_modelo]],
            groups = grupos_loso,
            size.max = maximo_grupo,
            verbose = TRUE
        )

        prediccion_cv <- resultado_cv$mean
        sd_predictiva_cv <- resultado_cv$sd
        observaciones_validas <- is.finite(prediccion_cv) &
            is.finite(sd_predictiva_cv) &
            is.finite(df$y)

        if (!any(observaciones_validas)) {
            stop("INLA no ha devuelto predicciones LOSO validas para ", id_modelo)
        }

        error_cv <- prediccion_cv[observaciones_validas] -
            df$y[observaciones_validas]
        limite_inferior <- prediccion_cv[observaciones_validas] -
            1.96 * sd_predictiva_cv[observaciones_validas]
        limite_superior <- prediccion_cv[observaciones_validas] +
            1.96 * sd_predictiva_cv[observaciones_validas]
        dentro_intervalo <- df$y[observaciones_validas] >= limite_inferior &
            df$y[observaciones_validas] <= limite_superior

        data.table(
            Modelo = id_modelo,
            Descripcion = ajustes[[id_modelo]]$metricas$Descripcion,
            RMSE_LOSO = sqrt(mean(error_cv^2)),
            COV95_LOSO = 100 * mean(dentro_intervalo),
            Observaciones_LOSO = sum(observaciones_validas)
        )
    }))

    tabla_seleccion[
        tabla_loso,
        on = .(Modelo, Descripcion),
        `:=`(
            RMSE_LOSO = i.RMSE_LOSO,
            COV95_LOSO = i.COV95_LOSO,
            Observaciones_LOSO = i.Observaciones_LOSO
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


# ==============================================================================
# 8. RESULTADOS
# ==============================================================================
for (id_modelo in names(tablas_coeficientes)) {
    cat("\n--- Tabla de coeficientes ", id_modelo, " ---\n", sep = "")
    print(tablas_coeficientes[[id_modelo]])
}

cat("\n--- Comparacion de modelos ---\n")
print(tabla_seleccion)
cat("El RMSE_INLA_orientativo es de ajuste y no de validacion cruzada.\n")
cat("El COV95_INLA_orientativo tampoco es validacion cruzada.\n")
cat("Los efectos RW2 se conservan en curvas_rw2_por_modelo y no se mezclan con los coeficientes lineales.\n")

if (CALCULAR_LOSO) {
    cat("\n--- RMSE y COV95 de validacion LOSO ---\n")
    print(tabla_loso)
    cat("COV95_LOSO cercano a 95 indica una cobertura predictiva adecuada.\n")
}

cat("\nResultados guardados en: ", DIR_OUT, "\n", sep = "")
