# ==============================================================================
# COMPARACION DE EFECTOS LINEALES Y RW2 CON SPDE - 2025
# ==============================================================================
# Objetivo:
#   Partir de un modelo INLA-SPDE con todas las covariables lineales y comparar,
#   para temperatura, radiacion solar y velocidad del viento, tres formas:
#   lineal original, transformacion parametrica y RW2. Las restantes
#   covariables permanecen lineales y el campo espacial no cambia.
#
# Comparacion principal:
#   Una unica tabla con WAIC, DIC y RMSE posterior de ajuste. Este RMSE se
#   obtiene del ajuste INLA ya calculado y se usa como cribado rapido.
#
# Diagnostico visual:
#   Curvas del efecto lineal, de la transformacion parametrica y del RW2,
#   con sus intervalos creibles posteriores del 95%.
#
# Validacion espacial opcional:ss
#   Si CALCULAR_LOSO = TRUE, se elimina una estacion completa cada vez para
#   los modelos indicados en MODELOS_LOSO. Sus resultados se guardan aparte.
# ==============================================================================


library(INLA)
library(data.table)
library(ggplot2)
library(here)
source(here("R", "modeling", "spde_config.R"))



# ==============================================================================
# 0. CONFIGURACION EDITABLE
# ==============================================================================
ANIO <- 2025
ESCALA <- "DIARIO"

RESPUESTA <- "LOG_NO2_DIARIO"
FECHA_INICIO <- as.Date("2025-01-01")
FECHA_FIN <- as.Date("2025-12-31")

# Introducimos las variables.
COVARIABLES <- c(
    "Temperatura",
    "Velocidad_Viento",
    "Radiacion_Solar",
    "Humedad_Relativa",
    "Presion_Barometrica",
    "intensidad",
    "Precipitaciones"
)

# Para cada variable se compara la forma lineal original, una transformacion
# parametrica (ya estandarizada en el maestro) y un efecto flexible RW2.
TRANSFORMACIONES <- c(
    Temperatura = "Temperatura_log",
    Radiacion_Solar = "Radiacion_Solar_log",
    Velocidad_Viento = "Velocidad_Viento_sqrt"
)
ETIQUETAS_TRANSFORMACIONES <- c(
    Temperatura = "Logaritmica",
    Radiacion_Solar = "Logaritmica",
    Velocidad_Viento = "Raiz cuadrada"
)

# Numero de grupos de los efectos RW2.
N_GRUPOS_RW2 <- c(
    Temperatura = 40L,
    Velocidad_Viento = 40L,
    Radiacion_Solar = 40L
)
ETIQUETAS_RW2 <- c(
    Temperatura = "Temperatura ",
    Velocidad_Viento = "Velocidad del viento ",
    Radiacion_Solar = "Radiación solar ",
    Presion_Barometrica = "Presion barometrica",
    intensidad = "Intensidad de trafico"
)
METODO_GRUPOS_RW2 <- "quantile" # Método utilizado para agrupar los valores de las covariables en N_GRUPOS_RW2. Puede ser "quantile" o "cut".

FAMILIA <- "gaussian"
MALLA <- "media" # "gruesa" | "media" | "fina"

CALCULAR_CPO <- FALSE # No se necesita para la tabla WAIC-DIC-RMSE.
CALCULAR_LOSO <- FALSE # Cambiar a TRUE solo para validar los candidatos finales.
# Valida solo los candidatos finales indicados despues del cribado rapido.
# Ejemplo: MODELOS_LOSO <- c("M1", "M2")
MODELOS_LOSO <- character()
GUARDAR_MODELOS <- TRUE
VERBOSE_INLA <- FALSE
NUM_THREADS <- 5L

DIR_OUT <- here("outputs", "modelo", "Modelo_1", paste0("modelo_", ANIO))
DIR_MODELO <- here(
    "data", "processed", "Modelos", "Modelo_1", paste0("modelo_", ANIO)
)
DIR_FIGURAS <- here(
    "outputs", "figures", "modelo", "Modelo_1", paste0("modelo_", ANIO)
)
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MODELO, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIGURAS, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. DATOS Y DISCRETIZACION PARA RW2
# ==============================================================================
ruta_datos <- here(
    "data", "processed", "Maestro", as.character(ANIO),
    sprintf("dataset_maestro_inla_%d_%s.rds", ANIO, ESCALA)
)
ruta_malla <- here(
    "data", "processed", "Malla", "NO2",
    sprintf("malla_spde_madrid_%s.rds", MALLA)
)
if (!file.exists(ruta_datos)) stop("No existe el maestro: ", ruta_datos)
if (!file.exists(ruta_malla)) stop("No existe la malla: ", ruta_malla)

df <- as.data.table(readRDS(ruta_datos))
columnas_necesarias <- unique(c(
    "ESTACION", "FECHA", "X_km", "Y_km", RESPUESTA,
    COVARIABLES, unname(TRANSFORMACIONES)
))
columnas_ausentes <- setdiff(columnas_necesarias, names(df))
if (length(columnas_ausentes) > 0L) {
    stop("Faltan columnas: ", paste(columnas_ausentes, collapse = ", "))
}

df[, FECHA := as.Date(FECHA)]
if (!is.null(FECHA_INICIO)) df <- df[FECHA >= FECHA_INICIO]
if (!is.null(FECHA_FIN)) df <- df[FECHA <= FECHA_FIN]
df <- df[complete.cases(df[, ..columnas_necesarias])]
if (nrow(df) == 0L) stop("No quedan observaciones completas.")

setorder(df, FECHA, ESTACION)
df[, y := get(RESPUESTA)]
variables_rw2 <- names(N_GRUPOS_RW2)
if (length(variables_rw2) == 0L || anyDuplicated(variables_rw2)) {
    stop("N_GRUPOS_RW2 debe tener nombres de variables unicos.")
}
if (any(!is.finite(N_GRUPOS_RW2) | N_GRUPOS_RW2 < 3L)) {
    stop("Cada efecto RW2 necesita al menos tres grupos.")
}

etiquetas_grafico_rw2 <- setNames(gsub("_", " ", variables_rw2), variables_rw2)
etiquetas_definidas <- intersect(variables_rw2, names(ETIQUETAS_RW2))
etiquetas_grafico_rw2[etiquetas_definidas] <- ETIQUETAS_RW2[etiquetas_definidas]

for (variable in variables_rw2) {
    nombre_agrupado <- paste0(variable, "_rw2")
    set(
        df,
        j = nombre_agrupado,
        value = inla.group(
            df[[variable]],
            n = N_GRUPOS_RW2[[variable]],
            method = METODO_GRUPOS_RW2
        )
    )
}

soporte_grupos_rw2 <- rbindlist(lapply(variables_rw2, function(variable) {
    nombre_agrupado <- paste0(variable, "_rw2")
    tabla <- data.table(nivel = df[[nombre_agrupado]])[, .N, by = nivel]
    tabla[, covariable := variable]
    setcolorder(tabla, c("covariable", "nivel", "N"))
    tabla
}))
setorder(soporte_grupos_rw2, covariable, nivel)

cat("\n--- Datos comunes a todos los modelos ---\n")
cat("Observaciones:", nrow(df), "| estaciones:", uniqueN(df$ESTACION), "\n")
grupos_obtenidos <- vapply(
    paste0(variables_rw2, "_rw2"),
    function(variable) uniqueN(df[[variable]]),
    integer(1)
)
names(grupos_obtenidos) <- variables_rw2
grupos_solicitados <- as.integer(N_GRUPOS_RW2[variables_rw2])
if (any(grupos_obtenidos < 3L)) {
    variables_sin_soporte <- variables_rw2[grupos_obtenidos < 3L]
    stop(
        "No hay al menos tres valores agrupados para: ",
        paste(variables_sin_soporte, collapse = ", "),
        ". No es posible estimar un RW2 fiable."
    )
}
if (any(grupos_obtenidos < grupos_solicitados)) {
    variables_con_empates <- variables_rw2[
        grupos_obtenidos < grupos_solicitados
    ]
    warning(
        "Algunas covariables tienen menos grupos distintos de los solicitados ",
        "por sus valores repetidos: ",
        paste(
            sprintf(
                "%s (%d de %d)",
                variables_con_empates,
                grupos_obtenidos[variables_con_empates],
                N_GRUPOS_RW2[variables_con_empates]
            ),
            collapse = ", "
        ),
        ". Revise especialmente esas curvas antes de mantener el RW2."
    )
}
cat(
    "Grupos RW2:",
    paste(sprintf("%s = %d", variables_rw2, grupos_obtenidos), collapse = " | "),
    "\n"
)


# ==============================================================================
# 2. MALLA, CAMPO ESPACIAL Y STACK COMUN
# ==============================================================================
# Todos los modelos comparten exactamente la misma malla, priors y matriz de
# proyeccion. De este modo, solo cambia la forma lineal/RW2 de una covariable.
mesh <- readRDS(ruta_malla)
spde <- crear_spde(mesh, escala = ESCALA)
indice_campo <- inla.spde.make.index(
    name = "campo_espacial",
    n.spde = spde$n.spde
)
coords <- as.matrix(df[, .(X_km, Y_km)])
A_campo <- inla.spde.make.A(mesh = mesh, loc = coords)

stopifnot(
    "A_campo y df tienen distinto numero de filas" =
        nrow(A_campo) == nrow(df),
    "A_campo y la malla tienen distinta dimension" =
        ncol(A_campo) == mesh$n
)

# Se construye un data.frame con los efectos fijos y los efectos aleatorios RW2 de las covariables.
# Este dataframe sirve para el stack de INLA y para extraer los efectos RW2 posteriores.
efectos_stack <- data.frame(Intercept = rep(1, nrow(df)))
for (variable in COVARIABLES) {
    efectos_stack[[variable]] <- df[[variable]]
}
for (variable_transformada in unname(TRANSFORMACIONES)) {
    efectos_stack[[variable_transformada]] <- df[[variable_transformada]]
}
for (variable in variables_rw2) {
    nombre_agrupado <- paste0(variable, "_rw2")
    efectos_stack[[nombre_agrupado]] <- df[[nombre_agrupado]]
}

# Se construye un stack espacial con los datos y las covariables.
inla_stack <- inla.stack(
    data = list(y_response = df$y),
    A = list(A_campo, 1),
    effects = list(indice_campo, efectos_stack),
    tag = "estimacion"
)

stopifnot(
    "Stack y datos tienen distinto numero de filas" =
        nrow(inla.stack.data(inla_stack)) == nrow(df)
)


# ==============================================================================
# 3. FORMULAS A COMPARAR
# ==============================================================================
construir_formula <- function(
  variable = NULL,
  tipo = "Lineal",
  variable_transformada = NULL
) {
    terminos_lineales <- if (is.null(variable)) {
        COVARIABLES
    } else {
        setdiff(COVARIABLES, variable)
    }
    termino_alternativo <- switch(tipo,
        Lineal = character(),
        Transformacion = variable_transformada,
        RW2 = sprintf(
            "f(%s_rw2, model = 'rw2', scale.model = TRUE)",
            variable
        ),
        stop("Tipo de efecto no reconocido: ", tipo)
    )
    as.formula(paste(
        "y_response ~ -1 + Intercept +",
        paste(
            c(
                terminos_lineales,
                termino_alternativo,
                "f(campo_espacial, model = spde)"
            ),
            collapse = " + "
        )
    ))
}

especificaciones_modelos <- rbindlist(list(
    data.table(
        Covariable = "Todas",
        Tipo = "Lineal",
        Forma = "Lineal original",
        Variable_modelo = NA_character_
    ),
    rbindlist(lapply(names(TRANSFORMACIONES), function(variable) {
        data.table(
            Covariable = variable,
            Tipo = c("Transformacion", "RW2"),
            Forma = c(ETIQUETAS_TRANSFORMACIONES[[variable]], "RW2"),
            Variable_modelo = c(
                TRANSFORMACIONES[[variable]],
                paste0(variable, "_rw2")
            )
        )
    }))
))
especificaciones_modelos[, Id := paste0("M", .I - 1L)]
especificaciones_modelos[, Modelo := fifelse(
    Tipo == "Lineal",
    paste0(Id, ": todas lineales originales"),
    sprintf("%s: %s - %s", Id, Covariable, Forma)
)]

formulas <- setNames(lapply(seq_len(nrow(especificaciones_modelos)), function(i) {
    if (especificaciones_modelos$Tipo[[i]] == "Lineal") {
        construir_formula()
    } else {
        construir_formula(
            variable = especificaciones_modelos$Covariable[[i]],
            tipo = especificaciones_modelos$Tipo[[i]],
            variable_transformada = if (
                especificaciones_modelos$Tipo[[i]] == "Transformacion"
            ) {
                especificaciones_modelos$Variable_modelo[[i]]
            } else {
                NULL
            }
        )
    }
}), especificaciones_modelos$Id)
etiquetas_modelos <- setNames(
    especificaciones_modelos$Modelo,
    especificaciones_modelos$Id
)


# ==============================================================================
# 4. AJUSTE DE LOS MODELOS
# ==============================================================================
ajustar_modelo <- function(
  formula,
  etiqueta,
  stack_modelo = inla_stack,
  calcular_cpo = CALCULAR_CPO,
  verbose = VERBOSE_INLA
) {
    cat("\nAjustando ", etiqueta, "...\n", sep = "")
    flush.console()
    inicio <- Sys.time()

    modelo <- inla(
        formula = formula,
        data = inla.stack.data(stack_modelo, spde = spde),
        family = FAMILIA,
        verbose = verbose,
        num.threads = NUM_THREADS,
        control.predictor = list(
            A = inla.stack.A(stack_modelo),
            compute = TRUE
        ),
        control.inla = list(
            strategy = "gaussian",
            int.strategy = "eb"
        ),
        inla.mode = "experimental",
        control.compute = list(
            cpo = calcular_cpo,
            dic = TRUE,
            waic = TRUE,
            openmp.strategy = "huge"
        )
    )

    minutos <- as.numeric(difftime(Sys.time(), inicio, units = "mins"))
    cat(sprintf("%s completado en %.2f minutos.\n", etiqueta, minutos))
    list(model = modelo, minutes = minutos)
}

ajustes <- lapply(names(formulas), function(id_modelo) {
    ajustar_modelo(formulas[[id_modelo]], etiquetas_modelos[[id_modelo]])
})
names(ajustes) <- names(formulas)
modelos <- lapply(ajustes, `[[`, "model")
tiempos <- vapply(ajustes, `[[`, numeric(1), "minutes")


# ==============================================================================
# 5. METRICAS COMPARABLES
# ==============================================================================
media_finita <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
}

varianza_residual_gaussiana <- function(modelo) {
    nombre <- grep(
        "^Precision for the Gaussian observations",
        names(modelo$marginals.hyperpar),
        value = TRUE
    )
    if (length(nombre) != 1L) {
        stop("No se ha podido identificar la precision residual gaussiana.")
    }
    inla.emarginal(
        function(precision) 1 / precision,
        modelo$marginals.hyperpar[[nombre]]
    )
}

resumir_spde <- function(modelo) {
    resultado <- inla.spde2.result(
        inla = modelo,
        name = "campo_espacial",
        spde = spde,
        do.transf = TRUE
    )
    rango <- inla.zmarginal(
        resultado$marginals.range.nominal[[1]],
        silent = TRUE
    )
    marginal_sd <- inla.tmarginal(
        sqrt,
        resultado$marginals.variance.nominal[[1]]
    )
    sd_espacial <- inla.zmarginal(marginal_sd, silent = TRUE)

    list(
        range_mean_km = rango$mean,
        range_q025_km = rango$quant0.025,
        range_q975_km = rango$quant0.975,
        spatial_sd_mean = sd_espacial$mean,
        spatial_sd_q025 = sd_espacial$quant0.025,
        spatial_sd_q975 = sd_espacial$quant0.975
    )
}

resumir_modelo <- function(modelo, id_modelo, minutos) {
    indices <- inla.stack.index(inla_stack, tag = "estimacion")$data
    fitted <- modelo$summary.fitted.values[indices, , drop = FALSE]
    media_predicha <- fitted[, "mean"]
    sd_predictiva <- sqrt(
        fitted[, "sd"]^2 + varianza_residual_gaussiana(modelo)
    )
    error <- media_predicha - df$y
    cov95 <- 100 * mean(
        df$y >= media_predicha - 1.96 * sd_predictiva &
            df$y <= media_predicha + 1.96 * sd_predictiva,
        na.rm = TRUE
    )

    cpo_medio <- NA_real_
    fallos_cpo <- NA_integer_
    if (CALCULAR_CPO) {
        cpo <- modelo$cpo$cpo[indices]
        cpo_valido <- is.finite(cpo) & cpo > 0
        cpo_medio <- media_finita(-log(cpo[cpo_valido]))
        fallos_cpo <- sum(!cpo_valido) +
            sum(modelo$cpo$failure[indices] != 0, na.rm = TRUE)
    }

    resumen_espacial <- resumir_spde(modelo)

    data.table(
        Id = id_modelo,
        WAIC = modelo$waic$waic,
        DIC = modelo$dic$dic,
        COV95 = cov95,
        CPO = cpo_medio,
        RMSE = sqrt(mean(error^2, na.rm = TRUE)),
        CPO_failures = fallos_cpo,
        p_eff_WAIC = modelo$waic$p.eff,
        p_eff_DIC = modelo$dic$p.eff,
        range_mean_km = resumen_espacial$range_mean_km,
        range_q025_km = resumen_espacial$range_q025_km,
        range_q975_km = resumen_espacial$range_q975_km,
        spatial_sd_mean = resumen_espacial$spatial_sd_mean,
        spatial_sd_q025 = resumen_espacial$spatial_sd_q025,
        spatial_sd_q975 = resumen_espacial$spatial_sd_q975,
        minutes = minutos
    )
}

metricas_completas <- rbindlist(lapply(names(modelos), function(id_modelo) {
    resumir_modelo(
        modelos[[id_modelo]],
        id_modelo,
        tiempos[[id_modelo]]
    )
}))

metricas_completas <- merge(
    especificaciones_modelos,
    metricas_completas,
    by = "Id",
    all.x = TRUE,
    sort = FALSE
)
metricas_completas[, orden_modelo := as.integer(sub("M", "", Id))]
setorder(metricas_completas, orden_modelo)
metricas_completas[, orden_modelo := NULL]

# RMSE posterior de ajuste: no necesita volver a ajustar el modelo y, por
# tanto, permite cribar rapidamente todos los candidatos. No es RMSE LOSO.
tabla_comparacion <- rbindlist(lapply(names(TRANSFORMACIONES), function(variable) {
    ids_variable <- c(
        "M0",
        especificaciones_modelos[Covariable == variable, Id]
    )
    resultado <- metricas_completas[match(ids_variable, Id)]
    resultado[, Covariable_comparada := variable]
    resultado[Id == "M0", Forma := "Lineal original"]
    resultado[, Delta_WAIC := WAIC - min(WAIC, na.rm = TRUE)]
    resultado[, .(
        Covariable = Covariable_comparada,
        Id,
        Forma,
        Modelo,
        WAIC,
        Delta_WAIC,
        DIC,
        RMSE,
        p_eff_WAIC,
        minutes
    )]
}))


# ==============================================================================
# 6. VALIDACION LOSO OPCIONAL PARA LOS CANDIDATOS FINALES
# ==============================================================================
tabla_comparacion_loso <- NULL
predicciones_loso <- NULL
metricas_loso_estacion <- NULL

if (CALCULAR_LOSO) {
    if (tolower(FAMILIA) != "gaussian") {
        stop("LOSO requiere FAMILIA = 'gaussian' en esta implementacion.")
    }
    modelos_loso_ausentes <- setdiff(MODELOS_LOSO, names(modelos))
    if (length(MODELOS_LOSO) == 0L || length(modelos_loso_ausentes) > 0L) {
        stop(
            "Revise MODELOS_LOSO. Identificadores disponibles: ",
            paste(names(modelos), collapse = ", ")
        )
    }

    indices_por_estacion <- split(
        seq_len(nrow(df)),
        factor(df$ESTACION, levels = unique(df$ESTACION))
    )
    grupos_loso <- unname(indices_por_estacion[as.character(df$ESTACION)])
    max_tamano_estacion <- max(lengths(indices_por_estacion))

    calcular_loso <- function(modelo, id_modelo) {
        cat("\nCalculando LOSO de ", id_modelo, "...\n", sep = "")
        resultado <- inla.group.cv(
            result = modelo,
            groups = grupos_loso,
            size.max = max_tamano_estacion,
            verbose = FALSE
        )
        if (length(resultado$mean) != nrow(df)) {
            stop("Numero inesperado de predicciones LOSO en ", id_modelo)
        }

        sd_predictiva <- sqrt(
            resultado$sd^2 + varianza_residual_gaussiana(modelo)
        )
        error <- resultado$mean - df$y
        cubierto <- df$y >= resultado$mean - 1.96 * sd_predictiva &
            df$y <= resultado$mean + 1.96 * sd_predictiva
        cpo_valido <- is.finite(resultado$cv) & resultado$cv > 0
        cpo_loso <- media_finita(-log(resultado$cv[cpo_valido]))

        predicciones <- data.table(
            Modelo = id_modelo,
            ESTACION = df$ESTACION,
            FECHA = df$FECHA,
            observed = df$y,
            predicted_mean = resultado$mean,
            predictive_sd = sd_predictiva,
            lower_95 = resultado$mean - 1.96 * sd_predictiva,
            upper_95 = resultado$mean + 1.96 * sd_predictiva,
            error = error,
            covered_95 = cubierto,
            predictive_density = resultado$cv
        )

        por_estacion <- predicciones[, .(
            observations = .N,
            RMSE = sqrt(mean(error^2, na.rm = TRUE)),
            COV95 = 100 * mean(covered_95, na.rm = TRUE),
            CPO = media_finita(-log(
                predictive_density[
                    is.finite(predictive_density) & predictive_density > 0
                ]
            ))
        ), by = .(Modelo, ESTACION)]

        global <- data.table(
            Modelo = id_modelo,
            WAIC = modelo$waic$waic,
            DIC = modelo$dic$dic,
            COV95 = 100 * mean(cubierto, na.rm = TRUE),
            CPO = cpo_loso,
            RMSE = sqrt(mean(error^2, na.rm = TRUE))
        )
        list(global = global, predictions = predicciones, station = por_estacion)
    }

    loso <- lapply(MODELOS_LOSO, function(id_modelo) {
        calcular_loso(modelos[[id_modelo]], id_modelo)
    })
    tabla_comparacion_loso <- rbindlist(lapply(loso, `[[`, "global"))
    tabla_comparacion_loso <- merge(
        especificaciones_modelos[, .(
            Modelo = Id,
            Descripcion = Modelo,
            Covariable,
            Forma
        )],
        tabla_comparacion_loso,
        by = "Modelo",
        all.y = TRUE,
        sort = FALSE
    )
    predicciones_loso <- rbindlist(lapply(loso, `[[`, "predictions"))
    metricas_loso_estacion <- rbindlist(lapply(loso, `[[`, "station"))
}


# ==============================================================================
# 7. EFECTOS NO LINEALES: DATOS Y GRAFICO
# ==============================================================================
extraer_coeficiente <- function(modelo, termino) {
    if (!termino %in% rownames(modelo$summary.fixed)) {
        stop("No se encontro el coeficiente fijo de ", termino, ".")
    }
    resumen <- modelo$summary.fixed[termino, , drop = FALSE]
    if (nrow(resumen) != 1L) {
        stop("No se encontro el coeficiente fijo de ", termino, ".")
    }
    c(
        mean = unname(resumen[1L, "mean"]),
        lower = unname(resumen[1L, "0.025quant"]),
        upper = unname(resumen[1L, "0.975quant"])
    )
}

curva_coeficiente <- function(
  soporte,
  z,
  coeficiente,
  variable,
  forma,
  modelo
) {
    extremo_1 <- coeficiente[["lower"]] * z
    extremo_2 <- coeficiente[["upper"]] * z
    data.table(
        Covariable = etiquetas_grafico_rw2[[variable]],
        variable = variable,
        Forma = forma,
        Modelo = modelo,
        x = soporte,
        effect_mean = coeficiente[["mean"]] * z,
        effect_lower = pmin(extremo_1, extremo_2),
        effect_upper = pmax(extremo_1, extremo_2)
    )
}

extraer_comparacion_efectos <- function(variable) {
    variable_transformada <- TRANSFORMACIONES[[variable]]
    variable_rw2 <- paste0(variable, "_rw2")
    variable_raw <- paste0(variable, "_raw")
    variable_eje <- if (variable_raw %in% names(df)) variable_raw else variable

    soporte <- df[, .(
        x = mean(get(variable_eje)),
        z_lineal = mean(get(variable)),
        z_transformada = mean(get(variable_transformada))
    ), by = .(nivel = get(variable_rw2))]
    soporte[, nivel := as.numeric(nivel)]
    setorder(soporte, nivel)

    # Se centran las contribuciones parametricas para compararlas con el RW2,
    # que INLA estima bajo una restriccion de suma cero.
    soporte[, z_lineal := z_lineal - mean(df[[variable]])]
    soporte[, z_transformada :=
        z_transformada - mean(df[[variable_transformada]])]

    id_transformacion <- especificaciones_modelos[
        Covariable == variable & Tipo == "Transformacion",
        Id
    ]
    id_rw2 <- especificaciones_modelos[
        Covariable == variable & Tipo == "RW2",
        Id
    ]

    curva_lineal <- curva_coeficiente(
        soporte = soporte$x,
        z = soporte$z_lineal,
        coeficiente = extraer_coeficiente(modelos$M0, variable),
        variable = variable,
        forma = "Lineal original",
        modelo = etiquetas_modelos[["M0"]]
    )
    curva_transformada <- curva_coeficiente(
        soporte = soporte$x,
        z = soporte$z_transformada,
        coeficiente = extraer_coeficiente(
            modelos[[id_transformacion]],
            variable_transformada
        ),
        variable = variable,
        forma = "Transformacion parametrica",
        modelo = etiquetas_modelos[[id_transformacion]]
    )

    resumen_rw2 <- as.data.table(
        modelos[[id_rw2]]$summary.random[[variable_rw2]]
    )
    if (!nrow(resumen_rw2)) {
        stop("No se encontro el efecto RW2 de ", variable, ".")
    }
    resumen_rw2[, nivel := as.numeric(ID)]
    curva_rw2 <- merge(
        soporte[, .(nivel, x)],
        resumen_rw2[, .(
            nivel,
            effect_mean = mean,
            effect_lower = `0.025quant`,
            effect_upper = `0.975quant`
        )],
        by = "nivel",
        all.x = TRUE,
        sort = FALSE
    )
    if (anyNA(curva_rw2[, .(effect_mean, effect_lower, effect_upper)])) {
        stop("No se pudo alinear el soporte agrupado y el resumen RW2 de ", variable, ".")
    }
    curva_rw2[, `:=`(
        Covariable = etiquetas_grafico_rw2[[variable]],
        variable = variable,
        Forma = "RW2",
        Modelo = etiquetas_modelos[[id_rw2]]
    )]
    curva_rw2 <- curva_rw2[, .(
        Covariable,
        variable,
        Forma,
        Modelo,
        x,
        effect_mean,
        effect_lower,
        effect_upper
    )]

    rbindlist(list(curva_lineal, curva_transformada, curva_rw2))
}

efectos_comparacion <- rbindlist(lapply(
    names(TRANSFORMACIONES),
    extraer_comparacion_efectos
))
efectos_comparacion[, Forma := factor(
    Forma,
    levels = c("Lineal original", "Transformacion parametrica", "RW2")
)]

capas_efectos <- function(datos) {
    ggplot(
        datos,
        aes(
            x = x,
            y = effect_mean,
            color = Forma,
            fill = Forma,
            linetype = Forma,
            group = Forma
        )
    ) +
        geom_ribbon(
            aes(ymin = effect_lower, ymax = effect_upper),
            color = NA,
            alpha = 0.10
        ) +
        geom_line(linewidth = 0.95) +
        geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
        scale_color_manual(values = c(
            "Lineal original" = "#D55E00",
            "Transformacion parametrica" = "#009E73",
            "RW2" = "#0072B2"
        )) +
        scale_fill_manual(values = c(
            "Lineal original" = "#D55E00",
            "Transformacion parametrica" = "#009E73",
            "RW2" = "#56B4E9"
        )) +
        scale_linetype_manual(values = c(
            "Lineal original" = "dashed",
            "Transformacion parametrica" = "dotdash",
            "RW2" = "solid"
        )) +
        guides(fill = "none", linetype = "none")
}

n_columnas_grafico <- 1L
n_filas_grafico <- length(variables_rw2)

grafico_efectos <- capas_efectos(efectos_comparacion) +
    facet_wrap(
        ~Covariable,
        scales = "free",
        ncol = n_columnas_grafico
    ) +
    labs(
        title = "Forma lineal, transformacion parametrica y RW2",
        subtitle = paste(
            "Log para temperatura y radiacion; raiz cuadrada para viento.",
            "Las bandas representan los",
            "intervalos creibles posteriores del 95%"
        ),
        x = "Valor de la covariable (escala original si esta disponible)",
        y = "Contribucion al predictor lineal",
        color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold")
    )

# Renderiza primero en un archivo temporal. Asi ragg no intenta escribir sobre
# un PNG abierto en el visor de RStudio o en otra aplicacion de Windows.
guardar_png_seguro <- function(grafico, ruta, width, height, dpi = 300) {
    dir.create(dirname(ruta), recursive = TRUE, showWarnings = FALSE)
    temporal <- tempfile(
        pattern = "grafico_",
        tmpdir = dirname(ruta),
        fileext = ".png"
    )
    on.exit(if (file.exists(temporal)) unlink(temporal), add = TRUE)

    ggsave(
        filename = temporal,
        plot = grafico,
        device = ragg::agg_png,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white"
    )

    copiado <- suppressWarnings(file.copy(temporal, ruta, overwrite = TRUE))
    if (!copiado) {
        ruta_alternativa <- paste0(
            tools::file_path_sans_ext(ruta),
            "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"
        )
        if (!file.copy(temporal, ruta_alternativa, overwrite = FALSE)) {
            stop("No se ha podido guardar el grafico en: ", dirname(ruta))
        }
        message(
            "El PNG anterior estaba abierto. El nuevo grafico se guardo en: ",
            ruta_alternativa
        )
        return(ruta_alternativa)
    }
    ruta
}

ruta_grafico_efectos <- guardar_png_seguro(
    grafico = grafico_efectos,
    ruta = file.path(
        DIR_FIGURAS,
        sprintf("efectos_transformacion_vs_rw2_%d.png", ANIO)
    ),
    width = if (n_columnas_grafico == 1L) 8 else 12,
    height = max(6, 3.2 * n_filas_grafico)
)

# Ademas del panel conjunto, se genera un PNG independiente por covariable.
rutas_graficos_individuales <- setNames(
    vapply(variables_rw2, function(variable) {
        variable_actual <- variable
        datos_variable <- efectos_comparacion[variable == variable_actual]
        grafico_variable <- capas_efectos(datos_variable) +
            labs(
                title = paste("Efecto de", etiquetas_grafico_rw2[[variable]]),
                subtitle = sprintf(
                    "Lineal original, %s y RW2; intervalos creibles del 95%%",
                    tolower(ETIQUETAS_TRANSFORMACIONES[[variable]])
                ),
                x = paste0(
                    etiquetas_grafico_rw2[[variable]],
                    " (escala original si esta disponible)"
                ),
                y = "Contribucion al predictor lineal",
                color = NULL
            ) +
            theme_minimal(base_size = 11) +
            theme(
                legend.position = "bottom",
                panel.grid.minor = element_blank()
            )

        guardar_png_seguro(
            grafico = grafico_variable,
            ruta = file.path(
                DIR_FIGURAS,
                sprintf("comparacion_formas_%s_%d.png", variable, ANIO)
            ),
            width = 8,
            height = 5.5
        )
    }, character(1)),
    variables_rw2
)


# ==============================================================================
# 8. SALIDAS
# ==============================================================================
fwrite(tabla_comparacion, file.path(DIR_OUT, "comparacion_modelos_2025.csv"))
fwrite(metricas_completas, file.path(DIR_OUT, "metricas_completas_2025.csv"))
fwrite(soporte_grupos_rw2, file.path(DIR_OUT, "soporte_grupos_rw_2025.csv"))
fwrite(
    efectos_comparacion,
    file.path(DIR_OUT, "efectos_lineal_transformacion_rw2_2025.csv")
)

extraer_resumen <- function(modelo, id_modelo, componente) {
    tabla <- as.data.table(modelo[[componente]], keep.rownames = "term")
    tabla[, Modelo := id_modelo]
    setcolorder(tabla, c("Modelo", "term"))
    tabla
}

fwrite(
    rbindlist(lapply(names(modelos), function(id_modelo) {
        extraer_resumen(
            modelos[[id_modelo]],
            etiquetas_modelos[[id_modelo]],
            "summary.fixed"
        )
    }), fill = TRUE),
    file.path(DIR_OUT, "efectos_fijos_modelos_2025.csv")
)
fwrite(
    rbindlist(lapply(names(modelos), function(id_modelo) {
        extraer_resumen(
            modelos[[id_modelo]],
            etiquetas_modelos[[id_modelo]],
            "summary.hyperpar"
        )
    }), fill = TRUE),
    file.path(DIR_OUT, "hiperparametros_modelos_2025.csv")
)

if (CALCULAR_LOSO) {
    fwrite(
        tabla_comparacion_loso,
        file.path(DIR_OUT, "comparacion_modelos_loso_2025.csv")
    )
    fwrite(
        predicciones_loso,
        file.path(DIR_OUT, "predicciones_loso_modelos_2025.csv")
    )
    fwrite(
        metricas_loso_estacion,
        file.path(DIR_OUT, "metricas_loso_por_estacion_modelos_2025.csv")
    )
}

configuracion <- list(
    year = ANIO,
    scale = ESCALA,
    start_date = FECHA_INICIO,
    end_date = FECHA_FIN,
    response = RESPUESTA,
    covariates = COVARIABLES,
    transformations = TRANSFORMACIONES,
    mesh = MALLA,
    mesh_vertices = mesh$n,
    prior_range = PRIORS_SPDE[[ESCALA]]$prior.range,
    prior_sigma = PRIORS_SPDE[[ESCALA]]$prior.sigma,
    rw_groups = N_GRUPOS_RW2,
    rw_grouping_method = METODO_GRUPOS_RW2,
    family = FAMILIA,
    calculate_cpo = CALCULAR_CPO,
    rmse_comparison_type = "posterior_fitted_screening",
    calculate_loso = CALCULAR_LOSO,
    loso_models = MODELOS_LOSO
)

if (GUARDAR_MODELOS) {
    saveRDS(
        list(
            models = modelos,
            formulas = formulas,
            mesh = mesh,
            spde = spde,
            stack = inla_stack,
            data = df,
            configuration = configuracion,
            comparison = tabla_comparacion,
            comparison_loso = tabla_comparacion_loso,
            comparison_effects = efectos_comparacion,
            plot_files = list(
                combined = ruta_grafico_efectos,
                individual = rutas_graficos_individuales
            )
        ),
        file.path(
            DIR_MODELO,
            sprintf(
                "prueba_transformaciones_%d_%s_spde_malla_%s.rds",
                ANIO,
                ESCALA,
                MALLA
            )
        )
    )
}


# ==============================================================================
# 9. TABLA ACADEMICA DE COMPARACION
# ==============================================================================
source(here("R", "utilities", "academic_quality_tables.R"))

tabla_para_presentar <- copy(tabla_comparacion)
metodo_tabla <- "RMSE posterior de ajuste (cribado rapido; no es validacion cruzada)"

formatear <- function(x, digitos = 3L) {
    ifelse(
        is.finite(x),
        formatC(x, format = "f", digits = digitos),
        "--"
    )
}

tabla_academica <- tabla_para_presentar[, .(
    Variable = gsub("_", " ", Covariable),
    Model = Id,
    Form = Forma,
    WAIC = formatear(WAIC, 1L),
    `Delta WAIC` = formatear(Delta_WAIC, 1L),
    DIC = formatear(DIC, 1L),
    RMSE = formatear(RMSE, 3L),
    `p eff` = formatear(p_eff_WAIC, 1L),
    Minutes = formatear(minutes, 2L)
)]

ruta_tabla_comparacion <- file.path(
    DIR_FIGURAS,
    sprintf("table_model_comparison_%d.png", ANIO)
)
ruta_tabla_comparacion <- booktabs_png(
    tabla_academica,
    ruta_tabla_comparacion,
    title = "Comparacion INLA-SPDE de formas funcionales",
    subtitle = sprintf(
        "%s; malla %s (%d vertices); RW2 groups (%s): %s",
        metodo_tabla,
        MALLA,
        mesh$n,
        METODO_GRUPOS_RW2,
        paste(
            sprintf("%s=%d", variables_rw2, N_GRUPOS_RW2[variables_rw2]),
            collapse = ", "
        )
    ),
    note = paste0(
        "Cada bloque compara la forma lineal original, la transformacion y RW2; ",
        "Delta WAIC se calcula dentro de cada variable. Valores menores son mejores. ",
        "El RMSE de esta tabla es optimista y se ",
        "usa solo para el cribado; confirme los candidatos finales activando ",
        "LOSO. Todos los modelos usan las mismas ",
        "observaciones, covariables, verosimilitud y opciones de integracion; ",
        "solo cambia un efecto cada vez y todos incluyen el mismo campo espacial SPDE."
    ),
    widths = c(1.65, 0.55, 1.55, 0.85, 0.90, 0.85, 0.85, 0.70, 0.75),
    align = c("left", "center", "left", rep("right", 6)),
    font_size = 7.8,
    row_height = 0.30
)


# ==============================================================================
# 10. RESUMEN EN CONSOLA
# ==============================================================================
cat("\n--- Comparacion completada ---\n")
print(tabla_comparacion)
cat("RMSE mostrado: posterior de ajuste (cribado rapido; no es LOSO).\n")
if (CALCULAR_LOSO) {
    cat("\n--- Validacion LOSO de candidatos ---\n")
    print(tabla_comparacion_loso)
    cat(
        "Tabla LOSO CSV: ",
        file.path(DIR_OUT, "comparacion_modelos_loso_2025.csv"),
        "\n",
        sep = ""
    )
}
cat("\nTabla CSV: ", file.path(DIR_OUT, "comparacion_modelos_2025.csv"), "\n", sep = "")
cat("Grafico lineal/transformacion/RW2: ", ruta_grafico_efectos, "\n", sep = "")
cat("Tabla PNG: ", ruta_tabla_comparacion, "\n", sep = "")
