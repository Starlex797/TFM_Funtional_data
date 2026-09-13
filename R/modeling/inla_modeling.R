# ==============================================================================
# FUNCIONES COMUNES PARA AJUSTAR MODELOS INLA
# ==============================================================================
# Este archivo concentra la parte operativa de los modelos INLA-SPDE:
#   - construccion del stack;
#   - llamada comun a inla();
#   - extraccion de coeficientes fijos y curvas RW2;
#   - ajuste y resumen de los modelos de seleccion;
#   - ajuste y evaluacion del hold-out espacial.
#
# Las decisiones propias de cada analisis (datos, formulas, malla, estaciones
# reservadas y nombres de modelos) permanecen en su script principal.
# ==============================================================================


#' Construye un stack INLA con uno o varios bloques de efectos.
#'
#' @param respuesta vector de respuesta. Puede contener NA para prediccion.
#' @param A lista de matrices de proyeccion.
#' @param effects lista de bloques de efectos asociada a A.
#' @param tag etiqueta utilizada para recuperar las filas del stack.
#' @param nombre_respuesta nombre que utiliza la respuesta en la formula.
#' @return objeto inla.stack.
crear_stack_inla <- function(
    respuesta,
    A,
    effects,
    tag = "estimacion",
    nombre_respuesta = "y_response"
) {
    if (length(A) != length(effects)) {
        stop("A y effects deben contener el mismo numero de bloques.")
    }
    datos_stack <- setNames(list(respuesta), nombre_respuesta)
    INLA::inla.stack(
        data = datos_stack,
        A = A,
        effects = effects,
        tag = tag
    )
}


#' Construye un stack con un campo SPDE y un bloque de efectos fijos.
#'
#' @param respuesta vector de respuesta. Puede contener NA para prediccion.
#' @param A_espacial matriz de proyeccion entre observaciones y malla.
#' @param indice_espacial indice creado con inla.spde.make.index().
#' @param efectos_fijos data.frame con intercepto y covariables.
#' @param tag etiqueta utilizada para recuperar las filas del stack.
#' @param nombre_respuesta nombre que utiliza la respuesta en la formula.
#' @return objeto inla.stack.
crear_stack_inla_spde <- function(
    respuesta,
    A_espacial,
    indice_espacial,
    efectos_fijos,
    tag = "estimacion",
    nombre_respuesta = "y_response"
) {
    n <- length(respuesta)
    if (nrow(A_espacial) != n || nrow(efectos_fijos) != n) {
        stop(
            "La respuesta, A_espacial y efectos_fijos deben tener ",
            "el mismo numero de filas."
        )
    }

    crear_stack_inla(
        respuesta = respuesta,
        A = list(A_espacial, 1),
        effects = list(indice_espacial, efectos_fijos),
        tag = tag,
        nombre_respuesta = nombre_respuesta
    )
}


#' Recupera los indices de las observaciones asociados a un tag del stack.
indices_stack_inla <- function(stack, tag) {
    INLA::inla.stack.index(stack, tag = tag)$data
}


#' Ejecuta una llamada comun a inla() para un modelo con stack y SPDE.
#'
#' @param calcular_criterios si TRUE calcula DIC y WAIC.
#' @return objeto devuelto por INLA::inla().
ejecutar_inla_spde <- function(
    formula_modelo,
    stack,
    spde = NULL,
    familia = "gaussian",
    verbose = FALSE,
    num_threads = 1L,
    calcular_criterios = TRUE,
    calcular_cpo = FALSE,
    strategy = "gaussian",
    int_strategy = "eb",
    inla_mode = NULL
) {
    datos_inla <- if (is.null(spde)) {
        INLA::inla.stack.data(stack)
    } else {
        INLA::inla.stack.data(stack, spde = spde)
    }

    argumentos <- list(
        formula = formula_modelo,
        data = datos_inla,
        family = familia,
        verbose = verbose,
        num.threads = num_threads,
        control.predictor = list(
            A = INLA::inla.stack.A(stack),
            compute = TRUE
        ),
        control.inla = list(
            strategy = strategy,
            int.strategy = int_strategy
        ),
        control.compute = list(
            dic = calcular_criterios,
            waic = calcular_criterios,
            cpo = calcular_cpo,
            openmp.strategy = "huge"
        )
    )
    if (!is.null(inla_mode)) argumentos$inla.mode <- inla_mode
    if (!is.null(datos_inla$link)) argumentos$control.predictor$link <- datos_inla$link
    if (!is.null(attr(stack, "control.fixed"))) {
        argumentos$control.fixed <- attr(stack, "control.fixed")
    }
    if (!is.null(attr(stack, "control.family"))) {
        argumentos$control.family <- attr(stack, "control.family")
    }
    resultado <- do.call(INLA::inla, argumentos)
    resultado$likelihood_no2 <- attr(stack, "likelihood_no2")
    resultado
}


#' Calcula la varianza residual de un modelo INLA gaussiano.
#'
#' Usa E(1 / precision) sobre la marginal posterior completa.
varianza_residual_gaussiana <- function(modelo, likelihood = modelo$likelihood_no2) {
    nombre <- grep(
        "^Precision for the Gaussian observations",
        names(modelo$marginals.hyperpar),
        value = TRUE
    )
    if (!is.null(likelihood) && length(nombre) > 1L) {
        # INLA enumera las precisiones por verosimilitud; primera columna = NO2.
        if (length(likelihood) != 1L || !likelihood %in% seq_along(nombre)) {
            stop("Indice de verosimilitud gaussiana invalido.")
        }
        nombre <- nombre[likelihood]
    }
    if (length(nombre) != 1L) {
        stop("No se ha podido identificar la precision residual gaussiana.")
    }

    INLA::inla.emarginal(
        function(precision) 1 / precision,
        modelo$marginals.hyperpar[[nombre]]
    )
}


#' Extrae una tabla academica de los efectos fijos de un modelo INLA.
extraer_coeficientes_fijos_inla <- function(
    modelo,
    excluir = "Intercept",
    id_modelo = NULL
) {
    resumen <- modelo$summary.fixed
    if (is.null(resumen) || nrow(resumen) == 0L) {
        return(data.table::data.table(
            Variable = character(),
            Coeficiente = numeric(),
            IC95 = character(),
            Significativa_95 = character()
        ))
    }

    tabla <- data.table::data.table(
        Variable = rownames(resumen),
        Coeficiente = resumen[["mean"]],
        IC95 = sprintf(
            "[%.4f, %.4f]",
            resumen[["0.025quant"]],
            resumen[["0.975quant"]]
        ),
        Significativa_95 = data.table::fifelse(
            resumen[["0.025quant"]] > 0 |
                resumen[["0.975quant"]] < 0,
            "Si",
            "No"
        )
    )
    tabla <- tabla[!Variable %in% excluir]
    if (!is.null(id_modelo)) {
        tabla[, Modelo := id_modelo]
        data.table::setcolorder(tabla, c("Modelo", setdiff(names(tabla), "Modelo")))
    }
    tabla
}


#' Ajusta un modelo INLA y registra el tiempo de ejecucion.
ajustar_modelo_cronometrado_inla <- function(
    formula_modelo,
    stack,
    id_modelo,
    spde = NULL,
    familia = "gaussian",
    verbose = FALSE,
    num_threads = 1L,
    calcular_cpo = TRUE,
    int_strategy = "eb"
) {
    cat("\nAjustando ", id_modelo, "...\n", sep = "")
    inicio <- Sys.time()
    modelo <- ejecutar_inla_spde(
        formula_modelo = formula_modelo,
        stack = stack,
        spde = spde,
        familia = familia,
        verbose = verbose,
        num_threads = num_threads,
        calcular_criterios = TRUE,
        calcular_cpo = calcular_cpo,
        int_strategy = int_strategy
    )
    minutos <- as.numeric(difftime(Sys.time(), inicio, units = "mins"))
    cat(sprintf("%s completado en %.2f minutos.\n", id_modelo, minutos))
    list(modelo = modelo, stack = stack, minutos = minutos)
}


#' Resume ajuste, complejidad y cobertura predictiva de un modelo INLA.
resumir_ajuste_inla <- function(
    ajuste,
    id_modelo,
    procedimiento,
    y_observada,
    tag = "estimacion",
    nivel = 0.95
) {
    modelo <- ajuste$modelo
    indices <- indices_stack_inla(ajuste$stack, tag)
    resumen <- modelo$summary.fitted.values[indices, , drop = FALSE]
    if (nrow(resumen) != length(y_observada)) {
        stop("El resumen ajustado y la respuesta tienen distinta longitud.")
    }

    media <- resumen[, "mean"]
    sd_predictiva <- sqrt(
        resumen[, "sd"]^2 + varianza_residual_gaussiana(modelo)
    )
    z <- stats::qnorm(1 - (1 - nivel) / 2)
    error <- media - y_observada
    dentro <- y_observada >= media - z * sd_predictiva &
        y_observada <= media + z * sd_predictiva

    cpo_score <- NA_real_
    cpo_failures <- NA_integer_
    if (!is.null(modelo$cpo$cpo)) {
        cpo <- modelo$cpo$cpo[indices]
        failure <- modelo$cpo$failure[indices]
        evaluable <- is.finite(y_observada)
        valid_cpo <- evaluable & is.finite(cpo) & cpo > 0
        if (any(valid_cpo)) {
            cpo_score <- mean(-log(cpo[valid_cpo]))
        }
        failed_cpo <- evaluable & (
            !is.finite(cpo) | cpo <= 0 |
                (!is.na(failure) & failure != 0)
        )
        cpo_failures <- sum(failed_cpo)
    }

    # En modelos conjuntos no sumar la verosimilitud meteorologica a la de NO2.
    # Los indices del tag contienen exclusivamente la misma respuesta evaluada.
    seleccion <- indices[is.finite(y_observada)]
    criterio_local <- function(bloque, local, total) {
        if (!is.null(modelo$likelihood_no2)) {
            valores <- bloque[[local]][seleccion]
            if (length(valores) != length(seleccion) || any(!is.finite(valores))) {
                stop("No se puede extraer el criterio NO2: ", local)
            }
            sum(valores)
        } else bloque[[total]]
    }
    data.table::data.table(
        Modelo = id_modelo,
        Procedimiento = procedimiento,
        WAIC = criterio_local(modelo$waic, "local.waic", "waic"),
        WAIC_p_eff = criterio_local(modelo$waic, "local.p.eff", "p.eff"),
        DIC = criterio_local(modelo$dic, "local.dic", "dic"),
        DIC_p_eff = criterio_local(modelo$dic, "local.p.eff", "p.eff"),
        CPO_mean_neg_log = cpo_score,
        CPO_failures = cpo_failures,
        RMSE_ajuste = sqrt(mean(error^2, na.rm = TRUE)),
        COV95_predictiva_ajuste = 100 * mean(dentro, na.rm = TRUE)
    )
}


#' Extrae las curvas RW2 presentes en un modelo INLA.
extraer_curvas_rw2_inla <- function(
    modelo,
    nombres_rw2,
    id_modelo,
    descripcion
) {
    presentes <- intersect(nombres_rw2, names(modelo$summary.random))
    if (length(presentes) == 0L) {
        return(data.table::data.table())
    }

    data.table::rbindlist(lapply(presentes, function(nombre_rw2) {
        resumen <- modelo$summary.random[[nombre_rw2]]
        resultado <- data.table::data.table(
            Modelo = id_modelo,
            Descripcion = descripcion,
            Variable = sub("_rw2$", "", nombre_rw2),
            Tipo_efecto = "RW2",
            Nivel_RW2 = resumen[["ID"]],
            Coeficiente = resumen[["mean"]],
            Q025 = resumen[["0.025quant"]],
            Q975 = resumen[["0.975quant"]]
        )
        resultado[, IC95 := sprintf("[%.4f, %.4f]", Q025, Q975)]
        resultado[, Significativa_95 := data.table::fifelse(
            Q025 > 0 | Q975 < 0,
            "Si (puntual)",
            "No (puntual)"
        )]
        resultado
    }), use.names = TRUE, fill = TRUE)
}


#' Ajusta y resume un modelo candidato de seleccion de variables.
ajustar_modelo_inla <- function(
    formula_modelo,
    id_modelo,
    descripcion,
    stack,
    spde,
    indices_observaciones,
    y_observada,
    nombres_rw2 = character(),
    familia = "gaussian",
    verbose = FALSE,
    num_threads = 1L
) {
    if (!exists("calcular_rmse_cov95", mode = "function")) {
        stop(
            "Debe cargar R/utilities/metricas_predictivas.R antes de ",
            "utilizar ajustar_modelo_inla()."
        )
    }

    cat("\nAjustando ", id_modelo, ": ", descripcion, "...\n", sep = "")
    modelo <- ejecutar_inla_spde(
        formula_modelo = formula_modelo,
        stack = stack,
        spde = spde,
        familia = familia,
        verbose = verbose,
        num_threads = num_threads,
        calcular_criterios = TRUE
    )

    resumen <- modelo$summary.fitted.values[
        indices_observaciones, ,
        drop = FALSE
    ]
    if (nrow(resumen) != length(y_observada)) {
        stop("El resumen ajustado y la respuesta tienen distinta longitud.")
    }

    metricas_ajuste <- calcular_rmse_cov95(
        y_obs = y_observada,
        media_pred = resumen[, "mean"],
        sd_pred = resumen[, "sd"],
        varianza_residual = varianza_residual_gaussiana(modelo)
    )

    list(
        modelo = modelo,
        metricas = data.table::data.table(
            Modelo = id_modelo,
            Descripcion = descripcion,
            WAIC = modelo$waic$waic,
            DIC = modelo$dic$dic,
            RMSE_INLA_orientativo = metricas_ajuste$RMSE,
            COV95_INLA_orientativo = metricas_ajuste$COV95
        ),
        coeficientes = extraer_coeficientes_fijos_inla(modelo),
        curvas_rw2 = extraer_curvas_rw2_inla(
            modelo = modelo,
            nombres_rw2 = nombres_rw2,
            id_modelo = id_modelo,
            descripcion = descripcion
        )
    )
}


#' Ajusta y evalua un modelo mediante un hold-out de estaciones completas.
ajustar_holdout_espacial_inla <- function(
    formula_modelo,
    id_modelo,
    descripcion,
    stack,
    indices_observaciones,
    datos,
    es_holdout,
    spde = NULL,
    familia = "gaussian",
    verbose = FALSE,
    num_threads = 1L,
    nivel = 0.95,
    inla_mode = NULL,
    idioma = c("es", "en"),
    int_strategy = "eb"
) {
    idioma <- match.arg(idioma)
    datos <- data.table::as.data.table(datos)
    columnas <- c("ESTACION", "NOM_TIPO", "FECHA", "y")
    faltan <- setdiff(columnas, names(datos))
    if (length(faltan) > 0L) {
        stop("Faltan columnas para el hold-out: ", paste(faltan, collapse = ", "))
    }
    if (length(es_holdout) != nrow(datos)) {
        stop("es_holdout debe tener una posicion por fila de datos.")
    }

    cat("\nAjustando hold-out de ", id_modelo, "...\n", sep = "")
    modelo <- ejecutar_inla_spde(
        formula_modelo = formula_modelo,
        stack = stack,
        spde = spde,
        familia = familia,
        verbose = verbose,
        num_threads = num_threads,
        calcular_criterios = FALSE,
        int_strategy = int_strategy,
        inla_mode = inla_mode
    )

    resumen <- modelo$summary.fitted.values[
        indices_observaciones, ,
        drop = FALSE
    ]
    if (nrow(resumen) != nrow(datos)) {
        stop("El resumen predictivo y los datos tienen distinta longitud.")
    }

    media <- resumen[, "mean"]
    sd_media <- resumen[, "sd"]
    var_residual <- varianza_residual_gaussiana(modelo)
    sd_predictiva <- sqrt(sd_media^2 + var_residual)
    z <- stats::qnorm(1 - (1 - nivel) / 2)

    predicciones <- datos[es_holdout, .(
        Modelo = id_modelo,
        ESTACION,
        NOM_TIPO = as.character(NOM_TIPO),
        FECHA,
        Observado = y
    )]
    predicciones[, `:=`(
        Predicho = media[es_holdout],
        SD_media = sd_media[es_holdout],
        SD_predictiva = sd_predictiva[es_holdout]
    )]
    predicciones[, `:=`(
        Limite_inferior_95 = Predicho - z * SD_predictiva,
        Limite_superior_95 = Predicho + z * SD_predictiva,
        Error = Predicho - Observado
    )]
    predicciones[, `:=`(
        Dentro_IC95 = Observado >= Limite_inferior_95 &
            Observado <= Limite_superior_95,
        Anchura_IC95 = Limite_superior_95 - Limite_inferior_95
    )]

    validas <- predicciones[
        is.finite(Observado) & is.finite(Predicho) &
            is.finite(SD_predictiva)
    ]
    if (nrow(validas) == 0L) {
        stop("No hay predicciones hold-out validas para ", id_modelo, ".")
    }

    global <- validas[, .(
        Modelo = id_modelo,
        Descripcion = descripcion,
        RMSE_HOLDOUT = sqrt(mean(Error^2)),
        MAE_HOLDOUT = mean(abs(Error)),
        Sesgo_HOLDOUT = mean(Error),
        COV95_HOLDOUT = 100 * mean(Dentro_IC95),
        Anchura_media_IC95_HOLDOUT = mean(Anchura_IC95),
        Observaciones_HOLDOUT = .N
    )]

    por_estacion <- validas[, .(
        Observaciones = .N,
        RMSE = sqrt(mean(Error^2)),
        MAE = mean(abs(Error)),
        Sesgo = mean(Error),
        COV95 = 100 * mean(Dentro_IC95),
        Anchura_media_IC95 = mean(Anchura_IC95)
    ), by = .(Modelo, ESTACION, NOM_TIPO)]

    if (idioma == "en") {
        data.table::setnames(
            global,
            c(
                "Modelo", "Descripcion", "Sesgo_HOLDOUT",
                "Anchura_media_IC95_HOLDOUT", "Observaciones_HOLDOUT"
            ),
            c(
                "Model", "Description", "Bias_HOLDOUT",
                "Mean_width_95_HOLDOUT", "Holdout_observations"
            )
        )
        data.table::setnames(
            por_estacion,
            c(
                "Modelo", "ESTACION", "NOM_TIPO", "Observaciones",
                "Sesgo", "Anchura_media_IC95"
            ),
            c(
                "Model", "Station", "Station_type", "Observations",
                "Bias", "Mean_width_95"
            )
        )
        data.table::setnames(
            predicciones,
            c(
                "Modelo", "ESTACION", "NOM_TIPO", "FECHA", "Observado",
                "Predicho", "SD_media", "SD_predictiva",
                "Limite_inferior_95", "Limite_superior_95", "Dentro_IC95",
                "Anchura_IC95"
            ),
            c(
                "Model", "Station", "Station_type", "Date", "Observed",
                "Predicted", "Mean_SD", "Predictive_SD", "Lower_95",
                "Upper_95", "Covered_95", "Width_95"
            )
        )
    }

    list(
        modelo = modelo,
        global = global,
        por_estacion = por_estacion,
        predicciones = predicciones
    )
}
