# ==============================================================================
# MODELO ESPACIAL INLA-SPDE - 2025
# ==============================================================================


library(INLA)
library(data.table)
library(here)


# ==============================================================================
# CONFIGURACION EDITABLE
# ==============================================================================
ANIO <- c(2019, 2020, 2021)
ESCALA <- "DIARIO"

# Periodo usado. Escribir NULL en ambos para utilizar todo el fichero.
FECHA_INICIO <- as.Date("2019-01-01")
FECHA_FIN <- as.Date("2021-12-31")

# Respuesta y covariables. Basta con editar este vector: la formula y el stack
# se actualizan automaticamente.
RESPUESTA <- "LOG_NO2_DIARIO"
COVARIABLES <- c(
    "intensidad",
    "Temperatura",
    "Velocidad_Viento",
)

# Malla elegida en el analisis de sensibilidad.
MALLA <- "media" # "gruesa" | "media" | "fina"

# Hiperparametros y opciones del ajuste.
PRIOR_RANGE <- c(9.3, 0.5)
PRIOR_SIGMA <- c(0.6, 0.02)
SPDE_ALPHA <- 2
SPDE_CONSTR <- FALSE
FAMILIA <- "gaussian"

# Discretizacion de las covariables con efecto no lineal en el segundo modelo.
# Se configura por variable para poder hacer analisis de sensibilidad despues.
N_GRUPOS_RW2 <- c(
    Temperatura = 15L,
    Velocidad_Viento = 15L
)
METODO_GRUPOS_RW2 <- "quantile"


CALCULAR_CPO <- TRUE
CALCULAR_LOSO <- FALSE # Cambiar a TRUE solo cuando se quiera ejecutar LOSO.
GUARDAR_MODELO <- TRUE

DIR_OUT <- here("outputs", "modelo", "Modelo_1", paste0("modelo_", ANIO))
DIR_MODELO <- here("data", "processed", "Modelos", "Modelo_1", paste0("modelo_", ANIO))
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MODELO, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. DATOS
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

df <- readRDS(ruta_datos)
setDT(df)

columnas_necesarias <- unique(c(
    "ESTACION", "FECHA", "X_km", "Y_km", RESPUESTA, COVARIABLES
))
columnas_ausentes <- setdiff(columnas_necesarias, names(df))
if (length(columnas_ausentes) > 0L) {
    stop("Faltan columnas en el maestro: ", paste(columnas_ausentes, collapse = ", "))
}

df[, FECHA := as.Date(FECHA)]
df <- df[FECHA >= FECHA_INICIO & FECHA <= FECHA_FIN]
df <- df[complete.cases(
    df[, c(RESPUESTA, COVARIABLES, "X_km", "Y_km"), with = FALSE]
)]
if (nrow(df) == 0L) stop("No quedan observaciones completas para ajustar el modelo.")

setorder(df, FECHA, ESTACION)
df[, y := get(RESPUESTA)]
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
coords <- as.matrix(df[, .(X_km, Y_km)])

# ==============================================================================
# 2. MALLA Y CAMPO ESPACIAL
# ==============================================================================
mesh <- readRDS(ruta_malla)
spde <- inla.spde2.pcmatern(
    mesh = mesh,
    alpha = SPDE_ALPHA,
    prior.range = PRIOR_RANGE,
    prior.sigma = PRIOR_SIGMA,
    constr = SPDE_CONSTR
)

indice_campo <- inla.spde.make.index(
    name = "campo_espacial",
    n.spde = spde$n.spde
)
A_campo <- inla.spde.make.A(mesh = mesh, loc = coords)

stopifnot(
    "A_campo y df tienen distinto numero de filas" = nrow(A_campo) == nrow(df),
    "A_campo y la malla tienen distinta dimension" = ncol(A_campo) == mesh$n
)


# ==============================================================================
# 3. AJUSTE
# ==============================================================================

inla_stack <- inla.stack(
    data = list(y_response = df$y),
    A = list(A_campo, 1),
    effects = list(
        indice_campo,
        data.frame(
            Intercept = rep(1, nrow(df)),
            as.data.frame(df[, ..COVARIABLES]),
            Temperatura_rw2 = df$Temperatura_rw2,
            Velocidad_Viento_rw2 = df$Velocidad_Viento_rw2
        )
    ),
    tag = "modelo_2025"
)


t0_modelo_lineal <- Sys.time()
cat("\nAjustando el modelo lineal INLA-SPDE...\n")
flush.console()
modelo_2025 <- inla(
    formula = y_response ~ -1 + Intercept +
        intensidad + Temperatura + Velocidad_Viento +
        f(campo_espacial, model = spde),
    data = inla.stack.data(inla_stack, spde = spde),
    family = FAMILIA, # verosimilitud de la respuesta
    verbose= TRUE,
    num.threads = 4, # Numero de hilos para paralelizar el ajuste
    control.predictor = list( # Controla como inla calcula las distribuciones a posterior del predictor final y sus fitted values
        A = inla.stack.A(inla_stack), # Matriz de proyeccion para el predictor final
        compute = TRUE # Para poder calcular los fitted values y sus marginals a posteriori
    ),
    control.inla = list(
        strategy = "gaussian", # Estrategia de aproximacion a posteriori
        int.strategy = "eb" # Estrategia de integracion numerica para los hiperparametros
    ),
    inla.mode = "experimental",# Modo experimental para usar la version mas reciente de inla

    control.compute = list(
        cpo = CALCULAR_CPO,
        dic = TRUE,
        waic = TRUE,
        openmp.strategy = "huge"
    )
)
minutos_modelo_lineal <- as.numeric(difftime(
    Sys.time(), t0_modelo_lineal, units = "mins"
))
cat(sprintf("Modelo lineal completado en %.2f minutos.\n", minutos_modelo_lineal))

t0_modelo_rw2 <- Sys.time()
cat("\nAjustando el modelo no lineal RW2 INLA-SPDE...\n")
flush.console()
modelo_2025_trata_varia <- inla(
    formula = y_response ~ -1 + Intercept +
        intensidad +
        f(Temperatura_rw2, model = "rw2", scale.model = TRUE) +
        f(Velocidad_Viento_rw2, model = "rw2", scale.model = TRUE) +
        f(campo_espacial, model = spde),
    data = inla.stack.data(inla_stack, spde = spde),
    family = FAMILIA, # verosimilitud de la respuesta
    verbose= TRUE,
    num.threads = 4, # Numero de hilos para paralelizar el ajuste
    control.predictor = list( # Controla como inla calcula las distribuciones a posterior del predictor final y sus fitted values
        A = inla.stack.A(inla_stack), # Matriz de proyeccion para el predictor final
        compute = TRUE # Para poder calcular los fitted values y sus marginals a posteriori
    ),
    control.inla = list(
        strategy = "gaussian", # Estrategia de aproximacion a posteriori
        int.strategy = "eb" # Estrategia de integracion numerica para los hiperparametros
    ),
    inla.mode = "experimental",# Modo experimental para usar la version mas reciente de inla
    control.compute = list(
        cpo = CALCULAR_CPO,
        dic = TRUE,
        waic = TRUE,
        openmp.strategy = "huge"
    )
)
minutos_modelo_rw2 <- as.numeric(difftime(
    Sys.time(), t0_modelo_rw2, units = "mins"
))
cat(sprintf("Modelo RW2 completado en %.2f minutos.\n", minutos_modelo_rw2))

# ==============================================================================
# 5. RESUMEN POSTERIOR Y VALIDACION
# ==============================================================================
resumir_modelo <- function(modelo, id_modelo, especificacion, minutos) {
    resultado_spde <- inla.spde2.result(
        inla = modelo,
        name = "campo_espacial",
        spde = spde,
        do.transf = TRUE
    )
    resumen_rango_modelo <- inla.zmarginal(
        resultado_spde$marginals.range.nominal[[1]],
        silent = TRUE
    )
    marginal_sigma_modelo <- inla.tmarginal(
        sqrt,
        resultado_spde$marginals.variance.nominal[[1]]
    )
    resumen_sigma_modelo <- inla.zmarginal(
        marginal_sigma_modelo,
        silent = TRUE
    )

    fila <- data.table(
        model = id_modelo,
        specification = especificacion,
        year = ANIO,
        scale = ESCALA,
        mesh = MALLA,
        vertices = mesh$n,
        observations = nrow(df),
        stations = uniqueN(df$ESTACION),
        DIC = modelo$dic$dic,
        DIC_p_eff = modelo$dic$p.eff,
        WAIC = modelo$waic$waic,
        WAIC_p_eff = modelo$waic$p.eff,
        range_mean_km = resumen_rango_modelo$mean,
        range_q025_km = resumen_rango_modelo$quant0.025,
        range_q975_km = resumen_rango_modelo$quant0.975,
        spatial_sd_mean = resumen_sigma_modelo$mean,
        spatial_sd_q025 = resumen_sigma_modelo$quant0.025,
        spatial_sd_q975 = resumen_sigma_modelo$quant0.975,
        minutes = minutos
    )

    if (CALCULAR_CPO) {
        cpo <- modelo$cpo$cpo
        cpo_valido <- is.finite(cpo) & cpo > 0
        fila[, `:=`(
            LCPO = mean(-log(cpo[cpo_valido])),
            PIT_mean = mean(modelo$cpo$pit, na.rm = TRUE),
            PIT_sd = sd(modelo$cpo$pit, na.rm = TRUE),
            CPO_failures = sum(modelo$cpo$failure != 0, na.rm = TRUE),
            CPO_invalid = sum(!cpo_valido)
        )]
    }
    fila
}

metricas_comparacion <- rbindlist(list(
    resumir_modelo(
        modelo_2025,
        "M1_lineal",
        "Temperatura y viento lineales",
        minutos_modelo_lineal
    ),
    resumir_modelo(
        modelo_2025_trata_varia,
        "M2_rw2",
        sprintf(
            "Temperatura y viento RW2 (%d y %d grupos %s)",
            N_GRUPOS_RW2[["Temperatura"]],
            N_GRUPOS_RW2[["Velocidad_Viento"]],
            METODO_GRUPOS_RW2
        ),
        minutos_modelo_rw2
    )
), use.names = TRUE, fill = TRUE)

# Se mantienen estos objetos con el modelo lineal para las tablas ya existentes.
spde_result <- inla.spde2.result(
    inla = modelo_2025,
    name = "campo_espacial",
    spde = spde,
    do.transf = TRUE
)
marginal_rango <- spde_result$marginals.range.nominal[[1]]
marginal_sigma <- inla.tmarginal(
    sqrt,
    spde_result$marginals.variance.nominal[[1]]
)
resumen_rango <- inla.zmarginal(marginal_rango, silent = TRUE)
resumen_sigma <- inla.zmarginal(marginal_sigma, silent = TRUE)
metricas <- copy(metricas_comparacion[model == "M1_lineal"])

# Validacion cruzada leave-one-station-out (LOSO). Para cada observacion, el
# grupo excluido contiene todas las observaciones de su estacion. A diferencia
# del CPO anterior, esta validacion condiciona la prediccion sin ninguna
# respuesta de la ubicacion evaluada. Para acelerar, los hiperparametros se fijan
# en la moda estimada por el modelo completo, como hace inla.group.cv.
predicciones_loso <- NULL
metricas_loso_estacion <- NULL
metricas_loso <- NULL
if (CALCULAR_LOSO) {
    if (tolower(FAMILIA) != "gaussian") {
        stop(
            "La validacion LOSO implementada requiere FAMILIA = 'gaussian' ",
            "para calcular predicciones y sus intervalos en la escala de la respuesta."
        )
    }

    indices_por_estacion <- split(
        seq_len(nrow(df)),
        factor(df$ESTACION, levels = unique(df$ESTACION))
    )
    grupos_loso <- unname(indices_por_estacion[as.character(df$ESTACION)])
    max_tamano_estacion <- max(lengths(indices_por_estacion))

    stopifnot(
        "Debe existir un grupo LOSO por observacion" =
            length(grupos_loso) == nrow(df),
        "Cada grupo LOSO debe contener la propia observacion" =
            all(vapply(
                seq_len(nrow(df)),
                function(i) i %in% grupos_loso[[i]],
                logical(1)
            ))
    )

    cat(sprintf(
        paste0(
            "\nCalculando LOSO de los dos modelos para %d estaciones ",
            "(hasta %d observaciones por grupo)...\n"
        ),
        length(indices_por_estacion), max_tamano_estacion
    ))
    flush.console()

    calcular_loso_modelo <- function(modelo, id_modelo) {
        cat("  ", id_modelo, "...\n", sep = "")
        t0_loso <- Sys.time()
        resultado_loso <- inla.group.cv(
            result = modelo,
            groups = grupos_loso,
            size.max = max_tamano_estacion,
            verbose = FALSE
        )
        minutos_loso <- as.numeric(difftime(
            Sys.time(), t0_loso, units = "mins"
        ))

        if (
            length(resultado_loso$mean) != nrow(df) ||
                length(resultado_loso$sd) != nrow(df)
        ) {
            stop(
                "INLA no ha devuelto una prediccion LOSO por observacion para ",
                id_modelo, "."
            )
        }

        # inla.group.cv devuelve la incertidumbre del predictor lineal. Se suma
        # la varianza residual para formar intervalos predictivos gaussianos.
        nombre_precision <- grep(
            "^Precision for the Gaussian observations",
            names(modelo$marginals.hyperpar),
            value = TRUE
        )
        if (length(nombre_precision) != 1L) {
            stop("No se ha podido identificar la precision residual de ", id_modelo)
        }
        varianza_residual <- inla.emarginal(
            function(precision) 1 / precision,
            modelo$marginals.hyperpar[[nombre_precision]]
        )
        sd_predictiva <- sqrt(resultado_loso$sd^2 + varianza_residual)
        error <- resultado_loso$mean - df$y
        limite_inferior <- resultado_loso$mean - 1.96 * sd_predictiva
        limite_superior <- resultado_loso$mean + 1.96 * sd_predictiva
        dentro_95 <- limite_inferior <= df$y & df$y <= limite_superior
        densidad <- resultado_loso$cv
        densidad_valida <- is.finite(densidad) & densidad > 0
        log_score_negativo <- rep(NA_real_, length(densidad))
        log_score_negativo[densidad_valida] <- -log(densidad[densidad_valida])

        predicciones <- df[, .(
            model = id_modelo,
            ESTACION,
            FECHA,
            observed = y
        )]
        predicciones[, `:=`(
            predicted_mean = resultado_loso$mean,
            predictor_sd = resultado_loso$sd,
            predictive_sd = sd_predictiva,
            lower_95 = limite_inferior,
            upper_95 = limite_superior,
            error = error,
            covered_95 = dentro_95,
            predictive_density = densidad,
            negative_log_score = log_score_negativo
        )]

        por_estacion <- predicciones[, .(
            observations = .N,
            RMSE = sqrt(mean(error^2, na.rm = TRUE)),
            MAE = mean(abs(error), na.rm = TRUE),
            Bias = mean(error, na.rm = TRUE),
            Coverage95 = 100 * mean(covered_95, na.rm = TRUE),
            MNLPD = mean(negative_log_score, na.rm = TRUE),
            density_failures = sum(is.na(negative_log_score))
        ), by = .(model, ESTACION)]

        globales <- data.table(
            model = id_modelo,
            RMSE_loso = sqrt(mean(error^2, na.rm = TRUE)),
            MAE_loso = mean(abs(error), na.rm = TRUE),
            Bias_loso = mean(error, na.rm = TRUE),
            Cov95_loso = 100 * mean(dentro_95, na.rm = TRUE),
            MNLPD_loso = mean(log_score_negativo, na.rm = TRUE),
            LOSO_density_failures = sum(is.na(log_score_negativo)),
            LOSO_minutes = minutos_loso
        )
        list(
            predictions = predicciones,
            by_station = por_estacion,
            global = globales
        )
    }

    resultados_loso <- list(
        calcular_loso_modelo(modelo_2025, "M1_lineal"),
        calcular_loso_modelo(modelo_2025_trata_varia, "M2_rw2")
    )
    predicciones_loso <- rbindlist(lapply(
        resultados_loso, `[[`, "predictions"
    ))
    metricas_loso_estacion <- rbindlist(lapply(
        resultados_loso, `[[`, "by_station"
    ))
    metricas_loso <- rbindlist(lapply(resultados_loso, `[[`, "global"))
    setorder(metricas_loso_estacion, model, ESTACION)

    metricas_comparacion <- merge(
        metricas_comparacion,
        metricas_loso,
        by = "model",
        all.x = TRUE,
        sort = FALSE
    )
    metricas <- copy(metricas_comparacion[model == "M1_lineal"])
}

# ==============================================================================
# 6. SALIDAS
# ==============================================================================
fwrite(metricas, file.path(DIR_OUT, "metricas_modelo_2025.csv"))
fwrite(
    metricas_comparacion,
    file.path(DIR_OUT, "comparacion_modelos_2025.csv")
)
if (CALCULAR_LOSO) {
    fwrite(
        predicciones_loso,
        file.path(DIR_OUT, "predicciones_loso_modelos_2025.csv")
    )
    fwrite(
        metricas_loso_estacion,
        file.path(DIR_OUT, "metricas_loso_por_estacion_modelos_2025.csv")
    )
}
efectos_fijos_comparacion <- rbindlist(list(
    as.data.table(
        modelo_2025$summary.fixed,
        keep.rownames = "term"
    )[, model := "M1_lineal"],
    as.data.table(
        modelo_2025_trata_varia$summary.fixed,
        keep.rownames = "term"
    )[, model := "M2_rw2"]
), use.names = TRUE, fill = TRUE)
setcolorder(efectos_fijos_comparacion, c("model", "term"))
fwrite(
    efectos_fijos_comparacion,
    file.path(DIR_OUT, "efectos_fijos_modelos_2025.csv")
)
hiperparametros_comparacion <- rbindlist(list(
    as.data.table(
        modelo_2025$summary.hyperpar,
        keep.rownames = "hyperparameter"
    )[, model := "M1_lineal"],
    as.data.table(
        modelo_2025_trata_varia$summary.hyperpar,
        keep.rownames = "hyperparameter"
    )[, model := "M2_rw2"]
), use.names = TRUE, fill = TRUE)
setcolorder(hiperparametros_comparacion, c("model", "hyperparameter"))
fwrite(
    hiperparametros_comparacion,
    file.path(DIR_OUT, "hiperparametros_modelos_2025.csv")
)
efectos_no_lineales <- rbindlist(list(
    as.data.table(
        modelo_2025_trata_varia$summary.random$Temperatura_rw2,
        keep.rownames = "level"
    )[, covariate := "Temperatura"],
    as.data.table(
        modelo_2025_trata_varia$summary.random$Velocidad_Viento_rw2,
        keep.rownames = "level"
    )[, covariate := "Velocidad_Viento"]
), use.names = TRUE, fill = TRUE)
setcolorder(efectos_no_lineales, c("covariate", "level", "ID"))
fwrite(
    efectos_no_lineales,
    file.path(DIR_OUT, "efectos_no_lineales_rw2_2025.csv")
)
configuracion <- list(
    year = ANIO,
    scale = ESCALA,
    start_date = FECHA_INICIO,
    end_date = FECHA_FIN,
    response = RESPUESTA,
    covariates = COVARIABLES,
    mesh = MALLA,
    prior_range = PRIOR_RANGE,
    prior_sigma = PRIOR_SIGMA,
    spde_alpha = SPDE_ALPHA,
    spde_constr = SPDE_CONSTR,
    family = FAMILIA,
    rw2_groups = N_GRUPOS_RW2,
    rw2_grouping_method = METODO_GRUPOS_RW2,
    calculate_cpo = CALCULAR_CPO,
    calculate_loso = CALCULAR_LOSO
)

if (GUARDAR_MODELO) {
    saveRDS(
        list(
            models = list(
                M1_lineal = modelo_2025,
                M2_rw2 = modelo_2025_trata_varia
            ),
            mesh = mesh,
            spde = spde,
            stack = inla_stack,
            data = df,
            configuration = configuracion,
            metrics = metricas_comparacion
        ),
        file.path(
            DIR_MODELO,
            sprintf("comparacion_modelos_%d_%s_malla_%s.rds", ANIO, ESCALA, MALLA)
        )
    )
}

cat("\n--- Comparacion de modelos completada ---\n")
print(metricas_comparacion)
cat("Resultados:", DIR_OUT, "\n")
if (GUARDAR_MODELO) cat("Modelo:", DIR_MODELO, "\n")

# ==============================================================================
# 7. TABLAS ACADEMICAS PARA EL TFM
# ==============================================================================
source(here("R", "utilities", "academic_quality_tables.R"))

DIR_TABLAS <- here(
    "outputs", "figures", "modelo", "Modelo_1", paste0("modelo_", ANIO)
)
dir.create(DIR_TABLAS, recursive = TRUE, showWarnings = FALSE)

formatear_numero <- function(x, decimales = 3L) {
    ifelse(
        is.finite(x),
        formatC(x, format = "f", digits = decimales, decimal.mark = "."),
        "--"
    )
}

# Efectos fijos: media, desviacion posterior e intervalo creible del 95%.
tabla_efectos_tfm <- as.data.table(
    modelo_2025$summary.fixed,
    keep.rownames = "termino_original"
)
etiquetas_terminos <- c(
    Intercept = "Intercept",
    intensidad = "Traffic intensity",
    Temperatura = "Temperature",
    Velocidad_Viento = "Wind speed"
)
tabla_efectos_tfm[, Term := fifelse(
    termino_original %chin% names(etiquetas_terminos),
    unname(etiquetas_terminos[termino_original]),
    gsub("_", " ", termino_original)
)]
tabla_efectos_tfm[, Mean := formatear_numero(mean)]
tabla_efectos_tfm[, SD := formatear_numero(sd)]
tabla_efectos_tfm[, `95% CrI` := sprintf(
    "[%s, %s]",
    formatear_numero(`0.025quant`),
    formatear_numero(`0.975quant`)
)]
tabla_efectos_tfm <- tabla_efectos_tfm[, .(Term, Mean, SD, `95% CrI`)]

ruta_tabla_efectos <- file.path(
    DIR_TABLAS,
    sprintf("table_fixed_effects_model_1_%d.png", ANIO)
)
booktabs_png(
    tabla_efectos_tfm,
    ruta_tabla_efectos,
    title = "Posterior estimates of the fixed effects",
    subtitle = sprintf(
        "Model 1: daily Gaussian INLA-SPDE model, Madrid, %d",
        ANIO
    ),
    note = paste0(
        "Posterior mean, standard deviation (SD), and 95% credible interval ",
        "(CrI). The response is log(1 + NO2), and continuous covariates are ",
        "standardized; coefficients therefore refer to a one-standard-deviation ",
        "increase in each covariate."
    ),
    widths = c(1.75, 0.80, 0.75, 1.55),
    align = c("left", "right", "right", "right"),
    font_size = 9,
    row_height = 0.28
)

# Parametros del campo espacial con la parametrizacion de los priors PC usada.
tabla_espacial_tfm <- data.table(
    Parameter = c("Effective range (km)", "Spatial SD"),
    `PC prior` = c(
        sprintf(
            "P(range < %.2f km) = %.2f",
            PRIOR_RANGE[1], PRIOR_RANGE[2]
        ),
        sprintf(
            "P(SD > %.2f) = %.2f",
            PRIOR_SIGMA[1], PRIOR_SIGMA[2]
        )
    ),
    Mean = c(
        formatear_numero(resumen_rango$mean),
        formatear_numero(resumen_sigma$mean)
    ),
    SD = c(
        formatear_numero(resumen_rango$sd),
        formatear_numero(resumen_sigma$sd)
    ),
    `95% CrI` = c(
        sprintf(
            "[%s, %s]",
            formatear_numero(resumen_rango$quant0.025),
            formatear_numero(resumen_rango$quant0.975)
        ),
        sprintf(
            "[%s, %s]",
            formatear_numero(resumen_sigma$quant0.025),
            formatear_numero(resumen_sigma$quant0.975)
        )
    )
)

ruta_tabla_espacial <- file.path(
    DIR_TABLAS,
    sprintf("table_spatial_field_model_1_%d.png", ANIO)
)
booktabs_png(
    tabla_espacial_tfm,
    ruta_tabla_espacial,
    title = "Posterior estimates of the spatial field",
    subtitle = sprintf(
        "Medium mesh (%d vertices) and Matern SPDE specification",
        mesh$n
    ),
    note = paste0(
        "The effective range is the distance at which spatial correlation ",
        "becomes approximately negligible. SD denotes the marginal standard ",
        "deviation of the latent spatial field. Priors are stated using the PC ",
        "prior probability interpretation."
    ),
    widths = c(1.65, 2.05, 0.72, 0.72, 1.52),
    align = c("left", "left", "right", "right", "right"),
    font_size = 9,
    row_height = 0.30
)

# Ajuste y capacidad predictiva. Las columnas no calculadas se muestran como --.
valor_metrica <- function(nombre) {
    if (nombre %chin% names(metricas_comparacion)) {
        metricas_comparacion[[nombre]]
    } else {
        rep(NA_real_, nrow(metricas_comparacion))
    }
}
escala_ingles <- c(HORARIO = "hourly", DIARIO = "daily", MENSUAL = "monthly")
pit_media <- valor_metrica("PIT_mean")
pit_sd <- valor_metrica("PIT_sd")
tabla_rendimiento_tfm <- data.table(
    Model = c("M1 linear", "M2 nonlinear RW2"),
    DIC = formatear_numero(valor_metrica("DIC"), 1L),
    WAIC = formatear_numero(valor_metrica("WAIC"), 1L),
    `Mean -log(CPO)` = formatear_numero(valor_metrica("LCPO")),
    `PIT mean (SD)` = fifelse(
        is.finite(pit_media) & is.finite(pit_sd),
        sprintf("%.3f (%.3f)", pit_media, pit_sd),
        "--"
    ),
    `CPO failures` = formatear_numero(valor_metrica("CPO_failures"), 0L)
)
anchos_rendimiento <- c(1.45, 0.82, 0.82, 1.22, 1.22, 1.02)
if (CALCULAR_LOSO) {
    tabla_rendimiento_tfm[, `:=`(
        `LOSO RMSE` = formatear_numero(valor_metrica("RMSE_loso")),
        `LOSO MAE` = formatear_numero(valor_metrica("MAE_loso")),
        `LOSO coverage` = paste0(
            formatear_numero(valor_metrica("Cov95_loso"), 1L),
            "%"
        )
    )]
    anchos_rendimiento <- c(anchos_rendimiento, 1.05, 1.05, 1.25)
}

ruta_tabla_rendimiento <- file.path(
    DIR_TABLAS,
    sprintf("table_model_comparison_%d.png", ANIO)
)
booktabs_png(
    tabla_rendimiento_tfm,
    ruta_tabla_rendimiento,
    title = "Linear versus nonlinear covariate effects",
    subtitle = sprintf(
        paste0(
            "%s data, %s observations, %d stations, %s mesh; ",
            "RW2 groups: %d temperature, %d wind (%s)"
        ),
        unname(escala_ingles[ESCALA]), format(nrow(df), big.mark = ","),
        uniqueN(df$ESTACION), MALLA,
        N_GRUPOS_RW2[["Temperatura"]],
        N_GRUPOS_RW2[["Velocidad_Viento"]],
        METODO_GRUPOS_RW2
    ),
    note = paste0(
        "Lower DIC, WAIC, and mean negative log-CPO indicate better performance ",
        "on the same observations. M1 uses linear temperature and wind effects; ",
        "M2 estimates both effects with scaled second-order random walks."
    ),
    widths = anchos_rendimiento,
    align = c("left", rep("right", ncol(tabla_rendimiento_tfm) - 1L)),
    font_size = 8.5,
    row_height = 0.30
)

cat("\nTablas academicas:\n")
cat("  Efectos fijos: ", ruta_tabla_efectos, "\n", sep = "")
cat("  Campo espacial: ", ruta_tabla_espacial, "\n", sep = "")
cat("  Rendimiento: ", ruta_tabla_rendimiento, "\n", sep = "")
