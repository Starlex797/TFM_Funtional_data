# Campo meteorologico espacial para el modelo conjunto NO2 + meteorologia.
#
# Las mediciones diarias z(s, t) informan un unico campo espacial X(s):
#   z(s, t) = alpha_meteo + X(s) + error_meteo(s, t).
# El campo se proyecta directamente a las estaciones de NO2 con la matriz A
# y entra en el modelo de NO2 como beta * X(s) mediante copy. No se necesita
# un efecto iid auxiliar ni una pseudo-observacion de restriccion. No hay
# grupos temporales, replicas diarias del SPDE ni estructura AR1.


cargar_meteorologia_espacial <- function(archivos, anios, variable) {
    if (length(archivos) != length(anios) || any(!file.exists(archivos))) {
        stop(
            "Se necesita un archivo meteorologico existente por cada anio: ",
            paste(archivos[!file.exists(archivos)], collapse = ", ")
        )
    }

    piezas <- lapply(seq_along(archivos), function(k) {
        d <- data.table::as.data.table(readRDS(archivos[k]))
        requeridas <- c("ESTACION", "FECHA", "X_km", "Y_km", variable)
        if (!all(requeridas %in% names(d))) {
            stop("Columnas ausentes en ", archivos[k])
        }

        d[, FECHA := as.Date(FECHA)]
        d <- d[as.integer(format(FECHA, "%Y")) == anios[k]]
        if (!nrow(d)) stop("Archivo sin datos del anio esperado: ", anios[k])

        estado <- paste0(variable, "_estado")
        if (!estado %in% names(d)) stop("Falta indicador de calidad: ", estado)

        # Solo las mediciones observadas informan la verosimilitud meteorologica.
        # X_km/Y_km estan en ETRS89 / UTM 30N y en kilometros, como las mallas.
        d <- d[
            get(estado) == "OK" & is.finite(get(variable)) &
                is.finite(X_km) & is.finite(Y_km)
        ]
        if (!nrow(d)) {
            stop("Sin mediciones validas de ", variable, " en ", anios[k])
        }

        d[, .(ESTACION, FECHA, X_km, Y_km, valor = get(variable))]
    })

    meteo <- data.table::rbindlist(piezas)
    if (anyDuplicated(meteo[, .(ESTACION, FECHA)])) {
        stop("Duplicados estacion-fecha en los archivos meteorologicos.")
    }

    centro <- mean(meteo$valor)
    escala <- stats::sd(meteo$valor)
    if (!is.finite(escala) || escala <= 0) {
        stop("Covariable meteorologica constante.")
    }

    meteo[, z := (valor - centro) / escala]
    data.table::setorder(meteo, FECHA, ESTACION)

    list(
        datos = meteo,
        centro = centro,
        escala = escala,
        variable = variable,
        archivos = normalizePath(archivos, winslash = "/")
    )
}


preparar_campo_meteorologico <- function(meteo, datos_no2, prior, malla) {
    requeridas_meteo <- c("ESTACION", "X_km", "Y_km", "z")
    requeridas_no2 <- c("ESTACION", "X_km", "Y_km")
    if (!all(requeridas_meteo %in% names(meteo))) {
        stop("Faltan columnas para construir el campo meteorologico.")
    }
    if (!all(requeridas_no2 %in% names(datos_no2))) {
        stop("Faltan columnas espaciales en los datos de NO2.")
    }

    comprobar_coordenadas <- function(d, etiqueta) {
        variacion <- d[, .(
            rango_x = max(X_km) - min(X_km),
            rango_y = max(Y_km) - min(Y_km)
        ), by = ESTACION]
        if (any(variacion$rango_x > 1e-7 | variacion$rango_y > 1e-7)) {
            stop("Una estacion cambia de coordenadas en ", etiqueta, ".")
        }
    }
    comprobar_coordenadas(meteo, "meteorologia")
    comprobar_coordenadas(datos_no2, "NO2")

    sitios_no2 <- unique(datos_no2[, .(ESTACION, X_km, Y_km)])
    data.table::setorder(sitios_no2, ESTACION)

    xy_meteo <- as.matrix(meteo[, .(X_km, Y_km)])
    xy_no2 <- as.matrix(datos_no2[, .(X_km, Y_km)])

    spde <- INLA::inla.spde2.pcmatern(
        malla,
        alpha = 2,
        prior.range = prior$prior.range,
        prior.sigma = prior$prior.sigma
    )
    indice <- INLA::inla.spde.make.index(
        name = "campo_meteo",
        n.spde = spde$n.spde
    )

    # Una fila por medicion diaria, pero solo n.spde columnas: todos los dias
    # observados en una estacion informan el mismo campo X(s).
    A_meteo <- INLA::inla.spde.make.A(
        mesh = malla,
        loc = xy_meteo
    )
    # La copia tiene la misma dimension que el SPDE original. A_no2 proyecta
    # beta * X desde los nodos de la malla a cada observacion de contaminacion.
    A_no2 <- INLA::inla.spde.make.A(
        mesh = malla,
        loc = xy_no2
    )

    fuera_meteo <- abs(Matrix::rowSums(A_meteo) - 1) > 1e-7
    fuera_no2 <- abs(Matrix::rowSums(A_no2) - 1) > 1e-7
    if (any(fuera_meteo) || any(fuera_no2)) {
        stop(
            "La malla meteorologica no cubre todas las localizaciones. ",
            "Estaciones meteorologicas fuera: ",
            paste(unique(meteo$ESTACION[fuera_meteo]), collapse = ", "),
            " | Estaciones NO2 fuera: ",
            paste(unique(datos_no2$ESTACION[fuera_no2]), collapse = ", ")
        )
    }

    list(
        malla = malla,
        spde = spde,
        indice = indice,
        indice_copia = data.frame(campo_meteo_copia = seq_len(spde$n.spde)),
        A_meteo = A_meteo,
        A_no2 = A_no2,
        z = meteo$z,
        sitios_no2 = sitios_no2,
        n_estaciones_meteo = data.table::uniqueN(meteo$ESTACION),
        n_estaciones_no2 = nrow(sitios_no2),
        nombre_campo = "campo_meteo",
        nombre_copia = "campo_meteo_copia"
    )
}


crear_stack_conjunto_meteo <- function(y, A_no2, indice_no2, efectos_fijos,
                                       campo, tag = "estimacion") {
    n_no2 <- length(y)
    n_meteo <- length(campo$z)

    if (nrow(A_no2) != n_no2 || nrow(efectos_fijos) != n_no2) {
        stop("Respuesta, A_no2 y efectos fijos tienen dimensiones incompatibles.")
    }
    if (nrow(campo$A_meteo) != n_meteo || nrow(campo$A_no2) != n_no2) {
        stop("Matrices del campo meteorologico con dimensiones incompatibles.")
    }

    # 1) Modelo de NO2: A_no2 proyecta directamente la copia beta * X(s).
    stack_no2 <- INLA::inla.stack(
        data = list(
            y_response = cbind(y, NA_real_),
            link = rep(1L, n_no2)
        ),
        A = list(A_no2, campo$A_no2, 1),
        effects = list(indice_no2, campo$indice_copia, efectos_fijos),
        tag = tag
    )

    # 2) Mediciones meteorologicas:
    # z(s,t) = Intercept_meteo + campo_meteo(s) + error_meteo(s,t).
    stack_meteo <- INLA::inla.stack(
        data = list(
            y_response = cbind(NA_real_, campo$z),
            link = rep(2L, n_meteo)
        ),
        A = list(1, campo$A_meteo),
        effects = list(
            data.frame(Intercept_meteo = rep(1, n_meteo)),
            campo$indice
        ),
        tag = "meteo"
    )

    stack <- INLA::inla.stack(stack_no2, stack_meteo)
    attr(stack, "likelihood_no2") <- 1L
    stack
}


extraer_campo_en_no2 <- function(ajuste, datos, climatologia, campo) {
    resumen <- ajuste$modelo$summary.random[[campo$nombre_campo]]
    if (is.null(resumen) || nrow(resumen) != campo$spde$n.spde) {
        stop("INLA no ha devuelto el campo meteorologico esperado.")
    }

    # La media se puede proyectar exactamente de forma lineal. No se publican
    # intervalos marginales aproximados porque requeririan la covarianza conjunta
    # completa del campo, no solo las desviaciones tipicas de cada nodo.
    campo_estandarizado <- as.numeric(campo$A_no2 %*% resumen[, "mean"])
    intercepto_meteo <- ajuste$modelo$summary.fixed["Intercept_meteo", "mean"]

    salida <- data.table::copy(datos[, .(ESTACION, FECHA, X_km, Y_km)])
    salida[, `:=`(
        Variable = climatologia$variable,
        Campo_espacial_estandarizado = campo_estandarizado,
        Media_meteorologica_estimada = climatologia$centro +
            climatologia$escala * (intercepto_meteo + campo_estandarizado)
    )]
    salida
}


extraer_beta_campo <- function(modelo, variable) {
    candidatos <- grep(
        "^Beta for campo_meteo_copia",
        rownames(modelo$summary.hyperpar),
        value = TRUE
    )
    if (length(candidatos) != 1L) {
        stop("INLA no ha devuelto el coeficiente del campo meteorologico.")
    }

    b <- modelo$summary.hyperpar[candidatos, ]
    data.table::data.table(
        Modelo = "M2",
        Variable = variable,
        Coeficiente = b[["mean"]],
        IC95 = sprintf("[%.4f, %.4f]", b[["0.025quant"]], b[["0.975quant"]]),
        Significativa_95 = if (
            b[["0.025quant"]] > 0 || b[["0.975quant"]] < 0
        ) {
            "Si"
        } else {
            "No"
        }
    )
}
