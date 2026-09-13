# ==============================================================================
# DIAGNOSTICO DE AUTOCORRELACION TEMPORAL EN LOS RESIDUOS
# ==============================================================================
# Sirve para decidir si un modelo puramente espacial necesita ademas una
# estructura temporal (AR1, RW1) antes de plantear el modelo espacio-temporal.
# La regla es la del proyecto: se ajusta primero el modelo espacial, se miran
# los residuos y el objetivo es RUIDO BLANCO; si queda autocorrelacion, se
# justifica el termino temporal.
#
# TRES DECISIONES METODOLOGICAS QUE CONDICIONAN TODO EL ARCHIVO
# -------------------------------------------------------------
# 1. Los datos son un PANEL (estacion x dia), no una sola serie. Aplicar acf()
#    al vector de residuos ordenado por (FECHA, ESTACION) mezclaria estaciones
#    distintas en el mismo desfase y el resultado no significaria nada. Por eso
#    cada estacion se trata como su propia serie.
#
# 2. Hay dias perdidos (complete.cases elimina filas). acf(), pacf() y
#    Box.test() suponen observaciones EQUIESPACIADAS: si se les pasa la serie
#    comprimida, el desfase k deja de ser "k dias". Por eso cada serie se
#    reindexa sobre una rejilla diaria completa con NA en los huecos. Los tres
#    procedimientos aceptan NA con na.action = na.pass, y Box.test() ya usa
#    internamente na.pass y n = numero de observaciones no NA.
#
# 3. Se calculan DOS diagnosticos porque responden a preguntas distintas y
#    llevan a terminos INLA distintos:
#      - modo "por_estacion": autocorrelacion DENTRO de cada estacion. Si
#        aparece, apunta a un campo espacio-temporal
#        f(campo, model = spde, group = tiempo, control.group = list(model = "ar1")).
#      - modo "media_diaria": autocorrelacion de la media diaria de los
#        residuos, es decir, la senal temporal COMPARTIDA por toda la ciudad.
#        Si aparece, apunta a un efecto temporal comun f(ID_TIEMPO, "ar1").
#    Lo habitual es mirar primero "media_diaria": es mucho mas barato de
#    modelar y suele explicar la mayor parte de la estructura.
# ==============================================================================


#' Residuos de un modelo INLA: observado menos media ajustada.
#'
#' @param modelo objeto devuelto por INLA::inla().
#' @param y_observada vector de respuesta observada.
#' @param indices_observaciones filas de summary.fitted.values que corresponden
#'   a las observaciones (ver indices_stack_inla()).
#' @return vector numerico de residuos.
residuos_inla <- function(modelo, y_observada, indices_observaciones) {
    ajustado <- modelo$summary.fitted.values[
        indices_observaciones, "mean"
    ]
    if (length(ajustado) != length(y_observada)) {
        stop("La respuesta y los valores ajustados tienen distinta longitud.")
    }
    y_observada - ajustado
}


#' Coloca una serie diaria sobre una rejilla completa, con NA en los huecos.
#'
#' Es el paso que hace que el desfase k signifique "k dias" y no "k
#' observaciones disponibles".
#'
#' @param fechas vector Date.
#' @param valores vector numerico de la misma longitud.
#' @param fecha_min,fecha_max extremos de la rejilla. Por defecto, los de fechas.
#' @return vector numerico indexado por dia consecutivo.
serie_diaria_regular <- function(fechas, valores,
                                 fecha_min = NULL, fecha_max = NULL) {
    fechas <- as.Date(fechas)
    if (length(fechas) != length(valores)) {
        stop("fechas y valores deben tener la misma longitud.")
    }
    if (is.null(fecha_min)) fecha_min <- min(fechas, na.rm = TRUE)
    if (is.null(fecha_max)) fecha_max <- max(fechas, na.rm = TRUE)

    rejilla <- seq(fecha_min, fecha_max, by = "day")
    serie <- rep(NA_real_, length(rejilla))
    posicion <- match(fechas, rejilla)
    conservar <- !is.na(posicion)
    serie[posicion[conservar]] <- valores[conservar]
    serie
}


#' ACF y PACF de una unica serie ya regularizada.
#'
#' pacf() resuelve las ecuaciones de Yule-Walker sobre la ACF estimada; con
#' muchos huecos esa matriz puede no ser definida positiva y la llamada falla.
#' En ese caso se devuelve NA en la columna PACF en lugar de interrumpir todo
#' el diagnostico.
#'
#' @param serie vector numerico, posiblemente con NA.
#' @param lag_max desfase maximo.
#' @return data.table con Lag, ACF, PACF, N_validas y Banda_95.
acf_pacf_serie <- function(serie, lag_max = 30L) {
    n_validas <- sum(!is.na(serie))
    if (n_validas <= lag_max + 1L) {
        return(data.table::data.table())
    }

    correlaciones <- stats::acf(
        serie,
        lag.max = lag_max,
        plot = FALSE,
        na.action = stats::na.pass
    )$acf[-1L]

    parciales <- tryCatch(
        as.numeric(stats::pacf(
            serie,
            lag.max = lag_max,
            plot = FALSE,
            na.action = stats::na.pass
        )$acf),
        error = function(e) rep(NA_real_, lag_max)
    )

    data.table::data.table(
        Lag = seq_len(lag_max),
        ACF = as.numeric(correlaciones),
        PACF = parciales,
        N_validas = n_validas,
        # Banda de ruido blanco: fuera de +/- 1.96/sqrt(n) el desfase es
        # incompatible con ausencia de autocorrelacion.
        Banda_95 = 1.96 / sqrt(n_validas)
    )
}


#' ACF y PACF de los residuos de un panel estacion x dia.
#'
#' @param residuos vector de residuos.
#' @param fechas vector Date de la misma longitud.
#' @param estaciones identificador de estacion de la misma longitud. Solo se
#'   utiliza cuando modo = "por_estacion".
#' @param modo "por_estacion" (una serie por estacion, luego se promedia) o
#'   "media_diaria" (media de los residuos de cada dia, una sola serie).
#' @param lag_max desfase maximo en dias.
#' @return lista con:
#'   - resumen: data.table con un desfase por fila. En "por_estacion" incluye
#'     la media entre estaciones y Prop_fuera_banda, la proporcion de
#'     estaciones cuya ACF supera su banda en ese desfase.
#'   - por_estacion: detalle por estacion (NULL en modo "media_diaria").
acf_pacf_residuos <- function(residuos, fechas, estaciones = NULL,
                              modo = c("por_estacion", "media_diaria"),
                              lag_max = 30L) {
    modo <- match.arg(modo)
    fechas <- as.Date(fechas)
    if (length(residuos) != length(fechas)) {
        stop("residuos y fechas deben tener la misma longitud.")
    }

    if (modo == "media_diaria") {
        panel <- data.table::data.table(FECHA = fechas, residuo = residuos)
        diaria <- panel[
            , .(residuo = mean(residuo, na.rm = TRUE)),
            by = FECHA
        ][order(FECHA)]

        resumen <- acf_pacf_serie(
            serie_diaria_regular(diaria$FECHA, diaria$residuo),
            lag_max = lag_max
        )
        if (nrow(resumen) == 0L) {
            stop("La serie diaria es demasiado corta para lag_max = ", lag_max)
        }
        resumen[, Modo := "media_diaria"]
        return(list(resumen = resumen, por_estacion = NULL))
    }

    if (is.null(estaciones)) {
        stop("El modo 'por_estacion' necesita el vector estaciones.")
    }
    if (length(estaciones) != length(residuos)) {
        stop("estaciones y residuos deben tener la misma longitud.")
    }

    # La rejilla es comun a todas las estaciones para que el desfase k sea el
    # mismo dia de calendario en todas ellas.
    fecha_min <- min(fechas, na.rm = TRUE)
    fecha_max <- max(fechas, na.rm = TRUE)

    panel <- data.table::data.table(
        ESTACION = as.character(estaciones),
        FECHA = fechas,
        residuo = residuos
    )
    data.table::setorder(panel, ESTACION, FECHA)

    detalle <- data.table::rbindlist(lapply(
        split(panel, by = "ESTACION", keep.by = TRUE),
        function(bloque) {
            tabla <- acf_pacf_serie(
                serie_diaria_regular(
                    bloque$FECHA, bloque$residuo,
                    fecha_min = fecha_min, fecha_max = fecha_max
                ),
                lag_max = lag_max
            )
            if (nrow(tabla) == 0L) return(data.table::data.table())
            tabla[, ESTACION := bloque$ESTACION[1L]]
            tabla
        }
    ), use.names = TRUE, fill = TRUE)

    if (nrow(detalle) == 0L) {
        stop("Ninguna estacion tiene serie suficiente para lag_max = ", lag_max)
    }

    resumen <- detalle[, .(
        ACF = mean(ACF, na.rm = TRUE),
        PACF = mean(PACF, na.rm = TRUE),
        # Cuantas estaciones, en proporcion, superan su propia banda en este
        # desfase. Es mas informativo que la media: una media pequena puede
        # esconder que todas las estaciones estan fuera de banda.
        Prop_fuera_banda = mean(abs(ACF) > Banda_95, na.rm = TRUE),
        Banda_95 = mean(Banda_95, na.rm = TRUE),
        N_estaciones = .N
    ), by = Lag][order(Lag)]
    resumen[, Modo := "por_estacion"]

    list(resumen = resumen, por_estacion = detalle[])
}


#' Prueba de Ljung-Box sobre los residuos de un panel estacion x dia.
#'
#' H0: los residuos son ruido blanco hasta el desfase indicado. Un p-valor
#' pequeno significa que queda autocorrelacion y que el termino temporal esta
#' justificado.
#'
#' En modo "por_estacion" se contrasta cada estacion por separado y los
#' p-valores se corrigen por comparaciones multiples: con ~24 estaciones, una
#' o dos podrian bajar de 0.05 solo por azar.
#'
#' @param residuos,fechas,estaciones igual que en acf_pacf_residuos().
#' @param modo "por_estacion" o "media_diaria".
#' @param lags vector de desfases a contrastar. 7 capta el ciclo semanal,
#'   30 el mes.
#' @param metodo_ajuste metodo de p.adjust() para el modo "por_estacion".
#' @return data.table con un contraste por fila.
ljung_box_residuos <- function(residuos, fechas, estaciones = NULL,
                               modo = c("por_estacion", "media_diaria"),
                               lags = c(1L, 7L, 14L, 30L),
                               metodo_ajuste = "holm") {
    modo <- match.arg(modo)
    fechas <- as.Date(fechas)
    if (length(residuos) != length(fechas)) {
        stop("residuos y fechas deben tener la misma longitud.")
    }

    contrastar <- function(serie, lag) {
        n_validas <- sum(!is.na(serie))
        if (n_validas <= lag + 1L) {
            return(list(estadistico = NA_real_, p_valor = NA_real_,
                        n = n_validas))
        }
        prueba <- stats::Box.test(serie, lag = lag, type = "Ljung-Box")
        list(
            estadistico = unname(prueba$statistic),
            p_valor = prueba$p.value,
            n = n_validas
        )
    }

    if (modo == "media_diaria") {
        panel <- data.table::data.table(FECHA = fechas, residuo = residuos)
        diaria <- panel[
            , .(residuo = mean(residuo, na.rm = TRUE)),
            by = FECHA
        ][order(FECHA)]
        serie <- serie_diaria_regular(diaria$FECHA, diaria$residuo)

        return(data.table::rbindlist(lapply(lags, function(lag) {
            resultado <- contrastar(serie, lag)
            data.table::data.table(
                Modo = "media_diaria",
                Lag = lag,
                Estadistico_Q = resultado$estadistico,
                P_valor = resultado$p_valor,
                N_validas = resultado$n,
                Autocorrelacion = data.table::fifelse(
                    is.na(resultado$p_valor), NA_character_,
                    data.table::fifelse(resultado$p_valor < 0.05, "Si", "No")
                )
            )
        })))
    }

    if (is.null(estaciones)) {
        stop("El modo 'por_estacion' necesita el vector estaciones.")
    }

    fecha_min <- min(fechas, na.rm = TRUE)
    fecha_max <- max(fechas, na.rm = TRUE)
    panel <- data.table::data.table(
        ESTACION = as.character(estaciones),
        FECHA = fechas,
        residuo = residuos
    )
    data.table::setorder(panel, ESTACION, FECHA)

    resultado <- data.table::rbindlist(lapply(
        split(panel, by = "ESTACION", keep.by = TRUE),
        function(bloque) {
            serie <- serie_diaria_regular(
                bloque$FECHA, bloque$residuo,
                fecha_min = fecha_min, fecha_max = fecha_max
            )
            data.table::rbindlist(lapply(lags, function(lag) {
                prueba <- contrastar(serie, lag)
                data.table::data.table(
                    Modo = "por_estacion",
                    ESTACION = bloque$ESTACION[1L],
                    Lag = lag,
                    Estadistico_Q = prueba$estadistico,
                    P_valor = prueba$p_valor,
                    N_validas = prueba$n
                )
            }))
        }
    ), use.names = TRUE, fill = TRUE)

    resultado[, P_valor_ajustado := stats::p.adjust(P_valor, method = metodo_ajuste),
        by = Lag
    ]
    resultado[, Autocorrelacion := data.table::fifelse(
        is.na(P_valor_ajustado), NA_character_,
        data.table::fifelse(P_valor_ajustado < 0.05, "Si", "No")
    )]
    data.table::setorder(resultado, Lag, ESTACION)
    resultado[]
}


#' Diagnostico temporal completo de los residuos de un modelo espacial.
#'
#' Ejecuta ACF, PACF y Ljung-Box en los dos modos y escribe por consola una
#' lectura orientativa. La decision final sobre el termino temporal es del
#' analista: esta funcion resume la evidencia, no elige el modelo.
#'
#' @param residuos,fechas,estaciones el panel de residuos.
#' @param lag_max desfase maximo de ACF y PACF.
#' @param lags_ljung desfases contrastados con Ljung-Box.
#' @param etiqueta nombre del modelo, solo para los mensajes.
#' @return lista con acf_por_estacion, acf_media_diaria, ljung_por_estacion y
#'   ljung_media_diaria.
diagnostico_temporal_residuos <- function(residuos, fechas, estaciones,
                                          lag_max = 30L,
                                          lags_ljung = c(1L, 7L, 14L, 30L),
                                          etiqueta = "modelo espacial") {
    acf_estacion <- acf_pacf_residuos(
        residuos, fechas, estaciones,
        modo = "por_estacion", lag_max = lag_max
    )
    acf_diaria <- acf_pacf_residuos(
        residuos, fechas,
        modo = "media_diaria", lag_max = lag_max
    )
    ljung_estacion <- ljung_box_residuos(
        residuos, fechas, estaciones,
        modo = "por_estacion", lags = lags_ljung
    )
    ljung_diaria <- ljung_box_residuos(
        residuos, fechas,
        modo = "media_diaria", lags = lags_ljung
    )

    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("DIAGNOSTICO TEMPORAL DE LOS RESIDUOS: ", etiqueta, "\n", sep = "")
    cat(strrep("=", 72), "\n", sep = "")

    cat("\nACF/PACF dentro de cada estacion (primeros 10 desfases):\n")
    print(acf_estacion$resumen[
        Lag <= 10L,
        .(Lag, ACF = round(ACF, 4), PACF = round(PACF, 4),
          Prop_fuera_banda = round(Prop_fuera_banda, 3))
    ])

    cat("\nACF/PACF de la media diaria (primeros 10 desfases):\n")
    print(acf_diaria$resumen[
        Lag <= 10L,
        .(Lag, ACF = round(ACF, 4), PACF = round(PACF, 4),
          Banda_95 = round(Banda_95, 4))
    ])

    cat("\nLjung-Box sobre la media diaria:\n")
    print(ljung_diaria[, .(Lag, Estadistico_Q = round(Estadistico_Q, 2),
                           P_valor = signif(P_valor, 4), Autocorrelacion)])

    cat("\nLjung-Box por estacion (resumen tras ajuste de p-valores):\n")
    print(ljung_estacion[!is.na(Autocorrelacion), .(
        Estaciones = .N,
        Con_autocorrelacion = sum(Autocorrelacion == "Si"),
        Proporcion = round(mean(Autocorrelacion == "Si"), 3)
    ), by = Lag][order(Lag)])

    cat("\nComo leerlo:\n")
    cat("  - ACF que decae de forma geometrica y PACF que se corta en el\n")
    cat("    desfase 1: estructura AR(1).\n")
    cat("  - Pico en el desfase 7 y sus multiplos: ciclo semanal. Se corrige\n")
    cat("    con un efecto de dia de la semana, no con un AR.\n")
    cat("  - ACF que decae muy despacio y casi lineal: tendencia o\n")
    cat("    estacionalidad sin modelar; un AR1 la taparia sin resolverla.\n")
    cat("  - Ljung-Box significativo en la media diaria: hay senal temporal\n")
    cat("    compartida por toda la ciudad; empieza por f(ID_TIEMPO, 'ar1').\n")

    list(
        acf_por_estacion = acf_estacion,
        acf_media_diaria = acf_diaria,
        ljung_por_estacion = ljung_estacion,
        ljung_media_diaria = ljung_diaria
    )
}


#' Grafico de ACF y PACF con la banda de ruido blanco.
#'
#' @param resumen data.table devuelto en $resumen por acf_pacf_residuos().
#' @param titulo titulo del grafico.
#' @return objeto ggplot.
graficar_acf_pacf <- function(resumen, titulo = "ACF y PACF de los residuos") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("graficar_acf_pacf() necesita el paquete ggplot2.")
    }

    largo <- data.table::melt(
        resumen[, .(Lag, ACF, PACF, Banda_95)],
        id.vars = c("Lag", "Banda_95"),
        variable.name = "Funcion",
        value.name = "Valor"
    )

    ggplot2::ggplot(largo, ggplot2::aes(x = Lag, y = Valor)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
        ggplot2::geom_hline(
            ggplot2::aes(yintercept = Banda_95),
            linetype = "dashed", colour = "steelblue"
        ) +
        ggplot2::geom_hline(
            ggplot2::aes(yintercept = -Banda_95),
            linetype = "dashed", colour = "steelblue"
        ) +
        ggplot2::geom_segment(
            ggplot2::aes(xend = Lag, yend = 0), linewidth = 0.6
        ) +
        ggplot2::facet_wrap(~Funcion, ncol = 1, scales = "free_y") +
        ggplot2::labs(
            title = titulo,
            x = "Desfase (dias)",
            y = NULL
        ) +
        ggplot2::theme_minimal(base_size = 11)
}
