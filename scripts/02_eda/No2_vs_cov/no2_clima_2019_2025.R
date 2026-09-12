# ==============================================================================
# RELACION ENTRE NO2 Y COVARIABLES CLIMATICAS: 2019 FRENTE A 2025
# ==============================================================================
# Bloques del analisis:
#   1. Escala mensual: todas las covariables climaticas frente a NO2.
#   2. Escala diaria: precipitacion y viento frente a NO2 en dias de episodio.
#   3. Escala horaria: perfil medio de precipitacion/viento y NO2 dentro de esos
#      mismos dias de episodio.
#
# Se usan siempre las variables originales (*_raw), no las estandarizadas.
# No se imputan valores ausentes. Las medias se calculan con los valores
# disponibles y se conserva el numero de observaciones utilizado.
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(here)
    library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. Configuracion
# ------------------------------------------------------------------------------

ANIOS <- c(2019L, 2025L)
UMBRAL_LLUVIA_MM <- 0.1
# Se selecciona aproximadamente el 25 % de dias mas ventosos de cada ano.
CUANTIL_DIA_VENTOSO <- 0.75
RETARDOS_DIARIOS <- 0:7
RETARDOS_HORARIOS <- 0:24

RUTA_DIARIA <- here(
    "data", "processed", "Maestro", "diario",
    "dataset_maestro_inla_2019_2025_DIARIO.rds"
)
RUTAS_HORARIAS <- c(
    `2019` = here(
        "data", "processed", "Maestro", "2019",
        "dataset_maestro_inla_2019_HORARIO.rds"
    ),
    `2025` = here(
        "data", "processed", "Maestro", "horario",
        "dataset_maestro_inla_2025_HORARIO.rds"
    )
)

DIR_SALIDA <- here("outputs", "figures", "no2_clima_2019_2025")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

CLIMAVARS <- c(
    "Temperatura_raw",
    "Humedad_Relativa_raw",
    "Precipitaciones_raw",
    "Presion_Barometrica_raw",
    "Radiacion_Solar_raw",
    "Velocidad_Viento_raw"
)
METADATOS_CLIMA <- data.table(
    Variable = CLIMAVARS,
    Etiqueta = c(
        "Temperatura",
        "Humedad relativa",
        "Precipitaciones",
        "Presi\u00f3n barom\u00e9trica",
        "Radiaci\u00f3n solar",
        "Velocidad del viento"
    ),
    Unidad = c("\u00b0C", "%", "mm/d\u00eda", "hPa", "W/m\u00b2", "m/s")
)

COLORES_ANIO <- c(`2019` = "#D55E00", `2025` = "#0072B2")
RELLENOS_ANIO <- c(`2019` = "#D9D9D9", `2025` = "#8C8C8C")
MESES_ABREV <- c(
    "Ene", "Feb", "Mar", "Abr", "May", "Jun",
    "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"
)

media_disponible <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
}

cor_segura <- function(x, y, metodo) {
    completos <- is.finite(x) & is.finite(y)
    if (sum(completos) < 3L) return(NA_real_)
    x <- x[completos]
    y <- y[completos]
    if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
    stats::cor(x, y, method = metodo)
}

cor_segura <- function(x, y, metodo = "spearman") {
    completos <- is.finite(x) & is.finite(y)
    if (sum(completos) < 3L) return(NA_real_)
    x <- x[completos]
    y <- y[completos]
    if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
    unname(cor(x, y, method = metodo))
}

comprobar_columnas <- function(dt, columnas, origen) {
    faltantes <- setdiff(columnas, names(dt))
    if (length(faltantes) > 0L) {
        stop(
            "Faltan columnas en ", origen, ": ",
            paste(faltantes, collapse = ", ")
        )
    }
}

if (!file.exists(RUTA_DIARIA)) stop("No existe el maestro diario: ", RUTA_DIARIA)
rutas_ausentes <- RUTAS_HORARIAS[!file.exists(RUTAS_HORARIAS)]
if (length(rutas_ausentes) > 0L) {
    stop(
        "No se encuentran los maestros horarios: ",
        paste(rutas_ausentes, collapse = ", ")
    )
}

# ------------------------------------------------------------------------------
# 2. Datos diarios: una observacion media de Madrid por fecha
# ------------------------------------------------------------------------------

datos_diarios <- as.data.table(readRDS(RUTA_DIARIA))
columnas_diarias <- c("FECHA", "DATO_DIARIO", CLIMAVARS)
comprobar_columnas(datos_diarios, columnas_diarias, "el maestro diario")
datos_diarios[, FECHA := as.Date(FECHA)]
datos_diarios[, Anio := as.integer(format(FECHA, "%Y"))]
datos_diarios <- datos_diarios[Anio %in% ANIOS]

# Se promedian primero las estaciones para que cada fecha pese una sola vez.
diario_madrid <- datos_diarios[, c(
    list(NO2 = media_disponible(DATO_DIARIO), N_NO2 = sum(!is.na(DATO_DIARIO))),
    lapply(.SD, media_disponible)
), by = .(Anio, FECHA), .SDcols = CLIMAVARS]

diario_madrid[, Mes := as.integer(format(FECHA, "%m"))]
diario_madrid[, Umbral_viento := quantile(
    Velocidad_Viento_raw,
    probs = CUANTIL_DIA_VENTOSO,
    na.rm = TRUE,
    names = FALSE
), by = Anio]
diario_madrid[, Dia_lluvioso := is.finite(Precipitaciones_raw) &
    Precipitaciones_raw >= UMBRAL_LLUVIA_MM]
diario_madrid[, Dia_ventoso := is.finite(Velocidad_Viento_raw) &
    Velocidad_Viento_raw >= Umbral_viento]

resumen_eventos <- diario_madrid[, .(
    Dias_totales = .N,
    Dias_lluviosos = sum(Dia_lluvioso),
    Dias_ventosos = sum(Dia_ventoso),
    Umbral_lluvia_mm = UMBRAL_LLUVIA_MM,
    Umbral_viento_ms = unique(Umbral_viento)
), by = Anio]

# ------------------------------------------------------------------------------
# 3. Escala mensual: todas las covariables frente a NO2
# ------------------------------------------------------------------------------

mensual <- diario_madrid[, c(
    list(NO2 = media_disponible(NO2), N_dias_NO2 = sum(!is.na(NO2))),
    lapply(.SD, media_disponible)
), by = .(Anio, Mes), .SDcols = CLIMAVARS]

mensual_largo <- melt(
    mensual,
    id.vars = c("Anio", "Mes", "NO2", "N_dias_NO2"),
    measure.vars = CLIMAVARS,
    variable.name = "Variable",
    value.name = "Valor_clima"
)
mensual_largo <- merge(
    mensual_largo,
    METADATOS_CLIMA,
    by = "Variable",
    all.x = TRUE,
    sort = FALSE
)
mensual_largo[, Anio := factor(Anio, levels = ANIOS)]
setorder(mensual_largo, Variable, Anio, Mes)

correlacion_mensual <- mensual_largo[
    is.finite(Valor_clima) & is.finite(NO2),
    .(
        N_meses = .N,
        Cor_Pearson = cor(Valor_clima, NO2, method = "pearson"),
        Cor_Spearman = cor(Valor_clima, NO2, method = "spearman")
    ),
    by = .(Anio, Variable, Etiqueta)
]

grafico_mensual_variable <- function(variable) {
    d <- mensual_largo[Variable == variable]
    etiqueta <- unique(d$Etiqueta)
    unidad <- unique(d$Unidad)
    max_clima <- max(d$Valor_clima, na.rm = TRUE)
    max_no2 <- max(d$NO2, na.rm = TRUE)
    if (!is.finite(max_clima) || !is.finite(max_no2) ||
            max_clima <= 0 || max_no2 <= 0) {
        stop("No se puede construir el eje doble para ", variable)
    }
    factor_eje <- max_clima / max_no2

    ggplot(d, aes(x = Mes)) +
        geom_col(
            aes(y = Valor_clima, fill = Anio, group = Anio),
            position = position_dodge(width = 0.82),
            width = 0.38,
            color = "grey35",
            linewidth = 0.18,
            alpha = 0.72,
            na.rm = TRUE
        ) +
        geom_line(
            aes(y = NO2 * factor_eje, color = Anio, group = Anio),
            linewidth = 0.95,
            na.rm = TRUE
        ) +
        geom_point(
            aes(y = NO2 * factor_eje, color = Anio, group = Anio),
            size = 1.4,
            na.rm = TRUE
        ) +
        scale_x_continuous(
            breaks = 1:12,
            labels = MESES_ABREV,
            expand = expansion(mult = c(0.02, 0.02))
        ) +
        scale_y_continuous(
            name = paste0(etiqueta, " (", unidad, ")"),
            expand = expansion(mult = c(0, 0.08)),
            sec.axis = sec_axis(
                transform = ~ . / factor_eje,
                name = "NO\u2082 medio (\u00b5g/m\u00b3)"
            )
        ) +
        scale_fill_manual(values = RELLENOS_ANIO, name = "Clima") +
        scale_color_manual(values = COLORES_ANIO, name = "NO\u2082") +
        labs(title = etiqueta, x = NULL) +
        theme_bw(base_size = 9.5) +
        theme(
            plot.title = element_text(face = "bold", size = 10.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank(),
            legend.position = "bottom"
        )
}

graficos_mensuales <- lapply(CLIMAVARS, grafico_mensual_variable)
panel_mensual <- wrap_plots(graficos_mensuales, ncol = 2, guides = "collect") +
    plot_annotation(
        title = "NO\u2082 y covariables clim\u00e1ticas a escala mensual",
        subtitle = paste0(
            "Madrid, 2019 frente a 2025 | barras = clima mensual medio; ",
            "l\u00edneas = NO\u2082 mensual medio"
        ),
        caption = paste0(
            "Cada panel utiliza su propia escala. Las medias mensuales se obtienen ",
            "a partir de las medias diarias de Madrid; no se imputan datos ausentes."
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 15),
            plot.subtitle = element_text(size = 10.5),
            plot.caption = element_text(hjust = 0, color = "grey35")
        )
    ) &
    theme(legend.position = "bottom")

# Relacion directa: cada punto es un mes y el numero identifica el mes.
dispersion_mensual <- ggplot(
    mensual_largo,
    aes(x = Valor_clima, y = NO2, color = Anio)
) +
    geom_point(size = 2.0, alpha = 0.82, na.rm = TRUE) +
    geom_text(
        aes(label = Mes),
        size = 2.5,
        nudge_y = 0.6,
        show.legend = FALSE,
        check_overlap = TRUE,
        na.rm = TRUE
    ) +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        linewidth = 0.75,
        na.rm = TRUE
    ) +
    facet_wrap(~ Etiqueta, scales = "free_x", ncol = 2) +
    scale_color_manual(values = COLORES_ANIO, name = "A\u00f1o") +
    labs(
        title = "Relaci\u00f3n mensual entre NO\u2082 y clima",
        subtitle = "Cada punto es un mes; la recta resume la asociaci\u00f3n dentro de cada a\u00f1o",
        x = "Valor clim\u00e1tico mensual",
        y = "NO\u2082 medio mensual (\u00b5g/m\u00b3)",
        caption = "Asociaci\u00f3n descriptiva; no implica causalidad."
    ) +
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

# ------------------------------------------------------------------------------
# 4. Escala diaria: precipitacion y viento en dias de episodio
# ------------------------------------------------------------------------------

diario_eventos <- rbindlist(list(
    diario_madrid[Dia_lluvioso == TRUE, .(
        Anio, FECHA, Variable = "Precipitaciones",
        Valor_clima = Precipitaciones_raw, NO2
    )],
    diario_madrid[Dia_ventoso == TRUE, .(
        Anio, FECHA, Variable = "Velocidad del viento",
        Valor_clima = Velocidad_Viento_raw, NO2
    )]
))
diario_eventos[, Anio := factor(Anio, levels = ANIOS)]

correlacion_diaria_eventos <- diario_eventos[
    is.finite(Valor_clima) & is.finite(NO2),
    .(
        N_dias = .N,
        Cor_Pearson = cor(Valor_clima, NO2, method = "pearson"),
        Cor_Spearman = cor(Valor_clima, NO2, method = "spearman")
    ),
    by = .(Anio, Variable)
]

dispersion_diaria <- ggplot(
    diario_eventos,
    aes(x = Valor_clima, y = NO2, color = Anio)
) +
    geom_point(size = 1.8, alpha = 0.58, na.rm = TRUE) +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        linewidth = 0.85,
        alpha = 0.14,
        na.rm = TRUE
    ) +
    facet_wrap(
        ~ Variable,
        scales = "free_x",
        nrow = 1,
        labeller = as_labeller(c(
            Precipitaciones = paste0(
                "D\u00edas lluviosos (\u2265 ", UMBRAL_LLUVIA_MM, " mm/d\u00eda)"
            ),
            `Velocidad del viento` = paste0(
                "D\u00edas ventosos (percentil ",
                CUANTIL_DIA_VENTOSO * 100, " de cada a\u00f1o)"
            )
        ))
    ) +
    scale_color_manual(values = COLORES_ANIO, name = "A\u00f1o") +
    labs(
        title = "NO\u2082, precipitaci\u00f3n y viento a escala diaria",
        subtitle = "Solamente se muestran los d\u00edas clasificados como episodios",
        x = "Valor clim\u00e1tico diario",
        y = "NO\u2082 medio diario (\u00b5g/m\u00b3)",
        caption = "Las bandas corresponden al intervalo de confianza de la recta descriptiva."
    ) +
    theme_bw(base_size = 11.5) +
    theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

# ------------------------------------------------------------------------------
# 5. Escala horaria dentro de los dias de episodio
# ------------------------------------------------------------------------------

normalizar_horario <- function(dt) {
    alternativas <- c(
        "Velocidad Viento_raw" = "Velocidad_Viento_raw",
        "Velocidad Viento" = "Velocidad_Viento"
    )
    for (origen in names(alternativas)) {
        destino <- alternativas[[origen]]
        if (origen %in% names(dt) && !destino %in% names(dt)) {
            setnames(dt, origen, destino)
        }
    }
    dt
}

datos_horarios <- rbindlist(
    lapply(names(RUTAS_HORARIAS), function(anio) {
        dt <- normalizar_horario(as.data.table(readRDS(RUTAS_HORARIAS[[anio]])))
        columnas <- c(
            "FECHA", "HORA", "DATO",
            "Precipitaciones_raw", "Velocidad_Viento_raw"
        )
        comprobar_columnas(dt, columnas, paste0("el maestro horario de ", anio))
        dt <- dt[, ..columnas]
        dt[, Anio := as.integer(anio)]
        dt
    }),
    use.names = TRUE
)

datos_horarios[, FECHA := as.Date(FECHA)]
datos_horarios[, Hora := as.integer(HORA) - 1L]

horario_madrid <- datos_horarios[, .(
    NO2 = media_disponible(DATO),
    Precipitaciones = media_disponible(Precipitaciones_raw),
    Viento = media_disponible(Velocidad_Viento_raw),
    N_NO2 = sum(!is.na(DATO))
), by = .(Anio, FECHA, Hora)]

clasificacion_dias <- diario_madrid[, .(
    Anio, FECHA, Dia_lluvioso, Dia_ventoso, Umbral_viento
)]
horario_madrid <- merge(
    horario_madrid,
    clasificacion_dias,
    by = c("Anio", "FECHA"),
    all.x = TRUE,
    sort = FALSE
)

perfil_horario_eventos <- rbindlist(list(
    horario_madrid[Dia_lluvioso == TRUE, .(
        Valor_clima = media_disponible(Precipitaciones),
        NO2 = media_disponible(NO2),
        N_horas_clima = sum(!is.na(Precipitaciones)),
        N_horas_NO2 = sum(!is.na(NO2))
    ), by = .(Anio, Hora)][, Variable := "Precipitaciones"],
    horario_madrid[Dia_ventoso == TRUE, .(
        Valor_clima = media_disponible(Viento),
        NO2 = media_disponible(NO2),
        N_horas_clima = sum(!is.na(Viento)),
        N_horas_NO2 = sum(!is.na(NO2))
    ), by = .(Anio, Hora)][, Variable := "Velocidad del viento"]
), use.names = TRUE)

perfil_horario_eventos[, Anio := factor(Anio, levels = ANIOS)]
setorder(perfil_horario_eventos, Variable, Anio, Hora)

grafico_horario_evento <- function(variable) {
    d <- perfil_horario_eventos[Variable == variable]
    max_clima <- max(d$Valor_clima, na.rm = TRUE)
    max_no2 <- max(d$NO2, na.rm = TRUE)
    factor_eje <- max_clima / max_no2
    unidad <- if (variable == "Precipitaciones") "mm" else "m/s"

    ggplot(d, aes(x = Hora)) +
        geom_col(
            aes(y = Valor_clima, fill = Anio, group = Anio),
            position = position_dodge(width = 0.82),
            width = 0.38,
            color = "grey35",
            linewidth = 0.18,
            alpha = 0.72,
            na.rm = TRUE
        ) +
        geom_line(
            aes(y = NO2 * factor_eje, color = Anio, group = Anio),
            linewidth = 1.0,
            na.rm = TRUE
        ) +
        geom_point(
            aes(y = NO2 * factor_eje, color = Anio, group = Anio),
            size = 1.25,
            na.rm = TRUE
        ) +
        scale_x_continuous(
            breaks = c(0, 5, 10, 15, 20, 23),
            limits = c(-0.55, 23.55),
            expand = expansion(mult = c(0, 0))
        ) +
        scale_y_continuous(
            name = paste0(variable, " media (", unidad, ")"),
            expand = expansion(mult = c(0, 0.08)),
            sec.axis = sec_axis(
                transform = ~ . / factor_eje,
                name = "NO\u2082 medio (\u00b5g/m\u00b3)"
            )
        ) +
        scale_fill_manual(values = RELLENOS_ANIO, name = "Clima") +
        scale_color_manual(values = COLORES_ANIO, name = "NO\u2082") +
        labs(title = variable, x = "Hora del d\u00eda") +
        theme_bw(base_size = 10.5) +
        theme(
            plot.title = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            legend.position = "bottom"
        )
}

panel_horario <- wrap_plots(
    lapply(c("Precipitaciones", "Velocidad del viento"), grafico_horario_evento),
    ncol = 2,
    guides = "collect"
) +
    plot_annotation(
        title = "Perfil horario de NO\u2082 durante episodios de lluvia y viento",
        subtitle = paste0(
            "Madrid, 2019 frente a 2025 | lluvia \u2265 ", UMBRAL_LLUVIA_MM,
            " mm/d\u00eda; viento = 25 % de d\u00edas m\u00e1s ventosos de cada a\u00f1o"
        ),
        caption = paste0(
            "Cada hora es la media de todos los d\u00edas de episodio del a\u00f1o. ",
            "Barras = clima; l\u00edneas = NO\u2082."
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10),
            plot.caption = element_text(hjust = 0, color = "grey35")
        )
    ) &
    theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# 6. Mapas de calor para explorar retardos temporales
# ------------------------------------------------------------------------------
# Un retardo positivo k compara X(t-k) con NO2(t). Antes de calcular la
# correlacion se elimina el ciclo medio mensual (datos diarios) o el ciclo
# medio mes-hora (datos horarios). Asi se reduce el riesgo de confundir un
# patron estacional compartido con una relacion retardada.

diario_anomalias <- copy(diario_madrid)
diario_anomalias[, NO2_anomalia := NO2 - media_disponible(NO2), by = .(Anio, Mes)]
for (variable in CLIMAVARS) {
    nombre_anomalia <- paste0(variable, "_anomalia")
    diario_anomalias[, (nombre_anomalia) :=
        get(variable) - media_disponible(get(variable)), by = .(Anio, Mes)]
}

calcular_retardos_diarios <- function(variable) {
    nombre_anomalia <- paste0(variable, "_anomalia")
    rbindlist(lapply(RETARDOS_DIARIOS, function(retardo) {
        objetivo <- diario_anomalias[, .(Anio, FECHA, NO2_anomalia)]
        predictor <- diario_anomalias[, .(
            Anio,
            FECHA = FECHA + retardo,
            Clima_anomalia = get(nombre_anomalia)
        )]
        unidos <- merge(
            objetivo, predictor,
            by = c("Anio", "FECHA"), all = FALSE, sort = FALSE
        )
        resultado <- unidos[, .(
            N_pares = sum(is.finite(NO2_anomalia) & is.finite(Clima_anomalia)),
            Cor_Pearson = cor_segura(Clima_anomalia, NO2_anomalia, "pearson"),
            Cor_Spearman = cor_segura(Clima_anomalia, NO2_anomalia, "spearman")
        ), by = Anio]
        resultado[, `:=`(Variable = variable, Retardo = retardo)]
        resultado
    }))
}

retardos_diarios <- rbindlist(lapply(CLIMAVARS, calcular_retardos_diarios))
retardos_diarios <- merge(
    retardos_diarios,
    METADATOS_CLIMA[, .(Variable, Etiqueta)],
    by = "Variable",
    all.x = TRUE,
    sort = FALSE
)
retardos_diarios[, Maximo_abs := {
    valores <- abs(Cor_Spearman)
    if (any(is.finite(valores))) valores == max(valores, na.rm = TRUE) else rep(FALSE, .N)
}, by = .(Anio, Variable)]
retardos_diarios[, Etiqueta := factor(
    Etiqueta,
    levels = rev(METADATOS_CLIMA$Etiqueta)
)]

limite_diario <- max(abs(retardos_diarios$Cor_Spearman), na.rm = TRUE)
mapa_retardos_diarios <- ggplot(
    retardos_diarios,
    aes(x = Retardo, y = Etiqueta, fill = Cor_Spearman)
) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(
        aes(label = ifelse(
            is.finite(Cor_Spearman),
            paste0(sprintf("%.2f", Cor_Spearman), ifelse(Maximo_abs, " *", "")),
            ""
        )),
        size = 2.7
    ) +
    facet_grid(Anio ~ .) +
    scale_x_continuous(breaks = RETARDOS_DIARIOS) +
    scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
        limits = c(-limite_diario, limite_diario),
        name = "Correlaci\u00f3n\nSpearman"
    ) +
    labs(
        title = "Retardos diarios entre el clima y el NO₂",
        subtitle = "Correlaci\u00f3n de anomal\u00edas respecto al promedio mensual; * marca el mayor |r|",
        x = "Retardo de la covariable (d\u00edas)",
        y = NULL,
        caption = "k > 0 significa X(t-k) frente a NO\u2082(t). Es una exploraci\u00f3n de candidatos, no evidencia causal."
    ) +
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right"
    )

horario_anomalias <- copy(horario_madrid)
horario_anomalias[, Mes := as.integer(format(FECHA, "%m"))]
horario_anomalias[, Fecha_hora := as.POSIXct(FECHA, tz = "UTC") + Hora * 3600]
horario_anomalias[, NO2_anomalia :=
    NO2 - media_disponible(NO2), by = .(Anio, Mes, Hora)]
for (variable in c("Precipitaciones", "Viento")) {
    nombre_anomalia <- paste0(variable, "_anomalia")
    horario_anomalias[, (nombre_anomalia) :=
        get(variable) - media_disponible(get(variable)), by = .(Anio, Mes, Hora)]
}

calcular_retardos_horarios <- function(variable, columna_evento, etiqueta) {
    nombre_anomalia <- paste0(variable, "_anomalia")
    rbindlist(lapply(RETARDOS_HORARIOS, function(retardo) {
        objetivo <- horario_anomalias[
            get(columna_evento) == TRUE,
            .(Anio, Fecha_hora, NO2_anomalia)
        ]
        predictor <- horario_anomalias[, .(
            Anio,
            Fecha_hora = Fecha_hora + retardo * 3600,
            Clima_anomalia = get(nombre_anomalia)
        )]
        unidos <- merge(
            objetivo, predictor,
            by = c("Anio", "Fecha_hora"), all = FALSE, sort = FALSE
        )
        resultado <- unidos[, .(
            N_pares = sum(is.finite(NO2_anomalia) & is.finite(Clima_anomalia)),
            Cor_Pearson = cor_segura(Clima_anomalia, NO2_anomalia, "pearson"),
            Cor_Spearman = cor_segura(Clima_anomalia, NO2_anomalia, "spearman")
        ), by = Anio]
        resultado[, `:=`(
            Variable = variable,
            Etiqueta = etiqueta,
            Retardo = retardo
        )]
        resultado
    }))
}

retardos_horarios <- rbindlist(list(
    calcular_retardos_horarios(
        "Precipitaciones", "Dia_lluvioso", "Precipitaciones"
    ),
    calcular_retardos_horarios(
        "Viento", "Dia_ventoso", "Velocidad del viento"
    )
))
retardos_horarios[, Maximo_abs := {
    valores <- abs(Cor_Spearman)
    if (any(is.finite(valores))) valores == max(valores, na.rm = TRUE) else rep(FALSE, .N)
}, by = .(Anio, Variable)]
retardos_horarios[, Etiqueta := factor(
    Etiqueta,
    levels = rev(c("Precipitaciones", "Velocidad del viento"))
)]

limite_horario <- max(abs(retardos_horarios$Cor_Spearman), na.rm = TRUE)
mapa_retardos_horarios <- ggplot(
    retardos_horarios,
    aes(x = Retardo, y = Etiqueta, fill = Cor_Spearman)
) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(
        aes(label = ifelse(
            is.finite(Cor_Spearman),
            paste0(sprintf("%.2f", Cor_Spearman), ifelse(Maximo_abs, " *", "")),
            ""
        )),
        size = 2.25
    ) +
    facet_grid(Anio ~ .) +
    scale_x_continuous(breaks = seq(0, 24, by = 2)) +
    scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
        limits = c(-limite_horario, limite_horario),
        name = "Correlaci\u00f3n\nSpearman"
    ) +
    labs(
        title = "Retardos horarios de lluvia y viento frente al NO₂",
        subtitle = "Anomal\u00edas respecto al promedio de cada mes y hora; * marca el mayor |r|",
        x = "Retardo de la covariable (horas)",
        y = NULL,
        caption = paste0(
            "k > 0 significa X(t-k) frente a NO\u2082(t). Lluvia se eval\u00faa en d\u00edas lluviosos y viento en los ",
            "d\u00edas ventosos seleccionados."
        )
    ) +
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right"
    )

retardos_candidatos <- rbindlist(list(
    retardos_diarios[Maximo_abs == TRUE, .(
        Escala = "diaria", Anio, Variable = as.character(Etiqueta),
        Retardo, Cor_Pearson, Cor_Spearman, N_pares
    )],
    retardos_horarios[Maximo_abs == TRUE, .(
        Escala = "horaria", Anio, Variable = as.character(Etiqueta),
        Retardo, Cor_Pearson, Cor_Spearman, N_pares
    )]
), use.names = TRUE)

# ------------------------------------------------------------------------------
# 7. Guardar figuras y tablas
# ------------------------------------------------------------------------------

ruta_panel_mensual <- file.path(DIR_SALIDA, "01_panel_mensual_no2_clima_2019_2025.png")
ruta_dispersion_mensual <- file.path(DIR_SALIDA, "02_relacion_mensual_no2_clima_2019_2025.png")
ruta_dispersion_diaria <- file.path(DIR_SALIDA, "03_relacion_diaria_no2_lluvia_viento_2019_2025.png")
ruta_panel_horario <- file.path(DIR_SALIDA, "04_perfil_horario_no2_lluvia_viento_2019_2025.png")
ruta_retardos_diarios <- file.path(DIR_SALIDA, "05_mapa_calor_retardos_diarios_no2_clima_2019_2025.png")
ruta_retardos_horarios <- file.path(DIR_SALIDA, "06_mapa_calor_retardos_horarios_no2_lluvia_viento_2019_2025.png")

ggsave(ruta_panel_mensual, panel_mensual, width = 15, height = 15, dpi = 300, bg = "white")
ggsave(ruta_dispersion_mensual, dispersion_mensual, width = 13, height = 10, dpi = 300, bg = "white")
ggsave(ruta_dispersion_diaria, dispersion_diaria, width = 13, height = 6.5, dpi = 300, bg = "white")
ggsave(ruta_panel_horario, panel_horario, width = 14, height = 7.5, dpi = 300, bg = "white")
ggsave(ruta_retardos_diarios, mapa_retardos_diarios, width = 13, height = 9, dpi = 300, bg = "white")
ggsave(ruta_retardos_horarios, mapa_retardos_horarios, width = 15, height = 6.8, dpi = 300, bg = "white")

fwrite(mensual_largo, file.path(DIR_SALIDA, "datos_mensuales_no2_clima.csv"))
fwrite(correlacion_mensual, file.path(DIR_SALIDA, "correlaciones_mensuales_no2_clima.csv"))
fwrite(diario_eventos, file.path(DIR_SALIDA, "datos_diarios_eventos.csv"))
fwrite(correlacion_diaria_eventos, file.path(DIR_SALIDA, "correlaciones_diarias_eventos.csv"))
fwrite(perfil_horario_eventos, file.path(DIR_SALIDA, "perfil_horario_eventos.csv"))
fwrite(resumen_eventos, file.path(DIR_SALIDA, "resumen_dias_evento.csv"))
fwrite(retardos_diarios, file.path(DIR_SALIDA, "correlaciones_retardos_diarios.csv"))
fwrite(retardos_horarios, file.path(DIR_SALIDA, "correlaciones_retardos_horarios.csv"))
fwrite(retardos_candidatos, file.path(DIR_SALIDA, "retardos_candidatos.csv"))

cat("\nAnalisis climatico terminado.\n")
cat("Anios utilizados: ", paste(ANIOS, collapse = ", "), "\n", sep = "")
print(resumen_eventos)
cat("\nFiguras guardadas en: ", DIR_SALIDA, "\n", sep = "")
