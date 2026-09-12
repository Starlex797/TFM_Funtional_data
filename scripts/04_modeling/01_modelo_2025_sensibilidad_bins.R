# ==============================================================================
# SENSIBILIDAD AL NUMERO DE GRUPOS PARA RW2 - 2025
# ==============================================================================
# Este script se ejecuta por separado de 01_modelo_2025.R. Para cada covariable
# sustituye solo su efecto lineal por RW2 y mantiene las demas lineales.
# Cuando se elijan los grupos, deben copiarse manualmente a N_GRUPOS_RW2 del
# script principal.
# ==============================================================================

library(INLA)
library(data.table)
library(ggplot2)
library(here)


# ==============================================================================
# 0. CONFIGURACION EDITABLE
# ==============================================================================
ANIO <- 2025
ESCALA <- "DIARIO"
RESPUESTA <- "LOG_NO2_DIARIO"
FECHA_INICIO <- as.Date("2025-01-01")
FECHA_FIN <- as.Date("2025-12-31")

COVARIABLES <- c(
    "Velocidad_Viento",
    "Radiacion_Solar",
    "Humedad_Relativa",
    "Presion_Barometrica",
    "intensidad",
    "Precipitaciones"
)

BINS_A_EVALUAR <- c(10L, 20L, 30L, 50L, 75L, 100L)
VARIABLES_SENSIBILIDAD <- COVARIABLES
FORMAS_SENSIBILIDAD <- "RW2"
METODO_GRUPOS <- "quantile"
TOLERANCIA_DELTA_WAIC <- 2

MALLA <- "media"
PRIOR_RANGE <- c(9.3, 0.5)
PRIOR_SIGMA <- c(0.6, 0.02)
SPDE_ALPHA <- 2
SPDE_CONSTR <- FALSE
FAMILIA <- "gaussian"
NUM_THREADS <- 4L
VERBOSE_INLA <- FALSE


DIR_OUT <- here(
    "outputs", "modelo", "Modelo_1", paste0("modelo_", ANIO),
    "sensibilidad_bins"
)
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
DIR_FIGURAS <- here(
    "outputs", "figures", "modelo", "Modelo_1", paste0("modelo_", ANIO),
    "sensibilidad_bins"
)
dir.create(DIR_FIGURAS, recursive = TRUE, showWarnings = FALSE)
ruta_resultados <- file.path(
    DIR_OUT,
    sprintf("sensibilidad_bins_%d.csv", ANIO)
)


# ==============================================================================
# 1. DATOS Y COMPONENTE ESPACIAL COMUN
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
    "ESTACION", "FECHA", "X_km", "Y_km", RESPUESTA, COVARIABLES
))
columnas_ausentes <- setdiff(columnas_necesarias, names(df))
if (length(columnas_ausentes) > 0L) {
    stop("Faltan columnas: ", paste(columnas_ausentes, collapse = ", "))
}

variables_ausentes <- setdiff(VARIABLES_SENSIBILIDAD, COVARIABLES)
if (length(variables_ausentes) > 0L) {
    stop(
        "VARIABLES_SENSIBILIDAD contiene variables no incluidas: ",
        paste(variables_ausentes, collapse = ", ")
    )
}
if (!all(FORMAS_SENSIBILIDAD %in% c("RW1", "RW2"))) {
    stop("FORMAS_SENSIBILIDAD solo puede contener RW1 y/o RW2.")
}
if (any(BINS_A_EVALUAR < 3L) || anyDuplicated(BINS_A_EVALUAR)) {
    stop("BINS_A_EVALUAR debe contener valores unicos de al menos 3.")
}

df[, FECHA := as.Date(FECHA)]
if (!is.null(FECHA_INICIO)) df <- df[FECHA >= FECHA_INICIO]
if (!is.null(FECHA_FIN)) df <- df[FECHA <= FECHA_FIN]
df <- df[complete.cases(df[, ..columnas_necesarias])]
if (nrow(df) == 0L) stop("No quedan observaciones completas.")

setorder(df, FECHA, ESTACION)
df[, y := get(RESPUESTA)]

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
coords <- as.matrix(df[, .(X_km, Y_km)])
A_campo <- inla.spde.make.A(mesh = mesh, loc = coords)

efectos_base <- data.frame(Intercept = rep(1, nrow(df)))
for (variable in COVARIABLES) {
    efectos_base[[variable]] <- df[[variable]]
}

config_id <- paste(
    ANIO,
    ESCALA,
    FECHA_INICIO,
    FECHA_FIN,
    RESPUESTA,
    MALLA,
    nrow(df),
    paste(COVARIABLES, collapse = "+"),
    METODO_GRUPOS,
    sep = "|"
)


# ==============================================================================
# 2. FUNCIONES DE AJUSTE Y GUARDADO
# ==============================================================================
construir_formula <- function(variable, forma) {
    terminos_lineales <- setdiff(COVARIABLES, variable)
    termino_suave <- sprintf(
        "f(x_grupo, model = '%s', scale.model = TRUE)",
        tolower(forma)
    )
    as.formula(paste(
        "y_response ~ -1 + Intercept +",
        paste(
            c(
                terminos_lineales,
                termino_suave,
                "f(campo_espacial, model = spde)"
            ),
            collapse = " + "
        )
    ))
}

media_finita <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
}

calcular_deltas <- function(tabla) {
    tabla <- copy(tabla)
    columnas_delta <- intersect(c("Delta_WAIC", "Delta_DIC", "Delta_CPO"), names(tabla))
    if (length(columnas_delta) > 0L) tabla[, (columnas_delta) := NULL]

    tabla[, Delta_WAIC := {
        minimo <- if (any(is.finite(WAIC))) min(WAIC, na.rm = TRUE) else NA_real_
        WAIC - minimo
    }, by = .(Covariable, Forma)]
    tabla[, Delta_DIC := {
        minimo <- if (any(is.finite(DIC))) min(DIC, na.rm = TRUE) else NA_real_
        DIC - minimo
    }, by = .(Covariable, Forma)]
    tabla[, Delta_CPO := {
        minimo <- if (any(is.finite(CPO))) min(CPO, na.rm = TRUE) else NA_real_
        CPO - minimo
    }, by = .(Covariable, Forma)]
    setorder(tabla, Covariable, Forma, bins_solicitados)
    tabla
}

guardar_resultados <- function(tabla) {
    tabla_salida <- calcular_deltas(tabla)
    fwrite(tabla_salida, ruta_resultados)
    tabla_salida
}


# ==============================================================================
# 3. BUCLE DE SENSIBILIDAD
# ==============================================================================
tabla_resultados <- data.table(
    config_id = character(),
    Covariable = character(),
    Forma = character(),
    bins_solicitados = integer(),
    bins_reales = integer(),
    WAIC = numeric(),
    DIC = numeric(),
    CPO = numeric(),
    p_efectivo = numeric(),
    CPO_fallos = integer(),
    minutos = numeric(),
    estado = character()
)
if (REANUDAR && file.exists(ruta_resultados)) {
    tabla_anterior <- fread(ruta_resultados)
    if ("config_id" %in% names(tabla_anterior)) {
        config_actual <- config_id
        tabla_resultados <- tabla_anterior[
            config_id == config_actual &
                Covariable %in% VARIABLES_SENSIBILIDAD &
                Forma %in% FORMAS_SENSIBILIDAD &
                bins_solicitados %in% BINS_A_EVALUAR
        ]
    }
}

numero_total <- length(VARIABLES_SENSIBILIDAD) *
    length(FORMAS_SENSIBILIDAD) * length(BINS_A_EVALUAR)
cat(
    "\nObservaciones: ", nrow(df),
    " | estaciones: ", uniqueN(df$ESTACION),
    "\nCombinaciones solicitadas: ", numero_total,
    "\nResultados: ", ruta_resultados, "\n",
    sep = ""
)

for (variable in VARIABLES_SENSIBILIDAD) {
    for (forma in FORMAS_SENSIBILIDAD) {
        for (nb in BINS_A_EVALUAR) {
            ya_calculado <- nrow(tabla_resultados[
                Covariable == variable &
                    Forma == forma &
                    bins_solicitados == nb &
                    estado == "OK"
            ]) > 0L
            if (ya_calculado) {
                cat("Omitido (ya calculado): ", variable, " ", forma, " n=", nb, "\n", sep = "")
                next
            }

            x_grupo <- inla.group(
                df[[variable]],
                n = nb,
                method = METODO_GRUPOS
            )
            bins_reales <- uniqueN(x_grupo)
            minimo_grupos <- if (forma == "RW1") 2L else 3L

            cat(
                "Ajustando ", variable, " ", forma,
                " | solicitados=", nb, " | reales=", bins_reales, "...\n",
                sep = ""
            )

            if (bins_reales < minimo_grupos) {
                resultado <- data.table(
                    config_id = config_id,
                    Covariable = variable,
                    Forma = forma,
                    bins_solicitados = nb,
                    bins_reales = bins_reales,
                    WAIC = NA_real_,
                    DIC = NA_real_,
                    CPO = NA_real_,
                    p_efectivo = NA_real_,
                    CPO_fallos = NA_integer_,
                    minutos = 0,
                    estado = "Soporte insuficiente"
                )
            } else {
                efectos <- efectos_base
                efectos$x_grupo <- x_grupo
                stack <- inla.stack(
                    data = list(y_response = df$y),
                    A = list(A_campo, 1),
                    effects = list(indice_campo, efectos),
                    tag = "estimacion"
                )

                inicio <- Sys.time()
                resultado <- tryCatch({
                    fit <- inla(
                        formula = construir_formula(variable, forma),
                        data = inla.stack.data(stack, spde = spde),
                        family = FAMILIA,
                        verbose = VERBOSE_INLA,
                        num.threads = NUM_THREADS,
                        control.predictor = list(
                            A = inla.stack.A(stack),
                            compute = TRUE
                        ),
                        control.inla = list(
                            strategy = "gaussian",
                            int.strategy = "eb"
                        ),
                        inla.mode = "experimental",
                        control.compute = list(
                            cpo = TRUE,
                            dic = TRUE,
                            waic = TRUE,
                            openmp.strategy = "huge"
                        )
                    )
                    indices <- inla.stack.index(stack, tag = "estimacion")$data
                    cpo <- fit$cpo$cpo[indices]
                    cpo_valido <- is.finite(cpo) & cpo > 0

                    data.table(
                        config_id = config_id,
                        Covariable = variable,
                        Forma = forma,
                        bins_solicitados = nb,
                        bins_reales = bins_reales,
                        WAIC = fit$waic$waic,
                        DIC = fit$dic$dic,
                        CPO = media_finita(-log(cpo[cpo_valido])),
                        p_efectivo = fit$waic$p.eff,
                        CPO_fallos = sum(!cpo_valido) +
                            sum(fit$cpo$failure[indices] != 0, na.rm = TRUE),
                        minutos = as.numeric(difftime(
                            Sys.time(), inicio, units = "mins"
                        )),
                        estado = "OK"
                    )
                }, error = function(error) {
                    data.table(
                        config_id = config_id,
                        Covariable = variable,
                        Forma = forma,
                        bins_solicitados = nb,
                        bins_reales = bins_reales,
                        WAIC = NA_real_,
                        DIC = NA_real_,
                        CPO = NA_real_,
                        p_efectivo = NA_real_,
                        CPO_fallos = NA_integer_,
                        minutos = as.numeric(difftime(
                            Sys.time(), inicio, units = "mins"
                        )),
                        estado = conditionMessage(error)
                    )
                })
            }

            tabla_resultados <- tabla_resultados[
                !(Covariable == variable &
                    Forma == forma &
                    bins_solicitados == nb)
            ]
            tabla_resultados <- rbindlist(
                list(tabla_resultados, resultado),
                fill = TRUE
            )
            tabla_resultados <- guardar_resultados(tabla_resultados)
            rm(x_grupo, resultado)
            if (exists("fit")) rm(fit)
            if (exists("stack")) rm(stack)
            gc(verbose = FALSE)
        }
    }
}

tabla_resultados <- guardar_resultados(tabla_resultados)


# ==============================================================================
# 4. GRAFICAS DE CODO POR COVARIABLE
# ==============================================================================
guardar_png_seguro <- function(grafico, ruta, width = 10, height = 7, dpi = 300) {
    temporal <- tempfile(
        pattern = "codo_",
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

    ruta_actual <- ruta
    copiado <- suppressWarnings(file.copy(temporal, ruta_actual, overwrite = TRUE))
    if (!copiado) {
        ruta_actual <- paste0(
            tools::file_path_sans_ext(ruta),
            "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"
        )
        if (!file.copy(temporal, ruta_actual, overwrite = FALSE)) {
            stop("No se ha podido guardar la grafica en: ", dirname(ruta))
        }
        message("El PNG anterior estaba abierto. Nuevo archivo: ", ruta_actual)
    }
    ruta_actual
}

datos_graficas <- tabla_resultados[
    Forma == "RW2" & estado == "OK" &
        is.finite(WAIC) & is.finite(DIC) & is.finite(CPO)
]
rutas_graficas_codo <- character()
bins_recomendados <- data.table()

if (nrow(datos_graficas) == 0L) {
    warning("No hay resultados RW2 validos para construir las graficas de codo.")
} else {
    # Regla de parsimonia: entre los ajustes practicamente equivalentes al
    # WAIC minimo (Delta_WAIC <= 2), se recomienda el de menor resolucion.
    bins_recomendados <- datos_graficas[
        Delta_WAIC <= TOLERANCIA_DELTA_WAIC
    ][order(Covariable, bins_reales, bins_solicitados), .SD[1L], by = Covariable]
    bins_recomendados <- bins_recomendados[, .(
        Covariable,
        bins_solicitados_recomendados = bins_solicitados,
        bins_reales_recomendados = bins_reales,
        WAIC,
        DIC,
        CPO,
        p_efectivo,
        Delta_WAIC
    )]
    bins_recomendados[, Estado := fifelse(
        bins_solicitados_recomendados == max(BINS_A_EVALUAR),
        "No estabilizado: el minimo esta en el limite evaluado",
        "Codo compatible con estabilizacion"
    )]
    fwrite(
        bins_recomendados,
        file.path(
            DIR_OUT,
            sprintf("bins_recomendados_rw2_%d.csv", ANIO)
        )
    )

    datos_graficas_largos <- melt(
        datos_graficas,
        id.vars = c(
            "Covariable", "bins_solicitados", "bins_reales"
        ),
        measure.vars = c("WAIC", "DIC", "CPO", "p_efectivo"),
        variable.name = "Metrica",
        value.name = "Valor"
    )
    datos_graficas_largos[, Metrica := factor(
        Metrica,
        levels = c("WAIC", "DIC", "CPO", "p_efectivo"),
        labels = c("WAIC", "DIC", "-log(CPO) medio", "Parametros efectivos")
    )]

    rutas_graficas_codo <- setNames(vapply(
        unique(datos_graficas_largos$Covariable),
        function(variable) {
            datos_variable <- datos_graficas_largos[Covariable == variable]
            recomendacion <- bins_recomendados[Covariable == variable]
            texto_recomendacion <- if (
                recomendacion$bins_solicitados_recomendados ==
                    max(BINS_A_EVALUAR)
            ) {
                sprintf(
                    paste(
                        "No se observa un codo: el minimo sigue en n=%d",
                        "(%d grupos reales)"
                    ),
                    recomendacion$bins_solicitados_recomendados,
                    recomendacion$bins_reales_recomendados
                )
            } else {
                sprintf(
                    paste(
                        "Candidato: n solicitado = %d, n real = %d",
                        "(linea magenta; Delta WAIC <= %.1f)"
                    ),
                    recomendacion$bins_solicitados_recomendados,
                    recomendacion$bins_reales_recomendados,
                    TOLERANCIA_DELTA_WAIC
                )
            }
            grafico <- ggplot(
                datos_variable,
                aes(x = bins_reales, y = Valor, group = 1)
            ) +
                geom_vline(
                    xintercept = recomendacion$bins_reales_recomendados,
                    color = "#CC79A7",
                    linewidth = 0.8,
                    linetype = "dotted"
                ) +
                geom_line(color = "#0072B2", linewidth = 0.9) +
                geom_point(color = "#0072B2", size = 2.3) +
                geom_text(
                    aes(label = paste0("n=", bins_solicitados)),
                    vjust = -0.7,
                    size = 2.8,
                    check_overlap = TRUE
                ) +
                facet_wrap(~Metrica, scales = "free_y", ncol = 2) +
                scale_x_continuous(breaks = sort(unique(datos_variable$bins_reales))) +
                labs(
                    title = paste(
                        "Sensibilidad RW2:",
                        gsub("_", " ", variable)
                    ),
                    subtitle = texto_recomendacion,
                    x = "Numero real de grupos",
                    y = NULL,
                    caption = paste(
                        "WAIC, DIC y -log(CPO): menor es mejor.",
                        "Los parametros efectivos describen la complejidad."
                    )
                ) +
                theme_minimal(base_size = 11) +
                theme(
                    panel.grid.minor = element_blank(),
                    strip.text = element_text(face = "bold"),
                    plot.caption = element_text(hjust = 0)
                )

            guardar_png_seguro(
                grafico,
                file.path(
                    DIR_FIGURAS,
                    sprintf("codo_rw2_%s_%d.png", variable, ANIO)
                )
            )
        },
        character(1)
    ), unique(datos_graficas_largos$Covariable))
}

cat("\n--- Sensibilidad completada ---\n")
print(tabla_resultados)
cat("\nTabla: ", ruta_resultados, "\n", sep = "")
if (length(rutas_graficas_codo) > 0L) {
    cat("Graficas de codo: ", DIR_FIGURAS, "\n", sep = "")
}
