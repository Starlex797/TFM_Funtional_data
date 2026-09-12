# ==============================================================================
# PERFIL HORARIO MEDIO DE TRAFICO Y NO2: 2019 FRENTE A 2025
# ==============================================================================
# El grafico reproduce la estructura de la figura de referencia:
#   - columnas: invierno y verano;
#   - filas: dias laborables y fines de semana;
#   - barras: intensidad media de trafico (eje izquierdo);
#   - lineas: concentracion media de NO2 (eje derecho).
#
# Se utilizan los maestros HORARIOS y las variables en su escala original:
#   DATO           -> NO2 en microgramos/m3
#   intensidad_raw -> intensidad de trafico en vehiculos/hora
#
# Convenciones:
#   HORA = 1 representa 00:00-00:59, por eso Hora = HORA - 1.
#   Laborable = lunes-viernes; fin de semana = sabado-domingo.
#   Invierno = diciembre-febrero; verano = junio-agosto.
#   No se imputan datos. Cada media se calcula con los valores disponibles.
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(here)
})

# ------------------------------------------------------------------------------
# 1. Configuracion
# ------------------------------------------------------------------------------

RUTAS_MAESTROS <- c(
    `2019` = here(
        "data", "processed", "Maestro", "2019",
        "dataset_maestro_inla_2019_HORARIO.rds"
    ),
    `2025` = here(
        "data", "processed", "Maestro", "horario",
        "dataset_maestro_inla_2025_HORARIO.rds"
    )
)

DIR_SALIDA <- here("outputs", "figures", "no2_trafico_2019_2025")
dir.create(DIR_SALIDA, recursive = TRUE, showWarnings = FALSE)

rutas_ausentes <- RUTAS_MAESTROS[!file.exists(RUTAS_MAESTROS)]
if (length(rutas_ausentes) > 0L) {
    stop(
        "No se encuentran los maestros horarios: ",
        paste(rutas_ausentes, collapse = ", ")
    )
}

# ------------------------------------------------------------------------------
# 2. Carga de los dos anos
# ------------------------------------------------------------------------------

datos <- rbindlist(
    lapply(names(RUTAS_MAESTROS), function(anio) {
        dt <- as.data.table(readRDS(RUTAS_MAESTROS[[anio]]))
        columnas <- c("FECHA", "HORA", "DATO", "intensidad_raw")
        faltantes <- setdiff(columnas, names(dt))

        if (length(faltantes) > 0L) {
            stop(
                "Faltan columnas en el maestro de ", anio, ": ",
                paste(faltantes, collapse = ", ")
            )
        }

        dt <- dt[, ..columnas]
        dt[, Anio := anio]
        dt
    }),
    use.names = TRUE
)

datos[, FECHA := as.Date(FECHA)]
datos[, Hora := as.integer(HORA) - 1L]
datos[, Mes := as.integer(format(FECHA, "%m"))]
datos[, Dia_semana := as.integer(format(FECHA, "%u"))]
datos[, Tipo_dia := fifelse(
    Dia_semana <= 5L,
    "D\u00edas laborables",
    "Fin de semana"
)]
datos[, Temporada := fcase(
    Mes %in% c(12L, 1L, 2L), "Invierno",
    Mes %in% c(6L, 7L, 8L), "Verano",
    default = NA_character_
)]

# Solamente se excluyen primavera y otono porque no forman parte de la
# comparacion solicitada. Los NA de trafico y NO2 permanecen en los datos.
datos_estaciones <- datos[!is.na(Temporada)]

media_disponible <- function(x) {
    if (all(is.na(x))) {
        return(NA_real_)
    }
    mean(x, na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# 3. Perfil horario medio
# ------------------------------------------------------------------------------

perfil_horario <- datos_estaciones[, .(
    Intensidad_media = media_disponible(intensidad_raw),
    NO2_medio = media_disponible(DATO),
    N_trafico = sum(!is.na(intensidad_raw)),
    N_NO2 = sum(!is.na(DATO))
), by = .(Anio, Temporada, Tipo_dia, Hora)]

perfil_horario[, Anio := factor(Anio, levels = c("2019", "2025"))]
perfil_horario[, Temporada := factor(
    Temporada,
    levels = c("Invierno", "Verano")
)]
perfil_horario[, Tipo_dia := factor(
    Tipo_dia,
    levels = c("D\u00edas laborables", "Fin de semana")
)]
setorder(perfil_horario, Tipo_dia, Temporada, Anio, Hora)

if (!any(is.finite(perfil_horario$Intensidad_media)) ||
        !any(is.finite(perfil_horario$NO2_medio))) {
    stop("No hay valores suficientes de trafico y NO2 para construir el perfil.")
}

# ggplot2 necesita una transformacion comun para representar el NO2 en el eje
# derecho. Se usa un unico factor para todos los paneles y anos, de forma que la
# comparacion visual entre ellos conserve la misma escala.
factor_eje <- max(perfil_horario$Intensidad_media, na.rm = TRUE) /
    max(perfil_horario$NO2_medio, na.rm = TRUE)

# ------------------------------------------------------------------------------
# 4. Grafico
# ------------------------------------------------------------------------------

colores_no2 <- c(`2019` = "#D55E00", `2025` = "#0072B2")
rellenos_trafico <- c(`2019` = "#D9D9D9", `2025` = "#8C8C8C")

grafico <- ggplot(perfil_horario, aes(x = Hora)) +
    geom_col(
        aes(y = Intensidad_media, fill = Anio, group = Anio),
        position = position_dodge(width = 0.82),
        width = 0.38,
        color = "grey35",
        linewidth = 0.18,
        alpha = 0.72,
        na.rm = TRUE
    ) +
    geom_line(
        aes(
            y = NO2_medio * factor_eje,
            color = Anio,
            group = Anio
        ),
        linewidth = 1.05,
        na.rm = TRUE
    ) +
    geom_point(
        aes(
            y = NO2_medio * factor_eje,
            color = Anio,
            group = Anio
        ),
        size = 1.35,
        na.rm = TRUE
    ) +
    facet_grid(Tipo_dia ~ Temporada) +
    scale_x_continuous(
        breaks = c(0, 5, 10, 15, 20, 23),
        limits = c(-0.55, 23.55),
        expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
        name = "Intensidad media de tr\u00e1fico (veh/h)",
        expand = expansion(mult = c(0, 0.08)),
        sec.axis = sec_axis(
            transform = ~ . / factor_eje,
            name = "NO\u2082 medio (\u00b5g/m\u00b3)"
        )
    ) +
    scale_fill_manual(
        values = rellenos_trafico,
        name = "Tr\u00e1fico"
    ) +
    scale_color_manual(
        values = colores_no2,
        name = "NO\u2082"
    ) +
    guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 2)
    ) +
    labs(
        title = "Tr\u00e1fico y NO\u2082 por tipo de d\u00eda y estaci\u00f3n del a\u00f1o",
        subtitle = paste0(
            "Perfil horario medio de Madrid: 2019 frente a 2025 | ",
            "barras = tr\u00e1fico; l\u00edneas = NO\u2082"
        ),
        x = "Hora del d\u00eda",
        caption = paste0(
            "Medias calculadas con los datos disponibles de todas las estaciones. ",
            "Laborables = lunes-viernes; no se excluyen festivos."
        )
    ) +
    theme_bw(base_size = 11.5) +
    theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10.5),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey94", color = "grey65"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey86", linewidth = 0.35),
        legend.position = "bottom",
        legend.box = "horizontal",
        axis.title.y.right = element_text(color = "grey20"),
        plot.caption = element_text(hjust = 0, color = "grey35")
    )

# ------------------------------------------------------------------------------
# 5. Salidas
# ------------------------------------------------------------------------------

ruta_png <- file.path(
    DIR_SALIDA,
    "perfil_horario_trafico_no2_invierno_verano_2019_2025.png"
)
ruta_csv <- file.path(
    DIR_SALIDA,
    "perfil_horario_trafico_no2_invierno_verano_2019_2025.csv"
)

ggsave(
    filename = ruta_png,
    plot = grafico,
    width = 13,
    height = 8.7,
    dpi = 300,
    bg = "white"
)
fwrite(perfil_horario, ruta_csv)

cat("\nGrafico guardado en: ", ruta_png, "\n", sep = "")
cat("Datos agregados guardados en: ", ruta_csv, "\n", sep = "")
cat("Factor comun entre los ejes: ", round(factor_eje, 3), "\n", sep = "")
