# ==============================================================================
# ANÁLISIS DESCRIPTIVO — NO₂ (Datos crudos)
# Madrid 2019-2025 · Escala diaria · Resumen anual y detalle de 2025
# ===============================================================================
# Proporciona la primera caracterización cuantitativa del NO₂.
# Permite comparar estaciones, detectar diferencias de nivel y variabilidad,
# revisar la cobertura de datos e identificar distribuciones asimétricas o
# valores extremos.
# ===============================================================================

# Finalidad. carga las concentraciones diarias de NO2, sin transformación logarítmica.
# Calcula por estación el número de días válidos, media y mediana, desviación típica, min y mac y asímetría y curtosis.
# Resultados se guardan en descriptivo_No2
# ==============================================================================
# Carga las concentraciones diarias de NO₂ de 2025, sin transformación logarítmica.
# Calcula por estación:número de días válidos;
# media y mediana;
# desviación típica;
# mínimo y máximo;
# asimetría y curtosis.

# Muestra los resultados en consola.
# Genera una tabla formateada y coloreada.
# Guarda los resultados por estacion de 2025 y la Tabla 3.1 anual en CSV y PNG.
# ==============================================================================
library(data.table)
library(here)
library(moments) # skewness(), kurtosis()
library(gt) # tabla formateada → PNG

# Directorio comun de resultados
carpeta_out <- here(
  "outputs", "figures", "EDA", "NO2", "Distribucion_global", "descriptivo_NO2"
)
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# BLOQUE 1: Estadísticas descriptivas básicas de NO2 por estación
# ==============================================================================

# Carga de datos diarios NO2 2025 (valores crudos, sin transformación log)
dt_no2 <- readRDS(here(
  "data", "processed", "contaminacion", "diario",
  "aire_madrid_2025_No2_trans_diarios1.rds"
))

# Estadísticas por estación sobre todas las observaciones diarias válidas
desc_no2 <- dt_no2[!is.na(DATO_DIARIO), .(
  N_dias    = .N,
  Media     = round(mean(DATO_DIARIO), 2),
  SD        = round(sd(DATO_DIARIO), 2),
  Mediana   = round(median(DATO_DIARIO), 2),
  Asimetria = round(skewness(DATO_DIARIO), 3),
  Curtosis  = round(kurtosis(DATO_DIARIO), 3),
  Min       = round(min(DATO_DIARIO), 2),
  Max       = round(max(DATO_DIARIO), 2)
), by = ESTACION]

setorder(desc_no2, ESTACION)

cat("=============================================================\n")
cat(" ESTADÍSTICAS DESCRIPTIVAS DE NO₂ — Madrid 2025 (µg/m³)\n")
cat(" Escala diaria · Datos crudos · Por estación\n")
cat("=============================================================\n\n")
print(desc_no2, digits = 4)

cat("\n-------------------------------------------------------------\n")
cat(sprintf(" Estaciones incluidas : %d\n", nrow(desc_no2)))
cat(sprintf(" Total días·estación  : %d\n", sum(desc_no2$N_dias)))
cat("-------------------------------------------------------------\n")

# ==============================================================================
# Tabla gt → PNG
# ==============================================================================

tbl_gt <- as.data.frame(desc_no2) |>
  gt() |>
  tab_header(
    title = html("<strong>Descriptive Statistics for NO<sub>2</sub></strong>"),
    subtitle = html(
      "Daily scale &middot; Raw data (&micro;g/m<sup>3</sup>) &middot; Madrid, 2025"
    )
  ) |>
  cols_label(
    ESTACION  = md("**Station**"),
    N_dias    = md("**Valid days**"),
    Media     = md("**Mean**"),
    SD        = md("**SD**"),
    Mediana   = md("**Median**"),
    Asimetria = md("**Skewness**"),
    Curtosis  = md("**Kurtosis**"),
    Min       = md("**Min**"),
    Max       = md("**Max**")
  ) |>
  tab_spanner(
    label   = md("**Central tendency**"),
    columns = c(Media, Mediana)
  ) |>
  tab_spanner(
    label   = md("**Dispersion**"),
    columns = c(SD, Min, Max)
  ) |>
  tab_spanner(
    label   = md("**Distribution shape**"),
    columns = c(Asimetria, Curtosis)
  ) |>
  data_color(
    columns = Media,
    method  = "numeric",
    palette = c("#d4e9f7", "#08519c")
  ) |>
  fmt_number(
    columns  = c(Media, SD, Mediana, Min, Max),
    decimals = 2
  ) |>
  fmt_number(
    columns  = c(Asimetria, Curtosis),
    decimals = 3
  ) |>
  tab_source_note(
    source_note = html(
      paste0(
        "<strong>Source.</strong> Madrid City Council Air Quality Monitoring ",
        "Network. Daily means are set to missing when at least 20% of hourly ",
        "measurements are missing."
      )
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#1a3a5c"),
      cell_text(color = "white", weight = "bold", size = px(15))
    ),
    locations = cells_title(groups = "title")
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#1a3a5c"),
      cell_text(color = "#d0e4f5", size = px(11))
    ),
    locations = cells_title(groups = "subtitle")
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#2c5f8a"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#2c5f8a"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_column_spanners()
  ) |>
  tab_options(
    table.font.names                = "Arial",
    table.font.size                 = 12,
    row.striping.include_table_body = TRUE,
    row.striping.background_color   = "#f7fafd",
    table.border.top.color          = "#1a3a5c",
    table.border.top.width          = px(3),
    table_body.hlines.color         = "#dde8f0",
    source_notes.font.size          = 9
  )

fwrite(
  desc_no2,
  file.path(carpeta_out, "descriptivo_NO2_diario_2025_por_estacion.csv")
)

gtsave(tbl_gt,
  filename = file.path(carpeta_out, "descriptivo_NO2_diario_2025_por_estacion.png"),
  zoom = 2, expand = 20
)

cat("\n\u2705 Archivos guardados en outputs/figures/EDA/NO2/Distribucion_global/descriptivo_NO2/:\n")
cat("   \u00b7 descriptivo_NO2_diario_2025_por_estacion.csv\n")
cat("   \u00b7 descriptivo_NO2_diario_2025_por_estacion.png\n")


# ==============================================================================
# BLOQUE 2: Tabla 3.1. Resumen anual de NO2, 2019-2025
# ==============================================================================
# La unidad de analisis es estacion-anio. Primero se obtiene la media anual de
# cada estacion a partir de sus dias validos; despues se resumen esas medias
# anuales y se cuentan las estaciones que superan cada umbral.

anios <- 2019:2025
rutas_diarias <- here(
  "data", "processed", "Contaminacion", "diario",
  sprintf("aire_madrid_%d_No2_trans_diarios1.rds", anios)
)

archivos_ausentes <- rutas_diarias[!file.exists(rutas_diarias)]
if (length(archivos_ausentes)) {
  stop(
    "No se encuentran los archivos diarios necesarios:\n  ",
    paste(archivos_ausentes, collapse = "\n  ")
  )
}

lista_no2_anual <- Map(function(ruta, anio) {
  datos <- as.data.table(readRDS(ruta))
  columnas_requeridas <- c("ESTACION", "FECHA", "DATO_DIARIO")
  columnas_ausentes <- setdiff(columnas_requeridas, names(datos))
  if (length(columnas_ausentes)) {
    stop(
      "El archivo de ", anio, " no contiene: ",
      paste(columnas_ausentes, collapse = ", ")
    )
  }
  datos[, Anio := as.integer(anio)]
  datos[, .(Anio, ESTACION, FECHA, DATO_DIARIO)]
}, rutas_diarias, anios)

dt_no2_historico <- rbindlist(lista_no2_anual, use.names = TRUE)

# Media anual por estacion: unidad comun de todas las columnas de la tabla.
medias_estacion_anio <- dt_no2_historico[, .(
  N_dias_validos = sum(!is.na(DATO_DIARIO)),
  Media_anual = if (all(is.na(DATO_DIARIO))) {
    NA_real_
  } else {
    mean(DATO_DIARIO, na.rm = TRUE)
  }
), by = .(Anio, ESTACION)]

resumen_anual <- medias_estacion_anio[!is.na(Media_anual), .(
  N = .N,
  Media = mean(Media_anual),
  Mediana = median(Media_anual),
  DE = sd(Media_anual),
  P95 = as.numeric(quantile(Media_anual, probs = 0.95, names = FALSE)),
  Estaciones_mayor_40 = sum(Media_anual > 40, na.rm = TRUE),
  Estaciones_mayor_20 = sum(Media_anual > 20, na.rm = TRUE)
), by = Anio]

setorder(resumen_anual, Anio)

if (!identical(resumen_anual$Anio, anios) || any(resumen_anual$N <= 0L)) {
  stop("El resumen anual no contiene exactamente los anos 2019-2025 esperados.")
}

# Archivo auxiliar con las estaciones que superan al menos el umbral de 20.
# Permite auditar los conteos presentados en la tabla principal.
detalle_umbrales <- copy(medias_estacion_anio[Media_anual > 20])
detalle_umbrales[, `:=`(
  Supera_40 = Media_anual > 40,
  Supera_20 = Media_anual > 20,
  Media_anual = round(Media_anual, 2)
)]
setorder(detalle_umbrales, Anio, -Media_anual, ESTACION)

resumen_anual_csv <- copy(resumen_anual)
columnas_decimales <- c("Media", "Mediana", "DE", "P95")
resumen_anual_csv[, (columnas_decimales) := lapply(.SD, round, 2),
  .SDcols = columnas_decimales
]

fwrite(
  resumen_anual_csv,
  file.path(carpeta_out, "tabla_3_1_resumen_anual_NO2_2019_2025.csv")
)
fwrite(
  detalle_umbrales,
  file.path(carpeta_out, "estaciones_sobre_umbrales_NO2_2019_2025.csv")
)

# Version de presentacion: todas las columnas se refieren a estaciones-anio.
tabla_anual_gt <- resumen_anual[, .(
  Anio, N, Media, Mediana, DE, P95,
  Estaciones_mayor_40, Estaciones_mayor_20
)] |>
  as.data.frame() |>
  gt() |>
  tab_header(
    title = html(
      "<strong>Table 3.1. Annual summary of NO<sub>2</sub> concentrations</strong>"
    ),
    subtitle = html(
      "Madrid, 2019&ndash;2025 &middot; Concentrations in &micro;g/m<sup>3</sup>"
    )
  ) |>
  cols_label(
    Anio = html("<strong>Year</strong>"),
    N = html("<strong>N</strong>"),
    Media = html("<strong>Mean</strong>"),
    Mediana = html("<strong>Median</strong>"),
    DE = html("<strong>SD</strong>"),
    P95 = html("<strong>P95</strong>"),
    Estaciones_mayor_40 = html("<strong>&gt; 40 &micro;g/m<sup>3</sup></strong>"),
    Estaciones_mayor_20 = html("<strong>&gt; 20 &micro;g/m<sup>3</sup></strong>")
  ) |>
  tab_spanner(
    label = html("<strong>Distribution of station annual means</strong>"),
    columns = c(Media, Mediana, DE, P95)
  ) |>
  tab_spanner(
    label = html(
      "<strong>Stations above the annual threshold</strong>"
    ),
    columns = c(Estaciones_mayor_40, Estaciones_mayor_20)
  ) |>
  fmt_integer(columns = Anio, use_seps = FALSE) |>
  fmt_integer(columns = N, use_seps = TRUE, sep_mark = ",") |>
  fmt_number(
    columns = c(Media, Mediana, DE, P95),
    decimals = 2, dec_mark = ".", sep_mark = ","
  ) |>
  fmt_integer(
    columns = c(Estaciones_mayor_40, Estaciones_mayor_20),
    use_seps = FALSE
  ) |>
  cols_align(align = "center", columns = everything()) |>
  cols_align(align = "right", columns = c(N, Media, Mediana, DE, P95)) |>
  tab_source_note(
    source_note = html(
      paste0(
        "<strong>Note.</strong> The unit of analysis is the station&ndash;year. ",
        "<em>N</em> is the number of stations with an estimable annual mean. ",
        "Mean, median, SD, and P95 summarize the distribution of station-level ",
        "annual means. Each annual mean is calculated from the available valid ",
        "daily observations; valid-day counts are reported in the detailed CSV file."
      )
    )
  ) |>
  tab_source_note(
    source_note = html(
      paste0(
        "<strong>Regulatory reference.</strong> 40 &micro;g/m<sup>3</sup> is the ",
        "annual limit value applicable to the study period; 20 ",
        "&micro;g/m<sup>3</sup> is the value adopted in Annex I to Directive (EU) ",
        "2024/2881, to be attained by 1 January 2030."
      )
    )
  ) |>
  tab_source_note(
    source_note = html(
      "<strong>Source.</strong> Madrid City Council Air Quality Monitoring Network."
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold", size = px(16)),
    locations = cells_title(groups = "title")
  ) |>
  tab_style(
    style = cell_text(size = px(11), color = "#333333"),
    locations = cells_title(groups = "subtitle")
  ) |>
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_borders(sides = "bottom", color = "#333333", weight = px(1.5))
    ),
    locations = list(cells_column_spanners(), cells_column_labels())
  ) |>
  tab_style(
    style = cell_borders(sides = "bottom", color = "#D9D9D9", weight = px(0.5)),
    locations = cells_body(rows = 1:(nrow(resumen_anual) - 1L))
  ) |>
  cols_width(
    Anio ~ px(72),
    N ~ px(95),
    c(Media, Mediana, DE, P95) ~ px(82),
    c(Estaciones_mayor_40, Estaciones_mayor_20) ~ px(120)
  ) |>
  tab_options(
    table.font.names = "Arial",
    table.font.size = 12,
    table.width = px(930),
    heading.align = "left",
    heading.border.bottom.color = "#111111",
    heading.border.bottom.width = px(2),
    table.border.top.color = "#111111",
    table.border.top.width = px(2),
    table.border.bottom.color = "#111111",
    table.border.bottom.width = px(2),
    data_row.padding = px(7),
    source_notes.font.size = 9,
    source_notes.padding = px(7)
  )

gtsave(
  tabla_anual_gt,
  filename = file.path(carpeta_out, "tabla_3_1_resumen_anual_NO2_2019_2025.png"),
  zoom = 2,
  expand = 20
)

cat("\nTABLA 3.1. RESUMEN ANUAL DE NO2, 2019-2025\n")
print(resumen_anual, digits = 4)
cat("\nArchivos anuales guardados:\n")
cat("   · tabla_3_1_resumen_anual_NO2_2019_2025.csv\n")
cat("   · tabla_3_1_resumen_anual_NO2_2019_2025.png\n")
cat("   · estaciones_sobre_umbrales_NO2_2019_2025.csv\n")
