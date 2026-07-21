# ==============================================================================
# ANÁLISIS DESCRIPTIVO — NO₂ (Datos crudos)
# Madrid 2025 · Escala diaria · Por estación de monitoreo
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
# Guarda:un CSV con las estadísticas;
# una imagen PNG de la tabla.
# ==============================================================================
library(data.table)
library(here)
library(moments) # skewness(), kurtosis()
library(gt) # tabla formateada → PNG


# ==============================================================================
# BLOQUE 1: Estadísticas descriptivas básicas de NO2 por estación
# ==============================================================================

# Carga de datos diarios NO2 2025 (valores crudos, sin transformación log)
dt_no2 <- readRDS(here(
  "data", "processed", "contaminacion", "diario",
  "aire_madrid_2025_No2_trans_diarios.rds"
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
    title    = md("**Estadísticas Descriptivas de NO\u2082**"),
    subtitle = md("Escala diaria \u00b7 Datos crudos (\u00b5g/m\u00b3) \u00b7 Madrid 2025")
  ) |>
  cols_label(
    ESTACION  = md("**Estaci\u00f3n**"),
    N_dias    = md("**N d\u00edas**"),
    Media     = md("**Media**"),
    SD        = md("**DE**"),
    Mediana   = md("**Mediana**"),
    Asimetria = md("**Asimetr\u00eda**"),
    Curtosis  = md("**Curtosis**"),
    Min       = md("**M\u00edn**"),
    Max       = md("**M\u00e1x**")
  ) |>
  tab_spanner(
    label   = md("**Tendencia central**"),
    columns = c(Media, Mediana)
  ) |>
  tab_spanner(
    label   = md("**Dispersi\u00f3n**"),
    columns = c(SD, Min, Max)
  ) |>
  tab_spanner(
    label   = md("**Forma**"),
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
    source_note = md("Fuente: Red de Monitoreo de Calidad del Aire de Madrid \u00b7 Media diaria de registros horarios v\u00e1lidos (umbral NA \u2264 20\u202f%)")
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

# Guardar outputs
carpeta_out <- here("outputs", "descriptivo_NO2")
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

fwrite(
  desc_no2,
  file.path(carpeta_out, "descriptivo_NO2_diario_2025_por_estacion.csv")
)

gtsave(tbl_gt,
  filename = file.path(carpeta_out, "descriptivo_NO2_diario_2025_por_estacion.png"),
  zoom = 2, expand = 20
)

cat("\n\u2705 Archivos guardados en outputs/descriptivo_NO2/:\n")
cat("   \u00b7 descriptivo_NO2_diario_2025_por_estacion.csv\n")
cat("   \u00b7 descriptivo_NO2_diario_2025_por_estacion.png\n")
