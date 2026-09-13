# ==============================================================================
# RENDER ACADEMIC TABLES FROM EXISTING VARIABLE-SELECTION RESULTS
# ==============================================================================
# This script does not fit INLA models. It reads the CSV files already produced
# by Seleccion_variables.R and creates English CSV and PNG tables for the TFM.
# ==============================================================================

library(data.table)
library(here)

source(here("R", "utilities", "academic_quality_tables.R"))

DIR_RESULTS <- here(
    "outputs", "modelo", "Modelo_1", "seleccion_variables", "2019_2021"
)
DIR_TABLES <- file.path(DIR_RESULTS, "academic_tables")
dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)

MODEL_DESCRIPTIONS <- c(
    M0 = "Full covariate model",
    M1 = "Without temperature and solar radiation",
    M2 = "M1 with nonlinear wind speed (RW2)",
    M3 = "Without relative humidity and solar radiation",
    M4 = "M3 with nonlinear temperature (RW2)",
    M5 = "Nonlinear temperature and wind speed (RW2)",
    M6 = "Without meteorological covariates",
    M7 = "Meteorological covariates only"
)

VARIABLE_LABELS <- c(
    intensidad = "Traffic intensity",
    Temperatura = "Temperature",
    Velocidad_Viento = "Wind speed",
    Velocidad_Viento_sqrt = "Square-root wind speed",
    Presion_Barometrica = "Barometric pressure",
    Humedad_Relativa = "Relative humidity",
    Precipitaciones = "Precipitation",
    Llueve = "Rain indicator",
    Radiacion_Solar = "Solar radiation",
    Tipo_Urbana_fondo = "Station type: urban background",
    Tipo_Urbana_trafico = "Station type: urban traffic",
    Intercept = "Intercept"
)

translate_variable <- function(x) {
    translated <- unname(VARIABLE_LABELS[x])
    translated[is.na(translated)] <- gsub("_", " ", x[is.na(translated)])
    translated
}

translate_station_type <- function(x) {
    translated <- as.character(x)
    translated[grepl("suburb", translated, ignore.case = TRUE)] <- "Suburban"
    translated[grepl("fondo", translated, ignore.case = TRUE)] <-
        "Urban background"
    translated[grepl("tr", translated, ignore.case = TRUE)] <- "Urban traffic"
    translated
}

format_number <- function(x, digits = 2L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "--")
}

required_file <- function(name) {
    path <- file.path(DIR_RESULTS, name)
    if (!file.exists(path)) stop("Required result file not found: ", path)
    path
}

# Overall model comparison -----------------------------------------------------
selection <- fread(required_file("tabla_comparacion_modelos_2019_2021.csv"))
selection[, Descripcion := unname(MODEL_DESCRIPTIONS[Modelo])]
model_comparison <- selection[, .(
    Model = Modelo,
    Description = Descripcion,
    WAIC = format_number(WAIC, 1L),
    `Delta WAIC` = format_number(Delta_WAIC_vs_M0, 1L),
    DIC = format_number(DIC, 1L),
    `Delta DIC` = format_number(Delta_DIC_vs_M0, 1L),
    `Fitted RMSE` = format_number(RMSE_INLA_orientativo, 3L),
    `Fitted COV95 (%)` = format_number(COV95_INLA_orientativo, 1L)
)]
fwrite(model_comparison, file.path(DIR_TABLES, "model_comparison_2019_2021.csv"))
path_model <- booktabs_png(
    model_comparison,
    file.path(DIR_TABLES, "table_model_comparison_2019_2021.png"),
    title = "Comparison of candidate INLA-SPDE models",
    subtitle = "Daily NO2 models for Madrid, 2019-2021",
    note = paste0(
        "Delta values are relative to M0; negative values favour the candidate. ",
        "Lower WAIC, DIC and RMSE values are preferred. Fitted RMSE and COV95 ",
        "are in-sample diagnostics, not spatial validation metrics."
    ),
    widths = c(0.55, 3.15, 0.78, 0.88, 0.78, 0.88, 0.90, 1.02),
    align = c("center", "left", rep("right", 6)),
    font_size = 7.8,
    row_height = 0.32
)

# Fixed effects: one table per model ------------------------------------------
coefficient_paths <- vapply(paste0("M", 0:7), function(id) {
    input <- fread(required_file(sprintf(
        "tabla_coeficientes_%s_2019_2021.csv",
        id
    )))
    output <- input[, .(
        Variable = translate_variable(Variable),
        `Posterior mean` = format_number(Coeficiente, 4L),
        `95% credible interval` = IC95,
        `Significant (95%)` = fifelse(Significativa_95 == "Si", "Yes", "No")
    )]
    fwrite(
        output,
        file.path(DIR_TABLES, sprintf("fixed_effects_%s_2019_2021.csv", id))
    )
    booktabs_png(
        output,
        file.path(
            DIR_TABLES,
            sprintf("table_fixed_effects_%s_2019_2021.png", id)
        ),
        title = sprintf("Fixed effects for model %s", id),
        subtitle = MODEL_DESCRIPTIONS[[id]],
        note = paste0(
            "A fixed effect is labelled significant when its 95% posterior credible ",
            "interval excludes zero. Station-type effects are differences from the ",
            "suburban reference category. RW2 effects have no single coefficient."
        ),
        widths = c(2.55, 1.15, 1.85, 1.25),
        align = c("left", "right", "center", "center"),
        font_size = 8.2,
        row_height = 0.31
    )
}, character(1))

# Spatial hold-out -------------------------------------------------------------
holdout <- fread(required_file("tabla_holdout_espacial_M5_2019_2021.csv"))
holdout_station <- fread(required_file(
    "tabla_holdout_por_estacion_M5_2019_2021.csv"
))
holdout[, Descripcion := unname(MODEL_DESCRIPTIONS[Modelo])]

holdout_global_english <- holdout[, .(
    Model = Modelo,
    Description = Descripcion,
    RMSE = format_number(RMSE_HOLDOUT, 3L),
    MAE = format_number(MAE_HOLDOUT, 3L),
    Bias = format_number(Sesgo_HOLDOUT, 3L),
    `COV95 (%)` = format_number(COV95_HOLDOUT, 1L),
    `Mean 95% width` = format_number(Anchura_media_IC95_HOLDOUT, 3L),
    Observations = as.character(Observaciones_HOLDOUT)
)]
holdout_station_english <- holdout_station[, .(
    Model = Modelo,
    Station = ESTACION,
    `Station type` = translate_station_type(NOM_TIPO),
    Observations = as.character(Observaciones),
    RMSE = format_number(RMSE, 3L),
    MAE = format_number(MAE, 3L),
    Bias = format_number(Sesgo, 3L),
    `COV95 (%)` = format_number(COV95, 1L),
    `Mean 95% width` = format_number(Anchura_media_IC95, 3L)
)]

fwrite(
    holdout_global_english,
    file.path(DIR_TABLES, "spatial_holdout_M5_2019_2021.csv")
)
fwrite(
    holdout_station_english,
    file.path(DIR_TABLES, "spatial_holdout_by_station_M5_2019_2021.csv")
)
path_holdout <- booktabs_png(
    holdout_global_english,
    file.path(DIR_TABLES, "table_spatial_holdout_M5_2019_2021.png"),
    title = "Spatial hold-out performance of model M5",
    subtitle = "Three monitoring stations excluded jointly",
    note = paste0(
        "Metrics use only excluded-station observations. Predictive COV95 includes ",
        "observation variance and must be interpreted with RMSE and interval width."
    ),
    widths = c(0.55, 2.75, 0.72, 0.70, 0.70, 0.78, 1.00, 0.90),
    align = c("center", "left", rep("right", 6)),
    font_size = 8.0,
    row_height = 0.32
)
path_holdout_station <- booktabs_png(
    holdout_station_english,
    file.path(
        DIR_TABLES,
        "table_spatial_holdout_by_station_M5_2019_2021.png"
    ),
    title = "Spatial hold-out performance by station",
    subtitle = "Model M5; daily NO2 observations, 2019-2021",
    note = paste0(
        "The three stations are excluded simultaneously during fitting. Coverage ",
        "must be assessed together with prediction error and interval width."
    ),
    widths = c(0.55, 1.45, 1.30, 0.82, 0.72, 0.68, 0.68, 0.78, 0.96),
    align = c("center", "left", "left", rep("right", 6)),
    font_size = 7.7,
    row_height = 0.34
)

all_png <- c(path_model, coefficient_paths, path_holdout, path_holdout_station)
if (!all(file.exists(all_png))) stop("Some academic PNG tables were not generated.")
cat("Academic PNG tables generated: ", length(all_png), "\n", sep = "")
cat("Directory: ", DIR_TABLES, "\n", sep = "")
