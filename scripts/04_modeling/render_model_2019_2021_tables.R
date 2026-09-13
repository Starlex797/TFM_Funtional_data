# ==============================================================================
# RENDER ACADEMIC TABLES FOR THE TWO MISALIGNMENT METHODS
# ==============================================================================
# Reuses saved INLA models and CSV results. It does not refit M1 or M2.
# ==============================================================================

library(data.table)
library(here)

source(here("R", "utilities", "academic_quality_tables.R"))

COVARIABLE_FIELD <- Sys.getenv(
    "INLA_COVARIABLE_FIELD", "Velocidad_Viento"
)
FIELD_OPTIONS <- c(
    Velocidad_Viento = "wind-speed",
    Presion_Barometrica = "pressure",
    Radiacion_Solar = "solar-radiation",
    Temperatura = "temperature"
)
if (!COVARIABLE_FIELD %in% names(FIELD_OPTIONS)) {
    stop(
        "INLA_COVARIABLE_FIELD must be one of: ",
        paste(names(FIELD_OPTIONS), collapse = ", ")
    )
}
DIR_RESULTS <- here(
    "outputs", "modelo", "Modelo_2019_2021",
    "desalineamiento_campo_espacial", COVARIABLE_FIELD, "2019_2021"
)
DIR_TABLES <- file.path(DIR_RESULTS, "academic_tables")
dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)

METHOD_LABELS <- c(
    M1 = "Assigned meteorological covariate",
    M2 = paste(
        "Joint spatial latent",
        FIELD_OPTIONS[[COVARIABLE_FIELD]],
        "field"
    )
)
VARIABLE_LABELS <- c(
    Velocidad_Viento = "Wind speed",
    Llueve = "Rain indicator",
    Presion_Barometrica = "Barometric pressure",
    Radiacion_Solar = "Solar radiation",
    Temperatura = "Temperature",
    intensidad = "Traffic intensity",
    Tipo_Urbana_fondo = "Station type: urban background",
    Tipo_Urbana_trafico = "Station type: urban traffic",
    Intercept = "Intercept"
)

format_number <- function(x, digits = 2L) {
    ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "--")
}
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
required_file <- function(path) {
    if (!file.exists(path)) stop("Required file not found: ", path)
    path
}

comparison_path <- required_file(file.path(
    DIR_RESULTS,
    "tabla_comparacion_modelos_2019_2021.csv"
))
comparison <- fread(comparison_path)
# Modelo_2019_2021.R ya calcula CPO, WAIC y DIC exclusivamente sobre las filas
# de NO2. No se recalculan aqui porque M2 contiene otras dos verosimilitudes.
required_columns <- c(
    "CPO_mean_neg_log", "CPO_failures", "Delta_CPO_vs_M1"
)
missing_columns <- setdiff(required_columns, names(comparison))
if (length(missing_columns)) {
    stop("Missing NO2 comparison columns: ", paste(missing_columns, collapse = ", "))
}
setorder(comparison, Modelo)

# Table 1: fit comparison, including CPO.
comparison_english <- comparison[, .(
    Model = Modelo,
    Method = unname(METHOD_LABELS[Modelo]),
    WAIC = format_number(WAIC, 1L),
    `Delta WAIC` = format_number(Delta_WAIC_vs_M1, 1L),
    DIC = format_number(DIC, 1L),
    `Delta DIC` = format_number(Delta_DIC_vs_M1, 1L),
    `Mean -log(CPO)` = format_number(CPO_mean_neg_log, 4L),
    `Delta CPO` = format_number(Delta_CPO_vs_M1, 4L),
    `WAIC p_eff` = format_number(WAIC_p_eff, 1L),
    `Fitted RMSE` = format_number(RMSE_ajuste, 3L),
    `Delta fitted RMSE` = format_number(Delta_RMSE_vs_M1, 3L)
)]
fwrite(comparison_english, file.path(DIR_TABLES, "misalignment_method_comparison.csv"))
path_comparison <- booktabs_png(
    comparison_english,
    file.path(DIR_TABLES, "table_misalignment_method_comparison.png"),
    title = "Comparison of spatial misalignment methods",
    subtitle = "Daily NO2 models for Madrid, 2019-2021",
    note = paste0(
        "Delta values are relative to M1; negative values favour M2. Lower WAIC, ",
        "DIC, mean -log(CPO) and RMSE values are preferred. Mean -log(CPO) is an ",
        "average leave-one-out logarithmic score. CPO failures: ",
        paste(sprintf("%s=%d", comparison$Modelo, comparison$CPO_failures), collapse = "; "),
        ". Fitted RMSE is an in-sample diagnostic."
    ),
    widths = c(.50, 3.10, .75, .82, .75, .82, .98, .82, .78, .88, 1.05),
    align = c("center", "left", rep("right", 9)),
    font_size = 7.3,
    row_height = .34
)

# Table 2: common spatial hold-out.
holdout <- fread(required_file(file.path(
    DIR_RESULTS,
    "tabla_holdout_espacial_2019_2021.csv"
)))
holdout_english <- comparison[, .(
    Model = Modelo,
    Method = unname(METHOD_LABELS[Modelo]),
    RMSE = format_number(RMSE_HOLDOUT, 3L),
    `Delta RMSE` = format_number(Delta_RMSE_HOLDOUT_vs_M1, 3L),
    MAE = format_number(MAE_HOLDOUT, 3L),
    Bias = format_number(Sesgo_HOLDOUT, 3L),
    `COV95 (%)` = format_number(COV95_HOLDOUT, 1L),
    `Mean 95% width` = format_number(Anchura_media_IC95_HOLDOUT, 3L),
    `Delta width` = format_number(Delta_Anchura_HOLDOUT_vs_M1, 3L),
    Observations = as.character(Observaciones_HOLDOUT)
)]
fwrite(holdout_english, file.path(DIR_TABLES, "spatial_holdout_comparison.csv"))
path_holdout <- booktabs_png(
    holdout_english,
    file.path(DIR_TABLES, "table_spatial_holdout_comparison.png"),
    title = "Spatial hold-out comparison of misalignment methods",
    subtitle = paste(
        "Stations excluded jointly: Plaza Castilla, Casa de Campo,",
        "Ensanche Vallecas and Villaverde Alto"
    ),
    note = paste0(
        "Metrics use only excluded-station observations. Negative deltas favour M2. ",
        "Predictive COV95 includes observation variance and must be interpreted with ",
        "RMSE and mean interval width."
    ),
    widths = c(.50, 3.10, .72, .82, .70, .70, .78, .98, .82, .88),
    align = c("center", "left", rep("right", 8)),
    font_size = 7.4,
    row_height = .34
)

# Table 3: hold-out by model and station.
holdout_station <- fread(required_file(file.path(
    DIR_RESULTS,
    "tabla_holdout_por_estacion_2019_2021.csv"
)))
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
fwrite(holdout_station_english, file.path(DIR_TABLES, "spatial_holdout_by_station.csv"))
path_holdout_station <- booktabs_png(
    holdout_station_english,
    file.path(DIR_TABLES, "table_spatial_holdout_by_station.png"),
    title = "Spatial hold-out performance by station",
    subtitle = "Two spatial misalignment methods; Madrid, 2019-2021",
    note = paste0(
        "The four stations are excluded simultaneously during fitting. Coverage ",
        "must be assessed together with prediction error and interval width."
    ),
    widths = c(.50, 1.48, 1.32, .82, .72, .68, .68, .78, .96),
    align = c("center", "left", "left", rep("right", 6)),
    font_size = 7.6,
    row_height = .32
)

# Tables 4-5: fixed coefficients for M1 and M2.
coefficients <- fread(required_file(file.path(
    DIR_RESULTS,
    "tabla_coeficientes_modelos_2019_2021.csv"
)))
coefficient_paths <- vapply(c("M1", "M2"), function(id) {
    table_model <- coefficients[Modelo == id, .(
        Variable = translate_variable(Variable),
        `Posterior mean` = format_number(Coeficiente, 4L),
        `95% credible interval` = IC95,
        `Significant (95%)` = fifelse(Significativa_95 == "Si", "Yes", "No")
    )]
    fwrite(table_model, file.path(DIR_TABLES, sprintf("fixed_effects_%s.csv", id)))
    booktabs_png(
        table_model,
        file.path(DIR_TABLES, sprintf("table_fixed_effects_%s.png", id)),
        title = sprintf("Fixed effects for model %s", id),
        subtitle = METHOD_LABELS[[id]],
        note = paste0(
            "A fixed effect is labelled significant when its 95% posterior credible ",
            "interval excludes zero. Station-type effects are differences from the ",
            "suburban reference category. In M2, the meteorological coefficient is ",
            "the estimated copy scaling parameter. RW2 curves are reported separately."
        ),
        widths = c(2.55, 1.15, 1.85, 1.25),
        align = c("left", "right", "center", "center"),
        font_size = 8.1,
        row_height = .31
    )
}, character(1))

all_png <- c(path_comparison, path_holdout, path_holdout_station, coefficient_paths)
if (!all(file.exists(all_png))) stop("Some academic PNG tables were not generated.")
cat("Academic PNG tables generated: ", length(all_png), "\n", sep = "")
cat("Directory: ", DIR_TABLES, "\n", sep = "")
