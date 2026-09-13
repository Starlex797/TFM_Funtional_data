# Temporary verification: real preparation, isolated table-layout fixtures.
expressions <- parse("scripts/04_modeling/Modelo/seleccion_variables/Seleccion_variables.R")
assigned <- function(e) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.symbol(e[[2]])) {
        as.character(e[[2]])
    } else ""
}
assignment_names <- vapply(expressions, assigned, character(1))
for (e in expressions[seq_len(match("ajustes", assignment_names) - 1L)]) eval(e)
stopifnot(identical(names(formulas_modelos), paste0("M", 0:13)))
variables <- lapply(formulas_modelos, all.vars)
available <- names(INLA::inla.stack.data(stack_modelo))
for (id in names(variables)) {
    missing <- setdiff(variables[[id]], c(available, "spde"))
    if (length(missing)) stop(id, ": missing variables: ", paste(missing, collapse = ", "))
}
removed <- list(M8 = "Temperatura_rw2", M9 = "Velocidad_Viento_rw2",
                M10 = "Presion_Barometrica", M11 = "intensidad", M12 = "Llueve",
                M13 = c("Tipo_Urbana_fondo", "Tipo_Urbana_trafico"))
for (id in names(removed)) {
    stopifnot(setequal(setdiff(variables$M5, variables[[id]]), removed[[id]]),
              length(setdiff(variables[[id]], variables$M5)) == 0L)
}
cat("PREPARATION_OK: 14 formulas, all variables in stack; M8-M13 removal terms verified.\n")
cat("CATALOG_PNG:", path_model_catalog, "\n")

# Fixtures are NOT model results and are written only to a temporary directory.
DIR_OUT <- tempfile("selection_table_layout_")
DIR_TABLAS_ACADEMICAS <- file.path(DIR_OUT, "academic_tables")
dir.create(DIR_TABLAS_ACADEMICAS, recursive = TRUE)
CALCULAR_HOLDOUT <- FALSE
ids <- names(formulas_modelos)
tabla_seleccion <- data.table(
    Modelo = ids, Descripcion = unname(descripcion_modelos[ids]),
    WAIC = seq_along(ids), DIC = seq_along(ids) + 1,
    RMSE_INLA_orientativo = seq_along(ids) / 10,
    Delta_WAIC_vs_M0 = seq_along(ids) - 1,
    Delta_DIC_vs_M0 = seq_along(ids) - 1,
    p_eff_WAIC = 1
)
first <- match("terminos_eliminados", assignment_names)
last <- match("tabla_holdout", assignment_names) - 1L
for (e in expressions[first:last]) eval(e)
stopifnot(identical(tabla_aportacion_vs_M5$Modelo, paste0("M", 8:13)),
          all(tabla_aportacion_vs_M5$Delta_WAIC_vs_M5 == 3:8))
tablas_coeficientes <- setNames(lapply(ids, function(id) {
    fixed <- setdiff(variables[[id]], c("y_response", "Intercept", "spde",
                                      "campo_espacial", "Temperatura_rw2", "Velocidad_Viento_rw2"))
    data.table(Variable = fixed, Coeficiente = 0, IC95 = "[-1.0000, 1.0000]",
               Significativa_95 = "No")
}), ids)
first <- match("VARIABLE_LABELS", assignment_names)
last <- match("expected_png", assignment_names)
for (e in expressions[first:last]) eval(e)
stopifnot(nrow(model_comparison_english) == 14L, ncol(model_comparison_english) == 8L,
          length(coefficient_table_paths) == 14L,
          nrow(variable_contribution_english) == 6L,
          ncol(variable_contribution_english) == 5L,
          all(file.exists(expected_png)),
          all(file.info(expected_png)$size > 0))
cat("RENDER_OK: comparison 14x8, removal table 6x5, 14 coefficient PNGs (4 columns).\n")
cat("LAYOUT_FIXTURES:", DIR_TABLAS_ACADEMICAS, "\n")
