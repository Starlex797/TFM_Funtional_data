# ==============================================================================
# FIGURA DE APOYO: la correlacion NO2-clima NO es estable
# Panel A: el coeficiente de Pearson cambia con la ESCALA (horaria/diaria/mensual)
# Panel B: y cambia con la ESTACION DEL ANO (invierno vs verano), a escala diaria
# Estacion representativa: Plaza Eliptica (2025). Correlaciones CRUDAS (para
# ilustrar justamente su inestabilidad, que motiva el enfoque por anomalias).
# Salida: outputs/analysis/limitaciones_correlacion/inestabilidad_correlacion.png
# ==============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(here)

EST <- "Plaza Elíptica"
COVS <- list(
  "Temperatura"   = "Temperatura_raw",
  "Humedad"       = "Humedad_Relativa_raw",
  "Precipitacion" = "Precipitaciones_raw",
  "Presion"       = "Presion Barométrica_raw",
  "Radiacion"     = "Radiación Solar_raw",
  "Viento"        = "Velocidad Viento_raw"
)
niv_cov <- names(COVS)

corp <- function(dt, no2, cc) {
  d <- dt[is.finite(get(no2)) & is.finite(get(cc))]
  if (nrow(d) < 5) return(NA_real_)
  cor(d[[no2]], d[[cc]])
}

h <- as.data.table(readRDS(here("data","processed","Maestro","horario",
        "dataset_maestro_inla_2025_HORARIO.rds")))[ESTACION == EST]
d <- as.data.table(readRDS(here("data","processed","Maestro","diario",
        "dataset_maestro_inla_2025_DIARIO.rds")))[ESTACION == EST]
m <- as.data.table(readRDS(here("data","processed","Maestro","mensual",
        "dataset_maestro_inla_2025_MENSUAL.rds")))[ESTACION == EST]

# ---- Panel A: por escala (todo el ano) ----
A <- rbindlist(lapply(names(COVS), function(lab) {
  cc <- COVS[[lab]]
  data.table(Covariable = lab,
    Escala = c("Horaria","Diaria","Mensual"),
    r = c(corp(h, "DATO", cc), corp(d, "DATO_DIARIO", cc), corp(m, "DATO_MENSUAL", cc)))
}))
A[, Covariable := factor(Covariable, levels = niv_cov)]
A[, Escala := factor(Escala, levels = c("Horaria","Diaria","Mensual"))]

# ---- Panel B: por estacion del ano (escala diaria) ----
d[, mes := as.integer(format(FECHA, "%m"))]
d[, Temporada := fifelse(mes %in% c(12,1,2), "Invierno",
                 fifelse(mes %in% c(6,7,8),  "Verano", NA_character_))]
B <- rbindlist(lapply(names(COVS), function(lab) {
  cc <- COVS[[lab]]
  data.table(Covariable = lab,
    Temporada = c("Invierno","Verano"),
    r = c(corp(d[Temporada=="Invierno"], "DATO_DIARIO", cc),
          corp(d[Temporada=="Verano"],   "DATO_DIARIO", cc)))
}))
B[, Covariable := factor(Covariable, levels = niv_cov)]
B[, Temporada := factor(Temporada, levels = c("Invierno","Verano"))]

tema <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "gray40", size = 9),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top", panel.grid.minor = element_blank())

pA <- ggplot(A, aes(Covariable, r, fill = Escala)) +
  geom_col(position = position_dodge(0.75), width = 0.7, color = "white") +
  geom_hline(yintercept = 0, color = "grey30") +
  scale_fill_manual(values = c(Horaria="#7fb3d5", Diaria="#f39c12", Mensual="#c0392b"),
                    name = NULL) +
  coord_cartesian(ylim = c(-0.8, 0.8)) +
  labs(title = "A. La correlacion cambia con la ESCALA",
       subtitle = "Pearson NO2-covariable (todo 2025) - Plaza Eliptica",
       x = NULL, y = "Correlacion de Pearson") + tema

pB <- ggplot(B, aes(Covariable, r, fill = Temporada)) +
  geom_col(position = position_dodge(0.7), width = 0.65, color = "white") +
  geom_hline(yintercept = 0, color = "grey30") +
  scale_fill_manual(values = c(Invierno="#2980b9", Verano="#e67e22"), name = NULL) +
  coord_cartesian(ylim = c(-0.8, 0.8)) +
  labs(title = "B. Y cambia con la ESTACION DEL ANO",
       subtitle = "Pearson NO2-covariable a escala diaria - Plaza Eliptica",
       x = NULL, y = "Correlacion de Pearson") + tema

fig <- pA + pB +
  plot_annotation(
    title = "La correlacion NO2-clima no es una constante: depende de la escala y de la epoca",
    subtitle = "Motiva medir sobre anomalias, por estacion y por escala (no un unico coeficiente crudo)",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(color = "gray40")))

dir_out <- here("outputs", "analysis", "limitaciones_correlacion")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(dir_out, "inestabilidad_correlacion.png"),
       fig, width = 14, height = 6.5, dpi = 200, bg = "white")
cat("Guardado: outputs/analysis/limitaciones_correlacion/inestabilidad_correlacion.png\n")
print(dcast(A, Covariable ~ Escala, value.var = "r"))
print(dcast(B, Covariable ~ Temporada, value.var = "r"))
