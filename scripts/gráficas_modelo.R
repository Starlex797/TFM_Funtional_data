# ==============================================================================
# PASO 6A: GRÁFICOS DE INFERENCIA (EFECTOS FIJOS E HIPERPARÁMETROS)
# ==============================================================================

library(INLA)
library(ggplot2)
library(data.table)
library(here)

# 1. Cargar el modelo entrenado
modelo_final <- readRDS(here("data", "processed", "modelo_final_no2_madrid.rds"))

# 2. EXTRAER LOS EFECTOS FIJOS (Covariables deterministas)
# Extraemos la media y los cuantiles del 2.5% y 97.5% (Intervalo de Credibilidad del 95%)
efectos_fijos <- as.data.frame(modelo_final$summary.fixed[, c("mean", "0.025quant", "0.975quant")])
efectos_fijos$Variable <- rownames(efectos_fijos)

# Limpiamos el nombre del intercepto y lo quitamos del gráfico para no distorsionar la escala
efectos_fijos <- efectos_fijos[efectos_fijos$Variable != "intercept", ]
setDT(efectos_fijos)

# 3. FOREST PLOT DE COVARIABLES
# Este gráfico te dirá qué aumenta o disminuye el NO2 en Madrid
pdf(here("output", "figures", "06A_forest_plot_covariables.pdf"), width = 7, height = 5)

ggplot(efectos_fijos, aes(x = reorder(Variable, mean), y = mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`), size = 0.8, color = "#2c3e50") +
  coord_flip() +
  labs(
    title = "Efectos del Tráfico y Clima sobre el Log(NO2)",
    subtitle = "Intervalos de Credibilidad del 95%",
    x = "Covariables",
    y = "Impacto Estimado (Beta)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = "lightgrey"),
    plot.title = element_text(face = "bold")
  )

dev.off()

# 4. EXTRAER Y DIBUJAR LOS HIPERPARÁMETROS DEL SPDE Y AR1
# Extraemos los marginales de la correlación temporal (rho), la varianza y el alcance (kappa)
marginal_rho   <- inla.smarginal(modelo_final$marginals.hyperpar$`GroupRho for campo_espacial`)
marginal_range <- inla.smarginal(modelo_final$marginals.hyperpar$`Range for campo_espacial`)

df_rho   <- data.frame(x = marginal_rho$x, y = marginal_rho$y, parametro = "Persistencia Temporal (AR1)")
df_range <- data.frame(x = marginal_range$x, y = marginal_range$y, parametro = "Alcance Espacial (km)")

pdf(here("output", "figures", "06B_densidad_hiperparametros.pdf"), width = 10, height = 4)

# Gráfico del AR1 temporal
p1 <- ggplot(df_rho, aes(x = x, y = y)) +
  geom_area(fill = "#3498db", alpha = 0.6) +
  geom_line(color = "#2980b9", size = 1) +
  labs(title = "Memoria Temporal del NO2 (Rho)", x = "Valor de Correlación AR(1)", y = "Densidad Posterior") +
  theme_minimal()

# Gráfico del Rango Espacial
p2 <- ggplot(df_range, aes(x = x, y = y)) +
  geom_area(fill = "#e67e22", alpha = 0.6) +
  geom_line(color = "#d35400", size = 1) +
  labs(title = "Radio de Contagio Espacial", x = "Distancia en Kilómetros", y = "Densidad Posterior") +
  theme_minimal()

library(gridExtra) # Usamos gridExtra como hizo Wuong
grid.arrange(p1, p2, ncol = 2)

dev.off()

