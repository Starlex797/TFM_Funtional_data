# ==============================================================================
# PASO 4: EL SÁNDWICH DE DATOS (CONSTRUCCIÓN DEL INLA STACK)
# ==============================================================================

library(INLA)
library(data.table)
library(here)

# 1. Cargar los objetos de los pasos anteriores 
dt_maestro   <- readRDS(here("data", "processed", "dataset_maestro_inla_2025.rds"))
malla_madrid <- readRDS(here("data", "processed", "malla_spde_madrid.rds"))
setDT(dt_maestro)

# Recordatorio crítico de ordenación para mantener la sincronía con la Matriz A
setorder(dt_maestro, ID_TIEMPO, ESTACION)

# 2. Empaquetar todo en el INLA Stack
# Vinculamos la variable respuesta, las matrices de proyección y los efectos fijos/aleatorios
stk_madrid <- inla.stack(
  tag = "estimacion_no2",                    # Etiqueta para identificar este bloque
  
  # Bloque A: Variable respuesta (lo que queremos predecir)
  data = list(y = dt_maestro$LOG_NO2_DIARIO), 
  
  # Bloque B: Las matrices de proyección correspondientes a cada efecto
  # A_espacial se aplica al campo aleatorio; el '1' indica una relación directa (1:1) para las covariables
  A = list(A_espacial, 1),                    
  
  # Bloque C: Los Efectos (Variables explicativas y campos latentes)
  effects = list(
    # Efecto 1: El campo espacio-temporal indexado y el intercepto base del modelo
    c(indice_s, list(intercept = 1)), 
    
    # Efecto 2: Covariables reales (Fijos y Aleatorios no espaciales)
    list(
      # --- Variables de Tráfico (Área de distrito) ---
      trafico_intensidad = dt_maestro$intensidad, 
      trafico_carga      = dt_maestro$carga,
      
      # --- Variables Meteorológicas (Media de Madrid) ---
      temperatura        = dt_maestro$Temperatura, 
      viento             = dt_maestro$`Velocidad Viento`,
      precipitacion      = dt_maestro$Precipitaciones,
      humedad            = dt_maestro$Humedad_Relativa,
      presion_barometrica  = dt_maestro$Presion_Barometrica,
      direccion_viento      = dt_maestro$Dir.Viento,
      
      # --- Efecto iid por Distrito ---
      id_distrito        = dt_maestro$ID_DISTRITO
    )
  ),
  compress = FALSE
)

# 3. Guardar el Stack preparado en tu disco
saveRDS(stk_madrid, here("data", "processed", "inla_stack_madrid_2025.rds"))

# ==============================================================================
# STACK PARA EL MODELO SOLO ESPACIAL
# Usa A_espacial_s e indice_s_solo definidos en Spde.R (sin grupo temporal).
# Las covariables y el efecto IID de distrito son idénticos al modelo anterior.

stk_madrid_s <- inla.stack(
  tag  = "estimacion_no2_espacial",
  data = list(y = dt_maestro$LOG_NO2_DIARIO),
  A    = list(A_espacial_s, 1),
  effects = list(
    c(indice_s_solo, list(intercept = 1)),
    list(
      trafico_intensidad = dt_maestro$intensidad,
      trafico_carga      = dt_maestro$carga,
      temperatura        = dt_maestro$Temperatura,
      viento             = dt_maestro$`Velocidad Viento`,
      precipitacion      = dt_maestro$Precipitaciones,
      humedad            = dt_maestro$Humedad_Relativa,
      presion_barometrica  = dt_maestro$Presion_Barometrica,
      direccion_viento      = dt_maestro$Dir.Viento,
      id_distrito        = dt_maestro$ID_DISTRITO
    )
  ),
  compress = FALSE
)

saveRDS(stk_madrid_s, here("data", "processed", "inla_stack_espacial_2025.rds"))

# ------------------------------------------------------------------------------
# CONTROL DE CALIDAD Y DIAGNÓSTICO
# ------------------------------------------------------------------------------
cat("✅ ¡Paso 4 completado con éxito!\n")
cat("- Stack espacio-temporal | obs:", inla.stack.ndata(stk_madrid),
    "| debe ser:", nrow(dt_maestro), "\n")
cat("- Stack solo espacial    | obs:", inla.stack.ndata(stk_madrid_s),
    "| debe ser:", nrow(dt_maestro), "\n")
