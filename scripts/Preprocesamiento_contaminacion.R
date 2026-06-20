#------------------------------------------------------------------------
# Goal: load the raw data of air pollution, clean it, and save it by year.
#------------------------------------------------------------------------
library(here)
library(data.table)
library(lubridate)

# 1). Charge the dictionaries and cleaning functions
source(here("R", "dictionaries.R"))
source(here("R", "cleaning_functions.R"))

# 2).Charge spatial locations  
dt_ubicaciones_aire <- fread(here("data", "raw", "Datos_contaminacion", "Estaciones", "datos.csv"))
setDT(dt_ubicaciones_aire)

# 3). LOOP in order to clean the data of the years of interest.
# Strategy: for each year we first look for a single annual CSV (e.g. calidad_2025.csv).
# If it does not exist we look for a monthly subfolder (e.g. Anio19/) and combine the
# files with combinar_meses_anual() before cleaning.
anios_analisis <- c(2019, 2023, 2024, 2025)
lista_historico_aire <- list()

# Base folder that contains both the annual CSVs and the monthly subfolders
carpeta_datos <- here("data", "raw", "Datos_contaminacion", "Datos")

for (anio in anios_analisis) {
  
  cat("\n── Año", anio, "──────────────────────────────────────────\n")
  
  # --- Option A: single annual CSV ------------------------------------------
  ruta_csv_anual <- file.path(carpeta_datos, paste0("calidad_", anio, ".csv"))
  
  # --- Option B: monthly subfolder (e.g. Anio19, Anio23, Anio24) ------------
  sufijo_anio    <- substr(as.character(anio), 3, 4)   # "19", "23", "24", "25"
  ruta_carpeta   <- file.path(carpeta_datos, paste0("Anio", sufijo_anio))
  
  dt_temporal <- tryCatch({
    
    if (file.exists(ruta_csv_anual)) {
      cat("  Fuente: archivo anual único ->", basename(ruta_csv_anual), "\n")
      fread(ruta_csv_anual, sep = ";")
      
    } else if (dir.exists(ruta_carpeta)) {
      cat("  Fuente: carpeta mensual     ->", basename(ruta_carpeta), "\n")
      combinar_meses_anual(ruta_carpeta)
      
    } else {
      cat("  ❌ No se encontró ni archivo anual ni carpeta mensual para", anio, "\n")
      cat("     Buscado en:\n")
      cat("       -", ruta_csv_anual, "\n")
      cat("       -", ruta_carpeta,   "\n")
      NULL
    }
    
  }, error = function(e) {
    cat("  ❌ Error al leer los datos:", conditionMessage(e), "\n")
    NULL
  })
  
  # Skip cleaning if loading failed
  if (is.null(dt_temporal)) next
  
  # Clean and store (keep only rows belonging to the target year)
  dt_limpio <- limpiar_aire_madrid(dt_temporal, dt_ubicaciones_aire)
  dt_limpio <- dt_limpio[year(FECHA) == anio]
  lista_historico_aire[[as.character(anio)]] <- dt_limpio
  cat("  ✅ Limpieza completada:", nrow(dt_limpio), "filas |",
      uniqueN(dt_limpio$ESTACION), "estaciones\n")
}

# 5. Guardado seguro en disco ---------------------------------------------
# Solo entramos a guardar si logramos limpiar algún año
if (length(lista_historico_aire) > 0) {
  for (anio in names(lista_historico_aire)) {
    dt_limpio <- lista_historico_aire[[anio]]
    
    # Guardamos cada año limpio de forma independiente en 'data/processed'
    saveRDS(dt_limpio, here("data", "processed", paste0("aire_madrid_", anio, "_limpio.rds")))
    cat("💾 Guardado con éxito el año:", anio, "en data/processed/\n")
  }
} else {
  cat("⚠️ No se ha guardado nada porque la lista de datos está vacía. Revisa las rutas de arriba.\n")
}


# ==============================================================================
# NO2: FILTRO HORARIO, AGREGACIÓN DIARIA Y LOG-TRANSFORMACIÓN (todos los años)
# ==============================================================================
# We work directly from lista_historico_aire so this block runs immediately
# after the loading loop without needing to re-read any file from disk.

for (anio in names(lista_historico_aire)) {
  
  cat("\n── NO2 procesando año", anio, "─────────────────────────────\n")
  
  dt_aire <- lista_historico_aire[[anio]]
  
  # 1. Filter to NO2 ----------------------------------------------------------
  dt_no2 <- dt_aire[MAGNITUD == "NO2"]
  cat("  Registros horarios de NO2:", nrow(dt_no2), "\n")
  
  if (nrow(dt_no2) == 0) {
    cat("  ⚠️  Sin datos de NO2 para", anio, "— se omite.\n")
    next
  }
  dt_no2[!is.na(DATO), LOG_NO2_HORARIO := log(DATO + 1)]
  
  # Save the raw hourly NO2 dataset
  saveRDS(dt_no2, here("data", "processed","contaminacion","horario",paste0("aire_madrid_", anio, "_No2_horarios.rds")))
  
  # 2. Hourly -> daily (daily mean; days with > 20 % NAs become NA) -----------
  dt_no2_diario <- agregar_a_diario(
    dt_no2,
    col_grupo    = "ESTACION",
    col_fecha    = "FECHA",
    col_valor    = "DATO",
    col_longitud = "LONGITUD",
    col_latitud  = "LATITUD",
    umbral_na    = 0.2
  )
  
  # 3. Log transformation: log(x + 1) to approximate normality for INLA ------
  dt_no2_diario[!is.na(DATO_DIARIO), LOG_NO2_DIARIO := log(DATO_DIARIO + 1)]
  
  # 4. Sequential time index (ID_TIEMPO) for the AR1 component in INLA -------
  setorder(dt_no2_diario, FECHA)
  dt_no2_diario[, ID_TIEMPO := .GRP, by = FECHA]
  
  # 5. Quality-control summary ------------------------------------------------
  cat("  📊 Resumen del dataset diario:\n")
  cat("     Estaciones × Días  :", nrow(dt_no2_diario), "\n")
  cat("     Días válidos        :", sum(!is.na(dt_no2_diario$LOG_NO2_DIARIO)), "\n")
  cat("     Días con NA (INLA)  :", sum( is.na(dt_no2_diario$LOG_NO2_DIARIO)), "\n")
  print(head(dt_no2_diario[, .(ESTACION, FECHA, DATO_DIARIO, LOG_NO2_DIARIO,
                               ID_TIEMPO, LONGITUD, LATITUD)], 5))
  
  # 6. Save final daily dataset -----------------------------------------------
  saveRDS(
    dt_no2_diario,
    here("data", "processed","contaminacion","diario", paste0("aire_madrid_", anio, "_No2_trans_diarios.rds"))
  )
  cat("  💾 Guardado: aire_madrid_", anio, "_log_No2_trans_diarios.rds\n", sep = "")
  
  # 7. Hourly -> monthly (daily mean per station-month; months with > 20 % NAs -> NA)
  # Usamos la función genérica convertir_resolucion() de cleaning_functions.R
  # que acepta columnas de grupo y valor libremente.
  dt_no2_mensual <- convertir_resolucion(
    dt          = dt_no2,
    a           = "mensual",
    col_fecha   = "FECHA",
    cols_grupo  = c("ESTACION", "LONGITUD", "LATITUD", "NOM_TIPO"),
    cols_valores = "DATO",
    umbral_na   = 0.2
  )
  
  # Opcional: log-transformación mensual consistente con la diaria
  dt_no2_mensual[!is.na(DATO), LOG_NO2_MENSUAL := log(DATO + 1)]
  setnames(dt_no2_mensual, "DATO", "DATO_MENSUAL")
  
  saveRDS(
    dt_no2_mensual,
    here("data", "processed","contaminacion","mensual", paste0("aire_madrid_", anio, "_log_No2_mensuales.rds"))
  )
  cat("  💾 Guardado: aire_madrid_", anio, "_log_No2_mensuales.rds\n", sep = "")
}

cat("\n✅ Pipeline NO2 completado para todos los años (horario, diario y mensual).\n")

## Checking all the cleaning data 

check_1<-readRDS(here("data", "processed","contaminacion","horario","aire_madrid_2025_No2_horarios.rds"))
View(check_1)

check_2<-readRDS(here("data", "processed","contaminacion","diario", "aire_madrid_2025_No2_trans_diarios.rds"))
View(check_2)

check_3<-readRDS(here("data", "processed","contaminacion","mensual", "aire_madrid_2025_log_No2_mensuales.rds"))
View(check_3)

