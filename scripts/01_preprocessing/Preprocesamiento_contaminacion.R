#------------------------------------------------------------------------
# Goal: load the raw data of air pollution, clean it, and save it by year.
#------------------------------------------------------------------------
library(here)
library(data.table)
library(lubridate)

# 1). Charge the dictionaries and cleaning functions
source(here("R", "utilities", "dictionaries.R"))
source(here("R", "cleaning", "cleaning_functions.R"))

# 2).Charge spatial locations  
dt_ubicaciones_aire <- fread(here("data", "raw", "Datos_contaminacion", "Estaciones", "datos.csv"))
setDT(dt_ubicaciones_aire)

# ------------------------------------------------------------------------
# Version de la salida
# ------------------------------------------------------------------------
# Se escribe en ficheros NUEVOS para no pisar los existentes: al cambiar una
# regla del preprocesamiento, los resultados anteriores siguen disponibles
# para comparar y cualquier cifra ya publicada sigue siendo rastreable.
# Poner "" para volver a escribir sobre los originales.
SUFIJO <- "1"

# Con TRUE el script se detiene si el fichero de destino ya existe, en lugar
# de sobrescribirlo sin avisar.
PROTEGER_EXISTENTES <- TRUE

# Helper unico para construir las rutas de salida: asi el sufijo no se puede
# olvidar en una de las cuatro llamadas a saveRDS.
ruta_salida <- function(..., nombre) {
  here("data", "processed", ..., paste0(nombre, SUFIJO, ".rds"))
}

guardar <- function(objeto, ruta) {
  if (PROTEGER_EXISTENTES && file.exists(ruta)) {
    stop("Ya existe el fichero de destino:
  ", ruta,
         "
Cambia SUFIJO o pon PROTEGER_EXISTENTES <- FALSE.")
  }
  # La carpeta destino puede no existir todavia (p. ej. al estrenar una
  # subcarpeta): saveRDS no la crea y fallaria.
  dir.create(dirname(ruta), recursive = TRUE, showWarnings = FALSE)
  saveRDS(objeto, ruta)
  cat("   guardado ->", basename(ruta), "
")
}

# 3). LOOP in order to clean the data of the years of interest.
# Strategy: for each year we first look for a single annual CSV (e.g. calidad_2025.csv).
# If it does not exist we look for a monthly subfolder (e.g. Anio19/) and combine the
# files with combinar_meses_anual() before cleaning.
anios_analisis <- c(2019,2020,2021,2022,2023,2024,2025)
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

  # --- Rejilla horaria completa --------------------------------------------
  # El proveedor no publica fila cuando la estacion esta caida el dia entero:
  # no escribe 24 registros vacios, sencillamente no escribe nada. Esos dias no
  # son NA, no existen, y desaparecerian del denominador sin dejar rastro.
  #
  # La rejilla se construye sobre los PARES estacion-magnitud OBSERVADOS, no
  # sobre el producto cartesiano de estaciones por magnitudes: ninguna estacion
  # mide los 8 contaminantes (solo 126 de los 192 pares posibles existen), y
  # cruzarlos todos inventaria ~578.000 filas de analizadores inexistentes.
  # Al estar el dato en formato LARGO, un contaminante que la estacion no mide
  # simplemente no genera filas, y no hace falta un estado "SIN_SENSOR".
  pares <- unique(dt_limpio[, .(ESTACION, MAGNITUD)])

  dias_anio <- seq(as.Date(paste0(anio, "-01-01")),
                   as.Date(paste0(anio, "-12-31")), by = "day")

  horas <- if (is.factor(dt_limpio$HORA)) {
    factor(levels(dt_limpio$HORA), levels = levels(dt_limpio$HORA))
  } else {
    sort(unique(dt_limpio$HORA))
  }

  rejilla <- pares[, CJ(FECHA = dias_anio, HORA = horas),
                   by = .(ESTACION, MAGNITUD)]

  # El right join conserva las filas presentes y crea las que faltan con NA.
  dt_limpio[, fila_ausente := FALSE]
  dt_limpio <- dt_limpio[rejilla, on = .(ESTACION, MAGNITUD, FECHA, HORA)]
  dt_limpio[is.na(fila_ausente), fila_ausente := TRUE]

  # Los atributos de estacion son constantes y se recuperan de sus otras filas.
  cols_est <- intersect(c("LONGITUD", "LATITUD", "NOM_TIPO"), names(dt_limpio))
  dt_limpio[, (cols_est) := lapply(.SD, function(x) x[!is.na(x)][1L]),
            by = ESTACION, .SDcols = cols_est]

  # Motivo por el que una celda esta vacia. En formato largo bastan tres
  # estados: no hay imputacion, y la ausencia de analizador no genera filas.
  dt_limpio[, ESTADO := factor(
    fifelse(fila_ausente,  "AUSENTE",     # la fila no existia en el crudo
    fifelse(is.na(DATO),   "FALLO",       # existia, invalidada por el flag "V"
                           "OK")),        # medida valida
    levels = c("OK", "FALLO", "AUSENTE")
  )]

  cat(sprintf("   Rejilla completada: +%d filas (%d dias estacion-magnitud)\n",
              sum(dt_limpio$fila_ausente),
              sum(dt_limpio$fila_ausente) / length(horas)))

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
    guardar(dt_limpio, ruta_salida("Contaminacion", "Con_8_magnitudes_primer_pass",
                                   nombre = paste0("aire_madrid_", anio, "_limpio")))
    cat("   ano", anio, "guardado en Contaminacion/limpio
")
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

  # De aqui en adelante se da por hecho que la tabla contiene UN SOLO
  # contaminante. Las agregaciones agrupan por estacion y fecha, sin mencionar
  # MAGNITUD: con dos magnitudes promediarian ambas juntas y devolverian un
  # valor plausible pero sin significado, sin lanzar ningun error.
  stopifnot(uniqueN(dt_no2$MAGNITUD) == 1L)

  dt_no2[!is.na(DATO), LOG_NO2_HORARIO := log(DATO + 1)]
  
  # Save the raw hourly NO2 dataset
  guardar(dt_no2, ruta_salida("Contaminacion", "horario",
                              nombre = paste0("aire_madrid_", anio, "_No2_horarios")))
  
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
  
  # 4. Ordenacion temporal ---------------------------------------------------
  # NO se crea aqui el indice temporal para el AR1 de INLA. El indice debe
  # construirse sobre la MUESTRA que entra al modelo, no sobre el dato
  # completo: inla.spde.make.A() recibe n.group = numero de periodos y
  # group = los valores del indice, y espera que estos vayan de 1 a n.group.
  #
  # Si el modelo filtra la muestra (un tramo de fechas, un holdout temporal),
  # un indice heredado del preprocesamiento tendria huecos: con los dias 267
  # a 300 quedaria n.group = 34 pero group = 267..300, y la matriz apuntaria
  # a grupos inexistentes SIN dar error.
  #
  # En cada script de modelado, despues de aplicar todos los filtros:
  #   setorder(dt, FECHA, ESTACION)
  #   dt[, ID_TIEMPO := .GRP, by = FECHA]
  #   stopifnot(identical(sort(unique(dt$ID_TIEMPO)), seq_len(uniqueN(dt$FECHA))))
  setorder(dt_no2_diario, FECHA)

  # 5. Quality-control summary ------------------------------------------------
  cat("  📊 Resumen del dataset diario:\n")
  cat("     Estaciones × Días  :", nrow(dt_no2_diario), "\n")
  cat("     Días válidos        :", sum(!is.na(dt_no2_diario$LOG_NO2_DIARIO)), "\n")
  cat("     Días con NA (INLA)  :", sum( is.na(dt_no2_diario$LOG_NO2_DIARIO)), "\n")
  print(head(dt_no2_diario[, .(ESTACION, FECHA, DATO_DIARIO, LOG_NO2_DIARIO,
                               LONGITUD, LATITUD)], 5))
  
  # 6. Save final daily dataset -----------------------------------------------
  guardar(
    dt_no2_diario,
    ruta_salida("Contaminacion", "diario",
                nombre = paste0("aire_madrid_", anio, "_No2_trans_diarios"))
  )
  cat("  💾 Guardado: aire_madrid_", anio, "_log_No2_trans_diarios.rds\n", sep = "")
  
  # 7. Diario -> mensual (media de las medias diarias; mes NA si faltan >= 20 %
  #    de los dias) ------------------------------------------------------------
  # Se parte del DIARIO, no del horario. Asi el valor mensual es la media de las
  # medias diarias y cada dia pesa lo mismo, con independencia de cuantas horas
  # validas tuviera. Promediando las ~720 horas del mes, un dia de 24 horas
  # aportaba mas que uno de 20, y el mensual no coincidia con el promedio de la
  # serie diaria (discrepancia media observada en 2019: 0,065 ug/m3).
  #
  # El umbral es PROPORCIONAL —sum(is.na(x)) / .N—, de modo que se adapta solo a
  # meses de 28, 29, 30 o 31 dias sin necesidad de un recuento fijo.
  #
  # Los dias que no superaron el umbral diario llegan aqui como NA y cuentan
  # como faltantes del mes, encadenando ambos criterios de forma coherente.
  dt_no2_mensual <- convertir_resolucion(
    dt           = dt_no2_diario,
    a            = "mensual",
    col_fecha    = "FECHA",
    cols_grupo   = c("ESTACION", "LONGITUD", "LATITUD"),
    cols_valores = "DATO_DIARIO",
    umbral_na    = 0.2
  )

  setnames(dt_no2_mensual, "DATO_DIARIO", "DATO_MENSUAL")

  # NOM_TIPO se pierde en la agregacion diaria; se recupera de la tabla horaria,
  # donde es constante por estacion.
  tipos <- unique(dt_no2[, .(ESTACION, NOM_TIPO)])
  dt_no2_mensual <- merge(dt_no2_mensual, tipos, by = "ESTACION", all.x = TRUE)

  # Log-transformacion mensual, consistente con la diaria: aplicada DESPUES de
  # promediar.
  dt_no2_mensual[!is.na(DATO_MENSUAL), LOG_NO2_MENSUAL := log(DATO_MENSUAL + 1)]
  
  guardar(
    dt_no2_mensual,
    ruta_salida("Contaminacion", "mensual",
                nombre = paste0("aire_madrid_", anio, "_log_No2_mensuales"))
  )
  cat("  💾 Guardado: aire_madrid_", anio, "_log_No2_mensuales.rds\n", sep = "")
}

cat("\n✅ Pipeline NO2 completado para todos los años (horario, diario y mensual).\n")

## Checking all the cleaning data 

check_1<-readRDS(here("data", "processed","Contaminacion","horario","aire_madrid_2025_No2_horarios.rds"))
View(check_1)

check_2<-readRDS(here("data", "processed","Contaminacion","diario", "aire_madrid_2025_No2_trans_diarios.rds"))
View(check_2)

check_3<-readRDS(here("data", "processed","Contaminacion","mensual", "aire_madrid_2025_log_No2_mensuales.rds"))
View(check_3)

check_4<-readRDS(here("data", "processed","aire_madrid_2025_limpio.rds"))
View(check_4)
