# ==============================================================================
# 33 — VARIABILIDAD ESPACIAL DE LAS VARIABLES CLIMÁTICAS DENTRO DE MADRID
# Escala diaria · 2019–2025 · red meteorológica municipal (SIN interpolar)
# ==============================================================================
#
# PREGUNTA: en un mismo día, ¿las estaciones meteorológicas de Madrid miden
# valores parecidos o distintos entre sí?
#
# PARA QUÉ: en el modelo espacio-temporal de NO2, una covariable cuyo campo es
# espacialmente homogéneo dentro del municipio puede entrar como UN valor de
# ciudad por día (media de la red) sin perder apenas información, mientras que
# una covariable con estructura espacial exige interpolarla a la localización de
# cada estación de NO2 (y esa interpolación aporta su propio error). Este script
# separa unas de otras, variable por variable.
#
# CUATRO BLOQUES (más un diagnóstico transversal, 6b):
#   1. Boxplot de la desviación típica ENTRE estaciones, un valor por día.
#   2. Series diarias de todas las estaciones superpuestas (una faceta por
#      variable) con la media de ciudad resaltada.
#   3. Mapas de un día representativo: valor crudo de cada estación, sin
#      interpolar (se quiere ver el dato, no una superficie inventada).
#   4. Comparación por ZONA GEOGRÁFICA (centro / corona / periferia abierta).
#      Ojo: es una clasificación geográfica de las estaciones de CLIMA; no tiene
#      nada que ver con la tipología tráfico/fondo/suburbana del NO2.
#  6b. Diagnóstico que separa el DESFASE FIJO de cada estación (altitud,
#      exposición del sensor) de la DIVERGENCIA que cambia día a día. Es lo que
#      decide si una covariable hay que interpolarla de verdad. Crítico en la
#      presión barométrica: las 8 estaciones con barómetro van de 586 m a 706 m
#      y ese desnivel son ~13 hPa de diferencia constante entre ellas, del mismo
#      orden que toda la variación temporal de la presión en Madrid.
#
# DATOS: data/processed/Clima/diario/meteo_madrid_<AÑO>_diario5.rds
#   Columnas: ESTACION, FECHA, LONGITUD, LATITUD, X_km, Y_km, las seis variables
#   y una columna <variable>_estado por variable
#   (OK / IMPUTADO / FALLO / AUSENTE / SIN_SENSOR).
#   Los valores diarios son los que publica el proveedor (D01..D31): no hay
#   agregación horaria->diaria, así que Precipitaciones ya es el acumulado del
#   día.
#
# DECISIÓN IMPORTANTE — solo se usan observaciones con estado OK. Los valores
# IMPUTADOS se construyen a partir de otras estaciones/variables, de modo que
# incluirlos REDUCIRÍA artificialmente la dispersión entre estaciones, que es
# justo lo que se quiere medir.
#
# COORDENADAS — se usan X_km / Y_km (ETRS89 UTM 30N, EPSG:25830, en km), NO
# LONGITUD/LATITUD: en el RDS estas últimas vienen corruptas para algunas
# estaciones (p. ej. Centro Mpal. De Acústica, E.D.A.R La China).
#
# Salidas: outputs/figures/EDA/variabilidad_espacial_clima/
#          outputs/tables/variabilidad_espacial_clima/
#
# Instalación de librerías, si faltara alguna:
#   install.packages(c("data.table", "ggplot2", "sf", "patchwork",
#                      "scales", "here"))
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
  library(patchwork)
  library(scales)
  library(here)
})

source(here("R", "utilities", "dictionaries.R"))


# ==============================================================================
# 1. CONFIGURACIÓN
# ==============================================================================

ANIOS <- 2019:2025 # periodo completo del TFM
ANIO_ZOOM <- 2025 # año que se dibuja "en detalle" (series legibles)
SOLO_OK <- TRUE # descartar IMPUTADO/FALLO/AUSENTE/SIN_SENSOR

# Nº mínimo de estaciones con dato para que un día cuente en la dispersión.
# Con 2 estaciones la sd está definida pero es puro ruido de un solo par.
MIN_ESTACIONES_DIA <- 3L

# Día para los mapas. NULL => se elige automáticamente el de mayor cobertura.
# Para fijar uno a mano: DIA_MAPA <- as.Date("2025-07-15")
DIA_MAPA <- NULL

# Umbrales de lectura del bloque 1 (sd espacial expresada en unidades de la sd
# total de la variable). Son una regla de trabajo declarada, no un test formal.
UMBRAL_HOMOGENEA <- 0.15 # por debajo: una media de ciudad resume bien
UMBRAL_HETEROGENEA <- 0.35 # por encima: conserva estructura espacial

VARIABLES <- c(
  "Temperatura", "Humedad_Relativa", "Precipitaciones",
  "Velocidad_Viento", "Radiacion_Solar", "Presion_Barometrica"
)

ETIQUETAS <- c(
  Temperatura         = "Temperatura (°C)",
  Humedad_Relativa    = "Humedad relativa (%)",
  Precipitaciones     = "Precipitación (mm/día)",
  Velocidad_Viento    = "Velocidad del viento (m/s)",
  Radiacion_Solar     = "Radiación solar (W/m²)",
  Presion_Barometrica = "Presión barométrica (hPa)"
)

# Paletas: una por variable, porque las unidades no son comparables entre mapas.
OPCION_COLOR <- c(
  Temperatura = "inferno", Humedad_Relativa = "viridis",
  Precipitaciones = "viridis", Velocidad_Viento = "viridis",
  Radiacion_Solar = "plasma", Presion_Barometrica = "cividis"
)

# Formato de la etiqueta numérica sobre cada punto de los mapas.
FORMATO <- c(
  Temperatura = "%.1f", Humedad_Relativa = "%.0f", Precipitaciones = "%.1f",
  Velocidad_Viento = "%.1f", Radiacion_Solar = "%.0f",
  Presion_Barometrica = "%.0f"
)

# ---- Zonas geográficas (bloque 4) --------------------------------------------
# Clasificación MANUAL por CODIGO_CORTO de la estación (no por nombre: así no
# dependemos de acentos ni de la codificación del fichero). Para reasignar una
# estación basta con cambiar su grupo aquí.
#
# Criterio: no es puramente geométrico. "Periferia abierta" son las estaciones
# en entorno no urbanizado o de borde (monte de El Pardo, parque forestal de la
# Casa de Campo, parque Juan Carlos I, borde de Ensanche de Vallecas), donde no
# cabe esperar isla de calor y el anemómetro está más expuesto. Casa de Campo
# está a solo 3,7 km de Sol y aun así entra ahí: por eso la agrupación es
# manual y no un simple anillo de distancia.
ZONAS_MANUAL <- c(
  # Centro urbano denso (almendra central, ~<2 km de Sol)
  "35"  = "Centro urbano", # Plaza del Carmen
  "110" = "Centro urbano", # J.M.D Centro
  "4"   = "Centro urbano", # Plaza España
  "109" = "Centro urbano", # J.M.D Chamberí
  "8"   = "Centro urbano", # Escuelas Aguirre
  # Corona urbana consolidada
  "38"  = "Corona urbana", # Cuatro Caminos
  "111" = "Corona urbana", # J.M.D Chamartín
  "112" = "Corona urbana", # J.M.D Vallecas 1
  "113" = "Corona urbana", # J.M.D Vallecas 2
  "114" = "Corona urbana", # Matadero 01
  "115" = "Corona urbana", # Matadero 02
  "56"  = "Corona urbana", # Plaza Elíptica
  "18"  = "Corona urbana", # Farolillo
  "36"  = "Corona urbana", # Moratalaz
  "102" = "Corona urbana", # J.M.D Moratalaz
  "104" = "Corona urbana", # E.D.A.R La China
  "103" = "Corona urbana", # J.M.D Villaverde
  "106" = "Corona urbana", # Centro Mpal. De Acústica
  "107" = "Corona urbana", # J.M.D Hortaleza
  "108" = "Corona urbana", # Peñagrande
  "39"  = "Corona urbana", # Barrio del Pilar
  "16"  = "Corona urbana", # Arturo Soria
  # Periferia abierta / forestal
  "58"  = "Periferia abierta", # El Pardo
  "24"  = "Periferia abierta", # Casa de Campo
  "59"  = "Periferia abierta", # Juan Carlos I
  "54"  = "Periferia abierta" # Ensanche de Vallecas
)

NIVELES_ZONA <- c("Centro urbano", "Corona urbana", "Periferia abierta")
COLOR_ZONA <- c(
  "Centro urbano"     = "#B2182B",
  "Corona urbana"     = "#E08214",
  "Periferia abierta" = "#2166AC"
)

# Respaldo automático para cualquier estación no listada arriba: anillos de
# distancia a Puerta del Sol (EPSG:25830, metros).
CENTRO_XY <- c(440300, 4474300)
RADIO_CENTRO_KM <- 3
RADIO_PERIFERIA_KM <- 8

DIR_FIG <- here("outputs", "figures", "EDA", "variabilidad_espacial_clima")
DIR_TAB <- here("outputs", "tables", "variabilidad_espacial_clima")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. UTILIDADES
# ==============================================================================

tema_base <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(colour = "grey35", size = base_size - 0.5),
      plot.caption = element_text(
        colour = "grey35", size = base_size - 1.5, hjust = 0, lineheight = 1.2
      ),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      strip.background = element_rect(fill = "grey96", colour = "grey85"),
      panel.grid.minor = element_blank()
    )
}

tema_mapa <- function(base_size = 10) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size, hjust = 0.5),
      plot.subtitle = element_text(
        size = base_size - 2, hjust = 0.5, colour = "grey40"
      ),
      legend.position = "right",
      legend.title = element_text(size = base_size - 2),
      legend.text = element_text(size = base_size - 2.5),
      legend.key.height = unit(0.75, "cm"),
      legend.key.width = unit(0.35, "cm")
    )
}

# Guarda a 300 dpi y además muestra la figura en el dispositivo activo.
guardar <- function(nombre, grafico, ancho, alto) {
  ruta <- file.path(DIR_FIG, nombre)
  ggsave(ruta, grafico, width = ancho, height = alto, dpi = 300, bg = "white")
  cat(sprintf("  [fig] %s\n", nombre))
  print(grafico)
  invisible(ruta)
}


# ==============================================================================
# 3. CARGA DE DATOS 2019–2025
# ==============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("VARIABILIDAD ESPACIAL DEL CLIMA EN MADRID — escala diaria\n")
cat(strrep("=", 78), "\n\n", sep = "")

# Cada año es un fichero independiente; se leen y se apilan. Si algún año no
# tiene alguna de las seis variables, se rellena con NA para poder apilar.
leer_anio <- function(anio) {
  ruta <- here(
    "data", "processed", "Clima", "diario",
    sprintf("meteo_madrid_%d_diario5.rds", anio)
  )
  if (!file.exists(ruta)) {
    cat(sprintf(
      "  [AVISO] %d: no existe %s — año omitido\n", anio, basename(ruta)
    ))
    return(NULL)
  }

  dt <- as.data.table(readRDS(ruta))
  dt[, FECHA := as.Date(FECHA)]

  # Se anula lo que no está marcado como OK (ver cabecera).
  if (SOLO_OK) {
    for (v in intersect(VARIABLES, names(dt))) {
      col_estado <- paste0(v, "_estado")
      if (col_estado %in% names(dt)) {
        dt[as.character(get(col_estado)) != "OK", (v) := NA_real_]
      }
    }
  }

  faltan <- setdiff(VARIABLES, names(dt))
  if (length(faltan)) {
    cat(sprintf(
      "  [AVISO] %d: sin columnas %s (se rellenan con NA)\n",
      anio, paste(faltan, collapse = ", ")
    ))
    dt[, (faltan) := NA_real_]
  }

  columnas_xy <- intersect(c("X_km", "Y_km"), names(dt))
  dt[, c("ESTACION", "FECHA", VARIABLES, columnas_xy), with = FALSE]
}

datos <- rbindlist(lapply(ANIOS, leer_anio), use.names = TRUE, fill = TRUE)
if (!nrow(datos)) stop("No se ha podido cargar ningún año.")

datos[, ESTACION := enc2utf8(as.character(ESTACION))]
setorder(datos, ESTACION, FECHA)

cat(sprintf(
  "\nCargadas %s filas · %d estaciones · %s a %s\n",
  format(nrow(datos), big.mark = "."), uniqueN(datos$ESTACION),
  format(min(datos$FECHA), "%d/%m/%Y"), format(max(datos$FECHA), "%d/%m/%Y")
))


# ==============================================================================
# 4. COORDENADAS Y LÍMITE DEL MUNICIPIO
# ==============================================================================

# X_km/Y_km ya vienen en EPSG:25830 (km). Se pasan a metros para trabajar en sf.
if (all(c("X_km", "Y_km") %in% names(datos))) {
  coords <- unique(
    datos[
      !is.na(X_km) & !is.na(Y_km),
      .(ESTACION, X_m = X_km * 1000, Y_m = Y_km * 1000)
    ],
    by = "ESTACION"
  )
} else {
  # Respaldo: catálogo oficial de estaciones (COORDENADA_X/Y_ETRS89, metros).
  cat("  [AVISO] El RDS no trae X_km/Y_km; se leen del catálogo de estaciones\n")
  catalogo <- fread(
    here(
      "data", "raw", "Datos metereologicos", "Estaciones_2019", "estaciones.csv"
    ),
    sep = ";", encoding = "Latin-1"
  )
  num_utm <- function(x) as.numeric(gsub(",", ".", as.character(x)))
  coords <- unique(
    catalogo[, .(
      ESTACION = unname(nombres_estaciones_clima[as.character(CODIGO_CORTO)]),
      X_m = num_utm(COORDENADA_X_ETRS89),
      Y_m = num_utm(COORDENADA_Y_ETRS89)
    )][!is.na(ESTACION) & is.finite(X_m)],
    by = "ESTACION"
  )
}

sin_coord <- setdiff(unique(datos$ESTACION), coords$ESTACION)
if (length(sin_coord)) {
  cat(sprintf(
    "  [AVISO] %d estación(es) sin coordenadas: %s\n          Entran en los bloques 1-2 pero no en mapas ni zonas.\n",
    length(sin_coord), paste(sin_coord, collapse = ", ")
  ))
}

# ALTITUD. Imprescindible para la presión barométrica: las estaciones con
# barómetro van de 586 m a 706 m, y a esa altura la presión cae ~0,11 hPa por
# metro, así que 120 m de desnivel son ~13 hPa de diferencia SISTEMÁTICA entre
# estaciones. Sin este dato no se puede distinguir la estructura espacial real
# del simple desnivel topográfico (ver bloque 6b).
catalogo_est <- fread(
  here(
    "data", "raw", "Datos metereologicos", "Estaciones_2019", "estaciones.csv"
  ),
  sep = ";", encoding = "Latin-1"
)
altitudes <- unique(
  catalogo_est[, .(
    ESTACION = unname(nombres_estaciones_clima[as.character(CODIGO_CORTO)]),
    ALTITUD  = as.numeric(ALTITUD)
  )][!is.na(ESTACION) & is.finite(ALTITUD)],
  by = "ESTACION"
)
coords <- merge(coords, altitudes, by = "ESTACION", all.x = TRUE)
if (coords[is.na(ALTITUD), .N]) {
  cat(sprintf(
    "  [AVISO] %d estación(es) sin altitud en el catálogo: %s\n",
    coords[is.na(ALTITUD), .N],
    paste(coords[is.na(ALTITUD), ESTACION], collapse = ", ")
  ))
}

# Límite municipal: unión de los distritos. Es el polígono real del municipio;
# solo si no existe la capa se recurre a la envolvente convexa de la red.
ruta_distritos <- here("data", "raw", "geometrias", "DISTRITOS.shp")
if (file.exists(ruta_distritos)) {
  distritos <- st_transform(st_read(ruta_distritos, quiet = TRUE), 25830)
  limite_madrid <- st_union(st_geometry(distritos))
} else {
  cat("  [AVISO] No se encuentra DISTRITOS.shp: el límite se aproxima con la\n")
  cat("          envolvente convexa de las estaciones (mapas menos fieles).\n")
  limite_madrid <- st_convex_hull(st_union(st_geometry(
    st_as_sf(coords, coords = c("X_m", "Y_m"), crs = 25830)
  )))
  distritos <- st_sf(geometry = limite_madrid)
}


# ==============================================================================
# 5. ZONAS GEOGRÁFICAS DE LAS ESTACIONES CLIMÁTICAS
# ==============================================================================
# Se traduce el diccionario código -> nombre para poder aplicar ZONAS_MANUAL.
# Las estaciones no listadas caen en el respaldo por distancia a Sol, y se avisa
# por consola para que puedan añadirse a mano si interesa.

zona_por_nombre <- setNames(
  unname(ZONAS_MANUAL),
  unname(nombres_estaciones_clima[names(ZONAS_MANUAL)])
)
zona_por_nombre <- zona_por_nombre[!is.na(names(zona_por_nombre))]

estaciones <- copy(coords)
estaciones[, dist_centro_km := sqrt(
  (X_m - CENTRO_XY[1])^2 + (Y_m - CENTRO_XY[2])^2
) / 1000]
estaciones[, Zona := zona_por_nombre[ESTACION]]
estaciones[, asignacion := fifelse(is.na(Zona), "automatica", "manual")]

# Respaldo automático por anillos de distancia.
estaciones[is.na(Zona), Zona := fcase(
  dist_centro_km < RADIO_CENTRO_KM, "Centro urbano",
  dist_centro_km <= RADIO_PERIFERIA_KM, "Corona urbana",
  default = "Periferia abierta"
)]
estaciones[, Zona := factor(Zona, levels = NIVELES_ZONA)]
setorder(estaciones, Zona, dist_centro_km)

n_auto <- estaciones[asignacion == "automatica", .N]
if (n_auto > 0) {
  cat(sprintf(
    "\n  [AVISO] %d estación(es) sin zona manual, asignadas por distancia a Sol:\n          %s\n",
    n_auto,
    paste(estaciones[asignacion == "automatica", ESTACION], collapse = ", ")
  ))
}

cat("\n--- Zonas geográficas de las estaciones climáticas ---\n")
print(estaciones[, .(
  ESTACION, Zona,
  dist_centro_km = round(dist_centro_km, 2), asignacion
)])
fwrite(
  estaciones[, .(
    ESTACION, Zona,
    dist_centro_km = round(dist_centro_km, 3), asignacion
  )],
  file.path(DIR_TAB, "zonas_estaciones.csv")
)


# ==============================================================================
# 6. FORMATO LARGO Y DISPERSIÓN ESPACIAL DIARIA
# ==============================================================================

largo <- melt(
  datos,
  id.vars = c("ESTACION", "FECHA"), measure.vars = VARIABLES,
  variable.name = "Variable", value.name = "valor", na.rm = TRUE
)
largo <- largo[is.finite(valor)]

# Escala de referencia de cada variable: su sd TOTAL sobre todo el registro
# (todas las estaciones y todos los días). Sirve para hacer comparables entre sí
# las dispersiones de variables con unidades distintas (ver bloque 7).
escala <- largo[, .(
  media_global = mean(valor),
  sd_global    = sd(valor),
  n_obs        = .N
), by = Variable]

# Dispersión ESPACIAL de cada día: sd entre las estaciones que midieron esa
# variable ese día. Los días con menos de MIN_ESTACIONES_DIA estaciones se
# excluyen del cálculo (no fallan: simplemente no entran).
dispersion <- largo[, .(
  n_est     = .N,
  media_dia = mean(valor),
  sd_dia    = sd(valor)
), by = .(Variable, FECHA)][n_est >= MIN_ESTACIONES_DIA & is.finite(sd_dia)]

dispersion <- merge(dispersion, escala, by = "Variable")

# sd espacial en unidades de la sd total de la variable: adimensional, luego
# comparable entre variables. Se interpreta como "qué fracción de la
# variabilidad total de la variable es variabilidad ENTRE estaciones".
dispersion[, sd_z := sd_dia / sd_global]

# Coeficiente de variación: se calcula como referencia, pero NO se usa para
# ordenar ni para la figura principal. Motivo: el CV exige un cero absoluto y
# una media que no se acerque a cero. La temperatura en °C tiene cero
# convencional (un día a 0,5 °C daría un CV disparatado) y la precipitación
# tiene media 0 en la mayoría de los días. El z-score no tiene ese problema.
dispersion[, cv_dia := fifelse(
  abs(media_dia) > 1e-9, sd_dia / abs(media_dia), NA_real_
)]

# Orden de menor a mayor dispersión mediana. Se aplica a TODAS las figuras para
# que los paneles se lean siempre en el mismo orden.
orden_variables <- dispersion[, .(m = median(sd_z)), by = Variable][
  order(m), as.character(Variable)
]

aplicar_orden <- function(x) {
  factor(
    ETIQUETAS[as.character(x)],
    levels = unname(ETIQUETAS[orden_variables])
  )
}
dispersion[, Variable_lab := aplicar_orden(Variable)]
largo[, Variable_lab := aplicar_orden(Variable)]

# Tabla resumen: es el número que se cita en la memoria.
#
# OJO CON LA MEDIANA EN LA PRECIPITACIÓN. En Madrid la mayoría de los días no
# llueve en ninguna estación: esos días tienen sd = 0 y arrastran la mediana a
# cero, con lo que la precipitación parecería la variable MÁS homogénea de las
# seis. Es un artefacto de resumir con la mediana una variable con exceso de
# ceros. Por eso se añaden tres columnas: el porcentaje de días con sd = 0, el
# percentil 90 y la mediana calculada solo sobre los días "activos" (sd > 0),
# que es la cifra que hay que mirar en precipitación.
resumen <- dispersion[, .(
  N_dias            = .N,
  Est_por_dia       = round(mean(n_est), 1),
  # Media global: es la referencia para juzgar si una sd espacial es grande o
  # pequeña. Sin ella, un "0,4 m/s de diferencia entre estaciones" no se puede
  # interpretar; sabiendo que la media de la ciudad ronda los 2 m/s, sí.
  media_total       = round(first(media_global), 3),
  sd_total          = round(first(sd_global), 3),
  sd_espacial_med   = round(median(sd_dia), 3),
  sd_espacial_p25   = round(quantile(sd_dia, 0.25), 3),
  sd_espacial_p75   = round(quantile(sd_dia, 0.75), 3),
  ratio_med         = round(median(sd_z), 3),
  ratio_p25         = round(quantile(sd_z, 0.25), 3),
  ratio_p75         = round(quantile(sd_z, 0.75), 3),
  ratio_p90         = round(quantile(sd_z, 0.90), 3),
  pct_dias_sd0      = round(100 * mean(sd_dia <= 0), 1),
  ratio_med_activo  = round(median(sd_z[sd_dia > 0]), 3),
  CV_med            = round(median(cv_dia, na.rm = TRUE), 3)
), by = Variable]
setorder(resumen, ratio_med)

# La clasificación usa la mediana sobre días activos cuando hay muchos días con
# dispersión nula, y la mediana normal en el resto. Así la precipitación se
# juzga por los días en que hay algo que repartir espacialmente.
resumen[, ratio_criterio := fifelse(
  pct_dias_sd0 > 25, ratio_med_activo, ratio_med
)]
resumen[, Lectura := fcase(
  ratio_criterio < UMBRAL_HOMOGENEA, "Homogenea: admite un valor de ciudad",
  ratio_criterio <= UMBRAL_HETEROGENEA, "Intermedia: conviene comprobar el efecto",
  default = "Heterogenea: conserva estructura espacial"
)]

fwrite(resumen, file.path(DIR_TAB, "dispersion_espacial_por_variable.csv"))
cat("\n--- Dispersión espacial diaria por variable (ordenada) ---\n")
print(resumen[, .(
  Variable, media_total, sd_total, sd_espacial_med, ratio_med,
  ratio_med_activo, pct_dias_sd0, N_dias, Est_por_dia, Lectura
)])


# ==============================================================================
# 6b. ¿DESFASE FIJO O ESTRUCTURA ESPACIAL? NIVEL FRENTE A ANOMALÍA
# ==============================================================================
# Que dos estaciones midan distinto puede deberse a dos cosas MUY diferentes:
#
#   (a) un DESFASE FIJO de la estación: siempre marca 13 hPa menos porque está
#       120 m más alta, o siempre menos viento porque está resguardada. Es un
#       efecto constante en el tiempo, no un campo espacial que evolucione.
#   (b) una DIVERGENCIA que cambia día a día: hoy llueve en el norte y no en el
#       sur. Esto sí es estructura espacio-temporal.
#
# Solo (b) justifica interpolar la covariable: un desfase fijo lo absorbe
# cualquier intercepto por estación (o el propio campo SPDE) sin necesidad de
# interpolar nada. Para separarlos se recalcula la dispersión sobre la ANOMALÍA
# de cada estación (valor menos la media histórica de esa estación): si al
# centrar cada estación la dispersión se desploma, lo que había era (a).
#
# Nota: el centrado usa la media de cada estación sobre su propio periodo
# disponible; con coberturas muy desiguales entre estaciones eso introduce algo
# de ruido, pero no cambia la lectura cualitativa.

largo[, valor_anom := valor - mean(valor), by = .(ESTACION, Variable)]

dispersion_anom <- largo[, .(
  n_est  = .N,
  sd_dia = sd(valor_anom)
), by = .(Variable, FECHA)][n_est >= MIN_ESTACIONES_DIA & is.finite(sd_dia)]
dispersion_anom <- merge(dispersion_anom, escala, by = "Variable")
dispersion_anom[, sd_z := sd_dia / sd_global]

descomposicion <- merge(
  resumen[, .(Variable,
    ratio_bruto = ratio_med,
    sd_espacial_bruta = sd_espacial_med
  )],
  dispersion_anom[, .(
    ratio_anomalia    = round(median(sd_z), 3),
    sd_espacial_anom  = round(median(sd_dia), 3)
  ), by = Variable],
  by = "Variable"
)
# Qué fracción de la discrepancia entre estaciones es puro desfase fijo.
# Se protege la división: la precipitación tiene ratio_bruto = 0 (mediana sobre
# días mayoritariamente secos) y ahí el porcentaje no está definido.
descomposicion[, pct_desfase_fijo := fifelse(
  ratio_bruto > 0, round(100 * (1 - ratio_anomalia / ratio_bruto), 1), NA_real_
)]
setorder(descomposicion, -pct_desfase_fijo, na.last = TRUE)
fwrite(descomposicion, file.path(DIR_TAB, "nivel_vs_anomalia.csv"))

cat("\n--- ¿Desfase fijo de estación o estructura espacial real? ---\n")
cat("    pct_desfase_fijo alto => las estaciones difieren SIEMPRE IGUAL:\n")
cat("    no hay campo que interpolar, basta un intercepto por estación.\n")
print(descomposicion)

# Diagnóstico de altitud: cuánto de la diferencia de NIVEL entre estaciones
# explica la topografía. Para la presión debería ser prácticamente todo, con
# una pendiente cercana a -0,11 hPa/m; si sale así, su aparente heterogeneidad
# espacial es un artefacto barométrico y no información meteorológica.
medias_estacion <- merge(
  largo[, .(media = mean(valor), n = .N), by = .(ESTACION, Variable)],
  coords[, .(ESTACION, ALTITUD)],
  by = "ESTACION"
)

altitud_efecto <- medias_estacion[
  !is.na(ALTITUD),
  {
    if (.N >= 4L && sd(ALTITUD) > 0) {
      ajuste <- lm(media ~ ALTITUD)
      .(
        N_estaciones = .N,
        pendiente    = round(unname(coef(ajuste)[2]), 4),
        R2           = round(summary(ajuste)$r.squared, 3),
        rango_alt_m  = round(diff(range(ALTITUD)), 0)
      )
    } else {
      .(
        N_estaciones = .N, pendiente = NA_real_, R2 = NA_real_,
        rango_alt_m = NA_real_
      )
    }
  },
  by = Variable
]
setorder(altitud_efecto, -R2, na.last = TRUE)
fwrite(altitud_efecto, file.path(DIR_TAB, "efecto_altitud.csv"))

cat("\n--- Efecto de la altitud sobre el nivel medio de cada estación ---\n")
cat("    (presión: R2 alto y pendiente ~ -0.11 hPa/m => desnivel, no clima)\n")
print(altitud_efecto)


# ==============================================================================
# 6c. RANGOS OBSERVADOS Y REFERENCIA EXTERNA: ¿ES MUCHO O ES POCO?
# ==============================================================================
# Una sd de 0,4 m/s no dice nada sin saber que la media de la red es 1,2 m/s.
# Este bloque sitúa cada variable en una escala interpretable: los percentiles
# diarios de la media de ciudad y el reparto de los días en clases con nombre.
#
# FUENTES DE LOS CORTES (declaradas; conviene citarlas en la memoria):
#   - Viento: escala Beaufort, expresada en medias diarias. Ojo: la red
#     municipal mide en emplazamientos urbanos resguardados y a poca altura, así
#     que lee sistemáticamente por debajo de un observatorio en campo abierto.
#   - Precipitación: clases de acumulado diario de uso habitual en AEMET.
#   - Temperatura, humedad y radiación: cortes convencionales anclados en los
#     valores normales 1991-2020 del observatorio Madrid-Retiro (AEMET).
#   - Presión: NO se usa un valor absoluto, porque depende de si la serie está
#     reducida a nivel del mar (~1013 hPa) o es presión de estación (~946 hPa a
#     640 m). Se clasifica la ANOMALÍA respecto a la mediana de la propia serie,
#     que funciona en los dos casos.
# Los valores "normal_ref" son orientativos: verifícalos contra la publicación
# de AEMET antes de citarlos en el texto.

REFERENCIAS <- list(
  Temperatura = list(
    cortes = c(5, 12, 20, 26),
    clases = c("muy frio (<5)", "frio (5-12)", "templado (12-20)",
               "calido (20-26)", "muy calido (>26)"),
    normal_ref = 15.6, fuente = "AEMET Madrid-Retiro, media anual 1991-2020"
  ),
  Humedad_Relativa = list(
    cortes = c(40, 60, 80),
    clases = c("seco (<40)", "normal (40-60)", "humedo (60-80)",
               "muy humedo (>80)"),
    normal_ref = 57, fuente = "AEMET Madrid-Retiro, media anual"
  ),
  Precipitaciones = list(
    cortes = c(0.1, 2, 10, 30),
    clases = c("seco (0)", "debil (0-2)", "moderada (2-10)",
               "fuerte (10-30)", "muy fuerte (>30)"),
    normal_ref = 415 / 365, fuente = "AEMET Madrid-Retiro, 415 mm/ano"
  ),
  Velocidad_Viento = list(
    cortes = c(0.5, 1.5, 3, 5),
    clases = c("calma (<0.5)", "ventolina (0.5-1.5)", "flojo (1.5-3)",
               "ventoso (3-5)", "muy ventoso (>5)"),
    normal_ref = 2.4, fuente = "Beaufort; observatorio abierto ~2-3 m/s"
  ),
  Radiacion_Solar = list(
    cortes = c(80, 150, 250, 320),
    clases = c("muy baja (<80)", "baja (80-150)", "media (150-250)",
               "alta (250-320)", "muy alta (>320)"),
    normal_ref = 200, fuente = "~4.8 kWh/m2/dia, tipico de Madrid"
  ),
  Presion_Barometrica = list(
    cortes = c(-8, -3, 3, 8), relativo = TRUE,
    clases = c("borrasca marcada", "baja", "normal", "alta",
               "anticiclon marcado"),
    normal_ref = NA_real_,
    fuente = "Anomalia (hPa) respecto a la mediana de la propia serie"
  )
)

# Media diaria de ciudad: es el valor que describe "cómo fue ese día en Madrid".
ciudad_dia <- largo[, .(valor = mean(valor)), by = .(Variable, FECHA)]

rangos <- ciudad_dia[, .(
  N_dias  = .N,
  minimo  = round(min(valor), 2),
  p01     = round(quantile(valor, 0.01), 2),
  p05     = round(quantile(valor, 0.05), 2),
  p25     = round(quantile(valor, 0.25), 2),
  mediana = round(median(valor), 2),
  p75     = round(quantile(valor, 0.75), 2),
  p95     = round(quantile(valor, 0.95), 2),
  p99     = round(quantile(valor, 0.99), 2),
  maximo  = round(max(valor), 2),
  media   = round(mean(valor), 2),
  sd      = round(sd(valor), 2)
), by = Variable]
rangos[, normal_ref := sapply(as.character(Variable), function(v) {
  if (!is.null(REFERENCIAS[[v]])) REFERENCIAS[[v]]$normal_ref else NA_real_
})]
rangos[, fuente_ref := sapply(as.character(Variable), function(v) {
  if (!is.null(REFERENCIAS[[v]])) REFERENCIAS[[v]]$fuente else NA_character_
})]

fwrite(rangos, file.path(DIR_TAB, "rangos_diarios_ciudad.csv"))
cat("\n--- Rangos diarios de la media de ciudad (¿es mucho o es poco?) ---\n")
print(rangos[, .(Variable, minimo, p05, mediana, p95, maximo, media, sd,
                 normal_ref)])

# Reparto de los días en clases con nombre: responde directamente a "¿hay mucho
# viento en Madrid o poco?" con un porcentaje de días, no con una media.
reparto <- rbindlist(lapply(names(REFERENCIAS), function(v) {
  d <- ciudad_dia[Variable == v]
  if (!nrow(d)) return(NULL)
  ref <- REFERENCIAS[[v]]
  x <- if (isTRUE(ref$relativo)) d$valor - median(d$valor) else d$valor
  clase <- cut(x, breaks = c(-Inf, ref$cortes, Inf), labels = ref$clases,
               right = FALSE)
  tabla <- as.data.table(table(Clase = clase))
  tabla[, `:=`(Variable = v, pct_dias = round(100 * N / sum(N), 1))]
  tabla[]
}))

if (nrow(reparto)) {
  fwrite(reparto[, .(Variable, Clase, N_dias = N, pct_dias)],
         file.path(DIR_TAB, "reparto_dias_por_clase.csv"))
  cat("\n--- Reparto de los días por clase (% de días del periodo) ---\n")
  print(dcast(reparto, Variable ~ Clase, value.var = "pct_dias", fill = 0))
}


# ==============================================================================
# 7. BLOQUE 1 — BOXPLOT DE LA DISPERSIÓN ENTRE ESTACIONES
# ==============================================================================
# DECISIÓN DE ESCALA: se estandariza. Comparar en crudo la sd de la temperatura
# (°C) con la del viento (m/s) o la presión (hPa) no significa nada, porque el
# tamaño de la caja lo fijaría la unidad y no la homogeneidad espacial. Se
# divide cada sd diaria por la sd TOTAL de su variable (z-score por variable),
# de modo que el eje pasa a ser "fracción de la variabilidad total de la
# variable que se observa entre estaciones el mismo día". Se descarta el
# coeficiente de variación por el motivo explicado en el bloque 6.

n_cero <- dispersion[sd_z <= 0, .N]
tope <- dispersion[, quantile(sd_z, 0.995)]
n_fuera <- dispersion[sd_z > tope, .N]

pie_box <- paste0(
  "Cada punto de la distribución es UN DÍA: la desviación típica entre las estaciones que midieron esa variable ese día\n",
  "(mínimo ", MIN_ESTACIONES_DIA, " estaciones; los días con menos se excluyen). Escala diaria, ",
  min(ANIOS), "–", max(ANIOS), ", solo observaciones con estado OK.\n",
  "Eje estandarizado: sd del día dividida por la sd total de esa misma variable en todo el registro (z-score por variable),\n",
  "porque las seis variables tienen unidades distintas y las cajas en crudo no serían comparables. Se descarta el coeficiente\n",
  "de variación porque la temperatura en °C no tiene cero absoluto y la precipitación tiene media cero la mayoría de los días.\n",
  "(a) escala lineal recortada al percentil 99,5 (", n_fuera, " días quedan fuera del recorte; las cajas se calculan con todos).\n",
  "(b) misma información en escala logarítmica, que separa los extremos bajos (", n_cero,
  " días con sd = 0 no son representables en log)."
)

p1a <- ggplot(dispersion, aes(Variable_lab, sd_z, fill = Variable_lab)) +
  geom_boxplot(
    outlier.alpha = 0.06, outlier.size = 0.5, width = 0.6, linewidth = 0.35
  ) +
  geom_hline(
    yintercept = c(UMBRAL_HOMOGENEA, UMBRAL_HETEROGENEA),
    linetype = "dashed", colour = "grey55", linewidth = 0.35
  ) +
  stat_summary(
    fun = median, geom = "text",
    aes(label = sprintf("%.2f", after_stat(y))),
    vjust = -0.9, size = 3, colour = "grey15"
  ) +
  coord_cartesian(ylim = c(0, tope)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    title = "(a) Dispersión espacial diaria — escala lineal",
    subtitle = paste0(
      "Líneas discontinuas: umbrales de lectura ", UMBRAL_HOMOGENEA, " y ",
      UMBRAL_HETEROGENEA, ". La cifra sobre cada caja es la mediana."
    ),
    x = NULL, y = "sd entre estaciones / sd total de la variable"
  ) +
  tema_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p1b <- ggplot(
  dispersion[sd_z > 0], aes(Variable_lab, sd_z, fill = Variable_lab)
) +
  geom_boxplot(
    outlier.alpha = 0.06, outlier.size = 0.5, width = 0.6, linewidth = 0.35
  ) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  annotation_logticks(sides = "l", linewidth = 0.25) +
  labs(
    title = "(b) La misma dispersión en escala logarítmica",
    subtitle = "Permite ver a la vez las variables casi homogéneas y las muy heterogéneas.",
    x = NULL, y = "sd entre estaciones / sd total (log10)"
  ) +
  tema_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

fig1 <- (p1a / p1b) +
  plot_annotation(
    title = "¿Cuánto difieren entre sí las estaciones meteorológicas de Madrid el mismo día?",
    subtitle = paste0(
      "Variables ordenadas de menor a mayor dispersión mediana · ",
      uniqueN(datos$ESTACION), " estaciones · ", min(ANIOS), "–", max(ANIOS)
    ),
    caption = pie_box,
    theme = tema_base(11)
  )
guardar("01_boxplot_dispersion_entre_estaciones.png", fig1, 11, 11)

# Versión en unidades originales: NO sirve para comparar variables entre sí
# (por eso va con escalas libres y en figura aparte), pero es la que da la
# magnitud física que se cita en el texto: "las estaciones difieren en X °C".
pie_box_crudo <- paste0(
  "Misma construcción que la figura 01, pero en las UNIDADES DE CADA VARIABLE y con escala libre por panel.\n",
  "Sirve para leer la magnitud física de la discrepancia entre estaciones (cuántos °C, cuántos mm, cuántos hPa),\n",
  "NO para comparar variables entre sí: las escalas de los paneles son distintas. Para eso, la figura 01."
)

fig1b <- ggplot(dispersion, aes(Variable_lab, sd_dia, fill = Variable_lab)) +
  geom_boxplot(outlier.alpha = 0.06, outlier.size = 0.5, linewidth = 0.35) +
  facet_wrap(~Variable_lab, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    title = "Dispersión espacial diaria en unidades originales",
    subtitle = paste0(
      "Desviación típica entre estaciones dentro de cada día · ",
      min(ANIOS), "–", max(ANIOS)
    ),
    x = NULL, y = "sd entre estaciones (unidades de la variable)",
    caption = pie_box_crudo
  ) +
  tema_base() +
  theme(axis.text.x = element_blank())
guardar("01b_boxplot_dispersion_unidades_originales.png", fig1b, 11, 7)

# ---- Figura 01c: desfase fijo frente a estructura espacial --------------------
# Es la figura que decide si merece la pena interpolar una covariable. Compara,
# para cada variable, la dispersión entre estaciones de los valores CRUDOS con
# la de las ANOMALÍAS por estación (bloque 6b). Si la caja se hunde al pasar a
# anomalías, las estaciones difieren siempre igual y no hay campo que interpolar.
comparacion_niveles <- rbind(
  dispersion[, .(Variable, sd_z, Version = "Valores crudos")],
  dispersion_anom[, .(Variable, sd_z, Version = "Anomalía por estación")]
)
comparacion_niveles[, Variable_lab := aplicar_orden(Variable)]
comparacion_niveles[, Version := factor(
  Version,
  levels = c("Valores crudos", "Anomalía por estación")
)]

fig1c <- ggplot(
  comparacion_niveles[sd_z > 0], aes(Variable_lab, sd_z, fill = Version)
) +
  geom_boxplot(
    outlier.alpha = 0.04, outlier.size = 0.4, linewidth = 0.35,
    position = position_dodge(width = 0.78), width = 0.7
  ) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  scale_fill_manual(
    values = c(
      "Valores crudos" = "#B2182B",
      "Anomalía por estación" = "#2166AC"
    ),
    name = NULL
  ) +
  labs(
    title = "¿Diferencia entre estaciones o desfase fijo de cada estación?",
    subtitle = paste0(
      "Rojo: dispersión entre estaciones de los valores crudos. ",
      "Azul: la misma dispersión tras restar a cada estación su media histórica."
    ),
    x = NULL, y = "sd entre estaciones / sd total (log10)",
    caption = paste0(
      "Que dos estaciones midan distinto puede ser un DESFASE FIJO (una está 120 m más alta, o resguardada del viento: siempre marca\n",
      "lo mismo de más o de menos) o una DIVERGENCIA que cambia cada día (hoy llueve en el norte y no en el sur). Solo la segunda\n",
      "justifica interpolar la covariable: un desfase constante lo absorbe un intercepto por estación o el propio campo SPDE.\n",
      "Caja azul muy por debajo de la roja = casi todo era desfase fijo. Cajas parecidas = hay estructura espacio-temporal real.\n",
      "Escala log10; se excluyen los días con dispersión nula. Cifras exactas en nivel_vs_anomalia.csv y efecto_altitud.csv."
    )
  ) +
  tema_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "top"
  )
guardar("01c_desfase_fijo_vs_estructura_espacial.png", fig1c, 11, 7)


# ==============================================================================
# 8. BLOQUE 2 — SERIES DIARIAS DE TODAS LAS ESTACIONES SUPERPUESTAS
# ==============================================================================
# Una línea fina y translúcida por estación, sin leyenda (con ~26 estaciones una
# leyenda de colores sería ilegible y no aporta: la pregunta es si las líneas se
# solapan, no cuál es cuál). Encima, en negro y más grueso, la media de ciudad.

media_ciudad <- largo[, .(valor = mean(valor)), by = .(Variable_lab, FECHA)]

panel_series <- function(datos_est, datos_med, titulo, subtitulo, pie,
                         alpha_linea, grosor) {
  ggplot() +
    geom_line(
      data = datos_est, aes(FECHA, valor, group = ESTACION),
      colour = "#2166AC", alpha = alpha_linea, linewidth = grosor, na.rm = TRUE
    ) +
    geom_line(
      data = datos_med, aes(FECHA, valor),
      colour = "black", linewidth = 0.45, na.rm = TRUE
    ) +
    facet_wrap(~Variable_lab, ncol = 2, scales = "free_y") +
    scale_x_date(date_labels = "%b\n%Y", expand = expansion(mult = 0.01)) +
    labs(
      title = titulo, subtitle = subtitulo, x = NULL, y = NULL, caption = pie
    ) +
    tema_base()
}

pie_series <- paste0(
  "Una línea azul translúcida por estación meteorológica (sin leyenda: con ", uniqueN(datos$ESTACION),
  " estaciones no sería legible, y la pregunta no es\n",
  "cuál es cuál sino si se solapan). En negro, la media diaria de la ciudad sobre las estaciones disponibles ese día.\n",
  "Escala diaria, solo observaciones con estado OK; una estación sin dato ese día no aporta línea en ese tramo. Escala libre por panel.\n",
  "LECTURA: líneas indistinguibles de la media negra = variable espacialmente homogénea; abanico de líneas separadas = heterogénea."
)

fig2 <- panel_series(
  largo, media_ciudad,
  "Series diarias superpuestas de todas las estaciones meteorológicas",
  paste0(
    "Madrid, ", min(ANIOS), "–", max(ANIOS),
    " · paneles ordenados de menor a mayor dispersión espacial"
  ),
  pie_series,
  alpha_linea = 0.12, grosor = 0.18
)
guardar("02_series_superpuestas_todas_estaciones.png", fig2, 13, 9)

# Zoom a un año: el periodo completo comprime siete años en el eje x y tapa el
# detalle diario. Un año se lee sin esfuerzo y muestra lo mismo.
largo_zoom <- largo[year(FECHA) == ANIO_ZOOM]
if (nrow(largo_zoom)) {
  fig2b <- panel_series(
    largo_zoom, media_ciudad[year(FECHA) == ANIO_ZOOM],
    sprintf("Series diarias superpuestas — detalle de %d", ANIO_ZOOM),
    "Mismo contenido que la figura 02, restringido a un año para poder ver el detalle diario",
    pie_series,
    alpha_linea = 0.30, grosor = 0.28
  )
  guardar(sprintf("02b_series_superpuestas_%d.png", ANIO_ZOOM), fig2b, 13, 9)
}


# ==============================================================================
# 9. BLOQUE 3 — MAPAS DE UN DÍA REPRESENTATIVO
# ==============================================================================
# Solo puntos coloreados por su valor observado: NO se interpola ninguna
# superficie, porque lo que se quiere ver es el dato crudo y su dispersión, no
# el suavizado que introduciría un kriging o un IDW.

largo_geo <- merge(
  largo, estaciones[, .(ESTACION, X_m, Y_m, Zona)],
  by = "ESTACION"
)

# Elección del día: entre los días en que las SEIS variables tienen al menos
# MIN_ESTACIONES_DIA estaciones, el de mayor número total de observaciones; si
# hay empate (lo habitual), el día central de los empatados. Es reproducible y
# no depende de ninguna fecha elegida a dedo.
elegir_dia_representativo <- function(d) {
  cob <- dcast(
    d[, .N, by = .(FECHA, Variable)], FECHA ~ Variable,
    value.var = "N", fill = 0L
  )
  vars_cob <- setdiff(names(cob), "FECHA")
  cob[, total := rowSums(.SD), .SDcols = vars_cob]
  cob[, minimo := do.call(pmin, .SD), .SDcols = vars_cob]

  candidatos <- cob[minimo >= MIN_ESTACIONES_DIA]
  if (!nrow(candidatos)) {
    cat("  [AVISO] Ningún día cubre las seis variables; se relaja el criterio.\n")
    candidatos <- cob
  }
  candidatos <- candidatos[total == max(total)]
  setorder(candidatos, FECHA)
  list(
    fecha = candidatos$FECHA[ceiling(nrow(candidatos) / 2)],
    n_empates = nrow(candidatos)
  )
}

if (is.null(DIA_MAPA)) {
  eleccion <- elegir_dia_representativo(largo_geo)
  DIA_MAPA <- eleccion$fecha
  cat(sprintf(
    "\nDía representativo elegido automáticamente: %s\n  (día central de los %d días con cobertura máxima en las seis variables)\n",
    format(DIA_MAPA, "%d/%m/%Y"), eleccion$n_empates
  ))
} else {
  cat(sprintf(
    "\nDía representativo fijado a mano: %s\n", format(DIA_MAPA, "%d/%m/%Y")
  ))
}

mapa_variable <- function(dia, variable) {
  d <- largo_geo[FECHA == dia & Variable == variable]
  if (nrow(d) < 1L) {
    return(NULL)
  }

  d_sf <- st_as_sf(d, coords = c("X_m", "Y_m"), crs = 25830)
  rango <- range(d$valor)
  sd_dia <- if (nrow(d) > 1L) sd(d$valor) else NA_real_

  ggplot() +
    geom_sf(
      data = distritos, fill = "grey97", colour = "grey86", linewidth = 0.2
    ) +
    geom_sf(
      data = st_sf(geometry = limite_madrid), fill = NA,
      colour = "grey45", linewidth = 0.45
    ) +
    geom_sf(data = d_sf, aes(colour = valor), size = 3.4, alpha = 0.95) +
    geom_sf_text(
      data = d_sf, aes(label = sprintf(FORMATO[[variable]], valor)),
      size = 2.1, nudge_y = 1100, colour = "grey20"
    ) +
    # Escala de color propia de cada variable: las unidades no son comparables.
    scale_colour_viridis_c(option = OPCION_COLOR[[variable]], name = NULL) +
    labs(
      title = ETIQUETAS[[variable]],
      subtitle = sprintf(
        "%d est. · rango %s–%s · sd = %s",
        nrow(d),
        sprintf(FORMATO[[variable]], rango[1]),
        sprintf(FORMATO[[variable]], rango[2]),
        ifelse(is.na(sd_dia), "-", sprintf("%.2f", sd_dia))
      )
    ) +
    tema_mapa()
}

rejilla_mapas <- function(dia, titulo, subtitulo, pie, nombre_archivo) {
  paneles <- Filter(
    Negate(is.null), lapply(orden_variables, function(v) mapa_variable(dia, v))
  )
  if (!length(paneles)) {
    cat(sprintf(
      "  [AVISO] Sin datos para el %s; no se genera %s\n",
      format(dia, "%d/%m/%Y"), nombre_archivo
    ))
    return(invisible(NULL))
  }
  fig <- wrap_plots(paneles, ncol = 3) +
    plot_annotation(
      title = titulo, subtitle = subtitulo, caption = pie,
      theme = tema_base(11)
    )
  guardar(nombre_archivo, fig, 14, 9)
}

pie_mapas <- paste0(
  "Valor DIARIO observado en cada estación el día indicado, sobre el límite del municipio y los distritos de Madrid.\n",
  "No se interpola ninguna superficie: son los puntos de medida en crudo, coloreados por su valor y con la cifra al lado.\n",
  "Cada panel lleva su propia escala de color porque las seis variables tienen unidades distintas; los paneles NO son\n",
  "comparables entre sí en color, solo en uniformidad. Día elegido automáticamente como el de mayor cobertura conjunta.\n",
  "LECTURA: color uniforme entre puntos = campo homogéneo ese día; gradiente centro-periferia o norte-sur = estructura espacial."
)

rejilla_mapas(
  DIA_MAPA,
  sprintf("Valores por estación el %s", format(DIA_MAPA, "%d de %B de %Y")),
  "Un mapa por variable, ordenados de menor a mayor dispersión espacial media · puntos crudos, sin interpolación",
  pie_mapas,
  sprintf("03_mapas_dia_%s.png", format(DIA_MAPA, "%Y%m%d"))
)

# Complemento: en un día de cobertura máxima la precipitación suele ser 0 en
# toda la ciudad, así que su mapa sale plano y no dice nada. Se añade el día
# más lluvioso (más estaciones con precipitación > 0), donde esa variable sí se
# puede leer. Va en figura aparte y etiquetada como complementaria.
lluvia <- largo_geo[
  Variable == "Precipitaciones",
  .(n_est = .N, n_con = sum(valor > 0), total = sum(valor)),
  by = FECHA
][n_est >= MIN_ESTACIONES_DIA & n_con > 0]

if (nrow(lluvia)) {
  setorder(lluvia, -n_con, -total)
  DIA_LLUVIA <- lluvia$FECHA[1]
  cat(sprintf(
    "Día lluvioso complementario: %s (%d de %d estaciones con precipitación)\n",
    format(DIA_LLUVIA, "%d/%m/%Y"), lluvia$n_con[1], lluvia$n_est[1]
  ))
  rejilla_mapas(
    DIA_LLUVIA,
    sprintf(
      "Día lluvioso complementario: %s", format(DIA_LLUVIA, "%d de %B de %Y")
    ),
    "En el día de cobertura máxima la precipitación es casi siempre nula en toda la ciudad; este día permite leerla",
    paste0(
      pie_mapas,
      "\nDía elegido como aquel con precipitación > 0 en el mayor número de estaciones (desempate: mayor acumulado)."
    ),
    sprintf("03b_mapas_dia_lluvioso_%s.png", format(DIA_LLUVIA, "%Y%m%d"))
  )
}


# ==============================================================================
# 10. BLOQUE 4 — COMPARACIÓN POR ZONA GEOGRÁFICA
# ==============================================================================
# Clasificación GEOGRÁFICA de las estaciones de clima (centro / corona /
# periferia abierta). No se usa aquí la tipología tráfico-fondo-suburbana: esa
# describe la exposición de las estaciones de NO2 y no aplica a la red
# meteorológica.

# Mapa de la agrupación, para que se vea qué estación está en cada zona.
est_sf <- st_as_sf(estaciones, coords = c("X_m", "Y_m"), crs = 25830)

fig4a <- ggplot() +
  geom_sf(
    data = distritos, fill = "grey97", colour = "grey86", linewidth = 0.2
  ) +
  geom_sf(
    data = st_sf(geometry = limite_madrid), fill = NA,
    colour = "grey45", linewidth = 0.5
  ) +
  geom_sf(data = est_sf, aes(colour = Zona, shape = asignacion), size = 3.2) +
  geom_sf_text(
    data = est_sf, aes(label = ESTACION),
    size = 2.1, nudge_y = 1200, colour = "grey25"
  ) +
  scale_colour_manual(values = COLOR_ZONA, name = "Zona") +
  scale_shape_manual(
    values = c("manual" = 16, "automatica" = 17), name = "Asignación"
  ) +
  labs(
    title = "Agrupación geográfica de las estaciones meteorológicas",
    subtitle = "Clasificación geográfica propia de la red de clima; no es la tipología tráfico/fondo/suburbana del NO2",
    caption = paste0(
      "'Periferia abierta' agrupa las estaciones en entorno no urbanizado o de borde (El Pardo, Casa de Campo, Juan Carlos I,\n",
      "Ensanche de Vallecas): el criterio no es solo la distancia al centro, sino el entorno. Casa de Campo está a 3,7 km de Sol\n",
      "y aun así entra ahí por ser parque forestal. Los puntos triangulares, si los hay, son estaciones asignadas automáticamente\n",
      "por anillos de distancia a Puerta del Sol (", RADIO_CENTRO_KM, " km y ", RADIO_PERIFERIA_KM,
      " km) por no figurar en la lista manual del script."
    )
  ) +
  tema_mapa(11) +
  theme(plot.caption = element_text(
    size = 7, hjust = 0, colour = "grey35", lineheight = 1.2
  ))
guardar("04a_mapa_zonas_estaciones.png", fig4a, 9, 9)

# Distribución de valores por zona, una faceta por variable.
pie_zonas <- paste0(
  "Distribución de TODOS los valores diarios ", min(ANIOS), "–", max(ANIOS),
  " agrupados por zona geográfica de la estación.\n",
  "Escala libre por panel (unidades distintas). Solo observaciones con estado OK; una estación sin dato un día no participa.\n",
  "LECTURA: cajas desplazadas entre zonas = la variable no es la misma en el centro que en la periferia (isla de calor,\n",
  "exposición al viento); cajas alineadas = el nivel no depende de la zona, aunque pueda haber ruido local."
)

fig4b <- ggplot(largo_geo, aes(Zona, valor, fill = Zona)) +
  geom_boxplot(outlier.alpha = 0.03, outlier.size = 0.3, linewidth = 0.35) +
  facet_wrap(~Variable_lab, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = COLOR_ZONA, guide = "none") +
  labs(
    title = "Distribución de cada variable climática por zona geográfica",
    subtitle = paste0(
      "Madrid, ", min(ANIOS), "–", max(ANIOS),
      " · escala diaria · paneles ordenados por dispersión espacial"
    ),
    x = NULL, y = NULL, caption = pie_zonas
  ) +
  tema_base() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 8))
guardar("04b_boxplot_por_zona.png", fig4b, 12, 8)

# Serie media por zona. Se dibuja el año de zoom en diario (legible) y todo el
# periodo en media mensual (un desfase sistemático se ve mucho mejor así que
# con siete años de líneas diarias solapadas).
media_zona_dia <- largo_geo[
  year(FECHA) == ANIO_ZOOM,
  .(valor = mean(valor)),
  by = .(Variable_lab, Zona, FECHA)
]
media_zona_mes <- largo_geo[
  , .(valor = mean(valor)),
  by = .(Variable_lab, Zona, MES = as.Date(format(FECHA, "%Y-%m-01")))
]

fig4c <- ggplot(media_zona_dia, aes(FECHA, valor, colour = Zona)) +
  geom_line(linewidth = 0.35, alpha = 0.85, na.rm = TRUE) +
  facet_wrap(~Variable_lab, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = COLOR_ZONA, name = NULL) +
  scale_x_date(date_labels = "%b", date_breaks = "2 months") +
  labs(
    title = sprintf("Serie diaria media por zona geográfica — %d", ANIO_ZOOM),
    subtitle = "Media de las estaciones de cada zona en cada día",
    x = NULL, y = NULL,
    caption = paste0(
      "Cada línea es la media diaria de las estaciones de esa zona (las que tienen dato ese día). Escala libre por panel.\n",
      "LECTURA: líneas superpuestas = la zona no cambia el valor; separación mantenida en el tiempo = diferencia sistemática\n",
      "entre centro y periferia, no ruido de un día suelto."
    )
  ) +
  tema_base() +
  theme(legend.position = "bottom")
guardar(sprintf("04c_series_medias_por_zona_%d.png", ANIO_ZOOM), fig4c, 12, 9)

fig4d <- ggplot(media_zona_mes, aes(MES, valor, colour = Zona)) +
  geom_line(linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = 0.6, na.rm = TRUE) +
  facet_wrap(~Variable_lab, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = COLOR_ZONA, name = NULL) +
  scale_x_date(date_labels = "%Y") +
  labs(
    title = sprintf(
      "Media mensual por zona geográfica — %d–%d", min(ANIOS), max(ANIOS)
    ),
    subtitle = "Agregado mensual de las medias por zona: aísla el desfase sistemático del ruido diario",
    x = NULL, y = NULL,
    caption = paste0(
      "Media mensual de los valores diarios de las estaciones de cada zona. Escala libre por panel.\n",
      "Complementa a la figura 04c: si el desfase entre zonas persiste mes a mes durante siete años, es estructural."
    )
  ) +
  tema_base() +
  theme(legend.position = "bottom")
guardar("04d_series_mensuales_por_zona.png", fig4d, 12, 9)

# Cuantificación: diferencia diaria pareada Periferia abierta - Centro urbano.
# Es la forma directa de contrastar la hipótesis de isla de calor: se comparan
# ambas zonas EL MISMO DÍA, con lo que el ciclo estacional se cancela.
resumen_zona <- largo_geo[, .(
  N       = .N,
  media   = round(mean(valor), 3),
  sd      = round(sd(valor), 3),
  mediana = round(median(valor), 3)
), by = .(Variable, Zona)]
setorder(resumen_zona, Variable, Zona)
fwrite(resumen_zona, file.path(DIR_TAB, "resumen_por_zona.csv"))

media_zona_amplia <- dcast(
  largo_geo[, .(valor = mean(valor)), by = .(Variable, Zona, FECHA)],
  Variable + FECHA ~ Zona,
  value.var = "valor"
)

if (all(c("Periferia abierta", "Centro urbano") %in% names(media_zona_amplia))) {
  media_zona_amplia[, dif := `Periferia abierta` - `Centro urbano`]
  diferencias <- media_zona_amplia[is.finite(dif), .(
    N_dias            = .N,
    dif_media         = round(mean(dif), 3),
    dif_mediana       = round(median(dif), 3),
    dif_p2.5          = round(quantile(dif, 0.025), 3),
    dif_p97.5         = round(quantile(dif, 0.975), 3),
    pct_dias_positiva = round(100 * mean(dif > 0), 1)
  ), by = Variable]
  setorder(diferencias, -N_dias)
  fwrite(diferencias, file.path(DIR_TAB, "diferencia_periferia_centro.csv"))

  cat("\n--- Diferencia diaria pareada: periferia abierta - centro urbano ---\n")
  cat("    (mismo día en ambas zonas, con lo que el ciclo estacional se cancela)\n")
  print(diferencias)
} else {
  cat("\n  [AVISO] Faltan zonas para la diferencia pareada centro/periferia.\n")
}


# ==============================================================================
# 11. CONCLUSIÓN POR VARIABLE
# ==============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("LECTURA VARIABLE A VARIABLE\n")
cat(strrep("=", 78), "\n", sep = "")
cat(sprintf(
  "Criterio declarado: ratio = mediana(sd entre estaciones del día) / sd total.\n  ratio < %.2f homogénea | %.2f-%.2f intermedia | > %.2f heterogénea\n  En las variables con más de un 25%% de días de dispersión nula se usa la\n  mediana sobre los días activos (relevante en precipitación).\n\n",
  UMBRAL_HOMOGENEA, UMBRAL_HOMOGENEA, UMBRAL_HETEROGENEA, UMBRAL_HETEROGENEA
))

final <- merge(
  resumen[, .(Variable, ratio_criterio, sd_espacial_med, Est_por_dia, Lectura)],
  descomposicion[, .(Variable, pct_desfase_fijo)],
  by = "Variable", all.x = TRUE
)
setorder(final, ratio_criterio)

for (i in seq_len(nrow(final))) {
  cat(sprintf(
    "  %-28s ratio = %.3f  (sd espacial mediana = %.2f · %.0f est./día)  -> %s\n",
    ETIQUETAS[[as.character(final$Variable[i])]],
    final$ratio_criterio[i], final$sd_espacial_med[i], final$Est_por_dia[i],
    final$Lectura[i]
  ))
  # Aviso cuando la dispersión resulta ser sobre todo un desfase fijo: la
  # variable parece heterogénea pero no hay campo espacial que interpolar.
  if (!is.na(final$pct_desfase_fijo[i]) && final$pct_desfase_fijo[i] > 60) {
    cat(sprintf(
      "  %-28s   OJO: el %.0f%% de esa dispersión es desfase fijo de estación;\n  %-28s   basta un intercepto por estación, no hace falta interpolar.\n",
      "", final$pct_desfase_fijo[i], ""
    ))
  }
}

cat(sprintf("\nFiguras: %s\nTablas : %s\n\n", DIR_FIG, DIR_TAB))
