# ==============================================================================
# GENERA UN HTML CON LAS TABLAS DE CALIDAD PARA COPIAR/PEGAR EN WORD
# ==============================================================================
# Abre el HTML resultante en el navegador, selecciona una tabla, Ctrl+C y pega
# en Word: se pega como tabla nativa (con bordes). Formato pensado para Word.
#
# Fuente: outputs/tables/calidad_datos/dias_ausentes_por_estacion/*.csv
# Salida: outputs/tables/calidad_datos/TABLAS_CALIDAD_WORD.html
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

DIR_IN  <- here("outputs", "tables", "calidad_datos", "dias_ausentes_por_estacion")
OUT     <- here("outputs", "tables", "calidad_datos", "TABLAS_CALIDAD_WORD.html")
ANIOS   <- 2019:2025

# Etiquetas legibles de columnas
ETI <- c("Temperatura"="Temperatura", "Humedad_Relativa"="Humedad relativa",
         "Precipitaciones"="Precipitaciones", "Presion Barométrica"="Presión barométrica",
         "Radiación Solar"="Radiación solar", "Velocidad Viento"="Velocidad viento",
         "NO2"="NO2", "ESTACION"="Estación", "Variable"="Variable", "total"="Total")

# ---- Helper: convierte un data.table en una <table> HTML lista para Word ------
tabla_html <- function(dt, titulo, nota = NULL, num_cols = NULL) {
  dt <- as.data.table(copy(dt))
  # NA -> "—"; decimales con coma (formato español)
  for (cc in names(dt)) {
    if (is.numeric(dt[[cc]])) {
      dt[[cc]] <- ifelse(is.na(dt[[cc]]), "—",
                         format(dt[[cc]], decimal.mark = ",", trim = TRUE))
    } else {
      dt[[cc]] <- ifelse(is.na(dt[[cc]]), "—", as.character(dt[[cc]]))
    }
  }
  encabezados <- sapply(names(dt), function(x) if (!is.na(ETI[x])) ETI[x] else x)

  th <- paste0("<th style='border:1px solid #888;padding:5px 9px;background:#20548f;",
               "color:#fff;text-align:center;'>", encabezados, "</th>", collapse = "")

  filas <- apply(dt, 1, function(fila) {
    celdas <- sapply(seq_along(fila), function(j) {
      val <- fila[j]
      alin <- if (j == 1) "left" else "center"
      paste0("<td style='border:1px solid #bbb;padding:4px 9px;text-align:", alin, ";'>",
             val, "</td>")
    })
    paste0("<tr>", paste(celdas, collapse = ""), "</tr>")
  })

  paste0(
    "<h3 style='font-family:Arial;color:#20548f;margin-bottom:2px;'>", titulo, "</h3>",
    if (!is.null(nota)) paste0("<p style='font-family:Arial;font-size:12px;color:#555;margin-top:0;'>", nota, "</p>") else "",
    "<table style='border-collapse:collapse;font-family:Arial;font-size:13px;margin-bottom:22px;'>",
    "<tr>", th, "</tr>", paste(filas, collapse = ""), "</table>"
  )
}

partes <- c(paste0(
  "<html><head><meta charset='utf-8'></head><body style='margin:24px;'>",
  "<h1 style='font-family:Arial;color:#20548f;'>Calidad de los datos — días ausentes por estación</h1>",
  "<p style='font-family:Arial;font-size:12px;color:#555;'>",
  "Un «día ausente» = una estación que, un día concreto, no tiene ninguna de sus 24 horas ",
  "(estación offline). «—» = la estación no mide esa variable. Para pegar en Word: abre este ",
  "archivo en el navegador, selecciona la tabla, Ctrl+C y pega en el documento.</p><hr>"
))

# ---- (A) Tablas anuales -------------------------------------------------------
media <- fread(file.path(DIR_IN, "matriz_media_por_estacion.csv"), header = TRUE)
total <- fread(file.path(DIR_IN, "matriz_total_dias_ausentes.csv"), header = TRUE)
# Ordenar filas de variables de forma logica
orden <- c("Temperatura","Humedad_Relativa","Precipitaciones",
           "Presion Barométrica","Radiación Solar","Velocidad Viento","NO2")
media <- media[match(orden, Variable)]
total <- total[match(orden, Variable)]

partes <- c(partes,
  tabla_html(media,
             "Tabla A1. Media de días ausentes POR ESTACIÓN (de 365/366 días)",
             "Promedio de días/año que estuvo caída una estación que mide esa variable."),
  tabla_html(total,
             "Tabla A2. Total de días-estación ausentes (suma de todas las estaciones)",
             "Suma sobre todas las estaciones de la red para cada variable y año.")
)

# ---- (B) Una tabla por año (estación x variable) ------------------------------
for (a in ANIOS) {
  f <- file.path(DIR_IN, sprintf("dias_ausentes_por_estacion_%d.csv", a))
  if (!file.exists(f)) next
  w <- fread(f)
  # Solo estaciones con algún dato ese año (total > 0 o alguna celda no NA)
  vcols <- setdiff(names(w), c("ESTACION","total"))
  w <- w[rowSums(!is.na(as.matrix(w[, ..vcols]))) > 0]
  setorder(w, -total)
  partes <- c(partes,
    tabla_html(w,
               sprintf("Tabla B%d. Días ausentes por estación en %d (de %d días)",
                       a - 2018, a, if (a %in% c(2020,2024)) 366 else 365),
               "Ordenado de más a menos días ausentes."))
}

partes <- c(partes, "</body></html>")
writeLines(paste(partes, collapse = "\n"), OUT)

cat("HTML generado en:\n  ", OUT, "\n")
cat("Ábrelo con doble clic (se abre en el navegador) y copia/pega las tablas en Word.\n")
