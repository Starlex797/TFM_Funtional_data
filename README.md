# TFM Functional Data

Proyecto de TFM para el análisis espacial y espacio-temporal de contaminación,
clima y tráfico en Madrid.

## Estructura

```text
TFM_Funtional_data/
├── data/
│   ├── raw/                 # Datos originales; no se modifican manualmente
│   ├── interim/             # Datos temporales entre etapas
│   └── processed/           # Datos limpios y conjuntos maestros
├── R/
│   ├── cleaning/            # Funciones de limpieza y unión
│   ├── interpolation/       # Funciones de interpolación
│   ├── spatial/             # Utilidades cartográficas y espaciales
│   └── utilities/           # Diccionarios y utilidades generales
├── scripts/
│   ├── 01_preprocessing/    # Preparación de datos
│   ├── 02_eda/              # Análisis exploratorio y visualizaciones
│   ├── 03_interpolation/    # Kriging e interpolación
│   ├── 04_modeling/         # Modelos INLA-SPDE
│   ├── 05_validation/       # Calidad de datos y validación
│   └── 06_simulation/       # Simulaciones
├── outputs/
│   ├── analysis/            # Resultados de análisis exploratorios
│   ├── figures/             # Figuras principales
│   ├── tables/              # Tablas y métricas
│   ├── models/              # Modelos exportados
│   ├── reports/             # Informes generados por los scripts
│   ├── logs/                # Registros de ejecución
│   └── archive/             # Resultados y sesiones históricas
├── report/
│   ├── source/              # Rmd, Markdown y LaTeX
│   └── rendered/            # PDF, HTML y DOCX
└── docs/                    # Metodología, esquemas y notas del proyecto
```

## Convenciones

- Las rutas del código se construyen desde la raíz con `here::here()`.
- `data/raw/` conserva los datos originales y `data/processed/` contiene los
  resultados reproducibles del preprocesamiento.
- `R/` contiene funciones reutilizables; `scripts/` contiene ejecuciones de
  análisis completas.
- Los resultados nuevos se guardan bajo `outputs/`. El contenido de
  `outputs/archive/` se conserva como histórico y no debe sobrescribirse.
- Los archivos fuente de los informes se editan en `report/source/`; sus
  versiones compiladas se guardan en `report/rendered/`.
