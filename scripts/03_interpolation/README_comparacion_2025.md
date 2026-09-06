# Sensibilidad de RMSE a k y eleccion manual

El analisis utiliza los archivos de 2025 con sufijo `5`, en las escalas
horaria, diaria y mensual. No elige automaticamente un k ni un metodo.

## Como usarlo

1. Ejecutar desde la raiz del proyecto:

   ```r
   source("scripts/03_interpolation/comparacion_metodos_interpolacion.R")
   ```

2. Abrir los graficos `RMSE_k_<variable>.png` en cada carpeta de escala.
   El eje horizontal es el numero de vecinos; el vertical, RMSE en las unidades
   de la variable. Cada figura muestra KNN, IDW p=1 e IDW p=2. Las lineas
   horizontales son Media y 1-NN. Un RMSE menor indica menor error.

3. Escribir los k elegidos en `K_ELEGIDO`, dentro de
   `scripts/03_interpolation/configuracion_comparacion_clima.R`.
   Hay una fila por variable y una columna por escala. El mismo k elegido
   se aplica a KNN y a ambas potencias IDW para esa variable y escala.

4. Actualizar las tablas y las marcas de k en las figuras, sin repetir las
   interpolaciones:

   ```r
   source("scripts/05_validation/tabla_comparacion_escalas.R")
   ```

Se genera una tabla **separada por escala**, en CSV, PNG y PDF. Si todavia no
se ha indicado un k, su fila aparece pendiente y los RMSE locales quedan vacios.
La tabla `sensibilidad_rmse_<escala>.csv` contiene todos los k representados,
para consultarlos antes de decidir. No se genera una tabla conjunta ni rankings.

## Calculo que se explica al profesor

Para cada ubicacion y cada instante:

1. Retirar la observacion objetivo y todos los sensores colocalizados.
2. Ordenar las otras ubicaciones por distancia.
3. Conservar los vecinos que tienen observaciones validas en ese instante.
4. Predecir usando la media de todos, el vecino mas cercano, la media de los k
   vecinos (KNN), o una media ponderada por distancia (IDW).
5. Calcular el error respecto a la observacion retirada.

En IDW los pesos son proporcionales a `1/distancia^p`, con p=1 o p=2.
Primero se calcula el error cuadratico medio (MSE) de cada ubicacion. El RMSE
mostrado es la raiz de la media de esos MSE: asi las ubicaciones tienen el
mismo peso aunque no todas tengan el mismo numero de observaciones.

Todos los puntos de una curva se calculan sobre **los mismos objetivos**.
El rango de k se limita para conservar al menos el 90% de los casos evaluables,
globalmente y en cada ubicacion. El umbral es visible en la configuracion.
Esto establece que valores pueden compararse con suficiente cobertura, pero
no escoge el k con menor error. Los CSV `cobertura_<variable>.csv` muestran
la disponibilidad incluso para los k que no se pueden incluir en la curva.

No se reduce k cuando faltan vecinos. Con k=1, KNN e IDW coinciden con 1-NN.
Si solo se puede representar k=1, los datos de esa muestra no permiten
comparar varios numeros de vecinos.

## Precipitacion

La curva y la tabla de precipitacion evaluan **solo objetivos con lluvia
observada mayor de 0.1 mm**. Los ceros de las estaciones vecinas se mantienen:
el interpolador debe reconstruir la lluvia aunque algunos vecinos esten secos.

Se anade la referencia `Siempre cero`. Su RMSE durante la lluvia permite
comprobar si interpolar mejora respecto a predecir cero en todos esos casos.
Esa referencia no es un metodo candidato.

Esta es una evaluacion condicionada a que llueva en el objetivo. No mide la
deteccion de lluvia en periodos secos. El criterio se aplica al valor de la
escala correspondiente: mm horarios, diarios o mensuales.

## Calidad de los datos

- Se usan exclusivamente los tres archivos `5`, sin modificarlos.
- Solo se validan valores finitos con estado `OK`.
- En diario y mensual se excluyen agregados que contengan alguna hora imputada,
  aunque el estado agregado sea `OK`. El fichero horario permite comprobarlo.
- Los agregados aceptados no se recalculan y pueden contener periodos incompletos.
- Las coordenadas con desplazamiento decimal se corrigen y contrastan con
  `X_km/Y_km`. Las coordenadas originales se guardan.
- Los sensores exactamente colocalizados se agrupan por ubicacion y se usa su
  media observada. La validacion retira la ubicacion completa.

En mensual solo hay doce periodos, y excluir agregados con horas imputadas
puede reducir mucho el rango de k. Comparar las escalas tambien implica comparar
redes efectivas distintas. Los acumulados de lluvia de distinta duracion no
tienen errores directamente comparables.

La preparacion esta aislada en `R/interpolation/preparacion_comparacion_clima.R`.
El calculo de vecinos, errores y graficos esta en
`R/interpolation/comparacion_clima_multiescala.R`, dividido en cinco pasos.
Las tablas manuales estan en `R/interpolation/resumen_comparacion_clima.R`.

Los scripts `simulacion_interpolacion_clima.R` y
`simulacion_interpolacion_clima_horaria.R` ejecutan el mismo procedimiento
solamente para diario y horario, respectivamente.

Las salidas se guardan bajo `outputs/interpolacion/rmse_k_manual_2025_5_<fecha>/`.
Las ejecuciones anteriores se conservan en sus carpetas. La comparacion actual
no contiene ensamble, seleccion automatica ni validacion anidada. Elegir k
mirando estas curvas es un analisis exploratorio: ese mismo RMSE no constituye
una evaluacion independiente del k elegido.
