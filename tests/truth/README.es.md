# Validación con verdad sintética

[English](README.md) · **Español**

La [suite diferencial](../differential/README.es.md) compara get_MNV con
`bcftools csq` y SnpEff. Esa es la comprobación más fuerte, porque esas
herramientas no son nuestras. Pero no puede validar justo aquello para lo que
existe get_MNV: **los MNV a nivel de codón**. Ningún anotador por variante
combina dos SNV de un codón como lo hace get_MNV, así que precisamente en esas
llamadas los oráculos discrepan por diseño y no pueden decir si el residuo es
correcto.

Esta suite cubre ese hueco, con un conjunto de datos cuya respuesta conocemos
porque lo construimos nosotros.

```bash
cargo build -p get_mnv --release
python3 tests/truth/run.py
```

## Cómo se deriva la verdad

La expectativa no puede salir de la lógica de get_MNV, o el test solo le da la
razón al código. Así que toma otra ruta hacia la misma respuesta:

```
mutar el genoma  ->  volver a extraer el CDS  ->  traducir  ->  comparar proteínas
```

get_MNV hace cirugía a nivel de codón, in situ. Esto vuelve a traducir el CDS
entero desde una tabla de codones escrita en `run.py`, con operaciones de cadena
sobre la referencia. Un desliz de implementación aparece como discrepancia.

## Qué cubre

Cuatro genes, generados de forma determinista: hebra más y hebra menos, de un
solo exón y con splicing. Para cada codón de cada uno:

- cada sustitución simple en cada posición (3 posiciones x 3 alternativos)
- cada par de sustituciones en el mismo codón
- el triple

Unos 1560 casos, cada uno como su propio VCF para que no se interfieran. De cada
uno se comprueban gen, cambio de aminoácido, tipo de cambio, codón de
referencia, codón alternativo y el descriptor `HGVS c.`.

Más 2130 indels: deleciones e inserciones de 1, 2 y 3 bases en **cada posición
de cada exón**. Dónde va una inserción es la decisión más sutil de todo el
camino de indels, así que no se muestrea: las bases deben ir en complemento
inverso en la hebra menos, deben caer entre los dos residuos correctos, y el
hueco por encima de la coordenada más alta del gen queda fuera del CDS en
*ambas* hebras, porque ahí es el extremo 3' de un transcrito de hebra más y el
5' de uno de hebra menos. Con anclas repartidas esos bordes quedarían sin probar. Su expectativa sale de la misma retraducción, aplicando el indel en el
espacio de coordenadas del CDS y no recortando el genoma, porque las bases de
alrededor conservan sus coordenadas. Se comprueban tipo de cambio, clase de
evento y ambos codones; la notación HGVS de proteína para indels es formato
propio de get_MNV, así que reproducirla aquí probaría ortografía, no biología.

La ejecución completa tarda unos seis segundos.

## Cuánto vale

Menos que un oráculo de terceros, y conviene decir por qué: es una segunda
implementación de la misma mano. Caza errores de implementación, no un
malentendido compartido sobre la biología. Las *etiquetas* de consecuencia
(`Synonymous`, `Stop gained`, ...) siguen la convención documentada de get_MNV,
así que en eso valida el código contra la documentación y no contra la
naturaleza. El aminoácido, los codones y la coordenada del CDS son la parte
independiente.

Y es la razón de que la suite diferencial siga existiendo a su lado. Ninguna
sustituye a la otra.

## Comprobar que todavía puede fallar

Un test que no puede fallar no vale nada. `GET_MNV_BIN` apunta la suite a otro
binario:

```bash
GET_MNV_BIN=/ruta/a/get_mnv/anterior python3 tests/truth/run.py
```

Contra el binario anterior al arreglo de orientación de codón reporta 240
discrepancias, todas en los genes de hebra menos. Contra el anterior al arreglo
del codón de inicio reporta 11 más, todas indels que eliminan el iniciador.

Los fallos se agrupan por categoría y se muestran ejemplos de cada grupo, en vez
de los primeros en orden, porque si no una categoría ruidosa entierra a otra: la
ejecución de arriba enseña tanto los 240 de orden como los 11 del codón de inicio.

## Una nota sobre escribir el oráculo

Su primera ejecución reportó 665 discrepancias. Todas eran un bug de este
fichero, no de get_MNV: una variable de bucle filtrada hacía que el VCF generado
llevara un alelo alternativo distinto de aquel para el que se había calculado la
expectativa. El oráculo hay que comprobarlo tanto como a lo que comprueba.
