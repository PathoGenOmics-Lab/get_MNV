# Conteo de lecturas contra bcftools mpileup

[English](README.md) · **Español**

get_MNV tiene tres suites que comprueban su anotación y una que comprueba su
fasing, y hasta esta toda comprobación sobre el *conteo de lecturas* estaba
escrita por la misma mano que el código: las expectativas de los escenarios y la
aritmética de la suite de fase. Las dos cazan errores de implementación. Ninguna
puede cazar un malentendido compartido sobre lo que vale una lectura, porque las
dos parten del mismo entendimiento.

Esta suite le da el mismo BAM a `bcftools mpileup`, que no lo ha escrito nadie
de aquí, y compara cuatro números por sitio:

| get_MNV | mpileup |
|---|---|
| `Total Reads` | `FORMAT/DP` |
| `SNP Reads` | `FORMAT/AD` de ese alelo |
| `SNP Forward Reads` | `FORMAT/ADF` |
| `SNP Reverse Reads` | `FORMAT/ADR` |

```bash
cargo build -p get_mnv --release
python3 tests/pileup/run.py
```

## Por qué esos cuatro

Son donde han estado de verdad los defectos de esta capa. El bug de atribución
de hebra arreglado en esta rama, donde un mate a cientos de bases de una
variante le prestaba su hebra, es una discrepancia directa con `ADF`/`ADR`, y se
habría cazado el mismo día en que se escribió en vez de sobrevivir hasta que a
alguien se le ocurrió mirar.

## Mantener comparables las dos herramientas

Una diferencia que viene de hacerles preguntas distintas no enseña nada, así que
los umbrales se le dan a las dos explícitamente, y BAQ se desactiva con `-B`
porque get_MNV no recalibra calidades de base. Todos los sitios de comparación
están en el relleno entre genes, así que cada uno es una fila de una sola
posición y ninguna agrupación por codones se interpone entre los dos recuentos.

Los casos cubren las dos hebras, un alelo que solo aparece en una, pares cuyos
dos mates leen sitios distintos, mates solapados, bases por debajo del umbral de
calidad, lecturas por debajo del umbral de mapeo, y varios sitios contados desde
una misma caché.

## Diferencias aceptadas

Como en la [suite diferencial](../differential/README.es.md), la regla no es que
las herramientas tengan que coincidir. Una diferencia entendida se escribe en el
caso que la produce, con su razón; una diferencia que nadie ha explicado hace
fallar la ejecución.

Hay una, en mates solapados. Las dos herramientas cuentan 16 moléculas y
coinciden en profundidad y soporte; reparten las hebras distinto. bcftools
elimina el solape quedándose con la base de un mate y poniendo la del otro a
cero, así que cada fragmento cae en exactamente un brazo y `ADF + ADR == AD`. En
cuál cae es en la práctica arbitrario: los dieciséis fragmentos son idénticos y
los divide en 5 forward y 11 reverse, un reparto que no lleva ninguna
información sobre hebra. get_MNV acredita los dos brazos, que es lo que de
verdad ocurre con esas moléculas, y significa que ninguna es evidencia de una
sola hebra. La convención difiere; ninguno de los dos recuentos está mal.

Con `--count-mates-separately` la discrepancia es total, porque bcftools no
tiene un modo equivalente. Ese caso está para fijar que el modo *por defecto* es
el que coincide con una herramienta de fuera.

## Comprobar que todavía puede fallar

Contra la build anterior al arreglo de hebra (`18e426c`):

```bash
GET_MNV=/ruta/a/get_mnv/anterior python3 tests/pileup/run.py
```

reporta el bug directamente:

```
03_mates_far_apart: reverse at 320 is 15 for get_MNV and 0 for mpileup
03_mates_far_apart: forward at 750 is 15 for get_MNV and 0 for mpileup
```

Quince moléculas reverse en una posición a la que ninguna lectura reverse llegó.

## Lo que no cubre

El conteo de indels y eventos complejos. La representación de indels de mpileup
es la suya (`IDV`, `IMF`, y los ALT candidatos que elige por su cuenta), así que
alinearla con los recuentos exactos por CIGAR de get_MNV compararía dos
definiciones distintas en vez de una implementación contra otra. La profundidad
y la hebra de las sustituciones son la parte donde las dos herramientas
responden de verdad a la misma pregunta.
