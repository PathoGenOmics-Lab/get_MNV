# Ligamiento: distinguir un haplotipo de una casualidad

Dos variantes en un mismo codón solo son una variante multinucleotídica si las
mismas moléculas llevan ambas. get_MNV reporta con cuánta fuerza van juntas,
como `D'` con su valor p, en las columnas `Haplotype LD` y `Haplotype LD p`
(`LD` / `LDP` en la salida VCF). Requiere `--bam`.

## Por qué un ratio de co-ocurrencia no basta

La medida obvia es la fracción de moléculas que llevan una variante y llevan
también la otra. Tiene un fallo que empeora cuanto más comunes son las variantes.

Coge dos sustituciones que están cada una en el 90% de las moléculas, sin nada
que las relacione. Por pura aritmética se encuentran juntas en el 81% de ellas, y
el ratio reporta `0.9000`. Leído como ligamiento, eso dice "estas dos casi
siempre van juntas". Lo que dice en realidad es "estas dos son comunes".

`MNV Phasing Support` es ese ratio y se sigue reportando, porque "con qué
frecuencia co-ocurren" es una pregunta legítima con una respuesta útil. Pero no
puede responder a la que hay debajo, que es si esa co-ocurrencia significa algo.

## Qué mide D'

`D'` es el exceso de co-ocurrencia sobre lo que predicen las dos frecuencias,
normalizado por el mayor exceso que esas frecuencias podrían haber producido.
Esa normalización es lo que lo hace comparable entre sitios, y deja el valor en
`[-1, 1]`.

| `D'` | Lectura |
|---|---|
| `+1` | Las variantes van juntas todo lo que sus frecuencias permiten. Un haplotipo. |
| `~0` | Co-ocurren exactamente lo que predice el azar. Dos variantes que comparten codón, que no es lo mismo que un MNV. |
| `-1` | Se excluyen: ambas presentes, nunca en la misma molécula. |

Los mismos datos de arriba, puntuados como toca:

| Moléculas con ambas | una | la otra | ninguna | Ratio | `D'` | p |
|---|---|---|---|---|---|---|
| 81 | 9 | 9 | 1 | 0.9000 | **0.0000** | 1.0 |
| 90 | 0 | 0 | 10 | 1.0000 | **1.0000** | 5.8e-14 |
| 3 | 0 | 0 | 97 | 1.0000 | **1.0000** | 6.2e-6 |
| 0 | 20 | 20 | 0 | 0.0000 | **-1.0000** | 1.5e-11 |
| 10 | 10 | 10 | 10 | 0.5000 | **0.0000** | 1.0 |

Las filas primera y quinta son la cuestión. Un ratio de 0.9 y uno de 0.5 parecen
ligamiento fuerte y moderado; los dos son exactamente lo que predice la
independencia.

## El lado negativo

Un `D'` cerca de `-1` es un hallazgo, no un fallo. Los dos alelos están
presentes en la muestra, a frecuencia apreciable, y ninguna molécula lleva ambos.

En un organismo haploide eso no es una variante con dos formas. Son **dos
linajes**, y a nivel de codón es la firma a nivel de lectura de una infección
mixta o de subpoblaciones distintas.

Un ratio de co-ocurrencia de `0.0000` sí lo reporta, y con honestidad: hubo
moléculas que cruzaron el sitio y ninguna llevaba ambas. Lo que no puede hacer es
decir cuánto se aleja eso del azar, así que se lee igual tanto si las dos se
excluyen como si una es simplemente demasiado rara para esperar una molécula
compartida. Un `D'` de `-1` con un p pequeño distingue los dos casos.

## El valor p

Un `D'` de `1.0` con cuatro moléculas y un `D'` de `1.0` con cuatrocientas son
el mismo número y no son la misma evidencia. `Haplotype LD p` es el test exacto
de Fisher a dos colas sobre la tabla 2x2 de recuentos de moléculas, así que dice
si el alejamiento de la independencia es más de lo que la profundidad podría
producir por azar.

Léelos juntos. Un `|D'|` grande con un p cerca de 1 es una muestra pequeña
hablando con seguridad de nada.

## Cuándo está ausente

`-` significa que la pregunta no se pudo responder, y eso es deliberadamente
distinto de una respuesta de cero:

- Sin `--bam`, o una fila de una sola variante: no hay nada que correlacionar.
- **Ninguna molécula observó las variantes juntas.** Habitual cuando un codón
  cruza un intrón y los fragmentos son más cortos que él.
- **Uno de los alelos está en todas las moléculas que observaron, o en
  ninguna.** Sin variación en ese locus no hay correlación que medir, y el
  denominador de `D'` es cero. Reportar `0.0` ahí afirmaría que se midió
  independencia y se encontró.

Ese último caso es más común de lo que parece. Si doce moléculas llevan tres
variantes y ocho llevan solo las dos primeras, esas dos primeras están en las
veinte moléculas que observaron la ventana: la fila sale con sus ocho moléculas,
y su ligamiento es `-`, porque dentro de esta ventana esas dos nunca varían.

## Los haplotipos de indel también lo llevan

El estadístico no es específico de codones. Una fila de haplotipo local de indel
también lo lleva, calculado sobre las variantes que esa fila afirma que van
juntas.

Esa fila antes solo tenía un recuento de lecturas, que no distingue una
asociación real de dos variantes comunes que se encuentran. Dos ventanas, cada
una con una inserción y una sustitución:

| Moléculas con ambas | solo inserción | solo sustitución | ninguna | Event Reads | `D'` | p |
|---|---|---|---|---|---|---|
| 5 | 15 | 10 | 20 | 5 | -0.1667 | 0.75 |
| 18 | 2 | 2 | 18 | 18 | 0.8000 | 5.3e-7 |

Las dos filas son reales: cinco y dieciocho moléculas llevan de verdad esas dos
variantes juntas. En la primera la inserción está en el 40% de las moléculas y la
sustitución en el 30%, y encontrarse en el 10% es lo que predice el azar. En la
segunda las dos están en la mitad y se encuentran en el 45% donde el azar predice
el 25%.

## Leerlo junto al recuento de lecturas

Responden a preguntas distintas y ninguno sustituye al otro.

- **`Event Reads` / `MNV Reads`**: cuántas moléculas *son* esa combinación.
- **`Haplotype LD`**: si sus variantes tienen algo que ver entre sí.

Un haplotipo con pocas moléculas puede estar perfectamente ligado, y uno con
muchas puede ser una casualidad. El primero es una especie minoritaria real que
merece reportarse; el segundo son dos variantes que resultan ser comunes en la
misma muestra.

## Detalles del cálculo

- La tabla se construye sobre **moléculas**, no lecturas: los dos mates de un
  fragmento paired-end son una molécula. Mira
  [Indels y haplotipos locales](indel-haplotypes.es.md) para qué cuenta como una.
- Solo se incluyen las moléculas que observaron **todas** las posiciones de la
  fila, en los dos lados de la tabla. Una molécula que paró entre dos posiciones
  no vio evidencia sobre el par.
- Con tres o más variantes se puntúa cada par y **decide el más débil**, porque
  la fila afirma que una molécula las lleva todas, y esa afirmación vale lo que
  su par menos ligado.
- `r²` se calcula a la vez y está disponible en la API. Se reporta `D'` porque
  la pregunta aquí es "están en las mismas moléculas", que es a lo que responde
  `D'`; `r²` refleja además cuánto se separan las dos frecuencias, que es otra
  cosa que saber.

## Cómo se comprueba esto

[`tests/phasing/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/tests/phasing)
construye alineamientos molécula a molécula, eligiendo qué lecturas llevan qué
variantes, y compara contra `D'` y un test exacto derivados por segunda vez en
la propia suite. Barre tablas 2x2 que cubren la independencia, las dos
direcciones de asociación, y las tablas degeneradas donde `D'` no tiene valor,
junto a los barridos de soporte de fase.

## Límites

- El estadístico es **dentro de una ventana**: un codón, o una ventana de
  haplotipo local. get_MNV no infiere fase de largo alcance, y una variante
  fuera de la ventana no forma parte de la tabla.
- Describe las moléculas que se secuenciaron. El ligamiento a nivel de lectura
  en un esquema de amplicones está acotado por el amplicón, y un par separado
  por más de un fragmento simplemente no tiene moléculas en común que contar.
- No dice cuál de dos linajes en competencia es cuál. Un `D'` cerca de `-1`
  reporta que existen dos; identificarlos es trabajo de una herramienta de
  fasing o de deconvolución de cepas.
