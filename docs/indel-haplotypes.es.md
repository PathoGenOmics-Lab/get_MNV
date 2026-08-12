# Indels y haplotipos locales

Cómo lee get_MNV un indel a partir de los alineamientos, qué cuenta, y qué
significa cada número de la salida. Para la compatibilidad con llamadores y las
reglas de límites, mira [Alcance y compatibilidad](indel-mnv-semantics.es.md);
esta página va de la lectura y del recuento.

## Qué es un indel, posicionalmente

Un indel de VCF está **anclado**. `POS=100 REF=A ALT=AGG` no cambia la posición
100: la `A` es relleno, y las `GG` se insertan entre la 100 y la 101. Una
deleción `POS=100 REF=AGG ALT=A` conserva la 100 y elimina la 101 y la 102.

Eso tiene una consecuencia que get_MNV tiene que respetar en todas partes: una
inserción no ocupa una base de referencia. Vive **entre** dos de ellas.

## Qué tiene que ver una lectura antes de poder opinar

Una lectura solo tiene derecho a opinar sobre un alelo si observó suficiente
referencia como para tener una opinión. Ese tramo no siempre es el del propio
alelo REF:

| Alelo | Tramo REF | Tramo que debe observar | Por qué |
|---|---|---|---|
| Sustitución | la base | la base | No hay nada fuera de ella. |
| Inserción | el ancla | ancla **+ 1** | Una lectura que termina en el ancla ha cubierto todo el tramo REF y no ha visto nada de la unión. Sin la base extra se lee como "aquí está la referencia", que es evidencia de que la inserción **falta** que esa lectura nunca tuvo. |
| Deleción | ancla y bases borradas | ese tramo **+ 1** | Dentro de su propio tramo, una lectura que borra más que este alelo es idéntica a una que borra exactamente eso. Dónde termina la deleción solo se ve desde fuera. |

Una sola función responde a esto para todos los caminos que observan un alelo
(`observation_ref_len`). Antes de existir, el conteo exacto, el descubrimiento
de haplotipos y el fasing indel-contra-sustitución lo decidían por separado, y
discrepaban. La peor combinación cancelaba soporte real: un mate que terminaba
en el ancla de una inserción se leía como "aquí hay referencia" y anulaba el
soporte que había aportado su pareja.

## Soporte exacto del evento

`Event Reads` cuenta moléculas que **reproducen el alelo exactamente** sobre su
tramo de referencia, CIGAR incluido. Una lectura debe producir la misma
secuencia local y contener los componentes de inserción y deleción esperados.
Esto importa en combinaciones de efecto neutro, donde una inserción más una
deleción pueden dar la misma secuencia que una sustitución simple bajo otro
alineamiento.

Dos reglas mantienen honesto ese recuento:

- Una lectura cuya deleción **se pasa de largo** del tramo del alelo lleva una
  deleción más larga y distinta, y no es soporte de esta.
- Una molécula cuyos dos mates **leen alelos distintos** sobre el locus no
  respalda nada. Uno de los dos se equivoca y no hay forma de saber cuál, que es
  la misma postura que toma el conteo de sustituciones ante un solape
  contradictorio. Sin esto, un solo mate con un artefacto de alineamiento
  bastaba para acreditar la molécula entera, y el solape de mates es justo donde
  aparecen los artefactos de realineamiento de indels.

`Event Depth` es la profundidad del locus: por defecto, toda molécula que
observe la base ancla con calidad suficiente. `--legacy-indel-depth` la
restringe a las que cruzan el alelo REF entero, lo que infravalora la
profundidad de una deleción multibase y sesga `Event Frequency` hacia arriba.

## Haplotipos locales

Cuando varias variantes caen cerca, la pregunta deja de ser qué hace cada una y
pasa a ser **cuáles comparten moléculas**. get_MNV lo responde leyendo las
moléculas, no enumerando posibilidades.

### Descubrimiento: lo que enseñan las moléculas

Las variantes cercanas se agrupan en una ventana por proximidad, pero la
proximidad solo *propone*. Cada variante candidata se juzga sobre su propio
tramo de referencia, así que un fragmento paired-end puede responder por una
variante que cubre su primer mate y por otra que cubre el segundo. Una molécula
que no observó todas las variantes de la ventana se aparta como parcial: lleva
lo que lleva, pero meterla en una combinación afirmaría que le faltan las
variantes que nadie miró.

Las moléculas que sí observaron la ventana entera se agrupan por la combinación
que enseñan. Eso es la salida del descubrimiento, y es una **distribución
conjunta**: cuántas moléculas son cada combinación.

De leer en vez de enumerar salen dos reglas:

- Una combinación que ninguna molécula lleva no se emite, **incluida la que es
  subconjunto de otra que sí se lleva**. Una molécula con A, B y C no es
  evidencia de una molécula con solo A y B.
- Dos combinaciones que conviven de verdad salen las dos, cada una con su
  recuento. Un esquema de enumerar y comprobar no puede distinguir eso de una
  sola molécula que lo lleva todo.

### Recuento: cuántas moléculas son este haplotipo

El soporte de una fila es el número de moléculas cuya combinación es
**exactamente** esa, y sale del descubrimiento. El conteo exacto del alelo no
puede responder a eso: empareja sobre el tramo de referencia de la propia
combinación, así que una molécula que lleva esa combinación *y algo más fuera
del tramo* también casa. Se queda con el trabajo que sí puede hacer, que es
confirmar que el alelo es reproducible, y dar la profundidad de la ventana.

La diferencia no es académica. Con doce moléculas que llevan tres variantes y
ocho que llevan solo las dos primeras:

| Event Components | Event Reads | Event Depth | Event Frequency |
|---|---|---|---|
| `SNV:28:G>T \| INS:29:+GCT \| SNV:30:T>A` | 12 | 20 | 0.6000 |
| `SNV:28:G>T \| INS:29:+GCT` | 8 | 20 | 0.4000 |

Contada por coincidencia de alelo, la segunda fila decía `20` al `1.0000`: el
soporte de la especie entera, para un haplotipo que llevan ocho moléculas. Ahora
los recuentos suman la profundidad, que es la imagen coherente de una ventana
con dos especies.

### Qué se emite

- Los propios registros de la entrada, anotados como siempre.
- Una fila `complex_indel` por cada combinación que enseñan las moléculas,
  cuando tiene más de una variante y al menos una es un indel.
- Nada más. `--phased-indel-min-reads` (por defecto `2`) y
  `--phased-indel-min-freq` (por defecto `0.0`) ponen el suelo: una lectura no
  es evidencia de un haplotipo, porque un solo error de secuenciación en una
  posición llamada acuña su propia combinación.

Una ventana que enseña más de 64 combinaciones distintas reporta las mejor
soportadas y dice en el log qué apartó, en vez de truncar en silencio.

## Leer una fila de indel

Para el ejemplo de arriba:

```text
Positions          28
Variant Type       INDEL
Change Type        In-frame Indel
AA Changes         Ala10delinsSerLeu
Event Class        complex_indel
Event Components   SNV:28:G>T | INS:29:+GCT | SNV:30:T>A
Event Reads        12
Event Depth        20
Event Frequency    0.6000
SO Term            inframe_insertion
HGVS g.            g.28_30delinsTCGCTA
```

`Positions` es el ancla del alelo compuesto, no una posición por componente; en
`Event Components` es donde se listan las partes. El cambio de aminoácido es el
efecto del alelo compuesto entero, así que un haplotipo de tres variantes tiene
una consecuencia proteica, no tres.

Una fila de MNV a nivel de codón que solapa un indel se conserva como contexto
posicional con `Change Type = Indel overlap` y `Unknown` en el residuo, porque
el efecto proteico del codón no está definido independientemente del indel. La
fila `complex_indel` de al lado es donde está la consecuencia real.

## Límites que conviene conocer

- **El conteo exacto necesita una observación contigua.** El descubrimiento
  trabaja a nivel de molécula, así que un haplotipo con sus variantes repartidas
  entre mates se encuentra; el alelo no lo puede reproducir ningún mate por
  separado, así que la fila no se puede respaldar y el umbral la descarta.
- **El alineamiento a la izquierda es trabajo de la entrada.** get_MNV empareja
  el alelo tal cual se le da. En un homopolímero o una repetición en tándem, un
  indel de entrada colocado distinto que el alineamiento del BAM da soporte
  exacto cero; un aviso lo dice, y `bcftools norm -f ref.fa` lo arregla.
- **Una inserción anclada en la última base de una feature queda fuera de ella**,
  porque la secuencia insertada cae más allá. En un gen de hebra menos esa misma
  coordenada es la *primera* base del transcrito, así que la inserción cae en el
  5' UTR. La misma regla, por la razón contraria.
- **Sin ensamblaje de novo.** get_MNV combina alelos que la entrada ya declara.
- **Los alelos simbólicos** (`<DEL>`, `<DUP>`) no tienen secuencia que una
  lectura pueda reproducir, así que la evidencia de lecturas no puede decir nada
  de ellos.

## Cómo se comprueba esto

- [`tests/scenarios/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/tests/scenarios)
  construye alineamientos a mano para cada regla de arriba, incluidas las que
  solo muerden en un borde: una lectura que acaba en el ancla de una inserción,
  mates que discrepan sobre un indel, una deleción que se pasa de la declarada,
  y un haplotipo anidado dentro de otro.
- [`tests/truth/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/tests/truth)
  barre una deleción y una inserción de 1 a 3 bases en **cada posición de cada
  exón** de cuatro genes, comprobando la consecuencia proteica contra una
  retraducción independiente.
- [`tests/pileup/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/tests/pileup)
  comprueba profundidad y hebra contra `bcftools mpileup`, que no lo ha escrito
  nadie de aquí.
