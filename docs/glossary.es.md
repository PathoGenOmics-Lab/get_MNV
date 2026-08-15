# Glosario

Las palabras que usa esta documentación, en el sentido en que las usa get_MNV.
Los términos en negrita aparecen como columnas o campos INFO en la salida; en
[Formatos de salida](output-formats.es.md) está dónde se escribe cada uno.

## El cambio en sí

**SNV** (variante de un solo nucleótido)
: Una base de la referencia sustituida por otra. La salida también la llama
  **SNP** en la columna `Variant Type`.

**MNV** (variante multinucleotídica)
: Dos o más sustituciones que caen dentro de un mismo codón y se leen juntas, así
  que el aminoácido sale del codón entero y no de cada base por separado. Es lo
  que get_MNV existe para encontrar.

**SNP/MNV**
: Una fila cuyo codón lleva varias sustituciones, con las dos lecturas: el
  aminoácido que daría cada sustitución por su cuenta y el que da el codón con
  todas. Una fila de tipo `MNV` a secas no tiene lectura por base que informar.

**Codón**
: Tres bases codificantes seguidas, leídas en la dirección del transcrito, que
  especifican un aminoácido. En la hebra menos el transcrito lee de coordenadas
  altas a bajas, así que la primera base de un codón es la de mayor coordenada de
  las tres.

**Indel**
: Una inserción o una deleción. Un VCF escribe las dos ancladas en una base que
  no cambia, así que `T>TGCT` inserta `GCT` después del ancla y `TGCT>T` borra las
  tres bases que la siguen.

**delins**
: Un cambio que borra e inserta a la vez, como `AT>GGG`. Ni sustitución pura ni
  indel puro.

**Indel complejo**
: Un alelo cuya descomposición lleva un indel junto con al menos una sustitución,
  así que no se puede describir como un solo evento simple.

**Ancla**
: La base que un registro VCF nombra para situar un indel: la base a la que sigue
  la secuencia insertada, o la anterior a las borradas. El indel no la cambia, y
  por eso get_MNV pregunta dónde *actúa* un registro y no dónde empieza.

**Clase de evento** y **componentes del evento**
: Cómo leyó get_MNV la pareja `REF`/`ALT`: la clase es la forma del alelo entero
  (`snp`, `mnv`, `insertion`, `deletion`, `delins`, `complex_indel`) y los
  componentes son los cambios individuales en los que se descompone, escritos como
  `SNV:110:A>C` o `INS:200:+GG`.

## Dónde cae el cambio

**Intergénico**
: Fuera de toda feature anotada. get_MNV escribe el marcador `intergenic` en la
  columna `Gene` para estas filas, que no es el nombre de un gen, y
  `--exclude-intergenic` quita solo estas.

**Desplazamiento de marco** y **en marco**
: Un indel cuya longitud no es múltiplo de tres desplaza el marco de lectura de
  todo lo que viene después en esa feature; uno cuya longitud sí lo es lo deja
  intacto y añade o quita residuos enteros.

**Donante de splice**, **aceptor de splice**, **región de splice**
: Las dos primeras bases de un intrón, las dos últimas, y las bases justo
  alrededor de una unión. Las dos primeras son los sitios esenciales y puntúan
  impacto HIGH; la región de alrededor puntúa LOW.

**NMD** (degradación mediada por mutaciones terminadoras)
: La vía por la que la célula destruye un transcrito que termina antes de tiempo.
  get_MNV la predice solo donde la parada es de verdad prematura, con la regla de
  que una parada a más de 50 nucleótidos de la última unión de exones la
  desencadena.

**Tabla de traducción**
: Con qué código genético del NCBI se leen los codones, elegido con
  `--translation-table`. Decide tanto los aminoácidos como qué codones terminan
  una proteína: `TGA` es parada en la tabla 11, codifica triptófano en la 2 y
  glicina en la 25.

## Qué dice la evidencia

**Haplotipo**
: Un conjunto de cambios que lleva la misma molécula. Las sustituciones de un
  codón forman un haplotipo solo si la misma molécula de ADN las lleva todas.

**Fasado**
: Colocar los cambios sobre moléculas. get_MNV nunca fasa por su cuenta: lee lo
  que declaró el llamador y cuenta lo que muestran las lecturas.

**Fase declarada**
: Lo que afirma el VCF de entrada, tomado de un `GT` separado por `|` y de su
  conjunto de fase `PS`. Un genotipo separado por `/` no está fasado y no afirma
  nada.

**Ligamiento**, **D'** y **r²**
: Cuánto viajan juntos dos cambios en las lecturas, más allá de lo que producirían
  por azar sus frecuencias individuales. Que dos variantes frecuentes coincidan a
  menudo no es prueba de un haplotipo; mira [Ligamiento](linkage.es.md).

**Soporte de lecturas** y **soporte de evento**
: Cuántas lecturas llevan un cambio. Las sustituciones se cuentan base a base con
  una caché por ventana; un indel se cuenta desde el CIGAR de cada lectura, y sus
  cuentas se escriben en las columnas `Event Reads`.

**Sesgo de hebra**
: Soporte que llega casi entero de una sola hebra, lo que suele apuntar a un
  artefacto y no a una variante. Se informa como un valor p de Fisher exacto con
  `--strand-bias-info` y se filtra con `--min-strand-bias-p`.

## Cómo se etiqueta la consecuencia

**Término SO**
: El nombre de la Sequence Ontology para la consecuencia, como `missense_variant`
  o `splice_donor_variant`. Una fila puede declarar más de uno a la vez, unidos
  con un ampersand, cuando una base es codificante y además está en región de
  splice.

**Impacto**
: El cajón de gravedad en el que cae el término SO: `HIGH`, `MODERATE`, `LOW` o
  `MODIFIER`. Cuando una fila declara dos consecuencias se queda con la más grave.

**Distancia de Grantham**
: Un número para lo químicamente distintos que son dos aminoácidos, que se informa
  para un missense genuino, de modo que un cambio conservador se distinga de uno
  drástico.
