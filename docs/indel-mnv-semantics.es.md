# Semántica de indels y MNV

get_MNV es un anotador de variantes y un generador de resúmenes de haplotipos
para llamadas de variantes ya existentes. No reemplaza a un llamador de
variantes: el VCF o el TSV de iVar de entrada sigue siendo el responsable de
decidir qué alelos existen. Cuando se proporciona un BAM, get_MNV usa la
evidencia de las lecturas para contar el soporte y para emitir haplotipos
combinados exactos de indel/SNV/MNV observados en las mismas lecturas.

## Compatibilidad con llamadores

Distintos llamadores pueden representar el mismo haplotipo local de maneras
diferentes. Por ejemplo, FreeBayes se basa en haplotipos y puede emitir SNPs,
indels, MNPs y eventos complejos como llamadas de haplotipo en bruto. Su propia
documentación recomienda descomponer esas llamadas tras el llamado si se
necesitan registros primitivos de SNP/indel.

bcftools `norm` se ocupa de la otra cara del problema: puede alinear a la
izquierda y normalizar indels con una referencia FASTA, separar sitios
multialélicos y atomizar variantes complejas o MNV en registros más simples. Ese
paso de normalización resulta útil antes de comparar conjuntos de llamadas
procedentes de distintos llamadores.

GATK HaplotypeCaller es un llamador, no un anotador. Recurre al ensamblaje local
y a parámetros de genotipado como la ploidía de la muestra, el número máximo de
haplotipos y la distancia de fusión de MNP. Por ahora, get_MNV no realiza este
descubrimiento de haplotipos de novo ni el modelado de verosimilitud de
genotipos.

## Qué aporta get_MNV

- Agrupación de SNVs por codón en filas SNP, MNV y mixtas `SNP/MNV`.
- Normalización de inserciones/deleciones de iVar en alelos anclados compatibles
  con VCF.
- Efectos proteicos de indels y alelos complejos a partir del feature
  seleccionado, la hebra, la fase del GFF y el offset proteico del transcrito.
- Soporte derivado del BAM para observaciones individuales de SNP/MNV y para
  eventos exactos de indel o indel complejo.
- Filas `complex_indel` adicionales solo cuando el haplotipo combinado completo
  de indel/SNV/MNV se observa en las lecturas del BAM.
- Haplotipos locales leídos de las moléculas. Los eventos cercanos se agrupan en
  una ventana por proximidad, pero la proximidad solo propone: son las lecturas
  que cruzan esa ventana las que deciden qué combinaciones se reportan, y cada
  una sale con el número de moléculas que la llevan. Una combinación que ninguna
  lectura lleva no se emite, incluida la que es subconjunto de otra que sí se
  observa, porque una molécula con A, B y C no es evidencia de una molécula con
  solo A y B. Dos combinaciones que conviven de verdad salen las dos, cada una
  con su recuento. No hay tope al número de variantes que una ventana admite.

## Reglas de límites

Las inserciones de VCF son eventos interbase anclados en la base de referencia
anterior. Una inserción se considera que solapa un feature solo cuando la
secuencia insertada cae entre dos bases de referencia dentro de ese feature. Una
inserción anclada en la última base del feature queda fuera de él.

Las deleciones de VCF están ancladas, pero su efecto biológico es el tramo de
referencia eliminado. Una deleción anclada justo antes de un CDS/gen afecta a ese
feature siempre que las bases eliminadas lo solapen.

Cuando una fila MNV solapa un indel, get_MNV conserva la fila MNV como contexto
posicional y marca su efecto en el aminoácido como `Unknown` con
`Change Type = Indel overlap`. Si el BAM respalda el evento combinado completo,
get_MNV emite por separado una fila exacta `complex_indel`.

El soporte exacto de indels complejos tiene en cuenta el CIGAR. Una lectura debe
producir la misma secuencia ALT local y contener los componentes de
inserción/deleción esperados. Esto importa en las combinaciones netamente
neutras, donde una inserción más una deleción pueden producir la misma secuencia
que un MNV simple bajo otro alineamiento.

## Notas sobre CDS eucariota

Para anotaciones eucariotas GFF/GTF, usa `--gff-features CDS` de modo que get_MNV
pueda aprovechar la fase del GFF y los identificadores de transcrito. Cuando las
filas CDS llevan `transcript_id` o `Parent`, get_MNV construye un modelo de CDS
empalmado para cada transcrito y evalúa la agrupación por codones, los efectos
del MNV en el aminoácido y el contexto de frameshift del indel sobre esa
secuencia codificante completa.

Esto implica que los codones repartidos entre uniones de exones se pueden anotar
como un único MNV a nivel de transcrito, y los SNP aguas abajo se marcan como
frameshift solo cuando el desplazamiento neto del indel codificante aguas
arriba permanece fuera de marco. Las filas CDS sin un identificador de transcrito
utilizable mantienen el comportamiento anterior por característica.

## Limitaciones actuales

- No se acepta un BCF de entrada directamente. Conviértelo antes con
  `bcftools view`.
- El alineamiento a la izquierda de indels y la normalización completa no se
  realizan de forma automática salvo que `--normalize-alleles` logre recortar un
  prefijo o sufijo compartido. Para comparaciones entre llamadores, normaliza
  primero las entradas con una herramienta que tenga en cuenta el FASTA, como
  `bcftools norm -f ref.fa`.
- Los genotipos, la ploidía, los conjuntos de fase y las verosimilitudes de
  genotipo no se vuelven a estimar. En datos eucariotas diploides/poliploides,
  interpreta con cautela los MNV heterocigotos sin fase a nivel de exón salvo que
  la evidencia del VCF o del BAM confirme que los alelos están en el mismo
  haplotipo.
- El ensamblaje local de novo queda fuera del alcance. get_MNV solo combina
  alelos ya presentes en la entrada y, cuando es posible, confirmados por el
  soporte de lecturas del BAM.
- El descubrimiento de haplotipos locales se limita a ventanas de eventos
  cercanos; la reconstrucción de haplotipos muy grandes debería seguir a cargo de
  un llamador dedicado.
- El ensamblaje local completo sigue fuera del alcance: get_MNV reanota los
  alelos presentes en el archivo VCF/iVar de entrada en lugar de descubrir nuevas
  variantes candidatas.

## Parámetros de ajuste

Estos flags opcionales modifican el comportamiento que tiene en cuenta los
indels. Todos adoptan por defecto el comportamiento histórico, de modo que los
pipelines existentes no se ven afectados salvo que se establezca un flag.

- `--frameshift-min-freq <0.0-1.0>` (por defecto `0.0`): frecuencia alélica
  mínima que un indel *aguas arriba* debe alcanzar para desplazar el marco de
  lectura de los codones SNV/MNV aguas abajo (el marcador `(fs)` y los tipos de
  cambio de frameshift). El valor por defecto propaga el efecto desde cualquier
  indel sin importar su frecuencia. En poblaciones intrahospedador o mixtas,
  elevar este umbral evita reetiquetar como frameshifted una sustitución de alta
  frecuencia aguas abajo a causa de un indel de baja frecuencia aguas arriba que
  casi con seguridad reside en otra molécula.
  Cuando se proporciona un BAM, la propagación además se **fasea por lecturas**:
  para cada indel aguas arriba y un SNV aguas abajo dentro del alcance de un read,
  se inspeccionan las lecturas que abarcan ambos loci, y el corrimiento de marco
  *no* se propaga a ese codón cuando las lecturas que llevan el SNV apenas llevan
  el indel (están en trans, en moléculas distintas). Es un refinamiento
  conservador que solo suprime (nunca añade propagación que el gate de frecuencia
  no haría), y solo aplica donde hay lecturas que abarcan ambos loci; más allá del
  alcance de un read sigue gobernando el gate de frecuencia.
- `--indel-anchor-depth` (desactivado por defecto): cuenta la profundidad del
  locus del indel (el denominador de `EDP`/`EFREQ`) a partir de las lecturas que
  observan la base de anclaje, en lugar de solo las que abarcan por completo el
  alelo REF. Reduce el subconteo de profundidad y el sesgo de `EFREQ` en las
  deleciones de varias bases, donde de otro modo las lecturas con solapamiento
  parcial quedan excluidas del denominador.
- `--count-mates-separately` (desactivado por defecto): cuenta los dos mates de
  un fragmento como dos observaciones en vez de una molécula. Por defecto son
  una. Un fragmento es una sola molécula de ADN secuenciada por los dos
  extremos, así que contar los mates por separado cuenta dos veces el solape, y
  trata una variante en cada mate como si no tuvieran relación cuando de hecho
  es la prueba de que van juntas: un codón partido por un intrón puede quedar
  sin respuesta lectura a lectura y resolverse con el par. Donde los mates
  solapan y discrepan sobre una base, uno de los dos se equivoca y no hay forma
  de saber cuál, así que se considera que la molécula no observó esa posición.
  Los datos single-end no se ven afectados.

  El *descubrimiento* de haplotipos locales también es a nivel de molécula:
  cada variante candidata se juzga sobre su propio tramo de referencia, así que
  un fragmento puede responder por una variante que cubre su primer mate y por
  otra que cubre el segundo. El *conteo* exacto del evento no lo es: una lectura
  tiene que reconstruir el alelo compuesto entero en una observación contigua, y
  eso no lo puede hacer ningún mate por separado. Un haplotipo local cuyas
  variantes caen en mates distintos se encuentra pero no se puede respaldar, y
  el umbral de soporte lo descarta en vez de reportarlo sin evidencia.
- `--phased-indel-min-reads <N>` (por defecto `2`) y
  `--phased-indel-min-freq <0.0-1.0>` (por defecto `0.0`): soporte mínimo del BAM
  que una fila de haplotipo indel/complejo con fase debe reunir para ser emitida.
  Una ventana densa en datos intra-host puede albergar varias combinaciones
  genuinas a baja frecuencia; elevar estos valores deja solo las bien
  soportadas.

Cuando un indel con cobertura de lecturas no produce ningún soporte de
CIGAR exacto, get_MNV emite una advertencia: por lo general, esto significa que
el indel de entrada no está alineado a la izquierda del mismo modo que el BAM
(algo habitual en homopolímeros y repeticiones en tándem). Normaliza antes la
entrada con una herramienta que tenga en cuenta el FASTA (`bcftools norm -f
ref.fa`) para que el alelo coincida con el alineamiento de las lecturas.

## Referencias

- FreeBayes: https://github.com/freebayes/freebayes
- bcftools norm: https://samtools.github.io/bcftools/bcftools.html#norm
- GATK HaplotypeCaller: https://gatk.broadinstitute.org/hc/en-us/articles/30332006386459-HaplotypeCaller
