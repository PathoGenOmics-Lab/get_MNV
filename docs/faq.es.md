# Preguntas frecuentes

## ¿Es get_MNV un llamador de variantes (variant caller)?

No. get_MNV **anota y resume llamadas de variantes existentes** (de un VCF o
un TSV de iVar) frente a una referencia. No llama variantes ni ensambla
lecturas: reinterpreta los alelos ya presentes en tu entrada, con conciencia de
codones.

## ¿Qué entradas necesito?

Un FASTA de referencia, un archivo de variantes (`--vcf` o `--tsv`) y una
anotación de genes (`--gff` o `--genes`). Un BAM (`--bam`) es opcional y solo se
necesita para el soporte de lecturas y los filtros basados en lecturas. Consulta
[Formatos de entrada](input-formats.es.md).

## ¿Cuándo debo usar `--gff` frente a `--genes`?

Usa `--genes` con un TSV sencillo de cuatro columnas (`name`, `start`, `end`,
`strand`) para genes de un solo exón / procariotas. Usa `--gff` con un archivo
GFF/GFF3 cuando necesites transcritos empalmados; añade `--gff-features CDS` para
que los codones se construyan a partir de los segmentos CDS unidos (incluyendo
codones que abarcan uniones de exones).

## ¿Qué `--translation-table` debo usar?

El valor por defecto es `11` (bacteriano), correcto para organismos como
*M. tuberculosis*. Para el código nuclear estándar usa `1`; otras tablas
admitidas son 2, 3, 4, 5, 6, 12 y 25. Consulta la
[Referencia de la CLI](cli-reference.es.md).

La tabla también decide qué codones pueden empezar una proteína. En la tabla 11
empiezan una GTG, TTG, CTG, ATT, ATC y ATA además de ATG, tal como los lista el
NCBI, y el primer codón de un gen se lee como metionina cuando es uno de ellos,
así que un gen que empieza por GTG se anota desde `Met1`, no desde `Val1`.

## ¿Por qué dos de mis SNP se convirtieron en un MNV?

Cuando dos o más SNV caen en el mismo codón, el efecto sobre el aminoácido
depende del codón **combinado**, no de los cambios individuales. get_MNV los
agrupa e informa una fila `MNV` (o `SNP/MNV`) con el codón y el aminoácido
reales. Consulta [Alcance y compatibilidad](indel-mnv-semantics.es.md).

## ¿Necesito un BAM?

Solo si quieres soporte de lecturas. Sin `--bam` igual obtienes la anotación
completa de codones/MNV e indels; con él también obtienes conteos de lecturas
por evento, frecuencias y métricas de hebra, y puedes aplicar filtros basados en
lecturas.

## ¿Por qué las frecuencias son distintas de las de mi VCF de entrada?

Los filtros basados en lecturas y el soporte informado se recalculan a partir de
`--bam`, no se toman de los `OFREQ`/`ODP` originales. Esto refleja lo que las
lecturas muestran realmente para cada SNV, haplotipo MNV y evento indel.

## ¿Puedo usar un archivo BCF como entrada?

Sí, desde la 1.1.6: pásalo a `--vcf` como cualquier VCF. Se decodifica a texto
VCF a medida que se lee, así que se anota exactamente igual que el VCF que
codifica. Las versiones anteriores necesitan convertirlo antes con
`bcftools view input.bcf -O v -o input.vcf`.

## Mi VCF tiene registros multialélicos: ¿qué ocurre?

get_MNV divide cada registro en alelos ALT independientes, con una fila anotada
por ALT, que es la misma salida que da `bcftools norm -m -`. Antes de la 1.1.6 la
ejecución se detenía salvo que se pasara `--split-multiallelic`; la opción se
sigue aceptando y no cambia nada.

## ¿Se admiten indels?

Sí. get_MNV descompone los alelos `REF/ALT` en componentes SNV, MNV, inserción,
deleción, delins y complejo, informa su efecto a nivel de proteína cuando se
solapan con una característica codificante y (con `--bam`) cuenta el soporte
exacto del evento indel. Con `--left-align` además alinea a la izquierda cada
indel sobre la referencia, como hace `bcftools norm -f ref.fa`, que es lo que
necesita un indel en un homopolímero o una repetición para encontrar sus
lecturas cuando la entrada lo escribió en otro sitio que el alineador.
`--normalize-alleles` solo recorta el contexto compartido y nunca mueve un indel.

## ¿Y las variantes intergénicas?

Se incluyen por defecto y se etiquetan como `intergenic`. Usa
`--exclude-intergenic` para descartar variantes fuera de los genes anotados.

## ¿Dónde va la salida?

Por defecto get_MNV escribe `<input_name>.MNV.tsv` en el directorio de trabajo actual, no junto a la entrada. Usa
`--convert` para VCF, `--both` para ambos, y las opciones de salida relacionadas
en la [Referencia de la CLI](cli-reference.es.md).

## La app de macOS no abre ("desarrollador no identificado")

La aplicación de escritorio no está firmada con un certificado de Apple
Developer. En el primer arranque, haz clic derecho en la app → **Abrir** → haz
clic en **Abrir** en el cuadro de diálogo. Consulta la guía de la
[GUI de escritorio](gui.es.md).

> [!NOTE]
> **Tu pregunta no está aquí**
>
> Pregúntala en [Q&A](https://github.com/PathoGenOmics-Lab/get_MNV/discussions/categories/q-a). Las preguntas sobre por qué una variante salió
> anotada de una forma concreta son bienvenidas ahí, y se responden mucho
> mejor con el registro y el comando que la produjeron.
