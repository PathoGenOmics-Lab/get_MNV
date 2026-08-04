# Uso

Esta página muestra los comandos más habituales y el significado de los argumentos principales.

## Comando básico

```bash
get_mnv \
  (--vcf <VCF_FILE> | --tsv <IVAR_TSV_FILE>) \
  --fasta <REFERENCE_FASTA> \
  (--gff <ANNOTATION_GFF> | --genes <ANNOTATION_TSV>)
```

Usa `--vcf` para entrada `.vcf` sin comprimir o `.vcf.gz` comprimida con BGZF, y
`--tsv` para el archivo `variants.tsv` que genera `ivar variants`. La entrada BCF
no se admite directamente; conviértela antes a VCF con `bcftools view`.

## Recetas habituales

### Entrada VCF

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3
```

### Entrada TSV de iVar

```bash
get_mnv \
  --tsv sample_variants.tsv \
  --fasta reference.fasta \
  --gff genes.gff3
```

### Añadir soporte de lecturas desde un BAM

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3
```

### Escribir tanto TSV como VCF

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both
```

### Analizar características CDS de un GFF

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --gff-features CDS
```

Usa `--gff-features CDS` cuando quieras una anotación proteica consciente del
codón a partir de las características CDS, sobre todo con archivos GFF/GTF
eucariotas. Las filas CDS con `transcript_id` o `Parent` se reconstruyen como
modelos de CDS de transcrito empalmado.

## Argumentos obligatorios

| Argumento | Significado |
|---|---|
| `--vcf <FILE>` | Llamadas de variantes en formato `.vcf` sin comprimir o `.vcf.gz`. |
| `--tsv <FILE>` | Llamadas `variants.tsv` de iVar. |
| `--fasta <FILE>` | FASTA de referencia con el que se llamaron las variantes. |
| `--gff <FILE>` | Anotación de genes en formato GFF/GFF3/GTF. |
| `--genes <FILE>` | TSV sencillo de anotación de genes. Úsalo en lugar de `--gff`. |

Debes proporcionar `--gff` o `--genes`.

## Argumentos de entrada

| Argumento | Por defecto | Significado |
|---|---:|---|
| `--bam <FILE>` | ninguno | BAM ordenado e indexado con el que se cuenta el soporte de lecturas. |
| `--sample <NAME>` | primera muestra | Muestra que se lee de un VCF multimuestra. Usa `all` para todas las muestras. |
| `--chrom <NAME>` | todos los contigs | Restringe la ejecución a un solo contig. |
| `--gff-features <LIST>` | `gene,pseudogene` | Tipos de característica que se analizan del GFF/GTF. |
| `--translation-table <N>` | `11` | Tabla del código genético del NCBI. Admitidas: `1,2,3,4,5,6,11,12,25`. |

## Argumentos de filtrado

| Argumento | Por defecto | Significado |
|---|---:|---|
| `--quality <N>` | `20` | Calidad Phred mínima de la base para el soporte de lecturas del BAM. |
| `--min-mapq <N>` | `0` | Calidad de mapeo mínima de las lecturas del BAM. |
| `--snp <N>` | `0` | Mínimo de lecturas que respaldan el SNP. |
| `--mnv <N>` | `0` | Mínimo de lecturas que respaldan el MNV. |
| `--min-snp-frequency <F>` | `0` | Frecuencia alélica mínima del SNP derivada del BAM (`0` a `1`). |
| `--min-mnv-frequency <F>` | `0` | Frecuencia mínima del haplotipo MNV derivada del BAM (`0` a `1`). |
| `--min-snp-strand <N>` | `0` | Mínimo de lecturas del SNP en cada hebra. |
| `--min-mnv-strand <N>` | `0` | Mínimo de lecturas del MNV en cada hebra. |
| `--min-strand-bias-p <P>` | `0` | p-valor exacto de Fisher mínimo para el filtrado por sesgo de hebra. |
| `--strict` | desactivado | Falla cuando faltan las métricas originales de profundidad o frecuencia. |

Los filtros de frecuencia requieren `--bam`, porque usan el soporte de lecturas
recalculado por get_MNV. No filtran según el valor `OFREQ` original de la
entrada. Por ejemplo, `--min-snp-frequency 0.05` conserva los registros SNP con
un 5% o más, y `--min-mnv-frequency 0.20` conserva los haplotipos MNV con un 20%
o más. Estos dos umbrales son independientes: en llamadas mixtas `SNP/MNV`, el
umbral del SNP no elimina un haplotipo MNV fuerte, y el umbral del MNV no elimina
las observaciones SNP que superan el umbral del SNP. Los filtros de recuento de
lecturas y de soporte por hebra siguen la misma regla: `--snp` y
`--min-snp-strand` se aplican a las observaciones SNP, mientras que `--mnv` y
`--min-mnv-strand` se aplican a los haplotipos MNV.

## Ajuste de la anotación de indels

Estos parámetros opcionales afinan la anotación de indels y frameshift. Los
valores por defecto reproducen el comportamiento histórico, así que los comandos
existentes no se ven afectados. Consulta
[indel-mnv-semantics.md](indel-mnv-semantics.es.md) para conocer el fundamento
biológico.

| Argumento | Por defecto | Significado |
|---|---:|---|
| `--frameshift-min-freq <F>` | `0.0` | Frecuencia alélica mínima que debe alcanzar un indel *aguas arriba* para marcar como frameshift los codones SNV/MNV aguas abajo. Con `0.0` la marca se propaga desde cualquier indel; súbela (p. ej. a `0.5`) para no reetiquetar sustituciones aguas abajo de alta frecuencia por culpa de un indel aguas arriba de baja frecuencia (datos intrahospedador). |
| `--indel-anchor-depth` | desactivado | Cuenta la profundidad del locus del indel (el denominador de EDP/EFREQ) a partir de las lecturas que observan la base de anclaje, en lugar de solo las que abarcan por completo el alelo REF. Reduce el subconteo de profundidad en deleciones multibase. Requiere `--bam`. |
| `--phased-indel-min-reads <N>` | `2` | Mínimo de lecturas del BAM necesarias para emitir una fila de haplotipo de indel en fase o complejo (indel+SNV). Una sola lectura no es evidencia de un haplotipo; pon `1` para emitir toda combinación que alguna lectura muestre. Requiere `--bam`. |
| `--phased-indel-min-freq <F>` | `0.0` | Frecuencia mínima derivada del BAM necesaria para emitir una fila de haplotipo de indel en fase o complejo. Requiere `--bam`. |

## Argumentos de salida

| Argumento | Significado |
|---|---|
| `--convert` | Escribe VCF en lugar de TSV. |
| `--both` | Escribe tanto TSV como VCF. |
| `--vcf-gz` | Escribe salida comprimida `.MNV.vcf.gz`. |
| `--index-vcf-gz` | Crea un índice Tabix para `.MNV.vcf.gz`. |
| `--bcf` | Escribe también salida BCF. Requiere salida VCF. |
| `--emit-filtered` | En la salida VCF, conserva los registros que no superan los filtros y los marca en `FILTER`. La salida TSV sigue omitiendo las filas rechazadas. |
| `--strand-bias-info` | Añade p-valores de sesgo de hebra a los campos INFO del VCF. |
| `--keep-original-info` | Conserva los campos INFO ajenos a get_MNV del VCF de entrada. Requiere salida VCF. |
| `--exclude-intergenic` | Omite las variantes que quedan fuera de las características anotadas. |
| `--summary-json <FILE>` | Escribe un resumen de la ejecución en JSON. |
| `--error-json <FILE>` | Escribe los detalles del error en JSON si la ejecución falla. |
| `--run-manifest <FILE>` | Escribe el comando, la versión, las entradas, las salidas y las sumas de verificación. |

## Argumentos de utilidad

| Argumento | Significado |
|---|---|
| `--dry-run` | Valida las entradas sin escribir archivos de salida. |
| `--threads <N>` | Número de hilos de trabajo. Por defecto: automático. |
| `--normalize-alleles` | Recorta el contexto REF/ALT compartido antes de procesar. |
| `--split-multiallelic` | Divide los registros VCF multialélicos dentro de get_MNV. Cada ALT se convierte en una fila de anotación independiente, incluidos los alts que comparten la misma posición del codón. |

## Notas

- Los nombres de los contigs deben coincidir exactamente entre el archivo de variantes, el FASTA, el GFF y el BAM.
- El análisis del TSV de iVar conserva las filas SNV e indel que pasan. La
  notación de indels como `+SEQ` o `-SEQ` se convierte en alelos `REF/ALT`
  anclados al estilo VCF usando la referencia FASTA.
- Si usas `--genes`, el TSV de anotación no tiene columna de contig. Para datos
  multicontig, usa preferiblemente `--gff`.
