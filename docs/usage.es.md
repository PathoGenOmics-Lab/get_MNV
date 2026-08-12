# Recetas habituales

Comandos listos para ejecutar de las tareas que más se piden. Cada opción, con
su valor por defecto y su significado, vive en un único sitio: la
[referencia de CLI](cli-reference.es.md).

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

## Entrada VCF

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3
```

## Entrada TSV de iVar

```bash
get_mnv \
  --tsv sample_variants.tsv \
  --fasta reference.fasta \
  --gff genes.gff3
```

## Añadir soporte de lecturas del BAM

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3
```

Sin BAM, get_MNV anota lo que reportó el llamador de variantes. Con BAM cuenta
las lecturas por su cuenta, y solo entonces puede distinguir un haplotipo real a
nivel de codón de dos sustituciones que nunca compartieron una molécula. Mira
[Ligamiento](linkage.es.md) para ver qué te da eso.

## Filtrar por el soporte recontado

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --min-snp-frequency 0.05 \
  --min-mnv-frequency 0.20
```

Los filtros de frecuencia y de número de lecturas usan el soporte que get_MNV
recalcula del BAM, no los valores `OFREQ`/`ODP` de la entrada, así que necesitan
`--bam`. Los umbrales de SNP y de MNV son independientes: un haplotipo MNV fuerte
sobrevive aunque sus sustituciones individuales queden por debajo del umbral de
SNP.

## Escribir TSV y VCF a la vez

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both
```

## Generar el informe interactivo

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --report sample.html
```

El informe es un único archivo HTML autocontenido, así que se abre sin servidor y
viaja como un solo adjunto. Necesita la salida TSV, que es la de por defecto; si
añades `--convert`, pide también `--both`.

Para una cohorte ya procesada con una muestra por ejecución, genera el informe a
partir de las salidas existentes en vez de volver a ejecutar el pipeline:

```bash
get_mnv --report-from run1.MNV.tsv run2.MNV.tsv --report cohort.html
```

## Analizar features CDS de un GFF

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --gff-features CDS
```

Usa `--gff-features CDS` cuando quieras anotación proteica consciente de codones
a partir de features CDS, sobre todo en archivos GFF/GTF eucariotas. Las filas
CDS con `transcript_id` o `Parent` se reconstruyen como modelos CDS de
transcrito con splicing.

## Notas

- Los nombres de contig deben coincidir exactamente entre el archivo de
  variantes, el FASTA, el GFF y el BAM.
- El parseo de TSV de iVar conserva las filas SNV e indel que pasan los filtros.
  La notación de indels tipo `+SEQ` o `-SEQ` se convierte a alelos `REF/ALT`
  anclados al estilo VCF usando el FASTA de referencia.
- Si usas `--genes`, el TSV de anotación no tiene columna de contig. Para datos
  con varios contigs, mejor `--gff`.
- El comportamiento con indels y frameshift tiene sus propios ajustes, y dos de
  sus valores por defecto se apartan a propósito del comportamiento original de
  la herramienta. El razonamiento está en
  [Alcance y compatibilidad](indel-mnv-semantics.es.md#parametros-de-ajuste).
