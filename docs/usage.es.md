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

Usa `--vcf` para entrada `.vcf` sin comprimir, `.vcf.gz` comprimida con BGZF o
`.bcf`, y `--tsv` para el archivo `variants.tsv` que genera `ivar variants`.

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
viaja como un solo adjunto. Necesita la salida TSV, que es la de por defecto.
`--convert` escribe el VCF *en lugar* del TSV, así que usa `--both` cuando
quieras informe y VCF a la vez.

Para una cohorte ya procesada con una muestra por ejecución, genera el informe a
partir de las salidas existentes en vez de volver a ejecutar el pipeline:

```bash
get_mnv --report-from run1.MNV.tsv run2.MNV.tsv --report cohort.html
```

## Marcar mutaciones de resistencia conocidas en *M. tuberculosis*

`example/who_catalogue_2023.tsv` es el catálogo de la OMS de mutaciones de
resistencia de *M. tuberculosis* (2023), sus grupos 1 y 2, en la forma que lee
`--catalogue`. Sobre la referencia y la tabla de genes del ejemplo:

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta example/MTB_ancestor.fas \
  --genes example/anot_genes.txt \
  --catalogue example/who_catalogue_2023.tsv
```

Una llamada que coincide lleva la entrada en su columna `Catalogue`. Un codón
con dos sustituciones se compara como el codón que forma, algo que la anotación
base a base no ve: rpoB `p.Ser450Phe` necesita dos cambios en el codón 450, y
su fila dice

```text
RIF: Assoc w R (rpoB p.Ser450Leu, per-SNV reading only); RIF: Assoc w R (rpoB p.Ser450Phe, codon reading only)
```

es decir, el codón es Ser450Phe, y un anotador que lee una base cada vez habría
dado Ser450Leu. Para `p.Ser441Gln` habría dado Ser441Pro y una parada
prematura. Para conservar otros grupos, o para colocar la tabla sobre otra
referencia, constrúyela a partir de los archivos de la OMS:

```bash
base="https://raw.githubusercontent.com/GTB-tbsequencing/mutation-catalogue-2023/main/Final%20Result%20Files"
curl -LO "$base/WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
curl -LO "$base/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
python3 tools/who_catalogue_to_tsv.py WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt \
  --grades 1,2,5 \
  --coordinates WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt \
  --fasta example/MTB_ancestor.fas \
  --unconverted who_unconverted.tsv \
  -o who_catalogue.tsv
```

Qué reescribe el conversor, qué no puede convertir y qué necesita una tabla de
genes o una referencia propias está en
[Formatos de entrada](input-formats.es.md#el-catálogo-de-la-oms-de-mutaciones-de-resistencia-de-m-tuberculosis).

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
- Si usas `--genes` con datos de varios contigs, nombra el contig de cada gen en
  la séptima columna del TSV de anotación, o usa `--gff`.
- El comportamiento con indels y frameshift tiene sus propios ajustes, y dos de
  sus valores por defecto se apartan a propósito del comportamiento original de
  la herramienta. El razonamiento está en
  [Alcance y compatibilidad](indel-mnv-semantics.es.md#parámetros-de-ajuste).
