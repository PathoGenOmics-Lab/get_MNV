# Formatos de entrada

get_MNV necesita tres archivos: llamadas de variantes, un FASTA de referencia y
anotación de genes. El archivo BAM es opcional.

## 1. Llamadas de variantes

Pasa las llamadas de variantes con:

```bash
--vcf <VCF_FILE>
# o
--tsv <IVAR_TSV_FILE>
```

Las entradas de llamadas de variantes admitidas son:

- VCF (`.vcf` o `.vcf.gz`)
- TSV de variantes de iVar (`.tsv`)

Usa `--vcf` para archivos `.vcf` simples o `.vcf.gz` comprimidos con BGZF, y
`--tsv` para archivos `variants.tsv` de iVar. La entrada BCF no se acepta
directamente; conviértela antes, por ejemplo con
`bcftools view input.bcf > input.vcf`.
Los comandos antiguos que pasan un TSV de iVar mediante `--vcf` se siguen
detectando de forma automática cuando la cabecera tiene las columnas estándar de
iVar.

### VCF

Usa un archivo VCF estándar que contenga llamadas SNV/MNV, indels o alelos
complejos.

Requisitos:

- Los nombres de contig del VCF deben coincidir con los del FASTA y los del
  GFF/GTF.
- Los alelos REF deben coincidir con la secuencia del FASTA.
- Los registros multialélicos deben dividirse previamente, o bien hay que
  ejecutar con `--split-multiallelic`.

get_MNV puede leer las métricas originales de profundidad y frecuencia desde
campos INFO o FORMAT habituales, como `DP`, `AF`, `FREQ`, `AD`, `AO` y `RO`.

Estos valores de frecuencia de entrada se conservan en el informe como `OFREQ`.
Los filtros de frecuencia de la línea de comandos (`--min-snp-frequency`,
`--min-mnv-frequency`) usan en su lugar el soporte de lecturas derivado del BAM,
por lo que requieren `--bam`.

### TSV de iVar

Usa el TSV generado por `ivar variants`.

Columnas obligatorias:

| Columna | Significado |
|---|---|
| `REGION` | Nombre del contig |
| `POS` | Posición en base 1 |
| `REF` | Base de referencia |
| `ALT` | Base alternativa |

Columnas opcionales que se usan cuando están presentes:

| Columna | Se usa como |
|---|---|
| `TOTAL_DP` | Profundidad original (`ODP`) |
| `ALT_FREQ` | Frecuencia original (`OFREQ`) |
| `REF_DP`, `ALT_DP` | Sirven para inferir profundidad y frecuencia si hace falta |
| `PASS` | Sirve para conservar las filas que pasan el filtro |

Filtrado:

- Si existe `PASS`, get_MNV conserva los valores con sentido verdadero, como
  `TRUE`, `PASS`, `1` o `YES`.
- Las filas donde `REF == ALT` se descartan.
- La notación de indels de iVar como `+SEQ` o `-SEQ` se convierte en alelos
  anclados al estilo VCF usando la referencia FASTA y, a continuación, se analiza
  con el mismo modelo de eventos de alelo que la entrada VCF.
- `ALT_FREQ` se incluye en el informe como frecuencia original (`OFREQ`), y es
  independiente de los filtros de frecuencia derivados del BAM.

## 2. FASTA de referencia

Pasa la referencia con:

```bash
--fasta reference.fasta
```

Requisitos:

- Los IDs de registro del FASTA deben coincidir con los nombres de contig de las
  variantes.
- Las bases deben ser bases de ADN IUPAC válidas.
- No se admiten nombres de contig duplicados.

## 3. Anotación de genes

Proporciona `--gff` o `--genes`.

### GFF/GFF3/GTF

Recomendado para la mayoría de los conjuntos de datos:

```bash
--gff genes.gff3
```

Por defecto, get_MNV analiza las features `gene,pseudogene`. Para la anotación de
codones de regiones codificantes de proteínas, usa:

```bash
--gff-features CDS
```

Detalles importantes:

- Las coordenadas se leen de las columnas 4 y 5.
- La hebra se lee de la columna 7.
- En las features `CDS`, se usa la fase de la columna 8 cuando está presente.
- En las filas `CDS` con `transcript_id` o `Parent`, get_MNV construye la
  secuencia CDS empalmada de cada transcrito. La agrupación de codones, los
  efectos sobre el aminoácido de los MNV y el contexto de frameshift de los
  indels se evalúan entonces sobre el CDS completo del transcrito.
- Si un GFF/GTF contiene varios transcritos para el mismo gen, una sola variante
  puede producir una línea de salida por cada transcrito solapante.

Los nombres de genes se leen de atributos habituales como `gene_name`, `gene`,
`Name`, `locus_tag`, `gene_id` e `ID`.

### Anotación TSV simple

Usa `--genes` para un archivo de anotación pequeño y sencillo:

```bash
--genes genes.tsv
```

Formato de cuatro columnas:

```text
GeneName	GeneStart	GeneEnd	Strand
```

Ejemplo:

```text
Rv0007_Rv0007	9914	10828	+
Rv0008c_Rv0008c	11874	12311	-
```

Formato opcional de cinco columnas con fase:

```text
GeneName	GeneStart	GeneEnd	Strand	Phase
```

La fase puede ser `0`, `1`, `2` o `.`. Si se omite la columna de fase, toma el
valor `0` por defecto.

Formato opcional de seis columnas con biotipo:

```text
GeneName	GeneStart	GeneEnd	Strand	Phase	Biotype
```

```text
Rv0007_Rv0007	9914	10828	+	0	protein_coding
mcr11_RVnc0013	1413094	1413224	-	0	ncRNA
```

Los valores aceptados son `protein_coding`, `coding`, `CDS` y `mRNA` para
features que se traducen, y `ncRNA`, `rRNA`, `tRNA`, `tmRNA`, `miRNA`, `snRNA`,
`snoRNA`, `misc_RNA`, `antisense_RNA`, `SRP_RNA`, `RNase_P_RNA`, `non_coding` y
`pseudogene` para las que no. Cualquier otro valor se rechaza con un error en
vez de adivinarse, porque adivinar mal inventa una proteína o esconde una real.

Una feature no codificante se reporta contra su gen como
`non_coding_transcript_exon_variant` (`MODIFIER`), sin cambio de aminoácido.
**Si se omite la columna de biotipo se asume que toda feature es codificante**,
que es lo que han hecho siempre los ficheros de cuatro y cinco columnas: un gen
de RNA se traduce entonces como si fuera una proteína. Declara el biotipo cuando
tu anotación contenga features no codificantes. La entrada GFF/GTF no se ve
afectada, porque selecciona features `CDS`.

Limitaciones de la anotación TSV:

- No tiene columna de contig.
- Para datos multi-contig, usa GFF/GTF o limita la ejecución con `--chrom`.

## 4. Lecturas BAM (opcional)

Pasa las lecturas BAM con:

```bash
--bam reads.bam
```

Cuando se proporciona un BAM, get_MNV calcula:

- El soporte de lecturas de los SNP
- El soporte de lecturas del haplotipo MNV
- La profundidad total y la frecuencia
- Los recuentos por hebra directa e inversa
- Estadísticas opcionales de sesgo de hebra

Requisitos:

- El BAM debe estar ordenado.
- El BAM debe estar indexado (`.bai`).
- Los nombres de contig del BAM deben coincidir con los del archivo de variantes
  y los del FASTA.
- Las lecturas duplicadas, secundarias y suplementarias se ignoran.

## Nombres de contig

Todos los archivos de entrada deben coincidir en los nombres de contig:

```text
variant contig == FASTA record ID == GFF sequence ID == BAM reference name
```

Por ejemplo, `chr1` y `1` son nombres diferentes. Renombra o normaliza las
entradas antes de ejecutar get_MNV.
