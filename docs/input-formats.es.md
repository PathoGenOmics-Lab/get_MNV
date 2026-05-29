# Formatos de entrada

get_MNV necesita tres archivos: llamadas de variantes, un FASTA de referencia y
anotación de genes. Un archivo BAM es opcional.

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

Usa `--vcf` para archivos `.vcf` simples o `.vcf.gz` comprimidos con BGZF y
`--tsv` para archivos `variants.tsv` de iVar. La entrada BCF no se acepta
directamente; conviértela primero, por ejemplo con
`bcftools view input.bcf > input.vcf`.
Los comandos antiguos que pasan un TSV de iVar a través de `--vcf` se siguen
detectando automáticamente cuando la cabecera tiene las columnas estándar de
iVar.

### VCF

Usa un archivo VCF estándar que contenga llamadas SNV/MNV, indels o alelos
complejos.

Requisitos:

- Los nombres de contig del VCF deben coincidir con los nombres de contig del
  FASTA y del GFF/GTF.
- Los alelos REF deben coincidir con la secuencia del FASTA.
- Los registros multialélicos deben dividirse previamente, o bien ejecutar con
  `--split-multiallelic`.

get_MNV puede leer las métricas originales de profundidad/frecuencia desde
campos INFO o FORMAT comunes, incluidos `DP`, `AF`, `FREQ`, `AD`, `AO` y `RO`.

Estos valores de frecuencia de entrada se conservan para el reporte como
`OFREQ`. Los filtros de frecuencia de la línea de comandos
(`--min-snp-frequency`, `--min-mnv-frequency`) usan en su lugar el soporte de
lecturas derivado del BAM, por lo que requieren `--bam`.

### TSV de iVar

Usa el TSV producido por `ivar variants`.

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
| `REF_DP`, `ALT_DP` | Se usan para inferir profundidad/frecuencia si es necesario |
| `PASS` | Se usa para conservar las filas que pasan |

Filtrado:

- Si existe `PASS`, get_MNV conserva los valores verdaderos como `TRUE`, `PASS`,
  `1` o `YES`.
- Las filas donde `REF == ALT` se omiten.
- La notación de indels de iVar como `+SEQ` o `-SEQ` se convierte a alelos
  anclados al estilo VCF usando la referencia FASTA, y luego se analiza con el
  mismo modelo de eventos de alelo que la entrada VCF.
- `ALT_FREQ` se reporta como frecuencia original (`OFREQ`). Es independiente de
  los filtros de frecuencia derivados del BAM.

## 2. FASTA de referencia

Pasa la referencia con:

```bash
--fasta reference.fasta
```

Requisitos:

- Los IDs de registro del FASTA deben coincidir con los nombres de contig de las
  variantes.
- Las bases deben ser bases de ADN IUPAC válidas.
- No se permiten nombres de contig duplicados.

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
- Para las features `CDS`, se usa la fase de la columna 8 cuando está presente.
- Para las filas `CDS` con `transcript_id` o `Parent`, get_MNV construye la
  secuencia CDS empalmada para cada transcrito. La agrupación de codones, los
  efectos de aminoácidos de los MNV y el contexto de frameshift de los indels se
  evalúan entonces sobre el CDS completo del transcrito.
- Si un GFF/GTF contiene múltiples transcritos para el mismo gen, una variante
  puede producir una línea de salida por cada transcrito solapante.

Los nombres de genes se leen de atributos comunes como `gene_name`, `gene`,
`Name`, `locus_tag`, `gene_id` e `ID`.

### Anotación TSV simple

Usa `--genes` para un archivo de anotación pequeño y simple:

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

La fase puede ser `0`, `1`, `2` o `.`. Si se omite la columna de fase, su valor
por defecto es `0`.

Limitaciones de la anotación TSV:

- No tiene columna de contig.
- Para datos multi-contig, usa GFF/GTF o restringe la ejecución con `--chrom`.

## 4. Lecturas BAM (opcional)

Pasa las lecturas BAM con:

```bash
--bam reads.bam
```

Cuando se proporciona un BAM, get_MNV calcula:

- Soporte de lecturas de SNP
- Soporte de lecturas del haplotipo MNV
- Profundidad total y frecuencia
- Recuentos de hebra directa/reversa
- Estadísticas opcionales de sesgo de hebra

Requisitos:

- El BAM debe estar ordenado.
- El BAM debe estar indexado (`.bai`).
- Los nombres de contig del BAM deben coincidir con el archivo de variantes y el
  FASTA.
- Las lecturas duplicadas, secundarias y suplementarias se ignoran.

## Nombres de contig

Todos los archivos de entrada deben coincidir en los nombres de contig:

```text
variant contig == FASTA record ID == GFF sequence ID == BAM reference name
```

Por ejemplo, `chr1` y `1` son nombres diferentes. Renombra o normaliza las
entradas antes de ejecutar get_MNV.
