# Resolución de problemas

Esta página enumera los errores comunes y las soluciones más rápidas.

## Los nombres de contig no coinciden

Ejemplo:

```text
[E002] Contig validation failed
```

Solución:

- Asegúrate de usar los mismos nombres de contig en el archivo de variantes, el
  FASTA, el GFF/GTF y el BAM.
- Los nombres distinguen mayúsculas de minúsculas.
- `chr1` y `1` son nombres distintos.

## El REF del VCF no coincide con el FASTA

Ejemplo:

```text
[E002] VCF REF/FASTA mismatch at <contig>:<pos>
```

Solución:

- Comprueba que el VCF se generó a partir del mismo FASTA de referencia.
- Comprueba que las coordenadas y los nombres de contig no se modificaron después
  de la llamada de variantes.

## No se detecta el TSV de iVar

Solución:

- Pasa el archivo con `--tsv sample_variants.tsv`.
- Comprueba que la cabecera del TSV contiene al menos `REGION`, `POS`, `REF` y
  `ALT`.
- Si el archivo es en realidad un TSV de anotación de genes, pásalo con
  `--genes`, no con `--tsv`.

## Bases inválidas

Ejemplo:

```text
[E002] Invalid base 'X' in REF/ALT allele
```

Solución:

- Elimina los alelos inválidos del archivo de variantes.
- get_MNV acepta bases IUPAC de ADN, pero no símbolos arbitrarios en REF/ALT.

## Registros VCF multialélicos

Ejemplo:

```text
Multiallelic VCF record is not supported
```

Solución:

```bash
get_mnv ... --split-multiallelic
```

o divide el VCF de antemano:

```bash
bcftools norm -m - input.vcf > split.vcf
```

`--split-multiallelic` emite una fila anotada por cada ALT, incluso cuando
varios alts comparten la misma posición de codón. Cada fila indica su propio
efecto sobre el aminoácido, su codón y el soporte de lecturas derivado del BAM.
Dividir el VCF de antemano con `bcftools norm -m -` produce la misma salida.

## Entrada BCF

Ejemplo:

```text
BCF input is not supported. Convert to VCF first
```

Solución:

```bash
bcftools view input.bcf > input.vcf
get_mnv --vcf input.vcf ...
```

`--bcf` es solo una opción de salida; no hace que el BCF sea válido como
entrada.

## Anotación TSV con múltiples contigs

Ejemplo:

```text
TSV annotation does not include contig names
```

Solución:

- Usa `--gff` para datos de múltiples contigs.
- O restringe la ejecución a un solo contig con `--chrom`.

## No se encuentra el nombre de la muestra

Ejemplo:

```text
Sample '<name>' not found in VCF header
```

Solución:

- Comprueba el nombre de la muestra en la cabecera del VCF.
- Omite `--sample` para usar la primera muestra.
- Usa `--sample all` solo cuando el VCF tiene columnas de muestra.

## El modo estricto falla

Ejemplo:

```text
--strict enabled, but original VCF metrics are missing
```

Solución:

- Desactiva `--strict`, o
- Asegúrate de que cada variante tenga métricas de profundidad y frecuencia que
  get_MNV pueda leer.

## Filtrado por frecuencia de alelos

Usa `--min-snp-frequency <F>` para los registros SNP y `--min-mnv-frequency <F>`
para los haplotipos MNV. Los valores son fracciones de `0` a `1`, así que `0.05`
significa el 5 %.

Estos filtros requieren `--bam` porque get_MNV los calcula a partir del soporte
de lecturas. No usan el valor `OFREQ` original del VCF/iVar.
Los umbrales de SNP y MNV son independientes. En llamadas mixtas `SNP/MNV`,
`--min-snp-frequency` filtra las observaciones de SNP y `--min-mnv-frequency`
filtra el haplotipo MNV en fase; un MNV fuerte no debería desaparecer solo
porque las observaciones individuales de SNP queden por debajo del umbral de SNP.

Soluciones comunes:

- Si ves un error `requires --bam`, añade un BAM ordenado e indexado o elimina
  los filtros de frecuencia.
- Si quieres filtrar por la frecuencia de alelos original del programa de llamada
  de variantes (`OFREQ`), prefiltra el VCF o el TSV de iVar antes de ejecutar
  get_MNV.
- Combina los filtros de frecuencia con filtros de soporte de lecturas como
  `--snp`, `--mnv`, `--min-snp-strand` y `--min-mnv-strand` para obtener llamadas
  más estrictas.

## El directorio de salida no admite escritura

Ejemplo:

```text
Read-only file system
```

Solución:

- Ejecuta el comando desde una carpeta con permiso de escritura, o
- En la GUI, elige un directorio de salida donde tengas permiso de escritura.

## Conflictos entre flags

| Error | Solución |
|---|---|
| `--index-vcf-gz requires --vcf-gz` | Añade `--vcf-gz`. |
| `--bcf requires --convert or --both` | Añade `--convert` o `--both`. |
| `--keep-original-info requires --convert or --both` | Añade `--convert` o `--both`. |
| `--min-strand-bias-p must be between 0 and 1` | Usa un valor de `0` a `1`. |

## Códigos de salida

| Código | Significado |
|---:|---|
| `0` | Éxito |
| `1` | Error genérico |
| `2` | Error de configuración |
| `3` | Error de validación de entrada |
| `10` | Error de lectura/escritura de archivo |
| `11` | Error de análisis de CSV/TSV |
| `12` | Error de análisis de BAM/VCF |
| `13` | Error de codificación UTF-8 |
| `14` | Error de análisis de enteros |
| `15` | Error de análisis de números en coma flotante |
| `16` | Error de JSON |
