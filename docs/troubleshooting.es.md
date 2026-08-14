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

## Salen menos filas que registros tiene el VCF

Ejemplo:

```text
Skipped 148 ALT alleles the selected sample's genotype does not carry
```

Causa: un registro VCF enumera todos los ALT vistos en ese sitio en el conjunto
de muestras, así que el `GT` de la muestra seleccionada decide cuáles le
corresponden. Un genotipo `0/0` no lleva ninguno y, en un registro multialélico,
`1/1` lleva solo el primer ALT.

Solución:

- Ninguna, si los genotipos son correctos: la ejecución anotó lo que tiene esa
  muestra.
- Elige la muestra que quieres con `--sample`, o escribe un archivo por muestra
  con `--sample all`.
- Si las llamadas no tienen un genotipo con sentido, quita el campo `GT` o ponlo
  a `./.`, que conserva todos los alelos porque desconocido no es ausente.

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

## Avisos

Un aviso nunca detiene la ejecución, así que es fácil pasarlo por alto. Estos
son los que cambian lo que acaba en la salida, y ninguno se ve en el propio TSV:
si falta una variante que esperabas, o su anotación es más gruesa de lo que
esperabas, la razón estaba en stderr.

### Falta una variante en la salida

| Aviso | Qué significa | Qué hacer |
|---|---|---|
| `falls in the phase-skipped region of CDS ... Variant skipped.` | El codón se mete en un exón vecino, y una sola fila del GFF no puede decir cuáles son las bases vecinas. La variante se descarta, no se reporta como intergénica. | Dale a las filas CDS un `transcript_id` o un `Parent` para que get_MNV pueda empalmar el transcrito, y ejecuta con `--gff-features CDS`. |

### La anotación es más gruesa de lo que permitiría el archivo

Cada uno de estos recae en la anotación por feature, así que los codones se
agrupan dentro de una fila CDS en vez de a lo largo del transcrito empalmado.
Una variante cerca de un borde de exón puede recibir entonces otro aminoácido.

| Aviso | Qué significa | Qué hacer |
|---|---|---|
| `has rows on both strands` | Dos filas de un mismo transcrito discrepan sobre la orientación, y con ellas no se puede construir un CDS empalmado. | Corrige la columna de hebra, o separa las filas en transcritos distintos. |
| `lists N duplicate CDS row(s)` | Las repeticiones se ignoran; contarlas metería las mismas bases dos veces en la secuencia codificante. | Elimina los duplicados del GFF. |
| `has a single CDS row with non-zero phase` / `starts with non-zero phase` | Las filas seleccionadas no parecen una secuencia codificante completa, así que no puede fiarse de la fase para colocar los codones. | Incluye en la ejecución todas las filas CDS del transcrito. |
| `contains CDS rows with non-zero phase, but --gff-features does not include 'CDS'` | La agrupación de codones está ignorando una fase que el archivo sí trae. | Añade `--gff-features CDS`. |

### Las métricas vienen de una muestra que no elegiste

| Aviso | Qué significa | Qué hacer |
|---|---|---|
| `Multi-sample VCF detected (N samples). Using first sample` | `ODP` y `OFREQ` se leen de la muestra que vaya primero en la cabecera, que rara vez es una elección deliberada. | Nombra una con `--sample`, o anótalas todas con `--sample all`. |

### El soporte de lecturas no es el que decía la entrada

| Aviso | Qué significa | Qué hacer |
|---|---|---|
| `has N spanning read(s) but 0 with exact CIGAR support` | Hay lecturas que cubren el sitio pero ninguna lleva ese indel exacto, que es lo que pasa cuando el alelo está escrito en otra posición dentro de un homopolímero o una repetición. | Alinea a la izquierda la entrada primero, por ejemplo con `bcftools norm -f ref.fa`. |
| `was declared X by the input VCF, but N of M reads spanning the site carry the whole haplotype` | El llamador y las lecturas discrepan sobre si esos cambios viajan juntos. | Mira el alineamiento antes de fiarte de ninguno de los dos. |
| `shows N distinct combinations on the reads; reporting the K best supported` | Una ventana con más haplotipos locales de los que get_MNV reporta. Los descartados son la cola con poco soporte. Consulta [Indels y haplotipos locales](indel-haplotypes.es.md). | Nada, salvo que te importe un haplotipo raro; entonces acota la ventana con `--chrom` e inspecciónala directamente. |

### Un flag no hizo nada

| Aviso | Qué significa |
|---|---|
| `--keep-original-info has no effect with iVar TSV input` | Un TSV de iVar no tiene columna INFO que conservar. |
| `--gff-features is ignored when using TSV annotation format (--genes)` | Los tipos de feature son una idea de GFF/GTF; el TSV de cuatro columnas no tiene ninguno. |
| `The MNV thresholds will not remove any SNP/MNV row` | Una fila SNP/MNV sobrevive cuando *cualquiera* de los dos lados supera su listón, así que con los umbrales de SNP a `0` el lado SNP siempre pasa. Sube `--snp` o `--min-snp-frequency` si querías filtrar esas filas. |

### Faltaba un programa externo opcional

| Aviso | Qué significa |
|---|---|
| `bcftools not found in PATH. Skipping BCF output.` | Se pidió `--bcf` y no se escribió ningún BCF. La ejecución continúa y sale con `0`; el TSV y el VCF no se ven afectados, y el resumen JSON no reporta ningún BCF. Instala bcftools, o quita `--bcf`. |
| `tabix not found in PATH. Skipping .tbi index` | Se pidió `--index-vcf-gz` y no se escribió ningún `.tbi`. El propio `.vcf.gz` está completo. Instala samtools/htslib, o indéxalo después con `tabix -p vcf <fichero>`. |

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
