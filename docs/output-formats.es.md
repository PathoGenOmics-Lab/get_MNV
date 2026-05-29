# Formatos de salida

get_MNV puede escribir archivos de metadatos TSV, VCF, BCF y JSON.

## Salida TSV por defecto

Nombre de archivo por defecto:

```text
<input_name>.MNV.tsv
```

Usa este formato para hojas de cálculo, análisis posteriores e inspección rápida.

Columnas principales:

| Columna | Significado |
|---|---|
| `Chromosome` | Nombre del contig |
| `Gene` | Nombre del gen o de la característica. Las variantes intergénicas se marcan como `intergenic`. |
| `Positions` | Una posición para los SNP, varias posiciones separadas por comas para los MNV. |
| `Reference Bases` | Bases de referencia en esas posiciones. |
| `Base Changes` | Bases alternativas. |
| `AA Changes` | Cambio de aminoácido tras combinar todos los SNV del codón. |
| `SNP AA Changes` | Cambio de aminoácido para cada SNV considerado por separado. |
| `Local AA Changes` | Numeración por característica para los modelos de característica heredados; idéntica a `AA Changes` cuando hay disponible un modelo de CDS de transcrito empalmado. |
| `Local SNP AA Changes` | Cambios de aminoácido por SNP en numeración local. |
| `Variant Type` | `SNP`, `MNV`, `SNP/MNV` o `INDEL`. |
| `Change Type` | Sinónimo, no sinónimo, codón de parada ganado/perdido, desconocido, etc. |
| `Reference Codon` | Codón original. |
| `SNP Codon` | Codón con sustituciones de SNP individuales. |
| `MNV Codon` | Codón con todas las sustituciones agrupadas. |
| `Event Class` | Clase canónica del evento de alelo: `snp`, `mnv`, `insertion`, `deletion`, `delins`, `complex_indel` o `symbolic`. |
| `Event Components` | Descomposición REF/ALT como `SNV:10:A>G`, `INS:10:+T` o `DEL:11-12:TG`. |

Columnas adicionales cuando se usa `--bam`:

| Columna | Significado |
|---|---|
| `SNP Reads` | Lecturas que respaldan cada SNV individual. |
| `SNP Forward/Reverse Reads` | Soporte de SNP específico de la hebra. |
| `MNV Reads` | Lecturas que respaldan el haplotipo MNV completo. |
| `MNV Forward/Reverse Reads` | Soporte de MNV específico de la hebra. |
| `Total Reads` | Profundidad en las posiciones de la variante. |
| `SNP Frequencies` | Frecuencias de SNP por posición. |
| `MNV Frequencies` | Frecuencia del haplotipo MNV. |
| `Event Reads` | Lecturas exactas que respaldan un evento indel/complejo. |
| `Event Forward/Reverse Reads` | Soporte exacto del evento específico de la hebra. |
| `Event Depth` | Lecturas con un alelo observado a lo largo del tramo del evento indel/complejo. |
| `Event Frequency` | Lecturas exactas del evento divididas por la profundidad del evento. |

El soporte exacto del evento tiene en cuenta el CIGAR. Una lectura debe reconstruir la misma
secuencia ALT local y, para los haplotipos complejos, contener los componentes de inserción
y deleción esperados. Esto evita que las combinaciones de inserción/deleción netamente neutras
se cuenten como soporte solo porque su secuencia parezca un MNV.

Las columnas de frecuencia se calculan a partir del soporte del BAM. `--min-snp-frequency` y
`--min-mnv-frequency` usan estos mismos valores derivados del BAM. Los filtros son
independientes: `--min-snp-frequency` se aplica a las observaciones de SNP individuales, y
`--min-mnv-frequency` se aplica a los haplotipos MNV fasados. En las llamadas mixtas `SNP/MNV`,
una fila o registro VCF se conserva cuando cualquiera de los componentes supera su propio
umbral activo.
Los filtros de recuento de lecturas y de soporte de hebra (`--snp`, `--mnv`, `--min-snp-strand`
y `--min-mnv-strand`) siguen el mismo comportamiento independiente de SNP/MNV.

Cuando un MNV a nivel de codón se solapa con un indel, la fila del MNV se conserva como una fila
de contexto posicional, pero su efecto a nivel de aminoácido se marca como `Unknown` con
`Change Type = Indel overlap`. Si las lecturas del BAM respaldan el evento combinado completo,
get_MNV emite una fila `complex_indel` exacta aparte con el REF/ALT combinado, los
componentes del evento y el soporte de lecturas del evento.

El solapamiento de indel sigue la semántica interbase de VCF. Las deleciones se solapan con una
característica por su tramo de referencia eliminado. Las inserciones se solapan con una característica
solo cuando la secuencia insertada cae entre dos bases de referencia dentro de esa característica, de
modo que una inserción anclada en la última base de la característica se reporta fuera de esa característica.

Ejemplo:

```text
Chromosome	Gene	Positions	Base Changes	AA Changes	Variant Type	Change Type
MTB_anc	Rv0095c	104838	T	Asp126Glu	SNP	Non-synonymous
MTB_anc	Rv0095c	104941,104942	T,G	Gly92Gln	SNP/MNV	Non-synonymous
```

## Salida VCF

Escribe VCF con:

```bash
--convert
```

o escribe tanto TSV como VCF con:

```bash
--both
```

Nombre de archivo por defecto:

```text
<input_name>.MNV.vcf
```

Usa `--vcf-gz` para salida comprimida:

```text
<input_name>.MNV.vcf.gz
```

Campos INFO comunes:

| Campo | Significado |
|---|---|
| `GENE` | Nombre del gen o de la característica |
| `AA` | Cambio de aminoácido |
| `CT` | Tipo de cambio |
| `TYPE` | Tipo de variante |
| `EC` | Clase canónica del evento de alelo |
| `COMP` | Componentes del evento REF/ALT |
| `ODP` | Profundidad original del archivo de variantes de entrada |
| `OFREQ` | Frecuencia alélica original del archivo de variantes de entrada |
| `SR`, `SRF`, `SRR` | Lecturas de SNP: total, forward, reverse |
| `MR`, `MRF`, `MRR` | Lecturas de MNV: total, forward, reverse |
| `DP` | Profundidad recalculada a partir del BAM |
| `FREQ` | Frecuencia recalculada a partir del BAM |
| `ER`, `ERF`, `ERR` | Lecturas exactas del evento indel/complejo: total, forward, reverse |
| `EDP` | Profundidad exacta del evento para alelos indel/complejos |
| `EFREQ` | Frecuencia exacta del evento para alelos indel/complejos |
| `SBP` | Valor p del sesgo de hebra de SNP |
| `MSBP` | Valor p del sesgo de hebra de MNV |

La cabecera del VCF registra la versión de get_MNV, la línea de comandos y los umbrales usados.
Cuando `--emit-filtered` está habilitado, los registros VCF por debajo de los umbrales de soporte de
lecturas, frecuencia, soporte de hebra o sesgo de hebra se escriben con etiquetas FILTER como
`LowSupport`, `LowFrequency`, `StrandSupport` o `StrandBias`; de lo contrario, se omiten.

## Salida BCF

Escribe BCF con:

```bash
--bcf
```

BCF requiere el modo de salida VCF, así que úsalo con `--convert` o `--both`.
Esto es solo conversión de salida; BCF no se acepta como formato de entrada.

Nombre de archivo por defecto:

```text
<input_name>.MNV.bcf
```

## Archivos JSON

### JSON de resumen

Escribe con:

```bash
--summary-json run.summary.json
```

Incluye:

- Sumas de comprobación de los archivos de entrada
- Recuentos de variantes por contig
- Recuentos globales de variantes
- Tiempos de ejecución
- Rutas de salida

### Manifiesto de la ejecución

Escribe con:

```bash
--run-manifest run.manifest.json
```

Incluye el resumen más:

- Línea de comandos
- Versión de la herramienta
- Sumas de comprobación de los archivos de salida
- Marca de tiempo

### JSON de errores

Escribe los errores como JSON con:

```bash
--error-json run.error.json
```

Esto es útil en pipelines automatizados.

## Notas

- Para los registros MNV, la profundidad y la frecuencia se calculan a partir de las lecturas que
  abarcan todas las posiciones del haplotipo agrupado.
- Las frecuencias se imprimen con 4 decimales.
- `--min-snp-frequency` y `--min-mnv-frequency` son valores de `0` a `1`
  y requieren `--bam`.
- Los filtros de frecuencia de SNP y MNV son independientes, de modo que un haplotipo MNV fuerte
  no se elimina por un umbral de frecuencia de SNP más estricto.
- Los filtros de soporte de lecturas y de soporte de hebra de SNP y MNV también son independientes.
- `--sample all` escribe un conjunto de salida por cada muestra del VCF.
- `--keep-original-info` conserva los campos INFO que no son de get_MNV del VCF de entrada.
