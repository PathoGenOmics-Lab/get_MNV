# Referencia de la CLI

Referencia completa de las opciones de línea de comandos de `get_mnv` (versión 1.1.5). Ejecuta
`get_mnv --help` para ver la misma lista en tu terminal.

## Sinopsis

```text
get_mnv [OPTIONS] --fasta <FASTA_FILE> <--vcf <VCF_FILE>|--tsv <TSV_FILE>>
```

Debes proporcionar una referencia (`--fasta`) y exactamente una fuente de variantes
(`--vcf` o `--tsv`), además de una anotación de genes (`--gff` o `--genes`). La
única excepción es `--report-from`, que construye un report a partir de salidas
TSV ya existentes y no necesita ninguna de las tres.

## Entrada

| Opción | Descripción |
|---|---|
| `-v, --vcf <FILE>` | Entrada de variantes en VCF plano o comprimido con BGZF (SNV/MNV e indels). |
| `--tsv <FILE>` | Entrada `variants.tsv` de iVar. |
| `-b, --bam <FILE>` | Lecturas alineadas opcionales para soporte de lecturas. Debe estar ordenado por coordenadas e indexado. |
| `-f, --fasta <FILE>` | FASTA de referencia (requerido). |
| `--sample <NAME>` | Muestra para las métricas FORMAT originales en un VCF multimuestra (por defecto: primera muestra; `all` para todas). |
| `--chrom <NAME>` | Restringe el procesamiento a un solo contig (por defecto: todos los contigs de la entrada). |

## Anotación

| Opción | Descripción |
|---|---|
| `--gff <FILE>` | Anotación de genes en formato GFF/GFF3. |
| `-g, --genes <FILE>` | Tabla de genes TSV simple: `gene,start,end,strand`. Úsala en lugar de `--gff`. |
| `--gff-features <LIST>` | Tipos de característica GFF separados por comas a analizar (por defecto: `gene,pseudogene`). Usa `CDS` para transcritos con splicing. |
| `--translation-table <N>` | Código genético NCBI (por defecto: `11`, bacteriano). Admitidos: 1, 2, 3, 4, 5, 6, 11, 12, 25. |
| `--exclude-intergenic` | Descarta las variantes fuera de los genes anotados. |

## Soporte de lecturas y calidad

Estas solo aplican cuando se proporciona `--bam`.

| Opción | Por defecto | Descripción |
|---|---|---|
| `-q, --quality <N>` | `20` | Calidad Phred de base mínima. |
| `--min-mapq <N>` | `0` | Calidad de mapeo mínima (MAPQ). |
| `-s, --snp <N>` | `0` | Lecturas mínimas que soportan el SNP. |
| `--min-snp-frequency <F>` | `0.0` | Frecuencia alélica mínima del SNP derivada del BAM (0.0–1.0). |
| `-m, --mnv <N>` | `0` | Lecturas mínimas que soportan el MNV. |
| `--min-mnv-frequency <F>` | `0.0` | Frecuencia mínima del haplotipo MNV derivada del BAM (0.0–1.0). |
| `--min-snp-strand <N>` | `0` | Lecturas mínimas que soportan el SNP en cada hebra. |
| `--min-mnv-strand <N>` | `0` | Lecturas mínimas que soportan el MNV en cada hebra. |
| `--min-strand-bias-p <F>` | `0.0` | Valor p mínimo del test exacto de Fisher aceptado para las métricas de sesgo de hebra. **Solo afecta a la salida VCF**: el escritor de TSV no tiene umbral de sesgo, así que esto nunca quita una fila del TSV. |

!!! note "Cómo se combinan los umbrales de SNP y MNV"
    Los filtros de frecuencia y recuento de lecturas usan el soporte recalculado
    a partir de `--bam`, no los `OFREQ`/`ODP` originales de la entrada.

    Una fila a nivel de codón (`SNP/MNV`) se conserva cuando **cualquiera** de
    los dos lados supera su listón: sus SNV individuales pasan los umbrales de
    SNP, **o** su haplotipo pasa los de MNV. Eso conserva un haplotipo bien
    soportado cuyos SNV individuales son débiles, y funciona también en el otro
    sentido: con los umbrales de SNP en su valor por defecto de `0`, el lado SNP
    pasa siempre, así que subir `--mnv` a solas no quita nada. Sube los dos, o
    ninguno. La ejecución avisa cuando has pedido un umbral MNV que no puede
    llegar a una fila a nivel de codón.

    Qué umbrales gobiernan cada fila:

    | Fila | La juzgan |
    |---|---|
    | `SNP` | los umbrales de SNP |
    | `MNV`, `SNP/MNV` | cualquiera de los dos lados, como arriba |
    | `INDEL` | los umbrales de **MNV**, medidos contra el soporte de evento del propio indel (`Event Reads`, `Event Forward/Reverse Reads`, `Event Depth`), no contra ninguna columna SNP |

    Así que un indel se filtra con `--mnv`, `--min-mnv-frequency` y
    `--min-mnv-strand`; `--snp` no lo toca nunca. Un indel intergénico no lleva
    soporte recalculado, así que no puede satisfacer un filtro basado en lecturas
    y se descarta en cuanto hay algún umbral MNV activo.

## Ajuste de indels

| Opción | Por defecto | Descripción |
|---|---|---|
| `--frameshift-min-freq <F>` | `0.5` | Frecuencia mínima que debe alcanzar un indel aguas arriba para marcar con frameshift los codones de SNV/MNV aguas abajo. Por defecto solo propaga desde un indel aguas arriba mayoritario; pon `0.0` para propagar desde cualquiera. Con `--bam` la frecuencia es la que get_MNV cuenta de las lecturas (el mismo número que reporta como `EFREQ`), no el `AF` que declaró el llamador; sin BAM se cae al `AF` declarado. Un indel cuya frecuencia no se conoce por ninguna de las dos vías siempre propaga, porque no hay nada que comparar. |
| `--legacy-indel-depth` | desactivado | Restringe la profundidad del locus del indel (el denominador de `EFREQ`) a las lecturas que abarcan el alelo REF entero. Por defecto se cuenta desde las que observan la base de anclaje, lo que evita subcontar profundidad en deleciones multibase; este flag recupera el denominador antiguo, más estrecho. |
| `--phased-indel-min-reads <N>` | `2` | Lecturas mínimas que soportan en el BAM para emitir una fila de haplotipo indel/complejo en fase. Una sola lectura no es evidencia de un haplotipo. |
| `--count-mates-separately` | apagado | Cuenta los dos mates de un fragmento como dos observaciones en vez de una molécula. |
| `--phased-indel-min-freq <F>` | `0.0` | Frecuencia mínima derivada del BAM para emitir una fila de haplotipo indel/complejo en fase. |
| `--normalize-alleles` | off | Recorta el prefijo/sufijo compartido de REF/ALT antes del procesamiento. |
| `--split-multiallelic` | off | Divide los registros VCF multialélicos en alelos ALT independientes en lugar de fallar. |

## Salida

| Opción | Descripción |
|---|---|
| `--convert` | Escribe salida VCF (`.MNV.vcf`) en lugar de TSV. |
| `--both` | Escribe TSV y VCF en una sola ejecución. |
| `--vcf-gz` | Escribe `.vcf.gz` comprimido con BGZF (modo de salida VCF). |
| `--index-vcf-gz` | Construye un índice Tabix `.tbi` (requiere `--vcf-gz`). |
| `--bcf` | Escribe también un BCF convertido a partir del VCF generado (requiere `--convert`/`--both`). |
| `--strand-bias-info` | Añade los valores p de sesgo de hebra del test exacto de Fisher al INFO del VCF (`SBP`/`MSBP`). |
| `--keep-original-info` | Conserva los campos INFO originales del VCF en la salida (requiere `--convert`/`--both`). |
| `--emit-filtered` | Emite los registros que no superan los umbrales con etiquetas `FILTER` en lugar de omitirlos. |

## Informe HTML interactivo

| Opción | Descripción |
|---|---|
| `--report <HTML_FILE>` | Escribe un informe HTML interactivo y autocontenido de las variantes llamadas. Necesita la salida TSV, que es la de por defecto. `--convert` escribe el VCF *en lugar* del TSV, y los dos flags son mutuamente excluyentes, así que usa `--both` en vez de `--convert` cuando quieras un informe junto a la salida VCF. Con `--sample all` el informe cubre todas las muestras. |
| `--report-from <TSV>...` | Construye el informe a partir de TSV de get_MNV ya existentes, sin ejecutar el pipeline, para cohortes procesadas muestra a muestra. Cada fichero es una muestra, etiquetada con su nombre de archivo. Requiere `--report` para la ruta de salida. |

Consulta [Formatos de salida](output-formats.es.md#informe-html-interactivo) para saber qué contiene el informe.

## Validación y metadatos

| Opción | Descripción |
|---|---|
| `--dry-run` | Valida las entradas e imprime un resumen por contig sin escribir salidas. |
| `--strict` | Falla si faltan en la entrada las métricas originales de profundidad/frecuencia (`ODP`/`OFREQ`). |
| `--summary-json <FILE>` | Escribe un resumen de la ejecución legible por máquina. |
| `--run-manifest <FILE>` | Escribe un manifiesto de reproducibilidad (entradas, salidas, checksums, metadatos de ejecución). |
| `--error-json <FILE>` | Escribe detalles de error estructurados en formato JSON cuando el comando falla. |
| `--threads <N>` | Número de hilos de trabajo (por defecto: automático de Rayon). |
| `-h, --help` | Imprime la ayuda. |
| `-V, --version` | Imprime la versión. |
