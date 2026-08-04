# End-to-end scenario tests

[English](README.md) · **Español**

Este directorio contiene un framework de pruebas en Python que construye
entradas sintéticas de FASTA, GFF, VCF (o TSV de iVar) y BAM a partir de
definiciones declarativas de escenarios, ejecuta `get_mnv` sobre cada una y
coteja el TSV de salida con los campos esperados predefinidos.

Es complementario a las pruebas unitarias y de integración de Rust de
[`tests/integration.rs`](../integration.rs). Las pruebas de Rust ejercitan la
biblioteca a nivel de función; estas pruebas de escenarios ejercitan el binario
completo sobre entradas realistas, lo que incluye:

- Soporte de lecturas calculado a partir del CIGAR del BAM
- Agrupación de codones a través de las uniones entre exones (modelos de CDS multiexón)
- Anotación de genes en hebra negativa
- Entrada TSV de iVar con notación `+SEQ` / `-SEQ`
- Entradas multicontig

## Quick start

Requisitos:

- Python 3.10 o posterior (solo se usa la biblioteca estándar)
- `samtools` en el `PATH` (o define la variable de entorno `SAMTOOLS`)
- Un binario `get_mnv` ya compilado: por defecto la suite busca
  `target/debug/get_mnv`, luego `target/release/get_mnv`, luego
  `dist/get_mnv` y, por último, `get_mnv` en el `PATH`. Cambia esta ruta con `GET_MNV=/path/to/get_mnv`.

```bash
# Compila get_mnv primero
cargo build

# Ejecuta los 36 escenarios
python3 tests/scenarios/run.py

# Ejecuta un subconjunto por prefijo de nombre
python3 tests/scenarios/run.py 18 22 27

# Usa un binario release en lugar de debug
GET_MNV=$(pwd)/target/release/get_mnv python3 tests/scenarios/run.py

# Usa el samtools de un entorno conda
SAMTOOLS=/path/to/conda/envs/bio/bin/samtools python3 tests/scenarios/run.py
```

Los archivos intermedios (FASTA, GFF, VCF, BAM, la salida de `get_mnv` y el log)
de cada escenario se escriben en `tests/scenarios/work/<scenario_name>/` y se
sobrescriben en cada ejecución. El directorio está en `.gitignore`.

## Visual overview (PDF)

[`plot_scenarios.py`](plot_scenarios.py) renderiza cada escenario en un único
PDF ilustrativo (`scenarios_overview.pdf`). Necesita `matplotlib` además de
`samtools`.

```bash
SAMTOOLS=/path/to/samtools python3 tests/scenarios/plot_scenarios.py
```

Las páginas se agrupan en **secciones** temáticas (SNV/MNV, clasificación de
indels, indels combinados con SNVs, hebra negativa, transcritos eucariotas
multiexón, formatos de entrada y casos límite), precedidas por una portada, una
leyenda de glifos y un índice. Cada página de escenario muestra:

- las **lecturas mapeadas** como un pequeño pileup de alineamiento, donde se
  marcan los SNVs, las inserciones, las deleciones y los saltos de intrón (CIGAR
  `N`), los genes/CDS solapantes y las variantes de entrada;
- una **pista de marco de lectura y aminoácidos**: tripletes de codón sombreados
  de forma alterna, junto al aminoácido de referencia de cada codón y con el AA
  alternativo en rojo cuando una variante vuelve no sinónimo al codón (teniendo
  en cuenta la hebra y la fase);
- dos tablas: **"Sin get_mnv"** (las llamadas en crudo del caller: posición,
  alelo y conjetura del marco a partir de la longitud) frente a la verdadera
  **"salida de get_mnv"**, con un resumen de una línea, **"Lo que añade
  get_mnv"**, entre ambas.

El PDF generado y cualquier vista previa `plot_*.png` están en `.gitignore`;
regenéralos desde el script.

## Mini-genome layout

Todos los escenarios comparten una referencia sintética (`framework.REFERENCE_SEQ`)
compuesta por dos contigs diseñados para mantener explícita la aritmética de codones:

### `chr_test` (1300 bp)

| Range       | Feature        | Content                                                   |
|-------------|----------------|-----------------------------------------------------------|
| 1–300       | `geneA` (+)    | `ATG` + 98 × `GCT` + `TAA` — Met, 98×Ala, stop           |
| 301–400     | filler         | `A` × 100                                                 |
| 401–700     | `geneB` (−)    | RC of the standard CDS: `TTA` + 98 × `AGC` + `CAT`        |
| 701–800     | filler         | `A` × 100                                                 |
| 801–900     | `geneC` exon 1 | `ATG` + 32 × `GCT` + `G` (codons 1–33 + base 1 of codon 34) |
| 901–1000    | `geneC` intron | `T` × 100                                                 |
| 1001–1200   | `geneC` exon 2 | `CT` + 65 × `GCT` + `TAA` (bases 2–3 of codon 34, codons 35–100) |
| 1201–1300   | filler         | `A` × 100                                                 |

### `chr_test2` (600 bp)

| Range  | Feature     | Content                                |
|--------|-------------|----------------------------------------|
| 1–300  | `geneD` (+) | `ATG` + 98 × `GCT` + `TAA` (same as geneA) |
| 301–600 | filler     | `A` × 300                              |

Resumen de la aritmética de codones:

- `geneA` codón N (+): posiciones `(3N-2, 3N-1, 3N)`
- `geneB` codón N (−): posiciones `(701-3N, 702-3N, 703-3N)`; base del mRNA = RC del genómico
- `geneC` codón 34 abarca la unión exón 1/intrón/exón 2: pos 900 (exón 1) + pos 1001-1002 (exón 2)
- Las columnas de codón se reportan en la orientación del transcrito, así que una expectativa de hebra menos es el codón del mRNA (`GCT`), no las bases genómicas (`AGC`)

## Validated scenarios

Actualmente hay 36 escenarios definidos en [`scenarios.py`](scenarios.py).
Cada uno declara las variantes de entrada, los grupos de lecturas del BAM (con
operaciones opcionales: sustitución SNV, inserción, deleción y salto de intrón) y
las filas TSV esperadas.

### Core SNP / MNV / SNP/MNV

| # | Name | What it validates |
|---|------|-------------------|
| 01 | `snp_simple` | Un único SNV `Ala10Thr` (geneA codón 10 GCT→ACT) con soporte de 20/20 lecturas |
| 02a | `snp_mnv_full_phasing` | Dos SNVs en el mismo codón y todas las lecturas portan ambos → `SNP/MNV`, `Ala10Ser`, 20 lecturas MNV, 0 lecturas solo-SNP |
| 02b | `vcf_mnp_decomposed` | Un registro MNP del VCF (`REF=GC ALT=TA`) se descompone en dos componentes SNV y se reporta como `SNP/MNV` con `event_class=mnv` |
| 03 | `snp_mnv_mixed` | Fase parcial: 10/10/10 lecturas (ambos / solo SNV1 / solo SNV2) → frecuencias 0.3333 / 0.3333 / 0.3333 |

### In-frame and frameshift indels

| # | Name | What it validates |
|---|------|-------------------|
| 04 | `ins_inframe_cds` | Inserción en marco de `GCT` (1 Ala) dentro del CDS, 20 lecturas, `INS:30:+GCT` |
| 05 | `del_frameshift_cds` | Deleción de 1 bp con frameshift, `Ala11Leufs`, `DEL:31:G` |
| 24 | `large_inframe_insertion` | Inserción en marco de 9 bp (`+GCTGCTGCT` = 3 Ala) |
| 25 | `large_inframe_deletion` | Deleción en marco de 6 bp (2 Ala consecutivas), `DEL:31-36:GCTGCT` |

### Complex haplotypes (indel + SNV/MNV in cis)

| # | Name | What it validates |
|---|------|-------------------|
| 06 | `indel_plus_snv_haplotype` | Indel + SNV a más de 3 bp de distancia: NO se fusionan en un `complex_indel` (quedan fuera de la ventana local de 3 bp) |
| 07 | `fs_del_plus_downstream_snv` | El frameshift se propaga a un SNV aguas abajo → sufijo `(fs)` en el AA Change y Change Type `Synonymous (frameshift)` |
| 08 | `inframe_ins_inside_codon_with_mnv` | Inserción en marco dentro de un codón + MNV en el mismo codón: emite 5 filas, que incluyen 2 `complex_indel`, la fila del MNV marcada como `Indel overlap` / `Unknown`, la inserción sola y el `complex_indel` ins+SNV |
| 09 | `fs_del_with_snv_overlap` | Deleción con frameshift que solapa el codón + SNV: la fila del SNV obtiene `Indel overlap` / `Unknown`; la fila de la deleción sola tiene `Event Reads = 0` (ninguna lectura porta solo la deleción); el `complex_indel` lleva el frameshift `Ala10Cysfs` |

### Negative-strand gene (`geneB`)

| # | Name | What it validates |
|---|------|-------------------|
| 10 | `minus_snp_simple` | Un único SNV en `geneB` (pos 673 C>T) — el codón 10 del mRNA `GCT`→`ACT` → `Ala10Thr` |
| 11 | `minus_mnv_same_codon` | MNV en el codón 10 de `geneB` (pos 671 A>T + pos 673 C>A en cis) — mRNA `GCT`→`TCA` → `Ala10Ser` |
| 12 | `minus_fs_del` | Deleción de 1 bp con frameshift dentro de `geneB` |

### Multi-exon CDS (`geneC`, with `--gff-features CDS`)

| # | Name | What it validates |
|---|------|-------------------|
| 13 | `multiexon_snp_exon2` | SNV en el exón 2 — confirma que la posición del codón se resuelve teniendo en cuenta el transcrito |
| 14 | `multiexon_junction_snp` | SNV en el codón que abarca la unión (base 1 en el exón 1, pos 900) — requiere el modelo empalmado |
| 15 | `multiexon_junction_mnv` | MNV con una base en el exón 1 (pos 900) y otra en el exón 2 (pos 1002) — adyacentes en el mRNA empalmado, pero a unos 100 bp de distancia en el genoma. Se valida mediante lecturas con N-CIGAR que abarcan el intrón |

### Operational edge cases

| # | Name | What it validates |
|---|------|-------------------|
| 16 | `no_bam_coverage` | El VCF declara un SNV en una posición donde el BAM no tiene lecturas solapantes — `get_mnv` aún emite la fila, con las columnas de soporte de lecturas vacías |
| 17 | `min_snp_frequency_filter` | Un SNV presente al 10% (2/20 lecturas) se descarta de la salida cuando se establece `--min-snp-frequency 0.5` |

### Amino-acid edge cases

| # | Name | What it validates |
|---|------|-------------------|
| 18 | `stop_gained_via_mnv` | MNV de tres SNVs en el codón 50 `GCT`→`TAA` → `Ala50Ter`, `Change Type = Stop gained` |
| 19 | `start_codon_altered` | SNV en el codón de inicio ATG en la pos 2 → `Met1Thr`. NOTA: `get_mnv` no tiene un Change Type `Start lost` dedicado; la fila se reporta como `Non-synonymous` |
| 20 | `stop_lost` | SNV en el codón de stop `TAA` (pos 298) → `Change Type = Stop lost` |

### Complex alleles

| # | Name | What it validates |
|---|------|-------------------|
| 21 | `intron_variant` | SNV en el intrón de geneC (pos 950), lejos de los sitios de splice → `intron_variant` contra geneC, no `intergenic` |
| 22 | `multiallelic_split` | VCF multialélico `pos 28 REF=G ALT=A,T` con `--split-multiallelic`: cada ALT produce ahora una fila de anotación independiente (prueba de regresión para la corrección del commit 9ea2aed) |
| 23 | `delins_subst_plus_del` | `REF=GCT ALT=GA`: un alelo compuesto de sustitución + deleción de 1 bp produce una fila INDEL con frameshift |

### Multi-contig and iVar TSV input

| # | Name | What it validates |
|---|------|-------------------|
| 26 | `multicontig` | Dos variantes en dos contigs distintos (`chr_test` y `chr_test2`) producen dos filas correctamente anotadas |
| 27 | `ivar_tsv_snv` | Entrada TSV de iVar con un SNV simple — la misma fila anotada que su equivalente en VCF |
| 28 | `ivar_tsv_insertion` | La notación `+GCT` del TSV de iVar se normaliza al alelo anclado estilo VCF y se anota como `INS:30:+GCT` |
| 29 | `ivar_tsv_deletion` | La notación `-G` del TSV de iVar se normaliza a `DEL:31:G` y produce la misma fila con frameshift |

### Indel refinements (indels branch)

Estos ejercitan los refinamientos en el manejo de indels que se añadieron en la
rama `indels` (commit 49d2f09): la detección de stop para indels en marco y el
control de propagación aguas abajo `--frameshift-min-freq`.

| # | Name | What it validates |
|---|------|-------------------|
| 30 | `stop_gained_inframe_ins` | Una inserción en marco de `TAA` tras la pos 30 forma un codón de stop prematuro → `Change Type = Stop gained`, `AA Changes = Ala10_Ala11ins*` (en lugar del genérico `In-frame Indel`). Lo gobierna `indel_stop_effect`, que compara el número de residuos de stop en la proteína ref frente a la alt |
| 31 | `stop_lost_inframe_del` | Deleción en marco de 3 bp del stop terminal `TAA` (pos 298-300, `DEL:298-300:TAA`) → `Change Type = Stop lost`, `AA Changes = *100del` |
| 32 | `fs_gate_default_propagates` | Deleción con frameshift aguas arriba de baja frecuencia (`AF=0.20`) + un SNV aguas abajo. Con el valor por defecto `--frameshift-min-freq 0.0` el frameshift se propaga → el SNV aguas abajo se etiqueta como `Synonymous (frameshift)` / `Ala13Ala (fs)` |
| 33 | `fs_gate_high_freq_suppressed` | Entradas idénticas a las del escenario 32, pero con `--frameshift-min-freq 0.5`. La deleción aguas arriba (`AF=0.20`) no supera el control, por lo que el frameshift **no** se propaga → el SNV aguas abajo es un simple `Synonymous` / `Ala13Ala` |

Los escenarios 32 y 33 forman un par A/B: el mismo VCF (que ahora lleva `AF` en
el `INFO`) y el mismo BAM, que se diferencian solo en el flag
`--frameshift-min-freq`. Juntos muestran que el control afecta únicamente a la
propagación aguas abajo del frameshift — la fila propia de la deleción aguas
arriba permanece como `Frameshift Indel` (`Ala11Leufs`) en ambos casos.

## Framework architecture

- [`framework.py`](framework.py) define los bloques de construcción:
  - `REFERENCE_SEQ`, `REFERENCE_SEQ2`, `GFF_GENE_ONLY`, `GFF_CDS_MULTIEXON`
  - Clases de datos: `Op`, `ReadGroup`, `VcfRecord`, `IvarRecord`, `ExpectedRow`, `Scenario`
  - Constructores: `write_fasta`, `write_gff`, `write_vcf`, `write_ivar_tsv`, `write_bam`
  - Driver: `run_get_mnv`, `parse_tsv`, `compare`, `run_scenario`
  - `VcfRecord(af=...)` emite un campo INFO `AF=` (omitido por defecto) para que
    los escenarios puedan ejercitar lógica condicionada por la frecuencia, como
    `--frameshift-min-freq`
  - `Scenario.extra_cli_args` pasa flags adicionales a `get_mnv`
- [`scenarios.py`](scenarios.py) declara los 36 escenarios.
- [`run.py`](run.py) es el driver de la CLI:
  ```bash
  python3 run.py             # all scenarios
  python3 run.py 18 22       # only scenarios starting with "18" or "22"
  ```

## Adding a new scenario

1. Elige un prefijo numérico único y un nombre descriptivo.
2. Decide en qué punto del minigenoma cae la variante (consulta la aritmética de
   codones de más arriba). Posiciones habituales:
   - `geneA` codón 10: pos 28-30
   - `geneA` codón 50: pos 148-150
   - `geneB` codón 10 (hebra −): pos 671-673
   - `geneC` codón de unión 34: pos 900 + pos 1001-1002
3. Construye los grupos de lecturas con las operaciones necesarias para dar
   soporte a la variante en el BAM. Los tipos de `Op` son:
   - `Op("snv", pos=P, seq="X")` — sustituye la base en la pos de ref `P` por `X`
   - `Op("ins", pos=P, seq="ABC")` — inserta `ABC` tras la pos de ref `P`
   - `Op("del", pos=P, length=N)` — borra `N` bases a partir de `P`
   - `Op("skip", pos=P, length=N)` — emite un CIGAR `N` (salto de intrón), para
     las lecturas empalmadas de los escenarios multiexón
4. Declara las filas TSV esperadas con los campos que quieras verificar
   (`positions`, `gene`, `base_changes`, `aa_changes`, `variant_type`,
   `change_type`, `event_class`, `event_components`, `*_reads`,
   `*_frequencies`, etc.). Los campos sin establecer no se comprueban.
5. Si quieres, establece `expected_row_count` para afirmar el número total
   exacto de filas producidas — resulta útil cuando la ventana de haplotipo
   local emite más de una fila (consulta los escenarios 08 y 09).
6. Añade el escenario a `ALL_SCENARIOS` y ejecuta `python3 run.py`.

El driver mostrará `PASS` o `FAIL` para cada escenario y, cuando falle,
imprimirá las filas reales producidas para que puedas ajustar la expectativa.

## Known behaviours documented by these tests

- `geneB` (hebra negativa) — las columnas `Reference Codon` y `MNV Codon`
  muestran la secuencia **genómica directa** (p. ej. `AGC`), no la del mRNA
  (`GCT`). El efecto en el AA sigue siendo correcto.
- Un CDS multiexón sin un atributo `Name=` en la feature de CDS reporta la
  columna del gen como el primer ID de CDS (p. ej. `cds-geneC-e1`). Los GFF
  estándar de NCBI/Ensembl incluyen `Name=`, así que esto solo afecta a los
  GFF hechos a mano.
- Un SNV/MNV que altera la Met iniciadora (posición proteica 1) se reporta con
  Change Type `Start lost` (escenario 19); `Met1` → stop permanece como
  `Stop gained`, y las sustituciones internas de Met no se ven afectadas.
- `--split-multiallelic`: cada ALT en la misma posición de codón emite ahora una
  fila de anotación independiente.
- Los indels en marco que **crean o eliminan un codón de stop** se reportan como
  `Stop gained` / `Stop lost` en lugar del genérico `In-frame Indel`. Los indels
  con frameshift conservan la etiqueta `Frameshift Indel` (casi siempre
  introducen un stop aguas abajo, por lo que la distinción no aportaría información).
- `--frameshift-min-freq F` controla la **propagación aguas abajo del frameshift**:
  un indel aguas arriba solo desplaza el marco de los codones aguas abajo si su
  frecuencia alélica reportada (`AF`/`FREQ`/`AD` del VCF, `ALT_FREQ` de iVar) es
  ≥ `F`. El valor por defecto `0.0` reproduce el comportamiento histórico de
  propagar siempre, y los indels sin una frecuencia conocida se propagan siempre.
  El control nunca cambia la clasificación propia del indel.
- Cuando no se proporciona `--gff-features` y el GFF contiene features `CDS`, se
  selecciona automáticamente el modelo de CDS que tiene en cuenta la fase y el
  empalme (escenario 34); pasa `--gff-features gene` para conservar el
  comportamiento heredado de gen completo.
- Una variante situada aguas abajo de un stop prematuro introducido por un
  frameshift se reporta con Change Type `Downstream of premature stop`, en lugar
  de llevar una anotación `(fs)` como si estuviera traducida (escenario 35). La
  propagación ordinaria del frameshift, en la que no se introduce ningún stop
  temprano, no cambia.
- La ventana de haplotipo local es de 3 bp (`LOCAL_HAPLOTYPE_JOIN_DISTANCE`) y
  acepta hasta 8 eventos (`MAX_LOCAL_HAPLOTYPE_VARIANTS`). Los eventos más
  alejados producen filas separadas en lugar de fusionarse en un `complex_indel`.
