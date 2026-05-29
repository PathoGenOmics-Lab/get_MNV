<p align="center">
  <img src="images/get_mnv.png" alt="logo de get_MNV" width="800" />
</p>

<div align="center">

[![License: AGPL v3](https://img.shields.io/badge/license-AGPL%20v3-%23af64d1?style=flat-square)](LICENSE)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/get_mnv.svg?style=flat-square&label=bioconda)](https://anaconda.org/bioconda/get_mnv)
[![Version](https://img.shields.io/badge/version-1.1.4-%23149389?style=flat-square)](https://github.com/PathoGenOmics-Lab/get_MNV/releases)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.13907423-%23ff0077?style=flat-square)](https://doi.org/10.5281/zenodo.13907423)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square)](https://github.com/PathoGenOmics-Lab)

**Detección de variantes multinucleotídicas (Multi-Nucleotide Variant) - anotación a nivel de codón a partir de VCF o TSV de iVar.**
**Rust puro · sin dependencias de C · multiplataforma (macOS, Linux, Windows)**

[Inicio rápido](#inicio-rápido) · [GUI](#gui-de-escritorio) · [Características](#características) · [Documentación](docs/) · [Citación](#citación)

[English](README.md) · **Español**

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__y Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>

---

## ¿Qué es get_MNV?

get_MNV detecta los casos en que dos o más SNV caen en el mismo codón y deben interpretarse de forma conjunta. Estos cambios combinados pueden producir un efecto sobre el aminoácido distinto del de cada SNV por separado.

La herramienta toma como entrada:

- Llamadas de variantes: VCF o `variants.tsv` de iVar
- Secuencia de referencia: FASTA
- Anotación de genes: GFF/GFF3/GTF o un archivo TSV simple
- Lecturas alineadas opcionales: BAM, que se usan para contar el soporte de eventos SNP, MNV e indel

Escribe las variantes anotadas en formato TSV, VCF o ambos.

<p align="center">
  <img src="images/get_mnv_aa.png" alt="reclasificación de aminoácidos por MNV" width="650" />
</p>

**Características principales:**

- Agrupa los SNV por codón y notifica llamadas SNP, MNV o SNP/MNV
- Recalcula los cambios de aminoácido a partir del haplotipo completo del codón
- Descompone los alelos `REF/ALT` de VCF/iVar en componentes de evento de tipo SNV, MNV, inserción, deleción,
  delins e indel complejo
- Lee llamadas de variantes en VCF y en TSV de iVar, incluida la notación de indels
  `+SEQ` y `-SEQ` de iVar
- Cuando se proporciona un BAM, usa sus lecturas para contar el soporte SNP/MNV, el soporte exacto de eventos
  indel y el sesgo de hebra
- Admite 9 tablas de código genético del NCBI
- Incluye una GUI de escritorio para analizar mediante arrastrar y soltar

## Instalación

### GUI de escritorio

Descarga la última release para tu plataforma:

| Plataforma | Descarga |
|---|---|
| 🍎 macOS (Apple Silicon) | [**get_MNV_1.1.4_aarch64.dmg**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🍎 macOS (Intel) | [**get_MNV_1.1.4_x64.dmg**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🐧 Linux | [**Página de releases**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🪟 Windows | [**Página de releases**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |

> [!NOTE]
> **Usuarios de macOS**: La aplicación no está firmada con un certificado de Apple Developer. La primera vez que la abras, haz clic derecho sobre ella → **Abrir** → haz clic en **Abrir** en el cuadro de diálogo. Consulta el [soporte de Apple](https://support.apple.com/en-us/HT202491) para más detalles.

Todas las releases están disponibles en la [página de releases](https://github.com/PathoGenOmics-Lab/get_MNV/releases).

### Línea de comandos

```bash
conda install -c bioconda get_mnv
```

o descarga un binario precompilado:

```bash
wget https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest/download/get_mnv
chmod +x get_mnv
./get_mnv --help
```

o compila desde el código fuente:

```bash
git clone https://github.com/PathoGenOmics-Lab/get_MNV.git
cd get_MNV
cargo install --path .
```

## Inicio rápido

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
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3
```

Usa `--tsv` para el archivo `variants.tsv` que genera `ivar variants`.

### Con soporte de lecturas de BAM

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3
```

### Salida TSV y VCF

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both \
  --summary-json run.summary.json \
  --run-manifest run.manifest.json
```

Ejecuta `get_mnv --help` para ver la lista completa de opciones.

## Argumentos comunes

| Argumento | Qué hace |
|---|---|
| `--vcf <FILE>` | Archivo de variantes de entrada en formato `.vcf` plano o `.vcf.gz` comprimido con BGZF. La entrada BCF debe convertirse antes a VCF. |
| `--tsv <FILE>` | Archivo `variants.tsv` de iVar de entrada. |
| `--bam <FILE>` | BAM ordenado e indexado, opcional, para el soporte de lecturas. |
| `--fasta <FILE>` | FASTA de referencia. Los nombres de los contigs deben coincidir con el archivo de variantes. |
| `--gff <FILE>` | Anotación de genes en formato GFF/GFF3/GTF. |
| `--genes <FILE>` | Anotación de genes en TSV simple. Úsalo en lugar de `--gff`. |
| `--gff-features <LIST>` | Tipos de feature a analizar, por ejemplo `CDS` o `gene,pseudogene`. |
| `--quality <N>` | Calidad Phred de base mínima para el soporte de lecturas del BAM. Por defecto: `20`. |
| `--min-mapq <N>` | Calidad de mapeo mínima de las lecturas al usar BAM. Por defecto: `0`. |
| `--snp <N>` | Número mínimo de lecturas que soportan el SNP. Por defecto: `0`. |
| `--min-snp-frequency <F>` | Frecuencia mínima de SNP derivada del BAM, de `0` a `1`. Por defecto: `0`. |
| `--min-snp-strand <N>` | Número mínimo de lecturas que soportan el SNP exigido en cada hebra. Por defecto: `0`. |
| `--mnv <N>` | Número mínimo de lecturas que soportan el MNV. Por defecto: `0`. |
| `--min-mnv-frequency <F>` | Frecuencia mínima del haplotipo MNV derivada del BAM, de `0` a `1`. Por defecto: `0`. |
| `--min-mnv-strand <N>` | Número mínimo de lecturas que soportan el MNV exigido en cada hebra. Por defecto: `0`. |
| `--both` | Escribe las salidas TSV y VCF a la vez. |
| `--summary-json <FILE>` | Escribe un resumen de la ejecución legible por máquina. |
| `--run-manifest <FILE>` | Escribe el comando, la versión, las entradas, las salidas y las sumas de comprobación. |

Los filtros de frecuencia usan el soporte de lecturas recalculado a partir de `--bam`, no el valor
`OFREQ` original de la entrada VCF/iVar. Usa valores como `0.05` para el 5% o `0.20`
para el 20%. Cuando se solicita salida VCF, los registros de baja frecuencia se omiten por
defecto o se marcan con `FILTER=LowFrequency` si `--emit-filtered` está activado.
Los filtros de frecuencia de SNP y de MNV son independientes: `--min-snp-frequency` se aplica a
las observaciones individuales de SNP, mientras que `--min-mnv-frequency` se aplica al haplotipo
MNV en fase. En las llamadas mixtas `SNP/MNV`, un haplotipo MNV sólido no se elimina
solo porque las observaciones individuales de SNP queden por debajo del umbral de SNP.
Los filtros de conteo de lecturas y de soporte por hebra son independientes del mismo modo:
`--snp` y `--min-snp-strand` se aplican a las observaciones de SNP, mientras que `--mnv` y
`--min-mnv-strand` se aplican al haplotipo MNV.

## Salidas

Por defecto, get_MNV escribe:

```text
<input_name>.MNV.tsv
```

Con `--convert` o `--both`, también escribe:

```text
<input_name>.MNV.vcf
```

Los campos de salida más importantes son:

| Columna | Significado |
|---|---|
| `Chromosome` | Nombre del contig |
| `Gene` | Nombre del gen o de la feature |
| `Positions` | Una posición para los SNP, varias posiciones para los MNV |
| `Base Changes` | Bases alternativas |
| `AA Changes` | Cambio de aminoácido tras combinar los SNV del codón |
| `Variant Type` | `SNP`, `MNV`, `SNP/MNV` o `INDEL` |
| `Change Type` | Sinónimo, no sinónimo, ganancia/pérdida de codón de stop, desconocido, etc. |

Cuando se proporciona un BAM, hay columnas adicionales que indican la profundidad de lecturas, el soporte SNP, el soporte MNV, la frecuencia y los conteos por hebra.

## Características

| Característica | Descripción |
|---|---|
| 🧬 Detección de MNV | Agrupa los SNV del mismo codón y los reclasifica como MNV |
| 🔬 Cambios de AA precisos | Calcula los cambios de aminoácido a partir del haplotipo completo del codón |
| 📊 Soporte de lecturas | Conteos de lecturas SNP/MNV basados en BAM, con métricas específicas por hebra |
| 🔍 Sesgo de hebra | Valores p de la prueba exacta de Fisher para el soporte de sesgo de hebra de SNP y MNV (`SBP`/`MSBP` en el campo INFO del VCF) |
| 📁 Múltiples salidas | TSV, VCF (plano/BGZF+Tabix), BCF, resumen JSON, manifiesto de ejecución |
| ⚡ Paralelo | Procesamiento de contigs multihilo con Rayon |
| 🧪 Códigos genéticos | 9 tablas de traducción del NCBI (1, 2, 3, 4, 5, 6, 11, 12, 25) |
| 🧩 Entrada flexible | Llamadas de variantes en VCF o TSV de iVar; anotaciones GFF3/GTF o TSV; VCF multicontig y multimuestra |
| ✅ Validación | Modo de simulación (dry-run), métricas estrictas, sumas de comprobación de las entradas, JSON de errores |
| 🖥️ GUI de escritorio | Aplicación nativa de Tauri con arrastrar y soltar, visor de pistas genómicas y modo oscuro |

## GUI de escritorio

La aplicación de escritorio ofrece el mismo flujo de trabajo de análisis en una interfaz visual:

- Suelta los archivos de variantes VCF o TSV de iVar
- Suelta los archivos FASTA, GFF/GTF/GFF3 y el BAM opcional
- Elige los parámetros comunes en el formulario
- Ejecuta una muestra o varias muestras emparejadas
- Inspecciona, filtra y exporta los resultados

```bash
bash scripts/dev.sh   # desarrollo
bash scripts/build_gui_bundle.sh  # paquete de producción .app + .dmg
```

## Ejemplo de salida

```
Chromosome  Gene      Positions       Base Changes  AA Changes  Variant Type  Change Type
MTB_anc     Rv0095c   104838          T             Asp126Glu   SNP           Non-synonymous
MTB_anc     Rv0095c   104941,104942   T,G           Gly92Gln    SNP/MNV       Non-synonymous
MTB_anc     esxL      1341102,1341103 T,C           Arg33Ser    MNV           Non-synonymous
```

**Tipos de variante:**
- **SNP**: cambio de un solo nucleótido, un SNV por codón
- **MNV**: varios SNV se representan como un único haplotipo de codón combinado
- **SNP/MNV**: fila a nivel de codón con el contexto del SNV individual y el del haplotipo MNV combinado a la vez; con BAM, las columnas de soporte distinguen la evidencia
- **INDEL**: inserción, deleción, delins o alelo complejo; se notifica con los componentes del evento, el soporte exacto del BAM cuando está disponible y el efecto codificante cuando solapa con una feature CDS/gen anotada

En [`example/`](example/README.es.md) tienes un conjunto de datos de *M. tuberculosis* listo para ejecutar (referencia, genes, VCF y un BAM
de demostración diminuto para el visor de lecturas).

## Documentación

| Documento | Descripción |
|---|---|
| [Uso](docs/usage.es.md) | Referencia completa de la CLI y ejemplos |
| [Formatos de entrada](docs/input-formats.es.md) | Especificaciones de VCF, FASTA, GFF, TSV, BAM |
| [Formatos de salida](docs/output-formats.es.md) | Detalles de la salida TSV, VCF, BCF, JSON |
| [Semántica de indels y MNV](docs/indel-mnv-semantics.es.md) | Cómo se representan los indels, los MNV, los límites y los haplotipos complejos |
| [Resolución de problemas](docs/troubleshooting.es.md) | Errores comunes y soluciones |
| [Benchmarking](docs/benchmarking.es.md) | Pruebas de rendimiento |
| [Changelog](CHANGELOG.md) | Historial de versiones |

## Para desarrolladores

La CLI y la biblioteca principales están en `src/`. La aplicación de escritorio usa Tauri en
`src-tauri/` y React/TypeScript en `frontend/`.

Comandos útiles:

```bash
cargo test --workspace
npm run build --prefix frontend
bash scripts/build_get_mnv.sh
bash scripts/build_gui_bundle.sh
```

### Pruebas de escenarios de extremo a extremo

`tests/scenarios/` contiene un arnés en Python que construye entradas sintéticas de FASTA,
GFF, VCF (o TSV de iVar) y BAM a partir de escenarios declarativos, ejecuta
el binario `get_mnv` compilado y comprueba cada salida TSV frente a las
filas esperadas. La suite cubre actualmente 30 escenarios, entre ellos
la agrupación SNP/MNV a nivel de codón, la propagación de frameshift, la emisión de
haplotipos complex_indel, la anotación de CDS en hebra negativa y multiexón,
la división de multialélicos y la entrada de TSV de iVar.

```bash
cargo build                                # produce target/debug/get_mnv
python3 tests/scenarios/run.py             # ejecuta los 30 escenarios
python3 tests/scenarios/run.py 22 27       # ejecuta un subconjunto por prefijo de nombre
```

Requiere `samtools` en el `PATH` (o `SAMTOOLS=/path/to/samtools`). Consulta
[tests/scenarios/README.md](tests/scenarios/README.es.md) para ver la lista completa
de casos validados, la disposición del mini-genoma y cómo añadir nuevos escenarios.

## Limitaciones

- Está diseñado para eventos pequeños de tipo SNV/MNV/indel frente a una secuencia de referencia
- Con `--gff-features CDS`, los registros GFF/GTF que aportan `transcript_id` o `Parent` se reconstruyen como modelos de CDS empalmados, lo que permite anotar codones en uniones de exones y el contexto de frameshift de los indels a nivel de transcrito.
- Las variantes heterocigotas eucariotas sin fase siguen requiriendo cuidado: get_MNV reanota los alelos del llamador de variantes, pero no reestima la ploidía, las verosimilitudes de genotipo ni la fase de largo alcance.
- Los registros VCF multialélicos requieren `--split-multiallelic` o una división previa (`bcftools norm -m -`)
- Los nombres de los contigs de las variantes deben coincidir exactamente con el FASTA y el GFF
- **Múltiples transcritos por gen**: al usar `--gff-features CDS` con un archivo GFF que contiene varios transcritos para el mismo gen, cada transcrito se anota de forma independiente, lo que produce una línea de salida por transcrito y por variante. Si quieres una sola línea por variante, filtra tu GFF para conservar solo el transcrito canónico antes de ejecutar get_MNV (por ejemplo, con [AGAT](https://github.com/NBISweden/AGAT) `agat_sp_keep_longest_isoform.pl` o una herramienta similar)

## Citación

Si usas get_MNV en tu investigación, cita por favor:

> Ruiz-Rodriguez P, Coscolla M. **get_MNV: Multi-Nucleotide Variant detection tool.** Zenodo. doi: [10.5281/zenodo.13907423](https://doi.org/10.5281/zenodo.13907423)

```bibtex
@software{ruiz-rodriguez_get_mnv_2026,
  title     = {get\_MNV: Multi-Nucleotide Variant detection tool},
  author    = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
  year      = {2026},
  doi       = {10.5281/zenodo.13907423},
  url       = {https://github.com/PathoGenOmics-Lab/get_MNV},
  version   = {1.1.4},
  license   = {AGPL-3.0}
}
```

## Licencia

[GNU Affero General Public License v3.0](LICENSE)

## Curiosidad

Haz clic para ver el logo imprimible en 3D:

<p align="center">
  <a href="https://www.printables.com/model/1030383-get_mnv-logo" target="_blank">
    <img src="https://media.printables.com/media/prints/1030383/images/7820375_62fee28c-1ef3-446a-9187-3a74e3912c09_7526c3fd-6f35-4ec1-ab2c-a8022ac592e9/thumbs/inside/1920x1440/jpg/img_3773.webp" height="350" alt="logo 3D de get_MNV">
  </a>
</p>

---

<h2 id="contributors" align="center">

✨ Colaboradores
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START -->
<div align="center">
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Design">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td>
  </tr>
</table>

Este proyecto sigue la especificación [all-contributors](https://github.com/all-contributors/all-contributors) ([clave de emojis](https://allcontributors.org/docs/en/emoji-key)).
</div>
<!-- ALL-CONTRIBUTORS-LIST:END -->
