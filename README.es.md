<p align="center">
  <img src="docs/assets/logo.svg#gh-light-mode-only" alt="logo de get_MNV" width="640">
  <img src="docs/assets/logo-dark.svg#gh-dark-mode-only" alt="logo de get_MNV" width="640">
</p>

<div align="center">

<a href="LICENSE"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/license-AGPL%20v3-%23af64d1?style=flat-square&amp;labelColor=21262d"><img alt="License: AGPL v3" src="https://img.shields.io/badge/license-AGPL%20v3-%23af64d1?style=flat-square"></picture></a>
<a href="https://anaconda.org/bioconda/get_mnv"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/conda/dn/bioconda/get_mnv.svg?style=flat-square&amp;label=bioconda&amp;labelColor=21262d"><img alt="Bioconda" src="https://img.shields.io/conda/dn/bioconda/get_mnv.svg?style=flat-square&amp;label=bioconda"></picture></a>
<a href="https://github.com/PathoGenOmics-Lab/get_MNV/releases"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/version-1.1.5-%23149389?style=flat-square&amp;labelColor=21262d"><img alt="Version" src="https://img.shields.io/badge/version-1.1.5-%23149389?style=flat-square"></picture></a>
<a href="https://doi.org/10.5281/zenodo.13907422"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.13907422-%23ff0077?style=flat-square&amp;labelColor=21262d"><img alt="DOI" src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.13907422-%23ff0077?style=flat-square"></picture></a>
<a href="https://pathogenomics-lab.github.io/get_MNV/es/"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/docs-en%20l%C3%ADnea-%230a7ea4?style=flat-square&amp;labelColor=21262d"><img alt="Documentación" src="https://img.shields.io/badge/docs-en%20l%C3%ADnea-%230a7ea4?style=flat-square"></picture></a>
<a href="https://github.com/PathoGenOmics-Lab/get_MNV/discussions/categories/q-a"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/pregunta-una%20duda-%23f5a623?style=flat-square&amp;labelColor=21262d"><img alt="Pregunta una duda" src="https://img.shields.io/badge/pregunta-una%20duda-%23f5a623?style=flat-square"></picture></a>
<a href="https://github.com/PathoGenOmics-Lab"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square&amp;labelColor=21262d"><img alt="PGO" src="https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square"></picture></a>

**Detección de variantes multinucleotídicas (Multi-Nucleotide Variant) - anotación a nivel de codón a partir de VCF o TSV de iVar.**
**Rust puro · sin dependencias de C · multiplataforma (macOS, Linux, Windows)**

[Documentación](https://pathogenomics-lab.github.io/get_MNV/es/) · [Instalación](#instalación) · [Inicio rápido](#inicio-rápido) · [Citación](#citación) · [Pregunta una duda](https://github.com/PathoGenOmics-Lab/get_MNV/discussions/categories/q-a)

[English](README.md) · **Español**

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__y Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>

---

## ¿Qué es get_MNV?

get_MNV detecta los casos en los que dos o más SNV caen en el mismo codón y
deben interpretarse de forma conjunta. Estos cambios combinados pueden producir
un efecto sobre el aminoácido distinto al de los SNV por separado.

Toma llamadas de variantes (VCF o `variants.tsv` de iVar), un FASTA de
referencia y una anotación de genes (GFF/GFF3/GTF o un TSV sencillo),
opcionalmente las lecturas alineadas, y escribe las variantes anotadas en TSV,
VCF o ambos, además de un informe HTML interactivo y autocontenido.

<p align="center">
  <img src="docs/assets/get_mnv_aa.png#gh-light-mode-only" alt="reclasificación de aminoácidos por MNV" width="650">
  <img src="docs/assets/get_mnv_aa-dark.png#gh-dark-mode-only" alt="reclasificación de aminoácidos por MNV" width="650">
</p>

- Agrupa los SNV por codón y reporta llamadas SNP, MNV o SNP/MNV
- Recalcula los cambios de aminoácido desde el haplotipo completo del codón
- Descompone los alelos `REF/ALT` en componentes SNV, MNV, inserción, deleción,
  delins e indel complejo
- Cuenta el soporte de SNP, MNV e indel exacto desde un BAM, con sesgo de hebra
- Admite 9 tablas de código genético del NCBI
- Incluye una GUI de escritorio con arrastrar y soltar y un visor de pistas genómicas

## Instalación

### GUI de escritorio

Descarga la última versión para tu plataforma:

| Plataforma | Descarga |
|---|---|
| 🍎 macOS (Apple Silicon) | [**get_MNV_1.1.5_aarch64.dmg**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🍎 macOS (Intel) | [**get_MNV_1.1.5_x64.dmg**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🐧 Linux | [**Página de releases**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🪟 Windows | [**Página de releases**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |

> [!NOTE]
> **Usuarios de macOS**: la aplicación no está firmada con un certificado de desarrollador de Apple. La primera vez, haz clic derecho sobre la app → **Abrir** → pulsa **Abrir** en el diálogo. Consulta el [soporte de Apple](https://support.apple.com/es-es/HT202491) para más detalles.

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

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both \
  --report run.html
```

`--bam` es opcional y es lo que convierte "estos dos cambios están en un codón"
en "estos dos cambios están en una misma molécula". Usa `--tsv` en vez de
`--vcf` para un `variants.tsv` de iVar, y `--genes` en vez de `--gff` para una
anotación de cuatro columnas. `get_mnv --help` lista todas las opciones, y la
[referencia de CLI](https://pathogenomics-lab.github.io/get_MNV/es/cli-reference/)
explica qué cambia cada una en la respuesta.

La salida tiene esta pinta:

```text
Chromosome  Gene      Positions       Base Changes  AA Changes  Variant Type  Change Type
MTB_anc     Rv0095c   104838          T             Asp126Glu   SNP           Non-synonymous
MTB_anc     Rv0095c   104941,104942   T,G           Gly92Gln    SNP/MNV       Non-synonymous
MTB_anc     esxL      1341102,1341103 T,C           Arg33Ser    SNP/MNV       Non-synonymous
```

En [`example/`](example/README.md) hay un conjunto de datos de *M. tuberculosis*
listo para ejecutar (referencia, genes, VCF y un BAM de demostración pequeño
para el visor de lecturas). El
[tutorial](https://pathogenomics-lab.github.io/get_MNV/es/getting-started/) lo
recorre entero.

## Documentación

**<https://pathogenomics-lab.github.io/get_MNV/es/>** es el manual, en español y
en inglés: cada opción con su valor por defecto, qué significa cada columna de
salida, y qué hace y qué no hace la herramienta. Las mismas páginas están en
[`docs/`](docs/) en este repositorio.

| Empieza por aquí | |
|---|---|
| [Tutorial de línea de comandos](https://pathogenomics-lab.github.io/get_MNV/es/getting-started/) | Una primera ejecución completa con los datos incluidos, explicando la salida línea a línea |
| [Recetas habituales](https://pathogenomics-lab.github.io/get_MNV/es/usage/) | Comandos listos para los trabajos de siempre |
| [Tutorial de la GUI](https://pathogenomics-lab.github.io/get_MNV/es/gui-tutorial/) | La misma ejecución en la aplicación, pantalla por pantalla |

| Referencia | |
|---|---|
| [Referencia de CLI](https://pathogenomics-lab.github.io/get_MNV/es/cli-reference/) | Todas las opciones, con su valor por defecto y qué cambian |
| [Formatos de entrada](https://pathogenomics-lab.github.io/get_MNV/es/input-formats/) | Cómo tienen que ser el VCF, el FASTA, la anotación y el BAM |
| [Formatos de salida](https://pathogenomics-lab.github.io/get_MNV/es/output-formats/) | Cada columna del TSV, clave INFO del VCF y campo JSON |
| [Informe de ejemplo](https://pathogenomics-lab.github.io/get_MNV/assets/example-report.html) | Un informe HTML real: ábrelo y trastea |

| Cómo funciona | |
|---|---|
| [Alcance y compatibilidad](https://pathogenomics-lab.github.io/get_MNV/es/indel-mnv-semantics/) | De qué se ocupa get_MNV, qué deja a tu llamador de variantes y dónde están sus límites |
| [Indels y haplotipos locales](https://pathogenomics-lab.github.io/get_MNV/es/indel-haplotypes/) | Cómo se lee un indel de los alineamientos y qué cuenta cada número |
| [Ligamiento](https://pathogenomics-lab.github.io/get_MNV/es/linkage/) | Distinguir un haplotipo real de dos variantes que solo comparten codón |
| [Resolución de problemas](https://pathogenomics-lab.github.io/get_MNV/es/troubleshooting/) | Los errores que detienen una ejecución, y qué te está diciendo cada aviso |

El historial de versiones está en [CHANGELOG.md](CHANGELOG.md).

Si el manual no lo responde, pregunta en **[Q&A](https://github.com/PathoGenOmics-Lab/get_MNV/discussions/categories/q-a)**. Las preguntas sobre por
qué una variante ha salido de una forma concreta van ahí, no al gestor de issues,
y una respuesta marcada como tal es la que encuentra el siguiente. Las issues son
para algo que se rompió, se negó a ejecutar o está claramente mal.

## Para desarrolladores

La CLI y la biblioteca principales están en `src/`. La aplicación de escritorio
usa Tauri en `src-tauri/` y React/TypeScript en `frontend/`.

```bash
cargo test --workspace
npm run build --prefix frontend
bash scripts/build_get_mnv.sh
bash scripts/build_gui_bundle.sh
```

`tests/scenarios/` es un arnés en Python que construye entradas sintéticas de
FASTA, GFF, VCF y BAM a partir de escenarios declarativos, ejecuta el binario
compilado y comprueba cada fila de la salida; necesita `samtools` en el `PATH`.
Consulta [tests/scenarios/README.es.md](tests/scenarios/README.es.md) para ver
los casos que cubre y cómo añadir uno.

```bash
cargo build                                # produce target/debug/get_mnv
python3 tests/scenarios/run.py             # ejecuta todos los escenarios
```

## Citación

Si usas get_MNV en tu investigación, cita por favor:

> Ruiz-Rodriguez P, Coscolla M. **get_MNV: Multi-Nucleotide Variant detection tool.** Zenodo. doi: [10.5281/zenodo.13907422](https://doi.org/10.5281/zenodo.13907422)

```bibtex
@software{ruiz-rodriguez_get_mnv_2026,
  title     = {get\_MNV: Multi-Nucleotide Variant detection tool},
  author    = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
  year      = {2026},
  doi       = {10.5281/zenodo.13907422},
  url       = {https://github.com/PathoGenOmics-Lab/get_MNV},
  version   = {1.1.5},
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
