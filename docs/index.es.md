# get_MNV

**Detección de variantes multinucleotídicas (MNV) a nivel de codón y anotación de indels desde VCF o iVar TSV.**
Rust puro · sin dependencias de C · multiplataforma (macOS, Linux, Windows).

get_MNV detecta los casos en los que dos o más SNV caen en el **mismo codón** y
deben interpretarse de forma conjunta. Estos cambios combinados pueden producir
un efecto sobre el aminoácido distinto al de los SNV por separado.

<div style="text-align: center;" markdown>
![Reclasificación de aminoácidos por MNV](assets/get_mnv_aa.png){ width="620" }
</div>

## Qué necesita

| Entrada | Formato |
|---|---|
| Llamadas de variantes | VCF o iVar `variants.tsv` |
| Secuencia de referencia | FASTA |
| Anotación de genes | GFF/GFF3/GTF o un TSV simple |
| Lecturas alineadas *(opcional)* | BAM, para contar el soporte de SNP/MNV/indel |

Escribe las variantes anotadas en TSV, VCF o ambos.

## Características principales

- Agrupa los SNV por codón y reporta llamadas **SNP**, **MNV** o **SNP/MNV**.
- Recalcula los cambios de aminoácido a partir del **haplotipo completo del codón**.
- Descompone los alelos `REF/ALT` en componentes SNV, MNV, inserción, deleción,
  delins e indels complejos.
- Lee VCF e iVar TSV, incluida la notación de indels `+SEQ` / `-SEQ` de iVar.
- Usa las lecturas del BAM (cuando se aportan) para el soporte de SNP/MNV, el
  soporte exacto de eventos indel y el sesgo de hebra.
- Admite 9 tablas de código genético del NCBI.
- Escribe un informe HTML interactivo y autocontenido:
  [ver un ejemplo](assets/example-report.html){ target=_blank }.
- Incluye una GUI de escritorio para análisis con arrastrar y soltar.

## Instalación

=== "GUI de escritorio"

    Descarga la última versión para tu plataforma desde la
    [página de Releases](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest).

    !!! note "Usuarios de macOS"
        La app no está firmada con un certificado de Apple Developer. En el
        primer arranque, haz clic derecho en la app → **Abrir** → pulsa
        **Abrir** en el diálogo.

=== "Línea de comandos"

    ```bash
    conda install -c bioconda get_mnv
    ```

    O compílala desde el código fuente:

    ```bash
    git clone https://github.com/PathoGenOmics-Lab/get_MNV.git
    cd get_MNV
    cargo install --path .
    ```

## Inicio rápido

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3
```

Hay un conjunto de datos de *M. tuberculosis* listo para usar (referencia, genes,
VCF y un BAM de demostración minúsculo para el visor de lecturas) en la carpeta
[`example/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/example)
del repositorio.

!!! tip "Siguientes pasos"
    - [Tutorial de línea de comandos](getting-started.es.md): una primera ejecución completa con los datos de ejemplo.
    - [Recetas habituales](usage.es.md): comandos listos para ejecutar de las tareas más frecuentes.
    - [Referencia de CLI](cli-reference.es.md): cada opción, con su valor por defecto.
    - [Formatos de salida](output-formats.es.md): qué significa cada columna y cada campo.
    - [Ligamiento](linkage.es.md): distinguir un haplotipo real de una coincidencia.
    - [Resolución de problemas](troubleshooting.es.md): soluciones a desajustes de entrada habituales.

## Cita

Si usas get_MNV, cítalo mediante su
[DOI de Zenodo](https://doi.org/10.5281/zenodo.13907423).

## Licencia

get_MNV se distribuye bajo la licencia **AGPL-3.0**.
