# Conjunto de datos de ejemplo

[English](README.md) · **Español**

Un pequeño conjunto de datos de *M. tuberculosis* para probar get_MNV de principio
a fin, incluido el visor de lecturas de la GUI.

| Archivo | Qué es |
|---|---|
| `MTB_ancestor.fas` (`.fai`) | Referencia, un único contig `MTB_anc` (4,411,532 bp) |
| `anot_genes.txt` | Tabla de genes: `name`, `start`, `end`, `strand` (separados por tabuladores) |
| `G35894.var.snp.vcf` | Llamadas VCF estilo VarScan para la muestra G35894 |
| `G35894.var.snp.MNV.tsv` | Salida de get_MNV precalculada para ese VCF (sin BAM) |
| `G35894.demo.bam` (`.bai`) | Alineamiento sintético minúsculo para el visor de lecturas (ver abajo) |
| `make_demo_bam.py` | Regenera `G35894.demo.bam` a partir de la referencia |

## Ejecuta la CLI

```bash
get_mnv \
  --vcf example/G35894.var.snp.vcf \
  --fasta example/MTB_ancestor.fas \
  --genes example/anot_genes.txt
```

Esto reproduce `G35894.var.snp.MNV.tsv`: 941 variantes (797 SNP, 10 SNP/MNV,
134 intergénicas) en 635 genes anotados.

## Visualiza las lecturas

`G35894.demo.bam` es un alineamiento **sintético y minúsculo** (24 lecturas, ~200 bp)
sobre un único MNV real de esta muestra, para que el soporte de lecturas y el pileup
de la GUI funcionen de fábrica:

> gen `Rv2036` (hebra +), posiciones genómicas `2282376` y `2282377`, codón `GTT → GCC`
> (Val93Ala). Todas las lecturas portan ambas bases ALT, repartidas a partes iguales entre las hebras.

Añade `--bam` para contar el soporte de lecturas a partir de él:

```bash
get_mnv \
  --vcf example/G35894.var.snp.vcf \
  --fasta example/MTB_ancestor.fas \
  --genes example/anot_genes.txt \
  --bam example/G35894.demo.bam --mnv 1
```

La fila `Rv2036` `2282376, 2282377` muestra entonces `MNV Reads=24` (12 forward /
12 reverse), `MNV Frequencies=1.0000`.

En la GUI de escritorio, carga el VCF, el FASTA, la tabla de genes y `G35894.demo.bam`,
ejecuta y luego abre esa fila `SNP/MNV` para ver las 24 lecturas en el pileup.

> El BAM de demostración solo cubre ese único locus, así que las lecturas se muestran
> únicamente ahí. Para un análisis real, proporciona tu propio BAM ordenado por coordenada e indexado.

## Regenerar el BAM de demostración

```bash
cd example && python3 make_demo_bam.py   # needs samtools on PATH (or $SAMTOOLS)
```
