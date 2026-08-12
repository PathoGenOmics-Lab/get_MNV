# Primeros pasos

Este breve tutorial ejecuta get_MNV de principio a fin sobre el ejemplo incluido
de *M. tuberculosis* y muestra cómo inspeccionar un MNV real a nivel de codón,
incluido su soporte de lecturas.

## 1. Obtén los datos de ejemplo

La carpeta [`example/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/example)
del repositorio incluye todo lo que necesitas:

| Archivo | Función |
|---|---|
| `MTB_ancestor.fas` | Referencia (un único contig `MTB_anc`) |
| `anot_genes.txt` | Tabla de genes (`name`, `start`, `end`, `strand`) |
| `G35894.var.snp.vcf` | Llamadas de variantes al estilo VarScan |
| `G35894.demo.bam` | Pequeño subconjunto de lecturas alineadas para el visor de lecturas |

```bash
git clone https://github.com/PathoGenOmics-Lab/get_MNV.git
cd get_MNV/example
```

## 2. Anota las variantes

```bash
get_mnv \
  --vcf G35894.var.snp.vcf \
  --fasta MTB_ancestor.fas \
  --genes anot_genes.txt
```

Esto escribe `G35894.var.snp.MNV.tsv` junto a la entrada. Para esta muestra
produce 941 variantes: 797 SNP, 10 SNP/MNV y 134 intergénicas.

!!! note "Código genético"
    get_MNV usa por defecto la tabla de traducción 11 del NCBI (bacteriana), que es correcta
    para *M. tuberculosis*. Usa `--translation-table` para otros organismos.

## 3. Lee la salida

Abre el TSV. Una fila destaca: el gen `Rv2036` lleva dos SNV en el mismo
codón:

| Positions | Reference Bases | Base Changes | Reference Codon | MNV Codon | AA Changes | Variant Type |
|---|---|---|---|---|---|---|
| `2282376, 2282377` | `T, T` | `C, C` | `GTT` | `GCC` | `Val93Ala` | `SNP/MNV` |

Anotado un SNV a la vez, cada cambio se lee de forma distinta; combinados, el codón
`GTT → GCC` es una única sustitución Val→Ala. Esa reclasificación consciente del
codón es exactamente para lo que sirve get_MNV. Consulta [Formatos de salida](output-formats.es.md) para
ver todas las columnas.

## 4. Añade soporte de lecturas (BAM)

Proporciona un BAM alineado e indexado para contar cuántas lecturas llevan realmente el
cambio combinado:

```bash
get_mnv \
  --vcf G35894.var.snp.vcf \
  --fasta MTB_ancestor.fas \
  --genes anot_genes.txt \
  --bam G35894.demo.bam --mnv 1
```

La fila de `Rv2036` ahora informa de `MNV Reads = 24` (12 directas / 12 inversas) y
`MNV Frequencies = 1.0000`: todas las lecturas que abarcan el codón llevan ambas
bases ALT.

!!! tip "Requisitos del BAM"
    El BAM debe estar ordenado por coordenada e indexado (`.bai`), y alineado a la
    misma referencia. El `G35894.demo.bam` incluido cubre solo ese único locus.

## 5. Visualiza las lecturas en la GUI

La aplicación de escritorio dibuja un pileup al estilo IGV para cualquier fila. Carga el VCF, el FASTA, la tabla
de genes y `G35894.demo.bam`, ejecuta el análisis y luego abre la fila de `Rv2036` para ver
las 24 lecturas con las dos bases ALT resaltadas. Consulta la
guía de la [GUI de escritorio](gui.es.md).

## Próximos pasos

- [Recetas habituales](usage.es.md): comandos listos para ejecutar de las tareas más frecuentes.
- [Referencia de CLI](cli-reference.es.md): cada opción, con su valor por defecto.
- [Formatos de entrada](input-formats.es.md) y [Formatos de salida](output-formats.es.md).
- [Alcance y compatibilidad](indel-mnv-semantics.es.md): reglas de límites y ajustes.
- [Ligamiento](linkage.es.md): distinguir un haplotipo real de una coincidencia.
