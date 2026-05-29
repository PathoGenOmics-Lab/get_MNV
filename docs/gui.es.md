# Interfaz gráfica de escritorio

get_MNV incluye una aplicación de escritorio nativa que ejecuta el mismo motor de
análisis que la línea de comandos, con entradas mediante arrastrar y soltar y un visor
interactivo de resultados.

## Instalación

Descarga la última versión para tu plataforma desde la
[página de Releases](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest).

!!! note "Usuarios de macOS"
    La app no está firmada con un certificado de Apple Developer. En el primer inicio,
    haz clic derecho en la app → **Abrir** → haz clic en **Abrir** en el cuadro de diálogo.

## Flujo de trabajo

1. **Añade las entradas.** Suelta tu archivo de variantes (VCF o TSV de iVar), el FASTA de
   referencia, una anotación génica (GFF/GFF3 o un TSV de genes) y, opcionalmente, un BAM
   ordenado por coordenada e indexado.
2. **Configura los parámetros.** El formulario expone las opciones habituales: código
   genético, umbrales de calidad y MAPQ, recuento de lecturas SNP/MNV, filtros de
   frecuencia y de hebra, y los ajustes finos de indels. Los valores por defecto coinciden
   con los de la CLI.
3. **Ejecuta.** Analiza una sola muestra, o varias muestras emparejadas en un mismo lote.
4. **Inspecciona y exporta.** Explora, ordena y filtra la tabla de resultados, y luego
   exporta a TSV o VCF.

## Visor de pistas genómicas

Al seleccionar una fila de variante se abre una vista de estilo IGV que alinea, columna a columna:

- una **regla** que marca las posiciones de la variante y la ventana mostrada;
- la secuencia de **referencia** y la **cobertura** por posición;
- **pistas de codones** que muestran el codón de referencia, los codones SNP individuales y el
  codón MNV combinado, con el cambio de aminoácido resultante;
- el **pileup de lecturas**, una fila por cada lectura que da soporte, con las bases ALT
  resaltadas y las lecturas coloreadas según el soporte que aportan (MNV / parcial /
  referencia).

Esto facilita confirmar visualmente que las lecturas portan el cambio de codón combinado,
y no solo los SNV individuales.

!!! tip "Pruébalo con el ejemplo"
    Sigue el tutorial de [Primeros pasos](getting-started.es.md) y carga el
    `example/G35894.demo.bam` incluido. Al abrir la fila de `Rv2036` se muestran las 24
    lecturas que portan ambas bases ALT del codón `GTT → GCC` (Val93Ala).

!!! note "Requisitos del visor de lecturas"
    El pileup necesita un BAM ordenado por coordenada e indexado (`.bai`) y un FASTA
    indexado (`.fai`); get_MNV crea el índice del FASTA automáticamente cuando es necesario.
    Las filas sin datos de lecturas siguen mostrando las pistas de codones y aminoácidos.

## Compilar desde el código fuente

La interfaz gráfica es una app de [Tauri](https://tauri.app/) (backend en Rust + frontend
web). Para ejecutarla desde una copia del repositorio:

```bash
bash scripts/dev.sh                # development
bash scripts/build_gui_bundle.sh   # production .app / .dmg bundle
```
