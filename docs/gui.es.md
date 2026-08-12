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
   frecuencia y de hebra, y los ajustes finos de indels.
3. **Ejecuta.** Analiza una sola muestra, o varias muestras emparejadas en un mismo lote.
4. **Inspecciona y exporta.** Explora, ordena y filtra la tabla de resultados, y luego
   exporta a TSV o VCF.

Para un recorrido con capturas, mira el
[tutorial de la GUI de escritorio](gui-tutorial.es.md).

### En qué se aparta el formulario de la CLI

La app ejecuta el mismo motor, y todo ajuste que no muestra cae en el valor por
defecto de la CLI. Cinco de los que sí muestra son a propósito más conservadores,
así que los mismos archivos pueden dar aquí menos filas que `get_mnv` sin flags:

| Campo del formulario | App | CLI | Efecto |
|---|---:|---:|---|
| Min MAPQ | `20` | `0` | las lecturas multimapeadas no se cuentan |
| Min SNP reads | `2` | `0` | un SNP con una sola lectura se descarta |
| Min MNV reads | `2` | `0` | un haplotipo con una sola lectura se descarta |
| Normalize alleles | activado | desactivado | se recorta el contexto REF/ALT compartido |
| Split multiallelic | activado | desactivado | los registros multialélicos se dividen en vez de rechazarse |

Pon esos valores a los de la CLI en el formulario para reproducir exactamente una
ejecución de línea de comandos. La lista la impone un test que lee los valores por
defecto del frontend y los compara con el parser de la propia CLI, así que una
sexta divergencia rompe el build.

## Visor de pistas genómicas

<div style="text-align: center;" markdown>
![Visor de pistas genómicas de get_MNV: pistas de codones y pileup de lecturas de un MNV](assets/gui-06-reads.png){ width="840" }
</div>

*El visor de pistas para el MNV de `Rv2036` (`GTT → GCC`, Val93Ala) en el ejemplo incluido: pistas de codones y el pileup de lecturas, con las bases ALT resaltadas en las 24 lecturas que dan soporte.*

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
