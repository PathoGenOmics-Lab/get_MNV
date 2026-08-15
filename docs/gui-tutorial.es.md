# Tutorial de la GUI de escritorio

Una ejecución completa de la aplicación de escritorio con los datos que vienen
con get_MNV, pantalla a pantalla: cargar los archivos, qué cambia cada
parámetro, cómo leer el resumen y cómo mirar las lecturas que hay detrás de una
llamada.

Si aún no has instalado la app, mira [GUI de escritorio](gui.es.md#instalacion).
El mismo recorrido en línea de comandos está en
[Tutorial de línea de comandos](getting-started.es.md).

!!! note "De dónde salen estas capturas"
    Están tomadas del demo de navegador que vive en `frontend/demo/`, que
    renderiza los componentes reales de la app contra datos de fixture en vez de
    llamar al backend de Rust. Los números, las lecturas y la tabla son la salida
    auténtica de `get_mnv` sobre los archivos de `example/`; lo único que se
    sustituye es la fontanería entre el formulario y el motor, para poder
    capturar las pantallas sin compilar un bundle de escritorio.

## Qué necesitas

Está todo en la carpeta [`example/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/example)
del repositorio:

| Archivo | Papel |
|---|---|
| `G35894.var.snp.vcf` | las llamadas de variantes |
| `MTB_ancestor.fas` | la referencia contra la que se llamaron |
| `anot_genes.txt` | una tabla de genes simple |
| `G35894.demo.bam` (+ `.bai`) | un alineamiento pequeño, para el soporte de lecturas |

## 1. Carga las entradas

![La pestaña Analysis con las cuatro entradas rellenas](assets/gui-02-inputs.png#only-light)
![La pestaña Analysis con las cuatro entradas rellenas](assets/gui-02-inputs-dark.png#only-dark)

Suelta cada archivo en su zona, o haz clic para buscarlo. Tres son obligatorios y
llevan un asterisco rojo; el contador marca **3/3 required + BAM** en cuanto
están puestos, y el botón de ejecutar se activa.

- **Variant calls** acepta un VCF plano o comprimido con BGZF, o un
  `variants.tsv` de iVar. La app detecta cuál es a partir del archivo.
- **FASTA reference** tiene que ser la referencia contra la que se llamaron las
  variantes. get_MNV escribe el índice `.fai` por su cuenta si falta.
- **Gene annotation** acepta un GFF/GFF3 o una tabla de genes simple.
  `anot_genes.txt` es una tabla de genes, así que no aparece el selector de
  features; carga un GFF y saldrá un panel **GFF features** donde eliges qué
  tipos leer (`gene,pseudogene` por defecto, `CDS` para transcritos con
  splicing).
- **BAM alignment** es opcional, y es lo que convierte una anotación en
  evidencia. Sin él, get_MNV reporta lo que dijo tu llamador. Con él, cuenta las
  lecturas por su cuenta y puede distinguir un haplotipo real a nivel de codón de
  dos sustituciones que nunca compartieron una molécula. Mira
  [Ligamiento](linkage.es.md).

## 2. Ajusta los parámetros

![La barra lateral de parámetros, con los mínimos de lecturas bajados a cero](assets/gui-03-parameters.png#only-light)
![La barra lateral de parámetros, con los mínimos de lecturas bajados a cero](assets/gui-03-parameters-dark.png#only-dark)

La barra lateral agrupa todos los ajustes que expone la app, y los cuatro chips
de preset de arriba los cambian en bloque. En cuanto tocas uno, el chip pasa a
**Custom**, que es lo que muestra la captura.

!!! warning "Por qué esta captura pone 0 y no el 2 de por defecto"
    El formulario viene con **Min SNP reads** y **Min MNV reads** en `2`, y esos
    umbrales solo actúan cuando hay un BAM cargado. `G35894.demo.bam` es un
    archivo de demostración diminuto que cubre un único locus, así que con el
    valor de por defecto todas las demás filas se quedan sin soporte y la salida
    baja de 941 filas a 1. Por eso aquí están los dos en `0`. Con un alineamiento
    real, déjalos como están. El [tutorial de línea de comandos](getting-started.es.md#5-filtra-por-ese-soporte)
    explica cómo se combinan los dos umbrales, que no es evidente.

Cinco de los valores por defecto del formulario son a propósito más estrictos que
los de la CLI; están listados en
[GUI de escritorio](gui.es.md#en-que-se-aparta-el-formulario-de-la-cli). Todo lo
que el formulario no muestra cae en el valor por defecto de la CLI.

## 3. Ejecuta

![La ejecución en marcha, mostrando la fase actual](assets/gui-04-running.png#only-light)
![La ejecución en marcha, mostrando la fase actual](assets/gui-04-running-dark.png#only-dark)

El botón va informando de la fase en la que está. Una ejecución sobre este
conjunto de datos tarda bastante menos de un segundo; la barra de progreso
importa en cohortes, donde puedes encolar varias muestras emparejadas y dejarlas
correr en un mismo lote.

## 4. Lee el resumen

![La pestaña Results: recuentos, desglose de variantes, tiempos y cifras por contig](assets/gui-05-results.png#only-light)
![La pestaña Results: recuentos, desglose de variantes, tiempos y cifras por contig](assets/gui-05-results-dark.png#only-dark)

La fila de arriba es la forma de la ejecución: **941 variantes** producidas a
partir de **950 registros del VCF** sobre **635 genes mapeados** en **1 contig**.
Registros y variantes no coinciden porque un codón que lleva dos sustituciones se
reporta una sola vez.

**Variant breakdown** es la parte que merece leerse despacio:

| Fila | Aquí | Qué cuenta |
|---|---:|---|
| SNP | 797 | una sustitución en un codón |
| MNV | 0 | un único registro del VCF que ya traía más de una base sustituida |
| SNP/MNV | 10 | registros separados del VCF que get_MNV encontró en el mismo codón |
| Indel | 0 | inserciones y deleciones |
| Intergenic | 134 | fuera de cualquier gen anotado |

`MNV` es `0` y `SNP/MNV` es `10` porque este llamador emitió un registro por
base. Esas diez filas son lo que un anotador por-SNV habría reportado como veinte
sustituciones independientes con los aminoácidos equivocados, y son la razón de
ser de la herramienta.

Debajo, **per-contig breakdown** repite las cifras por secuencia, y **output
files** muestra lo que se ha escrito, con un botón para abrirlo en tu gestor de
archivos.

## 5. Mira las lecturas

Busca `Rv2036` en la lista de loci y ábrelo. Es el codón que recorre
[Tutorial de línea de comandos](getting-started.es.md).

![El visor de tracks genómicos: regla, cobertura, referencia, tracks de codón y el pileup de lecturas](assets/gui-06-reads.png#only-light)
![El visor de tracks genómicos: regla, cobertura, referencia, tracks de codón y el pileup de lecturas](assets/gui-06-reads-dark.png#only-dark)

Los tracks se alinean columna a columna sobre las mismas coordenadas:

- la **regla**, con las dos posiciones de la variante marcadas en rojo
  (`2282376`, `2282377`);
- la **cobertura**, que aquí llega a 24x;
- la secuencia de **referencia**;
- los **tracks de codón**: referencia `GTT`, los codones SNP individuales `GCT` y
  `GTC`, y el **codón MNV combinado `GCC`**;
- el **pileup de lecturas**, una fila por lectura, cada una etiquetada con el
  soporte que aporta y la hebra de la que viene.

Las 24 lecturas están marcadas como `ALT`, 12 en cada hebra, y la fila reporta
`MNV Frequencies 1.0000`. Leídas por separado, las dos sustituciones dicen cosas
distintas: `GCT` sola es Val93Ala, `GTC` sola es Val93Val, un cambio silencioso.
Juntas dan `GCC`, que es Ala. Esa es una llamada que no se puede acertar base a
base.

!!! tip "Por qué los recuentos de SNP salen a cero"
    La fila muestra `SNP Reads 0, 0` al lado de `MNV Reads 24`. `SNP Reads`
    cuenta las lecturas que llevan una sustitución **sin** el haplotipo completo.
    Aquí todas llevan las dos, así que ninguna se cuenta como soporte suelto y
    las 24 se van a `MNV Reads`. Las dos columnas reparten la evidencia en vez de
    contarla dos veces. Mira
    [Formatos de salida](output-formats.es.md#salida-tsv-por-defecto).

## 6. Filtra y exporta

![La tabla de variantes con filtros por columna](assets/gui-07-table.png#only-light)
![La tabla de variantes con filtros por columna](assets/gui-07-table-dark.png#only-dark)

La tabla contiene todas las filas del TSV. Busca en todas las columnas a la vez,
o filtra de una en una: los desplegables toman un valor y las cajas de texto
filtran según escribes. Ordena haciendo clic en una cabecera, expande a pantalla
completa, y **Export** escribe la vista actual a TSV o VCF.

## La misma ejecución en línea de comandos

El estado del formulario en estas capturas son los valores por defecto de la app
con los dos mínimos de lecturas bajados, que en línea de comandos es:

```bash
get_mnv \
  --vcf example/G35894.var.snp.vcf \
  --fasta example/MTB_ancestor.fas \
  --genes example/anot_genes.txt \
  --bam example/G35894.demo.bam \
  --min-mapq 20 --normalize-alleles --split-multiallelic
```

`--snp` y `--mnv` ya valen `0` en línea de comandos, así que no hacen falta.
Añade `--snp 2 --mnv 2` para obtener lo que hace el formulario sin tocar nada, y
sobre este conjunto de datos la salida se queda en la única fila que tiene
soporte de lecturas.

## Por dónde seguir

- [Formatos de salida](output-formats.es.md) para saber qué significa cada
  columna.
- [Ligamiento](linkage.es.md) para ver cómo decide get_MNV que unas variantes
  viajan juntas de verdad.
- [Referencia de CLI](cli-reference.es.md) para las opciones que el formulario no
  expone.
