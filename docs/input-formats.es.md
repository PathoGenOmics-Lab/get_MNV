# Formatos de entrada

get_MNV necesita tres archivos: llamadas de variantes, un FASTA de referencia y
anotación de genes. El archivo BAM es opcional.

## 1. Llamadas de variantes

Pasa las llamadas de variantes con:

```bash
--vcf <VCF_FILE>
# o
--tsv <IVAR_TSV_FILE>
```

Las entradas de llamadas de variantes admitidas son:

- VCF (`.vcf` o `.vcf.gz`) y BCF (`.bcf`)
- TSV de variantes de iVar (`.tsv`)

Usa `--vcf` para archivos `.vcf` simples, `.vcf.gz` comprimidos con BGZF o
`.bcf`, y `--tsv` para archivos `variants.tsv` de iVar. Un BCF se decodifica a
texto VCF a medida que se lee, así que se anota exactamente igual que el VCF que
codifica; se reconoce por su extensión `.bcf`.
Los comandos antiguos que pasan un TSV de iVar mediante `--vcf` se siguen
detectando de forma automática cuando la cabecera tiene las columnas estándar de
iVar.

### VCF

Usa un archivo VCF estándar que contenga llamadas SNV/MNV, indels o alelos
complejos.

Requisitos:

- Los nombres de contig del VCF deben coincidir con los del FASTA y los del
  GFF/GTF.
- Los alelos REF deben coincidir con la secuencia del FASTA.
- Los registros multialélicos se dividen en una fila por alelo ALT.

get_MNV puede leer las métricas originales de profundidad y frecuencia desde
campos INFO o FORMAT habituales, como `DP`, `AF`, `FREQ`, `AD`, `AO` y `RO`.

Estos valores de frecuencia de entrada se conservan en el informe como `OFREQ`.
Los filtros de frecuencia de la línea de comandos (`--min-snp-frequency`,
`--min-mnv-frequency`) usan en su lugar el soporte de lecturas derivado del BAM,
por lo que requieren `--bam`.

### TSV de iVar

Usa el TSV generado por `ivar variants`.

Columnas obligatorias:

| Columna | Significado |
|---|---|
| `REGION` | Nombre del contig |
| `POS` | Posición en base 1 |
| `REF` | Base de referencia |
| `ALT` | Base alternativa |

Columnas opcionales que se usan cuando están presentes:

| Columna | Se usa como |
|---|---|
| `TOTAL_DP` | Profundidad original (`ODP`) |
| `ALT_FREQ` | Frecuencia original (`OFREQ`) |
| `REF_DP`, `ALT_DP` | Sirven para inferir profundidad y frecuencia si hace falta |
| `PASS` | Sirve para conservar las filas que pasan el filtro |

Filtrado:

- Si existe `PASS`, get_MNV conserva los valores con sentido verdadero, como
  `TRUE`, `PASS`, `1` o `YES`.
- Las filas donde `REF == ALT` se descartan.
- La notación de indels de iVar como `+SEQ` o `-SEQ` se convierte en alelos
  anclados al estilo VCF usando la referencia FASTA y, a continuación, se analiza
  con el mismo modelo de eventos de alelo que la entrada VCF.
- `ALT_FREQ` se incluye en el informe como frecuencia original (`OFREQ`), y es
  independiente de los filtros de frecuencia derivados del BAM.

## 2. FASTA de referencia

Pasa la referencia con:

```bash
--fasta reference.fasta
```

Requisitos:

- Los IDs de registro del FASTA deben coincidir con los nombres de contig de las
  variantes.
- Las bases deben ser bases de ADN IUPAC válidas.
- No se admiten nombres de contig duplicados.

## 3. Anotación de genes

Proporciona `--gff` o `--genes`.

### GFF/GFF3/GTF

Recomendado para la mayoría de los conjuntos de datos:

```bash
--gff genes.gff3
```

Por defecto, get_MNV analiza las features `gene,pseudogene`. Para la anotación de
codones de regiones codificantes de proteínas, usa:

```bash
--gff-features CDS
```

Detalles importantes:

- Las coordenadas se leen de las columnas 4 y 5.
- La hebra se lee de la columna 7.
- En las features `CDS`, se usa la fase de la columna 8 cuando está presente.
- En las filas `CDS` con `transcript_id` o `Parent`, get_MNV construye la
  secuencia CDS empalmada de cada transcrito. La agrupación de codones, los
  efectos sobre el aminoácido de los MNV y el contexto de frameshift de los
  indels se evalúan entonces sobre el CDS completo del transcrito.
- Si un GFF/GTF contiene varios transcritos para el mismo gen, una sola variante
  puede producir una línea de salida por cada transcrito solapante. Cada uno se
  anota por su cuenta, que es la respuesta honesta cuando la anotación ofrece
  varios. Si quieres una línea por variante, deja un solo transcrito por gen en
  el GFF antes de ejecutar, por ejemplo con
  `agat_sp_keep_longest_isoform.pl` de [AGAT](https://github.com/NBISweden/AGAT).

Los nombres de genes se leen de atributos habituales como `gene_name`, `gene`,
`Name`, `locus_tag`, `gene_id` e `ID`.

### Anotación TSV simple

Usa `--genes` para un archivo de anotación pequeño y sencillo:

```bash
--genes genes.tsv
```

Formato de cuatro columnas:

```text
GeneName	GeneStart	GeneEnd	Strand
```

Ejemplo:

```text
Rv0007_Rv0007	9914	10828	+
Rv0008c_Rv0008c	11874	12311	-
```

Formato opcional de cinco columnas con fase:

```text
GeneName	GeneStart	GeneEnd	Strand	Phase
```

La fase puede ser `0`, `1`, `2` o `.`. Si se omite la columna de fase, toma el
valor `0` por defecto.

Formato opcional de seis columnas con biotipo:

```text
GeneName	GeneStart	GeneEnd	Strand	Phase	Biotype
```

```text
Rv0007_Rv0007	9914	10828	+	0	protein_coding
mcr11_RVnc0013	1413094	1413224	-	0	ncRNA
```

Los valores aceptados son `protein_coding`, `coding`, `CDS` y `mRNA` para
features que se traducen, y `ncRNA`, `rRNA`, `tRNA`, `tmRNA`, `miRNA`, `snRNA`,
`snoRNA`, `misc_RNA`, `antisense_RNA`, `SRP_RNA`, `RNase_P_RNA`, `non_coding` y
`pseudogene` para las que no. Cualquier otro valor se rechaza con un error en
vez de adivinarse, porque adivinar mal inventa una proteína o esconde una real.

Una feature no codificante se reporta contra su gen como
`non_coding_transcript_exon_variant` (`MODIFIER`), sin cambio de aminoácido.
**Si se omite la columna de biotipo se asume que toda feature es codificante**,
que es lo que han hecho siempre los ficheros de cuatro y cinco columnas: un gen
de RNA se traduce entonces como si fuera una proteína. Declara el biotipo cuando
tu anotación contenga features no codificantes.

La entrada GFF y GTF lo dice por sí sola. Decide el atributo `gene_biotype`,
`biotype` o `transcript_biotype` del registro, leído con el mismo vocabulario que
la columna de arriba; y si no lo trae, decide el tipo de feature de la columna 3,
así que las filas `rRNA`, `tRNA`, `ncRNA` y `pseudogene` se leen como no
codificantes sin necesidad de decírselo. Una feature que no dice nada de sí misma
se asume codificante, que es lo que es una fila `gene` bacteriana. Ten en cuenta
que la lista por defecto de `--gff-features` es `gene,pseudogene`, así que los
pseudogenes entran por defecto y no se traducen.

Séptima columna opcional con el contig:

```text
GeneName	GeneStart	GeneEnd	Strand	Phase	Biotype	Contig
```

```text
Rv0007_Rv0007	9914	10828	+	0	protein_coding	NC_000962.3
geneD	1	300	+	0	protein_coding	plasmid1
```

Una fila que nombra un contig solo se aplica a las variantes de ese contig. Una
fila sin él, como en los formatos más cortos, se aplica al contig que se esté
anotando, que es todo lo que necesita un archivo de un solo contig. Por eso un
VCF con más de un contig necesita que todas las filas nombren su contig; si no,
la ejecución se detiene en vez de poner los genes de un contig sobre las
variantes de otro. Para usarla hay que escribir también las columnas de fase y
biotipo: `.` vale como fase desconocida.

## 4. Lecturas BAM o CRAM (opcional)

Pasa las lecturas alineadas con:

```bash
--bam reads.bam
--bam reads.cram
```

Cuando se proporcionan lecturas, get_MNV calcula:

- El soporte de lecturas de los SNP
- El soporte de lecturas del haplotipo MNV
- La profundidad total y la frecuencia
- Los recuentos por hebra directa e inversa
- Estadísticas opcionales de sesgo de hebra

Requisitos:

- El archivo debe estar ordenado por coordenadas.
- Debe estar indexado: `.bai` (o `.csi`) para un BAM, `.crai` para un CRAM.
  `samtools index` escribe cualquiera de los dos.
- Los nombres de contig deben coincidir con los del archivo de variantes y los
  del FASTA.
- Un CRAM guarda las lecturas como diferencias con una referencia, así que se
  decodifica con el `--fasta` de la ejecución. Tiene que ser la referencia con
  la que se escribió el CRAM; una distinta falla con un error de suma de
  comprobación en vez de dar bases equivocadas. Un `.fai` junto al FASTA permite
  leer la secuencia bajo demanda; sin él se carga entera.

Un CRAM da los mismos recuentos que el BAM del que sale: la batería de pruebas
ejecuta cada escenario de las dos maneras y compara las salidas byte a byte.
- Las lecturas duplicadas, secundarias, suplementarias y marcadas como QC-fail se ignoran, igual que las descarta `samtools mpileup` por defecto.

## 5. Catálogo de variantes conocidas (opcional)

`--catalogue conocidas.tsv` marca las llamadas que coinciden con una tabla de
variantes de las que ya sabes algo: mutaciones de resistencia, marcadores de
linaje, cambios patogénicos conocidos. Va separada por tabuladores, con una
cabecera que nombra una columna `gene` y otra `variant` y, de forma opcional,
`label`; las demás columnas se ignoran y las líneas que empiezan por `#` son
comentarios.

```text
gene	variant	label
rpoB	p.Ser450Leu	RIF: Assoc w R
fabG1	c.-15C>T	INH: Assoc w R
katG	LoF	INH: Assoc w R
```

Una variante se puede escribir como:

| Forma | Ejemplo | Coincide con |
|---|---|---|
| Proteína | `p.Ser450Leu`, `p.S450L`, `p.Gln10Ter`, `p.Ala11fs` | las filas del gen con ese cambio proteico. Un corrimiento de marco coincide por su primer residuo y su posición, siga lo que siga a `fs`. |
| Sustitución codificante | `c.1349C>T`, `c.-15C>T`, `c.*3G>A` | la base sustituida, colocada en el genoma con el modelo del gen (`c.-N` antes de la primera base codificante, `c.*N` después de la última), en cualquier fila que la lleve, también las intergénicas |
| Sustitución no codificante | `n.1401A>G`, `n.-119C>T` | lo mismo, en un gen que no se traduce, como un ARNr, numerado desde su primera base |
| Sustitución genómica | `g.761155C>T` | la base sustituida en esa posición |
| Otro HGVS | `c.516_518dup`, `g.1673425del` | las filas del gen cuyo `HGVS c.` se escribe exactamente así, o las filas de cualquier punto del contig del gen, también las intergénicas, cuyo `HGVS g.` se escribe así |
| `LoF` | `LoF` | cualquier corrimiento de marco, parada prematura o pérdida del inicio en el gen |

El gen se compara por nombre, o por cualquier parte separada por `_` de un
nombre compuesto, así que `rpoB` coincide con un gen llamado `rpoB_Rv0667`. Una
sustitución cuya base de referencia no es la que tiene el FASTA en esa posición
se escribió contra otra referencia; no coincide con nada, y la ejecución dice
cuántas hubo.

Una entrada proteica se compara con las dos lecturas de una fila de codón: el
residuo que da el codón entero y el que daría cada sustitución sola. Cuando solo
coincide una de ellas, la fila lo dice: `codon reading only` para un cambio del
catálogo que ninguna base sola produce, y que la anotación base a base no ve, y
`per-SNV reading only` para uno que esa anotación daría y que el codón no
produce. Las coincidencias van a la columna `Catalogue` (consulta
[Formatos de salida](output-formats.es.md)).

### El catálogo de la OMS de mutaciones de resistencia de *M. tuberculosis*

El catálogo de mutaciones del complejo *Mycobacterium tuberculosis* de la OMS
(segunda edición, 2023) clasifica unos 48.000 pares de variante y fármaco en
cinco grupos, de `1) Assoc w R` a `5) Not assoc w R`. `example/who_catalogue_2023.tsv`
recoge los grupos 1 y 2, las variantes asociadas a resistencia, listas para
`--catalogue`; en [Recetas habituales](usage.es.md#marcar-mutaciones-de-resistencia-conocidas-en-m-tuberculosis)
está la ejecución. Cada entrada lleva como etiqueta la abreviatura del fármaco
que usa la OMS y el grupo, así que una coincidencia se lee
`RIF: Assoc w R (rpoB p.Ser450Leu)`, y una variante clasificada para varios
fármacos es una entrada por fármaco. Otras tres columnas guardan el nombre del
fármaco, el grupo y el nombre de la variante de la OMS (`rpoB_p.Ser450Leu`),
para filtrar la tabla o volver a cruzarla con el catálogo; get_MNV las ignora.

`tools/who_catalogue_to_tsv.py` construye la tabla a partir de los archivos de
la OMS, así que se puede rehacer con otros grupos (`--grades 1,2,5` añade las
variantes que se ha visto que no se asocian, `--grades all` conserva todos) o
colocar sobre otra referencia. Solo necesita la biblioteca estándar de Python.
La mayoría de las entradas se escriben como get_MNV escribe el mismo cambio y se
conservan tal cual; las demás se reescriben:

| La OMS escribe | La tabla guarda |
|---|---|
| `p.Ser450Leu`, `p.Trp122*`, `p.Asp379fs`, `p.Met1?`, un `del`, `ins` o `delins` en fase | lo mismo |
| `c.-15C>T`, `c.102G>A`, `n.1401A>G`, `LoF` | lo mismo |
| `p.Phe433dup`, una duplicación en fase | `p.Phe433_Met434insPhe`, la inserción como la nombra get_MNV, con los residuos sacados de las demás entradas del catálogo |
| `p.Ter187Argext*?`, una parada perdida | `p.Ter187Arg`, el residuo que la sustituye |
| `c.421_423delTTGinsCTC`, misma longitud a la entrada y a la salida | `c.[421T>C;423G>C]`, las sustituciones que produce |
| `c.-5delG` o `n.2269_2270insT`, un indel fuera de una secuencia codificante | el `HGVS g.` que get_MNV da al cambio, `g.2289247del`, a partir de las coordenadas genómicas de la OMS y del FASTA de la ejecución (`--coordinates`, `--fasta`) |
| `deletion`, el gen entero | nada: la OMS no le da ninguna variante genómica, y un VCF de cambios dentro de los genes no tiene fila para un gen que a la muestra le falta |

Una entrada reescrita conserva en su etiqueta el nombre de la OMS:
`PZA: Assoc w R, WHO c.-5delG (pncA g.2289247del)`. `--unconverted` escribe lo
que no se ha convertido, con el motivo. En los grupos 1 y 2 son las cuatro
deleciones de gen entero, de katG, ethA, fgd1 y pncA, que tienen cada una una
entrada `LoF` para el mismo fármaco; sin `--coordinates` se añaden los tres
indels fuera de secuencias codificantes (pncA `c.-5delG`, eis `c.-8delC`, rrl
`n.2269_2270insT`).

Comparar cambios proteicos, y no registros genómicos, es lo que permite
encontrar un codón con dos sustituciones. 22 de las variantes de rpoB de los
grupos 1 y 2, entre ellas `p.Ser450Phe`, `p.Ser450Tyr`, `p.Ser441Gln` y
`p.His445Cys`, necesitan dos cambios en un mismo codón; leída base a base,
`p.Ser450Phe` es `p.Ser450Leu` más un cambio sinónimo. get_MNV marca cada una
con `codon reading only`.

La OMS también da, para cada variante clasificada, los cambios genómicos que la
producen: 26.231 en los grupos 1 y 2. Cada uno se anotó sobre el ancestro del
MTBC del ejemplo con la tabla incluida. 26.225 salen marcados con su propia
entrada, 198 de ellos solo por la lectura del codón, y 2 más con otra entrada
para el mismo fármaco, el `LoF` del gen. Los 4 restantes son 2 inserciones
justo antes de un codón de inicio, que la OMS cuenta como pérdidas de función y
get_MNV como intergénicas, y 2 cambios cuya consecuencia sobre el ancestro no
es la que tienen sobre H37Rv. Con la carpeta que guarda los dos
archivos de la OMS, `tests/who_catalogue/run.py` repite la comprobación.

Con tus propios datos:

- **Nombres de gen.** La OMS nombra los genes por sus símbolos y locus tags de
  H37Rv (`rpoB`, `katG`, `Rv0678`, `rrs`), y una entrada coincide con un gen de
  ese nombre o que lo tenga como una de sus partes separadas por `_`. Los 65
  genes del catálogo se encuentran así en `example/anot_genes.txt` y en la
  anotación de H37Rv del NCBI (GCF_000195955.2), cuyos genes get_MNV nombra
  `rpoB_Rv0667`. Con una anotación que los nombre de otra manera, cambia la
  columna `gene` de la tabla a los nombres de la anotación.
- **Referencia.** El catálogo está escrito contra H37Rv (NC_000962.3). Las
  entradas de proteína, `c.` y `n.` se numeran a lo largo de su gen, así que
  valen en cualquier referencia cuyos genes coincidan con los de H37Rv. Las
  entradas `g.` son posiciones de H37Rv: construye la tabla con tu propio FASTA
  en `--fasta` cuando esté en coordenadas de H37Rv, y deja fuera
  `--coordinates` cuando no lo esté. El ancestro del MTBC del ejemplo comparte
  las coordenadas de H37Rv, pero no los alelos que lleva el propio linaje de
  H37Rv: 38 entradas nombran un residuo o una base que el ancestro no tiene en
  ese punto, ninguna de los grupos 1 y 2 (`embB p.Glu378Ala`, un marcador de
  linaje, es una). Contra el ancestro no coinciden con nada, y una muestra con
  el alelo de H37Rv muestra el cambio inverso, `Ala378Glu`.

La tabla y los archivos de los que sale contienen información del
[catálogo de mutaciones de la OMS](https://github.com/GTB-tbsequencing/mutation-catalogue-2023),
que se distribuye bajo la
[licencia ODC Attribution](https://opendatacommons.org/licenses/by/1-0/); cítalo
como *Catalogue of mutations in Mycobacterium tuberculosis complex and their
association with drug resistance, second edition. Geneva: World Health
Organization; 2023*.

## Nombres de contig

Todos los archivos de entrada deben coincidir en los nombres de contig:

```text
variant contig == FASTA record ID == GFF sequence ID == BAM reference name
```

Por ejemplo, `chr1` y `1` son nombres diferentes. Renombra o normaliza las
entradas antes de ejecutar get_MNV.
