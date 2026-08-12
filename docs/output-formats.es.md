# Formatos de salida

get_MNV puede escribir archivos de metadatos TSV, VCF, BCF y JSON.

## Salida TSV por defecto

Nombre de archivo por defecto:

```text
<input_name>.MNV.tsv
```

Usa este formato para hojas de cálculo, análisis posteriores e inspección rápida.

Columnas principales:

| Columna | Significado |
|---|---|
| `Chromosome` | Nombre del contig |
| `Gene` | Nombre del gen o de la característica. Las variantes intergénicas se marcan como `intergenic`. |
| `Positions` | Una posición para los SNP y varias posiciones separadas por comas para los MNV, en orden genómico ascendente en ambas hebras. Las columnas por SNV (`Reference Bases`, `Base Changes`, `SNP AA Changes`, `SNP Codon`, `Event Components` y las de soporte de lecturas) son arrays paralelos en ese mismo orden. |
| `Reference Bases` | Bases de referencia en esas posiciones. |
| `Base Changes` | Bases alternativas. |
| `AA Changes` | Cambio de aminoácido tras combinar todos los SNV del codón. |
| `SNP AA Changes` | Cambio de aminoácido para cada SNV considerado por separado. |
| `Local AA Changes` | Numeración local al exón, que es lo que get_MNV reportaba antes de la 1.1.2. Idéntica a `AA Changes` cuando la característica no tiene contexto de transcrito (un TSV de genes, un CDS procariota o de un solo exón, y el primer exón de un modelo empalmado), y distinta en los exones posteriores de un transcrito empalmado, donde `AA Changes` cuenta desde el inicio de la proteína y esta desde el inicio del exón. |
| `Local SNP AA Changes` | Cambios de aminoácido por SNP en numeración local. |
| `Variant Type` | `SNP`, `MNV`, `SNP/MNV` o `INDEL`. |
| `Change Type` | Sinónimo, no sinónimo, codón de parada ganado/perdido, desconocido, etc. |
| `Reference Codon` | Codón original, en la orientación del propio transcrito: en un gen de hebra menos son las bases de la hebra codificante, no las genómicas, de modo que el codón siempre traduce al aminoácido que aparece a su lado. |
| `SNP Codon` | Codón con las sustituciones SNP individuales, misma orientación que `Reference Codon`. |
| `MNV Codon` | Codón con todas las sustituciones agrupadas, misma orientación que `Reference Codon`. |
| `Event Class` | Clase canónica del evento de alelo: `snp`, `mnv`, `insertion`, `deletion`, `delins`, `complex_indel` o `symbolic`. |
| `Event Components` | Descomposición REF/ALT como `SNV:10:A>G`, `INS:10:+T` o `DEL:11-12:TG`. |
| `SO Term` | Término de consecuencia Sequence Ontology (`missense_variant`, `synonymous_variant`, `stop_gained`, `start_lost`, `frameshift_variant`, `inframe_deletion`, `intergenic_variant`, …). Las variantes cerca de una unión exón-exón interna de un transcrito con splicing llevan además un término de splice: `splice_donor_variant` / `splice_acceptor_variant` (las dos bases intrónicas esenciales de cada extremo del intrón, `HIGH`) o `splice_region_variant` (las 3 primeras o las 3 últimas bases del exón, o las bases 3ª-8ª del intrón, `LOW`). Un cambio codificante exónico cerca de una unión se combina, p. ej. `missense_variant&splice_region_variant`. Una unión exige un intrón real: los segmentos CDS contiguos o solapados (una unión por deslizamiento ribosómico, como ORF1ab de SARS-CoV-2) forman una pauta de lectura continua y no llevan términos de splice. Una variante dentro de un intrón pero lejos de sus sitios de splice es `intron_variant` (`MODIFIER`), reportada contra su gen y no como intergénica, y una variante en una feature declarada no codificante es `non_coding_transcript_exon_variant` (`MODIFIER`), sin cambio de aminoácido. |
| `Impact` | Impacto previsto según convenciones SnpEff/VEP: `HIGH`, `MODERATE`, `LOW` o `MODIFIER`. Una consecuencia combinada splice/codificante conserva el impacto más severo. |
| `Grantham` | Distancia de Grantham y categoría de conservación de un cambio missense (p. ej. `177 (radical)`); `-` si es sinónimo, sin sentido o no codificante. |
| `MNV Consequence Shift` | Cómo compara el MNV combinado con sus SNV individuales: `MNV-gained` (más severo que cualquier SNV solo — lo que se pierden los anotadores por-SNV), `MNV-masked` (un SNV sin sentido rescatado por su vecino) o `Concordant`. `-` para SNV individuales. |
| `DBS Class` | Clase de sustitución de doblete (DBS) al estilo COSMIC para un MNV de dos sustituciones de una base adyacentes, p. ej. `CC>TT` (colapsado por complemento reverso, así que `GG>AA` se reporta como `CC>TT`). `-` para SNV individuales, indels y MNV no adyacentes o de 3 SNV. |
| `Declared Phase` | La fase que declaró el **VCF de entrada** para los alelos de esta fila, como `cis:12345` (veredicto y conjunto de fase `PS`) o solo `cis` si el llamador faseó sin conjunto. Se lee de un `GT` separado por `|`; un genotipo con `/` es el llamador diciendo que no resolvió la fase, y sale como `-`. Se añade `|contradicted-by-reads` cuando el BAM no le deja sitio: un cis declarado que ni una sola lectura completa lleva, o un trans declarado que llevan todas las lecturas informativas. Es la afirmación del llamador, no una observación; las columnas de fasing de al lado son la evidencia. `-` para filas de una sola posición y para entrada sin fasear. |
| `MNV Phasing Support` | Soporte de fasing (ligamiento) derivado del BAM: entre las lecturas que observan *todas* las posiciones del codón y llevan el SNV constituyente menos soportado, la fracción que además lleva el haplotipo MNV completo. `1.0000` = co-ocurrencia perfecta (un haplotipo real); valores bajos indican que los SNV caen mayoritariamente en moléculas distintas (una coincidencia en el mismo codón, no un MNV real). `0.0000` es un hallazgo: hubo lecturas que cruzaron el codón y ninguna llevaba ambos. `-` significa que no se pudo responder: sin `--bam`, un SNV suelto, o ninguna lectura que alcance de un extremo a otro del codón (habitual cuando un codón cruza un intrón y los fragmentos son más cortos que él). |
| `Haplotype LD` | Desequilibrio de ligamiento (`D'`) a nivel de lectura entre las variantes que esta fila afirma que van juntas, sobre las moléculas que las observaron. Cubre los dos tipos de fila con varias variantes: un MNV de codón y un haplotipo local de indel. Un ratio de co-ocurrencia no puede separar un haplotipo de una casualidad de frecuencias: dos sustituciones que están cada una en el 90% de las moléculas se encuentran juntas en el 81% por pura aritmética, y el ratio eso lo llama 0.9. `D'` mide el exceso sobre lo que predicen esas dos frecuencias, normalizado por el máximo que podría haber sido. `+1` = van juntas todo lo que sus frecuencias permiten, así que el MNV es un haplotipo real. `~0` = co-ocurren exactamente lo que predice el azar, así que solo comparten codón. `-1` = se excluyen: ambas presentes, nunca en la misma molécula, lo que en una población haploide son dos linajes en competencia y no una variante. Con tres o más variantes decide el par más débil, porque la fila afirma que una molécula las lleva todas. Responde a una pregunta distinta que el recuento de al lado: el recuento es cuántas moléculas *son* esa combinación, mientras que `D'` es si sus variantes co-ocurren más de lo que predicen sus propias frecuencias, así que un haplotipo con pocas moléculas puede estar perfectamente ligado y uno con muchas puede ser una casualidad. `-` cuando ninguna molécula observó las variantes juntas, o cuando una de ellas está en todas esas moléculas o en ninguna, que no deja nada que correlacionar. Solo aparece con `--bam`. |
| `Haplotype LD p` | Valor p exacto de Fisher a dos colas para esa tabla, para que un `D'` de 1.0 con cuatro moléculas no se lea igual que uno con cuatrocientas. `-` en los mismos casos que la columna anterior. Solo aparece con `--bam`. Mira [Ligamiento](linkage.es.md). |
| `MNV Phasing Reads` | Sobre cuántas lecturas se calculó `MNV Phasing Support`, para que un `1.0000` con 3 lecturas no se lea igual que uno con 300. Sale `-` exactamente cuando esa columna sale `-`, que no es el mismo criterio que el de las dos columnas de ligamiento de al lado. Solo aparece con `--bam`. |
| `Frameshift Phasing` | Lo que dijeron las lecturas sobre si este codón comparte moléculas con cada indel aguas arriba, como `trans:1234:0/18`: el veredicto, la posición del indel, y las lecturas en cis sobre las que podían responder. Varios se unen con ` | `. Un codón sin la marca de frameshift se ve igual tanto si las lecturas probaron que el indel está en otras moléculas como si nadie preguntó; `-` es el segundo caso. Solo aparece con `--bam`. |
| `NMD Prediction` | Predicción de decaimiento mediado por sin sentido (NMD) para un stop prematuro según la regla de los 50 nt: `NMD-triggering` cuando el PTC está a más de 50 nt aguas arriba de la última unión exón-exón, `NMD-escaping` cuando está en el último exón o a menos de 50 nt de esa unión. `-` para variantes sin stop prematuro y para transcritos sin unión exón-exón. Requiere un modelo de CDS con splicing (transcrito GFF/GTF) cuyos segmentos estén separados por intrones reales; un único segmento CDS, o segmentos unidos por deslizamiento ribosómico, no tiene unión. |
| `HGVS g.` | Descriptor HGVS genómico: `g.100A>G` para un SNV, el bracket de alelo `g.[28G>T;30T>A]` para un MNV, y `g.101_102del` / `g.100_101insTG` / `g.101delinsC` para indels. No hace 3'-shifting (usa la posición del alelo de entrada) y no lleva prefijo de accesión de referencia. |
| `HGVS c.` | Descriptor HGVS codificante para una sustitución codificante, numerado desde el inicio del CDS con bases de la hebra codificante: `c.30A>G` (SNV) o el bracket de alelo `c.[28G>A;30T>C]` (MNV). `-` para indels y variantes no codificantes; el cambio proteico (`p.`) en las columnas AA da la consecuencia de los indels. |

Columnas adicionales cuando se usa `--bam`:

| Columna | Significado |
|---|---|
| `SNP Reads` | Lecturas que respaldan cada SNV individual. |
| `SNP Forward/Reverse Reads` | Soporte de SNP por hebra. |
| `MNV Reads` | Lecturas que respaldan el haplotipo MNV completo. |
| `MNV Forward/Reverse Reads` | Soporte de MNV por hebra. |
| `Total Reads` | Profundidad en las posiciones de la variante. |
| `SNP Frequencies` | Frecuencias de SNP por posición. |
| `MNV Frequencies` | Frecuencia del haplotipo MNV. |
| `Event Reads` | Lecturas exactas que respaldan un evento indel/complejo. |
| `Event Forward/Reverse Reads` | Soporte exacto del evento por hebra. |
| `Event Depth` | Lecturas con un alelo observado a lo largo del tramo del evento indel/complejo. |
| `Event Frequency` | Lecturas exactas del evento divididas por la profundidad del evento. |

El soporte exacto del evento tiene en cuenta el CIGAR. Una lectura debe reconstruir la misma
secuencia ALT local y, en el caso de los haplotipos complejos, contener los componentes de inserción
y deleción esperados. Esto evita que las combinaciones de inserción/deleción de efecto neutro
se contabilicen como soporte solo porque su secuencia se parezca a un MNV.

Las columnas de frecuencia se calculan a partir del soporte de lecturas del BAM. `--min-snp-frequency` y
`--min-mnv-frequency` usan estos mismos valores derivados del BAM. Los filtros son
independientes: `--min-snp-frequency` se aplica a las observaciones de SNP individuales y
`--min-mnv-frequency` se aplica a los haplotipos MNV en fase. En las llamadas mixtas `SNP/MNV`,
una fila o registro VCF se conserva cuando cualquiera de los dos componentes supera su propio
umbral activo.
Los filtros de recuento de lecturas y de soporte de hebra (`--snp`, `--mnv`, `--min-snp-strand`
y `--min-mnv-strand`) siguen el mismo comportamiento independiente para SNP y MNV.

Cuando un MNV a nivel de codón se solapa con un indel, la fila del MNV se conserva como fila
de contexto posicional, pero su efecto a nivel de aminoácido se marca como `Unknown` con
`Change Type = Indel overlap`. Si las lecturas del BAM respaldan el evento combinado completo,
get_MNV emite una fila `complex_indel` exacta aparte con el REF/ALT combinado, los
componentes del evento y el soporte de lecturas del evento.

El solapamiento del indel sigue la semántica interbase de VCF. Las deleciones se solapan con una
característica a lo largo del tramo de referencia eliminado. Las inserciones solo se solapan con una
característica cuando la secuencia insertada cae entre dos bases de referencia situadas dentro de ella, de
modo que una inserción anclada en la última base de la característica se reporta fuera de esta.

Ejemplo:

```text
Chromosome	Gene	Positions	Base Changes	AA Changes	Variant Type	Change Type
MTB_anc	Rv0095c	104838	T	Asp126Glu	SNP	Non-synonymous
MTB_anc	Rv0095c	104941,104942	T,G	Gly92Gln	SNP/MNV	Non-synonymous
```

## Salida VCF

Escribe VCF con:

```bash
--convert
```

o escribe tanto TSV como VCF con:

```bash
--both
```

Nombre de archivo por defecto:

```text
<input_name>.MNV.vcf
```

Usa `--vcf-gz` para obtener salida comprimida:

```text
<input_name>.MNV.vcf.gz
```

Campos INFO comunes:

| Campo | Significado |
|---|---|
| `GENE` | Nombre del gen o de la característica |
| `AA` | Cambio de aminoácido |
| `CT` | Tipo de cambio |
| `TYPE` | Tipo de variante |
| `EC` | Clase canónica del evento de alelo |
| `COMP` | Componentes del evento REF/ALT |
| `ODP` | Profundidad original del archivo de variantes de entrada |
| `OFREQ` | Frecuencia alélica original del archivo de variantes de entrada |
| `SR`, `SRF`, `SRR` | Lecturas de SNP: total, forward, reverse |
| `MR`, `MRF`, `MRR` | Lecturas de MNV: total, forward, reverse |
| `DP` | Profundidad recalculada a partir del BAM |
| `FREQ` | Frecuencia recalculada a partir del BAM |
| `ER`, `ERF`, `ERR` | Lecturas exactas del evento indel/complejo: total, forward, reverse |
| `EDP` | Profundidad exacta del evento para alelos indel/complejos |
| `EFREQ` | Frecuencia exacta del evento para alelos indel/complejos |
| `SBP` | Valor p del sesgo de hebra de SNP |
| `MSBP` | Valor p del sesgo de hebra de MNV |
| `SO`, `IMPACT` | Término de consecuencia Sequence Ontology e impacto predicho |
| `GD` | Distancia de Grantham de un cambio missense |
| `MNVSHIFT` | Consecuencia del MNV combinado frente a sus SNV individuales |
| `DBS` | Clase de doblete estilo COSMIC para MNV de 2 SNV adyacentes (p. ej. `CC>TT`) |
| `MNVPS` | Soporte de fasing del MNV (de las lecturas que cruzan el codón y llevan el SNV limitante, la fracción que lleva el haplotipo completo) |
| `MNVPR` | Lecturas sobre las que se calculó ese ratio |
| `FSPH` | Fase con cada indel aguas arriba, como veredicto:posición:cis/informativas |
| `DPHASE` | Fase que declaró el VCF de entrada para esta fila, veredicto:conjunto |
| `LD` | Desequilibrio de ligamiento D-prime entre las sustituciones del codón |
| `LDP` | Valor p exacto de Fisher para esa tabla de ligamiento |
| `NMD` | Predicción de NMD para un stop prematuro (regla de los 50 nt) |
| `HGVSG` | Descriptor HGVS genómico (bracket de MNV con `;` percent-encoded) |
| `HGVSC` | Descriptor HGVS codificante para una sustitución (`;` percent-encoded) |

La cabecera del VCF registra la versión de get_MNV, la línea de comandos y los umbrales usados.
Cuando `--emit-filtered` está habilitado, los registros VCF que quedan por debajo de los umbrales de
soporte de lecturas, frecuencia, soporte de hebra o sesgo de hebra se escriben con etiquetas FILTER como
`LowSupport`, `LowFrequency`, `StrandSupport` o `StrandBias`; de lo contrario, se omiten.

## Salida BCF

Escribe BCF con:

```bash
--bcf
```

BCF requiere el modo de salida VCF, así que úsalo con `--convert` o `--both`.
Esto es solo conversión de salida; BCF no se acepta como formato de entrada.

Nombre de archivo por defecto:

```text
<input_name>.MNV.bcf
```

## Archivos JSON

### JSON de resumen

Escribe con:

```bash
--summary-json run.summary.json
```

Incluye:

- Sumas de comprobación de los archivos de entrada
- Recuentos de variantes por contig
- Recuentos globales de variantes
- Tiempos de ejecución
- Rutas de salida

### Manifiesto de la ejecución

Escribe con:

```bash
--run-manifest run.manifest.json
```

Incluye el resumen y, además:

- Línea de comandos
- Versión de la herramienta
- Sumas de comprobación de los archivos de salida
- Marca de tiempo

### JSON de errores

Escribe los errores como JSON con:

```bash
--error-json run.error.json
```

Esto resulta útil en pipelines automatizados.

## Informe HTML interactivo

`--report <FICHERO.html>` escribe un único fichero HTML autocontenido para explorar
las variantes llamadas. Incrusta sus datos y no carga scripts ni fuentes externas,
así que se abre sin conexión con un doble clic y se puede adjuntar en un correo o
archivar junto a los resultados.

[**Abrir el informe de ejemplo**](assets/example-report.html){ target=_blank }: 941
variantes llamadas sobre el conjunto de una muestra incluido en `example/`, o sea la
salida real del comando de abajo y no una maqueta. Se regenera con
`scripts/build_example_report.sh`. Las vistas de cohorte (las filas de la matriz, la
recurrencia de haplotipos) se llenan con `--sample all` o `--report-from` sobre varias
muestras.

```bash
# Informe de un run
get_mnv --vcf muestra.vcf --fasta ref.fasta --gff ref.gff --report muestra.html

# Un informe que cubre todas las muestras de un VCF multi-muestra
get_mnv --vcf cohorte.vcf --fasta ref.fasta --gff ref.gff --sample all --report cohorte.html

# Cohorte procesada muestra a muestra: agrega los TSV al final
get_mnv --report-from resultados/*.MNV.tsv --report cohorte.html
```

Con `--report-from` no se ejecuta el pipeline: el informe se construye a partir de
TSV de get_MNV ya existentes, que es la forma habitual en un flujo de Nextflow o
Snakemake que llama a get_MNV una vez por muestra. Cada fichero es una muestra,
etiquetada con su nombre de archivo sin el sufijo `.MNV.tsv` (`TB-001.MNV.tsv`
pasa a ser `TB-001`).

El informe contiene:

- **Una cifra principal** de las variantes mostradas, con tarjetas de apoyo para
  muestras, genes, filas MNV y variantes de alto impacto, y la **distribución de
  consecuencias** por término de Sequence Ontology. Todo sigue los filtros activos.
- **Una matriz de variantes**, dibujada como un pequeño navegador genómico: las
  muestras en el lateral, coordenadas genómicas reales arriba y color por la base
  alternativa. Rueda para hacer zoom alrededor del cursor, arrastrar para
  desplazarse, arrastrar sobre la tira del contig completo para saltar, y doble
  clic para volver al inicio. Arrastrar a lo largo de la tira del contig completo
  selecciona un rango; un clic simple sobre ella recentra. La caja de región acepta
  `inicio-fin`, `contig:inicio-fin`, una coordenada suelta o un nombre de gen, y lo
  encuadra. Sobre la regla, una pista de densidad muestra las llamadas por bin
  genómico en carriles independientes: uno por cada métrica (todas las llamadas,
  SNP, MNV, indel, impacto HIGH y número de muestras distintas con llamada). Cada
  carril tiene su propia escala y su propio máximo, así una clase rara se lee al
  lado de una común en vez de quedar aplastada en la base de un apilado, y el
  control `Tracks` alterna entre todos los carriles, un conjunto compacto y
  ninguno. Una pista de genes marca la extensión de los sitios llamados
  de cada gen (que es donde hubo llamadas, no el límite anotado del gen, que el
  informe no lleva). Al
  acercarse lo suficiente aparece la letra de la base dentro de cada celda. Las
  muestras se pueden ordenar por perfil compartido (los patrones idénticos quedan
  juntos), por número de variantes o por nombre, y cada etiqueta de muestra lleva
  una barra con cuántas llamadas tiene en la ventana visible, por eso no hay un
  gráfico aparte por muestra. El color significa una sola cosa: los tonos de
  nucleótido son de las celdas de la matriz, toda magnitud va en neutro, y el
  color de estado reservado marca el impacto HIGH. Las posiciones llamadas juntas en
  las mismas lecturas (MNV de codón, indels complejos faseados) llevan una marca
  sobre sus columnas, y aparece un selector de contig con datos multi-contig, ya
  que un eje de coordenadas continuo no puede abarcar varios contigs.

    Una celda solo puede ser "llamada ALT" o **"no llamada"**. La salida de get_MNV
    no distingue una base de referencia de una posición sin cobertura, así que una
    celda vacía nunca se presenta como referencia. Ese estado se lee como
    "referencia o sin cobertura".

- **Haplotipos respaldados por lecturas**: las combinaciones de alelos que get_MNV
  observó realmente en las mismas lecturas, es decir MNV de codón e indels complejos
  faseados localmente, ordenados por cuántas muestras llevan cada uno y con su
  soporte de fasing cuando se usó un BAM. El fasing de largo alcance entre sitios
  distantes no se infiere ni se muestra nunca.
- **Una tabla ordenable y filtrable** con todas las variantes, virtualizada para que
  decenas de miles de filas sigan siendo fluidas. Cada columna tiene su propio
  filtro: lista de casillas para muestra, contig, gen, tipo de variante,
  consecuencia e impacto (varios valores a la vez, con buscador cuando la lista es
  larga), coincidencia por texto para posición, cambio de bases y cambio de
  aminoácido, y rango mínimo/máximo para Grantham y frecuencia. Los filtros de
  columna se combinan entre sí y con el buscador libre, y gobiernan toda la página:
  las cifras principales, los gráficos, la matriz y el panel de haplotipos siguen
  la misma selección.
- **Un panel de detalle** de la variante seleccionada con su localización,
  consecuencia, descriptores HGVS `g.`/`c.`/`p.`, distancia de Grantham, clase DBS,
  predicción de NMD, codones, componentes del evento y soporte de lecturas.
- **Exportar TSV filtrado**, que descarga exactamente las filas mostradas.
- **Enlaces al repositorio y a la documentación** en la cabecera. Son hiperenlaces
  normales: no se descarga nada al abrir el informe, así que sigue siendo
  autocontenido sin conexión.

El informe sigue el tema claro u oscuro del sistema operativo y tiene su propio
selector. Como se construye desde el TSV, `--report` necesita salida TSV: funciona
con el modo de salida por defecto y con `--both`, pero no con `--convert` a solas
ni con `--dry-run`. Las columnas de soporte de lecturas (frecuencia, profundidad)
solo se rellenan si el run usó `--bam`.

El tamaño escala con el número de variantes. Los campos repetidos se codifican con
diccionario, así que una cohorte de decenas de miles de variantes se queda en pocos
megabytes.

## Notas

- Para los registros MNV, la profundidad y la frecuencia se calculan a partir de las lecturas que
  abarcan todas las posiciones del haplotipo agrupado.
- Las frecuencias se imprimen con 4 decimales.
- `--min-snp-frequency` y `--min-mnv-frequency` son valores de `0` a `1`
  y requieren `--bam`.
- Los filtros de frecuencia de SNP y MNV son independientes, de modo que un haplotipo MNV fuerte
  no se elimina por un umbral de frecuencia de SNP más estricto.
- Los filtros de soporte de lecturas y de soporte de hebra de SNP y MNV también son independientes.
- `--sample all` escribe un conjunto de salida por cada muestra del VCF.
- `--keep-original-info` conserva los campos INFO del VCF de entrada que no son de get_MNV.
