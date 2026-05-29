# Semántica de indels y MNV

get_MNV es un anotador de variantes y un resumidor de haplotipos para llamadas
de variantes ya existentes. No reemplaza a un llamador de variantes: el VCF o el
TSV de iVar de entrada sigue siendo el responsable de decidir qué alelos existen.
Cuando se proporciona un BAM, get_MNV usa la evidencia de las lecturas para
contar el soporte y para emitir haplotipos combinados exactos de indel/SNV/MNV
que se observan en las mismas lecturas.

## Compatibilidad con llamadores

Distintos llamadores pueden representar el mismo haplotipo local de maneras
diferentes. Por ejemplo, FreeBayes se basa en haplotipos y puede emitir SNPs,
indels, MNPs y eventos complejos como llamadas de haplotipo en bruto. Su propia
documentación recomienda descomponer esas llamadas después del llamado si se
necesitan registros primitivos de SNP/indel.

bcftools `norm` se ocupa del otro lado del problema: puede alinear a la izquierda
y normalizar indels con una referencia FASTA, separar sitios multialélicos y
atomizar variantes complejas o MNV en registros más simples. Ese paso de
normalización es útil antes de comparar conjuntos de llamadas de distintos
llamadores.

GATK HaplotypeCaller es un llamador, no un anotador. Usa ensamblaje local y
parámetros de genotipado como la ploidía de la muestra, el número máximo de
haplotipos y la distancia de fusión de MNP. get_MNV actualmente no realiza este
descubrimiento de haplotipos de novo ni el modelado de verosimilitud de
genotipos.

## Qué aporta get_MNV

- Agrupación de SNVs según el codón en filas SNP, MNV y mixtas `SNP/MNV`.
- Normalización de inserciones/deleciones de iVar en alelos anclados compatibles
  con VCF.
- Efectos proteicos de indels y alelos complejos usando el feature seleccionado,
  la hebra, la fase del GFF y el offset proteico del transcrito.
- Soporte derivado del BAM para observaciones individuales de SNP/MNV y para
  eventos exactos de indel o indel complejo.
- Filas adicionales `complex_indel` solo cuando el haplotipo combinado completo
  de indel/SNV/MNV se observa en las lecturas del BAM.
- Ventanas de haplotipo local acotadas que pueden combinar varios eventos
  cercanos, como haplotipos de inserción más deleción, cuando existe soporte
  exacto de lecturas.

## Reglas de límites

Las inserciones de VCF son eventos interbase anclados en la base de referencia
anterior. Una inserción se considera que solapa un feature únicamente cuando la
secuencia insertada cae entre dos bases de referencia dentro de ese feature. Una
inserción anclada en la última base del feature queda fuera de ese feature.

Las deleciones de VCF están ancladas, pero su efecto biológico es el tramo de
referencia eliminado. Una deleción anclada justo antes de un CDS/gen sigue
afectando a ese feature si las bases eliminadas lo solapan.

Cuando una fila MNV solapa un indel, get_MNV mantiene la fila MNV como contexto
posicional y marca su efecto en el aminoácido como `Unknown` con
`Change Type = Indel overlap`. Si el BAM respalda el evento combinado completo,
get_MNV emite una fila exacta `complex_indel` aparte.

El soporte exacto para indels complejos tiene en cuenta el CIGAR. Una lectura
debe producir la misma secuencia ALT local y contener los componentes esperados
de inserción/deleción. Esto importa para combinaciones netamente neutras, donde
una inserción más una deleción pueden producir la misma secuencia que un MNV
simple bajo otro alineamiento.

## Notas sobre CDS eucariota

Para anotaciones eucariotas GFF/GTF, usa `--gff-features CDS` para que get_MNV
pueda usar la fase del GFF y los identificadores de transcrito. Cuando las filas
CDS llevan `transcript_id` o `Parent`, get_MNV construye un modelo de CDS
empalmado para cada transcrito y evalúa la agrupación por codones, los efectos
del MNV en el aminoácido y el contexto de frameshift (cambio del marco de
lectura) del indel sobre esa secuencia codificante completa.

Esto significa que los codones repartidos entre uniones de exones se pueden
anotar como un único MNV a nivel de transcrito, y los SNPs corriente abajo se
marcan como frameshifted solo cuando el desplazamiento neto del indel codificante
corriente arriba permanece fuera de marco. Las filas CDS sin un identificador de
transcrito utilizable mantienen el comportamiento anterior por feature.

## Limitaciones actuales

- No se acepta un BCF de entrada directamente. Conviértelo primero con
  `bcftools view`.
- La alineación a la izquierda de indels y la normalización completa no se
  realizan automáticamente a menos que `--normalize-alleles` pueda recortar un
  prefijo o sufijo compartido. Para comparaciones entre llamadores, normaliza las
  entradas primero con una herramienta que tenga en cuenta el FASTA, como
  `bcftools norm -f ref.fa`.
- Los genotipos, la ploidía, los conjuntos de fase y las verosimilitudes de
  genotipo no se vuelven a estimar. Para datos eucariotas diploides/poliploides,
  interpreta con cautela los MNV heterocigotos sin fase a nivel de exón a menos
  que la evidencia del VCF o del BAM confirme que los alelos están en el mismo
  haplotipo.
- El ensamblaje local de novo queda fuera del alcance. get_MNV solo combina
  alelos ya presentes en la entrada y, cuando está disponible, confirmados por el
  soporte de lecturas del BAM.
- El descubrimiento de haplotipos locales está acotado a ventanas de eventos
  cercanos; la reconstrucción de haplotipos muy grandes debería seguir siendo
  manejada por un llamador dedicado.
- El ensamblaje local completo sigue fuera del alcance: get_MNV reanota los
  alelos presentes en el archivo VCF/iVar de entrada en lugar de descubrir nuevas
  variantes candidatas.

## Parámetros de ajuste

Estos flags opcionales cambian el comportamiento que tiene en cuenta los indels.
Todos toman por defecto el comportamiento histórico, así que los pipelines
existentes no se ven afectados a menos que se establezca un flag.

- `--frameshift-min-freq <0.0-1.0>` (por defecto `0.0`): frecuencia alélica
  mínima que un indel *corriente arriba* debe alcanzar antes de que desplace el
  marco de lectura de los codones SNV/MNV corriente abajo (el marcador `(fs)` y
  los change types de frameshift). El valor por defecto propaga desde cada indel
  sin importar la frecuencia. Para poblaciones intrahospedador / mixtas, subir
  esto evita reetiquetar una sustitución de alta frecuencia corriente abajo como
  frameshifted a causa de un indel de baja frecuencia corriente arriba que casi
  con seguridad está en otra molécula. La propagación del frameshift es
  posicional, no por fase de lecturas, así que es una protección tosca pero útil.
- `--indel-anchor-depth` (por defecto desactivado): cuenta la profundidad del
  locus del indel (el denominador de `EDP`/`EFREQ`) a partir de las lecturas que
  observan la base de anclaje, en lugar de solo las lecturas que abarcan por
  completo el alelo REF. Reduce el subconteo de profundidad y el sesgo de `EFREQ`
  para deleciones de múltiples bases, donde de otro modo las lecturas que solapan
  parcialmente quedan excluidas del denominador.
- `--phased-indel-min-reads <N>` (por defecto `1`) y
  `--phased-indel-min-freq <0.0-1.0>` (por defecto `0.0`): soporte mínimo del BAM
  que una fila de haplotipo indel/complejo con fase debe tener para ser emitida.
  Las ventanas locales pueden enumerar muchos sub-haplotipos solapados; subir
  estos valores suprime las filas de baja confianza provenientes de agrupaciones
  densas de variantes.

Cuando un indel que tiene cobertura de lecturas produce cero soporte de
CIGAR-exacto, get_MNV emite una advertencia: esto normalmente significa que el
indel de entrada no está alineado a la izquierda de la misma forma que el BAM
(común en homopolímeros/repeticiones en tándem). Normaliza la entrada primero con
una herramienta que tenga en cuenta el FASTA (`bcftools norm -f ref.fa`) para que
el alelo coincida con el alineamiento de las lecturas.

## Referencias

- FreeBayes: https://github.com/freebayes/freebayes
- bcftools norm: https://samtools.github.io/bcftools/bcftools.html#norm
- GATK HaplotypeCaller: https://gatk.broadinstitute.org/hc/en-us/articles/30332006386459-HaplotypeCaller
