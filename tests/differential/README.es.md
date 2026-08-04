# Suite de anotación diferencial

Los tests propios de get_MNV solo pueden comprobar que produce lo que nosotros
apuntamos como esperado. Cuando la expectativa está mal, el test se pone verde y
la herramienta se queda mal.

No es hipotético. Un CDS unido por deslizamiento ribosómico, la forma con la que
se anota ORF1ab de SARS-CoV-2, recibía `splice_donor_variant` (`HIGH`) sobre
bases codificantes normales. El test unitario que debía cubrir el splicing
construía un gen "de dos exones" con los exones pegados, sin intrón entre ellos.
Pasaba, y confirmaba el bug.

Esta suite nos saca del bucle comparando get_MNV con anotadores que no hemos
escrito nosotros.

```bash
cargo build -p get_mnv --release
python3 tests/differential/run.py
```

## Los dos oráculos

| Dataset | Oráculo | Por qué |
|---|---|---|
| `slippage` | `bcftools csq`, en vivo | Un *join* de deslizamiento: segmentos CDS contiguos, sin intrón, sin nada empalmado. La forma exacta que produjo el bug. |
| `spliced` | `bcftools csq`, en vivo | Un transcrito de dos exones con splicing real. El control positivo: sin él, "no emitir nunca términos de splice" también pasaría. |
| `mtb` | SnpEff 4.1l | El VCF de *M. tuberculosis* incluido ya trae campos `ANN=` de la ejecución original. Datos reales, anotador real, 950 llamadas. |

En los fixtures ambas herramientas leen la **misma** referencia y el **mismo**
GFF3, así que una discrepancia nunca se puede achacar a que vieran anotaciones
distintas. En MTB cada herramienta usa su propia anotación del mismo genoma, así
que las diferencias de modelo de gen forman parte de lo que registra la línea
base.

Los fixtures reflejan cómo se anotan estas cosas en la práctica. En concreto, el
fixture de deslizamiento **no tiene filas `exon`**, porque el GFF3 de NCBI para
ORF1ab tampoco las tiene: es un gen y un CDS de dos segmentos con una excepción
de deslizamiento ribosómico. Inventar exones ahí cambiaría lo que reporta
`bcftools csq` y convertiría el fixture en un test de algo que no existe.

## La regla no es "tienen que coincidir"

get_MNV discrepa a propósito. Reporta MNVs a nivel de codón que los anotadores
por variante separan, que es su razón de existir. En el dataset de MTB se ve
directamente: en la posición 2282377 SnpEff reporta `synonymous_variant` (`LOW`),
porque mira esa variante sola. get_MNV ve que 2282376 está en el mismo codón, así
que el codón pasa de verdad de `GTT` a `GCC` y el residuo cambia de verdad a Ala.
get_MNV tiene razón, y esa diferencia no debe romper CI.

Por eso cada discrepancia se contrasta con `baseline.tsv`, una lista versionada
de diferencias aceptadas **con el motivo escrito**. La suite falla solo cuando
aparece una diferencia que nadie ha explicado todavía.

Ese fichero merece leerse por sí solo: es el recuento honesto, posición a
posición, de en qué se diferencia get_MNV de las herramientas establecidas y por
qué.

### Cuando la suite falla

O get_MNV tiene un bug, o el comportamiento cambió a propósito. Mira la
diferencia antes de tocar nada. Si es intencionada, añade una fila a
`baseline.tsv` explicándola. `--update` reescribe el fichero entero con lo que se
observa ahora, con `TODO explain` en cada nota; úsalo para arrancar y luego
escribe las notas a mano. Una fila de línea base sin motivo es peor que no
tenerla.

La suite también avisa de entradas de la línea base que ya no observa. Eso suele
significar que un arreglo eliminó una diferencia, y la fila debe borrarse.

## Comprobar que la suite todavía puede fallar

Un test que no puede fallar no vale nada. `GET_MNV_BIN` apunta la suite a otro
binario, para confirmar que sigue cazando el bug para el que se escribió:

```bash
git worktree add /tmp/prefix <commit-anterior-al-arreglo>
cargo build --release --manifest-path /tmp/prefix/Cargo.toml
GET_MNV_BIN=/tmp/prefix/target/release/get_mnv python3 tests/differential/run.py
```

Contra el commit anterior al arreglo del deslizamiento, esto sale con código 1 y
apunta directamente a las posiciones del límite donde get_MNV emitía
`splice_donor` y `bcftools csq` no.

## Ficheros

| Ruta | Qué es |
|---|---|
| `run.py` | El arnés: ejecuta ambas herramientas, normaliza, compara y contrasta la línea base |
| `baseline.tsv` | Diferencias aceptadas, cada una con su motivo |
| `fixtures/` | Los genomas sintéticos pequeños, versionados para que la ejecución sea hermética |
| `make_fixtures.py` | Regenera `fixtures/`; solo se ejecuta para cambiarlos o ampliarlos |

`bcftools` debe estar en el `PATH` para los datasets de fixtures; si falta se
saltan con un mensaje claro y la comparación de MTB se sigue ejecutando.
