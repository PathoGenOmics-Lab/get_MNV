# Alineamiento a la izquierda de indels contra bcftools norm

[English](README.md) · **Español**

`--left-align` lleva cada inserción y deleción a la posición más a la izquierda
de su repetición que describe el mismo cambio, y recorta las bases que REF y ALT
comparten: la normalización que hace `bcftools norm -f ref.fa`. Los tests
unitarios fijan unos veinticinco registros a lo que bcftools escribe para
ellos, y un test de propiedades compara el resultado con una búsqueda por todas
las colocaciones del indel. Los dos los escribió la misma mano que el código.
Esta suite se lo pregunta a bcftools.

```bash
cargo build -p get_mnv --release
python3 tests/leftalign/run.py
```

## Qué genera

Cada registro cae en una repetición del genoma de *M. tuberculosis* incluido
(`example/MTB_ancestor.fas`), a 60 bases como mínimo de cualquier otro, a partir
de una semilla aleatoria fija:

- **Indels en repeticiones cortas**: homopolímeros y repeticiones en tándem de
  2 a 6 bases, borrando o insertando unidades enteras o cualquier otra longitud,
  en un sitio al azar de la repetición o en su extremo derecho, que es donde lo
  pone una descripción HGVS desplazada hacia 3'. Casi todas las inserciones
  continúan la repetición; algunas son bases al azar.
- **Indels de una unidad entera en minisatélites**, con unidades de 7 a 100
  bases, las únicas repeticiones de este genoma lo bastante largas para mover un
  registro cientos de bases.
- **Registros complejos**: una deleción y una inserción en un mismo registro,
  que bcftools recorta pero nunca mueve.
- **Los dos extremos del contig**, donde cambia la regla del ancla.

Cada uno se escribe de una de tres maneras: anclado en la base anterior al
cambio, en la posterior, o con contexto a los dos lados.

## Qué compara

El mismo VCF va a `bcftools norm -f` y a get_MNV con `--left-align`, una vez
como `.vcf` plano y otra como `.vcf.gz`, porque los dos los leen lectores
distintos. Cada registro lleva un identificador en INFO, que get_MNV conserva
con `--keep-original-info`, así que su VCF de salida se empareja con la entrada
registro a registro. Un registro pasa cuando get_MNV escribe el POS, REF y ALT
que escribe bcftools. Los registros pasan además por `--sample all` como una
cohorte de dos muestras, en los dos contenedores, porque una cohorte se lee de
memoria por un camino propio.

La ejecución cuenta además cuántos registros escribiría get_MNV de otra forma
sin `--left-align`, y falla si no hay ninguno, porque una comparación que nada
puede hacer fallar no demuestra nada.

## La diferencia deliberada

bcftools lee un código de ambigüedad de la referencia como su primera base (`R`
como `A`), mueve un indel a través de él y escribe esa base en REF, que entonces
no coincide con el FASTA. get_MNV detiene el indel en el código. El genoma de
ejemplo solo tiene A, C, G y T, así que la suite nunca se encuentra uno; el test
unitario `an_ambiguity_code_in_the_reference_stops_the_walk` fija ese
comportamiento.

## Comprobar que todavía puede fallar

Una compilación con la regla del inicio del contig rota (`start <= 1` en lugar
de `start == 0` en `src/io/left_align.rs`) informa de los dos registros que
llegan a la primera base:

```
[indels.vcf] edge_start_del (deletion, period 1, right): 2 TG>G -> bcftools 1 TT>T, get_MNV 2 TG>G
[indels.vcf] edge_start_ins (insertion, period 1, left): 2 T>TT -> bcftools 1 T>TT, get_MNV 2 T>TT
```

Apunta la suite a cualquier otro binario con `GET_MNV=/ruta/a/get_mnv`.
