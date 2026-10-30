# Conversor del catálogo de la OMS

[English](README.md) · **Español**

`tools/who_catalogue_to_tsv.py` convierte el catálogo de mutaciones del
complejo *M. tuberculosis* de la OMS (segunda edición, 2023) en una tabla para
`--catalogue` de get_MNV, y `example/who_catalogue_2023.tsv` es su salida para
los grupos 1 y 2. Un conversor puede equivocarse de una manera que su salida no
deja ver: una entrada que get_MNV lee pero que nunca coincide parece una más.
Así que esta batería pasa las entradas por las propias variantes genómicas de
la OMS y le pregunta a get_MNV si coinciden.

```bash
cargo build -p get_mnv
python3 tests/who_catalogue/run.py
```

## Qué comprueba

- **Conversión.** `fixtures/master_rows.tsv` guarda filas reales del archivo
  maestro, una por cada manera en que se reescribe una entrada: un cambio
  proteico, un cambio que necesita dos sustituciones en un codón (rpoB
  `p.Ser450Phe`, `p.Ser441Gln`), un corrimiento de marco, `LoF`, un inicio
  perdido, una parada perdida, una deleción y una duplicación en fase, una
  deleción-inserción de la misma longitud, cambios de promotor y de ARNr, un
  indel fuera de una secuencia codificante, una forma que get_MNV no escribe y
  una deleción del gen entero, que no se puede convertir. Tienen que dar
  exactamente las entradas esperadas, con y sin las coordenadas genómicas, y la
  deleción tiene que aparecer en la lista con su motivo.
- **Coincidencia.** `fixtures/coordinates_rows.tsv` guarda las variantes
  genómicas de la OMS para esas filas. Cada una va a un VCF sobre el ancestro
  del MTBC del ejemplo, como mucho una por gen y a 300 pb de cualquier otra,
  para que ninguna comparta codón con otra ni le corra el marco, y se anota con
  la tabla convertida. Cada una tiene que salir marcada con la entrada a la que
  la OMS la asocia, y los codones de dos cambios solo por la lectura del codón.
  Una no puede: la OMS escribe el cambio de `Rv0565c p.TrpLeu266*` sobre los
  codones 266 y 267, y get_MNV da cada codón en su propia fila, así que ninguna
  fila lleva las dos sustituciones que nombra el `g.` de la entrada. La batería
  espera ese fallo y falla si empieza a coincidir.
- **La tabla incluida.** `example/who_catalogue_2023.tsv` lleva el aviso de la
  licencia, coincide entrada a entrada con lo que escribe el conversor y se
  carga en get_MNV; la muestra de ejemplo G35894 no lleva ninguna variante de
  resistencia clasificada y no marca nada.

## Con los archivos de la OMS

Descarga los dos archivos del
[repositorio de la OMS](https://github.com/GTB-tbsequencing/mutation-catalogue-2023)
(`Final Result Files`) en una carpeta y pásala:

```bash
cargo build -p get_mnv --release
python3 tests/who_catalogue/run.py --who-dir ruta/a/la/carpeta
```

La batería reconstruye entonces la tabla incluida y la exige idéntica byte a
byte, y anota todas las variantes genómicas que la OMS asocia a los grupos 1 y
2: 26.231. 26.225 salen marcadas con su propia entrada, 198 solo por la lectura
del codón, y 2 más con la entrada `LoF` del gen para el mismo fármaco. Entre
las primeras están las deleciones que se llevan un codón de inicio, a menudo
con las bases de antes: get_MNV las da como pérdidas del inicio (`Met1?`), y
la OMS como pérdidas de función. Una de las 2 es una deleción que la OMS
escribe sobre la última base del codón de inicio de gid, gid `4407695`; deja
ATG delante, así que get_MNV lee el corrimiento desde `Ser2`, donde lo pone la
regla 3' de HGVS. Falla si quedan sin marcar para ninguno de sus fármacos más
de estas 4:

- 2 inserciones justo antes de un codón de inicio, tlyA `1917939 C>CAACG` y
  Rv2983 `3339117 C>CG`, que la OMS cuenta como pérdidas de función y get_MNV
  como intergénicas: el codón de inicio en sí queda intacto.
- 2 cambios cuya consecuencia depende de la base de al lado, gid
  `4407589 GC>CA` y tlyA `1917970 CT>TG`. Sobre H37Rv forman un codón de
  parada; sobre el ancestro del MTBC, que difiere de H37Rv justo al lado, dan
  `Ala205Trp` y `Leu11Trp`.

## Licencia de los fixtures

Los fixtures y la tabla incluida contienen información del catálogo de
mutaciones de la OMS, que se distribuye bajo la
[licencia ODC Attribution](https://opendatacommons.org/licenses/by/1-0/).
