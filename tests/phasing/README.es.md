# Validación de fasing contra una mezcla conocida

[English](README.md) · **Español**

La [suite de verdad sintética](../truth/README.es.md) valida la aritmética de
codones, pero no tiene BAM, así que no dice nada de la otra mitad de una llamada
MNV: si las sustituciones están en las *mismas moléculas*. Dos variantes en un
codón repartidas en dos moléculas distintas no son un MNV a nivel de codón, y
ninguna cantidad de aritmética de codones las distingue.

Esta suite cubre esa mitad. Construye alineamientos molécula a molécula,
eligiendo qué lecturas llevan qué variantes, y pregunta si get_MNV recupera la
mezcla que se le dio.

```bash
cargo build -p get_mnv --release
python3 tests/phasing/run.py
```

## Cómo se deriva la verdad

Aquí no hay ni un número copiado de una ejecución anterior. Un caso es una lista
de clases de molécula: cuántas llevan qué variantes, y si una lectura sacada de
esa molécula alcanza el codón entero. La expectativa sale solo de esos recuentos:

```
denominador = mínimo de moléculas que cruzan el codón y llevan un constituyente
numerador   = moléculas que cruzan el codón y los llevan todos
```

Si el código y la aritmética discrepan, uno de los dos está mal.

## Qué cubre

Cinco geometrías de codón, cada una barrida por las 21 mezclas que van de todo
en trans a todo en cis:

- un codón en hebra más
- un codón en hebra menos, que debe dar respuestas idénticas, porque el
  ligamiento es una propiedad de las moléculas y no de la hebra por la que se lee
  el gen
- un codón partido por un intrón, cruzado por lecturas con una `N` en el CIGAR
- las tres posiciones de un mismo codón, que ejercita la regla del constituyente
  menos soportado con más de dos
- el codón de hebra más otra vez, secuenciado como pares de mates donde cada uno
  alcanza un extremo del codón y no el otro, así que la respuesta solo existe a
  nivel de la molécula

Más los casos en los que la respuesta honesta no es un número:

- lecturas que llevan una variante y terminan antes de la otra, que no deben ni
  apoyar ni diluir el ratio
- ninguna lectura que cruce el codón, que es desconocido (`-`), no cero
- un haplotipo minoritario anidado dentro de una variante mayoritaria, que es
  soporte total sobre pocas lecturas, y la razón de que el recuento se publique
  al lado del ratio

120 mezclas, unos segundos.

## Comprobar que todavía puede fallar

Un test que no puede fallar no vale nada. Contra el binario anterior al arreglo
de fasing reporta 171 discrepancias:

```bash
GET_MNV=/ruta/a/get_mnv/anterior python3 tests/phasing/run.py
```

Merece la pena leer esos fallos, porque son las tres formas en que el ratio
antiguo se equivocaba a la vez. Dividía por las lecturas que llevaban un SNV
*solo*, excluyendo las que llevaban el haplotipo entero, así que sobreestimaba el
ligamiento en todas partes y se clavaba en `1.0000` a partir de la mezcla mitad y
mitad: un cara o cruz entre cis y trans se reportaba como haplotipo perfecto.
Cuando todas las lecturas llevaban el haplotipo completo no quedaban lecturas
sueltas por las que dividir, y un ligamiento impecable salía como `-`. Y cuando
ninguna lectura cruzaba el codón devolvía `0.0000`, que se lee como prueba de
trans, para una pregunta que nadie había respondido.

## Cuánto vale

Como la suite de verdad, es una segunda implementación de la misma mano, así que
caza errores de implementación y no un malentendido compartido. Lo que la hace
más fuerte que una expectativa escrita a mano es que la mezcla es la entrada: la
respuesta queda fijada antes de que get_MNV se ejecute, y un barrido no deja
sitio donde esconderse a un off-by-one ni a un ratio que satura.

Lo que no cubre: artefactos reales de alineamiento. Aquí toda lectura es un
registro limpio con MAPQ 60 y bases Q40, colocado según una geometría elegida.
El soft clipping, los pares solapados y el sesgo de referencia en los extremos de los amplicones
quedan fuera de su alcance.
