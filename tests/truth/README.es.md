# Validación con verdad sintética

[English](README.md) · **Español**

La [suite diferencial](../differential/README.es.md) compara get_MNV con
`bcftools csq` y SnpEff. Esa es la comprobación más fuerte, porque esas
herramientas no son nuestras. Pero no puede validar justo aquello para lo que
existe get_MNV: **los MNV a nivel de codón**. Ningún anotador por variante
combina dos SNV de un codón como lo hace get_MNV, así que precisamente en esas
llamadas los oráculos discrepan por diseño y no pueden decir si el residuo es
correcto.

Esta suite cubre ese hueco, con un conjunto de datos cuya respuesta conocemos
porque lo construimos nosotros.

```bash
cargo build -p get_mnv --release
python3 tests/truth/run.py
```

## Cómo se deriva la verdad

La expectativa no puede salir de la lógica de get_MNV, o el test solo le da la
razón al código. Así que toma otra ruta hacia la misma respuesta:

```
mutar el genoma  ->  volver a extraer el CDS  ->  traducir  ->  comparar proteínas
```

get_MNV hace cirugía a nivel de codón, in situ. Esto vuelve a traducir el CDS
entero desde una tabla de codones escrita en `run.py`, con operaciones de cadena
sobre la referencia. Un desliz de implementación aparece como discrepancia.

## Qué cubre

Cuatro genes, generados de forma determinista: hebra más y hebra menos, de un
solo exón y con splicing. Para cada codón de cada uno:

- cada sustitución simple en cada posición (3 posiciones x 3 alternativos)
- cada par de sustituciones en el mismo codón
- el triple

Unos 1560 casos, cada uno como su propio VCF para que no se interfieran. De cada
uno se comprueban gen, cambio de aminoácido, tipo de cambio, codón de
referencia, codón alternativo y el descriptor `HGVS c.`. Tarda unos cinco
segundos.

## Cuánto vale

Menos que un oráculo de terceros, y conviene decir por qué: es una segunda
implementación de la misma mano. Caza errores de implementación, no un
malentendido compartido sobre la biología. Las *etiquetas* de consecuencia
(`Synonymous`, `Stop gained`, ...) siguen la convención documentada de get_MNV,
así que en eso valida el código contra la documentación y no contra la
naturaleza. El aminoácido, los codones y la coordenada del CDS son la parte
independiente.

Y es la razón de que la suite diferencial siga existiendo a su lado. Ninguna
sustituye a la otra.

## Comprobar que todavía puede fallar

Un test que no puede fallar no vale nada. `GET_MNV_BIN` apunta la suite a otro
binario:

```bash
GET_MNV_BIN=/ruta/a/get_mnv/anterior python3 tests/truth/run.py
```

Contra el binario anterior al arreglo de orientación de codón reporta 240
discrepancias, todas en los genes de hebra menos.

## Una nota sobre escribir el oráculo

Su primera ejecución reportó 665 discrepancias. Todas eran un bug de este
fichero, no de get_MNV: una variable de bucle filtrada hacía que el VCF generado
llevara un alelo alternativo distinto de aquel para el que se había calculado la
expectativa. El oráculo hay que comprobarlo tanto como a lo que comprueba.
