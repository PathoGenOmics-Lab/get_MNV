# Benchmarking

get_MNV incluye un binario de benchmark integrado para medir el rendimiento de la anotación de variantes.

## Benchmark sintético

```bash
cargo run --release --bin bench_variants -- --warmup 5 --iters 30
```

Aumenta la escala del conjunto de datos sintético:

```bash
cargo run --release --bin bench_variants -- \
  --warmup 3 --iters 20 --threads 4 --synthetic-scale 4
```

## Benchmark con conjunto de datos

Ejecútalo sobre un conjunto de datos real y exporta los resultados a CSV:

```bash
cargo run --release --bin bench_variants -- \
  --dataset example \
  --warmup 5 --iters 20 \
  --threads 4 \
  --csv benchmark.csv \
  --max-avg-ms 200
```

Esto añade una fila por cada ejecución a `benchmark.csv`.

## Opciones del benchmark

| Argumento | Por defecto | Descripción |
|----------|---------|-------------|
| `--warmup <N>` | 5 | Iteraciones de calentamiento (no se miden) |
| `--iters <N>` | 30 | Iteraciones medidas |
| `--threads <N>` | 1 | Hilos de trabajo |
| `--dataset <DIR>` | none | Ruta al directorio del conjunto de datos (FASTA + VCF + anotación) |
| `--contig <NAME>` | all | Limita la ejecución a un solo contig |
| `--csv <FILE>` | none | Añade los resultados a un CSV |
| `--max-avg-ms <MS>` | none | Falla si la media supera el umbral (detección de regresiones) |
| `--synthetic-scale <N>` | 1 | Factor de escala para el número sintético de genes/SNP |

## Scripts reproducibles

Consulta `scripts/` para los ejecutores de benchmark por lotes:
- `scripts/run_reproducible_benchmark.sh`
- `scripts/reproduce_example_run.sh`
