# Banco de pruebas en navegador

[English](README.md) · **Español**

Ejecuta los componentes reales de la aplicación de escritorio en un navegador
normal, sin compilar el bundle de Tauri. Cada módulo `@tauri-apps/*` que importa
la app se redirige en `vite.demo.config.ts` a `tauri-mocks.ts`, que responde a
cada comando con datos de fixture. Aquí no se reimplementa ningún análisis: los
fixtures son salida auténtica del motor, y los componentes que los dibujan son
los mismos que lleva la app.

```bash
cd frontend
npm install        # solo la primera vez
npm run demo       # sirve en http://localhost:5180
```

Dos páginas:

| URL | Qué muestra |
|---|---|
| `/app.html` | la aplicación entera: entradas, parámetros, ejecución, resultados |
| `/` | el componente `BamViewer` por separado |

Añade `?reads=1` a la página del componente para que el pileup se desplace solo
hasta quedar a la vista.

## Fixtures

Los tres describen un mismo conjunto de datos, los archivos de `example/`, así
que los recuentos del resumen, la tabla y el pileup concuerdan entre sí.

- `sample_tsv.json`: el `TsvData` de 941 filas que sale de ejecutar `get_mnv`
  sobre `example/` con un BAM.
- `sample_bam_view.json`: una respuesta real de `get_bam_view` en el codón de
  `Rv2036` (`MTB_anc:2282376,2282377`, `GTT → GCC` = Val93Ala, 24 lecturas de
  soporte). Es el único locus que cubre el pequeño `G35894.demo.bam`, así que
  elige `Rv2036` en la lista de loci para ver lecturas que correspondan a la
  fila.
- El resumen de la ejecución está incrustado en `tauri-mocks.ts`, copiado de
  `--summary-json`.

Para regenerarlos:

```bash
# el pileup
cargo test -p get-mnv-gui --bins -- --ignored regenerate_demo_bam_view

# la tabla y el resumen
get_mnv --vcf example/G35894.var.snp.vcf --fasta example/MTB_ancestor.fas \
        --genes example/anot_genes.txt  --bam example/G35894.demo.bam \
        --summary-json summary.json
# luego convierte el .MNV.tsv a {headers, rows} y pega las cifras de
# summary.json en RUN_SUMMARY.
```

## Capturas

Las capturas de `docs/gui-tutorial.md` se toman de `/app.html`:

```bash
npm --prefix frontend run demo &
npm install --no-save puppeteer-core
node scripts/capture_gui_screenshots.mjs
```
