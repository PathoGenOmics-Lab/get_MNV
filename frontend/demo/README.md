# Browser demo harness

**English** · [Español](README.es.md)

Runs the desktop app's real components in a plain browser, without building the
Tauri bundle. Every `@tauri-apps/*` module the app imports is aliased in
`vite.demo.config.ts` to `tauri-mocks.ts`, which answers each command with
fixture data. Nothing here reimplements analysis: the fixtures are genuine
output of the engine, and the components rendering them are the ones the app
ships.

```bash
cd frontend
npm install        # first time only
npm run demo       # serves at http://localhost:5180
```

Two pages:

| URL | What it shows |
|---|---|
| `/app.html` | the whole application: inputs, parameters, run, results |
| `/` | the `BamViewer` component on its own |

Append `?reads=1` to the component page to auto-scroll the read pileup into view.

## Fixtures

All three describe one dataset, the files in `example/`, so the summary counts,
the table and the pileup agree with each other.

- `sample_tsv.json`: the 941-row `TsvData` from running `get_mnv` on
  `example/` with a BAM.
- `sample_bam_view.json`: a real `get_bam_view` response at the `Rv2036` codon
  (`MTB_anc:2282376,2282377`, `GTT → GCC` = Val93Ala, 24 supporting reads). It
  is the only locus the small `G35894.demo.bam` covers, so pick `Rv2036` in the
  locus list to see reads that belong to the row.
- The run summary is inlined in `tauri-mocks.ts`, copied from `--summary-json`.

To refresh them:

```bash
# the pileup
cargo test -p get-mnv-gui --bins -- --ignored regenerate_demo_bam_view

# the table and the summary
get_mnv --vcf example/G35894.var.snp.vcf --fasta example/MTB_ancestor.fas \
        --genes example/anot_genes.txt  --bam example/G35894.demo.bam \
        --summary-json summary.json
# then convert the .MNV.tsv to {headers, rows} and paste summary.json's
# figures into RUN_SUMMARY.
```

## Screenshots

The screenshots in `docs/gui-tutorial.md` are captured from `/app.html`:

```bash
npm --prefix frontend run demo &
npm install --no-save puppeteer-core
node scripts/capture_gui_screenshots.mjs
```
