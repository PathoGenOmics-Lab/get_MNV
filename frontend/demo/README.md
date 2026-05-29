# BamViewer demo harness

**English** · [Español](README.es.md)

Renders the real `BamViewer` component in a plain browser, without launching the
full Tauri desktop app. The Tauri `get_bam_view` command is mocked
(`tauri-core-mock.ts`) to return a fixture captured from the core engine, so the
read pileup, codon/AA tracks, coverage and reference rows render exactly as in
the app.

```bash
cd frontend
npm install        # first time only
npm run demo       # serves the demo at http://localhost:5180
```

Append `?reads=1` to the URL to auto-scroll the read pileup into view.

## Fixtures

- `sample_bam_view.json` — a real `get_bam_view` response for scenario `02a`
  (`chr_test` MNV `28 G>T` / `30 T>A`, codon `GCT→TCA` = Ala10Ser, 20 supporting
  reads), produced by calling the GUI backend command on the scenario BAM/FASTA.
- `sample_tsv.json` — the matching `TsvData` parsed from that scenario's
  `variants.MNV.tsv`.

To refresh the fixtures, call `get_bam_view` (src-tauri) on a scenario work dir
and serialize the response, then convert the scenario `.MNV.tsv` to
`{headers, rows}`.
