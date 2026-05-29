# Banco de pruebas de demostración de BamViewer

[English](README.md) · **Español**

Renderiza el componente real `BamViewer` en un navegador corriente, sin lanzar la
app de escritorio completa de Tauri. El comando `get_bam_view` de Tauri está mockeado
(`tauri-core-mock.ts`) para devolver un fixture capturado del motor central, de modo que el
pileup de lecturas, las pistas de codón/AA, la cobertura y las filas de referencia se renderizan exactamente como en
la app.

```bash
cd frontend
npm install        # first time only
npm run demo       # serves the demo at http://localhost:5180
```

Añade `?reads=1` a la URL para desplazar automáticamente el pileup de lecturas hasta la vista.

## Fixtures

- `sample_bam_view.json` — una respuesta real de `get_bam_view` para el escenario `02a`
  (`chr_test` MNV `28 G>T` / `30 T>A`, codón `GCT→TCA` = Ala10Ser, 20 lecturas de
  apoyo), producida llamando al comando del backend de la GUI sobre el BAM/FASTA del escenario.
- `sample_tsv.json` — el `TsvData` correspondiente, parseado a partir del
  `variants.MNV.tsv` de ese escenario.

Para actualizar los fixtures, llama a `get_bam_view` (src-tauri) sobre un directorio de trabajo de un escenario
y serializa la respuesta; luego convierte el `.MNV.tsv` del escenario a
`{headers, rows}`.
