# Banco de pruebas de demostración de BamViewer

[English](README.md) · **Español**

Renderiza el componente real `BamViewer` en un navegador normal, sin lanzar la
app de escritorio completa de Tauri. El comando `get_bam_view` de Tauri se
sustituye por un mock (`tauri-core-mock.ts`) que devuelve un fixture capturado
del motor principal, de modo que el pileup de lecturas, las pistas de codón/AA,
la cobertura y las filas de referencia se renderizan igual que en la app.

```bash
cd frontend
npm install        # first time only
npm run demo       # serves the demo at http://localhost:5180
```

Añade `?reads=1` a la URL para desplazar el pileup de lecturas hasta que quede a la vista.

## Fixtures

- `sample_bam_view.json` — una respuesta real de `get_bam_view` para el escenario `02a`
  (`chr_test` MNV `28 G>T` / `30 T>A`, codón `GCT→TCA` = Ala10Ser, 20 lecturas de
  soporte), generada al ejecutar el comando del backend de la GUI sobre el BAM/FASTA del escenario.
- `sample_tsv.json` — el `TsvData` correspondiente, obtenido del
  `variants.MNV.tsv` de ese escenario.

Para actualizar los fixtures, ejecuta `get_bam_view` (src-tauri) sobre el directorio de trabajo
de un escenario y serializa la respuesta; luego convierte el `.MNV.tsv` del escenario a
`{headers, rows}`.
