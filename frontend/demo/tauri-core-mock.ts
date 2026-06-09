// Standalone-demo mock of `@tauri-apps/api/core` so the BamViewer component can
// render in a plain browser (no Tauri runtime). Vite aliases the real module to
// this file via vite.demo.config.ts. Returns the real get_bam_view fixture that
// was dumped from scenario 02a by the get-mnv-gui fixture generator test.
import sampleBamView from "./sample_bam_view.json";

export function invoke<T = unknown>(cmd: string, _args?: unknown): Promise<T> {
  if (cmd === "get_bam_view") {
    return Promise.resolve(sampleBamView as unknown as T);
  }
  return Promise.reject(new Error(`demo: unmocked Tauri command '${cmd}'`));
}
