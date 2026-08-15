import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import path from "node:path";

// Standalone demo build/serve for the BamViewer component. Aliases the Tauri
// core API to a mock so the component renders in a plain browser with a real
// get_bam_view fixture. Run: npx vite --config vite.demo.config.ts
export default defineConfig({
  root: path.resolve(__dirname, "demo"),
  plugins: [react()],
  resolve: {
    alias: {
      // The BamViewer demo only needs `invoke`; the whole-app demo (app.html)
      // needs every Tauri module the app imports, so all of them point at the
      // shared mock. Nothing here reimplements analysis: see demo/tauri-mocks.ts.
      "@tauri-apps/api/core": path.resolve(__dirname, "demo/tauri-mocks.ts"),
      "@tauri-apps/api/event": path.resolve(__dirname, "demo/tauri-mocks.ts"),
      "@tauri-apps/api/webview": path.resolve(__dirname, "demo/tauri-mocks.ts"),
      "@tauri-apps/plugin-dialog": path.resolve(__dirname, "demo/tauri-mocks.ts"),
      "@tauri-apps/plugin-shell": path.resolve(__dirname, "demo/tauri-shell-mock.ts"),
    },
  },
  server: { port: 5180, strictPort: true },
});
