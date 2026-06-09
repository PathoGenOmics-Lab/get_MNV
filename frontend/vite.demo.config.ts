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
      "@tauri-apps/api/core": path.resolve(__dirname, "demo/tauri-core-mock.ts"),
    },
  },
  server: { port: 5180, strictPort: true },
});
