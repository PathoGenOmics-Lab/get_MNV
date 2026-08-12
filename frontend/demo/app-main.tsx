// Standalone demo entry for the whole desktop app: renders the real `App` with
// the Tauri runtime replaced by the fixtures in `tauri-mocks.ts`. Everything on
// screen is the actual component the desktop build ships; only the backend
// answers are canned.
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import App from "../src/App";
import "../src/index.css";
import "../src/App.css";

const root = document.getElementById("root");
if (!root) throw new Error("missing #root");

createRoot(root).render(
  <StrictMode>
    <App />
  </StrictMode>,
);
