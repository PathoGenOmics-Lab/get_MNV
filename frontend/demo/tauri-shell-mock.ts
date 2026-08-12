// `@tauri-apps/plugin-shell` exposes `open`; the shared mock names it `openPath`
// so it does not collide with the dialog plugin's own `open`.
export { openPath as open } from "./tauri-mocks";
