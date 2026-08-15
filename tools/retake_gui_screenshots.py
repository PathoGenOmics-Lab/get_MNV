#!/usr/bin/env python3
"""Retake the documentation screenshots, light and dark, at 2880x1840.

Every screenshot on the GUI pages comes from here, so the two themes are the
same window at the same size and the page does not jump when a reader switches.
What is on screen is the real application: the demo harness under `frontend/demo`
answers the Tauri commands with fixtures taken from a genuine run on `example/`,
and the components drawing them are the ones the desktop build ships.

    cd frontend && npm run demo          # serves the app at :5180
    get_mnv --vcf example/... --report /tmp/report.html
    GET_MNV_REPORT=/tmp/report.html python3 tools/retake_gui_screenshots.py

Only `gui-03` and its dark twin are cut from a clip; the rest are whole windows.
The report shot drives the page's own Theme button rather than stamping
`data-theme` on the root, because the variant matrix is drawn into a canvas that
only that handler tells to redraw.
"""

import base64
import json
import subprocess
import os
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

import websocket

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
PORT = 9333
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "docs/assets"
SCRATCH = Path(tempfile.mkdtemp(prefix="get_mnv_shots_"))
# Written by `get_mnv --report` on the example data; see the module docstring.
REPORT = Path(os.environ.get("GET_MNV_REPORT", SCRATCH / "report.html"))


class Chrome:
    def __init__(self, url: str):
        self.proc = subprocess.Popen(
            [CHROME, "--headless=new", f"--remote-debugging-port={PORT}",
             f"--user-data-dir={SCRATCH / 'chrome-profile'}",
             "--no-first-run", "--no-default-browser-check", "--disable-gpu",
             # Chrome refuses a DevTools websocket whose Origin it was not told
             # to expect, and the client library sends one.
             "--remote-allow-origins=*",
             "--hide-scrollbars", "--force-device-scale-factor=2", url],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        ws_url = None
        for _ in range(60):
            try:
                pages = json.loads(urllib.request.urlopen(
                    f"http://127.0.0.1:{PORT}/json").read())
                for p in pages:
                    if p.get("type") == "page" and p.get("webSocketDebuggerUrl"):
                        ws_url = p["webSocketDebuggerUrl"]
                        break
                if ws_url:
                    break
            except Exception:
                pass
            time.sleep(0.25)
        if not ws_url:
            sys.exit("chrome no levanto su puerto de depuracion")
        self.ws = websocket.create_connection(ws_url, timeout=30)
        self.n = 0
        self.send("Page.enable")
        self.send("Runtime.enable")

    def send(self, method, **params):
        self.n += 1
        self.ws.send(json.dumps({"id": self.n, "method": method, "params": params}))
        while True:
            msg = json.loads(self.ws.recv())
            if msg.get("id") == self.n:
                if "error" in msg:
                    raise RuntimeError(f"{method}: {msg['error']}")
                return msg.get("result", {})

    def js(self, expression, await_promise=True):
        r = self.send("Runtime.evaluate", expression=f"(async () => {{ {expression} }})()",
                      awaitPromise=await_promise, returnByValue=True)
        if r.get("exceptionDetails"):
            raise RuntimeError(r["exceptionDetails"].get("text", "js fallo"))
        return r.get("result", {}).get("value")

    def viewport(self, w, h):
        self.send("Emulation.setDeviceMetricsOverride",
                  width=w, height=h, deviceScaleFactor=2, mobile=False)

    def shot(self, name, clip=None):
        params = {"format": "png", "captureBeyondViewport": bool(clip)}
        if clip:
            # The device scale factor already doubles the capture; asking the
            # clip for 2 as well returns it at 4x.
            params["clip"] = {**clip, "scale": 1}
        data = self.send("Page.captureScreenshot", **params)["data"]
        path = OUT / name
        path.write_bytes(base64.b64decode(data))
        from PIL import Image
        print(f"  {name:28} {Image.open(path).size[0]}x{Image.open(path).size[1]}")

    def goto(self, url):
        self.send("Page.navigate", url=url)
        time.sleep(2.0)

    def close(self):
        try:
            self.ws.close()
        finally:
            self.proc.terminate()


ESPERAR_APP = """
  for (let i = 0; i < 80; i++) {
    if (document.querySelectorAll('.drop-zone').length >= 4) return true;
    await new Promise(r => setTimeout(r, 100));
  }
  throw new Error('la app no monto');
"""

PREPARAR = """
  const pausa = ms => new Promise(r => setTimeout(r, ms));
  if (document.documentElement.getAttribute('data-theme') !== 'dark') {
    document.querySelector('.theme-toggle').click();
    await pausa(250);
  }
  for (const zona of document.querySelectorAll('.drop-zone')) {
    if (!zona.querySelector('.drop-zone-clear')) { zona.click(); await pausa(300); }
  }
  await pausa(500);
  return document.querySelectorAll('.drop-zone-clear').length;
"""

# The caption for this one says the read-support floors are at zero, so they are
# set here rather than left at the preset's 2.
A_CERO = """
  const fijar = (el, v) => {
    const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
    setter.call(el, String(v));
    el.dispatchEvent(new Event('input', {bubbles: true}));
    el.dispatchEvent(new Event('change', {bubbles: true}));
  };
  const cabecera = [...document.querySelectorAll('.param-group-header')]
    .find(b => /read support/i.test(b.textContent));
  const campos = [...cabecera.parentElement.querySelectorAll('input[type=range]')];
  fijar(campos[0], 0);
  await new Promise(r => setTimeout(r, 150));
  fijar(campos[1], 0);
  await new Promise(r => setTimeout(r, 400));
  return [...cabecera.parentElement.querySelectorAll('input[type=text]')].map(i => i.value).join(',');
"""


def main():
    c = Chrome("http://localhost:5180/app.html")
    try:
        c.viewport(1440, 920)
        c.js(ESPERAR_APP)
        print("app lista:", c.js(PREPARAR), "entradas rellenas")
        # The light twin of this shot already has the floors at zero, so they go
        # down before it is taken and the pair differs only in the theme.
        print("umbrales:", c.js(A_CERO))
        c.shot("gui-02-inputs-dark.png")

        c.viewport(1440, 1240)
        time.sleep(0.6)
        caja = c.js("""
          const s = document.querySelector('.analysis-sidebar');
          s.scrollIntoView({block: 'start'});
          await new Promise(r => setTimeout(r, 300));
          const r = s.getBoundingClientRect();
          return {x: r.left, y: r.top + scrollY, width: r.width, height: r.height};
        """)
        # A few pixels of margin, or the panel's own left edge is shaved off.
        caja["x"] -= 12
        caja["width"] += 24
        recorte = {k: round(v) for k, v in caja.items()}
        c.shot("gui-03-parameters-dark.png", clip=recorte)
        # The light twin is taken from the same box in the same run. Any other
        # way and the two differ in height, which makes the page jump when the
        # reader switches theme.
        c.js("""
          document.querySelector('.theme-toggle').click();
          await new Promise(r => setTimeout(r, 500));
          return document.documentElement.getAttribute('data-theme');
        """)
        c.shot("gui-03-parameters.png", clip=recorte)
        c.js("""
          document.querySelector('.theme-toggle').click();
          await new Promise(r => setTimeout(r, 500));
          return document.documentElement.getAttribute('data-theme');
        """)
        c.viewport(1440, 920)
        time.sleep(0.4)

        c.js("document.querySelector('.run-button').click(); return true;")
        time.sleep(0.65)
        c.shot("gui-04-running-dark.png")

        time.sleep(2.2)
        c.js("""
          const r = [...document.querySelectorAll('.tab')].find(t => /results/i.test(t.textContent));
          if (!r.classList.contains('active')) r.click();
          await new Promise(x => setTimeout(x, 500));
          window.scrollTo(0, 0);
          return document.querySelector('.tab.active').textContent.trim();
        """)
        time.sleep(0.4)
        c.shot("gui-05-results-dark.png")

        c.js("""
          const s = [...document.querySelectorAll('section.step-section')]
            .find(x => /Genomic Track Viewer/.test(x.textContent));
          s.scrollIntoView({block: 'start'});
          await new Promise(r => setTimeout(r, 400));
          window.scrollBy(0, 70);
          await new Promise(r => setTimeout(r, 400));
          return Math.round(scrollY);
        """)
        c.shot("gui-06-reads-dark.png")

        c.js("""
          const s = [...document.querySelectorAll('section.step-section')]
            .find(x => /Variant Data/.test(x.textContent));
          s.scrollIntoView({block: 'start'});
          await new Promise(r => setTimeout(r, 500));
          return Math.round(scrollY);
        """)
        c.shot("gui-07-table-dark.png")

        # The report is a separate page. It stamps its own theme on the root, so
        # the dark one is asked for directly rather than by clicking through.
        c.goto(f"file://{REPORT}")
        c.viewport(1440, 920)
        # Driven through the report's own Theme button rather than by stamping
        # the attribute: the variant matrix is drawn into a canvas, and only the
        # button's handler tells it to redraw. Setting the attribute by hand
        # leaves the matrix painted for the theme it was drawn in.
        print("tema del informe:", c.js("""
          const boton = document.getElementById('themeBtn');
          for (let i = 0; i < 4; i++) {
            if (document.documentElement.getAttribute('data-theme') === 'dark') break;
            boton.click();
            await new Promise(r => setTimeout(r, 500));
          }
          await new Promise(r => setTimeout(r, 900));
          window.scrollTo(0, 0);
          return document.documentElement.getAttribute('data-theme');
        """))
        time.sleep(0.5)
        c.shot("cli-01-report-dark.png")
    finally:
        c.close()


if __name__ == "__main__":
    main()
