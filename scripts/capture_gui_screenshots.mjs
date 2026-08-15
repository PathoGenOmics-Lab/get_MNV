// Regenerates the screenshots in docs/gui-tutorial.md by driving the browser
// demo of the desktop app. Every component on screen is the one the desktop app
// ships; only the Tauri calls are mocked (see frontend/demo/tauri-mocks.ts), so
// the pages can be captured without building a desktop bundle.
//
// Usage, from the repository root:
//
//   npm --prefix frontend run demo &        # serves http://localhost:5180
//   npm install --no-save puppeteer-core    # not a project dependency
//   node scripts/capture_gui_screenshots.mjs
//
// CHROME can point at any Chromium build; the default is Chrome on macOS.
import puppeteer from "puppeteer-core";
import path from "node:path";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const CHROME =
  process.env.CHROME ||
  "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const URL = process.env.DEMO_URL || "http://localhost:5180/app.html";
const OUT = path.join(HERE, "..", "docs", "assets");

const wait = (ms) => new Promise((r) => setTimeout(r, ms));

async function shot(page, name, clip) {
  await page.screenshot({ path: path.join(OUT, name), ...(clip ? { clip } : {}) });
  console.log("  wrote", name);
}

const browser = await puppeteer.launch({
  executablePath: CHROME,
  headless: "new",
  defaultViewport: { width: 1440, height: 920, deviceScaleFactor: 2 },
  args: ["--hide-scrollbars"],
});
const page = await browser.newPage();
await page.goto(URL, { waitUntil: "networkidle0" });
await wait(600);

/** Types a value into the numeric box of a labelled parameter row. */
async function setParam(label, value) {
  const outcome = await page.evaluate(
    (label, value) => {
      const lbl = Array.from(document.querySelectorAll("*")).find(
        (e) => e.children.length === 0 && e.textContent.trim() === label,
      );
      if (!lbl) return "label not found";
      let row = lbl.parentElement;
      for (let i = 0; i < 4 && row; i++) {
        const input = row.querySelector('input[type="text"]');
        if (input) {
          const setter = Object.getOwnPropertyDescriptor(
            window.HTMLInputElement.prototype,
            "value",
          ).set;
          setter.call(input, String(value));
          input.dispatchEvent(new Event("input", { bubbles: true }));
          input.dispatchEvent(new Event("change", { bubbles: true }));
          return "set";
        }
        row = row.parentElement;
      }
      return "input not found";
    },
    label,
    value,
  );
  console.log(`  ${label} -> ${value}: ${outcome}`);
  if (outcome !== "set") throw new Error(`could not set ${label}`);
}

// The form ships with a 2-read floor on SNP and MNV support, and the bundled
// demo BAM covers a single locus, so that floor would leave a one-row table.
// Lower both so the form matches the fixture the demo answers with. The
// tutorial says so where it shows this screen.
await setParam("Min SNP reads", 0);
await setParam("Min MNV reads", 0);
await wait(300);

// Fill the inputs. Each click opens the mocked picker, which answers with the
// example-data path matching the filter it was given.
const zones = await page.$$("div.drop-zone");
if (zones.length !== 4) throw new Error(`expected 4 drop zones, saw ${zones.length}`);
for (const zone of zones) {
  await zone.click();
  await wait(250);
}
await wait(400);
await shot(page, "gui-02-inputs.png");

// The parameter sidebar on its own, so the knobs are legible.
const aside = await page.$("aside");
const box = aside && (await aside.boundingBox());
if (box) {
  await shot(page, "gui-03-parameters.png", {
    x: Math.floor(box.x),
    y: Math.floor(box.y),
    width: Math.ceil(box.width),
    height: Math.min(Math.ceil(box.height), 900),
  });
}

// Catch the run in flight; the mock walks three phases 450 ms apart.
await (await page.$("button.run-button")).click();
await wait(500);
await shot(page, "gui-04-running.png");

// Results.
await wait(1600);
for (const tab of await page.$$("nav button.tab")) {
  if ((await page.evaluate((e) => e.textContent.trim(), tab)) === "Results") {
    await tab.click();
    break;
  }
}
await wait(900);
await shot(page, "gui-05-results.png");

// The read pileup, at the one codon the example BAM covers.
const search = await page.$('input[placeholder*="Search gene"]');
if (search) {
  await search.click();
  await search.type("Rv2036");
  await wait(600);
}
const opened = await page.evaluate(() => {
  const hit = Array.from(document.querySelectorAll("*")).find(
    (e) => e.children.length === 0 && /^Rv2036/.test(e.textContent.trim()),
  );
  if (!hit) return false;
  (hit.closest("button, li, [role=button], div[class*=locus]") || hit).click();
  return true;
});
if (!opened) throw new Error("the Rv2036 locus was not in the list");
await wait(1400);

const scrollTo = (text, extra) =>
  page.evaluate(
    (t, dy) => {
      const heading = Array.from(document.querySelectorAll("*")).find(
        (e) => e.children.length === 0 && e.textContent.trim() === t,
      );
      heading?.scrollIntoView({ block: "start" });
      window.scrollBy(0, dy);
    },
    text,
    extra,
  );

await scrollTo("Genomic Track Viewer", 90);
await wait(700);
await shot(page, "gui-06-reads.png");

await scrollTo("Variant Data", -40);
await wait(700);
await shot(page, "gui-07-table.png");

await browser.close();
console.log("done");
