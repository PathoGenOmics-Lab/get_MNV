// Regenerates docs/assets/cli-01-report.png, the figure in the command-line
// tutorial, from a report get_mnv actually wrote.
//
// Usage, from the repository root:
//
//   cargo build -p get_mnv --release
//   mkdir -p /tmp/get_mnv_report && cd /tmp/get_mnv_report
//   <repo>/target/release/get_mnv \
//     --vcf <repo>/example/G35894.var.snp.vcf \
//     --fasta <repo>/example/MTB_ancestor.fas \
//     --genes <repo>/example/anot_genes.txt \
//     --bam <repo>/example/G35894.demo.bam \
//     --report sample.html
//   cd <repo>
//   npm install --no-save puppeteer-core
//   REPORT=/tmp/get_mnv_report/sample.html node scripts/capture_report_screenshot.mjs
import puppeteer from "puppeteer-core";
import path from "node:path";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const CHROME =
  process.env.CHROME ||
  "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const REPORT = process.env.REPORT;
if (!REPORT) throw new Error("set REPORT to the .html file get_mnv wrote");

const browser = await puppeteer.launch({
  executablePath: CHROME,
  headless: "new",
  defaultViewport: { width: 1440, height: 920, deviceScaleFactor: 2 },
  args: ["--hide-scrollbars", "--allow-file-access-from-files"],
});
const page = await browser.newPage();
await page.goto(`file://${path.resolve(REPORT)}`, { waitUntil: "networkidle0" });
await new Promise((r) => setTimeout(r, 1500));
const out = path.join(HERE, "..", "docs", "assets", "cli-01-report.png");
await page.screenshot({ path: out });
console.log("wrote", out);
await browser.close();
