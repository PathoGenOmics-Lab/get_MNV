// Behavioural checks for the JavaScript and TypeScript the project ships.
//
// The report's table, filters and matrix are built in the page, and the desktop
// form merges a preset in the browser. Nothing in the Rust suites can run either,
// so both were guarded by tests that read the source and checked a shape was
// still there. Those catch a deletion and nothing else: they cannot tell whether
// the code does what it says. This runs the real functions over real data.
//
// No test framework and no dependencies: node, the shipped sources, and a report
// this script generates from the bundled example. Exits non-zero on any failure,
// like the scenario runner beside it.

import { execFileSync } from "node:child_process";
import { mkdtempSync, readFileSync, writeFileSync, rmSync, readdirSync, statSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { createContext, runInContext } from "node:vm";

const REPO = join(dirname(fileURLToPath(import.meta.url)), "..", "..");

let failures = 0;
let checks = 0;

function check(name, got, want) {
  checks++;
  const ok = JSON.stringify(got) === JSON.stringify(want);
  if (!ok) failures++;
  const mark = ok ? "PASS" : "FAIL";
  console.log(`${mark}  ${name}`);
  if (!ok) console.log(`      got ${JSON.stringify(got)}\n      want ${JSON.stringify(want)}`);
}

/// The newest source file under a directory, which is what a binary has to be
/// newer than to be the binary these checks are about.
function newestSource(dir) {
  let newest = 0;
  for (const entry of readdirSync(dir, { withFileTypes: true })) {
    const path = join(dir, entry.name);
    newest = Math.max(newest, entry.isDirectory() ? newestSource(path) : statSync(path).mtimeMs);
  }
  return newest;
}

/// The binary these checks run against. Both profiles may exist, and picking
/// the wrong one means checking a build that no longer matches the source: the
/// report page is embedded at compile time, so a stale binary ships stale
/// JavaScript and every check here passes against code nobody can run.
function binary() {
  const built = [];
  for (const profile of ["debug", "release"]) {
    const path = join(REPO, "target", profile, "get_mnv");
    try {
      execFileSync(path, ["--version"], { stdio: "ignore" });
      built.push({ path, mtime: statSync(path).mtimeMs });
    } catch {
      /* not built under this profile */
    }
  }
  if (built.length === 0) {
    throw new Error("no get_mnv binary in target/debug or target/release; build it first");
  }
  built.sort((a, b) => b.mtime - a.mtime);
  const source = newestSource(join(REPO, "src"));
  if (built[0].mtime < source) {
    throw new Error(
      `${built[0].path} is older than src/; rebuild it, or these checks would ` +
        "be reading a report page that the source no longer produces");
  }
  return built[0].path;
}

// ---------------------------------------------------------------------------
// A cohort report with two overlapping genes, which is the shape the matrix and
// the filters have to get right.
// ---------------------------------------------------------------------------

function buildReport(work) {
  const bases = [];
  for (let i = 0; i < 900; i++) bases.push("ACGT"[(i * 7 + 3) % 4]);
  const sequence = bases.join("");
  writeFileSync(join(work, "ref.fasta"), `>chr1\n${sequence}\n`);
  // geneP and geneQ cover the same span in different frames, so a base belongs
  // to both; geneR sits apart.
  writeFileSync(
    join(work, "genes.gff3"),
    "##gff-version 3\n" +
      "chr1\tsyn\tCDS\t101\t400\t.\t+\t0\tID=c1;Name=geneP\n" +
      "chr1\tsyn\tCDS\t102\t401\t.\t+\t0\tID=c2;Name=geneQ\n" +
      "chr1\tsyn\tCDS\t501\t700\t.\t-\t0\tID=c3;Name=geneR\n",
  );
  const other = (b) => ({ A: "C", C: "G", G: "T", T: "A" })[b];
  // 372 and 373 share a codon in both genes, so they come back as one MNV; the
  // deletion at 620 gives the indel and HIGH lanes something to count.
  const sites = [150, 200, 260, 333, 372, 373, 520, 600, 620, 800];
  const records = sites.map((p, n) => {
    const gts = [0, 1, 2].map((k) => ["1/1", "0/1", "1/1"][(n + k) % 3]);
    const [ref, alt] = p === 620
      ? [bases[p - 1] + bases[p] + bases[p + 1], bases[p - 1]]
      : [bases[p - 1], other(bases[p - 1])];
    return `chr1\t${p}\t.\t${ref}\t${alt}\t100\tPASS\tDP=30\tGT:DP:AF\t` +
      gts.map((g) => `${g}:30:0.8`).join("\t");
  });
  writeFileSync(
    join(work, "cohort.vcf"),
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=900>\n" +
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n' +
      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">\n' +
      '##FORMAT=<ID=AF,Number=A,Type=Float,Description="f">\n' +
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n" +
      records.join("\n") + "\n",
  );
  execFileSync(binary(), [
    "-v", join(work, "cohort.vcf"),
    "-f", join(work, "ref.fasta"),
    "--gff", join(work, "genes.gff3"),
    "--sample", "all",
    "--report", join(work, "report.html"),
  ], { cwd: work, stdio: "ignore" });
  return join(work, "report.html");
}

/// The page's own functions, lifted out and given its own data.
function pageLogic(reportPath) {
  const page = readFileSync(reportPath, "utf8");
  const data = JSON.parse(page.match(/id="report-data"[^>]*>([\s\S]*?)<\/script>/)[1]);
  const body = [...page.matchAll(/<script>([\s\S]*?)<\/script>/g)].map((m) => m[1]).join("\n");

  const sandbox = { D: data.dict, ROWS: data.rows, console, clampView() {} };
  createContext(sandbox);
  // Everything from the row accessors down to the filter, which is pure.
  runInContext(body.slice(body.indexOf("var C = {"), body.indexOf("function recompute()")), sandbox);
  // The matrix builder, and the two constants it reads.
  for (const name of ["var NT", "var MAX_SITES"]) {
    const at = body.indexOf(name);
    runInContext(body.slice(at, body.indexOf("\n", at)), sandbox);
  }
  // The density lanes, whose predicates key on positions in the dictionaries.
  runInContext(
    body.slice(body.indexOf("var T_SNP ="), body.indexOf("function activeLanes()")), sandbox);
  const from = body.indexOf("var mtx = {");
  const to = body.indexOf("mtx.groups =", from);
  runInContext(body.slice(from, to) + "\n}", sandbox);
  return { data, sandbox };
}

function tsvRowsBySample(work) {
  const out = {};
  for (const file of readdirSync(work)) {
    const m = file.match(/^cohort\.sample_(.+)\.MNV\.tsv$/);
    if (!m) continue;
    const lines = readFileSync(join(work, file), "utf8").trim().split("\n");
    const header = lines[0].split("\t");
    out[m[1]] = lines.slice(1).map((line) => {
      const cells = line.split("\t");
      return Object.fromEntries(header.map((h, i) => [h, cells[i]]));
    });
  }
  return out;
}

// ---------------------------------------------------------------------------
// The report
// ---------------------------------------------------------------------------

function checkReport(work) {
  const reportPath = buildReport(work);
  const { sandbox } = pageLogic(reportPath);
  const { ROWS, passes, state, sortKey, COLS, geneOf, impOf, sampleOf } = sandbox;
  const perSample = tsvRowsBySample(work);
  const allRows = Object.values(perSample).flat();
  const viewSize = () => ROWS.filter((r) => passes(r)).length;

  check("report: with no filter the view holds every row of the TSVs", viewSize(), allRows.length);

  for (const sample of Object.keys(perSample).sort()) {
    state.col = { sample: { set: { [sample]: true } } };
    check(`report: filtering on sample ${sample} gives that sample's rows`,
      viewSize(), perSample[sample].length);
  }
  state.col = {};

  const byImpact = {};
  for (const row of allRows) byImpact[row.Impact] = (byImpact[row.Impact] || 0) + 1;
  for (const impact of Object.keys(byImpact).sort()) {
    state.col = { impact: { set: { [impact]: true } } };
    check(`report: filtering on impact ${impact} gives that many rows`, viewSize(), byImpact[impact]);
  }
  state.col = {};

  const byGene = {};
  for (const row of allRows) byGene[row.Gene] = (byGene[row.Gene] || 0) + 1;
  for (const gene of Object.keys(byGene).sort()) {
    state.q = gene.toLowerCase();
    check(`report: searching "${gene}" finds its rows`, viewSize(), byGene[gene]);
  }
  state.q = "";

  const posCol = COLS.find((c) => c.k === "pos");
  const order = ROWS.map((_, i) => i)
    .sort((a, b) => {
      const ka = sortKey(ROWS[a], posCol), kb = sortKey(ROWS[b], posCol);
      return ka < kb ? -1 : ka > kb ? 1 : a - b;
    });
  let ascending = true;
  for (let i = 1; i < order.length; i++) {
    if (sortKey(ROWS[order[i - 1]], posCol) > sortKey(ROWS[order[i]], posCol)) ascending = false;
  }
  check("report: sorting on position leaves the positions ascending", ascending, true);

  // The KPI headline is computed in the page from the rows it shows.
  const genesWithACall = new Set(
    ROWS.filter((r) => geneOf(r) && geneOf(r) !== "intergenic").map((r) => geneOf(r)));
  const genesInTsv = new Set(allRows.map((r) => r.Gene).filter((g) => g !== "intergenic"));
  check("report: the gene count matches the TSVs, excluding the intergenic placeholder",
    genesWithACall.size, genesInTsv.size);
  check("report: the HIGH count matches the TSVs",
    ROWS.filter((r) => impOf(r) === "HIGH").length,
    allRows.filter((r) => r.Impact === "HIGH").length);
  check("report: every row is attributed to a sample that has it",
    [...new Set(ROWS.map(sampleOf))].sort(), Object.keys(perSample).sort());

  // The matrix, over both orders of the table. A site belongs to every gene that
  // annotates it, and nothing about it may follow how the table is sorted.
  const genesPerPosition = {};
  for (const row of allRows) {
    for (const p of row.Positions.split(",").map((v) => v.trim())) {
      (genesPerPosition[p] = genesPerPosition[p] || new Set()).add(row.Gene);
    }
  }
  const matrixFor = (view) => {
    sandbox.view = view;
    runInContext("buildMatrix()", sandbox);
    return {
      sites: sandbox.mtx.sites.map((s) => `${s.pos}:${s.genes.join(",")}`),
      bars: sandbox.mtx.geneSpans.map((g) => `${g.gene}:${g.lo}-${g.hi}`),
    };
  };
  const forward = matrixFor(ROWS.map((_, i) => i));
  const reversed = matrixFor(ROWS.map((_, i) => i).reverse());
  check("report: the matrix does not change when the table is sorted the other way",
    forward, reversed);

  const asMatrix = Object.fromEntries(
    sandbox.mtx.sites.map((s) => [String(s.pos), [...s.genes].sort()]));
  const asTsv = Object.fromEntries(
    Object.entries(genesPerPosition)
      .filter(([p]) => p in asMatrix)
      .map(([p, g]) => [p, [...g].sort()]));
  check("report: a matrix site keeps every gene the TSVs give that base", asMatrix, asTsv);

  // The density lanes. Each predicate keys on a position in a dictionary, so
  // every lane is counted back against the bases the TSVs give it. The lanes
  // count changed bases, not rows, because the track is drawn over coordinates.
  sandbox.view = ROWS.map((_, i) => i);
  runInContext("buildMatrix()", sandbox);
  const changed = [];
  for (const [sample, rows] of Object.entries(perSample)) {
    for (const row of rows) {
      if (row.Chromosome !== sandbox.mtx.contig) continue;
      for (const _ of row.Positions.split(",")) {
        changed.push({ type: row["Variant Type"], impact: row.Impact, sample });
      }
    }
  }
  const laneWants = {
    all: changed.length,
    snp: changed.filter((c) => c.type === "SNP").length,
    mnv: changed.filter((c) => c.type === "MNV" || c.type === "SNP/MNV").length,
    indel: changed.filter((c) => c.type === "INDEL").length,
    high: changed.filter((c) => c.impact === "HIGH").length,
  };
  for (const lane of sandbox.LANES.filter((l) => !l.distinct)) {
    // A lane that counts nothing on both sides would agree without proving
    // anything, so the fixture is required to feed every one of them.
    check(`report: the "${lane.label}" lane has something to count`, laneWants[lane.key] > 0, true);
    check(`report: the "${lane.label}" lane counts the bases the TSVs give it`,
      sandbox.mtx.calls.filter((c) => lane.hit(c[1], c[2])).length, laneWants[lane.key]);
  }
  check("report: the samples lane sees every sample the TSVs name",
    new Set(sandbox.mtx.calls.map((c) => c[3])).size,
    new Set(changed.map((c) => c.sample)).size);
}

// ---------------------------------------------------------------------------
// The desktop form
// ---------------------------------------------------------------------------

/// The object literal starting at `open`, up to its matching brace.
function braceMatch(source, open) {
  let depth = 0;
  for (let i = open; i < source.length; i++) {
    const c = source[i];
    if (c === "{" || c === "[") depth++;
    else if (c === "}" || c === "]") {
      depth--;
      if (depth === 0) return source.slice(open, i + 1);
    }
  }
  throw new Error("unterminated object literal");
}

function literal(source, marker, opener) {
  const at = source.indexOf(marker);
  const eq = source.indexOf("=", at + marker.length - 1);
  const open = source.indexOf(opener, eq);
  let depth = 0;
  for (let i = open; i < source.length; i++) {
    const c = source[i];
    if (c === "{" || c === "[") depth++;
    else if (c === "}" || c === "]") {
      depth--;
      if (depth === 0) return source.slice(open, i + 1);
    }
  }
  throw new Error(`cannot slice ${marker}`);
}

function checkPresets() {
  const types = readFileSync(join(REPO, "frontend", "src", "types.ts"), "utf8");
  const form = readFileSync(
    join(REPO, "frontend", "src", "components", "ParameterForm.tsx"), "utf8");

  const sandbox = {};
  createContext(sandbox);
  runInContext("var DEFAULT_CONFIG = " + literal(types, "DEFAULT_CONFIG", "{"), sandbox);
  runInContext("var PRESETS = " + literal(types, "BUILT_IN_PRESETS", "["), sandbox);
  const { DEFAULT_CONFIG, PRESETS } = sandbox;

  // The merge the form actually performs. The object literal it hands to
  // onChange is sliced out of the component and evaluated with the same names in
  // scope, so this runs the form's own code rather than a copy of it: a rule
  // deleted there fails here.
  const call = form.indexOf("onChange({", form.indexOf("const preset = BUILT_IN_PRESETS.find"));
  const mergeLiteral = braceMatch(form, form.indexOf("{", call));
  const applyPreset = (config, preset) => {
    const scope = { DEFAULT_CONFIG, config, preset };
    createContext(scope);
    return runInContext(`(${mergeLiteral})`, scope);
  };

  const preserved = [...mergeLiteral.matchAll(/^\s+([a-zA-Z]+):/gm)].map((m) => m[1]);

  check("form: the preserved list covers the output formats",
    ["outputTsv", "outputVcf"].every((f) => preserved.includes(f)), true);

  const chose = { ...DEFAULT_CONFIG, outputTsv: false, outputVcf: true, vcfGz: true };
  for (const preset of PRESETS) {
    const after = applyPreset(chose, preset);
    check(`form: preset "${preset.name}" keeps a chosen VCF output`, after.outputVcf, true);
    check(`form: preset "${preset.name}" leaves no compression without a VCF`,
      after.vcfGz && !after.outputVcf, false);
  }

  // And a preset that states a format still imposes it.
  const stating = PRESETS.filter((p) => p.config.outputVcf !== undefined);
  check("form: at least one preset states an output format", stating.length > 0, true);
  for (const preset of stating) {
    const after = applyPreset({ ...DEFAULT_CONFIG, outputVcf: false }, preset);
    check(`form: preset "${preset.name}" imposes the format it states`,
      after.outputVcf, preset.config.outputVcf);
  }
}

// ---------------------------------------------------------------------------

const work = mkdtempSync(join(tmpdir(), "get_mnv_js_"));
try {
  checkReport(work);
  checkPresets();
} finally {
  rmSync(work, { recursive: true, force: true });
}

console.log(`\n${checks - failures}/${checks} checks passed`);
if (failures) {
  console.error(`${failures} check(s) failed`);
  process.exit(1);
}
