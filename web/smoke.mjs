/**
 * Run every danioprime web mode under Pyodide — the same CPython-on-WebAssembly
 * runtime the browser page uses.
 *
 * This is the only check that the shipped `index.html` can actually execute the
 * committed wheel. Everything else tests the Python under native CPython, which
 * is not the same interpreter build and not the same Python version.
 *
 * Usage, from the repository root:
 *
 *   npm install pyodide@$(node web/pyodide_version.mjs)
 *   node web/smoke.mjs
 *
 * Options via environment:
 *   PYODIDE_INDEX  path or URL to the Pyodide distribution
 *                  (default: ./node_modules/pyodide/)
 *   DANIOPRIME_WHEEL  wheel to load (default: the one beside index.html)
 */

import { readFileSync, readdirSync } from "node:fs";
import { join, resolve } from "node:path";
import { pathToFileURL } from "node:url";

const ROOT = new URL("..", import.meta.url).pathname;
const INDEX_URL = process.env.PYODIDE_INDEX || "./node_modules/pyodide/";
const EXAMPLES = join(ROOT, "examples");

function findWheel() {
  if (process.env.DANIOPRIME_WHEEL) return process.env.DANIOPRIME_WHEEL;
  const hit = readdirSync(ROOT).find((f) => /^danioprime-.*-py3-none-any\.whl$/.test(f));
  if (!hit) throw new Error("no danioprime wheel at the repository root");
  return join(ROOT, hit);
}

const read = (name) => readFileSync(join(EXAMPLES, name), "utf8");

const CASES = [
  {
    name: "genotype",
    payload: () => ({
      mode: "genotype", preset: "genotyping-gel", overrides: { max_pairs: "2" },
      template: read("kcnh6a_locus.fa"), avoid: read("kcnh6b_locus.fa"),
      protospacer: "AAGCCCTGTTATGTGCAACG", nuclease: "SpCas9",
      protected: "60", min_flank: "40", max_flank: "400", tails: "fluorescent",
    }),
    expect: (r) => {
      if (!r.pairs.length) return "no pairs";
      if (!r.pairs.every((p) => p.specific)) return "a reported pair cross-amplifies the panel";
      if (!r.order_sheet?.includes("M13F_universal")) return "fluorescent tails produced no companion oligo";
      if (!r.tables.some((t) => /allele sizes/i.test(t.title))) return "no allele-size table";
      return null;
    },
  },
  {
    name: "qpcr",
    payload: () => ({
      mode: "qpcr", preset: "qpcr", overrides: {},
      template: read("kcnh6a_transcript.fa"), avoid: read("kcnh6b_transcript.fa"),
      junctions: [210, 540], require_junction_spanning: true, junction_flank: "60",
    }),
    expect: (r) => {
      if (!r.pairs.length) return "no pairs";
      if (!r.pairs.some((p) => p.spans_junction != null)) return "no junction-spanning primer";
      if (!r.tables.some((t) => /protection/i.test(t.title))) return "no gDNA protection table";
      return null;
    },
  },
  {
    name: "dcaps",
    payload: () => ({
      mode: "dcaps", preset: "genotyping-gel",
      overrides: { min_product_size: "150", max_product_size: "450" },
      template: read("kcnh6a_locus.fa"), variant_position: "1768", variant_base: "C",
      cut_allele: "variant", max_introduced_mismatches: "2",
    }),
    expect: (r) => (r.tables?.[0]?.rows?.length ? null : "no dCAPS assays"),
  },
  {
    name: "design",
    payload: () => ({
      mode: "design", preset: "default", overrides: { max_pairs: "2" },
      template: read("kcnh6a_locus.fa"), avoid: "", targets: "1700-1800",
      excludes: "", tails: "none",
    }),
    expect: (r) => (r.pairs.length ? null : "no pairs"),
  },
  {
    name: "check",
    payload: () => ({
      mode: "check", preset: "default", overrides: {},
      // This reverse primer sits in the block conserved with the ohnolog, so a
      // correct run must flag it as priming on the panel.
      primers: "R_conserved=CTATCTGTATGCCGCGATG",
      avoid: read("kcnh6b_locus.fa"),
    }),
    expect: (r) => {
      const panel = r.tables.find((t) => /avoid-panel/.test(t.title));
      if (!panel) return "no panel table";
      if (panel.rows[0][2] !== "YES") return "failed to flag a primer that matches the ohnolog";
      return null;
    },
  },
];

// Resolve the runtime from PYODIDE_INDEX rather than by bare package name:
// a bare `import "pyodide"` resolves against this script's directory, which
// breaks whenever node_modules lives somewhere else (CI, a scratch dir).
const entry = INDEX_URL.startsWith("http")
  ? INDEX_URL + "pyodide.mjs"
  : pathToFileURL(resolve(process.cwd(), INDEX_URL, "pyodide.mjs")).href;
const { loadPyodide } = await import(entry);

const t0 = Date.now();
const py = await loadPyodide({ indexURL: INDEX_URL });
const wheel = findWheel();
py.FS.writeFile("/danioprime.whl", readFileSync(wheel));
// A wheel is a zip, so zipimport serves it straight from sys.path: a pure
// Python package with no dependencies needs no package manager at all.
py.runPython(`import sys; sys.path.insert(0, "/danioprime.whl")`);

const version = py.runPython(`import danioprime; danioprime.__version__`);
const pyver = py.runPython(`import sys; ".".join(map(str, sys.version_info[:3]))`);
console.log(`danioprime ${version} on CPython ${pyver} (WebAssembly), booted in ${((Date.now() - t0) / 1000).toFixed(1)}s`);
console.log(`wheel: ${wheel.split("/").pop()}\n`);

function call(payload) {
  py.globals.set("_payload", JSON.stringify(payload));
  return JSON.parse(
    py.runPython(`from danioprime.webapp import run_json\nrun_json(_payload)`)
  );
}

let failures = 0;
for (const c of CASES) {
  const t = Date.now();
  let result;
  try {
    result = call(c.payload());
  } catch (err) {
    console.error(`FAIL ${c.name}: the runtime raised ${err.message?.split("\n")[0]}`);
    failures++;
    continue;
  }
  const secs = ((Date.now() - t) / 1000).toFixed(1);
  if (!result.ok) {
    console.error(`FAIL ${c.name} (${secs}s): ${result.error}`);
    failures++;
    continue;
  }
  const why = c.expect(result);
  if (why) {
    console.error(`FAIL ${c.name} (${secs}s): ${why}`);
    failures++;
    continue;
  }
  const n = result.pairs.length || result.tables?.[0]?.rows?.length || 0;
  console.log(`ok   ${c.name.padEnd(9)} ${secs}s  ${n} result(s)`);
}

// An unusable input must come back as a readable message, never an exception:
// the page has no other way to tell the user what went wrong.
const bad = call({
  mode: "genotype",
  template: read("kcnh6a_locus.fa"),
  protospacer: "A".repeat(20),
});
if (bad.ok || !bad.error || bad.error.includes("Traceback")) {
  console.error("FAIL error path: an invalid protospacer did not produce a clean message");
  failures++;
} else {
  console.log(`ok   error path returns: "${bad.error.slice(0, 60)}…"`);
}

console.log(failures ? `\n${failures} failure(s)` : "\nall web modes pass under Pyodide");
process.exit(failures ? 1 : 0);
