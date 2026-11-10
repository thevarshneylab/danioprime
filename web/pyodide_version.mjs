/**
 * Print the Pyodide version that index.html is pinned to.
 *
 * One source of truth: CI installs whatever the page loads, so the smoke test
 * can never pass against a different runtime than visitors get.
 */
import { readFileSync } from "node:fs";

const html = readFileSync(new URL("../index.html", import.meta.url), "utf8");
const match = html.match(/\/v([0-9][^/"]*)\/full\//);
if (!match) {
  console.error("could not find a pinned Pyodide version in index.html");
  process.exit(1);
}
process.stdout.write(match[1] + "\n");
