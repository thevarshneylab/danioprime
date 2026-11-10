/**
 * Print the Pyodide version that index.html is pinned to.
 *
 * One source of truth: CI installs whatever the page loads, so the smoke test
 * can never pass against a different runtime than visitors get.
 */
import { readFileSync } from "node:fs";

// index.html is a sibling of this script in web/. It used to live at the repo
// root; when it moved, this path did not, and CI broke on a file-not-found
// that no local run reproduced because the script had last been exercised
// before the move. Hence the test that runs this script against the real file.
const html = readFileSync(new URL("./index.html", import.meta.url), "utf8");
const match = html.match(/\/v([0-9][^/"]*)\/full\//);
if (!match) {
  console.error("could not find a pinned Pyodide version in index.html");
  process.exit(1);
}
process.stdout.write(match[1] + "\n");
