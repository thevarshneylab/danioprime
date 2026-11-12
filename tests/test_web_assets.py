"""The browser build's static assets must stay in step with the package.

A stale wheel beside `index.html` means the website silently runs old code --
which is exactly what happened once during development, hence these checks.
The wheel's *contents* are verified in CI (where a fresh wheel can be built and
diffed); here we check everything that can be checked from the repo alone.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

import danioprime
from danioprime.webapp import MODES

ROOT = Path(__file__).resolve().parents[1]
#: Source of the page; docs/ is the build output GitHub Pages publishes.
INDEX = ROOT / "web" / "index.html"
DOCS = ROOT / "docs"
WHEEL = DOCS / f"danioprime-{danioprime.__version__}-py3-none-any.whl"


@pytest.fixture(scope="module")
def html() -> str:
    if not INDEX.exists():
        pytest.skip("index.html not present")
    return INDEX.read_text()


class TestWheelAsset:
    def test_the_wheel_for_this_version_is_committed(self):
        assert WHEEL.exists(), (
            f"docs/{WHEEL.name} is missing. Rebuild it with "
            f"`python web/check_wheel.py --fix && python web/build_site.py`."
        )

    def test_index_points_at_that_exact_wheel(self, html):
        match = re.search(r'const WHEEL\s*=\s*"\./([^"]+)"', html)
        assert match, "index.html does not declare a WHEEL constant"
        assert match.group(1) == WHEEL.name, (
            f"index.html loads {match.group(1)} but the committed wheel is {WHEEL.name}"
        )

    def test_wheel_is_pure_python(self):
        assert WHEEL.name.endswith("-py3-none-any.whl")

    def test_wheel_contains_the_web_entry_point(self):
        names = zipfile.ZipFile(WHEEL).namelist()
        assert "danioprime/webapp.py" in names
        assert "danioprime/__init__.py" in names

    def test_wheel_module_list_matches_the_source_tree(self):
        """Catches a wheel built before a module was added or removed."""
        in_wheel = {
            n[len("danioprime/") :]
            for n in zipfile.ZipFile(WHEEL).namelist()
            if n.startswith("danioprime/") and n.endswith(".py")
        }
        src = ROOT / "src" / "danioprime"
        # as_posix(), not str(): on Windows str() yields backslashes while a
        # wheel always stores forward slashes, and the mismatch would report a
        # perfectly fresh wheel as stale.
        on_disk = {p.relative_to(src).as_posix() for p in src.rglob("*.py")}
        missing = on_disk - in_wheel
        extra = in_wheel - on_disk
        assert not missing, f"the committed wheel is stale; it lacks: {sorted(missing)}"
        assert not extra, f"the committed wheel holds modules no longer in src/: {sorted(extra)}"

    def test_module_paths_are_compared_posix_style(self):
        """Regression guard: a backslash here reported a fresh wheel as stale.

        The first Windows CI run failed on exactly this, and a staleness check
        that cries wolf is worse than none -- people learn to ignore it.
        """
        src = ROOT / "src" / "danioprime"
        paths = {p.relative_to(src).as_posix() for p in src.rglob("*.py")}
        assert any("/" in p for p in paths), "expected nested modules in the package"
        assert all("\\" not in p for p in paths)

    def test_wheel_declares_no_runtime_dependencies(self):
        """The browser build works precisely because there is nothing to install."""
        meta = zipfile.ZipFile(WHEEL).read(
            f"danioprime-{danioprime.__version__}.dist-info/METADATA"
        ).decode()
        required = [
            line
            for line in meta.splitlines()
            if line.startswith("Requires-Dist:") and "extra ==" not in line
        ]
        assert not required, f"the browser build cannot install dependencies: {required}"


class TestIndexHtml:
    def test_pyodide_is_pinned_to_a_version(self, html):
        match = re.search(r'const PYODIDE\s*=\s*"https://[^"]*/v([0-9][^/"]*)/full/"', html)
        assert match, "the Pyodide URL should be pinned to an explicit version"
        assert re.fullmatch(r"[0-9.]+", match.group(1)), match.group(1)

    def test_every_run_button_names_a_real_mode(self, html):
        modes = set(re.findall(r'data-run="([a-z]+)"', html))
        assert modes, "no run buttons found"
        assert modes <= set(MODES), f"unknown modes in the page: {sorted(modes - set(MODES))}"

    def test_every_mode_is_reachable_from_the_page(self, html):
        modes = set(re.findall(r'data-run="([a-z]+)"', html))
        assert set(MODES) <= modes, f"modes with no UI: {sorted(set(MODES) - modes)}"

    def test_referenced_example_files_exist(self, html):
        names = set(re.findall(r'fetchExample\("([^"]+)"\)', html))
        assert names, "the page offers no examples"
        for name in names:
            assert (ROOT / "examples" / name).exists(), f"examples/{name} is referenced but missing"

    def test_every_referenced_example_is_bundled_into_the_page(self, html):
        """The demo must work without extra network requests."""
        from web_build import BUNDLED  # noqa: PLC0415

        names = set(re.findall(r'fetchExample\("([^"]+)"\)', html))
        assert names <= set(BUNDLED), f"not bundled by build_site.py: {sorted(names - set(BUNDLED))}"

    def test_every_demo_region_lies_inside_the_locus_it_is_set_with(self, html):
        """The 'load example' buttons must produce inputs that actually work.

        Each demo sets a region to interrogate alongside a template; a region
        past the end of that template would make the button produce an error
        instead of a design.
        """
        from danioprime.sequence import read_single_fasta  # noqa: PLC0415

        pairs = [
            ("g-targets", "demogenea_locus.fa"),
            ("x-targets", "myo7aa_locus.fa"),
        ]
        for field, fasta in pairs:
            match = re.search(rf'\$\("{field}"\)\.value = "(\d+)-(\d+)"', html)
            assert match, f"the demo does not set {field}"
            start, end = int(match.group(1)), int(match.group(2))
            locus = read_single_fasta(ROOT / "examples" / fasta)
            assert 0 < start < end <= len(locus.seq), (
                f"{field} is {start}-{end}, outside {fasta} ({len(locus.seq)} nt)"
            )

    def test_the_cut_site_inputs_are_gone(self, html):
        """They would send payload keys the backend no longer reads."""
        for stale in ('id="g-proto"', 'id="g-cut"', 'id="g-nuclease"', 'id="g-protected"',
                      'id="g-nopam"', 'id="g-minflank"', 'id="g-maxflank"'):
            assert stale not in html, f"{stale} is still in the page"

    def test_output_is_html_escaped(self, html):
        """Sequences and names are rendered into innerHTML, so they must be escaped."""
        assert re.search(r"const esc\s*=", html), "no escaping helper defined"
        # Every interpolation of a model field into markup has to go through
        # something that escapes: esc() itself, or one of the helpers that
        # escapes each piece it emits. Found by scanning the render functions
        # rather than by looking for one known spelling, so that rewriting the
        # markup cannot quietly drop the escaping -- and scoped to those
        # functions, because the clipboard and download paths build plain text,
        # where HTML-escaping a sequence would be the bug.
        bodies = re.findall(r"\nfunction (render\w*)\(.*?\n\}", html, re.DOTALL)
        assert bodies, "no render functions found to check"
        markup = "".join(
            m.group(0) for m in re.finditer(r"\nfunction render\w*\(.*?\n\}", html, re.DOTALL)
        )
        safe = ("esc(", "oligoSeq(", "tailChip(")
        for field in ("p.forward", "p.reverse", "p.name"):
            uses = re.findall(r"\$\{([^}]*\b" + re.escape(field) + r"\b[^}]*)\}", markup)
            assert uses, f"{field} is never rendered"
            for use in uses:
                assert use.lstrip().startswith(safe), f"{field} unescaped in ${{{use}}}"

    def test_colours_are_validated_before_reaching_a_style_attribute(self, html):
        """The palette is ours, but a value on its way into CSS gets checked."""
        assert re.search(r"const HEX\s*=\s*/\^#\[0-9a-fA-F\]", html), "no hex guard"
        for helper in ("function oligoSeq", "function tailChip"):
            body = html.split(helper, 1)[1].split("\n}", 1)[0]
            assert "HEX.test(" in body, f"{helper} does not validate its colours"

    def test_page_states_that_nothing_is_uploaded(self, html):
        assert "ever leaves this machine" in html

    def test_page_carries_the_validation_caveat(self, html):
        assert "Not validated against bench outcomes" in html

    def test_no_external_script_beyond_pyodide(self, html):
        srcs = re.findall(r'<script[^>]*\ssrc="([^"]+)"', html)
        remote = [s for s in srcs if not s.startswith("./")]
        assert remote == [], f"unexpected external scripts: {remote}"
        imports = re.findall(r'import\("?([^")]+)', html)
        for target in imports:
            assert "PYODIDE" in target or target.startswith("node:"), target


class TestPublishedSite:
    """`docs/` is a build output; CI must catch it drifting from its sources."""

    def test_docs_is_up_to_date(self):
        import subprocess  # noqa: PLC0415

        result = subprocess.run(
            [sys.executable, str(ROOT / "web" / "build_site.py"), "--check"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    def test_published_index_matches_the_source(self):
        assert (DOCS / "index.html").read_bytes() == INDEX.read_bytes()

    def test_nojekyll_is_present(self):
        """Without it Pages runs Jekyll, which rewrites the page."""
        assert (DOCS / ".nojekyll").exists()

    def test_page_contains_no_liquid_tags(self, html):
        """Belt and braces: `{{` or `{%` would be mangled if Jekyll ever ran."""
        assert "{{" not in html
        assert "{%" not in html

    def test_bundled_examples_match_the_dataset(self):
        import json  # noqa: PLC0415

        js = (DOCS / "examples.js").read_text()
        match = re.search(r"globalThis\.DP_EXAMPLES = (\{.*\});", js, re.S)
        assert match, "examples.js does not assign DP_EXAMPLES"
        data = json.loads(match.group(1))
        assert data, "no examples bundled"
        for name, text in data.items():
            assert text == (ROOT / "examples" / name).read_text(), f"{name} has drifted"

    def test_only_one_wheel_is_published(self):
        wheels = sorted(p.name for p in DOCS.glob("danioprime-*.whl"))
        assert wheels == [WHEEL.name], wheels


class TestSmokeTestCoverage:
    """The Pyodide smoke test must not fall behind the modes that exist.

    It gained `insitu` only after the mode shipped; this stops that recurring.
    """

    def test_every_mode_is_exercised_by_smoke_mjs(self):
        source = (ROOT / "web" / "smoke.mjs").read_text()
        covered = set(re.findall(r'name:\s*"([a-z]+)"', source))
        assert set(MODES) <= covered, f"not smoke-tested: {sorted(set(MODES) - covered)}"

    def test_smoke_test_checks_semantics_not_just_absence_of_errors(self):
        source = (ROOT / "web" / "smoke.mjs").read_text()
        assert source.count("expect:") >= len(MODES)


class TestWebScripts:
    """The helper scripts CI depends on must point at files that exist.

    `pyodide_version.mjs` kept reading the repo root after index.html moved
    into web/; it broke CI and no local run reproduced it, because the script
    had last been exercised before the move.
    """

    @pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
    def test_pyodide_version_script_reports_the_pinned_version(self):
        result = subprocess.run(
            ["node", str(ROOT / "web" / "pyodide_version.mjs")],
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        assert result.returncode == 0, result.stderr
        reported = result.stdout.strip()
        expected = re.search(r"/v([0-9][^/\"]*)/full/", INDEX.read_text()).group(1)
        assert reported == expected

    @pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
    def test_web_scripts_parse(self):
        for script in sorted((ROOT / "web").glob("*.mjs")):
            result = subprocess.run(
                ["node", "--check", str(script)], capture_output=True, text=True
            )
            assert result.returncode == 0, f"{script.name}: {result.stderr}"

    def test_scripts_reference_only_paths_that_exist(self):
        """Every repo-relative path literal in web/ must resolve."""
        for script in sorted((ROOT / "web").glob("*.mjs")):
            source = script.read_text()
            for rel in re.findall(r'new URL\("([^"]+)", import\.meta\.url\)', source):
                resolved = (script.parent / rel).resolve()
                assert resolved.exists(), f"{script.name} points at missing {rel}"


class TestLineEndings:
    """Generated artefacts are compared byte-for-byte, so their newlines are
    part of the contract.

    A Windows checkout translates LF to CRLF in the working tree unless told
    otherwise, which made `build_site.py --check` report docs/ as out of date
    on Windows CI when nothing had changed.
    """

    def test_gitattributes_pins_line_endings(self):
        path = ROOT / ".gitattributes"
        assert path.exists(), "no .gitattributes; Windows checkouts will use CRLF"
        text = path.read_text()
        assert "eol=lf" in text
        assert "*.whl binary" in text

    def test_generated_site_files_have_no_carriage_returns(self):
        for name in ("index.html", "examples.js", ".nojekyll"):
            data = (DOCS / name).read_bytes()
            assert b"\r\n" not in data, f"docs/{name} contains CRLF"

    def test_example_data_has_no_carriage_returns(self):
        for path in sorted((ROOT / "examples").glob("*")):
            if path.suffix in {".fa", ".gtf", ".vcf"}:
                assert b"\r\n" not in path.read_bytes(), f"{path.name} contains CRLF"


class TestOutputDeterminism:
    """Files the package writes must be byte-identical on every platform.

    Python's text mode translates "\\n" to the platform line ending, so on
    Windows the example generator produced CRLF while the committed data was
    LF, and CI's "example data is reproducible" check failed on a dataset that
    had not changed. Anything that writes a file users commit, diff or compare
    has to pin ``newline="\\n"``.

    Test code is exempt: it writes temporary files that are read back by the
    same process on the same platform.
    """

    SCOPE = ("src", "examples", "web")

    def _unpinned(self) -> list[str]:
        import ast

        offenders: list[str] = []
        for root in self.SCOPE:
            for path in sorted((ROOT / root).rglob("*.py")):
                for node in ast.walk(ast.parse(path.read_text())):
                    if not isinstance(node, ast.Call):
                        continue
                    func = node.func
                    name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
                    if name not in {"write_text", "open"}:
                        continue
                    if name == "open":
                        mode = None
                        if len(node.args) > 1 and isinstance(node.args[1], ast.Constant):
                            mode = node.args[1].value
                        for kw in node.keywords:
                            if kw.arg == "mode" and isinstance(kw.value, ast.Constant):
                                mode = kw.value.value
                        if not mode or "w" not in mode or "b" in mode:
                            continue
                    if not any(kw.arg == "newline" for kw in node.keywords):
                        rel = path.relative_to(ROOT)
                        offenders.append(f"{rel}:{node.lineno} {name}()")
        return offenders

    def test_package_text_writes_pin_the_newline(self):
        offenders = self._unpinned()
        assert not offenders, "text writes without newline='\\n':\n  " + "\n  ".join(offenders)

    def test_regenerating_the_example_data_changes_nothing(self):
        """The strongest form: run the generator and diff every output."""
        import subprocess

        examples = ROOT / "examples"
        before = {
            p.name: p.read_bytes()
            for p in sorted(examples.iterdir())
            if p.suffix in {".fa", ".gtf", ".vcf", ".md"}
        }
        result = subprocess.run(
            [sys.executable, "make_example_data.py"],
            cwd=examples,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        changed = [
            name
            for name, data in before.items()
            if (examples / name).read_bytes() != data
        ]
        assert not changed, f"regeneration is not deterministic: {changed}"
