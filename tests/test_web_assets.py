"""The browser build's static assets must stay in step with the package.

A stale wheel beside `index.html` means the website silently runs old code --
which is exactly what happened once during development, hence these checks.
The wheel's *contents* are verified in CI (where a fresh wheel can be built and
diffed); here we check everything that can be checked from the repo alone.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path

import pytest

import danioprime
from danioprime.webapp import MODES

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "index.html"
WHEEL = ROOT / f"danioprime-{danioprime.__version__}-py3-none-any.whl"


@pytest.fixture(scope="module")
def html() -> str:
    if not INDEX.exists():
        pytest.skip("index.html not present")
    return INDEX.read_text()


class TestWheelAsset:
    def test_the_wheel_for_this_version_is_committed(self):
        assert WHEEL.exists(), (
            f"{WHEEL.name} is missing. Rebuild it with "
            f"`pip wheel . --no-deps -w . ` after bumping the version."
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
        on_disk = {str(p.relative_to(src)) for p in src.rglob("*.py")}
        missing = on_disk - in_wheel
        extra = in_wheel - on_disk
        assert not missing, f"the committed wheel is stale; it lacks: {sorted(missing)}"
        assert not extra, f"the committed wheel holds modules no longer in src/: {sorted(extra)}"

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

    def test_demo_protospacer_is_present_in_the_demo_locus(self, html):
        """The 'load example' button must produce an input that actually works."""
        match = re.search(r'\$\("g-proto"\)\.value = "([ACGT]+)"', html)
        assert match, "the genotype demo does not set a protospacer"
        from danioprime.assays.genotyping import locate_protospacer
        from danioprime.sequence import read_single_fasta

        locus = read_single_fasta(ROOT / "examples" / "kcnh6a_locus.fa")
        sites = locate_protospacer(locus.seq, match.group(1))
        assert len(sites) == 1, (
            f"the demo protospacer matches {len(sites)} sites in the demo locus"
        )

    def test_output_is_html_escaped(self, html):
        """Sequences and names are rendered into innerHTML, so they must be escaped."""
        assert re.search(r"const esc\s*=", html), "no escaping helper defined"
        # Every interpolation of a model field in a template literal should be
        # wrapped in esc(...) -- spot-check the sequence fields, which are the
        # ones carrying user input.
        for field in ("p.forward", "p.reverse", "p.name"):
            assert f"esc({field})" in html, f"{field} is interpolated without escaping"

    def test_page_states_that_nothing_is_uploaded(self, html):
        assert "ever leaves this machine" in html

    def test_page_carries_the_validation_caveat(self, html):
        assert "Not validated against bench outcomes" in html

    def test_no_external_script_beyond_pyodide(self, html):
        srcs = re.findall(r'<script[^>]*\ssrc="([^"]+)"', html)
        assert srcs == [], f"unexpected external scripts: {srcs}"
        imports = re.findall(r'import\("?([^")]+)', html)
        for target in imports:
            assert "PYODIDE" in target or target.startswith("node:"), target
