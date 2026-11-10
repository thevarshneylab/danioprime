#!/usr/bin/env python3
"""Verify the wheel committed for the browser build matches the current source.

`index.html` loads a wheel committed beside it, so a stale wheel means the
website silently runs old code. This happened once during development, which is
why it is now checked.

Contents are compared rather than bytes: zip entry timestamps differ between
builds even when the sources are identical.

    python web/check_wheel.py          # check
    python web/check_wheel.py --fix    # rebuild the committed wheel

Exits non-zero when the committed wheel is missing, stale, or not pure Python.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src" / "danioprime"


def version() -> str:
    text = (SRC / "__init__.py").read_text()
    for line in text.splitlines():
        if line.startswith("__version__"):
            return line.split("=")[1].strip().strip('"').strip("'")
    raise SystemExit("could not read __version__ from src/danioprime/__init__.py")


def wheel_name(ver: str) -> str:
    return f"danioprime-{ver}-py3-none-any.whl"


def build(dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [sys.executable, "-m", "pip", "wheel", str(ROOT), "--no-deps", "-w", str(dest), "-q"],
        check=True,
    )
    wheels = list(dest.glob("danioprime-*.whl"))
    if len(wheels) != 1:
        raise SystemExit(f"expected one wheel in {dest}, found {len(wheels)}")
    return wheels[0]


def source_digests(wheel: Path) -> dict[str, str]:
    """sha256 of every package file inside the wheel, keyed by relative path."""
    out: dict[str, str] = {}
    with zipfile.ZipFile(wheel) as zf:
        for name in zf.namelist():
            if name.startswith("danioprime/") and not name.endswith("/"):
                out[name] = hashlib.sha256(zf.read(name)).hexdigest()
    return out


def disk_digests() -> dict[str, str]:
    out: dict[str, str] = {}
    for path in sorted(SRC.rglob("*")):
        if path.is_file() and path.suffix == ".py":
            rel = f"danioprime/{path.relative_to(SRC).as_posix()}"
            out[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fix", action="store_true", help="rebuild the committed wheel")
    args = ap.parse_args()

    ver = version()
    name = wheel_name(ver)
    committed = ROOT / name

    if args.fix:
        with tempfile.TemporaryDirectory() as tmp:
            fresh = build(Path(tmp))
            for old in ROOT.glob("danioprime-*-py3-none-any.whl"):
                old.unlink()
            shutil.copy2(fresh, committed)
        print(f"rebuilt {name}")
        return 0

    problems: list[str] = []
    if not committed.exists():
        print(f"FAIL {name} is not committed at the repo root")
        print("     fix: python web/check_wheel.py --fix")
        return 1

    stray = [p.name for p in ROOT.glob("danioprime-*.whl") if p.name != name]
    if stray:
        problems.append(f"wheels for other versions are still committed: {stray}")

    on_disk = disk_digests()
    in_wheel = source_digests(committed)
    py_in_wheel = {k: v for k, v in in_wheel.items() if k.endswith(".py")}

    missing = sorted(set(on_disk) - set(py_in_wheel))
    extra = sorted(set(py_in_wheel) - set(on_disk))
    changed = sorted(k for k in set(on_disk) & set(py_in_wheel) if on_disk[k] != py_in_wheel[k])

    if missing:
        problems.append(f"modules missing from the wheel: {missing}")
    if extra:
        problems.append(f"modules in the wheel but not in src/: {extra}")
    if changed:
        problems.append(f"modules whose contents differ: {changed}")

    if problems:
        print(f"FAIL the committed {name} does not match src/danioprime")
        for p in problems:
            print(f"     - {p}")
        print("     fix: python web/check_wheel.py --fix")
        return 1

    print(f"ok   {name} matches src/danioprime ({len(py_in_wheel)} modules)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
