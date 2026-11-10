"""Re-export the site builder for tests, which cannot import from ``web/``."""

from __future__ import annotations

import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "_danioprime_build_site", Path(__file__).resolve().parents[1] / "web" / "build_site.py"
)
_module = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_module)

BUNDLED = _module.BUNDLED
planned_files = _module.planned_files
