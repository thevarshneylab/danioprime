"""Whether to use colour in the terminal, and how.

Colour here is always redundant: everything it marks is also marked by letter
case, so nothing is lost when it is switched off, and it gets switched off a
lot. A report redirected into a file, piped into ``grep``, read by a
colour-blind user or rendered in the browser's report pane should contain no
escape sequences at all -- a file full of ``\\x1b[38;5;33m`` is worse than a
file with no colour in it, because it looks like corruption.

The rules follow the conventions other command-line tools have settled on, so
that existing habits and existing dotfiles work:

``--color never``
    Never. Nothing else is consulted.
``--color always``
    Always, even into a file or a pipe -- for ``less -R``, or for a build that
    captures output and renders it later.
``--color auto`` (the default)
    Only when the destination is a terminal. ``NO_COLOR`` (any value, per
    no-color.org) switches it off; ``FORCE_COLOR`` switches it on for a pipe;
    ``TERM=dumb`` switches it off. ``NO_COLOR`` wins over ``FORCE_COLOR``,
    because a user who has set it has said what they want about every program
    on the machine, while ``FORCE_COLOR`` is usually set by a CI harness that
    does not know what it is capturing.
"""

from __future__ import annotations

import os
import re
import sys
from typing import Any, Literal

ColourMode = Literal["auto", "always", "never"]

RESET = "\x1b[0m"

#: Matches the escape sequences this module emits, for tests and for anything
#: that needs the width of a painted string.
ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def strip_ansi(text: str) -> str:
    """*text* with colour escapes removed."""
    return ANSI_RE.sub("", text)


def _enable_windows_vt() -> bool:
    """Ask a Windows console to interpret escape sequences. True if it will.

    Windows Terminal and the modern console host both understand them, but only
    once the process has asked, and nothing asks on our behalf. Returns False
    when stdout is not a console at all, which is the right answer for ``auto``.
    """
    if sys.platform != "win32":  # pragma: no cover - exercised on Windows CI
        return True
    try:  # pragma: no cover - exercised on Windows CI
        import ctypes

        kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
        handle = kernel32.GetStdHandle(-11)  # STD_OUTPUT_HANDLE
        mode = ctypes.c_uint32()
        if not kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            return False
        ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x4
        return bool(
            kernel32.SetConsoleMode(handle, mode.value | ENABLE_VIRTUAL_TERMINAL_PROCESSING)
        )
    except Exception:  # pragma: no cover - any failure means "no colour"
        return False


def should_colour(
    mode: ColourMode | str = "auto",
    stream: Any | None = None,
    env: dict[str, str] | None = None,
) -> bool:
    """Resolve *mode* against the destination and the environment.

    *stream* is the destination; ``None`` means "not a terminal", which is what
    callers pass when the report is going to a file.
    """
    if mode == "never":
        return False
    if mode == "always":
        _enable_windows_vt()
        return True
    if mode != "auto":
        raise ValueError(f"unknown colour mode {mode!r}")

    env = os.environ if env is None else env
    if "NO_COLOR" in env:
        return False
    term = env.get("TERM", "")
    if term == "dumb":
        return False
    forced = env.get("FORCE_COLOR")
    if forced is not None and forced != "0":
        return _enable_windows_vt()
    if stream is None:
        return False
    try:
        if not stream.isatty():
            return False
    except Exception:
        return False
    if sys.platform != "win32" and not term:
        # A terminal with no TERM is unusual enough to be a harness pretending.
        return False
    return _enable_windows_vt()


def paint(text: str, ansi: int, *, colour: bool = True, bold: bool = False) -> str:
    """*text* in xterm-256 colour *ansi*, or unchanged when *colour* is false.

    Empty strings are returned unpainted so that a reset code never ends up
    sitting on its own in the output.
    """
    if not colour or not text:
        return text
    weight = "1;" if bold else ""
    return f"\x1b[{weight}38;5;{ansi}m{text}{RESET}"
