"""Pluggable thermodynamic backend.

``builtin`` (default) uses the pure-Python nearest-neighbour implementation in
:mod:`danioprime.thermo.nn` and :mod:`danioprime.thermo.structure`. It has no
compiled dependencies, so results are byte-for-byte reproducible across
platforms and CI.

``primer3`` delegates Tm, hairpin and dimer calculations to ``primer3-py``,
which wraps the same libprimer3 code that Primer3/Primer3Plus use. Prefer it
when you need numbers that match a Primer3 run exactly. Install with::

    pip install "danioprime[primer3]"

Both backends expose the same :class:`ThermoBackend` interface so the rest of
the package is agnostic to the choice.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from ..errors import DanioPrimeError
from . import nn, structure
from .nn import DEFAULT_CONDITIONS, SaltConditions


@dataclass(frozen=True, slots=True)
class StructureResult:
    """Backend-neutral summary of a secondary-structure calculation."""

    dg: float
    score: float
    tm: float = float("nan")
    three_prime_anchored: bool = False
    alignment: str = ""


class ThermoBackend(Protocol):
    """Interface implemented by every thermodynamic backend."""

    name: str

    def tm(self, seq: str) -> float: ...
    def delta_g(self, seq: str, temp_c: float = 37.0) -> float: ...
    def hairpin(self, seq: str) -> StructureResult: ...
    def self_dimer(self, seq: str) -> StructureResult: ...
    def hetero_dimer(self, a: str, b: str) -> StructureResult: ...
    def end_stability(self, seq: str) -> float: ...


class BuiltinBackend:
    """Pure-Python nearest-neighbour backend."""

    name = "builtin"

    def __init__(self, conditions: SaltConditions = DEFAULT_CONDITIONS, anneal_temp_c: float = 37.0):
        self.conditions = conditions
        self.anneal_temp_c = anneal_temp_c

    def tm(self, seq: str) -> float:
        return nn.tm(seq, self.conditions)

    def delta_g(self, seq: str, temp_c: float = 37.0) -> float:
        return nn.delta_g(seq, temp_c, self.conditions)

    def end_stability(self, seq: str) -> float:
        return nn.end_stability(seq, temp_c=self.anneal_temp_c)

    def _wrap(self, interaction: structure.Interaction) -> StructureResult:
        return StructureResult(
            dg=interaction.dg,
            score=interaction.score,
            three_prime_anchored=interaction.three_prime_anchored,
            alignment=interaction.alignment,
        )

    def hairpin(self, seq: str) -> StructureResult:
        return self._wrap(
            structure.hairpin(seq, temp_c=self.anneal_temp_c, conditions=self.conditions)
        )

    def self_dimer(self, seq: str) -> StructureResult:
        return self._wrap(
            structure.self_dimer(seq, temp_c=self.anneal_temp_c, conditions=self.conditions)
        )

    def hetero_dimer(self, a: str, b: str) -> StructureResult:
        return self._wrap(
            structure.hetero_dimer(a, b, temp_c=self.anneal_temp_c, conditions=self.conditions)
        )


class Primer3Backend:
    """Backend delegating to ``primer3-py`` (libprimer3)."""

    name = "primer3"

    def __init__(self, conditions: SaltConditions = DEFAULT_CONDITIONS, anneal_temp_c: float = 37.0):
        try:
            import primer3  # noqa: PLC0415
        except ImportError as exc:  # pragma: no cover - exercised only without the extra
            raise DanioPrimeError(
                "the 'primer3' backend requires primer3-py; install with "
                "`pip install \"danioprime[primer3]\"`"
            ) from exc
        self._p3 = primer3
        self.conditions = conditions
        self.anneal_temp_c = anneal_temp_c
        self._kwargs = dict(
            mv_conc=conditions.monovalent_mM,
            dv_conc=conditions.divalent_mM,
            dntp_conc=conditions.dntp_mM,
            dna_conc=conditions.primer_nM,
        )

    def tm(self, seq: str) -> float:
        return float(self._p3.calc_tm(seq, **self._kwargs))

    def delta_g(self, seq: str, temp_c: float = 37.0) -> float:
        # libprimer3 reports dG in cal/mol for the heterodimer of seq and its
        # reverse complement, which is the closest available analogue.
        from ..sequence import revcomp  # noqa: PLC0415

        res = self._p3.calc_heterodimer(seq, revcomp(seq), temp_c=temp_c, **self._kwargs)
        return float(res.dg) / 1000.0

    def end_stability(self, seq: str) -> float:
        return nn.end_stability(seq, temp_c=self.anneal_temp_c)

    def _wrap(self, res) -> StructureResult:
        if not getattr(res, "structure_found", False):
            return StructureResult(dg=0.0, score=0.0, tm=float("nan"))
        return StructureResult(
            dg=float(res.dg) / 1000.0,
            score=0.0,
            tm=float(res.tm),
            alignment=getattr(res, "ascii_structure", "") or "",
        )

    def hairpin(self, seq: str) -> StructureResult:
        return self._wrap(self._p3.calc_hairpin(seq, temp_c=self.anneal_temp_c, **self._kwargs))

    def self_dimer(self, seq: str) -> StructureResult:
        return self._wrap(self._p3.calc_homodimer(seq, temp_c=self.anneal_temp_c, **self._kwargs))

    def hetero_dimer(self, a: str, b: str) -> StructureResult:
        return self._wrap(
            self._p3.calc_heterodimer(a, b, temp_c=self.anneal_temp_c, **self._kwargs)
        )


_BACKENDS = {"builtin": BuiltinBackend, "primer3": Primer3Backend}

#: Names accepted by :func:`get_backend`.
BACKEND_NAMES = tuple(_BACKENDS)


def get_backend(
    name: str = "builtin",
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    anneal_temp_c: float = 37.0,
) -> ThermoBackend:
    """Instantiate a backend by name."""
    try:
        cls = _BACKENDS[name]
    except KeyError as exc:
        raise DanioPrimeError(
            f"unknown thermo backend {name!r}; choose from {sorted(_BACKENDS)}"
        ) from exc
    return cls(conditions=conditions, anneal_temp_c=anneal_temp_c)


def primer3_available() -> bool:
    """True if ``primer3-py`` can be imported."""
    try:
        import primer3  # noqa: F401, PLC0415

        return True
    except ImportError:
        return False
