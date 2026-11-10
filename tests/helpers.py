"""Deterministic sequence generators for the test suite."""

from __future__ import annotations

import random


def make_seq(n: int, seed: int) -> str:
    """A reproducible random nucleotide sequence."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq: str, rate: float, seed: int) -> str:
    """Substitute each base with probability *rate*, reproducibly."""
    rng = random.Random(seed)
    return "".join(
        rng.choice([b for b in "ACGT" if b != base]) if rng.random() < rate else base
        for base in seq
    )
