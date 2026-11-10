"""Secondary-structure heuristics: self-dimer, hetero-dimer and hairpin.

Two complementary measures are reported for each interaction:

``score``
    The classic Primer3-style complementarity score from the best ungapped
    antiparallel alignment, weighting a G:C pair 2 and an A:T pair 1. Cheap,
    dimensionless, and what most published primer-design rules of thumb are
    phrased in terms of.

``dg``
    Free energy (kcal/mol) of the longest contiguous base-paired stack in that
    alignment, from the nearest-neighbour model. Comparable between oligos of
    different length and the better quantity for hard filters.

A 3'-anchored interaction matters far more than an internal one, because only a
paired 3' end can be extended by the polymerase, so end-anchored interactions
are reported separately.
"""

from __future__ import annotations

from dataclasses import dataclass

from .nn import DEFAULT_CONDITIONS, NN_WC, SaltConditions, salt_corrected_entropy

_PAIR = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
_WEIGHT = {"G": 2.0, "C": 2.0, "A": 1.0, "T": 1.0}

#: Minimum unpaired bases in a hairpin loop.
MIN_LOOP = 3


def _pairs(x: str, y: str) -> bool:
    return (x, y) in _PAIR


def _stack_dg(pairs: list[tuple[str, str]], temp_c: float, conditions: SaltConditions) -> float:
    """Free energy of a contiguous run of base pairs.

    *pairs* is a list of (top, bottom) bases 5'->3' along the top strand.
    """
    if len(pairs) < 2:
        return 0.0
    top = "".join(p[0] for p in pairs)
    dh = 0.0
    ds = 0.0
    for i in range(len(top) - 1):
        step = top[i : i + 2]
        if step not in NN_WC:
            return 0.0
        sdh, sds = NN_WC[step]
        dh += sdh
        ds += sds
    ds = salt_corrected_entropy(ds, len(top), conditions)
    return dh - (temp_c + 273.15) * (ds / 1000.0)


@dataclass(frozen=True, slots=True)
class Interaction:
    """Result of a dimer or hairpin search."""

    score: float
    dg: float
    offset: int
    stack_length: int
    three_prime_anchored: bool
    alignment: str = ""

    @property
    def is_empty(self) -> bool:
        return self.stack_length == 0


EMPTY = Interaction(score=0.0, dg=0.0, offset=0, stack_length=0, three_prime_anchored=False)


def _render(top: str, bottom: str, top_off: int, bottom_off: int) -> str:
    width = max(top_off + len(top), bottom_off + len(bottom))
    t = " " * top_off + top + " " * (width - top_off - len(top))
    b = " " * bottom_off + bottom + " " * (width - bottom_off - len(bottom))
    bars = "".join(
        "|" if t[i] != " " and b[i] != " " and _pairs(t[i], b[i]) else " " for i in range(width)
    )
    return f"5'-{t}-3'\n   {bars}\n3'-{b}-5'"


def _align(
    a: str,
    b_anti: str,
    *,
    a_three_prime_index: int,
    b_three_prime_index: int | None,
    temp_c: float,
    conditions: SaltConditions,
    render: bool,
) -> Interaction:
    """Best ungapped antiparallel alignment of *a* against *b_anti*.

    *b_anti* is the partner written 3'->5' so that index *i* of *a* opposes
    index *i* of *b_anti* at zero offset.
    """
    best = EMPTY
    la, lb = len(a), len(b_anti)
    for offset in range(-(lb - 1), la):
        total = 0.0
        run: list[tuple[str, str]] = []
        best_run: list[tuple[str, str]] = []
        anchored = False
        for i in range(la):
            j = i - offset
            if j < 0 or j >= lb:
                if run and len(run) > len(best_run):
                    best_run = list(run)
                run = []
                continue
            if _pairs(a[i], b_anti[j]):
                total += _WEIGHT[a[i]]
                run.append((a[i], b_anti[j]))
                if i == a_three_prime_index or (
                    b_three_prime_index is not None and j == b_three_prime_index
                ):
                    anchored = True
            else:
                if len(run) > len(best_run):
                    best_run = list(run)
                run = []
        if len(run) > len(best_run):
            best_run = list(run)
        if not best_run:
            continue
        dg = _stack_dg(best_run, temp_c, conditions)
        # Rank primarily on stack free energy, then on raw complementarity.
        candidate_key = (dg, -total)
        incumbent_key = (best.dg, -best.score)
        if best.is_empty or candidate_key < incumbent_key:
            alignment = ""
            if render:
                top_off = max(0, -offset)
                bottom_off = max(0, offset)
                alignment = _render(a, b_anti, top_off, bottom_off)
            best = Interaction(
                score=total,
                dg=dg,
                offset=offset,
                stack_length=len(best_run),
                three_prime_anchored=anchored,
                alignment=alignment,
            )
    return best


def hetero_dimer(
    a: str,
    b: str,
    *,
    temp_c: float = 37.0,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    render: bool = False,
) -> Interaction:
    """Worst-case duplex between two different primers."""
    a = a.upper()
    b = b.upper()
    b_anti = b[::-1]
    return _align(
        a,
        b_anti,
        a_three_prime_index=len(a) - 1,
        b_three_prime_index=0,  # b's 3' end sits at index 0 of the reversed string
        temp_c=temp_c,
        conditions=conditions,
        render=render,
    )


def self_dimer(
    seq: str,
    *,
    temp_c: float = 37.0,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    render: bool = False,
) -> Interaction:
    """Worst-case duplex between two copies of the same primer."""
    return hetero_dimer(seq, seq, temp_c=temp_c, conditions=conditions, render=render)


def hairpin(
    seq: str,
    *,
    temp_c: float = 37.0,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    min_loop: int = MIN_LOOP,
    render: bool = False,
) -> Interaction:
    """Most stable intramolecular stem-loop.

    Searches every (i, j) seed with i < j, extends the stem inwards while bases
    remain complementary and the enclosed loop stays at or above *min_loop*.
    """
    s = seq.upper()
    n = len(s)
    best = EMPTY
    for i in range(n - min_loop - 1):
        for j in range(n - 1, i + min_loop, -1):
            if not _pairs(s[i], s[j]):
                continue
            pairs: list[tuple[str, str]] = []
            k = 0
            while (
                i + k < j - k
                and (j - k) - (i + k) - 1 >= min_loop
                and _pairs(s[i + k], s[j - k])
            ):
                pairs.append((s[i + k], s[j - k]))
                k += 1
            if len(pairs) < 2:
                continue
            dg = _stack_dg(pairs, temp_c, conditions)
            score = sum(_WEIGHT[p[0]] for p in pairs)
            anchored = (j + len(pairs) - 1) >= n - 1 or j == n - 1
            if best.is_empty or (dg, -score) < (best.dg, -best.score):
                alignment = ""
                if render:
                    stem_top = "".join(p[0] for p in pairs)
                    stem_bottom = "".join(p[1] for p in pairs)
                    loop = s[i + len(pairs) : j - len(pairs) + 1]
                    alignment = (
                        f"5'-{s[:i]}{stem_top} \n"
                        f"   {' ' * i}{'|' * len(pairs)}  loop={loop}\n"
                        f"3'-{s[j + 1:][::-1]}{stem_bottom} "
                    )
                best = Interaction(
                    score=score,
                    dg=dg,
                    offset=i,
                    stack_length=len(pairs),
                    three_prime_anchored=anchored,
                    alignment=alignment,
                )
    return best


def three_prime_dimer_score(a: str, b: str, *, window: int = 5) -> float:
    """Complementarity between the 3' termini of two primers.

    Only the terminal *window* bases of each primer are considered, in the
    register where the two 3' ends face each other. This catches the primer-dimer
    geometry that actually produces artefact product.
    """
    a_tail = a.upper()[-window:]
    b_tail = b.upper()[-window:]
    best = 0.0
    for offset in range(-(len(b_tail) - 1), len(a_tail)):
        total = 0.0
        for i, base in enumerate(a_tail):
            j = i - offset
            if 0 <= j < len(b_tail):
                partner = b_tail[::-1][j]
                if _pairs(base, partner):
                    total += _WEIGHT[base]
        best = max(best, total)
    return best
