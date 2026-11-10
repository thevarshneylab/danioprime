"""Assembly and scoring of primer pairs."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

from ..config import DesignParams
from ..thermo.structure import hetero_dimer, three_prime_dimer_score
from .candidate import Primer, TemplateContext
from .regions import Interval, IntervalSet


@dataclass(slots=True)
class PrimerPair:
    """A forward/reverse pair and its predicted product."""

    forward: Primer
    reverse: Primer
    product_size: int
    tm_difference: float = 0.0
    hetero_dimer_dg: float = 0.0
    three_prime_dimer_score: float = 0.0
    spans_targets: bool = True
    paralog_amplicons: list[dict[str, Any]] = field(default_factory=list)
    offtarget_amplicons: list[dict[str, Any]] = field(default_factory=list)
    penalty: float = 0.0
    penalty_terms: dict[str, float] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    name: str = ""

    @property
    def amplicon_interval(self) -> Interval:
        return Interval(self.forward.start, self.reverse.end)

    def amplicon_seq(self, template: str) -> str:
        return template[self.forward.start : self.reverse.end]

    @property
    def is_specific(self) -> bool:
        """True when no paralog or off-target amplicon is predicted."""
        return not self.paralog_amplicons and not self.offtarget_amplicons

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "product_size": self.product_size,
            "amplicon_1based": self.amplicon_interval.as_one_based(),
            "tm_difference": round(self.tm_difference, 2),
            "hetero_dimer_dg": round(self.hetero_dimer_dg, 2),
            "three_prime_dimer_score": self.three_prime_dimer_score,
            "spans_targets": self.spans_targets,
            "specific": self.is_specific,
            "paralog_amplicons": self.paralog_amplicons,
            "offtarget_amplicons": self.offtarget_amplicons,
            "penalty": round(self.penalty, 3),
            "penalty_terms": {k: round(v, 3) for k, v in self.penalty_terms.items()},
            "notes": list(self.notes),
            "forward": self.forward.to_dict(),
            "reverse": self.reverse.to_dict(),
        }


def pair_is_geometrically_valid(
    fwd: Primer, rev: Primer, params: DesignParams
) -> tuple[bool, str]:
    """Check orientation, overlap and product size."""
    if fwd.orientation != "F" or rev.orientation != "R":
        return False, "orientation"
    if rev.end <= fwd.start:
        return False, "inverted"
    if fwd.end > rev.start:
        return False, "primers overlap"
    size = rev.end - fwd.start
    if not (params.min_product_size <= size <= params.max_product_size):
        return False, "product size"
    return True, ""


def score_pair(pair: PrimerPair, ctx: TemplateContext) -> PrimerPair:
    """Compute the pair-level penalty on top of the two primer penalties."""
    p = ctx.params
    w = p.weights
    terms: dict[str, float] = {
        "forward_primer": pair.forward.penalty,
        "reverse_primer": pair.reverse.penalty,
    }

    terms["tm_difference"] = w.tm_pair_difference * pair.tm_difference
    size_dev = abs(pair.product_size - p.opt_product_size)
    terms["product_size"] = w.product_size_deviation * size_dev

    excess = max(0.0, p.max_hetero_dimer_dg - pair.hetero_dimer_dg)
    if excess:
        terms["hetero_dimer"] = w.hetero_dimer_dg * excess
    tp_excess = max(0.0, pair.three_prime_dimer_score - p.max_three_prime_dimer_score)
    if tp_excess:
        terms["three_prime_dimer"] = w.three_prime_dimer * tp_excess

    if pair.paralog_amplicons:
        terms["paralog_amplicon"] = w.paralog_amplicon * len(pair.paralog_amplicons)
    if pair.offtarget_amplicons:
        terms["offtarget_amplicon"] = w.offtarget_amplicon * len(pair.offtarget_amplicons)

    pair.penalty_terms = {k: v for k, v in terms.items() if v}
    pair.penalty = sum(pair.penalty_terms.values())
    return pair


def build_pairs(
    forwards: Iterable[Primer],
    reverses: Iterable[Primer],
    ctx: TemplateContext,
    *,
    targets: IntervalSet | None = None,
    max_pairs: int | None = None,
) -> list[PrimerPair]:
    """Form, filter, score and diversify primer pairs.

    Candidates are considered in order of combined single-primer penalty, so
    the cheap thermodynamic work is done on the most promising combinations
    first. Pairs whose amplicon fails to span every target interval are
    discarded, which is what keeps a genotyping pair from sitting beside the
    lesion instead of across it.
    """
    p = ctx.params
    max_pairs = p.max_pairs if max_pairs is None else max_pairs
    fwd_list = sorted(forwards, key=lambda pr: pr.penalty)
    rev_list = sorted(reverses, key=lambda pr: pr.penalty)

    scored: list[PrimerPair] = []
    for fwd in fwd_list:
        for rev in rev_list:
            ok, _reason = pair_is_geometrically_valid(fwd, rev, p)
            if not ok:
                continue
            tm_diff = abs(fwd.tm - rev.tm)
            if tm_diff > p.max_pair_tm_difference:
                continue
            amplicon = Interval(fwd.start, rev.end)
            if targets and not all(amplicon.contains_interval(t) for t in targets):
                continue
            hd = hetero_dimer(
                fwd.seq, rev.seq, temp_c=p.anneal_temp_c, conditions=p.conditions
            )
            tp = three_prime_dimer_score(fwd.seq, rev.seq)
            if hd.dg < p.max_hetero_dimer_dg - 4.0:
                continue
            pair = PrimerPair(
                forward=fwd,
                reverse=rev,
                product_size=rev.end - fwd.start,
                tm_difference=tm_diff,
                hetero_dimer_dg=hd.dg,
                three_prime_dimer_score=tp,
                spans_targets=True,
            )
            score_pair(pair, ctx)
            scored.append(pair)

    scored.sort(key=lambda pr: pr.penalty)
    return diversify(scored, window=p.pair_diversity_window, limit=max_pairs)


def diversify(pairs: list[PrimerPair], *, window: int, limit: int) -> list[PrimerPair]:
    """Keep the best pairs whose 3' ends are at least *window* apart.

    Without this, the top 20 pairs are typically the same primer shifted by one
    base, which gives the user no real alternatives if the first pair fails at
    the bench.
    """
    chosen: list[PrimerPair] = []
    for pair in pairs:
        clash = False
        for kept in chosen:
            if (
                abs(pair.forward.three_prime_pos - kept.forward.three_prime_pos) < window
                and abs(pair.reverse.three_prime_pos - kept.reverse.three_prime_pos) < window
            ):
                clash = True
                break
        if clash:
            continue
        chosen.append(pair)
        if len(chosen) >= limit:
            break
    return chosen


def name_pairs(pairs: list[PrimerPair], prefix: str) -> list[PrimerPair]:
    """Assign stable, orderable names like ``kcnh6a_p1``."""
    for i, pair in enumerate(pairs, 1):
        pair.name = f"{prefix}_p{i}"
    return pairs
