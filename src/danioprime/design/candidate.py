"""Primer candidate representation, property evaluation and enumeration."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any, Literal

from ..config import DesignParams
from ..errors import DesignError
from ..sequence import (
    DNA_BASES,
    gc_percent,
    max_dinucleotide_repeat,
    max_homopolymer,
    revcomp,
)
from ..thermo.backend import ThermoBackend
from .regions import Interval, IntervalSet

Orientation = Literal["F", "R"]


@dataclass(slots=True)
class Primer:
    """A candidate primer placed on a template.

    Coordinates are always template coordinates, half-open and zero-based,
    regardless of orientation. ``seq`` is the oligo as it would be ordered
    (5'->3'), so for a reverse primer it is the reverse complement of the
    template footprint.
    """

    seq: str
    start: int
    end: int
    orientation: Orientation
    tm: float = 0.0
    gc: float = 0.0
    hairpin_dg: float = 0.0
    self_dimer_dg: float = 0.0
    end_stability: float = 0.0
    gc_clamp: int = 0
    poly_x: int = 0
    dinucleotide_repeat: int = 0
    ambiguous: int = 0
    softmask_fraction: float = 0.0
    variants_overlapping: int = 0
    variants_three_prime: int = 0
    variant_positions: list[int] = field(default_factory=list)
    spans_junction: int | None = None
    offtarget_sites: int = 0
    paralog_sites: int = 0
    tail: str = ""
    tail_name: str = ""
    penalty: float = 0.0
    penalty_terms: dict[str, float] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def footprint(self) -> Interval:
        """Template interval covered by the primer."""
        return Interval(self.start, self.end)

    @property
    def three_prime_pos(self) -> int:
        """Template position of the 3'-terminal base."""
        return self.end - 1 if self.orientation == "F" else self.start

    @property
    def five_prime_pos(self) -> int:
        """Template position of the 5'-terminal base."""
        return self.start if self.orientation == "F" else self.end - 1

    @property
    def ordered_seq(self) -> str:
        """Full oligo to order, including any 5' tail."""
        return self.tail + self.seq

    def three_prime_window(self, n: int) -> Interval:
        """Template interval covered by the *n* 3'-terminal bases."""
        n = max(1, min(n, self.length))
        if self.orientation == "F":
            return Interval(self.end - n, self.end)
        return Interval(self.start, self.start + n)

    def to_dict(self) -> dict[str, Any]:
        return {
            "seq": self.seq,
            "ordered_seq": self.ordered_seq,
            "orientation": self.orientation,
            "start": self.start,
            "end": self.end,
            "position_1based": self.footprint.as_one_based(),
            "length": self.length,
            "tm": round(self.tm, 2),
            "gc_percent": round(self.gc, 1),
            "hairpin_dg": round(self.hairpin_dg, 2),
            "self_dimer_dg": round(self.self_dimer_dg, 2),
            "end_stability_dg": round(self.end_stability, 2),
            "gc_clamp": self.gc_clamp,
            "poly_x": self.poly_x,
            "dinucleotide_repeat": self.dinucleotide_repeat,
            "softmask_fraction": round(self.softmask_fraction, 3),
            "variants_overlapping": self.variants_overlapping,
            "variants_three_prime": self.variants_three_prime,
            "variant_positions_1based": [p + 1 for p in self.variant_positions],
            "spans_junction": self.spans_junction,
            "offtarget_sites": self.offtarget_sites,
            "paralog_sites": self.paralog_sites,
            "tail_name": self.tail_name,
            "penalty": round(self.penalty, 3),
            "penalty_terms": {k: round(v, 3) for k, v in self.penalty_terms.items()},
            "notes": list(self.notes),
        }


@dataclass(slots=True)
class TemplateContext:
    """Everything a primer is evaluated against."""

    seq: str
    params: DesignParams
    backend: ThermoBackend
    name: str = "template"
    excluded: IntervalSet = field(default_factory=IntervalSet)
    variants: IntervalSet = field(default_factory=IntervalSet)
    softmask: IntervalSet = field(default_factory=IntervalSet)
    #: Template offsets of exon-exon junctions (the position of the first base
    #: of the downstream exon in a spliced transcript).
    junctions: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if not self.seq:
            raise DesignError("template sequence is empty")

    def __len__(self) -> int:
        return len(self.seq)


def count_gc_clamp(seq: str) -> int:
    """Number of consecutive G/C bases at the 3' end."""
    n = 0
    for base in reversed(seq.upper()):
        if base in "GC":
            n += 1
        else:
            break
    return n


def evaluate_composition(primer: Primer, ctx: TemplateContext) -> Primer:
    """Fill in the cheap, sequence-only properties.

    Separated from the thermodynamic stage so that enumeration can discard
    obviously unusable placements before paying for a structure search.
    """
    seq = primer.seq
    params = ctx.params

    primer.ambiguous = sum(1 for b in seq if b not in DNA_BASES)
    primer.gc = gc_percent(seq)
    primer.gc_clamp = count_gc_clamp(seq)
    primer.poly_x = max_homopolymer(seq)
    primer.dinucleotide_repeat = max_dinucleotide_repeat(seq)

    fp = primer.footprint
    if ctx.softmask:
        primer.softmask_fraction = ctx.softmask.overlap_length(fp) / primer.length
    if ctx.variants:
        hits = ctx.variants.positions_in(fp)
        primer.variant_positions = hits
        primer.variants_overlapping = len(hits)
        tp_window = primer.three_prime_window(params.variant_three_prime_window)
        primer.variants_three_prime = sum(1 for pos in hits if pos in tp_window)

    if ctx.junctions:
        for j in ctx.junctions:
            # A primer only resists genomic-DNA amplification if it has a real
            # foothold on both exons; fewer than 4 nt on one side will still
            # anneal to the unspliced locus.
            if fp.start + 4 <= j <= fp.end - 4:
                primer.spans_junction = j
                break
    return primer


def evaluate_thermo(primer: Primer, ctx: TemplateContext, *, structure: bool = True) -> Primer:
    """Fill in Tm and, unless *structure* is False, the structure properties."""
    backend = ctx.backend
    if primer.ambiguous and not ctx.params.allow_ambiguous:
        primer.notes.append("contains ambiguity codes")
        primer.tm = float("nan")
        return primer
    # Ambiguity codes are resolved to A for the thermodynamic estimate; the
    # penalty model flags their presence separately.
    probe = (
        primer.seq
        if not primer.ambiguous
        else "".join(b if b in DNA_BASES else "A" for b in primer.seq)
    )
    primer.tm = backend.tm(probe)
    if structure:
        primer.hairpin_dg = backend.hairpin(probe).dg
        primer.self_dimer_dg = backend.self_dimer(probe).dg
        primer.end_stability = backend.end_stability(probe)
    return primer


def evaluate_primer(primer: Primer, ctx: TemplateContext) -> Primer:
    """Fill in every property of *primer* in place (convenience wrapper)."""
    evaluate_composition(primer, ctx)
    return evaluate_thermo(primer, ctx)


def check_composition(primer: Primer, ctx: TemplateContext) -> tuple[bool, str]:
    """Disqualifying constraints that need no thermodynamics."""
    p = ctx.params
    if primer.ambiguous and not p.allow_ambiguous:
        return False, "ambiguous base"
    if not (p.min_length <= primer.length <= p.max_length):
        return False, "length"
    if primer.gc_clamp < p.require_gc_clamp:
        return False, "gc-clamp-missing"
    if primer.gc_clamp > p.max_gc_clamp:
        return False, "gc-clamp-excessive"
    if not (p.min_gc <= primer.gc <= p.max_gc):
        return False, "gc"
    if primer.poly_x > p.max_poly_x:
        return False, "poly-x"
    if primer.dinucleotide_repeat > p.max_dinucleotide_repeat:
        return False, "dinucleotide-repeat"
    if primer.softmask_fraction > p.max_softmask_fraction:
        return False, "repeat-masked"
    if ctx.excluded.overlaps(primer.footprint):
        return False, "excluded-region"
    if p.reject_variant_overlap and primer.variants_overlapping:
        return False, "variant-overlap"
    if primer.variants_three_prime:
        return False, "variant-at-3prime"
    return True, ""


def check_tm(primer: Primer, ctx: TemplateContext) -> tuple[bool, str]:
    """Melting-temperature window."""
    p = ctx.params
    if primer.tm != primer.tm:  # NaN
        return False, "tm-undefined"
    if not (p.min_tm <= primer.tm <= p.max_tm):
        return False, "tm"
    return True, ""


def check_structure(primer: Primer, ctx: TemplateContext) -> tuple[bool, str]:
    """Hairpin, self-dimer and 3'-end stability limits."""
    p = ctx.params
    if primer.hairpin_dg < p.max_hairpin_dg:
        return False, "hairpin"
    if primer.self_dimer_dg < p.max_self_dimer_dg:
        return False, "self-dimer"
    if primer.end_stability < p.max_end_stability:
        return False, "end-stability"
    return True, ""


def passes_hard_filters(primer: Primer, ctx: TemplateContext) -> tuple[bool, str]:
    """Apply every disqualifying constraint. Returns ``(ok, reason)``."""
    for check in (check_composition, check_tm, check_structure):
        ok, reason = check(primer, ctx)
        if not ok:
            return False, reason
    return True, ""


def iter_placements(
    ctx: TemplateContext,
    orientation: Orientation,
    search: IntervalSet | None = None,
) -> Iterator[tuple[int, int]]:
    """Yield ``(start, end)`` footprints to consider.

    When *search* is given, only footprints lying wholly inside one of its
    intervals are produced. Restricting the search space is the main lever on
    runtime: a 50 kb locus unrestricted yields ~500k candidates per strand,
    whereas the two 300 bp windows used for genotyping yield ~6k.
    """
    p = ctx.params
    n = len(ctx)
    windows = list(search) if search else [Interval(0, n)]
    seen: set[tuple[int, int]] = set()
    for window in windows:
        lo = max(0, window.start)
        hi = min(n, window.end)
        for start in range(lo, hi):
            for length in range(p.min_length, p.max_length + 1):
                end = start + length
                if end > hi:
                    break
                key = (start, end)
                if key in seen:
                    continue
                seen.add(key)
                yield key


@dataclass(slots=True)
class CandidateSet:
    """Candidates for one strand, plus a breakdown of what was filtered out."""

    primers: list[Primer] = field(default_factory=list)
    rejections: dict[str, int] = field(default_factory=dict)
    truncated: bool = False

    def __iter__(self) -> Iterator[Primer]:
        return iter(self.primers)

    def __len__(self) -> int:
        return len(self.primers)

    def __bool__(self) -> bool:
        return bool(self.primers)

    def reason_summary(self) -> str:
        """Human-readable summary of why candidates were discarded."""
        if not self.rejections:
            return "no candidates were filtered"
        ordered = sorted(self.rejections.items(), key=lambda kv: -kv[1])
        return ", ".join(f"{reason}: {count}" for reason, count in ordered)


def enumerate_candidates(
    ctx: TemplateContext,
    orientation: Orientation,
    search: IntervalSet | None = None,
    *,
    score: bool = True,
) -> CandidateSet:
    """Build the filtered, scored candidate list for one strand.

    Filters are applied in increasing order of cost -- composition, then Tm,
    then secondary structure -- so the expensive hairpin and dimer searches
    only ever run on placements that already satisfy everything cheaper.
    """
    from .scoring import score_primer  # local import avoids a circular dependency

    result = CandidateSet()
    for start, end in iter_placements(ctx, orientation, search):
        footprint_seq = ctx.seq[start:end]
        seq = footprint_seq if orientation == "F" else revcomp(footprint_seq)
        primer = Primer(seq=seq, start=start, end=end, orientation=orientation)

        evaluate_composition(primer, ctx)
        ok, reason = check_composition(primer, ctx)
        if not ok:
            result.rejections[reason] = result.rejections.get(reason, 0) + 1
            continue

        evaluate_thermo(primer, ctx, structure=False)
        ok, reason = check_tm(primer, ctx)
        if not ok:
            result.rejections[reason] = result.rejections.get(reason, 0) + 1
            continue

        probe = (
            primer.seq
            if not primer.ambiguous
            else "".join(b if b in DNA_BASES else "A" for b in primer.seq)
        )
        primer.hairpin_dg = ctx.backend.hairpin(probe).dg
        primer.self_dimer_dg = ctx.backend.self_dimer(probe).dg
        primer.end_stability = ctx.backend.end_stability(probe)
        ok, reason = check_structure(primer, ctx)
        if not ok:
            result.rejections[reason] = result.rejections.get(reason, 0) + 1
            continue

        if score:
            score_primer(primer, ctx)
        result.primers.append(primer)
        if len(result.primers) >= ctx.params.max_candidates_per_side:
            result.truncated = True
            break
    result.primers.sort(key=lambda pr: pr.penalty)
    return result
