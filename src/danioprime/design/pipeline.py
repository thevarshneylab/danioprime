"""The design driver shared by every assay mode.

One place that knows the order of operations:

1. enumerate candidates on each strand within the allowed search windows;
2. screen each surviving candidate against the "must not amplify" panel, so
   paralog-priming primers are demoted before pairs are built;
3. build, score and diversify pairs;
4. predict panel and off-target amplicons for the pairs that remain;
5. re-score and re-rank.

Screening primers before pairing and only predicting amplicons afterwards keeps
the expensive specificity work proportional to the number of pairs reported
rather than to the square of the candidate count.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from ..config import DesignParams
from ..errors import DesignError
from ..sequence import SeqRecord
from ..specificity.paralog import (
    ParalogPanel,
    assess_pair,
    assess_primer,
    summarise_panel_amplicons,
)
from ..specificity.scan import PrimingRules
from ..thermo.backend import get_backend
from .candidate import CandidateSet, Primer, TemplateContext, enumerate_candidates
from .pairing import PrimerPair, build_pairs, diversify, name_pairs, score_pair
from .regions import Interval, IntervalSet
from .scoring import score_primer


@dataclass(slots=True)
class DesignRequest:
    """Everything needed for one design run."""

    template: SeqRecord
    params: DesignParams
    #: Windows in which forward primers may sit; empty means the whole template.
    forward_search: IntervalSet = field(default_factory=IntervalSet)
    reverse_search: IntervalSet = field(default_factory=IntervalSet)
    #: The amplicon must span every one of these intervals.
    targets: IntervalSet = field(default_factory=IntervalSet)
    excluded: IntervalSet = field(default_factory=IntervalSet)
    variants: IntervalSet = field(default_factory=IntervalSet)
    softmask: IntervalSet = field(default_factory=IntervalSet)
    junctions: tuple[int, ...] = ()
    panel: ParalogPanel | None = None
    #: Require that a primer straddles an exon-exon junction: ``"forward"``,
    #: ``"reverse"``, or ``"either"`` (at least one of the two). Enforced as a
    #: hard constraint on the candidate pools rather than via the scoring
    #: bonus, because a bonus can still be outranked and silently dropped.
    require_spanning: str | None = None
    name: str = "design"
    #: Free-form notes carried into the report (assay geometry, cut site, ...).
    context: dict[str, Any] = field(default_factory=dict)

    def priming_rules(self) -> PrimingRules:
        p = self.params
        return PrimingRules(
            clean_three_prime=min(3, p.priming_three_prime_window),
            three_prime_window=p.priming_three_prime_window,
            max_mismatches=p.priming_max_mismatches,
        )


@dataclass(slots=True)
class DesignResult:
    """Output of a design run, ready for the report writers."""

    request: DesignRequest
    pairs: list[PrimerPair] = field(default_factory=list)
    forward_candidates: CandidateSet = field(default_factory=CandidateSet)
    reverse_candidates: CandidateSet = field(default_factory=CandidateSet)
    discrimination: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    stats: dict[str, Any] = field(default_factory=dict)

    @property
    def template(self) -> SeqRecord:
        return self.request.template

    @property
    def best(self) -> PrimerPair | None:
        return self.pairs[0] if self.pairs else None

    def diagnosis(self) -> str:
        """Why no pair was found, phrased so the user knows what to relax."""
        if self.pairs:
            return ""
        f, r = self.forward_candidates, self.reverse_candidates
        if not f and not r:
            return (
                "no candidate primer passed the single-primer filters on either strand. "
                f"Filters that removed candidates -- forward: {f.reason_summary()}; "
                f"reverse: {r.reason_summary()}. "
                "Widening --min-tm/--max-tm or relaxing --require-gc-clamp usually helps."
            )
        if not f:
            return f"no usable forward primer (filtered by: {f.reason_summary()})"
        if not r:
            return f"no usable reverse primer (filtered by: {r.reason_summary()})"
        return (
            f"{len(f)} forward and {len(r)} reverse candidates passed, but no pair satisfied "
            "the product-size window, pair Tm difference, target-spanning requirement or "
            "specificity filter. Try widening --max-product-size or --max-pair-tm-difference."
        )


def _screen_primer_against_panel(
    primer: Primer, ctx: TemplateContext, request: DesignRequest
) -> dict[str, Any] | None:
    """Annotate a primer with panel cross-reactivity and re-score it."""
    panel = request.panel
    if not panel:
        return None
    report = assess_primer(
        primer.seq,
        panel,
        primer_name=f"{primer.orientation}@{primer.start + 1}",
        rules=request.priming_rules(),
        conditions=ctx.params.conditions,
    )
    primer.paralog_sites = len(report.priming_sites)
    if report.primes_on_panel:
        primer.notes.append(report.verdict())
    elif report.best_site is not None and report.discriminates_at_3p:
        offsets = report.best_site.mismatch_offsets_from_3p
        if offsets:
            primer.notes.append(
                f"discriminates {report.best_site.subject} by a mismatch "
                f"{min(offsets) + 1} nt from the 3' end"
            )
    score_primer(primer, ctx)
    return report.to_dict()


def _spanning(primers: Sequence[Primer]) -> list[Primer]:
    return [p for p in primers if p.spans_junction is not None]


def _compatible_partners(
    anchor: Primer, partners: Sequence[Primer], params: DesignParams, limit: int
) -> list[Primer]:
    """The best *limit* partners that could actually pair with *anchor*.

    Selecting partners by absolute penalty alone biases the pool towards
    Tm == opt_tm, which then cannot satisfy ``max_pair_tm_difference`` against
    an anchor whose own Tm is off-optimum -- the exact situation a forced
    junction-spanning primer creates. Filtering on compatibility with the
    anchor first fixes that without enlarging the search.
    """
    window = params.max_pair_tm_difference
    out: list[Primer] = []
    for candidate in partners:
        if abs(candidate.tm - anchor.tm) > window:
            continue
        out.append(candidate)
        if len(out) >= limit:
            break
    return out


def _build_pair_pool(
    fwd: CandidateSet,
    rev: CandidateSet,
    ctx: TemplateContext,
    request: DesignRequest,
    screen_top: int,
) -> list[PrimerPair]:
    """Form pairs, honouring any junction-spanning requirement.

    With ``require_spanning="either"`` both assignments are enumerated --
    spanning forward against a compatible reverse, and the mirror image -- so a
    spanning candidate is never lost to the pool cutoff just because its
    thermodynamics rank it below the others.
    """
    params = ctx.params
    limit = max(params.max_pairs * 3, params.max_pairs)
    need = request.require_spanning
    merged: list[PrimerPair] = []
    seen: set[tuple[int, int, int, int]] = set()

    def collect(f_pool: Sequence[Primer], r_pool: Sequence[Primer]) -> None:
        if not f_pool or not r_pool:
            return
        for pair in build_pairs(
            f_pool, r_pool, ctx, targets=request.targets or None, max_pairs=limit
        ):
            key = (
                pair.forward.start,
                pair.forward.end,
                pair.reverse.start,
                pair.reverse.end,
            )
            if key in seen:
                continue
            seen.add(key)
            merged.append(pair)

    anchored = False
    if need in {"forward", "either"}:
        for anchor in _spanning(fwd.primers)[:screen_top]:
            anchored = True
            collect([anchor], _compatible_partners(anchor, rev.primers, params, screen_top))
    if need in {"reverse", "either"}:
        for anchor in _spanning(rev.primers)[:screen_top]:
            anchored = True
            collect(_compatible_partners(anchor, fwd.primers, params, screen_top), [anchor])

    if not anchored:
        collect(list(fwd.primers[:screen_top]), list(rev.primers[:screen_top]))

    merged.sort(key=lambda pr: pr.penalty)
    return diversify(merged, window=params.pair_diversity_window, limit=limit)


def _windows_implied_by_targets(
    request: DesignRequest, length: int
) -> tuple[IntervalSet | None, IntervalSet | None]:
    """Search windows that the target and product-size constraints already imply.

    A valid pair spans every target and is at most ``max_product_size`` long, so
    a forward primer further upstream of the first target than that -- or a
    reverse primer further downstream of the last -- cannot appear in any pair
    that would be accepted. Narrowing the enumeration to those windows
    therefore removes no possible result.

    It is not an optimisation for its own sake. Candidate enumeration is capped
    (``max_candidates_per_side``), and on hitting the cap it keeps the
    candidates it found first, which are the ones nearest the 5' end. On a
    multi-kilobase template with a target in the middle, an unaimed enumeration
    spends its whole budget at the wrong end and then reports that no pair
    satisfied the constraints -- which is true, and useless. This is the
    replacement for the flanking windows that the old cut-site geometry
    supplied as a side effect of anchoring on a coordinate.
    """
    if not request.targets:
        return None, None
    starts = [interval.start for interval in request.targets]
    ends = [interval.end for interval in request.targets]
    span = max(ends) - min(starts)
    reach = max(0, request.params.max_product_size - span)
    forward = (
        request.forward_search
        or IntervalSet([Interval(max(0, min(starts) - reach), min(starts))])
    )
    reverse = (
        request.reverse_search
        or IntervalSet([Interval(max(ends), min(length, max(ends) + reach))])
    )
    return forward, reverse


def run_design(request: DesignRequest, *, screen_top: int = 60) -> DesignResult:
    """Execute a design request."""
    params = request.params
    params.validate()
    template = request.template
    if len(template.seq) < params.min_product_size:
        raise DesignError(
            f"template {template.name} is {len(template.seq)} nt, shorter than the minimum "
            f"product size ({params.min_product_size} nt)"
        )

    backend = get_backend(params.thermo_backend, params.conditions, params.anneal_temp_c)
    ctx = TemplateContext(
        seq=template.seq,
        params=params,
        backend=backend,
        name=template.name,
        excluded=request.excluded,
        variants=request.variants,
        softmask=request.softmask,
        junctions=request.junctions,
    )

    if request.targets:
        # A pair that spans the targets needs room for the span plus both
        # primers. When that exceeds the product-size window, no pair can
        # satisfy both constraints -- so say which two constraints conflict
        # rather than enumerating thousands of candidates and reporting that
        # none of them paired.
        starts = [interval.start for interval in request.targets]
        ends = [interval.end for interval in request.targets]
        span = max(ends) - min(starts)
        needed = span + 2 * params.min_length
        if needed > params.max_product_size:
            raise DesignError(
                f"the region to interrogate spans {span} nt, so the shortest product that "
                f"could span it with a primer at each end is {needed} nt -- longer than "
                f"--max-product-size ({params.max_product_size} nt). Either narrow the "
                f"target to the part that has to be covered, or raise the product-size "
                f"window; note that a readout chosen for small products, such as "
                f"fragment analysis, may not be able to accommodate a region this wide."
            )

    f_window, r_window = _windows_implied_by_targets(request, len(template.seq))
    fwd = enumerate_candidates(ctx, "F", f_window or request.forward_search or None)
    rev = enumerate_candidates(ctx, "R", r_window or request.reverse_search or None)

    result = DesignResult(
        request=request,
        forward_candidates=fwd,
        reverse_candidates=rev,
    )
    result.stats = {
        "template_length": len(template.seq),
        "forward_candidates": len(fwd),
        "reverse_candidates": len(rev),
        "forward_rejections": dict(fwd.rejections),
        "reverse_rejections": dict(rev.rejections),
        "masked_variant_positions": request.variants.total_length,
        "excluded_bases": request.excluded.total_length,
    }
    if fwd.truncated or rev.truncated:
        result.warnings.append(
            "candidate enumeration hit --max-candidates-per-side; "
            "results are the first candidates found, not an exhaustive search"
        )

    # Screen only the most promising candidates against the panel: the panel
    # scan is the expensive step and the rest will not reach the output anyway.
    if request.panel:
        for primer in list(fwd.primers)[:screen_top] + list(rev.primers)[:screen_top]:
            _screen_primer_against_panel(primer, ctx, request)
        if params.reject_paralog_priming:
            for side, candidates in (("forward", fwd), ("reverse", rev)):
                dropped = sum(1 for p in candidates.primers[:screen_top] if p.paralog_sites)
                candidates.primers = [
                    p for p in candidates.primers[:screen_top] if not p.paralog_sites
                ]
                if dropped:
                    candidates.rejections["primes-on-avoid-panel"] = dropped
                if not candidates.primers:
                    result_note = (
                        f"every screened {side} candidate can prime on the avoid-panel; "
                        "drop --strict-specificity to allow primers that bind the panel "
                        "but cannot produce a product with their partner"
                    )
                    request.context.setdefault("notes", []).append(result_note)
        fwd.primers.sort(key=lambda pr: pr.penalty)
        rev.primers.sort(key=lambda pr: pr.penalty)

    pairs = _build_pair_pool(fwd, rev, ctx, request, screen_top)

    if request.panel and pairs:
        kept: list[PrimerPair] = []
        for pair in pairs:
            amplicons, reports = assess_pair(
                pair.forward.seq,
                pair.reverse.seq,
                request.panel,
                min_size=params.min_product_size,
                max_size=max(params.max_product_size, 3000),
                forward_name="F",
                reverse_name="R",
                rules=request.priming_rules(),
                conditions=params.conditions,
            )
            pair.paralog_amplicons = summarise_panel_amplicons(amplicons)
            if pair.paralog_amplicons:
                names = sorted({a["subject"] for a in pair.paralog_amplicons})
                pair.notes.append(
                    "predicted to co-amplify " + ", ".join(names)
                )
                if params.reject_paralog_amplicon:
                    continue
            else:
                pair.notes.append("no product predicted on the avoid-panel")
            result.discrimination.setdefault(
                "pairs", {}
            )[f"{pair.forward.start + 1}-{pair.reverse.end}"] = {
                k: v.to_dict() for k, v in reports.items()
            }
            score_pair(pair, ctx)
            kept.append(pair)
        pairs = kept
        pairs.sort(key=lambda pr: pr.penalty)

    result.pairs = name_pairs(pairs[: params.max_pairs], request.name)

    if not result.pairs:
        result.warnings.append(result.diagnosis())
    return result


def windows_around(
    point: int,
    template_length: int,
    *,
    inner: int,
    outer: int,
) -> tuple[IntervalSet, IntervalSet]:
    """Search windows flanking *point*.

    Returns ``(forward_windows, reverse_windows)``: the forward window lies
    ``inner..outer`` bases 5' of *point*, the reverse window the same distance
    3' of it. Used by every assay that designs across a feature.
    """
    if inner < 0 or outer <= inner:
        raise DesignError(f"need 0 <= inner < outer, got inner={inner}, outer={outer}")
    f_start = max(0, point - outer)
    f_end = max(0, point - inner)
    r_start = min(template_length, point + inner)
    r_end = min(template_length, point + outer)
    forward = IntervalSet([Interval(f_start, f_end)] if f_end > f_start else [])
    reverse = IntervalSet([Interval(r_start, r_end)] if r_end > r_start else [])
    return forward, reverse


def merge_search(primary: IntervalSet, fallback: Sequence[Interval]) -> IntervalSet:
    """Use *primary* when non-empty, else *fallback*."""
    return primary if primary else IntervalSet(fallback)
