"""Additive penalty model.

Each criterion contributes a non-negative penalty (except the deliberate
junction-spanning bonus); the total is the sort key. Keeping the terms in a
dict rather than collapsing them to a single float means the report can say
*why* a primer ranked where it did, which matters when a user wants to relax
one constraint rather than guess.

Structure penalties are applied only to the excess beyond the configured
threshold, so a primer that is comfortably within limits is not penalised for
having any structure at all.
"""

from __future__ import annotations

from .candidate import Primer, TemplateContext


def _excess_below(value: float, limit: float) -> float:
    """How far *value* falls below *limit* (both are negative dG values)."""
    return max(0.0, limit - value)


def _excess_above(value: float, limit: float) -> float:
    return max(0.0, value - limit)


def score_primer(primer: Primer, ctx: TemplateContext) -> Primer:
    """Compute and attach the penalty breakdown for a single primer."""
    p = ctx.params
    w = p.weights
    terms: dict[str, float] = {}

    terms["tm"] = w.tm_deviation * abs(primer.tm - p.opt_tm)
    terms["gc"] = w.gc_deviation * abs(primer.gc - p.opt_gc)
    terms["length"] = w.length_deviation * abs(primer.length - p.opt_length)

    hairpin = _excess_below(primer.hairpin_dg, p.max_hairpin_dg)
    if hairpin:
        terms["hairpin"] = w.hairpin_dg * hairpin
    selfd = _excess_below(primer.self_dimer_dg, p.max_self_dimer_dg)
    if selfd:
        terms["self_dimer"] = w.self_dimer_dg * selfd
    endstab = _excess_below(primer.end_stability, p.max_end_stability)
    if endstab:
        terms["end_stability"] = w.end_stability * endstab

    if primer.gc_clamp == 0 and p.require_gc_clamp == 0:
        # Not disqualifying when the clamp is not required, but still worse.
        terms["gc_clamp"] = w.gc_clamp_missing

    polyx = _excess_above(primer.poly_x, max(2, p.max_poly_x - 1))
    if polyx:
        terms["poly_x"] = w.poly_x * polyx
    dinuc = _excess_above(primer.dinucleotide_repeat, max(2, p.max_dinucleotide_repeat - 1))
    if dinuc:
        terms["dinucleotide_repeat"] = w.dinucleotide_repeat * dinuc

    if primer.ambiguous:
        terms["ambiguous"] = w.ambiguous_base * primer.ambiguous
    if primer.softmask_fraction:
        terms["softmask"] = w.softmask_overlap * primer.softmask_fraction
    if primer.variants_overlapping:
        terms["variants"] = w.variant_overlap * primer.variants_overlapping
    if primer.variants_three_prime:
        terms["variants_3p"] = w.variant_three_prime * primer.variants_three_prime
    if primer.offtarget_sites:
        terms["offtarget_sites"] = w.offtarget_site * primer.offtarget_sites
    if primer.paralog_sites:
        terms["paralog_sites"] = w.paralog_priming_site * primer.paralog_sites
    if primer.spans_junction is not None:
        terms["junction_bonus"] = w.junction_spanning_bonus

    primer.penalty_terms = {k: v for k, v in terms.items() if v}
    primer.penalty = sum(primer.penalty_terms.values())
    return primer


def rescore_primer(primer: Primer, ctx: TemplateContext) -> Primer:
    """Recompute the penalty after specificity results have been attached."""
    return score_primer(primer, ctx)
