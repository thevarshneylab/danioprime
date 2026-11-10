"""The penalty model: every term, and the behaviour that depends on them."""

from __future__ import annotations

import pytest

from danioprime.config import DesignParams, PenaltyWeights
from danioprime.design.candidate import Primer, TemplateContext
from danioprime.design.scoring import score_primer
from danioprime.thermo.backend import get_backend


@pytest.fixture
def ctx():
    params = DesignParams()
    return TemplateContext(
        seq="ACGT" * 200,
        params=params,
        backend=get_backend("builtin", params.conditions, params.anneal_temp_c),
    )


def primer_with(**kwargs) -> Primer:
    """A primer sitting exactly at every optimum, then nudged by *kwargs*."""
    base = {
        "seq": "A" * 20,
        "start": 0,
        "end": 20,
        "orientation": "F",
        "tm": 60.0,
        "gc": 50.0,
        "gc_clamp": 1,
        "poly_x": 2,
        "dinucleotide_repeat": 1,
        "hairpin_dg": 0.0,
        "self_dimer_dg": 0.0,
        "end_stability": 0.0,
    }
    base.update(kwargs)
    return Primer(**base)


class TestIndividualTerms:
    def test_a_primer_at_every_optimum_scores_zero(self, ctx):
        p = score_primer(primer_with(), ctx)
        assert p.penalty == pytest.approx(0.0)
        assert p.penalty_terms == {}

    def test_tm_deviation_is_symmetric_and_linear(self, ctx):
        hot = score_primer(primer_with(tm=63.0), ctx).penalty_terms["tm"]
        cold = score_primer(primer_with(tm=57.0), ctx).penalty_terms["tm"]
        assert hot == pytest.approx(cold)
        assert hot == pytest.approx(3.0 * ctx.params.weights.tm_deviation)

    def test_gc_and_length_deviation(self, ctx):
        p = score_primer(primer_with(gc=60.0, seq="A" * 24, end=24), ctx)
        assert p.penalty_terms["gc"] == pytest.approx(10.0 * ctx.params.weights.gc_deviation)
        assert p.penalty_terms["length"] == pytest.approx(
            4.0 * ctx.params.weights.length_deviation
        )

    def test_structure_penalties_apply_only_to_the_excess(self, ctx):
        limit = ctx.params.max_hairpin_dg  # -3.0
        inside = score_primer(primer_with(hairpin_dg=limit + 0.5), ctx)
        assert "hairpin" not in inside.penalty_terms
        outside = score_primer(primer_with(hairpin_dg=limit - 2.0), ctx)
        assert outside.penalty_terms["hairpin"] == pytest.approx(
            2.0 * ctx.params.weights.hairpin_dg
        )

    def test_self_dimer_and_end_stability_excess(self, ctx):
        p = score_primer(
            primer_with(
                self_dimer_dg=ctx.params.max_self_dimer_dg - 1.0,
                end_stability=ctx.params.max_end_stability - 2.0,
            ),
            ctx,
        )
        assert p.penalty_terms["self_dimer"] > 0
        assert p.penalty_terms["end_stability"] > 0

    def test_poly_x_and_repeats_are_penalised_before_the_hard_limit(self, ctx):
        # max_poly_x is 4, so a run of 4 is allowed but already costed.
        p = score_primer(primer_with(poly_x=4, dinucleotide_repeat=4), ctx)
        assert p.penalty_terms["poly_x"] > 0
        assert p.penalty_terms["dinucleotide_repeat"] > 0

    def test_variants_cost_more_at_the_three_prime_end(self, ctx):
        internal = score_primer(primer_with(variants_overlapping=1), ctx)
        terminal = score_primer(
            primer_with(variants_overlapping=1, variants_three_prime=1), ctx
        )
        assert terminal.penalty > internal.penalty
        assert (
            terminal.penalty_terms["variants_3p"] > internal.penalty_terms["variants"]
        )

    def test_panel_and_offtarget_sites_are_costed(self, ctx):
        p = score_primer(primer_with(paralog_sites=2, offtarget_sites=3), ctx)
        assert p.penalty_terms["paralog_sites"] == pytest.approx(
            2 * ctx.params.weights.paralog_priming_site
        )
        assert p.penalty_terms["offtarget_sites"] == pytest.approx(
            3 * ctx.params.weights.offtarget_site
        )

    def test_junction_spanning_is_a_bonus_not_a_penalty(self, ctx):
        plain = score_primer(primer_with(), ctx).penalty
        spanning = score_primer(primer_with(spans_junction=120), ctx)
        assert spanning.penalty_terms["junction_bonus"] < 0
        assert spanning.penalty < plain

    def test_ambiguous_bases_and_softmask(self, ctx):
        p = score_primer(primer_with(ambiguous=2, softmask_fraction=0.25), ctx)
        assert p.penalty_terms["ambiguous"] == pytest.approx(
            2 * ctx.params.weights.ambiguous_base
        )
        assert p.penalty_terms["softmask"] > 0

    def test_zero_terms_are_omitted_from_the_breakdown(self, ctx):
        p = score_primer(primer_with(tm=61.0), ctx)
        assert set(p.penalty_terms) == {"tm"}
        assert all(v != 0 for v in p.penalty_terms.values())

    def test_penalty_is_the_sum_of_its_terms(self, ctx):
        p = score_primer(
            primer_with(tm=62.0, gc=58.0, seq="A" * 25, end=25, paralog_sites=1), ctx
        )
        assert p.penalty == pytest.approx(sum(p.penalty_terms.values()))


class TestWeightCustomisation:
    def test_zeroing_a_weight_removes_its_influence(self, ctx):
        params = ctx.params.replace(weights=PenaltyWeights(tm_deviation=0.0))
        indifferent = TemplateContext(
            seq=ctx.seq, params=params, backend=ctx.backend
        )
        p = score_primer(primer_with(tm=70.0), indifferent)
        assert "tm" not in p.penalty_terms

    def test_raising_a_weight_changes_the_ranking(self, locus):
        from danioprime.design.pipeline import DesignRequest, run_design

        base = DesignParams().replace(max_pairs=3)
        size_obsessed = base.replace(
            weights=PenaltyWeights(product_size_deviation=5.0), opt_product_size=200
        )
        default_result = run_design(DesignRequest(template=locus, params=base))
        tuned_result = run_design(DesignRequest(template=locus, params=size_obsessed))
        assert default_result.pairs and tuned_result.pairs
        assert abs(tuned_result.pairs[0].product_size - 200) <= abs(
            default_result.pairs[0].product_size - 200
        )
