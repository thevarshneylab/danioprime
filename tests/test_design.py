"""Candidate enumeration, filtering, pairing and the design pipeline."""

from __future__ import annotations

import pytest

from danioprime.config import DesignParams, get_preset
from danioprime.design.candidate import (
    Primer,
    TemplateContext,
    check_composition,
    count_gc_clamp,
    enumerate_candidates,
    evaluate_primer,
)
from danioprime.design.pairing import build_pairs, diversify, pair_is_geometrically_valid
from danioprime.design.pipeline import DesignRequest, run_design, windows_around
from danioprime.design.regions import Interval, IntervalSet
from danioprime.errors import DesignError
from danioprime.sequence import SeqRecord, revcomp
from danioprime.thermo.backend import get_backend


def ctx_for(seq: str, params: DesignParams) -> TemplateContext:
    return TemplateContext(
        seq=seq,
        params=params,
        backend=get_backend(params.thermo_backend, params.conditions, params.anneal_temp_c),
    )


class TestPrimerModel:
    def test_gc_clamp_counting(self):
        assert count_gc_clamp("AAAAGC") == 2
        assert count_gc_clamp("AAAAAA") == 0
        assert count_gc_clamp("GCGCGC") == 6

    def test_forward_and_reverse_coordinates_and_sequence(self, locus):
        fwd = Primer(seq=locus.seq[100:120], start=100, end=120, orientation="F")
        rev = Primer(seq=revcomp(locus.seq[100:120]), start=100, end=120, orientation="R")
        assert fwd.three_prime_pos == 119
        assert rev.three_prime_pos == 100
        assert fwd.five_prime_pos == 100
        assert rev.five_prime_pos == 119
        assert rev.seq == revcomp(locus.seq[100:120])

    def test_three_prime_window_follows_orientation(self):
        fwd = Primer(seq="A" * 20, start=100, end=120, orientation="F")
        rev = Primer(seq="A" * 20, start=100, end=120, orientation="R")
        assert fwd.three_prime_window(5) == Interval(115, 120)
        assert rev.three_prime_window(5) == Interval(100, 105)

    def test_ordered_sequence_includes_tail(self):
        p = Primer(seq="ACGTACGTACGTACGTACGT", start=0, end=20, orientation="F")
        p.tail = "TGTAAAACGACGGCCAGT"
        assert p.ordered_seq.startswith("TGTAAAACGACGGCCAGT")
        assert p.ordered_seq.endswith("ACGTACGTACGTACGTACGT")


class TestFilters:
    def test_homopolymer_is_rejected(self, params):
        # Chosen to pass the clamp and GC filters so the poly-X rule is what
        # fires: 45% GC, a 3' CC clamp, and a run of seven A.
        seq = "AAAAAAAGCGCGCGTATAGC"
        ctx = ctx_for(seq + "ACGT" * 15, params)
        p = Primer(seq=seq, start=0, end=20, orientation="F")
        evaluate_primer(p, ctx)
        ok, reason = check_composition(p, ctx)
        assert not ok
        assert reason == "poly-x"

    def test_missing_gc_clamp_is_rejected_before_cheaper_checks(self, params):
        ctx = ctx_for("A" * 60, params)
        p = Primer(seq="A" * 20, start=0, end=20, orientation="F")
        evaluate_primer(p, ctx)
        ok, reason = check_composition(p, ctx)
        assert not ok
        assert reason == "gc-clamp-missing"

    def test_excluded_region_blocks_a_primer(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        ctx.excluded = IntervalSet([Interval(100, 200)])
        p = Primer(seq=locus.seq[150:170], start=150, end=170, orientation="F")
        evaluate_primer(p, ctx)
        ok, reason = check_composition(p, ctx)
        assert not ok
        assert reason == "excluded-region"

    def test_variant_at_three_prime_end_is_always_disqualifying(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        ctx.variants = IntervalSet([Interval(168, 169)])
        p = Primer(seq=locus.seq[150:170], start=150, end=170, orientation="F")
        evaluate_primer(p, ctx)
        ok, reason = check_composition(p, ctx)
        assert not ok
        assert reason == "variant-at-3prime"

    def test_internal_variant_is_penalised_not_rejected_by_default(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        ctx.variants = IntervalSet([Interval(155, 156)])
        p = Primer(seq=locus.seq[150:170], start=150, end=170, orientation="F")
        evaluate_primer(p, ctx)
        ok, _ = check_composition(p, ctx)
        assert ok
        assert p.variants_overlapping == 1

    def test_reject_variant_overlap_makes_it_disqualifying(self, locus, params):
        ctx = ctx_for(locus.seq, params.replace(reject_variant_overlap=True))
        ctx.variants = IntervalSet([Interval(155, 156)])
        p = Primer(seq=locus.seq[150:170], start=150, end=170, orientation="F")
        evaluate_primer(p, ctx)
        ok, reason = check_composition(p, ctx)
        assert not ok
        assert reason == "variant-overlap"

    def test_ambiguous_bases_are_rejected_by_default(self, params):
        ctx = ctx_for("ACGTN" * 12, params)
        p = Primer(seq="ACGTNACGTNACGTNACGTN", start=0, end=20, orientation="F")
        evaluate_primer(p, ctx)
        ok, reason = check_composition(p, ctx)
        assert not ok
        assert reason == "ambiguous base"


class TestEnumeration:
    def test_candidates_are_scored_and_sorted(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        result = enumerate_candidates(ctx, "F")
        assert result.primers
        penalties = [p.penalty for p in result.primers]
        assert penalties == sorted(penalties)

    def test_every_candidate_satisfies_the_constraints(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        for p in enumerate_candidates(ctx, "F").primers[:100]:
            assert params.min_length <= p.length <= params.max_length
            assert params.min_tm <= p.tm <= params.max_tm
            assert params.min_gc <= p.gc <= params.max_gc
            assert p.gc_clamp >= params.require_gc_clamp
            assert p.hairpin_dg >= params.max_hairpin_dg

    def test_reverse_primers_are_reverse_complements_of_the_template(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        for p in enumerate_candidates(ctx, "R").primers[:20]:
            assert p.seq == revcomp(locus.seq[p.start : p.end])

    def test_search_windows_are_respected(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        window = IntervalSet([Interval(500, 700)])
        for p in enumerate_candidates(ctx, "F", window).primers:
            assert 500 <= p.start and p.end <= 700

    def test_rejection_reasons_are_recorded(self, params):
        ctx = ctx_for("AT" * 400, params)
        result = enumerate_candidates(ctx, "F")
        assert not result.primers
        assert result.rejections
        assert "no candidates were filtered" not in result.reason_summary()

    def test_truncation_is_flagged(self, locus, params):
        ctx = ctx_for(locus.seq, params.replace(max_candidates_per_side=5))
        result = enumerate_candidates(ctx, "F")
        assert result.truncated
        assert len(result.primers) == 5


class TestPairing:
    def test_geometry_rules(self, params):
        fwd = Primer(seq="A" * 20, start=100, end=120, orientation="F")
        rev_ok = Primer(seq="A" * 20, start=400, end=420, orientation="R")
        rev_overlap = Primer(seq="A" * 20, start=110, end=130, orientation="R")
        rev_inverted = Primer(seq="A" * 20, start=10, end=30, orientation="R")
        assert pair_is_geometrically_valid(fwd, rev_ok, params)[0]
        assert pair_is_geometrically_valid(fwd, rev_overlap, params)[1] == "primers overlap"
        assert pair_is_geometrically_valid(fwd, rev_inverted, params)[1] == "inverted"
        assert pair_is_geometrically_valid(rev_ok, fwd, params)[1] == "orientation"

    def test_product_size_window_is_enforced(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        f = enumerate_candidates(ctx, "F").primers[:40]
        r = enumerate_candidates(ctx, "R").primers[:40]
        for pair in build_pairs(f, r, ctx):
            assert params.min_product_size <= pair.product_size <= params.max_product_size

    def test_pair_tm_difference_is_enforced(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        f = enumerate_candidates(ctx, "F").primers[:40]
        r = enumerate_candidates(ctx, "R").primers[:40]
        for pair in build_pairs(f, r, ctx):
            assert pair.tm_difference <= params.max_pair_tm_difference

    def test_diversify_spreads_three_prime_ends(self):
        from danioprime.design.pairing import PrimerPair

        pairs = []
        for shift in range(6):
            f = Primer(seq="A" * 20, start=100 + shift, end=120 + shift, orientation="F")
            r = Primer(seq="A" * 20, start=400, end=420, orientation="R")
            pairs.append(PrimerPair(forward=f, reverse=r, product_size=320 - shift))
        kept = diversify(pairs, window=8, limit=10)
        assert len(kept) == 1

    def test_amplicon_sequence_matches_the_template(self, locus, params):
        ctx = ctx_for(locus.seq, params)
        f = enumerate_candidates(ctx, "F").primers[:20]
        r = enumerate_candidates(ctx, "R").primers[:20]
        pair = build_pairs(f, r, ctx)[0]
        amplicon = pair.amplicon_seq(locus.seq)
        assert amplicon.startswith(pair.forward.seq)
        assert amplicon.endswith(revcomp(pair.reverse.seq))
        assert len(amplicon) == pair.product_size


class TestPipeline:
    def test_target_must_be_spanned(self, locus, params):
        target = Interval(900, 960)
        result = run_design(
            DesignRequest(
                template=locus, params=params, targets=IntervalSet([target]), name="t"
            )
        )
        assert result.pairs
        for pair in result.pairs:
            assert pair.forward.start < target.start
            assert pair.reverse.end > target.end

    def test_pairs_are_named_and_sorted(self, locus, params):
        result = run_design(DesignRequest(template=locus, params=params, name="gene"))
        assert [p.name for p in result.pairs] == [
            f"gene_p{i}" for i in range(1, len(result.pairs) + 1)
        ]
        assert [p.penalty for p in result.pairs] == sorted(p.penalty for p in result.pairs)

    def test_impossible_constraints_produce_a_diagnosis(self, locus):
        params = get_preset("genotyping-gel").replace(min_tm=85.0, opt_tm=88.0, max_tm=90.0)
        result = run_design(DesignRequest(template=locus, params=params))
        assert not result.pairs
        assert result.warnings
        assert "tm" in result.diagnosis().lower()

    def test_template_shorter_than_the_product_is_an_error(self, params):
        with pytest.raises(DesignError, match="shorter than the minimum"):
            run_design(
                DesignRequest(template=SeqRecord(name="tiny", seq="ACGT" * 10), params=params)
            )

    def test_statistics_are_reported(self, locus, params):
        result = run_design(DesignRequest(template=locus, params=params))
        assert result.stats["template_length"] == len(locus.seq)
        assert result.stats["forward_candidates"] > 0

    def test_windows_around_a_point(self):
        fwd, rev = windows_around(1000, 3000, inner=50, outer=400)
        assert fwd.intervals == [Interval(600, 950)]
        assert rev.intervals == [Interval(1050, 1400)]

    def test_windows_around_validates_its_arguments(self):
        with pytest.raises(DesignError, match="inner < outer"):
            windows_around(100, 1000, inner=400, outer=50)


class TestConfig:
    def test_every_preset_validates(self):
        from danioprime.config import PRESETS

        for name in PRESETS:
            get_preset(name).validate()

    def test_json_round_trip(self, tmp_path):
        params = get_preset("qpcr")
        path = tmp_path / "p.json"
        params.to_json(path)
        assert DesignParams.from_json(path) == params

    def test_unknown_parameter_is_rejected(self):
        from danioprime.errors import DanioPrimeError

        with pytest.raises(DanioPrimeError, match="unknown design parameters"):
            DesignParams().replace(nonsense=1)

    def test_contradictory_parameters_are_caught(self):
        from danioprime.errors import DanioPrimeError

        with pytest.raises(DanioPrimeError, match="Tm must satisfy"):
            DesignParams().replace(min_tm=70.0, opt_tm=60.0, max_tm=65.0).validate()
        with pytest.raises(DanioPrimeError, match="product size"):
            DesignParams().replace(min_product_size=900, max_product_size=100).validate()

    def test_unknown_preset_lists_the_alternatives(self):
        from danioprime.errors import DanioPrimeError

        with pytest.raises(DanioPrimeError, match="qpcr"):
            get_preset("nope")
