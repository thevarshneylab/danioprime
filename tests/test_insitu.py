"""Riboprobe templates and probe-level cross-hybridisation.

The distinction this module exists for: primer specificity is decided at a 3'
end, probe specificity along the whole amplicon. A pair can pass the first and
fail the second, and these tests assert exactly that case.
"""

from __future__ import annotations

import random

import pytest
from helpers import make_seq, mutate

from danioprime.assays.insitu import (
    CROSS_HYB_IDENTITY,
    _longest_window_above,
    assess_probe,
    design_insitu,
    find_conserved_blocks,
)
from danioprime.assays.tails import RIBOPROBE_SCHEMES
from danioprime.config import get_preset
from danioprime.design.regions import Interval, IntervalSet
from danioprime.errors import DesignError
from danioprime.sequence import SeqRecord, revcomp
from danioprime.specificity.paralog import ParalogPanel


@pytest.fixture
def ohnolog_transcripts():
    """A transcript whose middle is conserved and whose ends have diverged.

    This is the shape that breaks a primer-only analysis: primers placed in the
    diverged ends discriminate perfectly, while an amplicon spanning the middle
    yields a probe that hybridises to both genes.
    """
    left = make_seq(700, 11)
    middle = make_seq(900, 22)
    right = make_seq(700, 33)
    target = left + middle + right
    paralog = mutate(left, 0.35, 1) + mutate(middle, 0.04, 2) + mutate(right, 0.35, 3)
    return (
        SeqRecord(name="kcnh6a", seq=target),
        ParalogPanel.from_records([SeqRecord(name="kcnh6b", seq=paralog)]),
        (700, 1600),
    )


class TestLongestWindow:
    """The window finder had two real defects; both are pinned here."""

    @pytest.mark.parametrize("seed", range(25))
    def test_agrees_with_brute_force(self, seed):
        rng = random.Random(seed)
        n = rng.randint(1, 40)
        matches = [rng.random() < rng.choice([0.3, 0.6, 0.9]) for _ in range(n)]
        threshold = rng.choice([0.0, 0.5, 0.75, 0.8, 0.9, 1.0])
        got = _longest_window_above(matches, threshold)
        best = 0
        for i in range(n):
            for j in range(i + 1, n + 1):
                window = matches[i:j]
                if sum(window) / len(window) >= threshold - 1e-12:
                    best = max(best, j - i)
        assert (got[1] - got[0] if got else 0) == best

    def test_floating_point_case_that_used_to_fail(self):
        """Prefix sums of 0.2 and -0.8 collide by an ulp; integers fixed it."""
        matches = [True] * 3 + [False] * 2 + [True] * 3 + [False, True, False, False, True]
        got = _longest_window_above(matches, 0.8)
        assert got is not None
        assert got[1] - got[0] == 5

    def test_empty_and_all_mismatch(self):
        assert _longest_window_above([], 0.8) is None
        assert _longest_window_above([False] * 10, 0.8) is None

    def test_all_match_returns_everything(self):
        assert _longest_window_above([True] * 10, 0.8) == (0, 10)


class TestConservedBlocks:
    def test_reported_blocks_always_meet_the_threshold(self, ohnolog_transcripts):
        """The original seed-extension took an unvalidated span between two
        chance seeds and reported 588 nt at 64% identity."""
        target, panel, _ = ohnolog_transcripts
        for start in (0, 300, 700, 1200):
            for block in find_conserved_blocks(target.seq[start : start + 600], panel):
                assert block.identity >= CROSS_HYB_IDENTITY - 1e-9
                assert block.length > 0

    def test_a_conserved_region_is_found(self, ohnolog_transcripts):
        target, panel, (lo, hi) = ohnolog_transcripts
        blocks = find_conserved_blocks(target.seq[lo:hi], panel)
        assert blocks
        assert max(b.length for b in blocks) > 500

    def test_a_diverged_region_yields_no_long_block(self, ohnolog_transcripts):
        target, panel, _ = ohnolog_transcripts
        blocks = find_conserved_blocks(target.seq[0:700], panel)
        assert all(b.length < 200 for b in blocks)

    def test_block_coordinates_are_in_probe_space(self, ohnolog_transcripts):
        target, panel, (lo, hi) = ohnolog_transcripts
        probe = target.seq[lo:hi]
        for block in find_conserved_blocks(probe, panel):
            assert 0 <= block.probe_start < block.probe_end <= len(probe)

    def test_reverse_complement_similarity_is_detected(self):
        """A probe is single-stranded RNA; an inverted paralog still binds it."""
        core = make_seq(400, 77)
        panel = ParalogPanel.from_records([SeqRecord(name="inverted", seq=revcomp(core))])
        blocks = find_conserved_blocks(core, panel)
        assert blocks
        assert max(b.length for b in blocks) > 300

    def test_no_panel_means_no_blocks(self):
        assert find_conserved_blocks(make_seq(500, 1), ParalogPanel()) == []

    def test_probe_shorter_than_the_seed(self):
        panel = ParalogPanel.from_records([SeqRecord(name="x", seq=make_seq(300, 2))])
        assert find_conserved_blocks("ACGT", panel) == []


class TestAssessProbe:
    def test_conserved_probe_is_flagged(self, ohnolog_transcripts):
        target, panel, (lo, hi) = ohnolog_transcripts
        result = assess_probe(target.seq[lo:hi], panel)
        assert result.will_cross_hybridise
        assert "cross-hybridise" in result.verdict()
        assert "probe body" in result.verdict()

    def test_diverged_probe_is_cleared(self, ohnolog_transcripts):
        target, panel, _ = ohnolog_transcripts
        result = assess_probe(target.seq[0:700], panel)
        assert not result.will_cross_hybridise

    def test_no_panel_is_always_clean(self):
        result = assess_probe(make_seq(600, 5), None)
        assert not result.will_cross_hybridise
        assert "no similarity" in result.verdict()

    def test_thresholds_are_honoured(self, ohnolog_transcripts):
        target, panel, (lo, hi) = ohnolog_transcripts
        probe = target.seq[lo:hi]
        assert assess_probe(probe, panel, min_length=50).will_cross_hybridise
        assert not assess_probe(probe, panel, min_length=10_000).will_cross_hybridise
        assert not assess_probe(probe, panel, min_identity=0.999).will_cross_hybridise

    def test_serialises(self, ohnolog_transcripts):
        import json

        target, panel, (lo, hi) = ohnolog_transcripts
        json.dumps(assess_probe(target.seq[lo:hi], panel).to_dict())


class TestDesignInSitu:
    @pytest.fixture
    def params(self):
        return get_preset("insitu-probe").replace(
            max_pairs=3, min_product_size=450, opt_product_size=550, max_product_size=650
        )

    def test_antisense_promoter_lands_on_the_reverse_primer(self, ohnolog_transcripts, params):
        target, panel, _ = ohnolog_transcripts
        design = design_insitu(target, params, scheme="riboprobe-t7", panel=panel)
        assert design.pairs
        assert design.antisense
        for oligo in design.antisense.oligos:
            if oligo.role == "reverse":
                assert oligo.sequence.startswith("TAATACGACTCACTATAGGG")
            elif oligo.role == "forward":
                assert not oligo.sequence.startswith("TAATACGACTCACTATAGGG")

    def test_sense_control_is_emitted_as_a_separate_plan(self, ohnolog_transcripts, params):
        target, panel, _ = ohnolog_transcripts
        design = design_insitu(target, params, scheme="riboprobe-t3", panel=panel)
        assert design.sense is not None
        forward = [o for o in design.sense.oligos if o.role == "forward"]
        assert all(o.sequence.startswith("AATTAACCCTCACTAAAGGG") for o in forward)
        assert any("second PCR" in w for w in design.result.warnings)

    def test_sense_control_can_be_suppressed(self, ohnolog_transcripts, params):
        target, panel, _ = ohnolog_transcripts
        design = design_insitu(
            target, params, scheme="riboprobe-t7", panel=panel, sense_control=False
        )
        assert design.sense is None

    def test_dual_scheme_needs_no_separate_sense_plan(self, ohnolog_transcripts, params):
        target, panel, _ = ohnolog_transcripts
        design = design_insitu(target, params, scheme="riboprobe-dual-t7-sp6", panel=panel)
        assert design.sense is None
        assert design.antisense
        fwd = next(o for o in design.antisense.oligos if o.role == "forward")
        rev = next(o for o in design.antisense.oligos if o.role == "reverse")
        assert fwd.tail == "SP6"
        assert rev.tail == "T7"

    def test_cross_hybridising_probes_are_flagged_and_demoted(self, ohnolog_transcripts, params):
        target, panel, _ = ohnolog_transcripts
        design = design_insitu(target, params, scheme="riboprobe-t7", panel=panel)
        assert design.cross_hyb
        assert all(p.name in design.cross_hyb for p in design.pairs)
        flags = [design.cross_hyb[p.name].will_cross_hybridise for p in design.pairs]
        # Clean templates must be ranked ahead of cross-hybridising ones.
        assert flags == sorted(flags)

    def test_a_pair_can_pass_primer_specificity_and_fail_as_a_probe(
        self, ohnolog_transcripts, params
    ):
        """The reason this module exists, asserted directly."""
        target, panel, (lo, hi) = ohnolog_transcripts
        design = design_insitu(
            target,
            params,
            scheme="riboprobe-t7",
            panel=panel,
            targets=IntervalSet([Interval(lo + 100, lo + 300)]),
        )
        assert design.pairs
        offenders = [
            p
            for p in design.pairs
            if p.is_specific and design.cross_hyb[p.name].will_cross_hybridise
        ]
        assert offenders, "expected a pair whose primers are specific but whose probe is not"

    def test_steering_into_a_divergent_region_gives_a_clean_probe(
        self, ohnolog_transcripts, params
    ):
        target, panel, _ = ohnolog_transcripts
        design = design_insitu(
            target,
            params,
            scheme="riboprobe-t7",
            panel=panel,
            targets=IntervalSet([Interval(200, 400)]),
        )
        assert design.pairs
        assert not design.cross_hyb[design.pairs[0].name].will_cross_hybridise

    def test_warns_when_every_template_cross_hybridises(self, ohnolog_transcripts, params):
        target, panel, (lo, hi) = ohnolog_transcripts
        design = design_insitu(
            target,
            params,
            scheme="riboprobe-t7",
            panel=panel,
            targets=IntervalSet([Interval(lo + 100, lo + 400)]),
        )
        assert any("3' UTR" in w for w in design.result.warnings)

    def test_short_probe_warning(self, ohnolog_transcripts):
        target, panel, _ = ohnolog_transcripts
        params = get_preset("insitu-probe").replace(
            max_pairs=2, min_product_size=150, opt_product_size=200, max_product_size=250
        )
        design = design_insitu(target, params, scheme="riboprobe-t7")
        assert any("weak signal" in w for w in design.result.warnings)

    def test_a_scheme_without_a_promoter_is_refused(self, ohnolog_transcripts, params):
        target, _panel, _ = ohnolog_transcripts
        with pytest.raises(DesignError, match="attaches no RNA-polymerase promoter"):
            design_insitu(target, params, scheme="fluorescent")

    def test_unknown_scheme_lists_the_riboprobe_options(self, ohnolog_transcripts, params):
        target, _panel, _ = ohnolog_transcripts
        with pytest.raises(DesignError, match="riboprobe-t7"):
            design_insitu(target, params, scheme="nonsense")

    @pytest.mark.parametrize("scheme", sorted(RIBOPROBE_SCHEMES))
    def test_every_riboprobe_scheme_designs(self, ohnolog_transcripts, params, scheme):
        target, _panel, _ = ohnolog_transcripts
        design = design_insitu(target, params, scheme=scheme)
        assert design.pairs
        assert design.antisense
        tails = {o.tail for o in design.antisense.oligos if o.tail}
        assert tails, f"{scheme} attached no tail"

    def test_serialises(self, ohnolog_transcripts, params):
        import json

        target, panel, _ = ohnolog_transcripts
        design = design_insitu(target, params, scheme="riboprobe-t7", panel=panel)
        json.dumps(design.to_dict())


class TestThresholdCalibration:
    """The defaults are calibrated; this pins the calibration.

    The first version used 80% over 50 nt, which fired on more than half of
    all diverged-paralog backgrounds purely by chance. A cross-hybridisation
    warning that is wrong half the time is worse than none.
    """

    @staticmethod
    def _chance_trips(background: float, identity: float, min_length: int, trials: int) -> int:
        rng = random.Random(4242)
        trips = 0
        for _ in range(trials):
            matches = [rng.random() < background for _ in range(700)]
            window = _longest_window_above(matches, identity)
            if window and (window[1] - window[0]) >= min_length:
                trips += 1
        return trips

    def test_defaults_are_quiet_against_a_diverged_paralog(self):
        from danioprime.assays.insitu import CROSS_HYB_MIN_LENGTH

        trips = self._chance_trips(0.65, CROSS_HYB_IDENTITY, CROSS_HYB_MIN_LENGTH, 60)
        assert trips <= 1, f"{trips}/60 chance trips at the default threshold"

    def test_defaults_are_silent_against_unrelated_sequence(self):
        from danioprime.assays.insitu import CROSS_HYB_MIN_LENGTH

        assert self._chance_trips(0.25, CROSS_HYB_IDENTITY, CROSS_HYB_MIN_LENGTH, 40) == 0

    def test_the_rejected_lenient_threshold_really_was_noisy(self):
        """Documents why the defaults are what they are."""
        assert self._chance_trips(0.65, 0.80, 50, 60) > 20

    def test_real_conservation_still_fires(self):
        """Calibration must not have made the check useless."""
        core = make_seq(800, 909)
        panel = ParalogPanel.from_records(
            [SeqRecord(name="ohnolog", seq=mutate(core, 0.05, 8))]
        )
        assert assess_probe(core, panel).will_cross_hybridise
