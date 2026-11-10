"""Thermodynamics: nearest-neighbour sums, Tm behaviour, structure detection."""

from __future__ import annotations

import math

import pytest

from danioprime.errors import SequenceError
from danioprime.thermo import nn, structure
from danioprime.thermo.backend import BuiltinBackend, get_backend, primer3_available
from danioprime.thermo.nn import SaltConditions


class TestNearestNeighbour:
    def test_polya_enthalpy_entropy_matches_hand_calculation(self):
        # 9 AA/TT stacks plus two A-T initiations:
        #   dH = 9 * -7.6 + 2 * 2.3 = -63.8
        #   dS = 9 * -21.3 + 2 * 4.1 = -183.5
        dh, ds = nn.duplex_enthalpy_entropy("A" * 10)
        assert dh == pytest.approx(-63.8, abs=1e-9)
        assert ds == pytest.approx(-183.5, abs=1e-9)

    def test_gc_palindrome_includes_symmetry_correction(self):
        # 5 GC + 4 CG stacks, two G-C initiations, plus the -1.4 symmetry term.
        dh, ds = nn.duplex_enthalpy_entropy("GCGCGCGCGC")
        assert dh == pytest.approx(5 * -9.8 + 4 * -10.6 + 2 * 0.1, abs=1e-9)
        assert ds == pytest.approx(5 * -24.4 + 4 * -27.2 + 2 * -2.8 - 1.4, abs=1e-9)

    def test_complement_pairs_have_identical_parameters(self):
        for seq in ("ACGTTGCATGCCAATGGCTA", "GGGCCCAAATTT", "ATATATATAT"):
            from danioprime.sequence import revcomp

            assert nn.duplex_enthalpy_entropy(seq) == pytest.approx(
                nn.duplex_enthalpy_entropy(revcomp(seq))
            )

    def test_self_complementarity_detection(self):
        assert nn.is_self_complementary("GAATTC")
        assert nn.is_self_complementary("GCGCGCGC")
        assert not nn.is_self_complementary("GAATTA")
        assert not nn.is_self_complementary("ACGTA")

    def test_requires_unambiguous_bases(self):
        with pytest.raises(SequenceError, match="unambiguous"):
            nn.tm("ACGTN")
        with pytest.raises(SequenceError, match="at least 2"):
            nn.tm("A")


class TestMeltingTemperature:
    def test_gc_rich_melts_higher_than_at_rich_at_equal_length(self):
        assert nn.tm("GCGCGCGCGCGCGCGCGCGC") > nn.tm("ATATATATATATATATATAT")

    def test_tm_trends_upward_with_length(self):
        # Not strictly monotonic, and that is correct rather than a defect:
        # appending a weak A-T pair to a duplex pulls the dH/dS ratio towards
        # the A-T asymptote, so Tm can fall by a few tenths of a degree. The
        # trend over several bases is what has to hold.
        base = "ACGTTGCATGCCAATGGCTAGCCTGAAC"
        assert nn.tm(base[:12]) < nn.tm(base[:18]) < nn.tm(base[:24])

    def test_appending_gc_raises_tm_more_than_appending_at(self):
        stem = "ACGTTGCATGCCAATGGCT"
        assert nn.tm(stem + "GC") > nn.tm(stem + "AT")

    def test_salt_raises_tm(self):
        seq = "ACGTTGCATGCCAATGGCTA"
        low = nn.tm(seq, SaltConditions(monovalent_mM=20, divalent_mM=0, dntp_mM=0))
        high = nn.tm(seq, SaltConditions(monovalent_mM=200, divalent_mM=0, dntp_mM=0))
        assert high > low

    def test_divalent_folds_into_equivalent_monovalent(self):
        cond = SaltConditions(monovalent_mM=50, divalent_mM=1.5, dntp_mM=0.6)
        # 50 + 120 * sqrt(0.9)
        assert cond.equivalent_monovalent_mM() == pytest.approx(50 + 120 * math.sqrt(0.9))

    def test_dntp_chelation_lowers_effective_salt(self):
        free = SaltConditions(monovalent_mM=50, divalent_mM=2.0, dntp_mM=0.0)
        chelated = SaltConditions(monovalent_mM=50, divalent_mM=2.0, dntp_mM=1.8)
        assert chelated.equivalent_monovalent_mM() < free.equivalent_monovalent_mM()

    def test_dmso_lowers_tm(self):
        seq = "ACGTTGCATGCCAATGGCTA"
        assert nn.tm(seq, SaltConditions(dmso_pct=5.0)) < nn.tm(seq)

    def test_higher_primer_concentration_raises_tm(self):
        seq = "ACGTTGCATGCCAATGGCTA"
        assert nn.tm(seq, SaltConditions(primer_nM=1000)) > nn.tm(
            seq, SaltConditions(primer_nM=50)
        )

    def test_typical_20mer_lands_in_the_pcr_range(self):
        # Sanity bound rather than a reference value: a 50% GC 20-mer under
        # standard PCR conditions should anneal somewhere in the 55-70 C band.
        assert 55.0 < nn.tm("ACGTTGCATGCCAATGGCTA") < 70.0


class TestFreeEnergy:
    def test_delta_g_is_negative_for_a_stable_duplex(self):
        assert nn.delta_g("ACGTTGCATGCCAATGGCTA") < -15.0

    def test_gc_rich_three_prime_end_is_more_stable(self):
        stable = nn.end_stability("AAAAAAAAAAAAAAAGGCGC")
        unstable = nn.end_stability("GGCGCAAAAAAAAAAAAAAA")
        assert stable < unstable


class TestImperfectDuplex:
    def test_perfect_register_reproduces_tm(self):
        seq = "ACGTTGCATGCCAATGGCTA"
        assert nn.imperfect_duplex_tm(seq, seq) == pytest.approx(nn.tm(seq))

    def test_internal_mismatch_lowers_tm(self):
        primer = "ACGTTGCATGCCAATGGCTA"
        mismatched = "ACGTTGCATACCAATGGCTA"
        assert nn.imperfect_duplex_tm(primer, mismatched) < nn.tm(primer)

    def test_more_mismatches_lower_tm_further(self):
        primer = "ACGTTGCATGCCAATGGCTA"
        one = nn.imperfect_duplex_tm(primer, "ACGTTGCATACCAATGGCTA")
        two = nn.imperfect_duplex_tm(primer, "ACGTTACATACCAATGGCTA")
        assert two < one

    def test_register_length_must_match(self):
        with pytest.raises(SequenceError, match="same length"):
            nn.imperfect_duplex_tm("ACGTACGTAC", "ACGTACGT")


class TestStructure:
    def test_palindrome_self_dimerises_strongly(self):
        result = structure.self_dimer("GCGCGCGCGC")
        assert result.stack_length == 10
        assert result.dg < -15.0

    def test_polya_has_no_self_dimer(self):
        result = structure.self_dimer("A" * 14)
        assert result.stack_length == 0
        assert result.dg == 0.0

    def test_designed_hairpin_is_found(self):
        result = structure.hairpin("GGGGGGAAATTTCCCCCC")
        assert result.stack_length >= 6
        assert result.dg < -7.0

    def test_loop_minimum_is_respected(self):
        # Complementary arms separated by only 2 nt cannot form a stem.
        result = structure.hairpin("GGGGCCGGGG"[:10], min_loop=5)
        assert result.dg > -8.0

    def test_complementary_primers_hetero_dimerise(self):
        from danioprime.sequence import revcomp

        a = "ACGTTGCATGCCAATGGCTA"
        result = structure.hetero_dimer(a, revcomp(a))
        assert result.dg < -20.0

    def test_three_prime_dimer_score_detects_facing_ends(self):
        # 3' ends fully complementary over the window.
        high = structure.three_prime_dimer_score("AAAAAAAAAAGGGGG", "AAAAAAAAAACCCCC")
        low = structure.three_prime_dimer_score("AAAAAAAAAAGGGGG", "AAAAAAAAAAGGGGG")
        assert high > low


class TestBackends:
    def test_builtin_backend_is_default_and_consistent(self):
        backend = get_backend("builtin")
        assert isinstance(backend, BuiltinBackend)
        seq = "ACGTTGCATGCCAATGGCTA"
        assert backend.tm(seq) == pytest.approx(nn.tm(seq))

    def test_unknown_backend_is_rejected(self):
        from danioprime.errors import DanioPrimeError

        with pytest.raises(DanioPrimeError, match="unknown thermo backend"):
            get_backend("nonsense")

    @pytest.mark.skipif(not primer3_available(), reason="primer3-py not installed")
    def test_primer3_backend_agrees_within_a_few_degrees(self):
        seq = "ACGTTGCATGCCAATGGCTA"
        builtin = get_backend("builtin").tm(seq)
        p3 = get_backend("primer3").tm(seq)
        assert abs(builtin - p3) < 3.0


class TestPrimer3BackendGuard:
    def test_absent_primer3_gives_an_actionable_error(self, monkeypatch):
        """The error must name the install command, not just fail to import."""
        import builtins

        from danioprime.errors import DanioPrimeError
        from danioprime.thermo.backend import Primer3Backend

        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "primer3":
                raise ImportError("no module named primer3")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        with pytest.raises(DanioPrimeError, match=r"danioprime\[primer3\]"):
            Primer3Backend()

    def test_backend_names_are_exposed(self):
        from danioprime.thermo.backend import BACKEND_NAMES

        assert set(BACKEND_NAMES) == {"builtin", "primer3"}
