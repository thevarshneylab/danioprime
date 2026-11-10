"""Binding-site detection, priming rules, ohnolog discrimination."""

from __future__ import annotations

import pytest
from helpers import make_seq

from danioprime.config import get_preset
from danioprime.design.pipeline import DesignRequest, run_design
from danioprime.sequence import SeqRecord, revcomp
from danioprime.specificity.external import _md_mismatch_offsets, available_aligners
from danioprime.specificity.index import KmerIndex
from danioprime.specificity.paralog import ParalogPanel, assess_pair, assess_primer
from danioprime.specificity.scan import (
    PrimingRules,
    find_binding_sites,
    find_binding_sites_seeded,
    predict_amplicons,
)

PRIMER = "ACGTTGCATGCCAATGGCTA"


@pytest.fixture
def planted() -> SeqRecord:
    """A subject with a perfect plus site, a perfect minus site, and a decoy."""
    bg = make_seq(6000, 404)
    seq = (
        bg[:1000]
        + PRIMER
        + bg[1000:3000]
        + PRIMER[:-1] + ("G" if PRIMER[-1] != "G" else "T")  # 3'-terminal mismatch
        + bg[3000:5000]
        + revcomp(PRIMER)
        + bg[5000:]
    )
    return SeqRecord(name="subject", seq=seq)


class TestBindingSites:
    def test_perfect_sites_are_found_on_both_strands(self, planted):
        sites = find_binding_sites(PRIMER, planted, primer_name="P")
        perfect = [s for s in sites if s.is_perfect]
        assert len(perfect) == 2
        assert {s.strand for s in perfect} == {"+", "-"}
        assert all(s.can_prime for s in perfect)

    def test_plus_site_coordinates_and_three_prime_end(self, planted):
        sites = [s for s in find_binding_sites(PRIMER, planted) if s.is_perfect and s.strand == "+"]
        site = sites[0]
        assert planted.seq[site.start : site.end] == PRIMER
        assert site.three_prime_pos == site.end - 1

    def test_minus_site_three_prime_end_points_the_other_way(self, planted):
        sites = [s for s in find_binding_sites(PRIMER, planted) if s.is_perfect and s.strand == "-"]
        site = sites[0]
        assert revcomp(planted.seq[site.start : site.end]) == PRIMER
        assert site.three_prime_pos == site.start

    def test_three_prime_terminal_mismatch_blocks_priming(self, planted):
        sites = find_binding_sites(PRIMER, planted, max_mismatches=2)
        decoys = [s for s in sites if s.mismatches == 1 and 0 in s.mismatch_offsets_from_3p]
        assert decoys
        for site in decoys:
            assert not site.can_prime
            assert "3' end" in site.note

    def test_internal_mismatches_still_prime(self):
        # Three mismatches, all 12 or more bases from the 3' end: the duplex is
        # weakened but the polymerase has a correctly paired 3'-OH, so this
        # site primes. A BLAST-identity filter would have discarded it.
        mutated = list(PRIMER)
        for i in (1, 4, 7):
            mutated[i] = "A" if PRIMER[i] != "A" else "C"
        subject = "TTTT" + "".join(mutated) + "TTTT"
        sites = find_binding_sites(PRIMER, subject, max_mismatches=6)
        priming = [s for s in sites if s.can_prime]
        assert priming
        assert priming[0].mismatches == 3
        assert priming[0].three_prime_mismatches == 0
        assert min(priming[0].mismatch_offsets_from_3p) >= 5

    def test_priming_only_filters_the_output(self, planted):
        all_sites = find_binding_sites(PRIMER, planted, max_mismatches=3)
        only = find_binding_sites(PRIMER, planted, max_mismatches=3, priming_only=True)
        assert len(only) < len(all_sites)
        assert all(s.can_prime for s in only)

    def test_mismatch_budget_is_respected(self, planted):
        for site in find_binding_sites(PRIMER, planted, max_mismatches=2):
            assert site.mismatches <= 2

    def test_rules_can_be_tightened(self, planted):
        strict = PrimingRules(clean_three_prime=10, max_mismatches=0)
        sites = find_binding_sites(PRIMER, planted, rules=strict, max_mismatches=3)
        decoys = [s for s in sites if s.mismatches == 1]
        assert all(not s.can_prime for s in decoys)


class TestSeededSearch:
    def test_seeded_agrees_with_exhaustive_on_priming_sites(self, planted):
        index = KmerIndex.build([planted], k=12)
        seeded = find_binding_sites_seeded(PRIMER, index)
        exhaustive = [s for s in find_binding_sites(PRIMER, planted) if s.can_prime]
        assert sorted((s.start, s.strand) for s in seeded) == sorted(
            (s.start, s.strand) for s in exhaustive
        )

    def test_primer_shorter_than_the_seed_is_an_error(self, planted):
        index = KmerIndex.build([planted], k=12)
        with pytest.raises(ValueError, match="shorter than the seed"):
            find_binding_sites_seeded("ACGTACGT", index)

    def test_index_skips_windows_containing_n(self):
        rec = SeqRecord(name="n", seq="ACGTACGTNACGTACGT")
        index = KmerIndex.build([rec], k=8)
        assert index.lookup("ACGTACGT")
        assert all("N" not in kmer for kmer in index.table)

    def test_index_stats(self, planted):
        index = KmerIndex.build([planted], k=12)
        stats = index.stats()
        assert stats["subjects"] == 1
        assert stats["seed_length"] == 12
        assert stats["total_length"] == len(planted.seq)


class TestAmpliconPrediction:
    def test_two_facing_sites_predict_a_product(self, planted):
        f = find_binding_sites(PRIMER, planted, priming_only=True, primer_name="F")
        r = find_binding_sites(revcomp(PRIMER), planted, priming_only=True, primer_name="R")
        amplicons = predict_amplicons({"F": f, "R": r}, min_size=50, max_size=10000)
        assert amplicons
        assert all(a.size <= 10000 for a in amplicons)

    def test_self_pair_products_are_detected_and_labelled(self, planted):
        sites = find_binding_sites(PRIMER, planted, priming_only=True, primer_name="P")
        amplicons = predict_amplicons({"P": sites}, min_size=50, max_size=10000)
        assert amplicons
        assert any(a.is_same_primer for a in amplicons)
        assert any("self-pair" in a.note for a in amplicons)

    def test_size_window_is_enforced(self, planted):
        sites = find_binding_sites(PRIMER, planted, priming_only=True, primer_name="P")
        amplicons = predict_amplicons({"P": sites}, min_size=50, max_size=100)
        assert amplicons == []

    def test_non_priming_sites_are_ignored_by_default(self, planted):
        sites = find_binding_sites(PRIMER, planted, max_mismatches=3, primer_name="P")
        strict = predict_amplicons({"P": sites}, min_size=50, max_size=10000)
        loose = predict_amplicons(
            {"P": sites}, min_size=50, max_size=10000, require_priming=False
        )
        assert len(loose) >= len(strict)


class TestOhnologDiscrimination:
    def test_conserved_primer_primes_on_the_paralog(self, ohnolog_pair, panel):
        target, _paralog, block = ohnolog_pair
        primer = target.seq[block[0] + 40 : block[0] + 60]
        report = assess_primer(primer, panel, primer_name="conserved")
        assert report.primes_on_panel
        assert report.best_identity_mismatches == 0
        assert "primes on" in report.verdict()

    def test_divergent_primer_is_discriminated(self, ohnolog_pair, panel):
        target, _paralog, block = ohnolog_pair
        # Scan outside the conserved block for a primer whose 3' end differs.
        found = None
        for i in range(0, block[0] - 30):
            report = assess_primer(target.seq[i : i + 20], panel, primer_name="d")
            if report.best_site and not report.primes_on_panel and report.discriminates_at_3p:
                found = report
                break
        assert found is not None
        assert "3' end" in found.verdict()

    def test_pair_in_the_conserved_block_predicts_a_paralog_product(self, ohnolog_pair, panel):
        target, _paralog, block = ohnolog_pair
        fwd = target.seq[block[0] + 5 : block[0] + 25]
        rev = revcomp(target.seq[block[1] - 25 : block[1] - 5])
        amplicons, reports = assess_pair(fwd, rev, panel, min_size=50, max_size=1000)
        assert amplicons
        assert all(r.primes_on_panel for r in reports.values())

    def test_pipeline_discards_cross_amplifying_pairs(self, ohnolog_pair, panel):
        target, _paralog, block = ohnolog_pair
        params = get_preset("genotyping-gel").replace(max_pairs=5)
        result = run_design(
            DesignRequest(template=target, params=params, panel=panel, name="t")
        )
        assert result.pairs
        for pair in result.pairs:
            assert pair.is_specific
            assert any("avoid-panel" in n for n in pair.notes)

    def test_allowing_paralog_amplicons_reports_them_instead(self, ohnolog_pair):
        from danioprime.design.regions import Interval, IntervalSet

        target, paralog, block = ohnolog_pair
        panel = ParalogPanel.from_records([paralog])
        params = get_preset("genotyping-gel").replace(
            max_pairs=5, reject_paralog_amplicon=False, min_product_size=120
        )
        # Force the design into the conserved block.
        result = run_design(
            DesignRequest(
                template=target,
                params=params,
                forward_search=IntervalSet([Interval(block[0], block[0] + 60)]),
                reverse_search=IntervalSet([Interval(block[1] - 60, block[1])]),
                panel=panel,
                name="t",
            )
        )
        assert result.pairs
        assert any(p.paralog_amplicons for p in result.pairs)

    def test_panel_strategy_switches_on_size(self, ohnolog_pair):
        _target, paralog, _block = ohnolog_pair
        small = ParalogPanel.from_records([paralog])
        assert small.resolved_strategy() == "exhaustive"
        forced = ParalogPanel.from_records([paralog], strategy="seeded")
        assert forced.resolved_strategy() == "seeded"

    def test_empty_panel_is_falsy(self):
        assert not ParalogPanel()


class TestExternalTools:
    def test_availability_probe_does_not_raise(self):
        assert set(available_aligners()) == {"blastn", "makeblastdb", "bowtie2", "bowtie2-build"}

    def test_md_tag_parsing(self):
        assert _md_mismatch_offsets("10A9", "20M", 20, "+") == [10]
        assert _md_mismatch_offsets("0T19", "20M", 20, "-") == [19]
        assert _md_mismatch_offsets("20", "20M", 20, "+") == []

    def test_gapped_and_deleted_alignments_are_rejected(self):
        assert _md_mismatch_offsets("10A9", "5M2I13M", 20, "+") is None
        assert _md_mismatch_offsets("10^AC9", "20M", 20, "+") is None

    def test_missing_tool_raises_a_helpful_error(self):
        from danioprime.errors import ExternalToolError
        from danioprime.specificity.external import blastn_sites, which

        if which("blastn"):
            pytest.skip("blastn is installed")
        with pytest.raises(ExternalToolError, match="not installed"):
            blastn_sites({"p": PRIMER}, "nonexistent_db")


class TestExternalParsers:
    """Parse canned aligner output, so the parsers are tested without the tools.

    The minus-strand cases matter most: ``blastn`` writes ``sseq`` in the
    query's alignment register while SAM stores ``SEQ`` in the reference's, and
    getting either convention backwards turns a perfect hit into a wall of
    mismatches.
    """

    FIELDS = "qseqid sseqid sstart send sstrand nident length qstart qend qseq sseq mismatch"

    @staticmethod
    def _blast_line(**kw) -> str:
        order = TestExternalParsers.FIELDS.split()
        return "\t".join(str(kw[f]) for f in order)

    def _parse(self, lines, primers=None, **kwargs):
        from danioprime.specificity.external import _parse_blast
        from danioprime.specificity.scan import DEFAULT_RULES
        from danioprime.thermo.nn import DEFAULT_CONDITIONS

        return _parse_blast(
            "\n".join(lines),
            primers or {"p": PRIMER},
            rules=kwargs.pop("rules", DEFAULT_RULES),
            conditions=DEFAULT_CONDITIONS,
            max_mismatches=kwargs.pop("max_mismatches", 6),
        )

    def test_perfect_plus_strand_hit(self):
        line = self._blast_line(
            qseqid="p", sseqid="chr5", sstart=1001, send=1020, sstrand="plus",
            nident=20, length=20, qstart=1, qend=20, qseq=PRIMER, sseq=PRIMER, mismatch=0,
        )
        sites = self._parse([line])["p"]
        assert len(sites) == 1
        site = sites[0]
        assert site.strand == "+"
        assert (site.start, site.end) == (1000, 1020)
        assert site.mismatches == 0
        assert site.can_prime

    def test_perfect_minus_strand_hit_is_not_reported_as_mismatched(self):
        # blastn reports sstart > send for a minus hit, and sseq equals the
        # primer when the match is perfect.
        line = self._blast_line(
            qseqid="p", sseqid="chr5", sstart=2020, send=2001, sstrand="minus",
            nident=20, length=20, qstart=1, qend=20, qseq=PRIMER, sseq=PRIMER, mismatch=0,
        )
        sites = self._parse([line])["p"]
        assert len(sites) == 1
        site = sites[0]
        assert site.strand == "-"
        assert (site.start, site.end) == (2000, 2020)
        assert site.mismatches == 0
        assert site.can_prime

    def test_three_prime_mismatch_is_located_correctly(self):
        subject = PRIMER[:-1] + ("G" if PRIMER[-1] != "G" else "T")
        line = self._blast_line(
            qseqid="p", sseqid="chr5", sstart=1, send=20, sstrand="plus",
            nident=19, length=20, qstart=1, qend=20, qseq=PRIMER, sseq=subject, mismatch=1,
        )
        site = self._parse([line])["p"][0]
        assert site.mismatch_offsets_from_3p == (0,)
        assert not site.can_prime
        assert "3' end" in site.note

    def test_five_prime_mismatch_still_primes(self):
        subject = ("G" if PRIMER[0] != "G" else "T") + PRIMER[1:]
        line = self._blast_line(
            qseqid="p", sseqid="chr5", sstart=1, send=20, sstrand="plus",
            nident=19, length=20, qstart=1, qend=20, qseq=PRIMER, sseq=subject, mismatch=1,
        )
        site = self._parse([line])["p"][0]
        assert site.mismatch_offsets_from_3p == (19,)
        assert site.can_prime

    def test_partial_alignments_are_discarded(self):
        # A hit covering only the primer's 5' half cannot be extended from.
        line = self._blast_line(
            qseqid="p", sseqid="chr5", sstart=1, send=10, sstrand="plus",
            nident=10, length=10, qstart=1, qend=10, qseq=PRIMER[:10],
            sseq=PRIMER[:10], mismatch=0,
        )
        assert self._parse([line])["p"] == []

    def test_unknown_query_and_blank_lines_are_ignored(self):
        line = self._blast_line(
            qseqid="other", sseqid="chr5", sstart=1, send=20, sstrand="plus",
            nident=20, length=20, qstart=1, qend=20, qseq=PRIMER, sseq=PRIMER, mismatch=0,
        )
        result = self._parse([line, "", "too\tfew\tfields"])
        assert result == {"p": []}

    def test_mismatch_budget_is_applied(self):
        subject = "A" * 8 + PRIMER[8:]
        line = self._blast_line(
            qseqid="p", sseqid="chr5", sstart=1, send=20, sstrand="plus",
            nident=14, length=20, qstart=1, qend=20, qseq=PRIMER, sseq=subject, mismatch=6,
        )
        assert self._parse([line], max_mismatches=2)["p"] == []
        assert self._parse([line], max_mismatches=8)["p"]

    def test_sam_plus_and_minus_offsets_are_mirrored(self):
        from danioprime.specificity.external import _parse_sam
        from danioprime.specificity.scan import DEFAULT_RULES
        from danioprime.thermo.nn import DEFAULT_CONDITIONS

        # MD "19A0" is a mismatch at reference-order position 19.
        plus = "\t".join(
            ["p", "0", "chr5", "1001", "42", "20M", "*", "0", "0", PRIMER, "I" * 20, "MD:Z:19A0"]
        )
        minus = "\t".join(
            ["p", "16", "chr5", "2001", "42", "20M", "*", "0", "0", PRIMER, "I" * 20, "MD:Z:19A0"]
        )
        sites = _parse_sam(
            "\n".join(["@HD\tVN:1.6", plus, minus]),
            {"p": PRIMER},
            rules=DEFAULT_RULES,
            conditions=DEFAULT_CONDITIONS,
            max_mismatches=6,
        )["p"]
        by_strand = {s.strand: s for s in sites}
        assert set(by_strand) == {"+", "-"}
        # SEQ is stored in reference orientation, so the same MD tag describes a
        # 5'-end mismatch on the plus strand and a 3'-end mismatch on the minus.
        assert by_strand["+"].mismatch_offsets_from_3p == (0,)
        assert by_strand["-"].mismatch_offsets_from_3p == (19,)
        assert not by_strand["+"].can_prime
        assert by_strand["-"].can_prime

    def test_sam_unmapped_and_gapped_records_are_skipped(self):
        from danioprime.specificity.external import _parse_sam
        from danioprime.specificity.scan import DEFAULT_RULES
        from danioprime.thermo.nn import DEFAULT_CONDITIONS

        unmapped = "\t".join(["p", "4", "*", "0", "0", "*", "*", "0", "0", PRIMER, "*"])
        gapped = "\t".join(
            ["p", "0", "chr5", "1", "42", "5M2I13M", "*", "0", "0", PRIMER, "I" * 20, "MD:Z:18"]
        )
        no_md = "\t".join(
            ["p", "0", "chr5", "1", "42", "20M", "*", "0", "0", PRIMER, "I" * 20]
        )
        sites = _parse_sam(
            "\n".join([unmapped, gapped, no_md]),
            {"p": PRIMER},
            rules=DEFAULT_RULES,
            conditions=DEFAULT_CONDITIONS,
            max_mismatches=6,
        )["p"]
        assert sites == []
