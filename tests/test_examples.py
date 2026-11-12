"""The shipped example dataset must actually demonstrate what the README claims.

These tests are the regression guard on the headline feature: without the
ohnolog panel the design lands on sequence shared with the paralog, and with it
the design moves. If that stops being true the example is no longer a
demonstration of anything.
"""

from __future__ import annotations

import pytest

from danioprime.assays.genotyping import design_genotyping, locate_protospacer
from danioprime.assays.qpcr import design_qpcr
from danioprime.config import get_preset
from danioprime.sequence import read_single_fasta
from danioprime.specificity.paralog import ParalogPanel, assess_primer

PROTOSPACER = "AAGCCCTGTTATGTGCAACG"


@pytest.fixture
def locus(example_dir):
    return read_single_fasta(example_dir / "demogenea_locus.fa")


@pytest.fixture
def ohnolog_panel(example_dir):
    return ParalogPanel.from_fasta(str(example_dir / "demogeneb_locus.fa"))


class TestExampleData:
    def test_protospacer_is_present_with_a_valid_pam(self, locus):
        sites = locate_protospacer(locus.seq, PROTOSPACER)
        assert len(sites) == 1
        assert sites[0].pam.endswith("GG")
        assert sites[0].position == 1767

    def test_gtf_matches_the_transcript(self, example_dir):
        from danioprime.annotation.gtf import read_annotation

        txs = read_annotation(example_dir / "demogenea.gtf")
        tx = txs["demogenea-201"]
        transcript = read_single_fasta(example_dir / "demogenea_transcript.fa")
        assert tx.spliced_length == len(transcript.seq)
        assert tx.junction_offsets() == (210, 540)

    def test_the_ohnolog_shares_a_conserved_block(self, locus, ohnolog_panel):
        # A primer inside the conserved block must be unable to discriminate.
        conserved_primer = locus.seq[1830:1850]
        report = assess_primer(conserved_primer, ohnolog_panel, primer_name="conserved")
        assert report.primes_on_panel

    def test_unrelated_sequence_has_no_panel_site(self, locus, ohnolog_panel):
        # Upstream of exon 1 the two loci are unrelated. Exon 1 itself is only
        # ~6% divergent, so individual 20-mers there can still match by chance;
        # that is the point of checking rather than assuming.
        report = assess_primer(locus.seq[100:120], ohnolog_panel, primer_name="unique")
        assert not report.primes_on_panel
        # Chance alignments at 5-6 mismatches are expected in a few kilobases of
        # unrelated sequence; what matters is that none of them could extend.
        assert all(not site.can_prime for site in report.all_sites)
        assert report.best_site is None or report.best_site.mismatches >= 4


class TestExampleWorkflow:
    def test_ohnolog_panel_changes_the_genotyping_design(self, locus, ohnolog_panel):
        params = get_preset("genotyping-gel").replace(max_pairs=5)
        site = locate_protospacer(locus.seq, PROTOSPACER)[0]

        naive = design_genotyping(locus, site, params)
        guarded = design_genotyping(locus, site, params, panel=ohnolog_panel)

        assert naive.pairs and guarded.pairs

        # The guarantee the panel buys: no reported pair is predicted to
        # produce a product on the ohnolog. A single primer with a panel
        # binding site is tolerated, because one primer gives linear
        # extension rather than a band.
        for pair in guarded.pairs:
            assert not pair.paralog_amplicons
            assert pair.is_specific

        # The naive design must contain at least one primer that does prime on
        # the ohnolog, or this example is not exercising the feature at all.
        naive_hits = [
            primer
            for pair in naive.pairs
            for primer in (pair.forward, pair.reverse)
            if assess_primer(primer.seq, ohnolog_panel).primes_on_panel
        ]
        assert naive_hits

        # Ranking must also have moved: the panel-aware run should not simply
        # return the naive top pair unchanged.
        assert (naive.pairs[0].forward.seq, naive.pairs[0].reverse.seq) != (
            guarded.pairs[0].forward.seq,
            guarded.pairs[0].reverse.seq,
        )

    def test_strict_specificity_removes_panel_binding_primers(self, locus, ohnolog_panel):
        params = get_preset("genotyping-gel").replace(max_pairs=5, reject_paralog_priming=True)
        site = locate_protospacer(locus.seq, PROTOSPACER)[0]
        strict = design_genotyping(locus, site, params, panel=ohnolog_panel)
        assert strict.pairs
        for pair in strict.pairs:
            for primer in (pair.forward, pair.reverse):
                assert not assess_primer(primer.seq, ohnolog_panel).primes_on_panel

    def test_qpcr_on_the_example_transcript_finds_protected_pairs(self, example_dir):
        transcript = read_single_fasta(example_dir / "demogenea_transcript.fa")
        ohnolog = ParalogPanel.from_fasta(str(example_dir / "demogeneb_transcript.fa"))
        design = design_qpcr(
            transcript,
            get_preset("qpcr").replace(max_pairs=5),
            junctions=(210, 540),
            require_junction_spanning=True,
            panel=ohnolog,
        )
        assert design.pairs
        assert any(row["protection"] == "strong" for row in design.assessment)

    def test_strain_variants_mask_real_positions(self, example_dir, locus):
        from danioprime.variants.vcf import load_mask

        mask, variants = load_mask(
            example_dir / "strain_variants.vcf",
            chrom="demogenea_locus",
            samples=["TU"],
        )
        assert variants
        assert mask.total_length == len(variants)
        for variant in variants:
            assert locus.seq[variant.pos] == variant.ref

    def test_heterozygous_only_recovers_fixed_positions(self, example_dir):
        from danioprime.variants.vcf import load_mask

        all_tu, _ = load_mask(
            example_dir / "strain_variants.vcf", chrom="demogenea_locus", samples=["TU"]
        )
        het_tu, _ = load_mask(
            example_dir / "strain_variants.vcf",
            chrom="demogenea_locus",
            samples=["TU"],
            heterozygous_only=True,
        )
        assert het_tu.total_length < all_tu.total_length


class TestInSituExample:
    """The dataset must demonstrate the in situ point, or the docs are bluffing.

    Ohnolog coding sequence is usually too conserved for any probe to
    discriminate, while the 3' UTR diverges fast enough that a specific probe
    can be placed there. If that contrast ever stops holding in the example
    data, the README's in situ section becomes a claim with nothing behind it.
    """

    @pytest.fixture
    def transcripts(self, example_dir):
        return (
            read_single_fasta(example_dir / "demogenea_transcript.fa"),
            ParalogPanel.from_fasta(str(example_dir / "demogeneb_transcript.fa")),
        )

    #: Coding sequence ends here in the example transcript.
    CDS_END = 800

    def test_a_probe_over_the_coding_region_cross_hybridises(self, transcripts):
        from danioprime.assays.insitu import assess_probe

        tx, panel = transcripts
        result = assess_probe(tx.seq[100 : self.CDS_END - 100], panel)
        assert result.will_cross_hybridise
        assert "probe body" in result.verdict()

    def test_a_probe_over_the_three_prime_utr_is_clean(self, transcripts):
        from danioprime.assays.insitu import assess_probe

        tx, panel = transcripts
        assert len(tx.seq) > self.CDS_END + 200, "the example transcript has no 3' UTR"
        result = assess_probe(tx.seq[self.CDS_END + 20 :], panel)
        assert not result.will_cross_hybridise

    def test_the_demo_design_reproduces_both_outcomes(self, transcripts):
        from danioprime.assays.insitu import design_insitu
        from danioprime.design.regions import Interval, IntervalSet

        tx, panel = transcripts
        params = get_preset("insitu-probe").replace(
            max_pairs=3, min_product_size=250, opt_product_size=290, max_product_size=320
        )
        over_cds = design_insitu(
            tx, params, scheme="riboprobe-t3", panel=panel,
            targets=IntervalSet([Interval(300, 500)]),
        )
        over_utr = design_insitu(
            tx, params, scheme="riboprobe-t3", panel=panel,
            targets=IntervalSet([Interval(849, 1049)]),
        )
        assert over_cds.pairs and over_utr.pairs
        assert over_cds.cross_hyb[over_cds.pairs[0].name].will_cross_hybridise
        assert not over_utr.cross_hyb[over_utr.pairs[0].name].will_cross_hybridise

    def test_the_antisense_promoter_is_on_the_reverse_primer(self, transcripts):
        from danioprime.assays.insitu import design_insitu

        tx, panel = transcripts
        params = get_preset("insitu-probe").replace(
            max_pairs=1, min_product_size=250, opt_product_size=290, max_product_size=320
        )
        design = design_insitu(tx, params, scheme="riboprobe-t3", panel=panel)
        rev = next(o for o in design.antisense.oligos if o.role == "reverse")
        assert rev.sequence.startswith("AATTAACCCTCACTAAAGGG")
