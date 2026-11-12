"""Assay modes: genotyping geometry, qPCR protection, dCAPS digests, tails."""

from __future__ import annotations

from dataclasses import replace

import pytest
from helpers import make_seq

from danioprime.annotation.gtf import (
    Exon,
    Transcript,
    canonical_transcript,
    find_transcripts,
    guess_ohnolog_symbol,
    parse_annotation,
)
from danioprime.assays.dcaps import design_dcaps
from danioprime.assays.enzymes import (
    DEFAULT_ENZYMES,
    concretise,
    digest_fragments,
    find_sites,
    load_enzymes,
    site_matches,
)
from danioprime.assays.genotyping import (
    design_genotyping,
    design_junction_pcr,
    predict_allele_sizes,
)
from danioprime.assays.qpcr import design_qpcr, intron_sizes_from_transcript
from danioprime.assays.tails import (
    RIBOPROBE_SCHEMES,
    TAIL_SCHEMES,
    TAILS,
    apply_scheme,
    order_sheet,
    reaction_plan,
)
from danioprime.config import get_preset
from danioprime.design.regions import Interval, IntervalSet
from danioprime.errors import AnnotationError, DesignError
from danioprime.sequence import SeqRecord, revcomp


class TestGenotyping:
    """No cut site anywhere: see the module docstring in assays/genotyping.py.

    What is tested instead is the contract that replaced it -- a region the
    amplicon must span and no primer may sit in -- plus the readout check,
    which was always the sound half of the old geometry.
    """

    @pytest.fixture
    def design(self, locus, params):
        return design_genotyping(
            locus, params, targets=IntervalSet([Interval(955, 1045)])
        )

    def test_the_amplicon_spans_the_region(self, design):
        assert design.pairs
        for pair in design.pairs:
            assert pair.forward.start < 955
            assert pair.reverse.end > 1045

    def test_no_primer_sits_in_the_region(self, design):
        """A primer inside the lesion does not amplify the mutant allele."""
        for pair in design.pairs:
            assert pair.forward.end <= 955
            assert pair.reverse.start >= 1045

    def test_a_design_with_no_region_at_all_still_works(self, locus, params):
        """The right answer when the extent of the lesion is unknown."""
        design = design_genotyping(locus, params)
        assert design.pairs
        assert design.allele_sizes

    def test_a_region_outside_the_template_is_an_error(self, locus, params):
        with pytest.raises(DesignError, match="outside the template"):
            design_genotyping(
                locus, params, targets=IntervalSet([Interval(99999, 100050)])
            )

    def test_a_region_wider_than_the_product_window_is_diagnosed(self, locus, params):
        """Naming the two constraints that conflict, not "no pair was found".

        This is the trade-off the old point-anchored geometry hid: a cut
        coordinate always fitted inside any product-size window, so nobody had
        to notice that a readout chosen for short products cannot interrogate a
        wide region.
        """
        narrow = replace(params, max_product_size=200, opt_product_size=180,
                         min_product_size=100)
        with pytest.raises(DesignError, match="--max-product-size"):
            design_genotyping(
                locus, narrow, targets=IntervalSet([Interval(500, 900)])
            )

    def test_enumeration_is_aimed_at_the_region(self, locus, params):
        """Otherwise the candidate cap is spent at the 5' end of the template.

        With no windows, enumeration walks the template from the start, hits
        max_candidates_per_side, and keeps candidates that provably cannot pair
        with a target in the middle -- reporting no pairs on a template that
        has plenty.
        """
        starved = replace(params, max_candidates_per_side=60)
        design = design_genotyping(
            locus, starved, targets=IntervalSet([Interval(955, 1045)])
        )
        assert design.pairs, "the candidate budget was spent in the wrong place"
        reach = starved.max_product_size
        for candidate in design.result.forward_candidates:
            assert candidate.start >= 955 - reach

    def test_allele_sizes_are_predicted(self, design):
        rows = design.allele_sizes
        assert rows[0]["allele"] == "wild type"
        wt = rows[0]["size"]
        for row in rows[1:]:
            assert row["size"] == wt + row["indel"]

    def test_small_indels_are_flagged_as_unresolvable(self, design):
        assert any("will not resolve" in w for w in design.result.warnings)

    def test_allele_size_prediction_drops_impossible_alleles(self):
        rows = predict_allele_sizes(20, indels=(-50, 5))
        assert [r["indel"] for r in rows] == [0, 5]

    def test_nothing_named_after_a_cut_site_survives(self):
        """A guard against it coming back by accident."""
        import danioprime.assays.genotyping as module

        leaked = [n for n in dir(module) if "cut" in n.lower() or "protospacer" in n.lower()]
        assert leaked == [], f"cut-site machinery is back: {leaked}"


class TestJunctionPcr:
    def test_both_junctions_are_designed_across_the_boundary(self, locus, params):
        cassette = SeqRecord(name="cassette", seq=make_seq(900, 77))
        designs = design_junction_pcr(
            locus, 1000, cassette, params, genomic_flank=500, cassette_depth=400
        )
        assert set(designs) == {"5prime", "3prime"}
        five = designs["5prime"]
        assert five.pairs
        for pair in five.pairs:
            # boundary sits at 500 in the junction template
            assert pair.forward.end <= 500 < pair.reverse.start + pair.reverse.length

    def test_homology_arm_keeps_primers_outside_it(self, locus, params):
        cassette = SeqRecord(name="cassette", seq=make_seq(900, 78))
        designs = design_junction_pcr(
            locus, 1000, cassette, params, genomic_flank=500, homology_arm=60
        )
        for pair in designs["5prime"].pairs:
            assert pair.forward.end <= 500 - 60

    def test_insertion_point_is_validated(self, locus, params):
        cassette = SeqRecord(name="cassette", seq=make_seq(900, 79))
        with pytest.raises(DesignError, match="outside the genomic template"):
            design_junction_pcr(locus, 99999, cassette, params)


class TestQpcr:
    @pytest.fixture
    def transcript_model(self):
        return Transcript(
            transcript_id="T1",
            chrom="chr1",
            strand="+",
            exons=[Exon(1000, 1200, 1), Exon(5000, 5180, 2), Exon(9000, 9220, 3)],
            gene_name="demogenea",
        )

    def test_junction_offsets_and_intron_sizes(self, transcript_model):
        assert transcript_model.junction_offsets() == (200, 380)
        assert intron_sizes_from_transcript(transcript_model) == {200: 3800, 380: 3820}

    def test_minus_strand_intron_sizes(self):
        tx = Transcript(
            transcript_id="M",
            chrom="chr1",
            strand="-",
            exons=[Exon(9000, 9220), Exon(5000, 5180), Exon(1000, 1200)],
        )
        assert tx.sorted_exons()[0].start == 9000
        assert intron_sizes_from_transcript(tx) == {220: 3820, 400: 3800}

    def test_junction_spanning_is_enforced_when_reachable(self, transcript_model):
        seq = SeqRecord(name="T1", seq=make_seq(600, 123))
        params = get_preset("qpcr").replace(
            max_pairs=5, min_tm=55.0, max_tm=65.0, min_gc=30.0, max_gc=70.0, require_gc_clamp=0
        )
        design = design_qpcr(
            seq, params, transcript_model=transcript_model, require_junction_spanning=True
        )
        assert design.pairs
        for pair in design.pairs:
            assert (
                pair.forward.spans_junction is not None
                or pair.reverse.spans_junction is not None
            )
        assert all(row["protection"] == "strong" for row in design.assessment)

    def test_unreachable_requirement_warns_and_ranks_by_protection(self, transcript_model):
        seq = SeqRecord(name="T1", seq=make_seq(600, 123))
        design = design_qpcr(
            SeqRecord(name="T1", seq=seq.seq),
            get_preset("qpcr").replace(max_pairs=4),
            transcript_model=transcript_model,
            require_junction_spanning=True,
        )
        assert any("no primer could be placed" in w for w in design.result.warnings)
        strengths = [row["protection"] for row in design.assessment]
        rank = {"strong": 0, "moderate": 1, "weak": 2, "none": 3}
        assert [rank[s] for s in strengths] == sorted(rank[s] for s in strengths)

    def test_pairs_within_one_exon_report_no_protection(self, transcript_model):
        seq = SeqRecord(name="T1", seq=make_seq(600, 123))
        design = design_qpcr(
            seq, get_preset("qpcr").replace(max_pairs=6), junctions=(200, 380)
        )
        assert any(row["protection"] == "none" for row in design.assessment)

    def test_model_length_must_match_the_sequence(self, transcript_model):
        with pytest.raises(DesignError, match="must correspond"):
            design_qpcr(
                SeqRecord(name="T1", seq=make_seq(500, 1)),
                get_preset("qpcr"),
                transcript_model=transcript_model,
            )

    def test_junction_offsets_are_validated(self):
        with pytest.raises(DesignError, match="outside the transcript"):
            design_qpcr(
                SeqRecord(name="T", seq=make_seq(400, 2)),
                get_preset("qpcr"),
                junctions=(5000,),
            )

    def test_requirement_without_junctions_is_an_error(self):
        with pytest.raises(DesignError, match="needs junction offsets"):
            design_qpcr(
                SeqRecord(name="T", seq=make_seq(400, 2)),
                get_preset("qpcr"),
                require_junction_spanning=True,
            )


class TestEnzymes:
    def test_site_finding_and_digestion_are_consistent(self):
        ecori = DEFAULT_ENZYMES["EcoRI"]
        seq = "AAAGAATTCTTTGAATTCAAA"
        assert find_sites(ecori, seq) == [3, 12]
        frags = digest_fragments(ecori, seq)
        assert sum(frags) == len(seq)
        assert len(frags) == 3

    def test_degenerate_sites_match_iupac(self):
        styi = DEFAULT_ENZYMES["StyI"]
        assert find_sites(styi, "TTCCAAGGTT") == [2]
        assert find_sites(styi, "TTCCTTGGTT") == [2]
        assert find_sites(styi, "TTCCGCGGTT") == []

    def test_no_site_gives_a_single_fragment(self):
        assert digest_fragments(DEFAULT_ENZYMES["NotI"], "ACGT" * 20) == [80]

    def test_site_matches_rejects_length_mismatch(self):
        assert not site_matches("GAATTC", "GAATT")

    def test_concretise_picks_a_representative(self):
        assert concretise("CCWWGG") == "CCAAGG"
        assert concretise("GAATTC") == "GAATTC"

    def test_custom_table_loading(self, tmp_path):
        path = tmp_path / "e.tsv"
        path.write_text("# name site cut\nMyI\tGGCC\t2\nDegen\tGGWCC\t2\n")
        table = load_enzymes(str(path))
        assert table["MyI"].display == "GG^CC"
        assert table["Degen"].is_degenerate

    def test_bad_custom_table_is_reported(self, tmp_path):
        path = tmp_path / "bad.tsv"
        path.write_text("MyI\tGGXX\t2\n")
        with pytest.raises(AnnotationError, match="bad IUPAC"):
            load_enzymes(str(path))
        path.write_text("MyI\tGGCC\t99\n")
        with pytest.raises(AnnotationError, match="outside the site"):
            load_enzymes(str(path))


class TestDcaps:
    @pytest.fixture
    def setup(self):
        seq = make_seq(1200, 55)
        pos = 600
        ref = seq[pos]
        alt = next(b for b in "ACGT" if b != ref)
        params = get_preset("genotyping-gel").replace(
            min_product_size=150, opt_product_size=250, max_product_size=450
        )
        return SeqRecord(name="locus", seq=seq), pos, ref, alt, params

    def test_assays_are_produced_and_discriminate(self, setup):
        template, pos, _ref, alt, params = setup
        assays = design_dcaps(template, pos, alt, params, max_assays=5)
        assert assays
        for assay in assays:
            assert assay.discriminating
            assert sum(assay.reference_fragments) == assay.product_size
            assert sum(assay.variant_fragments) == assay.product_size

    def test_the_engineered_site_is_genuinely_allele_specific(self, setup):
        template, pos, _ref, alt, params = setup
        assay = design_dcaps(template, pos, alt, params, max_assays=1)[0]
        d = assay.dcaps_primer
        if d.orientation == "F":
            start, end = d.start, assay.partner.end
        else:
            start, end = assay.partner.start, d.end
        amplicon = list(template.seq[start:end])
        engineered = d.seq if d.orientation == "F" else revcomp(d.seq)
        off = d.start - start
        amplicon[off : off + len(engineered)] = list(engineered)
        ref_amplicon = "".join(amplicon)
        vi = pos - start
        var_amplicon = ref_amplicon[:vi] + alt + ref_amplicon[vi + 1 :]
        ref_sites = find_sites(d.enzyme, ref_amplicon)
        var_sites = find_sites(d.enzyme, var_amplicon)
        cut, uncut = (
            (var_sites, ref_sites) if d.cut_allele == "variant" else (ref_sites, var_sites)
        )
        assert len(cut) == len(uncut) + 1

    def test_the_variant_lies_outside_the_dcaps_primer(self, setup):
        template, pos, _ref, alt, params = setup
        for assay in design_dcaps(template, pos, alt, params, max_assays=5):
            d = assay.dcaps_primer
            assert not (d.start <= pos < d.end)

    def test_mismatch_is_kept_away_from_the_three_prime_terminus(self, setup):
        template, pos, _ref, alt, params = setup
        for assay in design_dcaps(template, pos, alt, params, max_assays=5):
            assert min(assay.dcaps_primer.mismatch_offsets_from_3p) >= 1

    def test_mismatch_budget_is_respected(self, setup):
        template, pos, _ref, alt, params = setup
        for assay in design_dcaps(
            template, pos, alt, params, max_assays=5, max_introduced_mismatches=1
        ):
            assert assay.dcaps_primer.mismatches <= 1

    def test_diagnostic_bands_exclude_shared_fragments(self, setup):
        template, pos, _ref, alt, params = setup
        for assay in design_dcaps(template, pos, alt, params, max_assays=5):
            ref_only, var_only = assay.diagnostic_bands
            assert ref_only or var_only
            shared = set(assay.reference_fragments) & set(assay.variant_fragments)
            assert not (set(ref_only) & set(var_only) & shared)

    def test_identical_base_is_rejected(self, setup):
        template, pos, ref, _alt, params = setup
        with pytest.raises(DesignError, match="nothing to discriminate"):
            design_dcaps(template, pos, ref, params)

    def test_invalid_base_is_rejected(self, setup):
        template, pos, _ref, _alt, params = setup
        with pytest.raises(DesignError, match="must be A/C/G/T"):
            design_dcaps(template, pos, "N", params)

    def test_position_is_validated(self, setup):
        template, _pos, _ref, alt, params = setup
        with pytest.raises(DesignError, match="outside template"):
            design_dcaps(template, 99999, alt, params)

    def test_tiny_enzyme_table_fails_with_guidance(self, setup):
        template, pos, _ref, alt, params = setup
        only_notI = {"NotI": DEFAULT_ENZYMES["NotI"]}
        with pytest.raises(DesignError, match="no enzyme in the table"):
            design_dcaps(
                template, pos, alt, params, enzymes=only_notI, max_introduced_mismatches=1
            )


class TestTails:
    @pytest.fixture
    def pairs(self, locus, params):
        from danioprime.design.pipeline import DesignRequest, run_design

        return run_design(
            DesignRequest(template=locus, params=params.replace(max_pairs=2), name="g")
        ).pairs

    def test_fluorescent_scheme_tails_both_primers(self, pairs):
        pair = apply_scheme(pairs[0], "fluorescent")
        assert pair.forward.ordered_seq.startswith("TGTAAAACGACGGCCAGT")
        assert pair.reverse.ordered_seq.startswith("GTTTCTT")
        assert pair.forward.tail_name == "M13F"

    def test_annealing_sequence_is_unchanged_by_a_tail(self, pairs):
        before = pairs[0].forward.seq
        apply_scheme(pairs[0], "fluorescent")
        assert pairs[0].forward.seq == before

    def test_companion_oligo_is_listed_once(self, pairs):
        items = order_sheet(pairs, scheme="fluorescent")
        companions = [i for i in items if i.role == "companion"]
        assert len(companions) == 1
        assert companions[0].sequence == "TGTAAAACGACGGCCAGT"

    def test_no_scheme_leaves_sequences_bare(self, pairs):
        items = order_sheet(pairs, scheme="none")
        assert all(not i.tail for i in items)
        assert all(i.role in {"forward", "reverse"} for i in items)

    def test_reaction_plan_explains_the_three_primer_system(self, pairs):
        plan = reaction_plan(pairs, scheme="fluorescent")
        assert any("Three-primer" in n for n in plan.notes)

    def test_antisense_scheme_puts_the_promoter_on_the_reverse_primer(self, pairs):
        plan = reaction_plan(pairs, scheme="riboprobe-t7")
        reverse = [o for o in plan.oligos if o.role == "reverse"]
        forward = [o for o in plan.oligos if o.role == "forward"]
        assert all(o.sequence.startswith(TAILS["T7"].seq) for o in reverse)
        assert not any(o.sequence.startswith(TAILS["T7"].seq) for o in forward)
        assert any("ANTISENSE" in n for n in plan.notes)

    def test_sense_scheme_puts_the_promoter_on_the_forward_primer(self, pairs):
        """Orientation is the thing people get backwards, so assert both ways."""
        plan = reaction_plan(pairs, scheme="riboprobe-t7-sense")
        forward = [o for o in plan.oligos if o.role == "forward"]
        reverse = [o for o in plan.oligos if o.role == "reverse"]
        assert all(o.sequence.startswith(TAILS["T7"].seq) for o in forward)
        assert not any(o.sequence.startswith(TAILS["T7"].seq) for o in reverse)
        assert any("SENSE" in n and "control" in n for n in plan.notes)

    @pytest.mark.parametrize("promoter", ["T7", "T3", "SP6"])
    def test_all_three_polymerases_have_an_antisense_scheme(self, pairs, promoter):
        plan = reaction_plan(pairs, scheme=f"riboprobe-{promoter.lower()}")
        reverse = [o for o in plan.oligos if o.role == "reverse"]
        assert reverse
        assert all(o.sequence.startswith(TAILS[promoter].seq) for o in reverse)

    @pytest.mark.parametrize("promoter", ["T7", "T3", "SP6"])
    def test_all_three_polymerases_have_a_sense_scheme(self, pairs, promoter):
        plan = reaction_plan(pairs, scheme=f"riboprobe-{promoter.lower()}-sense")
        forward = [o for o in plan.oligos if o.role == "forward"]
        assert forward
        assert all(o.sequence.startswith(TAILS[promoter].seq) for o in forward)

    def test_dual_scheme_uses_two_different_promoters(self, pairs):
        """One PCR product that templates both probes needs distinct promoters."""
        plan = reaction_plan(pairs, scheme="riboprobe-dual-t7-sp6")
        fwd = next(o for o in plan.oligos if o.role == "forward")
        rev = next(o for o in plan.oligos if o.role == "reverse")
        assert fwd.sequence.startswith(TAILS["SP6"].seq)
        assert rev.sequence.startswith(TAILS["T7"].seq)
        assert fwd.tail != rev.tail
        assert any("one PCR product" in n.lower() or "both probes" in n for n in plan.notes)

    def test_promoter_sequences_match_the_published_consensus(self):
        """Wrong by one base means wasted oligos and a failed transcription.

        Consensus sequences as published in Thermo Fisher's in vitro
        transcription guide (MAN0011982); in each the final G of the core is
        the +1 base of the transcript.
        """
        assert TAILS["T7"].seq == "TAATACGACTCACTATAGGG"
        assert TAILS["T3"].seq == "AATTAACCCTCACTAAAGGG"
        assert TAILS["SP6"].seq == "ATTTAGGTGACACTATAGAA"
        for name in ("T7", "T3", "SP6"):
            assert TAILS[name].is_promoter
            assert TAILS[name].source

    def test_every_registered_promoter_is_reachable_from_a_scheme(self):
        """T3 was registered as a tail but had no scheme, so no UI could offer it."""
        used = {
            tail
            for scheme in RIBOPROBE_SCHEMES
            for tail in TAIL_SCHEMES[scheme].values()
        }
        promoters = {n for n, t in TAILS.items() if t.is_promoter}
        assert promoters <= used, f"promoters with no scheme: {sorted(promoters - used)}"

    def test_riboprobe_notes_mention_cleanup_and_leader(self, pairs):
        plan = reaction_plan(pairs, scheme="riboprobe-t7")
        joined = " ".join(plan.notes)
        assert "5' of the promoter" in joined
        assert "purif" in joined or "clean" in joined

    def test_unknown_scheme_is_rejected(self, pairs):
        with pytest.raises(DesignError, match="unknown tail scheme"):
            apply_scheme(pairs[0], "nope")

    def test_every_scheme_is_described(self):
        from danioprime.assays.tails import SCHEME_DESCRIPTIONS

        assert set(TAIL_SCHEMES) == set(SCHEME_DESCRIPTIONS)


class TestAnnotation:
    GTF = (
        'chr5\tx\texon\t1001\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
        'gene_name "demogenea"; transcript_biotype "protein_coding";\n'
        'chr5\tx\texon\t2001\t2150\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
        'gene_name "demogenea"; transcript_biotype "protein_coding";\n'
        'chr9\tx\texon\t500\t700\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; '
        'gene_name "demogeneb";\n'
    )

    def test_gtf_parsing(self):
        txs = parse_annotation(self.GTF.splitlines())
        assert set(txs) == {"T1", "T2"}
        assert txs["T1"].spliced_length == 350
        assert txs["T1"].junction_offsets() == (200,)

    def test_gff3_attributes_are_understood(self):
        gff = "chr5\t.\texon\t10\t20\t.\t+\t.\tID=exon1;Parent=transcript:T9;Name=foo\n"
        assert "T9" in parse_annotation(gff.splitlines())

    def test_lookup_by_symbol_and_id(self):
        txs = parse_annotation(self.GTF.splitlines())
        assert [t.transcript_id for t in find_transcripts(txs, "demogenea")] == ["T1"]
        assert [t.transcript_id for t in find_transcripts(txs, "T2")] == ["T2"]

    def test_canonical_prefers_protein_coding_then_length(self):
        txs = parse_annotation(self.GTF.splitlines())
        assert canonical_transcript(txs.values()).transcript_id == "T1"

    def test_ohnolog_symbol_inference(self):
        assert guess_ohnolog_symbol("demogenea") == "demogeneb"
        assert guess_ohnolog_symbol("demogeneb") == "demogenea"
        assert guess_ohnolog_symbol("tp53") is None
        assert guess_ohnolog_symbol("a") is None

    def test_annotation_without_exons_is_an_error(self, tmp_path):
        from danioprime.annotation.gtf import read_annotation

        path = tmp_path / "x.gtf"
        path.write_text("chr1\tx\tgene\t1\t100\t.\t+\t.\tgene_id \"g\";\n")
        with pytest.raises(AnnotationError, match="no exon features"):
            read_annotation(path)


class TestTranscriptSequenceAssembly:
    """Building spliced and genomic sequence from a GTF plus an indexed FASTA."""

    @pytest.fixture
    def reference(self, tmp_path):
        from danioprime.sequence import write_fasta

        seq = make_seq(10000, 4242)
        path = tmp_path / "ref.fa"
        write_fasta([SeqRecord(name="chr1", seq=seq)], path, wrap=60)
        return path, seq

    def test_plus_strand_splicing_concatenates_exons_in_order(self, reference):
        from danioprime.sequence import IndexedFasta

        path, seq = reference
        tx = Transcript(
            transcript_id="T",
            chrom="chr1",
            strand="+",
            exons=[Exon(100, 200), Exon(500, 650), Exon(900, 1000)],
        )
        with IndexedFasta(path) as fa:
            spliced = tx.spliced_sequence(fa)
        assert spliced.seq == seq[100:200] + seq[500:650] + seq[900:1000]
        assert len(spliced.seq) == tx.spliced_length
        assert tx.junction_offsets() == (100, 250)

    def test_minus_strand_splicing_reverse_complements_and_reorders(self, reference):
        from danioprime.sequence import IndexedFasta

        path, seq = reference
        tx = Transcript(
            transcript_id="T",
            chrom="chr1",
            strand="-",
            exons=[Exon(100, 200), Exon(500, 650), Exon(900, 1000)],
        )
        with IndexedFasta(path) as fa:
            spliced = tx.spliced_sequence(fa)
        expected = revcomp(seq[900:1000]) + revcomp(seq[500:650]) + revcomp(seq[100:200])
        assert spliced.seq == expected

    def test_genomic_record_includes_flanks(self, reference):
        from danioprime.sequence import IndexedFasta

        path, seq = reference
        tx = Transcript(
            transcript_id="T", chrom="chr1", strand="+", exons=[Exon(500, 600)]
        )
        with IndexedFasta(path) as fa:
            locus = tx.genomic_record(fa, flank=200)
        assert locus.seq == seq[300:800]

    def test_exon_junction_windows(self, reference):
        from danioprime.annotation.gtf import iter_exon_junction_windows

        tx = Transcript(
            transcript_id="T",
            chrom="chr1",
            strand="+",
            exons=[Exon(0, 100), Exon(200, 300)],
        )
        windows = list(iter_exon_junction_windows(tx, flank=30))
        assert windows == [(100, 70, 130)]
