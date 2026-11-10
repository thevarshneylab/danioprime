"""Sequence primitives, FASTA I/O and indexed random access."""

from __future__ import annotations

import pytest

from danioprime.errors import SequenceError
from danioprime.sequence import (
    IndexedFasta,
    SeqRecord,
    gc_percent,
    iupac_matches,
    load_sequence,
    max_dinucleotide_repeat,
    max_homopolymer,
    normalise,
    parse_fasta,
    read_fasta,
    revcomp,
    soft_mask_runs,
    write_fasta,
)


class TestBasics:
    def test_revcomp_handles_iupac_and_case(self):
        assert revcomp("ACGT") == "ACGT"
        assert revcomp("AAAA") == "TTTT"
        assert revcomp("ACGTRYN") == "NRYACGT"
        assert revcomp("acgt") == "acgt"

    def test_revcomp_is_an_involution(self):
        for seq in ("ACGTTGCATGCC", "N", "RYSWKM"):
            assert revcomp(revcomp(seq)) == seq

    def test_gc_percent_ignores_ambiguous_bases(self):
        assert gc_percent("GGCC") == 100.0
        assert gc_percent("ATAT") == 0.0
        assert gc_percent("GGCCATAT") == 50.0
        assert gc_percent("GGCCNNNN") == 100.0
        assert gc_percent("") == 0.0

    def test_homopolymer_and_dinucleotide_runs(self):
        assert max_homopolymer("AACCCCGT") == 4
        assert max_homopolymer("") == 0
        assert max_dinucleotide_repeat("GGATATATATCC") == 4
        assert max_dinucleotide_repeat("ACGT") == 1

    def test_normalise_rejects_junk_and_converts_u(self):
        assert normalise("ac gt\n") == "ACGT"
        assert normalise("ACGU") == "ACGT"
        with pytest.raises(SequenceError, match="non-IUPAC"):
            normalise("ACGTX")
        with pytest.raises(SequenceError, match="ambiguity"):
            normalise("ACGTN", allow_ambiguous=False)

    def test_iupac_matching(self):
        assert iupac_matches("N", "A")
        assert iupac_matches("R", "G")
        assert not iupac_matches("R", "C")
        with pytest.raises(SequenceError):
            iupac_matches("X", "A")

    def test_softmask_runs(self):
        assert soft_mask_runs("ACGTacgtACGT") == [(4, 8)]
        assert soft_mask_runs("acgtACGTacgt") == [(0, 4), (8, 12)]
        assert soft_mask_runs("ACGT") == []
        assert soft_mask_runs("acgt", min_length=5) == []


class TestFasta:
    def test_parse_multiple_records(self):
        text = [">a desc here", "ACGT", "ACGT", ">b", "TTTT"]
        records = list(parse_fasta(text))
        assert [r.name for r in records] == ["a", "b"]
        assert records[0].seq == "ACGTACGT"
        assert records[0].description == "desc here"

    def test_softmask_preserved_but_sequence_uppercased(self):
        records = list(parse_fasta([">a", "acgtACGT"]))
        assert records[0].seq == "ACGTACGT"
        assert records[0].softmask == [(0, 4)]

    def test_sequence_before_header_is_an_error(self):
        with pytest.raises(SequenceError, match="before any"):
            list(parse_fasta(["ACGT"]))

    def test_round_trip(self, tmp_path):
        records = [SeqRecord(name="x", seq="ACGT" * 30, description="d")]
        path = tmp_path / "x.fa"
        write_fasta(records, path, wrap=25)
        back = read_fasta(path)
        assert back[0].seq == records[0].seq
        assert back[0].description == "d"

    def test_subseq_bounds_are_checked(self):
        rec = SeqRecord(name="x", seq="ACGTACGT")
        assert rec.subseq(2, 5) == "GTA"
        with pytest.raises(SequenceError, match="out of range"):
            rec.subseq(0, 99)

    def test_load_sequence_accepts_path_or_literal(self, tmp_path):
        path = tmp_path / "y.fa"
        write_fasta([SeqRecord(name="y", seq="ACGTACGT")], path)
        assert load_sequence(path).name == "y"
        assert load_sequence("ACGTACGT").seq == "ACGTACGT"
        with pytest.raises(SequenceError, match="neither an existing file"):
            load_sequence("not a sequence!")

    def test_empty_file_is_an_error(self, tmp_path):
        path = tmp_path / "empty.fa"
        path.write_text("")
        with pytest.raises(SequenceError, match="no FASTA records"):
            read_fasta(path)


class TestIndexedFasta:
    @pytest.fixture
    def reference(self, tmp_path):
        import random

        rng = random.Random(5)
        records = [
            SeqRecord(name=f"chr{i}", seq="".join(rng.choice("ACGT") for _ in range(3000)))
            for i in (1, 2)
        ]
        path = tmp_path / "ref.fa"
        write_fasta(records, path, wrap=60)
        return path, records

    def test_index_is_written_and_slices_match(self, reference):
        path, records = reference
        with IndexedFasta(path) as fa:
            assert fa.contigs == ["chr1", "chr2"]
            assert fa.length("chr1") == 3000
            for start, end in [(0, 10), (59, 61), (119, 125), (2990, 3000), (0, 3000)]:
                assert fa.fetch("chr1", start, end) == records[0].seq[start:end]
            assert fa.fetch("chr2", 100, 110) == records[1].seq[100:110]
        assert path.with_name(path.name + ".fai").exists()

    def test_out_of_range_is_clamped(self, reference):
        path, records = reference
        with IndexedFasta(path) as fa:
            assert fa.fetch("chr1", 2995, 4000) == records[0].seq[2995:]
            assert fa.fetch("chr1", 10, 5) == ""

    def test_unknown_contig_names_the_alternatives(self, reference):
        path, _ = reference
        with IndexedFasta(path) as fa:
            with pytest.raises(SequenceError, match="chr1"):
                fa.fetch("chrZ", 0, 10)

    def test_record_helper_labels_coordinates(self, reference):
        path, _ = reference
        with IndexedFasta(path) as fa:
            rec = fa.record("chr1", 99, 200)
            assert rec.name == "chr1:100-200"
            assert len(rec.seq) == 101

    def test_missing_file_is_reported(self, tmp_path):
        with pytest.raises(SequenceError, match="not found"):
            IndexedFasta(tmp_path / "nope.fa")

    def test_gzip_is_refused_with_an_explanation(self, tmp_path):
        path = tmp_path / "r.fa.gz"
        path.write_bytes(b"")
        with pytest.raises(SequenceError, match="gzipped"):
            IndexedFasta(path)
