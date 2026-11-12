"""The one piece of real zebrafish sequence in this repository.

``examples/myo7aa_locus.fa`` is 1,000 bp of actual *Danio rerio* chromosome 18
from danRer11, and ``examples/myo7aa.gtf`` is the real RefSeq exon that sits in
it. Everything else in ``examples/`` is synthetic and says so.

Real data in a repository needs guarding differently from generated data. There
is no generator to re-run, so a corrupted file would simply sit there looking
plausible -- and the thing that makes it valuable, that it is genuinely this
locus, is exactly what a diff cannot tell you. These tests therefore check the
sequence against the annotation rather than against a stored copy of itself:

The gene is on the minus strand, so every intron's canonical GT..AG appears on
the plus strand as AC at the high end and CT at the low end. If a base were
dropped, a line reflowed, the window shifted, or the exon coordinates edited,
the dinucleotides would no longer land where the GTF says the boundaries are.
That is a check the data has to pass on its own terms, and synthetic sequence
pasted in behind the same header would fail it.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from danioprime.sequence import read_fasta

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
LOCUS = EXAMPLES / "myo7aa_locus.fa"
GTF = EXAMPLES / "myo7aa.gtf"

#: The window the sequence was taken from, 1-based inclusive.
WINDOW_START = 3405001
WINDOW_END = 3406000

#: UCSC refGene NM_152983, exon 24 of 48, 1-based inclusive.
EXON_START = 3405591
EXON_END = 3405767


@pytest.fixture(scope="module")
def record():
    [rec] = read_fasta(LOCUS)
    return rec


@pytest.fixture(scope="module")
def raw() -> str:
    """The sequence with its case intact; soft-masking is part of the data."""
    text = LOCUS.read_text()
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


class TestTheSequence:
    def test_it_is_the_window_it_claims_to_be(self, raw):
        assert len(raw) == WINDOW_END - WINDOW_START + 1 == 1000

    def test_it_is_dna(self, raw):
        assert set(raw) <= set("ACGTacgt"), "no Ns, no IUPAC codes, no stray characters"

    def test_the_soft_masking_survived(self, raw):
        """Lower case is UCSC's repeat masking, and it is the point.

        Zebrafish is repeat-dense; a primer landing in one is a real hazard
        that synthetic sequence cannot demonstrate. Upper-casing the file on
        the way through an editor would throw that away silently.
        """
        masked = sum(c.islower() for c in raw)
        assert masked > 0, "the repeat masking has been flattened"
        assert 0.2 < masked / len(raw) < 0.7, f"{masked / len(raw):.0%} masked looks wrong"

    def test_it_looks_like_zebrafish_rather_than_a_random_string(self, raw):
        """AT-rich, as this genome is. A uniform random sequence would be ~50%."""
        gc = sum(raw.upper().count(b) for b in "GC") / len(raw)
        assert 0.33 < gc < 0.47, f"GC {gc:.1%} is not in the zebrafish range"


class TestTheAnnotationAgreesWithTheSequence:
    """The check that cannot be faked by pasting in a different sequence."""

    @staticmethod
    def _at(raw: str, genomic_1based: int, length: int = 2) -> str:
        index = genomic_1based - WINDOW_START
        return raw[index : index + length].upper()

    def test_the_gtf_exon_maps_back_to_the_refseq_coordinates(self):
        [line] = [
            ln for ln in GTF.read_text().splitlines() if ln and not ln.startswith("#")
        ]
        fields = line.split("\t")
        assert fields[0] == "myo7aa_locus"
        assert fields[2] == "exon"
        assert fields[6] == "-", "myo7aa is on the minus strand"
        local_start, local_end = int(fields[3]), int(fields[4])
        assert local_start + WINDOW_START - 1 == EXON_START
        assert local_end + WINDOW_START - 1 == EXON_END
        assert 'transcript_id "NM_152983"' in line
        assert 'gene_name "myo7aa"' in line
        assert 'exon_number "24"' in line

    def test_the_five_prime_splice_site_is_where_the_annotation_says(self, raw):
        """A minus-strand intron's GT reads AC on the plus strand."""
        assert self._at(raw, EXON_START - 2) == "AC"

    def test_the_three_prime_splice_site_is_too(self, raw):
        """And its AG reads CT."""
        assert self._at(raw, EXON_END + 1) == "CT"

    def test_the_exon_is_the_length_refseq_says(self, raw):
        assert EXON_END - EXON_START + 1 == 177
        exon = raw[EXON_START - WINDOW_START : EXON_END - WINDOW_START + 1]
        assert len(exon) == 177

    def test_a_one_base_shift_would_be_caught(self, raw):
        """The splice-site check has to be sharp enough to be worth having."""
        shifted = raw[1:] + "A"
        assert shifted[EXON_START - 2 - WINDOW_START :][:2].upper() != "AC"


class TestTheProvenanceIsRecorded:
    """A real locus with no provenance is indistinguishable from a made-up one."""

    @pytest.mark.parametrize(
        "token",
        [
            "danRer11",                 # the assembly, which the variant guard reads
            "chr18:3405001-3406000",    # where it came from
            "NM_152983",                # RefSeq
            "ENSDARG00000099871",       # Ensembl
            "ZDB-GENE-020709-1",        # ZFIN
            "MINUS strand",             # why the splice sites read AC/CT
            "UCSC Genome Browser API",  # the source
            "2026-11-12",               # when
        ],
    )
    def test_the_header_records_it(self, record, token):
        assert token in record.description

    def test_it_does_not_claim_to_be_synthetic(self, record):
        assert "SYNTHETIC" not in record.description.upper()
        assert "real Danio rerio sequence" in record.description

    def test_every_other_example_locus_says_it_is_synthetic(self):
        """The inverse guard: nothing else may look real by omission."""
        for path in sorted(EXAMPLES.glob("*.fa")):
            if path.name.startswith("myo7aa"):
                continue
            [rec] = read_fasta(path)
            assert "SYNTHETIC" in rec.description.upper(), f"{path.name} does not say so"


class TestItWorksAsATemplate:
    def test_a_design_finds_pairs_flanking_the_real_exon(self):
        from danioprime.config import get_preset
        from danioprime.design.pipeline import DesignRequest, run_design
        from danioprime.design.regions import Interval, IntervalSet

        [rec] = read_fasta(LOCUS)
        midpoint = (EXON_START + EXON_END) // 2 - WINDOW_START
        result = run_design(
            DesignRequest(
                template=rec,
                params=get_preset("genotyping-gel"),
                targets=IntervalSet([Interval(midpoint, midpoint + 1)]),
            )
        )
        assert result.pairs, "no pairs on the real locus"
        best = result.pairs[0]
        assert best.forward.end <= midpoint <= best.reverse.start, "amplicon spans the exon"

    def test_the_cli_runs_on_it(self, capsys):
        from danioprime.cli import main

        assert main(["design", "--template", str(LOCUS), "--target", "679",
                     "--preset", "genotyping-gel", "--show", "1"]) == 0
        out = capsys.readouterr().out
        assert "myo7aa_locus" in out
        assert re.search(r"\d+ bp  penalty", out)
