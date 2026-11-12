"""Fetching a real locus.

The network is the one thing these tests cannot exercise, so it is pushed
behind a two-method transport and everything else is driven from recorded
response shapes. What is actually being tested: that the URLs are built to each
service's conventions (including the coordinate systems, which differ), that a
response is validated before it can become a FASTA file, and that a dry run
touches nothing.

The transport itself is tested separately against a real local HTTP server, so
the only untested link in the chain is the remote service's own behaviour.
"""

from __future__ import annotations

import json
import re

import pytest
from rangeserver import RangeServer

from danioprime.cli import main
from danioprime.errors import AnnotationError, SequenceError
from danioprime.genome.remote import (
    EnsemblRest,
    FetchPlan,
    UcscApi,
    UrllibTransport,
    _add_chr,
    _strip_chr,
    fetch_locus,
    fetch_transcript,
    resolve_provider,
)
from danioprime.sequence import read_fasta

# A real myo7aa lookup is shaped like this. The sequence is not myo7aa -- there
# is no real zebrafish sequence anywhere in this repository -- but every field
# the code reads is in the place Ensembl puts it.
MYO7AA_LOOKUP = {
    "id": "ENSDARG00000099871",
    "display_name": "myo7aa",
    "description": "myosin VIIAa [Source:ZFIN;Acc:ZDB-GENE-020709-1]",
    "species": "danio_rerio",
    "assembly_name": "GRCz11",
    "seq_region_name": "18",
    "start": 1000,
    "end": 1199,
    "strand": -1,
    "biotype": "protein_coding",
    "object_type": "Gene",
    "Transcript": [
        {"id": "ENSDART00000172076", "display_name": "myo7aa-202", "is_canonical": 0},
        {"id": "ENSDART00000148103", "display_name": "myo7aa-201", "is_canonical": 1},
    ],
}


class FakeTransport:
    """Answers by substring match, and records what it was asked for.

    A route may be a string or a callable taking the URL. The callable form is
    what ``region_fasta`` uses to answer with exactly the span that was asked
    for, so a test cannot accidentally pass because the fixture and the request
    disagreed about length -- the thing the code under test is checking.
    """

    def __init__(self, routes):
        self.routes = routes
        self.calls: list[str] = []

    def get(self, url: str, *, timeout: float = 30.0) -> str:
        self.calls.append(url)
        for fragment, body in self.routes.items():
            if fragment in url:
                return body(url) if callable(body) else body
        raise AssertionError(f"unexpected request: {url}")


def fasta(seq: str, name: str = "chromosome:GRCz11:18:1000:1199:1") -> str:
    lines = [seq[i : i + 60] for i in range(0, len(seq), 60)]
    return "\n".join([f">{name}", *lines, ""])


def region_fasta(url: str) -> str:
    """A FASTA of exactly the length the Ensembl region URL asks for."""
    chrom, start, end = re.search(r"/([^/:]+):(\d+)\.\.(\d+):", url).groups()
    span = int(end) - int(start) + 1
    return fasta(("ACGT" * (span // 4 + 1))[:span], name=f"{chrom}:{start}:{end}")


@pytest.fixture
def ensembl_transport():
    return FakeTransport(
        {
            "/lookup/symbol/": json.dumps(MYO7AA_LOOKUP),
            "/sequence/region/": region_fasta,
            "/sequence/id/": fasta("TTGCA" * 400, name="ENSDART00000148103.1 cdna"),
        }
    )


class TestUrlConventions:
    """The two services disagree about prefixes, separators and base.

    Each disagreement is a silent wrong answer rather than an error, so each
    one is pinned.
    """

    def test_ensembl_region_is_one_based_inclusive_with_dots(self):
        url = EnsemblRest().region_url("18", 1000, 1199, 1)
        assert "/sequence/region/danio_rerio/18:1000..1199:1" in url
        assert "content-type=text/x-fasta" in url

    def test_ucsc_region_is_zero_based_half_open(self):
        """The off-by-one that would shift every coordinate downstream."""
        url = UcscApi(genome="danRer11").region_url("18", 1000, 1199, 1)
        assert "start=999" in url and "end=1199" in url
        assert "genome=danRer11" in url

    def test_the_two_services_describe_the_same_200_bases(self):
        start, end = 1000, 1199
        ensembl = EnsemblRest().region_url("18", start, end)
        ucsc = UcscApi().region_url("chr18", start, end)
        assert f"{start}..{end}" in ensembl
        assert f"start={start - 1}" in ucsc and f"end={end}" in ucsc
        assert end - start + 1 == 200

    def test_ensembl_drops_the_chr_prefix(self):
        assert "/18:" in EnsemblRest().region_url("chr18", 1, 2)

    def test_ucsc_adds_it(self):
        assert "chrom=chr18" in UcscApi().region_url("18", 1, 2)

    @pytest.mark.parametrize("name", ["chrM", "M", "MT", "chrm"])
    def test_the_mitochondrion_is_renamed_not_just_reprefixed(self, name):
        assert _strip_chr(name) == "MT"
        assert _add_chr(name) == "chrM"

    @pytest.mark.parametrize("name", ["KZ115986.1", "CHR_ALT", "18", "chr18"])
    def test_a_scaffold_name_survives_the_round_trip(self, name):
        """Stripping three characters from anything that starts with 'chr'
        would mangle a scaffold called CHR_ALT."""
        assert _strip_chr(_add_chr(name)) == _strip_chr(name)

    def test_minus_strand_is_requested_as_minus_one(self):
        assert EnsemblRest().region_url("18", 1, 2, -1).split("?")[0].endswith(":-1")

    def test_transcript_url_asks_for_cdna_by_default(self):
        url = EnsemblRest().transcript_url("ENSDART00000148103")
        assert "/sequence/id/ENSDART00000148103" in url and "type=cdna" in url


class TestFetchBySymbol:
    def test_coordinates_and_provenance_come_from_the_server(self, ensembl_transport):
        locus = fetch_locus("myo7aa", flank=0, transport=ensembl_transport)
        assert locus.gene == "myo7aa"
        assert locus.gene_id == "ENSDARG00000099871"
        assert locus.assembly == "GRCz11", "never assumed; read from the response"
        assert locus.chrom == "18"
        assert (locus.start, locus.end) == (1000, 1199)
        assert locus.description == "myosin VIIAa"
        assert locus.retrieved, "a template with no date cannot be retraced"

    def test_flank_widens_the_request_on_both_sides(self, ensembl_transport):
        fetch_locus("myo7aa", flank=400, transport=ensembl_transport)
        assert "18:600..1599:1" in ensembl_transport.calls[-1]

    def test_flank_is_clamped_at_the_start_of_the_contig(self):
        transport = FakeTransport(
            {
                "/lookup/symbol/": json.dumps({**MYO7AA_LOOKUP, "start": 50, "end": 249}),
                "/sequence/region/": region_fasta,
            }
        )
        locus = fetch_locus("myo7aa", flank=2000, transport=transport)
        assert locus.start == 1, "position 0 does not exist in a 1-based assembly"
        assert "18:1..2249:1" in transport.calls[-1]

    def test_a_minus_strand_gene_is_still_fetched_in_assembly_orientation(
        self, ensembl_transport
    ):
        """myo7aa is on the minus strand in this fixture.

        Reverse-complementing the template here would renumber every
        coordinate the run then reports, for no benefit: the design code picks
        primer orientation itself.
        """
        assert MYO7AA_LOOKUP["strand"] == -1
        locus = fetch_locus("myo7aa", transport=ensembl_transport)
        assert locus.strand == 1
        assert ensembl_transport.calls[-1].split("?")[0].endswith(":1")

    def test_an_unknown_symbol_says_what_to_try(self):
        transport = FakeTransport({"/lookup/symbol/": json.dumps({"error": "not found"})})
        with pytest.raises(AnnotationError, match="ohnolog"):
            fetch_locus("myo7a", transport=transport)

    def test_a_non_json_lookup_is_an_error_not_a_crash(self):
        transport = FakeTransport({"/lookup/symbol/": "<html>503</html>"})
        with pytest.raises(AnnotationError, match="did not return JSON"):
            fetch_locus("myo7aa", transport=transport)


class TestFetchByRegion:
    def test_no_lookup_is_made(self):
        transport = FakeTransport({"/sequence/region/": fasta("ACGT" * 50)})
        locus = fetch_locus("ignored", chrom="18", start=1, end=200, transport=transport)
        assert len(transport.calls) == 1
        assert locus.gene == ""
        assert locus.name == "18_1_200"

    def test_ucsc_returns_json_with_a_dna_field(self):
        transport = FakeTransport({"getData/sequence": json.dumps({"dna": "ACGT" * 50})})
        locus = fetch_locus(
            "x", source="ucsc", chrom="chr18", start=1, end=200, transport=transport
        )
        assert locus.assembly == "danRer11"
        assert locus.seq.startswith("ACGT")

    def test_ucsc_minus_strand_is_reverse_complemented_locally(self):
        """UCSC has no strand parameter, so the flip has to happen here."""
        transport = FakeTransport({"getData/sequence": json.dumps({"dna": "AAAACCCC"})})
        locus = fetch_locus(
            "x", source="ucsc", chrom="chr18", start=1, end=8, strand=-1, transport=transport
        )
        assert locus.seq == "GGGGTTTT"

    def test_an_end_before_the_start_is_refused(self):
        with pytest.raises(AnnotationError, match="before start"):
            fetch_locus("x", chrom="18", start=500, end=100, transport=FakeTransport({}))

    def test_a_half_specified_region_is_refused(self):
        with pytest.raises(AnnotationError, match="both a start and an end"):
            fetch_locus("x", chrom="18", start=5, transport=FakeTransport({}))


class TestValidation:
    """The expensive failure is a file that looks fine and is the wrong DNA."""

    def test_a_short_response_is_refused(self):
        transport = FakeTransport({"/sequence/region/": fasta("ACGT")})
        with pytest.raises(AnnotationError, match="but 200 bp was requested"):
            fetch_locus("x", chrom="18", start=1, end=200, transport=transport)

    def test_an_error_page_is_not_mistaken_for_sequence(self):
        transport = FakeTransport({"/sequence/region/": "<html><body>Bad Gateway</body></html>"})
        with pytest.raises(SequenceError, match="not DNA"):
            fetch_locus("x", chrom="18", start=1, end=200, transport=transport)

    def test_an_empty_response_is_refused(self):
        transport = FakeTransport({"/sequence/region/": ">header only\n"})
        with pytest.raises(AnnotationError, match="no sequence"):
            fetch_locus("x", chrom="18", start=1, end=200, transport=transport)

    def test_assembly_gaps_and_soft_masking_are_allowed(self):
        """N is real and lower case is soft-masked repeat; neither is an error."""
        transport = FakeTransport({"/sequence/region/": fasta("acgtN" * 40)})
        locus = fetch_locus("x", chrom="18", start=1, end=200, transport=transport)
        assert len(locus.seq) == 200
        assert locus.to_record().seq.isupper(), "the record is upper case"

    def test_whitespace_in_the_body_is_not_counted_as_sequence(self):
        transport = FakeTransport({"/sequence/region/": fasta("ACGT" * 50)})
        locus = fetch_locus("x", chrom="18", start=1, end=200, transport=transport)
        assert len(locus.seq) == 200


class TestTranscript:
    def test_the_canonical_transcript_is_chosen(self, ensembl_transport):
        locus = fetch_transcript("myo7aa", transport=ensembl_transport)
        assert "ENSDART00000148103" in ensembl_transport.calls[-1], "the is_canonical one"
        assert locus.is_transcript
        assert len(locus.seq) == 2000

    def test_an_explicit_transcript_id_skips_the_lookup(self, ensembl_transport):
        fetch_transcript("ENSDART00000148103", transport=ensembl_transport)
        assert len(ensembl_transport.calls) == 1
        assert "/sequence/id/" in ensembl_transport.calls[0]

    def test_a_gene_with_no_transcripts_is_an_error(self):
        transport = FakeTransport(
            {"/lookup/symbol/": json.dumps({**MYO7AA_LOOKUP, "Transcript": []})}
        )
        with pytest.raises(AnnotationError, match="no transcripts"):
            fetch_transcript("myo7aa", transport=transport)

    def test_the_header_does_not_claim_genomic_coordinates(self, ensembl_transport):
        """A spliced sequence has none, and implying it does would be a lie."""
        header = fetch_transcript("myo7aa", transport=ensembl_transport).header()
        assert "spliced transcript" in header
        assert ":1-" not in header

    def test_no_length_check_applies_to_a_transcript(self, ensembl_transport):
        """Its length is what the exons add up to, which is not known in advance."""
        locus = fetch_transcript("myo7aa", transport=ensembl_transport)
        assert locus.end == len(locus.seq)


class TestDryRun:
    def test_nothing_is_requested_for_a_symbol(self):
        transport = FakeTransport({})  # any call raises
        plan = FetchPlan()
        locus = fetch_locus("myo7aa", transport=transport, plan=plan)
        assert transport.calls == []
        assert len(plan.urls) == 1 and "/lookup/symbol/danio_rerio/myo7aa" in plan.urls[0]
        assert locus.seq == ""

    def test_nothing_is_requested_for_a_region(self):
        transport = FakeTransport({})
        plan = FetchPlan()
        fetch_locus("x", chrom="18", start=1, end=200, transport=transport, plan=plan)
        assert transport.calls == []
        assert "/sequence/region/" in plan.urls[0]

    def test_nothing_is_requested_for_a_transcript(self):
        transport = FakeTransport({})
        plan = FetchPlan()
        fetch_transcript("ENSDART00000148103", transport=transport, plan=plan)
        assert transport.calls == []


class TestHeaderAndRecord:
    def test_the_header_carries_everything_needed_to_retrace_the_fetch(
        self, ensembl_transport
    ):
        locus = fetch_locus("myo7aa", flank=500, transport=ensembl_transport)
        header = locus.header()
        for expected in ("GRCz11", "18:500-1699", "ENSDARG00000099871", "myosin VIIAa"):
            assert expected in header
        assert "Ensembl REST" in header
        assert locus.retrieved in header

    def test_the_record_round_trips_through_fasta(self, ensembl_transport, tmp_path):
        from danioprime.sequence import write_fasta

        locus = fetch_locus("myo7aa", transport=ensembl_transport)
        path = tmp_path / "myo7aa.fa"
        write_fasta([locus.to_record()], path)
        [back] = read_fasta(path)
        assert back.name == "myo7aa"
        assert back.seq == locus.seq.upper()
        assert "GRCz11" in back.description

    def test_to_dict_is_json_serialisable(self, ensembl_transport):
        locus = fetch_locus("myo7aa", transport=ensembl_transport)
        assert json.loads(json.dumps(locus.to_dict()))["assembly"] == "GRCz11"


class TestProviderResolution:
    def test_unknown_source(self):
        with pytest.raises(AnnotationError, match="unknown source"):
            resolve_provider("genbank")

    def test_ucsc_cannot_resolve_symbols_and_says_so(self):
        with pytest.raises(AnnotationError, match="cannot look up a gene symbol"):
            resolve_provider("ucsc").lookup_url("myo7aa")


class TestTransportAgainstARealServer:
    """The urllib transport, over HTTP, with no mocking.

    Everything above fakes the network; this does not. It leaves the remote
    service's own behaviour as the only thing not covered by a test.
    """

    def test_it_fetches_a_body(self, tmp_path):
        (tmp_path / "seq.fa").write_text(fasta("ACGT" * 50), newline="\n")
        with RangeServer(tmp_path) as base:
            body = UrllibTransport().get(f"{base}/seq.fa")
        assert "ACGT" in body

    def test_a_region_fetch_works_end_to_end_over_http(self, tmp_path, monkeypatch):
        """Real sockets, real FASTA parsing, real length check."""
        (tmp_path / "region").write_text(fasta("ACGT" * 50), newline="\n")
        with RangeServer(tmp_path) as base:
            monkeypatch.setattr(
                EnsemblRest, "region_url", lambda self, c, s, e, strand=1: f"{base}/region"
            )
            locus = fetch_locus(
                "x", chrom="18", start=1, end=200, transport=UrllibTransport()
            )
        assert len(locus.seq) == 200
        assert locus.urls[0].startswith(base)
        assert locus.retrieved

    def test_a_truncated_response_is_caught_over_http_too(self, tmp_path, monkeypatch):
        (tmp_path / "region").write_text(fasta("ACGT" * 5), newline="\n")
        with RangeServer(tmp_path) as base:
            monkeypatch.setattr(
                EnsemblRest, "region_url", lambda self, c, s, e, strand=1: f"{base}/region"
            )
            with pytest.raises(AnnotationError, match="but 200 bp was requested"):
                fetch_locus("x", chrom="18", start=1, end=200, transport=UrllibTransport())

    def test_a_404_names_the_url(self, tmp_path):
        with RangeServer(tmp_path) as base:
            with pytest.raises(AnnotationError, match="404"):
                UrllibTransport().get(f"{base}/missing.fa")


class TestThroughTheCli:
    def test_dry_run_prints_urls_and_writes_nothing(self, capsys, tmp_path):
        out = tmp_path / "x.fa"
        assert main(["fetch", "myo7aa", "--dry-run", "--out", str(out)]) == 0
        captured = capsys.readouterr().out
        assert "rest.ensembl.org/lookup/symbol/danio_rerio/myo7aa" in captured
        assert "Nothing was fetched" in captured
        assert not out.exists(), "--dry-run must not create or clobber --out"

    def test_a_region_dry_run_shows_the_resolved_span(self, capsys):
        assert main(["fetch", "chr18:1,000,000-1,002,000", "--dry-run"]) == 0
        captured = capsys.readouterr().out
        assert "18:1000000..1002000:1" in captured
        assert "chr18:1,000,000-1,002,000" in captured

    def test_ucsc_dry_run_uses_ucsc_conventions(self, capsys):
        assert main(["fetch", "18:1000-2000", "--source", "ucsc", "--dry-run"]) == 0
        captured = capsys.readouterr().out
        assert "chrom=chr18" in captured and "start=999" in captured
        assert "danRer11" in captured

    def test_a_bare_interval_without_a_contig_is_rejected(self, capsys):
        """Otherwise it would be sent off as a gene symbol."""
        assert main(["fetch", "1000-2000", "--dry-run"]) == 2
        assert "needs a contig" in capsys.readouterr().err

    def test_ucsc_with_a_symbol_explains_the_alternative(self, capsys):
        assert main(["fetch", "myo7aa", "--source", "ucsc", "--dry-run"]) == 2
        assert "--source ensembl" in capsys.readouterr().err


class TestWholeChainAgainstAStub:
    """Fetch, write, design -- over real HTTP, nothing mocked.

    This is the test that would have caught a wrong URL, a coordinate that was
    1 out, a FASTA the design pipeline cannot read, or an assembly that never
    reached the header. The only thing it does not establish is whether the
    real Ensembl answers in the shape ``tests/stub_ensembl.py`` records.
    """

    def test_fetch_then_genotype(self, tmp_path, capsys):
        from stub_ensembl import StubEnsembl

        locus_path = tmp_path / "myo7aa_locus.fa"
        with StubEnsembl() as base:
            assert (
                main(
                    [
                        "fetch", "myo7aa",
                        "--base-url", base,
                        "--flank", "2000",
                        "--out", str(locus_path),
                    ]
                )
                == 0
            )
        written = capsys.readouterr().out
        assert "GRCz11" in written
        assert "--assembly GRCz11" in written, "the next command is spelled out"

        [record] = read_fasta(locus_path)
        assert record.name == "myo7aa"
        # The stub's gene is 3001..5000; 2 kb of flank on each side.
        assert len(record.seq) == 2000 + 2 * 2000
        assert "GRCz11 18:1001-7000" in record.description
        assert "ENSDARG00000099871" in record.description

        assert (
            main(
                [
                    "genotype",
                    "--template", str(locus_path),
                    "--cut-site", "3000",
                    "--tails", "fluorescent",
                    "--show", "1",
                ]
            )
            == 0
        )
        report = capsys.readouterr().out
        assert "myo7aa" in report, "the design is named after the fetched record"
        assert "tgtaaaacgacggccagt" in report, "and tailed"

    def test_fetch_a_transcript_then_design_a_probe(self, tmp_path, capsys):
        from stub_ensembl import StubEnsembl

        cdna = tmp_path / "myo7aa_cdna.fa"
        with StubEnsembl() as base:
            assert main(["fetch", "myo7aa", "--transcript", "--base-url", base,
                         "--out", str(cdna)]) == 0
        capsys.readouterr()
        [record] = read_fasta(cdna)
        assert "spliced transcript" in record.description
        assert "ENSDART00000148103" in record.description, "the canonical one"

        assert main(["insitu", "--template", str(cdna), "--scheme", "riboprobe-t7",
                     "--show", "1"]) == 0
        assert "taatacgactcactataggg" in capsys.readouterr().out, "T7, lower case"

    def test_an_unknown_symbol_is_reported_not_guessed_at(self, tmp_path, capsys):
        from stub_ensembl import StubEnsembl

        with StubEnsembl() as base:
            assert main(["fetch", "nosuchgene", "--base-url", base]) == 2
        assert "nosuchgene" in capsys.readouterr().err

    def test_json_output_records_the_provenance(self, capsys):
        from stub_ensembl import StubEnsembl

        with StubEnsembl() as base:
            assert main(["fetch", "myo7aa", "--base-url", base, "--format", "json"]) == 0
        payload = json.loads(capsys.readouterr().out)
        assert payload["assembly"] == "GRCz11"
        assert payload["gene_id"] == "ENSDARG00000099871"
        assert payload["length"] == payload["end"] - payload["start"] + 1
        assert payload["urls"] and all(u.startswith("http") for u in payload["urls"])

    def test_the_ohnolog_can_be_fetched_as_a_panel(self, tmp_path, capsys):
        """myo7aa has an ohnolog, which is what --avoid is for."""
        from stub_ensembl import StubEnsembl

        a, b = tmp_path / "a.fa", tmp_path / "b.fa"
        with StubEnsembl() as base:
            for symbol, path in (("myo7aa", a), ("myo7ab", b)):
                assert main(["fetch", symbol, "--base-url", base, "--flank", "500",
                             "--out", str(path)]) == 0
        capsys.readouterr()
        assert main(["genotype", "--template", str(a), "--cut-site", "1500",
                     "--avoid", str(b), "--show", "1"]) == 0
        assert "myo7aa" in capsys.readouterr().out
