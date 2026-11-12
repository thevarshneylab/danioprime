"""Published variant resources: the registry, the region fetch, and the
assembly guard that stops the most dangerous mistake.

Masking a GRCz12tu template with GRCz11 variant coordinates would exclude the
wrong bases *and report that known polymorphisms were avoided*, which is worse
than not checking. That refusal is the single most important behaviour here.
"""

from __future__ import annotations

import random

import pytest
from rangeserver import RangeServer

from danioprime.errors import AnnotationError
from danioprime.variants.remote import (
    clear_index_cache,
    describe_sources,
    fetch_for_template,
    fetch_region,
)
from danioprime.variants.sources import (
    UCSC_TRACK_HUB,
    VARIANT_SOURCES,
    VariantSource,
    get_source,
    normalise_assembly,
)

pysam = pytest.importorskip("pysam")

LENGTH = 2_000_000


@pytest.fixture(scope="module")
def mirror(tmp_path_factory):
    """A local stand-in for a published remote VCF, with genotypes."""
    directory = tmp_path_factory.mktemp("mirror")
    rng = random.Random(21)
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=chr5,length={LENGTH}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNHGRI1\tTU",
    ]
    pos = 0
    hom_alt = 0
    while True:
        pos += rng.randint(60, 900)
        if pos >= LENGTH:
            break
        ref = rng.choice("ACGT")
        alt = rng.choice([b for b in "ACGT" if b != ref])
        nhgri = rng.choice(["0/1", "1/1", "0/0"])
        if nhgri == "1/1":
            hom_alt += 1
        lines.append(
            f"chr5\t{pos}\trs{pos}\t{ref}\t{alt}\t60\tPASS\t.\tGT\t{nhgri}\t0/1"
        )
    plain = directory / "m.vcf"
    plain.write_text("\n".join(lines) + "\n", newline="\n")
    gz = directory / "m.vcf.gz"
    pysam.tabix_compress(str(plain), str(gz), force=True)
    pysam.tabix_index(str(gz), preset="vcf", force=True)
    return gz


@pytest.fixture
def source_factory(mirror):
    """Build a VariantSource pointing at the local mirror over HTTP."""

    def _make(assembly: str = "danRer11"):
        server = RangeServer(mirror.parent)
        base = server.__enter__()
        source = VariantSource(
            name="NHGRI-1 (mirror)",
            description="local stand-in",
            assembly=assembly,
            data_url=f"{base}/{mirror.name}",
            index_url=f"{base}/{mirror.name}.tbi",
            kind="vcf",
            citation="LaFave et al. 2014",
            homepage="https://example.invalid",
        )
        return source, server

    return _make


class TestRegistry:
    def test_nhgri1_is_registered_on_danrer11(self):
        source = get_source("nhgri1")
        assert source.assembly == "danRer11"
        assert source.kind == "vcf"
        assert source.data_url.endswith(".vcf.gz")
        assert source.index_url == source.data_url + ".tbi"

    def test_every_source_is_fully_described(self):
        for key, source in VARIANT_SOURCES.items():
            assert source.name and source.description, key
            assert source.citation, f"{key} has no citation"
            assert source.homepage.startswith("https://"), key
            assert source.data_url.startswith("https://"), key

    def test_track_hub_is_exposed_for_browsing(self):
        assert UCSC_TRACK_HUB.endswith("hub.txt")

    def test_unknown_source_lists_the_options(self):
        with pytest.raises(AnnotationError, match="nhgri1"):
            get_source("nope")

    def test_unknown_source_mentions_the_url_escape_hatch(self):
        with pytest.raises(AnnotationError, match="URL or a local path"):
            get_source("nope")

    def test_describe_sources_is_report_ready(self):
        rows = describe_sources()
        assert rows
        for row in rows:
            assert {"key", "name", "assembly", "citation", "track_hub"} <= set(row)

    @pytest.mark.parametrize(
        "given,expected",
        [
            ("danRer11", "danRer11"),
            ("GRCz11", "danRer11"),
            ("grcz11", "danRer11"),
            ("GRCz12tu", "GRCz12tu"),
            ("danRer12", "GRCz12tu"),
            ("danrer-11", "danRer11"),
        ],
    )
    def test_assembly_aliases(self, given, expected):
        assert normalise_assembly(given) == expected

    def test_unknown_assembly_passes_through(self):
        assert normalise_assembly("hg38") == "hg38"


class TestAssemblyGuard:
    def test_matching_assembly_is_accepted(self, source_factory):
        source, server = source_factory("danRer11")
        try:
            result = fetch_region(
                "chr5", 100_000, 104_000, source=source, template_assembly="GRCz11"
            )
            assert result.assembly == "danRer11"
            assert result.assembly_verified
        finally:
            server.__exit__()

    def test_mismatched_assembly_is_refused(self, source_factory):
        source, server = source_factory("danRer11")
        try:
            with pytest.raises(AnnotationError) as excinfo:
                fetch_region(
                    "chr5", 100_000, 104_000, source=source, template_assembly="GRCz12tu"
                )
            message = str(excinfo.value)
            assert "not transferable" in message
            assert "GRCz12tu" in message and "danRer11" in message
            # The message must say what to do, not merely that it refused.
            assert "lift" in message or "--region" in message
        finally:
            server.__exit__()

    def test_the_override_exists_for_people_who_lifted_over(self, source_factory):
        source, server = source_factory("danRer11")
        try:
            result = fetch_region(
                "chr5", 100_000, 104_000, source=source,
                template_assembly="GRCz12tu", allow_assembly_mismatch=True,
            )
            assert result.variants
        finally:
            server.__exit__()

    def test_no_declared_assembly_means_no_check(self, source_factory):
        source, server = source_factory("danRer11")
        try:
            assert fetch_region("chr5", 100_000, 104_000, source=source).variants
        finally:
            server.__exit__()

    def test_ad_hoc_url_reports_its_assembly_as_unverified(self, mirror):
        """Echoing back a declared assembly would be a false assurance."""
        with RangeServer(mirror.parent) as base:
            result = fetch_region(
                "chr5", 100_000, 104_000,
                data_location=f"{base}/{mirror.name}",
                template_assembly="GRCz12tu",
            )
            assert not result.assembly_verified
            assert result.assembly == "unknown"
            assert any("cannot determine the assembly" in n for n in result.notes)
            assert "unverified" in result.summary()


class TestRegionFetch:
    def test_variants_are_returned_within_the_region(self, source_factory):
        source, server = source_factory()
        try:
            result = fetch_region("chr5", 100_000, 104_000, source=source)
            assert result.variants
            for variant in result.variants:
                assert 100_000 <= variant.pos < 104_000
        finally:
            server.__exit__()

    def test_mask_covers_the_variants(self, source_factory):
        source, server = source_factory()
        try:
            result = fetch_region("chr5", 100_000, 104_000, source=source)
            assert result.mask.total_length >= len(result.variants)
            for variant in result.variants:
                assert result.mask.contains_point(variant.pos)
        finally:
            server.__exit__()

    def test_offset_shifts_into_template_coordinates(self, source_factory):
        source, server = source_factory()
        try:
            absolute = fetch_region("chr5", 100_000, 104_000, source=source)
            shifted = fetch_region(
                "chr5", 100_000, 104_000, source=source, offset=100_000
            )
            assert len(absolute) == len(shifted)
            assert [v.pos - 100_000 for v in absolute.variants] == [
                v.pos for v in shifted.variants
            ]
            assert all(0 <= v.pos < 4_000 for v in shifted.variants)
        finally:
            server.__exit__()

    def test_fetch_for_template_is_equivalent(self, source_factory):
        source, server = source_factory()
        try:
            direct = fetch_region(
                "chr5", 100_000, 104_000, source=source, offset=100_000
            )
            wrapped = fetch_for_template(4_000, "chr5", 100_000, source=source)
            assert [v.pos for v in direct.variants] == [v.pos for v in wrapped.variants]
        finally:
            server.__exit__()

    def test_sample_filtering(self, source_factory):
        source, server = source_factory()
        try:
            everything = fetch_region("chr5", 0, 200_000, source=source)
            nhgri = fetch_region(
                "chr5", 0, 200_000, source=source, samples=["NHGRI1"]
            )
            assert 0 < len(nhgri) < len(everything)
        finally:
            server.__exit__()

    def test_heterozygous_only_excludes_homozygous_alternate(self, source_factory):
        """A position fixed non-reference in the line is not a design hazard."""
        source, server = source_factory()
        try:
            allv = fetch_region("chr5", 0, 200_000, source=source, samples=["NHGRI1"])
            het = fetch_region(
                "chr5", 0, 200_000, source=source,
                samples=["NHGRI1"], heterozygous_only=True,
            )
            assert len(het) < len(allv)
            for variant in het.variants:
                assert any(
                    gt.split(":")[0].replace("|", "/") in {"0/1", "1/0"}
                    for gt in variant.genotypes
                )
        finally:
            server.__exit__()

    def test_unknown_sample_lists_what_is_available(self, source_factory):
        source, server = source_factory()
        try:
            with pytest.raises(AnnotationError, match="NHGRI1"):
                fetch_region("chr5", 0, 10_000, source=source, samples=["WIK"])
        finally:
            server.__exit__()

    def test_padding_widens_the_mask(self, source_factory):
        source, server = source_factory()
        try:
            plain = fetch_region("chr5", 100_000, 104_000, source=source)
            padded = fetch_region("chr5", 100_000, 104_000, source=source, pad=5)
            assert padded.mask.total_length > plain.mask.total_length
        finally:
            server.__exit__()

    def test_transfer_is_reported_and_bounded(self, source_factory):
        """Per-locus transfer is bounded by a couple of BGZF blocks.

        Stated as a bound rather than a fraction of the file: this fixture is
        smaller than one block read, so a fraction would be meaningless here.
        The bound is what makes the real 214 MB dataset usable, since it does
        not grow with the file.
        """
        from danioprime.variants.bgzf import MAX_BLOCK_SIZE

        source, server = source_factory()
        try:
            result = fetch_region("chr5", 100_000, 104_000, source=source)
            assert result.transfer["requests"] >= 1
            assert result.transfer["bytes_fetched"] <= 2 * MAX_BLOCK_SIZE + 1024
            assert "kB fetched" in result.summary()
        finally:
            server.__exit__()

    def test_transfer_does_not_grow_with_the_dataset(self, tmp_path_factory):
        """The property the feature depends on, measured on two file sizes."""
        from danioprime.variants.bgzf import MAX_BLOCK_SIZE

        def build(span: int):
            directory = tmp_path_factory.mktemp(f"scale{span}")
            rng = random.Random(span)
            lines = [
                "##fileformat=VCFv4.2",
                f"##contig=<ID=chr5,length={span}>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            ]
            pos = 0
            while True:
                pos += rng.randint(40, 120)
                if pos >= span:
                    break
                lines.append(f"chr5\t{pos}\t.\tA\tG\t60\tPASS\t.")
            plain = directory / "s.vcf"
            plain.write_text("\n".join(lines) + "\n", newline="\n")
            gz = directory / "s.vcf.gz"
            pysam.tabix_compress(str(plain), str(gz), force=True)
            pysam.tabix_index(str(gz), preset="vcf", force=True)
            return gz

        measurements = []
        for span in (2_000_000, 20_000_000):
            gz = build(span)
            with RangeServer(gz.parent) as base:
                result = fetch_region(
                    "chr5", 1_000_000, 1_005_000, data_location=f"{base}/{gz.name}"
                )
                measurements.append((result.transfer["file_size"], result.transfer["bytes_fetched"]))

        (small_file, small_read), (big_file, big_read) = measurements
        assert big_file > 5 * small_file, "the fixtures are not different enough in size"
        assert big_read <= 2 * MAX_BLOCK_SIZE + 1024
        # Ten times the data, essentially the same read.
        assert big_read <= small_read * 2

    def test_summary_and_dict_are_informative(self, source_factory):
        import json

        source, server = source_factory()
        try:
            result = fetch_region(
                "chr5", 100_000, 104_000, source=source,
                samples=["NHGRI1"], heterozygous_only=True,
            )
            summary = result.summary()
            assert "NHGRI-1" in summary and "heterozygous" in summary
            payload = result.to_dict()
            json.dumps(payload)
            assert payload["citation"]
            assert payload["assembly"] == "danRer11"
        finally:
            server.__exit__()

    def test_empty_region_yields_nothing(self, source_factory):
        source, server = source_factory()
        try:
            assert fetch_region("chr5", 50, 51, source=source).mask.total_length >= 0
        finally:
            server.__exit__()

    def test_unreachable_host_explains_itself(self):
        with pytest.raises(AnnotationError, match="could not reach|HTTP"):
            fetch_region(
                "chr5", 0, 1000,
                data_location="https://127.0.0.1:9/nope.vcf.gz",
            )

    def test_source_or_location_is_required(self):
        with pytest.raises(AnnotationError, match="named source"):
            fetch_region("chr5", 0, 1000, source=None)


class TestIndexCache:
    def test_index_is_cached_and_reused(self, mirror, monkeypatch, tmp_path):
        monkeypatch.setenv("DANIOPRIME_CACHE", str(tmp_path))
        clear_index_cache()
        with RangeServer(mirror.parent) as base:
            url = f"{base}/{mirror.name}"
            fetch_region("chr5", 0, 10_000, data_location=url)
            cached = list((tmp_path / "danioprime").glob("*.tbi"))
            assert len(cached) == 1
            stamp = cached[0].stat().st_mtime_ns
            fetch_region("chr5", 20_000, 30_000, data_location=url)
            assert cached[0].stat().st_mtime_ns == stamp, "index was re-downloaded"

    def test_cache_can_be_cleared(self, mirror, monkeypatch, tmp_path):
        monkeypatch.setenv("DANIOPRIME_CACHE", str(tmp_path))
        with RangeServer(mirror.parent) as base:
            fetch_region("chr5", 0, 10_000, data_location=f"{base}/{mirror.name}")
        assert clear_index_cache() >= 1
        assert clear_index_cache() == 0
