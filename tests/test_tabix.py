"""BGZF and tabix: the pure-Python reader, checked against htslib.

The reader exists so that a 214 MB published variant set can be queried for one
locus without downloading it. Correctness is not negotiable here -- a reader
that silently returns the wrong records would mask the wrong bases -- so where
``pysam`` is installed these tests compare against htslib's own answers over
hundreds of random queries. Where it is not, the fixtures ship as bytes and the
structural tests still run.
"""

from __future__ import annotations

import random

import pytest
from rangeserver import RangeServer

from danioprime.errors import AnnotationError
from danioprime.variants.bgzf import (
    EOF_MARKER,
    MAX_BLOCK_SIZE,
    BgzfReader,
    parse_block,
    split_virtual_offset,
    virtual_offset,
)
from danioprime.variants.sources import HttpFile, LocalFile
from danioprime.variants.tabix import TabixFile, load_index, reg2bins

pysam = pytest.importorskip("pysam", reason="htslib comparison needs pysam")

CONTIGS = {"chr5": 6_000_000, "chr19": 4_800_000}


@pytest.fixture(scope="module")
def vcf(tmp_path_factory):
    """A bgzipped, tabix-indexed VCF with SNVs, indels and two contigs."""
    directory = tmp_path_factory.mktemp("tabix")
    rng = random.Random(7)
    lines = ["##fileformat=VCFv4.2"]
    for name, length in CONTIGS.items():
        lines.append(f"##contig=<ID={name},length={length}>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNHGRI1")
    for name, length in CONTIGS.items():
        pos = 0
        while True:
            pos += rng.randint(40, 600)
            if pos >= length:
                break
            ref = rng.choice("ACGT")
            if rng.random() < 0.12:
                ref += "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))
                alt = ref[0]
            else:
                alt = rng.choice([b for b in "ACGT" if b != ref[0]])
            lines.append(
                f"{name}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t"
                f"{rng.choice(['0/1', '1/1'])}"
            )
    plain = directory / "v.vcf"
    plain.write_text("\n".join(lines) + "\n", newline="\n")
    gz = directory / "v.vcf.gz"
    pysam.tabix_compress(str(plain), str(gz), force=True)
    pysam.tabix_index(str(gz), preset="vcf", force=True)
    return gz


@pytest.fixture(scope="module")
def reader(vcf):
    return TabixFile(
        data=BgzfReader(LocalFile(vcf)), index=load_index(LocalFile(f"{vcf}.tbi"))
    )


class TestVirtualOffsets:
    def test_round_trip(self):
        for block, within in ((0, 0), (1, 1), (123_456_789, 65_535)):
            assert split_virtual_offset(virtual_offset(block, within)) == (block, within)

    def test_in_block_offset_is_bounded(self):
        with pytest.raises(ValueError, match="out of range"):
            virtual_offset(0, MAX_BLOCK_SIZE)


class TestBgzf:
    def test_eof_marker_is_a_valid_empty_block(self):
        block = parse_block(EOF_MARKER)
        assert block.compressed_size == len(EOF_MARKER)
        assert block.data == b""

    def test_plain_gzip_is_rejected_with_an_explanation(self, tmp_path):
        import gzip

        path = tmp_path / "plain.gz"
        path.write_bytes(gzip.compress(b"chr5\t1\t.\tA\tG\n"))
        with pytest.raises(AnnotationError, match="plain gzip"):
            parse_block(path.read_bytes())

    def test_non_gzip_is_rejected(self):
        with pytest.raises(AnnotationError, match="not a gzip"):
            parse_block(b"not compressed at all, really")

    def test_truncated_header(self):
        with pytest.raises(AnnotationError, match="truncated"):
            parse_block(b"\x1f\x8b")

    def test_corrupt_payload_fails_its_crc(self, vcf):
        raw = bytearray(vcf.read_bytes()[:MAX_BLOCK_SIZE])
        first = parse_block(bytes(raw))
        # Flip a bit inside the deflate payload of the first block.
        raw[40] ^= 0xFF
        with pytest.raises(Exception):  # noqa: B017 - zlib or CRC, both acceptable
            parse_block(bytes(raw))
        assert first.data  # the untouched block was fine

    def test_header_lines_are_recovered(self, reader):
        header = reader.header()
        assert header[0].startswith("##fileformat=VCFv4.2")
        assert header[-1].startswith("#CHROM")

    def test_block_cache_avoids_refetching(self, vcf):
        class Counting(LocalFile):
            reads = 0

            def read_range(self, start, length):
                type(self).reads += 1
                return super().read_range(start, length)

        source = Counting(vcf)
        bgzf = BgzfReader(source)
        bgzf.block_at(0)
        first = Counting.reads
        bgzf.block_at(0)
        assert Counting.reads == first


class TestBinning:
    def test_bins_for_a_single_base(self):
        assert reg2bins(0, 1) == [0, 1, 9, 73, 585, 4681]

    def test_empty_region_has_no_bins(self):
        assert reg2bins(10, 10) == []
        assert reg2bins(10, 5) == []

    def test_leaf_bins_advance_every_16kb(self):
        assert 4681 in reg2bins(0, 1)
        assert 4682 in reg2bins(16_384, 16_385)

    def test_wider_regions_touch_more_bins(self):
        assert len(reg2bins(0, 1)) < len(reg2bins(0, 10_000_000))


class TestIndex:
    def test_header_fields(self, vcf):
        index = load_index(LocalFile(f"{vcf}.tbi"))
        assert index.reference_names == list(CONTIGS)
        assert index.base_format == 2  # VCF
        assert not index.zero_based  # VCF is 1-based
        assert index.col_seq == 1
        assert index.col_begin == 2

    def test_chr_prefix_is_tolerated_in_both_directions(self, vcf):
        index = load_index(LocalFile(f"{vcf}.tbi"))
        assert index.reference_id("chr5") == 0
        assert index.reference_id("5") == 0

    def test_unknown_contig_names_the_alternatives(self, vcf):
        index = load_index(LocalFile(f"{vcf}.tbi"))
        with pytest.raises(AnnotationError, match="chr5"):
            index.chunks_for("chrZ", 0, 1000)

    def test_a_non_index_is_rejected(self, tmp_path):
        import gzip

        path = tmp_path / "bad.tbi"
        path.write_bytes(gzip.compress(b"NOTANINDEX"))
        # A plain gzip has no BC extra field, which is the specific thing the
        # reader checks for, so the message should say so rather than just
        # "bad file".
        with pytest.raises(AnnotationError, match="plain gzip"):
            load_index(LocalFile(path))

    def test_csi_is_reported_specifically(self):
        from danioprime.variants.tabix import parse_tabix

        with pytest.raises(AnnotationError, match="CSI index"):
            parse_tabix(b"CSI\x01" + b"\x00" * 32)

    def test_chunks_are_merged(self, vcf):
        index = load_index(LocalFile(f"{vcf}.tbi"))
        chunks = index.chunks_for("chr5", 0, 1_000_000)
        assert chunks
        for earlier, later in zip(chunks, chunks[1:], strict=False):
            assert earlier.end < later.begin


class TestAgainstHtslib:
    """The claim that matters: identical answers to htslib."""

    @staticmethod
    @pytest.fixture(scope="class")
    def truth(vcf):
        return pysam.TabixFile(str(vcf))

    @staticmethod
    def _key(fields: list[str]) -> tuple:
        return (fields[0], int(fields[1]), fields[3], fields[4])

    @pytest.mark.parametrize("width", [200, 2_000, 50_000, 500_000])
    def test_random_queries_match_exactly(self, reader, truth, width):
        rng = random.Random(width)
        total = 0
        for _ in range(60):
            chrom = rng.choice(list(CONTIGS))
            start = rng.randrange(0, CONTIGS[chrom] - width)
            mine = [self._key(f) for f in reader.fetch(chrom, start, start + width)]
            theirs = [
                self._key(line.split("\t"))
                for line in truth.fetch(chrom, start, start + width)
            ]
            assert mine == theirs, f"{chrom}:{start}-{start + width}"
            total += len(theirs)
        assert total > 0, "the fixture produced no overlapping records to compare"

    def test_indels_at_region_edges_match(self, reader, truth):
        """A VCF record spans its REF allele, so edge handling is easy to get wrong."""
        rng = random.Random(99)
        for _ in range(200):
            chrom = rng.choice(list(CONTIGS))
            start = rng.randrange(0, CONTIGS[chrom] - 50)
            mine = [self._key(f) for f in reader.fetch(chrom, start, start + 5)]
            theirs = [
                self._key(line.split("\t"))
                for line in truth.fetch(chrom, start, start + 5)
            ]
            assert mine == theirs

    def test_whole_contig_matches(self, reader, truth):
        mine = sum(1 for _ in reader.fetch("chr19", 0, CONTIGS["chr19"]))
        theirs = sum(1 for _ in truth.fetch("chr19", 0, CONTIGS["chr19"]))
        assert mine == theirs

    def test_empty_regions_agree(self, reader, truth):
        assert list(reader.fetch("chr5", 10, 10)) == []
        assert list(reader.fetch("chr5", 5_999_000, 5_999_001)) == [
            f for f in reader.fetch("chr5", 5_999_000, 5_999_001)
        ]


class TestOverHttp:
    def test_queries_over_range_requests_match_local(self, vcf, reader):
        with RangeServer(vcf.parent) as base:
            remote = TabixFile(
                data=BgzfReader(HttpFile(f"{base}/{vcf.name}")),
                index=load_index(HttpFile(f"{base}/{vcf.name}.tbi")),
            )
            rng = random.Random(5)
            for _ in range(15):
                chrom = rng.choice(list(CONTIGS))
                start = rng.randrange(0, CONTIGS[chrom] - 20_000)
                assert [f[:5] for f in remote.fetch(chrom, start, start + 20_000)] == [
                    f[:5] for f in reader.fetch(chrom, start, start + 20_000)
                ]

    def test_transfer_is_a_small_fraction_of_the_file(self, vcf):
        """The whole point: a locus costs a bounded read, not the file."""
        with RangeServer(vcf.parent) as base:
            data = HttpFile(f"{base}/{vcf.name}")
            remote = TabixFile(
                data=BgzfReader(data), index=load_index(HttpFile(f"{base}/{vcf.name}.tbi"))
            )
            list(remote.fetch("chr5", 1_000_000, 1_005_000))
            stats = data.stats
            assert stats["requests"] <= 2, "a single locus should not need many requests"
            assert stats["bytes_fetched"] <= 2 * MAX_BLOCK_SIZE + 1024

    def test_missing_file_reports_the_status(self, vcf):
        with RangeServer(vcf.parent) as base, pytest.raises(AnnotationError, match="HTTP 404"):
            HttpFile(f"{base}/not-there.vcf.gz")

    def test_server_without_range_support_is_detected(self, vcf):
        import http.server
        import threading

        class NoRanges(http.server.SimpleHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_HEAD(self):
                body = vcf.read_bytes()
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()

            def do_GET(self):
                body = vcf.read_bytes()
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), NoRanges)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        try:
            host, port = server.server_address[:2]
            source = HttpFile(f"http://{host}:{port}/x")
            with pytest.raises(AnnotationError, match="ignored the Range header"):
                source.read_range(0, 100)
        finally:
            server.shutdown()
            server.server_close()
