"""VCF parsing, strain filtering and mask construction."""

from __future__ import annotations

import gzip

import pytest

from danioprime.design.regions import Interval
from danioprime.errors import AnnotationError
from danioprime.variants.vcf import load_mask, read_vcf, summarise, variants_to_intervals

VCF = """##fileformat=VCFv4.2
##contig=<ID=chr5,length=100000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTU\tAB
chr5\t101\tsnp1\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
chr5\t205\tsnp2\tC\tT\t50\tPASS\t.\tGT\t1/1\t0/1
chr5\t310\tins1\tT\tTAAG\t50\tPASS\t.\tGT\t0/0\t0/1
chr5\t400\tsnp3\tG\tA\t50\tPASS\t.\tGT\t0/0\t0/0
chr5\t500\tdel1\tGATC\tG\t50\tPASS\t.\tGT\t0/1\t0/0
chr6\t600\tsnp4\tG\tC\t50\tPASS\t.\tGT\t0/1\t0/1
"""


@pytest.fixture
def vcf_path(tmp_path):
    path = tmp_path / "v.vcf"
    path.write_text(VCF)
    return path


@pytest.fixture
def vcf_gz(tmp_path):
    path = tmp_path / "v.vcf.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(VCF)
    return path


class TestParsing:
    def test_reads_every_record(self, vcf_path):
        variants = read_vcf(vcf_path)
        assert len(variants) == 6
        assert [v.identifier for v in variants][:2] == ["snp1", "snp2"]

    def test_gzip_is_transparent(self, vcf_path, vcf_gz):
        assert len(read_vcf(vcf_gz)) == len(read_vcf(vcf_path))

    def test_contig_filter(self, vcf_path):
        assert len(read_vcf(vcf_path, chrom="chr5")) == 5
        assert len(read_vcf(vcf_path, chrom="chr6")) == 1

    def test_positions_are_converted_to_zero_based(self, vcf_path):
        snp1 = read_vcf(vcf_path)[0]
        assert snp1.pos == 100
        assert snp1.interval == Interval(100, 101)

    def test_indel_footprint_covers_the_ref_allele(self, vcf_path):
        del1 = next(v for v in read_vcf(vcf_path) if v.identifier == "del1")
        assert del1.is_indel
        assert del1.interval == Interval(499, 503)

    def test_variant_classification(self, vcf_path):
        counts = summarise(read_vcf(vcf_path))
        assert counts["snv"] == 4
        assert counts["indel"] == 2

    def test_offset_shifts_into_template_coordinates(self, vcf_path):
        variants = read_vcf(vcf_path, chrom="chr5", offset=100)
        assert variants[0].pos == 0

    def test_records_before_the_offset_are_dropped(self, vcf_path):
        variants = read_vcf(vcf_path, chrom="chr5", offset=300)
        assert all(v.pos >= 0 for v in variants)
        assert "snp1" not in [v.identifier for v in variants]

    def test_malformed_line_is_reported(self, tmp_path):
        path = tmp_path / "bad.vcf"
        path.write_text("#CHROM\tPOS\tID\tREF\tALT\nchr1\tnotanumber\tx\tA\tG\n")
        with pytest.raises(AnnotationError, match="bad POS"):
            read_vcf(path)


class TestStrainFiltering:
    def test_only_variants_carried_by_the_strain_are_kept(self, vcf_path):
        tu = [v.identifier for v in read_vcf(vcf_path, chrom="chr5", samples=["TU"])]
        assert tu == ["snp1", "snp2", "del1"]
        ab = [v.identifier for v in read_vcf(vcf_path, chrom="chr5", samples=["AB"])]
        assert ab == ["snp2", "ins1"]

    def test_homozygous_alternate_positions_are_not_hazards(self, vcf_path):
        # snp2 is 1/1 in TU: the line matches the alternate allele consistently,
        # so a primer designed against it works -- it is not a design hazard.
        het = [
            v.identifier
            for v in read_vcf(vcf_path, chrom="chr5", samples=["TU"], heterozygous_only=True)
        ]
        assert "snp2" not in het
        assert "snp1" in het

    def test_multiple_strains_are_unioned(self, vcf_path):
        both = [v.identifier for v in read_vcf(vcf_path, chrom="chr5", samples=["TU", "AB"])]
        assert set(both) == {"snp1", "snp2", "ins1", "del1"}

    def test_unknown_sample_lists_the_alternatives(self, vcf_path):
        with pytest.raises(AnnotationError, match="TU"):
            read_vcf(vcf_path, samples=["WIK"])

    def test_sample_filtering_without_genotypes_is_an_error(self, tmp_path):
        # Caught at the header, which lets the message say what is available.
        path = tmp_path / "nogt.vcf"
        path.write_text("#CHROM\tPOS\tID\tREF\tALT\nchr1\t10\tx\tA\tG\n")
        with pytest.raises(AnnotationError, match="available: none"):
            read_vcf(path, samples=["TU"])


class TestMask:
    def test_intervals_are_merged_and_padded(self, vcf_path):
        variants = read_vcf(vcf_path, chrom="chr5")
        plain = variants_to_intervals(variants)
        padded = variants_to_intervals(variants, pad=5)
        assert padded.total_length > plain.total_length

    def test_kind_filter(self, vcf_path):
        variants = read_vcf(vcf_path, chrom="chr5")
        only_indels = variants_to_intervals(variants, kinds=["indel"])
        assert only_indels.total_length == 1 + 4  # ins1 (1 bp REF) + del1 (4 bp REF)

    def test_load_mask_handles_vcf(self, vcf_path):
        mask, variants = load_mask(vcf_path, chrom="chr5", samples=["AB"])
        assert len(variants) == 2
        assert mask.total_length == 2

    def test_load_mask_handles_bed(self, tmp_path):
        path = tmp_path / "m.bed"
        path.write_text("chr5\t100\t110\n")
        mask, variants = load_mask(path, chrom="chr5")
        assert variants == []
        assert mask.intervals == [Interval(100, 110)]

    def test_bed_offset_is_applied(self, tmp_path):
        path = tmp_path / "m.bed"
        path.write_text("chr5\t100\t110\n")
        mask, _ = load_mask(path, chrom="chr5", offset=50)
        assert mask.intervals == [Interval(50, 60)]
