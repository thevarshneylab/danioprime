"""Interval arithmetic and coordinate parsing."""

from __future__ import annotations

import pytest

from danioprime.design.regions import Interval, IntervalSet, parse_interval, read_bed
from danioprime.errors import RegionError


class TestInterval:
    def test_length_and_membership(self):
        iv = Interval(10, 20)
        assert len(iv) == 10
        assert 10 in iv
        assert 19 in iv
        assert 20 not in iv

    def test_rejects_inverted_or_negative(self):
        with pytest.raises(RegionError, match="precedes start"):
            Interval(20, 10)
        with pytest.raises(RegionError, match=">= 0"):
            Interval(-1, 10)

    def test_overlap_and_containment(self):
        a, b = Interval(10, 20), Interval(15, 25)
        assert a.overlaps(b)
        assert not a.overlaps(Interval(20, 30))  # half-open: touching is not overlapping
        assert a.contains_interval(Interval(12, 18))
        assert not a.contains_interval(b)

    def test_intersect_and_pad_and_clamp(self):
        assert Interval(10, 20).intersect(Interval(15, 25)) == Interval(15, 20)
        assert Interval(10, 20).intersect(Interval(30, 40)) is None
        assert Interval(10, 20).pad(5) == Interval(5, 25)
        assert Interval(2, 20).pad(5) == Interval(0, 25)
        assert Interval(10, 20).clamp(12, 15) == Interval(12, 15)

    def test_one_based_rendering(self):
        assert Interval(0, 10).as_one_based() == "1-10"


class TestParsing:
    def test_range_forms(self):
        assert parse_interval("chr5:1,000-1200") == ("chr5", Interval(999, 1200))
        assert parse_interval("1000-1200") == (None, Interval(999, 1200))
        assert parse_interval("1000..1200")[1] == Interval(999, 1200)

    def test_length_and_point_forms(self):
        assert parse_interval("100+50")[1] == Interval(99, 149)
        assert parse_interval("250")[1] == Interval(249, 250)

    def test_zero_based_mode(self):
        assert parse_interval("0-10", one_based=False)[1] == Interval(0, 10)

    def test_rejects_garbage_and_zero(self):
        with pytest.raises(RegionError, match="could not parse"):
            parse_interval("not-a-range!")
        with pytest.raises(RegionError, match="start at 1"):
            parse_interval("0-10")


class TestIntervalSet:
    def test_merges_overlapping_and_sorts(self):
        s = IntervalSet.from_specs(["10-20", "15-30", "100-110"])
        assert s.intervals == [Interval(9, 30), Interval(99, 110)]
        assert s.total_length == 32

    def test_drops_empty_intervals(self):
        assert len(IntervalSet([Interval(5, 5), Interval(10, 20)])) == 1

    def test_overlap_queries(self):
        s = IntervalSet([Interval(10, 20), Interval(100, 110)])
        assert s.overlaps(Interval(15, 17))
        assert not s.overlaps(Interval(50, 60))
        assert s.overlap_length(Interval(15, 105)) == 5 + 5
        assert s.contains_point(10)
        assert not s.contains_point(20)
        assert s.contains_interval(Interval(12, 18))
        assert not s.contains_interval(Interval(15, 105))

    def test_complement(self):
        s = IntervalSet([Interval(10, 20)])
        assert s.complement(0, 30).intervals == [Interval(0, 10), Interval(20, 30)]
        assert IntervalSet().complement(0, 10).intervals == [Interval(0, 10)]

    def test_positions_in(self):
        s = IntervalSet([Interval(10, 13)])
        assert s.positions_in(Interval(11, 20)) == [11, 12]

    def test_union_and_add(self):
        a = IntervalSet([Interval(0, 10)])
        b = IntervalSet([Interval(5, 15)])
        assert a.union(b).intervals == [Interval(0, 15)]
        assert a.add(Interval(20, 25)).total_length == 15

    def test_bed_reading(self, tmp_path):
        path = tmp_path / "x.bed"
        path.write_text("track name=x\nchr1\t10\t20\tfeat\nchr2\t0\t5\n#comment\n")
        assert read_bed(str(path)).intervals == [Interval(0, 5), Interval(10, 20)]
        assert read_bed(str(path), chrom="chr1").intervals == [Interval(10, 20)]

    def test_bed_with_too_few_columns(self, tmp_path):
        path = tmp_path / "bad.bed"
        path.write_text("chr1\t10\n")
        with pytest.raises(RegionError, match="at least 3 columns"):
            read_bed(str(path))
