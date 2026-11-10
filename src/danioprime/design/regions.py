"""Half-open interval arithmetic for design constraints.

Three region roles, matching Primer3's vocabulary:

``include``
    Primers must fall entirely inside one of these intervals. Empty means the
    whole template is available.
``exclude``
    Primers must not overlap any of these intervals. Used for known variants,
    repeats, the CRISPR cut window, and vector junctions.
``target``
    The amplicon must span every target interval. Used for the lesion being
    genotyped, so a pair that sits beside the edit instead of across it is
    rejected.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from ..errors import RegionError

_SPEC = re.compile(
    r"^\s*(?:(?P<chrom>[^\s:]+):)?(?P<start>[\d,]+)\s*(?:-|\.\.)\s*(?P<end>[\d,]+)\s*$"
)
_LEN_SPEC = re.compile(r"^\s*(?:(?P<chrom>[^\s:]+):)?(?P<start>[\d,]+)\s*\+\s*(?P<length>\d+)\s*$")
_POINT = re.compile(r"^\s*(?:(?P<chrom>[^\s:]+):)?(?P<pos>[\d,]+)\s*$")


@dataclass(frozen=True, order=True, slots=True)
class Interval:
    """Half-open, zero-based interval ``[start, end)``."""

    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise RegionError(f"interval start must be >= 0, got {self.start}")
        if self.end < self.start:
            raise RegionError(f"interval end {self.end} precedes start {self.start}")

    def __len__(self) -> int:
        return self.end - self.start

    def __contains__(self, pos: int) -> bool:
        return self.start <= pos < self.end

    def overlaps(self, other: Interval) -> bool:
        return self.start < other.end and other.start < self.end

    def contains_interval(self, other: Interval) -> bool:
        return self.start <= other.start and other.end <= self.end

    def intersect(self, other: Interval) -> Interval | None:
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        return Interval(start, end) if start < end else None

    def clamp(self, lo: int, hi: int) -> Interval:
        return Interval(max(self.start, lo), min(self.end, hi))

    def pad(self, left: int, right: int | None = None) -> Interval:
        right = left if right is None else right
        return Interval(max(0, self.start - left), self.end + right)

    def as_one_based(self) -> str:
        """Render as an inclusive 1-based string, the convention in papers."""
        return f"{self.start + 1}-{self.end}"


def parse_interval(spec: str, *, one_based: bool = True) -> tuple[str | None, Interval]:
    """Parse ``chr1:1000-1200``, ``1000-1200``, ``1000+200`` or a bare point.

    Coordinates are interpreted as inclusive 1-based (the genome-browser
    convention) unless *one_based* is False, and returned as half-open
    zero-based.
    """
    for pattern in (_SPEC, _LEN_SPEC, _POINT):
        m = pattern.match(spec)
        if not m:
            continue
        groups = m.groupdict()
        chrom = groups.get("chrom")
        if "end" in groups and groups.get("end") is not None:
            start = int(groups["start"].replace(",", ""))
            end = int(groups["end"].replace(",", ""))
        elif groups.get("length") is not None:
            start = int(groups["start"].replace(",", ""))
            end = start + int(groups["length"]) - (1 if one_based else 0)
        else:
            start = end = int(groups["pos"].replace(",", ""))
        if one_based:
            if start < 1:
                raise RegionError(f"1-based coordinates start at 1, got {start} in {spec!r}")
            return chrom, Interval(start - 1, end)
        return chrom, Interval(start, end)
    raise RegionError(
        f"could not parse interval {spec!r}; expected 'start-end', 'start+length' or 'pos'"
    )


class IntervalSet:
    """A normalised, merged, sorted collection of intervals."""

    __slots__ = ("_items",)

    def __init__(self, intervals: Iterable[Interval] = ()):
        self._items = self._merge(intervals)

    @staticmethod
    def _merge(intervals: Iterable[Interval]) -> list[Interval]:
        items = sorted((i for i in intervals if len(i) > 0), key=lambda i: (i.start, i.end))
        merged: list[Interval] = []
        for item in items:
            if merged and item.start <= merged[-1].end:
                last = merged[-1]
                merged[-1] = Interval(last.start, max(last.end, item.end))
            else:
                merged.append(item)
        return merged

    @classmethod
    def from_specs(cls, specs: Iterable[str], *, one_based: bool = True) -> IntervalSet:
        return cls(parse_interval(s, one_based=one_based)[1] for s in specs)

    def __iter__(self) -> Iterator[Interval]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __bool__(self) -> bool:
        return bool(self._items)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        inner = ", ".join(f"{i.start}-{i.end}" for i in self._items)
        return f"IntervalSet([{inner}])"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, IntervalSet) and self._items == other._items

    @property
    def intervals(self) -> list[Interval]:
        return list(self._items)

    @property
    def total_length(self) -> int:
        return sum(len(i) for i in self._items)

    def overlaps(self, query: Interval) -> bool:
        return any(i.overlaps(query) for i in self._items)

    def overlapping(self, query: Interval) -> list[Interval]:
        return [i for i in self._items if i.overlaps(query)]

    def overlap_length(self, query: Interval) -> int:
        total = 0
        for i in self._items:
            hit = i.intersect(query)
            if hit:
                total += len(hit)
        return total

    def contains_point(self, pos: int) -> bool:
        return any(pos in i for i in self._items)

    def contains_interval(self, query: Interval) -> bool:
        """True if *query* lies wholly inside a single member interval."""
        return any(i.contains_interval(query) for i in self._items)

    def union(self, other: IntervalSet) -> IntervalSet:
        return IntervalSet(self._items + other._items)

    def add(self, interval: Interval) -> IntervalSet:
        return IntervalSet([*self._items, interval])

    def complement(self, start: int, end: int) -> IntervalSet:
        """Intervals within ``[start, end)`` not covered by this set."""
        gaps: list[Interval] = []
        cursor = start
        for item in self._items:
            if item.end <= start or item.start >= end:
                continue
            if item.start > cursor:
                gaps.append(Interval(cursor, min(item.start, end)))
            cursor = max(cursor, item.end)
            if cursor >= end:
                break
        if cursor < end:
            gaps.append(Interval(cursor, end))
        return IntervalSet(gaps)

    def positions_in(self, query: Interval) -> list[int]:
        """Every covered position inside *query*."""
        out: list[int] = []
        for item in self._items:
            hit = item.intersect(query)
            if hit:
                out.extend(range(hit.start, hit.end))
        return sorted(out)


def read_bed(path: str, *, chrom: str | None = None) -> IntervalSet:
    """Read a BED file into an :class:`IntervalSet` (BED is already 0-based)."""
    intervals: list[Interval] = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 3:
                raise RegionError(f"{path}:{lineno}: BED needs at least 3 columns")
            if chrom is not None and parts[0] != chrom:
                continue
            try:
                intervals.append(Interval(int(parts[1]), int(parts[2])))
            except ValueError as exc:
                raise RegionError(f"{path}:{lineno}: bad coordinates: {exc}") from exc
    return IntervalSet(intervals)
