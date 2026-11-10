"""Strain-variant masking.

Wild-type zebrafish lines are not isogenic. TU, AB, TL and WIK each carry
substantial standing variation, and a primer sitting on a SNP that is
polymorphic in the line you are working with simply fails -- often silently, as
one allele amplifying and the other not, which looks exactly like a genotyping
result and is the worst kind of failure.

Variants are therefore treated as first-class design constraints:

* a variant anywhere in the primer footprint is penalised;
* a variant within the 3'-terminal window is disqualifying, because that is
  where a mismatch stops extension outright;
* ``--reject-variant-overlap`` makes any overlap disqualifying.

Input may be VCF (plain or gzipped) or BED. For VCF, per-sample genotypes can
be used to keep only variants that are actually segregating in the line being
worked on, which matters because a position fixed as homozygous non-reference
in your line is not a problem -- the primer just needs to match the line, not
the reference.
"""

from __future__ import annotations

import gzip
import io
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO

from ..design.regions import Interval, IntervalSet
from ..errors import AnnotationError


@dataclass(frozen=True, slots=True)
class Variant:
    """A single VCF record, reduced to what primer design needs."""

    chrom: str
    pos: int  # zero-based start of REF
    ref: str
    alt: tuple[str, ...]
    identifier: str = "."
    genotypes: tuple[str, ...] = ()

    @property
    def end(self) -> int:
        """Half-open end of the REF allele footprint."""
        return self.pos + max(1, len(self.ref))

    @property
    def interval(self) -> Interval:
        return Interval(self.pos, self.end)

    @property
    def is_indel(self) -> bool:
        return any(len(a) != len(self.ref) for a in self.alt)

    @property
    def kind(self) -> str:
        if self.is_indel:
            return "indel"
        if len(self.ref) == 1 and all(len(a) == 1 for a in self.alt):
            return "snv"
        return "mnv"


def _open_text(path: str | Path) -> IO[str]:
    path = Path(path)
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return path.open("r")


_HOM_REF = {"0/0", "0|0", "0", "./.", ".|.", "."}


def _genotype_is_variant(gt_field: str) -> bool:
    """True if a VCF GT sub-field carries a non-reference allele."""
    gt = gt_field.split(":", 1)[0]
    return gt not in _HOM_REF and any(a not in {"0", "."} for a in gt.replace("|", "/").split("/"))


def _genotype_is_homozygous_alt(gt_field: str) -> bool:
    gt = gt_field.split(":", 1)[0].replace("|", "/")
    alleles = [a for a in gt.split("/") if a != "."]
    return bool(alleles) and len(set(alleles)) == 1 and alleles[0] not in {"0", "."}


def parse_vcf(
    handle: Iterable[str],
    *,
    chrom: str | None = None,
    samples: Iterable[str] | None = None,
    heterozygous_only: bool = False,
    offset: int = 0,
) -> Iterator[Variant]:
    """Stream variants from VCF lines.

    Parameters
    ----------
    chrom:
        Keep only records on this contig.
    samples:
        Keep only records carrying a non-reference genotype in at least one of
        these samples. Requires genotype columns.
    heterozygous_only:
        Keep only records that are heterozygous in a selected sample. A position
        homozygous for the alternate allele in your line is not a design hazard,
        so excluding it recovers primer sites that naive masking throws away.
    offset:
        Subtract this many bases from VCF positions, for when the template is a
        slice of the contig (``--template-offset``).
    """
    wanted = list(samples) if samples else []
    sample_names: list[str] = []
    for lineno, line in enumerate(handle, 1):
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            cols = line.split("\t")
            sample_names = cols[9:] if len(cols) > 9 else []
            missing = [s for s in wanted if s not in sample_names]
            if missing:
                raise AnnotationError(
                    f"samples not present in VCF: {missing}; available: {sample_names or 'none'}"
                )
            continue
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            raise AnnotationError(f"VCF line {lineno}: expected at least 5 columns")
        rec_chrom, pos_s, vid, ref, alt_s = f[0], f[1], f[2], f[3].upper(), f[4].upper()
        if chrom is not None and rec_chrom != chrom:
            continue
        try:
            pos1 = int(pos_s)
        except ValueError as exc:
            raise AnnotationError(f"VCF line {lineno}: bad POS {pos_s!r}") from exc
        start = pos1 - 1 - offset
        if start < 0:
            continue
        alts = tuple(a for a in alt_s.split(",") if a not in {".", ""})
        gts = tuple(f[9:]) if len(f) > 9 else ()

        if wanted:
            if not sample_names:
                raise AnnotationError("sample filtering requested but VCF has no genotype columns")
            idxs = [sample_names.index(s) for s in wanted]
            relevant = [gts[i] for i in idxs if i < len(gts)]
            if not any(_genotype_is_variant(g) for g in relevant):
                continue
            if heterozygous_only and not any(
                _genotype_is_variant(g) and not _genotype_is_homozygous_alt(g) for g in relevant
            ):
                continue

        yield Variant(
            chrom=rec_chrom,
            pos=start,
            ref=ref,
            alt=alts,
            identifier=vid,
            genotypes=gts,
        )


def read_vcf(path: str | Path, **kwargs) -> list[Variant]:
    """Read every (filtered) variant from a VCF file."""
    with _open_text(path) as fh:
        return list(parse_vcf(fh, **kwargs))


def variants_to_intervals(
    variants: Iterable[Variant], *, pad: int = 0, kinds: Iterable[str] | None = None
) -> IntervalSet:
    """Convert variants to a maskable interval set.

    *pad* widens each variant, which is worth doing for indels: an indel shifts
    everything downstream of it, so a primer immediately 3' of a polymorphic
    indel may anneal but place its 3' end at an unexpected position.
    """
    keep = set(kinds) if kinds else None
    intervals: list[Interval] = []
    for v in variants:
        if keep and v.kind not in keep:
            continue
        iv = v.interval
        intervals.append(iv.pad(pad) if pad else iv)
    return IntervalSet(intervals)


def load_mask(
    path: str | Path,
    *,
    chrom: str | None = None,
    samples: Iterable[str] | None = None,
    heterozygous_only: bool = False,
    offset: int = 0,
    pad: int = 0,
) -> tuple[IntervalSet, list[Variant]]:
    """Load a VCF or BED file into a mask, returning the intervals and variants.

    BED input yields no :class:`Variant` objects -- it is treated purely as an
    interval mask, which is the right behaviour for repeat tracks or
    hand-curated "avoid this" regions.
    """
    p = Path(path)
    suffixes = {s.lower() for s in p.suffixes}
    if ".bed" in suffixes:
        from ..design.regions import read_bed  # local import avoids a cycle

        mask = read_bed(str(p), chrom=chrom)
        if offset:
            mask = IntervalSet(
                Interval(max(0, i.start - offset), max(0, i.end - offset)) for i in mask
            )
        return mask, []
    variants = read_vcf(
        p,
        chrom=chrom,
        samples=samples,
        heterozygous_only=heterozygous_only,
        offset=offset,
    )
    return variants_to_intervals(variants, pad=pad), variants


def summarise(variants: Iterable[Variant]) -> dict[str, int]:
    """Counts by variant class, for the run log."""
    out: dict[str, int] = {}
    for v in variants:
        out[v.kind] = out.get(v.kind, 0) + 1
    return out
