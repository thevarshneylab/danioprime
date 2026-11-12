"""Region-scoped variant lookup from a tabix-indexed VCF or BED.

This is what makes a published strain-variant resource usable during primer
design. The NHGRI-1 zebrafish VCF is 214 MB; a design needs the few kilobases
around one locus. With the tabix index, that costs about 66 kB of transfer
regardless of how large the file is, so there is no reason to download it.

The index itself (766 kB for NHGRI-1) is cached on disk between runs, because
it is fetched whole and does not change.

Assembly safety
---------------
Published zebrafish variant sets are on danRer11 (GRCz11). ZFIN has moved to
GRCz12tu, and the coordinates are **not** transferable -- shifts of tens to
hundreds of kilobases are normal. Masking a GRCz12tu template with GRCz11
variant coordinates would exclude the wrong bases while reporting that the
primers avoid known polymorphisms, which is worse than not checking at all.
:func:`fetch_region` therefore refuses a mismatch instead of warning about it.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..design.regions import Interval, IntervalSet
from ..errors import AnnotationError
from .bgzf import BgzfReader
from .sources import (
    LocalFile,
    VariantSource,
    cache_dir,
    get_source,
    normalise_assembly,
    open_source,
)
from .tabix import TabixFile, load_index
from .vcf import Variant, variants_to_intervals

#: Columns of a VCF record.
_CHROM, _POS, _ID, _REF, _ALT = 0, 1, 2, 3, 4
_FORMAT_COL = 8


@dataclass(slots=True)
class RegionVariants:
    """Variants found in one region, plus provenance worth reporting."""

    variants: list[Variant] = field(default_factory=list)
    mask: IntervalSet = field(default_factory=IntervalSet)
    source_name: str = ""
    assembly: str = ""
    region: str = ""
    citation: str = ""
    samples: tuple[str, ...] = ()
    heterozygous_only: bool = False
    transfer: dict[str, int] = field(default_factory=dict)
    #: False for an ad-hoc URL, whose assembly danioprime cannot determine.
    assembly_verified: bool = True
    notes: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.variants)

    def summary(self) -> str:
        assembly = self.assembly if self.assembly_verified else f"{self.assembly}, unverified"
        bits = [
            f"{len(self.variants)} variant(s) from {self.source_name} "
            f"({assembly}) in {self.region}"
        ]
        if self.samples:
            bits.append(f"restricted to sample(s) {', '.join(self.samples)}")
        if self.heterozygous_only:
            bits.append("heterozygous positions only")
        if self.transfer.get("bytes_fetched"):
            bits.append(
                f"{self.transfer['bytes_fetched'] / 1000:.0f} kB fetched in "
                f"{self.transfer.get('requests', 0)} request(s)"
            )
        return "; ".join(bits)

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source_name,
            "assembly": self.assembly,
            "assembly_verified": self.assembly_verified,
            "notes": list(self.notes),
            "region": self.region,
            "citation": self.citation,
            "variant_count": len(self.variants),
            "masked_bases": self.mask.total_length,
            "samples": list(self.samples),
            "heterozygous_only": self.heterozygous_only,
            "transfer": dict(self.transfer),
        }


def _cached_index(source: VariantSource | None, index_location: str) -> Path | None:
    """Local path for a cached remote index, downloading it if absent."""
    if not index_location.startswith(("http://", "https://")):
        return None
    digest = hashlib.sha256(index_location.encode()).hexdigest()[:16]
    name = (source.name.replace(" ", "_") + "." if source else "") + digest + ".tbi"
    path = cache_dir() / name
    if path.exists() and path.stat().st_size > 0:
        return path
    remote = open_source(index_location)
    data = remote.read_range(0, remote.size())
    tmp = path.with_suffix(".tbi.part")
    tmp.write_bytes(data)
    tmp.replace(path)
    return path


def open_tabix(
    data_location: str,
    index_location: str | None = None,
    *,
    source: VariantSource | None = None,
    cache_index: bool = True,
) -> tuple[TabixFile, Any]:
    """Open a bgzipped, tabix-indexed file, local or remote.

    Returns the :class:`TabixFile` and the underlying data byte source, so the
    caller can report how much was transferred.
    """
    index_location = index_location or f"{data_location}.tbi"
    if cache_index:
        cached = _cached_index(source, index_location)
        index_src = LocalFile(cached) if cached else open_source(index_location)
    else:
        index_src = open_source(index_location)
    index = load_index(index_src)
    data_src = open_source(data_location)
    return TabixFile(data=BgzfReader(data_src), index=index), data_src


_HOM_REF = {"0/0", "0|0", "0", "./.", ".|.", "."}


def _genotype_fields(header: list[str]) -> list[str]:
    for line in reversed(header):
        if line.startswith("#CHROM"):
            columns = line.lstrip("#").split("\t")
            return columns[9:] if len(columns) > 9 else []
    return []


def _is_variant(gt: str) -> bool:
    call = gt.split(":", 1)[0]
    return call not in _HOM_REF and any(
        a not in {"0", "."} for a in call.replace("|", "/").split("/")
    )


def _is_homozygous_alt(gt: str) -> bool:
    call = gt.split(":", 1)[0].replace("|", "/")
    alleles = [a for a in call.split("/") if a != "."]
    return bool(alleles) and len(set(alleles)) == 1 and alleles[0] not in {"0", "."}


def fetch_region(
    chrom: str,
    start: int,
    end: int,
    *,
    source: str | VariantSource | None = "nhgri1",
    data_location: str | None = None,
    index_location: str | None = None,
    template_assembly: str | None = None,
    samples: Iterable[str] | None = None,
    heterozygous_only: bool = False,
    pad: int = 0,
    offset: int = 0,
    allow_assembly_mismatch: bool = False,
    cache_index: bool = True,
) -> RegionVariants:
    """Variants overlapping ``chrom:[start, end)``, as a maskable interval set.

    Coordinates in and out are zero-based half-open. *offset* shifts results
    into template coordinates, exactly as :func:`danioprime.variants.load_mask`
    does, so a template cut from a contig lines up.

    *template_assembly* is checked against the source's assembly and a mismatch
    raises. Pass *allow_assembly_mismatch* only if you have already lifted the
    coordinates over yourself.
    """
    resolved: VariantSource | None
    notes: list[str] = []
    assembly_verified = True
    if data_location is not None:
        resolved = source if isinstance(source, VariantSource) else None
        data_url = data_location
        index_url = index_location
        source_name = resolved.name if resolved else data_location.rsplit("/", 1)[-1]
        citation = resolved.citation if resolved else ""
        if resolved:
            assembly = resolved.assembly
        else:
            # An arbitrary URL carries no assembly declaration that danioprime
            # can read, so claiming one would be a false assurance. Report it
            # as unknown and say so, rather than echoing back whatever the user
            # declared as though it had been checked.
            assembly = "unknown"
            assembly_verified = False
            notes.append(
                "danioprime cannot determine the assembly of an ad-hoc URL, so the "
                "coordinate check was skipped. Make sure this dataset is on the same "
                "assembly as your template: GRCz11 and GRCz12 coordinates differ by "
                "tens to hundreds of kilobases."
            )
    else:
        if source is None:
            raise AnnotationError("supply a named source, or a data location")
        resolved = source if isinstance(source, VariantSource) else get_source(str(source))
        data_url, index_url = resolved.data_url, resolved.index_url
        assembly, source_name, citation = resolved.assembly, resolved.name, resolved.citation

    if template_assembly and assembly_verified:
        wanted = normalise_assembly(template_assembly)
        have = normalise_assembly(assembly)
        if wanted != have and not allow_assembly_mismatch:
            raise AnnotationError(
                f"{source_name} is on {have}, but the template coordinates were "
                f"declared as {wanted}. These assemblies' coordinates are not "
                f"transferable -- shifts of tens to hundreds of kilobases are "
                f"normal between GRCz11 and GRCz12 -- so masking would exclude "
                f"the wrong bases while reporting that known variants were "
                f"avoided.\n"
                f"Either supply a {have} template (for example with --region "
                f"against a {have} reference), lift your coordinates over to "
                f"{have} first, or pass --allow-assembly-mismatch if you have "
                f"already done so."
            )

    tabix, data_src = open_tabix(
        data_url, index_url, source=resolved, cache_index=cache_index
    )
    wanted_samples = tuple(samples or ())
    sample_names = _genotype_fields(tabix.header()) if wanted_samples else []
    if wanted_samples:
        missing = [s for s in wanted_samples if s not in sample_names]
        if missing:
            raise AnnotationError(
                f"sample(s) {missing} are not in {source_name}; it has "
                f"{sample_names or 'no genotype columns'}"
            )
        indices = [sample_names.index(s) for s in wanted_samples]
    else:
        indices = []

    found: list[Variant] = []
    for fields in tabix.fetch(chrom, start, end):
        if len(fields) <= _ALT:
            continue
        try:
            pos0 = int(fields[_POS]) - 1 - offset
        except ValueError:
            continue
        if pos0 < 0:
            continue
        ref = fields[_REF].upper()
        alts = tuple(a for a in fields[_ALT].upper().split(",") if a not in {".", ""})
        genotypes = tuple(fields[_FORMAT_COL + 1 :]) if len(fields) > _FORMAT_COL else ()

        if indices:
            relevant = [genotypes[i] for i in indices if i < len(genotypes)]
            if not any(_is_variant(g) for g in relevant):
                continue
            if heterozygous_only and not any(
                _is_variant(g) and not _is_homozygous_alt(g) for g in relevant
            ):
                continue
        elif heterozygous_only and genotypes:
            if not any(_is_variant(g) and not _is_homozygous_alt(g) for g in genotypes):
                continue

        found.append(
            Variant(
                chrom=fields[_CHROM],
                pos=pos0,
                ref=ref,
                alt=alts,
                identifier=fields[_ID] if len(fields) > _ID else ".",
                genotypes=genotypes,
            )
        )

    return RegionVariants(
        variants=found,
        mask=variants_to_intervals(found, pad=pad),
        source_name=source_name,
        assembly=assembly,
        region=f"{chrom}:{start + 1}-{end}",
        citation=citation,
        samples=wanted_samples,
        heterozygous_only=heterozygous_only,
        transfer=getattr(data_src, "stats", {}) or {},
        assembly_verified=assembly_verified,
        notes=notes,
    )


def fetch_for_template(
    length: int,
    chrom: str,
    template_start: int,
    **kwargs: Any,
) -> RegionVariants:
    """Variants over a template that was cut from *chrom* at *template_start*.

    Convenience wrapper: converts the template's genomic span into a query and
    shifts the results into template coordinates in one step.
    """
    kwargs.setdefault("offset", template_start)
    return fetch_region(chrom, template_start, template_start + length, **kwargs)


def describe_sources() -> list[dict[str, str]]:
    """Registry contents, for ``danioprime sources``."""
    from .sources import UCSC_TRACK_HUB, VARIANT_SOURCES

    rows = [
        {
            "key": key,
            "name": src.name,
            "assembly": src.assembly,
            "kind": src.kind,
            "description": src.description,
            "citation": src.citation,
            "homepage": src.homepage,
            "data_url": src.data_url,
        }
        for key, src in sorted(VARIANT_SOURCES.items())
    ]
    for row in rows:
        row["track_hub"] = UCSC_TRACK_HUB
    return rows


def clear_index_cache() -> int:
    """Delete cached indexes. Returns how many files were removed."""
    removed = 0
    for path in cache_dir().glob("*.tbi"):
        path.unlink()
        removed += 1
    return removed


__all__ = [
    "Interval",
    "RegionVariants",
    "clear_index_cache",
    "describe_sources",
    "fetch_for_template",
    "fetch_region",
    "open_tabix",
]
