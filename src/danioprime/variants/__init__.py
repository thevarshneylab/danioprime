"""Strain-variant handling: local files and published remote resources."""

from .remote import (
    RegionVariants,
    clear_index_cache,
    describe_sources,
    fetch_for_template,
    fetch_region,
    open_tabix,
)
from .sources import (
    ASSEMBLY_ALIASES,
    UCSC_TRACK_HUB,
    VARIANT_SOURCES,
    HttpFile,
    LocalFile,
    VariantSource,
    get_source,
    normalise_assembly,
    open_source,
)
from .tabix import TabixFile, TabixIndex, load_index, reg2bins
from .vcf import Variant, load_mask, parse_vcf, read_vcf, summarise, variants_to_intervals

__all__ = [
    "ASSEMBLY_ALIASES",
    "UCSC_TRACK_HUB",
    "VARIANT_SOURCES",
    "HttpFile",
    "LocalFile",
    "RegionVariants",
    "TabixFile",
    "TabixIndex",
    "Variant",
    "VariantSource",
    "clear_index_cache",
    "describe_sources",
    "fetch_for_template",
    "fetch_region",
    "get_source",
    "load_index",
    "load_mask",
    "normalise_assembly",
    "open_source",
    "open_tabix",
    "parse_vcf",
    "read_vcf",
    "reg2bins",
    "summarise",
    "variants_to_intervals",
]
