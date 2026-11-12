"""Fetching a real locus by gene name."""

from .remote import (
    PROVIDERS,
    EnsemblRest,
    Locus,
    Provider,
    UcscApi,
    UrllibTransport,
    fetch_locus,
    fetch_transcript,
    resolve_provider,
)

__all__ = [
    "PROVIDERS",
    "EnsemblRest",
    "Locus",
    "Provider",
    "UcscApi",
    "UrllibTransport",
    "fetch_locus",
    "fetch_transcript",
    "resolve_provider",
]
