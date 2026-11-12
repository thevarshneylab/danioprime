"""Pull a real locus out of Ensembl or UCSC by gene name.

The examples in this repository used to be synthetic. Synthetic data earns its
place in the test suite -- it is deterministic, it runs offline, and it can be
built to contain exactly the hard cases -- but it makes a poor example, because
the first thing anyone does with an example is substitute their own gene, and
nothing in a mock locus tells them whether the result they then get is
believable.

The fix is not to ship a real locus in the repository. A FASTA file of genomic
sequence in a git repository is a copy that silently goes stale, pinned to one
assembly, for one gene, for whoever happened to write the example. The fix is
for the tool to go and get the locus, so the worked example is a command rather
than a file, and the same command works for the reader's own gene.

Two rules shape everything here.

**Report the assembly; never assume it.** Published zebrafish variant sets are
on danRer11 (GRCz11), ZFIN has moved to GRCz12tu, and their coordinates differ
by tens to hundreds of kilobases. So the assembly recorded in the FASTA header
is the one the server reported for the sequence it actually returned, and that
is what the downstream assembly guard reads. Nothing here hard-codes an
assembly name, because the day Ensembl moves is the day a hard-coded one starts
lying.

**Validate, then fail loudly.** A fetch that returns the wrong span, or truncates,
or hands back an HTML error page, must not produce a FASTA file. Every response
is checked against what was asked for -- exact length, nucleotides only -- and a
mismatch raises. The failure mode this rules out is the expensive one: a locus
that looks fine, designs fine, and is the wrong piece of DNA.
"""

from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from datetime import date
from typing import Any, Protocol

from ..errors import AnnotationError, SequenceError
from ..sequence import SeqRecord, revcomp

#: Seconds to wait on any single request.
DEFAULT_TIMEOUT = 30.0

USER_AGENT = "danioprime (+https://github.com/thevarshneylab/danioprime)"

#: Everything a genomic sequence may legitimately contain. ``N`` is expected --
#: assembly gaps are real -- and lower case marks soft-masked repeat.
_ALLOWED = re.compile(r"\A[ACGTNacgtn]+\Z")


def _strip_chr(chrom: str) -> str:
    """``chr18`` -> ``18``, for Ensembl, which does not use the prefix.

    Worth doing rather than leaving to the user: Ensembl answers a prefixed
    name with an empty sequence or a 400, and both are a confusing way to
    discover a naming convention. The mitochondrion is renamed too, since the
    two projects disagree about that one by name rather than by prefix.
    """
    if chrom.lower() in {"chrm", "m"}:
        return "MT"
    return chrom[3:] if chrom[:3].lower() == "chr" and len(chrom) > 3 else chrom


def _add_chr(chrom: str) -> str:
    """``18`` -> ``chr18``, for UCSC, which does use the prefix."""
    if _strip_chr(chrom) == "MT":
        return "chrM"
    return chrom if chrom[:3].lower() == "chr" else f"chr{chrom}"


class Transport(Protocol):
    """Fetches a URL and returns the body as text.

    An interface rather than a direct call to urllib so the providers can be
    driven from recorded responses in the tests. The network is the one part of
    this module that cannot be tested offline; keeping it behind two methods
    means it is also the only part.
    """

    def get(self, url: str, *, timeout: float = DEFAULT_TIMEOUT) -> str: ...


class UrllibTransport:
    """The real thing."""

    __slots__ = ()

    def get(self, url: str, *, timeout: float = DEFAULT_TIMEOUT) -> str:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                charset = response.headers.get_content_charset() or "utf-8"
                return response.read().decode(charset, errors="replace")
        except urllib.error.HTTPError as exc:  # pragma: no cover - needs a server
            detail = ""
            try:
                detail = exc.read().decode("utf-8", errors="replace")[:400].strip()
            except Exception:
                pass
            raise AnnotationError(
                f"{exc.code} {exc.reason} from {url}" + (f"\n{detail}" if detail else "")
            ) from exc
        except urllib.error.URLError as exc:  # pragma: no cover - needs a network
            raise AnnotationError(
                f"could not reach {urllib.parse.urlsplit(url).netloc}: {exc.reason}. "
                "A locus fetch needs network access; on a machine without it, "
                "download the FASTA elsewhere and pass it with --template."
            ) from exc


@dataclass(frozen=True, slots=True)
class Locus:
    """A stretch of sequence, and where it came from.

    The provenance is not decoration. A template whose assembly, coordinates
    and retrieval date are unknown cannot be checked against a variant set, a
    GTF, or a primer ordered from it six months ago.
    """

    seq: str
    assembly: str
    chrom: str
    #: 1-based inclusive, the genome-browser convention.
    start: int
    end: int
    strand: int = 1
    gene: str = ""
    gene_id: str = ""
    description: str = ""
    source: str = ""
    urls: tuple[str, ...] = ()
    retrieved: str = ""
    flank: int = 0
    #: A spliced transcript has no genomic coordinates to report.
    is_transcript: bool = False

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.start:,}-{self.end:,}"

    @property
    def name(self) -> str:
        return self.gene or f"{self.chrom}_{self.start}_{self.end}"

    def header(self) -> str:
        """The FASTA description line: everything needed to retrace this."""
        if self.is_transcript:
            parts = [f"{self.chrom} spliced transcript", f"{len(self.seq):,} nt"]
            if self.assembly:
                parts.append(self.assembly)
        else:
            parts = [
                f"{self.assembly} {self.chrom}:{self.start}-{self.end}",
                f"strand {'+' if self.strand >= 0 else '-'}",
                f"{len(self.seq):,} bp",
            ]
        if self.flank:
            parts.append(f"gene +/- {self.flank:,} bp")
        if self.gene_id:
            parts.append(self.gene_id)
        if self.description:
            parts.append(self.description)
        if self.source:
            parts.append(f"via {self.source}")
        if self.retrieved:
            parts.append(f"retrieved {self.retrieved}")
        return " | ".join(parts)

    def to_record(self) -> SeqRecord:
        return SeqRecord(name=self.name, seq=self.seq.upper(), description=self.header())

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "gene": self.gene,
            "gene_id": self.gene_id,
            "assembly": self.assembly,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "length": len(self.seq),
            "flank": self.flank,
            "source": self.source,
            "urls": list(self.urls),
            "retrieved": self.retrieved,
            "description": self.description,
        }


def _check_sequence(seq: str, *, expected: int | None, url: str) -> str:
    """Reject anything that is not the sequence that was asked for."""
    seq = "".join(seq.split())
    if not seq:
        raise AnnotationError(f"{url} returned no sequence")
    if not _ALLOWED.match(seq):
        bad = sorted({c for c in seq if c not in "ACGTNacgtn"})[:8]
        raise SequenceError(
            f"{url} returned something that is not DNA (unexpected characters: "
            f"{', '.join(repr(c) for c in bad)}). This is usually an error page "
            "rather than a sequence."
        )
    if expected is not None and len(seq) != expected:
        raise AnnotationError(
            f"{url} returned {len(seq):,} bp but {expected:,} bp was requested. "
            "Refusing to write a template whose coordinates do not describe it: "
            "every position downstream -- variant masking, the GTF, the primer "
            "coordinates you would order from -- would be shifted."
        )
    return seq


def _parse_fasta_body(text: str) -> str:
    """Sequence out of a FASTA response, headers dropped."""
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


class EnsemblRest:
    """Ensembl's REST service: resolves gene symbols and returns sequence.

    Chosen as the default because it is the only one of the two that can turn
    ``myo7aa`` into coordinates, and because it reports the assembly it is
    serving, which is the fact this tool most needs and most often lacks.
    """

    name = "Ensembl REST"
    homepage = "https://rest.ensembl.org"

    def __init__(self, base: str = "https://rest.ensembl.org", species: str = "danio_rerio"):
        self.base = base.rstrip("/")
        self.species = species

    def lookup_url(self, symbol: str) -> str:
        quoted = urllib.parse.quote(symbol, safe="")
        return (
            f"{self.base}/lookup/symbol/{self.species}/{quoted}"
            "?content-type=application/json;expand=1"
        )

    def region_url(self, chrom: str, start: int, end: int, strand: int = 1) -> str:
        """*start* and *end* are 1-based inclusive, as Ensembl expects."""
        region = f"{_strip_chr(chrom)}:{start}..{end}:{1 if strand >= 0 else -1}"
        return (
            f"{self.base}/sequence/region/{self.species}/{urllib.parse.quote(region, safe=':.-')}"
            "?content-type=text/x-fasta"
        )

    def transcript_url(self, transcript_id: str, kind: str = "cdna") -> str:
        quoted = urllib.parse.quote(transcript_id, safe="")
        return f"{self.base}/sequence/id/{quoted}?type={kind};content-type=text/x-fasta"

    def lookup(self, symbol: str, transport: Transport, timeout: float) -> dict[str, Any]:
        url = self.lookup_url(symbol)
        body = transport.get(url, timeout=timeout)
        try:
            data = json.loads(body)
        except ValueError as exc:
            raise AnnotationError(f"{url} did not return JSON") from exc
        if not isinstance(data, dict) or "seq_region_name" not in data:
            raise AnnotationError(
                f"Ensembl has no gene {symbol!r} in {self.species}. Check the "
                "symbol (zebrafish ohnologs take an 'a'/'b' suffix, so myo7aa "
                "rather than myo7a), or pass a region instead."
            )
        return data

    def sequence(
        self, chrom: str, start: int, end: int, strand: int, transport: Transport, timeout: float
    ) -> tuple[str, str]:
        url = self.region_url(chrom, start, end, strand)
        return _parse_fasta_body(transport.get(url, timeout=timeout)), url


class UcscApi:
    """UCSC's REST API: sequence by coordinates, for a named assembly.

    No symbol lookup, but it serves assemblies by name -- which is how you get
    danRer11 sequence on purpose, after Ensembl has moved on.
    """

    name = "UCSC REST API"
    homepage = "https://api.genome.ucsc.edu"

    def __init__(self, base: str = "https://api.genome.ucsc.edu", genome: str = "danRer11"):
        self.base = base.rstrip("/")
        self.genome = genome

    _no_symbols = (
        "UCSC cannot look up a gene symbol. Either fetch by region "
        "(--source ucsc chr18:1,000,000-1,010,000), or resolve the symbol "
        "with Ensembl first (--source ensembl)."
    )

    def region_url(self, chrom: str, start: int, end: int, strand: int = 1) -> str:
        """*start* and *end* are 1-based inclusive; UCSC wants 0-based half-open.

        The conversion is the whole reason this method exists rather than a
        format string at the call site: an off-by-one here shifts every primer
        coordinate the run goes on to report.
        """
        query = urllib.parse.urlencode(
            {"genome": self.genome, "chrom": _add_chr(chrom), "start": start - 1, "end": end}
        )
        return f"{self.base}/getData/sequence?{query}"

    def lookup_url(self, symbol: str) -> str:
        del symbol
        raise AnnotationError(self._no_symbols)

    def lookup(self, symbol: str, transport: Transport, timeout: float) -> dict[str, Any]:
        del symbol, transport, timeout
        raise AnnotationError(self._no_symbols)

    def sequence(
        self, chrom: str, start: int, end: int, strand: int, transport: Transport, timeout: float
    ) -> tuple[str, str]:
        url = self.region_url(chrom, start, end, strand)
        body = transport.get(url, timeout=timeout)
        try:
            data = json.loads(body)
        except ValueError as exc:
            raise AnnotationError(f"{url} did not return JSON") from exc
        if "dna" not in data:
            raise AnnotationError(f"{url} returned no sequence: {str(data)[:200]}")
        seq = str(data["dna"])
        if strand < 0:
            # UCSC has no strand parameter, so the flip happens here.
            seq = revcomp(seq.upper())
        return seq, url


class Provider(Protocol):  # pragma: no cover - typing only
    name: str

    def lookup(self, symbol: str, transport: Transport, timeout: float) -> dict[str, Any]: ...
    def sequence(
        self, chrom: str, start: int, end: int, strand: int, transport: Transport, timeout: float
    ) -> tuple[str, str]: ...


PROVIDERS = ("ensembl", "ucsc")


def resolve_provider(
    name: str = "ensembl",
    *,
    species: str = "danio_rerio",
    genome: str = "danRer11",
    base_url: str = "",
) -> Any:
    """Build a provider. *base_url* overrides the public service.

    Useful for a lab running a local Ensembl REST mirror, for one of the public
    mirrors, and for the tests, which run the whole fetch-to-design chain
    against a stub on localhost rather than mocking it out.
    """
    if name == "ensembl":
        return EnsemblRest(base=base_url or "https://rest.ensembl.org", species=species)
    if name == "ucsc":
        return UcscApi(base=base_url or "https://api.genome.ucsc.edu", genome=genome)
    raise AnnotationError(f"unknown source {name!r}; available: {', '.join(PROVIDERS)}")


@dataclass
class FetchPlan:
    """What a fetch is about to do, without doing it.

    Exposed so ``--dry-run`` can show the exact URLs. A user who is about to
    put a locus into a figure legend should be able to see where it came from
    before anything is written, and a user behind a restrictive proxy should be
    able to find out which hosts to allow without guessing.
    """

    urls: list[str] = field(default_factory=list)

    def add(self, url: str) -> str:
        self.urls.append(url)
        return url


def fetch_locus(
    target: str,
    *,
    source: str = "ensembl",
    species: str = "danio_rerio",
    genome: str = "danRer11",
    flank: int = 0,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    strand: int | None = None,
    transport: Transport | None = None,
    timeout: float = DEFAULT_TIMEOUT,
    plan: FetchPlan | None = None,
    base_url: str = "",
) -> Locus:
    """Fetch a gene by symbol, or an explicit region.

    Pass *target* as a gene symbol, or give *chrom*, *start* and *end* (1-based
    inclusive) for a region. *flank* extends a symbol lookup on both sides,
    which is what makes the result usable as a design template: a primer has to
    sit somewhere, and the ends of a gene are not it.
    """
    provider = resolve_provider(source, species=species, genome=genome, base_url=base_url)
    transport = UrllibTransport() if transport is None else transport
    urls: list[str] = []
    gene = gene_id = description = ""
    assembly = ""

    if chrom is None and plan is not None:
        # A symbol's coordinates are only knowable once the lookup has
        # answered, so a dry run shows the request it would make and stops
        # rather than pretending to know the second URL.
        plan.add(provider.lookup_url(target))
        return Locus(
            seq="",
            assembly="",
            chrom="",
            start=0,
            end=0,
            gene=target,
            source=provider.name,
            urls=tuple(plan.urls),
            flank=flank,
        )
    if chrom is None:
        data = provider.lookup(target, transport, timeout)
        urls.append(provider.lookup_url(target))
        gene = str(data.get("display_name") or target)
        gene_id = str(data.get("id") or "")
        description = str(data.get("description") or "").split(" [Source")[0].strip()
        assembly = str(data.get("assembly_name") or "")
        chrom = str(data["seq_region_name"])
        start = int(data["start"]) - flank
        end = int(data["end"]) + flank
        if strand is None:
            # Deliberately *not* data["strand"]: a design template is read in
            # assembly orientation, and flipping it here would silently
            # renumber every coordinate the run reports. A reverse-strand gene
            # is still designed against the plus strand; the tool picks primer
            # orientation itself.
            strand = 1
        if start < 1:
            start = 1
    if start is None or end is None:
        raise AnnotationError("a region fetch needs both a start and an end")
    strand = 1 if strand is None else strand
    if end < start:
        raise AnnotationError(f"end {end} is before start {start}")

    if plan is not None:
        for url in urls:
            plan.add(url)
        plan.add(provider.region_url(chrom, start, end, strand))
        return Locus(
            seq="",
            assembly=assembly or (genome if source == "ucsc" else ""),
            chrom=chrom,
            start=start,
            end=end,
            strand=strand,
            gene=gene,
            gene_id=gene_id,
            description=description,
            source=provider.name,
            urls=tuple(plan.urls),
            flank=flank,
        )

    seq, seq_url = provider.sequence(chrom, start, end, strand, transport, timeout)
    urls.append(seq_url)
    seq = _check_sequence(seq, expected=end - start + 1, url=seq_url)
    if not assembly:
        assembly = genome if source == "ucsc" else "unknown"
    return Locus(
        seq=seq,
        assembly=assembly,
        chrom=chrom,
        start=start,
        end=end,
        strand=strand,
        gene=gene,
        gene_id=gene_id,
        description=description,
        source=provider.name,
        urls=tuple(urls),
        retrieved=date.today().isoformat(),
        flank=flank,
    )


def fetch_transcript(
    symbol_or_id: str,
    *,
    species: str = "danio_rerio",
    kind: str = "cdna",
    transport: Transport | None = None,
    timeout: float = DEFAULT_TIMEOUT,
    plan: FetchPlan | None = None,
    base_url: str = "",
) -> Locus:
    """Fetch a spliced transcript, for qPCR and in situ probe design.

    A gene symbol resolves to its canonical transcript; an ENSDART id is used
    as given. The result carries no genomic coordinates, because a spliced
    sequence has none -- which is also why it is the right template for an
    assay that must not amplify genomic DNA.
    """
    provider = EnsemblRest(base=base_url or "https://rest.ensembl.org", species=species)
    transport = UrllibTransport() if transport is None else transport
    urls: list[str] = []
    transcript_id = symbol_or_id
    gene = gene_id = description = ""
    assembly = ""

    if not re.match(r"\AENS\w*T\d+", symbol_or_id):
        data = provider.lookup(symbol_or_id, transport, timeout)
        urls.append(provider.lookup_url(symbol_or_id))
        gene = str(data.get("display_name") or symbol_or_id)
        gene_id = str(data.get("id") or "")
        description = str(data.get("description") or "").split(" [Source")[0].strip()
        assembly = str(data.get("assembly_name") or "")
        transcripts = data.get("Transcript") or []
        if not transcripts:
            raise AnnotationError(f"Ensembl lists no transcripts for {symbol_or_id!r}")
        chosen = next((t for t in transcripts if t.get("is_canonical")), transcripts[0])
        transcript_id = str(chosen.get("id"))

    url = provider.transcript_url(transcript_id, kind)
    if plan is not None:
        for existing in urls:
            plan.add(existing)
        plan.add(url)
        return Locus(
            seq="",
            assembly=assembly,
            chrom=transcript_id,
            start=1,
            end=0,
            gene=gene or transcript_id,
            gene_id=gene_id,
            description=description,
            source=provider.name,
            urls=tuple(plan.urls),
            is_transcript=True,
        )

    seq = _parse_fasta_body(transport.get(url, timeout=timeout))
    urls.append(url)
    seq = _check_sequence(seq, expected=None, url=url)
    return Locus(
        seq=seq,
        assembly=assembly or "unknown",
        chrom=transcript_id,
        start=1,
        end=len(seq),
        gene=gene or transcript_id,
        gene_id=gene_id,
        description=f"{kind} of {transcript_id}"
        + (f"; {description}" if description else ""),
        source=provider.name,
        urls=tuple(urls),
        retrieved=date.today().isoformat(),
        is_transcript=True,
    )


__all__ = [
    "DEFAULT_TIMEOUT",
    "PROVIDERS",
    "EnsemblRest",
    "FetchPlan",
    "Locus",
    "Provider",
    "Transport",
    "UcscApi",
    "UrllibTransport",
    "fetch_locus",
    "fetch_transcript",
    "resolve_provider",
]
