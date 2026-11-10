"""Transcript structure from GTF/GFF3.

Only what primer design needs: exon coordinates per transcript, strand, and the
gene each transcript belongs to. Gene and transcript identifiers are read from
both GTF attribute syntax (``gene_id "x";``) and GFF3 syntax (``ID=x;``), so
Ensembl, RefSeq and GENCODE releases all parse without conversion.

Zebrafish annotation note: ZFIN gene symbols distinguish ohnologs with an
``a``/``b`` suffix (``kcnh6a``/``kcnh6b``). :func:`guess_ohnolog_symbol` uses
that convention to suggest the paralog a user probably needs to avoid, which is
a suggestion to confirm against ZFIN rather than an orthology assertion -- the
suffix is a naming convention, not a claim about synteny.
"""

from __future__ import annotations

import gzip
import io
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from ..errors import AnnotationError
from ..sequence import IndexedFasta, SeqRecord, revcomp

_GTF_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
_GFF_ATTR = re.compile(r"([^=;]+)=([^;]*)")


@dataclass(slots=True)
class Exon:
    """One exon, half-open zero-based."""

    start: int
    end: int
    number: int | None = None

    def __len__(self) -> int:
        return self.end - self.start


@dataclass(slots=True)
class Transcript:
    """A transcript and its exon structure."""

    transcript_id: str
    chrom: str
    strand: str
    exons: list[Exon] = field(default_factory=list)
    gene_id: str = ""
    gene_name: str = ""
    biotype: str = ""

    def sorted_exons(self) -> list[Exon]:
        """Exons in transcription order (5'->3' along the transcript)."""
        return sorted(self.exons, key=lambda e: e.start, reverse=self.strand == "-")

    @property
    def spliced_length(self) -> int:
        return sum(len(e) for e in self.exons)

    @property
    def span(self) -> tuple[int, int]:
        return min(e.start for e in self.exons), max(e.end for e in self.exons)

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    def junction_offsets(self) -> tuple[int, ...]:
        """Offsets in the spliced transcript at which each exon-exon junction falls.

        An offset *j* means that transcript bases ``[0, j)`` come from upstream
        exons and base *j* is the first base of the next exon, so a primer
        covering positions either side of *j* cannot anneal to genomic DNA.
        """
        offsets: list[int] = []
        cursor = 0
        exons = self.sorted_exons()
        for exon in exons[:-1]:
            cursor += len(exon)
            offsets.append(cursor)
        return tuple(offsets)

    def spliced_sequence(self, fasta: IndexedFasta) -> SeqRecord:
        """Assemble the mature transcript sequence from a reference FASTA."""
        parts: list[str] = []
        for exon in self.sorted_exons():
            chunk = fasta.fetch(self.chrom, exon.start, exon.end)
            parts.append(revcomp(chunk) if self.strand == "-" else chunk)
        return SeqRecord(
            name=self.transcript_id,
            seq="".join(parts),
            description=(
                f"{self.gene_name or self.gene_id} {self.chrom}:{self.span[0] + 1}-{self.span[1]} "
                f"strand={self.strand} exons={self.exon_count}"
            ),
        )

    def genomic_record(self, fasta: IndexedFasta, *, flank: int = 0) -> SeqRecord:
        """The unspliced locus, with optional flanking sequence."""
        lo, hi = self.span
        lo = max(0, lo - flank)
        hi = hi + flank
        seq = fasta.fetch(self.chrom, lo, hi)
        if self.strand == "-":
            seq = revcomp(seq)
        return SeqRecord(
            name=f"{self.transcript_id}_locus",
            seq=seq,
            description=f"{self.chrom}:{lo + 1}-{hi} strand={self.strand}",
        )


def _open_text(path: str | Path) -> IO[str]:
    path = Path(path)
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return path.open("r")


def _parse_attributes(field_text: str) -> dict[str, str]:
    attrs = {k: v for k, v in _GTF_ATTR.findall(field_text)}
    if attrs:
        return attrs
    for key, value in _GFF_ATTR.findall(field_text):
        attrs[key.strip()] = value.strip()
    return attrs


_TRANSCRIPT_KEYS = ("transcript_id", "Parent", "ID", "transcript")
_GENE_KEYS = ("gene_id", "gene", "geneID")
_NAME_KEYS = ("gene_name", "Name", "gene_symbol")
_BIOTYPE_KEYS = ("transcript_biotype", "gene_biotype", "biotype", "transcript_type")


def _first(attrs: dict[str, str], keys: Iterable[str]) -> str:
    for key in keys:
        value = attrs.get(key)
        if value:
            return value.split(",")[0]
    return ""


def parse_annotation(
    handle: Iterable[str],
    *,
    features: tuple[str, ...] = ("exon",),
    chrom: str | None = None,
) -> dict[str, Transcript]:
    """Build transcripts from GTF/GFF lines."""
    transcripts: dict[str, Transcript] = {}
    for lineno, line in enumerate(handle, 1):
        if not line or line.startswith("#"):
            continue
        f = line.rstrip("\n\r").split("\t")
        if len(f) < 9:
            continue
        seqname, _source, feature, start_s, end_s, _score, strand, _frame, attr_text = f[:9]
        if feature not in features:
            continue
        if chrom is not None and seqname != chrom:
            continue
        try:
            start = int(start_s) - 1
            end = int(end_s)
        except ValueError as exc:
            raise AnnotationError(f"line {lineno}: bad coordinates {start_s!r}-{end_s!r}") from exc
        attrs = _parse_attributes(attr_text)
        tid = _first(attrs, _TRANSCRIPT_KEYS)
        if not tid:
            continue
        tid = tid.removeprefix("transcript:")
        tx = transcripts.get(tid)
        if tx is None:
            tx = Transcript(
                transcript_id=tid,
                chrom=seqname,
                strand=strand if strand in "+-" else "+",
                gene_id=_first(attrs, _GENE_KEYS).removeprefix("gene:"),
                gene_name=_first(attrs, _NAME_KEYS),
                biotype=_first(attrs, _BIOTYPE_KEYS),
            )
            transcripts[tid] = tx
        number = attrs.get("exon_number")
        tx.exons.append(Exon(start=start, end=end, number=int(number) if number else None))
    return transcripts


def read_annotation(path: str | Path, **kwargs) -> dict[str, Transcript]:
    """Read transcripts from a GTF/GFF file (optionally gzipped)."""
    with _open_text(path) as fh:
        txs = parse_annotation(fh, **kwargs)
    if not txs:
        raise AnnotationError(f"no exon features with a transcript identifier found in {path}")
    return txs


def find_transcripts(
    transcripts: dict[str, Transcript], query: str
) -> list[Transcript]:
    """Look up by transcript ID, gene ID or gene symbol (case-insensitive)."""
    q = query.strip().lower()
    exact = [t for t in transcripts.values() if t.transcript_id.lower() == q]
    if exact:
        return exact
    by_gene = [
        t
        for t in transcripts.values()
        if t.gene_id.lower() == q or (t.gene_name and t.gene_name.lower() == q)
    ]
    if by_gene:
        return sorted(by_gene, key=lambda t: (-t.spliced_length, t.transcript_id))
    partial = [
        t
        for t in transcripts.values()
        if q in t.transcript_id.lower() or q in (t.gene_name or "").lower()
    ]
    return sorted(partial, key=lambda t: (-t.spliced_length, t.transcript_id))


def canonical_transcript(transcripts: Iterable[Transcript]) -> Transcript:
    """Pick the longest protein-coding transcript, else simply the longest.

    A heuristic stand-in for a curated canonical set: use the real canonical
    flag from your annotation release when you have one.
    """
    items = list(transcripts)
    if not items:
        raise AnnotationError("no transcripts to choose from")
    coding = [t for t in items if "protein_coding" in (t.biotype or "")]
    pool = coding or items
    return max(pool, key=lambda t: (t.spliced_length, t.exon_count))


def group_by_gene(transcripts: dict[str, Transcript]) -> dict[str, list[Transcript]]:
    out: dict[str, list[Transcript]] = defaultdict(list)
    for tx in transcripts.values():
        out[tx.gene_name or tx.gene_id or tx.transcript_id].append(tx)
    return dict(out)


def guess_ohnolog_symbol(symbol: str) -> str | None:
    """Suggest the ohnolog partner implied by the ZFIN ``a``/``b`` convention.

    ``kcnh6a`` -> ``kcnh6b`` and vice versa. Returns None when the symbol does
    not carry the suffix. This is a naming-convention inference and must be
    confirmed against ZFIN or an orthology database before it is relied on:
    plenty of genes end in 'a' or 'b' without being ohnologs, and plenty of
    ohnolog pairs are not named this way.
    """
    s = symbol.strip()
    if len(s) < 2:
        return None
    stem, last = s[:-1], s[-1]
    if last in "ab" and len(stem) >= 2:
        return stem + ("b" if last == "a" else "a")
    return None


def iter_exon_junction_windows(
    tx: Transcript, *, flank: int = 30
) -> Iterator[tuple[int, int, int]]:
    """Yield ``(junction_offset, window_start, window_end)`` in transcript coordinates.

    Useful for restricting a qPCR primer search to the neighbourhood of
    junctions when a junction-spanning primer is required.
    """
    length = tx.spliced_length
    for offset in tx.junction_offsets():
        yield offset, max(0, offset - flank), min(length, offset + flank)
