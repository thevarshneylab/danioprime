"""Nucleotide sequence primitives: IUPAC handling, FASTA I/O, composition stats.

Sequences are handled as plain upper-case ASCII strings. Soft-masked (lower-case)
input from reference genomes is preserved as a separate mask so that repeat
content can be penalised without losing the underlying bases.
"""

from __future__ import annotations

import gzip
import io
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from .errors import SequenceError

#: Unambiguous DNA bases.
DNA_BASES = frozenset("ACGT")

#: IUPAC ambiguity codes mapped to the set of bases they represent.
IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)

_NON_IUPAC = re.compile(r"[^ACGTURYSWKMBDHVN]")


def normalise(seq: str, *, allow_ambiguous: bool = True) -> str:
    """Upper-case *seq*, strip whitespace, and validate the alphabet."""
    cleaned = "".join(seq.split()).upper().replace("U", "T")
    bad = _NON_IUPAC.findall(cleaned)
    if bad:
        raise SequenceError(f"sequence contains non-IUPAC characters: {sorted(set(bad))}")
    if not allow_ambiguous:
        amb = set(cleaned) - DNA_BASES
        if amb:
            raise SequenceError(f"sequence contains ambiguity codes: {sorted(amb)}")
    return cleaned


def revcomp(seq: str) -> str:
    """Reverse complement, preserving case and IUPAC codes."""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    """Fraction of G+C among unambiguous bases. Returns 0.0 for an empty count."""
    counted = [b for b in seq.upper() if b in DNA_BASES]
    if not counted:
        return 0.0
    gc = sum(1 for b in counted if b in "GC")
    return gc / len(counted)


def gc_percent(seq: str) -> float:
    """G+C content as a percentage."""
    return 100.0 * gc_fraction(seq)


def max_homopolymer(seq: str) -> int:
    """Length of the longest single-base run (``poly-X``)."""
    best = run = 0
    prev = ""
    for base in seq.upper():
        run = run + 1 if base == prev else 1
        prev = base
        best = max(best, run)
    return best


def max_dinucleotide_repeat(seq: str) -> int:
    """Number of tandem copies in the longest dinucleotide repeat (e.g. ``ATATAT`` -> 3)."""
    s = seq.upper()
    best = 0
    for i in range(len(s) - 1):
        unit = s[i : i + 2]
        if unit[0] == unit[1]:
            continue
        copies = 1
        j = i + 2
        while s[j : j + 2] == unit:
            copies += 1
            j += 2
        best = max(best, copies)
    return best


def ambiguous_positions(seq: str) -> list[int]:
    """Zero-based offsets of bases that are not A/C/G/T."""
    return [i for i, b in enumerate(seq.upper()) if b not in DNA_BASES]


def iupac_matches(code: str, base: str) -> bool:
    """True if unambiguous *base* is compatible with IUPAC *code*."""
    allowed = IUPAC.get(code.upper())
    if allowed is None:
        raise SequenceError(f"unknown IUPAC code: {code!r}")
    return base.upper() in allowed


def soft_mask_runs(raw: str, *, min_length: int = 1) -> list[tuple[int, int]]:
    """Half-open intervals of lower-case (soft-masked) bases in *raw*.

    Reference FASTA from Ensembl/UCSC marks interspersed repeats in lower case;
    those intervals are useful as a repeat penalty for primer placement.
    """
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for i, ch in enumerate(raw):
        if ch.islower():
            if start is None:
                start = i
        elif start is not None:
            if i - start >= min_length:
                runs.append((start, i))
            start = None
    if start is not None and len(raw) - start >= min_length:
        runs.append((start, len(raw)))
    return runs


def annotate_tail(tail: str, annealing: str) -> str:
    """Join a 5' tail and an annealing sequence, the tail in lower case.

    The convention is the one oligo order forms and papers already use: bases
    that anneal are upper case, bases that were bolted on are not. Note that
    lower case means something different here than it does in a template --
    there it marks interspersed repeat (see :func:`soft_mask_runs`), here it
    marks sequence that is not in the template at all. The two never meet: one
    describes a locus, the other an oligo, and nothing reads an oligo back as
    if it were a locus.
    """
    return tail.lower() + annealing.upper()


@dataclass(slots=True)
class SeqRecord:
    """A single FASTA record.

    Attributes
    ----------
    name:
        First whitespace-delimited token of the header.
    description:
        Remainder of the header line.
    seq:
        Upper-case, validated sequence.
    softmask:
        Half-open intervals that were lower-case in the source file.
    """

    name: str
    seq: str
    description: str = ""
    softmask: list[tuple[int, int]] = field(default_factory=list)

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self.seq)

    def subseq(self, start: int, end: int) -> str:
        """Half-open slice with bounds checking."""
        if start < 0 or end > len(self.seq) or start > end:
            raise SequenceError(f"slice [{start},{end}) out of range for {self.name} ({len(self.seq)} nt)")
        return self.seq[start:end]


def _open_text(path: str | Path) -> IO[str]:
    path = Path(path)
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return path.open("r")


def parse_fasta(handle: Iterable[str]) -> Iterator[SeqRecord]:
    """Stream FASTA records from any iterable of lines."""
    header: str | None = None
    chunks: list[str] = []

    def build(hdr: str, parts: list[str]) -> SeqRecord:
        raw = "".join(parts)
        name, _, desc = hdr.partition(" ")
        return SeqRecord(
            name=name,
            seq=normalise(raw),
            description=desc.strip(),
            softmask=soft_mask_runs("".join(raw.split())),
        )

    for line in handle:
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield build(header, chunks)
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise SequenceError("FASTA data encountered before any '>' header")
            chunks.append(line)
    if header is not None:
        yield build(header, chunks)


def read_fasta(path: str | Path) -> list[SeqRecord]:
    """Read every record from a (optionally gzipped) FASTA file."""
    with _open_text(path) as fh:
        records = list(parse_fasta(fh))
    if not records:
        raise SequenceError(f"no FASTA records found in {path}")
    return records


def read_single_fasta(path: str | Path) -> SeqRecord:
    """Read a FASTA file that is expected to hold exactly one record."""
    records = read_fasta(path)
    if len(records) != 1:
        raise SequenceError(f"expected 1 record in {path}, found {len(records)}")
    return records[0]


#: Strings longer than this are never probed as filesystem paths: most systems
#: cap a path component at 255 bytes and ``Path.exists()`` raises rather than
#: returning False, so a pasted sequence would crash the lookup.
MAX_PATH_PROBE = 255


def load_sequence(spec: str | Path) -> SeqRecord:
    """Accept either a FASTA path or a bare nucleotide string.

    A path that exists wins, so a file named ``ACGT`` is read as a file. When
    the argument is too long to be a path, or probing the filesystem fails for
    any OS-level reason, it is treated as a literal sequence.
    """
    text = str(spec)
    compact = "".join(text.split()).upper()
    looks_like_sequence = bool(compact) and not _NON_IUPAC.search(compact)

    if len(text) <= MAX_PATH_PROBE:
        try:
            candidate = Path(text)
            if candidate.exists():
                return read_single_fasta(candidate)
        except OSError:
            pass

    if looks_like_sequence:
        return SeqRecord(name="input", seq=normalise(compact))
    raise SequenceError(f"{text!r} is neither an existing file nor a nucleotide sequence")


def write_fasta(records: Iterable[SeqRecord], path: str | Path, *, wrap: int = 60) -> None:
    """Write records to a FASTA file."""
    # newline="\n" so output is byte-identical on every platform: text mode
    # would emit CRLF on Windows, and the example dataset is regenerated and
    # diffed in CI.
    with Path(path).open("w", newline="\n") as fh:
        for rec in records:
            header = f">{rec.name}" + (f" {rec.description}" if rec.description else "")
            fh.write(header + "\n")
            for i in range(0, len(rec.seq), wrap):
                fh.write(rec.seq[i : i + wrap] + "\n")


class IndexedFasta:
    """Random access to a FASTA file via its ``.fai`` index.

    Reference genomes are too large to hold in memory, but exon extraction only
    needs a few kilobases at a time. The index is the standard ``samtools faidx``
    format; it is written on first use if absent, which takes one linear pass
    and only works on a FASTA with fixed-width lines (as every reference release
    is).
    """

    __slots__ = ("path", "_index", "_fh")

    def __init__(self, path: str | Path, *, build_index: bool = True):
        self.path = Path(path)
        if self.path.suffix == ".gz":
            raise SequenceError(
                f"{self.path.name} is gzipped; indexed access needs an uncompressed or "
                "bgzipped-with-.gzi reference. Decompress it, or pass sequences directly."
            )
        if not self.path.exists():
            raise SequenceError(f"reference FASTA not found: {self.path}")
        fai = self.path.with_name(self.path.name + ".fai")
        if not fai.exists():
            if not build_index:
                raise SequenceError(f"missing index {fai}; run `samtools faidx {self.path}`")
            self._write_index(fai)
        self._index = self._read_index(fai)
        self._fh: IO[bytes] | None = None

    def _write_index(self, fai: Path) -> None:
        records: list[tuple[str, int, int, int, int]] = []
        with self.path.open("rb") as fh:
            name: str | None = None
            offset = 0
            length = 0
            line_bases = 0
            line_width = 0
            pos = 0
            for raw in fh:
                if raw.startswith(b">"):
                    if name is not None:
                        records.append((name, length, offset, line_bases, line_width))
                    name = raw[1:].split()[0].decode() if len(raw) > 1 else ""
                    offset = pos + len(raw)
                    length = 0
                    line_bases = 0
                    line_width = 0
                else:
                    stripped = raw.rstrip(b"\r\n")
                    if line_bases == 0:
                        line_bases = len(stripped)
                        line_width = len(raw)
                    length += len(stripped)
                pos += len(raw)
            if name is not None:
                records.append((name, length, offset, line_bases, line_width))
        with fai.open("w", newline="\n") as out:
            for rec in records:
                out.write("\t".join(str(x) for x in rec) + "\n")

    @staticmethod
    def _read_index(fai: Path) -> dict[str, tuple[int, int, int, int]]:
        index: dict[str, tuple[int, int, int, int]] = {}
        with fai.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                index[parts[0]] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
        return index

    @property
    def contigs(self) -> list[str]:
        return list(self._index)

    def length(self, contig: str) -> int:
        try:
            return self._index[contig][0]
        except KeyError as exc:
            raise SequenceError(
                f"contig {contig!r} not in {self.path.name}; first few are {self.contigs[:5]}"
            ) from exc

    def fetch(self, contig: str, start: int, end: int) -> str:
        """Half-open, zero-based slice of *contig*."""
        if contig not in self._index:
            raise SequenceError(
                f"contig {contig!r} not in {self.path.name}; first few are {self.contigs[:5]}"
            )
        clen, offset, line_bases, line_width = self._index[contig]
        start = max(0, start)
        end = min(clen, end)
        if start >= end:
            return ""
        if self._fh is None:
            self._fh = self.path.open("rb")
        first_line, first_col = divmod(start, line_bases)
        byte_start = offset + first_line * line_width + first_col
        last_line, last_col = divmod(end - 1, line_bases)
        byte_end = offset + last_line * line_width + last_col + 1
        self._fh.seek(byte_start)
        raw = self._fh.read(byte_end - byte_start)
        return normalise(raw.decode("ascii", errors="replace").replace("\n", "").replace("\r", ""))

    def record(self, contig: str, start: int, end: int, *, name: str | None = None) -> SeqRecord:
        seq = self.fetch(contig, start, end)
        return SeqRecord(
            name=name or f"{contig}:{start + 1}-{end}",
            seq=seq,
            description=f"{self.path.name} {contig}:{start + 1}-{end}",
        )

    def close(self) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    def __enter__(self) -> IndexedFasta:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
