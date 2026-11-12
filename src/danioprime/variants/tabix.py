"""Tabix indexes, in pure Python.

A ``.tbi`` is a BGZF-compressed binary index over a position-sorted, bgzipped
tab-delimited file. It lets a reader turn "chromosome 5, 1.2 Mb to 1.3 Mb" into
a short list of byte ranges, so a 214 MB VCF can be queried for one locus over
HTTP without downloading it.

The scheme is the UCSC binning index. The genome is cut into a hierarchy of
bins -- one bin covering 512 Mb, then 8 of 64 Mb, 64 of 8 Mb, 512 of 1 Mb, 4096
of 128 kb, 32768 of 16 kb -- and every record is registered in the smallest bin
that wholly contains it. A query collects the bins at every level that overlap
the region, which is why a handful of bin lookups suffices. Each bin holds
*chunks*: pairs of BGZF virtual offsets. A **linear index** then gives, for
each 16 kb window, the smallest virtual offset of any record overlapping it,
which prunes chunks that provably end before the query starts.

Reference: Li H. (2011) Tabix: fast retrieval of sequence features from generic
TAB-delimited files. *Bioinformatics* 27:718-719, and the SAM/BAM/tabix
specification.
"""

from __future__ import annotations

import struct
from collections.abc import Iterator
from dataclasses import dataclass, field

from ..errors import AnnotationError
from .bgzf import BgzfReader, ByteSource, parse_block

TBI_MAGIC = b"TBI\x01"
CSI_MAGIC = b"CSI\x01"

#: Size of the smallest bin, and so of each linear-index window.
MIN_SHIFT = 14  # 16 kb
#: Number of levels in the bin hierarchy.
DEPTH = 5

#: Format codes in the index header.
FORMAT_GENERIC = 0
FORMAT_SAM = 1
FORMAT_VCF = 2
#: Set in the format field when coordinates are zero-based half-open (BED).
FLAG_ZERO_BASED = 0x10000


def reg2bins(begin: int, end: int) -> list[int]:
    """Every bin that could hold a record overlapping ``[begin, end)``.

    Zero-based half-open input. Mirrors the reference implementation's
    ``reg2bins`` for a 5-level, 16 kb-leaf hierarchy.
    """
    if end <= begin:
        return []
    end -= 1
    bins = [0]
    for level_shift, offset in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        bins.extend(range(offset + (begin >> level_shift), offset + (end >> level_shift) + 1))
    return bins


@dataclass(slots=True)
class Chunk:
    """A half-open range of BGZF virtual offsets holding candidate records."""

    begin: int
    end: int


@dataclass(slots=True)
class RefIndex:
    """Index for one reference sequence."""

    bins: dict[int, list[Chunk]] = field(default_factory=dict)
    #: Smallest virtual offset of any record overlapping each 16 kb window.
    intervals: list[int] = field(default_factory=list)

    def linear_floor(self, begin: int) -> int:
        """Lowest virtual offset worth reading for a query starting at *begin*."""
        window = begin >> MIN_SHIFT
        if not self.intervals:
            return 0
        if window >= len(self.intervals):
            return self.intervals[-1]
        return self.intervals[window]


@dataclass(slots=True)
class TabixIndex:
    """A parsed ``.tbi``."""

    format: int
    col_seq: int
    col_begin: int
    col_end: int
    meta_char: str
    skip_lines: int
    reference_names: list[str] = field(default_factory=list)
    refs: list[RefIndex] = field(default_factory=list)

    @property
    def zero_based(self) -> bool:
        return bool(self.format & FLAG_ZERO_BASED)

    @property
    def base_format(self) -> int:
        return self.format & 0xFFFF

    def reference_id(self, name: str) -> int | None:
        """Index of *name*, tolerating a missing or extra ``chr`` prefix."""
        try:
            return self.reference_names.index(name)
        except ValueError:
            pass
        alternatives = (
            name[3:] if name.lower().startswith("chr") else f"chr{name}",
        )
        for alt in alternatives:
            try:
                return self.reference_names.index(alt)
            except ValueError:
                continue
        return None

    def chunks_for(self, name: str, begin: int, end: int) -> list[Chunk]:
        """Merged, pruned byte ranges to read for a region.

        *begin* and *end* are zero-based half-open.
        """
        ref_id = self.reference_id(name)
        if ref_id is None:
            raise AnnotationError(
                f"{name!r} is not in this index. It names "
                f"{len(self.reference_names)} sequence(s), the first few being "
                f"{self.reference_names[:5]}. Chromosome naming differs between "
                f"sources: Ensembl uses '5', UCSC uses 'chr5'."
            )
        ref = self.refs[ref_id]
        floor = ref.linear_floor(begin)
        candidates: list[Chunk] = []
        for bin_id in reg2bins(begin, end):
            for chunk in ref.bins.get(bin_id, ()):
                # The linear index proves nothing before `floor` can overlap.
                if chunk.end > floor:
                    candidates.append(Chunk(max(chunk.begin, floor), chunk.end))
        if not candidates:
            return []
        candidates.sort(key=lambda c: c.begin)
        merged = [candidates[0]]
        for chunk in candidates[1:]:
            last = merged[-1]
            if chunk.begin <= last.end:
                last.end = max(last.end, chunk.end)
            else:
                merged.append(chunk)
        return merged


def _read_u32(buf: bytes, pos: int) -> tuple[int, int]:
    return struct.unpack_from("<I", buf, pos)[0], pos + 4


def parse_tabix(raw: bytes) -> TabixIndex:
    """Parse a decompressed ``.tbi`` payload."""
    if raw[:4] == CSI_MAGIC:
        raise AnnotationError(
            "this is a CSI index, not a tabix index. CSI is used for references "
            "longer than 512 Mb; danioprime reads .tbi only. Re-index with "
            "`tabix -p vcf` if the reference allows it."
        )
    if raw[:4] != TBI_MAGIC:
        raise AnnotationError(
            f"not a tabix index (magic {raw[:4]!r}); expected a .tbi file"
        )
    pos = 4
    n_ref, pos = _read_u32(raw, pos)
    fmt, pos = _read_u32(raw, pos)
    col_seq, pos = _read_u32(raw, pos)
    col_begin, pos = _read_u32(raw, pos)
    col_end, pos = _read_u32(raw, pos)
    meta, pos = _read_u32(raw, pos)
    skip, pos = _read_u32(raw, pos)
    l_nm, pos = _read_u32(raw, pos)

    names_blob = raw[pos : pos + l_nm]
    pos += l_nm
    names = [n.decode("ascii") for n in names_blob.split(b"\x00") if n]

    refs: list[RefIndex] = []
    for _ in range(n_ref):
        ref = RefIndex()
        n_bin, pos = _read_u32(raw, pos)
        for _ in range(n_bin):
            bin_id, pos = _read_u32(raw, pos)
            n_chunk, pos = _read_u32(raw, pos)
            chunks = []
            for _ in range(n_chunk):
                begin, end = struct.unpack_from("<QQ", raw, pos)
                pos += 16
                chunks.append(Chunk(begin, end))
            ref.bins[bin_id] = chunks
        n_intv, pos = _read_u32(raw, pos)
        if n_intv:
            ref.intervals = list(struct.unpack_from(f"<{n_intv}Q", raw, pos))
            pos += 8 * n_intv
        refs.append(ref)

    if len(names) != n_ref:
        raise AnnotationError(
            f"index declares {n_ref} references but names {len(names)}"
        )
    return TabixIndex(
        format=fmt,
        col_seq=col_seq,
        col_begin=col_begin,
        col_end=col_end,
        meta_char=chr(meta) if meta else "#",
        skip_lines=skip,
        reference_names=names,
        refs=refs,
    )


def load_index(source: ByteSource) -> TabixIndex:
    """Read and parse a ``.tbi`` from a byte source, inflating every block."""
    total = source.size()
    raw = source.read_range(0, total)
    chunks: list[bytes] = []
    offset = 0
    while offset < len(raw):
        block = parse_block(raw[offset:], offset)
        if not block.compressed_size:
            break
        chunks.append(block.data)
        offset += block.compressed_size
    if not chunks:
        raise AnnotationError("tabix index is empty or not BGZF-compressed")
    return parse_tabix(b"".join(chunks))


@dataclass(slots=True)
class TabixFile:
    """A bgzipped, tabix-indexed tab-delimited file, queryable by region."""

    data: BgzfReader
    index: TabixIndex

    def header(self) -> list[str]:
        meta = self.index.meta_char.encode()
        return [line.decode("utf-8", "replace") for line in self.data.read_header_lines(meta)]

    def fetch(self, name: str, begin: int, end: int) -> Iterator[list[str]]:
        """Records overlapping ``[begin, end)``, zero-based half-open.

        Chunks are a superset of the answer -- a bin can hold records that miss
        the query -- so every candidate is re-checked against its own
        coordinates here rather than trusted.
        """
        if end <= begin:
            return
        meta = self.index.meta_char
        col_seq = self.index.col_seq - 1
        col_begin = self.index.col_begin - 1
        col_end = self.index.col_end - 1
        zero_based = self.index.zero_based
        is_vcf = self.index.base_format == FORMAT_VCF

        for chunk in self.index.chunks_for(name, begin, end):
            for raw_line in self.data.iter_lines(chunk.begin, chunk.end):
                line = raw_line.decode("utf-8", "replace")
                if not line or line.startswith(meta):
                    continue
                fields = line.split("\t")
                if len(fields) <= max(col_seq, col_begin):
                    continue
                if fields[col_seq] != name and self.index.reference_id(
                    fields[col_seq]
                ) != self.index.reference_id(name):
                    continue
                try:
                    start = int(fields[col_begin])
                except ValueError:
                    continue
                if not zero_based:
                    start -= 1
                if col_end > 0 and col_end < len(fields) and not is_vcf:
                    try:
                        stop = int(fields[col_end])
                    except ValueError:
                        stop = start + 1
                elif is_vcf:
                    # A VCF record spans the length of its REF allele.
                    stop = start + max(1, len(fields[3])) if len(fields) > 3 else start + 1
                else:
                    stop = start + 1
                if start < end and begin < stop:
                    yield fields
