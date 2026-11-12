"""BGZF: block-gzip random access, in pure Python.

A ``.vcf.gz`` produced by ``bgzip`` is not an ordinary gzip file. It is a
sequence of independently inflatable gzip members, each holding at most 64 KiB
of uncompressed data and carrying an extra field that records its own
compressed length. That is what makes random access possible: with an index
giving the offset of the block containing a record, a reader can fetch and
inflate that block alone.

Positions in BGZF are **virtual offsets**: a 64-bit value whose high 48 bits
are the byte offset of the block in the file and whose low 16 bits are the
offset of the datum within that block once inflated. Tabix indexes store
chunks as pairs of virtual offsets.

Why implement this rather than depend on pysam: the NHGRI-1 zebrafish variant
VCF is 214 MB and the whole point is to read a few kilobases of it without
downloading the rest. pysam would do that, but it is a compiled extension, and
this package has no required dependencies -- which is also what lets the same
code run in the browser build. The format is small and fully specified, so the
cost of implementing it is bounded. htslib remains the reference: the test
suite checks this reader against files written by ``bgzip`` and queries
answered by pysam.

Reference: SAM/BAM specification, section on BGZF compression.
"""

from __future__ import annotations

import struct
import zlib
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Protocol

from ..errors import AnnotationError

#: Gzip member header, fixed part: magic, CM, FLG, MTIME, XFL, OS, XLEN.
_HEADER = struct.Struct("<4BI2BH")
_HEADER_SIZE = _HEADER.size  # 12
#: Largest a single BGZF block can be on disk.
MAX_BLOCK_SIZE = 65_536
#: The 28-byte empty block that marks a well-formed end of file.
EOF_MARKER = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def virtual_offset(block_start: int, within: int) -> int:
    """Pack a (block offset, in-block offset) pair into a virtual offset."""
    if not 0 <= within < MAX_BLOCK_SIZE:
        raise ValueError(f"in-block offset {within} out of range")
    return (block_start << 16) | within


def split_virtual_offset(voffset: int) -> tuple[int, int]:
    """Unpack a virtual offset into (block offset, in-block offset)."""
    return voffset >> 16, voffset & 0xFFFF


class ByteSource(Protocol):
    """Anything that can serve byte ranges: a local file or an HTTP server."""

    def read_range(self, start: int, length: int) -> bytes: ...
    def size(self) -> int: ...


@dataclass(slots=True)
class Block:
    """One inflated BGZF block and where it came from."""

    block_start: int
    compressed_size: int
    data: bytes

    @property
    def next_block_start(self) -> int:
        return self.block_start + self.compressed_size


def parse_block(raw: bytes, block_start: int = 0) -> Block:
    """Inflate a single BGZF block from *raw*, which must begin at its header."""
    if len(raw) < _HEADER_SIZE:
        raise AnnotationError("truncated BGZF block header")
    m1, m2, cm, flg, _mtime, _xfl, _os, xlen = _HEADER.unpack_from(raw, 0)
    if (m1, m2) != (0x1F, 0x8B):
        raise AnnotationError(
            f"not a gzip/BGZF block at offset {block_start} "
            f"(magic {m1:#04x}{m2:#04x}); is the file bgzip-compressed?"
        )
    if cm != 8:
        raise AnnotationError(f"unsupported compression method {cm}")
    if not flg & 0x04:
        raise AnnotationError(
            "gzip member has no extra field, so this is plain gzip rather than "
            "BGZF. Random access needs a bgzip-compressed file; recompress with "
            "`bgzip` and index with `tabix`."
        )

    # Walk the extra subfields looking for BC, which holds the block size - 1.
    extra = raw[_HEADER_SIZE : _HEADER_SIZE + xlen]
    if len(extra) < xlen:
        raise AnnotationError("truncated BGZF extra field")
    bsize: int | None = None
    pos = 0
    while pos + 4 <= len(extra):
        si1, si2, slen = struct.unpack_from("<2BH", extra, pos)
        payload = extra[pos + 4 : pos + 4 + slen]
        if (si1, si2) == (0x42, 0x43) and slen == 2:
            bsize = struct.unpack("<H", payload)[0] + 1
        pos += 4 + slen
    if bsize is None:
        raise AnnotationError("BGZF block is missing its BC size field")
    if len(raw) < bsize:
        raise AnnotationError(
            f"need {bsize} bytes for the block at {block_start}, have {len(raw)}"
        )

    deflate_start = _HEADER_SIZE + xlen
    deflate_end = bsize - 8  # trailing CRC32 and ISIZE
    payload = raw[deflate_start:deflate_end]
    crc, isize = struct.unpack("<II", raw[deflate_end:bsize])
    data = zlib.decompress(payload, -zlib.MAX_WBITS)
    if len(data) != isize:
        raise AnnotationError(
            f"BGZF block at {block_start} inflated to {len(data)} bytes, "
            f"header says {isize}"
        )
    if zlib.crc32(data) != crc:
        raise AnnotationError(f"BGZF block at {block_start} failed its CRC check")
    return Block(block_start=block_start, compressed_size=bsize, data=data)


class BgzfReader:
    """Random access to a BGZF file over any :class:`ByteSource`.

    Blocks are cached, because a tabix query typically lands several chunks in
    the same block and refetching it over HTTP would be wasteful.
    """

    __slots__ = ("_source", "_cache", "_cache_order", "_max_cached")

    def __init__(self, source: ByteSource, *, max_cached_blocks: int = 64):
        self._source = source
        self._cache: dict[int, Block] = {}
        self._cache_order: list[int] = []
        self._max_cached = max_cached_blocks

    def block_at(self, block_start: int) -> Block:
        """The block beginning at *block_start*, inflated."""
        cached = self._cache.get(block_start)
        if cached is not None:
            return cached
        # A block is at most 64 KiB on disk, so one ranged read always suffices.
        raw = self._source.read_range(block_start, MAX_BLOCK_SIZE)
        if not raw:
            raise AnnotationError(f"no data at offset {block_start}")
        block = parse_block(raw, block_start)
        self._remember(block)
        return block

    def _remember(self, block: Block) -> None:
        self._cache[block.block_start] = block
        self._cache_order.append(block.block_start)
        if len(self._cache_order) > self._max_cached:
            self._cache.pop(self._cache_order.pop(0), None)

    def read_chunk(self, begin: int, end: int) -> bytes:
        """Uncompressed bytes between two virtual offsets.

        The whole byte span of the chunk is fetched in **one** ranged read and
        the blocks decoded from that buffer, rather than one request per block.
        Over HTTP that is the difference between one round trip and dozens: a
        chunk spanning twenty blocks used to cost twenty requests and fetch
        64 KiB for each, most of it discarded.
        """
        if end < begin:
            raise ValueError("chunk end precedes its begin")
        begin_block, begin_in = split_virtual_offset(begin)
        end_block, end_in = split_virtual_offset(end)
        if begin_block == end_block:
            data = self.block_at(begin_block).data
            return data[begin_in:end_in]

        # One read covering every block in the chunk. The last block must be
        # complete, so extend by the maximum block size; that is the only
        # over-read, and it is bounded by 64 KiB regardless of chunk length.
        span_end = min(end_block + MAX_BLOCK_SIZE, self._source.size())
        buffer = self._source.read_range(begin_block, span_end - begin_block)

        out: list[bytes] = []
        offset = begin_block
        while offset <= end_block:
            local = offset - begin_block
            if local >= len(buffer):
                break
            cached = self._cache.get(offset)
            block = cached if cached is not None else parse_block(buffer[local:], offset)
            if cached is None:
                self._remember(block)
            start = begin_in if offset == begin_block else 0
            stop = end_in if offset == end_block else len(block.data)
            out.append(block.data[start:stop])
            if not block.compressed_size:
                break
            offset = block.next_block_start
        return b"".join(out)

    def iter_lines(self, begin: int, end: int) -> Iterator[bytes]:
        """Whole lines within a virtual-offset chunk.

        The final line may be truncated when the chunk boundary falls mid-line;
        tabix chunks start on a record boundary, so only the tail is at risk and
        callers re-filter by position anyway.
        """
        data = self.read_chunk(begin, end)
        for line in data.split(b"\n"):
            if line:
                yield line

    def read_header_lines(self, prefix: bytes = b"#", limit: int = 4_000_000) -> list[bytes]:
        """Leading comment lines, read from the start of the file."""
        lines: list[bytes] = []
        offset = 0
        consumed = 0
        pending = b""
        while consumed < limit:
            block = self.block_at(offset)
            if not block.data:
                break
            pending += block.data
            consumed += len(block.data)
            *complete, pending = pending.split(b"\n")
            for line in complete:
                if line.startswith(prefix):
                    lines.append(line)
                else:
                    return lines
            offset = block.next_block_start
            if offset >= self._source.size():
                break
        return lines
