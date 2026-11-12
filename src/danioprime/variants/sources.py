"""Byte sources, and the named remote variant datasets.

Two concerns live here. First, the :class:`~danioprime.variants.bgzf.ByteSource`
implementations: a local file, and an HTTP server that honours range requests.
Second, a small registry of published zebrafish variant resources, so a user
does not have to remember a URL or which assembly it is on.

The assembly is the part that bites. Every dataset here is on **danRer11 /
GRCz11**, because that is the assembly the zebrafish community has been using.
ZFIN has since moved to GRCz12tu, whose coordinates are *not* transferable:
shifts of tens to hundreds of kilobases are normal. Applying GRCz11 variant
coordinates to a GRCz12tu template would mask the wrong bases and report
confident nonsense, so the design code refuses to do it rather than warning.
"""

from __future__ import annotations

import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from ..errors import AnnotationError

#: How long to wait on any single HTTP request, in seconds.
DEFAULT_TIMEOUT = 30.0

USER_AGENT = "danioprime (+https://github.com/thevarshneylab/danioprime)"


class LocalFile:
    """A :class:`ByteSource` over a file on disk."""

    __slots__ = ("path", "_size")

    def __init__(self, path: str | Path):
        self.path = Path(path)
        if not self.path.exists():
            raise AnnotationError(f"file not found: {self.path}")
        self._size = self.path.stat().st_size

    def read_range(self, start: int, length: int) -> bytes:
        with self.path.open("rb") as fh:
            fh.seek(start)
            return fh.read(length)

    def size(self) -> int:
        return self._size

    def __repr__(self) -> str:  # pragma: no cover
        return f"LocalFile({self.path})"


class HttpFile:
    """A :class:`ByteSource` over HTTP, using Range requests.

    This is what makes a 214 MB remote VCF usable: a tabix query needs a few
    tens of kilobytes of it. The server must support ``Range``; most do, and
    one that does not is detected on the first read rather than silently
    returning the whole file.
    """

    __slots__ = ("url", "_size", "_timeout", "_bytes_fetched", "_requests")

    def __init__(self, url: str, *, timeout: float = DEFAULT_TIMEOUT):
        self.url = url
        self._timeout = timeout
        self._bytes_fetched = 0
        self._requests = 0
        self._size = self._probe_size()

    def _open(self, headers: dict[str, str]):
        request = urllib.request.Request(  # noqa: S310 - https URLs from a fixed registry or the user
            self.url, headers={"User-Agent": USER_AGENT, **headers}
        )
        return urllib.request.urlopen(request, timeout=self._timeout)  # noqa: S310

    def _probe_size(self) -> int:
        try:
            request = urllib.request.Request(  # noqa: S310
                self.url, method="HEAD", headers={"User-Agent": USER_AGENT}
            )
            with urllib.request.urlopen(request, timeout=self._timeout) as response:  # noqa: S310
                length = response.headers.get("Content-Length")
                accepts = (response.headers.get("Accept-Ranges") or "").lower()
            if length is None:
                raise AnnotationError(
                    f"{self.url} did not report a Content-Length, so it cannot be "
                    f"read in ranges. Download the file and pass it locally."
                )
            if accepts == "none":
                raise AnnotationError(
                    f"{self.url} declares it does not accept range requests. "
                    f"Download the file and pass it locally instead."
                )
            return int(length)
        except urllib.error.HTTPError as exc:
            raise AnnotationError(
                f"{self.url} returned HTTP {exc.code} ({exc.reason}). Check the URL; "
                f"published data locations do move."
            ) from exc
        except urllib.error.URLError as exc:
            raise AnnotationError(
                f"could not reach {self.url}: {exc.reason}. "
                f"This needs outbound network access."
            ) from exc

    def read_range(self, start: int, length: int) -> bytes:
        if length <= 0:
            return b""
        stop = min(start + length, self._size) - 1
        if stop < start:
            return b""
        try:
            with self._open({"Range": f"bytes={start}-{stop}"}) as response:
                if response.status == 200 and self._size > (stop - start + 1):
                    raise AnnotationError(
                        f"{self.url} ignored the Range header and returned the whole "
                        f"file. Download it and pass it locally instead."
                    )
                data = response.read()
        except urllib.error.HTTPError as exc:
            raise AnnotationError(
                f"range request to {self.url} failed with HTTP {exc.code}"
            ) from exc
        except urllib.error.URLError as exc:
            raise AnnotationError(f"range request to {self.url} failed: {exc.reason}") from exc
        self._requests += 1
        self._bytes_fetched += len(data)
        return data

    def size(self) -> int:
        return self._size

    @property
    def stats(self) -> dict[str, int]:
        """Requests made and bytes transferred, for reporting."""
        return {
            "requests": self._requests,
            "bytes_fetched": self._bytes_fetched,
            "file_size": self._size,
        }

    def __repr__(self) -> str:  # pragma: no cover
        return f"HttpFile({self.url})"


def open_source(location: str | Path, *, timeout: float = DEFAULT_TIMEOUT):
    """Open a local path or an http(s) URL as a byte source."""
    text = str(location)
    if text.startswith(("http://", "https://")):
        return HttpFile(text, timeout=timeout)
    return LocalFile(text)


# --------------------------------------------------------------------------
# Published zebrafish variant resources
# --------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class VariantSource:
    """A published variant dataset that can be queried by region."""

    name: str
    description: str
    assembly: str
    data_url: str
    index_url: str
    kind: str  # "vcf" or "bed"
    citation: str
    homepage: str
    #: True when the intervals mark sequence *confirmed to match* the
    #: reference, rather than sequence that differs from it.
    invariant: bool = False

    @property
    def strain(self) -> str:
        return self.name


_BURGESS = "https://research.nhgri.nih.gov/manuscripts/Burgess/zebrafish"
_NHGRI1_TRACKS = f"{_BURGESS}/downloads/NHGRI-1/danRer11/danRer11Tracks"

#: Named sources, keyed by the value accepted on the command line.
#:
#: URLs are as published in the primary literature and on the Burgess lab
#: download page. Data locations do move; a failure names the homepage so the
#: current link can be checked, and ``--variants`` always accepts a URL or a
#: local file directly, so a stale entry here never blocks the capability.
VARIANT_SOURCES: dict[str, VariantSource] = {
    "nhgri1": VariantSource(
        name="NHGRI-1",
        description=(
            "Variants in the NHGRI-1 line relative to the reference, from founder "
            "parents sequenced to ~50x. The line exists precisely so that the "
            "sequence of a given region, with all its variants, is known."
        ),
        assembly="danRer11",
        data_url=f"{_NHGRI1_TRACKS}/NHGRI1.danRer11.variant.vcf.gz",
        index_url=f"{_NHGRI1_TRACKS}/NHGRI1.danRer11.variant.vcf.gz.tbi",
        kind="vcf",
        citation=(
            "LaFave MC, Varshney GK, Vemulapalli M, Mullikin JC, Burgess SM (2014) "
            "A defined zebrafish line for high-throughput genetics and genomics: "
            "NHGRI-1. Genetics 198:167-170. doi:10.1534/genetics.114.166769"
        ),
        homepage=f"{_BURGESS}/download.shtml",
    ),
}

#: UCSC track hub carrying the same data, for browsing rather than querying.
UCSC_TRACK_HUB = f"{_BURGESS}/downloads/NHGRI-1/hub.txt"

#: Assemblies whose coordinates are mutually incompatible. Used to refuse a
#: cross-assembly mask rather than silently masking the wrong bases.
ASSEMBLY_ALIASES: dict[str, str] = {
    "danrer11": "danRer11",
    "grcz11": "danRer11",
    "danrer10": "danRer10",
    "grcz10": "danRer10",
    "grcz12tu": "GRCz12tu",
    "danrer12": "GRCz12tu",
    "grcz12ab": "GRCz12ab",
}


def normalise_assembly(name: str) -> str:
    """Canonical name for an assembly, accepting UCSC or GRC spellings."""
    key = (name or "").strip().lower().replace("_", "").replace("-", "")
    return ASSEMBLY_ALIASES.get(key, name.strip())


def get_source(name: str) -> VariantSource:
    try:
        return VARIANT_SOURCES[name.strip().lower()]
    except KeyError as exc:
        raise AnnotationError(
            f"unknown variant source {name!r}; available: "
            f"{', '.join(sorted(VARIANT_SOURCES))}. You can also pass a URL or a "
            f"local path to any bgzip-compressed, tabix-indexed VCF or BED."
        ) from exc


def cache_dir() -> Path:
    """Where downloaded indexes are kept between runs."""
    base = os.environ.get("DANIOPRIME_CACHE") or os.environ.get("XDG_CACHE_HOME")
    root = Path(base) if base else Path.home() / ".cache"
    path = root / "danioprime"
    path.mkdir(parents=True, exist_ok=True)
    return path
