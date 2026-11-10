"""Restriction enzyme recognition sites for CAPS and dCAPS design.

A curated subset of common, inexpensive, widely stocked enzymes. Sites are
written 5'->3' on the top strand in IUPAC notation; ``cut_offset`` is the number
of bases from the start of the recognition sequence to the top-strand cut.

Scope and caveats
-----------------
* Type IIS enzymes that cut at a distance from their recognition sequence
  (FokI, BsaI, BsmBI, ...) are **not** included: dCAPS fragment prediction
  assumes the cut falls inside the recognition site.
* Sites with long internal N runs (BstXI, PflMI) are excluded for the same
  reason.
* Always confirm the site, buffer and incubation temperature against your
  supplier's current catalogue before ordering. Isoschizomers differ in
  methylation sensitivity, which this table does not model -- relevant because
  Dam/Dcm methylation of a PCR product is not an issue but of plasmid template
  is.

Supply your own table with ``--enzymes FILE`` (TSV: ``name  site  cut_offset``)
to use an enzyme not listed here.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..errors import AnnotationError
from ..sequence import IUPAC


@dataclass(frozen=True, slots=True)
class Enzyme:
    """A restriction enzyme and where it cuts."""

    name: str
    site: str
    cut_offset: int

    @property
    def length(self) -> int:
        return len(self.site)

    @property
    def is_degenerate(self) -> bool:
        return any(b not in "ACGT" for b in self.site)

    @property
    def display(self) -> str:
        return f"{self.site[: self.cut_offset]}^{self.site[self.cut_offset :]}"


#: name, recognition site (IUPAC, top strand), top-strand cut offset
_TABLE: tuple[tuple[str, str, int], ...] = (
    ("AluI", "AGCT", 2),
    ("ApaI", "GGGCCC", 5),
    ("AseI", "ATTAAT", 2),
    ("AvrII", "CCTAGG", 1),
    ("BamHI", "GGATCC", 1),
    ("BclI", "TGATCA", 1),
    ("BglII", "AGATCT", 1),
    ("BspHI", "TCATGA", 1),
    ("BsrGI", "TGTACA", 1),
    ("BstEII", "GGTNACC", 1),
    ("ClaI", "ATCGAT", 2),
    ("DdeI", "CTNAG", 1),
    ("DpnI", "GATC", 2),
    ("DraI", "TTTAAA", 3),
    ("EagI", "CGGCCG", 1),
    ("EcoRI", "GAATTC", 1),
    ("EcoRV", "GATATC", 3),
    ("Fnu4HI", "GCNGC", 2),
    ("FspI", "TGCGCA", 3),
    ("HaeIII", "GGCC", 2),
    ("HhaI", "GCGC", 3),
    ("HincII", "GTYRAC", 3),
    ("HindIII", "AAGCTT", 1),
    ("HinfI", "GANTC", 1),
    ("HpaI", "GTTAAC", 3),
    ("KpnI", "GGTACC", 5),
    ("MboI", "GATC", 0),
    ("MfeI", "CAATTG", 1),
    ("MluI", "ACGCGT", 1),
    ("MscI", "TGGCCA", 3),
    ("MseI", "TTAA", 1),
    ("MspI", "CCGG", 1),
    ("NcoI", "CCATGG", 1),
    ("NdeI", "CATATG", 2),
    ("NheI", "GCTAGC", 1),
    ("NlaIII", "CATG", 4),
    ("NotI", "GCGGCCGC", 2),
    ("NruI", "TCGCGA", 3),
    ("NsiI", "ATGCAT", 5),
    ("NspI", "RCATGY", 5),
    ("PciI", "ACATGT", 1),
    ("PstI", "CTGCAG", 5),
    ("PvuI", "CGATCG", 4),
    ("PvuII", "CAGCTG", 3),
    ("RsaI", "GTAC", 2),
    ("SacI", "GAGCTC", 5),
    ("SalI", "GTCGAC", 1),
    ("ScaI", "AGTACT", 3),
    ("SmaI", "CCCGGG", 3),
    ("SnaBI", "TACGTA", 3),
    ("SpeI", "ACTAGT", 1),
    ("SphI", "GCATGC", 5),
    ("SspI", "AATATT", 3),
    ("StuI", "AGGCCT", 3),
    ("StyI", "CCWWGG", 1),
    ("TaqI", "TCGA", 1),
    ("Tsp509I", "AATT", 0),
    ("XbaI", "TCTAGA", 1),
    ("XhoI", "CTCGAG", 1),
    ("XmaI", "CCCGGG", 1),
)

#: Default enzyme set, keyed by name.
DEFAULT_ENZYMES: dict[str, Enzyme] = {
    name: Enzyme(name=name, site=site, cut_offset=offset) for name, site, offset in _TABLE
}


def load_enzymes(path: str) -> dict[str, Enzyme]:
    """Read a TSV of ``name  site  cut_offset`` into an enzyme table."""
    table: dict[str, Enzyme] = {}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                raise AnnotationError(f"{path}:{lineno}: need at least name and site")
            name, site = parts[0], parts[1].upper()
            offset = int(parts[2]) if len(parts) > 2 else 0
            bad = [b for b in site if b not in IUPAC]
            if bad:
                raise AnnotationError(f"{path}:{lineno}: bad IUPAC code(s) {bad} in {site!r}")
            if not 0 <= offset <= len(site):
                raise AnnotationError(
                    f"{path}:{lineno}: cut offset {offset} outside the site (0-{len(site)})"
                )
            table[name] = Enzyme(name=name, site=site, cut_offset=offset)
    if not table:
        raise AnnotationError(f"no enzymes read from {path}")
    return table


def site_matches(site: str, seq: str) -> bool:
    """IUPAC-aware comparison of a recognition site against a sequence."""
    if len(site) != len(seq):
        return False
    for code, base in zip(site, seq, strict=True):
        allowed = IUPAC.get(code)
        if allowed is None or base not in allowed:
            return False
    return True


def find_sites(enzyme: Enzyme, seq: str) -> list[int]:
    """Start offsets of every occurrence of *enzyme*'s site in *seq*."""
    k = enzyme.length
    return [i for i in range(len(seq) - k + 1) if site_matches(enzyme.site, seq[i : i + k])]


def digest_fragments(enzyme: Enzyme, seq: str) -> list[int]:
    """Fragment lengths from a complete digest of linear *seq*."""
    cuts = sorted(i + enzyme.cut_offset for i in find_sites(enzyme, seq))
    cuts = [c for c in cuts if 0 < c < len(seq)]
    if not cuts:
        return [len(seq)]
    sizes: list[int] = []
    prev = 0
    for cut in cuts:
        sizes.append(cut - prev)
        prev = cut
    sizes.append(len(seq) - prev)
    return sizes


def concretise(site: str) -> str:
    """A representative unambiguous sequence for a degenerate site.

    Used when a dCAPS primer must *write* a degenerate position: the first base
    of the IUPAC set is chosen, which is arbitrary but deterministic.
    """
    return "".join(sorted(IUPAC[b])[0] if b not in "ACGT" else b for b in site.upper())
