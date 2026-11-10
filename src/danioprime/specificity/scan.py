"""Finding the places a primer can actually prime.

The question that matters is not "does this primer have a similar sequence
elsewhere" but "can a polymerase extend from it elsewhere". Those are different
questions, and conflating them is why BLAST output is such a poor proxy for
primer specificity:

* A primer with 6 internal mismatches and a perfect 3' end will prime.
* A primer with one mismatch at its 3'-terminal base will not, however well the
  rest of it matches, because the polymerase needs a correctly paired 3'-OH.

So every candidate site is judged on three criteria, reported separately:
mismatches in the 3'-terminal window, total mismatches, and the estimated Tm of
the imperfect duplex.

Two scanning strategies are available. ``exhaustive`` slides the primer across
every position of the subject and is exact; with 3'-first comparison and early
abort it handles sequence sets up to a few megabases comfortably, which covers
paralog panels, transcriptomes and single loci. ``seeded`` uses a k-mer index
and is faster on larger subjects at the cost of sensitivity set by the seed
length. For a whole 1.4 Gb genome, delegate to ``blastn`` or ``bowtie2`` via
:mod:`danioprime.specificity.external`.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

from ..sequence import DNA_BASES, SeqRecord, revcomp
from ..thermo.nn import DEFAULT_CONDITIONS, SaltConditions, imperfect_duplex_tm

if TYPE_CHECKING:  # pragma: no cover
    from .index import KmerIndex

Strand = Literal["+", "-"]

#: Default number of 3'-terminal bases that must match perfectly for a site to
#: be considered extensible. A mismatch within three bases of the 3' end
#: reduces extension efficiency by orders of magnitude.
DEFAULT_CLEAN_3P = 3


@dataclass(slots=True)
class BindingSite:
    """A position where a primer could anneal to a subject sequence.

    ``start``/``end`` are half-open subject coordinates of the footprint.
    ``strand`` is ``+`` when the primer sequence matches the subject's top
    strand (so it extends rightwards) and ``-`` when its reverse complement
    does (extending leftwards).
    """

    subject: str
    start: int
    end: int
    strand: Strand
    mismatches: int
    mismatch_offsets_from_3p: tuple[int, ...]
    three_prime_mismatches: int
    tm_estimate: float
    primer_name: str = ""
    can_prime: bool = False
    note: str = ""

    @property
    def three_prime_pos(self) -> int:
        """Subject position of the primer's 3'-terminal base."""
        return self.end - 1 if self.strand == "+" else self.start

    @property
    def is_perfect(self) -> bool:
        return self.mismatches == 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "subject": self.subject,
            "primer": self.primer_name,
            "start": self.start,
            "end": self.end,
            "position_1based": f"{self.start + 1}-{self.end}",
            "strand": self.strand,
            "mismatches": self.mismatches,
            "mismatch_offsets_from_3p": list(self.mismatch_offsets_from_3p),
            "three_prime_mismatches": self.three_prime_mismatches,
            "tm_estimate": round(self.tm_estimate, 2),
            "can_prime": self.can_prime,
            "note": self.note,
        }


@dataclass(frozen=True, slots=True)
class PrimingRules:
    """When an imperfect site is treated as able to prime."""

    #: Number of 3'-terminal bases that must match perfectly.
    clean_three_prime: int = DEFAULT_CLEAN_3P
    #: Window in which 3'-proximal mismatches are counted for reporting.
    three_prime_window: int = 5
    #: Maximum total mismatches still considered extensible.
    max_mismatches: int = 5
    #: Minimum estimated duplex Tm (degrees C) for the site to prime.
    min_tm: float = 45.0

    def judge(self, mismatches: int, offsets_from_3p: Iterable[int], tm: float) -> tuple[bool, str]:
        offsets = tuple(offsets_from_3p)
        if any(o < self.clean_three_prime for o in offsets):
            return False, f"mismatch within {self.clean_three_prime} nt of the 3' end"
        if mismatches > self.max_mismatches:
            return False, f"{mismatches} mismatches exceeds limit"
        if tm < self.min_tm:
            return False, f"duplex Tm {tm:.1f}C below {self.min_tm:.1f}C"
        return True, ""


DEFAULT_RULES = PrimingRules()


def _scan_one_register(
    query: str,
    subject: str,
    subject_name: str,
    strand: Strand,
    *,
    max_mismatches: int,
    rules: PrimingRules,
    conditions: SaltConditions,
    primer_seq: str,
    primer_name: str,
) -> list[BindingSite]:
    """Slide *query* along *subject*, reporting windows within the mismatch budget.

    *query* is the string whose literal match defines this strand: the primer
    itself for ``+``, its reverse complement for ``-``. Comparison runs from the
    primer's 3' end inwards so that the common case -- a 3' mismatch -- aborts
    after one character.
    """
    sites: list[BindingSite] = []
    lq = len(query)
    ls = len(subject)
    if lq == 0 or ls < lq:
        return sites

    # Index order: iterate query positions nearest the primer's 3' end first.
    if strand == "+":
        order = range(lq - 1, -1, -1)

        def offset_from_3p(i: int) -> int:
            return lq - 1 - i
    else:
        order = range(lq)

        def offset_from_3p(i: int) -> int:
            return i

    budget = max_mismatches
    for pos in range(ls - lq + 1):
        mism = 0
        offsets: list[int] = []
        abort = False
        for i in order:
            s_base = subject[pos + i]
            if s_base == query[i]:
                continue
            if s_base not in DNA_BASES:
                # Treat N in the subject as a mismatch but flag it.
                pass
            mism += 1
            if mism > budget:
                abort = True
                break
            offsets.append(offset_from_3p(i))
        if abort:
            continue
        window = subject[pos : pos + lq]
        template_register = window if strand == "+" else revcomp(window)
        try:
            tm = imperfect_duplex_tm(primer_seq, template_register, conditions)
        except Exception:  # noqa: BLE001 - ambiguous bases in the subject
            tm = float("nan")
        offsets_sorted = tuple(sorted(offsets))
        tp = sum(1 for o in offsets_sorted if o < rules.three_prime_window)
        can_prime, note = rules.judge(mism, offsets_sorted, tm if tm == tm else -999.0)
        sites.append(
            BindingSite(
                subject=subject_name,
                start=pos,
                end=pos + lq,
                strand=strand,
                mismatches=mism,
                mismatch_offsets_from_3p=offsets_sorted,
                three_prime_mismatches=tp,
                tm_estimate=tm,
                primer_name=primer_name,
                can_prime=can_prime,
                note=note,
            )
        )
    return sites


def find_binding_sites(
    primer_seq: str,
    subject: str | SeqRecord,
    *,
    subject_name: str = "subject",
    max_mismatches: int = 5,
    rules: PrimingRules = DEFAULT_RULES,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    primer_name: str = "",
    priming_only: bool = False,
) -> list[BindingSite]:
    """All sites in *subject* (both strands) where *primer_seq* anneals within budget.

    Set *priming_only* to discard sites the polymerase could not extend from;
    that is usually what you want for a specificity verdict, and keeping the
    rest is for diagnostics.
    """
    if isinstance(subject, SeqRecord):
        subject_name = subject.name
        subject_seq = subject.seq
    else:
        subject_seq = subject.upper()
    primer_seq = primer_seq.upper()

    sites = _scan_one_register(
        primer_seq,
        subject_seq,
        subject_name,
        "+",
        max_mismatches=max_mismatches,
        rules=rules,
        conditions=conditions,
        primer_seq=primer_seq,
        primer_name=primer_name,
    )
    sites += _scan_one_register(
        revcomp(primer_seq),
        subject_seq,
        subject_name,
        "-",
        max_mismatches=max_mismatches,
        rules=rules,
        conditions=conditions,
        primer_seq=primer_seq,
        primer_name=primer_name,
    )
    if priming_only:
        sites = [s for s in sites if s.can_prime]
    sites.sort(key=lambda s: (s.mismatches, s.start))
    return sites


@dataclass(slots=True)
class PredictedAmplicon:
    """A product predicted from two priming sites on the same subject."""

    subject: str
    forward_site: BindingSite
    reverse_site: BindingSite
    size: int
    forward_primer: str = ""
    reverse_primer: str = ""
    note: str = ""
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def is_same_primer(self) -> bool:
        return self.forward_primer == self.reverse_primer

    def to_dict(self) -> dict[str, Any]:
        return {
            "subject": self.subject,
            "size": self.size,
            "forward_primer": self.forward_primer,
            "reverse_primer": self.reverse_primer,
            "forward_site": self.forward_site.to_dict(),
            "reverse_site": self.reverse_site.to_dict(),
            "total_mismatches": self.forward_site.mismatches + self.reverse_site.mismatches,
            "same_primer_both_ends": self.is_same_primer,
            "note": self.note,
            **self.extra,
        }


def predict_amplicons(
    sites_by_primer: dict[str, list[BindingSite]],
    *,
    min_size: int,
    max_size: int,
    require_priming: bool = True,
) -> list[PredictedAmplicon]:
    """Pair up priming sites into predicted products.

    Every combination of primers is considered, including a primer with itself,
    because self-pairs (F/F and R/R products from inverted repeats or paralogous
    arrays) are a real source of artefact bands and a pure "is my primer
    unique" check misses them entirely.
    """
    by_subject: dict[str, dict[str, list[BindingSite]]] = {}
    for name, sites in sites_by_primer.items():
        for site in sites:
            if require_priming and not site.can_prime:
                continue
            by_subject.setdefault(site.subject, {}).setdefault(name, []).append(site)

    out: list[PredictedAmplicon] = []
    for subject, per_primer in by_subject.items():
        names = sorted(per_primer)
        for f_name in names:
            for r_name in names:
                for f_site in per_primer[f_name]:
                    if f_site.strand != "+":
                        continue
                    for r_site in per_primer[r_name]:
                        if r_site.strand != "-":
                            continue
                        if r_site.start < f_site.end:
                            continue
                        size = r_site.end - f_site.start
                        if not (min_size <= size <= max_size):
                            continue
                        note = (
                            "self-pair product (same oligo at both ends)"
                            if f_name == r_name
                            else ""
                        )
                        out.append(
                            PredictedAmplicon(
                                subject=subject,
                                forward_site=f_site,
                                reverse_site=r_site,
                                size=size,
                                forward_primer=f_name,
                                reverse_primer=r_name,
                                note=note,
                            )
                        )
    out.sort(key=lambda a: (a.forward_site.mismatches + a.reverse_site.mismatches, a.size))
    return out


def find_binding_sites_seeded(
    primer_seq: str,
    index: KmerIndex,
    *,
    max_mismatches: int = 5,
    rules: PrimingRules = DEFAULT_RULES,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    primer_name: str = "",
    priming_only: bool = True,
) -> list[BindingSite]:
    """Seeded search for priming sites across an indexed subject set.

    Anchors on the primer's 3'-terminal k-mer (and that k-mer's reverse
    complement, for minus-strand sites), then extends 5' counting mismatches.
    See :mod:`danioprime.specificity.index` for the sensitivity contract.
    """
    primer_seq = primer_seq.upper()
    k = index.k
    lp = len(primer_seq)
    if lp < k:
        raise ValueError(f"primer ({lp} nt) is shorter than the seed length ({k})")

    sites: list[BindingSite] = []

    def consider(sid: int, start: int, strand: Strand) -> None:
        record = index.subject(sid)
        seq = record.seq
        if start < 0 or start + lp > len(seq):
            return
        window = seq[start : start + lp]
        query = primer_seq if strand == "+" else revcomp(primer_seq)
        mism = 0
        offsets: list[int] = []
        for i in range(lp):
            if window[i] != query[i]:
                mism += 1
                if mism > max_mismatches:
                    return
                offsets.append(lp - 1 - i if strand == "+" else i)
        template_register = window if strand == "+" else revcomp(window)
        try:
            tm = imperfect_duplex_tm(primer_seq, template_register, conditions)
        except Exception:  # noqa: BLE001
            tm = float("nan")
        offsets_sorted = tuple(sorted(offsets))
        tp = sum(1 for o in offsets_sorted if o < rules.three_prime_window)
        can_prime, note = rules.judge(mism, offsets_sorted, tm if tm == tm else -999.0)
        if priming_only and not can_prime:
            return
        sites.append(
            BindingSite(
                subject=record.name,
                start=start,
                end=start + lp,
                strand=strand,
                mismatches=mism,
                mismatch_offsets_from_3p=offsets_sorted,
                three_prime_mismatches=tp,
                tm_estimate=tm,
                primer_name=primer_name,
                can_prime=can_prime,
                note=note,
            )
        )

    # Plus strand: the primer's 3'-terminal k-mer is the last k bases.
    for sid, pos in index.lookup(primer_seq[-k:]):
        consider(sid, pos - (lp - k), "+")
    # Minus strand: the primer anneals to the top strand, so the subject
    # carries the reverse complement and the seed sits at the window's start.
    for sid, pos in index.lookup(revcomp(primer_seq)[:k]):
        consider(sid, pos, "-")

    seen: set[tuple[str, int, str]] = set()
    unique: list[BindingSite] = []
    for site in sites:
        key = (site.subject, site.start, site.strand)
        if key in seen:
            continue
        seen.add(key)
        unique.append(site)
    unique.sort(key=lambda s: (s.mismatches, s.subject, s.start))
    return unique
