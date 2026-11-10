"""5' tails and the oligos that go with them.

Several zebrafish workflows add fixed sequence to a primer's 5' end. Tails do
not change where a primer anneals, so they are applied after design -- but they
do change what you order, what it costs, and in the three-primer fluorescent
system, how many oligos the reaction needs. Keeping them out of the design and
applying them at the ordering step avoids the common mistake of letting a
tail's own GC content distort the Tm calculation for the annealing portion.

The fluorescent-PCR system (Varshney et al. 2015, *Genome Res.* 25:1030) is the
one most often got wrong: it needs **three** oligos, not two. An M13F-tailed
forward primer, a pig-tailed reverse primer, and a separately ordered
fluorophore-labelled M13F primer that is shared across every assay in the lab.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any, Literal

from ..design.pairing import PrimerPair
from ..errors import DesignError
from ..sequence import gc_percent, normalise

Which = Literal["forward", "reverse", "both"]


@dataclass(frozen=True, slots=True)
class Tail:
    """A fixed 5' extension."""

    name: str
    seq: str
    purpose: str
    #: Companion oligo that must also be ordered, if any.
    companion: str = ""
    companion_note: str = ""

    def __len__(self) -> int:
        return len(self.seq)


TAILS: dict[str, Tail] = {
    "M13F": Tail(
        name="M13F",
        seq="TGTAAAACGACGGCCAGT",
        purpose="universal forward tail for fluorescent PCR and shared-primer sequencing",
        companion="TGTAAAACGACGGCCAGT",
        companion_note=(
            "order this separately with a 5' fluorophore (commonly 6-FAM) and include it "
            "in the reaction as a third primer; one labelled oligo then serves every "
            "M13F-tailed assay in the lab"
        ),
    ),
    "M13R": Tail(
        name="M13R",
        seq="CAGGAAACAGCTATGAC",
        purpose="universal reverse tail for shared-primer sequencing",
        companion="CAGGAAACAGCTATGAC",
        companion_note="order separately if a labelled or sequencing primer is needed",
    ),
    "pigtail": Tail(
        name="pigtail",
        seq="GTTTCTT",
        purpose=(
            "forces complete 3' A addition so fragment-analysis peaks are a single "
            "size rather than a +/-1 nt doublet (Brownstein et al. 1996)"
        ),
    ),
    "T7": Tail(
        name="T7",
        seq="TAATACGACTCACTATAGGG",
        purpose="T7 promoter for in vitro transcription of a riboprobe or sgRNA template",
    ),
    "T3": Tail(
        name="T3",
        seq="AATTAACCCTCACTAAAGGG",
        purpose="T3 promoter for in vitro transcription",
    ),
    "SP6": Tail(
        name="SP6",
        seq="ATTTAGGTGACACTATAG",
        purpose="SP6 promoter for in vitro transcription",
    ),
}

#: Tail combinations for named assay styles.
TAIL_SCHEMES: dict[str, dict[str, str]] = {
    "fluorescent": {"forward": "M13F", "reverse": "pigtail"},
    "m13-sequencing": {"forward": "M13F", "reverse": "M13R"},
    "riboprobe-t7": {"reverse": "T7"},
    "riboprobe-sp6": {"reverse": "SP6"},
    "none": {},
}

SCHEME_DESCRIPTIONS: dict[str, str] = {
    "fluorescent": "M13F-tailed forward + pig-tailed reverse; needs a labelled M13F oligo too",
    "m13-sequencing": "M13F/M13R tails so one pair of sequencing primers serves every amplicon",
    "riboprobe-t7": "T7 promoter on the reverse primer: transcribes an antisense probe",
    "riboprobe-sp6": "SP6 promoter on the reverse primer: transcribes an antisense probe",
    "none": "no tails",
}


def register_tail(name: str, seq: str, purpose: str = "user-supplied") -> Tail:
    """Add a custom tail to the registry."""
    tail = Tail(name=name, seq=normalise(seq, allow_ambiguous=False), purpose=purpose)
    TAILS[name] = tail
    return tail


def get_tail(name: str) -> Tail:
    try:
        return TAILS[name]
    except KeyError as exc:
        raise DesignError(f"unknown tail {name!r}; available: {sorted(TAILS)}") from exc


@dataclass(slots=True)
class OrderItem:
    """One oligo to order."""

    name: str
    sequence: str
    length: int
    role: str
    tail: str = ""
    annealing_seq: str = ""
    gc_percent: float = 0.0
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "sequence": self.sequence,
            "length": self.length,
            "role": self.role,
            "tail": self.tail,
            "annealing_sequence": self.annealing_seq,
            "gc_percent": round(self.gc_percent, 1),
            "note": self.note,
        }


def apply_scheme(pair: PrimerPair, scheme: str = "none") -> PrimerPair:
    """Attach the tails of a named scheme to a pair, in place."""
    if scheme not in TAIL_SCHEMES:
        raise DesignError(
            f"unknown tail scheme {scheme!r}; available: {sorted(TAIL_SCHEMES)}"
        )
    mapping = TAIL_SCHEMES[scheme]
    for role, attr in (("forward", "forward"), ("reverse", "reverse")):
        tail_name = mapping.get(role)
        primer = getattr(pair, attr)
        if tail_name:
            tail = get_tail(tail_name)
            primer.tail = tail.seq
            primer.tail_name = tail.name
        else:
            primer.tail = ""
            primer.tail_name = ""
    return pair


def order_sheet(
    pairs: Iterable[PrimerPair],
    *,
    scheme: str = "none",
    include_companions: bool = True,
) -> list[OrderItem]:
    """Build the list of oligos to order for a set of pairs.

    Companion oligos (the labelled M13F primer, say) are emitted once for the
    whole sheet rather than once per assay, because that is how they are used
    and how they should be costed.
    """
    items: list[OrderItem] = []
    companions: dict[str, Tail] = {}

    for pair in pairs:
        apply_scheme(pair, scheme)
        for role, primer in (("forward", pair.forward), ("reverse", pair.reverse)):
            suffix = "F" if role == "forward" else "R"
            items.append(
                OrderItem(
                    name=f"{pair.name}_{suffix}" if pair.name else f"{role}",
                    sequence=primer.ordered_seq,
                    length=len(primer.ordered_seq),
                    role=role,
                    tail=primer.tail_name,
                    annealing_seq=primer.seq,
                    gc_percent=gc_percent(primer.seq),
                    note=(
                        f"{primer.tail_name} tail; Tm {primer.tm:.1f}C refers to the "
                        f"{len(primer.seq)} nt annealing portion only"
                        if primer.tail_name
                        else f"Tm {primer.tm:.1f}C"
                    ),
                )
            )
            if primer.tail_name:
                tail = get_tail(primer.tail_name)
                if tail.companion:
                    companions[tail.name] = tail

    if include_companions:
        for tail in companions.values():
            items.append(
                OrderItem(
                    name=f"{tail.name}_universal",
                    sequence=tail.companion,
                    length=len(tail.companion),
                    role="companion",
                    tail=tail.name,
                    gc_percent=gc_percent(tail.companion),
                    note=tail.companion_note,
                )
            )
    return items


@dataclass(slots=True)
class ReactionPlan:
    """Human-readable description of how to set the reaction up."""

    scheme: str
    oligos: list[OrderItem] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "scheme": self.scheme,
            "description": SCHEME_DESCRIPTIONS.get(self.scheme, ""),
            "oligos": [o.to_dict() for o in self.oligos],
            "notes": list(self.notes),
        }


def reaction_plan(pairs: Iterable[PrimerPair], scheme: str = "none") -> ReactionPlan:
    """Order sheet plus setup notes for a named scheme."""
    pairs = list(pairs)
    items = order_sheet(pairs, scheme=scheme)
    plan = ReactionPlan(scheme=scheme, oligos=items)
    if scheme == "fluorescent":
        plan.notes.append(
            "Three-primer reaction: use the M13F-tailed forward primer at reduced "
            "concentration relative to the labelled M13F and pig-tailed reverse "
            "primers, so that the tail is incorporated and then extended by the "
            "labelled oligo."
        )
        plan.notes.append(
            "Reported product sizes exclude the tails. Add 18 nt (M13F) and 7 nt "
            "(pig-tail) to predict the fragment-analysis peak position."
        )
    if scheme in {"riboprobe-t7", "riboprobe-sp6"}:
        promoter = "T7" if scheme.endswith("t7") else "SP6"
        plan.notes.append(
            f"The {promoter} promoter is on the reverse primer, so in vitro "
            "transcription from the purified product yields an antisense probe. "
            "Put the promoter on the forward primer instead for a sense control."
        )
    if scheme == "m13-sequencing":
        plan.notes.append(
            "Sequence the product with the universal M13F/M13R primers rather than "
            "the locus-specific ones; the read then starts in known sequence, which "
            "gives cleaner base calls over the first ~30 bases."
        )
    return plan
