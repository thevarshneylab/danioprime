"""Ohnolog-aware specificity: keeping a primer off its own paralog.

The teleost-specific whole-genome duplication left zebrafish with thousands of
retained ohnolog pairs, conventionally named ``geneA``/``geneB``. Coding
sequence identity between them is often high enough that a primer pair designed
against one gene model amplifies both, and because the paralog is a real gene
rather than a repeat, nothing in a standard primer-design workflow flags it.

This module answers two questions a generic tool does not ask:

1. Can either primer of a pair prime on the paralog at all?
2. If both can, is a product of amplifiable size predicted?

It also rewards the design trick experimentalists use by hand: placing the
primer's 3'-terminal bases over a position where the two ohnologs differ, so
that discrimination comes from extension failure rather than from a Tm
difference that a permissive annealing temperature will wash out.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any

from ..sequence import SeqRecord, read_fasta
from ..thermo.nn import DEFAULT_CONDITIONS, SaltConditions
from .index import DEFAULT_SEED, KmerIndex
from .scan import (
    DEFAULT_RULES,
    BindingSite,
    PredictedAmplicon,
    PrimingRules,
    find_binding_sites,
    find_binding_sites_seeded,
    predict_amplicons,
)

#: Subject sets larger than this (bases) switch from exhaustive to seeded
#: scanning unless the caller overrides the strategy.
EXHAUSTIVE_LIMIT = 3_000_000


@dataclass(slots=True)
class ParalogPanel:
    """Sequences a primer must *not* amplify.

    Typically the ohnolog of the target gene, but the same machinery serves any
    "avoid this sequence" requirement: a closely related family member, a
    transgene backbone, a pseudogene, or the unedited allele.
    """

    records: list[SeqRecord] = field(default_factory=list)
    strategy: str = "auto"
    seed_length: int = DEFAULT_SEED
    _index: KmerIndex | None = field(default=None, repr=False)

    @classmethod
    def from_fasta(cls, paths: str | Iterable[str], **kwargs: Any) -> ParalogPanel:
        """Load one or more FASTA files into a panel."""
        if isinstance(paths, str):
            paths = [paths]
        records: list[SeqRecord] = []
        for path in paths:
            records.extend(read_fasta(path))
        return cls(records=records, **kwargs)

    @classmethod
    def from_records(cls, records: Sequence[SeqRecord], **kwargs: Any) -> ParalogPanel:
        return cls(records=list(records), **kwargs)

    def __bool__(self) -> bool:
        return bool(self.records)

    def __len__(self) -> int:
        return len(self.records)

    @property
    def names(self) -> list[str]:
        return [r.name for r in self.records]

    @property
    def total_length(self) -> int:
        return sum(len(r.seq) for r in self.records)

    def resolved_strategy(self) -> str:
        if self.strategy != "auto":
            return self.strategy
        return "exhaustive" if self.total_length <= EXHAUSTIVE_LIMIT else "seeded"

    def index(self) -> KmerIndex:
        if self._index is None:
            self._index = KmerIndex.build(self.records, k=self.seed_length)
        return self._index

    def scan(
        self,
        primer_seq: str,
        *,
        primer_name: str = "",
        rules: PrimingRules = DEFAULT_RULES,
        conditions: SaltConditions = DEFAULT_CONDITIONS,
        max_mismatches: int = 5,
        priming_only: bool = False,
    ) -> list[BindingSite]:
        """All binding sites for one primer across the panel."""
        if self.resolved_strategy() == "seeded":
            sites = find_binding_sites_seeded(
                primer_seq,
                self.index(),
                max_mismatches=max_mismatches,
                rules=rules,
                conditions=conditions,
                primer_name=primer_name,
                priming_only=priming_only,
            )
        else:
            sites = []
            for record in self.records:
                sites.extend(
                    find_binding_sites(
                        primer_seq,
                        record,
                        max_mismatches=max_mismatches,
                        rules=rules,
                        conditions=conditions,
                        primer_name=primer_name,
                        priming_only=priming_only,
                    )
                )
        return sites


@dataclass(slots=True)
class Discrimination:
    """How well one primer distinguishes its target from the panel."""

    primer_name: str
    primer_seq: str
    best_site: BindingSite | None
    priming_sites: list[BindingSite] = field(default_factory=list)
    all_sites: list[BindingSite] = field(default_factory=list)

    @property
    def primes_on_panel(self) -> bool:
        return bool(self.priming_sites)

    @property
    def best_identity_mismatches(self) -> int | None:
        return self.best_site.mismatches if self.best_site else None

    @property
    def discriminates_at_3p(self) -> bool:
        """True when the closest panel match carries a 3'-proximal mismatch.

        This is the robust form of discrimination: the polymerase cannot extend,
        irrespective of annealing temperature.
        """
        if self.best_site is None:
            return True
        return self.best_site.three_prime_mismatches > 0

    def verdict(self) -> str:
        if self.best_site is None:
            return "no detectable binding site on the panel"
        if self.primes_on_panel:
            return (
                f"primes on {self.priming_sites[0].subject} "
                f"({self.priming_sites[0].mismatches} mismatches)"
            )
        if self.discriminates_at_3p:
            offsets = self.best_site.mismatch_offsets_from_3p
            nearest = min(offsets) + 1 if offsets else None
            return f"discriminated by a mismatch {nearest} nt from the 3' end"
        return "binds the panel but is not predicted to extend"

    def to_dict(self) -> dict[str, Any]:
        return {
            "primer": self.primer_name,
            "primes_on_panel": self.primes_on_panel,
            "discriminates_at_3p": self.discriminates_at_3p,
            "best_match_mismatches": self.best_identity_mismatches,
            "verdict": self.verdict(),
            "priming_sites": [s.to_dict() for s in self.priming_sites],
            "closest_site": self.best_site.to_dict() if self.best_site else None,
        }


def assess_primer(
    primer_seq: str,
    panel: ParalogPanel,
    *,
    primer_name: str = "",
    rules: PrimingRules = DEFAULT_RULES,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    max_mismatches: int = 6,
) -> Discrimination:
    """Characterise one primer against the panel."""
    sites = panel.scan(
        primer_seq,
        primer_name=primer_name,
        rules=rules,
        conditions=conditions,
        max_mismatches=max_mismatches,
        priming_only=False,
    )
    priming = [s for s in sites if s.can_prime]
    best = min(sites, key=lambda s: (s.mismatches, -s.tm_estimate)) if sites else None
    return Discrimination(
        primer_name=primer_name or primer_seq,
        primer_seq=primer_seq,
        best_site=best,
        priming_sites=priming,
        all_sites=sites,
    )


def assess_pair(
    forward_seq: str,
    reverse_seq: str,
    panel: ParalogPanel,
    *,
    min_size: int,
    max_size: int,
    forward_name: str = "F",
    reverse_name: str = "R",
    rules: PrimingRules = DEFAULT_RULES,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    max_mismatches: int = 6,
) -> tuple[list[PredictedAmplicon], dict[str, Discrimination]]:
    """Predict paralog products for a primer pair.

    Returns the predicted amplicons and the per-primer discrimination reports.
    Self-pairs (the same oligo priming both ends) are included, because a
    conserved inverted arrangement between ohnologs can produce exactly that.
    """
    fwd = assess_primer(
        forward_seq,
        panel,
        primer_name=forward_name,
        rules=rules,
        conditions=conditions,
        max_mismatches=max_mismatches,
    )
    rev = assess_primer(
        reverse_seq,
        panel,
        primer_name=reverse_name,
        rules=rules,
        conditions=conditions,
        max_mismatches=max_mismatches,
    )
    amplicons = predict_amplicons(
        {forward_name: fwd.all_sites, reverse_name: rev.all_sites},
        min_size=min_size,
        max_size=max_size,
        require_priming=True,
    )
    return amplicons, {forward_name: fwd, reverse_name: rev}


def summarise_panel_amplicons(amplicons: Sequence[PredictedAmplicon]) -> list[dict[str, Any]]:
    """Condense predicted products into report-ready dictionaries."""
    out: list[dict[str, Any]] = []
    for amp in amplicons:
        out.append(
            {
                "subject": amp.subject,
                "size": amp.size,
                "forward_primer": amp.forward_primer,
                "reverse_primer": amp.reverse_primer,
                "forward_mismatches": amp.forward_site.mismatches,
                "reverse_mismatches": amp.reverse_site.mismatches,
                "forward_position_1based": f"{amp.forward_site.start + 1}-{amp.forward_site.end}",
                "reverse_position_1based": f"{amp.reverse_site.start + 1}-{amp.reverse_site.end}",
                "note": amp.note,
            }
        )
    return out
