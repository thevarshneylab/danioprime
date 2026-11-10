"""K-mer index for seeded binding-site search.

Rationale for seeding on the primer's **3'-terminal** k-mer rather than any
k-mer: a site that can be extended by the polymerase necessarily has a
perfectly matched 3' end. Anchoring the seed there therefore loses nothing that
could prime efficiently, while cutting the search space by 4**k.

Sensitivity contract
--------------------
With seed length *k*, seeded search finds **every** site whose 3'-terminal *k*
bases match the primer exactly. It misses sites carrying a mismatch inside that
terminal window but outside the "clean 3' end" window used by
:class:`~danioprime.specificity.scan.PrimingRules`. Those sites are, by
construction, poor primers -- but if you need a guarantee, use exhaustive
scanning (exact, ~300 kb/s/primer) or an external aligner for genome scale.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field

from ..sequence import DNA_BASES, SeqRecord

DEFAULT_SEED = 12


@dataclass(slots=True)
class KmerIndex:
    """Positions of every k-mer across a set of subject sequences."""

    k: int = DEFAULT_SEED
    #: kmer -> list of (subject_id, position)
    table: dict[str, list[tuple[int, int]]] = field(default_factory=lambda: defaultdict(list))
    subjects: list[SeqRecord] = field(default_factory=list)

    @property
    def total_length(self) -> int:
        return sum(len(s.seq) for s in self.subjects)

    def add(self, record: SeqRecord) -> None:
        """Index one record. Windows containing ambiguity codes are skipped."""
        sid = len(self.subjects)
        self.subjects.append(record)
        seq = record.seq
        k = self.k
        n = len(seq)
        if n < k:
            return
        # Track the position of the most recent non-ACGT base so windows
        # containing one can be skipped without rescanning.
        last_bad = -1
        for i, base in enumerate(seq):
            if base not in DNA_BASES:
                last_bad = i
            start = i - k + 1
            if start < 0:
                continue
            if last_bad >= start:
                continue
            self.table[seq[start : i + 1]].append((sid, start))

    @classmethod
    def build(cls, records: Iterable[SeqRecord], k: int = DEFAULT_SEED) -> KmerIndex:
        idx = cls(k=k)
        for record in records:
            idx.add(record)
        return idx

    def lookup(self, kmer: str) -> list[tuple[int, int]]:
        """Positions of *kmer*; empty if absent or wrong length."""
        if len(kmer) != self.k:
            return []
        return self.table.get(kmer.upper(), [])

    def subject(self, sid: int) -> SeqRecord:
        return self.subjects[sid]

    def iter_subjects(self) -> Iterator[tuple[int, SeqRecord]]:
        return iter(enumerate(self.subjects))

    def stats(self) -> dict[str, int | float]:
        """Index size and seed specificity, useful for choosing *k*."""
        if not self.table:
            return {"subjects": len(self.subjects), "distinct_kmers": 0, "total_length": 0}
        counts = [len(v) for v in self.table.values()]
        return {
            "subjects": len(self.subjects),
            "total_length": self.total_length,
            "seed_length": self.k,
            "distinct_kmers": len(counts),
            "max_kmer_occurrences": max(counts),
            "mean_kmer_occurrences": round(sum(counts) / len(counts), 3),
        }
