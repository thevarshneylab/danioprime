"""Riboprobe templates for in situ hybridisation.

Two things make this different from every other mode in the package, and both
are ways a probe can fail that a PCR-shaped analysis will not catch.

**Promoter orientation decides which strand you get.** A promoter on the
*reverse* primer transcribes the antisense strand -- the probe that detects
mRNA. On the *forward* primer it transcribes the sense strand, which is the
negative control. Getting this backwards yields a probe that cannot hybridise
and a control that can, and the experiment looks merely "failed". Putting
different promoters at the two ends makes one PCR product serve as template for
both.

**Cross-hybridisation is not primer specificity.** Everywhere else in this
package, specificity is decided at the primer's 3' end, because that is what a
polymerase needs. A riboprobe has no 3' end that matters: it is hundreds of
bases of RNA, and it will hybridise to anything sufficiently similar along its
whole length. So a pair whose *primers* discriminate the ohnolog perfectly can
still produce a probe that stains both genes -- the primers sit in divergent
sequence while the body of the amplicon is conserved. This module therefore
analyses the amplicon itself, not just its ends.

The published expression pattern depends on getting that right, which is why
the check is on by default whenever an avoid-panel is supplied.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any

from ..config import DesignParams
from ..design.pipeline import DesignRequest, DesignResult, run_design
from ..design.regions import IntervalSet
from ..errors import DesignError
from ..sequence import SeqRecord
from ..specificity.index import KmerIndex
from ..specificity.paralog import ParalogPanel
from .tails import RIBOPROBE_SCHEMES, TAIL_SCHEMES, ReactionPlan, reaction_plan

#: Seed length for finding conserved blocks between a probe and the panel.
BLOCK_SEED = 12

#: A contiguous block at or above this identity, and at least
#: :data:`CROSS_HYB_MIN_LENGTH` long, is treated as able to cross-hybridise
#: under standard in situ stringency.
#:
#: These defaults are calibrated rather than guessed. The failure mode to avoid
#: is a chance window: between two genuinely homologous but diverged ohnologs
#: (say 65% identity) a locally better-conserved patch turns up by luck, and a
#: lenient threshold reports it as a risk. Measured over 100 simulated
#: 700 nt backgrounds, the proportion of runs tripped purely by chance was:
#:
#: ===============  ==========  ==============  ================
#: identity         min length  65% background  25% background
#: ===============  ==========  ==============  ================
#: 0.80             50          57%             0%
#: 0.85             50          15%             0%
#: 0.80             100         2%              0%
#: **0.85**         **100**     **0%**          **0%**
#: ===============  ==========  ==============  ================
#:
#: So 85% over 100 nt is the first combination that is quiet against a diverged
#: paralog while still firing on real conservation. It is also the biologically
#: defensible one: 100 nt at 85% identity will hold a riboprobe through a
#: standard high-stringency wash, whereas 50 nt at 80% will not reliably.
#: Both are tunable, and a lower threshold is reasonable if you would rather
#: investigate false alarms than miss a real one.
CROSS_HYB_IDENTITY = 0.85

#: Blocks shorter than this are unlikely to hold a probe under wash stringency.
CROSS_HYB_MIN_LENGTH = 100


@dataclass(slots=True)
class ConservedBlock:
    """An ungapped stretch shared between a probe and a panel sequence."""

    subject: str
    probe_start: int
    probe_end: int
    subject_start: int
    identity: float

    @property
    def length(self) -> int:
        return self.probe_end - self.probe_start

    def to_dict(self) -> dict[str, Any]:
        return {
            "subject": self.subject,
            "probe_region_1based": f"{self.probe_start + 1}-{self.probe_end}",
            "subject_start_1based": self.subject_start + 1,
            "length": self.length,
            "identity_percent": round(100 * self.identity, 1),
        }


@dataclass(slots=True)
class CrossHybridisation:
    """Whether a probe is predicted to stain more than its own transcript."""

    blocks: list[ConservedBlock] = field(default_factory=list)
    min_identity: float = CROSS_HYB_IDENTITY
    min_length: int = CROSS_HYB_MIN_LENGTH
    probe_length: int = 0

    @property
    def risky_blocks(self) -> list[ConservedBlock]:
        return [
            b
            for b in self.blocks
            if b.length >= self.min_length and b.identity >= self.min_identity
        ]

    @property
    def will_cross_hybridise(self) -> bool:
        return bool(self.risky_blocks)

    @property
    def worst(self) -> ConservedBlock | None:
        if not self.blocks:
            return None
        return max(self.blocks, key=lambda b: (b.length * b.identity, b.length))

    def verdict(self) -> str:
        if not self.blocks:
            return "no similarity to the avoid-panel detected along the probe"
        worst = self.worst
        assert worst is not None
        if self.will_cross_hybridise:
            pct = 100 * worst.length / max(1, self.probe_length)
            return (
                f"predicted to cross-hybridise: {worst.length} nt at "
                f"{100 * worst.identity:.0f}% identity to {worst.subject} "
                f"({pct:.0f}% of the probe). Primer-level discrimination does not "
                f"prevent this -- the probe body is what hybridises."
            )
        return (
            f"longest shared block is {worst.length} nt at "
            f"{100 * worst.identity:.0f}% identity to {worst.subject}, below the "
            f"{self.min_length} nt / {100 * self.min_identity:.0f}% threshold"
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "will_cross_hybridise": self.will_cross_hybridise,
            "verdict": self.verdict(),
            "probe_length": self.probe_length,
            "threshold": {
                "min_length": self.min_length,
                "min_identity_percent": round(100 * self.min_identity, 1),
            },
            "blocks": [b.to_dict() for b in sorted(self.blocks, key=lambda x: -x.length)[:10]],
        }


def _identity(a: str, b: str) -> float:
    if not a:
        return 0.0
    return sum(1 for x, y in zip(a, b, strict=True) if x == y) / len(a)


def _longest_window_above(matches: list[bool], threshold: float) -> tuple[int, int] | None:
    """Longest window of *matches* whose mean is at least *threshold*.

    Linear time. Subtracting the threshold turns "mean at least t" into "sum at
    least zero", and the longest such window is found by scanning a stack of
    decreasing prefix minima from the right.

    This replaces a seed-extension loop that had a real defect: it took the
    span between the first and last seed on a diagonal as the starting block
    without checking it, so two chance 12-mer matches 600 nt apart were
    reported as one 600 nt block at 64% identity.
    """
    n = len(matches)
    if n == 0:
        return None
    # Integer arithmetic throughout. In floating point the prefix sums of
    # (match - threshold) accumulate error, and two sums that are nominally
    # equal can differ by an ulp -- which silently loses the longest window,
    # because the comparison that should accept it fails by a hair.
    scale = 10_000
    step = round(scale * threshold)
    prefix = [0] * (n + 1)
    for i, m in enumerate(matches):
        prefix[i + 1] = prefix[i] + ((scale if m else 0) - step)

    # Indices of strictly decreasing prefix minima: any window must start here.
    stack: list[int] = []
    for i in range(n + 1):
        if not stack or prefix[i] < prefix[stack[-1]]:
            stack.append(i)

    best: tuple[int, int] | None = None
    best_len = 0
    for j in range(n, -1, -1):
        while stack and prefix[j] >= prefix[stack[-1]]:
            i = stack.pop()
            if j - i > best_len:
                best_len = j - i
                best = (i, j)
    return best if best_len > 0 else None


def find_conserved_blocks(
    probe: str,
    panel: ParalogPanel,
    *,
    seed: int = BLOCK_SEED,
    min_identity: float = CROSS_HYB_IDENTITY,
) -> list[ConservedBlock]:
    """Ungapped stretches the probe shares with the panel, on either strand.

    Exact k-mers seed candidate alignment diagonals; on each diagonal the
    longest window whose identity reaches *min_identity* is then reported. Only
    ungapped similarity is modelled, so an indel splits one block into two.
    That understates the risk slightly, which is the safe direction for a check
    whose job is to warn.
    """
    from ..sequence import revcomp  # noqa: PLC0415

    if not panel or len(probe) < seed:
        return []

    index = KmerIndex.build(panel.records, k=seed)
    blocks: list[ConservedBlock] = []

    for is_forward, strand_seq in ((True, probe), (False, revcomp(probe))):
        diagonals: set[tuple[int, int]] = set()
        for i in range(len(strand_seq) - seed + 1):
            for sid, pos in index.lookup(strand_seq[i : i + seed]):
                diagonals.add((sid, pos - i))

        for sid, offset in diagonals:
            subject = index.subject(sid).seq
            # Region of the probe that overlaps the subject on this diagonal.
            lo = max(0, -offset)
            hi = min(len(strand_seq), len(subject) - offset)
            if hi - lo < seed:
                continue
            matches = [strand_seq[i] == subject[i + offset] for i in range(lo, hi)]
            window = _longest_window_above(matches, min_identity)
            if window is None:
                continue
            w_lo, w_hi = lo + window[0], lo + window[1]
            if w_hi - w_lo < seed:
                continue
            a = strand_seq[w_lo:w_hi]
            b = subject[w_lo + offset : w_hi + offset]
            identity = _identity(a, b)
            p_start, p_end = (
                (w_lo, w_hi) if is_forward else (len(probe) - w_hi, len(probe) - w_lo)
            )
            blocks.append(
                ConservedBlock(
                    subject=index.subject(sid).name,
                    probe_start=p_start,
                    probe_end=p_end,
                    subject_start=w_lo + offset,
                    identity=identity,
                )
            )

    # Collapse blocks contained within a longer one on the same subject.
    blocks.sort(key=lambda b: (-b.length, b.probe_start))
    kept: list[ConservedBlock] = []
    for block in blocks:
        if any(
            k.subject == block.subject
            and k.probe_start <= block.probe_start
            and block.probe_end <= k.probe_end
            for k in kept
        ):
            continue
        kept.append(block)
    return kept


def assess_probe(
    probe: str,
    panel: ParalogPanel | None,
    *,
    seed: int = BLOCK_SEED,
    min_identity: float = CROSS_HYB_IDENTITY,
    min_length: int = CROSS_HYB_MIN_LENGTH,
) -> CrossHybridisation:
    """Decide whether a probe sequence will also stain the avoid-panel."""
    result = CrossHybridisation(
        min_identity=min_identity, min_length=min_length, probe_length=len(probe)
    )
    if panel:
        result.blocks = find_conserved_blocks(
            probe, panel, seed=seed, min_identity=min_identity
        )
    return result


@dataclass(slots=True)
class InSituDesign:
    """A probe-template design with its probe-level specificity assessment."""

    result: DesignResult
    scheme: str
    antisense: ReactionPlan | None = None
    sense: ReactionPlan | None = None
    cross_hyb: dict[str, CrossHybridisation] = field(default_factory=dict)

    @property
    def pairs(self):  # noqa: ANN201
        return self.result.pairs

    def to_dict(self) -> dict[str, Any]:
        return {
            "assay": "insitu",
            "scheme": self.scheme,
            "antisense": self.antisense.to_dict() if self.antisense else None,
            "sense": self.sense.to_dict() if self.sense else None,
            "cross_hybridisation": {k: v.to_dict() for k, v in self.cross_hyb.items()},
        }


#: Which scheme transcribes the sense control for a given antisense scheme.
_SENSE_PARTNER = {
    "riboprobe-t7": "riboprobe-t7-sense",
    "riboprobe-t3": "riboprobe-t3-sense",
    "riboprobe-sp6": "riboprobe-sp6-sense",
}


def design_insitu(
    template: SeqRecord,
    params: DesignParams,
    *,
    scheme: str = "riboprobe-t7",
    panel: ParalogPanel | None = None,
    variants: IntervalSet | None = None,
    targets: IntervalSet | None = None,
    sense_control: bool = True,
    min_identity: float = CROSS_HYB_IDENTITY,
    min_block: int = CROSS_HYB_MIN_LENGTH,
    name: str | None = None,
) -> InSituDesign:
    """Design a template for riboprobe synthesis.

    *scheme* names the promoter arrangement (see
    :data:`~danioprime.assays.tails.TAIL_SCHEMES`). A dual-promoter scheme
    produces one template for both probes; otherwise, with *sense_control*, a
    second ordering plan is emitted for the sense control.

    Every reported pair's amplicon is checked against *panel* for stretches
    long and similar enough to cross-hybridise, which is a different and
    stricter question than whether the primers discriminate.
    """
    if scheme not in TAIL_SCHEMES:
        raise DesignError(
            f"unknown tail scheme {scheme!r}; riboprobe schemes are "
            f"{', '.join(RIBOPROBE_SCHEMES)}"
        )
    if scheme not in RIBOPROBE_SCHEMES and scheme != "none":
        raise DesignError(
            f"{scheme!r} attaches no RNA-polymerase promoter, so the product cannot be "
            f"transcribed. Choose one of: {', '.join(RIBOPROBE_SCHEMES)}"
        )

    result = run_design(
        DesignRequest(
            template=template,
            params=params,
            targets=targets or IntervalSet(),
            variants=variants or IntervalSet(),
            panel=panel,
            name=name or f"{template.name}_probe",
            context={"assay": "insitu", "scheme": scheme},
        )
    )

    design = InSituDesign(result=result, scheme=scheme)

    if result.pairs:
        mapping = TAIL_SCHEMES[scheme]
        is_dual = bool(mapping.get("forward")) and bool(mapping.get("reverse"))
        if sense_control and not is_dual and scheme in _SENSE_PARTNER:
            # The sense control needs its own PCR with the promoter at the
            # other end. Build it on a copy, because reaction_plan attaches the
            # tails to the primers it is given.
            design.sense = reaction_plan(
                copy.deepcopy(result.pairs), scheme=_SENSE_PARTNER[scheme]
            )
            result.warnings.append(
                f"The sense control needs a second PCR with the "
                f"'{_SENSE_PARTNER[scheme]}' primers listed separately. A dual-promoter "
                f"scheme such as riboprobe-dual-t7-sp6 would give both probes from one "
                f"product instead."
            )
        # Applied last so the reported pairs carry the antisense tails.
        design.antisense = reaction_plan(result.pairs, scheme=scheme)

    if panel:
        for pair in result.pairs:
            probe = pair.amplicon_seq(template.seq)
            assessment = assess_probe(
                probe, panel, min_identity=min_identity, min_length=min_block
            )
            design.cross_hyb[pair.name] = assessment
            if assessment.will_cross_hybridise:
                pair.notes.append(f"probe: {assessment.verdict()}")
            else:
                pair.notes.append(f"probe: {assessment.verdict()}")
        clean = [p for p in result.pairs if not design.cross_hyb[p.name].will_cross_hybridise]
        if clean and len(clean) < len(result.pairs):
            result.warnings.append(
                f"{len(result.pairs) - len(clean)} of {len(result.pairs)} probe templates "
                f"are predicted to cross-hybridise with the panel; the clean ones are "
                f"listed first."
            )
            order = {p.name: 0 if p in clean else 1 for p in result.pairs}
            result.pairs.sort(key=lambda p: (order[p.name], p.penalty))
        elif not clean and result.pairs:
            result.warnings.append(
                "every probe template shares a long, highly similar block with the "
                "avoid-panel. Target a divergent region (often the 3' UTR, which "
                "diverges faster than coding sequence between ohnologs) with --target, "
                "or accept that the probe reports both genes."
            )

    if result.pairs:
        shortest = min(p.product_size for p in result.pairs)
        if shortest < 300:
            result.warnings.append(
                f"the shortest template is {shortest} bp; riboprobes below about 300 nt "
                f"give weak signal. Raise --min-product-size unless you specifically "
                f"need a short probe."
            )
    return design
