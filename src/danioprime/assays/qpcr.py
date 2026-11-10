"""RT-qPCR primer design across exon-exon junctions.

The specific hazard in RT-qPCR is amplifying contaminating genomic DNA, which
inflates every measurement and is invisible unless a no-RT control is run.
Two independent defences, in order of strength:

1. **A junction-spanning primer.** Its 3' half anneals to one exon and its 5'
   half to the next, so on genomic DNA it has no contiguous binding site at
   all. This works regardless of intron size and is the only defence that
   holds for an intronless gene's retained-intron neighbours.
2. **An intron-spanning amplicon.** Both primers sit in different exons, so a
   genomic template yields a longer product. Depends on the intron being large
   enough that the genomic product fails to amplify under qPCR cycling -- a
   70 bp intron gives no protection.

Zebrafish adds a third consideration: with an ohnolog present, the amplicon
must also discriminate the paralog, or the measurement reports the sum of two
genes. That check comes from the shared specificity engine, with the ohnolog's
*transcript* sequence as the avoid-panel.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ..annotation.gtf import Transcript
from ..config import DesignParams
from ..design.pipeline import DesignRequest, DesignResult, run_design
from ..design.regions import Interval, IntervalSet
from ..errors import DesignError
from ..sequence import SeqRecord
from ..specificity.paralog import ParalogPanel


@dataclass(slots=True)
class QpcrDesign:
    """A qPCR design plus its genomic-DNA protection assessment."""

    result: DesignResult
    junctions: tuple[int, ...]
    require_junction_spanning: bool
    intron_sizes: dict[int, int] = field(default_factory=dict)
    assessment: list[dict[str, Any]] = field(default_factory=list)

    @property
    def pairs(self):  # noqa: ANN201
        return self.result.pairs

    def to_dict(self) -> dict[str, Any]:
        return {
            "junction_offsets": list(self.junctions),
            "require_junction_spanning": self.require_junction_spanning,
            "intron_sizes": {str(k): v for k, v in self.intron_sizes.items()},
            "gdna_protection": self.assessment,
        }


def intron_sizes_from_transcript(tx: Transcript) -> dict[int, int]:
    """Intron length keyed by the transcript-coordinate junction it creates."""
    exons = tx.sorted_exons()
    sizes: dict[int, int] = {}
    cursor = 0
    for upstream, downstream in zip(exons, exons[1:], strict=False):
        cursor += len(upstream)
        if tx.strand == "+":
            gap = downstream.start - upstream.end
        else:
            gap = upstream.start - downstream.end
        sizes[cursor] = max(0, gap)
    return sizes


def assess_gdna_protection(
    design: QpcrDesign, *, min_protective_intron: int = 1000
) -> list[dict[str, Any]]:
    """For each reported pair, say how it is protected against genomic DNA."""
    rows: list[dict[str, Any]] = []
    for pair in design.result.pairs:
        spanning = [
            p.orientation
            for p in (pair.forward, pair.reverse)
            if p.spans_junction is not None
        ]
        crossed = [
            j for j in design.junctions if pair.forward.start < j < pair.reverse.end
        ]
        protective = [
            j for j in crossed if design.intron_sizes.get(j, 0) >= min_protective_intron
        ]
        if spanning:
            mechanism = f"junction-spanning primer ({', '.join(spanning)})"
            strength = "strong"
        elif protective:
            biggest = max(design.intron_sizes.get(j, 0) for j in protective)
            mechanism = f"amplicon spans an intron of {biggest} bp"
            strength = "moderate"
        elif crossed:
            biggest = max((design.intron_sizes.get(j, 0) for j in crossed), default=0)
            mechanism = (
                f"amplicon spans a junction but the intron is only {biggest} bp"
                if biggest
                else "amplicon spans a junction of unknown intron size"
            )
            strength = "weak"
        else:
            mechanism = "no junction crossed: genomic DNA will amplify identically"
            strength = "none"
        rows.append(
            {
                "pair": pair.name,
                "protection": strength,
                "mechanism": mechanism,
                "junctions_crossed": crossed,
                "junction_spanning_primers": spanning,
            }
        )
        if strength in {"none", "weak"}:
            pair.notes.append(f"gDNA protection {strength}: {mechanism}")
    design.assessment = rows
    return rows


def design_qpcr(
    transcript: SeqRecord,
    params: DesignParams,
    *,
    junctions: tuple[int, ...] = (),
    transcript_model: Transcript | None = None,
    require_junction_spanning: bool = False,
    junction_flank: int = 60,
    panel: ParalogPanel | None = None,
    variants: IntervalSet | None = None,
    min_protective_intron: int = 1000,
    name: str | None = None,
) -> QpcrDesign:
    """Design a qPCR pair on a spliced transcript.

    *junctions* are offsets in the transcript at which exon-exon boundaries
    fall; pass *transcript_model* instead to take them (and the intron sizes)
    from a parsed GTF record.

    With *require_junction_spanning*, the primer search is confined to windows
    around the junctions, which both enforces the requirement and makes the
    search much cheaper than scoring the whole transcript and discarding
    almost all of it.
    """
    if transcript_model is not None:
        junctions = transcript_model.junction_offsets()
        intron_sizes = intron_sizes_from_transcript(transcript_model)
        if transcript_model.spliced_length != len(transcript.seq):
            raise DesignError(
                f"transcript model says {transcript_model.spliced_length} nt but the "
                f"supplied sequence is {len(transcript.seq)} nt; they must correspond"
            )
    else:
        intron_sizes = {}

    n = len(transcript.seq)
    bad = [j for j in junctions if not 0 < j < n]
    if bad:
        raise DesignError(f"junction offsets outside the transcript (1-{n - 1}): {bad}")

    forward_search = IntervalSet()
    reverse_search = IntervalSet()
    if require_junction_spanning:
        if not junctions:
            raise DesignError(
                "--require-junction-spanning needs junction offsets; supply --junctions "
                "or a GTF with --gtf/--transcript"
            )
        # Either primer may be the junction-spanning one, so both strands are
        # searched over the same space: the junction neighbourhoods plus enough
        # room on each side for a partner within the product-size window. Any
        # valid junction-spanning pair must lie inside that, so this is a pure
        # performance restriction, not an extra constraint.
        reach = junction_flank + params.max_product_size
        windows = IntervalSet(
            [Interval(max(0, j - reach), min(n, j + reach)) for j in junctions]
        )
        forward_search = windows
        reverse_search = windows

    request = DesignRequest(
        template=transcript,
        params=params,
        forward_search=forward_search,
        reverse_search=reverse_search,
        variants=variants or IntervalSet(),
        junctions=tuple(junctions),
        panel=panel,
        require_spanning="either" if require_junction_spanning else None,
        name=name or f"{transcript.name}_qpcr",
        context={
            "assay": "qpcr",
            "junctions": list(junctions),
            "require_junction_spanning": require_junction_spanning,
        },
    )
    result = run_design(request)

    if require_junction_spanning:
        spanning = [
            pair
            for pair in result.pairs
            if pair.forward.spans_junction is not None or pair.reverse.spans_junction is not None
        ]
        if spanning:
            result.pairs = spanning
        elif result.pairs:
            blocked = {
                **result.forward_candidates.rejections,
                **result.reverse_candidates.rejections,
            }
            worst = max(blocked, key=blocked.get) if blocked else "the active filters"
            result.warnings.append(
                "no primer could be placed across an exon-exon junction under the current "
                f"constraints (most candidates there were removed by: {worst}). Falling back "
                "to intron-spanning pairs, ranked with the better-protected ones first -- "
                "check the gDNA protection column before ordering, and consider widening "
                "--min-tm/--max-tm or --junction-flank."
            )

    design = QpcrDesign(
        result=result,
        junctions=tuple(junctions),
        require_junction_spanning=require_junction_spanning,
        intron_sizes=intron_sizes,
    )
    assess_gdna_protection(design, min_protective_intron=min_protective_intron)

    if require_junction_spanning and not any(
        p.forward.spans_junction is not None or p.reverse.spans_junction is not None
        for p in result.pairs
    ):
        rank = {"strong": 0, "moderate": 1, "weak": 2, "none": 3}
        strength = {row["pair"]: rank.get(row["protection"], 3) for row in design.assessment}
        result.pairs.sort(key=lambda pr: (strength.get(pr.name, 3), pr.penalty))
        design.assessment.sort(key=lambda row: rank.get(row["protection"], 3))
    return design
