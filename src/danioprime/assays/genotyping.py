"""Genotyping an edited locus.

**There is deliberately no cut site here.** An earlier version anchored the
whole design on one predicted cut coordinate: it centred the amplicon there,
carved out a protected window of +/- n bases, and searched for primers in two
symmetric flanking windows. That was the wrong criterion, for a reason that
only shows up in the fish rather than in the output.

A predicted cut coordinate is not where the lesion is. The blunt cut is, but
what you genotype is the repair: NHEJ routinely gives deletions of tens to
hundreds of bases, insertions, and complex rearrangements, and an F0 mosaic
carries a spectrum of them at once. Sizing a protected window from "the largest
indel you expect" is a guess that holds right up until the allele that matters
is bigger than the guess -- and then a primer sits inside the deletion, the
mutant allele does not amplify, and the fish reads as wild type. A design
anchored on the cut optimises confidently for a coordinate whose relevance it
cannot check, and fails silently in the direction that costs a generation.

So genotyping here is the ordinary design problem plus a readout check. Primer
placement is driven by primer quality, specificity against the ohnolog, and
strain variants -- and, where the user says so, by ``targets``: the amplicon
must span this region. That is a statement about what is to be interrogated,
which the user knows, rather than about where a nuclease is predicted to cut,
which the tool was only ever inferring.

What remains from the old geometry is the part that was always sound:
:func:`predict_allele_sizes` asks whether the readout can resolve the alleles
you care about. That needs the product size and nothing else, and it is where
the real constraint lives -- a 400 bp product cannot show you a 2 bp indel on
agarose however carefully it was placed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ..config import DesignParams
from ..design.pipeline import DesignRequest, DesignResult, run_design
from ..design.regions import Interval, IntervalSet
from ..errors import DesignError
from ..sequence import SeqRecord
from ..specificity.paralog import ParalogPanel

#: Indel sizes product sizes are predicted for when the caller names none.
#: Spread rather than exhaustive: the point is to find the size below which the
#: chosen readout stops separating alleles, not to enumerate outcomes.
DEFAULT_EXPECTED_INDELS: tuple[int, ...] = (-10, -5, -2, -1, 1, 2, 5, 10)


@dataclass(slots=True)
class GenotypingDesign:
    """Result of a genotyping design, with predicted allele sizes."""

    result: DesignResult
    expected_indels: tuple[int, ...] = ()
    allele_sizes: list[dict[str, Any]] = field(default_factory=list)

    @property
    def pairs(self):  # noqa: ANN201 - mirrors DesignResult
        return self.result.pairs

    def to_dict(self) -> dict[str, Any]:
        return {"allele_sizes": self.allele_sizes}


def predict_allele_sizes(
    wt_size: int, indels: tuple[int, ...] = DEFAULT_EXPECTED_INDELS
) -> list[dict[str, Any]]:
    """Product sizes for a wild-type and a set of hypothetical indel alleles.

    Also flags which of them a standard 2% agarose gel will not separate from
    wild type, because that determines whether the assay needs PAGE, fragment
    analysis or a restriction digest instead.
    """
    rows: list[dict[str, Any]] = [
        {"allele": "wild type", "indel": 0, "size": wt_size, "resolvable_on_gel": True}
    ]
    for indel in sorted(indels):
        if indel == 0:
            continue
        size = wt_size + indel
        if size <= 0:
            continue
        # 2% agarose resolves roughly 10 bp around a few hundred bp; below that
        # a size shift is not reliably callable by eye.
        resolvable = abs(indel) >= 10 or (abs(indel) / wt_size) > 0.03
        rows.append(
            {
                "allele": f"{indel:+d} bp",
                "indel": indel,
                "size": size,
                "resolvable_on_gel": resolvable,
            }
        )
    return rows


def design_genotyping(
    template: SeqRecord,
    params: DesignParams,
    *,
    targets: IntervalSet | None = None,
    panel: ParalogPanel | None = None,
    variants: IntervalSet | None = None,
    excluded: IntervalSet | None = None,
    expected_indels: tuple[int, ...] = DEFAULT_EXPECTED_INDELS,
    name: str | None = None,
) -> GenotypingDesign:
    """Design a pair for genotyping, and say what its readout can resolve.

    *targets*, when given, are regions the amplicon must span and no primer may
    overlap -- the exon you disrupted, say. Give none and this is a plain design
    over the whole template, which is the right answer when you do not know
    where the lesion is: a wide amplicon that certainly contains it beats a
    narrow one placed around a coordinate that might be wrong.
    """
    n = len(template.seq)
    targets = targets or IntervalSet()
    for interval in targets:
        if not (0 <= interval.start < interval.end <= n):
            raise DesignError(
                f"target {interval.start + 1}-{interval.end} lies outside the template (1-{n})"
            )

    request = DesignRequest(
        template=template,
        params=params,
        targets=targets,
        # A region worth interrogating is a region no primer may sit in: a
        # primer inside the lesion does not amplify the mutant allele.
        excluded=(excluded or IntervalSet()).union(targets),
        variants=variants or IntervalSet(),
        panel=panel,
        name=name or f"{template.name}_gt",
        context={"assay": "genotyping"},
    )
    result = run_design(request)

    design = GenotypingDesign(result=result, expected_indels=expected_indels)
    if result.best:
        design.allele_sizes = predict_allele_sizes(result.best.product_size, expected_indels)
        unresolvable = [r["allele"] for r in design.allele_sizes if not r["resolvable_on_gel"]]
        if unresolvable:
            result.warnings.append(
                f"a {result.best.product_size} bp product will not resolve "
                f"{', '.join(unresolvable)} alleles on a standard gel; consider the "
                "'genotyping-fluorescent' preset, a dCAPS assay, or PAGE"
            )
    return design


def design_junction_pcr(
    genomic: SeqRecord,
    insertion_point: int,
    cassette: SeqRecord,
    params: DesignParams,
    *,
    genomic_flank: int = 500,
    cassette_depth: int = 400,
    homology_arm: int = 0,
    panel: ParalogPanel | None = None,
    variants: IntervalSet | None = None,
    name: str | None = None,
) -> dict[str, DesignResult]:
    """Design the two junction PCRs that verify a targeted knock-in.

    One primer must sit in genomic sequence outside the homology arm and the
    other inside the cassette; a pair with both primers in the cassette
    amplifies the donor plasmid just as happily as the integrated allele, which
    is the classic false positive in knock-in screening.

    Returns the ``5'`` and ``3'`` junction designs separately, since each needs
    its own pair.
    """
    n = len(genomic.seq)
    if not 0 <= insertion_point <= n:
        raise DesignError(f"insertion point {insertion_point + 1} lies outside the genomic template")

    out: dict[str, DesignResult] = {}

    # 5' junction: forward primer in genomic sequence upstream of the arm,
    # reverse primer reading back out of the cassette's 5' end.
    five_template = SeqRecord(
        name=f"{genomic.name}_5p_junction",
        seq=genomic.seq[max(0, insertion_point - genomic_flank) : insertion_point]
        + cassette.seq[:cassette_depth],
    )
    boundary_5 = min(insertion_point, genomic_flank)
    out["5prime"] = run_design(
        DesignRequest(
            template=five_template,
            params=params,
            forward_search=IntervalSet([Interval(0, max(1, boundary_5 - homology_arm))]),
            reverse_search=IntervalSet([Interval(boundary_5, len(five_template.seq))]),
            targets=IntervalSet([Interval(max(0, boundary_5 - 1), min(len(five_template.seq), boundary_5 + 1))]),
            panel=panel,
            variants=variants or IntervalSet(),
            name=name and f"{name}_5p" or f"{genomic.name}_5p",
            context={"assay": "junction-pcr", "junction": "5prime", "boundary": boundary_5},
        )
    )

    # 3' junction: forward primer inside the cassette, reverse primer in
    # genomic sequence downstream of the arm.
    three_template = SeqRecord(
        name=f"{genomic.name}_3p_junction",
        seq=cassette.seq[-cassette_depth:]
        + genomic.seq[insertion_point : insertion_point + genomic_flank],
    )
    boundary_3 = min(cassette_depth, len(cassette.seq))
    out["3prime"] = run_design(
        DesignRequest(
            template=three_template,
            params=params,
            forward_search=IntervalSet([Interval(0, boundary_3)]),
            reverse_search=IntervalSet(
                [Interval(min(len(three_template.seq), boundary_3 + homology_arm), len(three_template.seq))]
            ),
            targets=IntervalSet([Interval(max(0, boundary_3 - 1), min(len(three_template.seq), boundary_3 + 1))]),
            panel=panel,
            variants=variants or IntervalSet(),
            name=name and f"{name}_3p" or f"{genomic.name}_3p",
            context={"assay": "junction-pcr", "junction": "3prime", "boundary": boundary_3},
        )
    )
    return out
