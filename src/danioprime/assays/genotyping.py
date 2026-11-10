"""Genotyping an edited locus.

Designing a genotyping assay is a geometry problem before it is a
thermodynamics problem, and it is the geometry that generic tools leave to the
user. Three constraints have to hold simultaneously:

* the amplicon must span the lesion, or the assay reports nothing;
* neither primer may overlap the lesion, or the mutant allele will not amplify
  -- which looks like a homozygote and is the failure mode that costs a
  generation of fish;
* the primers must sit far enough out that the readout works: enough margin for
  a Sanger read to clear the primer, or a short enough product for a few-base
  indel to resolve on a gel.

The protected window is sized from the largest indel the user expects rather
than the one they got, because the assay has to work before the alleles are
known.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from ..config import DesignParams
from ..design.pipeline import DesignRequest, DesignResult, run_design
from ..design.regions import Interval, IntervalSet
from ..errors import DesignError
from ..sequence import SeqRecord, revcomp
from ..specificity.paralog import ParalogPanel

Nuclease = Literal["SpCas9", "SaCas9", "Cas12a", "manual"]

#: PAM and cut geometry per nuclease. ``cut_offset`` is measured from the 5' end
#: of the protospacer on the strand the protospacer is written on.
NUCLEASES: dict[str, dict[str, Any]] = {
    "SpCas9": {
        "pam": "NGG",
        "pam_side": 3,
        "spacer_length": 20,
        "cut_offset": 17,
        "blunt": True,
        "description": "SpCas9: 20 nt spacer, 3' NGG PAM, blunt cut 3 bp from the PAM",
    },
    "SaCas9": {
        "pam": "NNGRRT",
        "pam_side": 3,
        "spacer_length": 21,
        "cut_offset": 18,
        "blunt": True,
        "description": "SaCas9: 21 nt spacer, 3' NNGRRT PAM",
    },
    "Cas12a": {
        "pam": "TTTV",
        "pam_side": 5,
        "spacer_length": 23,
        "cut_offset": 18,
        "blunt": False,
        "description": "Cas12a/Cpf1: 5' TTTV PAM, staggered cut leaving 5' overhangs",
    },
}


def _iupac_match(pattern: str, seq: str) -> bool:
    from ..sequence import IUPAC

    if len(pattern) != len(seq):
        return False
    return all(base in IUPAC.get(code, frozenset()) for code, base in zip(pattern, seq, strict=True))


@dataclass(slots=True)
class CutSite:
    """A nuclease target located on a template."""

    position: int
    strand: Literal["+", "-"]
    protospacer: str = ""
    pam: str = ""
    nuclease: str = "manual"
    spacer_start: int = -1
    spacer_end: int = -1

    @property
    def description(self) -> str:
        if not self.protospacer:
            return f"manual cut site at {self.position + 1}"
        return (
            f"{self.nuclease} cut at {self.position + 1} "
            f"({self.strand} strand, spacer {self.spacer_start + 1}-{self.spacer_end}, "
            f"PAM {self.pam})"
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "position_1based": self.position + 1,
            "strand": self.strand,
            "nuclease": self.nuclease,
            "protospacer": self.protospacer,
            "pam": self.pam,
            "spacer_1based": (
                f"{self.spacer_start + 1}-{self.spacer_end}" if self.spacer_start >= 0 else None
            ),
            "description": self.description,
        }


def locate_protospacer(
    template: str,
    protospacer: str,
    *,
    nuclease: Nuclease = "SpCas9",
    require_pam: bool = True,
) -> list[CutSite]:
    """Find a protospacer on either strand and compute the cut position.

    The protospacer is given without its PAM, as ordered. Both strands are
    searched, so the sequence as copied from a design tool will be found
    regardless of which strand of the locus the user happens to have.
    """
    if nuclease not in NUCLEASES:
        raise DesignError(f"unknown nuclease {nuclease!r}; choose from {sorted(NUCLEASES)}")
    spec = NUCLEASES[nuclease]
    template = template.upper()
    spacer = protospacer.upper().replace(" ", "")
    if not spacer:
        raise DesignError("protospacer is empty")

    pam_pattern: str = spec["pam"]
    pam_side: int = spec["pam_side"]
    cut_offset: int = spec["cut_offset"]
    sites: list[CutSite] = []

    for strand in ("+", "-"):
        needle = spacer if strand == "+" else revcomp(spacer)
        start = 0
        while True:
            idx = template.find(needle, start)
            if idx < 0:
                break
            start = idx + 1
            end = idx + len(needle)
            if strand == "+":
                if pam_side == 3:
                    pam_seq = template[end : end + len(pam_pattern)]
                    pam_ok = _iupac_match(pam_pattern, pam_seq)
                    cut = idx + cut_offset
                else:
                    pam_seq = template[max(0, idx - len(pam_pattern)) : idx]
                    pam_ok = _iupac_match(pam_pattern, pam_seq)
                    cut = idx + cut_offset
            else:
                # Protospacer lies on the minus strand; PAM and cut mirror.
                if pam_side == 3:
                    pam_window = template[max(0, idx - len(pam_pattern)) : idx]
                    pam_seq = revcomp(pam_window)
                    pam_ok = _iupac_match(pam_pattern, pam_seq)
                    cut = end - cut_offset
                else:
                    pam_window = template[end : end + len(pam_pattern)]
                    pam_seq = revcomp(pam_window)
                    pam_ok = _iupac_match(pam_pattern, pam_seq)
                    cut = end - cut_offset
            if require_pam and not pam_ok:
                continue
            sites.append(
                CutSite(
                    position=max(0, min(len(template), cut)),
                    strand=strand,
                    protospacer=needle,
                    pam=pam_seq,
                    nuclease=nuclease,
                    spacer_start=idx,
                    spacer_end=end,
                )
            )
    return sites


@dataclass(slots=True)
class GenotypingGeometry:
    """How far from the lesion the primers are allowed to sit."""

    #: Bases each side of the cut that no primer may touch. Should exceed the
    #: largest indel expected, since a primer inside the deleted span will not
    #: amplify the mutant allele at all.
    protected: int = 40
    #: Minimum distance from the cut to the nearest primer base.
    min_flank: int = 40
    #: Maximum distance from the cut to the furthest primer base.
    max_flank: int = 400

    def validate(self) -> None:
        if self.protected < 0:
            raise DesignError("protected window cannot be negative")
        if self.min_flank < self.protected:
            self.min_flank = self.protected
        if self.max_flank <= self.min_flank:
            raise DesignError(
                f"max_flank ({self.max_flank}) must exceed min_flank ({self.min_flank})"
            )


@dataclass(slots=True)
class GenotypingDesign:
    """Result of a genotyping design, with predicted allele sizes."""

    result: DesignResult
    cut_site: CutSite
    geometry: GenotypingGeometry
    expected_indels: tuple[int, ...] = ()
    allele_sizes: list[dict[str, Any]] = field(default_factory=list)

    @property
    def pairs(self):  # noqa: ANN201 - mirrors DesignResult
        return self.result.pairs

    def to_dict(self) -> dict[str, Any]:
        return {
            "cut_site": self.cut_site.to_dict(),
            "geometry": {
                "protected_window": self.geometry.protected,
                "min_flank": self.geometry.min_flank,
                "max_flank": self.geometry.max_flank,
            },
            "allele_sizes": self.allele_sizes,
        }


def predict_allele_sizes(
    wt_size: int, indels: tuple[int, ...] = (-10, -5, -2, -1, 1, 2, 5, 10)
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
    cut_site: CutSite | int,
    params: DesignParams,
    *,
    geometry: GenotypingGeometry | None = None,
    panel: ParalogPanel | None = None,
    variants: IntervalSet | None = None,
    excluded: IntervalSet | None = None,
    expected_indels: tuple[int, ...] = (-10, -5, -2, -1, 1, 2, 5, 10),
    name: str | None = None,
) -> GenotypingDesign:
    """Design a pair that brackets an edited site.

    The protected window is both excluded from primer placement and required to
    be spanned by the amplicon, which is what distinguishes this from simply
    designing primers "near" the cut.
    """
    site = CutSite(position=cut_site, strand="+") if isinstance(cut_site, int) else cut_site
    geom = geometry or GenotypingGeometry()
    geom.validate()

    n = len(template.seq)
    if not 0 <= site.position <= n:
        raise DesignError(f"cut site {site.position + 1} lies outside the template (1-{n})")

    protected = Interval(
        max(0, site.position - geom.protected), min(n, site.position + geom.protected)
    )
    f_lo = max(0, site.position - geom.max_flank)
    f_hi = max(0, site.position - geom.min_flank)
    r_lo = min(n, site.position + geom.min_flank)
    r_hi = min(n, site.position + geom.max_flank)

    if f_hi - f_lo < params.min_length:
        raise DesignError(
            f"only {max(0, f_hi - f_lo)} nt of template available 5' of the cut site; "
            f"need at least {params.min_length} nt. Supply more flanking sequence or "
            f"reduce --min-flank (currently {geom.min_flank})."
        )
    if r_hi - r_lo < params.min_length:
        raise DesignError(
            f"only {max(0, r_hi - r_lo)} nt of template available 3' of the cut site; "
            f"need at least {params.min_length} nt. Supply more flanking sequence or "
            f"reduce --min-flank (currently {geom.min_flank})."
        )

    request = DesignRequest(
        template=template,
        params=params,
        forward_search=IntervalSet([Interval(f_lo, f_hi)]),
        reverse_search=IntervalSet([Interval(r_lo, r_hi)]),
        targets=IntervalSet([protected]),
        excluded=(excluded or IntervalSet()).union(IntervalSet([protected])),
        variants=variants or IntervalSet(),
        panel=panel,
        name=name or f"{template.name}_gt",
        context={"assay": "genotyping", "cut_site": site.to_dict()},
    )
    result = run_design(request)

    design = GenotypingDesign(
        result=result, cut_site=site, geometry=geom, expected_indels=expected_indels
    )
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
