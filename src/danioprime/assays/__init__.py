"""Assay-specific design modes."""

from .dcaps import DcapsAssay, DcapsPrimer, design_dcaps
from .enzymes import DEFAULT_ENZYMES, Enzyme, digest_fragments, find_sites, load_enzymes
from .genotyping import (
    NUCLEASES,
    CutSite,
    GenotypingDesign,
    GenotypingGeometry,
    design_genotyping,
    design_junction_pcr,
    locate_protospacer,
    predict_allele_sizes,
)
from .insitu import (
    ConservedBlock,
    CrossHybridisation,
    InSituDesign,
    assess_probe,
    design_insitu,
    find_conserved_blocks,
)
from .qpcr import QpcrDesign, assess_gdna_protection, design_qpcr, intron_sizes_from_transcript
from .tails import (
    RIBOPROBE_SCHEMES,
    SCHEME_DESCRIPTIONS,
    TAIL_SCHEMES,
    TAILS,
    OrderItem,
    ReactionPlan,
    Tail,
    apply_scheme,
    order_sheet,
    reaction_plan,
)

__all__ = [
    "DEFAULT_ENZYMES",
    "RIBOPROBE_SCHEMES",
    "NUCLEASES",
    "SCHEME_DESCRIPTIONS",
    "TAILS",
    "TAIL_SCHEMES",
    "ConservedBlock",
    "CrossHybridisation",
    "CutSite",
    "DcapsAssay",
    "DcapsPrimer",
    "Enzyme",
    "InSituDesign",
    "GenotypingDesign",
    "GenotypingGeometry",
    "OrderItem",
    "QpcrDesign",
    "ReactionPlan",
    "Tail",
    "apply_scheme",
    "assess_gdna_protection",
    "assess_probe",
    "design_dcaps",
    "design_genotyping",
    "design_insitu",
    "design_junction_pcr",
    "design_qpcr",
    "digest_fragments",
    "find_conserved_blocks",
    "find_sites",
    "intron_sizes_from_transcript",
    "load_enzymes",
    "locate_protospacer",
    "order_sheet",
    "predict_allele_sizes",
    "reaction_plan",
]
