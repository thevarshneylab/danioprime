"""Specificity: where a primer can actually prime, and what it will amplify."""

from .index import DEFAULT_SEED, KmerIndex
from .paralog import (
    Discrimination,
    ParalogPanel,
    assess_pair,
    assess_primer,
    summarise_panel_amplicons,
)
from .scan import (
    BindingSite,
    PredictedAmplicon,
    PrimingRules,
    find_binding_sites,
    find_binding_sites_seeded,
    predict_amplicons,
)

__all__ = [
    "DEFAULT_SEED",
    "BindingSite",
    "Discrimination",
    "KmerIndex",
    "ParalogPanel",
    "PredictedAmplicon",
    "PrimingRules",
    "assess_pair",
    "assess_primer",
    "find_binding_sites",
    "find_binding_sites_seeded",
    "predict_amplicons",
    "summarise_panel_amplicons",
]
