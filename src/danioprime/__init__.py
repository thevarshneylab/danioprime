"""danioprime: zebrafish-aware PCR primer design.

What this package adds over a generic primer designer:

* **Ohnolog discrimination.** The teleost genome duplication left zebrafish with
  thousands of retained paralog pairs. A primer pair designed against one gene
  model frequently amplifies its ``a``/``b`` partner, and nothing in a standard
  workflow flags it because the paralog is a real gene, not a repeat.
* **Strain-variant masking.** TU, AB, TL and WIK are not isogenic. A primer on a
  polymorphic SNP fails on one allele and looks like a genotype.
* **Assay geometry.** Genotyping must span the lesion without touching it; qPCR
  wants a junction-spanning primer; small indels need a dCAPS digest. These are
  the steps normally done by hand.

Quick start::

    from danioprime import DesignParams, DesignRequest, run_design, read_single_fasta

    template = read_single_fasta("locus.fa")
    result = run_design(DesignRequest(template=template, params=DesignParams()))
    for pair in result.pairs:
        print(pair.name, pair.product_size, pair.forward.seq, pair.reverse.seq)
"""

from __future__ import annotations

__version__ = "0.1.0"

from .config import PRESETS, DesignParams, PenaltyWeights, get_preset
from .design.candidate import Primer
from .design.pairing import PrimerPair
from .design.pipeline import DesignRequest, DesignResult, run_design
from .design.regions import Interval, IntervalSet
from .errors import (
    AnnotationError,
    DanioPrimeError,
    DesignError,
    ExternalToolError,
    RegionError,
    SequenceError,
)
from .sequence import SeqRecord, gc_percent, read_fasta, read_single_fasta, revcomp
from .specificity.paralog import ParalogPanel, assess_pair, assess_primer
from .thermo.nn import SaltConditions, delta_g, tm

__all__ = [
    "PRESETS",
    "AnnotationError",
    "DanioPrimeError",
    "DesignError",
    "DesignParams",
    "DesignRequest",
    "DesignResult",
    "ExternalToolError",
    "Interval",
    "IntervalSet",
    "ParalogPanel",
    "PenaltyWeights",
    "Primer",
    "PrimerPair",
    "RegionError",
    "SaltConditions",
    "SeqRecord",
    "SequenceError",
    "__version__",
    "assess_pair",
    "assess_primer",
    "delta_g",
    "gc_percent",
    "get_preset",
    "read_fasta",
    "read_single_fasta",
    "revcomp",
    "run_design",
    "tm",
]
