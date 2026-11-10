"""Candidate generation, scoring, pairing and the design driver."""

from .candidate import (
    CandidateSet,
    Orientation,
    Primer,
    TemplateContext,
    enumerate_candidates,
    evaluate_primer,
)
from .pairing import PrimerPair, build_pairs, diversify, name_pairs, score_pair
from .pipeline import DesignRequest, DesignResult, run_design, windows_around
from .regions import Interval, IntervalSet, parse_interval, read_bed
from .scoring import score_primer

__all__ = [
    "CandidateSet",
    "DesignRequest",
    "DesignResult",
    "Interval",
    "IntervalSet",
    "Orientation",
    "Primer",
    "PrimerPair",
    "TemplateContext",
    "build_pairs",
    "diversify",
    "enumerate_candidates",
    "evaluate_primer",
    "name_pairs",
    "parse_interval",
    "read_bed",
    "run_design",
    "score_pair",
    "score_primer",
    "windows_around",
]
