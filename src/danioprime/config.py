"""Design parameters, penalty weights and zebrafish assay presets.

Every tunable lives in :class:`DesignParams`. Presets are named bundles matched
to assays zebrafish labs actually run; they are starting points, not
prescriptions, and any field can be overridden from the CLI or a JSON file.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path
from typing import Any

from .errors import DanioPrimeError
from .thermo.nn import SaltConditions


@dataclass(frozen=True, slots=True)
class PenaltyWeights:
    """Weights of the additive penalty model.

    Penalties are in arbitrary units; only their ratios matter. The defaults
    reproduce Primer3's qualitative ordering (Tm deviation dominates, then
    size, then structure) with extra terms for the zebrafish-specific
    criteria that have no Primer3 equivalent.
    """

    tm_deviation: float = 1.0
    tm_pair_difference: float = 1.0
    gc_deviation: float = 0.2
    length_deviation: float = 0.2
    product_size_deviation: float = 0.05
    hairpin_dg: float = 0.5
    self_dimer_dg: float = 0.5
    hetero_dimer_dg: float = 0.75
    three_prime_dimer: float = 1.0
    end_stability: float = 0.3
    gc_clamp_missing: float = 1.0
    poly_x: float = 0.5
    dinucleotide_repeat: float = 0.5
    ambiguous_base: float = 3.0
    softmask_overlap: float = 1.0
    #: Per variant overlapping the primer footprint.
    variant_overlap: float = 4.0
    #: Extra penalty per variant in the 3'-terminal window.
    variant_three_prime: float = 10.0
    #: Per off-target binding site predicted to prime.
    offtarget_site: float = 2.0
    #: Per predicted off-target amplicon in the size window.
    offtarget_amplicon: float = 25.0
    #: Per paralog predicted to be co-amplified.
    paralog_amplicon: float = 40.0
    #: Applied when a primer binds a paralog well enough to extend.
    paralog_priming_site: float = 5.0
    #: Reward (negative penalty) for a primer that straddles an exon junction.
    junction_spanning_bonus: float = -5.0


@dataclass(frozen=True, slots=True)
class DesignParams:
    """Complete specification of a design run."""

    # --- primer geometry -------------------------------------------------
    min_length: int = 18
    opt_length: int = 20
    max_length: int = 27

    # --- thermodynamics --------------------------------------------------
    min_tm: float = 57.0
    opt_tm: float = 60.0
    max_tm: float = 63.0
    max_pair_tm_difference: float = 2.0
    min_gc: float = 35.0
    opt_gc: float = 50.0
    max_gc: float = 65.0

    # --- product ---------------------------------------------------------
    min_product_size: int = 150
    opt_product_size: int = 400
    max_product_size: int = 800

    # --- structure limits (kcal/mol; less negative is safer) --------------
    max_hairpin_dg: float = -3.0
    max_self_dimer_dg: float = -6.0
    max_hetero_dimer_dg: float = -8.0
    max_three_prime_dimer_score: float = 6.0
    max_end_stability: float = -9.0

    # --- composition -----------------------------------------------------
    require_gc_clamp: int = 1
    max_gc_clamp: int = 5
    max_poly_x: int = 4
    max_dinucleotide_repeat: int = 4
    allow_ambiguous: bool = False
    max_softmask_fraction: float = 0.5

    # --- zebrafish-specific ----------------------------------------------
    #: Primers overlapping a known variant are discarded outright when True,
    #: otherwise they are penalised and reported.
    reject_variant_overlap: bool = False
    #: Size of the 3' window in which a variant is considered disqualifying.
    variant_three_prime_window: int = 5
    #: A primer is treated as able to prime an off-target/paralog site when it
    #: has no mismatch in this many terminal 3' bases ...
    priming_three_prime_window: int = 5
    #: ... and no more than this many mismatches overall.
    priming_max_mismatches: int = 5
    #: ... and the imperfect-duplex Tm is at least this far below the on-target
    #: Tm to be considered safe.
    paralog_min_tm_margin: float = 5.0
    #: Discard pairs predicted to amplify a paralog at all.
    reject_paralog_amplicon: bool = True
    #: Discard any primer that can prime on the avoid-panel, even when no
    #: product is predicted. Off by default: a single priming primer gives
    #: linear extension, not exponential amplification, so it does not produce
    #: a band and rejecting it can leave no design at all. Worth turning on for
    #: multiplex reactions and for qPCR, where primer consumed on the wrong
    #: template still distorts the measurement.
    reject_paralog_priming: bool = False

    # --- search behaviour -------------------------------------------------
    max_candidates_per_side: int = 2000
    max_pairs: int = 20
    #: Minimum distance between the 3' ends of two reported pairs, so the
    #: output is a set of genuinely different options rather than 20 primers
    #: shifted by one base.
    pair_diversity_window: int = 8

    # --- reaction conditions ---------------------------------------------
    conditions: SaltConditions = SaltConditions()
    anneal_temp_c: float = 37.0
    thermo_backend: str = "builtin"

    weights: PenaltyWeights = PenaltyWeights()

    def validate(self) -> None:
        """Raise :class:`DanioPrimeError` on contradictory settings."""
        problems: list[str] = []
        if not self.min_length <= self.opt_length <= self.max_length:
            problems.append("primer length must satisfy min <= opt <= max")
        if self.min_length < 8:
            problems.append("min_length below 8 nt is not usable for PCR")
        if not self.min_tm <= self.opt_tm <= self.max_tm:
            problems.append("Tm must satisfy min <= opt <= max")
        if not self.min_gc <= self.opt_gc <= self.max_gc:
            problems.append("GC%% must satisfy min <= opt <= max")
        if not self.min_product_size <= self.max_product_size:
            problems.append("min_product_size exceeds max_product_size")
        if not (self.min_product_size <= self.opt_product_size <= self.max_product_size):
            problems.append("opt_product_size lies outside the product size range")
        if self.min_product_size < 2 * self.min_length:
            problems.append("min_product_size is smaller than two primers")
        if self.max_pairs < 1:
            problems.append("max_pairs must be at least 1")
        if problems:
            raise DanioPrimeError("; ".join(problems))

    def replace(self, **kwargs: Any) -> DesignParams:
        """Return a copy with *kwargs* overridden."""
        known = {f.name for f in fields(self)}
        unknown = set(kwargs) - known
        if unknown:
            raise DanioPrimeError(f"unknown design parameters: {sorted(unknown)}")
        return replace(self, **kwargs)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2) + "\n")

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DesignParams:
        data = dict(data)
        conditions = data.pop("conditions", None)
        weights = data.pop("weights", None)
        base = cls()
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise DanioPrimeError(f"unknown design parameters: {sorted(unknown)}")
        params = replace(base, **data)
        if conditions:
            params = replace(params, conditions=SaltConditions(**conditions))
        if weights:
            params = replace(params, weights=PenaltyWeights(**weights))
        return params

    @classmethod
    def from_json(cls, path: str | Path) -> DesignParams:
        try:
            data = json.loads(Path(path).read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise DanioPrimeError(f"could not read parameter file {path}: {exc}") from exc
        return cls.from_dict(data)


#: Named parameter bundles. Keys are the values accepted by ``--preset``.
PRESETS: dict[str, dict[str, Any]] = {
    "genotyping-gel": dict(
        min_product_size=180,
        opt_product_size=320,
        max_product_size=600,
        min_tm=57.0,
        opt_tm=60.0,
        max_tm=63.0,
        max_poly_x=4,
    ),
    "genotyping-fluorescent": dict(
        # Fragment analysis needs short amplicons that size-separate cleanly,
        # and the forward primer carries an M13F tail (see assays.tails).
        min_product_size=90,
        opt_product_size=180,
        max_product_size=300,
        min_tm=58.0,
        opt_tm=60.0,
        max_tm=62.0,
        max_pair_tm_difference=1.5,
    ),
    "qpcr": dict(
        min_product_size=70,
        opt_product_size=110,
        max_product_size=160,
        min_length=18,
        opt_length=20,
        max_length=24,
        min_tm=59.0,
        opt_tm=60.0,
        max_tm=62.0,
        max_pair_tm_difference=1.0,
        min_gc=40.0,
        max_gc=60.0,
        max_poly_x=3,
        max_hairpin_dg=-2.0,
        max_self_dimer_dg=-5.0,
        max_hetero_dimer_dg=-6.0,
        require_gc_clamp=1,
    ),
    "sanger": dict(
        # Amplicon long enough that the read clears the primer and any indel.
        min_product_size=500,
        opt_product_size=800,
        max_product_size=1200,
        min_tm=56.0,
        opt_tm=59.0,
        max_tm=63.0,
    ),
    "amplicon-seq": dict(
        # Illumina-length amplicons for MiSeq indel quantification.
        min_product_size=180,
        opt_product_size=250,
        max_product_size=280,
        min_tm=58.0,
        opt_tm=60.0,
        max_tm=63.0,
        max_pair_tm_difference=1.5,
    ),
    "insitu-probe": dict(
        # Template for a riboprobe: long product, relaxed structure limits.
        min_product_size=500,
        opt_product_size=900,
        max_product_size=1500,
        min_tm=55.0,
        opt_tm=59.0,
        max_tm=64.0,
        max_hairpin_dg=-5.0,
    ),
}

PRESET_DESCRIPTIONS: dict[str, str] = {
    "genotyping-gel": "Indel genotyping resolved on an agarose/PAGE gel (180-600 bp).",
    "genotyping-fluorescent": "Short amplicons for fluorescent PCR / fragment analysis (90-300 bp).",
    "qpcr": "RT-qPCR: 70-160 bp, tight Tm window, strict structure limits.",
    "sanger": "Sanger sequencing across an edited site (500-1200 bp).",
    "amplicon-seq": "Short-read amplicon sequencing of an edited locus (180-280 bp).",
    "insitu-probe": "Riboprobe template for in situ (500-1500 bp); see the `insitu` mode.",
}


def get_preset(name: str) -> DesignParams:
    """Build :class:`DesignParams` from a named preset."""
    try:
        overrides = PRESETS[name]
    except KeyError as exc:
        raise DanioPrimeError(
            f"unknown preset {name!r}; choose from {sorted(PRESETS)}"
        ) from exc
    params = DesignParams().replace(**overrides)
    params.validate()
    return params
