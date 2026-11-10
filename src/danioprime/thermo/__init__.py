"""Thermodynamic calculations for primer design."""

from .backend import (
    BACKEND_NAMES,
    BuiltinBackend,
    Primer3Backend,
    StructureResult,
    ThermoBackend,
    get_backend,
    primer3_available,
)
from .nn import (
    DEFAULT_CONDITIONS,
    SaltConditions,
    delta_g,
    duplex_enthalpy_entropy,
    end_stability,
    imperfect_duplex_tm,
    is_self_complementary,
    tm,
)
from .structure import (
    Interaction,
    hairpin,
    hetero_dimer,
    self_dimer,
    three_prime_dimer_score,
)

__all__ = [
    "BACKEND_NAMES",
    "DEFAULT_CONDITIONS",
    "BuiltinBackend",
    "Interaction",
    "Primer3Backend",
    "SaltConditions",
    "StructureResult",
    "ThermoBackend",
    "delta_g",
    "duplex_enthalpy_entropy",
    "end_stability",
    "get_backend",
    "hairpin",
    "hetero_dimer",
    "imperfect_duplex_tm",
    "is_self_complementary",
    "primer3_available",
    "self_dimer",
    "three_prime_dimer_score",
    "tm",
]
