"""Exception hierarchy for danioprime."""


class DanioPrimeError(Exception):
    """Base class for all danioprime errors."""


class SequenceError(DanioPrimeError):
    """Malformed or unusable nucleotide sequence."""


class RegionError(DanioPrimeError):
    """Invalid or contradictory interval specification."""


class DesignError(DanioPrimeError):
    """Design could not be completed under the supplied constraints."""


class AnnotationError(DanioPrimeError):
    """Malformed annotation input (GTF/GFF/BED/VCF)."""


class ExternalToolError(DanioPrimeError):
    """An optional external executable failed or is unavailable."""
