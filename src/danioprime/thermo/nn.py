"""Nearest-neighbour duplex thermodynamics.

Watson-Crick stacking parameters are the unified set of SantaLucia & Hicks
(2004), *Annu. Rev. Biophys. Biomol. Struct.* 33:415-440, in kcal/mol (dH) and
cal/(mol K) (dS). Internal single mismatch parameters are from the Allawi &
SantaLucia series (1997-1998) and are used only to estimate the stability of
imperfect primer-template duplexes when assessing paralog cross-reactivity.

Salt correction follows SantaLucia & Hicks (2004) eq. 6, applied to dS:

    dS_salt = dS + 0.368 * (N - 1) * ln([Na+]_eq)

Divalent cations and dNTPs are folded into an equivalent monovalent
concentration with the von Ahsen approximation used by Primer3:

    [Na+]_eq = [Na+] + [K+] + [Tris]/2 + 120 * sqrt([Mg2+] - [dNTP])
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..errors import SequenceError

#: Gas constant, cal/(mol K).
R_CAL = 1.9872

#: Watson-Crick nearest-neighbour stacks: (dH kcal/mol, dS cal/(mol K)).
NN_WC: dict[str, tuple[float, float]] = {
    "AA": (-7.6, -21.3),
    "TT": (-7.6, -21.3),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

#: Helix initiation terms keyed by the identity of the terminal base pair.
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)

#: Entropic penalty applied to self-complementary (palindromic) duplexes.
SYMMETRY_DS = -1.4

#: Internal single-mismatch stacks, keyed "XY/ZW" where XY is the top-strand
#: dinucleotide 5'->3' and ZW the bottom-strand dinucleotide 3'->5'.
#: Values from Allawi & SantaLucia 1997-1998 (NAR/Biochemistry series).
NN_MISMATCH: dict[str, tuple[float, float]] = {
    # A.A
    "AA/TA": (1.2, 1.7), "CA/GA": (-0.9, -4.2), "GA/CA": (-2.9, -9.8),
    "TA/AA": (4.7, 12.9), "AA/AT": (2.3, 4.6), "AC/AG": (5.3, 14.6),
    "AG/AC": (1.9, 3.7), "AT/AA": (3.4, 8.0),
    # C.C
    "AC/TC": (0.0, -4.4), "CC/GC": (-1.5, -7.2), "GC/CC": (3.6, 8.9),
    "TC/AC": (6.1, 16.4), "CA/CT": (1.9, 3.7), "CC/CG": (-0.7, -3.8),
    "CG/CC": (-0.6, -0.6), "CT/CA": (-0.7, -3.8),
    # G.G
    "AG/TG": (-3.1, -9.5), "CG/GG": (-4.9, -15.3), "GG/CG": (-6.0, -15.8),
    "TG/AG": (1.6, 3.6), "GA/GT": (-0.1, -1.8), "GC/GG": (-4.3, -11.1),
    "GG/GC": (-4.9, -15.3), "GT/GA": (-1.1, -2.1),
    # T.T
    "AT/TT": (-2.7, -10.8), "CT/GT": (-5.0, -15.8), "GT/CT": (-2.2, -8.4),
    "TT/AT": (0.2, -1.5), "TA/TT": (-1.3, -5.3), "TC/TG": (-0.1, -1.7),
    "TG/TC": (-1.9, -6.3), "TT/TA": (-2.7, -10.8),
    # G.T wobbles
    "AG/TT": (1.0, 0.9), "AT/TG": (-2.5, -8.3), "CG/GT": (-4.1, -11.7),
    "CT/GG": (-2.8, -8.0), "GG/CT": (3.3, 10.4), "GT/CG": (-4.4, -12.3),
    "GG/TC": (5.8, 16.3), "GT/TG": (4.1, 9.5), "TG/AT": (-0.1, -1.7),
    "TT/AG": (-1.3, -5.3), "TG/GT": (-1.4, -6.2), "TT/GA": (-1.3, -5.3),
    # G.A / A.C / C.T mixed, abbreviated set
    "AA/TC": (2.3, 4.6), "AC/TA": (5.3, 14.6), "CA/GC": (1.9, 3.7),
    "CC/GA": (0.6, -0.6), "GA/CC": (5.2, 14.2), "GC/CA": (-0.7, -3.8),
    "TA/AC": (3.4, 8.0), "TC/AA": (7.6, 20.2),
    "AA/TG": (-0.6, -2.3), "AG/TA": (-0.7, -2.3), "CA/GG": (-0.7, -2.3),
    "CG/GA": (-4.0, -13.2), "GA/CG": (-0.6, -1.0), "GG/CA": (0.5, 3.2),
    "TA/AG": (0.7, 0.7), "TG/AA": (3.0, 7.4),
}

#: Terminal mismatch / dangling-end approximation: a blunt mismatch at the end
#: of a duplex contributes essentially no stacking stability.
TERMINAL_MISMATCH = (0.0, 0.0)


@dataclass(frozen=True, slots=True)
class SaltConditions:
    """Ionic and oligo conditions of the reaction (all concentrations in mM
    except *primer_conc* / *dntp_conc*, which are nM and mM respectively)."""

    monovalent_mM: float = 50.0
    divalent_mM: float = 1.5
    dntp_mM: float = 0.6
    primer_nM: float = 250.0
    dmso_pct: float = 0.0
    formamide_M: float = 0.0

    def equivalent_monovalent_mM(self) -> float:
        """Fold divalent cations and dNTPs into an equivalent [Na+]."""
        div = self.divalent_mM
        dntp = self.dntp_mM
        if div <= 0:
            return max(self.monovalent_mM, 1e-6)
        if div < dntp:
            div = dntp
        free = div - dntp
        return max(self.monovalent_mM + 120.0 * math.sqrt(free) if free > 0 else self.monovalent_mM, 1e-6)


DEFAULT_CONDITIONS = SaltConditions()


def _validate(seq: str) -> str:
    s = seq.upper()
    if len(s) < 2:
        raise SequenceError("nearest-neighbour model needs at least 2 nt")
    bad = set(s) - set("ACGT")
    if bad:
        raise SequenceError(
            f"nearest-neighbour model requires unambiguous bases; found {sorted(bad)}"
        )
    return s


def is_self_complementary(seq: str) -> bool:
    """True if *seq* is its own reverse complement."""
    s = seq.upper()
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    try:
        return s == "".join(comp[b] for b in reversed(s))
    except KeyError:
        return False


def duplex_enthalpy_entropy(seq: str) -> tuple[float, float]:
    """Return (dH kcal/mol, dS cal/(mol K)) for a perfectly matched duplex."""
    s = _validate(seq)
    dh = 0.0
    ds = 0.0
    for i in range(len(s) - 1):
        step = s[i : i + 2]
        sdh, sds = NN_WC[step]
        dh += sdh
        ds += sds
    for terminal in (s[0], s[-1]):
        idh, ids = INIT_GC if terminal in "GC" else INIT_AT
        dh += idh
        ds += ids
    if is_self_complementary(s):
        ds += SYMMETRY_DS
    return dh, ds


def salt_corrected_entropy(ds: float, length: int, conditions: SaltConditions) -> float:
    """Apply the SantaLucia (2004) monovalent entropy correction."""
    na = conditions.equivalent_monovalent_mM() / 1000.0
    return ds + 0.368 * (length - 1) * math.log(na)


def tm(seq: str, conditions: SaltConditions = DEFAULT_CONDITIONS) -> float:
    """Melting temperature in degrees Celsius for a primer/template duplex.

    Uses total strand concentration ``CT = [primer]`` (primer in large excess
    over template, the usual PCR regime) and doubles it for self-complementary
    oligos as required by the two-state model.
    """
    s = _validate(seq)
    dh, ds = duplex_enthalpy_entropy(s)
    ds = salt_corrected_entropy(ds, len(s), conditions)
    ct = conditions.primer_nM * 1e-9
    if is_self_complementary(s):
        ct *= 4.0
    if ct <= 0:
        raise ValueError("primer concentration must be positive")
    denom = ds + R_CAL * math.log(ct)
    if denom == 0:
        raise ValueError("degenerate thermodynamic denominator")
    kelvin = (dh * 1000.0) / denom
    celsius = kelvin - 273.15
    # Co-solvent corrections (Chester & Marshak 1993; Blake & Delcourt 1996).
    celsius -= 0.75 * conditions.dmso_pct
    celsius -= 0.65 * conditions.formamide_M
    return celsius


def delta_g(seq: str, temp_c: float = 37.0, conditions: SaltConditions = DEFAULT_CONDITIONS) -> float:
    """Duplex free energy in kcal/mol at *temp_c*, salt corrected."""
    dh, ds = duplex_enthalpy_entropy(seq)
    ds = salt_corrected_entropy(ds, len(seq), conditions)
    t = temp_c + 273.15
    return dh - t * (ds / 1000.0)


def end_stability(seq: str, n: int = 5, temp_c: float = 37.0) -> float:
    """Free energy (kcal/mol) of the 3'-terminal *n* bases.

    Strongly negative values indicate a 3' end that is prone to non-specific
    priming; Primer3 exposes the same quantity as ``PRIMER_MAX_END_STABILITY``.
    """
    s = _validate(seq)
    tail = s[-n:] if len(s) >= n else s
    dh = 0.0
    ds = 0.0
    for i in range(len(tail) - 1):
        sdh, sds = NN_WC[tail[i : i + 2]]
        dh += sdh
        ds += sds
    t = temp_c + 273.15
    return dh - t * (ds / 1000.0)


def mismatch_stack(top: str, bottom: str) -> tuple[float, float]:
    """Stacking parameters for a dinucleotide step containing a mismatch.

    *top* is the 5'->3' dinucleotide of the primer; *bottom* is the opposing
    3'->5' dinucleotide of the template. Unparameterised combinations fall back
    to a mildly destabilising default.
    """
    key = f"{top.upper()}/{bottom.upper()}"
    if key in NN_MISMATCH:
        return NN_MISMATCH[key]
    return (0.0, 2.0)


_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def imperfect_duplex_tm(
    primer: str,
    template_rc: str,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
) -> float:
    """Estimate Tm of a primer annealed to an imperfectly matched site.

    *template_rc* must be the template strand written 5'->3' in the same
    register as *primer*, so that position *i* of each string is a base pair
    (matched when ``template_rc[i] == primer[i]``). Mismatched steps draw on
    the internal-mismatch tables; terminal mismatches contribute no stacking.

    This is a coarse estimate, adequate for ranking paralog cross-reactivity,
    not a substitute for a full partition-function calculation.
    """
    p = primer.upper()
    t = template_rc.upper()
    if len(p) != len(t):
        raise SequenceError("primer and template registers must be the same length")
    if len(p) < 2:
        raise SequenceError("need at least 2 nt")

    dh = 0.0
    ds = 0.0
    for i in range(len(p) - 1):
        p_step = p[i : i + 2]
        t_step = t[i : i + 2]
        left_match = p[i] == t[i]
        right_match = p[i + 1] == t[i + 1]
        if left_match and right_match:
            sdh, sds = NN_WC[p_step]
        elif i == 0 and not left_match:
            sdh, sds = TERMINAL_MISMATCH
        elif i == len(p) - 2 and not right_match:
            sdh, sds = TERMINAL_MISMATCH
        else:
            # Bottom strand written 3'->5' opposite the primer step.
            bottom = "".join(_COMP.get(b, "N") for b in t_step)
            sdh, sds = mismatch_stack(p_step, bottom)
        dh += sdh
        ds += sds

    for idx in (0, -1):
        if p[idx] == t[idx]:
            idh, ids = INIT_GC if p[idx] in "GC" else INIT_AT
        else:
            idh, ids = INIT_AT
        dh += idh
        ds += ids

    ds = salt_corrected_entropy(ds, len(p), conditions)
    ct = conditions.primer_nM * 1e-9
    denom = ds + R_CAL * math.log(ct)
    if dh >= 0 or denom >= 0:
        return -273.15  # duplex not predicted to form
    return (dh * 1000.0) / denom - 273.15
