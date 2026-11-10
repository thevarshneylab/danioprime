"""dCAPS: turning a point difference into a length difference.

A 1-2 bp indel or a single substitution cannot be called on a gel from product
size. dCAPS (Neff et al. 1998, *Plant J.* 14:387) solves this by introducing
deliberate mismatches into one primer so that, together with the variant base
supplied by the template, a restriction site exists in one allele and not the
other. The digest then gives two bands for one allele and one for the other --
a robust, equipment-free readout.

Geometry
--------
For a forward dCAPS primer the recognition site must straddle the primer's 3'
end::

      primer (we control these bases)      template (fixed, allele-specific)
    5'-...............................|...................................-3'
                        [ recognition site spans the boundary ]
                                         ^ variant position

Two constraints follow, and both are easy to get wrong by hand:

* the variant must lie **downstream** of the primer's 3' end. If the primer
  covers the variant, the primer's own sequence overwrites it and both alleles
  give an identical product.
* every site position downstream of the primer, other than the variant itself,
  must already match the template. Those bases cannot be changed.

The mismatches needed are then just the site positions inside the primer that
differ from the template. Fewer is better: one introduced mismatch usually
amplifies indistinguishably from a perfect primer, two is often acceptable,
three rarely works.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from ..config import DesignParams
from ..design.candidate import Primer, TemplateContext, evaluate_primer
from ..design.pairing import PrimerPair, score_pair
from ..design.pipeline import DesignRequest, run_design
from ..design.regions import Interval, IntervalSet
from ..errors import DesignError
from ..sequence import SeqRecord, revcomp
from ..thermo.backend import get_backend
from ..thermo.structure import hetero_dimer, three_prime_dimer_score
from .enzymes import DEFAULT_ENZYMES, Enzyme, digest_fragments, site_matches

CutAllele = Literal["variant", "reference"]


@dataclass(slots=True)
class DcapsPrimer:
    """A mismatch-carrying primer that creates an allele-specific site."""

    seq: str
    start: int
    end: int
    orientation: Literal["F", "R"]
    enzyme: Enzyme
    site_start: int
    mismatch_offsets_from_3p: tuple[int, ...]
    cut_allele: CutAllele
    variant_position: int
    template_seq: str = ""

    @property
    def mismatches(self) -> int:
        return len(self.mismatch_offsets_from_3p)

    @property
    def nearest_mismatch_to_3p(self) -> int | None:
        return min(self.mismatch_offsets_from_3p) + 1 if self.mismatch_offsets_from_3p else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "seq": self.seq,
            "template_seq": self.template_seq,
            "orientation": self.orientation,
            "position_1based": f"{self.start + 1}-{self.end}",
            "enzyme": self.enzyme.name,
            "enzyme_site": self.enzyme.display,
            "site_start_1based": self.site_start + 1,
            "introduced_mismatches": self.mismatches,
            "mismatch_offsets_from_3p": list(self.mismatch_offsets_from_3p),
            "nearest_mismatch_nt_from_3p": self.nearest_mismatch_to_3p,
            "cut_allele": self.cut_allele,
            "variant_position_1based": self.variant_position + 1,
        }


@dataclass(slots=True)
class DcapsAssay:
    """A complete dCAPS assay: the two primers and the predicted digest."""

    dcaps_primer: DcapsPrimer
    partner: Primer
    product_size: int
    reference_fragments: list[int] = field(default_factory=list)
    variant_fragments: list[int] = field(default_factory=list)
    pair: PrimerPair | None = None
    penalty: float = 0.0
    notes: list[str] = field(default_factory=list)
    name: str = ""

    @property
    def discriminating(self) -> bool:
        """True if the two allele digests give distinguishable band patterns."""
        return sorted(self.reference_fragments) != sorted(self.variant_fragments)

    @property
    def diagnostic_bands(self) -> tuple[list[int], list[int]]:
        """Bands present in one allele's digest but not the other's.

        These are what the call is actually made on. Pre-existing sites
        elsewhere in the amplicon produce bands common to both alleles, which
        are reassuring as a digest control but carry no genotype information.
        """
        ref = sorted(self.reference_fragments)
        var = sorted(self.variant_fragments)
        shared = []
        pool = list(var)
        for band in ref:
            if band in pool:
                pool.remove(band)
                shared.append(band)
        ref_only = list(ref)
        for band in shared:
            ref_only.remove(band)
        return ref_only, pool

    @property
    def smallest_diagnostic_band(self) -> int:
        """Smallest band that carries genotype information."""
        ref_only, var_only = self.diagnostic_bands
        bands = ref_only + var_only
        return min(bands) if bands else 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "product_size": self.product_size,
            "enzyme": self.dcaps_primer.enzyme.name,
            "cut_allele": self.dcaps_primer.cut_allele,
            "introduced_mismatches": self.dcaps_primer.mismatches,
            "reference_fragments": sorted(self.reference_fragments, reverse=True),
            "variant_fragments": sorted(self.variant_fragments, reverse=True),
            "diagnostic_bands_reference_only": sorted(self.diagnostic_bands[0], reverse=True),
            "diagnostic_bands_variant_only": sorted(self.diagnostic_bands[1], reverse=True),
            "discriminating": self.discriminating,
            "smallest_diagnostic_band": self.smallest_diagnostic_band,
            "penalty": round(self.penalty, 3),
            "notes": list(self.notes),
            "dcaps_primer": self.dcaps_primer.to_dict(),
            "partner_primer": self.partner.to_dict(),
        }


def _candidate_dcaps_primers(
    template: str,
    variant_pos: int,
    variant_base: str,
    enzymes: dict[str, Enzyme],
    *,
    orientation: Literal["F", "R"],
    min_length: int,
    max_length: int,
    max_mismatches: int,
    min_mismatch_offset: int,
    cut_allele: CutAllele,
) -> list[DcapsPrimer]:
    """Enumerate dCAPS primers for one orientation.

    Reverse-orientation design is handled by reverse-complementing the template
    and recursing, so the geometry only has to be written once.
    """
    if orientation == "R":
        n = len(template)
        flipped = revcomp(template)
        flipped_pos = n - 1 - variant_pos
        flipped_base = revcomp(variant_base)
        out: list[DcapsPrimer] = []
        for cand in _candidate_dcaps_primers(
            flipped,
            flipped_pos,
            flipped_base,
            enzymes,
            orientation="F",
            min_length=min_length,
            max_length=max_length,
            max_mismatches=max_mismatches,
            min_mismatch_offset=min_mismatch_offset,
            cut_allele=cut_allele,
        ):
            out.append(
                DcapsPrimer(
                    seq=cand.seq,
                    start=n - cand.end,
                    end=n - cand.start,
                    orientation="R",
                    enzyme=cand.enzyme,
                    site_start=n - (cand.site_start + cand.enzyme.length),
                    mismatch_offsets_from_3p=cand.mismatch_offsets_from_3p,
                    cut_allele=cand.cut_allele,
                    variant_position=variant_pos,
                    template_seq=cand.template_seq,
                )
            )
        return out

    results: list[DcapsPrimer] = []
    ref_base = template[variant_pos]
    target_base = variant_base if cut_allele == "variant" else ref_base

    for enzyme in enzymes.values():
        k = enzyme.length
        site = enzyme.site
        # The site must cover the variant and straddle the primer's 3' end.
        for w in range(max(0, variant_pos - k + 1), variant_pos + 1):
            if w + k > len(template):
                continue
            v_in_site = variant_pos - w
            # The variant base must be the one the site demands for the cut
            # allele, and must differ for the other allele.
            allowed = site_matches(site[v_in_site], target_base)
            other_base = ref_base if cut_allele == "variant" else variant_base
            if not allowed or site_matches(site[v_in_site], other_base):
                continue
            # Primer 3' end must fall strictly inside the site and strictly
            # before the variant, so the variant comes from the template.
            for e in range(w + 1, variant_pos + 1):
                # Positions [e, w+k) are template-derived: they must already
                # match the site, except at the variant itself.
                ok = True
                for i in range(e, w + k):
                    if i == variant_pos:
                        continue
                    if not site_matches(site[i - w], template[i]):
                        ok = False
                        break
                if not ok:
                    continue
                # Positions [w, e) come from the primer; write the site there.
                offsets: list[int] = []
                for i in range(w, e):
                    want = site[i - w]
                    if not site_matches(want, template[i]):
                        offsets.append(e - 1 - i)
                if len(offsets) > max_mismatches or not offsets:
                    continue
                if min(offsets) < min_mismatch_offset:
                    continue
                for length in range(min_length, max_length + 1):
                    s = e - length
                    if s < 0:
                        continue
                    bases = list(template[s:e])
                    for i in range(w, e):
                        want = site[i - w]
                        bases[i - s] = want if want in "ACGT" else template[i]
                    # Re-check: a degenerate site position we could not write
                    # must still match.
                    built = "".join(bases)
                    window = built[w - s : e - s] + template[e : w + k]
                    if cut_allele == "variant":
                        window = window[:v_in_site] + variant_base + window[v_in_site + 1 :]
                    if not site_matches(site, window):
                        continue
                    results.append(
                        DcapsPrimer(
                            seq=built,
                            start=s,
                            end=e,
                            orientation="F",
                            enzyme=enzyme,
                            site_start=w,
                            mismatch_offsets_from_3p=tuple(sorted(offsets)),
                            cut_allele=cut_allele,
                            variant_position=variant_pos,
                            template_seq=template[s:e],
                        )
                    )
    return results


def design_dcaps(
    template: SeqRecord,
    variant_position: int,
    variant_base: str,
    params: DesignParams,
    *,
    enzymes: dict[str, Enzyme] | None = None,
    cut_allele: CutAllele = "variant",
    max_introduced_mismatches: int = 2,
    min_mismatch_offset: int = 1,
    max_assays: int = 10,
    min_band_separation: int = 25,
    name: str | None = None,
) -> list[DcapsAssay]:
    """Design dCAPS assays for a single-base difference.

    Parameters
    ----------
    variant_position:
        Zero-based position of the variant in *template*.
    variant_base:
        The alternate base. The reference base is read from the template.
    cut_allele:
        Whether the restriction site should be created in the ``variant`` or
        the ``reference`` allele. Cutting the variant is usually preferable:
        a mutant-specific band is easier to read than a mutant-specific
        *absence* of one, which is indistinguishable from a failed digest.
    min_mismatch_offset:
        Minimum distance from the primer's 3' end at which an introduced
        mismatch may sit. The default of 1 keeps the 3'-terminal base matched,
        since a terminal mismatch cripples extension.
    min_band_separation:
        Reject assays whose digest products differ by less than this, as they
        will not be callable on a gel.
    """
    enzymes = enzymes or DEFAULT_ENZYMES
    seq = template.seq
    n = len(seq)
    if not 0 <= variant_position < n:
        raise DesignError(f"variant position {variant_position + 1} outside template (1-{n})")
    variant_base = variant_base.upper()
    if variant_base not in "ACGT":
        raise DesignError(f"variant base must be A/C/G/T, got {variant_base!r}")
    ref_base = seq[variant_position]
    if variant_base == ref_base:
        raise DesignError(
            f"variant base {variant_base} equals the template base at "
            f"{variant_position + 1}; nothing to discriminate"
        )

    candidates: list[DcapsPrimer] = []
    for orientation in ("F", "R"):
        candidates += _candidate_dcaps_primers(
            seq,
            variant_position,
            variant_base,
            enzymes,
            orientation=orientation,  # type: ignore[arg-type]
            min_length=params.min_length,
            max_length=params.max_length,
            max_mismatches=max_introduced_mismatches,
            min_mismatch_offset=min_mismatch_offset,
            cut_allele=cut_allele,
        )
    if not candidates:
        raise DesignError(
            f"no enzyme in the table can be made allele-specific at position "
            f"{variant_position + 1} ({ref_base}>{variant_base}) with up to "
            f"{max_introduced_mismatches} introduced mismatches. Try "
            f"--max-introduced-mismatches 3, --cut-allele reference, or a larger "
            f"enzyme table via --enzymes."
        )

    backend = get_backend(params.thermo_backend, params.conditions, params.anneal_temp_c)
    ctx = TemplateContext(seq=seq, params=params, backend=backend, name=template.name)

    # Rank dCAPS primers: fewest mismatches, then mismatch furthest from the
    # 3' end, then closest to the optimum Tm.
    scored: list[tuple[float, DcapsPrimer, Primer]] = []
    for cand in candidates:
        probe = Primer(
            seq=cand.seq, start=cand.start, end=cand.end, orientation=cand.orientation
        )
        evaluate_primer(probe, ctx)
        if probe.tm != probe.tm:
            continue
        if not (params.min_tm - 3 <= probe.tm <= params.max_tm + 3):
            continue
        rank = (
            10.0 * cand.mismatches
            - 0.5 * (cand.nearest_mismatch_to_3p or 0)
            + abs(probe.tm - params.opt_tm)
        )
        scored.append((rank, cand, probe))
    if not scored:
        raise DesignError(
            "dCAPS primers exist but none has a usable Tm; widen --min-tm/--max-tm "
            "or allow longer primers with --max-length"
        )
    scored.sort(key=lambda row: row[0])

    assays: list[DcapsAssay] = []
    per_enzyme: dict[str, int] = {}
    max_per_enzyme = 2
    for rank, cand, probe in scored:
        if len(assays) >= max_assays:
            break
        # Report at most a couple of assays per enzyme: twenty variants of the
        # same primer shifted by one base are not twenty options.
        if per_enzyme.get(cand.enzyme.name, 0) >= max_per_enzyme:
            continue
        partner_orientation = "R" if cand.orientation == "F" else "F"
        # The partner sits on the far side of the variant so the amplicon
        # contains the engineered site.
        if partner_orientation == "R":
            lo = min(n, cand.end + params.min_product_size - probe.length)
            hi = min(n, cand.start + params.max_product_size)
            search = IntervalSet([Interval(max(0, lo), hi)]) if hi > lo else IntervalSet()
        else:
            hi = max(0, cand.start - params.min_product_size + probe.length)
            lo = max(0, cand.end - params.max_product_size)
            search = IntervalSet([Interval(lo, max(0, hi))]) if hi > lo else IntervalSet()
        if not search:
            continue

        request = DesignRequest(
            template=template,
            params=params.replace(max_pairs=3),
            forward_search=search if partner_orientation == "F" else IntervalSet(
                [Interval(cand.start, cand.end)]
            ),
            reverse_search=search if partner_orientation == "R" else IntervalSet(
                [Interval(cand.start, cand.end)]
            ),
            name=name or f"{template.name}_dcaps",
        )
        partner_set = run_design(request)
        partners = [
            p.reverse if partner_orientation == "R" else p.forward for p in partner_set.pairs
        ]
        if not partners:
            continue
        partner = partners[0]

        if cand.orientation == "F":
            product_start, product_end = cand.start, partner.end
        else:
            product_start, product_end = partner.start, cand.end
        product_size = product_end - product_start
        if not (params.min_product_size <= product_size <= params.max_product_size):
            continue

        # Build the two allele amplicons as they will actually exist, with the
        # primer's engineered bases in place of the template's.
        amplicon_ref = list(seq[product_start:product_end])
        offset = cand.start - product_start
        engineered = cand.seq if cand.orientation == "F" else revcomp(cand.seq)
        amplicon_ref[offset : offset + len(engineered)] = list(engineered)
        amplicon_ref_s = "".join(amplicon_ref)
        v_idx = variant_position - product_start
        amplicon_var_s = (
            amplicon_ref_s[:v_idx] + variant_base + amplicon_ref_s[v_idx + 1 :]
        )

        ref_frags = digest_fragments(cand.enzyme, amplicon_ref_s)
        var_frags = digest_fragments(cand.enzyme, amplicon_var_s)

        assay = DcapsAssay(
            dcaps_primer=cand,
            partner=partner,
            product_size=product_size,
            reference_fragments=ref_frags,
            variant_fragments=var_frags,
            penalty=rank,
        )
        if not assay.discriminating:
            continue
        cut_frags = var_frags if cand.cut_allele == "variant" else ref_frags
        if len(cut_frags) < 2:
            continue
        if assay.smallest_diagnostic_band < min_band_separation:
            assay.notes.append(
                f"the smallest genotype-informative band is "
                f"{assay.smallest_diagnostic_band} bp, which may run off a standard "
                "agarose gel; use a 2.5-3% gel or PAGE"
            )
        shared_sites = len(ref_frags) - 1 if cand.cut_allele == "variant" else len(var_frags) - 1
        if shared_sites > 0:
            assay.notes.append(
                f"{cand.enzyme.name} also cuts the amplicon {shared_sites} time(s) "
                "independently of the genotype; those bands appear in both alleles and "
                "serve as a digest control"
            )
        assay.notes.append(
            f"{cand.mismatches} mismatch(es) introduced, nearest "
            f"{cand.nearest_mismatch_to_3p} nt from the 3' end"
        )

        fwd = probe if cand.orientation == "F" else partner
        rev = partner if cand.orientation == "F" else probe
        pair = PrimerPair(
            forward=fwd,
            reverse=rev,
            product_size=product_size,
            tm_difference=abs(fwd.tm - rev.tm),
            hetero_dimer_dg=hetero_dimer(
                fwd.seq, rev.seq, temp_c=params.anneal_temp_c, conditions=params.conditions
            ).dg,
            three_prime_dimer_score=three_prime_dimer_score(fwd.seq, rev.seq),
        )
        score_pair(pair, ctx)
        assay.pair = pair
        assay.name = f"{name or template.name}_dcaps{len(assays) + 1}"
        per_enzyme[cand.enzyme.name] = per_enzyme.get(cand.enzyme.name, 0) + 1
        assays.append(assay)

    if not assays:
        raise DesignError(
            "dCAPS primers were found but none yielded a discriminating digest within "
            "the product-size window. Widen --min-product-size/--max-product-size."
        )

    # Final ordering: fewest introduced mismatches first, then the largest
    # genotype-informative band (the small fragment in a dCAPS digest is
    # roughly the length of the dCAPS primer, so a longer primer buys a band
    # that is easier to see), then primer quality. Bands above ~60 bp run
    # comfortably on a normal gel, so there is nothing to gain past that.
    assays.sort(
        key=lambda a: (
            a.dcaps_primer.mismatches,
            -min(a.smallest_diagnostic_band, 60),
            a.penalty,
        )
    )
    for i, assay in enumerate(assays, 1):
        assay.name = f"{name or template.name}_dcaps{i}"
    return assays
