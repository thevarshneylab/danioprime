"""Command-line interface.

Subcommands mirror the assays rather than the machinery, because that is how
the work arrives: "I need to genotype this crispant", not "I need a primer with
these thermodynamic properties".

    danioprime genotype  --template locus.fa --target 591-767 --avoid ohnolog.fa
    danioprime qpcr      --transcript tx.fa --gtf ann.gtf.gz --transcript-id ENSDART...
    danioprime dcaps     --template locus.fa --variant 601 --variant-base A
    danioprime junction  --template locus.fa --cassette pGTag.fa --insert-at 1200
    danioprime design    --template any.fa --target 700-760
    danioprime check     --primers F.seq,R.seq --avoid ohnolog.fa
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from . import __version__
from .annotation.gtf import (
    canonical_transcript,
    find_transcripts,
    guess_ohnolog_symbol,
    read_annotation,
)
from .assays.dcaps import design_dcaps
from .assays.enzymes import DEFAULT_ENZYMES, load_enzymes
from .assays.genotyping import (
    DEFAULT_EXPECTED_INDELS,
    design_genotyping,
    design_junction_pcr,
)
from .assays.insitu import design_insitu
from .assays.qpcr import design_qpcr
from .assays.tails import (
    RIBOPROBE_SCHEMES,
    SCHEME_DESCRIPTIONS,
    TAIL_SCHEMES,
    reaction_plan,
)
from .config import PRESET_DESCRIPTIONS, DesignParams, get_preset
from .design.pipeline import DesignRequest, DesignResult, run_design
from .design.regions import IntervalSet, parse_interval
from .errors import DanioPrimeError, RegionError
from .genome.remote import PROVIDERS
from .report.colour import should_colour
from .report.writers import (
    emit,
    format_text,
    render_order_sequence,
    summarise_iterable,
    write_bed,
    write_json,
    write_order_sheet,
    write_tsv,
)
from .sequence import (
    IndexedFasta,
    SeqRecord,
    load_sequence,
    read_fasta,
    soft_mask_runs,
    write_fasta,
)
from .specificity.external import available_aligners
from .specificity.paralog import ParalogPanel, assess_pair, assess_primer
from .thermo.backend import BACKEND_NAMES, primer3_available
from .thermo.nn import SaltConditions
from .variants.remote import clear_index_cache, describe_sources, fetch_region
from .variants.sources import VARIANT_SOURCES
from .variants.vcf import load_mask, summarise


# --------------------------------------------------------------------------
# argument groups shared between subcommands
# --------------------------------------------------------------------------
def _add_template_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("template")
    g.add_argument("--template", "-t", help="FASTA file, or a bare nucleotide sequence")
    g.add_argument(
        "--region",
        help="extract this region from --reference instead, e.g. chr5:1,000,000-1,010,000",
    )
    g.add_argument("--reference", help="indexed reference FASTA to extract --region from")
    g.add_argument(
        "--template-offset",
        type=int,
        default=0,
        help="genomic coordinate of template base 1, so VCF and BED input line up "
        "(set automatically when --region is used)",
    )
    g.add_argument(
        "--respect-softmask",
        action="store_true",
        help="penalise primers overlapping lower-case (repeat-masked) template bases",
    )


def _add_param_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("design parameters")
    g.add_argument(
        "--preset",
        choices=sorted(PRESET_DESCRIPTIONS),
        help="start from a named parameter bundle (see `danioprime presets`)",
    )
    g.add_argument("--params", help="JSON file of design parameters (overrides --preset)")
    g.add_argument("--min-length", type=int, help="minimum primer length")
    g.add_argument("--opt-length", type=int, help="optimal primer length")
    g.add_argument("--max-length", type=int, help="maximum primer length")
    g.add_argument("--min-tm", type=float, help="minimum primer Tm (C)")
    g.add_argument("--opt-tm", type=float, help="optimal primer Tm (C)")
    g.add_argument("--max-tm", type=float, help="maximum primer Tm (C)")
    g.add_argument("--max-pair-tm-difference", type=float, help="maximum Tm difference in a pair")
    g.add_argument("--min-gc", type=float, help="minimum GC%%")
    g.add_argument("--max-gc", type=float, help="maximum GC%%")
    g.add_argument("--min-product-size", type=int, help="minimum amplicon size")
    g.add_argument("--opt-product-size", type=int, help="optimal amplicon size")
    g.add_argument("--max-product-size", type=int, help="maximum amplicon size")
    g.add_argument("--require-gc-clamp", type=int, help="required consecutive 3' G/C bases")
    g.add_argument("--max-poly-x", type=int, help="longest allowed single-base run")
    g.add_argument("--max-pairs", type=int, help="number of pairs to report")
    g.add_argument(
        "--thermo-backend",
        choices=sorted(BACKEND_NAMES),
        help="thermodynamic engine (builtin is dependency-free and reproducible)",
    )
    g.add_argument("--monovalent", type=float, help="monovalent cation concentration (mM)")
    g.add_argument("--divalent", type=float, help="Mg2+ concentration (mM)")
    g.add_argument("--dntp", type=float, help="total dNTP concentration (mM)")
    g.add_argument("--primer-conc", type=float, help="primer concentration (nM)")
    g.add_argument("--dmso", type=float, help="DMSO (%% v/v), lowers Tm")


def _add_specificity_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("specificity")
    g.add_argument(
        "--avoid",
        action="append",
        default=[],
        metavar="FASTA",
        help="sequence the primers must not amplify: the ohnolog, a family member, a "
        "transgene. Repeatable.",
    )
    g.add_argument(
        "--avoid-strategy",
        choices=("auto", "exhaustive", "seeded"),
        default="auto",
        help="scan strategy for the avoid-panel (auto picks by panel size)",
    )
    g.add_argument("--seed-length", type=int, default=12, help="k-mer seed length for seeded scans")
    g.add_argument(
        "--allow-paralog-amplicon",
        action="store_true",
        help="report pairs predicted to co-amplify the panel instead of discarding them",
    )
    g.add_argument(
        "--strict-specificity",
        action="store_true",
        help="reject any primer that can prime on the avoid-panel, even when no product "
        "is predicted. Worth it for multiplex and qPCR; may leave no design.",
    )
    g.add_argument(
        "--priming-3p-window",
        type=int,
        default=5,
        help="3'-terminal window in which mismatches are counted",
    )


def _add_variant_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("strain variants")
    g.add_argument("--variants", help="local VCF or BED of known variants to avoid")
    g.add_argument(
        "--variant-source",
        choices=sorted(VARIANT_SOURCES),
        help="published variant set to query by region instead of a local file "
        "(see `danioprime sources`). Only the locus is transferred, not the "
        "whole dataset.",
    )
    g.add_argument(
        "--variants-url",
        metavar="URL",
        help="any bgzip-compressed, tabix-indexed VCF or BED, local or over https; "
        "its index is assumed to be the same path plus '.tbi'",
    )
    g.add_argument(
        "--assembly",
        help="assembly your template coordinates are on (danRer11/GRCz11, "
        "GRCz12tu, ...). Required with --variant-source: a mismatch is refused, "
        "because GRCz11 and GRCz12 coordinates are not transferable.",
    )
    g.add_argument(
        "--allow-assembly-mismatch",
        action="store_true",
        help="apply a variant set from another assembly anyway. Only correct if "
        "you have already lifted the coordinates over yourself.",
    )
    g.add_argument(
        "--no-index-cache",
        action="store_true",
        help="re-download the tabix index instead of reusing the cached copy",
    )
    g.add_argument("--variant-chrom", help="restrict VCF records to this contig")
    g.add_argument(
        "--strain",
        action="append",
        default=[],
        help="VCF sample name(s) whose variants matter, e.g. --strain TU. Repeatable.",
    )
    g.add_argument(
        "--heterozygous-only",
        action="store_true",
        help="mask only positions heterozygous in the named strain; positions fixed "
        "non-reference in your line are not a design hazard",
    )
    g.add_argument("--variant-pad", type=int, default=0, help="widen each variant by N bases")
    g.add_argument(
        "--reject-variant-overlap",
        action="store_true",
        help="discard any primer overlapping a variant, not merely penalise it",
    )


def _add_output_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("output")
    g.add_argument("--out", "-o", help="write the main report here instead of stdout")
    g.add_argument(
        "--format",
        "-f",
        choices=("text", "tsv", "json"),
        default="text",
        help="format of the main report",
    )
    g.add_argument("--json-out", help="also write the full JSON result here")
    g.add_argument("--tsv-out", help="also write the pair table here")
    g.add_argument("--bed-out", help="also write primer footprints as BED here")
    g.add_argument("--order-out", help="also write an oligo order sheet here")
    g.add_argument(
        "--tails",
        choices=sorted(TAIL_SCHEMES),
        default="none",
        help="5' tail scheme to apply when ordering (see `danioprime tails`)",
    )
    g.add_argument("--show", type=int, default=5, help="pairs to show in the text report")
    g.add_argument("--verbose", "-v", action="store_true", help="include penalty breakdowns")
    g.add_argument(
        "--color",
        "--colour",
        dest="color",
        choices=("auto", "always", "never"),
        default="auto",
        help=(
            "colour 5' tails in the text report; auto means only when stdout is a "
            "terminal (NO_COLOR is honoured). Tails are written in lower case either way"
        ),
    )


# --------------------------------------------------------------------------
# argument resolution
# --------------------------------------------------------------------------
_PARAM_FLAGS = (
    "min_length", "opt_length", "max_length", "min_tm", "opt_tm", "max_tm",
    "max_pair_tm_difference", "min_gc", "max_gc", "min_product_size",
    "opt_product_size", "max_product_size", "require_gc_clamp", "max_poly_x",
    "max_pairs", "thermo_backend",
)


def resolve_params(args: argparse.Namespace) -> DesignParams:
    """Build :class:`DesignParams` from preset, file and explicit flags."""
    if getattr(args, "params", None):
        params = DesignParams.from_json(args.params)
    elif getattr(args, "preset", None):
        params = get_preset(args.preset)
    else:
        params = DesignParams()

    overrides: dict[str, Any] = {}
    for flag in _PARAM_FLAGS:
        value = getattr(args, flag, None)
        if value is not None:
            overrides[flag] = value
    if getattr(args, "reject_variant_overlap", False):
        overrides["reject_variant_overlap"] = True
    if getattr(args, "allow_paralog_amplicon", False):
        overrides["reject_paralog_amplicon"] = False
    if getattr(args, "strict_specificity", False):
        overrides["reject_paralog_priming"] = True
    if getattr(args, "priming_3p_window", None):
        overrides["priming_three_prime_window"] = args.priming_3p_window

    cond_fields = {
        "monovalent_mM": getattr(args, "monovalent", None),
        "divalent_mM": getattr(args, "divalent", None),
        "dntp_mM": getattr(args, "dntp", None),
        "primer_nM": getattr(args, "primer_conc", None),
        "dmso_pct": getattr(args, "dmso", None),
    }
    cond_set = {k: v for k, v in cond_fields.items() if v is not None}
    if cond_set:
        base = params.conditions
        overrides["conditions"] = SaltConditions(
            monovalent_mM=cond_set.get("monovalent_mM", base.monovalent_mM),
            divalent_mM=cond_set.get("divalent_mM", base.divalent_mM),
            dntp_mM=cond_set.get("dntp_mM", base.dntp_mM),
            primer_nM=cond_set.get("primer_nM", base.primer_nM),
            dmso_pct=cond_set.get("dmso_pct", base.dmso_pct),
        )

    params = params.replace(**overrides) if overrides else params
    params.validate()
    if params.thermo_backend == "primer3" and not primer3_available():
        raise DanioPrimeError(
            "--thermo-backend primer3 needs primer3-py: pip install 'danioprime[primer3]'"
        )
    return params


def resolve_template(args: argparse.Namespace) -> tuple[SeqRecord, int]:
    """Load the template and the genomic offset of its first base."""
    offset = getattr(args, "template_offset", 0) or 0
    if getattr(args, "region", None):
        if not getattr(args, "reference", None):
            raise DanioPrimeError("--region needs --reference")
        chrom, interval = parse_interval(args.region)
        if chrom is None:
            raise DanioPrimeError(
                f"--region {args.region!r} needs a contig, e.g. chr5:{args.region}"
            )
        with IndexedFasta(args.reference) as fa:
            record = fa.record(chrom, interval.start, interval.end)
        return record, interval.start
    if not getattr(args, "template", None):
        raise DanioPrimeError("supply --template (FASTA or sequence) or --region with --reference")
    return load_sequence(args.template), offset


def resolve_panel(args: argparse.Namespace) -> ParalogPanel | None:
    """Load the avoid-panel, if any."""
    paths = getattr(args, "avoid", None) or []
    if not paths:
        return None
    records: list[SeqRecord] = []
    for path in paths:
        records.extend(read_fasta(path))
    return ParalogPanel(
        records=records,
        strategy=getattr(args, "avoid_strategy", "auto"),
        seed_length=getattr(args, "seed_length", 12),
    )


def resolve_remote_variants(
    args: argparse.Namespace, offset: int, template_length: int
) -> Any:
    """Query a published variant resource for the template's own interval."""
    source = getattr(args, "variant_source", None)
    url = getattr(args, "variants_url", None)
    if not source and not url:
        return None

    chrom = getattr(args, "variant_chrom", None)
    if not chrom and getattr(args, "region", None):
        parsed, _ = parse_interval(args.region)
        chrom = parsed
    if not chrom:
        raise DanioPrimeError(
            "a remote variant source needs to know which contig your template "
            "came from. Either use --region chr5:START-END with --reference, or "
            "pass --variant-chrom together with --template-offset."
        )
    if offset == 0 and not getattr(args, "region", None):
        sys.stderr.write(
            "note: --template-offset is 0, so the template is assumed to start at "
            f"the first base of {chrom}. If it is a slice, pass its genomic start "
            "or the wrong variants will be fetched.\n"
        )

    return fetch_region(
        chrom,
        offset,
        offset + template_length,
        source=source if source else None,
        data_location=url,
        template_assembly=getattr(args, "assembly", None),
        samples=getattr(args, "strain", None) or None,
        heterozygous_only=getattr(args, "heterozygous_only", False),
        pad=getattr(args, "variant_pad", 0),
        offset=offset,
        allow_assembly_mismatch=getattr(args, "allow_assembly_mismatch", False),
        cache_index=not getattr(args, "no_index_cache", False),
    )


def resolve_variants(args: argparse.Namespace, offset: int) -> tuple[IntervalSet, list[Any]]:
    """Load the variant mask, shifted into template coordinates."""
    path = getattr(args, "variants", None)
    if not path:
        return IntervalSet(), []
    mask, variants = load_mask(
        path,
        chrom=getattr(args, "variant_chrom", None),
        samples=getattr(args, "strain", None) or None,
        heterozygous_only=getattr(args, "heterozygous_only", False),
        offset=offset,
        pad=getattr(args, "variant_pad", 0),
    )
    return mask, variants


def resolve_all_variants(
    args: argparse.Namespace, offset: int, template_length: int
) -> tuple[IntervalSet, list[Any], Any]:
    """Local file and remote source combined into one mask.

    Both may be given: a published strain set plus a lab's own calls is a
    reasonable thing to want.
    """
    mask, local = resolve_variants(args, offset)
    remote = resolve_remote_variants(args, offset, template_length)
    if remote is not None:
        mask = mask.union(remote.mask)
        local = [*local, *remote.variants]
        sys.stderr.write(f"note: {remote.summary()}\n")
        for note in remote.notes:
            sys.stderr.write(f"warning: {note}\n")
    return mask, local, remote


def variant_sections(
    variant_list: list[Any], remote: Any
) -> list[tuple[str, Sequence[str]]]:
    """Report lines describing what was masked and where it came from."""
    sections: list[tuple[str, Sequence[str]]] = []
    if variant_list:
        lines = [f"{kind}: {count}" for kind, count in summarise(variant_list).items()]
        if remote is not None:
            assembly = (
                remote.assembly
                if remote.assembly_verified
                else f"{remote.assembly} (unverified)"
            )
            lines.append(f"source: {remote.source_name} ({assembly})")
            if remote.citation:
                lines.append(f"cite: {remote.citation}")
        sections.append(("Variant mask:", lines))
    return sections


def softmask_for(record: SeqRecord, args: argparse.Namespace) -> IntervalSet:
    if not getattr(args, "respect_softmask", False):
        return IntervalSet()
    from .design.regions import Interval

    runs = record.softmask or soft_mask_runs(record.seq)
    return IntervalSet(Interval(a, b) for a, b in runs)


# --------------------------------------------------------------------------
# shared output handling
# --------------------------------------------------------------------------
def write_outputs(
    result: DesignResult,
    args: argparse.Namespace,
    *,
    offset: int = 0,
    extra: dict[str, Any] | None = None,
    extra_sections: Sequence[tuple[str, Sequence[str]]] = (),
) -> None:
    """Emit the main report plus any additional requested files."""
    plan = None
    if result.pairs and getattr(args, "tails", "none") != "none":
        plan = reaction_plan(result.pairs, scheme=args.tails)
        if extra is None:
            extra = {}
        extra["reaction_plan"] = plan.to_dict()

    fmt = getattr(args, "format", "text")
    if fmt == "text":
        # The report is only coloured when it is going to a terminal, so a
        # report written to a file never carries escape sequences.
        colour = should_colour(
            getattr(args, "color", "auto"),
            None if args.out else sys.stdout,
        )
        body = format_text(
            result,
            verbose=getattr(args, "verbose", False),
            show=getattr(args, "show", 5),
            extra_sections=extra_sections,
            colour=colour,
        )
        if plan:
            lines = [f"scheme: {plan.scheme} - {SCHEME_DESCRIPTIONS.get(plan.scheme, '')}"]
            lines += [
                f"{o.name}: {render_order_sequence(o, colour=colour)} "
                f"({o.length} nt, {o.role})"
                for o in plan.oligos
            ]
            lines += [f"note: {n}" for n in plan.notes]
            lines.append(
                "order from the upper-case column of --order-out; the lower-case "
                "bases above are the tail"
            )
            body += "\n\n" + "Ordering:\n" + "\n".join(f"  {line}" for line in lines)
        emit(body, args.out)
    elif fmt == "tsv":
        emit(write_tsv(result.pairs), args.out)
    else:
        emit(write_json(result, extra=extra), args.out)

    if getattr(args, "json_out", None):
        write_json(result, args.json_out, extra=extra)
    if getattr(args, "tsv_out", None):
        write_tsv(result.pairs, args.tsv_out)
    if getattr(args, "bed_out", None):
        chrom = "template"
        if getattr(args, "region", None):
            parsed_chrom, _ = parse_interval(args.region)
            chrom = parsed_chrom or chrom
        write_bed(result.pairs, args.bed_out, chrom=chrom, offset=offset)
    if getattr(args, "order_out", None):
        plan = plan or reaction_plan(result.pairs, scheme=getattr(args, "tails", "none"))
        write_order_sheet(plan.oligos, args.order_out)


# --------------------------------------------------------------------------
# subcommands
# --------------------------------------------------------------------------
def cmd_design(args: argparse.Namespace) -> int:
    params = resolve_params(args)
    template, offset = resolve_template(args)
    panel = resolve_panel(args)
    variants, variant_list, remote = resolve_all_variants(args, offset, len(template.seq))

    targets = IntervalSet.from_specs(args.target) if args.target else IntervalSet()
    excluded = IntervalSet.from_specs(args.exclude) if args.exclude else IntervalSet()
    fwd = IntervalSet.from_specs(args.forward_region) if args.forward_region else IntervalSet()
    rev = IntervalSet.from_specs(args.reverse_region) if args.reverse_region else IntervalSet()

    request = DesignRequest(
        template=template,
        params=params,
        forward_search=fwd,
        reverse_search=rev,
        targets=targets,
        excluded=excluded,
        variants=variants,
        softmask=softmask_for(template, args),
        panel=panel,
        name=args.name or template.name,
        context={"assay": "generic", "template_offset": offset},
    )
    result = run_design(request)
    write_outputs(
        result, args, offset=offset, extra_sections=variant_sections(variant_list, remote)
    )
    return 0 if result.pairs else 1


def cmd_genotype(args: argparse.Namespace) -> int:
    params = resolve_params(args)
    template, offset = resolve_template(args)
    panel = resolve_panel(args)
    variants, variant_list, remote = resolve_all_variants(args, offset, len(template.seq))

    design = design_genotyping(
        template,
        params,
        targets=IntervalSet.from_specs(args.target) if args.target else IntervalSet(),
        excluded=IntervalSet.from_specs(args.exclude) if args.exclude else IntervalSet(),
        panel=panel,
        variants=variants,
        expected_indels=(
            tuple(args.expected_indel) if args.expected_indel else DEFAULT_EXPECTED_INDELS
        ),
        name=args.name,
    )
    sections: list[tuple[str, Sequence[str]]] = []
    if design.allele_sizes:
        table = summarise_iterable(
            [
                {
                    "allele": r["allele"],
                    "size": r["size"],
                    "gel-resolvable": "yes" if r["resolvable_on_gel"] else "no",
                }
                for r in design.allele_sizes
            ],
            ("allele", "size", "gel-resolvable"),
        )
        sections.append(("Predicted allele sizes (best pair):", table.splitlines()))
    sections.extend(variant_sections(variant_list, remote))
    write_outputs(
        design.result, args, offset=offset, extra=design.to_dict(), extra_sections=sections
    )
    return 0 if design.pairs else 1


def cmd_qpcr(args: argparse.Namespace) -> int:
    params = resolve_params(args)
    if not getattr(args, "preset", None) and not getattr(args, "params", None):
        params = get_preset("qpcr").replace(
            **{
                k: v
                for k, v in ((f, getattr(args, f, None)) for f in _PARAM_FLAGS)
                if v is not None
            }
        )
        params.validate()

    model = None
    junctions: tuple[int, ...] = tuple(args.junction) if args.junction else ()

    if args.gtf:
        transcripts = read_annotation(args.gtf, chrom=args.variant_chrom)
        if args.transcript_id:
            hits = find_transcripts(transcripts, args.transcript_id)
            if not hits:
                raise DanioPrimeError(
                    f"no transcript matching {args.transcript_id!r} in {args.gtf}"
                )
            model = hits[0] if len(hits) == 1 else canonical_transcript(hits)
        else:
            model = canonical_transcript(transcripts.values())

    if args.transcript:
        template = load_sequence(args.transcript)
        offset = 0
    elif model is not None and args.reference:
        with IndexedFasta(args.reference) as fa:
            template = model.spliced_sequence(fa)
        offset = 0
    else:
        template, offset = resolve_template(args)

    panel = resolve_panel(args)
    if panel is None and model is not None and model.gene_name:
        guess = guess_ohnolog_symbol(model.gene_name)
        if guess:
            sys.stderr.write(
                f"note: {model.gene_name} looks like an ohnolog of {guess} under the ZFIN "
                f"naming convention. Pass --avoid {guess}.fa (its transcript sequence) to "
                "check that this assay does not measure both genes. Confirm the pairing "
                "against ZFIN first -- the suffix is a convention, not an orthology call.\n"
            )

    variants, _variant_list, _remote = resolve_all_variants(args, offset, len(template.seq))
    design = design_qpcr(
        template,
        params,
        junctions=junctions,
        transcript_model=model,
        require_junction_spanning=args.require_junction_spanning,
        junction_flank=args.junction_flank,
        panel=panel,
        variants=variants,
        min_protective_intron=args.min_protective_intron,
        name=args.name,
    )
    sections: list[tuple[str, Sequence[str]]] = []
    if design.assessment:
        table = summarise_iterable(
            [
                {"pair": r["pair"], "gDNA protection": r["protection"], "mechanism": r["mechanism"]}
                for r in design.assessment
            ],
            ("pair", "gDNA protection", "mechanism"),
        )
        sections.append(("Genomic-DNA protection:", table.splitlines()))
    write_outputs(design.result, args, extra=design.to_dict(), extra_sections=sections)
    return 0 if design.pairs else 1


def cmd_dcaps(args: argparse.Namespace) -> int:
    params = resolve_params(args)
    template, offset = resolve_template(args)
    enzymes = load_enzymes(args.enzymes) if args.enzymes else DEFAULT_ENZYMES

    assays = design_dcaps(
        template,
        args.variant - 1,
        args.variant_base,
        params,
        enzymes=enzymes,
        cut_allele=args.cut_allele,
        max_introduced_mismatches=args.max_introduced_mismatches,
        min_mismatch_offset=args.min_mismatch_offset,
        max_assays=args.max_assays,
        min_band_separation=args.min_band_separation,
        name=args.name,
    )

    payload = {
        "assay": "dcaps",
        "template": {"name": template.name, "length": len(template.seq)},
        "variant": {
            "position_1based": args.variant,
            "reference_base": template.seq[args.variant - 1],
            "variant_base": args.variant_base.upper(),
            "cut_allele": args.cut_allele,
        },
        "assays": [a.to_dict() for a in assays],
    }

    if args.format == "json":

        emit(json.dumps(payload, indent=2), args.out)
    elif args.format == "tsv":
        rows = [
            {
                "assay": a.name,
                "enzyme": a.dcaps_primer.enzyme.name,
                "site": a.dcaps_primer.enzyme.display,
                "dcaps_primer": a.dcaps_primer.seq,
                "orientation": a.dcaps_primer.orientation,
                "mismatches": a.dcaps_primer.mismatches,
                "partner_primer": a.partner.seq,
                "product_size": a.product_size,
                "ref_fragments": "/".join(str(x) for x in sorted(a.reference_fragments, reverse=True)),
                "var_fragments": "/".join(str(x) for x in sorted(a.variant_fragments, reverse=True)),
                "cut_allele": a.dcaps_primer.cut_allele,
            }
            for a in assays
        ]
        emit(
            summarise_iterable(
                rows,
                (
                    "assay", "enzyme", "site", "dcaps_primer", "orientation", "mismatches",
                    "partner_primer", "product_size", "ref_fragments", "var_fragments",
                    "cut_allele",
                ),
            ),
            args.out,
        )
    else:
        ref = template.seq[args.variant - 1]
        lines = [
            f"danioprime: dCAPS design on {template.name} ({len(template.seq):,} nt)",
            f"  variant: position {args.variant} {ref}>{args.variant_base.upper()}, "
            f"creating a site in the {args.cut_allele} allele",
            f"  enzymes screened: {len(enzymes)}",
            "",
            f"{len(assays)} assay(s):",
        ]
        for a in assays:
            d = a.dcaps_primer
            lines.append("")
            lines.append(
                f"  {a.name}  {d.enzyme.name} ({d.enzyme.display})  "
                f"{d.mismatches} introduced mismatch(es)  product {a.product_size} bp"
            )
            lines.append(f"    dCAPS primer ({d.orientation})  {d.seq}")
            lines.append(f"    template was              {d.template_seq}")
            lines.append(
                f"    mismatch {d.nearest_mismatch_to_3p} nt from the 3' end "
                f"(offsets {list(d.mismatch_offsets_from_3p)})"
            )
            lines.append(f"    partner primer ({a.partner.orientation})  {a.partner.seq}")
            ref_only, var_only = a.diagnostic_bands
            lines.append(
                f"    digest: reference {sorted(a.reference_fragments, reverse=True)} bp  "
                f"variant {sorted(a.variant_fragments, reverse=True)} bp"
            )
            lines.append(
                f"    diagnostic bands: reference-only {sorted(ref_only, reverse=True)}, "
                f"variant-only {sorted(var_only, reverse=True)}"
            )
            for note in a.notes:
                lines.append(f"    - {note}")
        emit("\n".join(lines), args.out)

    if args.json_out:

        Path(args.json_out).write_text(json.dumps(payload, indent=2) + "\n", newline="\n")
    return 0 if assays else 1


def cmd_junction(args: argparse.Namespace) -> int:
    params = resolve_params(args)
    template, offset = resolve_template(args)
    cassette = load_sequence(args.cassette)
    panel = resolve_panel(args)
    variants, _variant_list, _remote = resolve_all_variants(args, offset, len(template.seq))

    designs = design_junction_pcr(
        template,
        args.insert_at - 1,
        cassette,
        params,
        genomic_flank=args.genomic_flank,
        cassette_depth=args.cassette_depth,
        homology_arm=args.homology_arm,
        panel=panel,
        variants=variants,
        name=args.name,
    )
    chunks: list[str] = []
    ok = False
    for which, result in designs.items():
        chunks.append(f"=== {which} junction ===")
        chunks.append(
            format_text(result, verbose=args.verbose, show=args.show)
        )
        ok = ok or bool(result.pairs)
        chunks.append("")
    emit("\n".join(chunks), args.out)
    if args.json_out:

        from .report.writers import result_to_dict

        Path(args.json_out).write_text(
            json.dumps({k: result_to_dict(v) for k, v in designs.items()}, indent=2) + "\n",
            newline="\n",
        )
    if not ok:
        sys.stderr.write(
            "note: a junction PCR needs one primer in genomic sequence and one in the "
            "cassette. If no pair was found, check that --insert-at is right and that "
            "--genomic-flank covers enough sequence outside the homology arm.\n"
        )
    return 0 if ok else 1


def cmd_insitu(args: argparse.Namespace) -> int:
    params = resolve_params(args)
    if not getattr(args, "preset", None) and not getattr(args, "params", None):
        params = get_preset("insitu-probe").replace(
            **{
                k: v
                for k, v in ((f, getattr(args, f, None)) for f in _PARAM_FLAGS)
                if v is not None
            }
        )
        params.validate()
    template, offset = resolve_template(args)
    panel = resolve_panel(args)
    variants, _variant_list, _remote = resolve_all_variants(args, offset, len(template.seq))
    targets = IntervalSet.from_specs(args.target) if args.target else IntervalSet()

    design = design_insitu(
        template,
        params,
        scheme=args.scheme,
        panel=panel,
        variants=variants,
        targets=targets,
        sense_control=not args.no_sense_control,
        min_identity=args.cross_hyb_identity / 100.0,
        min_block=args.cross_hyb_length,
        name=args.name,
    )

    sections: list[tuple[str, Sequence[str]]] = []
    if design.cross_hyb:
        rows = [
            {
                "pair": name,
                "cross-hybridises": "YES" if ch.will_cross_hybridise else "no",
                "longest shared block": (
                    f"{ch.worst.length} nt @ {100 * ch.worst.identity:.0f}% "
                    f"({ch.worst.subject})"
                    if ch.worst
                    else "none"
                ),
            }
            for name, ch in design.cross_hyb.items()
        ]
        sections.append(
            (
                "Probe cross-hybridisation (whole-amplicon, not primer 3' ends):",
                summarise_iterable(
                    rows, ("pair", "cross-hybridises", "longest shared block")
                ).splitlines(),
            )
        )
    if design.sense:
        lines = [f"{o.name}: {o.sequence}" for o in design.sense.oligos]
        sections.append(("Sense control — order these as a second PCR:", lines))

    # The antisense scheme is already attached to the reported pairs.
    args.tails = design.scheme
    write_outputs(design.result, args, offset=offset, extra=design.to_dict(), extra_sections=sections)
    return 0 if design.pairs else 1


def cmd_check(args: argparse.Namespace) -> int:
    params = resolve_params(args)
    panel = resolve_panel(args)
    from .design.candidate import Primer, TemplateContext, evaluate_primer
    from .thermo.backend import get_backend
    from .thermo.structure import hetero_dimer, three_prime_dimer_score

    seqs: list[tuple[str, str]] = []
    for raw in args.primers:
        for piece in raw.split(","):
            piece = piece.strip()
            if not piece:
                continue
            if "=" in piece:
                name, _, seq = piece.partition("=")
            else:
                name, seq = f"primer{len(seqs) + 1}", piece
            candidate = Path(seq)
            if candidate.exists():
                for rec in read_fasta(candidate):
                    seqs.append((rec.name, rec.seq))
            else:
                seqs.append((name.strip(), seq.strip().upper()))
    if not seqs:
        raise DanioPrimeError("no primer sequences supplied to --primers")

    backend = get_backend(params.thermo_backend, params.conditions, params.anneal_temp_c)
    lines = [f"danioprime: evaluating {len(seqs)} primer(s)", ""]
    rows: list[dict[str, Any]] = []
    evaluated: list[tuple[str, str, Primer]] = []
    for name, seq in seqs:
        ctx = TemplateContext(seq=seq, params=params, backend=backend, name=name)
        primer = Primer(seq=seq, start=0, end=len(seq), orientation="F")
        evaluate_primer(primer, ctx)
        evaluated.append((name, seq, primer))
        rows.append(
            {
                "primer": name,
                "len": primer.length,
                "Tm": f"{primer.tm:.1f}",
                "GC%": f"{primer.gc:.0f}",
                "hairpin": f"{primer.hairpin_dg:.1f}",
                "self-dimer": f"{primer.self_dimer_dg:.1f}",
                "3'clamp": primer.gc_clamp,
                "polyX": primer.poly_x,
            }
        )
    lines.append(
        summarise_iterable(
            rows, ("primer", "len", "Tm", "GC%", "hairpin", "self-dimer", "3'clamp", "polyX")
        )
    )

    if len(evaluated) > 1:
        lines.append("")
        lines.append("Pairwise interactions (dG kcal/mol / 3'-end score):")
        pair_rows = []
        for i in range(len(evaluated)):
            for j in range(i + 1, len(evaluated)):
                a_name, a_seq, _ = evaluated[i]
                b_name, b_seq, _ = evaluated[j]
                hd = hetero_dimer(
                    a_seq, b_seq, temp_c=params.anneal_temp_c, conditions=params.conditions
                )
                pair_rows.append(
                    {
                        "pair": f"{a_name} x {b_name}",
                        "dimer dG": f"{hd.dg:.1f}",
                        "3' score": f"{three_prime_dimer_score(a_seq, b_seq):.0f}",
                        "verdict": "watch" if hd.dg < params.max_hetero_dimer_dg else "ok",
                    }
                )
        lines.append(summarise_iterable(pair_rows, ("pair", "dimer dG", "3' score", "verdict")))

    if panel:
        lines.append("")
        lines.append(f"Against the avoid-panel ({len(panel)} sequence(s)):")
        panel_rows = []
        for name, seq, _ in evaluated:
            report = assess_primer(
                seq, panel, primer_name=name, conditions=params.conditions
            )
            panel_rows.append(
                {
                    "primer": name,
                    "closest match": (
                        f"{report.best_site.subject} ({report.best_site.mismatches} mm)"
                        if report.best_site
                        else "none"
                    ),
                    "primes?": "YES" if report.primes_on_panel else "no",
                    "verdict": report.verdict(),
                }
            )
        lines.append(summarise_iterable(panel_rows, ("primer", "closest match", "primes?", "verdict")))

        if len(evaluated) >= 2:
            amplicons, _ = assess_pair(
                evaluated[0][1],
                evaluated[1][1],
                panel,
                min_size=params.min_product_size,
                max_size=max(params.max_product_size, 5000),
                forward_name=evaluated[0][0],
                reverse_name=evaluated[1][0],
                conditions=params.conditions,
            )
            lines.append("")
            if amplicons:
                lines.append("Predicted products on the avoid-panel:")
                lines.append(
                    summarise_iterable(
                        [
                            {
                                "subject": a.subject,
                                "size": a.size,
                                "F mm": a.forward_site.mismatches,
                                "R mm": a.reverse_site.mismatches,
                                "note": a.note,
                            }
                            for a in amplicons
                        ],
                        ("subject", "size", "F mm", "R mm", "note"),
                    )
                )
            else:
                lines.append("No product predicted on the avoid-panel.")

    emit("\n".join(lines), args.out)
    return 0


def cmd_fetch(args: argparse.Namespace) -> int:
    """Pull a real locus or transcript out of Ensembl or UCSC."""
    from .genome.remote import FetchPlan, fetch_locus, fetch_transcript

    plan = FetchPlan() if args.dry_run else None
    kwargs: dict[str, Any] = {
        "species": args.species,
        "transport": None,
        "timeout": args.timeout,
        "plan": plan,
        "base_url": args.base_url,
    }

    if args.transcript:
        locus = fetch_transcript(args.target, kind=args.transcript_type, **kwargs)
    else:
        chrom = start = end = None
        parsed_chrom, interval = (None, None)
        try:
            parsed_chrom, interval = parse_interval(args.target)
        except RegionError:
            pass
        if interval is not None:
            if parsed_chrom is None:
                raise RegionError(
                    f"{args.target!r} needs a contig to be fetched as a region, "
                    f"e.g. chr18:{args.target}"
                )
            # parse_interval hands back half-open zero-based; Ensembl and this
            # module both speak 1-based inclusive.
            chrom, start, end = parsed_chrom, interval.start + 1, interval.end
        locus = fetch_locus(
            args.target,
            source=args.source,
            genome=args.assembly,
            flank=args.flank,
            chrom=chrom,
            start=start,
            end=end,
            **kwargs,
        )

    if plan is not None:
        lines = ["Would request:", *(f"  {url}" for url in plan.urls)]
        if locus.chrom and not locus.is_transcript:
            lines += ["", f"Region: {locus.region} on {locus.assembly or 'the server default'}"]
        lines.append("")
        lines.append("Nothing was fetched and nothing was written (--dry-run).")
        # Deliberately stdout, not --out: --out names where the FASTA goes, and
        # a dry run that overwrote it with a list of URLs would be a poor way
        # to find out what --dry-run means.
        emit("\n".join(lines), None)
        return 0

    record = locus.to_record()
    if args.format == "json":
        emit(json.dumps(locus.to_dict(), indent=2), args.out)
        return 0

    if args.out:
        write_fasta([record], args.out)
        summary = [
            f"wrote {args.out}: {record.name}, {len(locus.seq):,} bp",
            f"  {locus.header()}",
        ]
        if not locus.is_transcript:
            summary += [
                "",
                "Pass this to a design command with the assembly it is on, so that a",
                "variant set on another assembly is refused rather than misapplied:",
                "",
                f"  danioprime genotype --template {args.out} \\",
                f"      --assembly {locus.assembly} --variant-source nhgri1",
            ]
        print("\n".join(summary))
        return 0

    wrapped = "\n".join(record.seq[i : i + 60] for i in range(0, len(record.seq), 60))
    emit(f">{record.name} {record.description}\n{wrapped}", None)
    return 0


def cmd_sources(args: argparse.Namespace) -> int:
    if getattr(args, "clear_cache", False):
        removed = clear_index_cache()
        emit(f"removed {removed} cached index file(s)", args.out)
        return 0

    rows = describe_sources()
    lines: list[str] = ["Published variant resources that can be queried by region:", ""]
    for row in rows:
        lines.append(f"  {row['key']}  —  {row['name']}  [{row['assembly']}, {row['kind']}]")
        lines.append(f"      {row['description']}")
        lines.append(f"      data:  {row['data_url']}")
        lines.append(f"      page:  {row['homepage']}")
        lines.append(f"      hub:   {row['track_hub']}  (attach in the UCSC browser)")
        lines.append(f"      cite:  {row['citation']}")
        lines.append("")
    lines += [
        "Only the queried locus is transferred, via the tabix index: about 66 kB",
        "per locus regardless of how large the dataset is. The index is cached.",
        "",
        "Every source here is on danRer11 (GRCz11). ZFIN now uses GRCz12tu, whose",
        "coordinates are NOT transferable. Declare your template's assembly with",
        "--assembly; a mismatch is refused rather than silently masking the wrong",
        "bases. Example:",
        "",
        "  danioprime genotype --region chr5:1,000,000-1,004,000 \\",
        "      --reference danRer11.fa --assembly danRer11 \\",
        "      --variant-source nhgri1 --strain NHGRI1 --heterozygous-only",
        "",
        "--variants-url takes any bgzip+tabix VCF or BED, so a source missing from",
        "this list is never a blocker.",
    ]
    emit("\n".join(lines), args.out)
    return 0


def cmd_presets(args: argparse.Namespace) -> int:
    rows = []
    for name, desc in sorted(PRESET_DESCRIPTIONS.items()):
        p = get_preset(name)
        rows.append(
            {
                "preset": name,
                "product": f"{p.min_product_size}-{p.max_product_size}",
                "Tm": f"{p.min_tm:.0f}-{p.max_tm:.0f}",
                "description": desc,
            }
        )
    emit(summarise_iterable(rows, ("preset", "product", "Tm", "description")), args.out)
    return 0


def cmd_enzymes(args: argparse.Namespace) -> int:
    table = load_enzymes(args.enzymes) if args.enzymes else DEFAULT_ENZYMES
    rows = [
        {"enzyme": e.name, "site": e.site, "cut": e.display, "degenerate": "yes" if e.is_degenerate else ""}
        for e in sorted(table.values(), key=lambda x: x.name)
    ]
    emit(summarise_iterable(rows, ("enzyme", "site", "cut", "degenerate")), args.out)
    return 0


def cmd_tails(args: argparse.Namespace) -> int:
    from .assays.tails import TAILS, tail_style
    from .report.colour import paint

    rows = [
        {"scheme": name, "forward": mapping.get("forward", "-"), "reverse": mapping.get("reverse", "-"),
         "description": SCHEME_DESCRIPTIONS.get(name, "")}
        for name, mapping in sorted(TAIL_SCHEMES.items())
    ]
    text = ["Schemes:", summarise_iterable(rows, ("scheme", "forward", "reverse", "description")), ""]
    text.append("Tails:")
    tails = sorted(TAILS.values(), key=lambda x: x.name)
    text.append(
        summarise_iterable(
            # Lower case because that is how a tail is written everywhere else
            # in the output; the table is the legend for the reports.
            [
                {"tail": t.name, "sequence": t.seq.lower(), "length": len(t), "purpose": t.purpose}
                for t in tails
            ],
            ("tail", "sequence", "length", "purpose"),
        )
    )
    # Painted after the table is laid out, so no escape sequence ever lands
    # inside a column whose width was measured without it.
    colour = should_colour(getattr(args, "color", "auto"), None if args.out else sys.stdout)
    text.append("")
    swatches = " ".join(paint(t.name, tail_style(t.name).ansi, colour=colour) for t in tails)
    text.append(f"Colours in the text report: {swatches}")
    if not colour:
        text.append(
            "(not shown: stdout is not a terminal, or colour is off. The lower case "
            "carries the same information.)"
        )
    emit("\n".join(text), args.out)
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    aligners = available_aligners()
    lines = [
        f"danioprime {__version__}",
        f"python {sys.version.split()[0]}",
        "",
        "Thermodynamic backends:",
        "  builtin  available (pure Python, no dependencies, reproducible)",
        f"  primer3  {'available' if primer3_available() else 'not installed (pip install danioprime[primer3])'}",
        "",
        "Optional external aligners for genome-scale specificity:",
    ]
    for name, path in aligners.items():
        lines.append(f"  {name:<14} {path or 'not found on PATH'}")
    lines += [
        "",
        f"Built-in enzymes: {len(DEFAULT_ENZYMES)}",
        f"Presets: {', '.join(sorted(PRESET_DESCRIPTIONS))}",
        f"Tail schemes: {', '.join(sorted(TAIL_SCHEMES))}",
    ]
    emit("\n".join(lines), args.out)
    return 0


# --------------------------------------------------------------------------
# parser construction
# --------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="danioprime",
        description=(
            "Zebrafish-aware PCR primer design: ohnolog discrimination, strain-variant "
            "masking and assay-specific geometry."
        ),
        epilog="Run `danioprime <subcommand> --help` for the options of each mode.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"danioprime {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<subcommand>")

    # --- design -----------------------------------------------------------
    p = sub.add_parser("design", help="general-purpose pair design on any template")
    _add_template_args(p)
    p.add_argument("--target", action="append", default=[], help="amplicon must span this region (repeatable)")
    p.add_argument("--exclude", action="append", default=[], help="primers must avoid this region (repeatable)")
    p.add_argument("--forward-region", action="append", default=[], help="restrict forward primers here")
    p.add_argument("--reverse-region", action="append", default=[], help="restrict reverse primers here")
    p.add_argument("--name", help="prefix for pair names")
    _add_param_args(p)
    _add_specificity_args(p)
    _add_variant_args(p)
    _add_output_args(p)
    p.set_defaults(func=cmd_design)

    # --- genotype ---------------------------------------------------------
    p = sub.add_parser("genotype", help="genotype an edited locus (CRISPR indel or point mutation)")
    _add_template_args(p)
    g = p.add_argument_group("region to interrogate")
    g.add_argument(
        "--target", action="append", default=[],
        help="the amplicon must span this region and no primer may sit in it, e.g. the exon "
             "you disrupted (repeatable). Omit it to design over the whole template, which is "
             "the right answer when the extent of the lesion is unknown",
    )
    g.add_argument("--exclude", action="append", default=[],
                   help="primers must avoid this region (repeatable)")
    g.add_argument("--expected-indel", type=int, action="append", default=[],
                   help="indel size to predict product sizes for (repeatable)")
    p.add_argument("--name", help="prefix for pair names")
    _add_param_args(p)
    _add_specificity_args(p)
    _add_variant_args(p)
    _add_output_args(p)
    p.set_defaults(func=cmd_genotype)

    # --- qpcr -------------------------------------------------------------
    p = sub.add_parser("qpcr", help="RT-qPCR pair, optionally spanning an exon junction")
    _add_template_args(p)
    g = p.add_argument_group("transcript")
    g.add_argument("--transcript", help="FASTA of the spliced transcript")
    g.add_argument("--gtf", help="GTF/GFF annotation to take exon structure from")
    g.add_argument("--transcript-id", help="transcript ID, gene ID or gene symbol to design against")
    g.add_argument("--junction", type=int, action="append", default=[],
                   help="exon-exon junction offset in the transcript (repeatable)")
    g.add_argument("--require-junction-spanning", action="store_true",
                   help="require that one primer straddles a junction (the strongest "
                        "protection against genomic-DNA amplification)")
    g.add_argument("--junction-flank", type=int, default=60,
                   help="how far either side of a junction to search")
    g.add_argument("--min-protective-intron", type=int, default=1000,
                   help="intron size below which intron-spanning is not counted as protection")
    p.add_argument("--name", help="prefix for pair names")
    _add_param_args(p)
    _add_specificity_args(p)
    _add_variant_args(p)
    _add_output_args(p)
    p.set_defaults(func=cmd_qpcr)

    # --- dcaps ------------------------------------------------------------
    p = sub.add_parser("dcaps", help="dCAPS assay for a variant too small to size on a gel")
    _add_template_args(p)
    g = p.add_argument_group("variant")
    g.add_argument("--variant", type=int, required=True, help="1-based position of the variant")
    g.add_argument("--variant-base", required=True, help="the alternate base (A/C/G/T)")
    g.add_argument("--cut-allele", choices=("variant", "reference"), default="variant",
                   help="which allele the engineered site should cut")
    g.add_argument("--enzymes", help="TSV of custom enzymes (name, site, cut offset)")
    g.add_argument("--max-introduced-mismatches", type=int, default=2,
                   help="mismatches the dCAPS primer may introduce")
    g.add_argument("--min-mismatch-offset", type=int, default=1,
                   help="keep mismatches at least this far from the 3' end")
    g.add_argument("--max-assays", type=int, default=10, help="assays to report")
    g.add_argument("--min-band-separation", type=int, default=25,
                   help="warn below this diagnostic band size")
    p.add_argument("--name", help="prefix for assay names")
    _add_param_args(p)
    _add_variant_args(p)
    _add_output_args(p)
    p.set_defaults(func=cmd_dcaps)

    # --- junction ---------------------------------------------------------
    p = sub.add_parser("junction", help="junction PCRs verifying a targeted knock-in")
    _add_template_args(p)
    g = p.add_argument_group("knock-in")
    g.add_argument("--cassette", required=True, help="FASTA of the inserted cassette")
    g.add_argument("--insert-at", type=int, required=True,
                   help="1-based genomic position of the insertion point")
    g.add_argument("--genomic-flank", type=int, default=500, help="genomic sequence each side to use")
    g.add_argument("--cassette-depth", type=int, default=400, help="cassette sequence each end to use")
    g.add_argument("--homology-arm", type=int, default=0,
                   help="homology arm length: primers are kept outside it, so the assay "
                        "cannot read the donor plasmid instead of the integrated allele")
    p.add_argument("--name", help="prefix for pair names")
    _add_param_args(p)
    _add_specificity_args(p)
    _add_variant_args(p)
    _add_output_args(p)
    p.set_defaults(func=cmd_junction)

    # --- insitu -----------------------------------------------------------
    p = sub.add_parser(
        "insitu", help="riboprobe template for in situ hybridisation (T7/T3/SP6)"
    )
    _add_template_args(p)
    g = p.add_argument_group("probe")
    g.add_argument(
        "--scheme",
        choices=sorted(RIBOPROBE_SCHEMES),
        default="riboprobe-t7",
        help="promoter arrangement. A promoter on the reverse primer transcribes the "
        "antisense probe; on the forward primer, the sense control. A dual scheme "
        "puts different promoters at each end so one product makes both.",
    )
    g.add_argument(
        "--no-sense-control",
        action="store_true",
        help="do not emit a second primer set for the sense control",
    )
    g.add_argument(
        "--target", action="append", default=[],
        help="confine the probe to this region (repeatable); useful for steering it "
        "into a divergent stretch such as the 3' UTR",
    )
    g.add_argument(
        "--cross-hyb-identity", type=float, default=85.0,
        help="percent identity at which a shared block is treated as able to "
        "cross-hybridise (default 85; calibrated so a diverged paralog does not "
        "trip it by chance)",
    )
    g.add_argument(
        "--cross-hyb-length", type=int, default=100,
        help="shortest shared block considered a cross-hybridisation risk (default 100)",
    )
    p.add_argument("--name", help="prefix for pair names")
    _add_param_args(p)
    _add_specificity_args(p)
    _add_variant_args(p)
    _add_output_args(p)
    p.set_defaults(func=cmd_insitu)

    # --- check ------------------------------------------------------------
    p = sub.add_parser("check", help="evaluate primers you already have")
    p.add_argument("--primers", action="append", required=True,
                   help="comma-separated sequences, name=sequence pairs, or a FASTA path")
    _add_param_args(p)
    _add_specificity_args(p)
    p.add_argument("--out", "-o", help="write the report here instead of stdout")
    p.set_defaults(func=cmd_check)

    # --- fetching real sequence -------------------------------------------
    p = sub.add_parser(
        "fetch",
        help="fetch a real locus or transcript by gene symbol or region",
        description=(
            "Fetch sequence from Ensembl or UCSC so a design runs against the real "
            "locus. The FASTA header records the assembly the server reported, the "
            "coordinates, the source URL and the date, which is what lets a variant "
            "set on another assembly be refused rather than misapplied."
        ),
    )
    p.add_argument(
        "target",
        help="gene symbol (myo7aa), a region (chr18:1,000,000-1,010,000), or an Ensembl ID",
    )
    p.add_argument("--species", default="danio_rerio", help="Ensembl species name")
    p.add_argument(
        "--source", choices=PROVIDERS, default="ensembl",
        help="ensembl resolves gene symbols; ucsc serves a named assembly by coordinates",
    )
    p.add_argument(
        "--assembly", default="danRer11",
        help="UCSC genome to fetch from (--source ucsc); ignored for Ensembl, which "
             "reports its own",
    )
    p.add_argument(
        "--flank", type=int, default=2000,
        help="bases to add on each side of a gene, so primers have somewhere to sit",
    )
    p.add_argument(
        "--transcript", action="store_true",
        help="fetch the canonical spliced transcript instead of the genomic locus "
             "(for qPCR and in situ probes)",
    )
    p.add_argument(
        "--transcript-type", default="cdna", choices=("cdna", "cds", "protein"),
        help="what --transcript returns",
    )
    p.add_argument("--timeout", type=float, default=30.0)
    p.add_argument(
        "--base-url", default="",
        help="override the service root, e.g. a local Ensembl REST mirror",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="print the URLs that would be requested, fetch nothing, write nothing",
    )
    p.add_argument("--format", choices=("fasta", "json"), default="fasta")
    p.add_argument("--out", "-o", help="write the FASTA here instead of stdout")
    p.set_defaults(func=cmd_fetch)

    # --- informational ----------------------------------------------------
    p = sub.add_parser(
        "sources", help="list published variant resources that can be queried by region"
    )
    p.add_argument(
        "--clear-cache", action="store_true", help="delete cached tabix indexes"
    )
    p.add_argument("--out", "-o")
    p.set_defaults(func=cmd_sources)

    p = sub.add_parser("presets", help="list design presets")
    p.add_argument("--out", "-o")
    p.set_defaults(func=cmd_presets)

    p = sub.add_parser("enzymes", help="list restriction enzymes available for dCAPS")
    p.add_argument("--enzymes", help="TSV of custom enzymes to list instead")
    p.add_argument("--out", "-o")
    p.set_defaults(func=cmd_enzymes)

    p = sub.add_parser("tails", help="list 5' tail schemes, their sequences and their colours")
    p.add_argument("--out", "-o")
    p.add_argument(
        "--color", "--colour", dest="color", choices=("auto", "always", "never"), default="auto"
    )
    p.set_defaults(func=cmd_tails)

    p = sub.add_parser("info", help="show backends, optional tools and built-in tables")
    p.add_argument("--out", "-o")
    p.set_defaults(func=cmd_info)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 2
    try:
        return int(args.func(args))
    except DanioPrimeError as exc:
        sys.stderr.write(f"danioprime: {exc}\n")
        return 2
    except BrokenPipeError:  # pragma: no cover - piping into head
        return 0
    except KeyboardInterrupt:  # pragma: no cover
        sys.stderr.write("interrupted\n")
        return 130


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
