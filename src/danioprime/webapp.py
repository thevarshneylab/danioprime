"""A single entry point for the browser build.

The web UI runs this package unmodified under Pyodide (CPython on WebAssembly),
which works because `danioprime` is pure Python with no required dependencies:
the wheel is a zip, so `sys.path.insert(0, "danioprime.whl")` is the whole
install step -- no package manager, no compilation, no server.

All of the glue between a form's worth of strings and the design API lives here
rather than in JavaScript, for one reason: this is the layer most likely to
contain bugs, and here it is reachable by the normal test suite. `index.html`
only has to hand :func:`run` a dict and render what comes back.

Everything is one function with a JSON-shaped input and a JSON-shaped output,
so the same call works from the browser, a notebook, or a test.
"""

from __future__ import annotations

import json
import traceback
from typing import Any

from .annotation.gtf import canonical_transcript, find_transcripts, parse_annotation
from .assays.dcaps import design_dcaps
from .assays.genotyping import (
    NUCLEASES,
    CutSite,
    GenotypingGeometry,
    design_genotyping,
    locate_protospacer,
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
from .design.regions import IntervalSet
from .errors import DanioPrimeError
from .report.writers import format_text, result_to_dict, write_order_sheet, write_tsv
from .sequence import SeqRecord, normalise, parse_fasta
from .specificity.paralog import ParalogPanel, assess_pair, assess_primer
from .thermo.nn import SaltConditions

MODES = ("genotype", "qpcr", "dcaps", "insitu", "design", "check")

#: Fields of :class:`DesignParams` the UI is allowed to override directly.
TUNABLE = (
    "min_length", "opt_length", "max_length",
    "min_tm", "opt_tm", "max_tm", "max_pair_tm_difference",
    "min_gc", "max_gc",
    "min_product_size", "opt_product_size", "max_product_size",
    "require_gc_clamp", "max_poly_x", "max_pairs",
)

_INT_FIELDS = frozenset(
    {
        "min_length", "opt_length", "max_length", "min_product_size",
        "opt_product_size", "max_product_size", "require_gc_clamp",
        "max_poly_x", "max_pairs",
    }
)


def parse_sequences(text: str, *, default_name: str = "input") -> list[SeqRecord]:
    """Read pasted text as FASTA if it has a header, else as a bare sequence.

    Users paste both, and silently mangling one of them would be worse than
    either guessing or asking: a bare sequence with a stray '>' line is FASTA,
    and anything else is one sequence.
    """
    text = (text or "").strip()
    if not text:
        return []
    if text.lstrip().startswith(">"):
        records = list(parse_fasta(text.splitlines()))
        if not records:
            raise DanioPrimeError("the text starts with '>' but holds no FASTA records")
        # A header with no bases under it parses happily into a zero-length
        # record, which then fails much later with a confusing message. Catch it
        # here, where the name of the offending record is still to hand.
        empty = [r.name or "(unnamed)" for r in records if not r.seq]
        if empty:
            raise DanioPrimeError(
                "these FASTA records have a header but no sequence: " + ", ".join(empty)
            )
        return records
    return [SeqRecord(name=default_name, seq=normalise(text))]


def parse_one_sequence(text: str, *, what: str, default_name: str = "input") -> SeqRecord:
    records = parse_sequences(text, default_name=default_name)
    if not records:
        raise DanioPrimeError(f"no {what} supplied")
    if len(records) > 1:
        raise DanioPrimeError(
            f"{what} holds {len(records)} sequences; paste only one "
            f"(the first is '{records[0].name}')"
        )
    return records[0]


def build_params(payload: dict[str, Any]) -> DesignParams:
    """Preset plus whatever overrides the form supplied."""
    preset = payload.get("preset") or ""
    params = get_preset(preset) if preset and preset != "default" else DesignParams()

    overrides: dict[str, Any] = {}
    for key, value in (payload.get("overrides") or {}).items():
        if key not in TUNABLE:
            raise DanioPrimeError(f"parameter {key!r} cannot be set from the web form")
        if value in (None, ""):
            continue
        overrides[key] = int(value) if key in _INT_FIELDS else float(value)

    if payload.get("strict_specificity"):
        overrides["reject_paralog_priming"] = True
    if payload.get("allow_paralog_amplicon"):
        overrides["reject_paralog_amplicon"] = False
    if payload.get("reject_variant_overlap"):
        overrides["reject_variant_overlap"] = True

    cond = payload.get("conditions") or {}
    if cond:
        base = params.conditions
        overrides["conditions"] = SaltConditions(
            monovalent_mM=float(cond.get("monovalent", base.monovalent_mM)),
            divalent_mM=float(cond.get("divalent", base.divalent_mM)),
            dntp_mM=float(cond.get("dntp", base.dntp_mM)),
            primer_nM=float(cond.get("primer", base.primer_nM)),
            dmso_pct=float(cond.get("dmso", base.dmso_pct)),
        )

    params = params.replace(**overrides) if overrides else params
    params.validate()
    return params


def build_panel(payload: dict[str, Any]) -> ParalogPanel | None:
    records = parse_sequences(payload.get("avoid") or "", default_name="avoid")
    if not records:
        return None
    for i, rec in enumerate(records, 1):
        if rec.name == "avoid" and len(records) == 1:
            rec.name = "avoid-panel"
        elif not rec.name:
            rec.name = f"avoid{i}"
    return ParalogPanel(records=records)


def _pair_rows(result: DesignResult) -> list[dict[str, Any]]:
    """Flat, display-ready rows; the UI should not have to walk nested dicts."""
    rows: list[dict[str, Any]] = []
    for pair in result.pairs:
        rows.append(
            {
                "name": pair.name,
                "forward": pair.forward.ordered_seq,
                "reverse": pair.reverse.ordered_seq,
                "forward_len": len(pair.forward.ordered_seq),
                "reverse_len": len(pair.reverse.ordered_seq),
                "product_size": pair.product_size,
                "forward_tm": round(pair.forward.tm, 1),
                "reverse_tm": round(pair.reverse.tm, 1),
                "tm_diff": round(pair.tm_difference, 1),
                "forward_gc": round(pair.forward.gc),
                "reverse_gc": round(pair.reverse.gc),
                "forward_pos": pair.forward.footprint.as_one_based(),
                "reverse_pos": pair.reverse.footprint.as_one_based(),
                "hetero_dimer_dg": round(pair.hetero_dimer_dg, 1),
                "specific": pair.is_specific,
                "penalty": round(pair.penalty, 2),
                "penalty_terms": {k: round(v, 2) for k, v in pair.penalty_terms.items()},
                "notes": list(pair.notes)
                + [f"forward: {n}" for n in pair.forward.notes]
                + [f"reverse: {n}" for n in pair.reverse.notes],
                "spans_junction": (
                    pair.forward.spans_junction
                    if pair.forward.spans_junction is not None
                    else pair.reverse.spans_junction
                ),
            }
        )
    return rows


def _package_result(
    result: DesignResult,
    *,
    tails: str = "none",
    extra: dict[str, Any] | None = None,
    tables: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    plan = None
    if result.pairs and tails and tails != "none":
        plan = reaction_plan(result.pairs, scheme=tails)
    payload: dict[str, Any] = {
        "ok": True,
        "mode": result.request.context.get("assay", "design"),
        "report": format_text(result, verbose=True, show=len(result.pairs) or 1),
        "pairs": _pair_rows(result),
        "warnings": list(result.warnings),
        "stats": {
            "template": result.template.name,
            "template_length": len(result.template.seq),
            "forward_candidates": result.stats.get("forward_candidates", 0),
            "reverse_candidates": result.stats.get("reverse_candidates", 0),
        },
        "diagnosis": result.diagnosis(),
        "tsv": write_tsv(result.pairs),
        "json": json.dumps(result_to_dict(result, extra=extra), indent=2),
        "tables": tables or [],
    }
    if plan:
        payload["order_sheet"] = write_order_sheet(plan.oligos)
        payload["order_notes"] = list(plan.notes)
        payload["order_scheme"] = f"{plan.scheme} - {SCHEME_DESCRIPTIONS.get(plan.scheme, '')}"
    if extra:
        payload["extra"] = extra
    return payload


def _mode_genotype(payload: dict[str, Any]) -> dict[str, Any]:
    template = parse_one_sequence(payload.get("template"), what="template", default_name="locus")
    params = build_params(payload)
    panel = build_panel(payload)

    protospacer = (payload.get("protospacer") or "").strip()
    cut_site_value = payload.get("cut_site")
    nuclease = payload.get("nuclease") or "SpCas9"

    if protospacer:
        sites = locate_protospacer(
            template.seq, protospacer, nuclease=nuclease, require_pam=not payload.get("no_pam")
        )
        if not sites:
            pam = NUCLEASES[nuclease]["pam"]
            raise DanioPrimeError(
                f"that protospacer is not in the template on either strand with a valid "
                f"{pam} PAM. Check the sequence (paste it without the PAM), or tick "
                f"'locate without PAM', or give a cut position instead."
            )
        if len(sites) > 1:
            where = ", ".join(str(s.position + 1) for s in sites)
            raise DanioPrimeError(
                f"that protospacer occurs {len(sites)} times in the template (cut sites "
                f"{where}); give a cut position instead to say which one you mean"
            )
        site = sites[0]
    elif cut_site_value not in (None, ""):
        site = CutSite(position=int(cut_site_value) - 1, strand="+")
    else:
        raise DanioPrimeError("supply a protospacer or a cut position")

    geometry = GenotypingGeometry(
        protected=int(payload.get("protected") or 40),
        min_flank=int(payload.get("min_flank") or 40),
        max_flank=int(payload.get("max_flank") or 400),
    )
    design = design_genotyping(
        template, site, params, geometry=geometry, panel=panel, name=payload.get("name") or None
    )
    tables = []
    if design.allele_sizes:
        tables.append(
            {
                "title": "Predicted allele sizes (best pair)",
                "columns": ["allele", "size (bp)", "resolvable on a gel"],
                "rows": [
                    [r["allele"], r["size"], "yes" if r["resolvable_on_gel"] else "no"]
                    for r in design.allele_sizes
                ],
            }
        )
    return _package_result(
        design.result,
        tails=payload.get("tails") or "none",
        extra=design.to_dict(),
        tables=tables,
    )


def _mode_qpcr(payload: dict[str, Any]) -> dict[str, Any]:
    transcript = parse_one_sequence(
        payload.get("template"), what="transcript", default_name="transcript"
    )
    params = build_params(payload)
    panel = build_panel(payload)

    model = None
    gtf_text = (payload.get("gtf") or "").strip()
    if gtf_text:
        transcripts = parse_annotation(gtf_text.splitlines())
        if not transcripts:
            raise DanioPrimeError("no exon features with a transcript identifier in that GTF")
        wanted = (payload.get("transcript_id") or "").strip()
        if wanted:
            hits = find_transcripts(transcripts, wanted)
            if not hits:
                raise DanioPrimeError(
                    f"no transcript matching {wanted!r}; the GTF has: "
                    + ", ".join(sorted(transcripts)[:8])
                )
            model = hits[0] if len(hits) == 1 else canonical_transcript(hits)
        else:
            model = canonical_transcript(transcripts.values())

    junctions = tuple(
        int(j) for j in (payload.get("junctions") or []) if str(j).strip() not in ("", "None")
    )
    design = design_qpcr(
        transcript,
        params,
        junctions=junctions,
        transcript_model=model,
        require_junction_spanning=bool(payload.get("require_junction_spanning")),
        junction_flank=int(payload.get("junction_flank") or 60),
        panel=panel,
        name=payload.get("name") or None,
    )
    tables = []
    if design.assessment:
        tables.append(
            {
                "title": "Genomic-DNA protection",
                "columns": ["pair", "protection", "mechanism"],
                "rows": [[r["pair"], r["protection"], r["mechanism"]] for r in design.assessment],
            }
        )
    return _package_result(
        design.result,
        tails=payload.get("tails") or "none",
        extra=design.to_dict(),
        tables=tables,
    )


def _mode_dcaps(payload: dict[str, Any]) -> dict[str, Any]:
    template = parse_one_sequence(payload.get("template"), what="template", default_name="locus")
    params = build_params(payload)
    position = payload.get("variant_position")
    base = (payload.get("variant_base") or "").strip().upper()
    if position in (None, ""):
        raise DanioPrimeError("supply the variant position (1-based)")
    if not base:
        raise DanioPrimeError("supply the variant base (A, C, G or T)")

    assays = design_dcaps(
        template,
        int(position) - 1,
        base,
        params,
        cut_allele=payload.get("cut_allele") or "variant",
        max_introduced_mismatches=int(payload.get("max_introduced_mismatches") or 2),
        min_mismatch_offset=int(payload.get("min_mismatch_offset") or 1),
        max_assays=int(payload.get("max_assays") or 8),
        name=payload.get("name") or None,
    )

    ref = template.seq[int(position) - 1]
    lines = [
        f"dCAPS design on {template.name} ({len(template.seq):,} nt)",
        f"  variant: position {position} {ref}>{base}, creating a site in the "
        f"{payload.get('cut_allele') or 'variant'} allele",
        "",
        f"{len(assays)} assay(s):",
    ]
    rows = []
    for a in assays:
        d = a.dcaps_primer
        ref_only, var_only = a.diagnostic_bands
        lines += [
            "",
            f"  {a.name}  {d.enzyme.name} ({d.enzyme.display})  "
            f"{d.mismatches} introduced mismatch(es)  product {a.product_size} bp",
            f"    dCAPS primer ({d.orientation})  {d.seq}",
            f"    template was              {d.template_seq}",
            f"    mismatch {d.nearest_mismatch_to_3p} nt from the 3' end",
            f"    partner primer ({a.partner.orientation})  {a.partner.seq}",
            f"    digest: reference {sorted(a.reference_fragments, reverse=True)} bp  "
            f"variant {sorted(a.variant_fragments, reverse=True)} bp",
            f"    diagnostic bands: reference-only {sorted(ref_only, reverse=True)}, "
            f"variant-only {sorted(var_only, reverse=True)}",
        ]
        lines += [f"    - {n}" for n in a.notes]
        rows.append(
            [
                a.name,
                d.enzyme.name,
                d.enzyme.display,
                d.orientation,
                d.seq,
                d.mismatches,
                a.partner.seq,
                a.product_size,
                "/".join(str(x) for x in sorted(a.reference_fragments, reverse=True)),
                "/".join(str(x) for x in sorted(a.variant_fragments, reverse=True)),
            ]
        )

    payload_out: dict[str, Any] = {
        "ok": True,
        "mode": "dcaps",
        "report": "\n".join(lines),
        "pairs": [],
        "warnings": [],
        "stats": {"template": template.name, "template_length": len(template.seq)},
        "diagnosis": "",
        "tables": [
            {
                "title": "dCAPS assays",
                "columns": [
                    "assay", "enzyme", "site", "orient", "dCAPS primer", "mismatches",
                    "partner primer", "product", "reference bands", "variant bands",
                ],
                "rows": rows,
            }
        ],
        "json": json.dumps(
            {
                "assay": "dcaps",
                "template": {"name": template.name, "length": len(template.seq)},
                "variant": {
                    "position_1based": int(position),
                    "reference_base": ref,
                    "variant_base": base,
                },
                "assays": [a.to_dict() for a in assays],
            },
            indent=2,
        ),
    }
    header = [
        "assay", "enzyme", "site", "orientation", "dcaps_primer", "mismatches",
        "partner_primer", "product_size", "reference_fragments", "variant_fragments",
    ]
    payload_out["tsv"] = "\n".join(
        ["\t".join(header)] + ["\t".join(str(c) for c in r) for r in rows]
    )
    return payload_out


def _mode_insitu(payload: dict[str, Any]) -> dict[str, Any]:
    template = parse_one_sequence(
        payload.get("template"), what="transcript", default_name="transcript"
    )
    params = build_params(payload)
    panel = build_panel(payload)

    def specs(key: str) -> IntervalSet:
        raw = payload.get(key) or []
        if isinstance(raw, str):
            raw = [piece.strip() for piece in raw.replace(",", " ").split() if piece.strip()]
        return IntervalSet.from_specs(raw) if raw else IntervalSet()

    scheme = payload.get("scheme") or "riboprobe-t7"
    design = design_insitu(
        template,
        params,
        scheme=scheme,
        panel=panel,
        targets=specs("targets"),
        sense_control=not payload.get("no_sense_control"),
        min_identity=float(payload.get("cross_hyb_identity") or 85.0) / 100.0,
        min_block=int(payload.get("cross_hyb_length") or 100),
        name=payload.get("name") or None,
    )

    out = _package_result(design.result, tails=scheme, extra=design.to_dict())
    if design.cross_hyb:
        out["tables"].append(
            {
                "title": "Probe cross-hybridisation (whole amplicon, not primer 3' ends)",
                "columns": ["pair", "cross-hybridises", "longest shared block", "verdict"],
                "rows": [
                    [
                        name,
                        "YES" if ch.will_cross_hybridise else "no",
                        (
                            f"{ch.worst.length} nt @ {100 * ch.worst.identity:.0f}% "
                            f"({ch.worst.subject})"
                            if ch.worst
                            else "none"
                        ),
                        ch.verdict(),
                    ]
                    for name, ch in design.cross_hyb.items()
                ],
            }
        )
    if design.sense:
        out["tables"].append(
            {
                "title": "Sense control — order as a separate PCR",
                "columns": ["oligo", "sequence", "length"],
                "rows": [[o.name, o.sequence, o.length] for o in design.sense.oligos],
            }
        )
    return out


def _mode_design(payload: dict[str, Any]) -> dict[str, Any]:
    template = parse_one_sequence(payload.get("template"), what="template", default_name="template")
    params = build_params(payload)
    panel = build_panel(payload)

    def specs(key: str) -> IntervalSet:
        raw = payload.get(key) or []
        if isinstance(raw, str):
            raw = [piece.strip() for piece in raw.replace(",", " ").split() if piece.strip()]
        return IntervalSet.from_specs(raw) if raw else IntervalSet()

    result = run_design(
        DesignRequest(
            template=template,
            params=params,
            targets=specs("targets"),
            excluded=specs("excludes"),
            panel=panel,
            name=payload.get("name") or template.name,
            context={"assay": "generic"},
        )
    )
    return _package_result(result, tails=payload.get("tails") or "none")


def _mode_check(payload: dict[str, Any]) -> dict[str, Any]:
    from .design.candidate import Primer, TemplateContext, evaluate_primer
    from .thermo.backend import get_backend
    from .thermo.structure import hetero_dimer, three_prime_dimer_score

    params = build_params(payload)
    panel = build_panel(payload)
    text = (payload.get("primers") or "").strip()
    if not text:
        raise DanioPrimeError("paste at least one primer")

    entries: list[tuple[str, str]] = []
    if text.lstrip().startswith(">"):
        for rec in parse_sequences(text):
            entries.append((rec.name, rec.seq))
    else:
        for piece in text.replace(",", "\n").splitlines():
            piece = piece.strip()
            if not piece:
                continue
            if "=" in piece:
                name, _, seq = piece.partition("=")
                entries.append((name.strip(), normalise(seq.strip())))
            else:
                entries.append((f"primer{len(entries) + 1}", normalise(piece)))
    if not entries:
        raise DanioPrimeError("paste at least one primer")

    backend = get_backend(params.thermo_backend, params.conditions, params.anneal_temp_c)
    rows = []
    evaluated = []
    for name, seq in entries:
        ctx = TemplateContext(seq=seq, params=params, backend=backend, name=name)
        primer = Primer(seq=seq, start=0, end=len(seq), orientation="F")
        evaluate_primer(primer, ctx)
        evaluated.append((name, seq, primer))
        rows.append(
            [
                name, seq, primer.length,
                f"{primer.tm:.1f}" if primer.tm == primer.tm else "n/a",
                f"{primer.gc:.0f}",
                f"{primer.hairpin_dg:.1f}",
                f"{primer.self_dimer_dg:.1f}",
                primer.gc_clamp, primer.poly_x,
            ]
        )

    tables = [
        {
            "title": "Primer properties",
            "columns": ["primer", "sequence", "length", "Tm", "GC%", "hairpin dG",
                        "self-dimer dG", "3' clamp", "poly-X"],
            "rows": rows,
        }
    ]
    lines = [f"Evaluated {len(evaluated)} primer(s)."]

    if len(evaluated) > 1:
        pair_rows = []
        for i in range(len(evaluated)):
            for j in range(i + 1, len(evaluated)):
                an, a_seq, _ = evaluated[i]
                bn, b_seq, _ = evaluated[j]
                hd = hetero_dimer(
                    a_seq, b_seq, temp_c=params.anneal_temp_c, conditions=params.conditions
                )
                pair_rows.append(
                    [
                        f"{an} x {bn}",
                        f"{hd.dg:.1f}",
                        f"{three_prime_dimer_score(a_seq, b_seq):.0f}",
                        "watch" if hd.dg < params.max_hetero_dimer_dg else "ok",
                    ]
                )
        tables.append(
            {
                "title": "Pairwise interactions",
                "columns": ["pair", "dimer dG (kcal/mol)", "3'-end score", "verdict"],
                "rows": pair_rows,
            }
        )

    if panel:
        panel_rows = []
        for name, seq, _ in evaluated:
            report = assess_primer(seq, panel, primer_name=name, conditions=params.conditions)
            panel_rows.append(
                [
                    name,
                    (
                        f"{report.best_site.subject} ({report.best_site.mismatches} mm)"
                        if report.best_site
                        else "none"
                    ),
                    "YES" if report.primes_on_panel else "no",
                    report.verdict(),
                ]
            )
        tables.append(
            {
                "title": f"Against the avoid-panel ({len(panel)} sequence(s))",
                "columns": ["primer", "closest match", "primes?", "verdict"],
                "rows": panel_rows,
            }
        )
        if len(evaluated) >= 2:
            amplicons, _ = assess_pair(
                evaluated[0][1], evaluated[1][1], panel,
                min_size=params.min_product_size,
                max_size=max(params.max_product_size, 5000),
                forward_name=evaluated[0][0], reverse_name=evaluated[1][0],
                conditions=params.conditions,
            )
            if amplicons:
                tables.append(
                    {
                        "title": "Predicted products on the avoid-panel",
                        "columns": ["subject", "size", "forward mm", "reverse mm", "note"],
                        "rows": [
                            [
                                a.subject, a.size,
                                a.forward_site.mismatches, a.reverse_site.mismatches,
                                a.note,
                            ]
                            for a in amplicons
                        ],
                    }
                )
                lines.append("A product IS predicted on the avoid-panel.")
            else:
                lines.append("No product predicted on the avoid-panel.")

    return {
        "ok": True,
        "mode": "check",
        "report": "\n".join(lines),
        "pairs": [],
        "warnings": [],
        "stats": {"primers": len(evaluated)},
        "diagnosis": "",
        "tables": tables,
        "tsv": "\n".join(
            ["\t".join(tables[0]["columns"])]
            + ["\t".join(str(c) for c in r) for r in tables[0]["rows"]]
        ),
        "json": json.dumps(
            {"mode": "check", "tables": tables}, indent=2
        ),
    }


_HANDLERS = {
    "genotype": _mode_genotype,
    "qpcr": _mode_qpcr,
    "dcaps": _mode_dcaps,
    "insitu": _mode_insitu,
    "design": _mode_design,
    "check": _mode_check,
}


def run(payload: dict[str, Any] | str) -> dict[str, Any]:
    """Run one design from a JSON-shaped payload and return a JSON-shaped result.

    Never raises: a failure comes back as ``{"ok": False, "error": ...}`` with a
    message meant for the person at the keyboard, because a traceback in a
    browser console helps nobody. The traceback is kept under ``detail`` for bug
    reports.
    """
    if isinstance(payload, str):
        payload = json.loads(payload)
    try:
        mode = payload.get("mode") or "design"
        handler = _HANDLERS.get(mode)
        if handler is None:
            return {
                "ok": False,
                "error": f"unknown mode {mode!r}; choose from {', '.join(MODES)}",
            }
        return handler(payload)
    except DanioPrimeError as exc:
        return {"ok": False, "error": str(exc)}
    except (ValueError, TypeError, KeyError) as exc:
        return {
            "ok": False,
            "error": f"could not read the form input: {exc}",
            "detail": traceback.format_exc(),
        }
    except Exception as exc:  # noqa: BLE001 - last resort, must not escape to JS
        return {
            "ok": False,
            "error": f"unexpected error: {type(exc).__name__}: {exc}",
            "detail": traceback.format_exc(),
        }


def run_json(payload: dict[str, Any] | str) -> str:
    """:func:`run`, serialised. Convenient across the Python/JS boundary."""
    return json.dumps(run(payload))


def ui_metadata() -> dict[str, Any]:
    """Everything the form needs to build itself, so nothing is hard-coded in HTML."""
    from . import __version__

    return {
        "version": __version__,
        "modes": list(MODES),
        "presets": [
            {"name": name, "description": desc}
            for name, desc in sorted(PRESET_DESCRIPTIONS.items())
        ],
        "tail_schemes": [
            {"name": name, "description": SCHEME_DESCRIPTIONS.get(name, "")}
            for name in sorted(TAIL_SCHEMES)
        ],
        "riboprobe_schemes": [
            {"name": name, "description": SCHEME_DESCRIPTIONS.get(name, "")}
            for name in sorted(RIBOPROBE_SCHEMES)
        ],
        "nucleases": [
            {"name": name, "description": spec["description"]}
            for name, spec in sorted(NUCLEASES.items())
        ],
        "tunable": list(TUNABLE),
        "defaults": {k: getattr(DesignParams(), k) for k in TUNABLE},
    }
