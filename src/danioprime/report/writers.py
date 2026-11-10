"""Output formats.

Four of them, because primer design output is consumed in four different ways:

``text``
    What you read in the terminal to decide whether the design is any good.
``tsv``
    What you paste into a lab notebook or a spreadsheet.
``json``
    Everything, including the penalty breakdown and every predicted off-target
    site, for programmatic use and for reproducing a decision later.
``bed``
    Primer footprints to load into IGV or the UCSC browser alongside the locus,
    which is how you catch a primer sitting in a repeat that the numbers did
    not convey.
"""

from __future__ import annotations

import csv
import json
import sys
from collections.abc import Iterable, Sequence
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, TextIO

from ..design.pairing import PrimerPair
from ..design.pipeline import DesignResult

PAIR_COLUMNS: tuple[str, ...] = (
    "pair",
    "forward_seq",
    "reverse_seq",
    "product_size",
    "forward_tm",
    "reverse_tm",
    "tm_diff",
    "forward_gc",
    "reverse_gc",
    "forward_pos",
    "reverse_pos",
    "hetero_dimer_dg",
    "specific",
    "paralog_amplicons",
    "penalty",
    "notes",
)


def _pair_row(pair: PrimerPair) -> dict[str, Any]:
    return {
        "pair": pair.name,
        "forward_seq": pair.forward.ordered_seq,
        "reverse_seq": pair.reverse.ordered_seq,
        "product_size": pair.product_size,
        "forward_tm": f"{pair.forward.tm:.1f}",
        "reverse_tm": f"{pair.reverse.tm:.1f}",
        "tm_diff": f"{pair.tm_difference:.1f}",
        "forward_gc": f"{pair.forward.gc:.0f}",
        "reverse_gc": f"{pair.reverse.gc:.0f}",
        "forward_pos": pair.forward.footprint.as_one_based(),
        "reverse_pos": pair.reverse.footprint.as_one_based(),
        "hetero_dimer_dg": f"{pair.hetero_dimer_dg:.1f}",
        "specific": "yes" if pair.is_specific else "no",
        "paralog_amplicons": ";".join(
            f"{a['subject']}:{a['size']}bp" for a in pair.paralog_amplicons
        ),
        "penalty": f"{pair.penalty:.2f}",
        "notes": "; ".join(pair.notes),
    }


def write_tsv(pairs: Sequence[PrimerPair], path: str | Path | None = None) -> str:
    """Write (or return) the pair table as TSV."""
    rows = [_pair_row(p) for p in pairs]
    if path is None:
        import io

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=PAIR_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
        return buf.getvalue()
    with Path(path).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PAIR_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


def _jsonable(obj: Any) -> Any:
    if is_dataclass(obj) and not isinstance(obj, type):
        return _jsonable(asdict(obj))
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, float):
        return round(obj, 4) if obj == obj else None
    return obj


def result_to_dict(result: DesignResult, *, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    """Full structured representation of a design run."""
    payload: dict[str, Any] = {
        "template": {
            "name": result.template.name,
            "length": len(result.template.seq),
            "description": result.template.description,
        },
        "assay": result.request.context.get("assay", "generic"),
        "context": _jsonable(result.request.context),
        "parameters": _jsonable(result.request.params.to_dict()),
        "statistics": _jsonable(result.stats),
        "warnings": list(result.warnings),
        "pairs": [p.to_dict() for p in result.pairs],
    }
    if result.discrimination:
        payload["specificity"] = _jsonable(result.discrimination)
    if extra:
        payload.update(_jsonable(extra))
    return payload


def write_json(
    result: DesignResult, path: str | Path | None = None, *, extra: dict[str, Any] | None = None
) -> str:
    """Write (or return) the complete result as JSON."""
    text = json.dumps(result_to_dict(result, extra=extra), indent=2) + "\n"
    if path is None:
        return text
    Path(path).write_text(text, newline="\n")
    return str(path)


def write_bed(
    pairs: Sequence[PrimerPair],
    path: str | Path | None = None,
    *,
    chrom: str = "template",
    offset: int = 0,
) -> str:
    """Primer footprints as BED.

    Pass ``--template-offset`` as *offset* and the contig name as *chrom* when
    the template was cut out of a genome, and the output lands in the right
    place in a browser.
    """
    lines = ["track name=danioprime description=\"primer footprints\" itemRgb=On"]
    for pair in pairs:
        for primer, colour in ((pair.forward, "27,120,55"), (pair.reverse, "166,31,55")):
            strand = "+" if primer.orientation == "F" else "-"
            lines.append(
                "\t".join(
                    [
                        chrom,
                        str(primer.start + offset),
                        str(primer.end + offset),
                        f"{pair.name}_{primer.orientation}",
                        str(min(1000, int(primer.penalty * 100))),
                        strand,
                        str(primer.start + offset),
                        str(primer.end + offset),
                        colour,
                    ]
                )
            )
        lines.append(
            "\t".join(
                [
                    chrom,
                    str(pair.forward.start + offset),
                    str(pair.reverse.end + offset),
                    f"{pair.name}_amplicon",
                    "0",
                    "+",
                    str(pair.forward.start + offset),
                    str(pair.reverse.end + offset),
                    "120,120,120",
                ]
            )
        )
    text = "\n".join(lines) + "\n"
    if path is None:
        return text
    Path(path).write_text(text, newline="\n")
    return str(path)


def write_order_sheet(items: Sequence[Any], path: str | Path | None = None) -> str:
    """Write (or return) an oligo order sheet as TSV."""
    columns = ("name", "sequence", "length", "role", "tail", "gc_percent", "note")
    rows = [i.to_dict() if hasattr(i, "to_dict") else dict(i) for i in items]
    if path is None:
        import io

        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
        return buf.getvalue()
    with Path(path).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


def _bar(value: float, width: int = 12, scale: float = 10.0) -> str:
    filled = max(0, min(width, int(round(value / scale * width))))
    return "#" * filled + "." * (width - filled)


def format_text(
    result: DesignResult,
    *,
    verbose: bool = False,
    show: int = 5,
    extra_sections: Sequence[tuple[str, Sequence[str]]] = (),
) -> str:
    """Human-readable summary for the terminal."""
    out: list[str] = []
    tpl = result.template
    assay = result.request.context.get("assay", "generic")
    out.append(f"danioprime: {assay} design on {tpl.name} ({len(tpl.seq):,} nt)")

    ctx = result.request.context
    if "cut_site" in ctx:
        out.append(f"  cut site: {ctx['cut_site'].get('description', '')}")
    if result.request.panel:
        out.append(
            f"  avoid-panel: {len(result.request.panel)} sequence(s), "
            f"{result.request.panel.total_length:,} nt "
            f"({result.request.panel.resolved_strategy()} scan)"
        )
    if result.request.variants:
        out.append(
            f"  variant mask: {result.request.variants.total_length} position(s) excluded"
        )

    stats = result.stats
    out.append(
        f"  candidates: {stats.get('forward_candidates', 0)} forward, "
        f"{stats.get('reverse_candidates', 0)} reverse"
    )

    for title, lines in extra_sections:
        out.append("")
        out.append(title)
        out.extend(f"  {line}" for line in lines)

    if not result.pairs:
        out.append("")
        out.append("No primer pair satisfied the constraints.")
        if result.warnings:
            for w in result.warnings:
                out.append(f"  {w}")
        return "\n".join(out)

    out.append("")
    out.append(f"Top {min(show, len(result.pairs))} of {len(result.pairs)} pair(s):")
    for pair in result.pairs[:show]:
        out.append("")
        flag = "" if pair.is_specific else "  [!] cross-amplification predicted"
        out.append(f"  {pair.name}  {pair.product_size} bp  penalty {pair.penalty:.2f}{flag}")
        for label, primer in (("F", pair.forward), ("R", pair.reverse)):
            tail = f" (+{primer.tail_name})" if primer.tail_name else ""
            out.append(
                f"    {label}  {primer.ordered_seq}{tail}"
            )
            out.append(
                f"       {primer.length} nt  Tm {primer.tm:.1f}C  GC {primer.gc:.0f}%  "
                f"pos {primer.footprint.as_one_based()}  hairpin {primer.hairpin_dg:.1f}"
            )
            if primer.notes:
                for note in primer.notes:
                    out.append(f"       - {note}")
        out.append(
            f"    pair  dTm {pair.tm_difference:.1f}C  "
            f"hetero-dimer {pair.hetero_dimer_dg:.1f} kcal/mol  "
            f"3' dimer score {pair.three_prime_dimer_score:.0f}"
        )
        for note in pair.notes:
            out.append(f"    - {note}")
        if verbose and pair.penalty_terms:
            out.append("    penalty breakdown:")
            for term, value in sorted(pair.penalty_terms.items(), key=lambda kv: -abs(kv[1])):
                out.append(f"       {term:<22} {value:7.2f}  {_bar(abs(value))}")

    if result.warnings:
        out.append("")
        out.append("Warnings:")
        for w in result.warnings:
            out.append(f"  - {w}")
    return "\n".join(out)


def emit(text: str, path: str | Path | None, stream: TextIO | None = None) -> None:
    """Write *text* to *path*, or to *stream* (default stdout) when path is None."""
    if path is None:
        (stream or sys.stdout).write(text if text.endswith("\n") else text + "\n")
    else:
        Path(path).write_text(
            text if text.endswith("\n") else text + "\n", newline="\n"
        )


def summarise_iterable(rows: Iterable[dict[str, Any]], columns: Sequence[str]) -> str:
    """Render a list of dicts as a fixed-width table for the text report."""
    rows = list(rows)
    if not rows:
        return "(none)"
    widths = {c: max(len(c), *(len(str(r.get(c, ""))) for r in rows)) for c in columns}
    lines = ["  ".join(c.ljust(widths[c]) for c in columns)]
    lines.append("  ".join("-" * widths[c] for c in columns))
    for row in rows:
        lines.append("  ".join(str(row.get(c, "")).ljust(widths[c]) for c in columns))
    return "\n".join(lines)
