"""Optional external aligners for genome-scale specificity checks.

The built-in scanners are exact and dependency-free but scale to roughly a few
megabases per primer. A whole zebrafish genome is ~1.4 Gb, so genome-wide
checks delegate to a real aligner when one is on ``PATH``:

``blastn``
    Short-query BLAST with the parameters NCBI Primer-BLAST uses
    (``-task blastn-short``, ``-word_size 7``). Needs a BLAST database built
    once with ``makeblastdb``.
``bowtie2``
    End-to-end alignment reporting many alignments per read. Fast, needs an
    index built once with ``bowtie2-build``.

Both are used purely as *candidate site finders*. Their hits are re-judged with
the same :class:`~danioprime.specificity.scan.PrimingRules` as the built-in
scanner, so a primer's verdict does not depend on which tool found the site.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from ..errors import ExternalToolError
from ..thermo.nn import DEFAULT_CONDITIONS, SaltConditions, imperfect_duplex_tm
from .scan import DEFAULT_RULES, BindingSite, PrimingRules


def which(tool: str) -> str | None:
    """Path to *tool* if it is callable, else None."""
    return shutil.which(tool)


def available_aligners() -> dict[str, str | None]:
    """Report which optional aligners are installed."""
    return {name: which(name) for name in ("blastn", "makeblastdb", "bowtie2", "bowtie2-build")}


@dataclass(slots=True)
class ExternalHit:
    """One alignment reported by an external tool, before priming rules."""

    subject: str
    start: int
    end: int
    strand: str
    query_name: str
    identity: float = 0.0
    aligned_length: int = 0
    subject_seq: str = ""


def _run(cmd: list[str], *, stdin: str | None = None, timeout: int = 1800) -> str:
    try:
        proc = subprocess.run(  # noqa: S603 - argv list, no shell
            cmd,
            input=stdin,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except FileNotFoundError as exc:
        raise ExternalToolError(f"{cmd[0]} is not installed or not on PATH") from exc
    except subprocess.TimeoutExpired as exc:
        raise ExternalToolError(f"{cmd[0]} timed out after {timeout}s") from exc
    if proc.returncode != 0:
        tail = (proc.stderr or "").strip().splitlines()[-5:]
        raise ExternalToolError(f"{cmd[0]} failed ({proc.returncode}): {' | '.join(tail)}")
    return proc.stdout


def blastn_sites(
    primers: dict[str, str],
    database: str | Path,
    *,
    word_size: int = 7,
    evalue: float = 1000.0,
    max_targets: int = 500,
    threads: int = 1,
    rules: PrimingRules = DEFAULT_RULES,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    max_mismatches: int = 6,
) -> dict[str, list[BindingSite]]:
    """Find binding sites with ``blastn -task blastn-short``.

    *database* is a BLAST database prefix, e.g. the output of::

        makeblastdb -in GRCz12tu.fa -dbtype nucl -out GRCz12tu

    The aligned subject sequence is requested in the tabular output so that
    mismatch positions relative to the primer's 3' end can be recovered without
    a second pass over the reference.
    """
    if not which("blastn"):
        raise ExternalToolError("blastn is not installed; `conda install -c bioconda blast`")
    fasta = "".join(f">{name}\n{seq}\n" for name, seq in primers.items())
    fields = "qseqid sseqid sstart send sstrand nident length qstart qend qseq sseq mismatch"
    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-db", str(database),
        "-word_size", str(word_size),
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_targets),
        "-num_threads", str(threads),
        "-dust", "no",
        "-soft_masking", "false",
        "-outfmt", f"6 {fields}",
    ]
    out = _run(cmd, stdin=fasta)
    return _parse_blast(out, primers, rules=rules, conditions=conditions, max_mismatches=max_mismatches)


def _parse_blast(
    text: str,
    primers: dict[str, str],
    *,
    rules: PrimingRules,
    conditions: SaltConditions,
    max_mismatches: int,
) -> dict[str, list[BindingSite]]:
    sites: dict[str, list[BindingSite]] = {name: [] for name in primers}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        qid, sid, sstart, send, sstrand, _nident, _length, qstart, qend, _qseq, sseq, _mm = parts[:12]
        primer_seq = primers.get(qid)
        if primer_seq is None:
            continue
        lp = len(primer_seq)
        # Only full-length alignments can be judged for priming; a partial
        # alignment that stops short of the 3' end cannot extend.
        if int(qstart) != 1 or int(qend) != lp:
            continue
        s1, s2 = int(sstart), int(send)
        lo, hi = min(s1, s2) - 1, max(s1, s2)
        strand = "+" if sstrand.lower().startswith("plus") else "-"
        # blastn writes sseq in the *query's* alignment register: for a
        # minus-strand hit it is already the strand the primer anneals to, read
        # 5'->3' alongside the primer. So position i of sseq pairs with position
        # i of the primer on both strands, and no complementing is needed --
        # complementing it here would report a perfect minus-strand hit as
        # twenty mismatches.
        aligned_subject = sseq.replace("-", "")
        if len(aligned_subject) != lp:
            continue
        template_register = aligned_subject
        offsets = [
            lp - 1 - i for i in range(lp) if aligned_subject[i] != primer_seq[i]
        ]
        mism = len(offsets)
        if mism > max_mismatches:
            continue
        try:
            tm = imperfect_duplex_tm(primer_seq, template_register, conditions)
        except Exception:  # noqa: BLE001
            tm = float("nan")
        offsets_sorted = tuple(sorted(offsets))
        can_prime, note = rules.judge(mism, offsets_sorted, tm if tm == tm else -999.0)
        sites[qid].append(
            BindingSite(
                subject=sid,
                start=lo,
                end=hi,
                strand=strand,
                mismatches=mism,
                mismatch_offsets_from_3p=offsets_sorted,
                three_prime_mismatches=sum(1 for o in offsets_sorted if o < rules.three_prime_window),
                tm_estimate=tm,
                primer_name=qid,
                can_prime=can_prime,
                note=note,
            )
        )
    for name in sites:
        sites[name].sort(key=lambda s: (s.mismatches, s.subject, s.start))
    return sites


def bowtie2_sites(
    primers: dict[str, str],
    index_prefix: str | Path,
    *,
    max_alignments: int = 100,
    threads: int = 1,
    reference_fasta: str | Path | None = None,
    rules: PrimingRules = DEFAULT_RULES,
    conditions: SaltConditions = DEFAULT_CONDITIONS,
    max_mismatches: int = 6,
) -> dict[str, list[BindingSite]]:
    """Find binding sites with ``bowtie2 -k`` end-to-end alignment.

    Bowtie2's SAM output does not carry the reference window, so mismatch
    positions are reconstructed from the MD tag. Alignments with insertions or
    deletions are discarded, since a gapped primer-template duplex is not a
    priming geometry this model covers.
    """
    if not which("bowtie2"):
        raise ExternalToolError("bowtie2 is not installed; `conda install -c bioconda bowtie2`")
    with tempfile.TemporaryDirectory(prefix="danioprime_") as tmp:
        reads = Path(tmp) / "primers.fa"
        reads.write_text("".join(f">{n}\n{s}\n" for n, s in primers.items()))
        cmd = [
            "bowtie2",
            "-f",
            "-x", str(index_prefix),
            "-U", str(reads),
            "-k", str(max_alignments),
            "--end-to-end",
            "--very-sensitive",
            "-N", "1",
            "-L", "10",
            "--no-unal",
            "-p", str(threads),
        ]
        out = _run(cmd)
    _ = reference_fasta  # retained for API symmetry; MD tags make it unnecessary
    return _parse_sam(out, primers, rules=rules, conditions=conditions, max_mismatches=max_mismatches)


def _md_mismatch_offsets(md: str, cigar: str, length: int, strand: str) -> list[int] | None:
    """Positions (0-based along the read) that mismatch the reference.

    Returns None when the alignment is gapped or the MD tag is unusable.
    """
    if any(op in cigar for op in ("I", "D", "N", "S", "H", "P")):
        return None
    offsets: list[int] = []
    pos = 0
    num = ""
    i = 0
    while i < len(md):
        ch = md[i]
        if ch.isdigit():
            num += ch
            i += 1
            continue
        if num:
            pos += int(num)
            num = ""
        if ch == "^":
            return None  # deletion from the reference
        offsets.append(pos)
        pos += 1
        i += 1
    if pos > length:
        return None
    # Read coordinates are reference-forward; flip for minus-strand alignments
    # so that offsets are expressed along the primer as ordered.
    if strand == "-":
        offsets = [length - 1 - o for o in offsets]
    return offsets


def _parse_sam(
    text: str,
    primers: dict[str, str],
    *,
    rules: PrimingRules,
    conditions: SaltConditions,
    max_mismatches: int,
) -> dict[str, list[BindingSite]]:
    sites: dict[str, list[BindingSite]] = {name: [] for name in primers}
    for line in text.splitlines():
        if line.startswith("@") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        qname, flag, rname, pos, _mapq, cigar = f[0], int(f[1]), f[2], int(f[3]), f[4], f[5]
        if flag & 0x4 or rname == "*":
            continue
        primer_seq = primers.get(qname)
        if primer_seq is None:
            continue
        lp = len(primer_seq)
        strand = "-" if flag & 0x10 else "+"
        md = next((t.split(":", 2)[2] for t in f[11:] if t.startswith("MD:Z:")), None)
        if md is None:
            continue
        read_offsets = _md_mismatch_offsets(md, cigar, lp, strand)
        if read_offsets is None:
            continue
        mism = len(read_offsets)
        if mism > max_mismatches:
            continue
        # Convert read offsets to distance-from-3'-end of the primer.
        offsets = sorted(lp - 1 - o for o in read_offsets)
        start = pos - 1
        tm = float("nan")
        offsets_t = tuple(offsets)
        can_prime, note = rules.judge(mism, offsets_t, 60.0 if mism <= 2 else 50.0)
        sites[qname].append(
            BindingSite(
                subject=rname,
                start=start,
                end=start + lp,
                strand=strand,
                mismatches=mism,
                mismatch_offsets_from_3p=offsets_t,
                three_prime_mismatches=sum(1 for o in offsets_t if o < rules.three_prime_window),
                tm_estimate=tm,
                primer_name=qname,
                can_prime=can_prime,
                note=note or "Tm not estimated (bowtie2 backend)",
            )
        )
        _ = conditions
    for name in sites:
        sites[name].sort(key=lambda s: (s.mismatches, s.subject, s.start))
    return sites
