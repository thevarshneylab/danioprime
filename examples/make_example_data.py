#!/usr/bin/env python3
"""Generate the example dataset.

Synthetic, but built to have the properties that make zebrafish design hard, so
the examples exercise real behaviour rather than a happy path:

* ``kcnh6a`` and ``kcnh6b`` are a mock ohnolog pair. Their coding exons are 94%
  identical with a block of near-perfect conservation, which is exactly where a
  naive design lands and the reason the pair gets co-amplified.
* Their 3' UTRs are unrelated, as real ohnolog UTRs largely are: UTR sequence
  diverges far faster than coding sequence. That matters for in situ probes,
  where the coding region is often too conserved for any probe to discriminate
  and the UTR is the only place a specific probe can be placed.
* The intergenic and intronic sequence diverges freely, so primers placed there
  discriminate easily. The contrast between the two is the point.
* A SpCas9 protospacer with a real NGG PAM sits in exon 2 of ``kcnh6a``.
* A VCF carries TU- and AB-specific variants, including one sitting in an
  otherwise attractive primer site, and one homozygous in TU (which is *not* a
  hazard and should be recovered by ``--heterozygous-only``).

Nothing here is real zebrafish sequence; the gene names are placeholders chosen
because the ``a``/``b`` suffix convention is what the ohnolog logic keys on.
"""

from __future__ import annotations

import random
from pathlib import Path

HERE = Path(__file__).parent
RNG = random.Random(20260211)


def rand_seq(n: int) -> str:
    return "".join(RNG.choice("ACGT") for _ in range(n))


def mutate(seq: str, rate: float) -> str:
    out = []
    for base in seq:
        if RNG.random() < rate:
            out.append(RNG.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out)


def write_fasta(path: Path, records: list[tuple[str, str, str]], wrap: int = 70) -> None:
    with path.open("w") as fh:
        for name, desc, seq in records:
            fh.write(f">{name}" + (f" {desc}" if desc else "") + "\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i : i + wrap] + "\n")


def main() -> None:
    # --- exons shared between the two ohnologs ---------------------------
    exon1 = rand_seq(210)
    exon2 = rand_seq(330)
    exon3_cds = rand_seq(260)
    # A 120 bp block of near-perfect conservation inside exon 2: the trap.
    conserved = exon2[100:220]

    # --- kcnh6a locus: exon1 - intron1 - exon2 - intron2 - exon3 ---------
    a_intron1 = rand_seq(900)
    a_intron2 = rand_seq(1100)
    a_upstream = rand_seq(600)
    a_downstream = rand_seq(600)

    # Plant a SpCas9 protospacer with an NGG PAM inside exon 2 of kcnh6a,
    # outside the conserved block so the genotyping example is not forced to
    # sit on paralogous sequence.
    protospacer = rand_seq(20)
    exon2_a = exon2[:40] + protospacer + "TGG" + exon2[63:]
    assert len(exon2_a) == len(exon2)

    # The last exon carries a 3' UTR. The two genes' UTRs are generated
    # independently rather than mutated from a common ancestor, which is the
    # point: a probe over the CDS cannot discriminate these ohnologs, a probe
    # over the UTR can. Drawn after the protospacer so that adding them did not
    # shift the coordinates already published in the README and demo.
    utr3_a = rand_seq(320)
    utr3_b = rand_seq(320)
    exon3 = exon3_cds + utr3_a

    a_locus = (
        a_upstream + exon1 + a_intron1 + exon2_a + a_intron2 + exon3 + a_downstream
    )
    a_exon_starts = [
        len(a_upstream),
        len(a_upstream) + len(exon1) + len(a_intron1),
        len(a_upstream) + len(exon1) + len(a_intron1) + len(exon2_a) + len(a_intron2),
    ]
    a_exon_lens = [len(exon1), len(exon2_a), len(exon3)]

    # --- kcnh6b locus: diverged except for the conserved exon block ------
    b_exon1 = mutate(exon1, 0.06)
    b_exon2 = mutate(exon2, 0.06)
    # restore the conserved block almost exactly (2 substitutions over 120 bp)
    b_conserved = mutate(conserved, 0.017)
    b_exon2 = b_exon2[:100] + b_conserved + b_exon2[220:]
    b_exon3 = mutate(exon3_cds, 0.06) + utr3_b
    b_locus = (
        rand_seq(500) + b_exon1 + rand_seq(820) + b_exon2 + rand_seq(1000) + b_exon3 + rand_seq(500)
    )

    write_fasta(
        HERE / "kcnh6a_locus.fa",
        [("kcnh6a_locus", "mock chr5 locus, 3 exons, SpCas9 protospacer in exon 2", a_locus)],
    )
    write_fasta(
        HERE / "kcnh6b_locus.fa",
        [("kcnh6b_locus", "mock ohnolog of kcnh6a: 94% identical exons", b_locus)],
    )

    # --- spliced transcripts --------------------------------------------
    a_tx = exon1 + exon2_a + exon3
    b_tx = b_exon1 + b_exon2 + b_exon3
    cds_end = len(exon1) + len(exon2_a) + len(exon3_cds)
    write_fasta(
        HERE / "kcnh6a_transcript.fa",
        [(
            "kcnh6a-201",
            f"spliced transcript, junctions at {len(exon1)} and {len(exon1) + len(exon2_a)}, "
            f"3' UTR from {cds_end + 1}",
            a_tx,
        )],
    )
    write_fasta(
        HERE / "kcnh6b_transcript.fa",
        [("kcnh6b-201", "ohnolog transcript: avoid-panel for qPCR and in situ", b_tx)],
    )

    # --- GTF for the kcnh6a locus (coordinates within the locus FASTA) ---
    with (HERE / "kcnh6a.gtf").open("w") as fh:
        fh.write("#!genome-build mock-GRCz12tu-like\n")
        fh.write("# coordinates are within kcnh6a_locus.fa, not a real genome\n")
        for i, (start, length) in enumerate(zip(a_exon_starts, a_exon_lens, strict=True), 1):
            attrs = (
                'gene_id "ENSDARGMOCK001"; transcript_id "kcnh6a-201"; '
                f'gene_name "kcnh6a"; exon_number "{i}"; '
                'transcript_biotype "protein_coding";'
            )
            fh.write(
                f"kcnh6a_locus\tmock\texon\t{start + 1}\t{start + length}\t.\t+\t.\t{attrs}\n"
            )

    # --- strain variants -------------------------------------------------
    # Place variants deterministically at positions that fall in otherwise
    # attractive primer territory flanking the cut site.
    cut = a_exon_starts[1] + 40 + 17  # protospacer starts at exon2+40
    positions = [cut - 180, cut - 120, cut + 150, cut + 240, cut + 600]
    with (HERE / "strain_variants.vcf").open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=danioprime example data (synthetic)\n")
        fh.write(f"##contig=<ID=kcnh6a_locus,length={len(a_locus)}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTU\tAB\n")
        genotypes = [("0/1", "0/0"), ("0/1", "0/1"), ("1/1", "0/0"), ("0/0", "0/1"), ("0/1", "0/0")]
        for i, (pos, (tu, ab)) in enumerate(zip(positions, genotypes, strict=True), 1):
            ref = a_locus[pos]
            alt = RNG.choice([b for b in "ACGT" if b != ref])
            fh.write(
                f"kcnh6a_locus\t{pos + 1}\tex_snp{i}\t{ref}\t{alt}\t60\tPASS\t.\tGT\t{tu}\t{ab}\n"
            )

    # --- a knock-in cassette --------------------------------------------
    cassette = (
        "ATGGTGAGCAAGGGCGAGGAG"  # mock reporter start
        + rand_seq(680)
        + "TTACTTGTACAGCTCGTCCATGCC"  # mock reporter end
    )
    write_fasta(HERE / "reporter_cassette.fa", [("mScarlet_cassette", "mock knock-in cassette", cassette)])

    # --- summary ---------------------------------------------------------
    summary = HERE / "DATA.md"
    summary.write_text(
        f"""# Example dataset (synthetic)

Regenerate with `python make_example_data.py`. Deterministic: seed 20260211.

| file | contents |
|---|---|
| `kcnh6a_locus.fa` | {len(a_locus)} bp mock locus, 3 exons, SpCas9 protospacer in exon 2 |
| `kcnh6b_locus.fa` | {len(b_locus)} bp mock ohnolog, exons ~94% identical to kcnh6a |
| `kcnh6a_transcript.fa` | {len(a_tx)} bp spliced transcript |
| `kcnh6b_transcript.fa` | {len(b_tx)} bp ohnolog transcript (avoid-panel for qPCR) |
| `kcnh6a.gtf` | exon structure in `kcnh6a_locus.fa` coordinates |
| `strain_variants.vcf` | 5 variants with TU/AB genotypes |
| `reporter_cassette.fa` | {len(cassette)} bp mock knock-in cassette |

Key coordinates in `kcnh6a_locus.fa` (1-based):

- exon 1: {a_exon_starts[0] + 1}-{a_exon_starts[0] + a_exon_lens[0]}
- exon 2: {a_exon_starts[1] + 1}-{a_exon_starts[1] + a_exon_lens[1]}
- exon 3: {a_exon_starts[2] + 1}-{a_exon_starts[2] + a_exon_lens[2]}
- protospacer: {a_exon_starts[1] + 41}-{a_exon_starts[1] + 60} (PAM {a_exon_starts[1] + 61}-{a_exon_starts[1] + 63})
- SpCas9 cut site: {cut + 1}
- conserved block shared with kcnh6b: transcript {len(exon1) + 101}-{len(exon1) + 220}

Transcript junctions: {len(exon1)} and {len(exon1) + len(exon2_a)}.

Coding sequence ends at transcript position {cds_end}; the 3' UTR
({cds_end + 1}-{len(a_tx)}) is unrelated between the two genes. A probe over the
CDS cannot discriminate the ohnologs; a probe over the UTR can. That is the
contrast `danioprime insitu` is there to surface.

Nothing here is real zebrafish sequence. The `a`/`b` names are placeholders
chosen because the ZFIN suffix convention is what the ohnolog logic keys on.
"""
    )

    print(f"protospacer: {protospacer}")
    print(f"cut site (1-based): {cut + 1}")
    print(f"locus lengths: kcnh6a {len(a_locus)}, kcnh6b {len(b_locus)}")
    print(f"transcript junctions: {len(exon1)}, {len(exon1) + len(exon2_a)}")
    print(f"conserved block in transcript: {len(exon1) + 101}-{len(exon1) + 220}")
    print(f"CDS ends at transcript {cds_end}; divergent 3' UTR {cds_end + 1}-{len(a_tx)}")


if __name__ == "__main__":
    main()
