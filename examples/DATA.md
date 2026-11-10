# Example dataset (synthetic)

Regenerate with `python make_example_data.py`. Deterministic: seed 20260211.

| file | contents |
|---|---|
| `kcnh6a_locus.fa` | 4000 bp mock locus, 3 exons, SpCas9 protospacer in exon 2 |
| `kcnh6b_locus.fa` | 3620 bp mock ohnolog, exons ~94% identical to kcnh6a |
| `kcnh6a_transcript.fa` | 800 bp spliced transcript |
| `kcnh6b_transcript.fa` | 800 bp ohnolog transcript (avoid-panel for qPCR) |
| `kcnh6a.gtf` | exon structure in `kcnh6a_locus.fa` coordinates |
| `strain_variants.vcf` | 5 variants with TU/AB genotypes |
| `reporter_cassette.fa` | 725 bp mock knock-in cassette |

Key coordinates in `kcnh6a_locus.fa` (1-based):

- exon 1: 601-810
- exon 2: 1711-2040
- exon 3: 3141-3400
- protospacer: 1751-1770 (PAM 1771-1773)
- SpCas9 cut site: 1768
- conserved block shared with kcnh6b: transcript 311-430

Transcript junctions: 210 and 540.

Nothing here is real zebrafish sequence. The `a`/`b` names are placeholders
chosen because the ZFIN suffix convention is what the ohnolog logic keys on.
