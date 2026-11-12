# Example dataset (synthetic)

**None of this is real zebrafish sequence.** The names are deliberately not
those of real genes: a file called `kcnh6a_locus.fa` holding invented bases is
an impersonation, and someone will eventually BLAST it. The `a`/`b` suffix is
kept because it is what the ohnolog logic keys on.

This dataset exists so the tests, `demo.sh` and the browser demo run offline and
deterministically, and it is built to contain the cases that are hard rather
than a happy path. To design against a real gene, fetch one:

```bash
danioprime fetch myo7aa --flank 2000 --out myo7aa_locus.fa
bash myo7aa.sh          # the whole workflow, on the real locus
```

Regenerate with `python make_example_data.py`. Deterministic: seed 20260211.

| file | contents |
|---|---|
| `demogenea_locus.fa` | 4320 bp synthetic locus, 3 exons, SpCas9 protospacer in exon 2 |
| `demogeneb_locus.fa` | 3940 bp synthetic ohnolog, exons ~94% identical to demogenea |
| `demogenea_transcript.fa` | 1120 bp spliced transcript |
| `demogeneb_transcript.fa` | 1120 bp ohnolog transcript (avoid-panel for qPCR) |
| `demogenea.gtf` | exon structure in `demogenea_locus.fa` coordinates |
| `strain_variants.vcf` | 5 variants with TU/AB genotypes |
| `reporter_cassette.fa` | 725 bp synthetic knock-in cassette |

Key coordinates in `demogenea_locus.fa` (1-based):

- exon 1: 601-810
- exon 2: 1711-2040
- exon 3: 3141-3720
- protospacer: 1751-1770 (PAM 1771-1773)
- SpCas9 cut site: 1768
- conserved block shared with demogeneb: transcript 311-430

Transcript junctions: 210 and 540.

Coding sequence ends at transcript position 800; the 3' UTR
(801-1120) is unrelated between the two genes. A probe over the
CDS cannot discriminate the ohnologs; a probe over the UTR can. That is the
contrast `danioprime insitu` is there to surface.

Nothing here is real zebrafish sequence. The `a`/`b` names are placeholders
chosen because the ZFIN suffix convention is what the ohnolog logic keys on.
