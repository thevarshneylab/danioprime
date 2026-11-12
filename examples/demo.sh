#!/usr/bin/env bash
# Every danioprime mode, run against the example dataset.
#
# The headline comparison is the first two blocks: the same genotyping design
# with and without the ohnolog panel. Without it, the reverse primer lands in a
# conserved exon block and matches the paralog perfectly.
set -euo pipefail
cd "$(dirname "$0")"

PROTOSPACER=AAGCCCTGTTATGTGCAACG
hr() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }

hr "1. Genotyping, with no knowledge of the ohnolog"
danioprime genotype \
    --template demogenea_locus.fa \
    --protospacer "$PROTOSPACER" \
    --preset genotyping-gel \
    --show 1

hr "2. The same design, now avoiding the ohnolog"
danioprime genotype \
    --template demogenea_locus.fa \
    --protospacer "$PROTOSPACER" \
    --avoid demogeneb_locus.fa \
    --preset genotyping-gel \
    --show 1

hr "3. Why pair 1 above was a trap: check its reverse primer against the ohnolog"
danioprime check \
    --primers "R_from_design_1=CTATCTGTATGCCGCGATG" \
    --avoid demogeneb_locus.fa

hr "4. Masking variants that segregate in the TU line"
danioprime genotype \
    --template demogenea_locus.fa \
    --protospacer "$PROTOSPACER" \
    --variants strain_variants.vcf \
    --variant-chrom demogenea_locus \
    --strain TU --heterozygous-only \
    --reject-variant-overlap \
    --preset genotyping-gel \
    --show 1

hr "5. Fluorescent-PCR genotyping, with the order sheet"
danioprime genotype \
    --template demogenea_locus.fa \
    --protospacer "$PROTOSPACER" \
    --avoid demogeneb_locus.fa \
    --preset genotyping-fluorescent \
    --tails fluorescent \
    --show 1

hr "6. RT-qPCR: junction-spanning, and not measuring the ohnolog as well"
danioprime qpcr \
    --transcript demogenea_transcript.fa \
    --gtf demogenea.gtf --transcript-id demogenea \
    --avoid demogeneb_transcript.fa \
    --require-junction-spanning \
    --show 2

hr "7. dCAPS for a point mutation a gel cannot size"
danioprime dcaps \
    --template demogenea_locus.fa \
    --variant 1768 --variant-base C \
    --preset genotyping-gel \
    --min-product-size 150 --max-product-size 450 \
    --max-assays 2

hr "8. In situ probe: dual promoters, and a probe-level specificity check"
# Note what this shows: the primers discriminate the ohnolog perfectly, yet the
# PROBE made from them would stain both genes. Priming is decided at a 3' end;
# hybridisation is decided along the whole molecule.
danioprime insitu \
    --template demogenea_transcript.fa \
    --avoid demogeneb_transcript.fa \
    --scheme riboprobe-dual-t7-sp6 \
    --min-product-size 400 --opt-product-size 500 --max-product-size 700 \
    --show 1

hr "9. The same gene, probed over the divergent 3' UTR instead"
# Ohnolog coding sequence is usually too conserved for any probe to
# discriminate. The 3' UTR diverges far faster, and is where a specific probe
# can actually be placed. CDS ends at transcript position 800.
danioprime insitu \
    --template demogenea_transcript.fa \
    --avoid demogeneb_transcript.fa \
    --scheme riboprobe-t3 \
    --target 850-1050 \
    --min-product-size 250 --opt-product-size 290 --max-product-size 320 \
    --show 1

hr "10. Junction PCRs verifying a knock-in"
danioprime junction \
    --template demogenea_locus.fa \
    --cassette reporter_cassette.fa \
    --insert-at 1768 \
    --homology-arm 24 \
    --preset genotyping-gel \
    --show 1

hr "11. Machine-readable output"
danioprime genotype \
    --template demogenea_locus.fa \
    --protospacer "$PROTOSPACER" \
    --avoid demogeneb_locus.fa \
    --preset genotyping-gel \
    --format json --out /tmp/danioprime_demo.json
python3 - <<'PY'
import json
d = json.load(open("/tmp/danioprime_demo.json"))
best = d["pairs"][0]
print(f"best pair: {best['name']}  {best['product_size']} bp")
print(f"  forward: {best['forward']['seq']}  Tm {best['forward']['tm']}")
print(f"  reverse: {best['reverse']['seq']}  Tm {best['reverse']['tm']}")
print(f"  why this rank: {best['penalty_terms']}")
print(f"  cut site: {d['cut_site']['description']}")
PY

hr "9. Real zebrafish sequence: 1 kb inside myo7aa on danRer11"
# Everything above ran on synthetic sequence. This block does not: it is
# chr18:3,405,001-3,406,000, and the target is the real RefSeq exon in it.
# Note the repeat-masked (lower case) stretches in the template -- 41% of this
# window -- which is the most common reason a zebrafish primer fails.
danioprime design \
    --template myo7aa_locus.fa \
    --target 679 \
    --preset genotyping-gel \
    --tails fluorescent \
    --show 2

printf '\n\033[1mDone.\033[0m See docs/zebrafish-notes.md for the reasoning behind the presets.\n'
printf 'Blocks 1-8 ran on synthetic sequence; block 9 ran on 1 kb of real chr18.\n'
printf 'For the whole workflow on the full myo7aa gene and its cDNA, run\n'
printf 'bash myo7aa.sh -- it fetches them from Ensembl, so it needs network access.\n'
