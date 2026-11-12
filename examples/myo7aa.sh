#!/usr/bin/env bash
# A real gene, end to end: fetch myo7aa, then design against it.
#
# Unlike demo.sh, this script needs network access, because it starts by
# fetching sequence from Ensembl. That is the point: the rest of the example
# dataset in this directory is synthetic, which is right for tests and wrong
# for showing you what a design on your own gene looks like.
#
# myo7aa is the mariner gene -- zebrafish myosin VIIAa, the orthologue of human
# MYO7A, where mutations cause Usher syndrome type 1B. Ernest et al. (2000)
# showed mariner mutants have defective hair cells, making it the first fish
# model of human hereditary deafness. It is also an ohnolog, so it exercises
# the paralog logic, and the CRISPR screen whose M13-tailed fluorescent-PCR
# genotyping this tool implements -- Varshney et al. 2015, Genome Res
# 25:1030-42 -- is among the publications ZFIN lists for it.
#
# Usage:  bash myo7aa.sh [OUTPUT_DIR]
set -euo pipefail

OUT="${1:-myo7aa_out}"
mkdir -p "$OUT"
run() { printf '\n\033[1m$ %s\033[0m\n' "$*"; "$@"; }

# Set DANIOPRIME_BASE_URL to point at an Ensembl mirror, a local REST install,
# or the stub this repository's tests use. Unset, the public service is used.
BASE=()
[ -n "${DANIOPRIME_BASE_URL:-}" ] && BASE=(--base-url "$DANIOPRIME_BASE_URL")

# What would be requested, before anything is requested. Worth running first on
# a restricted network: it names the host to allow.
run danioprime fetch myo7aa "${BASE[@]}" --dry-run

# The genomic locus, plus 2 kb on each side so primers have somewhere to sit.
# The FASTA header records the assembly Ensembl served, which is what the
# variant guard reads later -- published zebrafish variant sets are on
# danRer11/GRCz11 and ZFIN has moved to GRCz12tu, so this is not a detail.
run danioprime fetch myo7aa "${BASE[@]}" --flank 2000 --out "$OUT/myo7aa_locus.fa"

# The spliced transcript, for assays that must not amplify genomic DNA.
run danioprime fetch myo7aa "${BASE[@]}" --transcript --out "$OUT/myo7aa_cdna.fa"

# The ohnolog, as the panel to discriminate against. If Ensembl has no myo7ab
# for your release, drop the --avoid flags below; the designs still run.
if danioprime fetch myo7ab "${BASE[@]}" --flank 2000 --out "$OUT/myo7ab_locus.fa" 2>/dev/null; then
  AVOID_LOCUS=(--avoid "$OUT/myo7ab_locus.fa")
else
  echo "note: no myo7ab in this Ensembl release; running without a paralog panel"
  AVOID_LOCUS=()
fi

# The header is ">myo7aa <assembly> <chrom>:<start>-<end> | ...", so the
# assembly is the second field. Reported rather than assumed, because this is
# the fact that decides whether a variant set may be applied at all.
ASSEMBLY=$(head -1 "$OUT/myo7aa_locus.fa" | awk '{print $2}')
LENGTH=$(awk '!/^>/ {n += length($0)} END {print n}' \
  "$OUT/myo7aa_locus.fa")
MIDPOINT=$((LENGTH / 2))
echo "note: ${LENGTH} bp on ${ASSEMBLY:-an unreported assembly}; targeting the midpoint"

# Genotyping primers for a CRISPR allele. The region to interrogate is the
# thing you know -- the exon you targeted -- not a predicted cut coordinate:
# a repaired lesion is not confined to the cut, so primers placed a fixed
# distance from it can land inside a larger deletion and fail to amplify the
# mutant allele at all. Give --target the exon's span in template coordinates.
run danioprime genotype \
  --template "$OUT/myo7aa_locus.fa" \
  --preset genotyping-gel \
  --target "$((MIDPOINT - 100))-$((MIDPOINT + 100))" \
  --tails fluorescent \
  "${AVOID_LOCUS[@]}" \
  --order-out "$OUT/myo7aa_order_sheet.tsv" \
  --show 2

# An in situ probe on the transcript, with T7 on the reverse primer so the
# product transcribes as the antisense strand.
run danioprime insitu \
  --template "$OUT/myo7aa_cdna.fa" \
  --scheme riboprobe-t7 \
  --show 2

# qPCR across the transcript.
# qPCR across the transcript. Pass --gtf with exon coordinates to get the
# junction-spanning check: `fetch` does not retrieve exon structure yet, so
# without one every pair is reported as offering no gDNA protection -- which is
# true, not a bug.
run danioprime qpcr \
  --transcript "$OUT/myo7aa_cdna.fa" \
  --show 2

echo
echo "Wrote: $OUT/"
echo "Order from the 'sequence' column of $OUT/myo7aa_order_sheet.tsv (upper case);"
echo "the lower-case bases in the report are the 5' tails."
