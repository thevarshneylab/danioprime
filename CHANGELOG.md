# Changelog

All notable changes to this project are documented here, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — initial release

First working version. Nothing in it has been validated against wet-lab
outcomes yet; see the Validation section of the README.

### Design engine

- Nearest-neighbour thermodynamics (SantaLucia & Hicks 2004) with salt,
  divalent, dNTP-chelation and co-solvent corrections. Pure Python, no required
  dependencies, byte-for-byte reproducible.
- Optional `primer3-py` backend for numbers matching a Primer3 run exactly.
- Hairpin, self-dimer and hetero-dimer detection reporting both a classic
  complementarity score and a stack free energy.
- Candidate enumeration with cost-staged filters, interval-based target/exclude
  region logic, and pair diversification by 3'-end separation.
- Additive penalty model kept as labelled terms, so a ranking can be explained.

### Zebrafish-specific

- Ohnolog discrimination: avoid-panel scanning that judges sites on 3'-window
  mismatches, total mismatches and imperfect-duplex Tm rather than identity,
  with predicted cross-amplicons including self-pair products.
- Strain-variant masking from VCF or BED, with per-sample genotype filtering and
  `--heterozygous-only` to recover sites fixed non-reference in a line.
- `guess_ohnolog_symbol` hint from the ZFIN `a`/`b` convention, flagged as a
  prompt to confirm rather than an orthology call.

### Assay modes

- `genotype` — protospacer location for SpCas9/SaCas9/Cas12a on either strand,
  protected-window geometry, predicted allele sizes with gel-resolvability flags.
- `qpcr` — junction-spanning primers as a hard pool constraint, with a graded
  genomic-DNA protection verdict per pair.
- `dcaps` — engineered allele-specific restriction sites over 60 enzymes, with
  predicted digests and diagnostic bands separated from shared bands.
- `junction` — knock-in verification PCRs anchored outside the homology arm.
- `design` — general-purpose mode. `check` — evaluate existing primers.

### Interface

- `text`, `tsv`, `json` and `bed` output; oligo order sheets with tail schemes
  (M13F, pig-tail, T7/T3/SP6) and companion oligos listed once.
- Six presets covering gel and fluorescent genotyping, qPCR, Sanger,
  amplicon-seq and riboprobe templates.
- Indexed FASTA random access with automatic `.fai` creation; GTF/GFF3 transcript
  parsing.
- Optional `blastn` / `bowtie2` backends for genome-scale specificity, re-judged
  with the same priming rules as the built-in scanner.
