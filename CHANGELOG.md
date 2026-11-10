# Changelog

All notable changes to this project are documented here, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — in situ hybridisation

- `insitu` mode (CLI, Python API and browser) for riboprobe templates.
- **T3 is now reachable.** It was registered as a tail but had no scheme, so no
  interface could offer it. All of T7, T3 and SP6 now have schemes in both
  orientations, plus two dual-promoter schemes that put different promoters at
  each end so one PCR product templates both the antisense probe and its sense
  control.
- **Probe-level cross-hybridisation analysis**, which is a different question
  from primer specificity: priming is decided at a 3' end, hybridisation along
  the whole molecule. A pair whose primers discriminate the ohnolog perfectly
  can still yield a probe that stains both genes, and only this check sees it.
  Thresholds (85% identity over at least 100 nt) are calibrated against chance
  windows rather than guessed; see `docs/algorithm.md`.
- The example transcripts now carry independently generated 3' UTRs, as real
  ohnolog UTRs largely are, so the dataset demonstrates both the problem (a CDS
  probe cross-hybridises) and the fix (a UTR probe does not).

### Fixed

- **SP6 promoter was two bases short.** It was the 18-mer universal primer
  (`ATTTAGGTGACACTATAG`) rather than the in vitro transcription consensus
  (`ATTTAGGTGACACTATAGAA`). T7 and T3 were already correct. All three are now
  checked against the published consensus by a test, because being wrong by one
  base means wasted oligos and a failed transcription.
- Two defects in the conserved-block finder, both found by testing rather than
  by review: it reported the unchecked span between two chance seeds on a
  diagonal (600 nt at 64% identity was presented as one block), and its
  longest-window scan was in floating point, where prefix sums that are
  nominally equal differ by an ulp and the longest valid window is silently
  lost. Now integer arithmetic, verified against brute force on 2,000 random
  cases.

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
