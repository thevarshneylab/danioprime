# Changelog

All notable changes to this project are documented here, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — a tail you can see

- A 5' tail is now written in **lower case** wherever a human reads it, and
  drawn in a **colour of its own**: M13F blue, M13R cyan, the pig-tail slate,
  and T7, T3 and SP6 in green, orange and magenta. The promoters take the three
  most separable colours on purpose, since transcribing from the wrong one
  gives a probe that cannot hybridise.
- The palette lives with the tail definitions, so a tail is the same colour in
  the terminal, in the browser and in `danioprime tails`, which now prints it.
- Case is the primary channel and colour the secondary one: everything the
  colour says, the case says too, so a redirected report, a monochrome terminal
  and a colour-blind reader lose nothing.
- `--color auto|always|never` (`--colour` also accepted), honouring `NO_COLOR`,
  `FORCE_COLOR` and `TERM=dumb`. `auto` colours only a terminal, so a report
  written with `--out` and the report string the browser embeds never contain
  escape sequences. Windows consoles are switched into VT mode when asked.
- A tailed primer now reports both lengths — `38 nt, 20 anneal` — because the
  oligo ordered is longer than the primer that anneals, and the Tm, GC and
  position given are properties of the annealing portion only.
- `Primer.annotated_seq` alongside `ordered_seq`; `annotated_seq.upper()` is
  `ordered_seq`, so nothing here can change a base.

### Changed

- The pair table's `forward_seq` and `reverse_seq`, and the sequences in the
  browser's pair cards, now carry the tail in lower case. The bases are
  unchanged; only the case differs.
- The order sheet gains a `sequence_annotated` column. `sequence` stays upper
  case throughout and remains the column to order from: it is what gets pasted
  into a vendor's form, and some bulk-upload parsers give lower case a meaning
  of their own.

### Added — published variant sets, queried by region

- `--variant-source nhgri1` reads the **NHGRI-1** zebrafish variant set
  (Burgess lab, NHGRI; LaFave et al. 2014) over the network, and
  `--variants-url` accepts any bgzip-compressed, tabix-indexed VCF or BED.
  `danioprime sources` lists what is available with citations and the
  `ZebrafishGenomics` UCSC track hub URL.
- Pure-Python **BGZF and tabix** readers, so a 214 MB dataset can be queried
  for one locus at a cost of about 66 kB that does not grow with the file. No
  compiled dependency, which keeps the browser build working. Validated against
  htslib: identical answers to pysam over hundreds of random regions, including
  indels at region boundaries, locally and over HTTP range requests.
- Each chunk is fetched in a single ranged read rather than one request per
  block, with a bounded one-block tail over-read. Tabix indexes are cached on
  disk between runs.
- An **assembly guard**. Published zebrafish variant sets are on danRer11
  (GRCz11) while ZFIN has moved to GRCz12tu, and the coordinates differ by tens
  to hundreds of kilobases. A declared mismatch is *refused*, because masking
  the wrong bases while reporting that known variants were avoided is worse
  than not checking; `--allow-assembly-mismatch` exists for coordinates already
  lifted over. An ad-hoc URL has no readable assembly, so it is reported as
  unverified rather than echoing back the declaration.
- `--strain` and `--heterozygous-only` apply to remote sources as they already
  did to local files, so a position fixed non-reference in the line is not
  masked.

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
