# Changelog

All notable changes to this project are documented here, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Removed — the cut site is no longer a design criterion

- `genotype` no longer takes `--protospacer`, `--cut-site`, `--nuclease`,
  `--no-pam`, `--protected`, `--min-flank` or `--max-flank`, and the nuclease
  geometry table, `CutSite` and `locate_protospacer` are gone with them.
  **Removed rather than deprecated**, so a stale command line fails loudly
  instead of looking as though it still constrains the design.
- The criterion is now `--target`: the region to interrogate. The amplicon must
  span it and no primer may sit inside it. Give none and the design runs across
  the whole template.
- Why: a predicted cut coordinate is not where the lesion is. NHEJ gives
  deletions of tens to hundreds of bases, and an F0 mosaic carries a spectrum
  at once. Sizing a protected window from "the largest indel you expect" holds
  until the allele that matters is bigger than the guess — then a primer sits
  inside the deletion, the mutant allele does not amplify, and the animal reads
  as wild type. The old geometry optimised confidently for a coordinate whose
  relevance it could not check, and failed silently in the costly direction.
- `predict_allele_sizes` and the gel-resolvability warning stay: they concern
  the readout, need only the product size, and were always the sound half.

### Fixed — enumeration is aimed at the target

- A valid pair spans the target and is at most `--max-product-size` long, so
  candidate enumeration is now restricted to the windows that implies. This
  removes no possible result, and it matters because enumeration is capped and
  keeps what it finds first: unaimed, it spent its whole budget near the 5' end
  of a long template and then reported that no pair satisfied the constraints.
  The old cut-site windows had been hiding this.
- A region wider than the product-size window can accommodate now raises an
  error naming both constraints, instead of failing as "no pair found". A cut
  coordinate always fitted inside any product-size window, so this conflict —
  a short fragment-analysis product cannot span a wide exon — was invisible.

### Added — a real zebrafish locus, checked in

- `examples/myo7aa_locus.fa` is **real sequence**: 1,000 bp of *Danio rerio*
  chromosome 18 on danRer11, `chr18:3,405,001-3,406,000`, inside *myo7aa* — the
  *mariner* gene. `examples/myo7aa.gtf` carries the real RefSeq exon in it
  (NM_152983, exon 24 of 48). Retrieved from the UCSC Genome Browser API; the
  provenance is in the FASTA header.
- It earns its place by showing two things synthetic sequence cannot: this
  window is 40% GC, because the genome is AT-rich, and **41% of it is
  repeat-masked** — the most common reason a zebrafish primer that looks
  perfect does not work.
- `tests/test_real_locus.py` checks the sequence *against its annotation*
  rather than against a stored copy. *myo7aa* is on the minus strand, so each
  intron's canonical GT..AG reads as AC and CT on the plus strand, and those
  dinucleotides must land exactly where the GTF puts the boundaries. A dropped
  base, a shifted window, an edited coordinate, or synthetic sequence pasted in
  behind the same header all fail it.
- The browser demo's **design** tab now loads this locus, with the target over
  the real exon. The other tabs stay on synthetic data, which they need for a
  planted protospacer, a variant, or a spliced transcript; every synthetic
  header says SYNTHETIC and a test asserts it.
- It is a 1 kb window, not the whole 96 kb gene, so it has one exon and no
  spliced transcript. For those, `fetch`.

### Added — real sequence, by gene name

- **`danioprime fetch`** retrieves a locus or a spliced transcript from Ensembl
  (or by coordinates from UCSC), so a design starts from a gene name rather
  than from hunting for a FASTA:

  ```bash
  danioprime fetch myo7aa --flank 2000 --out myo7aa_locus.fa
  danioprime fetch myo7aa --transcript --out myo7aa_cdna.fa
  ```

- The FASTA header records **the assembly the server reported** — never one
  assumed here — with the coordinates, the source URL and the date. That is
  what lets the variant guard refuse a set on another assembly instead of
  misapplying it, and what makes a primer ordered six months ago retraceable.
- Every response is validated before it can become a file: exact span length,
  nucleotides only. A truncated response or an HTML error page raises rather
  than producing a locus that looks fine and is the wrong piece of DNA.
- `--dry-run` prints the URLs it would request and fetches nothing, which also
  names the hosts to allow through a restrictive proxy. `--base-url` points at
  an Ensembl mirror or a local REST install. Chromosome naming is translated
  per service (`chr18` ↔ `18`, `chrM` ↔ `MT`) and so is UCSC's 0-based
  half-open coordinate system.
- **`examples/myo7aa.sh`** runs the whole workflow on the real gene: locus,
  transcript, ohnolog panel, genotyping, in situ probe, qPCR. *myo7aa* is the
  *mariner* gene — the orthologue of human *MYO7A*, and the first fish model of
  human hereditary deafness (Ernest et al. 2000).

### Changed — the example data no longer impersonates real genes

- The synthetic dataset was named after `kcnh6a`/`kcnh6b`, which are real
  zebrafish genes. Invented bases in a file named after a real gene is an
  impersonation that survives until someone BLASTs it, so the pair is now
  `demogenea`/`demogeneb` and every header says SYNTHETIC. The `a`/`b` suffix
  stays, because that is what the ohnolog logic keys on.
- The dataset keeps its job — offline, deterministic tests and a browser demo
  that needs no network — and the README now says plainly that it is not a
  model of what your input looks like.
- The browser demo says the same, and shows the `fetch` command for real work.
  The page cannot fetch for itself: it is the real package compiled to
  WebAssembly, which has no sockets. That is also why nothing pasted there
  leaves the machine.

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
