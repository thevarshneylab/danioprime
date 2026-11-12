# danioprime

[![CI](https://github.com/thevarshneylab/danioprime/actions/workflows/ci.yml/badge.svg)](https://github.com/thevarshneylab/danioprime/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org)
[![Licence](https://img.shields.io/badge/licence-MIT-green)](LICENSE)

Zebrafish-aware PCR primer design: **ohnolog discrimination**, **strain-variant masking**, and the **assay geometry** that is normally done by hand.

Primer3 solves the thermodynamics problem well, and `danioprime` does not try to replace it. What a generic designer cannot know is that the gene you are working on probably has a near-identical paralog, that your wild-type line is not isogenic, and that a genotyping amplicon has to span the lesion without touching it. Those are the three things this tool is for.

**[Use it in your browser →](https://thevarshneylab.github.io/danioprime/)** — no install, and nothing you paste is uploaded.

Or locally:

```bash
pip install danioprime
```

Pure Python, no required dependencies. Optional extras: `pip install "danioprime[primer3]"` to use libprimer3 for the thermodynamics instead of the built-in engine.

---

## Start with a real gene

`danioprime fetch` pulls the locus from Ensembl, so the first step is a gene name rather than hunting for a FASTA:

```bash
danioprime fetch myo7aa --flank 2000 --out myo7aa_locus.fa
danioprime fetch myo7aa --transcript --out myo7aa_cdna.fa
```

```
wrote myo7aa_locus.fa: myo7aa, <N> bp
  <assembly> <chr>:<start>-<end> | strand + | <N> bp | gene +/- 2,000 bp |
  ENSDARG00000099871 | myosin VIIAa | via Ensembl REST | retrieved 2026-11-12
```

(The assembly, contig and coordinates are whatever Ensembl serves for the release you fetch from, which is why they are not written out here.)

The header is the point. It records **the assembly the server reported** — never one assumed here — along with the coordinates, the source and the date, because a template whose assembly is unknown cannot be checked against a variant set, a GTF, or a primer you ordered six months ago. Pass it on and a variant set on another assembly is refused rather than misapplied:

```bash
danioprime genotype --template myo7aa_locus.fa --assembly GRCz11 \
    --protospacer <your guide> --preset genotyping-fluorescent \
    --tails fluorescent --variant-source nhgri1
```

`--dry-run` prints the URLs it would request and fetches nothing, which is the quickest way to find out what to allow through a restrictive proxy. `--source ucsc` fetches by coordinates from a named assembly instead, for when you need danRer11 specifically; `--base-url` points at a mirror or a local Ensembl REST install.

`examples/myo7aa.sh` runs the whole thing — locus, transcript, ohnolog panel, genotyping, in situ probe, qPCR — against the real gene. *myo7aa* is the **mariner** gene: zebrafish myosin VIIAa, orthologue of human *MYO7A*, where mutations cause Usher syndrome type 1B. Ernest et al. (2000) showed *mariner* mutants have defective sensory hair cells, making it the first fish model of human hereditary deafness. It is a reasonable example here because it has an ohnolog, it is large enough to be realistic, and the high-throughput CRISPR screen whose M13-tailed fluorescent-PCR genotyping this tool implements — Varshney et al. (2015), *Genome Res* 25:1030–42 — is among the publications ZFIN lists for it.

### One real locus, checked in

`examples/myo7aa_locus.fa` is **real zebrafish sequence**: 1,000 bp of chromosome 18 on danRer11, `chr18:3,405,001–3,406,000`, inside *myo7aa*. `examples/myo7aa.gtf` carries the real RefSeq exon that sits in it.

```bash
danioprime design --template examples/myo7aa_locus.fa --target 679 \
    --preset genotyping-gel --tails fluorescent
```

It is worth having checked in because it shows two things synthetic sequence cannot. It is 40% GC, because this genome is AT-rich — so the Tm window a human-tuned tool would give you is the wrong one. And **41% of it is repeat-masked**, in lower case, which is the single most common reason a zebrafish primer that looks perfect does not work.

The provenance is in the header, and the tests check the sequence *against the annotation* rather than against a stored copy of itself: *myo7aa* is on the minus strand, so each intron's canonical GT…AG reads as AC and CT on the plus strand, and those dinucleotides have to land exactly where the GTF says the boundaries are. Drop a base, shift the window, edit the coordinates, or paste in a synthetic sequence behind the same header, and that check fails.

It is a 1 kb window rather than the whole 96 kb gene, so it suits `design`, `genotype` and `dcaps` but has only one exon — not enough for a spliced transcript. For the full locus or the cDNA, `fetch` it.

> **Everything else in `examples/` is synthetic**, and named so that it cannot be mistaken for a real gene (`demogenea`, `demogeneb`). It exists so the tests and the browser demo run offline and deterministically, and it is built to contain the hard cases — a 94%-identical ohnolog, a conserved exon block, strain variants in an otherwise attractive primer site. Every one of those headers says SYNTHETIC, and a test asserts it.

---

## Why a zebrafish-specific tool

### 1. The teleost genome duplication

Zebrafish retain thousands of ohnolog pairs from the teleost-specific whole-genome duplication, conventionally named `geneA`/`geneB`. Coding-sequence identity between them is often high enough that a pair designed against one gene model amplifies both — and because the paralog is a real gene rather than a repeat, nothing in a standard workflow flags it.

`danioprime` asks the question directly:

```bash
danioprime genotype --template demogenea_locus.fa \
                    --protospacer AAGCCCTGTTATGTGCAACG \
                    --avoid demogeneb_locus.fa
```

```
  demogenea_locus_gt_p1  314 bp  penalty 1.93
    F  GTGGCTCATAACCTGTGATAG
       21 nt  Tm 59.9C  GC 48%  pos 1674-1694  hairpin -2.4
    R  GAGAGGCGATGGTGTATTATG
       21 nt  Tm 60.0C  GC 48%  pos 1967-1987  hairpin -0.7
       - discriminates demogeneb_locus by a mismatch 1 nt from the 3' end
    - no product predicted on the avoid-panel
```

Without `--avoid`, the same run picks a reverse primer at 1877-1895 — a **perfect** match to the ohnolog, because it sits in a conserved exon block. The tool rewards the trick experimentalists use by hand: putting the primer's 3'-terminal bases over a position where the two ohnologs differ, so discrimination comes from extension failure rather than from a Tm difference that a permissive annealing temperature washes out.

The `--avoid` panel is not ohnolog-specific. Any "must not amplify" sequence works: a family member, a pseudogene, a transgene backbone, the unedited allele.

### 2. Priming is not sequence identity

The question that matters is not *"does this primer have a similar sequence elsewhere"* but *"can a polymerase extend from it elsewhere"*. Those differ, and conflating them is why BLAST output is a poor proxy for specificity:

- A primer with 6 internal mismatches and a perfect 3' end **will** prime.
- A primer with one mismatch at its 3'-terminal base **will not**, however well the rest matches, because the polymerase needs a correctly paired 3'-OH.

Every candidate site is judged on three criteria, reported separately: mismatches in the 3'-terminal window, total mismatches, and the estimated Tm of the imperfect duplex.

### 3. Strain polymorphism

TU, AB, TL and WIK are not isogenic. A primer sitting on a SNP that segregates in your line fails on one allele and looks exactly like a genotype — the worst kind of failure, because it is silent.

```bash
danioprime genotype --template demogenea_locus.fa --cut-site 1768 \
                    --variants sanger_zfish.vcf.gz --variant-chrom chr5 \
                    --template-offset 12400000 \
                    --strain TU --heterozygous-only
```

A variant anywhere in the footprint is penalised; a variant within the 3'-terminal window is **always** disqualifying. `--heterozygous-only` is the detail that matters in practice: a position fixed homozygous-alternate in your line is *not* a hazard — the primer just has to match the line, not the reference — so excluding those recovers primer sites that naive masking throws away.

#### Published variant sets, queried by region

You do not need a local VCF. `danioprime sources` lists published resources that can be read over the network, and only the locus is transferred:

```bash
danioprime genotype \
    --region chr5:1,000,000-1,006,000 --reference danRer11.fa \
    --assembly danRer11 \
    --variant-source nhgri1 --strain NHGRI1 --heterozygous-only
```

```
note: 1,284 variant(s) from NHGRI-1 (danRer11) in chr5:1000001-1006000;
      66 kB fetched in 1 request(s)
```

`nhgri1` is the **NHGRI-1** line from the Burgess lab at NHGRI ([LaFave et al. 2014](https://doi.org/10.1534/genetics.114.166769)), the strain built so that "the exact sequence (with all variants) of a particular genomic region" is known. The VCF is 214 MB; using its tabix index, a locus costs about **66 kB of transfer regardless of dataset size**, so there is nothing to download and nothing to keep in sync. The index itself is cached after first use.

`--variants-url` takes any bgzip-compressed, tabix-indexed VCF or BED, local or remote, so a resource missing from the registry is never a blocker. The same data is also browsable as the `ZebrafishGenomics` UCSC track hub, which `danioprime sources` prints.

#### The assembly trap, and why this refuses rather than warns

Every published zebrafish variant set is on **danRer11 / GRCz11**. ZFIN has moved to GRCz12tu, and the coordinates are not transferable — shifts of tens to hundreds of kilobases are normal.

Masking a GRCz12tu template with GRCz11 coordinates would exclude the wrong bases *while reporting that known polymorphisms were avoided*. That is worse than not checking at all, so `danioprime` **refuses** the combination:

```
danioprime: NHGRI-1 is on danRer11, but the template coordinates were declared
as GRCz12tu. These assemblies' coordinates are not transferable [...] so masking
would exclude the wrong bases while reporting that known variants were avoided.
Either supply a danRer11 template, lift your coordinates over to danRer11
first, or pass --allow-assembly-mismatch if you have already done so.
```

An ad-hoc `--variants-url` carries no assembly that can be read, so its assembly is reported as `unknown (unverified)` with a warning, rather than echoing back whatever you declared as though it had been checked.

### 4. Assay geometry

| Mode | What it enforces |
|---|---|
| `genotype` | Amplicon spans the lesion; **neither primer touches it**; primers far enough out for the readout; predicted allele sizes with a flag for indels a gel cannot resolve |
| `qpcr` | A junction-spanning primer where possible, with an explicit genomic-DNA protection verdict per pair |
| `dcaps` | Engineered mismatches that create an allele-specific restriction site, with predicted digests for both alleles |
| `junction` | One primer in genomic sequence outside the homology arm, one in the cassette — so the assay cannot read the donor plasmid instead of the integrated allele |
| `insitu` | Riboprobe template with a T7, T3 or SP6 promoter in either orientation, checked for **probe-level** cross-hybridisation |
| `design` | General-purpose target/exclude-region design |
| `check` | Evaluate primers you already have, including against an avoid-panel |

---

## Worked examples

The `examples/` directory ships a synthetic dataset built to have the properties that make zebrafish design hard — a mock ohnolog pair with a conserved exon block, a real NGG protospacer, and a VCF with TU/AB genotypes. `examples/demo.sh` runs every mode.

### Genotype a crispant

```bash
danioprime genotype \
    --template examples/demogenea_locus.fa \
    --protospacer AAGCCCTGTTATGTGCAACG \
    --avoid examples/demogeneb_locus.fa \
    --preset genotyping-fluorescent \
    --tails fluorescent \
    --protected 60 \
    --order-out oligos.tsv
```

`--protected 60` is the parameter people get wrong: it must exceed the **largest indel you expect**, not the one you got, because the assay has to work before the alleles are known. A primer inside a deleted span will not amplify the mutant allele at all.

`--tails fluorescent` emits the **three** oligos the Varshney fluorescent-PCR system actually needs — an M13F-tailed forward, a pig-tailed reverse, and one labelled M13F primer listed once for the whole sheet, since that oligo is shared across every assay in the lab.

A tail is printed in lower case, and in a colour of its own:

```
  demogenea_locus_gt_p1  327 bp  penalty 1.36
    F  tgtaaaacgacggccagtCGCTACCGACTCAAGTAAAG  (+M13F)
       38 nt, 20 anneal  Tm 60.1C  GC 50%  pos 1569-1588  hairpin -2.0
    R  gtttcttCTATCTGTATGCCGCGATG  (+pigtail)
       26 nt, 19 anneal  Tm 60.1C  GC 53%  pos 1877-1895  hairpin -1.8
```

Each tail keeps the same colour in the terminal, in the browser and in `danioprime tails`, which prints the palette. The colour makes a T3-tailed primer impossible to mistake for a T7 one at a glance; the lower case says the same thing without it, so a redirected report, a monochrome terminal or a reader who cannot tell the two greens apart loses nothing. `--color never` switches the colour off, `NO_COLOR` is honoured, and a report written with `--out` never contains escape sequences. Both numbers are given for a tailed primer because the oligo you order is longer than the primer that anneals, and the Tm, GC and position are properties of the part that anneals.

In the order sheet the `sequence` column is upper case throughout — that is the column to paste into a vendor's form — and `sequence_annotated` beside it carries the same bases with the tail in lower case, for checking by eye that the right tail went on the right primer.

### RT-qPCR that cannot read genomic DNA or the ohnolog

```bash
danioprime qpcr \
    --transcript examples/demogenea_transcript.fa \
    --gtf examples/demogenea.gtf --transcript-id demogenea \
    --avoid examples/demogeneb_transcript.fa \
    --require-junction-spanning
```

```
Genomic-DNA protection:
  pair                gDNA protection  mechanism
  ------------------  ---------------  ----------------------------
  demogenea-201_qpcr_p1  strong           junction-spanning primer (R)
```

Protection is graded rather than asserted. `strong` means a primer straddles a junction, so it has no contiguous binding site on genomic DNA at all. `moderate` means the amplicon spans an intron large enough that the genomic product fails under qPCR cycling. `weak` means it spans a junction whose intron is too small to protect anything — a 70 bp intron gives no protection. `none` means genomic DNA amplifies identically, and you are relying entirely on your DNase step.

If a junction-spanning primer is impossible under your constraints, the tool says so and names the filter that blocked it rather than silently returning unprotected pairs.

### dCAPS for a 1 bp indel or point mutation

```bash
danioprime dcaps --template examples/demogenea_locus.fa \
                 --variant 1768 --variant-base C \
                 --min-product-size 150 --max-product-size 450
```

```
  demogenea_locus_dcaps1  RsaI (GT^AC)  1 introduced mismatch(es)  product 320 bp
    dCAPS primer (F)  CCTAAAGCCCTGTTATGTGTA
    template was      CCTAAAGCCCTGTTATGTGCA
    mismatch 2 nt from the 3' end (offsets [1])
    digest: reference [184, 136] bp  variant [164, 136, 20] bp
    diagnostic bands: reference-only [184], variant-only [164, 20]
    - RsaI also cuts the amplicon 1 time(s) independently of the genotype;
      those bands appear in both alleles and serve as a digest control
```

Diagnostic bands are separated from bands common to both alleles, because only the former carry genotype information. The default keeps the 3'-terminal base matched, since a terminal mismatch cripples extension.

### An in situ probe that does not stain the ohnolog too

```bash
danioprime insitu \
    --template examples/demogenea_transcript.fa \
    --avoid examples/demogeneb_transcript.fa \
    --scheme riboprobe-dual-t7-sp6
```

Two things here that a primer-shaped analysis gets wrong.

**Promoter orientation decides which strand you get.** On the *reverse* primer a promoter transcribes the antisense strand — the probe that detects mRNA. On the *forward* primer it transcribes the sense strand, the negative control. Get it backwards and the probe cannot hybridise while the control can, and the experiment just looks "failed". T7, T3 and SP6 are all available in both orientations; a **dual** scheme puts different promoters at the two ends so one PCR product templates both probes and no second reaction is needed.

**Cross-hybridisation is not primer specificity.** Everywhere else in this package, specificity is decided at a primer's 3' end, because that is what a polymerase needs. A riboprobe has no 3' end that matters: it is hundreds of bases of RNA and will bind anything similar enough along its whole length. So this mode analyses the amplicon body, not its ends:

```
Probe cross-hybridisation (whole-amplicon, not primer 3' ends):
  pair                 cross-hybridises  longest shared block
  -------------------  ----------------  -------------------------
  demogenea-201_probe_p1  YES               496 nt @ 93% (demogeneb-201)
```

That pair's primers are *perfectly* specific — the reverse primer discriminates `demogeneb` by a 3'-terminal mismatch and no product is predicted on the panel. The probe made from it would nonetheless stain both genes. A tool that only checked primers would have called it clean, and the published expression pattern would belong to two genes.

The thresholds (85% identity over ≥100 nt) are calibrated rather than guessed — see [`docs/algorithm.md`](docs/algorithm.md#in-situ-probes).

### Check primers you already have

```bash
danioprime check --primers "F=CGCTACCGACTCAAGTAAAG,R=CTATCTGTATGCCGCGATG" \
                 --avoid examples/demogeneb_locus.fa
```

```
primer  closest match        primes?  verdict
------  -------------------  -------  -------------------------------------
F       none                 no       no detectable binding site on the panel
R       demogeneb_locus (0 mm)  YES      primes on demogeneb_locus (0 mismatches)
```

---

## Running it in a browser

<https://thevarshneylab.github.io/danioprime/> runs this package — not a reimplementation of it — as CPython 3.14 compiled to WebAssembly, via [Pyodide](https://pyodide.org). All five modes are available.

Two consequences worth knowing:

- **Nothing is uploaded.** There is no server. The page fetches a Python runtime and the `danioprime` wheel, then does every calculation in your tab. Unpublished sequence is as safe there as in a local terminal.
- **It works at all because the package has no dependencies.** A wheel is a zip, so `sys.path.insert(0, "danioprime.whl")` is the entire install step — no package manager, no compilation. That property was a design goal, not a lucky accident.

The first visit downloads ~12 MB of runtime, then caches it. A design takes a few seconds; the page is unresponsive while Python works, which is honest about what is happening rather than hiding it behind a spinner that implies otherwise.

What the browser build cannot do: genome-scale off-target screening, which needs `blastn` or `bowtie2`. In the browser the avoid-panel is whatever you paste. For a whole-genome check, use the CLI.

The site is published from `main` / `/docs`, which is a **build output**: edit `web/index.html`, then run `python web/build_site.py`. Publishing from `docs/` rather than the repository root keeps the served surface to four deliberate files instead of exposing the whole tree, and `docs/.nojekyll` makes Pages serve them verbatim rather than running the HTML through Jekyll.

CI verifies the browser build on every push — that the committed wheel matches `src/`, and that all five modes plus the error path still run under Pyodide. The wheel going stale is a real failure mode: it happened once during development, and the website would have silently served old code.

## Python API

```python
from danioprime import DesignParams, DesignRequest, ParalogPanel, read_single_fasta, run_design

result = run_design(DesignRequest(
    template=read_single_fasta("demogenea_locus.fa"),
    params=DesignParams().replace(min_product_size=200, max_product_size=600),
    panel=ParalogPanel.from_fasta("demogeneb_locus.fa"),
    name="demogenea",
))

for pair in result.pairs:
    print(pair.name, pair.product_size, pair.forward.seq, pair.reverse.seq)
    print("  why this rank:", pair.penalty_terms)
```

Every penalty is kept as a labelled term rather than collapsed to one number, so the report can say *why* a pair ranked where it did — which is what you need when the answer is to relax one constraint rather than guess.

---

## Specificity at genome scale

The built-in scanners are exact and dependency-free, and handle sequence sets up to a few megabases: paralog panels, transcriptomes, single loci, plasmids. A whole zebrafish genome is ~1.4 Gb, so genome-wide checks delegate to a real aligner when one is on `PATH`:

```bash
makeblastdb -in GRCz12tu.fa -dbtype nucl -out GRCz12tu
danioprime info   # reports which optional aligners were found
```

Both `blastn` (with Primer-BLAST's short-query parameters) and `bowtie2` are used purely as **candidate site finders**. Their hits are re-judged with the same priming rules as the built-in scanner, so a primer's verdict does not depend on which tool found the site.

Scanning strategies and their guarantees are documented in [`docs/algorithm.md`](docs/algorithm.md).

---

## What this tool does not do

Being clear about the edges:

- **It does not fetch sequence for you.** No Ensembl or ZFIN API calls. Supply FASTA, or a `--region` with an indexed `--reference`. This keeps runs reproducible and offline, but it does mean you fetch the locus yourself. (Variant *data* is the exception: `--variant-source` reads published sets over the network, since downloading 214 MB to read 6 kb would be absurd.)
- **Ohnolog pairing is your call.** `guess_ohnolog_symbol` infers `kcnh6a` → `kcnh6b` from the ZFIN suffix convention and the `qpcr` mode prints it as a hint, but it is a naming convention, not an orthology assertion. Plenty of genes end in `a`/`b` without being ohnologs. Confirm against ZFIN or an orthology database, then pass the sequence with `--avoid`.
- **The thermodynamic model is two-state nearest-neighbour**, not a partition function. Mismatched-duplex Tm is a coarse estimate adequate for ranking cross-reactivity, not a substitute for a full calculation.
- **No liftOver.** It will refuse a cross-assembly variant mask but will not convert coordinates for you; use UCSC chain files and spot-check the result.
- **The browser build cannot reach remote variant sets**, because Pyodide has no sockets. Paste a VCF there, or use the CLI.
- **No multiplex optimisation.** Pairs are scored individually; it will not solve a 24-plex panel for you.
- **The enzyme table is a curated common subset** (60 enzymes). Type IIS enzymes that cut at a distance are deliberately excluded, since dCAPS fragment prediction assumes the cut falls inside the recognition site. Supply your own with `--enzymes`.
- **It has not been validated against bench results at scale.** The algorithms are tested for correctness and the predictions are internally consistent and independently re-derived in the test suite, but "these primers worked in N of M loci" is the evidence that would make this publishable, and it does not exist yet. See [Validation](#validation).

---

## Validation

What is currently verified, by 287 tests covering 93% of the package:

- Nearest-neighbour enthalpy and entropy sums match hand calculations from the published SantaLucia & Hicks (2004) tables, term by term.
- Salt, divalent, dNTP-chelation, DMSO and primer-concentration effects move Tm in the right direction and magnitude.
- Seeded and exhaustive binding-site search return **identical** priming-site sets on planted targets.
- dCAPS digests are re-derived independently of the design code: the engineered amplicon is rebuilt, re-digested, and the cut allele is confirmed to gain exactly one site.
- Genotyping geometry is asserted directly: every reported amplicon spans the protected window and neither primer enters it.
- Aligner output parsers are tested against canned `blastn` tabular and SAM records, including the minus-strand register conventions that differ between the two formats — writing those tests found a real bug in which a perfect minus-strand BLAST hit was discarded as fully mismatched.
- The shipped example is a regression guard on the headline claim — without the ohnolog panel the design lands on sequence shared with the paralog; with it, the design moves.

What is **not** yet verified: agreement with wet-lab success rates. If you use this and record which primers worked, that data is the most valuable contribution you could make — see [`CONTRIBUTING.md`](CONTRIBUTING.md).

---

## Reproducibility

The default `builtin` thermodynamic backend is pure Python with no compiled dependencies, so results are byte-for-byte identical across platforms and CI runs. Every run can emit its complete parameter set:

```bash
danioprime design --template locus.fa --preset qpcr --format json --out run.json
python -c "import json;from danioprime import DesignParams;DesignParams.from_dict(json.load(open('run.json'))['parameters'])"
```

Switch to `--thermo-backend primer3` when you need numbers that match a Primer3 run exactly.

**The two backends are not interchangeable mid-project.** They agree on Tm to within a couple of degrees — CI asserts this on every push — but a few tenths of a degree is enough to reshuffle a ranking in which hundreds of candidates sit within one penalty unit of each other, so they routinely return *different* top pairs. Both pairs are valid and pass every filter; they are just not the same pair. Pick a backend at the start of a project and record it alongside the primers.

---

## Repository layout

```
src/danioprime/      the package
  thermo/            nearest-neighbour model, structure, backends
  design/            candidates, pairing, scoring, the shared pipeline
  specificity/       priming rules, scanners, ohnolog panels, aligners
  assays/            genotyping, qPCR, dCAPS, enzymes, tails
  annotation/        GTF/GFF transcripts
  variants/          VCF/BED masking
  report/            text, TSV, JSON, BED writers
  webapp.py          the single entry point the browser build calls
web/index.html       source of the browser app (self-contained)
web/build_site.py    assembles docs/ ; --check fails CI when it drifts
web/check_wheel.py   fails CI when the published wheel drifts from src/
web/smoke.mjs        runs every web mode under real Pyodide
docs/                published by GitHub Pages: the built site + the wheel
examples/            synthetic mock-ohnolog dataset + demo.sh
```

## Documentation

- [`docs/algorithm.md`](docs/algorithm.md) — the thermodynamic model, the priming rules, scanning strategies and their sensitivity guarantees, the penalty model
- [`docs/zebrafish-notes.md`](docs/zebrafish-notes.md) — ohnologs, strain variation, GRCz12tu coordinates, and the practical notes behind the presets
- [`examples/DATA.md`](examples/DATA.md) — what the example dataset contains and how to regenerate it
- [`web/`](web/) — the browser build's checks: `check_wheel.py` and `smoke.mjs`

## Citation

If this is useful in published work, cite the repository and the underlying methods it builds on: SantaLucia & Hicks (2004) for the thermodynamics, Neff et al. (1998) for dCAPS, Brownstein et al. (1996) for the pig-tail, and Varshney et al. (2015) for the fluorescent-PCR genotyping system.

## Licence

MIT — see [`LICENSE`](LICENSE).
