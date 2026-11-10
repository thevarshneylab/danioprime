# Algorithm and model reference

## 1. Thermodynamics

### Nearest-neighbour model

Watson-Crick stacking parameters are the unified set of SantaLucia & Hicks (2004), *Annu. Rev. Biophys. Biomol. Struct.* 33:415-440 — ΔH in kcal/mol, ΔS in cal/(mol·K). Helix initiation uses the terminal-base-pair terms (G·C: 0.1, −2.8; A·T: 2.3, 4.1) and self-complementary oligos take the −1.4 cal/(mol·K) symmetry correction.

```
Tm = ΔH·1000 / (ΔS_salt + R·ln C_T) − 273.15
```

with `R = 1.9872 cal/(mol·K)` and `C_T = [primer]` — primer in large excess over template, the usual PCR regime — quadrupled for self-complementary oligos as the two-state model requires.

### Salt correction

```
ΔS_salt = ΔS + 0.368·(N−1)·ln[Na⁺]_eq
```

Divalent cations and dNTPs fold into an equivalent monovalent concentration using the von Ahsen approximation that Primer3 uses:

```
[Na⁺]_eq = [Na⁺] + [K⁺] + [Tris]/2 + 120·√([Mg²⁺] − [dNTP])
```

The dNTP term matters more than people expect: at 1.5 mM Mg²⁺ and 0.6 mM dNTPs the free Mg²⁺ is 0.9 mM, and the chelated fraction is not available to stabilise the duplex. Co-solvent corrections are linear: −0.75 °C per % DMSO, −0.65 °C per molar formamide.

### Backend agreement, and why it does not mean identical output

CI runs the full suite with `primer3-py` installed and asserts the two backends'
Tm estimates agree within 3 °C for a typical 20-mer, and separately runs a real
design through each. They nonetheless report different top pairs: with hundreds
of candidates inside one penalty unit of each other, a few tenths of a degree
is enough to reorder them. Agreement on the *physics* does not imply agreement
on the *ranking*, and only the former is a correctness claim.

### A note on Tm monotonicity

Tm is **not** strictly monotonic in primer length, and this is correct rather than a defect. Appending a weak A·T pair to a GC-rich duplex pulls the ΔH/ΔS ratio towards the A·T asymptote, so the computed Tm can fall by a few tenths of a degree. The test suite asserts the trend over several bases rather than step-by-step monotonicity.

### Imperfect duplexes

Mismatched primer-template duplexes use internal single-mismatch stacks from the Allawi & SantaLucia series (1997-1998). Terminal mismatches contribute no stacking stability. The result is a **coarse estimate** adequate for ranking cross-reactivity, not a substitute for a partition-function calculation.

One consequence is worth stating plainly: because terminal mismatches are modelled as merely non-stabilising, a primer with two 3'-terminal mismatches can compute a *higher* Tm than one with a single disruptive internal mismatch. That is thermodynamically defensible and biologically misleading, which is exactly why 3'-end integrity is handled as a **separate, categorical** criterion rather than being folded into Tm.

### Secondary structure

Two measures per interaction:

- **`score`** — classic Primer3-style complementarity from the best ungapped antiparallel alignment, G:C weighted 2 and A:T weighted 1. Dimensionless, and what published rules of thumb are phrased in.
- **`dg`** — free energy of the longest contiguous base-paired stack in that alignment. Comparable across oligo lengths, and the better basis for a hard filter.

Hairpins search every `(i, j)` seed and extend inwards while bases pair and the loop stays at or above 3 nt. 3'-anchored interactions are reported separately, because only a paired 3' end can be extended.

---

## 2. Priming rules

A site is treated as able to prime when **all three** hold:

| Criterion | Default | Rationale |
|---|---|---|
| No mismatch within *k* nt of the 3' end | *k* = 3 | The polymerase needs a correctly paired 3'-OH. A mismatch here reduces extension by orders of magnitude, independent of Tm. |
| Total mismatches ≤ *m* | *m* = 5 | Beyond this the duplex does not form at annealing temperature. |
| Estimated duplex Tm ≥ *T* | *T* = 45 °C | Catches sites that satisfy both mismatch rules but are too unstable to anneal. |

All three are configurable via `PrimingRules`. The 3'-window criterion is the one that distinguishes this from an identity filter, and it cuts both ways: it **admits** sites with many internal mismatches and a clean 3' end (which do prime, and which a BLAST identity cut-off discards) and **rejects** near-perfect sites with a terminal mismatch (which do not prime, and which BLAST reports as excellent hits).

### Amplicon prediction

Predicted products pair up priming sites in the correct relative orientation within the size window. Every primer combination is considered **including a primer with itself**, because F/F and R/R products from inverted repeats or paralogous arrays are a real source of artefact bands that a "is my primer unique" check misses entirely.

---

## 3. Scanning strategies

| Strategy | Exactness | Throughput | Use for |
|---|---|---|---|
| `exhaustive` | Exact within the mismatch budget | ~300 kb/s/primer | Paralog panels, single loci, transcriptomes, plasmids (≤ ~3 Mb) |
| `seeded` | See contract below | ~4,000× faster per query after a one-off index build | Large panels, whole transcriptomes |
| `blastn` / `bowtie2` | Tool-dependent | Genome scale | Whole-genome off-target checks |

`auto` picks `exhaustive` below 3 Mb of panel sequence and `seeded` above it.

### Seeded sensitivity contract

The seed is the primer's **3'-terminal k-mer**, not an arbitrary one. A site that can be extended necessarily has a perfectly matched 3' end, so anchoring there loses nothing that could prime efficiently.

With seed length *k*, seeded search finds **every** site whose 3'-terminal *k* bases match exactly. It **misses** sites carrying a mismatch inside that terminal window but outside the "clean 3' end" window. Those sites are, by construction, poor primers — but if you need a guarantee, use `exhaustive` or an external aligner.

Why not seed on any k-mer with full pigeonhole sensitivity? For a 20-mer at *m* mismatches, guaranteed sensitivity needs `k ≤ ⌊L/(m+1)⌋` — that is `k = 6` at *m* = 2, which returns a hit every ~4 kb and makes the index useless. The 3'-anchored choice buys speed by discarding exactly the sites the biology says do not matter.

### External aligners

`blastn` runs with `-task blastn-short -word_size 7 -dust no -soft_masking false`, matching NCBI Primer-BLAST's short-query configuration, and requests the aligned subject sequence so mismatch positions can be recovered without a second reference pass. Partial alignments that stop short of the primer's 3' end are discarded, since they cannot extend.

`bowtie2` runs `--end-to-end --very-sensitive -N 1 -L 10 -k <n>`; mismatch positions are reconstructed from the MD tag, and gapped alignments are discarded because a gapped primer-template duplex is not a priming geometry this model covers.

In both cases hits are re-judged with the same `PrimingRules` as the built-in scanner.

---

## 4. Penalty model

An additive weighted sum, kept as **labelled terms** rather than a single float so the report can explain a ranking.

Primer-level terms: Tm deviation from optimum, GC deviation, length deviation, hairpin ΔG excess, self-dimer ΔG excess, 3'-end stability excess, missing GC clamp, poly-X excess, dinucleotide-repeat excess, ambiguous bases, soft-mask overlap, variant overlap, 3'-proximal variants, off-target sites, paralog priming sites, and a **negative** term (a bonus) for straddling an exon junction.

Pair-level terms: both primer penalties, pair Tm difference, product-size deviation, hetero-dimer ΔG excess, 3'-end dimer score excess, predicted paralog amplicons, predicted off-target amplicons.

Structure penalties apply only to the **excess beyond the threshold**, so a primer comfortably within limits is not penalised for having any structure at all.

### Hard filters vs penalties

Disqualifying (primer never reported): length, Tm window, GC window, GC clamp, poly-X, dinucleotide repeats, hairpin/self-dimer/end-stability limits, repeat-mask fraction, excluded regions, **any variant in the 3' window**, ambiguous bases (unless allowed).

Penalised but reported: internal variants, soft-masked overlap, panel priming sites, off-target sites, suboptimal everything else.

Filters are applied in **increasing order of cost** — composition, then Tm, then secondary structure — so the expensive hairpin and dimer searches only run on placements that already satisfy everything cheaper. On a 1.2 kb template this is about six times faster than evaluating everything up front, with identical results.

### Why a junction bonus is not enough

A scoring bonus can still be outranked and silently dropped by the pool cutoff. When `--require-junction-spanning` is set, junction spanning becomes a **hard constraint on the candidate pools**: both assignments (spanning forward with any reverse, any forward with spanning reverse) are enumerated separately and merged.

That exposes a second subtlety. Selecting the partner pool by absolute penalty biases it towards Tm ≈ optimum, which then cannot satisfy `max_pair_tm_difference` against an anchor whose own Tm is off-optimum — exactly what a forced junction-spanning primer creates. Partners are therefore filtered for **compatibility with the anchor** before truncation.

### Diversification

Reported pairs must have 3' ends at least `pair_diversity_window` (default 8 nt) apart. Without this, the top 20 pairs are typically the same primer shifted by one base, which gives no real alternative when the first fails at the bench.

---

## 5. Assay geometry

### Genotyping

Three constraints must hold simultaneously:

1. the amplicon spans the lesion, or the assay reports nothing;
2. neither primer overlaps the lesion, or the mutant allele does not amplify — which looks like a homozygote;
3. primers sit far enough out for the readout (Sanger read clearing the primer, or a short enough product for a gel).

The protected window is both **excluded** from primer placement and **required to be spanned** by the amplicon. That pairing is what distinguishes this from designing primers "near" the cut.

Nuclease geometry for cut-site calculation:

| Nuclease | Spacer | PAM | Cut |
|---|---|---|---|
| SpCas9 | 20 nt | 3′ NGG | Blunt, 3 bp from the PAM |
| SaCas9 | 21 nt | 3′ NNGRRT | Blunt |
| Cas12a | 23 nt | 5′ TTTV | Staggered, 5′ overhangs |

Both strands are searched, so a protospacer copied from a design tool is found regardless of which strand of the locus you happen to have. PAM validation is on by default; `--no-pam` locates the sequence without it.

### qPCR genomic-DNA protection

Graded, not asserted:

- **strong** — a primer straddles a junction, so it has no contiguous binding site on genomic DNA. Works regardless of intron size.
- **moderate** — the amplicon spans an intron ≥ `--min-protective-intron` (default 1 kb).
- **weak** — spans a junction whose intron is too small to protect anything.
- **none** — genomic DNA amplifies identically.

A primer counts as junction-spanning only with ≥ 4 nt on each side; fewer and it still anneals to the unspliced locus.

### In situ probes

Two constraints that have no analogue elsewhere in the package.

**Promoter orientation.** A promoter on the reverse primer transcribes the antisense strand (the probe); on the forward primer, the sense strand (the control). Different promoters at the two ends make one PCR product serve as template for both. All of T7, T3 and SP6 are available in either orientation.

Promoter sequences are the consensus forms from Thermo Fisher's in vitro transcription guide (MAN0011982), in each of which the final G of the core is the +1 base of the transcript:

| Polymerase | Sequence |
|---|---|
| T7 | `TAATACGACTCACTATAGGG` |
| T3 | `AATTAACCCTCACTAAAGGG` |
| SP6 | `ATTTAGGTGACACTATAGAA` |

**Probe cross-hybridisation is a different question from primer specificity.** Priming is decided at a 3' end; hybridisation is decided along the whole molecule. The check therefore scans the amplicon, not its ends: exact k-mers seed alignment diagonals against the avoid-panel, and on each diagonal the *longest window whose identity reaches the threshold* is reported. That window is found in linear time by subtracting the threshold from each position and locating the longest subarray with a non-negative sum.

Two defects were found and fixed while building this, both worth recording because they are easy to reintroduce:

- An earlier version took the span between the first and last seed on a diagonal as the starting block without checking it, so two chance 12-mer matches 600 nt apart were reported as a single 600 nt block at 64% identity.
- The subarray scan was written in floating point. Prefix sums of `0.2` and `-0.8` accumulate error, so two nominally equal sums differ by an ulp and the comparison that should accept the longest window fails. It now uses integers throughout.

#### Threshold calibration

The failure mode to avoid is a chance window: between genuinely homologous but diverged ohnologs, a locally better-conserved patch turns up by luck. Measured over 100 simulated 700 nt backgrounds, the proportion of runs tripped purely by chance:

| identity | min length | 65% background (diverged paralog) | 25% background (unrelated) |
|---|---|---|---|
| 0.80 | 50 | 57% | 0% |
| 0.85 | 50 | 15% | 0% |
| 0.80 | 100 | 2% | 0% |
| **0.85** | **100** | **0%** | **0%** |

85% over 100 nt is the first combination quiet against a diverged paralog while still firing on real conservation, and it is the biologically defensible one: 100 nt at 85% identity holds a riboprobe through a standard high-stringency wash, 50 nt at 80% does not reliably. Both are tunable with `--cross-hyb-identity` and `--cross-hyb-length`.

Only ungapped similarity is modelled, so an indel splits one block into two. That understates the risk slightly, which is the safe direction for a check whose job is to warn.

### dCAPS

The recognition site must straddle the dCAPS primer's 3' end:

```
    primer (bases we control)        template (fixed, allele-specific)
  5'-...........................|...................................-3'
                    [ recognition site spans the boundary ]
                                     ^ variant position
```

Two constraints follow, and both are easy to get wrong by hand:

- the variant must lie **downstream** of the primer's 3' end. If the primer covers the variant, the primer's own sequence overwrites it and both alleles give an identical product;
- every site position downstream of the primer, other than the variant, must **already** match the template — those bases cannot be changed.

Required mismatches are then the site positions inside the primer that differ from the template. Assays are ranked by fewest mismatches, then largest diagnostic band, then primer quality; the 3'-terminal base is kept matched by default.

---

## 6. Complexity and runtime

| Operation | Cost |
|---|---|
| Candidate enumeration | O(L · (max_len − min_len)) placements; staged filters make the structure search proportional to survivors |
| Hairpin search | O(n³) worst case for n ≤ `max_length` (≈ 27), so bounded and small |
| Pair assembly | O(F · R) over the screened pools, bounded by `screen_top` (60) |
| Exhaustive scan | O(subject · primer) with 3'-first early abort; ~300 kb/s |
| Seeded scan | O(1) index lookup plus O(hits · primer) |
| Panel amplicon prediction | Runs on reported pairs only, not on the candidate square |

A 2 kb locus with a 3.6 kb avoid-panel completes in a few seconds on one core.

---

## References

- SantaLucia J. Jr. & Hicks D. (2004) The thermodynamics of DNA structural motifs. *Annu. Rev. Biophys. Biomol. Struct.* 33:415-440.
- Allawi H.T. & SantaLucia J. Jr. (1997-1998) Nearest-neighbour thermodynamics of internal mismatches. *Biochemistry* / *Nucleic Acids Res.* series.
- Owczarzy R. et al. (2004) Effects of sodium ions on DNA duplex oligomers. *Biochemistry* 43:3537-3554.
- Untergasser A. et al. (2012) Primer3 — new capabilities and interfaces. *Nucleic Acids Res.* 40:e115.
- Neff M.M. et al. (1998) dCAPS, a simple technique for the genetic analysis of single nucleotide polymorphisms. *Plant J.* 14:387-392.
- Brownstein M.J. et al. (1996) Modulation of non-templated nucleotide addition by Taq DNA polymerase. *BioTechniques* 20:1004-1010.
- Varshney G.K. et al. (2015) High-throughput gene targeting and phenotyping in zebrafish using CRISPR/Cas9. *Genome Res.* 25:1030-1042.
