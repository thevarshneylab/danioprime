# Zebrafish notes

Background for the design decisions in `danioprime`, and the practical details behind the presets. Where a number is an estimate that varies between studies, it says so.

---

## Ohnologs and the teleost genome duplication

The teleost lineage underwent a whole-genome duplication (the TGD, commonly placed in the ray-finned fish stem lineage several hundred million years ago). Most duplicates were subsequently lost, but a substantial minority — frequently quoted as roughly a fifth of the duplicated pairs, with estimates varying by method and study — were retained. These retained duplicate pairs are **ohnologs**, and they are why so many zebrafish genes have two copies where humans have one.

### Why this breaks primer design

Retained ohnologs are often strongly conserved in coding sequence, because both copies remained under selection. Identity in exons is routinely high enough that a 20-mer drawn from one gene matches the other perfectly. Three consequences:

1. **Genotyping** can amplify both loci, giving a band pattern that does not track the allele you are scoring.
2. **RT-qPCR** silently reports the *sum* of two genes' expression. This is the most damaging case, because the result looks entirely normal.
3. **In situ probes** cross-hybridise, and the expression pattern you publish belongs to two genes.

Intronic and intergenic sequence diverges much faster, which is the practical escape route: primers placed in introns or UTRs usually discriminate easily, while primers in conserved exon blocks usually cannot. For a cDNA assay that option is unavailable, which is why the `qpcr` mode treats the ohnolog transcript as an explicit avoid-panel.

### The naming convention, and its limits

ZFIN distinguishes ohnologs with an `a`/`b` suffix: `kcnh6a` / `kcnh6b`. `danioprime` uses this to *suggest* the paralog you probably need to avoid — the `qpcr` mode prints a hint when it sees such a symbol.

**Treat it as a prompt, not an answer.** The suffix is a nomenclature convention, not an orthology assertion. Many gene symbols end in `a` or `b` for unrelated reasons, some ohnolog pairs are not named this way, and some `a`/`b` pairs are tandem duplicates rather than TGD ohnologs. Confirm the pairing against ZFIN, the Alliance of Genome Resources, or an Ensembl gene tree, then supply the sequence yourself with `--avoid`.

### Discriminating deliberately

The robust way to separate ohnologs is **not** a Tm difference. Annealing temperature in a real thermocycler is permissive enough that a 3-4 °C difference often fails to discriminate, and it degrades further in a gradient or on a different machine.

The robust way is to place the primer's 3'-terminal bases over a position where the two genes differ. The polymerase then cannot extend on the paralog at all, regardless of annealing temperature. `danioprime` rewards exactly this: the penalty model gives credit when the closest panel match carries a 3'-proximal mismatch, and the report tells you how many bases from the 3' end that mismatch sits.

When you need the discrimination to be absolute, use `--strict-specificity`, which rejects any primer that can prime on the panel even when no product is predicted. Be aware this can leave no design at all.

---

## Strain variation

Wild-type zebrafish lines are **not** isogenic, and this surprises people coming from mouse work. TU (Tübingen), AB, TL (Tupfel long-fin) and WIK all carry substantial standing variation, both between and within lines, and individual facility stocks drift.

### The failure mode

A primer sitting on a position that segregates in your line amplifies one allele and not the other. On a genotyping gel that is indistinguishable from a homozygote. It is the worst class of failure because it produces a clean, plausible, wrong result — and because it is heritable in your stock, it reproduces.

### How this tool handles it

- A variant anywhere in the primer footprint is **penalised**.
- A variant within the 3'-terminal window is **always disqualifying**, because that is where a mismatch stops extension outright.
- `--reject-variant-overlap` makes any overlap disqualifying.

### `--heterozygous-only` is the useful flag

A position homozygous for the alternate allele in your line is **not** a design hazard. Your primer has to match your fish, not the reference. Masking such positions throws away perfectly good primer sites for no benefit.

```bash
--variants sanger_zfish.vcf.gz --strain TU --heterozygous-only
```

This keeps only positions actually segregating in TU. In practice it recovers a useful fraction of the masked sequence.

### Getting variant data

The Sanger Zebrafish Genomes Project publishes variant calls across strains and additional *Danio* species. Note that its data-sharing terms ask you to contact the project before publishing genome-scale analyses of their data — read the policy on the project page. For routine primer design you are consuming a handful of positions, not performing a genome-wide analysis, but check if in doubt.

If your lab has resequenced its own stocks, that VCF is better than any public resource, because it describes the fish in your tanks.

---

## Reference assemblies and coordinates

GRCz11 (May 2017) was the reference for many years, at roughly 1.37 Gb. It has now been superseded: complete *de novo* assemblies from homozygous TU and AB material produced **GRCz12tu** and **GRCz12ab**, each around 1.4 Gb — roughly 100 Mb more sequence than GRCz11, including a large amount of previously unassembled material — and both have been elevated to reference status. ZFIN has moved to GRCz12tu for genomic data.

### What this means for you

**Coordinates are not transferable between GRCz11 and GRCz12tu.** Shifts of tens to hundreds of kilobases are normal. Positions in an old lab notebook, an old VCF, or an old BED file will silently point somewhere else.

Practical guidance:

- Prefer **identifiers over coordinates** wherever possible: gene symbols, ZFIN IDs, Ensembl transcript IDs. `danioprime` takes transcripts by ID from a GTF for exactly this reason.
- When you must convert, use liftOver chain files for the specific assembly pair, and **spot-check** the result by confirming the sequence at the new coordinates is what you expect.
- Keep the assembly name in your file names and primer records. `kcnh6a_F_GRCz12tu.txt` costs nothing now and saves a day later.
- Use `--template-offset` so a VCF or BED in genome coordinates lines up with a template you cut out of a contig, and `--bed-out` to push primer footprints back into genome coordinates for browser inspection.

---

## Notes on the presets

The presets are starting points shaped by what the readout actually requires. Any field can be overridden.

### `genotyping-gel` (180-600 bp)

Standard indel genotyping. The thing to notice is the **resolution limit**: on a 300 bp product, a 2% agarose gel does not reliably separate a few-base indel from wild type by eye. `danioprime` warns when the predicted allele sizes are not gel-resolvable, and then the options are a shorter product, PAGE, fragment analysis, or a dCAPS digest.

### `genotyping-fluorescent` (90-300 bp)

For fragment analysis. Pair it with `--tails fluorescent`, which emits the **three** oligos the system needs: an M13F-tailed forward primer, a pig-tailed reverse primer, and one fluorophore-labelled M13F primer. That last oligo is shared across every M13F-tailed assay in the lab, so it is listed once per order sheet rather than once per assay, and it should be costed that way.

The pig-tail (`GTTTCTT`) forces complete non-templated 3' A addition, so peaks come out as a single size rather than a ±1 nt doublet. Without it, scoring a 1 bp indel is guesswork.

Reported product sizes exclude the tails — add 18 nt and 7 nt to predict peak positions.

### `qpcr` (70-160 bp, tight Tm window)

Short amplicons, a 1 °C pair Tm difference, strict structure limits. Combine with `--require-junction-spanning` and an `--avoid` panel carrying the ohnolog **transcript**. Both matter: the junction-spanning primer protects against genomic DNA, the panel protects against measuring two genes as one.

### `sanger` (500-1200 bp)

Long enough that the read clears the primer and any indel before quality degrades. For a crispant with a possible large deletion, raise `--protected` well above your largest expected indel.

### `amplicon-seq` (180-280 bp)

Sized for short-read sequencing of an edited locus. Add your platform adapters as tails or in a second PCR; they are not part of the annealing portion and should not enter the Tm calculation.

### `insitu-probe` (500-1500 bp)

Template for riboprobe synthesis, with relaxed structure limits because the product is transcribed rather than amplified efficiently. Use `--tails riboprobe-t7` to put a T7 promoter on the reverse primer, giving an antisense probe; put it on the forward primer for a sense control. **Check the ohnolog here too** — a cross-hybridising probe produces a published expression pattern belonging to two genes.

---

## Practical geometry for CRISPR work

A few numbers worth having in mind, which the defaults reflect:

- **Protected window.** Set `--protected` above the **largest indel you expect**, not the one you got. The assay has to work before the alleles are known, and a primer inside a deleted span does not amplify the mutant allele at all. Large deletions — hundreds of base pairs or more — are common enough with a single guide, and routine with two.
- **Primer distance from the cut.** Far enough that a Sanger read has cleared the primer and settled before reaching the lesion; close enough that the product size suits the readout. The `--min-flank` / `--max-flank` defaults (40 / 400) are a reasonable compromise for gel genotyping.
- **Knock-in junction PCR.** One primer **must** sit in genomic sequence outside the homology arm. A pair with both primers in the cassette amplifies the donor plasmid just as happily as the integrated allele — the classic false positive in knock-in screening. `--homology-arm` keeps primers clear of the arm.
- **Founder screening.** G0 injected animals are mosaic, so a genotyping assay there reports a mixture. Design for the F1, and treat G0 results as presence/absence of editing rather than a genotype.

---

## Further reading

- ZFIN (<https://zfin.org>) — nomenclature, genome browsers, alleles, the authority on which genes are ohnologs
- Alliance of Genome Resources — orthology and comparative data
- Sanger Zebrafish Genomes Project — variant calls across strains, and the data-use policy
- Ensembl / UCSC — browsers, GTF annotation and liftOver chain files
