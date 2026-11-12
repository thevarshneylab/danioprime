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

**NHGRI-1 (Burgess lab, NHGRI)** is the most directly useful set for this purpose, and `danioprime` can query it by region without any download:

```bash
danioprime sources                       # what is available, and how to cite it
danioprime genotype --region chr5:1,000,000-1,006,000 \
    --reference danRer11.fa --assembly danRer11 \
    --variant-source nhgri1 --strain NHGRI1 --heterozygous-only
```

NHGRI-1 was derived from founder parents sequenced to ~50×, and the stated purpose of the resource is that a researcher can know the exact sequence, with all variants, of a given region — which is exactly the question primer design asks. The same data is browsable as the `ZebrafishGenomics` UCSC track hub.

There is a second file on that site worth knowing about: an **invariant** BED, marking the portion of the genome *confirmed to match* the reference. That is a stronger guarantee than the absence of a called variant, because absence can also mean no coverage or unmappable sequence. If you want positive assurance rather than the absence of a warning, design inside the invariant regions; pass the BED as `--variants` inverted, or intersect it yourself and pass the complement with `--exclude`.

**Sanger Zebrafish Genomes Project** publishes calls across strains and additional *Danio* species. Its data-sharing terms ask you to contact the project before publishing genome-scale analyses — read the policy on the project page. Routine primer design consumes a handful of positions rather than performing a genome-wide analysis, but check if in doubt.

**Your own lab's resequencing** beats any public resource, because it describes the fish in your tanks. Pass it with `--variants`, or bgzip and tabix it and pass `--variants-url` with a local path to get region-scoped reads.

#### All of this is on danRer11

Every published set above is on danRer11 (GRCz11). If you are working in GRCz12tu coordinates, they do not apply, and `danioprime` will refuse the combination rather than mask the wrong bases. Either work in danRer11 for the variant check, or lift the variant coordinates over and pass `--allow-assembly-mismatch` to say you have done so.

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

Long enough that the read clears the primer and any indel before quality degrades. For a crispant with a possible large deletion, widen the product-size window rather than narrowing the region: the read has to reach past a lesion whose size you do not yet know.

### `amplicon-seq` (180-280 bp)

Sized for short-read sequencing of an edited locus. Add your platform adapters as tails or in a second PCR; they are not part of the annealing portion and should not enter the Tm calculation.

### `insitu-probe` (500-1500 bp)

Template for riboprobe synthesis, with relaxed structure limits because the product is transcribed rather than amplified efficiently. Use the dedicated `insitu` mode rather than `design`, because two things here differ from every other assay.

**Orientation decides which strand you get.** A promoter on the reverse primer transcribes the antisense strand — the probe. On the forward primer it transcribes the sense strand — the negative control. Backwards, the probe cannot hybridise and the control can, and the result merely looks like a failed experiment. T7, T3 and SP6 are offered in both orientations, and the dual schemes put different promoters at the two ends so one product makes both probes.

**The ohnolog check has to be different here.** Everywhere else, specificity is decided at a primer's 3' end. A riboprobe has no 3' end that matters: it is hundreds of bases of RNA and binds anything similar enough along its length. So a pair whose primers discriminate the ohnolog perfectly can still yield a probe that stains both genes — the primers sit in divergent sequence while the body of the amplicon is conserved. The `insitu` mode scans the amplicon itself and reports the longest block it shares with the panel.

In practice ohnolog coding sequence is often too conserved for *any* probe over the ORF to discriminate. The usual escape is the **3' UTR**, which diverges considerably faster than coding sequence; steer the probe there with `--target`. If no divergent region exists, the honest options are to accept that the probe reports both genes and say so, or to switch to a method with single-base discrimination.

---

## Practical geometry for CRISPR work

A few numbers worth having in mind, which the defaults reflect:

- **Give the region, not the cut.** `--target` takes the span you need interrogated — the exon you disrupted, the window you edited. There is deliberately no cut-site input: large deletions, hundreds of base pairs or more, are common with a single guide and routine with two, so a primer placed a fixed distance from the predicted cut can land inside the deletion and silently fail to amplify the mutant allele. `docs/algorithm.md` sets out the reasoning.
- **Widen the region when you are unsure, or give none.** A design across the whole template produces a wide amplicon that certainly contains the lesion, which beats a narrow one placed around a coordinate that may be wrong. Narrow it only when the readout forces you to — fragment analysis wants a short product, and a short product cannot span a wide region.
- **Knock-in junction PCR.** One primer **must** sit in genomic sequence outside the homology arm. A pair with both primers in the cassette amplifies the donor plasmid just as happily as the integrated allele — the classic false positive in knock-in screening. `--homology-arm` keeps primers clear of the arm.
- **Founder screening.** G0 injected animals are mosaic, so a genotyping assay there reports a mixture. Design for the F1, and treat G0 results as presence/absence of editing rather than a genotype.

---

## Further reading

- ZFIN (<https://zfin.org>) — nomenclature, genome browsers, alleles, the authority on which genes are ohnologs
- Alliance of Genome Resources — orthology and comparative data
- Sanger Zebrafish Genomes Project — variant calls across strains, and the data-use policy
- Ensembl / UCSC — browsers, GTF annotation and liftOver chain files
