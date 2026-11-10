# Contributing

## The most valuable contribution

**Bench outcomes.** The algorithms here are tested for correctness, but nobody has yet shown that `danioprime` designs succeed more often than Primer3 plus manual checking. That evidence does not exist, and it is the thing that would make this tool worth trusting.

If you design primers with this and record what happened, please open an issue with:

- the command you ran (and the `--format json` output if you kept it);
- whether each pair amplified, and whether the assay discriminated what it was meant to;
- for genotyping, whether the allele calls held up against sequencing;
- for qPCR, whether the no-RT control was clean.

Failures are more informative than successes. A pair that the tool rated highly and that failed at the bench points at a real modelling gap.

## Development setup

```bash
git clone https://github.com/thevarshneylab/danioprime
cd danioprime
pip install -e ".[dev]"
pytest
ruff check .
```

Tests must pass with no required dependencies installed. Anything needing `primer3-py`, `blastn` or `bowtie2` must be guarded with `pytest.mark.skipif`.

## Conventions

- **No required runtime dependencies.** The default thermodynamic backend stays pure Python so results are reproducible across platforms and CI. Optional extras are fine behind a guarded import with an error message that names the install command.
- **Coordinates are zero-based and half-open internally**, converted to 1-based only at the edges (CLI input, report output). Every function that takes a position says which convention it uses.
- **Comments explain why, not what.** If a threshold has a biological reason, state the reason.
- **Errors tell the user what to do.** `"no usable forward primer"` is not enough; name the filter that removed the candidates and the flag that would relax it.
- **New biology needs a citation** in the docstring, and a test that an independent calculation agrees with.

## Adding an assay mode

1. Put the geometry in `src/danioprime/assays/`, building a `DesignRequest` and calling `run_design`. Do not reimplement candidate enumeration.
2. Assert the geometry in tests — "every reported amplicon spans the protected window" is the kind of claim that must be checked directly, not assumed from the search windows.
3. Add a CLI subcommand with a `--help` string that explains the *assay*, not the parameters.
4. Document the constraint the mode enforces in `docs/algorithm.md` section 5.

## Adding enzymes

The table in `src/danioprime/assays/enzymes.py` is a curated common subset. Additions should be widely stocked, cut **inside** the recognition site (type IIS enzymes that cut at a distance are out of scope for dCAPS fragment prediction), and carry a correct top-strand cut offset. Users can always supply their own table with `--enzymes`.

## Reporting a bug

Include the command, the version (`danioprime --version`), and a minimal template that reproduces it — a synthetic sequence is ideal, since real locus sequence is rarely necessary to show the problem.
