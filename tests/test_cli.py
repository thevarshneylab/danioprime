"""End-to-end CLI behaviour, driven through main() rather than a subprocess."""

from __future__ import annotations

import json

import pytest
from helpers import make_seq, mutate

from danioprime.cli import main
from danioprime.sequence import SeqRecord, write_fasta


@pytest.fixture
def locus_file(tmp_path):
    path = tmp_path / "locus.fa"
    write_fasta([SeqRecord(name="locus", seq=make_seq(3000, 2026))], path)
    return path


@pytest.fixture
def ohnolog_file(tmp_path, locus_file):
    from danioprime.sequence import read_single_fasta

    target = read_single_fasta(locus_file).seq
    paralog = mutate(target, 0.10, 31)
    paralog = paralog[:1400] + target[1400:1600] + paralog[1600:]
    path = tmp_path / "ohnolog.fa"
    write_fasta([SeqRecord(name="ohnolog", seq=paralog)], path)
    return path


class TestInformationalCommands:
    def test_no_subcommand_prints_help(self, capsys):
        assert main([]) == 2
        assert "subcommand" in capsys.readouterr().out

    def test_presets(self, capsys):
        assert main(["presets"]) == 0
        out = capsys.readouterr().out
        assert "qpcr" in out
        assert "genotyping-gel" in out

    def test_enzymes(self, capsys):
        assert main(["enzymes"]) == 0
        out = capsys.readouterr().out
        assert "EcoRI" in out
        assert "G^AATTC" in out

    def test_tails(self, capsys):
        assert main(["tails"]) == 0
        out = capsys.readouterr().out
        assert "fluorescent" in out
        assert "TGTAAAACGACGGCCAGT" in out

    def test_info(self, capsys):
        assert main(["info"]) == 0
        out = capsys.readouterr().out
        assert "builtin" in out
        assert "danioprime" in out


class TestDesign:
    def test_basic_design(self, locus_file, capsys):
        assert main(["design", "--template", str(locus_file), "--preset", "genotyping-gel"]) == 0
        out = capsys.readouterr().out
        assert "bp" in out
        assert "penalty" in out

    def test_target_is_spanned(self, locus_file, tmp_path):
        out = tmp_path / "r.json"
        assert (
            main(
                [
                    "design",
                    "--template", str(locus_file),
                    "--target", "1400-1460",
                    "--format", "json",
                    "--out", str(out),
                ]
            )
            == 0
        )
        data = json.loads(out.read_text())
        assert data["pairs"]
        for pair in data["pairs"]:
            assert pair["forward"]["start"] < 1399
            assert pair["reverse"]["end"] > 1460

    def test_exclude_region_is_honoured(self, locus_file, tmp_path):
        out = tmp_path / "r.json"
        main(
            [
                "design",
                "--template", str(locus_file),
                "--exclude", "1-1500",
                "--format", "json",
                "--out", str(out),
            ]
        )
        data = json.loads(out.read_text())
        for pair in data["pairs"]:
            assert pair["forward"]["start"] >= 1500

    def test_tsv_and_bed_and_order_outputs(self, locus_file, tmp_path):
        tsv, bed, order = tmp_path / "p.tsv", tmp_path / "p.bed", tmp_path / "o.tsv"
        assert (
            main(
                [
                    "design",
                    "--template", str(locus_file),
                    "--tsv-out", str(tsv),
                    "--bed-out", str(bed),
                    "--order-out", str(order),
                    "--tails", "fluorescent",
                    "--out", str(tmp_path / "main.txt"),
                ]
            )
            == 0
        )
        assert "forward_seq" in tsv.read_text()
        assert "track name=danioprime" in bed.read_text()
        assert "M13F_universal" in order.read_text()

    def test_impossible_constraints_exit_nonzero(self, locus_file, capsys):
        code = main(
            ["design", "--template", str(locus_file), "--min-tm", "85", "--opt-tm", "87", "--max-tm", "90"]
        )
        assert code == 1
        assert "No primer pair" in capsys.readouterr().out

    def test_contradictory_parameters_are_reported(self, locus_file, capsys):
        code = main(["design", "--template", str(locus_file), "--min-tm", "70", "--max-tm", "60"])
        assert code == 2
        assert "Tm must satisfy" in capsys.readouterr().err

    def test_bare_sequence_is_accepted_as_a_template(self, capsys):
        seq = make_seq(900, 17)
        code = main(["design", "--template", seq, "--min-product-size", "300", "--max-product-size", "800"])
        assert code == 0
        assert "bp" in capsys.readouterr().out

    def test_missing_template_is_reported(self, capsys):
        assert main(["design"]) == 2
        assert "supply --template" in capsys.readouterr().err

    def test_avoid_panel_changes_the_result(self, locus_file, ohnolog_file, tmp_path):
        plain, guarded = tmp_path / "a.json", tmp_path / "b.json"
        args = ["design", "--template", str(locus_file), "--target", "1480-1520", "--format", "json"]
        main([*args, "--out", str(plain)])
        main([*args, "--avoid", str(ohnolog_file), "--out", str(guarded)])
        a = json.loads(plain.read_text())
        b = json.loads(guarded.read_text())
        assert b["pairs"]
        assert all(p["specific"] for p in b["pairs"])
        assert a["pairs"][0]["forward"]["seq"] != b["pairs"][0]["forward"]["seq"] or all(
            p["specific"] for p in a["pairs"]
        )


class TestGenotype:
    @pytest.fixture
    def with_protospacer(self, tmp_path):
        bg = make_seq(3000, 606)
        spacer = make_seq(20, 707)
        seq = bg[:1500] + spacer + "TGG" + bg[1523:]
        path = tmp_path / "g.fa"
        write_fasta([SeqRecord(name="g", seq=seq)], path)
        return path, spacer

    def test_protospacer_drives_the_design(self, with_protospacer, capsys):
        path, spacer = with_protospacer
        assert main(["genotype", "--template", str(path), "--protospacer", spacer]) == 0
        out = capsys.readouterr().out
        assert "SpCas9 cut at 1518" in out
        assert "Predicted allele sizes" in out

    def test_cut_site_may_be_given_directly(self, locus_file, capsys):
        assert main(["genotype", "--template", str(locus_file), "--cut-site", "1500"]) == 0
        assert "cut site" in capsys.readouterr().out

    def test_missing_target_specification_is_reported(self, locus_file, capsys):
        assert main(["genotype", "--template", str(locus_file)]) == 2
        assert "--protospacer or --cut-site" in capsys.readouterr().err

    def test_absent_protospacer_is_explained(self, locus_file, capsys):
        code = main(["genotype", "--template", str(locus_file), "--protospacer", "A" * 20])
        assert code == 2
        err = capsys.readouterr().err
        assert "not found" in err
        assert "--no-pam" in err

    def test_protected_window_is_respected(self, with_protospacer, tmp_path):
        path, spacer = with_protospacer
        out = tmp_path / "r.json"
        main(
            [
                "genotype",
                "--template", str(path),
                "--protospacer", spacer,
                "--protected", "60",
                "--format", "json",
                "--out", str(out),
            ]
        )
        data = json.loads(out.read_text())
        cut = data["cut_site"]["position_1based"] - 1
        for pair in data["pairs"]:
            assert pair["forward"]["end"] <= cut - 60
            assert pair["reverse"]["start"] >= cut + 60


class TestQpcr:
    @pytest.fixture
    def transcript_files(self, tmp_path):
        exon1, exon2, exon3 = make_seq(210, 1), make_seq(330, 2), make_seq(260, 3)
        tx = exon1 + exon2 + exon3
        tx_path = tmp_path / "tx.fa"
        write_fasta([SeqRecord(name="tx", seq=tx)], tx_path)
        gtf = tmp_path / "a.gtf"
        rows = []
        starts = [0, 1000, 3000]
        for i, (start, exon) in enumerate(zip(starts, (exon1, exon2, exon3), strict=True), 1):
            rows.append(
                f'chr1\tx\texon\t{start + 1}\t{start + len(exon)}\t.\t+\t.\t'
                f'gene_id "G"; transcript_id "tx"; gene_name "kcnh6a"; exon_number "{i}";\n'
            )
        gtf.write_text("".join(rows))
        return tx_path, gtf

    def test_junctions_from_gtf(self, transcript_files, capsys):
        tx, gtf = transcript_files
        code = main(
            ["qpcr", "--transcript", str(tx), "--gtf", str(gtf), "--transcript-id", "tx"]
        )
        assert code == 0
        out = capsys.readouterr().out
        assert "Genomic-DNA protection" in out

    def test_explicit_junctions(self, transcript_files, capsys):
        tx, _gtf = transcript_files
        code = main(["qpcr", "--transcript", str(tx), "--junction", "210", "--junction", "540"])
        assert code == 0
        assert "protection" in capsys.readouterr().out

    def test_ohnolog_hint_is_emitted(self, transcript_files, capsys):
        tx, gtf = transcript_files
        main(["qpcr", "--transcript", str(tx), "--gtf", str(gtf), "--transcript-id", "tx"])
        assert "kcnh6b" in capsys.readouterr().err

    def test_unknown_transcript_is_reported(self, transcript_files, capsys):
        tx, gtf = transcript_files
        code = main(
            ["qpcr", "--transcript", str(tx), "--gtf", str(gtf), "--transcript-id", "nope"]
        )
        assert code == 2
        assert "no transcript matching" in capsys.readouterr().err


class TestDcaps:
    def test_assays_are_reported(self, locus_file, capsys):
        from danioprime.sequence import read_single_fasta

        seq = read_single_fasta(locus_file).seq
        alt = next(b for b in "ACGT" if b != seq[1499])
        code = main(
            [
                "dcaps",
                "--template", str(locus_file),
                "--variant", "1500",
                "--variant-base", alt,
                "--min-product-size", "150",
                "--max-product-size", "450",
                "--max-assays", "3",
            ]
        )
        assert code == 0
        out = capsys.readouterr().out
        assert "dCAPS primer" in out
        assert "diagnostic bands" in out

    def test_json_output_round_trips(self, locus_file, tmp_path):
        from danioprime.sequence import read_single_fasta

        seq = read_single_fasta(locus_file).seq
        alt = next(b for b in "ACGT" if b != seq[1499])
        out = tmp_path / "d.json"
        main(
            [
                "dcaps",
                "--template", str(locus_file),
                "--variant", "1500",
                "--variant-base", alt,
                "--min-product-size", "150",
                "--max-product-size", "450",
                "--format", "json",
                "--out", str(out),
            ]
        )
        data = json.loads(out.read_text())
        assert data["assay"] == "dcaps"
        assert data["assays"]
        assert data["variant"]["variant_base"] == alt

    def test_identical_base_is_reported(self, locus_file, capsys):
        from danioprime.sequence import read_single_fasta

        ref = read_single_fasta(locus_file).seq[1499]
        code = main(
            ["dcaps", "--template", str(locus_file), "--variant", "1500", "--variant-base", ref]
        )
        assert code == 2
        assert "nothing to discriminate" in capsys.readouterr().err


class TestCheck:
    def test_evaluates_supplied_primers(self, capsys):
        code = main(["check", "--primers", "F=ACGTTGCATGCCAATGGCTA,R=TAGCCATTGGCATGCAACGT"])
        assert code == 0
        out = capsys.readouterr().out
        assert "Tm" in out
        assert "Pairwise interactions" in out

    def test_reports_panel_cross_reactivity(self, tmp_path, capsys):
        target = make_seq(2000, 808)
        path = tmp_path / "p.fa"
        write_fasta([SeqRecord(name="paralog", seq=target)], path)
        primer = target[500:520]
        code = main(["check", "--primers", primer, "--avoid", str(path)])
        assert code == 0
        out = capsys.readouterr().out
        assert "YES" in out

    def test_fasta_input_is_accepted(self, tmp_path, capsys):
        path = tmp_path / "primers.fa"
        write_fasta(
            [
                SeqRecord(name="F1", seq="ACGTTGCATGCCAATGGCTA"),
                SeqRecord(name="R1", seq="TAGCCATTGGCATGCAACGT"),
            ],
            path,
        )
        assert main(["check", "--primers", str(path)]) == 0
        assert "F1" in capsys.readouterr().out

    def test_empty_primer_list_is_reported(self, capsys):
        assert main(["check", "--primers", " "]) == 2
        assert "no primer sequences" in capsys.readouterr().err


class TestJunction:
    def test_both_junctions_are_designed(self, locus_file, tmp_path, capsys):
        cassette = tmp_path / "c.fa"
        write_fasta([SeqRecord(name="cassette", seq=make_seq(900, 44))], cassette)
        code = main(
            [
                "junction",
                "--template", str(locus_file),
                "--cassette", str(cassette),
                "--insert-at", "1500",
                "--preset", "genotyping-gel",
            ]
        )
        assert code == 0
        out = capsys.readouterr().out
        assert "5prime junction" in out
        assert "3prime junction" in out


class TestRemoteVariantSources:
    """The CLI surface for published variant resources."""

    @pytest.fixture
    def mirror(self, tmp_path_factory):
        pysam = pytest.importorskip("pysam")
        import random

        directory = tmp_path_factory.mktemp("cli_mirror")
        rng = random.Random(8)
        lines = [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr5,length=2000000>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNHGRI1",
        ]
        pos = 0
        while True:
            pos += rng.randint(60, 400)
            if pos >= 2_000_000:
                break
            ref = rng.choice("ACGT")
            alt = rng.choice([b for b in "ACGT" if b != ref])
            lines.append(
                f"chr5\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t{rng.choice(['0/1', '1/1'])}"
            )
        plain = directory / "m.vcf"
        plain.write_text("\n".join(lines) + "\n", newline="\n")
        gz = directory / "m.vcf.gz"
        pysam.tabix_compress(str(plain), str(gz), force=True)
        pysam.tabix_index(str(gz), preset="vcf", force=True)
        return gz

    @pytest.fixture
    def reference(self, tmp_path):
        from danioprime.sequence import SeqRecord, write_fasta

        path = tmp_path / "danRer11.fa"
        write_fasta([SeqRecord(name="chr5", seq=make_seq(1_010_000, 4))], path, wrap=60)
        return path

    def test_sources_lists_nhgri1_with_its_citation(self, capsys):
        assert main(["sources"]) == 0
        out = capsys.readouterr().out
        assert "nhgri1" in out
        assert "danRer11" in out
        assert "LaFave" in out
        assert "hub.txt" in out

    def test_sources_warns_about_the_assembly(self, capsys):
        main(["sources"])
        out = capsys.readouterr().out
        assert "GRCz12tu" in out
        assert "not transferable" in out.replace("NOT", "not")

    def test_remote_mask_reduces_the_candidate_pool(self, mirror, reference, capsys):
        from rangeserver import RangeServer

        args = [
            "genotype", "--region", "chr5:1,000,000-1,006,000",
            "--reference", str(reference), "--cut-site", "3000",
            "--preset", "genotyping-gel", "--show", "0",
        ]
        assert main(args) == 0
        plain = capsys.readouterr().out

        with RangeServer(mirror.parent) as base:
            assert (
                main([*args, "--variants-url", f"{base}/{mirror.name}",
                      "--reject-variant-overlap"])
                == 0
            )
        captured = capsys.readouterr()
        masked = captured.out
        assert "Variant mask" in masked
        assert "variant(s) from" in captured.err

        def forwards(text: str) -> int:
            line = next(line for line in text.splitlines() if "candidates:" in line)
            return int(line.split()[1])

        assert forwards(masked) < forwards(plain)

    def test_assembly_mismatch_is_refused_through_the_cli(self, mirror, reference, capsys):
        from rangeserver import RangeServer

        from danioprime.variants.sources import VARIANT_SOURCES, VariantSource

        with RangeServer(mirror.parent) as base:
            patched = dict(VARIANT_SOURCES)
            patched["testsrc"] = VariantSource(
                name="test", description="d", assembly="danRer11",
                data_url=f"{base}/{mirror.name}",
                index_url=f"{base}/{mirror.name}.tbi",
                kind="vcf", citation="c", homepage="https://example.invalid",
            )
            VARIANT_SOURCES.update(patched)
            try:
                code = main([
                    "genotype", "--region", "chr5:1,000,000-1,006,000",
                    "--reference", str(reference), "--cut-site", "3000",
                    "--variant-source", "testsrc", "--assembly", "GRCz12tu",
                ])
            finally:
                VARIANT_SOURCES.pop("testsrc", None)
        assert code == 2
        err = capsys.readouterr().err
        assert "not transferable" in err

    def test_ad_hoc_url_warns_that_the_assembly_is_unchecked(self, mirror, reference, capsys):
        from rangeserver import RangeServer

        with RangeServer(mirror.parent) as base:
            main([
                "genotype", "--region", "chr5:1,000,000-1,006,000",
                "--reference", str(reference), "--cut-site", "3000",
                "--show", "0", "--variants-url", f"{base}/{mirror.name}",
            ])
        err = capsys.readouterr().err
        assert "cannot determine the assembly" in err

    def test_remote_source_without_a_contig_is_reported(self, capsys):
        code = main([
            "design", "--template", make_seq(2000, 5),
            "--variant-source", "nhgri1", "--assembly", "danRer11",
        ])
        assert code == 2
        assert "which contig" in capsys.readouterr().err

    def test_strain_and_heterozygous_filters_reach_the_remote_source(
        self, mirror, reference, capsys
    ):
        from rangeserver import RangeServer

        base_args = [
            "genotype", "--region", "chr5:1,000,000-1,006,000",
            "--reference", str(reference), "--cut-site", "3000", "--show", "0",
        ]
        with RangeServer(mirror.parent) as base:
            url = f"{base}/{mirror.name}"
            main([*base_args, "--variants-url", url])
            everything = capsys.readouterr().err
            main([*base_args, "--variants-url", url, "--strain", "NHGRI1",
                  "--heterozygous-only"])
            het = capsys.readouterr().err

        def count(text: str) -> int:
            line = next(line for line in text.splitlines() if "variant(s) from" in line)
            return int(line.split()[1])

        assert count(het) < count(everything)
        assert "heterozygous positions only" in het
