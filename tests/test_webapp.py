"""The browser entry point.

This is the layer the web UI calls, so it is tested here rather than left to
JavaScript. The contract that matters: `run` never raises, and every failure
comes back as a readable message rather than a traceback.
"""

from __future__ import annotations

import json

import pytest
from helpers import make_seq, mutate

from danioprime.webapp import (
    MODES,
    TUNABLE,
    build_params,
    parse_one_sequence,
    parse_sequences,
    run,
    run_json,
    ui_metadata,
)


@pytest.fixture
def locus_text():
    return ">locus mock\n" + make_seq(3000, 2026)


@pytest.fixture
def with_protospacer():
    bg = make_seq(3000, 606)
    spacer = make_seq(20, 707)
    return ">g\n" + bg[:1500] + spacer + "TGG" + bg[1523:], spacer


@pytest.fixture
def ohnolog_text(locus_text):
    target = locus_text.split("\n", 1)[1]
    paralog = mutate(target, 0.10, 31)
    paralog = paralog[:1400] + target[1400:1600] + paralog[1600:]
    return ">ohnolog\n" + paralog


class TestInputParsing:
    def test_bare_sequence(self):
        records = parse_sequences("acgt acgt\nACGT")
        assert len(records) == 1
        assert records[0].seq == "ACGTACGTACGT"
        assert records[0].name == "input"

    def test_fasta_with_several_records(self):
        records = parse_sequences(">a\nACGTACGT\n>b desc\nTTTTAAAA")
        assert [r.name for r in records] == ["a", "b"]
        assert records[1].description == "desc"

    def test_empty_text_gives_no_records(self):
        assert parse_sequences("") == []
        assert parse_sequences("   \n  ") == []

    def test_one_sequence_helper_rejects_several(self):
        with pytest.raises(Exception, match="holds 2 sequences"):
            parse_one_sequence(">a\nACGT\n>b\nTTTT", what="template")

    def test_one_sequence_helper_rejects_none(self):
        with pytest.raises(Exception, match="no template supplied"):
            parse_one_sequence("", what="template")

    def test_header_without_a_sequence_is_reported_by_name(self):
        # A header with no bases parses into a zero-length record, which would
        # otherwise fail much later with a confusing message.
        with pytest.raises(Exception, match="header but no sequence"):
            parse_sequences(">")
        with pytest.raises(Exception, match=r"no sequence: a\b"):
            parse_sequences(">a\n>b\nACGT")


class TestParameterBuilding:
    def test_preset_is_applied(self):
        assert build_params({"preset": "qpcr"}).min_product_size == 70

    def test_default_when_no_preset(self):
        assert build_params({}).opt_tm == 60.0

    def test_overrides_are_typed_correctly(self):
        params = build_params({"overrides": {"max_pairs": "7", "min_tm": "58.5"}})
        assert params.max_pairs == 7
        assert isinstance(params.max_pairs, int)
        assert params.min_tm == 58.5

    def test_blank_overrides_are_ignored(self):
        params = build_params({"preset": "qpcr", "overrides": {"min_tm": "", "max_pairs": None}})
        assert params.min_tm == 59.0

    def test_only_whitelisted_parameters_may_be_set(self):
        with pytest.raises(Exception, match="cannot be set from the web form"):
            build_params({"overrides": {"thermo_backend": "primer3"}})

    def test_every_tunable_is_a_real_parameter(self):
        from dataclasses import fields

        from danioprime.config import DesignParams

        names = {f.name for f in fields(DesignParams)}
        assert set(TUNABLE) <= names

    def test_specificity_switches(self):
        params = build_params({"strict_specificity": True, "allow_paralog_amplicon": True})
        assert params.reject_paralog_priming
        assert not params.reject_paralog_amplicon

    def test_reaction_conditions(self):
        params = build_params({"conditions": {"divalent": "3.0", "dmso": "5"}})
        assert params.conditions.divalent_mM == 3.0
        assert params.conditions.dmso_pct == 5.0

    def test_contradictory_overrides_are_rejected(self):
        with pytest.raises(Exception, match="Tm must satisfy"):
            build_params({"overrides": {"min_tm": 70, "max_tm": 60}})


class TestGenotypeMode:
    def test_protospacer_design(self, with_protospacer):
        text, spacer = with_protospacer
        r = run({"mode": "genotype", "template": text, "protospacer": spacer})
        assert r["ok"]
        assert r["pairs"]
        assert r["mode"] == "genotyping"
        assert any("allele sizes" in t["title"] for t in r["tables"])
        assert "cut site" in r["report"]

    def test_cut_position_design(self, locus_text):
        r = run({"mode": "genotype", "template": locus_text, "cut_site": 1500})
        assert r["ok"]
        assert r["pairs"]

    def test_protected_window_is_honoured(self, with_protospacer):
        text, spacer = with_protospacer
        r = run(
            {
                "mode": "genotype",
                "template": text,
                "protospacer": spacer,
                "protected": 70,
                "overrides": {"max_pairs": 5},
            }
        )
        cut = r["extra"]["cut_site"]["position_1based"] - 1
        for pair in r["pairs"]:
            f_end = int(pair["forward_pos"].split("-")[1])
            r_start = int(pair["reverse_pos"].split("-")[0])
            assert f_end <= cut - 70
            assert r_start - 1 >= cut + 70

    def test_absent_protospacer_gives_actionable_advice(self, locus_text):
        r = run({"mode": "genotype", "template": locus_text, "protospacer": "A" * 20})
        assert not r["ok"]
        assert "PAM" in r["error"]
        assert "cut position" in r["error"]

    def test_no_target_specified(self, locus_text):
        r = run({"mode": "genotype", "template": locus_text})
        assert not r["ok"]
        assert "protospacer or a cut position" in r["error"]

    def test_ohnolog_panel_is_applied(self, locus_text, ohnolog_text):
        r = run(
            {
                "mode": "genotype",
                "template": locus_text,
                "avoid": ohnolog_text,
                "cut_site": 1500,
                "overrides": {"max_pairs": 4},
            }
        )
        assert r["ok"]
        assert all(p["specific"] for p in r["pairs"])
        assert any("avoid-panel" in n for p in r["pairs"] for n in p["notes"])

    def test_tails_produce_an_order_sheet(self, locus_text):
        r = run(
            {"mode": "genotype", "template": locus_text, "cut_site": 1500, "tails": "fluorescent"}
        )
        assert "M13F_universal" in r["order_sheet"]
        # Rows carry the annotated form: the tail in lower case, so the page
        # and anything pasted out of it show where the tail stops.
        assert r["pairs"][0]["forward"].startswith("TGTAAAACGACGGCCAGT".lower())
        assert any("Three-primer" in n for n in r["order_notes"])


class TestQpcrMode:
    @pytest.fixture
    def transcript(self):
        e1, e2, e3 = make_seq(210, 1), make_seq(330, 2), make_seq(260, 3)
        gtf = []
        for i, (start, exon) in enumerate(zip([0, 1000, 3000], (e1, e2, e3), strict=True), 1):
            gtf.append(
                f'chr1\tx\texon\t{start + 1}\t{start + len(exon)}\t.\t+\t.\t'
                f'gene_id "G"; transcript_id "tx"; gene_name "demogenea"; exon_number "{i}";'
            )
        return ">tx\n" + e1 + e2 + e3, "\n".join(gtf)

    def test_explicit_junctions(self, transcript):
        text, _gtf = transcript
        r = run({"mode": "qpcr", "template": text, "junctions": [210, 540]})
        assert r["ok"]
        assert any("protection" in t["title"].lower() for t in r["tables"])

    def test_junctions_from_gtf(self, transcript):
        text, gtf = transcript
        r = run({"mode": "qpcr", "template": text, "gtf": gtf, "transcript_id": "demogenea"})
        assert r["ok"]
        assert r["extra"]["junction_offsets"] == [210, 540]

    def test_unknown_transcript_lists_alternatives(self, transcript):
        text, gtf = transcript
        r = run({"mode": "qpcr", "template": text, "gtf": gtf, "transcript_id": "nope"})
        assert not r["ok"]
        assert "no transcript matching" in r["error"]

    def test_junction_spanning_requirement(self, transcript):
        text, _gtf = transcript
        r = run(
            {
                "mode": "qpcr",
                "template": text,
                "junctions": [210, 540],
                "require_junction_spanning": True,
                "overrides": {
                    "min_tm": 55, "max_tm": 65, "min_gc": 30, "max_gc": 70,
                    "require_gc_clamp": 0, "max_pairs": 4,
                },
            }
        )
        assert r["ok"]
        assert r["pairs"]
        assert all(p["spans_junction"] is not None for p in r["pairs"])

    def test_malformed_gtf_is_reported(self, transcript):
        text, _gtf = transcript
        r = run({"mode": "qpcr", "template": text, "gtf": "not a gtf at all"})
        assert not r["ok"]
        assert "exon features" in r["error"]


class TestDcapsMode:
    @pytest.fixture
    def setup(self):
        seq = make_seq(1200, 55)
        ref = seq[599]
        alt = next(b for b in "ACGT" if b != ref)
        return ">locus\n" + seq, 600, alt

    def test_assays_are_returned(self, setup):
        text, pos, alt = setup
        r = run(
            {
                "mode": "dcaps",
                "template": text,
                "variant_position": pos,
                "variant_base": alt,
                "preset": "genotyping-gel",
                "overrides": {"min_product_size": 150, "max_product_size": 450},
                "max_assays": 3,
            }
        )
        assert r["ok"]
        table = r["tables"][0]
        assert table["rows"]
        assert "dCAPS primer" in r["report"]
        assert len(table["rows"]) <= 3

    def test_missing_variant_details(self, setup):
        text, _pos, _alt = setup
        r = run({"mode": "dcaps", "template": text, "variant_base": "A"})
        assert not r["ok"]
        assert "variant position" in r["error"]
        r = run({"mode": "dcaps", "template": text, "variant_position": 600})
        assert not r["ok"]
        assert "variant base" in r["error"]

    def test_identical_base_is_explained(self, setup):
        text, pos, _alt = setup
        ref = text.split("\n", 1)[1][pos - 1]
        r = run({"mode": "dcaps", "template": text, "variant_position": pos, "variant_base": ref})
        assert not r["ok"]
        assert "nothing to discriminate" in r["error"]

    def test_json_payload_is_valid(self, setup):
        text, pos, alt = setup
        r = run(
            {
                "mode": "dcaps", "template": text, "variant_position": pos, "variant_base": alt,
                "overrides": {"min_product_size": 150, "max_product_size": 450},
            }
        )
        data = json.loads(r["json"])
        assert data["variant"]["variant_base"] == alt
        assert data["assays"]


class TestInSituMode:
    @pytest.fixture
    def transcripts(self):
        from helpers import make_seq, mutate

        left, mid, right = make_seq(700, 11), make_seq(900, 22), make_seq(700, 33)
        target = left + mid + right
        paralog = mutate(left, 0.35, 1) + mutate(mid, 0.04, 2) + mutate(right, 0.35, 3)
        return ">demogenea\n" + target, ">demogeneb\n" + paralog

    def _payload(self, transcripts, **kw):
        a, b = transcripts
        base = {
            "mode": "insitu", "template": a, "avoid": b, "preset": "insitu-probe",
            "overrides": {
                "min_product_size": "450", "opt_product_size": "550",
                "max_product_size": "650", "max_pairs": "3",
            },
        }
        base.update(kw)
        return base

    def test_default_scheme_designs_an_antisense_probe(self, transcripts):
        r = run(self._payload(transcripts))
        assert r["ok"]
        assert r["pairs"]
        assert r["pairs"][0]["reverse"].startswith("TAATACGACTCACTATAGGG".lower())

    @pytest.mark.parametrize(
        "scheme,tail",
        [
            ("riboprobe-t7", "TAATACGACTCACTATAGGG"),
            ("riboprobe-t3", "AATTAACCCTCACTAAAGGG"),
            ("riboprobe-sp6", "ATTTAGGTGACACTATAGAA"),
        ],
    )
    def test_all_three_polymerases_are_reachable(self, transcripts, scheme, tail):
        r = run(self._payload(transcripts, scheme=scheme))
        assert r["ok"], r.get("error")
        assert r["pairs"][0]["reverse"].startswith(tail.lower())

    def test_sense_schemes_tail_the_forward_primer(self, transcripts):
        r = run(self._payload(transcripts, scheme="riboprobe-t3-sense"))
        assert r["ok"]
        assert r["pairs"][0]["forward"].startswith("AATTAACCCTCACTAAAGGG".lower())
        assert not r["pairs"][0]["reverse"].startswith("AATTAACCCTCACTAAAGGG".lower())

    def test_dual_scheme_tails_both_ends_differently(self, transcripts):
        r = run(self._payload(transcripts, scheme="riboprobe-dual-t7-sp6"))
        assert r["ok"]
        assert r["pairs"][0]["forward"].startswith("ATTTAGGTGACACTATAGAA".lower())
        assert r["pairs"][0]["reverse"].startswith("TAATACGACTCACTATAGGG".lower())

    def test_cross_hybridisation_table_is_returned(self, transcripts):
        r = run(self._payload(transcripts, targets="900-1100"))
        assert r["ok"]
        table = next(t for t in r["tables"] if "cross-hybridisation" in t["title"].lower())
        assert any(row[1] == "YES" for row in table["rows"])

    def test_sense_control_table_is_returned(self, transcripts):
        r = run(self._payload(transcripts, scheme="riboprobe-t7"))
        assert any("Sense control" in t["title"] for t in r["tables"])

    def test_sense_control_can_be_suppressed(self, transcripts):
        r = run(self._payload(transcripts, scheme="riboprobe-t7", no_sense_control=True))
        assert not any("Sense control" in t["title"] for t in r["tables"])

    def test_a_non_promoter_scheme_is_refused_readably(self, transcripts):
        r = run(self._payload(transcripts, scheme="fluorescent"))
        assert not r["ok"]
        assert "promoter" in r["error"]

    def test_thresholds_are_passed_through(self, transcripts):
        strict = run(self._payload(transcripts, cross_hyb_identity="99", targets="900-1100"))
        lenient = run(self._payload(transcripts, cross_hyb_identity="60", cross_hyb_length="20",
                                    targets="900-1100"))
        assert strict["ok"] and lenient["ok"]
        strict_t = next(t for t in strict["tables"] if "cross-hyb" in t["title"].lower())
        lenient_t = next(t for t in lenient["tables"] if "cross-hyb" in t["title"].lower())
        assert sum(r[1] == "YES" for r in lenient_t["rows"]) >= sum(
            r[1] == "YES" for r in strict_t["rows"]
        )


class TestDesignMode:
    def test_basic(self, locus_text):
        r = run({"mode": "design", "template": locus_text, "preset": "genotyping-gel"})
        assert r["ok"]
        assert r["pairs"]
        assert r["stats"]["template_length"] == 3000

    def test_targets_as_a_string(self, locus_text):
        r = run({"mode": "design", "template": locus_text, "targets": "1400-1460"})
        assert r["ok"]
        for pair in r["pairs"]:
            assert int(pair["forward_pos"].split("-")[0]) < 1400
            assert int(pair["reverse_pos"].split("-")[1]) > 1460

    def test_targets_as_a_list(self, locus_text):
        r = run({"mode": "design", "template": locus_text, "targets": ["1400-1460"]})
        assert r["ok"]

    def test_excludes(self, locus_text):
        r = run({"mode": "design", "template": locus_text, "excludes": "1-1500"})
        assert r["ok"]
        for pair in r["pairs"]:
            assert int(pair["forward_pos"].split("-")[0]) >= 1500

    def test_impossible_constraints_report_a_diagnosis(self, locus_text):
        r = run(
            {
                "mode": "design",
                "template": locus_text,
                "overrides": {"min_tm": 85, "opt_tm": 87, "max_tm": 90},
            }
        )
        assert r["ok"]
        assert r["pairs"] == []
        assert r["diagnosis"]
        assert "tm" in r["diagnosis"].lower()

    def test_template_too_short(self):
        r = run({"mode": "design", "template": "ACGT" * 10})
        assert not r["ok"]
        assert "shorter than the minimum" in r["error"]


class TestCheckMode:
    def test_named_primers(self):
        r = run({"mode": "check", "primers": "F=ACGTTGCATGCCAATGGCTA\nR=TAGCCATTGGCATGCAACGT"})
        assert r["ok"]
        titles = [t["title"] for t in r["tables"]]
        assert "Primer properties" in titles
        assert "Pairwise interactions" in titles
        assert r["tables"][0]["rows"][0][0] == "F"

    def test_unnamed_and_comma_separated(self):
        r = run({"mode": "check", "primers": "ACGTTGCATGCCAATGGCTA,TAGCCATTGGCATGCAACGT"})
        assert r["ok"]
        assert len(r["tables"][0]["rows"]) == 2

    def test_fasta_input(self):
        r = run({"mode": "check", "primers": ">F1\nACGTTGCATGCCAATGGCTA\n>R1\nTAGCCATTGGCATGCAACGT"})
        assert r["ok"]
        assert r["tables"][0]["rows"][0][0] == "F1"

    def test_panel_cross_reactivity_is_reported(self):
        target = make_seq(2000, 808)
        primer = target[500:520]
        r = run({"mode": "check", "primers": primer, "avoid": ">paralog\n" + target})
        assert r["ok"]
        panel_table = next(t for t in r["tables"] if "avoid-panel" in t["title"])
        assert panel_table["rows"][0][2] == "YES"

    def test_empty_input(self):
        r = run({"mode": "check", "primers": "  "})
        assert not r["ok"]
        assert "at least one primer" in r["error"]

    def test_invalid_sequence_is_reported_readably(self):
        r = run({"mode": "check", "primers": "F=NOTASEQUENCE!!"})
        assert not r["ok"]
        assert "error" in r


class TestContract:
    def test_run_never_raises_on_nonsense(self):
        for payload in (
            {},
            {"mode": "nope"},
            {"mode": "genotype"},
            {"mode": "design", "template": "???"},
            {"mode": "qpcr", "template": ">a\nACGT", "junctions": ["x"]},
            {"mode": "dcaps", "template": ">a\n" + make_seq(500, 9), "variant_position": "abc"},
            {"mode": "design", "template": ">a\n" + make_seq(900, 4), "targets": "not-a-range"},
        ):
            result = run(payload)
            assert isinstance(result, dict)
            assert result["ok"] is False or result["ok"] is True
            if not result["ok"]:
                assert result["error"]
                assert "Traceback" not in result["error"]

    def test_unknown_mode_lists_the_valid_ones(self):
        r = run({"mode": "wat"})
        assert not r["ok"]
        for mode in MODES:
            assert mode in r["error"]

    def test_accepts_a_json_string(self, locus_text):
        r = run(json.dumps({"mode": "design", "template": locus_text}))
        assert r["ok"]

    def test_run_json_returns_valid_json(self, locus_text):
        text = run_json({"mode": "design", "template": locus_text})
        data = json.loads(text)
        assert data["ok"]
        assert data["pairs"]

    def test_every_successful_result_has_the_full_shape(self, locus_text):
        r = run({"mode": "design", "template": locus_text})
        for key in ("ok", "mode", "report", "pairs", "warnings", "stats", "tsv", "json", "tables"):
            assert key in r, key
        json.loads(r["json"])
        assert r["tsv"].startswith("pair\t")

    def test_results_are_json_serialisable(self, locus_text, with_protospacer):
        text, spacer = with_protospacer
        for payload in (
            {"mode": "design", "template": locus_text},
            {"mode": "genotype", "template": text, "protospacer": spacer},
            {"mode": "check", "primers": "ACGTTGCATGCCAATGGCTA"},
        ):
            json.dumps(run(payload))


class TestUiMetadata:
    def test_riboprobe_schemes_are_exposed(self):
        meta = ui_metadata()
        names = {s["name"] for s in meta["riboprobe_schemes"]}
        assert {"riboprobe-t7", "riboprobe-t3", "riboprobe-sp6"} <= names
        assert all(s["description"] for s in meta["riboprobe_schemes"])

    def test_shape(self):
        meta = ui_metadata()
        assert meta["version"]
        assert meta["modes"] == list(MODES)
        assert {p["name"] for p in meta["presets"]} >= {"qpcr", "genotyping-gel"}
        assert any(n["name"] == "SpCas9" for n in meta["nucleases"])
        assert set(meta["defaults"]) == set(TUNABLE)

    def test_metadata_is_json_serialisable(self):
        json.dumps(ui_metadata())
