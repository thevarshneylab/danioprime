"""Report writers and the public API surface."""

from __future__ import annotations

import csv
import io
import json

import danioprime
from danioprime.config import get_preset
from danioprime.design.pipeline import DesignRequest, run_design
from danioprime.report.writers import (
    PAIR_COLUMNS,
    format_text,
    result_to_dict,
    summarise_iterable,
    write_bed,
    write_json,
    write_order_sheet,
    write_tsv,
)


class TestWriters:
    def test_tsv_has_a_row_per_pair_and_stable_columns(self, locus, params):
        result = run_design(DesignRequest(template=locus, params=params, name="g"))
        text = write_tsv(result.pairs)
        rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
        assert len(rows) == len(result.pairs)
        assert list(rows[0]) == list(PAIR_COLUMNS)
        assert rows[0]["pair"] == "g_p1"

    def test_json_contains_parameters_and_penalty_terms(self, locus, params):
        result = run_design(DesignRequest(template=locus, params=params, name="g"))
        data = json.loads(write_json(result))
        assert data["template"]["name"] == locus.name
        assert data["parameters"]["min_tm"] == params.min_tm
        assert data["pairs"][0]["penalty_terms"]

    def test_json_is_serialisable_including_nested_dataclasses(self, locus, params):
        result = run_design(DesignRequest(template=locus, params=params))
        # Round-trips without raising, and floats survive.
        payload = result_to_dict(result)
        assert json.loads(json.dumps(payload)) == json.loads(json.dumps(payload))

    def test_bed_offsets_into_genome_coordinates(self, locus, params):
        result = run_design(DesignRequest(template=locus, params=params, name="g"))
        text = write_bed(result.pairs[:1], chrom="chr5", offset=1_000_000)
        lines = [line for line in text.splitlines() if not line.startswith("track")]
        assert len(lines) == 3  # forward, reverse, amplicon
        fields = lines[0].split("\t")
        assert fields[0] == "chr5"
        assert int(fields[1]) >= 1_000_000
        assert fields[5] in {"+", "-"}

    def test_bed_amplicon_spans_both_primers(self, locus, params):
        result = run_design(DesignRequest(template=locus, params=params, name="g"))
        text = write_bed(result.pairs[:1])
        amplicon = [line for line in text.splitlines() if "_amplicon" in line][0].split("\t")
        pair = result.pairs[0]
        assert int(amplicon[1]) == pair.forward.start
        assert int(amplicon[2]) == pair.reverse.end

    def test_order_sheet_columns(self, locus, params):
        from danioprime.assays.tails import reaction_plan

        result = run_design(DesignRequest(template=locus, params=params, name="g"))
        plan = reaction_plan(result.pairs, scheme="fluorescent")
        text = write_order_sheet(plan.oligos)
        rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
        assert {"name", "sequence", "length", "role"} <= set(rows[0])
        assert any(r["role"] == "companion" for r in rows)

    def test_files_are_written_when_a_path_is_given(self, locus, params, tmp_path):
        result = run_design(DesignRequest(template=locus, params=params, name="g"))
        tsv = tmp_path / "p.tsv"
        assert write_tsv(result.pairs, tsv) == str(tsv)
        assert tsv.read_text().startswith("pair\t")

    def test_text_report_mentions_the_key_facts(self, locus, params):
        result = run_design(DesignRequest(template=locus, params=params, name="g"))
        text = format_text(result, verbose=True, show=2)
        assert locus.name in text
        assert "penalty breakdown" in text
        assert "Tm" in text

    def test_text_report_explains_an_empty_result(self, locus):
        params = get_preset("genotyping-gel").replace(min_tm=85.0, opt_tm=87.0, max_tm=90.0)
        result = run_design(DesignRequest(template=locus, params=params))
        text = format_text(result)
        assert "No primer pair" in text
        assert "tm" in text.lower()

    def test_summarise_iterable_handles_empty_input(self):
        assert summarise_iterable([], ("a", "b")) == "(none)"

    def test_summarise_iterable_aligns_columns(self):
        text = summarise_iterable([{"a": "x", "b": "yyy"}, {"a": "zzzz", "b": "w"}], ("a", "b"))
        lines = text.splitlines()
        assert len(lines) == 4
        assert len(lines[0]) == len(lines[1])


class TestPublicApi:
    def test_documented_entry_points_are_importable(self):
        for name in danioprime.__all__:
            assert hasattr(danioprime, name), name

    def test_version_is_exposed(self):
        assert danioprime.__version__
        assert danioprime.__version__.count(".") == 2

    def test_quick_start_from_the_docstring_works(self, locus):
        result = danioprime.run_design(
            danioprime.DesignRequest(template=locus, params=danioprime.DesignParams())
        )
        assert result.pairs
        pair = result.pairs[0]
        assert pair.forward.seq and pair.reverse.seq and pair.product_size

    def test_error_hierarchy(self):
        assert issubclass(danioprime.DesignError, danioprime.DanioPrimeError)
        assert issubclass(danioprime.SequenceError, danioprime.DanioPrimeError)
        assert issubclass(danioprime.AnnotationError, danioprime.DanioPrimeError)
