"""How a tailed oligo is displayed: lower case for the tail, plus a colour.

The two carry the same information on purpose. These tests hold both halves of
that: the case must always be there, and the colour must never be the only
thing that is -- and must never leak into output that is not a terminal.
"""

from __future__ import annotations

import csv
import io
import re

import pytest
from helpers import make_seq

from danioprime.assays.tails import (
    FALLBACK_STYLES,
    TAIL_SCHEMES,
    TAILS,
    order_sheet,
    register_tail,
    tail_style,
)
from danioprime.cli import main
from danioprime.design.candidate import Primer
from danioprime.report.colour import ANSI_RE, paint, should_colour, strip_ansi
from danioprime.report.writers import (
    TAIL_LEGEND,
    format_text,
    render_oligo,
    render_order_sequence,
    write_order_sheet,
    write_tsv,
)
from danioprime.sequence import SeqRecord, annotate_tail, write_fasta

HEX = re.compile(r"^#[0-9a-f]{6}$")


def tailed_primer(tail_name: str = "M13F") -> Primer:
    tail = TAILS[tail_name]
    return Primer(
        seq="ACGTACGTACGTACGTACGT",
        start=0,
        end=20,
        orientation="F",
        tm=60.0,
        gc=50.0,
        tail=tail.seq,
        tail_name=tail.name,
    )


class TestAnnotation:
    def test_tail_is_lower_and_annealing_portion_is_upper(self):
        p = tailed_primer()
        assert p.annotated_seq == "tgtaaaacgacggccagtACGTACGTACGTACGTACGT"

    @pytest.mark.parametrize("name", sorted(TAILS))
    def test_only_the_case_differs_from_the_ordered_oligo(self, name):
        """The bases must be identical; nothing here may edit a sequence."""
        p = tailed_primer(name)
        assert p.annotated_seq.upper() == p.ordered_seq
        assert len(p.annotated_seq) == len(p.ordered_seq)

    @pytest.mark.parametrize("name", sorted(TAILS))
    def test_the_lower_case_run_is_exactly_the_tail(self, name):
        """What a reader infers from the case has to be the truth.

        Were a tail ever to end in a base that is also lower case in the
        annealing portion, or were the join off by one, the reader would be
        told the wrong boundary -- so check the boundary, not just the string.
        """
        p = tailed_primer(name)
        n = len(p.tail)
        assert p.annotated_seq[:n].islower()
        assert p.annotated_seq[n:].isupper()
        assert p.annotated_seq[:n].upper() == p.tail

    def test_an_untailed_primer_is_unchanged(self):
        p = Primer(seq="ACGTACGTACGTACGTACGT", start=0, end=20, orientation="F")
        assert p.annotated_seq == p.ordered_seq == "ACGTACGTACGTACGTACGT"
        assert p.tail_length == 0

    def test_annotate_tail_normalises_mixed_input(self):
        assert annotate_tail("TGTaaa", "acgT") == "tgtaaaACGT"

    def test_to_dict_carries_both_forms(self):
        d = tailed_primer().to_dict()
        assert d["ordered_seq"] == d["annotated_seq"].upper()


class TestPalette:
    @pytest.mark.parametrize("name", sorted(TAILS))
    def test_every_registered_tail_has_a_valid_style(self, name):
        style = tail_style(name)
        assert TAILS[name].style is not None, "registered tails get a chosen colour"
        assert HEX.match(style.light) and HEX.match(style.dark)
        assert 0 <= style.ansi <= 255

    def test_registered_colours_are_all_different(self):
        """Two tails the same colour would be worse than no colour at all."""
        styles = [tail_style(n) for n in TAILS]
        for field in ("light", "dark", "ansi"):
            values = [getattr(s, field) for s in styles]
            assert len(set(values)) == len(values), f"duplicate {field}"

    def test_the_three_promoters_differ(self):
        """The costly confusion: the wrong polymerase makes an unusable probe."""
        promoters = [n for n, t in TAILS.items() if t.is_promoter]
        assert len(promoters) >= 3
        assert len({tail_style(n).ansi for n in promoters}) == len(promoters)

    def test_an_unknown_name_gets_a_colour_rather_than_an_error(self):
        style = tail_style("no-such-tail")
        assert style in FALLBACK_STYLES

    def test_a_custom_tail_keeps_its_colour_across_processes(self):
        """A hash, not registration order, so the colour is reproducible."""
        before = tail_style("my-weird-tail")
        register_tail("my-weird-tail", "ACGTACGT")
        try:
            assert tail_style("my-weird-tail") == before
        finally:
            TAILS.pop("my-weird-tail", None)

    def test_fallback_styles_are_valid_and_distinct(self):
        assert len({s.ansi for s in FALLBACK_STYLES}) == len(FALLBACK_STYLES)
        for s in FALLBACK_STYLES:
            assert HEX.match(s.light) and HEX.match(s.dark)


class TestShouldColour:
    class Tty:
        def isatty(self):
            return True

    class NotTty:
        def isatty(self):
            return False

    def test_never_ignores_everything_else(self):
        assert not should_colour("never", self.Tty(), {"TERM": "xterm", "FORCE_COLOR": "1"})

    def test_always_even_without_a_terminal(self):
        assert should_colour("always", None, {"NO_COLOR": "1"})

    def test_auto_on_a_terminal(self):
        assert should_colour("auto", self.Tty(), {"TERM": "xterm-256color"})

    def test_auto_off_when_redirected(self):
        assert not should_colour("auto", self.NotTty(), {"TERM": "xterm"})
        assert not should_colour("auto", None, {"TERM": "xterm"})

    def test_no_color_wins_over_force_color(self):
        env = {"TERM": "xterm", "NO_COLOR": "", "FORCE_COLOR": "1"}
        assert not should_colour("auto", self.Tty(), env)

    def test_force_color_covers_a_pipe(self):
        assert should_colour("auto", None, {"TERM": "xterm", "FORCE_COLOR": "1"})
        assert not should_colour("auto", None, {"TERM": "xterm", "FORCE_COLOR": "0"})

    def test_dumb_terminal(self):
        assert not should_colour("auto", self.Tty(), {"TERM": "dumb"})

    def test_a_stream_that_raises_is_treated_as_not_a_terminal(self):
        class Broken:
            def isatty(self):
                raise OSError("closed")

        assert not should_colour("auto", Broken(), {"TERM": "xterm"})

    def test_an_unknown_mode_is_an_error_not_a_guess(self):
        with pytest.raises(ValueError):
            should_colour("maybe", None, {})


class TestPaint:
    def test_round_trip(self):
        painted = paint("acgt", 33)
        assert painted != "acgt"
        assert strip_ansi(painted) == "acgt"

    def test_nothing_happens_when_colour_is_off(self):
        assert paint("acgt", 33, colour=False) == "acgt"

    def test_an_empty_string_gets_no_escapes(self):
        """Otherwise a stray reset code lands in the output on its own."""
        assert paint("", 33) == ""


class TestRenderOligo:
    def test_plain_form_has_no_escapes(self):
        text = render_oligo(tailed_primer())
        assert not ANSI_RE.search(text)
        assert text.startswith("tgtaaaacgacggccagt")
        assert "(+M13F)" in text

    def test_colour_adds_only_colour(self):
        primer = tailed_primer()
        assert strip_ansi(render_oligo(primer, colour=True)) == render_oligo(primer)

    def test_the_tail_and_its_name_share_one_colour(self):
        """The name is the legend; a different colour would break the link."""
        text = render_oligo(tailed_primer("T7"), colour=True)
        codes = set(ANSI_RE.findall(text))
        ansi = tail_style("T7").ansi
        assert f"\x1b[38;5;{ansi}m" in codes
        assert len(codes) == 2  # the colour and the reset

    def test_an_untailed_primer_is_never_painted(self):
        primer = Primer(seq="ACGTACGTACGTACGTACGT", start=0, end=20, orientation="F")
        assert render_oligo(primer, colour=True) == "ACGTACGTACGTACGTACGT"

    def test_a_tail_with_no_name_still_renders(self):
        primer = Primer(
            seq="ACGTACGTACGTACGTACGT", start=0, end=20, orientation="F", tail="GGGG"
        )
        assert render_oligo(primer, colour=True).startswith("\x1b[")
        assert strip_ansi(render_oligo(primer, colour=True)) == "ggggACGTACGTACGTACGTACGT"


@pytest.fixture
def tailed_result(locus, params):
    from danioprime.assays.tails import apply_scheme
    from danioprime.design.pipeline import DesignRequest, run_design

    result = run_design(DesignRequest(template=locus, params=params, name="g"))
    assert result.pairs
    for pair in result.pairs:
        apply_scheme(pair, "fluorescent")
    return result


class TestTextReport:
    def test_no_escapes_by_default(self, tailed_result):
        """The browser shows this string in a <pre>; escapes would be garbage."""
        assert not ANSI_RE.search(format_text(tailed_result))

    def test_colour_adds_only_colour(self, tailed_result):
        coloured = format_text(tailed_result, colour=True)
        assert ANSI_RE.search(coloured)
        assert strip_ansi(coloured) == format_text(tailed_result)

    def test_the_legend_appears_once(self, tailed_result):
        assert format_text(tailed_result, show=3).count(TAIL_LEGEND) == 1

    def test_no_legend_without_tails(self, locus, params):
        from danioprime.design.pipeline import DesignRequest, run_design

        result = run_design(DesignRequest(template=locus, params=params))
        assert TAIL_LEGEND not in format_text(result)

    def test_both_lengths_are_given_for_a_tailed_primer(self, tailed_result):
        """The oligo on screen is longer than the primer that anneals."""
        text = format_text(tailed_result, show=1)
        assert "38 nt, 20 anneal" in text or re.search(r"\d+ nt, \d+ anneal", text)
        pair = tailed_result.pairs[0]
        total = len(pair.forward.ordered_seq)
        assert f"{total} nt, {pair.forward.length} anneal" in text

    def test_the_legend_is_short_enough_to_read(self):
        assert len(TAIL_LEGEND) <= 90


class TestTables:
    def test_pair_tsv_marks_the_tail(self, tailed_result):
        rows = list(csv.DictReader(io.StringIO(write_tsv(tailed_result.pairs)), delimiter="\t"))
        pair = tailed_result.pairs[0]
        assert rows[0]["forward_seq"] == pair.forward.annotated_seq
        assert rows[0]["forward_seq"].upper() == pair.forward.ordered_seq

    def test_order_sheet_keeps_an_upper_case_column_to_order_from(self, tailed_result):
        items = order_sheet(tailed_result.pairs, scheme="fluorescent")
        text = write_order_sheet(items)
        rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
        assert rows, "the sheet must not be empty"
        for row in rows:
            assert row["sequence"] == row["sequence"].upper(), "the ordering column"
            assert row["sequence_annotated"].upper() == row["sequence"], "same bases"

    def test_the_companion_oligo_reads_as_tail_throughout(self, tailed_result):
        items = order_sheet(tailed_result.pairs, scheme="fluorescent")
        companions = [i for i in items if i.role == "companion"]
        assert companions, "the labelled M13F oligo is the one people forget"
        for item in companions:
            assert item.annotated == item.sequence.lower()

    def test_render_order_sequence_colours_only_the_tail(self, tailed_result):
        items = order_sheet(tailed_result.pairs, scheme="fluorescent")
        for item in items:
            plain = render_order_sequence(item)
            assert not ANSI_RE.search(plain)
            assert strip_ansi(render_order_sequence(item, colour=True)) == plain
            if item.role != "companion":
                head = len(plain) - len(item.annealing_seq)
                painted = render_order_sequence(item, colour=True)
                assert painted.startswith("\x1b[")
                assert plain[head:] in painted.split("\x1b[0m")[-1]


class TestThroughTheCli:
    @pytest.fixture
    def locus(self, tmp_path):
        path = tmp_path / "locus.fa"
        write_fasta([SeqRecord(name="locus", seq=make_seq(2500, 11))], path)
        return path

    def test_captured_stdout_is_never_coloured(self, locus, capsys):
        """pytest's capture is not a tty, which is the common CI case."""
        assert main(["design", "--template", str(locus), "--tails", "fluorescent"]) == 0
        out = capsys.readouterr().out
        assert not ANSI_RE.search(out)
        assert "tgtaaaacgacggccagt" in out

    def test_a_report_written_to_a_file_is_never_coloured(self, locus, tmp_path):
        out = tmp_path / "report.txt"
        assert main(
            ["design", "--template", str(locus), "--tails", "fluorescent", "--out", str(out)]
        ) == 0
        assert not ANSI_RE.search(out.read_text())

    def test_color_always_is_honoured(self, locus, capsys):
        assert main(
            ["design", "--template", str(locus), "--tails", "fluorescent", "--color", "always"]
        ) == 0
        assert ANSI_RE.search(capsys.readouterr().out)

    def test_color_never_wins(self, locus, capsys):
        assert main(
            ["design", "--template", str(locus), "--tails", "fluorescent", "--color", "never"]
        ) == 0
        assert not ANSI_RE.search(capsys.readouterr().out)

    def test_colour_spelling_is_accepted_too(self, locus, capsys):
        assert main(
            ["design", "--template", str(locus), "--tails", "fluorescent", "--colour", "always"]
        ) == 0
        assert ANSI_RE.search(capsys.readouterr().out)

    def test_tails_command_lists_the_palette(self, capsys):
        assert main(["tails", "--color", "always"]) == 0
        out = capsys.readouterr().out
        assert "tgtaaaacgacggccagt" in out, "sequences read as tails"
        for name in TAILS:
            assert f"\x1b[38;5;{tail_style(name).ansi}m{name}\x1b[0m" in out

    def test_tails_command_explains_a_missing_palette(self, capsys):
        assert main(["tails"]) == 0
        out = capsys.readouterr().out
        assert not ANSI_RE.search(out)
        assert "lower case" in out.lower()


class TestWebPayload:
    def test_rows_carry_what_the_page_needs(self):
        from danioprime.webapp import run

        result = run(
            {
                "mode": "design",
                "template": make_seq(1500, 13),
                "tails": "fluorescent",
            }
        )
        assert result["ok"], result.get("error")
        row = result["pairs"][0]
        assert row["forward_tail"] == "M13F"
        assert row["reverse_tail"] == "pigtail"
        for prefix in ("forward", "reverse"):
            seq, n = row[prefix], row[f"{prefix}_tail_len"]
            assert n > 0
            assert seq[:n].islower() and seq[n:].isupper()
            assert HEX.match(row[f"{prefix}_tail_light"])
            assert HEX.match(row[f"{prefix}_tail_dark"])
            assert row[f"{prefix}_len"] == len(seq)
            assert row[f"{prefix}_anneal_len"] == len(seq) - n

    def test_untailed_rows_say_so_without_colour_keys(self):
        from danioprime.webapp import run

        result = run({"mode": "design", "template": make_seq(1500, 13)})
        row = result["pairs"][0]
        assert row["forward_tail"] == "" and row["forward_tail_len"] == 0
        assert "forward_tail_light" not in row

    def test_the_embedded_report_is_plain_text(self):
        from danioprime.webapp import run

        result = run(
            {"mode": "design", "template": make_seq(1500, 13), "tails": "fluorescent"}
        )
        assert not ANSI_RE.search(result["report"])


class TestSchemesAreAllRenderable:
    @pytest.mark.parametrize("scheme", sorted(TAIL_SCHEMES))
    def test_every_scheme_renders_and_keeps_its_bases(self, scheme, tailed_result):
        from danioprime.assays.tails import apply_scheme

        pair = tailed_result.pairs[0]
        apply_scheme(pair, scheme)
        for primer in (pair.forward, pair.reverse):
            assert primer.annotated_seq.upper() == primer.ordered_seq
            assert strip_ansi(render_oligo(primer, colour=True)) == render_oligo(primer)
