"""Shared fixtures.

Templates are generated from a fixed seed rather than shipped as files so the
tests stay readable, and a helper builds the ohnolog scenario that most of the
specificity tests need.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers import make_seq, mutate

from danioprime.config import DesignParams, get_preset
from danioprime.sequence import SeqRecord, write_fasta
from danioprime.specificity.paralog import ParalogPanel

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture
def locus() -> SeqRecord:
    return SeqRecord(name="locus", seq=make_seq(2000, 1234))


@pytest.fixture
def params() -> DesignParams:
    return get_preset("genotyping-gel").replace(max_pairs=5)


@pytest.fixture
def ohnolog_pair() -> tuple[SeqRecord, SeqRecord, tuple[int, int]]:
    """A target and a paralog sharing one perfectly conserved block.

    Returns ``(target, paralog, (block_start, block_end))`` with the block in
    target coordinates. Inside the block a primer cannot discriminate; outside
    it, it can.
    """
    target_seq = make_seq(2000, 99)
    block = (800, 1000)
    paralog_seq = mutate(target_seq, 0.10, 7)
    paralog_seq = (
        paralog_seq[: block[0]] + target_seq[block[0] : block[1]] + paralog_seq[block[1] :]
    )
    return (
        SeqRecord(name="target", seq=target_seq),
        SeqRecord(name="paralog", seq=paralog_seq),
        block,
    )


@pytest.fixture
def panel(ohnolog_pair) -> ParalogPanel:
    _target, paralog, _block = ohnolog_pair
    return ParalogPanel.from_records([paralog])


@pytest.fixture
def example_dir() -> Path:
    if not (EXAMPLES / "kcnh6a_locus.fa").exists():
        pytest.skip("example data not generated; run examples/make_example_data.py")
    return EXAMPLES


@pytest.fixture
def tmp_fasta(tmp_path) -> callable:
    def _write(name: str, records: list[SeqRecord]) -> Path:
        path = tmp_path / name
        write_fasta(records, path)
        return path

    return _write


