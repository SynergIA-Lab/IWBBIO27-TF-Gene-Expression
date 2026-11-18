from pathlib import Path

import pytest

from regbio.io import read_clinical, read_expression, read_regulon


def test_read_expression_contract(tmp_path: Path):
    path = tmp_path / "expr.tsv"
    path.write_text("gene\tS1\nG1\t1.0\n", encoding="utf-8")
    expr = read_expression(path)
    assert expr.shape == (1, 1)
    assert list(expr.index) == ["G1"]


def test_read_clinical_requires_columns(tmp_path: Path):
    path = tmp_path / "clinical.tsv"
    path.write_text("sample_id\ttime\nS1\t1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="event"):
        read_clinical(path)


def test_read_regulon_adds_weight(tmp_path: Path):
    path = tmp_path / "regulon.tsv"
    path.write_text("source\ttarget\nTF1\tG1\n", encoding="utf-8")
    regulon = read_regulon(path)
    assert regulon.loc[0, "weight"] == 1.0
