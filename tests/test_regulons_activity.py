import subprocess
import sys

import pandas as pd
import pytest

from regbio.activity import tf_expression_control, ulm_activity, viper_activity
from regbio.regulons import filter_regulons, load_collectri, load_dorothea, tf_control_set


def test_tf_expression_uses_filtered_collectri_tfs():
    regulon = pd.DataFrame({
        "source": ["TF1", "TF1", "TF2", "TF3"],
        "target": ["G1", "G2", "G1", "G9"],
        "weight": [1.0, 1.0, 1.0, 1.0],
    })
    filtered = filter_regulons(regulon, common_genes={"G1", "G2"}, min_targets=2)
    assert set(filtered["source"]) == {"TF1"}
    assert tf_control_set(filtered) == ["TF1"]


def test_activity_and_tf_control_share_tf_names():
    expr = pd.DataFrame({"S1": [1.0, 2.0, 3.0], "S2": [2.0, 3.0, 4.0]}, index=["TF1", "G1", "G2"])
    regulon = pd.DataFrame({"source": ["TF1", "TF1"], "target": ["G1", "G2"], "weight": [1.0, 1.0]})
    activity = ulm_activity(expr, regulon, tmin=1)
    control = tf_expression_control(expr, ["TF1"])
    assert list(activity.index) == list(control.index) == ["TF1"]


def test_ulm_activity_matches_direct_decoupler_call():
    dc = pytest.importorskip("decoupler")
    expr = pd.read_csv("tests/fixtures/expression_tcga.tsv", sep="\t", index_col=0)
    regulon = pd.read_csv("tests/fixtures/collectri_fixture.tsv", sep="\t")
    expected, _ = dc.mt.ulm(expr.T, regulon, tmin=1, verbose=False)
    observed = ulm_activity(expr, regulon, tmin=1)
    pd.testing.assert_frame_equal(observed, expected.T)


def test_viper_activity_matches_direct_decoupler_call():
    dc = pytest.importorskip("decoupler")
    expr = pd.read_csv("tests/fixtures/expression_tcga.tsv", sep="\t", index_col=0)
    regulon = pd.read_csv("tests/fixtures/collectri_fixture.tsv", sep="\t")
    expected, _ = dc.mt.viper(expr.T, regulon, tmin=1, verbose=False)
    observed = viper_activity(expr, regulon, tmin=1)
    pd.testing.assert_frame_equal(observed, expected.T)


def test_build_features_tf_control_matches_collectri_tf_set(tmp_path):
    pytest.importorskip("decoupler")
    expr = pd.read_csv("tests/fixtures/expression_tcga.tsv", sep="\t", index_col=0).drop(index="TF3")
    expr.to_csv(tmp_path / "expr.tsv", sep="\t")
    prefix = tmp_path / "f"
    subprocess.run(
        [sys.executable, "scripts/build_features.py", "--expression", str(tmp_path / "expr.tsv"),
         "--collectri", "tests/fixtures/collectri_fixture.tsv", "--dorothea", "tests/fixtures/dorothea_fixture.tsv",
         "--min-targets", "4", "--out-prefix", str(prefix), "--target-counts", str(tmp_path / "c.tsv")],
        check=True,
    )
    ulm = pd.read_csv(f"{prefix}.collectri_ulm.tsv", sep="\t", index_col=0)
    tfx = pd.read_csv(f"{prefix}.tf_expression.tsv", sep="\t", index_col=0)
    assert sorted(ulm.index) == sorted(tfx.index)
    assert "TF3" not in ulm.index  # TF whose own expression is unmeasured is dropped from both
    assert not (tmp_path / "f.aracne_activity.tsv").exists()  # no silent ARACNe placeholder


def test_dorothea_loader_requests_confidence_abc(monkeypatch):
    calls = {}

    def fake_dorothea(**kwargs):
        calls.update(kwargs)
        return pd.DataFrame({"source": ["TF1"], "target": ["G1"], "weight": [1.0]})

    import decoupler as dc

    monkeypatch.setattr(dc.op, "dorothea", fake_dorothea)
    loaded = load_dorothea()
    assert calls["levels"] == ["A", "B", "C"]
    assert list(loaded.columns) == ["source", "target", "weight"]


def test_collectri_loader_returns_core_columns(monkeypatch):
    def fake_collectri(**kwargs):
        return pd.DataFrame({"source": ["TF1"], "target": ["G1"], "weight": [1.0], "extra": ["x"]})

    import decoupler as dc

    monkeypatch.setattr(dc.op, "collectri", fake_collectri)
    loaded = load_collectri()
    assert list(loaded.columns) == ["source", "target", "weight"]
