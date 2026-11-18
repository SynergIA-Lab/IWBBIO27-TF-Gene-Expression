import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from regbio.harmonize import filter_low_expression, subset_samples
from regbio.survival import clinical_matrix, oof_risk


def test_filter_low_expression_on_training_only():
    expr = pd.DataFrame({"S1": [0.0, 5.0, 2.0], "S2": [0.5, 6.0, 0.0], "S3": [0.0, 7.0, 0.0]}, index=["low", "high", "sparse"])
    kept = filter_low_expression(expr, min_value=1.0, min_fraction=0.5)
    assert list(kept.index) == ["high"]


def test_subset_samples_drops_samples_without_clinical():
    expr = pd.DataFrame(np.ones((2, 3)), index=["G1", "G2"], columns=["tumor1", "normal1", "tumor2"])
    assert list(subset_samples(expr, {"tumor1", "tumor2"}).columns) == ["tumor1", "tumor2"]
    with pytest.raises(ValueError):
        subset_samples(expr, {"other"})


def test_clinical_matrix_ordinal_stage_and_complete_cases():
    clin = pd.DataFrame({"age": [60, None, 70, 55], "stage": ["I", "II", None, "IV"]}, index=list("ABCD"))
    m = clinical_matrix(clin)
    assert list(m.index) == ["A", "D"] and m["stage"].tolist() == [1, 4]


def test_oof_risk_covers_every_training_patient():
    pytest.importorskip("sksurv")
    rng = np.random.default_rng(0)
    n = 80
    X = pd.DataFrame(rng.normal(size=(n, 3)), columns=["F1", "F2", "F3"], index=[f"S{i}" for i in range(n)])
    t = rng.exponential(1000 / np.exp(X["F1"].to_numpy()))
    clin = pd.DataFrame({"time": np.maximum(t, 1), "event": rng.integers(0, 2, size=n)}, index=X.index)
    risk = oof_risk(X, clin, ["F1", "F2"], l1_ratio=0.05, alpha=0.01)
    assert risk.notna().all() and len(risk) == n


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_real_workflow_dag_builds():
    res = subprocess.run(["snakemake", "real", "-n", "--quiet", "--configfile", "tests/real_workflow_test.yaml"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-2000:]


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_real_workflow_report_labels(tmp_path):
    """Stability rows carry the analysis; clinical comparisons report complete-case events."""
    import yaml

    # validation cohort with some missing stage -> complete cases differ from the full cohort
    geo = pd.read_csv("tests/fixtures/clinical_geo.tsv", sep="\t")
    geo.loc[geo.index[:6], "stage"] = None
    geo_path = tmp_path / "clinical_geo_missing_stage.tsv"
    geo.to_csv(geo_path, sep="\t", index=False)
    cfg = yaml.safe_load(open("tests/real_workflow_test.yaml"))
    cfg["real"]["cancers"]["FIX"]["cohorts"]["FIX_GEO"]["clinical"]["RFS"] = str(geo_path)
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    res = subprocess.run(["snakemake", "real", "--cores", "2", "--forceall", "--quiet", "--configfile", str(cfg_path)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-3000:]

    stab = pd.read_csv("outputs/real/FIX_PFI/gene_expression_stability.tsv", sep="\t")
    assert stab["analysis"].tolist() == ["FIX_PFI"]
    metrics = pd.read_csv("outputs/real/FIX_PFI/FIX_GEO_metrics.tsv", sep="\t")
    comps = pd.read_csv("outputs/real/FIX_PFI/FIX_GEO_comparisons.tsv", sep="\t")
    clin_events = int(metrics.loc[metrics["model"] == "clinical", "n_events"].iloc[0])
    added = comps[comps["comparison"].str.endswith("- clinical")]
    assert not added.empty and (added["n_events"] == clin_events).all()
    assert clin_events < int(metrics.loc[metrics["model"] == "gene_expression", "n_events"].iloc[0])
