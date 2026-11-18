import numpy as np
import pandas as pd
import pytest

from regbio.survival import fit_survival_models, make_survival_y, model_plan, predict_risk_scores, validate_penalty_cv_source


def _informative(n=80, seed=0):
    rng = np.random.default_rng(seed)
    X = pd.DataFrame(rng.normal(size=(n, 3)), columns=["F1", "F2", "F3"], index=[f"S{i}" for i in range(n)])
    t = rng.exponential(1000 / np.exp(X["F1"].to_numpy()))
    c = rng.uniform(100, 2500, size=n)
    clinical = pd.DataFrame({"time": np.maximum(np.minimum(t, c), 1), "event": (t <= c).astype(int)}, index=X.index)
    return X, clinical


def test_survival_target_fields():
    clinical = pd.DataFrame({"event": [1, 0], "time": [10.0, 20.0]})
    y = make_survival_y(clinical)
    assert y.dtype.names == ("event", "time")


def test_penalty_cv_rejects_validation_source():
    with pytest.raises(ValueError, match="TCGA"):
        validate_penalty_cv_source("GEO")


def test_primary_model_plan_uses_low_l1():
    plan = model_plan(primary_l1_ratio=0.05)
    assert plan["primary"]["l1_ratio"] <= 0.10
    assert plan["lasso_sensitivity"]["l1_ratio"] >= 0.90


def test_informative_data_uses_cv_not_fallback(tmp_path):
    pytest.importorskip("sksurv")
    X, clinical = _informative()
    artifact = fit_survival_models(X, clinical, ["F1", "F2", "F3"], "TCGA", tmp_path / "m.joblib", cv_splits=5, rsf_n_estimators=20)
    assert artifact["primary"]["alpha_selection"] == "cv"
    assert artifact["lasso_sensitivity"]["alpha_selection"] == "cv"
    assert artifact["primary"]["n_informative_folds"] >= 2
    assert (tmp_path / "m.joblib").exists()


def test_ridge_primary_keeps_full_panel(tmp_path):
    pytest.importorskip("sksurv")
    X, clinical = _informative()
    artifact = fit_survival_models(X, clinical, ["F1", "F2", "F3"], "TCGA", tmp_path / "m.joblib", alpha_grid=[0.01], rsf_n_estimators=20)
    assert np.count_nonzero(artifact["primary"]["estimator"].coef_) == 3


def test_noninformative_folds_fall_back_and_are_recorded(tmp_path):
    pytest.importorskip("sksurv")
    X = pd.DataFrame({"F1": [0.0, 0.2, 0.4, 0.6], "F2": [1.0, 0.9, 0.7, 0.4]}, index=[f"S{i}" for i in range(4)])
    clinical = pd.DataFrame({"event": [1, 0, 1, 1], "time": [100, 200, 150, 250]}, index=X.index)
    artifact = fit_survival_models(X, clinical, ["F1", "F2"], "TCGA", tmp_path / "s.joblib", cv_splits=2, alpha_grid=[0.01, 0.1], rsf_n_estimators=10)
    assert artifact["primary"]["alpha_selection"] == "fallback_no_informative_folds"
    assert artifact["primary"]["alpha"] in {0.01, 0.1}


def test_fit_rejects_validation_penalty_source(tmp_path):
    pytest.importorskip("sksurv")
    X, clinical = _informative(20)
    with pytest.raises(ValueError, match="TCGA"):
        fit_survival_models(X, clinical, ["F1"], "GEO", tmp_path / "bad.joblib")


def test_predict_risk_scores_from_frozen_artifact(tmp_path):
    pytest.importorskip("sksurv")
    X, clinical = _informative()
    artifact = fit_survival_models(X, clinical, ["F1", "F2"], "TCGA", tmp_path / "m.joblib", rsf_n_estimators=20)
    risks = predict_risk_scores(artifact, X)
    assert set(risks) == {"primary", "lasso_sensitivity", "rsf"}
    assert all(len(v) == len(X) for v in risks.values())
