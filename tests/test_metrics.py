import numpy as np
import pandas as pd
import pytest

from regbio.metrics import (
    comparison_status,
    cumulative_dynamic_auc_scores,
    ipcw_cindex,
    paired_bootstrap_cindex_difference,
    valid_time_horizons,
)
from regbio.survival import make_survival_y


def _cohort(n, seed):
    rng = np.random.default_rng(seed)
    signal = rng.normal(size=n)
    t = rng.exponential(1000 / np.exp(signal))
    c = rng.uniform(100, 2000, size=n)
    clin = pd.DataFrame({"time": np.maximum(np.minimum(t, c), 1), "event": (t <= c).astype(int)})
    return signal, make_survival_y(clin)


def test_invalid_horizon_removed():
    times = valid_time_horizons([365, 5000], train_max_event_time=2000, validation_max_event_time=1500)
    assert times == [365]


def test_paired_bootstrap_identical_risks_gives_zero():
    pytest.importorskip("sksurv")
    _, y_train = _cohort(150, 1)
    signal, y_valid = _cohort(120, 2)
    res = paired_bootstrap_cindex_difference(y_train, y_valid, signal, signal.copy(), n_boot=50, seed=1)
    assert res["difference"] == 0.0
    assert res["ci_lower"] == 0.0 and res["ci_upper"] == 0.0


def test_paired_bootstrap_detects_informative_vs_noise():
    pytest.importorskip("sksurv")
    _, y_train = _cohort(200, 3)
    signal, y_valid = _cohort(200, 4)
    noise = np.random.default_rng(9).normal(size=len(signal))
    res = paired_bootstrap_cindex_difference(y_train, y_valid, signal, noise, n_boot=200, seed=1)
    assert res["cindex_a"] > res["cindex_b"]
    assert res["ci_lower"] > 0
    assert res["difference"] == pytest.approx(res["cindex_a"] - res["cindex_b"])


def test_primary_comparison_status():
    primary = [["collectri_ulm", "gene_expression"], ["collectri_ulm", "tf_expression"]]
    assert comparison_status("collectri_ulm", "gene_expression", primary) == "confirmatory"
    assert comparison_status("dorothea_viper", "gene_expression", primary) == "exploratory"


def test_ipcw_cindex_matches_direct_sksurv_call():
    pytest.importorskip("sksurv")
    from sksurv.metrics import concordance_index_ipcw

    train = make_survival_y(pd.DataFrame({"event": [1, 1, 0, 1, 0], "time": [2, 4, 5, 6, 7]}))
    valid = make_survival_y(pd.DataFrame({"event": [1, 0, 1], "time": [3, 4, 5]}))
    risk = np.array([0.8, 0.2, 0.6])
    assert ipcw_cindex(train, valid, risk) == concordance_index_ipcw(train, valid, risk)[0]


def test_dynamic_auc_matches_direct_sksurv_call():
    pytest.importorskip("sksurv")
    from sksurv.metrics import cumulative_dynamic_auc

    train = make_survival_y(pd.DataFrame({"event": [1, 1, 0, 1, 0], "time": [2, 4, 5, 6, 7]}))
    valid = make_survival_y(pd.DataFrame({"event": [1, 0, 1], "time": [3, 4, 5]}))
    risk = np.array([0.8, 0.2, 0.6])
    observed = cumulative_dynamic_auc_scores(train, valid, risk, [3])
    expected_auc, expected_mean = cumulative_dynamic_auc(train, valid, risk, [3])
    assert np.allclose(observed["auc"], expected_auc)
    assert observed["mean_auc"] == expected_mean


def test_holm_adjust_matches_definition():
    from regbio.metrics import holm_adjust

    assert holm_adjust([0.25, 0.035]) == pytest.approx([0.25, 0.07])
    assert holm_adjust([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])


def test_inverse_variance_meta_fixed_effect_and_i2():
    from regbio.metrics import inverse_variance_meta

    m = inverse_variance_meta([0.004, 0.044], [0.031, 0.030])
    assert m["fe_estimate"] == pytest.approx(0.0247, abs=1e-4)
    assert m["fe_p"] == pytest.approx(0.253, abs=1e-3)
    same = inverse_variance_meta([0.05, 0.05, 0.05], [0.02, 0.03, 0.04])
    assert same["i2"] == 0.0 and same["tau2"] == 0.0 and same["fe_estimate"] == pytest.approx(0.05)
