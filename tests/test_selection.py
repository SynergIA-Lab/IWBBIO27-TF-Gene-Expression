import numpy as np
import pandas as pd
import pytest

from regbio.metrics import kuncheva_index, mean_pairwise_kuncheva
from regbio.selection import assert_no_validation_labels, cox_score_statistics, select_top_k_univariate


def _brute_force_score(x, time, event):
    U = 0.0
    info = 0.0
    for i in np.where(event)[0]:
        risk = time >= time[i]
        m = x[risk].mean()
        U += x[i] - m
        info += ((x[risk] - m) ** 2).mean()
    return U * U / info


def test_score_statistic_matches_brute_force_with_ties():
    rng = np.random.default_rng(3)
    n = 40
    X = pd.DataFrame(rng.normal(size=(n, 4)), columns=list("ABCD"))
    y = pd.DataFrame({"time": rng.integers(1, 15, size=n).astype(float), "event": rng.integers(0, 2, size=n)})
    observed = cox_score_statistics(X, y)
    for col in X.columns:
        expected = _brute_force_score(X[col].to_numpy(), y["time"].to_numpy(), y["event"].astype(bool).to_numpy())
        assert observed[col] == pytest.approx(expected, rel=1e-9)


def test_selection_ranks_prognostic_feature_first():
    rng = np.random.default_rng(5)
    n = 200
    signal = rng.normal(size=n)
    X = pd.DataFrame({"signal": signal, "noise1": rng.normal(size=n), "noise2": rng.normal(size=n)})
    t = rng.exponential(1.0 / np.exp(1.2 * signal))
    y = pd.DataFrame({"time": t, "event": np.ones(n, dtype=int)})
    assert select_top_k_univariate(X, y, k=1) == ["signal"]


def test_fixed_k_selection_returns_exact_k():
    X = pd.DataFrame({f"F{i}": np.arange(8) * (i + 1) % 5 for i in range(6)}, dtype=float)
    y = pd.DataFrame({"time": [8, 7, 6, 5, 4, 3, 2, 1], "event": [1, 1, 0, 1, 1, 0, 1, 1]})
    selected = select_top_k_univariate(X, y, k=3)
    assert len(selected) == 3


def test_kuncheva_identical_sets_is_one():
    assert kuncheva_index({"A", "B"}, {"A", "B"}, universe_size=5) == 1.0


def test_mean_pairwise_kuncheva():
    panels = [{"A", "B"}, {"A", "B"}, {"A", "B"}]
    assert mean_pairwise_kuncheva(panels, universe_size=10) == 1.0


def test_validation_labels_are_rejected_for_selection():
    metadata = {"label_source": "validation"}
    with pytest.raises(ValueError, match="validation"):
        assert_no_validation_labels(metadata)
