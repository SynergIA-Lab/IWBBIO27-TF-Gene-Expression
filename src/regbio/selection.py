from __future__ import annotations

import numpy as np
import pandas as pd


def assert_no_validation_labels(metadata: dict[str, str]) -> None:
    if metadata.get("label_source") == "validation" or metadata.get("dataset_role") == "validation":
        raise ValueError("validation labels cannot be used for feature selection")


def cox_score_statistics(X: pd.DataFrame, y: pd.DataFrame) -> pd.Series:
    """Univariate Cox score-test statistic (Breslow ties) for every column of X.

    For each feature, the score test of H0: beta = 0 in a univariate Cox model is
    U^2 / I, with U = sum over events of (x_i - mean of x in the risk set) and
    I = sum over events of the variance of x in the risk set. It is asymptotically
    chi-square(1) and equals the log-rank statistic for a binary covariate.
    Computed for all features at once with reverse cumulative sums.
    """
    time = y["time"].astype(float).to_numpy()
    event = y["event"].astype(bool).to_numpy()
    if not event.any():
        raise ValueError("score statistics require at least one event")
    order = np.argsort(time, kind="mergesort")
    t_sorted = time[order]
    x = X.astype(float).to_numpy()[order]
    e = event[order]

    s0 = np.arange(len(t_sorted), 0, -1, dtype=float)  # risk-set size from each index onwards
    s1 = np.cumsum(x[::-1], axis=0)[::-1]
    s2 = np.cumsum((x * x)[::-1], axis=0)[::-1]
    # risk set of a subject = everyone with time >= its time -> first index of its tie group
    first = np.searchsorted(t_sorted, t_sorted, side="left")[e]
    mean = s1[first] / s0[first, None]
    var = s2[first] / s0[first, None] - mean * mean
    U = (x[e] - mean).sum(axis=0)
    info = np.clip(var, 0.0, None).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = np.where(info > 1e-12, U * U / info, 0.0)
    return pd.Series(stat, index=X.columns, name="cox_score_chi2")


def select_top_k_univariate(X: pd.DataFrame, y: pd.DataFrame, k: int) -> list[str]:
    if k <= 0:
        raise ValueError("k must be positive")
    if k > X.shape[1]:
        raise ValueError("k cannot exceed number of features")
    scores = cox_score_statistics(X, y)
    ranked = sorted(scores.index, key=lambda feature: (-scores[feature], feature))
    return ranked[:k]
