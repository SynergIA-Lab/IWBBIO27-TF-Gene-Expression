from __future__ import annotations

from itertools import combinations

import numpy as np


def kuncheva_index(a: set[str], b: set[str], universe_size: int) -> float:
    if len(a) != len(b):
        raise ValueError("Kuncheva index requires sets of equal size")
    k = len(a)
    if k == 0:
        raise ValueError("Kuncheva index requires non-empty sets")
    if universe_size < k:
        raise ValueError("universe_size must be at least the set size")
    if universe_size == k:
        return 1.0 if a == b else 0.0
    r = len(a & b)
    return float((r * universe_size - k * k) / (k * (universe_size - k)))


def mean_pairwise_kuncheva(panels: list[set[str]], universe_size: int) -> float:
    if len(panels) < 2:
        raise ValueError("stability needs at least two panels")
    values = [kuncheva_index(a, b, universe_size) for a, b in combinations(panels, 2)]
    return float(np.mean(values))


def valid_time_horizons(horizons: list[int], train_max_event_time: float, validation_max_event_time: float) -> list[int]:
    upper = min(train_max_event_time, validation_max_event_time)
    return [int(t) for t in horizons if 0 < t < upper]


def default_tau(y_train: np.ndarray, y_valid: np.ndarray) -> float:
    """Truncation time for IPCW: inside the follow-up support of both cohorts."""
    return float(min(y_train["time"].max(), y_valid["time"].max()) * 0.999)


def comparison_status(left: str, right: str, primary: list[list[str]]) -> str:
    pair = [left, right]
    reverse = [right, left]
    return "confirmatory" if pair in primary or reverse in primary else "exploratory"


def ipcw_cindex(y_train: np.ndarray, y_valid: np.ndarray, risk_scores: np.ndarray, tau: float | None = None) -> float:
    from sksurv.metrics import concordance_index_ipcw

    kwargs = {"tau": tau} if tau is not None else {}
    return float(concordance_index_ipcw(y_train, y_valid, risk_scores, **kwargs)[0])


def paired_bootstrap_cindex_difference(
    y_train: np.ndarray,
    y_valid: np.ndarray,
    risk_a: np.ndarray,
    risk_b: np.ndarray,
    n_boot: int = 1000,
    seed: int = 1,
    tau: float | None = None,
) -> dict[str, float | int]:
    """Paired bootstrap of C-index(a) - C-index(b) on the same validation patients.

    Patients are resampled with replacement; both models are scored on the same
    resample, so the interval reflects the paired difference in discrimination.
    """
    risk_a = np.asarray(risk_a, dtype=float)
    risk_b = np.asarray(risk_b, dtype=float)
    if not (len(risk_a) == len(risk_b) == len(y_valid)):
        raise ValueError("risk vectors and validation outcomes must have equal length")
    if tau is None:
        tau = default_tau(y_train, y_valid)
    c_a = ipcw_cindex(y_train, y_valid, risk_a, tau)
    c_b = ipcw_cindex(y_train, y_valid, risk_b, tau)
    rng = np.random.default_rng(seed)
    n = len(y_valid)
    diffs: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y_valid[idx]
        if yb["event"].sum() < 2:
            continue
        try:
            diffs.append(ipcw_cindex(y_train, yb, risk_a[idx], tau) - ipcw_cindex(y_train, yb, risk_b[idx], tau))
        except (ValueError, ZeroDivisionError):
            continue
    arr = np.asarray(diffs, dtype=float)
    if arr.size == 0:
        raise ValueError("no valid bootstrap resamples")
    se = float(arr.std(ddof=1)) if arr.size > 1 else float("nan")
    return {
        "cindex_a": c_a,
        "cindex_b": c_b,
        "difference": c_a - c_b,
        "ci_lower": float(np.percentile(arr, 2.5)),
        "ci_upper": float(np.percentile(arr, 97.5)),
        "bootstrap_se": se,
        "p_two_sided": float(min(1.0, 2 * min((arr <= 0).mean(), (arr >= 0).mean()))),
        "n_boot_valid": int(arr.size),
        "tau": float(tau),
    }


def cumulative_dynamic_auc_scores(
    y_train: np.ndarray,
    y_valid: np.ndarray,
    risk_scores: np.ndarray,
    times: list[int] | np.ndarray,
) -> dict[str, np.ndarray | float]:
    from sksurv.metrics import cumulative_dynamic_auc

    auc, mean_auc = cumulative_dynamic_auc(y_train, y_valid, risk_scores, np.asarray(times, dtype=float))
    return {"times": np.asarray(times, dtype=float), "auc": auc, "mean_auc": float(mean_auc)}


def bootstrap_cindex_ci(y: np.ndarray, risk: np.ndarray, tau: float, n_boot: int = 1000, seed: int = 1) -> tuple[float, float]:
    """Percentile 95% CI of the IPCW C-index (censoring estimated on the same cohort)."""
    risk = np.asarray(risk, dtype=float)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), size=len(y))
        yb = y[idx]
        if yb["event"].sum() < 2:
            continue
        try:
            vals.append(ipcw_cindex(yb, yb, risk[idx], tau))
        except (ValueError, ZeroDivisionError):
            continue
    arr = np.asarray(vals)
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def holm_adjust(pvalues: list[float]) -> list[float]:
    """Holm step-down adjusted p-values (family-wise error control)."""
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[i]))
        adj[i] = running
    return adj.tolist()


def inverse_variance_meta(estimates: list[float], ses: list[float]) -> dict[str, float]:
    """Fixed-effect and DerSimonian-Laird random-effects pooling, with I^2."""
    from scipy.stats import norm

    y = np.asarray(estimates, dtype=float)
    se = np.asarray(ses, dtype=float)
    w = 1.0 / se**2
    fe = float(np.sum(w * y) / np.sum(w))
    fe_se = float(np.sqrt(1.0 / np.sum(w)))
    k = len(y)
    q = float(np.sum(w * (y - fe) ** 2))
    tau2 = max(0.0, (q - (k - 1)) / (np.sum(w) - np.sum(w**2) / np.sum(w))) if k > 1 else 0.0
    wr = 1.0 / (se**2 + tau2)
    re = float(np.sum(wr * y) / np.sum(wr))
    re_se = float(np.sqrt(1.0 / np.sum(wr)))
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0

    def p(est, s):
        return float(2 * norm.sf(abs(est / s)))

    return {
        "k_cohorts": k,
        "fe_estimate": fe, "fe_se": fe_se, "fe_ci_lower": fe - 1.96 * fe_se, "fe_ci_upper": fe + 1.96 * fe_se, "fe_p": p(fe, fe_se),
        "re_estimate": re, "re_se": re_se, "re_ci_lower": re - 1.96 * re_se, "re_ci_upper": re + 1.96 * re_se, "re_p": p(re, re_se),
        "tau2": float(tau2), "i2": float(i2), "q": q,
    }
