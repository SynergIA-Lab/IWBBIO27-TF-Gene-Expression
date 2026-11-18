from __future__ import annotations

import importlib.util
import warnings
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

DEFAULT_ALPHA_GRID = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]


def make_survival_y(clinical: pd.DataFrame) -> np.ndarray:
    event = clinical["event"].astype(bool).to_numpy()
    time = clinical["time"].astype(float).to_numpy()
    if importlib.util.find_spec("sksurv") is not None:
        from sksurv.util import Surv

        return Surv.from_arrays(event=event, time=time)
    return np.array(list(zip(event, time)), dtype=[("event", "?"), ("time", "<f8")])


def validate_penalty_cv_source(source: str) -> None:
    if source not in {"TCGA", "training"}:
        raise ValueError("penalty selection must use TCGA training data only")


def model_plan(primary_l1_ratio: float = 0.05, lasso_l1_ratio: float = 0.95) -> dict[str, dict[str, float]]:
    return {
        "primary": {"model": "coxnet", "l1_ratio": float(primary_l1_ratio)},
        "lasso_sensitivity": {"model": "coxnet", "l1_ratio": float(lasso_l1_ratio)},
    }


def require_sksurv() -> None:
    if importlib.util.find_spec("sksurv") is None:
        raise RuntimeError("scikit-survival is required for real survival model fitting")


def _cv_folds(y: np.ndarray, cv_splits: int, random_state: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """Folds stratified by event status so every fold keeps events when possible."""
    from sklearn.model_selection import StratifiedKFold

    event = y["event"].astype(int)
    n_splits = int(min(cv_splits, event.sum(), len(event) - event.sum()))
    if n_splits < 2:
        return []
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    return list(skf.split(np.zeros(len(event)), event))


def _score_alpha(X: pd.DataFrame, y: np.ndarray, alpha: float, l1_ratio: float, folds) -> list[float]:
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.metrics import concordance_index_censored

    scores: list[float] = []
    for train_idx, test_idx in folds:
        y_test = y[test_idx]
        if y_test["event"].sum() < 2:
            continue  # non-informative fold: too few events to rank
        model = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=[alpha], fit_baseline_model=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                model.fit(X.iloc[train_idx], y[train_idx])
                risk = model.predict(X.iloc[test_idx])
                score = concordance_index_censored(y_test["event"], y_test["time"], risk)[0]
            except (ValueError, ArithmeticError, np.linalg.LinAlgError):
                continue  # e.g. no comparable pairs in this fold
        if np.isfinite(score):
            scores.append(float(score))
    return scores


def _choose_alpha(X: pd.DataFrame, y: np.ndarray, alpha_grid: list[float], l1_ratio: float, cv_splits: int, random_state: int) -> dict[str, Any]:
    folds = _cv_folds(y, cv_splits, random_state)
    per_alpha = {alpha: _score_alpha(X, y, alpha, l1_ratio, folds) for alpha in alpha_grid}
    valid = {alpha: float(np.mean(s)) for alpha, s in per_alpha.items() if s}
    if not valid:
        # Deterministic fallback, recorded in the artifact. Must never happen with real cohorts.
        alpha = sorted(alpha_grid)[len(alpha_grid) // 2]
        return {"alpha": float(alpha), "alpha_selection": "fallback_no_informative_folds", "cv_scores": {}, "n_informative_folds": 0}
    best = max(valid, key=lambda a: (valid[a], a))  # ties -> stronger penalty
    return {
        "alpha": float(best),
        "alpha_selection": "cv",
        "cv_scores": {str(k): v for k, v in valid.items()},
        "n_informative_folds": max(len(s) for s in per_alpha.values()),
    }


def fit_survival_models(
    X: pd.DataFrame,
    clinical: pd.DataFrame,
    selected_features: list[str],
    label_source: str,
    output_path: str | Path,
    cv_splits: int = 5,
    alpha_grid: list[float] | None = None,
    primary_l1_ratio: float = 0.05,
    lasso_l1_ratio: float = 0.95,
    random_state: int = 20261026,
    rsf_n_estimators: int = 500,
    rsf_min_samples_leaf: int = 15,
) -> dict[str, Any]:
    require_sksurv()
    validate_penalty_cv_source(label_source)
    from sksurv.ensemble import RandomSurvivalForest
    from sksurv.linear_model import CoxnetSurvivalAnalysis

    if alpha_grid is None:
        alpha_grid = list(DEFAULT_ALPHA_GRID)
    missing = [feature for feature in selected_features if feature not in X.columns]
    if missing:
        raise ValueError(f"selected features missing from X: {missing}")
    X_selected = X.loc[:, selected_features].astype(float)
    y = make_survival_y(clinical)

    primary_choice = _choose_alpha(X_selected, y, alpha_grid, primary_l1_ratio, cv_splits, random_state)
    lasso_choice = _choose_alpha(X_selected, y, alpha_grid, lasso_l1_ratio, cv_splits, random_state)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        primary = CoxnetSurvivalAnalysis(l1_ratio=primary_l1_ratio, alphas=[primary_choice["alpha"]], fit_baseline_model=True)
        primary.fit(X_selected, y)
        lasso = CoxnetSurvivalAnalysis(l1_ratio=lasso_l1_ratio, alphas=[lasso_choice["alpha"]], fit_baseline_model=True)
        lasso.fit(X_selected, y)
    leaf = int(max(1, min(rsf_min_samples_leaf, len(X_selected) // 4)))
    rsf = RandomSurvivalForest(n_estimators=rsf_n_estimators, min_samples_leaf=leaf, random_state=random_state, n_jobs=-1)
    rsf.fit(X_selected, y)

    artifact: dict[str, Any] = {
        "selected_features": list(selected_features),
        "label_source": label_source,
        "primary": {"estimator": primary, "l1_ratio": primary_l1_ratio, **primary_choice},
        "lasso_sensitivity": {
            "estimator": lasso,
            "l1_ratio": lasso_l1_ratio,
            "n_nonzero": int(np.count_nonzero(lasso.coef_)),
            **lasso_choice,
        },
        "rsf": {"estimator": rsf, "n_estimators": rsf_n_estimators, "min_samples_leaf": leaf},
    }
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(artifact, output_path)
    return artifact


def load_model_artifact(path: str | Path) -> dict[str, Any]:
    return joblib.load(path)


def predict_risk_scores(artifact: dict[str, Any], X: pd.DataFrame) -> dict[str, np.ndarray]:
    features = artifact["selected_features"]
    X_selected = X.loc[:, features].astype(float)
    return {
        "primary": np.asarray(artifact["primary"]["estimator"].predict(X_selected), dtype=float),
        "lasso_sensitivity": np.asarray(artifact["lasso_sensitivity"]["estimator"].predict(X_selected), dtype=float),
        "rsf": np.asarray(artifact["rsf"]["estimator"].predict(X_selected), dtype=float),
    }


# ---------------------------------------------------------------- clinical baseline

STAGE_ORDINAL = {"I": 1, "II": 2, "III": 3, "IV": 4}


def clinical_matrix(clinical: pd.DataFrame) -> pd.DataFrame:
    """age and ordinal stage (I-IV -> 1-4); rows with either missing are dropped."""
    df = pd.DataFrame({
        "age": pd.to_numeric(clinical.get("age"), errors="coerce"),
        "stage": clinical.get("stage").map(lambda v: STAGE_ORDINAL.get(str(v).upper(), np.nan) if pd.notna(v) else np.nan),
    }, index=clinical.index)
    return df.dropna()


def oof_risk(X: pd.DataFrame, clinical: pd.DataFrame, features: list[str], l1_ratio: float, alpha: float, cv_splits: int = 5, random_state: int = 20261026) -> pd.Series:
    """Out-of-fold linear predictor of the fixed-panel Coxnet on the training cohort.

    Used to fit the clinical+biomarker model without the optimism of in-sample risk scores.
    """
    require_sksurv()
    from sksurv.linear_model import CoxnetSurvivalAnalysis

    Xs = X.loc[:, features].astype(float)
    y = make_survival_y(clinical.loc[Xs.index])
    out = pd.Series(np.nan, index=Xs.index)
    for train_idx, test_idx in _cv_folds(y, cv_splits, random_state):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=[alpha]).fit(Xs.iloc[train_idx], y[train_idx])
        out.iloc[test_idx] = m.predict(Xs.iloc[test_idx])
    return out


def fit_cox(features: pd.DataFrame, clinical: pd.DataFrame, alpha: float = 1e-3):
    """Low-dimensional Cox model (clinical baseline or clinical + biomarker risk)."""
    require_sksurv()
    from sksurv.linear_model import CoxPHSurvivalAnalysis

    y = make_survival_y(clinical.loc[features.index])
    return CoxPHSurvivalAnalysis(alpha=alpha).fit(features.astype(float), y)
