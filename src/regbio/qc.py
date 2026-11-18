from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from regbio.config import validate_endpoint_alignment


def characterize_expression(expression: pd.DataFrame) -> dict[str, Any]:
    return {
        "n_features": int(expression.shape[0]),
        "n_samples": int(expression.shape[1]),
        "missing_values": int(expression.isna().sum().sum()),
    }


def characterize_survival(clinical: pd.DataFrame, min_confirmatory_events: int = 50) -> dict[str, Any]:
    if (clinical["time"] <= 0).any():
        raise ValueError("survival time must be positive")
    events = int(clinical["event"].sum())
    if events == 0:
        raise ValueError("zero events in clinical table")
    n = int(clinical.shape[0])
    return {
        "n_patients": n,
        "n_events": events,
        "n_censored": n - events,
        "censoring_fraction": float((n - events) / n),
        "analysis_status": "confirmatory_eligible" if events >= min_confirmatory_events else "exploratory",
    }


def endpoint_status(training_endpoint: str, validation_endpoint: str) -> str:
    return validate_endpoint_alignment(training_endpoint, validation_endpoint)["status"]


def summarize_composition(clinical: pd.DataFrame) -> dict[str, int]:
    summary: dict[str, int] = {}
    for column in ("stage", "subtype"):
        if column in clinical.columns:
            counts = clinical[column].fillna("missing").astype(str).value_counts()
            for value, count in counts.items():
                summary[f"{column}:{value}"] = int(count)
    return summary


def validate_sample_alignment(expression: pd.DataFrame, clinical: pd.DataFrame) -> None:
    expr_samples = set(expression.columns.astype(str))
    clinical_samples = set(clinical["sample_id"].astype(str))
    missing_expr = clinical_samples - expr_samples
    missing_clinical = expr_samples - clinical_samples
    if missing_expr or missing_clinical:
        raise ValueError(
            "sample mismatch between expression and clinical tables: "
            f"missing expression={sorted(missing_expr)}, missing clinical={sorted(missing_clinical)}"
        )


def write_qc_summary(summary: dict[str, Any], path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(path, sep="\t", index=False)
