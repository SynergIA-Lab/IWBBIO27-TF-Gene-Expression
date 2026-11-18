from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


DEFAULTS: dict[str, Any] = {
    "feature_panel_size": 25,
    "harmonization": "cohort_zscore",
    "min_regulon_targets": 5,
    "min_confirmatory_events": 50,
    "random_seed": 20261026,
    "time_horizons_days": [365, 1095, 1825],
    "primary_l1_ratio": 0.05,
    "lasso_l1_ratio": 0.95,
    "primary_comparisons": [
        ["collectri_ulm", "gene_expression"],
        ["collectri_ulm", "tf_expression"],
    ],
    "representations": [
        "gene_expression",
        "tf_expression",
        "collectri_ulm",
        "dorothea_viper",
        "aracne_activity",
    ],
}

PROGRESSION_ENDPOINTS = {"PFI", "RFS", "DFS", "DMFS"}


def _read_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML file must contain a mapping: {path}")
    return data


def load_analysis_config(path: str | Path) -> dict[str, Any]:
    cfg = dict(DEFAULTS)
    cfg.update(_read_yaml(path))
    if cfg["harmonization"] not in {"cohort_zscore", "sample_rank"}:
        raise ValueError("harmonization must be 'cohort_zscore' or 'sample_rank'")
    if cfg["feature_panel_size"] <= 0:
        raise ValueError("feature_panel_size must be positive")
    if cfg["min_confirmatory_events"] <= 0:
        raise ValueError("min_confirmatory_events must be positive")
    return cfg


def load_cohort_config(path: str | Path) -> dict[str, Any]:
    data = _read_yaml(path)
    cohorts = data.get("cohorts")
    if not isinstance(cohorts, dict):
        raise ValueError("cohort config must contain a 'cohorts' mapping")
    required = {"cancer", "role", "platform_class", "endpoint", "source"}
    for name, cohort in cohorts.items():
        if not isinstance(cohort, dict):
            raise ValueError(f"cohort {name} must be a mapping")
        missing = required - set(cohort)
        if missing:
            raise ValueError(f"cohort {name} missing required fields: {', '.join(sorted(missing))}")
    return data


def validate_endpoint_alignment(training_endpoint: str, validation_endpoint: str) -> dict[str, str]:
    train = training_endpoint.upper()
    valid = validation_endpoint.upper()
    if train == "OS" and valid == "OS":
        return {"status": "confirmatory_eligible", "reason": "OS paired with OS"}
    if train in PROGRESSION_ENDPOINTS and valid in PROGRESSION_ENDPOINTS:
        return {"status": "confirmatory_eligible", "reason": "progression-related endpoints paired"}
    return {"status": "exploratory", "reason": f"endpoint mismatch: {train} vs {valid}"}
