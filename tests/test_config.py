from pathlib import Path

import pytest

from regbio.config import load_analysis_config, load_cohort_config, validate_endpoint_alignment


def test_analysis_config_defaults(tmp_path: Path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text("feature_panel_size: 25\n", encoding="utf-8")
    loaded = load_analysis_config(cfg)
    assert loaded["feature_panel_size"] == 25
    assert loaded["harmonization"] == "cohort_zscore"
    assert loaded["min_confirmatory_events"] == 50
    assert loaded["primary_comparisons"] == [
        ["collectri_ulm", "gene_expression"],
        ["collectri_ulm", "tf_expression"],
    ]


def test_project_config_loads():
    loaded = load_analysis_config("config/config.yaml")
    assert loaded["cv_splits"] >= 2
    load_cohort_config("config/cohorts.yaml")


def test_cohort_config_requires_endpoint(tmp_path: Path):
    cohorts = tmp_path / "cohorts.yaml"
    cohorts.write_text("cohorts:\n  TCGA_BRCA:\n    cancer: BRCA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="endpoint"):
        load_cohort_config(cohorts)


def test_endpoint_alignment_marks_os_to_pfi_exploratory():
    result = validate_endpoint_alignment("OS", "PFI")
    assert result["status"] == "exploratory"
