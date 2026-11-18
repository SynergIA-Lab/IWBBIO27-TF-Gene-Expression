import pandas as pd
import pytest

from regbio.qc import characterize_survival, endpoint_status, summarize_composition, validate_sample_alignment


def test_zero_events_rejected():
    clinical = pd.DataFrame({"sample_id": ["S1", "S2"], "time": [10, 20], "event": [0, 0]})
    with pytest.raises(ValueError, match="zero events"):
        characterize_survival(clinical, min_confirmatory_events=50)


def test_low_events_are_exploratory():
    clinical = pd.DataFrame({
        "sample_id": [f"S{i}" for i in range(60)],
        "time": list(range(1, 61)),
        "event": [1] * 20 + [0] * 40,
    })
    summary = characterize_survival(clinical, min_confirmatory_events=50)
    assert summary["analysis_status"] == "exploratory"


def test_endpoint_status_os_to_rfs_is_exploratory():
    assert endpoint_status("OS", "RFS") == "exploratory"


def test_stage_summary_counts_values():
    clinical = pd.DataFrame({"stage": ["I", "I", "II"], "subtype": ["A", "B", "A"]})
    summary = summarize_composition(clinical)
    assert summary["stage:I"] == 2


def test_sample_mismatch_rejected():
    expr = pd.DataFrame({"S1": [1.0]}, index=["G1"])
    clinical = pd.DataFrame({"sample_id": ["S2"], "time": [1], "event": [1]})
    with pytest.raises(ValueError, match="sample mismatch"):
        validate_sample_alignment(expr, clinical)
