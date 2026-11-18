import pandas as pd

from regbio.report import SECTIONS, build_markdown_report


def test_report_contains_required_sections():
    report = build_markdown_report()
    for section in SECTIONS:
        assert f"## {section}" in report


def test_fixture_report_is_flagged_and_splits_comparisons():
    comps = pd.DataFrame({
        "comparison": ["collectri_ulm - gene_expression", "dorothea_viper - tf_expression"],
        "comparison_status": ["confirmatory", "exploratory"],
        "difference": [0.05, -0.01],
    })
    report = build_markdown_report(comparisons=comps, is_fixture=True)
    assert "FIXTURE RUN" in report
    conf = report.split("## Confirmatory Comparisons")[1].split("## Exploratory Analyses")[0]
    assert "collectri_ulm - gene_expression" in conf and "dorothea_viper" not in conf
