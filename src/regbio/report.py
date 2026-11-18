from __future__ import annotations

from pathlib import Path

import pandas as pd


SECTIONS = [
    "Data",
    "Inference",
    "Primary Analysis (pooled, Holm)",
    "Confirmatory Comparisons",
    "Exploratory Analyses",
    "Interpretation",
    "Limitations",
]


def _table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        cells = [f"{v:.3f}" if isinstance(v, float) else str(v) for v in row.tolist()]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def build_markdown_report(
    qc: pd.DataFrame | None = None,
    evaluation: pd.DataFrame | None = None,
    comparisons: pd.DataFrame | None = None,
    stability: pd.DataFrame | None = None,
    meta: pd.DataFrame | None = None,
    is_fixture: bool = True,
) -> str:
    lines = ["# Regulatory Biomarker Transfer Summary", ""]
    for section in SECTIONS:
        lines.extend([f"## {section}", ""])
        if section == "Data" and qc is not None:
            keep = [c for c in ["cohort", "platform_class", "endpoint", "n_features", "n_patients", "n_events", "censoring_fraction", "analysis_status"] if c in qc.columns]
            lines.append(_table(qc[keep]))
        elif section == "Inference" and evaluation is not None:
            lines.append(_table(evaluation))
            if stability is not None:
                lines.extend(["", "Panel stability (mean pairwise Kuncheva over TCGA bootstrap resamples):", "", _table(stability)])
        elif section == "Primary Analysis (pooled, Holm)" and meta is not None and not meta.empty:
            keep = [c for c in ["comparison", "scope", "role", "k_cohorts", "fe_estimate", "fe_ci_lower", "fe_ci_upper", "fe_p", "fe_p_holm", "re_estimate", "re_p", "i2"] if c in meta.columns]
            lines.append("Fixed-effect inverse-variance pooling of paired C-index differences; Holm across the pre-specified comparisons (role = primary). Per-cohort rows below are descriptive.")
            lines.extend(["", _table(meta[keep])])
        elif section == "Confirmatory Comparisons" and comparisons is not None:
            lines.append("Per-cohort results for the pre-specified comparisons (descriptive; inference is the pooled primary analysis).")
            lines.append("")
            conf = comparisons[comparisons["comparison_status"] == "confirmatory"]
            lines.append(_table(conf) if not conf.empty else "No confirmatory comparisons available.")
        elif section == "Exploratory Analyses" and comparisons is not None:
            expl = comparisons[comparisons["comparison_status"] == "exploratory"]
            lines.append(_table(expl) if not expl.empty else "No exploratory comparisons.")
        elif section == "Limitations":
            if is_fixture:
                lines.append("FIXTURE RUN: synthetic data, results test the machinery only and do not support biological conclusions.")
            lines.append("Cohort-wise z-scoring, endpoint mismatches, low event counts, ARACNe leakage flags and any 'fallback' penalty selection must be reviewed before making claims.")
        else:
            lines.append("No results for this section.")
        lines.append("")
    return "\n".join(lines)


def write_report(markdown: str, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(markdown, encoding="utf-8")
