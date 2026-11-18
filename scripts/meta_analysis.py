"""Pre-specified confirmatory analysis (design amendment 2026-10-28, option C).

For each pre-specified comparison, pool the paired C-index differences of every
confirmatory-eligible validation cohort with inverse-variance weights (bootstrap SE).
Primary: fixed effect, Holm-adjusted across the pre-specified comparisons.
Sensitivity: DerSimonian-Laird random effects + I^2, and pooling within endpoint class.
Per-cohort results remain descriptive.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.config import PROGRESSION_ENDPOINTS, load_analysis_config
from regbio.io import write_table
from regbio.metrics import holm_adjust, inverse_variance_meta


def _endpoint_class(row: pd.Series) -> str:
    return "OS" if str(row["endpoint_alignment"]).startswith("OS") else "progression"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--comparisons", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    cfg = load_analysis_config(args.config)

    comps = pd.concat([pd.read_csv(p, sep="\t") for p in args.comparisons], ignore_index=True)
    conf = comps[comps["comparison_status"] == "confirmatory"].copy()
    conf["endpoint_class"] = conf.apply(_endpoint_class, axis=1)
    rows = []
    for left, right in cfg["primary_comparisons"]:
        label = f"{left} - {right}"
        sub = conf[conf["comparison"] == label]
        if sub.empty:
            continue
        rows.append({"comparison": label, "scope": "all_cohorts", "role": "primary",
                     "cohorts": ",".join(sub["validation_cohort"]),
                     **inverse_variance_meta(sub["difference"].tolist(), sub["bootstrap_se"].tolist())})
        for cls, g in sub.groupby("endpoint_class"):
            rows.append({"comparison": label, "scope": f"endpoint_{cls}", "role": "sensitivity",
                         "cohorts": ",".join(g["validation_cohort"]),
                         **inverse_variance_meta(g["difference"].tolist(), g["bootstrap_se"].tolist())})
    out = pd.DataFrame(rows)
    if not out.empty:
        primary = out["role"] == "primary"
        out["fe_p_holm"] = pd.NA
        out.loc[primary, "fe_p_holm"] = holm_adjust(out.loc[primary, "fe_p"].tolist())
    write_table(out, args.output)
    print(out[["comparison", "scope", "k_cohorts", "fe_estimate", "fe_ci_lower", "fe_ci_upper", "fe_p"] + (["fe_p_holm"] if "fe_p_holm" in out else [])].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
