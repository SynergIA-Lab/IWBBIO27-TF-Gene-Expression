from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.io import read_clinical, read_expression
from regbio.survival import fit_survival_models


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--label-source", default="TCGA")
    parser.add_argument("--primary-l1-ratio", type=float, default=0.05)
    parser.add_argument("--lasso-l1-ratio", type=float, default=0.95)
    parser.add_argument("--cv-splits", type=int, default=5)
    parser.add_argument("--rsf-n-estimators", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    features = read_expression(args.features).T
    clinical = read_clinical(args.clinical).set_index("sample_id")
    common = [s for s in features.index if s in clinical.index]
    features, clinical = features.loc[common], clinical.loc[common]
    panel = pd.read_csv(args.panel, sep="\t")
    artifact = fit_survival_models(
        features,
        clinical,
        selected_features=panel["feature"].tolist(),
        label_source=args.label_source,
        output_path=args.output,
        cv_splits=args.cv_splits,
        primary_l1_ratio=args.primary_l1_ratio,
        lasso_l1_ratio=args.lasso_l1_ratio,
        random_state=args.seed,
        rsf_n_estimators=args.rsf_n_estimators,
    )
    for name in ("primary", "lasso_sensitivity"):
        if artifact[name]["alpha_selection"] != "cv":
            print(f"WARNING: {name} penalty used {artifact[name]['alpha_selection']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
