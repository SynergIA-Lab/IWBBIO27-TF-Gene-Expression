"""Kuncheva stability of the fixed top-k panel under bootstrap resampling of TCGA."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.io import read_clinical, read_expression, write_table
from regbio.metrics import mean_pairwise_kuncheva
from regbio.selection import assert_no_validation_labels, select_top_k_univariate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--analysis", default="fixture")
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--n-resamples", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--label-source", default="TCGA")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    assert_no_validation_labels({"label_source": args.label_source})
    X = read_expression(args.features).T
    clinical = read_clinical(args.clinical).set_index("sample_id")
    common = [s for s in X.index if s in clinical.index]
    X, clinical = X.loc[common], clinical.loc[common]
    k = min(args.k, X.shape[1])
    rng = np.random.default_rng(args.seed)
    panels = []
    for _ in range(args.n_resamples):
        idx = rng.integers(0, len(X), size=len(X))
        yb = clinical.iloc[idx]
        if yb["event"].sum() < 2:
            continue
        panels.append(set(select_top_k_univariate(X.iloc[idx].reset_index(drop=True), yb.reset_index(drop=True), k)))
    write_table(
        pd.DataFrame([{
            "analysis": args.analysis,
            "representation": args.representation,
            "k": k,
            "universe_size": X.shape[1],
            "n_resamples": len(panels),
            "kuncheva_mean": mean_pairwise_kuncheva(panels, X.shape[1]),
        }]),
        args.output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
