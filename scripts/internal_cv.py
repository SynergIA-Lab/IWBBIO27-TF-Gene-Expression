"""Internal (TCGA) performance of one representation by repeated, fully nested cross-validation.

Every step that looks at survival labels is repeated inside each outer training fold, so the
held-out fold never influences its own prediction:

  1. top-k panel by univariate Cox score test        (training fold only)
  2. Coxnet penalty by inner stratified CV            (training fold only)
  3. primary Coxnet fit                               (training fold only)
  4. IPCW C-index on the held-out fold, censoring estimated on the training fold,
     truncated at the same horizon tau as the external validation (cindex_tau_days).

Reusing the panel chosen on the whole of TCGA would leak the held-out patients into the
selection and inflate the internal C-index; this script avoids that on purpose.

Outputs
  --output        one summary row: mean / SD / median of the fold C-indices
  --folds-output  optional, one row per (repeat, fold)
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.config import load_analysis_config
from regbio.io import read_clinical, read_expression, write_table
from regbio.metrics import ipcw_cindex
from regbio.selection import assert_no_validation_labels, select_top_k_univariate
from regbio.survival import DEFAULT_ALPHA_GRID, _choose_alpha, _cv_folds, make_survival_y, require_sksurv


def nested_cv(X: pd.DataFrame, clinical: pd.DataFrame, k: int, l1_ratio: float, outer_splits: int,
              inner_splits: int, repeats: int, seed: int, tau_max: float) -> pd.DataFrame:
    require_sksurv()
    from sksurv.linear_model import CoxnetSurvivalAnalysis

    y = make_survival_y(clinical)
    rows = []
    for rep in range(repeats):
        for fold, (tr, te) in enumerate(_cv_folds(y, outer_splits, seed + rep)):
            Xtr, Xte, ytr, yte = X.iloc[tr], X.iloc[te], y[tr], y[te]
            if yte["event"].sum() < 2:
                continue
            panel = select_top_k_univariate(Xtr, clinical.iloc[tr], min(k, X.shape[1]))
            choice = _choose_alpha(Xtr[panel], ytr, list(DEFAULT_ALPHA_GRID), l1_ratio, inner_splits, seed + rep)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=[choice["alpha"]]).fit(Xtr[panel], ytr)
            risk = model.predict(Xte[panel])
            # tau inside the follow-up of both folds (IPCW needs censoring support up to tau)
            tau = float(min(tau_max, ytr["time"].max() * 0.999, yte["time"].max() * 0.999))
            rows.append({"repeat": rep, "fold": fold, "n_test": len(te), "events_test": int(yte["event"].sum()),
                         "alpha": choice["alpha"], "alpha_selection": choice["alpha_selection"], "tau": tau,
                         "ipcw_cindex": ipcw_cindex(ytr, yte, risk, tau)})
    if not rows:
        raise ValueError("no informative outer folds (too few events)")
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--features", required=True, help="TCGA matrix of one representation (features x samples)")
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--analysis", required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--label-source", default="TCGA")
    parser.add_argument("--output", required=True)
    parser.add_argument("--folds-output")
    args = parser.parse_args()

    assert_no_validation_labels({"label_source": args.label_source})
    cfg = load_analysis_config(args.config)
    X = read_expression(args.features).T
    clinical = read_clinical(args.clinical).set_index("sample_id")
    common = [s for s in X.index if s in clinical.index]
    X, clinical = X.loc[common].astype(float), clinical.loc[common]

    folds = nested_cv(X, clinical, args.k, cfg["primary_l1_ratio"], cfg.get("cv_splits", 5), cfg.get("cv_splits", 5),
                      args.repeats, cfg["random_seed"], float(cfg.get("cindex_tau_days", 1825)))
    c = folds["ipcw_cindex"]
    write_table(pd.DataFrame([{
        "analysis": args.analysis, "representation": args.representation, "k": args.k,
        "n_patients": len(X), "n_events": int(clinical["event"].sum()),
        "repeats": args.repeats, "n_folds": len(folds),
        "internal_cindex_mean": float(c.mean()), "internal_cindex_sd": float(c.std(ddof=1)),
        "internal_cindex_median": float(c.median()),
        "internal_cindex_q025": float(np.percentile(c, 2.5)), "internal_cindex_q975": float(np.percentile(c, 97.5)),
    }]), args.output)
    if args.folds_output:
        write_table(folds.assign(analysis=args.analysis, representation=args.representation), args.folds_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
