"""Apply frozen TCGA models to one validation cohort and compare representations.

Outputs
  metrics:      one row per model (each representation, the clinical baseline and each
                clinical+representation model) with IPCW C-index and mean dynamic AUC.
  comparisons:  paired bootstrap of C-index differences on the same validation patients.
                Representation pairs use the primary model; pre-specified pairs are
                confirmatory only if endpoints align and the cohort has enough events.
                'clinical+rep - clinical' rows measure added value (exploratory).

The IPCW censoring distribution is estimated on the validation cohort itself (its
follow-up differs from TCGA) and C-indices are truncated at a common horizon tau.
"""
from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.config import load_analysis_config, validate_endpoint_alignment
from regbio.io import read_clinical, read_expression, write_table
from regbio.metrics import bootstrap_cindex_ci, comparison_status, cumulative_dynamic_auc_scores, ipcw_cindex, paired_bootstrap_cindex_difference
from regbio.survival import clinical_matrix, fit_cox, load_model_artifact, make_survival_y, oof_risk, predict_risk_scores


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--representations", nargs="+", required=True)
    parser.add_argument("--models", nargs="+", required=True, help="one artifact per representation, same order")
    parser.add_argument("--train-features", nargs="+", required=True, help="TCGA matrix per representation, same order")
    parser.add_argument("--valid-features", nargs="+", required=True, help="validation matrix per representation, same order")
    parser.add_argument("--train-clinical", required=True)
    parser.add_argument("--valid-clinical", required=True)
    parser.add_argument("--train-endpoint", required=True)
    parser.add_argument("--valid-endpoint", required=True)
    parser.add_argument("--valid-qc", help="cohort QC table; its analysis_status decides confirmatory eligibility")
    parser.add_argument("--analysis", default="fixture")
    parser.add_argument("--cohort", default="fixture_geo")
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--metrics-output", required=True)
    parser.add_argument("--comparisons-output", required=True)
    parser.add_argument("--risks-output", help="per-patient primary-model risk scores + TCGA-median risk groups")
    args = parser.parse_args()
    n = len(args.representations)
    if not (len(args.models) == len(args.train_features) == len(args.valid_features) == n):
        raise ValueError("--representations, --models, --train-features and --valid-features must have the same length")

    cfg = load_analysis_config(args.config)
    seed = cfg["random_seed"]
    train_clin = read_clinical(args.train_clinical).set_index("sample_id")
    valid_clin = read_clinical(args.valid_clinical).set_index("sample_id")
    endpoint = validate_endpoint_alignment(args.train_endpoint, args.valid_endpoint)
    status = pd.read_csv(args.valid_qc, sep="\t")["analysis_status"][0] if args.valid_qc else "confirmatory_eligible"
    eligible = endpoint["status"] == "confirmatory_eligible" and status == "confirmatory_eligible"

    # validation patients = those with clinical data and features in every representation
    valid_X = {rep: read_expression(p).T for rep, p in zip(args.representations, args.valid_features)}
    patients = [s for s in valid_clin.index if all(s in X.index for X in valid_X.values())]
    valid_clin = valid_clin.loc[patients]
    y_valid = make_survival_y(valid_clin)
    tau = float(min(cfg.get("cindex_tau_days", 1825), valid_clin["time"].max() * 0.999))
    horizons = [t for t in cfg["time_horizons_days"] if t < valid_clin.loc[valid_clin["event"] == 1, "time"].max()]

    def cidx(y, risk):
        return ipcw_cindex(y, y, risk, tau)  # censoring estimated on the validation cohort

    def pair(y, a, b):
        return paired_bootstrap_cindex_difference(y, y, a, b, n_boot=args.n_boot, seed=seed, tau=tau)

    risks, artifacts = {}, {}
    for rep, model_path in zip(args.representations, args.models):
        artifacts[rep] = load_model_artifact(model_path)
        risks[rep] = predict_risk_scores(artifacts[rep], valid_X[rep].loc[patients])

    base = {"analysis": args.analysis, "validation_cohort": args.cohort, "n_patients": len(patients), "n_events": int(valid_clin["event"].sum()), "tau_days": tau}
    metric_rows = []
    for rep, r in risks.items():
        row = {**base, "model": rep}
        for m in ("primary", "lasso_sensitivity", "rsf"):
            row[f"{m}_ipcw_cindex"] = cidx(y_valid, r[m])
        row["primary_ci_lower"], row["primary_ci_upper"] = bootstrap_cindex_ci(y_valid, r["primary"], tau, args.n_boot, seed)
        row["primary_mean_dynamic_auc"] = cumulative_dynamic_auc_scores(y_valid, y_valid, r["primary"], horizons)["mean_auc"] if horizons else np.nan
        row["primary_alpha_selection"] = artifacts[rep]["primary"]["alpha_selection"]
        metric_rows.append(row)

    comp_rows = []
    for left, right in combinations(args.representations, 2):
        for p in cfg["primary_comparisons"]:
            if set(p) == {left, right}:
                left, right = p
        st = comparison_status(left, right, cfg["primary_comparisons"])
        if st == "confirmatory" and not eligible:
            st = "exploratory"
        comp_rows.append({**base, "comparison": f"{left} - {right}", "comparison_status": st, "endpoint_alignment": endpoint["reason"],
                          **pair(y_valid, risks[left]["primary"], risks[right]["primary"])})

    # clinical baseline and clinical + biomarker (complete cases for age and stage)
    tr_c = clinical_matrix(train_clin)
    va_c = clinical_matrix(valid_clin)
    if len(tr_c) >= 30 and len(va_c) >= 30 and va_c["stage"].nunique() > 1:
        cc = list(va_c.index)
        y_cc = make_survival_y(valid_clin.loc[cc])
        clin_model = fit_cox(tr_c, train_clin)
        clin_risk = clin_model.predict(va_c)
        lo, hi = bootstrap_cindex_ci(y_cc, clin_risk, tau, args.n_boot, seed)
        metric_rows.append({**base, "n_patients": len(cc), "n_events": int(valid_clin.loc[cc, "event"].sum()), "model": "clinical", "primary_ipcw_cindex": cidx(y_cc, clin_risk), "primary_ci_lower": lo, "primary_ci_upper": hi})
        for rep, train_path in zip(args.representations, args.train_features):
            art = artifacts[rep]
            Xtr = read_expression(train_path).T
            common = [s for s in tr_c.index if s in Xtr.index]
            oof = oof_risk(Xtr.loc[common], train_clin, art["selected_features"], art["primary"]["l1_ratio"], art["primary"]["alpha"], cfg.get("cv_splits", 5), seed)
            tr_comb = tr_c.loc[common].assign(biomarker=oof.loc[common].values).dropna()
            comb = fit_cox(tr_comb, train_clin)
            va_comb = va_c.assign(biomarker=pd.Series(risks[rep]["primary"], index=patients).loc[cc].values)
            comb_risk = comb.predict(va_comb)
            lo, hi = bootstrap_cindex_ci(y_cc, comb_risk, tau, args.n_boot, seed)
            metric_rows.append({**base, "n_patients": len(cc), "n_events": int(valid_clin.loc[cc, "event"].sum()), "model": f"clinical+{rep}", "primary_ipcw_cindex": cidx(y_cc, comb_risk), "primary_ci_lower": lo, "primary_ci_upper": hi})
            comp_rows.append({**base, "n_patients": len(cc), "n_events": int(valid_clin.loc[cc, "event"].sum()), "comparison": f"clinical+{rep} - clinical", "comparison_status": "exploratory", "endpoint_alignment": endpoint["reason"], **pair(y_cc, comb_risk, clin_risk)})
    else:
        print("clinical baseline skipped: not enough complete age/stage data or a single stage", file=sys.stderr)

    if args.risks_output:
        # Risk groups use a cutoff frozen on TCGA (median of the in-sample training risk),
        # never re-estimated on the validation cohort.
        risk_df = valid_clin[["time", "event"]].copy()
        for rep, train_path in zip(args.representations, args.train_features):
            Xtr = read_expression(train_path).T
            cutoff = float(np.median(predict_risk_scores(artifacts[rep], Xtr)["primary"]))
            risk_df[f"{rep}_risk"] = risks[rep]["primary"]
            risk_df[f"{rep}_high"] = (risks[rep]["primary"] > cutoff).astype(int)
        risk_df.index.name = "sample_id"
        risk_df.reset_index().assign(analysis=args.analysis, validation_cohort=args.cohort).to_csv(args.risks_output, sep="\t", index=False)

    write_table(pd.DataFrame(metric_rows), args.metrics_output)
    write_table(pd.DataFrame(comp_rows), args.comparisons_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
