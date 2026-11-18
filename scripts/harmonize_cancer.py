"""Harmonize all cohorts of one cancer at once.

1. map each cohort to gene symbols (duplicate probes/genes averaged);
2. keep only samples present in that cohort's clinical table(s) (e.g. drops normal tissue);
3. drop lowly expressed genes, judged on the training cohort only;
4. lock the common gene universe across ALL cohorts of the cancer;
5. cohort-wise z-score (or per-sample ranks), each cohort on its own patients.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.harmonize import aggregate_duplicates, cohort_zscore, filter_low_expression, lock_common_genes, sample_rank, subset_samples
from regbio.io import read_clinical, read_expression, write_expression


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", nargs=4, action="append", metavar=("NAME", "EXPRESSION", "GENE_MAP", "CLINICALS"),
                        required=True, help="CLINICALS is a comma-separated list of clinical tables")
    parser.add_argument("--training", required=True, help="cohort NAME used for the expression filter")
    parser.add_argument("--min-expression", type=float, default=1.0)
    parser.add_argument("--min-fraction", type=float, default=0.2)
    parser.add_argument("--method", choices=["cohort_zscore", "sample_rank"], default="cohort_zscore")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    cohorts, log = {}, {"cohorts": {}}
    for name, expr_path, map_path, clinicals in args.cohort:
        expr = aggregate_duplicates(read_expression(expr_path), pd.read_csv(map_path, sep="\t"))
        samples = set().union(*(set(read_clinical(p)["sample_id"]) for p in clinicals.split(",")))
        expr = subset_samples(expr, samples)
        cohorts[name] = expr
        log["cohorts"][name] = {"genes_mapped": int(expr.shape[0]), "samples": int(expr.shape[1])}
    train = filter_low_expression(cohorts[args.training], args.min_expression, args.min_fraction)
    log["training_genes_after_expression_filter"] = int(train.shape[0])
    cohorts[args.training] = train
    locked = lock_common_genes(cohorts)
    log["common_gene_universe"] = int(next(iter(locked.values())).shape[0])
    transformed = cohort_zscore(locked) if args.method == "cohort_zscore" else sample_rank(locked)
    out = Path(args.outdir)
    for name, df in transformed.items():
        write_expression(df, out / f"{name}.expression.tsv")
    (out / "harmonization_log.json").write_text(json.dumps(log, indent=2))
    print(json.dumps(log, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
