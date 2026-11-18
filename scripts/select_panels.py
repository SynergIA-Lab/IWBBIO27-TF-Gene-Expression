from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.io import read_clinical, read_expression, write_table
from regbio.selection import assert_no_validation_labels, cox_score_statistics, select_top_k_univariate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--label-source", default="TCGA")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    assert_no_validation_labels({"label_source": args.label_source})
    features = read_expression(args.features).T
    clinical = read_clinical(args.clinical).set_index("sample_id")
    common = [s for s in features.index if s in clinical.index]
    features, clinical = features.loc[common], clinical.loc[common]
    k = min(args.k, features.shape[1])
    selected = select_top_k_univariate(features, clinical, k)
    scores = cox_score_statistics(features, clinical)
    write_table(
        pd.DataFrame({"feature": selected, "rank": range(1, len(selected) + 1), "cox_score_chi2": scores.loc[selected].values}),
        args.output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
