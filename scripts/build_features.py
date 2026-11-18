from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.activity import tf_expression_control, ulm_activity, viper_activity
from regbio.io import read_expression, read_regulon, write_expression, write_table
from regbio.regulons import filter_regulons, target_counts, tf_control_set


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True, help="harmonized, common-gene-locked matrix (genes x samples)")
    parser.add_argument("--collectri", required=True)
    parser.add_argument("--dorothea", required=True)
    parser.add_argument("--aracne", help="cancer-specific ARACNe regulon; aracne_activity is only built when given")
    parser.add_argument("--min-targets", type=int, default=5)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--target-counts", required=True)
    args = parser.parse_args()

    expression = read_expression(args.expression)
    common = set(expression.index)

    collectri = filter_regulons(read_regulon(args.collectri), common, args.min_targets)
    # Primary control: the TF set must be identical for collectri_ulm and tf_expression,
    # so keep only TFs whose own expression is measured in the locked universe.
    collectri = collectri[collectri["source"].isin(common)].reset_index(drop=True)
    tf_set = tf_control_set(collectri)
    if not tf_set:
        raise ValueError("no CollecTRI TF passes the common-gene and min_regulon_targets filters")

    write_expression(expression, f"{args.out_prefix}.gene_expression.tsv")
    write_expression(tf_expression_control(expression, tf_set), f"{args.out_prefix}.tf_expression.tsv")
    write_expression(ulm_activity(expression, collectri, tmin=args.min_targets), f"{args.out_prefix}.collectri_ulm.tsv")

    dorothea = filter_regulons(read_regulon(args.dorothea), common, args.min_targets)
    write_expression(viper_activity(expression, dorothea, tmin=args.min_targets), f"{args.out_prefix}.dorothea_viper.tsv")

    if args.aracne:
        aracne = filter_regulons(read_regulon(args.aracne), common, args.min_targets)
        write_expression(viper_activity(expression, aracne, tmin=args.min_targets), f"{args.out_prefix}.aracne_activity.tsv")

    counts = pd.concat(
        [target_counts(collectri).assign(regulon="collectri"), target_counts(dorothea).assign(regulon="dorothea")],
        ignore_index=True,
    )
    write_table(counts, args.target_counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
