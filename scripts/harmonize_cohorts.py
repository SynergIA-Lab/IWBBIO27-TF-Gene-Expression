from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.harmonize import aggregate_duplicates, cohort_zscore, lock_common_genes, sample_rank
from regbio.io import read_expression, write_expression


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tcga-expression", required=True)
    parser.add_argument("--geo-expression", required=True)
    parser.add_argument("--gene-map", help="one map for both cohorts (fixtures)")
    parser.add_argument("--tcga-gene-map")
    parser.add_argument("--geo-gene-map")
    parser.add_argument("--method", choices=["cohort_zscore", "sample_rank"], default="cohort_zscore")
    parser.add_argument("--tcga-output", required=True)
    parser.add_argument("--geo-output", required=True)
    args = parser.parse_args()

    tcga_map = args.tcga_gene_map or args.gene_map
    geo_map = args.geo_gene_map or args.gene_map
    if not (tcga_map and geo_map):
        parser.error("give --gene-map, or both --tcga-gene-map and --geo-gene-map")
    cohorts = {
        "tcga": aggregate_duplicates(read_expression(args.tcga_expression), pd.read_csv(tcga_map, sep="\t")),
        "geo": aggregate_duplicates(read_expression(args.geo_expression), pd.read_csv(geo_map, sep="\t")),
    }
    print(f"common gene universe: {len(set(cohorts['tcga'].index) & set(cohorts['geo'].index))} genes")
    locked = lock_common_genes(cohorts)
    transformed = cohort_zscore(locked) if args.method == "cohort_zscore" else sample_rank(locked)
    write_expression(transformed["tcga"], args.tcga_output)
    write_expression(transformed["geo"], args.geo_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
