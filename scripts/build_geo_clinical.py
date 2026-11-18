"""Apply config/clinical_mappings.yaml to a GEO phenotype table.

Writes clinical.tsv and an expression matrix restricted to the samples kept, so that
expression columns and clinical sample_ids match exactly (required by the QC step).
Refuses to run on DRAFT mappings unless --allow-draft is given.
Example:
  python scripts/build_geo_clinical.py --gse GSE31210 --rawdir data/raw/GSE31210
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.acquire import apply_clinical_mapping


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse", required=True)
    parser.add_argument("--rawdir", required=True)
    parser.add_argument("--mappings", default="config/clinical_mappings.yaml")
    parser.add_argument("--allow-draft", action="store_true")
    args = parser.parse_args()
    raw = Path(args.rawdir)

    mapping = yaml.safe_load(Path(args.mappings).read_text())[args.gse]
    if str(mapping.get("status", "")).upper().startswith("DRAFT") and not args.allow_draft:
        raise SystemExit(f"{args.gse} mapping is DRAFT: verify it against characteristics_overview.tsv, set status: VERIFIED, then rerun")
    # dtype=str: a characteristic split over several columns contains NaN, so pandas would
    # read '0'/'1' as floats ('0.0'/'1.0') and no filter or event label would match.
    pheno = pd.read_csv(raw / "pheno_raw.tsv", sep="\t", index_col=0, dtype=str)
    clinical = apply_clinical_mapping(pheno, mapping)
    n_dropped = clinical.attrs.get("n_dropped", 0)
    if clinical.empty:
        raise SystemExit(f"{args.gse}: no patient left after sample_filter/time/event; check the values in clinical_mappings.yaml against characteristics_overview.tsv")

    expr = pd.read_csv(raw / "expression_probes.tsv", sep="\t", index_col=0)
    keep = [s for s in clinical["sample_id"] if s in expr.columns]
    clinical = clinical.set_index("sample_id").loc[keep].reset_index()
    expr[keep].to_csv(raw / "expression.tsv", sep="\t")
    clinical.to_csv(raw / "clinical.tsv", sep="\t", index=False)
    print(
        f"{args.gse}: {len(clinical)} patients, {int(clinical['event'].sum())} {mapping['endpoint']} events, "
        f"{n_dropped} dropped for missing/unknown time or event; "
        f"stage counts: {clinical['stage'].value_counts(dropna=False).to_dict()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
