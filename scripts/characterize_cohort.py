from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.io import read_clinical, read_expression
from regbio.qc import characterize_expression, characterize_survival, summarize_composition, validate_sample_alignment, write_qc_summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--endpoint", required=True)
    parser.add_argument("--min-events", type=int, default=50)
    parser.add_argument("--subset-to-clinical", action="store_true", help="restrict expression to patients with clinical data")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    expression = read_expression(args.expression)
    clinical = read_clinical(args.clinical)
    if args.subset_to_clinical:
        common = [s for s in expression.columns if s in set(clinical["sample_id"])]
        expression = expression[common]
        clinical = clinical[clinical["sample_id"].isin(common)]
    validate_sample_alignment(expression, clinical)
    summary = {
        "cohort": args.cohort,
        "platform_class": args.platform,
        "endpoint": args.endpoint,
        **characterize_expression(expression),
        **characterize_survival(clinical, args.min_events),
        **summarize_composition(clinical),
    }
    write_qc_summary(summary, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
