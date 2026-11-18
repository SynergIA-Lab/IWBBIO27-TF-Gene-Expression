from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.report import build_markdown_report, write_report


def _read_many(paths: list[str] | None) -> pd.DataFrame | None:
    if not paths:
        return None
    return pd.concat([pd.read_csv(p, sep="\t") for p in paths], ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc", nargs="*")
    parser.add_argument("--metrics", nargs="*")
    parser.add_argument("--comparisons", nargs="*")
    parser.add_argument("--stability", nargs="*")
    parser.add_argument("--meta")
    parser.add_argument("--real-data", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    md = build_markdown_report(
        qc=_read_many(args.qc),
        evaluation=_read_many(args.metrics),
        comparisons=_read_many(args.comparisons),
        stability=_read_many(args.stability),
        meta=pd.read_csv(args.meta, sep="\t") if args.meta else None,
        is_fixture=not args.real_data,
    )
    write_report(md, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
