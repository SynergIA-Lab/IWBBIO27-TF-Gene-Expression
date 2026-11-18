"""Download a GEO series and write expression, probe->gene map and raw phenotype tables.

Outputs in --outdir:
  expression_probes.tsv        probes x samples (log2 applied if values look raw; see acquisition_log.json)
  gene_map.tsv                 probe -> gene symbol (multi-gene probes dropped)
  pheno_raw.tsv                GEOparse phenotype table (all characteristics)
  characteristics_overview.tsv one row per phenotype column with examples -> use it to fill config/clinical_mappings.yaml
  acquisition_log.json
Example:
  python scripts/download_geo.py --gse GSE31210 --outdir data/raw/GSE31210
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.acquire import characteristics_overview, gpl_gene_map, needs_log2


def soft_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz"


def fetch_soft(gse: str, outdir: Path, attempts: int = 5, base_wait: float = 10.0) -> Path:
    """Download the SOFT family file ourselves (resumable cache, retries with backoff).

    GEOparse's own downloader hides network errors behind a misleading
    'FileNotFoundError ... ID could be incorrect' message and does not retry.
    """
    dest = outdir / f"{gse}_family.soft.gz"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"cached {dest.name}")
        return dest
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(1, attempts + 1):
        try:
            print(f"downloading {soft_url(gse)}", flush=True)
            with urllib.request.urlopen(soft_url(gse), timeout=300) as resp, tmp.open("wb") as out:
                shutil.copyfileobj(resp, out, length=1 << 20)
            with gzip.open(tmp, "rb") as check:  # reject truncated archives
                while check.read(1 << 24):
                    pass
            tmp.rename(dest)
            return dest
        except (urllib.error.URLError, ConnectionError, TimeoutError, OSError, EOFError) as exc:
            if attempt == attempts:
                raise SystemExit(f"could not download {soft_url(gse)}: {exc}")
            wait = base_wait * 2 ** (attempt - 1)
            print(f"  network error ({type(exc).__name__}: {exc}); retry {attempt}/{attempts - 1} in {wait:.0f}s", flush=True)
            time.sleep(wait)
    return dest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--symbol-column", help="override the platform column holding gene symbols")
    args = parser.parse_args()
    import GEOparse

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    soft = fetch_soft(args.gse, outdir)
    gse = GEOparse.get_GEO(filepath=str(soft), silent=True)
    if len(gse.gpls) != 1:
        raise ValueError(f"{args.gse} has {len(gse.gpls)} platforms ({list(gse.gpls)}); handle them separately")
    gpl_name, gpl = next(iter(gse.gpls.items()))

    expr = gse.pivot_samples("VALUE").apply(pd.to_numeric, errors="coerce")
    expr.index = expr.index.astype(str)
    log2_applied = needs_log2(expr)
    if log2_applied:
        expr = np.log2(expr.clip(lower=1.0))
    expr.index.name = "feature_id"

    gene_map, symbol_col = gpl_gene_map(gpl.table, args.symbol_column)
    pheno = gse.phenotype_data
    expr.to_csv(outdir / "expression_probes.tsv", sep="\t")
    gene_map.to_csv(outdir / "gene_map.tsv", sep="\t", index=False)
    pheno.to_csv(outdir / "pheno_raw.tsv", sep="\t")
    characteristics_overview(pheno).to_csv(outdir / "characteristics_overview.tsv", sep="\t", index=False)
    log = {
        "gse": args.gse,
        "platform": gpl_name,
        "n_probes": int(expr.shape[0]),
        "n_samples": int(expr.shape[1]),
        "symbol_column": symbol_col,
        "n_probes_mapped": int(gene_map.shape[0]),
        "n_genes_mapped": int(gene_map["gene_id"].nunique()),
        "log2_applied": log2_applied,
    }
    (outdir / "acquisition_log.json").write_text(json.dumps(log, indent=2))
    print(json.dumps(log, indent=2))
    print(f"Now inspect {outdir / 'characteristics_overview.tsv'} and fill config/clinical_mappings.yaml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
