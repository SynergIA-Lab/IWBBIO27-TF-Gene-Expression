"""Download a GEO series whose expression lives in a supplementary matrix (e.g. RNA-seq).

Use it when download_geo.py cannot: series with several platforms or without per-sample
VALUE tables, such as GSE96058 (SCAN-B, log2(FPKM+0.1) per gene symbol).

Phenotypes come from every '<GSE>*_series_matrix.txt.gz' file of the series; the
expression matrix is the supplementary file matching --suppl-pattern, whose columns are
sample titles and are mapped to GSM accessions. Writes the same outputs as
download_geo.py, so build_geo_clinical.py and the workflow work unchanged:
  expression_probes.tsv, gene_map.tsv (identity: rows are gene symbols), pheno_raw.tsv,
  characteristics_overview.tsv, acquisition_log.json
Example:
  python scripts/download_geo_suppl.py --gse GSE96058 --outdir data/raw/GSE96058 \
      --suppl-pattern "gene_expression.*transformed.*csv" --drop-columns-regex "repl"
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.acquire import characteristics_overview, links_in_listing, map_columns_to_gsm, needs_log2, parse_series_matrix_pheno


def series_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}"


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"cached {dest.name}")
        return dest
    print(f"downloading {url}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    return dest


def listing(url: str) -> str:
    with urllib.request.urlopen(url + "/", timeout=60) as resp:
        return resp.read().decode("utf-8", "replace")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--suppl-pattern", required=True, help="regex selecting the expression file in suppl/")
    parser.add_argument("--drop-columns-regex", help="regex of expression columns to drop (e.g. technical replicates)")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    base = series_url(args.gse)

    matrices = links_in_listing(listing(base + "/matrix"), r"_series_matrix\.txt\.gz$")
    if not matrices:
        raise SystemExit(f"no series matrix files listed at {base}/matrix/")
    phenos = []
    for name in matrices:
        path = fetch(f"{base}/matrix/{name}", outdir / name)
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            phenos.append(parse_series_matrix_pheno(handle))
    pheno = pd.concat(phenos, axis=0, sort=False)
    pheno = pheno[~pheno.index.duplicated(keep="first")]

    suppl = links_in_listing(listing(base + "/suppl"), args.suppl_pattern)
    if len(suppl) != 1:
        raise SystemExit(f"--suppl-pattern must match exactly one file in {base}/suppl/, matched: {suppl}")
    expr_path = fetch(f"{base}/suppl/{suppl[0]}", outdir / suppl[0])
    sep = "\t" if ".tsv" in suppl[0] or ".txt" in suppl[0] else ","
    expr = pd.read_csv(expr_path, sep=sep, index_col=0)
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.groupby(level=0).mean()  # guard against repeated symbols
    expr, stats = map_columns_to_gsm(expr, pheno, args.drop_columns_regex)
    log2_applied = needs_log2(expr)
    if log2_applied:
        expr = np.log2(expr.clip(lower=0) + 1.0)
    expr.index.name = "feature_id"

    expr.to_csv(outdir / "expression_probes.tsv", sep="\t")
    pd.DataFrame({"feature_id": expr.index, "gene_id": expr.index}).to_csv(outdir / "gene_map.tsv", sep="\t", index=False)
    pheno.to_csv(outdir / "pheno_raw.tsv", sep="\t")
    characteristics_overview(pheno).to_csv(outdir / "characteristics_overview.tsv", sep="\t", index=False)
    log = {"gse": args.gse, "series_matrices": matrices, "expression_file": suppl[0], "n_genes": int(expr.shape[0]),
           "n_pheno_samples": int(len(pheno)), **stats, "log2_applied": log2_applied}
    (outdir / "acquisition_log.json").write_text(json.dumps(log, indent=2))
    print(json.dumps(log, indent=2))
    print(f"Now inspect {outdir / 'characteristics_overview.tsv'} and fill config/clinical_mappings.yaml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
