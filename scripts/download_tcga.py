"""Download TCGA STAR counts (open access) from GDC and build expression + clinical tables.

Outputs in --outdir:
  expression.tsv   log2(TPM+1), protein-coding genes x patients (one primary tumour per patient)
  gene_map.tsv     identity map (genes are already symbols)
  clinical_<ENDPOINT>.tsv  from TCGA-CDR: sample_id, time, event, endpoint, age, stage
  manifest.tsv     GDC files used
Example:
  python scripts/download_tcga.py --project TCGA-LUAD --cancer LUAD --endpoint PFI --outdir data/raw/TCGA_LUAD
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.acquire import (
    GDC_DATA_ENDPOINT,
    GDC_FILES_ENDPOINT,
    cdr_clinical_table,
    gdc_file_table,
    gdc_star_counts_filters,
    log2p1,
    one_sample_per_patient,
    parse_star_counts,
)

# GDC PanCanAtlas page, file "TCGA-CDR-SupplementalTableS1.xlsx". If this id stops working,
# download the file manually from https://gdc.cancer.gov/about-data/publications/pancanatlas
# and pass --cdr-file.
TCGA_CDR_URL = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"


def _post_json(url: str, payload: dict, timeout: int = 600) -> bytes:
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def query_files(project: str) -> pd.DataFrame:
    payload = {
        "filters": gdc_star_counts_filters(project),
        "fields": "file_id,file_name,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type",
        "format": "JSON",
        "size": "5000",
    }
    hits = json.loads(_post_json(GDC_FILES_ENDPOINT, payload))["data"]["hits"]
    return gdc_file_table(hits)


RETRYABLE = (urllib.error.URLError, ConnectionError, TimeoutError, OSError, tarfile.TarError, EOFError)


def with_retries(fn, *args, attempts: int = 5, base_wait: float = 5.0, sleep=time.sleep):
    """Call fn(*args), retrying on network/archive errors with exponential backoff."""
    for attempt in range(1, attempts + 1):
        try:
            return fn(*args)
        except RETRYABLE as exc:
            if attempt == attempts:
                raise
            wait = base_wait * 2 ** (attempt - 1)
            print(f"  network error ({type(exc).__name__}); retry {attempt}/{attempts - 1} in {wait:.0f}s", flush=True)
            sleep(wait)


def download_one(file_id: str, timeout: int = 300) -> dict[str, bytes]:
    """Single file via a plain GET (lighter for the GDC server than a POST bundle)."""
    with urllib.request.urlopen(f"{GDC_DATA_ENDPOINT}/{file_id}", timeout=timeout) as resp:
        return {file_id: resp.read()}


def download_batch(file_ids: list[str]) -> dict[str, bytes]:
    if len(file_ids) == 1:
        return download_one(file_ids[0])
    raw = _post_json(GDC_DATA_ENDPOINT, {"ids": file_ids})
    out: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith(".tsv"):
                out[member.name.split("/")[0]] = tar.extractfile(member).read()
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--cancer", required=True, help="TCGA-CDR 'type' value, e.g. LUAD")
    parser.add_argument("--endpoint", default="PFI")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cdr-file", help="local TCGA-CDR-SupplementalTableS1.xlsx (skips download)")
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--passes", type=int, default=5, help="extra passes over files that failed")
    parser.add_argument("--pass-wait", type=int, default=60, help="seconds to wait between passes")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = one_sample_per_patient(with_retries(query_files, args.project))
    files.to_csv(outdir / "manifest.tsv", sep="\t", index=False)
    print(f"{len(files)} primary-tumour files (one per patient)")

    # Every file is cached on disk, so an interrupted run resumes where it stopped.
    cache = outdir / "star_counts"
    cache.mkdir(exist_ok=True)
    ids = files["file_id"].tolist()
    # GDC often drops connections under load. A failing file never aborts the run: it is
    # skipped and retried in a later pass; everything already cached is kept.
    for n_pass in range(1, args.passes + 1):
        pending = [f for f in ids if not (cache / f"{f}.tsv").exists()]
        if not pending:
            break
        if n_pass > 1:
            print(f"pass {n_pass}/{args.passes}: {len(pending)} files left; waiting {args.pass_wait}s", flush=True)
            time.sleep(args.pass_wait)
        else:
            print(f"{len(ids) - len(pending)} already cached, {len(pending)} to download")
        for i in range(0, len(pending), args.batch_size):
            batch = pending[i:i + args.batch_size]
            try:
                got = with_retries(download_batch, batch, attempts=3)
            except RETRYABLE:
                print("  batch keeps failing; downloading its files one by one", flush=True)
                got = {}
                for file_id in batch:
                    try:
                        got.update(with_retries(download_one, file_id, attempts=3))
                    except RETRYABLE:
                        print(f"  skipping {file_id} for now", flush=True)
            for file_id, content in got.items():
                (cache / f"{file_id}.tsv").write_bytes(content)
            print(f"downloaded {min(i + args.batch_size, len(pending))}/{len(pending)}", flush=True)
    missing = [f for f in ids if not (cache / f"{f}.tsv").exists()]
    if missing:
        raise SystemExit(f"{len(missing)} files still missing; rerun the same command to resume")

    series = {}
    for file_id, patient in zip(files["file_id"], files["patient_barcode"]):
        series[patient] = parse_star_counts(cache / f"{file_id}.tsv")
    expression = log2p1(pd.DataFrame(series))
    expression.index.name = "gene"

    cdr_path = Path(args.cdr_file) if args.cdr_file else outdir / "TCGA-CDR-SupplementalTableS1.xlsx"
    if not cdr_path.exists():
        urllib.request.urlretrieve(TCGA_CDR_URL, cdr_path)
    xls = pd.ExcelFile(cdr_path)
    sheet = "TCGA-CDR" if "TCGA-CDR" in xls.sheet_names else xls.sheet_names[0]
    clinical = cdr_clinical_table(xls.parse(sheet), args.cancer, args.endpoint)

    # Expression keeps every patient; one clinical table per endpoint (clinical_PFI.tsv,
    # clinical_OS.tsv, ...). Samples are aligned later, per analysis, by the workflow.
    # Re-running with another --endpoint reuses the cached STAR files (no download).
    expression.to_csv(outdir / "expression.tsv", sep="\t")
    pd.DataFrame({"feature_id": expression.index, "gene_id": expression.index}).to_csv(outdir / "gene_map.tsv", sep="\t", index=False)
    clinical = clinical[clinical["sample_id"].isin(expression.columns)].reset_index(drop=True)
    clinical.to_csv(outdir / f"clinical_{args.endpoint.upper()}.tsv", sep="\t", index=False)
    print(
        f"final: {expression.shape[0]} genes x {expression.shape[1]} patients with expression; "
        f"{len(clinical)} with {args.endpoint.upper()} data, {int(clinical['event'].sum())} events"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
