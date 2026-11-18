"""Pure parsing helpers for public-data acquisition (TCGA via GDC, GEO via GEOparse).

Network access lives in the CLI scripts; everything here is testable offline.
"""
from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- TCGA / GDC

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"


def gdc_star_counts_filters(project: str, sample_type: str = "Primary Tumor") -> dict[str, Any]:
    def _in(field: str, values: list[str]) -> dict[str, Any]:
        return {"op": "in", "content": {"field": field, "value": values}}

    return {
        "op": "and",
        "content": [
            _in("cases.project.project_id", [project]),
            _in("data_type", ["Gene Expression Quantification"]),
            _in("analysis.workflow_type", ["STAR - Counts"]),
            _in("cases.samples.sample_type", [sample_type]),
            _in("access", ["open"]),
        ],
    }


def gdc_file_table(hits: list[dict[str, Any]]) -> pd.DataFrame:
    """Flatten GDC /files hits into file_id, file_name, sample barcode, patient barcode."""
    rows = []
    for hit in hits:
        for case in hit.get("cases", []):
            for sample in case.get("samples", []):
                rows.append({
                    "file_id": hit["file_id"],
                    "file_name": hit["file_name"],
                    "sample_barcode": sample["submitter_id"],
                    "patient_barcode": case["submitter_id"],
                    "sample_type": sample.get("sample_type"),
                })
    return pd.DataFrame(rows)


def one_sample_per_patient(files: pd.DataFrame) -> pd.DataFrame:
    """Keep one primary-tumour file per patient: lowest sample/vial code (01A before 01B), then file_id."""
    files = files.sort_values(["patient_barcode", "sample_barcode", "file_id"])
    return files.drop_duplicates("patient_barcode", keep="first").reset_index(drop=True)


def parse_star_counts(text_or_path, value_column: str = "tpm_unstranded", protein_coding_only: bool = True) -> pd.Series:
    """Parse a GDC STAR augmented gene counts file into a gene_name-indexed series."""
    df = pd.read_csv(text_or_path, sep="\t", comment="#")
    df = df[~df["gene_id"].astype(str).str.startswith("N_")]
    if protein_coding_only:
        df = df[df["gene_type"] == "protein_coding"]
    df = df.dropna(subset=["gene_name"])
    return df.groupby("gene_name")[value_column].mean()


def log2p1(matrix: pd.DataFrame) -> pd.DataFrame:
    return np.log2(matrix.clip(lower=0) + 1.0)


def _stage_group(value: Any) -> Any:
    if not isinstance(value, str):
        return pd.NA
    text = value.strip()
    match = re.search(r"\b(?:stage\s*)?(IV|III|II|I)[ABC]?\d?\b", text, flags=re.IGNORECASE)
    if match:
        return match.group(1).upper()
    arabic = re.search(r"(?:^|stage\s*|\b)([1-4])[ABC]?\d?\b", text, flags=re.IGNORECASE)
    return {"1": "I", "2": "II", "3": "III", "4": "IV"}[arabic.group(1)] if arabic else pd.NA


def cdr_clinical_table(cdr: pd.DataFrame, cancer: str, endpoint: str = "PFI") -> pd.DataFrame:
    """Build the pipeline clinical table from the TCGA-CDR supplementary sheet.

    TCGA-CDR columns used: bcr_patient_barcode, type, <endpoint>, <endpoint>.time,
    age_at_initial_pathologic_diagnosis, ajcc_pathologic_tumor_stage.
    """
    endpoint = endpoint.upper()
    needed = ["bcr_patient_barcode", "type", endpoint, f"{endpoint}.time"]
    missing = [c for c in needed if c not in cdr.columns]
    if missing:
        raise ValueError(f"TCGA-CDR table missing columns: {missing}")
    sub = cdr[cdr["type"] == cancer].copy()
    out = pd.DataFrame({
        "sample_id": sub["bcr_patient_barcode"].astype(str),
        "time": pd.to_numeric(sub[f"{endpoint}.time"], errors="coerce"),
        "event": pd.to_numeric(sub[endpoint], errors="coerce"),
        "endpoint": endpoint,
        "age": pd.to_numeric(sub.get("age_at_initial_pathologic_diagnosis"), errors="coerce"),
        "stage": sub.get("ajcc_pathologic_tumor_stage", pd.Series(pd.NA, index=sub.index)).map(_stage_group),
    })
    out = out.dropna(subset=["time", "event"])
    out = out[out["time"] > 0]
    out["time"] = out["time"].astype(float)
    out["event"] = out["event"].astype(int)
    return out.reset_index(drop=True)


# ---------------------------------------------------------------- GEO

SYMBOL_COLUMNS = ["Gene Symbol", "GENE_SYMBOL", "GeneSymbol", "Gene symbol", "Symbol", "SYMBOL", "gene_symbol"]


def _first_symbol(value: Any) -> Any:
    if not isinstance(value, str) or not value.strip() or value.strip() in {"---", "-"}:
        return pd.NA
    if "///" in value:  # probe maps to several genes -> ambiguous, drop
        return pd.NA
    return value.strip()


def _symbol_from_gene_assignment(value: Any) -> Any:
    # Affymetrix "gene_assignment": "NM_000000 // SYMBOL // description // ... /// NM_... // OTHER"
    if not isinstance(value, str) or "//" not in value:
        return pd.NA
    if "///" in value:
        symbols = {chunk.split("//")[1].strip() for chunk in value.split("///") if chunk.count("//") >= 1}
        return symbols.pop() if len(symbols) == 1 else pd.NA
    return value.split("//")[1].strip() or pd.NA


def gpl_gene_map(gpl_table: pd.DataFrame, symbol_column: str | None = None) -> tuple[pd.DataFrame, str]:
    """Probe -> gene symbol map from a GEO platform table. Multi-gene probes are dropped."""
    if symbol_column is None:
        symbol_column = next((c for c in SYMBOL_COLUMNS if c in gpl_table.columns), None)
        if symbol_column is None and "gene_assignment" in gpl_table.columns:
            symbol_column = "gene_assignment"
    if symbol_column is None or symbol_column not in gpl_table.columns:
        raise ValueError(f"no gene symbol column found; available: {list(gpl_table.columns)}")
    parse = _symbol_from_gene_assignment if symbol_column == "gene_assignment" else _first_symbol
    gm = pd.DataFrame({"feature_id": gpl_table["ID"].astype(str), "gene_id": gpl_table[symbol_column].map(parse)})
    return gm.dropna().reset_index(drop=True), symbol_column


def needs_log2(matrix: pd.DataFrame) -> bool:
    """Heuristic: raw (MAS5/intensity) arrays have values in the thousands; log2 data stays < ~20."""
    values = matrix.to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    return bool(values.size and np.nanpercentile(values, 99) > 100)


def characteristics_overview(pheno: pd.DataFrame, max_examples: int = 5) -> pd.DataFrame:
    """One row per phenotype column with distinct-value count and examples, to write the mapping."""
    rows = []
    for col in pheno.columns:
        values = pheno[col].dropna().astype(str)
        uniq = values.unique()
        rows.append({
            "column": col,
            "n_distinct": len(uniq),
            "examples": " | ".join(uniq[:max_examples]),
        })
    return pd.DataFrame(rows)


def _strip_key(value: Any) -> Any:
    # GEO characteristics are often "key: value"; GEOparse may already split them.
    if isinstance(value, str) and ": " in value:
        return value.split(": ", 1)[1].strip()
    return value


def resolve_columns(pheno: pd.DataFrame, name: str | list[str]) -> list[str]:
    """All columns holding a characteristic.

    `name` is an exact column name or a characteristic key (or a list of them). GEOparse
    names characteristics 'characteristics_ch1.<i>.<key>', and when samples list their
    characteristics in different orders the same key is split over several <i> columns
    (e.g. '.10.relapse', '.11.relapse', '.12.relapse'). Matching is case-insensitive.
    """
    names = [name] if isinstance(name, str) else list(name)
    hits: list[str] = []
    for n in names:
        if n in pheno.columns:
            hits.append(n)
            continue
        hits.extend(c for c in pheno.columns if c.lower().endswith("." + n.lower()) and c not in hits)
    if not hits:
        raise ValueError(f"no column found for {names}")
    return hits


def resolve_column(pheno: pd.DataFrame, name: str) -> str:
    """Single-column variant kept for callers that need exactly one column."""
    hits = resolve_columns(pheno, name)
    if len(hits) != 1:
        raise ValueError(f"column '{name}' not found uniquely; candidates: {hits}")
    return hits[0]


def coalesce_characteristic(pheno: pd.DataFrame, name: str | list[str]) -> pd.Series:
    """One value per sample: the first non-missing value across all matching columns.

    Raises if a sample has conflicting non-missing values for the same characteristic.
    """
    cols = resolve_columns(pheno, name)
    values = pheno[cols].apply(lambda col: col.map(_strip_key))
    for idx, row in values.iterrows():
        distinct = {str(v).strip() for v in row.dropna()}
        if len(distinct) > 1:
            raise ValueError(f"sample {idx} has conflicting values for {name}: {sorted(distinct)}")
    return values.bfill(axis=1).iloc[:, 0]


def breast_stage_from_tnm(t: Any, n: Any, m: Any) -> Any:
    """AJCC 7th anatomic stage group (I-IV) for breast cancer from T, N, M categories.

    T1 (incl. T1mi/a/b/c) N0 -> I; T0-1 N1 or T2 N0-1 or T3 N0 -> II; T0-2 N2, T3 N1-2,
    T4 any N, any T N3 -> III; M1 -> IV. Unknown or Tis/TX/NX values give NA.
    """
    def _num(value: Any, prefix: str) -> Any:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return None
        hit = re.match(rf"^\s*{prefix}?\s*([0-4])", str(value), flags=re.IGNORECASE)
        return int(hit.group(1)) if hit else None

    tt, nn, mm = _num(t, "T"), _num(n, "N"), _num(m, "M")
    if mm == 1:
        return "IV"
    if tt is None or nn is None or mm not in (0, None):
        return pd.NA
    if tt == 4 or nn == 3 or nn == 2 or (tt == 3 and nn >= 1):
        return "III"
    if (tt <= 1 and nn == 1) or (tt == 2 and nn <= 1) or (tt == 3 and nn == 0):
        return "II"
    if tt <= 1 and nn == 0:
        return "I" if tt == 1 else pd.NA  # T0 N0 is not a tumour stage
    return pd.NA


def t_from_size_mm(size: Any) -> Any:
    """AJCC T category from invasive tumour size in mm (T4 needs clinical extension, not size)."""
    value = pd.to_numeric(size, errors="coerce")
    if pd.isna(value) or value <= 0:
        return None
    return "1" if value <= 20 else ("2" if value <= 50 else "3")


def _stage_from_tnm_spec(df: pd.DataFrame, spec: dict[str, Any]) -> pd.Series:
    """Breast stage group from a mapping spec.

    spec: t (T column) or t_size_mm (size column in mm), n (N column), optional n_map
    (raw N label -> N category, e.g. lymph node groups), optional m (M column; when
    absent M0 is assumed, so the cohort must not include metastatic disease at diagnosis).
    """
    idx = df.index
    if "t" in spec:
        t = coalesce_characteristic(df, spec["t"])
    else:
        t = coalesce_characteristic(df, spec["t_size_mm"]).map(t_from_size_mm)
    n = coalesce_characteristic(df, spec["n"])
    if spec.get("n_map"):
        n_map = {str(k).lower(): str(v) for k, v in spec["n_map"].items()}
        n = n.map(lambda v: n_map.get(str(v).strip().lower()) if pd.notna(v) else None)
    m = coalesce_characteristic(df, spec["m"]) if spec.get("m") else pd.Series("0", index=idx)
    return pd.Series([breast_stage_from_tnm(a, b, c) for a, b, c in zip(t, n, m)], index=idx)


def apply_clinical_mapping(pheno: pd.DataFrame, mapping: dict[str, Any]) -> pd.DataFrame:
    """Build the curated clinical table from GEO phenotype data and an explicit mapping.

    mapping keys: endpoint, time_column, time_unit (days|months|years), event_column,
    event_values (list of values meaning event=1), censor_values (values meaning 0),
    optional age_column, stage_column, sample_filter {characteristic: [allowed values]},
    event_time_column (time used for patients with an event, e.g. time to metastasis,
    while time_column holds follow-up for the censored ones), and stage_from_tnm
    {t, n, m} (breast AJCC stage group derived from T/N/M when no stage is published).
    Column entries are characteristic keys (or lists of keys) and are coalesced across
    GEOparse's split columns. Rows whose event value is in neither list are dropped.
    """
    df = pheno.copy()
    for col, allowed in (mapping.get("sample_filter") or {}).items():
        values = coalesce_characteristic(df, col).astype(str).str.strip().str.lower()
        df = df[values.isin([str(a).lower() for a in allowed])]
    factor = {"days": 1.0, "months": 30.4375, "years": 365.25}[mapping.get("time_unit", "days")]
    time = pd.to_numeric(coalesce_characteristic(df, mapping["time_column"]), errors="coerce") * factor
    raw_event = coalesce_characteristic(df, mapping["event_column"]).astype(str).str.strip().str.lower()
    pos = {str(v).lower() for v in mapping["event_values"]}
    neg = {str(v).lower() for v in mapping["censor_values"]}
    event = raw_event.map(lambda v: 1 if v in pos else (0 if v in neg else np.nan))
    if mapping.get("event_time_column"):
        event_time = pd.to_numeric(coalesce_characteristic(df, mapping["event_time_column"]), errors="coerce") * factor
        time = time.where(event != 1, event_time)  # events without an event time -> NaN -> dropped
    age = pd.to_numeric(coalesce_characteristic(df, mapping["age_column"]), errors="coerce") if mapping.get("age_column") else pd.Series(np.nan, index=df.index)
    if mapping.get("stage_from_tnm"):
        stage = _stage_from_tnm_spec(df, mapping["stage_from_tnm"])
    elif mapping.get("stage_column"):
        stage = coalesce_characteristic(df, mapping["stage_column"]).map(_stage_group)
    else:
        stage = pd.Series(pd.NA, index=df.index)
    out = pd.DataFrame({
        "sample_id": df.index.astype(str),
        "time": time.values,
        "event": event.values,
        "endpoint": mapping["endpoint"],
        "age": age.values,
        "stage": stage.values,
    })
    n_before = len(out)
    out = out.dropna(subset=["time", "event"])
    out = out[out["time"] > 0]
    out["event"] = out["event"].astype(int)
    out.attrs["n_dropped"] = n_before - len(out)
    return out.reset_index(drop=True)


# ---------------------------------------------------------------- GEO RNA-seq (supplementary matrices)

def parse_series_matrix_pheno(lines) -> pd.DataFrame:
    """Phenotype table from one GEO series-matrix file (only the '!Sample_' header lines).

    Columns follow GEOparse naming so the same clinical mappings work: 'title',
    'geo_accession', ..., and 'characteristics_ch1.<i>.<key>' where <i> is the index of
    the characteristics line. Samples that list different keys on the same line end up
    in different columns, which coalesce_characteristic merges.
    """
    import csv

    fields: dict[str, list[str]] = {}
    char_lines: list[list[str]] = []
    for raw in lines:
        line = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else raw
        if line.startswith("!series_matrix_table_begin"):
            break
        if not line.startswith("!Sample_"):
            continue
        row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
        key, values = row[0][len("!Sample_"):], row[1:]
        if key == "characteristics_ch1":
            char_lines.append(values)
        elif key not in fields:
            fields[key] = values
    if "geo_accession" not in fields:
        raise ValueError("no !Sample_geo_accession line found")
    index = fields["geo_accession"]
    out = pd.DataFrame({k: v for k, v in fields.items() if len(v) == len(index)}, index=index)
    for i, values in enumerate(char_lines):
        for gsm, text in zip(index, values):
            if not text:
                continue
            key, _, value = text.partition(": ")
            out.loc[gsm, f"characteristics_ch1.{i}.{key.strip()}"] = value.strip() if _ else text.strip()
    out.index.name = "sample_id"
    return out


def map_columns_to_gsm(expr: pd.DataFrame, pheno: pd.DataFrame, drop_regex: str | None = None) -> tuple[pd.DataFrame, dict[str, int]]:
    """Rename supplementary-matrix columns (sample titles) to GSM accessions.

    Columns matching drop_regex (e.g. technical replicates) are removed first; columns
    without a matching title are dropped and counted. Raises if nothing matches.
    """
    cols = list(expr.columns)
    dropped_regex = [c for c in cols if drop_regex and re.search(drop_regex, str(c), flags=re.IGNORECASE)]
    cols = [c for c in cols if c not in dropped_regex]
    title_to_gsm = {str(t).strip(): gsm for gsm, t in pheno["title"].items()}
    if len(title_to_gsm) != len(pheno):
        raise ValueError("sample titles are not unique; cannot map columns to GSM")
    mapped = {c: title_to_gsm[str(c).strip()] for c in cols if str(c).strip() in title_to_gsm}
    if not mapped:
        raise ValueError(f"no expression column matches a sample title; columns look like {cols[:5]}, titles like {list(title_to_gsm)[:5]}")
    out = expr[list(mapped)].rename(columns=mapped)
    stats = {"n_columns": len(expr.columns), "n_dropped_regex": len(dropped_regex), "n_unmatched": len(cols) - len(mapped), "n_kept": len(mapped)}
    return out, stats


def links_in_listing(html: str, pattern: str) -> list[str]:
    """File names in an NCBI FTP-over-HTTPS directory listing that match a regex."""
    names = re.findall(r'href="([^"/?]+)"', html)
    return sorted({n for n in names if re.search(pattern, n)})
