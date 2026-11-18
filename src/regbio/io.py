from __future__ import annotations

from pathlib import Path

import pandas as pd

REQUIRED_CLINICAL_COLUMNS = {"sample_id", "time", "event"}
REQUIRED_REGULON_COLUMNS = {"source", "target"}


def read_expression(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    if df.empty:
        raise ValueError(f"expression matrix is empty: {path}")
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.apply(pd.to_numeric, errors="raise")


def write_expression(df: pd.DataFrame, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t")


def read_clinical(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    missing = REQUIRED_CLINICAL_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"clinical table missing required columns: {', '.join(sorted(missing))}")
    df = df.copy()
    df["sample_id"] = df["sample_id"].astype(str)
    df["time"] = pd.to_numeric(df["time"], errors="raise")
    df["event"] = pd.to_numeric(df["event"], errors="raise").astype(int)
    return df


def read_regulon(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    missing = REQUIRED_REGULON_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"regulon missing required columns: {', '.join(sorted(missing))}")
    if "weight" not in df.columns:
        df = df.assign(weight=1.0)
    return df[["source", "target", "weight"]].copy()


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
