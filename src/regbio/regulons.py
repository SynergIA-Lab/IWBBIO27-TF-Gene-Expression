from __future__ import annotations

from pathlib import Path

import pandas as pd

from regbio.io import read_regulon


def filter_regulons(regulon: pd.DataFrame, common_genes: set[str], min_targets: int) -> pd.DataFrame:
    filtered = regulon[regulon["target"].astype(str).isin(common_genes)].copy()
    counts = filtered.groupby("source")["target"].nunique()
    keep = counts[counts >= min_targets].index
    return filtered[filtered["source"].isin(keep)].reset_index(drop=True)


def tf_control_set(filtered_collectri: pd.DataFrame) -> list[str]:
    return sorted(filtered_collectri["source"].astype(str).unique())


def target_counts(filtered_regulon: pd.DataFrame) -> pd.DataFrame:
    if filtered_regulon.empty:
        return pd.DataFrame(columns=["source", "n_targets"])
    return (
        filtered_regulon.groupby("source")["target"]
        .nunique()
        .reset_index(name="n_targets")
        .sort_values("source")
    )


def load_regulon_file(path: str | Path) -> pd.DataFrame:
    return read_regulon(path)


def _core_regulon_columns(regulon: pd.DataFrame) -> pd.DataFrame:
    if "weight" not in regulon.columns:
        regulon = regulon.assign(weight=1.0)
    return regulon[["source", "target", "weight"]].copy()


def load_collectri(organism: str = "human", license: str = "academic") -> pd.DataFrame:
    import decoupler as dc

    return _core_regulon_columns(dc.op.collectri(organism=organism, license=license, verbose=False))


def load_dorothea(organism: str = "human", levels: list[str] | None = None, license: str = "academic") -> pd.DataFrame:
    import decoupler as dc

    if levels is None:
        levels = ["A", "B", "C"]
    return _core_regulon_columns(dc.op.dorothea(organism=organism, levels=levels, license=license, verbose=False))
