from __future__ import annotations

import pandas as pd


def aggregate_duplicates(expression: pd.DataFrame, gene_map: pd.DataFrame) -> pd.DataFrame:
    required = {"feature_id", "gene_id"}
    missing = required - set(gene_map.columns)
    if missing:
        raise ValueError(f"gene_map missing required columns: {', '.join(sorted(missing))}")
    mapped = expression.copy()
    mapper = gene_map.set_index("feature_id")["gene_id"].astype(str)
    mapped = mapped.loc[mapped.index.intersection(mapper.index)]
    mapped.index = mapper.loc[mapped.index].values
    if mapped.empty:
        raise ValueError("no expression features overlap gene_map")
    return mapped.groupby(level=0).mean().sort_index()


def lock_common_genes(cohorts: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    if not cohorts:
        raise ValueError("at least one cohort is required")
    common = set.intersection(*(set(df.index.astype(str)) for df in cohorts.values()))
    if not common:
        raise ValueError("no common genes across cohorts")
    ordered = sorted(common)
    return {name: df.loc[ordered].copy() for name, df in cohorts.items()}


def cohort_zscore(cohorts: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    transformed: dict[str, pd.DataFrame] = {}
    for name, df in cohorts.items():
        means = df.mean(axis=1)
        stds = df.std(axis=1, ddof=1).replace(0, pd.NA)
        z = df.sub(means, axis=0).div(stds, axis=0).fillna(0.0)
        transformed[name] = z
    return transformed


def sample_rank(cohorts: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    transformed: dict[str, pd.DataFrame] = {}
    for name, df in cohorts.items():
        ranks = df.rank(axis=0, method="average", pct=True)
        transformed[name] = ranks.fillna(0.0)
    return transformed


def filter_low_expression(expression: pd.DataFrame, min_value: float, min_fraction: float) -> pd.DataFrame:
    """Keep genes with value >= min_value in at least min_fraction of samples (training cohort only)."""
    keep = (expression >= min_value).mean(axis=1) >= min_fraction
    return expression.loc[keep]


def subset_samples(expression: pd.DataFrame, samples: set[str]) -> pd.DataFrame:
    cols = [c for c in expression.columns if c in samples]
    if not cols:
        raise ValueError("no expression samples overlap the clinical tables")
    return expression[cols]
