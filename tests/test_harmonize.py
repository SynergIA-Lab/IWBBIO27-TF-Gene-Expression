import numpy as np
import pandas as pd

from regbio.harmonize import aggregate_duplicates, cohort_zscore, lock_common_genes, sample_rank


def test_duplicate_genes_are_averaged():
    expr = pd.DataFrame({"S1": [1.0, 3.0], "S2": [2.0, 4.0]}, index=["P1", "P2"])
    gene_map = pd.DataFrame({"feature_id": ["P1", "P2"], "gene_id": ["GENE1", "GENE1"]})
    result = aggregate_duplicates(expr, gene_map)
    assert result.loc["GENE1", "S1"] == 2.0


def test_common_genes_identical_across_cohorts():
    a = pd.DataFrame(index=["A", "B", "C"])
    b = pd.DataFrame(index=["B", "C", "D"])
    locked = lock_common_genes({"tcga": a, "geo": b})
    assert list(locked["tcga"].index) == ["B", "C"]
    assert list(locked["geo"].index) == ["B", "C"]


def test_cohort_zscore_centers_each_cohort_independently():
    tcga = pd.DataFrame({"T1": [1.0, 10.0], "T2": [3.0, 20.0]}, index=["G1", "G2"])
    geo = pd.DataFrame({"G1": [100.0, 5.0], "G2": [200.0, 7.0]}, index=["G1", "G2"])
    z = cohort_zscore({"tcga": tcga, "geo": geo})
    assert np.allclose(z["tcga"].mean(axis=1).values, 0.0)
    assert np.allclose(z["geo"].mean(axis=1).values, 0.0)


def test_sample_rank_scales_columns_to_percentiles():
    expr = pd.DataFrame({"S1": [1.0, 2.0, 3.0]}, index=["A", "B", "C"])
    ranked = sample_rank({"x": expr})["x"]
    assert ranked.loc["C", "S1"] == 1.0
