"""Generate deterministic synthetic fixtures (run once; outputs are committed).

Two cohorts on different measurement scales (RNA-seq-like vs microarray-like) share a
regulatory structure: TF1 drives its targets and higher TF1 activity raises hazard.
This lets the fixture workflow exercise the real code paths (CV with informative folds,
IPCW C-index, paired bootstrap) instead of only fallbacks.
"""
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "tests" / "fixtures"
RNG = np.random.default_rng(20261026)
TFS = [f"TF{i}" for i in range(1, 7)]
GENES = [f"G{i}" for i in range(1, 35)]


def regulons() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, rows_d = [], []
    for i, tf in enumerate(TFS):
        targets = GENES[i * 5:(i * 5) + 6]
        for j, g in enumerate(targets):
            w = -1.0 if j == 5 else 1.0
            rows.append((tf, g, w))
            if j < 5:
                rows_d.append((tf, g, w))
        rows_d.append((tf, GENES[(i * 5 + 7) % len(GENES)], 1.0))  # DoRothEA differs slightly
    col = pd.DataFrame(rows, columns=["source", "target", "weight"])
    dor = pd.DataFrame(rows_d, columns=["source", "target", "weight"])
    return col, dor


def cohort(n: int, prefix: str, scale: float, shift: float, collectri: pd.DataFrame, followup: float):
    latent = RNG.normal(size=(n, len(TFS)))
    expr = pd.DataFrame(RNG.normal(scale=0.6, size=(n, len(TFS) + len(GENES))), columns=TFS + GENES)
    for k, tf in enumerate(TFS):
        expr[tf] += 0.8 * latent[:, k]
        for _, r in collectri[collectri["source"] == tf].iterrows():
            expr[r["target"]] += r["weight"] * 1.0 * latent[:, k]
    expr = expr * scale + shift
    hazard = np.exp(0.9 * latent[:, 0] - 0.4 * latent[:, 2])
    t_event = RNG.exponential(1500 / hazard)
    t_cens = RNG.uniform(200, followup, size=n)
    time = np.maximum(np.minimum(t_event, t_cens).round(), 1)
    event = (t_event <= t_cens).astype(int)
    samples = [f"{prefix}{i:03d}" for i in range(n)]
    expr.index = samples
    clinical = pd.DataFrame({
        "sample_id": samples,
        "time": time.astype(int),
        "event": event,
        "endpoint": "PFI" if prefix == "T" else "RFS",
        "age": RNG.integers(40, 85, size=n),
        "stage": RNG.choice(["I", "II", "III"], size=n, p=[0.4, 0.4, 0.2]),
        "subtype": RNG.choice(["A", "B"], size=n),
    })
    out = expr.T
    out.index.name = "gene"
    return out, clinical


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    col, dor = regulons()
    tcga_e, tcga_c = cohort(80, "T", scale=1.0, shift=6.0, collectri=col, followup=4000)
    geo_e, geo_c = cohort(60, "M", scale=2.5, shift=40.0, collectri=col, followup=3000)
    tcga_e.to_csv(OUT / "expression_tcga.tsv", sep="\t")
    geo_e.to_csv(OUT / "expression_geo.tsv", sep="\t")
    tcga_c.to_csv(OUT / "clinical_tcga.tsv", sep="\t", index=False)
    geo_c.to_csv(OUT / "clinical_geo.tsv", sep="\t", index=False)
    col.to_csv(OUT / "collectri_fixture.tsv", sep="\t", index=False)
    dor.to_csv(OUT / "dorothea_fixture.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": TFS + GENES, "gene_id": TFS + GENES}).to_csv(OUT / "gene_map.tsv", sep="\t", index=False)
    print("events: TCGA", int(tcga_c.event.sum()), "/", len(tcga_c), "| GEO", int(geo_c.event.sum()), "/", len(geo_c))


if __name__ == "__main__":
    main()
