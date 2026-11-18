"""Transfer gap: internal TCGA C-index (nested CV) versus external C-index (frozen model on GEO).

Inputs
  --internal  *_internal_cv.tsv files from scripts/internal_cv.py (one per analysis x representation)
  --metrics   *_metrics.tsv files from scripts/compare_representations.py (one per validation cohort)
Outputs
  --output    one row per analysis x validation cohort x representation:
              internal mean C-index, external C-index (95% CI), drop = external - internal
  --figure-dir  optional; writes fig_performance_drop.pdf/.png (dumbbell: hollow = TCGA, filled = GEO)

Descriptive only: the internal and external C-indices come from different patients, so their
difference has no paired test. Confirmatory claims stay with the pre-specified paired comparisons.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from regbio.io import write_table


def build_table(internal: pd.DataFrame, metrics: pd.DataFrame) -> pd.DataFrame:
    ext = metrics.rename(columns={"model": "representation", "primary_ipcw_cindex": "external_cindex",
                                  "primary_ci_lower": "external_ci_lower", "primary_ci_upper": "external_ci_upper",
                                  "n_patients": "external_n_patients", "n_events": "external_n_events"})
    ext = ext[ext["representation"].isin(internal["representation"].unique())]
    cols = ["analysis", "validation_cohort", "representation", "external_n_patients", "external_n_events",
            "external_cindex", "external_ci_lower", "external_ci_upper"]
    out = internal[["analysis", "representation", "n_patients", "n_events", "n_folds", "internal_cindex_mean",
                    "internal_cindex_sd"]].rename(columns={"n_patients": "tcga_n_patients", "n_events": "tcga_n_events"})
    out = out.merge(ext[cols], on=["analysis", "representation"], how="inner")
    out["drop"] = out["external_cindex"] - out["internal_cindex_mean"]
    return out


def fig_performance_drop(table: pd.DataFrame, reps: list[str], outdir: Path) -> None:
    import matplotlib.pyplot as plt

    from make_figures import INK2, REP_STYLE, SURFACE, WIDTH_IN, _grid_positions, _save, _xgrid

    groups = dict(list(table.groupby(["analysis", "validation_cohort"], sort=False)))
    nrows, ncols, pos = _grid_positions([a for a, _ in groups])
    fig, axes = plt.subplots(nrows, ncols, figsize=(WIDTH_IN, (0.17 * len(reps) + 0.6) * nrows + 0.55),
                             sharey=True, sharex=True, squeeze=False)
    used = set()
    for (analysis, cohort), g in groups.items():
        r_, c_ = pos[analysis]
        used.add((r_, c_))
        ax = axes[r_][c_]
        g = g.set_index("representation")
        for i, rep in enumerate(reps):
            if rep not in g.index:
                continue
            _, color, marker = REP_STYLE[rep]
            y = len(reps) - i
            a, b = g.loc[rep, "internal_cindex_mean"], g.loc[rep, "external_cindex"]
            ax.plot([a, b], [y, y], color=color, lw=1.0, solid_capstyle="round")
            ax.plot(a, y, marker=marker, ms=4.2, mfc=SURFACE, mec=color, mew=0.9, ls="none")
            ax.plot(b, y, marker=marker, ms=4.2, color=color, mec=SURFACE, mew=0.6, ls="none")
        ax.axvline(0.5, color=INK2, lw=0.7, ls=(0, (3, 2)))
        ax.set_yticks(range(len(reps), 0, -1))
        ax.set_yticklabels([REP_STYLE[r][0] for r in reps])
        ax.set_ylim(0.4, len(reps) + 0.6)
        ax.set_title(f"{analysis.replace('_', ' ')} → {cohort}\nn = {int(g['tcga_n_patients'].max())} → {int(g['external_n_patients'].max())}",
                     loc="left", fontsize=6.2, linespacing=1.3)
        _xgrid(ax, 0.1)
    # common x-range with a margin, so no marker sits on the panel edge
    lo = min(table["internal_cindex_mean"].min(), table["external_cindex"].min(), 0.5)
    hi = max(table["internal_cindex_mean"].max(), table["external_cindex"].max())
    axes[0][0].set_xlim(lo - 0.02, hi + 0.02)
    for r_ in range(nrows):
        for c_ in range(ncols):
            if (r_, c_) not in used:
                axes[r_][c_].axis("off")
    fig.supxlabel("IPCW C-index: hollow = TCGA internal (nested CV) → filled = external GEO; dashed = 0.5", fontsize=7)
    fig.tight_layout(h_pad=0.8, w_pad=0.6)
    _save(fig, outdir, "fig_performance_drop")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--internal", nargs="+", required=True)
    parser.add_argument("--metrics", nargs="+", required=True)
    parser.add_argument("--representations", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--figure-dir")
    args = parser.parse_args()

    internal = pd.concat([pd.read_csv(p, sep="\t") for p in args.internal], ignore_index=True)
    metrics = pd.concat([pd.read_csv(p, sep="\t") for p in args.metrics], ignore_index=True)
    table = build_table(internal, metrics)
    write_table(table, args.output)
    if args.figure_dir:
        fig_performance_drop(table, args.representations, Path(args.figure_dir))
    summary = table.groupby("representation")["drop"].agg(["mean", "min", "max"]).reindex(args.representations)
    print("External minus internal C-index (mean, min, max over validation cohorts):")
    print(summary.round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
