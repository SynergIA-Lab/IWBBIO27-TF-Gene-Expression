"""Paper figures from the real-data outputs (vector PDF + 300 dpi PNG).

fig_design                study design, filled in from the run (cohorts, sizes, events, settings)
fig_confirmatory_forest   paired C-index differences per cohort + pooled fixed effect (main result)
fig_external_cindex       external IPCW C-index (95% CI) of every representation and the clinical model
fig_stability             Kuncheva stability of the top-k panel, all analyses in one panel
fig_added_value           C-index gain of clinical+representation over the clinical model
fig_km                    Kaplan-Meier curves of TCGA-median CollecTRI risk groups per validation cohort
fig_km_gene               the same for the gene-expression panel (supplementary)

Multi-cohort panels are laid out as a grid: one column per cancer, one row per
training endpoint (PFI, OS). Sized for the LNBI text width (12.2 cm). Identity is
never colour-alone: every representation has its own marker and is labelled on the
axis or in a legend, so figures survive greyscale printing.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import MaxNLocator, MultipleLocator

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# Reference categorical palette (validated: CVD-safe adjacent pairs; contrast relief via markers + axis labels)
INK, INK2, GRID, SURFACE, MUTED = "#0b0b0b", "#52514e", "#e2e1dc", "#ffffff", "#f4f3ef"
REP_STYLE = {
    "gene_expression": ("Gene expression", "#2a78d6", "o"),
    "tf_expression": ("TF expression", "#eb6834", "s"),
    "collectri_ulm": ("CollecTRI + ULM", "#1baf7a", "D"),
    "dorothea_viper": ("DoRothEA + VIPER", "#eda100", "^"),
    "clinical": ("Clinical (age + stage)", INK2, "X"),
}
REP_FEATURE = {"gene_expression": "genes", "tf_expression": "TFs", "collectri_ulm": "TF activities", "dorothea_viper": "TF activities"}
WIDTH_IN = 12.2 / 2.54

plt.rcParams.update({
    "font.size": 7.5, "axes.titlesize": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 6.5, "ytick.labelsize": 7,
    "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0, "pdf.fonttype": 42, "font.family": "DejaVu Sans",
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def _save(fig, outdir: Path, name: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _xgrid(ax, step: float | None = None) -> None:
    ax.xaxis.set_major_locator(MultipleLocator(step) if step else MaxNLocator(nbins=5))
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def _cancer(analysis: str) -> str:
    return analysis.split("_", 1)[0]


def _kind(analysis: str) -> str:
    return analysis.split("_", 1)[1] if "_" in analysis else analysis


def _grid_positions(analyses: list[str]) -> tuple[int, int, dict[str, tuple[int, int]]]:
    """Columns = cancers, rows = training endpoints, when that tiles exactly; else one row."""
    cancers = list(dict.fromkeys(_cancer(a) for a in analyses))
    kinds = list(dict.fromkeys(_kind(a) for a in analyses))
    if 1 < len(cancers) and len(cancers) * len(kinds) == len(analyses):
        return len(kinds), len(cancers), {a: (kinds.index(_kind(a)), cancers.index(_cancer(a))) for a in analyses}
    ncols = min(3, len(analyses))
    nrows = math.ceil(len(analyses) / ncols)
    return nrows, ncols, {a: divmod(i, ncols) for i, a in enumerate(analyses)}


def _cohort_label(row) -> str:
    return f"{row['validation_cohort']}\n({row['analysis'].replace('_', ' ')})"


# ------------------------------------------------------------------ study design

def _box(ax, x, y, w, h, text, edge=INK2, face=SURFACE, weight="normal", size=6.3, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012", fc=face, ec=edge, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, fontweight=weight, linespacing=1.35)


def _arrow(ax, y_from, y_to, x=0.3, label=None):
    ax.annotate("", xy=(x, y_to), xytext=(x, y_from), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.8, mutation_scale=7))
    if label:
        ax.text(x + 0.015, (y_from + y_to) / 2, label, va="center", ha="left", fontsize=5.8, color=INK2, style="italic")


def fig_design(cfg: dict, qc: pd.DataFrame | None, metrics: pd.DataFrame, stab: pd.DataFrame, reps: list[str], outdir: Path) -> None:
    analyses = list(dict.fromkeys(metrics["analysis"]))
    cancers = list(dict.fromkeys(_cancer(a) for a in analyses))
    real = cfg.get("real", {})
    k = real.get("panel_size", cfg.get("feature_panel_size", 25))
    tau_y = cfg.get("cindex_tau_days", 1825) / 365.25
    n_res = cfg.get("n_stability_resamples", 50)
    n_val = metrics.drop_duplicates(["analysis", "validation_cohort"]).shape[0]

    fig, ax = plt.subplots(figsize=(WIDTH_IN, 4.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.04, 1)
    ax.axis("off")
    gap = 0.015
    colw = (1 - gap * (len(cancers) - 1)) / len(cancers)

    def stage(y, text):
        ax.text(0, y, text, fontsize=6.8, fontweight="bold", color=INK, va="bottom")

    # 1. training cohorts
    stage(0.955, "1  Training (TCGA, RNA-seq)")
    for i, c in enumerate(cancers):
        lines = [f"TCGA-{c}"]
        if qc is not None:
            sub = qc[qc["cohort"] == f"TCGA_{c}"]
            if not sub.empty:
                ns = sub.set_index("endpoint")["n_patients"].astype(int)
                if ns.nunique() == 1:
                    lines.append(f"n = {ns.iloc[0]}")
                else:
                    lines.append("n = " + " / ".join(f"{n} {e}" for e, n in ns.items()))
                lines.append(" · ".join(f"{r.endpoint} {int(r.n_events)} ev" for r in sub.itertuples()))
        _box(ax, i * (colw + gap), 0.83, colw, 0.115, "\n".join(lines))
    _arrow(ax, 0.828, 0.79, label="same patients, cohort-wise z-score, common gene set")

    # 2. representations
    stage(0.765, "2  Four representations of each tumour")
    rw = (1 - gap * (len(reps) - 1)) / len(reps)
    for i, rep in enumerate(reps):
        name, color, _ = REP_STYLE[rep]
        sizes = stab[stab["representation"] == rep]["universe_size"]
        size = f"{int(sizes.min()):,}–{int(sizes.max()):,}" if len(sizes) and sizes.min() != sizes.max() else (f"{int(sizes.iloc[0]):,}" if len(sizes) else "")
        _box(ax, i * (rw + gap), 0.655, rw, 0.1, f"{name}\n{size} {REP_FEATURE.get(rep, '')}", edge=color, lw=1.4, size=5.3)
    _arrow(ax, 0.653, 0.615)

    # 3. model
    stage(0.59, "3  Biomarker panel and model, fitted on TCGA only and frozen")
    _box(ax, 0, 0.475, 1, 0.105,
         f"Top-{k} features by univariate Cox score test  →  ridge-like Coxnet\n"
         f"(l1 ratio {cfg.get('primary_l1_ratio', 0.05)}, α by {cfg.get('cv_splits', 5)}-fold CV)  ·  sensitivity: LASSO, random survival forest\n"
         f"panel stability: {n_res} bootstrap resamples, Kuncheva index", face=MUTED, edge=GRID, size=6)
    _arrow(ax, 0.473, 0.435, label="risk scores computed in each external cohort, no refitting")

    # 4. validation cohorts, one column per cancer
    stage(0.41, f"4  External validation ({n_val} GEO cohorts, other platforms and centres)")
    plat = {} if qc is None else dict(zip(qc["cohort"], qc["platform_class"]))
    endp = {} if qc is None else dict(zip(qc["cohort"], qc["endpoint"]))
    vals = metrics.drop_duplicates(["analysis", "validation_cohort"])
    for i, c in enumerate(cancers):
        sub = vals[vals["analysis"].map(_cancer) == c].reset_index(drop=True)
        h = 0.3 / max(1, len(sub)) - 0.012
        for j, r in sub.iterrows():
            y = 0.395 - (j + 1) * (h + 0.012)
            p = {"microarray": "array", "rnaseq": "RNA-seq"}.get(plat.get(r["validation_cohort"], ""), "")
            _box(ax, i * (colw + gap), y, colw, h,
                 f"{r['validation_cohort']}  ({p})\n{endp.get(r['validation_cohort'], '')} · "
                 f"n = {int(r['n_patients'])}, {int(r['n_events'])} ev\nmodel: TCGA {_kind(r['analysis'])}")
    _arrow(ax, 0.083, 0.06)

    # 5. inference
    _box(ax, 0, -0.03, 1, 0.085,
         f"5  IPCW C-index (τ = {tau_y:.0f} y)  ·  paired bootstrap ΔC per cohort\n"
         f"fixed-effect pooling across cohorts  ·  Holm over {len(cfg.get('primary_comparisons', []))} pre-specified comparisons",
         face=MUTED, edge=GRID, weight="bold", size=5.9)
    _save(fig, outdir, "fig_design")


# ------------------------------------------------------------------ results

def fig_confirmatory_forest(comps: pd.DataFrame, meta: pd.DataFrame, primary: list[list[str]], outdir: Path) -> None:
    labels = [f"{a} - {b}" for a, b in primary]
    n = comps["validation_cohort"].nunique()
    fig, axes = plt.subplots(1, len(labels), figsize=(WIDTH_IN, 0.26 * (n + 2) + 0.8), sharey=True, sharex=True)
    axes = np.atleast_1d(axes)
    for ax, label in zip(axes, labels):
        sub = comps[(comps["comparison"] == label) & (comps["comparison_status"] == "confirmatory")].reset_index(drop=True)
        ys = np.arange(len(sub))[::-1] + 1.4
        d_, lo_, hi_ = (sub[c].to_numpy(dtype=float) for c in ("difference", "ci_lower", "ci_upper"))
        ax.errorbar(d_, ys, xerr=[d_ - lo_, hi_ - d_], fmt="o", color=INK, ms=3.5, elinewidth=0.9, capsize=0)
        m = meta[(meta["comparison"] == label) & (meta["role"] == "primary")]
        summary = ""
        if not m.empty:
            m = m.iloc[0]
            d = [m["fe_ci_lower"], m["fe_estimate"], m["fe_ci_upper"], m["fe_estimate"]]
            ax.fill(d, [0.4, 0.62, 0.4, 0.18], color=REP_STYLE["collectri_ulm"][1], lw=0)
            i2 = f", I² = {100 * float(m['i2']):.0f}%" if "i2" in m and pd.notna(m["i2"]) else ""
            summary = (f"\npooled Δ = {m['fe_estimate']:+.3f} [{m['fe_ci_lower']:+.3f}, {m['fe_ci_upper']:+.3f}]"
                       f"\nHolm p = {float(m['fe_p_holm']):.3f}{i2}")
        ax.axvline(0, color=INK2, lw=0.7, ls=(0, (3, 2)))
        ax.axhline(0.95, color=GRID, lw=0.6)
        ax.set_yticks(list(ys) + [0.4])
        ax.set_yticklabels([_cohort_label(r) for _, r in sub.iterrows()] + ["Pooled (FE)"], fontsize=6.3)
        ax.set_ylim(-0.1, ys.max() + 0.6)
        _, right = label.split(" - ")
        ax.set_title(f"CollecTRI vs {REP_STYLE[right][0]}{summary}", loc="left", fontsize=7, linespacing=1.4)
        _xgrid(ax, 0.05)
    fig.supxlabel("Δ C-index, CollecTRI + ULM minus comparator (95% CI); > 0 favours CollecTRI", fontsize=7)
    fig.tight_layout()
    _save(fig, outdir, "fig_confirmatory_forest")


def _grid_dot_panels(df: pd.DataFrame, models: list[str], value: str, lo: str, hi: str, xlabel: str, ref: float,
                     outdir: Path, name: str, title_fn, tick_step: float | None) -> None:
    groups = dict(list(df.groupby(["analysis", "validation_cohort"], sort=False)))
    analyses = [a for a, _ in groups]
    nrows, ncols, pos = _grid_positions(analyses)
    fig, axes = plt.subplots(nrows, ncols, figsize=(WIDTH_IN, (0.135 * len(models) + 0.5) * nrows + 0.15),
                             sharey=True, sharex=True, squeeze=False)
    used = set()
    for (analysis, cohort), g in groups.items():
        r_, c_ = pos[analysis]
        used.add((r_, c_))
        ax = axes[r_][c_]
        g = g.set_index("model")
        for i, mdl in enumerate(models):
            if mdl not in g.index:
                continue
            rep = mdl.replace("clinical+", "")
            _, color, marker = REP_STYLE.get(rep, (mdl, INK2, "o"))
            r = g.loc[mdl]
            y = len(models) - i
            if lo in g.columns and pd.notna(r.get(lo)):
                ax.plot([r[lo], r[hi]], [y, y], color=color, lw=1.0, solid_capstyle="round")
            ax.plot(r[value], y, marker=marker, color=color, ms=4.2, mec=SURFACE, mew=0.6, ls="none")
        ax.axvline(ref, color=INK2, lw=0.7, ls=(0, (3, 2)))
        ax.set_yticks(range(len(models), 0, -1))
        ax.set_yticklabels([REP_STYLE.get(m.replace("clinical+", ""), (m,))[0] for m in models])
        ax.set_ylim(0.4, len(models) + 0.6)
        ax.set_title(title_fn(analysis, cohort, g), loc="left", fontsize=6.2, linespacing=1.3)
        _xgrid(ax, tick_step)
    for r_ in range(nrows):
        for c_ in range(ncols):
            if (r_, c_) not in used:
                axes[r_][c_].axis("off")
    fig.supxlabel(xlabel, fontsize=7)
    fig.tight_layout(h_pad=0.4, w_pad=0.6)
    _save(fig, outdir, name)


def _panel_title(a, c, g):
    return f"{c} · {a.replace('_', ' ')}\nn = {int(g['n_patients'].max())}, {int(g['n_events'].max())} events"


def fig_external_cindex(metrics: pd.DataFrame, reps: list[str], outdir: Path) -> None:
    models = reps + ["clinical"]
    df = metrics[metrics["model"].isin(models)]
    _grid_dot_panels(df, models, "primary_ipcw_cindex", "primary_ci_lower", "primary_ci_upper",
                     "External IPCW C-index (95% CI); dashed line = no discrimination", 0.5,
                     outdir, "fig_external_cindex", _panel_title, tick_step=0.1)


def fig_added_value(comps: pd.DataFrame, reps: list[str], outdir: Path) -> None:
    sub = comps[comps["comparison"].str.endswith(" - clinical")].copy()
    if sub.empty:
        return
    sub["model"] = sub["comparison"].str.replace(" - clinical", "", regex=False)
    models = [f"clinical+{r}" for r in reps]
    _grid_dot_panels(sub, models, "difference", "ci_lower", "ci_upper",
                     "Δ C-index of clinical + panel over clinical model alone (95% CI); > 0 = added value", 0.0,
                     outdir, "fig_added_value", _panel_title, tick_step=0.1)


def fig_stability(stab: pd.DataFrame, reps: list[str], outdir: Path) -> None:
    """All analyses in one panel: one row per analysis, one marker per representation."""
    analyses = list(dict.fromkeys(stab["analysis"]))
    fig, ax = plt.subplots(figsize=(WIDTH_IN, 0.2 * len(analyses) + 0.75))
    for i, a in enumerate(analyses):
        y = len(analyses) - i
        g = stab[stab["analysis"] == a].set_index("representation")
        vals = [g.loc[r, "kuncheva_mean"] for r in reps if r in g.index]
        ax.plot([min(vals), max(vals)], [y, y], color=GRID, lw=2.2, solid_capstyle="round", zorder=1)
        for j, rep in enumerate(reps):
            if rep not in g.index:
                continue
            _, color, marker = REP_STYLE[rep]
            dy = (j - (len(reps) - 1) / 2) * 0.09  # small dodge so equal values stay visible
            ax.plot(g.loc[rep, "kuncheva_mean"], y + dy, marker=marker, color=color, ms=5, mec=SURFACE, mew=0.7, ls="none", zorder=2)
    ax.set_yticks(range(len(analyses), 0, -1))
    ax.set_yticklabels([a.replace("_", " ") for a in analyses])
    ax.set_ylim(0.4, len(analyses) + 0.6)
    ax.set_xlim(0, float(stab["kuncheva_mean"].max()) * 1.12)
    _xgrid(ax, 0.1)
    handles = [Line2D([], [], marker=REP_STYLE[r][2], color=REP_STYLE[r][1], mec=SURFACE, ls="none", ms=5,
                      label=f"{REP_STYLE[r][0]}") for r in reps]
    ax.legend(handles=handles, ncol=len(reps), loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False,
              handletextpad=0.2, columnspacing=1.0)
    k = int(stab["k"].iloc[0]) if "k" in stab else None
    ax.set_xlabel(f"Mean pairwise Kuncheva index of the top-{k} panel over TCGA bootstraps (higher = more stable)", fontsize=7)
    fig.tight_layout()
    _save(fig, outdir, "fig_stability")


def fig_km(risks: pd.DataFrame, rep: str, outdir: Path, name: str) -> None:
    from sksurv.compare import compare_survival
    from sksurv.nonparametric import kaplan_meier_estimator
    from sksurv.util import Surv

    groups = dict(list(risks.groupby(["analysis", "validation_cohort"], sort=False)))
    nrows, ncols, pos = _grid_positions([a for a, _ in groups])
    label, color, _ = REP_STYLE[rep]
    fig, axes = plt.subplots(nrows, ncols, figsize=(WIDTH_IN, 1.95 * nrows + 0.45), sharey=True, squeeze=False)
    used = set()
    for (analysis, cohort), g in groups.items():
        r_, c_ = pos[analysis]
        used.add((r_, c_))
        ax = axes[r_][c_]
        years = g["time"] / 365.25
        ns = {}
        for grp, style in ((0, dict(color=color, lw=1.1, alpha=0.55, ls=(0, (4, 2)))), (1, dict(color=color, lw=1.3))):
            mask = g[f"{rep}_high"].to_numpy() == grp
            ns[grp] = int(mask.sum())
            if mask.sum() < 2:
                continue
            t, s_ = kaplan_meier_estimator(g["event"].astype(bool)[mask], years[mask])
            ax.step(t, s_, where="post", **style)
        p_txt = ""
        try:
            _, p = compare_survival(Surv.from_arrays(event=g["event"].astype(bool), time=years), g[f"{rep}_high"].to_numpy())
            p_txt = f"log-rank p = {p:.1e}" if p < 1e-3 else f"log-rank p = {p:.3f}"
        except ValueError:
            pass
        ax.set_ylim(0, 1.02)
        ax.xaxis.set_major_locator(MultipleLocator(2))
        # no in-panel legend: curves can occupy any corner; groups are explained once below
        ax.set_title(f"{cohort} · {analysis.replace('_', ' ')}\nn low/high = {ns.get(0, 0)}/{ns.get(1, 0)}\n{p_txt}",
                     loc="left", fontsize=6, linespacing=1.3)
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        if c_ == 0:
            ax.set_ylabel("Event-free probability")
    for r_ in range(nrows):
        for c_ in range(ncols):
            if (r_, c_) not in used:
                axes[r_][c_].axis("off")
    handles = [Line2D([], [], color=color, lw=1.1, alpha=0.55, ls=(0, (4, 2)), label="Low risk"),
               Line2D([], [], color=color, lw=1.3, label="High risk")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.045))
    fig.supxlabel(f"Years  ·  {label} panel, risk groups split at the TCGA median (frozen cutoff)", fontsize=7, y=0.035)
    fig.tight_layout(h_pad=0.8, w_pad=0.6, rect=(0, 0.03, 1, 1))
    _save(fig, outdir, name)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--metrics", nargs="+", required=True)
    parser.add_argument("--comparisons", nargs="+", required=True)
    parser.add_argument("--meta", required=True)
    parser.add_argument("--stability", nargs="+", required=True)
    parser.add_argument("--risks", nargs="*", default=[], help="per-patient risk tables; Kaplan-Meier figures are skipped without them")
    parser.add_argument("--qc", nargs="*", default=[], help="cohort QC tables; fill training/validation sizes in fig_design")
    parser.add_argument("--representations", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    from regbio.config import load_analysis_config

    cfg = load_analysis_config(args.config)
    read = lambda paths: pd.concat([pd.read_csv(p, sep="\t") for p in paths], ignore_index=True)
    metrics, comps, stab = read(args.metrics), read(args.comparisons), read(args.stability)
    meta = pd.read_csv(args.meta, sep="\t")
    qc = read(args.qc).drop_duplicates(["cohort", "endpoint"]) if args.qc else None
    out = Path(args.outdir)
    fig_design(cfg, qc, metrics, stab, args.representations, out)
    fig_confirmatory_forest(comps, meta, cfg["primary_comparisons"], out)
    fig_external_cindex(metrics, args.representations, out)
    fig_added_value(comps, args.representations, out)
    fig_stability(stab, args.representations, out)
    if args.risks:
        risks = read(args.risks)
        for rep, name in (("collectri_ulm", "fig_km"), ("gene_expression", "fig_km_gene")):
            if f"{rep}_high" in risks.columns:
                fig_km(risks, rep, out, name)
    print(f"figures written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
