# Real-data workflow. Target: `snakemake --cores 4 real`
R = config.get("real", {})
REPS = R.get("representations", [])
CANCERS = R.get("cancers", {})
if R.get("active_cancers"):  # restrict the run (config files are deep-merged by Snakemake)
    CANCERS = {c: CANCERS[c] for c in R["active_cancers"]}
ANALYSES = {a: dict(cancer=c, **s) for c, cs in CANCERS.items() for a, s in cs.get("analyses", {}).items()}
PANEL = R.get("panel_size", config["feature_panel_size"])


def _training(analysis):
    a = ANALYSES[analysis]
    return CANCERS[a["cancer"]]["training"]


def _clinical(cancer, cohort, endpoint):
    return CANCERS[cancer]["cohorts"][cohort]["clinical"][endpoint]


def _train_clinical(w):
    a = ANALYSES[w.analysis]
    return _clinical(a["cancer"], _training(w.analysis), a["endpoint"])


def _train_features(w):
    a = ANALYSES[w.analysis]
    return f"data/processed/{a['cancer']}/features/{_training(w.analysis)}.{w.rep}.tsv"


DEFAULT_REGULONS = {"collectri": "resources/regulons/collectri.tsv", "dorothea": "resources/regulons/dorothea_abc.tsv"}
REGULONS = {**DEFAULT_REGULONS, **R.get("regulons", {})}


# Only the default resource paths are produced by this rule. Regulon files given at any
# other path are treated as plain inputs, so --forceall can never delete them.
rule fetch_regulons:
    output:
        collectri=DEFAULT_REGULONS["collectri"],
        dorothea=DEFAULT_REGULONS["dorothea"]
    shell:
        "python scripts/fetch_regulons.py --collectri {output.collectri} --dorothea {output.dorothea}"


for _cancer, _spec in CANCERS.items():
    _cohorts = _spec["cohorts"]
    rule:
        name: f"harmonize_{_cancer}"
        input:
            [f for c in _cohorts.values() for f in [c["expression"], c["gene_map"], *c["clinical"].values()]]
        output:
            expand(f"data/processed/{_cancer}/harmonized/{{cohort}}.expression.tsv", cohort=list(_cohorts)),
            log=f"data/processed/{_cancer}/harmonized/harmonization_log.json"
        params:
            cohorts=" ".join(
                f"--cohort {n} {c['expression']} {c['gene_map']} {','.join(c['clinical'].values())}" for n, c in _cohorts.items()
            ),
            training=_spec["training"],
            outdir=f"data/processed/{_cancer}/harmonized"
        shell:
            "python scripts/harmonize_cancer.py {params.cohorts} --training {params.training} "
            "--min-expression {R[min_expression]} --min-fraction {R[min_fraction]} "
            "--method {config[harmonization]} --outdir {params.outdir}"


rule real_features:
    input:
        expression="data/processed/{cancer}/harmonized/{cohort}.expression.tsv",
        collectri=REGULONS["collectri"],
        dorothea=REGULONS["dorothea"]
    output:
        expand("data/processed/{{cancer}}/features/{{cohort}}.{rep}.tsv", rep=REPS),
        counts="outputs/real/{cancer}/regulon_target_counts_{cohort}.tsv"
    shell:
        "python scripts/build_features.py --expression {input.expression} --collectri {input.collectri} "
        "--dorothea {input.dorothea} --min-targets {R[min_regulon_targets]} "
        "--out-prefix data/processed/{wildcards.cancer}/features/{wildcards.cohort} --target-counts {output.counts}"


rule real_qc:
    input:
        expression=lambda w: f"data/processed/{ANALYSES[w.analysis]['cancer']}/harmonized/{w.cohort}.expression.tsv",
        clinical=lambda w: _clinical(ANALYSES[w.analysis]["cancer"], w.cohort,
                                     ANALYSES[w.analysis]["endpoint"] if w.cohort == _training(w.analysis) else ANALYSES[w.analysis]["validation"][w.cohort])
    output:
        "outputs/real/{analysis}/qc_{cohort}.tsv"
    params:
        platform=lambda w: CANCERS[ANALYSES[w.analysis]["cancer"]]["cohorts"][w.cohort]["platform"],
        endpoint=lambda w: ANALYSES[w.analysis]["endpoint"] if w.cohort == _training(w.analysis) else ANALYSES[w.analysis]["validation"][w.cohort]
    shell:
        "python scripts/characterize_cohort.py --expression {input.expression} --clinical {input.clinical} "
        "--cohort {wildcards.cohort} --platform {params.platform} --endpoint {params.endpoint} "
        "--min-events {config[min_confirmatory_events]} --subset-to-clinical --output {output}"


rule real_panel:
    input:
        features=_train_features,
        clinical=_train_clinical
    output:
        "outputs/real/{analysis}/{rep}_panel.tsv"
    shell:
        "python scripts/select_panels.py --features {input.features} --clinical {input.clinical} --k {PANEL} --label-source TCGA --output {output}"


rule real_stability:
    input:
        features=_train_features,
        clinical=_train_clinical
    output:
        "outputs/real/{analysis}/{rep}_stability.tsv"
    shell:
        "python scripts/panel_stability.py --features {input.features} --clinical {input.clinical} --representation {wildcards.rep} "
        "--analysis {wildcards.analysis} --k {PANEL} --n-resamples {config[n_stability_resamples]} --seed {config[random_seed]} --output {output}"


rule real_models:
    input:
        features=_train_features,
        clinical=_train_clinical,
        panel="outputs/real/{analysis}/{rep}_panel.tsv"
    output:
        "outputs/real/{analysis}/{rep}_models.joblib"
    threads: 4
    shell:
        "python scripts/fit_models.py --features {input.features} --clinical {input.clinical} --panel {input.panel} "
        "--label-source TCGA --primary-l1-ratio {config[primary_l1_ratio]} --lasso-l1-ratio {config[lasso_l1_ratio]} "
        "--cv-splits {config[cv_splits]} --rsf-n-estimators {config[rsf_n_estimators]} --seed {config[random_seed]} --output {output}"


rule real_compare:
    input:
        models=lambda w: expand("outputs/real/{analysis}/{rep}_models.joblib", analysis=w.analysis, rep=REPS),
        train=lambda w: [f"data/processed/{ANALYSES[w.analysis]['cancer']}/features/{_training(w.analysis)}.{r}.tsv" for r in REPS],
        valid=lambda w: [f"data/processed/{ANALYSES[w.analysis]['cancer']}/features/{w.cohort}.{r}.tsv" for r in REPS],
        train_clinical=_train_clinical,
        valid_clinical=lambda w: _clinical(ANALYSES[w.analysis]["cancer"], w.cohort, ANALYSES[w.analysis]["validation"][w.cohort]),
        valid_qc="outputs/real/{analysis}/qc_{cohort}.tsv"
    output:
        metrics="outputs/real/{analysis}/{cohort}_metrics.tsv",
        comparisons="outputs/real/{analysis}/{cohort}_comparisons.tsv",
        risks="outputs/real/{analysis}/{cohort}_risks.tsv"
    params:
        reps=" ".join(REPS),
        train_ep=lambda w: ANALYSES[w.analysis]["endpoint"],
        valid_ep=lambda w: ANALYSES[w.analysis]["validation"][w.cohort]
    shell:
        "python scripts/compare_representations.py --config config/config.yaml --representations {params.reps} "
        "--models {input.models} --train-features {input.train} --valid-features {input.valid} "
        "--train-clinical {input.train_clinical} --valid-clinical {input.valid_clinical} "
        "--train-endpoint {params.train_ep} --valid-endpoint {params.valid_ep} --valid-qc {input.valid_qc} "
        "--analysis {wildcards.analysis} --cohort {wildcards.cohort} --n-boot {config[n_bootstrap]} "
        "--metrics-output {output.metrics} --comparisons-output {output.comparisons} --risks-output {output.risks}"


def _real_targets(kind):
    out = []
    for a, s in ANALYSES.items():
        if kind == "qc":
            out += [f"outputs/real/{a}/qc_{c}.tsv" for c in [_training(a), *s["validation"]]]
        elif kind in ("metrics", "comparisons", "risks"):
            out += [f"outputs/real/{a}/{c}_{kind}.tsv" for c in s["validation"]]
        elif kind == "stability":
            out += [f"outputs/real/{a}/{r}_stability.tsv" for r in REPS]
    return out


rule real_meta:
    input:
        _real_targets("comparisons")
    output:
        "outputs/real/meta_analysis.tsv"
    shell:
        "python scripts/meta_analysis.py --config config/config.yaml --comparisons {input} --output {output}"


FIGURES = ["fig_design", "fig_confirmatory_forest", "fig_external_cindex", "fig_added_value", "fig_stability", "fig_km", "fig_km_gene"]


rule real_figures:
    input:
        metrics=_real_targets("metrics"),
        comparisons=_real_targets("comparisons"),
        stability=_real_targets("stability"),
        risks=_real_targets("risks"),
        qc=_real_targets("qc"),
        meta="outputs/real/meta_analysis.tsv"
    output:
        expand("outputs/figures/{fig}.{ext}", fig=FIGURES, ext=["pdf", "png"])
    params:
        reps=" ".join(REPS)
    shell:
        "python scripts/make_figures.py --config config/config.yaml --metrics {input.metrics} --comparisons {input.comparisons} "
        "--meta {input.meta} --stability {input.stability} --risks {input.risks} --qc {input.qc} --representations {params.reps} --outdir outputs/figures"


rule real_report:
    input:
        qc=_real_targets("qc"),
        metrics=_real_targets("metrics"),
        comparisons=_real_targets("comparisons"),
        stability=_real_targets("stability"),
        meta="outputs/real/meta_analysis.tsv"
    output:
        "outputs/reports/real_summary.md"
    shell:
        "python scripts/build_report.py --real-data --qc {input.qc} --metrics {input.metrics} "
        "--comparisons {input.comparisons} --stability {input.stability} --meta {input.meta} --output {output}"


rule real:
    input:
        "outputs/reports/real_summary.md",
        expand("outputs/figures/{fig}.pdf", fig=FIGURES)


# ------------------------------------------------------------- transfer gap (TCGA vs GEO)
# Separate target, so adding it never re-runs the analyses above:
#   snakemake --cores 4 real_drop

rule real_internal_cv:
    input:
        features=_train_features,
        clinical=_train_clinical
    output:
        summary="outputs/real/{analysis}/{rep}_internal_cv.tsv",
        folds="outputs/real/{analysis}/{rep}_internal_cv_folds.tsv"
    shell:
        "python scripts/internal_cv.py --config config/config.yaml --features {input.features} --clinical {input.clinical} "
        "--representation {wildcards.rep} --analysis {wildcards.analysis} --k {PANEL} "
        "--repeats {config[internal_cv_repeats]} --output {output.summary} --folds-output {output.folds}"


rule real_drop:
    input:
        internal=[f"outputs/real/{a}/{r}_internal_cv.tsv" for a in ANALYSES for r in REPS],
        metrics=_real_targets("metrics")
    output:
        table="outputs/real/performance_drop.tsv",
        fig=expand("outputs/figures/fig_performance_drop.{ext}", ext=["pdf", "png"])
    params:
        reps=" ".join(REPS)
    shell:
        "python scripts/performance_drop.py --internal {input.internal} --metrics {input.metrics} "
        "--representations {params.reps} --output {output.table} --figure-dir outputs/figures"
