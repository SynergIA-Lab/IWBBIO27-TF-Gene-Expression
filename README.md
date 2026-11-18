# Regulatory Biomarker Transfer Pipeline

Compares prognostic biomarkers built from gene expression, TF expression and
regulatory-network-inferred TF activity (CollecTRI+ULM, DoRothEA+VIPER, optional
ARACNe+VIPER) when models trained on TCGA are transferred to independent GEO cohorts.

## Milestones

1. Fixture MVP with real engines (done): decoupler for activity, scikit-survival for
   models and metrics, synthetic fixtures on two measurement scales.
2. Conference-scale execution: LUAD first (TCGA + GSE31210 + GSE72094), then BRCA and COAD.

## Run

```bash
python scripts/check_dependencies.py
python -m pytest -q            # 71 tests
snakemake -n                   # dry run
snakemake --cores 4 --forceall # fixture workflow -> outputs/reports/fixture_summary.md
```

`scripts/make_fixtures.py` regenerates the synthetic fixtures deterministically.

## Phase 3: get the LUAD data

```bash
# TCGA-LUAD: STAR counts from GDC (open access) + TCGA-CDR outcomes
python scripts/download_tcga.py --project TCGA-LUAD --cancer LUAD --endpoint PFI --outdir data/raw/TCGA_LUAD

# GEO cohorts: expression, probe->gene map, raw phenotypes
python scripts/download_geo.py --gse GSE31210 --outdir data/raw/GSE31210
python scripts/download_geo.py --gse GSE72094 --outdir data/raw/GSE72094

# Inspect data/raw/<GSE>/characteristics_overview.tsv, fix config/clinical_mappings.yaml,
# set status: VERIFIED, then build the curated clinical tables:
python scripts/build_geo_clinical.py --gse GSE31210 --rawdir data/raw/GSE31210
python scripts/build_geo_clinical.py --gse GSE72094 --rawdir data/raw/GSE72094
```

Check each `acquisition_log.json` (platform, symbol column, genes mapped, log2 applied).
The clinical mappings ship as DRAFT and the builder refuses to run until they are verified.

## Phase 4: run LUAD

```bash
python scripts/download_tcga.py --project TCGA-LUAD --cancer LUAD --endpoint OS --outdir data/raw/TCGA_LUAD  # also PFI
snakemake real -n                 # dry run: lists every step
snakemake real --cores 4          # -> outputs/reports/real_summary.md + outputs/figures/*.pdf|png
```

Analyses are declared in `config/config.yaml` under `real:` (LUAD_PFI validated on GSE31210
RFS; LUAD_OS validated on GSE72094 OS). The first run downloads CollecTRI and DoRothEA
(`resources/regulons/`). `tests/real_workflow_test.yaml` runs the same rules on fixtures:
`snakemake real --cores 4 --configfile tests/real_workflow_test.yaml`.

## Method notes

- Harmonization: common gene universe locked across cohorts, then cohort-wise z-score
  (default) or per-sample ranks. No TCGA scaling parameter is ever applied to GEO.
- Panel: top-k features by univariate Cox score test (Breslow ties), TCGA only.
  Stability: mean pairwise Kuncheva index of the top-k panel over TCGA bootstrap resamples.
- Models: primary ridge-like Coxnet (l1_ratio 0.05), LASSO sensitivity (0.95), RSF secondary.
  Penalty chosen by event-stratified CV on TCGA only; the artifact records
  `alpha_selection` = `cv` or `fallback_no_informative_folds` (the latter must never
  appear with real cohorts).
- `tf_expression` uses exactly the TF set of `collectri_ulm`: CollecTRI TFs passing the
  common-gene and `min_regulon_targets` filters whose own expression is measured.
- Confirmatory inference (amendment 2026-10-28, `docs/superpowers/specs/2026-10-28-design-amendment-confirmatory-inference.md`):
  fixed-effect pooling of the paired C-index differences of all eligible cohorts for the two
  pre-specified comparisons, Holm-adjusted; per-cohort results are descriptive.
- Evaluation: frozen TCGA models applied to validation cohorts; IPCW C-index (censoring
  estimated on the validation cohort, truncated at `cindex_tau_days`, 5 years) and
  cumulative/dynamic AUC; paired bootstrap of the C-index difference between
  representations on the same patients. Only the pre-specified comparisons are
  confirmatory, and only if endpoints align and the cohort has >= min_confirmatory_events.
- Clinical baseline: Cox on age + ordinal stage; clinical+biomarker adds the out-of-fold
  TCGA risk score of each representation (complete cases); added value is reported as
  exploratory paired comparisons.
- Real-data harmonization (`harmonize_cancer.py`): all cohorts of a cancer at once, samples
  restricted to those with clinical data, low-expression filter on TCGA only, one common
  gene universe, then cohort-wise z-score.
- `aracne_activity` is built only when a cancer-specific ARACNe regulon is supplied;
  there is no placeholder.
- VIPER is used through decoupler; check licensing before non-academic use.

## Limitations

Fixture results are not biological evidence. Cohort-wise z-scoring can shift scores when
validation cohorts differ in stage or subtype composition. Real runs must inspect cohort
QC, endpoint alignment, event counts and confirmatory/exploratory labels before claims.
