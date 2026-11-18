# Design Amendment: Confirmatory Inference And Figures

Amends: docs/superpowers/specs/2026-10-26-regulatory-biomarker-pipeline-design.md

Status: adopted on 2026-10-28, after the LUAD analysis (TCGA-LUAD -> GSE31210, TCGA-LUAD -> GSE72094) had been run and inspected, and before any BRCA or COAD validation result exists. The paper must state this timing.

Reason: with three cancers and about six validation cohorts, testing the two pre-specified comparisons separately in every cohort gives up to twelve tests. Without correction the family-wise false-positive rate is uncontrolled; with Holm over twelve per-cohort tests the per-cohort power is too low to detect a consistent effect.

Primary confirmatory analysis (replaces per-cohort confirmatory tests):

- For each pre-specified comparison (`collectri_ulm - gene_expression`, `collectri_ulm - tf_expression`), the paired C-index differences of all confirmatory-eligible validation cohorts (aligned endpoints and at least `min_confirmatory_events`) are pooled by fixed-effect inverse-variance meta-analysis, using the paired-bootstrap standard error of each cohort.
- The two pooled p-values are adjusted with Holm at a family-wise alpha of 0.05. These two tests are the only confirmatory claims of the study.
- Per-cohort paired differences are reported with 95% bootstrap intervals as descriptive, supporting evidence.

Pooling across endpoints: the Evaluation section forbade pooling metrics across incompatible endpoints. That rule protects absolute C-indices. The pooled quantity here is the within-patient difference between two representations, a relative measure of discrimination that is comparable across endpoint classes. Pooling RFS/PFI and OS cohorts is therefore allowed for this quantity only, and is checked by a sensitivity analysis that pools within endpoint class (progression-related vs OS).

Sensitivity analyses: DerSimonian-Laird random effects with I^2 heterogeneity; pooling within endpoint class.

Stated limitation: intervals reflect sampling variability of validation patients conditional on the frozen TCGA models; training variability is not propagated.

Additional evaluation choices fixed in the same amendment:

- IPCW C-index: the censoring distribution is estimated on each validation cohort (its follow-up differs from TCGA), and C-indices are truncated at a common horizon `cindex_tau_days` = 1825 (5 years).
- Clinical baseline: Cox model on age and ordinal stage (I-IV -> 1-4), complete cases. The clinical-plus-biomarker model adds the out-of-fold TCGA risk score of each representation; added value is reported as exploratory paired comparisons.
- Risk groups for Kaplan-Meier figures use a cutoff frozen on TCGA (median of the training risk score), never re-estimated on validation cohorts.

Paper figures produced by the workflow (`outputs/figures/`): confirmatory forest plot (per-cohort differences and pooled estimate), external C-index with 95% CI per representation and clinical model, panel selection stability (Kuncheva), added value over the clinical model, and Kaplan-Meier curves of TCGA-median risk groups for gene expression and CollecTRI + ULM.

LUAD values at the time of adoption (for transparency, two cohorts only): pooled difference vs gene expression +0.025 (p = 0.25); vs TF expression +0.055 (p = 0.036, Holm-adjusted 0.071).
