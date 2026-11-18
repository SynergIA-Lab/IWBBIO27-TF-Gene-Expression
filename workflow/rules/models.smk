rule fixture_panel:
    input:
        features="data/processed/fixture/features/tcga.{rep}.tsv",
        clinical="tests/fixtures/clinical_tcga.tsv"
    output:
        "outputs/tables/models/fixture_{rep}_panel.tsv"
    shell:
        "python scripts/select_panels.py --features {input.features} --clinical {input.clinical} --k {FX[panel_size]} --label-source TCGA --output {output}"

rule fixture_stability:
    input:
        features="data/processed/fixture/features/tcga.{rep}.tsv",
        clinical="tests/fixtures/clinical_tcga.tsv"
    output:
        "outputs/tables/models/fixture_{rep}_stability.tsv"
    shell:
        "python scripts/panel_stability.py --features {input.features} --clinical {input.clinical} --representation {wildcards.rep} --analysis fixture --k {FX[panel_size]} --n-resamples {FX[n_stability_resamples]} --seed {config[random_seed]} --output {output}"

rule fixture_models:
    input:
        features="data/processed/fixture/features/tcga.{rep}.tsv",
        clinical="tests/fixtures/clinical_tcga.tsv",
        panel="outputs/tables/models/fixture_{rep}_panel.tsv"
    output:
        "outputs/tables/models/fixture_{rep}_models.joblib"
    shell:
        "python scripts/fit_models.py --features {input.features} --clinical {input.clinical} --panel {input.panel} --label-source TCGA --primary-l1-ratio {config[primary_l1_ratio]} --lasso-l1-ratio {config[lasso_l1_ratio]} --cv-splits {config[cv_splits]} --rsf-n-estimators {FX[rsf_n_estimators]} --seed {config[random_seed]} --output {output}"
