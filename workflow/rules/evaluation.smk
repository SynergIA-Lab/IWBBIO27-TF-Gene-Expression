rule fixture_compare:
    input:
        models=expand("outputs/tables/models/fixture_{rep}_models.joblib", rep=FIXTURE_REPS),
        train=expand("data/processed/fixture/features/tcga.{rep}.tsv", rep=FIXTURE_REPS),
        valid=expand("data/processed/fixture/features/geo.{rep}.tsv", rep=FIXTURE_REPS),
        train_clinical="tests/fixtures/clinical_tcga.tsv",
        valid_clinical="tests/fixtures/clinical_geo.tsv",
        qc_geo="outputs/tables/qc/fixture_geo.tsv"
    output:
        metrics="outputs/tables/evaluation/fixture_metrics.tsv",
        comparisons="outputs/tables/evaluation/fixture_comparisons.tsv"
    params:
        reps=" ".join(FIXTURE_REPS)
    shell:
        "python scripts/compare_representations.py --config config/config.yaml "
        "--representations {params.reps} --models {input.models} --train-features {input.train} --valid-features {input.valid} "
        "--train-clinical {input.train_clinical} --valid-clinical {input.valid_clinical} "
        "--train-endpoint PFI --valid-endpoint RFS "
        "--valid-qc {input.qc_geo} "
        "--n-boot {FX[n_bootstrap]} --metrics-output {output.metrics} --comparisons-output {output.comparisons}"
