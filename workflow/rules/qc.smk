rule fixture_qc_tcga:
    input:
        expression="tests/fixtures/expression_tcga.tsv",
        clinical="tests/fixtures/clinical_tcga.tsv"
    output:
        "outputs/tables/qc/fixture_tcga.tsv"
    shell:
        "python scripts/characterize_cohort.py --expression {input.expression} --clinical {input.clinical} --cohort fixture_tcga --platform rnaseq --endpoint PFI --min-events {FX[min_confirmatory_events]} --output {output}"

rule fixture_qc_geo:
    input:
        expression="tests/fixtures/expression_geo.tsv",
        clinical="tests/fixtures/clinical_geo.tsv"
    output:
        "outputs/tables/qc/fixture_geo.tsv"
    shell:
        "python scripts/characterize_cohort.py --expression {input.expression} --clinical {input.clinical} --cohort fixture_geo --platform microarray --endpoint RFS --min-events {FX[min_confirmatory_events]} --output {output}"
