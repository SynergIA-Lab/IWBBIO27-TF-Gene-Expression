rule fixture_features:
    input:
        expression="data/processed/fixture/{cohort}_expression.tsv",
        collectri="tests/fixtures/collectri_fixture.tsv",
        dorothea="tests/fixtures/dorothea_fixture.tsv"
    output:
        expand("data/processed/fixture/features/{{cohort}}.{rep}.tsv", rep=FIXTURE_REPS),
        counts="outputs/tables/qc/fixture_{cohort}_regulon_target_counts.tsv"
    wildcard_constraints:
        cohort="tcga|geo"
    shell:
        "python scripts/build_features.py --expression {input.expression} --collectri {input.collectri} --dorothea {input.dorothea} --min-targets {FX[min_regulon_targets]} --out-prefix data/processed/fixture/features/{wildcards.cohort} --target-counts {output.counts}"
