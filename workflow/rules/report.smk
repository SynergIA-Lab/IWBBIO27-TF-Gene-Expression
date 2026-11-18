rule fixture_report:
    input:
        qc=["outputs/tables/qc/fixture_tcga.tsv", "outputs/tables/qc/fixture_geo.tsv"],
        metrics="outputs/tables/evaluation/fixture_metrics.tsv",
        comparisons="outputs/tables/evaluation/fixture_comparisons.tsv",
        stability=expand("outputs/tables/models/fixture_{rep}_stability.tsv", rep=FIXTURE_REPS)
    output:
        "outputs/reports/fixture_summary.md"
    shell:
        "python scripts/build_report.py --qc {input.qc} --metrics {input.metrics} --comparisons {input.comparisons} --stability {input.stability} --output {output}"
