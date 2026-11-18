configfile: "config/config.yaml"

FX = config["fixture"]
FIXTURE_REPS = FX["representations"]

include: "workflow/rules/data.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/features.smk"
include: "workflow/rules/models.smk"
include: "workflow/rules/evaluation.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/real.smk"

rule all:
    input:
        "outputs/reports/fixture_summary.md",
        "outputs/tables/evaluation/fixture_metrics.tsv",
        "outputs/tables/evaluation/fixture_comparisons.tsv"
