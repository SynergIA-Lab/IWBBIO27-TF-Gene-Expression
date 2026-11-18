rule fixture_harmonized:
    input:
        tcga="tests/fixtures/expression_tcga.tsv",
        geo="tests/fixtures/expression_geo.tsv",
        gene_map="tests/fixtures/gene_map.tsv"
    output:
        tcga="data/processed/fixture/tcga_expression.tsv",
        geo="data/processed/fixture/geo_expression.tsv"
    shell:
        "python scripts/harmonize_cohorts.py --tcga-expression {input.tcga} --geo-expression {input.geo} --gene-map {input.gene_map} --method {config[harmonization]} --tcga-output {output.tcga} --geo-output {output.geo}"

rule download_tcga:
    output:
        touch("data/raw/.tcga_download_requested")
    shell:
        "python -c 'raise SystemExit(\"TCGA download is opt-in; configure GDC inputs before running full conference targets\")'"

rule download_geo:
    output:
        touch("data/raw/.geo_download_requested")
    shell:
        "python -c 'raise SystemExit(\"GEO download is opt-in; configure GEOparse inputs before running full conference targets\")'"
