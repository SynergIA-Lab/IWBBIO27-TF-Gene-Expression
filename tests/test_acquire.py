import io
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from regbio.acquire import (
    apply_clinical_mapping,
    cdr_clinical_table,
    characteristics_overview,
    gdc_file_table,
    gpl_gene_map,
    needs_log2,
    one_sample_per_patient,
    parse_star_counts,
    resolve_column,
)

STAR = """# gene-model: GENCODE v36
gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded
N_unmapped\t\t\t100\t100\t100\t\t\t
N_multimapping\t\t\t50\t50\t50\t\t\t
ENSG01.1\tTP53\tprotein_coding\t10\t5\t5\t12.5\t1\t1
ENSG02.1\tLINC1\tlncRNA\t10\t5\t5\t3.0\t1\t1
ENSG03.1\tPAR1\tprotein_coding\t10\t5\t5\t2.0\t1\t1
ENSG03.2_PAR_Y\tPAR1\tprotein_coding\t10\t5\t5\t4.0\t1\t1
"""


def test_parse_star_counts_keeps_protein_coding_and_averages_duplicates():
    s = parse_star_counts(io.StringIO(STAR))
    assert set(s.index) == {"TP53", "PAR1"}
    assert s["TP53"] == 12.5 and s["PAR1"] == 3.0


def test_gdc_file_table_and_one_sample_per_patient():
    hits = [
        {"file_id": "f2", "file_name": "b.tsv", "cases": [{"submitter_id": "TCGA-AA-0001", "samples": [{"submitter_id": "TCGA-AA-0001-01B", "sample_type": "Primary Tumor"}]}]},
        {"file_id": "f1", "file_name": "a.tsv", "cases": [{"submitter_id": "TCGA-AA-0001", "samples": [{"submitter_id": "TCGA-AA-0001-01A", "sample_type": "Primary Tumor"}]}]},
        {"file_id": "f3", "file_name": "c.tsv", "cases": [{"submitter_id": "TCGA-AA-0002", "samples": [{"submitter_id": "TCGA-AA-0002-01A", "sample_type": "Primary Tumor"}]}]},
    ]
    kept = one_sample_per_patient(gdc_file_table(hits))
    assert list(kept["file_id"]) == ["f1", "f3"]


def test_cdr_clinical_table_filters_cancer_and_parses_stage():
    cdr = pd.DataFrame({
        "bcr_patient_barcode": ["P1", "P2", "P3", "P4"],
        "type": ["LUAD", "LUAD", "BRCA", "LUAD"],
        "PFI": [1, 0, 1, 1],
        "PFI.time": [300, 1200, 50, "#N/A"],
        "age_at_initial_pathologic_diagnosis": [65, 70, 50, 60],
        "ajcc_pathologic_tumor_stage": ["Stage IB", "Stage IIIA", "Stage I", "[Not Available]"],
    })
    out = cdr_clinical_table(cdr, "LUAD", "PFI")
    assert list(out["sample_id"]) == ["P1", "P2"]  # BRCA excluded, P4 has no time
    assert list(out["stage"]) == ["I", "III"]
    assert list(out["event"]) == [1, 0]


def test_gpl_gene_map_drops_ambiguous_probes():
    gpl = pd.DataFrame({"ID": ["p1", "p2", "p3", "p4"], "Gene Symbol": ["TP53", "A /// B", "---", "EGFR"]})
    gm, col = gpl_gene_map(gpl)
    assert col == "Gene Symbol"
    assert dict(zip(gm["feature_id"], gm["gene_id"])) == {"p1": "TP53", "p4": "EGFR"}


def test_gpl_gene_map_parses_gene_assignment():
    gpl = pd.DataFrame({"ID": ["p1", "p2"], "gene_assignment": ["NM_1 // KRAS // desc // 12p // 3845", "---"]})
    gm, col = gpl_gene_map(gpl)
    assert col == "gene_assignment" and gm.iloc[0]["gene_id"] == "KRAS" and len(gm) == 1


def test_gpl_gene_map_errors_without_symbol_column():
    with pytest.raises(ValueError, match="no gene symbol column"):
        gpl_gene_map(pd.DataFrame({"ID": ["p1"], "SEQ": ["ACGT"]}))


def test_needs_log2():
    assert needs_log2(pd.DataFrame({"a": [5.0, 12000.0, 300.0] * 50}))
    assert not needs_log2(pd.DataFrame({"a": [5.0, 12.0, 3.0] * 50}))


def _pheno():
    return pd.DataFrame({
        "characteristics_ch1.0.tissue": ["primary lung tumor", "primary lung tumor", "normal lung", "primary lung tumor", "primary lung tumor"],
        "characteristics_ch1.3.months to relapse": ["12", "30", "", "6", "8"],
        "characteristics_ch1.4.relapse": ["relapsed", "not relapsed", "", "Relapsed", "unknown"],
        "characteristics_ch1.5.age": ["61", "70", "55", "48", "50"],
        "characteristics_ch1.6.stage": ["IA", "Stage 2B", "", "IB", "II"],
    }, index=["GSM1", "GSM2", "GSM3", "GSM4", "GSM5"])


def test_resolve_column_by_key_suffix():
    assert resolve_column(_pheno(), "relapse") == "characteristics_ch1.4.relapse"
    with pytest.raises(ValueError):
        resolve_column(_pheno(), "nonexistent")


def test_apply_clinical_mapping_filters_converts_and_drops_unknown_events():
    mapping = {
        "endpoint": "RFS",
        "time_column": "months to relapse",
        "time_unit": "months",
        "event_column": "relapse",
        "event_values": ["relapsed"],
        "censor_values": ["not relapsed"],
        "age_column": "age",
        "stage_column": "stage",
        "sample_filter": {"tissue": ["primary lung tumor"]},
    }
    out = apply_clinical_mapping(_pheno(), mapping)
    assert list(out["sample_id"]) == ["GSM1", "GSM2", "GSM4"]  # GSM3 normal, GSM5 unknown event
    assert out["time"].tolist() == pytest.approx([12 * 30.4375, 30 * 30.4375, 6 * 30.4375])
    assert out["event"].tolist() == [1, 0, 1]  # case-insensitive event values
    assert out["stage"].tolist() == ["I", "II", "I"]  # arabic "Stage 2B" handled
    assert out.attrs["n_dropped"] == 1


def test_apply_clinical_mapping_strips_key_prefixes():
    pheno = pd.DataFrame({"c.time": ["time: 100", "time: 200"], "c.status": ["status: dead", "status: alive"]}, index=["A", "B"])
    out = apply_clinical_mapping(pheno, {"endpoint": "OS", "time_column": "time", "event_column": "status", "event_values": ["dead"], "censor_values": ["alive"]})
    assert out["time"].tolist() == [100.0, 200.0] and out["event"].tolist() == [1, 0]


def test_characteristics_overview_lists_examples():
    ov = characteristics_overview(_pheno())
    assert "characteristics_ch1.4.relapse" in set(ov["column"])


def test_build_geo_clinical_refuses_draft_mapping(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    _pheno().to_csv(raw / "pheno_raw.tsv", sep="\t")
    pd.DataFrame({"GSM1": [1.0]}, index=["p1"]).to_csv(raw / "expression_probes.tsv", sep="\t")
    maps = tmp_path / "m.yaml"
    maps.write_text("GSEX:\n  status: DRAFT_VERIFY\n  endpoint: RFS\n  time_column: x\n  event_column: y\n  event_values: [a]\n  censor_values: [b]\n")
    res = subprocess.run([sys.executable, "scripts/build_geo_clinical.py", "--gse", "GSEX", "--rawdir", str(raw), "--mappings", str(maps)], capture_output=True, text=True)
    assert res.returncode != 0 and "DRAFT" in (res.stderr + res.stdout)


def test_build_geo_clinical_aligns_expression(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    _pheno().to_csv(raw / "pheno_raw.tsv", sep="\t")
    pd.DataFrame(np.ones((2, 5)), index=["p1", "p2"], columns=["GSM1", "GSM2", "GSM3", "GSM4", "GSM5"]).to_csv(raw / "expression_probes.tsv", sep="\t")
    maps = tmp_path / "m.yaml"
    maps.write_text(
        "GSEX:\n  status: VERIFIED\n  endpoint: RFS\n  time_column: months to relapse\n  time_unit: months\n"
        "  event_column: relapse\n  event_values: [relapsed]\n  censor_values: [not relapsed]\n"
        "  sample_filter:\n    tissue: [primary lung tumor]\n"
    )
    subprocess.run([sys.executable, "scripts/build_geo_clinical.py", "--gse", "GSEX", "--rawdir", str(raw), "--mappings", str(maps)], check=True, capture_output=True)
    expr = pd.read_csv(raw / "expression.tsv", sep="\t", index_col=0)
    clin = pd.read_csv(raw / "clinical.tsv", sep="\t")
    assert list(expr.columns) == list(clin["sample_id"]) == ["GSM1", "GSM2", "GSM4"]


def test_download_retries_then_succeeds():
    import importlib.util
    import urllib.error

    spec = importlib.util.spec_from_file_location("download_tcga", "scripts/download_tcga.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    calls = {"n": 0}

    def flaky(x):
        calls["n"] += 1
        if calls["n"] < 3:
            raise urllib.error.URLError("RemoteDisconnected")
        return x * 2

    assert mod.with_retries(flaky, 21, attempts=5, sleep=lambda s: None) == 42
    assert calls["n"] == 3


def test_download_retries_give_up():
    import importlib.util
    import urllib.error

    spec = importlib.util.spec_from_file_location("download_tcga", "scripts/download_tcga.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    def always_fails():
        raise urllib.error.URLError("down")

    with pytest.raises(urllib.error.URLError):
        mod.with_retries(always_fails, attempts=3, sleep=lambda s: None)


def _gse31210_like():
    # Same key split over different column positions, as GEOparse produces for GSE31210.
    return pd.DataFrame({
        "characteristics_ch1.0.tissue": ["primary lung tumor", "primary lung tumor", "normal lung", "primary lung tumor", "primary lung tumor"],
        "characteristics_ch1.1.age (years)": ["55", None, "60", "70", None],
        "characteristics_ch1.1.age": [None, "64", None, None, "48"],
        "characteristics_ch1.5.pathological stage": ["IA", "II", None, "IB", "IA"],
        "characteristics_ch1.10.relapse": ["relapsed", None, None, None, None],
        "characteristics_ch1.11.relapse": [None, "not relapsed", None, "relapsed", None],
        "characteristics_ch1.12.relapse": [None, None, None, None, "not relapsed"],
        "characteristics_ch1.11.days before relapse/censor": ["253", None, None, None, None],
        "characteristics_ch1.12.days before relapse/censor": [None, "2601", None, "438", None],
        "characteristics_ch1.13.days before relapse/censor": [None, None, None, None, "1638"],
        "characteristics_ch1.16.exclude for prognosis analysis due to incomplete resection or adjuvant therapy": ["none", "none", None, "exclude", None],
        "characteristics_ch1.18.exclude for prognosis analysis due to incomplete resection or adjuvant therapy": [None, None, None, None, "none"],
    }, index=["S1", "S2", "N1", "S3", "S4"])


def test_split_characteristics_are_coalesced_like_gse31210():
    import yaml

    mapping = yaml.safe_load(open("config/clinical_mappings.yaml"))["GSE31210"]
    out = apply_clinical_mapping(_gse31210_like(), mapping)
    assert list(out["sample_id"]) == ["S1", "S2", "S4"]  # normal and excluded samples removed
    assert out["time"].tolist() == [253.0, 2601.0, 1638.0]
    assert out["event"].tolist() == [1, 0, 0]
    assert out["age"].tolist() == [55.0, 64.0, 48.0]
    assert out["stage"].tolist() == ["I", "II", "I"]


def test_conflicting_split_values_raise():
    from regbio.acquire import coalesce_characteristic

    pheno = pd.DataFrame({"c.1.months before relapse/censor": ["10"], "c.2.months before relapse/censor": ["20"]}, index=["S1"])
    with pytest.raises(ValueError, match="conflicting"):
        coalesce_characteristic(pheno, "months before relapse/censor")


def test_gse72094_mapping_drops_na_vital_status():
    import yaml

    mapping = yaml.safe_load(open("config/clinical_mappings.yaml"))["GSE72094"]
    pheno = pd.DataFrame({
        "characteristics_ch1.5.age_at_diagnosis": ["70", "59", "61"],
        "characteristics_ch1.9.vital_status": ["Dead", "Alive", "NA"],
        "characteristics_ch1.10.survival_time_in_days": ["895", "1249", "NA"],
        "characteristics_ch1.11.Stage": ["2B", "1A", "NA"],
    }, index=["A", "B", "C"])
    out = apply_clinical_mapping(pheno, mapping)
    assert list(out["sample_id"]) == ["A", "B"]
    assert out["event"].tolist() == [1, 0]
    assert out["stage"].tolist() == ["II", "I"]


SERIES_MATRIX = """!Series_title\t"x"
!Sample_title\t"F1"\t"F2"\t"F1repl"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_characteristics_ch1\t"age at diagnosis: 60"\t"age at diagnosis: 71"\t"age at diagnosis: 60"
!Sample_characteristics_ch1\t"overall survival days: 1200"\t"lymph node status: NodeNegative"\t"overall survival days: 1200"
!Sample_characteristics_ch1\t"overall survival event: 1"\t"overall survival days: 900"\t"overall survival event: 1"
!Sample_characteristics_ch1\t""\t"overall survival event: 0"\t""
!series_matrix_table_begin
"ID_REF"\t"GSM1"
"""


def test_parse_series_matrix_pheno_and_mapping_roundtrip():
    from regbio.acquire import apply_clinical_mapping, map_columns_to_gsm, parse_series_matrix_pheno

    pheno = parse_series_matrix_pheno(SERIES_MATRIX.splitlines(True))
    assert list(pheno.index) == ["GSM1", "GSM2", "GSM3"]
    assert pheno.loc["GSM2", "title"] == "F2"
    clin = apply_clinical_mapping(pheno, {
        "endpoint": "OS", "time_column": "overall survival days", "time_unit": "days",
        "event_column": "overall survival event", "event_values": ["1"], "censor_values": ["0"],
        "age_column": "age at diagnosis",
    }).set_index("sample_id")
    assert clin.loc["GSM2", "time"] == 900 and clin.loc["GSM2", "event"] == 0
    assert clin.loc["GSM1", "event"] == 1 and clin.loc["GSM1", "age"] == 60

    expr = pd.DataFrame({"F1": [1.0], "F2": [2.0], "F1repl": [1.1], "F9": [0.0]}, index=["TP53"])
    out, stats = map_columns_to_gsm(expr, pheno, "repl")
    assert list(out.columns) == ["GSM1", "GSM2"]
    assert stats == {"n_columns": 4, "n_dropped_regex": 1, "n_unmatched": 1, "n_kept": 2}


def test_links_in_listing_filters_names():
    from regbio.acquire import links_in_listing

    html = '<a href="../">..</a><a href="GSE1-GPL1_series_matrix.txt.gz">a</a><a href="GSE1_RAW.tar">b</a>'
    assert links_in_listing(html, r"_series_matrix\.txt\.gz$") == ["GSE1-GPL1_series_matrix.txt.gz"]


def test_breast_stage_from_tnm():
    from regbio.acquire import breast_stage_from_tnm as st

    assert st("1c", "0", "0") == "I"
    assert st("1", "1", "0") == "II" and st("2", "0", "0") == "II" and st("3", "0", "0") == "II"
    assert st("2", "2", "0") == "III" and st("3", "1", "0") == "III" and st("4", "0", "0") == "III" and st("1", "3", "0") == "III"
    assert st("2", "1", "1") == "IV"
    assert pd.isna(st("NA", "1", "0")) and pd.isna(st("2", None, "0"))


def test_gse20685_style_dmfs_uses_event_time_for_events_only():
    from regbio.acquire import apply_clinical_mapping

    pheno = pd.DataFrame({
        "characteristics_ch1.2.age at diagnosis": ["50", "60", "45", "70"],
        "characteristics_ch1.4.follow_up_duration (years)": ["7", "8", "6", "5"],
        "characteristics_ch1.5.event_metastasis": ["0", "1", "1", "1"],
        "characteristics_ch1.6.time_to_metastasis (years)": [None, "2", None, "0"],
        "characteristics_ch1.7.t_stage": ["1c", "2", "3", "2"],
        "characteristics_ch1.8.n_stage": ["0", "1", "2", "0"],
        "characteristics_ch1.9.m_stage": ["0", "0", "0", "1"],
    }, index=["A", "B", "C", "D"])
    clin = apply_clinical_mapping(pheno, {
        "endpoint": "DMFS", "time_column": "follow_up_duration (years)", "time_unit": "years",
        "event_time_column": "time_to_metastasis (years)",
        "event_column": "event_metastasis", "event_values": ["1"], "censor_values": ["0"],
        "age_column": "age at diagnosis", "sample_filter": {"m_stage": ["0"]},
        "stage_from_tnm": {"t": "t_stage", "n": "n_stage", "m": "m_stage"},
    }).set_index("sample_id")
    assert list(clin.index) == ["A", "B"]          # D is M1; C is an event without event time
    assert clin.loc["A", "time"] == 7 * 365.25 and clin.loc["A", "event"] == 0
    assert clin.loc["B", "time"] == 2 * 365.25 and clin.loc["B", "event"] == 1
    assert clin.loc["A", "stage"] == "I" and clin.loc["B", "stage"] == "II"


def test_build_geo_clinical_keeps_numeric_codes_as_text(tmp_path):
    # split characteristic columns contain NaN; read back as float they become '0.0'/'1.0'
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.DataFrame({
        "characteristics_ch1.4.follow_up_duration (years)": ["7", "8", "6"],
        "characteristics_ch1.5.event_metastasis": ["0", "1", "0"],
        "characteristics_ch1.6.time_to_metastasis (years)": [None, "2", None],
        "characteristics_ch1.9.m_stage": ["0", None, "1"],
        "characteristics_ch1.10.m_stage": [None, "0", None],
    }, index=["GSM1", "GSM2", "GSM3"]).to_csv(raw / "pheno_raw.tsv", sep="\t")
    pd.DataFrame(np.ones((1, 3)), index=["g"], columns=["GSM1", "GSM2", "GSM3"]).to_csv(raw / "expression_probes.tsv", sep="\t")
    maps = tmp_path / "m.yaml"
    maps.write_text(
        "GSEX:\n  status: VERIFIED\n  endpoint: DMFS\n  time_column: follow_up_duration (years)\n  time_unit: years\n"
        "  event_time_column: time_to_metastasis (years)\n"
        "  event_column: event_metastasis\n  event_values: ['1']\n  censor_values: ['0']\n"
        "  sample_filter:\n    m_stage: ['0']\n"
    )
    subprocess.run([sys.executable, "scripts/build_geo_clinical.py", "--gse", "GSEX", "--rawdir", str(raw), "--mappings", str(maps)], check=True, capture_output=True)
    clin = pd.read_csv(raw / "clinical.tsv", sep="\t")
    assert list(clin["sample_id"]) == ["GSM1", "GSM2"] and list(clin["event"]) == [0, 1]


def test_stage_from_size_and_node_group_like_gse96058():
    from regbio.acquire import apply_clinical_mapping

    pheno = pd.DataFrame({
        "characteristics_ch1.21.overall survival days": ["2000", "1500", "900", "800"],
        "characteristics_ch1.22.overall survival event": ["0", "1", "0", "1"],
        "characteristics_ch1.3.tumor size": ["15", "30", "60", "NA"],
        "characteristics_ch1.4.lymph node group": ["NodeNegative", "1to3", "4toX", "1to3"],
    }, index=["A", "B", "C", "D"])
    clin = apply_clinical_mapping(pheno, {
        "endpoint": "OS", "time_column": "overall survival days", "time_unit": "days",
        "event_column": "overall survival event", "event_values": ["1"], "censor_values": ["0"],
        "stage_from_tnm": {"t_size_mm": "tumor size", "n": "lymph node group",
                           "n_map": {"NodeNegative": "0", "SubMicroMet": "0", "1to3": "1", "4toX": "2"}},
    }).set_index("sample_id")
    assert clin.loc["A", "stage"] == "I" and clin.loc["B", "stage"] == "II" and clin.loc["C", "stage"] == "III"
    assert pd.isna(clin.loc["D", "stage"])
