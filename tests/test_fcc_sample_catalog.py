from pathlib import Path
import pytest
from modules.fcc_sample_catalog import load_catalog, resolve_records

REPO = Path(__file__).resolve().parents[1]

def test_catalog_keeps_p8c_and_p8o_distinct():
    data = load_catalog(REPO / "configs/analysis/fcc_sample_catalog_v1.yaml")
    assert data["samples"]["P8C"]["expected_events"] == 100000
    assert data["samples"]["P8O"]["expected_events"] == 18000
    assert data["comparisons"]["whizard_kkmcee_2k"]["W_source_id"] == "000242385"

def test_resolve_records(tmp_path, monkeypatch):
    manifest = tmp_path / "m.csv"; manifest.write_text("id,rec,direct,ancestor\n1,a,b,c\n2,d,e,f\n")
    sample = {"manifest": str(manifest), "columns": {"source_id":"id","source_rec":"rec","direct":"direct","ancestor":"ancestor"}, "include_source_ids":["2"]}
    assert resolve_records(sample) == [{"source_id":"2","source_rec":"d","direct":"e","ancestor":"f"}]
    monkeypatch.delenv("FCC_TAU_KKMCEE_MATERIAL", raising=False)
    with pytest.raises(ValueError, match="unresolved"):
        resolve_records({"single_source":{"source_id":"7","source_rec":"$FCC_TAU_KKMCEE_MATERIAL/x"}})

def test_generic_catalog_accepts_named_campaign_samples(tmp_path):
    catalog = tmp_path / "generic.yaml"
    catalog.write_text("schema_version: fcc_sample_catalog_v1\ncatalog_scope: generic\nsamples:\n  A: {}\n  B: {}\n")
    assert set(load_catalog(catalog)["samples"]) == {"A", "B"}
