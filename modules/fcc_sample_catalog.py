"""Configuration and manifest normalization for maintained FCC samples."""
from __future__ import annotations
import csv
import os
from pathlib import Path
import yaml

CATALOG_SCHEMA = "fcc_sample_catalog_v1"

def expand_path(value: str) -> Path:
    expanded = os.path.expandvars(value)
    if "$" in expanded:
        raise ValueError(f"unresolved environment variable in path: {value}")
    return Path(expanded)

def load_catalog(path: Path) -> dict:
    data = yaml.safe_load(Path(path).read_text())
    if data.get("schema_version") != CATALOG_SCHEMA:
        raise ValueError("unsupported sample catalog")
    required = {"W", "P8C", "P8O", "KKMCee"}
    if set(data.get("samples", {})) != required:
        raise ValueError("catalog must distinguish W, P8C, P8O, and KKMCee")
    return data

def resolve_records(sample: dict) -> list[dict]:
    if "single_source" in sample:
        return [{key: (str(expand_path(value)) if key != "source_id" else str(value)) for key, value in sample["single_source"].items()}]
    manifest = expand_path(sample["manifest"])
    with manifest.open(newline="") as handle:
        source = list(csv.DictReader(handle))
    columns = sample["columns"]
    records = [{target: row[column] for target, column in columns.items()} for row in source]
    include = set(map(str, sample.get("include_source_ids", [])))
    if include:
        records = [row for row in records if row["source_id"] in include]
        if {row["source_id"] for row in records} != include:
            raise ValueError("configured source IDs are not all present in manifest")
    if len({row["source_id"] for row in records}) != len(records):
        raise ValueError("duplicate source_file_id")
    return records
