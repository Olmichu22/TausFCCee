"""Data-only consumer for FCC-tau-workflow association products.

No Python code is imported from the workflow repository.  A versioned YAML or
JSON manifest supplies product paths and provenance; Parquet schemas are
validated before rows are exposed to analysis code.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import yaml


WORKFLOW_ASSOCIATION_CONTRACT = "fcc_tau_association_v1"
PRODUCT_MANIFEST_SCHEMA = "fcc_tau_workflow_product_manifest_v1"
EVENT_KEY_FIELDS = ("sample", "source_file_id", "event_in_file")
DIRECT_REQUIRED_COLUMNS = {
    "source_file_id", "event_in_file", "pfo_index", "truthlink_status",
    "assigned_mc_index",
}
ANCESTOR_REQUIRED_COLUMNS = {
    "sample", "source_file_id", "event_in_file", "pfo_index",
    "ancestor_status", "ancestor_mc_index",
}


def load_product_manifest(path: str | Path) -> dict:
    """Load and validate a workflow product manifest without opening products."""
    source = Path(path)
    if source.suffix.lower() == ".json":
        data = json.loads(source.read_text())
    else:
        data = yaml.safe_load(source.read_text())
    if not isinstance(data, dict):
        raise ValueError("workflow manifest must contain a mapping")
    if data.get("schema_version") != PRODUCT_MANIFEST_SCHEMA:
        raise ValueError("unsupported workflow product manifest schema")
    if data.get("association_contract") != WORKFLOW_ASSOCIATION_CONTRACT:
        raise ValueError("unsupported workflow association contract")
    if not isinstance(data.get("sample"), str) or not data["sample"]:
        raise ValueError("workflow manifest requires a non-empty sample")
    products = data.get("products")
    if not isinstance(products, list) or not products:
        raise ValueError("workflow manifest requires at least one product entry")
    required = {
        "source_file_id", "source_rec", "direct_assignment",
        "ancestor_assignment", "truth_definition_version", "input_provenance",
    }
    seen: set[str] = set()
    for entry in products:
        if not isinstance(entry, dict) or required - set(entry):
            raise ValueError(f"incomplete workflow product entry: {entry!r}")
        source_file_id = str(entry["source_file_id"])
        if source_file_id in seen:
            raise ValueError(f"duplicate source_file_id: {source_file_id}")
        seen.add(source_file_id)
    return data


def authoritative_event_key(sample: str, source_file_id: str, event_in_file: int) -> tuple[str, str, int]:
    """Return the cross-repository composite event identity."""
    if not sample or not str(source_file_id) or int(event_in_file) < 0:
        raise ValueError("invalid event identity")
    return str(sample), str(source_file_id), int(event_in_file)


def validate_association_columns(columns: Iterable[str], level: str) -> None:
    required = DIRECT_REQUIRED_COLUMNS if level == "direct" else ANCESTOR_REQUIRED_COLUMNS if level == "ancestor" else None
    if required is None:
        raise ValueError(f"unknown association level: {level}")
    missing = required - set(columns)
    if missing:
        raise ValueError(f"{level} association product missing columns: {sorted(missing)}")


def read_association_table(path: str | Path, level: str):
    """Read a Parquet product after checking the Stage-2 physical schema."""
    import pyarrow.parquet as pq

    table = pq.read_table(path)
    validate_association_columns(table.column_names, level)
    return table


def truthlink_representative(rows: list[dict]) -> dict:
    """Invert exported PFO assignments using the frozen representative rule.

    This consumes assignment rows; it does not reimplement L_direct or
    L_ancestor. Exact terminal ties remain ambiguous and no PFO-index tie-break
    is introduced.
    """
    if not rows:
        return {"status": "unmatched", "pfo_index": None, "multiplicity": 0}
    if len(rows) == 1:
        return {"status": "assigned", "pfo_index": rows[0]["pfo_index"], "multiplicity": 1}
    track = [row for row in rows if int(row["track_permille"]) > 0]
    if track:
        best_t = max(int(row["track_permille"]) for row in track)
        stage = [row for row in track if int(row["track_permille"]) == best_t]
        best_c = max(int(row["cluster_permille"]) for row in stage)
        winners = [row for row in stage if int(row["cluster_permille"]) == best_c]
    else:
        best_c = max(int(row["cluster_permille"]) for row in rows)
        winners = [row for row in rows if int(row["cluster_permille"]) == best_c]
    if len(winners) != 1:
        return {
            "status": "ambiguous_multiple_pfo", "pfo_index": None,
            "multiplicity": len(rows), "tied_pfo_count": len(winners),
        }
    return {"status": "assigned", "pfo_index": winners[0]["pfo_index"], "multiplicity": len(rows)}
