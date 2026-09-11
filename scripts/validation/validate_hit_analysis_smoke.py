#!/usr/bin/env python3
"""Validate a tiny G run and its data-only interoperability with L products."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow.parquet as pq
import yaml

from modules.fcc_workflow_interface import load_product_manifest, read_association_table


def stable_source_id(raw_file_ordinal: int, inputs: list[dict]) -> str:
    if not 0 <= int(raw_file_ordinal) < len(inputs):
        raise ValueError(f"G source-file ordinal outside manifest: {raw_file_ordinal}")
    return str(inputs[int(raw_file_ordinal)]["source_file_id"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hit-manifest", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--workflow-manifest", type=Path, required=True)
    parser.add_argument("--golden-g-pfo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise FileExistsError(args.output)

    hit_manifest = yaml.safe_load(args.hit_manifest.read_text())
    inputs = hit_manifest["inputs"]
    g_table = pq.read_table(args.result_dir / "association_results_full_dR.parquet")
    g_rows = g_table.to_pylist()
    mapped = []
    for row in g_rows:
        source_id = stable_source_id(int(row["source_file_id"]), inputs)
        expected_name = Path(inputs[int(row["source_file_id"])]["path"]).name
        if row["source_file"] != expected_name:
            raise ValueError("G source filename does not match input manifest")
        mapped.append({**row, "stable_source_file_id": source_id})

    workflow = load_product_manifest(args.workflow_manifest)
    product = workflow["products"][0]
    direct = read_association_table(product["direct_assignment"], "direct").to_pylist()
    ancestor = read_association_table(product["ancestor_assignment"], "ancestor").to_pylist()
    direct_keys = [(str(row["source_file_id"]), int(row["event_in_file"]), int(row["pfo_index"])) for row in direct]
    ancestor_keys = [(str(row["source_file_id"]), int(row["event_in_file"]), int(row["pfo_index"])) for row in ancestor]
    matched = sorted(
        (int(row["reco"]), int(row["gen"]), int(row["Gen_pid"]))
        for row in mapped if int(row["reco"]) >= 0 and int(row["gen"]) >= 0
    )
    golden = pq.read_table(args.golden_g_pfo, filters=[("event_in_file", "=", 0)]).to_pylist()
    golden_matched = sorted(
        (int(row["pfo_index"]), int(row["G_mc_index"]), int(row["G_mc_pdg"]))
        for row in golden if row["G_status"] == "assigned"
    )
    stable_ids = {row["stable_source_file_id"] for row in mapped}
    event_keys = {(row["stable_source_file_id"], int(row["event_in_file"])) for row in mapped}
    reco_indices = [int(row["reco"]) for row in mapped if int(row["reco"]) >= 0]
    pfo_keys_compatible = {
        (stable_source_id, event, pfo) for stable_source_id, event, pfo in direct_keys
    } == {
        (row["stable_source_file_id"], int(row["event_in_file"]), int(row["reco"]))
        for row in mapped if int(row["reco"]) >= 0
    }
    result = {
        "schema_version": "fcc_tau_g_hit_analysis_smoke_v1",
        "status": "PASS",
        "events": len(event_keys),
        "g_rows": len(mapped),
        "selected_truth_rows": sum(int(row["gen"]) >= 0 for row in mapped),
        "matched_rows": len(matched),
        "unmatched_truth_rows": sum(int(row["gen"]) >= 0 and int(row["reco"]) < 0 for row in mapped),
        "fake_pfo_rows": sum(int(row["gen"]) < 0 and int(row["reco"]) >= 0 for row in mapped),
        "pfo_count": len(reco_indices),
        "raw_source_file_ordinals": sorted({int(row["source_file_id"]) for row in mapped}),
        "stable_source_file_ids": sorted(stable_ids),
        "dedup_mode": "reco",
        "duplicate_reco_indices": len(reco_indices) - len(set(reco_indices)),
        "direct_rows": len(direct),
        "ancestor_rows": len(ancestor),
        "duplicate_direct_keys": len(direct_keys) - len(set(direct_keys)),
        "duplicate_ancestor_keys": len(ancestor_keys) - len(set(ancestor_keys)),
        "direct_ancestor_keys_equal": set(direct_keys) == set(ancestor_keys),
        "g_l_pfo_keys_compatible_after_manifest_mapping": pfo_keys_compatible,
        "historical_g_pfo_mapping_equal": matched == golden_matched,
        "historical_g_pfo_mapping_rows": len(golden_matched),
        "g_full_regression_fixture_status": "G_GOLDEN_FULL_COUNT_FIXTURE_STILL_MISSING",
        "note": "The independent historical product validates event-0 PFO-to-MC mapping, not the complete raw G row set.",
    }
    passed = (
        result["events"] == 1 and result["g_rows"] > 0 and result["pfo_count"] == len(direct)
        and stable_ids == {str(product["source_file_id"])}
        and result["duplicate_reco_indices"] == 0
        and result["duplicate_direct_keys"] == 0 and result["duplicate_ancestor_keys"] == 0
        and result["direct_ancestor_keys_equal"] and pfo_keys_compatible
        and result["historical_g_pfo_mapping_equal"]
    )
    result["status"] = "PASS" if passed else "FAIL"
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
