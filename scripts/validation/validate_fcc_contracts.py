#!/usr/bin/env python3
"""Static consistency checks for maintained FCC analysis contracts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from modules.fcc_truth_definitions import (  # noqa: E402
    MIN_SELECTED_MOMENTUM_GEV,
    NEUTRINO_ABS_PDGS,
    SELECTED_TRUTH_VERSION,
)
from modules.fcc_workflow_interface import (  # noqa: E402
    PRODUCT_MANIFEST_SCHEMA,
    WORKFLOW_ASSOCIATION_CONTRACT,
    load_product_manifest,
    read_association_table,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--product-manifest", type=Path)
    args = parser.parse_args()
    g = yaml.safe_load((REPO / "configs/analysis/geometric_association_v1.yaml").read_text())
    assert g["schema_version"] == "geometric_association_v1"
    assert g["distance"]["threshold"] == 0.1
    assert g["distance"]["comparison"] == "strict_less_than"
    assert g["truth_selection"]["version"] == SELECTED_TRUTH_VERSION
    assert set(g["truth_selection"]["excluded_abs_pdg"]) == set(NEUTRINO_ABS_PDGS)
    assert float(g["truth_selection"]["minimum_momentum_GeV"]) == MIN_SELECTED_MOMENTUM_GEV
    interface = yaml.safe_load((REPO / "configs/interfaces/fcc-tau-workflow-v1.yaml").read_text())
    assert interface["accepted_association_contract"] == WORKFLOW_ASSOCIATION_CONTRACT
    assert interface["accepted_product_manifest_schema"] == PRODUCT_MANIFEST_SCHEMA
    load_product_manifest(REPO / "configs/interfaces/example_workflow_product_manifest.yaml")
    if args.product_manifest:
        product_manifest = load_product_manifest(args.product_manifest)
        for product in product_manifest["products"]:
            if not Path(product["source_rec"]).is_file():
                raise FileNotFoundError(product["source_rec"])
            direct = read_association_table(product["direct_assignment"], "direct")
            ancestor = read_association_table(product["ancestor_assignment"], "ancestor")
            source_id = str(product["source_file_id"])
            if set(map(str, direct.column("source_file_id").to_pylist())) != {source_id}:
                raise ValueError("direct source_file_id mismatch")
            if set(map(str, ancestor.column("source_file_id").to_pylist())) != {source_id}:
                raise ValueError("ancestor source_file_id mismatch")
    freeze = REPO / "docs/scientific/freezes/20260901"
    for path in freeze.glob("*.json"):
        json.loads(path.read_text())
    print("fcc_contract_validation=PASS")


if __name__ == "__main__":
    main()
