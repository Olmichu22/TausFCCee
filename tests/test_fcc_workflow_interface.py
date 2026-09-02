import tempfile
from pathlib import Path
import unittest

import yaml

from modules.fcc_workflow_interface import (
    authoritative_event_key,
    load_product_manifest,
    truthlink_representative,
    validate_association_columns,
)


class WorkflowInterfaceTest(unittest.TestCase):
    def manifest(self):
        return {
            "schema_version": "fcc_tau_workflow_product_manifest_v1",
            "association_contract": "fcc_tau_association_v1",
            "sample": "W",
            "products": [{
                "source_file_id": "033851393", "source_rec": "/data/a.root",
                "direct_assignment": "/data/a_direct.parquet",
                "ancestor_assignment": "/data/a_ancestor.parquet",
                "truth_definition_version": "selected_truth_v1",
                "input_provenance": "/data/a.json",
            }],
        }

    def test_manifest_schema_and_composite_event_key(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "manifest.yaml"
            path.write_text(yaml.safe_dump(self.manifest()))
            self.assertEqual(load_product_manifest(path)["sample"], "W")
        self.assertEqual(authoritative_event_key("W", "033851393", 1999), ("W", "033851393", 1999))

    def test_duplicate_source_id_rejected(self):
        data = self.manifest(); data["products"].append(dict(data["products"][0]))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "manifest.yaml"; path.write_text(yaml.safe_dump(data))
            with self.assertRaisesRegex(ValueError, "duplicate"):
                load_product_manifest(path)

    def test_physical_product_schema(self):
        validate_association_columns(
            ["source_file_id", "event_in_file", "pfo_index", "truthlink_status", "assigned_mc_index"],
            "direct",
        )
        with self.assertRaisesRegex(ValueError, "missing columns"):
            validate_association_columns(["source_file_id"], "ancestor")

    def test_representative_rule_has_no_index_tiebreak(self):
        rows = [
            {"pfo_index": 8, "track_permille": 0, "cluster_permille": 10},
            {"pfo_index": 2, "track_permille": 0, "cluster_permille": 10},
        ]
        self.assertEqual(truthlink_representative(rows)["status"], "ambiguous_multiple_pfo")
        rows[1]["track_permille"] = 1
        self.assertEqual(truthlink_representative(rows)["pfo_index"], 2)


if __name__ == "__main__":
    unittest.main()
