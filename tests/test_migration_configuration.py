import tempfile
import unittest
from pathlib import Path

import yaml

from scripts.analysis import produce_migration_matrices as migration


class MigrationConfigurationTest(unittest.TestCase):
    def test_path_fields_remain_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = {
                "schema_version": "fcc_migration_matrix_inputs_v1",
                "output_root": str(root / "out"),
                "report_path": str(root / "report.md"),
                "g_products": {sample: str(root / f"{sample}.parquet") for sample in ("W", "P8C", "P8O")},
                "workflow_manifests": {sample: str(root / f"{sample}.csv") for sample in ("W", "P8C", "P8O")},
                "validation_inputs": {
                    "fiducial_tables": str(root / "fiducial"),
                    "lancestor_tables": str(root / "lancestor"),
                },
                "provenance": {
                    "historical_definition_source": str(root / "definition.py"),
                    "historical_comparison_reference": str(root / "comparison"),
                },
            }
            path = root / "config.yaml"
            path.write_text(yaml.safe_dump(config))
            migration.configure(path)
        self.assertIsInstance(migration.HIST_SCRIPT, Path)
        self.assertIsInstance(migration.HIST_COMPARISON, Path)
        self.assertEqual(set(migration.G_PATHS), {"W", "P8C", "P8O"})


if __name__ == "__main__":
    unittest.main()
