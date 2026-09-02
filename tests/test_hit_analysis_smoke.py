import unittest

from scripts.validation.validate_hit_analysis_smoke import stable_source_id


class HitAnalysisSmokeIdentityTest(unittest.TestCase):
    def test_manifest_ordinal_maps_to_stable_source_id(self):
        inputs = [
            {"source_file_id": "033851393", "path": "/data/a.root"},
            {"source_file_id": "090000001", "path": "/data/b.root"},
        ]
        self.assertEqual(stable_source_id(0, inputs), "033851393")
        self.assertEqual(stable_source_id(1, inputs), "090000001")

    def test_out_of_range_ordinal_is_rejected(self):
        with self.assertRaises(ValueError):
            stable_source_id(1, [{"source_file_id": "0", "path": "/data/a.root"}])


if __name__ == "__main__":
    unittest.main()
