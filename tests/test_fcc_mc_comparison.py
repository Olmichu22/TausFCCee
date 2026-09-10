from collections import Counter
from pathlib import Path
import tempfile
import unittest

import numpy as np

from modules.fcc_mc_comparison import (
    classify_terminal_tau, hist_with_flow, load_comparison,
    strict_outcome, truthlink_representative,
)


ROOT = Path(__file__).resolve().parents[1]


class ComparisonContractTest(unittest.TestCase):
    def test_both_comparison_modes_share_one_contract(self):
        path = ROOT / "configs/analysis/fcc_mc_comparisons_v1.yaml"
        p8 = load_comparison(path, "whizard_p8o_18k")
        kk = load_comparison(path, "whizard_kkmcee_2k")
        stable = load_comparison(path, "p8o_p8h_stable10k")
        self.assertEqual(p8["contract"], kk["contract"])
        self.assertEqual([x["internal_name"] for x in p8["samples"]], ["W", "P8O"])
        self.assertEqual([x["internal_name"] for x in kk["samples"]], ["W", "KKMCee"])
        self.assertEqual(kk["samples"][0]["include_source_ids"], ["000242385"])
        self.assertEqual(kk["samples"][1]["expected_events"], 2000)
        self.assertEqual(stable["association_truth_scope"], "inclusive")
        self.assertEqual([x["internal_name"] for x in stable["samples"]], ["P8O", "P8H"])

    def test_representative_rule_is_track_then_cluster_without_index_tiebreak(self):
        rows = [
            {"pfo_index": 9, "track_permille": 400, "cluster_permille": 800},
            {"pfo_index": 2, "track_permille": 500, "cluster_permille": 10},
        ]
        self.assertEqual(truthlink_representative(rows)["pfo_index"], 2)
        tied = [dict(rows[1]), {"pfo_index": 1, "track_permille": 500, "cluster_permille": 10}]
        self.assertEqual(strict_outcome(tied), ("ambiguous_multiple_pfo", None))

    def test_cluster_only_rule_and_flow_normalization(self):
        rows = [{"pfo_index": 4, "track_permille": 0, "cluster_permille": 2},
                {"pfo_index": 1, "track_permille": 0, "cluster_permille": 5}]
        self.assertEqual(truthlink_representative(rows)["pfo_index"], 1)
        histogram = hist_with_flow([-2, .2, .8, 3], np.array([0., .5, 1.]), 4)
        self.assertEqual(sum(row["count"] for row in histogram), 4)
        self.assertAlmostEqual(sum(row["fraction"] for row in histogram), 1.)


if __name__ == "__main__":
    unittest.main()
