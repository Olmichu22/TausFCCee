import unittest

from modules.fcc_truth_definitions import (
    descriptive_photon_origin,
    has_tau_origin,
    photon_origin_v1,
    reconstructed_pid_category,
    selected_truth_values,
)


class FrozenTruthDefinitionsTest(unittest.TestCase):
    def test_selected_truth(self):
        self.assertTrue(selected_truth_values(1, 22, (1.0, 0.0, 0.0)))
        self.assertFalse(selected_truth_values(0, 22, (1.0, 0.0, 0.0)))
        self.assertFalse(selected_truth_values(1, 14, (1.0, 0.0, 0.0)))
        self.assertFalse(selected_truth_values(1, 22, (0.0, 0.0, 0.0)))
        self.assertTrue(selected_truth_values(1, 22, (1.0e-10, 0.0, 0.0)))

    def test_pid_mapping_and_sentinel(self):
        expected = {
            -11: "electron", 13: "muon", 22: "photon", -211: "charged_pion",
            310: "K0S", 2112: "neutron", 3122: "Lambda",
        }
        self.assertEqual({pdg: reconstructed_pid_category(pdg) for pdg in expected}, expected)
        with self.assertRaisesRegex(ValueError, "unsupported"):
            reconstructed_pid_category(321)
        with self.assertRaisesRegex(ValueError, "sentinel"):
            reconstructed_pid_category(999)

    def test_tau_origin_uses_recursive_stored_parents(self):
        pdgs = [15, 111, 22, 22]
        parents = [[], [0], [1], []]
        self.assertTrue(has_tau_origin(pdgs, parents, 2))
        self.assertFalse(has_tau_origin(pdgs, parents, 3))

    def test_tau_origin_cycle_is_fatal(self):
        with self.assertRaisesRegex(ValueError, "cycle"):
            has_tau_origin([22, 15], [[1], [0]], 0)

    def test_photon_origin_v1_frozen_labels(self):
        desc = descriptive_photon_origin(
            tau_origin=True, pi0_ancestor=False, direct_tau_daughter=True,
            immediate_parent_pdgs=[15],
        )
        self.assertEqual(desc, "tau_direct_daughter")
        self.assertEqual(photon_origin_v1(desc, False), "tau_direct_daughter_unresolved")
        self.assertEqual(photon_origin_v1("parentless_non_tau", False), "parentless_unresolved")
        self.assertEqual(photon_origin_v1("electron_parent_non_tau", True), "explicit_ISR")
        self.assertEqual(photon_origin_v1("multiple_parent_non_tau", False), "resolved_non_tau_other")


if __name__ == "__main__":
    unittest.main()
