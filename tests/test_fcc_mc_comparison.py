from collections import Counter
from pathlib import Path
import tempfile
import unittest

import numpy as np

from modules.fcc_mc_comparison import (
    PID_CATEGORIES, SampleData, angle3d_mrad, binned_efficiency, classify_terminal_tau,
    clamped_acos_mrad, fiducial_outcomes, fiducial_provenance, hist_with_flow,
    load_comparison, nominal_truth_species, normalize_performance_config,
    performance_residuals, selected_truth_inventory, strict_outcome,
    truth_fiducial_accepts, truthlink_representative, wrap_delta_phi,
)
from modules.fcc_mc_comparison_outputs import (
    _fiducial_efficiency_tables, _nominal_performance_tables, _nominal_pid_tables,
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
        self.assertEqual(p8["performance"]["association_method"], "Lancestor")
        self.assertIsNone(p8["performance"]["truth_p_min"])
        self.assertIsNone(p8["performance"]["truth_theta_min_deg"])
        self.assertIsNone(p8["performance"]["truth_theta_max_deg"])

    def test_representative_rule_is_track_then_cluster_without_index_tiebreak(self):
        rows = [
            {"pfo_index": 9, "track_permille": 400, "cluster_permille": 800},
            {"pfo_index": 2, "track_permille": 500, "cluster_permille": 10},
        ]
        self.assertEqual(truthlink_representative(rows)["pfo_index"], 2)
        tied = [dict(rows[1]), {"pfo_index": 1, "track_permille": 500, "cluster_permille": 10}]
        self.assertEqual(strict_outcome(tied), ("ambiguous_multiple_pfo", None))


    def test_nominal_species_policy_and_other_inventory(self):
        self.assertEqual(nominal_truth_species(321), "charged_kaon")
        self.assertEqual(nominal_truth_species(-321), "charged_kaon")
        self.assertEqual(nominal_truth_species(130), "K0L")
        self.assertIsNone(nominal_truth_species(111))
        self.assertIsNone(nominal_truth_species(310))
        rows = [{"truth_pdg": pdg} for pdg in (11, -11, 321, -321, 130, 111, 310, 2112, -2112)]
        inventory = selected_truth_inventory(rows)
        lookup = {(row["category"], row["species"], row["pdg"]): row["count"] for row in inventory}
        self.assertEqual(lookup[("nominal_species", "electron", "")], 2)
        self.assertEqual(lookup[("nominal_species", "charged_kaon", "")], 2)
        self.assertEqual(lookup[("nominal_species", "K0L", "")], 1)
        for pdg in (-2112, 111, 310, 2112):
            self.assertEqual(lookup[("other_selected_truth", "other_selected_truth", pdg)], 1)
        self.assertEqual(sum(row["count"] for row in inventory), len(rows))

    def test_wrapped_and_ordinary_signed_delta_phi(self):
        epsilon = 1e-6
        self.assertAlmostEqual(wrap_delta_phi(-np.pi + epsilon, np.pi - epsilon), 2 * epsilon)
        self.assertAlmostEqual(wrap_delta_phi(0.2, 0.5), -0.3)
        self.assertEqual(wrap_delta_phi(-np.pi, 0.0), np.pi)

    def test_exact_opening_angle_and_clamp(self):
        self.assertAlmostEqual(angle3d_mrad((1, 0, 0), (0, 1, 0)), np.pi / 2 * 1000)
        self.assertEqual(angle3d_mrad((1, 2, 3), (1, 2, 3)), 0.0)
        self.assertEqual(clamped_acos_mrad(1.0 + 1e-12), 0.0)
        self.assertAlmostEqual(clamped_acos_mrad(-1.0 - 1e-12), np.pi * 1000)

    def test_additive_residuals_preserve_dp_and_dtheta_and_add_energy(self):
        pair = {
            "truth_p": 2.0, "reco_p": 3.0,
            "truth_theta": 10.0, "reco_theta": 11.0,
            "truth_phi": np.pi - 1e-6, "reco_phi": -np.pi + 1e-6,
            "truth_energy": 4.0, "reco_energy": 5.0,
            "truth_px": 1.0, "truth_py": 0.0, "truth_pz": 0.0,
            "reco_px": 0.0, "reco_py": 1.0, "reco_pz": 0.0,
        }
        residuals = performance_residuals(pair)
        self.assertEqual(residuals["dp_over_p"], 0.5)
        self.assertAlmostEqual(residuals["dtheta_mrad"], np.pi / 180 * 1000)
        self.assertAlmostEqual(residuals["dphi_mrad"], 0.002)
        self.assertAlmostEqual(residuals["angle3d_mrad"], np.pi / 2 * 1000)
        self.assertEqual(residuals["de_over_e"], 0.25)

    def test_nominal_performance_summary_distinguishes_signed_and_unsigned(self):
        pairs = []
        for signed, angle in ((-1.0, 1.0), (1.0, 3.0)):
            pairs.append({
                "association_method": "Lancestor", "truth_species": "photon",
                "dp_over_p": signed, "dtheta_mrad": signed,
                "dphi_mrad": signed, "de_over_e": signed,
                "angle3d_mrad": angle,
            })
        data = SampleData(
            "X", "X", 1, [{"truth_pdg": 22}], [], pairs, [], Counter(), Counter(), [],
            {"adapter": "workflow_products"},
        )
        inventory, summaries = _nominal_performance_tables([data])
        self.assertEqual(next(row for row in inventory if row["species"] == "photon")["count"], 1)
        signed = next(row for row in summaries
                      if row["species"] == "photon" and row["residual"] == "dp_over_p")
        unsigned = next(row for row in summaries
                        if row["species"] == "photon" and row["residual"] == "angle3d_mrad")
        self.assertAlmostEqual(signed["q16"], -0.68)
        self.assertAlmostEqual(signed["q84"], 0.68)
        self.assertAlmostEqual(signed["central68_halfwidth"], 0.68)
        self.assertEqual(unsigned["median"], 2.0)
        self.assertAlmostEqual(unsigned["q68"], 2.36)
        self.assertAlmostEqual(unsigned["q95"], 2.9)
        self.assertEqual(unsigned["q16"], "")

    @staticmethod
    def _truth_outcome_row(p, theta, phi, outcome="association_unmatched", **reco):
        return {"p": p, "theta": theta, "phi": phi, "outcome": outcome, **reco}

    def test_default_fiducial_is_identity_and_uses_truth_only(self):
        rows = [
            self._truth_outcome_row(0.1, 0.0, -1.0, reco_p=1000.0, reco_theta=90.0),
            self._truth_outcome_row(2.0, 180.0, 1.0, reco_p=0.0, reco_theta=-999.0),
        ]
        self.assertEqual(fiducial_outcomes(rows), rows)
        self.assertTrue(truth_fiducial_accepts(rows[0], None))
        self.assertEqual(normalize_performance_config()["truth_p_min"], None)

    def test_fiducial_p_and_theta_cuts_are_inclusive(self):
        rows = [
            self._truth_outcome_row(0.9, 9.9, 0.0, reco_p=100.0, reco_theta=90.0),
            self._truth_outcome_row(1.0, 10.0, 0.0, reco_p=0.0, reco_theta=-10.0),
            self._truth_outcome_row(2.0, 170.0, 0.0, reco_p=0.0, reco_theta=999.0),
            self._truth_outcome_row(3.0, 170.1, 0.0, reco_p=0.0, reco_theta=90.0),
        ]
        config = {"truth_p_min": 1.0, "truth_theta_min_deg": 10.0,
                  "truth_theta_max_deg": 170.0}
        self.assertEqual(fiducial_outcomes(rows, config), rows[1:3])

    def test_binned_efficiency_counts_subset_and_zero_denominator(self):
        rows = [
            self._truth_outcome_row(0.5, 20.0, 0.0, "associated_unique"),
            self._truth_outcome_row(0.6, 20.0, 0.0, "ambiguous_multiple_pfo"),
            self._truth_outcome_row(1.5, 20.0, 0.0, "association_unmatched"),
        ]
        binned = binned_efficiency(rows, "p", (0.0, 1.0, 2.0, 3.0))
        regular = [row for row in binned if row["bin_kind"] == "regular"]
        self.assertEqual((regular[0]["N_truth"], regular[0]["N_associated_unique"]), (2, 1))
        self.assertEqual((regular[1]["N_truth"], regular[1]["N_associated_unique"]), (1, 0))
        self.assertEqual(regular[2]["N_truth"], 0)
        self.assertEqual(regular[2]["efficiency_percent"], "")
        self.assertLessEqual(sum(row["N_associated_unique"] for row in binned),
                             sum(row["N_truth"] for row in binned))

    def test_binned_theta_and_periodic_phi_account_once(self):
        theta_rows = [
            self._truth_outcome_row(1.0, 10.0, 0.0, "associated_unique"),
            self._truth_outcome_row(1.0, 20.0, 0.0, "association_unmatched"),
        ]
        theta = [row for row in binned_efficiency(theta_rows, "theta", (0.0, 15.0, 30.0))
                 if row["bin_kind"] == "regular"]
        self.assertEqual([(row["N_truth"], row["N_associated_unique"]) for row in theta],
                         [(1, 1), (1, 0)])
        epsilon = 1e-9
        phi_rows = [self._truth_outcome_row(1.0, 90.0, value, "associated_unique")
                    for value in (-np.pi, -np.pi + epsilon, np.pi - epsilon, np.pi)]
        phi = binned_efficiency(phi_rows, "phi", (-np.pi, 0.0, np.pi))
        regular = [row for row in phi if row["bin_kind"] == "regular"]
        self.assertEqual([row["N_truth"] for row in regular], [1, 3])
        self.assertEqual(sum(row["N_truth"] for row in phi), 4)
        self.assertEqual(phi[0]["N_truth"], 0)
        self.assertEqual(phi[-1]["N_truth"], 0)

    def test_fiducial_outputs_serialize_provenance_and_ambiguity_is_failure(self):
        outcomes = [
            {**self._truth_outcome_row(1.0, 10.0, 0.0, "associated_unique"),
             "truth_definition": "Lancestor", "truth_species": "electron", "tau_ancestor": True},
            {**self._truth_outcome_row(2.0, 20.0, 0.1, "ambiguous_multiple_pfo"),
             "truth_definition": "Lancestor", "truth_species": "electron", "tau_ancestor": True},
        ]
        data = SampleData("X", "X", 1, [], outcomes, [], [], Counter(), Counter(), [],
                          {"adapter": "workflow_products"})
        config = {"truth_p_min": 1.0, "truth_theta_min_deg": 10.0,
                  "truth_theta_max_deg": 170.0, "association_method": "Lancestor",
                  "efficiency_bins": {"p": (0.0, 3.0), "theta": (0.0, 180.0),
                                      "phi": (-np.pi, 0.0, np.pi)}}
        integrated, differential = _fiducial_efficiency_tables(data, True, config)
        electron = next(row for row in integrated if row["species"] == "electron")
        self.assertEqual((electron["denominator_count"], electron["numerator_count"]), (2, 1))
        self.assertEqual(electron["truth_p_min"], 1.0)
        self.assertEqual(electron["truth_theta_min_deg"], 10.0)
        self.assertEqual(electron["truth_theta_max_deg"], 170.0)
        self.assertEqual(sum(row["numerator_count"] for row in differential["p"]
                             if row["species"] == "electron"), 1)
        self.assertEqual(fiducial_provenance({}), {
            "truth_p_min": "none", "truth_theta_min_deg": "none",
            "truth_theta_max_deg": "none",
        })

    def test_charged_kaon_truth_is_excluded_from_pid_evaluation(self):
        pair = {"association_method": "Lancestor", "truth_species": "charged_kaon",
                "tau_ancestor": False, "p": 2.0, "theta": 90.0,
                "reco_category": "charged_pion"}
        data = SampleData("X", "X", 1, [], [], [pair], [], Counter(), Counter(), [],
                          {"adapter": "workflow_products"})
        confusion, summary = _nominal_pid_tables(data, None, {})
        self.assertFalse(any(row["truth_species"] == "charged_kaon" for row in confusion))
        self.assertFalse(any(row["truth_species"] == "charged_kaon" for row in summary))
        self.assertNotIn("charged_kaon", PID_CATEGORIES)

    def test_cluster_only_rule_and_flow_normalization(self):
        rows = [{"pfo_index": 4, "track_permille": 0, "cluster_permille": 2},
                {"pfo_index": 1, "track_permille": 0, "cluster_permille": 5}]
        self.assertEqual(truthlink_representative(rows)["pfo_index"], 1)
        histogram = hist_with_flow([-2, .2, .8, 3], np.array([0., .5, 1.]), 4)
        self.assertEqual(sum(row["count"] for row in histogram), 4)
        self.assertAlmostEqual(sum(row["fraction"] for row in histogram), 1.)


if __name__ == "__main__":
    unittest.main()
