# FCC-ee tau photon/pion reconstruction: consolidated W–P8O report

Date: 2026-09-01  
Status: **SCIENTIFIC ANALYSIS FROZEN AT 2026-09-01**  
Scientific status: **PASS_WITH_WARNINGS**

## Executive summary

This note freezes the validated FCC-ee tau photon/pion diagnostic programme. Its central result is that the large global W/P8O photon difference is driven primarily by sample/generator truth composition and truth-to-PFO association/reconstructability, not by a large change in conditional reconstructed photon PID. Tau topology, visible tau-decay truth, tau-decay photon kinematics, reconstructed PFO composition, and conditional PID are closely compatible between W and P8O.

The decisive factorization is

```text
origin / truth composition
  -> phase space
  -> unique truth-to-PFO association
  -> conditional reconstructed PID
  -> truth-normalized global migration outcome
```

P8O contains 49,120 photons with positive stored genealogy evidence for `explicit_ISR`; this population is exceptionally soft and beam-collinear. W contains 25,879 parentless non-tau photons with no stored photon-to-beam ancestry. They remain `parentless_unresolved` and may only be described as having partially similar or ISR-like kinematics. Direct tau daughters remain `tau_direct_daughter_unresolved`: the stored records do not support a positive-evidence, cross-generator FSR classification.

## 1. Samples, event provenance, and frozen scope

| Label | Generated sample | Reconstruction | Events relevant here |
|---|---|---|---:|
| W | WHIZARD | our reconstruction | 1,992,000 available; first deterministic 18,000 used in equal-statistics audits |
| P8C | PYTHIA8 | collaborator reconstruction | migration-matrix baseline only |
| P8O | PYTHIA8 | our reconstruction | 18,000 available and used |

W versus P8O isolates different generated samples under the same reconstruction chain. The detailed audits use the same deterministic event manifest logic and 18,000 events per sample. The manifest SHA-256 is `b2f81449579d33ec228b65c20a7927f8930ef19368aedd97244085c6e587cf3b`.

The sole inherited topology warning is P8O event `out_sim_edm4hep_1:838`, which has four terminal tau representatives. It is retained. All other 17,999 P8O events and all 18,000 W events have two terminal representatives.

## 2. Frozen scientific definitions

Selected truth is frozen as already implemented: `generatorStatus == 1`, neutrinos excluded, non-zero truth momentum, and the established species mapping. The reconstructed PID map is `|PDG|=11` electron, `13` muon, `22` photon, `211` charged pion, `310` K0S, `2112` neutron, and `3122` Lambda.

Association definitions are complementary diagnostics:

- **G:** geometric selected-truth association.
- **L_direct:** immediate detector-level MC contributor, using the frozen representative-PFO inversion.
- **L_ancestor:** nearest unique selected generator-level ancestor, with the same inversion.

`tau_origin=True` exactly when recursively stored ancestry contains `abs(PDG)==15`. G/L_direct/L_ancestor differences are **association-definition dependence**, not changes in the detector or detector efficiency.

The frozen `photon_origin_v1` categories are `tau_decay_pi0`, `tau_direct_daughter_unresolved`, `tau_other_descendant`, `explicit_ISR`, `parentless_unresolved`, `electron_parent_non_tau_unresolved`, `positron_parent_non_tau_unresolved`, `photon_chain_non_tau_unresolved`, and `resolved_non_tau_other`.

## 3. Why global migration matrices are not sufficient

A truth-row-normalized migration matrix combines three physically distinct effects: the composition and phase space of selected truth particles, whether a unique representative PFO exists, and the PFO PID conditional on that association. An unmatched truth object is not a PID failure. Consequently, the raw photon diagonal cannot by itself support a claim that one reconstruction reconstructs photons better globally.

### Global migration headline

| Association | Sample | photon→photon | photon→unmatched |
|---|---|---:|---:|
| G | W | 51.249% | 46.610% |
| G | P8C | 37.053% | 61.109% |
| G | P8O | 36.846% | 61.379% |
| L_ancestor | W | 51.149% | 45.068% |
| L_ancestor | P8C | 39.526% | 57.496% |
| L_ancestor | P8O | 39.355% | 57.786% |

The dominant difference is diagonal ↔ unmatched, not transfer to a particular wrong PID. Charged-pion diagonals remain high and comparatively stable.

![W G inclusive migration matrix](../artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/inclusive/migration_matrix_W_G_inclusive.png)

![P8O G inclusive migration matrix](../artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/inclusive/migration_matrix_P8O_G_inclusive.png)

## 4. Tau topology and visible tau-decay truth

The tau-pair topology and principal decay fractions agree closely.

| Tau decay category | W | P8O |
|---|---:|---:|
| 1-prong with pi0 | 36.797% | 36.892% |
| electron | 18.336% | 17.688% |
| muon | 17.583% | 17.493% |
| 1-prong no pi0 | 12.256% | 12.144% |
| 3-prong no pi0 | 9.744% | 9.991% |
| 3-prong with pi0 | 5.153% | 5.664% |

Representative tau means are `p = 44.960/45.007 GeV`, `pT = 33.142/33.279 GeV`, and `theta = 89.961/89.983 deg` for W/P8O. The stored-tau-copy multiplicity diagnostic is deliberately excluded from the scientific figure set.

![Tau decay channels](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/01_tau_decay_channels.png)

![Tau momentum](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/03_tau_p.png)

![Tau transverse momentum](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/04_tau_pt.png)

![Tau angular distributions](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/05_tau_theta.png)

Both tau decay chains are included and neutrinos are excluded. Selected visible truth particles from tau decays occur at 5.350/event W and 5.334/event P8O. Their principal composition is:

| Species | W | P8O |
|---|---:|---:|
| photon | 50.329% | 49.478% |
| charged pion | 34.701% | 35.444% |
| electron | 7.237% | 7.242% |
| muon | 6.593% | 6.575% |

Tau-origin photons constitute 65.190% of selected photons in the W 18k subset (65.481% in the full W reference) and 49.163% in P8O.

![Tau-decay truth composition](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/07_tau_decay_truth_composition.png)

## 5. The non-tau photon population

| Quantity | W non-tau photons | P8O non-tau photons |
|---|---:|---:|
| count | 25,879 | 49,121 |
| photons/event | 1.438 | 2.729 |
| P8O/W count ratio | — | 1.898 |
| mean E | 0.217 GeV | 0.100 GeV |
| median E | 3.17e-4 GeV | 1.41e-6 GeV |
| median pT | 2.86e-5 GeV | 9.88e-9 GeV |
| `|cos(theta)| > 0.98` | 50.191% | 83.748% |
| `|cos(theta)| > 0.999` | 31.62% | 78.57% |

The figures retain physical E and pT in GeV, logarithmic x-axes, logarithmically spaced bins, identical W/P8O bins and ranges, and explicit `/ log bin` normalization. No variable is transformed to `log10(X)`.

![Non-tau photon energy](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/16_non_tau_photon_energy_logx.png)

![Non-tau photon transverse momentum](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/17_non_tau_photon_pt_logx.png)

![Non-tau photon angle](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/19_non_tau_photon_costheta.png)

## 6. Reconstructed PFO composition

The large truth-level difference is strongly attenuated in reconstructed PFO content.

| Reconstructed PID | W / event | P8O / event |
|---|---:|---:|
| electron | 0.5403 | 0.5333 |
| muon | 0.3265 | 0.3228 |
| photon | 2.8321 | 2.8820 |
| charged pion | 1.8621 | 1.8847 |
| K0S | 0.0152 | 0.0148 |
| neutron | 0.4409 | 0.4573 |
| Lambda | 0.0182 | 0.0193 |

Despite 1.898 times more non-tau truth photons in P8O, photon-PFO yield differs by only 1.76%. Charged/neutral PFO populations and the main reconstructed PID populations are very similar.

![PFO PID composition](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/22_pfo_pid_composition.png)

![PFO PID yield per event](../artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/plots/primary/23_pfo_pid_yield.png)

## 7. Conditional PID audit

Conditional PID is `P(correct reconstructed PID | unique associated PFO exists)`. Unmatched and ambiguous particles are excluded from this denominator and are not PID failures.

| Truth species | W associated | W accuracy | P8O associated | P8O accuracy |
|---|---:|---:|---:|---:|
| electron | 6,725 | 89.844% | 6,662 | 89.433% |
| muon | 6,254 | 88.647% | 6,185 | 88.747% |
| photon | 40,478 | 93.438% | 40,627 | 93.600% |
| charged pion | 32,447 | 92.267% | 32,908 | 91.965% |

Mis-ID populations are often lower-pT and, particularly for pions, more forward/backward. The W/P8O conditional accuracies agree within 0.42 percentage points.

![W conditional PID matrix](../artifacts/pid_audit_W_P8O_equalstats_20260901/primary/01_conditional_pid_matrix_W_Lancestor.png)

![P8O conditional PID matrix](../artifacts/pid_audit_W_P8O_equalstats_20260901/primary/02_conditional_pid_matrix_P8O_Lancestor.png)

![Conditional correct and mis-ID fractions](../artifacts/pid_audit_W_P8O_equalstats_20260901/primary/03_pid_correct_misid_fractions_Lancestor.png)

![Pion pT correct versus mis-ID](../artifacts/pid_audit_W_P8O_equalstats_20260901/primary/05_truth_pion_pt_correct_vs_misid.png)

![Photon kinematics correct versus mis-ID](../artifacts/pid_audit_W_P8O_equalstats_20260901/primary/09_truth_photon_energy_pt_correct_vs_misid.png)

## 8. Truth-to-PFO association/reconstruction

Primary efficiency is `epsilon_unique = N(unique associated PFO) / N(selected truth)` under L_ancestor.

| Species | W selected | W epsilon_unique | P8O selected | P8O epsilon_unique |
|---|---:|---:|---:|---:|
| electron | 6,969 | 96.499% | 6,959 | 95.732% |
| muon | 6,349 | 98.504% | 6,313 | 97.972% |
| photon | 74,344 | 54.447% | 96,625 | 42.046% |
| charged pion | 33,416 | 97.100% | 34,030 | 96.703% |

Electrons, muons, and charged pions usually obtain a unique PFO. Photons are exceptional.

![PFO association efficiency by species](../artifacts/pfo_reconstruction_efficiency_W_P8O_equalstats_20260901/primary/02_global_epsilon_by_species.png)

### Tau-origin versus non-tau photons

| Photon population | W epsilon_unique | P8O epsilon_unique |
|---|---:|---:|
| tau-origin | 81.013% | 81.873% |
| non-tau | 4.695% | 3.530% |
| non-tau, E < 0.1 GeV | 0.313% | 0.230% |
| non-tau, pT < 0.1 GeV | 1.632% | 1.243% |
| non-tau, `|cos(theta)| > 0.98` | 2.048% | 0.887% |

The P8O non-tau excess predominantly occupies extremely low-efficiency phase space.

![Tau-photon association efficiency](../artifacts/pfo_reconstruction_efficiency_W_P8O_equalstats_20260901/primary/09_tau_photon_efficiency_energy.png)

![Non-tau photon energy efficiency](../artifacts/pfo_reconstruction_efficiency_W_P8O_equalstats_20260901/primary/14_non_tau_photon_efficiency_energy.png)

![Non-tau photon angular efficiency](../artifacts/pfo_reconstruction_efficiency_W_P8O_equalstats_20260901/primary/17_non_tau_photon_efficiency_costheta.png)

Charged pions form the control: association is about 97%, conditional PID about 92%, and W/P8O differences are small. Their problematic region is low pT and strongly forward/backward phase space.

## 9. Photon-origin discovery

The event records encode the two dominant non-tau populations differently. In P8O, 49,120 selected photons satisfy a positive-evidence rule tracing same-sign electron-copy chains to parentless incoming e-/e+ objects with approximately `E=45.595 GeV` and `|cos(theta)|=0.9999755`. They are classified `explicit_ISR`.

In W, 25,879 non-tau photons are parentless in the stored record. Beam e-/e+ objects exist, but no photon ancestry relationship is stored. These photons remain `parentless_unresolved`; their kinematics may be described as partially similar or ISR-like, never as proof of physical equivalence.

| photon_origin_v1 | W count | P8O count |
|---|---:|---:|
| parentless_unresolved | 25,879 | 0 |
| explicit_ISR | 0 | 49,120 |
| electron_parent_non_tau_unresolved | 0 | 1 |
| tau_decay_pi0 | 39,585 | 39,727 |
| tau_direct_daughter_unresolved | 8,640 | 7,697 |
| tau_other_descendant | 240 | 80 |

![Photon-origin composition](../artifacts/photon_origin_discovery_W_P8O_equalstats_20260901/plots/01_genealogy_class_composition.png)

![Initial electron-chain evidence](../artifacts/photon_origin_discovery_W_P8O_equalstats_20260901/plots/13_initial_electron_chain_summary.png)

Direct tau daughters usually appear beside continuing tau copies: 8,640 W and 7,697 P8O. The topology mixes radiation/bookkeeping and terminal decay configurations and does not supply a positive-evidence cross-generator FSR definition. The frozen category is therefore `tau_direct_daughter_unresolved`, not FSR.

## 10. Reconstruction and PID by photon origin

| Origin | W/P8O association | W/P8O conditional photon PID | W/P8O full photon outcome |
|---|---:|---:|---:|
| parentless / explicit ISR | 4.695 / 3.530% | 94.239 / 94.348% | 4.424 / 3.330% |
| tau_decay_pi0 | 91.350 / 90.757% | 93.393 / 93.491% | 85.316 / 84.845% |
| tau_direct_daughter_unresolved | 33.854 / 36.001% | 93.573 / 94.551% | 31.678 / 34.040% |

Once a unique PFO exists, photon PID is stable at approximately 93–95% across samples and origins. The large variation occurs before PID, in truth population, phase space, and association.

![Association efficiency by photon origin](../artifacts/photon_origin_reco_audit_W_P8O_equalstats_20260901/primary/07_association_efficiency_by_origin.png)

![Conditional photon PID by origin](../artifacts/photon_origin_reco_audit_W_P8O_equalstats_20260901/primary/13_conditional_photon_pid_accuracy.png)

![Full origin-to-association-to-PID chain](../artifacts/photon_origin_reco_audit_W_P8O_equalstats_20260901/primary/16_full_chain_outcomes.png)

G/L_direct/L_ancestor remain complementary association semantics. For example, W parentless efficiency is 4.475/3.938/4.695% and P8O explicit ISR efficiency is 3.361/3.013/3.530%. These changes must not be called detector-efficiency variation.

## 11. Frozen interpretation

The large global W/P8O photon difference is primarily driven by different non-tau truth-photon populations and their phase-space composition. P8O contains an explicit-ISR population that is substantially softer and more beam-collinear than the unresolved W parentless population. Tau-decay photons have closely matching truth kinematics and PFO association efficiencies. Conditional photon PID, given a unique PFO, is also stable across samples and origins. Thus the global migration-matrix difference is dominated by generator/sample truth composition plus truth-to-PFO association/reconstructability, not by a large difference in reconstructed photon PID.

## 12. Explicit limitations

1. W `parentless_unresolved` photons cannot be called ISR.
2. P8O `explicit_ISR` is supported by stored generator genealogy.
3. `tau_direct_daughter_unresolved` photons cannot currently be called FSR.
4. G/L_direct/L_ancestor are association semantics, not detector definitions.
5. This programme is not a generator-physics validation of WHIZARD versus PYTHIA8.
6. Raw all-photon outcomes cannot support a global claim that one reconstruction reconstructs photons better.
7. Detailed equal-statistics audits use 18,000 events per sample.

These are limitations of stored event-record information, not pending software bugs.

## 13. Freeze declaration and provenance

**STATUS: SCIENTIFIC ANALYSIS FROZEN AT 2026-09-01.**

Frozen scope: W/P8C/P8O migration comparison; equal-statistics W/P8O tau topology; conditional PID; PFO association efficiency; photon-origin discovery; and photon-origin reconstruction/PID factorization. The individual reports, machine-readable tables, validation summaries, scripts, and figures remain authoritative sources. Their sizes and SHA-256 hashes are recorded under `artifacts/fcc_tau_w_p8o_consolidated_freeze_20260901/`.

No simulation, reconstruction, generator, linker, HitAnalysis, or scientific audit was rerun to prepare this document. No scientific definition was changed and no derived output was written to Lustre.
