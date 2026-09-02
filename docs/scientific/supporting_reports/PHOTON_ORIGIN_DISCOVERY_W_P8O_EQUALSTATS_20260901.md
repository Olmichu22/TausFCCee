# Photon-origin discovery: W versus P8O (2026-09-01)

**Scientific status: PASS_WITH_WARNINGS.** Integrity and coverage pass, but cross-generator physical equivalence is limited by different event-record bookkeeping.

## 1. Scope

Read-only event-record/genealogy audit of the exact frozen 18,000-event manifests. It uses only MCParticle truth fields and the frozen selected-photon keys; it performs no reconstruction, PID, PFO-efficiency, matching, or association analysis.

## 2. Frozen selected-photon sample

Manifest SHA-256: `b2f81449579d33ec228b65c20a7927f8930ef19368aedd97244085c6e587cf3b`. Selected photons reproduce exactly: W 74,344, P8O 96,625. Tau-origin counts reproduce exactly: W 48,465, P8O 47,504.

## 3. Event-record integrity

All 36,000 events passed relation-range/index and full stored parent/daughter graph cycle checks. Kinematics are finite and truth keys are unique. MCParticle provenance is the `MCParticles` collection in each manifest REC, with `_MCParticles_parents` and `_MCParticles_daughters` relation arrays.

## 4. W photon genealogy

All 25,879 W non-tau photons have zero stored parents and generatorStatus 1. Their stored production vertices and times are identically zero. They may have zero or simulation-created daughters—typically conversion-like e+/e- descendants—but no stored ancestor or incoming-electron relation exists. Collection order, simulatorStatus, endpoints, and beam-collinear kinematics do not supply positive origin evidence. They therefore remain `parentless_unresolved`; kinematic similarity cannot restore missing genealogy.

## 5. P8O photon genealogy

All 49,121 P8O non-tau photons have exactly one parent. Immediate parents are e- (PDG +11): 23,030 (46.884%), e+ (PDG -11): 23,100 (47.027%), and photon: 2,991 (6.089%). Electron-parent statuses are 4 for 18,000 e- and 18,000 e+, 41 for 5,029 e- and 5,100 e+, plus one status-44 e-. Photon-parent statuses are 43 (2,709) or 44 (282). No fact was hard-coded.

## 6. Tau-origin photons

Descriptive classes preserve pi0 descendants, direct tau daughters, and other tau descendants. Direct tau daughters number 8,640 W and 7,697 P8O. In W, 8,475 have a status-2 parent tau with one tau-copy daughter and 165 have two tau-copy daughters. In P8O, 6,186 have a status-23 parent with one tau-copy daughter, while 740/716 have status 52/51 parents embedded between tau copies; only 55 terminal status-2 cases lack a tau-copy daughter and show ordinary decay siblings. Median E_gamma/E_tau is 4.21e-4 W and 5.13e-4 P8O; median angular distance is 0.531 and 0.492 rad. This mixture of tau-copy radiation/bookkeeping and terminal decay configurations is not a mutually exclusive proof of tau FSR.

## 7. Non-tau photons

W parentlessness and P8O explicit parent chains are generator-bookkeeping differences. Both populations are predominantly soft and beam-directed, but P8O is more extreme: E<0.1 GeV contains 90.46% P8O versus 82.77% W; pT<0.1 GeV contains 97.27% versus 96.14%; |cos(theta)|>0.98 contains 83.75% versus 50.19%. These comparisons are descriptive and do not classify W photons.

## 8. Initial-state electron-chain representation

Representative W events do contain recognizable parentless, beam-energy e-/e+ objects (generatorStatus 2, E=45.605129 GeV, |cos(theta)|=0.9998875), but the W non-tau photons carry no ancestry relation to them; they therefore provide event context, not usable photon-origin evidence. Positive stored photon-to-beam-chain evidence is found for W 0 and P8O 49,120 selected photons. The P8O chains terminate at explicit parentless generatorStatus-4 e-/e+ objects with E=45.595117 GeV and |cos(theta)|=0.9999755. The rule requires a unique e+/e- ancestor and a same-sign electron-copy chain ending at that stored beam object. This covers all 2,991 photon-parent chains: their status-43/44 parent photons trace in two to four steps to a status-41 electron and then the beam chain. The sole exception is a 1.22e-6 GeV photon emitted by a status-44, 0.108 GeV secondary electron reached through an intervening photon/electron chain; it correctly remains `electron_parent_non_tau_unresolved`.

## 9. FSR-like genealogy patterns

No `explicit_FSR` category is frozen. Direct tau daughters frequently sit beside a continuing tau copy, but some are terminal decay daughters; status, sibling, angle, and energy-fraction patterns overlap. The one non-tau secondary-electron case is also insufficient to define a general final-state-radiation category. All retain descriptive/unresolved labels.

## 10. ISR-like genealogy patterns

`explicit_ISR` is used only where positive stored incoming-electron-chain evidence passes the stated genealogy/status/beam-object rule. It is not assigned from parent PDG, parentlessness, tau status, or angle alone.

## 11. Parentless photons

W parentless non-tau photons contain no ancestry relation to an incoming e+/e- object. They remain unresolved even when soft and beam-collinear.

## 12. Proposed photon_origin_v1

Precedence is: `tau_decay_pi0`; `tau_direct_daughter_unresolved`; `tau_other_descendant`; positive-evidence `explicit_ISR`; `parentless_unresolved`; unresolved electron/positron/photon-parent classes; `resolved_non_tau_other`. The CSV rule table contains exact Boolean conditions and limitations.

Coverage:

| category | W count | W fraction | P8O count | P8O fraction |
|---|---:|---:|---:|---:|
| electron_parent_non_tau_unresolved | 0 | 0.000% | 1 | 0.001% |
| explicit_ISR | 0 | 0.000% | 49,120 | 50.836% |
| parentless_unresolved | 25,879 | 34.810% | 0 | 0.000% |
| tau_decay_pi0 | 39,585 | 53.246% | 39,727 | 41.115% |
| tau_direct_daughter_unresolved | 8,640 | 11.622% | 7,697 | 7.966% |
| tau_other_descendant | 240 | 0.323% | 80 | 0.083% |

## 13. Limitations

The common classification is exhaustive because unresolved categories are explicit. Generator bookkeeping is not assumed equivalent: a W parentless photon cannot be promoted to ISR merely because a P8O photon with similar kinematics has an explicit beam chain.

## 14. Conclusions

The event record supports descriptive tau-decay ancestry in both samples and positive `explicit_ISR` only for the 49,120 P8O photons satisfying the explicit incoming-electron-chain rule. The identical rule finds no W cases because W omits the relevant ancestry, so cross-generator category fractions are not physically symmetric. The record does not support a sufficiently unambiguous `explicit_FSR` category in this version. All remaining physically underdetermined populations retain unresolved labels.
