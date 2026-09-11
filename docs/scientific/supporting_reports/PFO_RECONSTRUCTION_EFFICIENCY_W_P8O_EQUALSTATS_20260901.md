# PFO reconstruction/association efficiency: W versus P8O (2026-09-01)

## 1. Scope

This read-only audit uses the exact frozen 18,000-event W and P8O manifests. It measures truth-to-PFO association/reconstruction. It does **not** measure PID correctness. No REC file was reopened and no matching was rerun.

## 2. Definition of reconstruction/association efficiency

Primary `L_ancestor` efficiency is `epsilon_unique = N(associated_unique) / N(selected truth)`. Unmatched and ambiguous outcomes remain separate; ambiguous is excluded from the numerator. The three outcome fractions close to 100%.

Binomial uncertainties are the unweighted Wald one-standard-deviation value `sqrt(epsilon*(1-epsilon)/N_truth)`. Zero-denominator bins are omitted and bins with `N_truth < 25` are marked low-statistics.

## 3. Absolute truth populations

Absolute plots are event-normalized first. Physical E/p/pT variables use geometrically spaced bins and logarithmic x-axes; no `log10(X)` transform is used.

## 4. Species-level efficiencies

| Species | W selected | W epsilon | P8O selected | P8O epsilon | P8O-W |
|---|---:|---:|---:|---:|---:|
| electron | 6,969 | 96.499% | 6,959 | 95.732% | -0.767 pp |
| muon | 6,349 | 98.504% | 6,313 | 97.972% | -0.531 pp |
| photon | 74,344 | 54.447% | 96,625 | 42.046% | -12.401 pp |
| charged pion | 33,416 | 97.100% | 34,030 | 96.703% | -0.397 pp |

## 5. Tau-ancestor versus non-tau-ancestor

The origin flag is the frozen recursively stored `abs(PDG)==15` ancestry definition. No absent ancestry was re-inferred from REC. For W electrons/muons, the frozen full-inventory equality `N_all=N_tau_origin` supplies the already-established flag coverage; persisted compact flags are used elsewhere.

## 6. Photon efficiency

| Origin | W selected | W epsilon | P8O selected | P8O epsilon |
|---|---:|---:|---:|---:|
| all photons | 74,344 | 54.447% | 96,625 | 42.046% |
| tau-origin | 48,465 | 81.013% | 47,504 | 81.873% |
| non-tau-origin | 25,879 | 4.695% | 49,121 | 3.530% |

The P8O/W non-tau photon multiplicity ratio is 1.898 (49,121/25,879), but these photons have only 3.530%/4.695% PFO association efficiency. The excess is overwhelmingly in low-efficiency phase space: 44,437 P8O versus 21,419 W non-tau photons have E<0.1 GeV, with efficiencies 0.230% and 0.313%; 47,778 versus 24,880 have pT<0.1 GeV, with efficiencies 1.243% and 1.632%. Median non-tau photon energies are 1.41e-6 GeV (P8O) and 3.17e-4 GeV (W), and median pT values are 9.88e-9 and 2.86e-5 GeV. Thus the large P8O excess does occupy extremely soft, very-low-efficiency phase space. Near |cos(theta)|>0.98, non-tau efficiencies fall to 0.887% P8O and 2.048% W.

## 7. Charged-pion efficiency

Charged-pion efficiency is 97.100% W and 96.703% P8O globally, but falls to 85.584%/82.474% for pT<1 GeV and 69.342%/68.432% for |cos(theta)|>0.98. Among the strict association losses, ambiguity (542 W, 621 P8O) is slightly larger than unmatched (427, 501). The separate PID audit found 2,509/2,644 associated-but-misidentified pions, so the larger pion loss after truth selection is PID after reconstruction, not absence of a unique PFO.

## 8. Electron/muon efficiency

Electron efficiency is 96.499% W and 95.732% P8O; muon efficiency is 98.504% and 97.972%. Both degrade at low pT and forward/backward angles. For pT<1 GeV, electron efficiency is 68.196%/63.188% and muon efficiency 86.120%/80.357%; at |cos(theta)|>0.98 the corresponding values are 71.751%/74.725% and 59.487%/56.853%. Electron loss is unmatched-dominated. Muon loss is slightly ambiguity-dominated in W (51 versus 44) and unmatched-dominated in P8O (71 versus 57). Conditional PID losses are larger than association losses for both species.

## 9. Charged versus neutral truth

| Charge group | W epsilon | P8O epsilon |
|---|---:|---:|
| charged | 97.201% | 96.730% |
| neutral | 54.447% | 42.046% |

Charged truth is reconstructed much more efficiently than neutral truth in this frozen four-species selection. The neutral category is exactly the photon category here, so photons drive the difference. Tau ancestry barely changes the charged aggregate (97.201% W and 96.742% P8O for tau-origin charged truth; only six P8O charged particles and none in W are non-tau), but it changes neutral efficiency from 81.013%/81.873% for tau-origin photons to 4.695%/3.530% for non-tau photons.

## 10. Association-definition robustness

G and L_direct are reported only as association-definition dependence. Differences from L_ancestor must not be interpreted as changes in detector efficiency. W/P8O efficiencies for G, L_direct, and L_ancestor are respectively: electron 93.184/91.737%, 95.006/93.907%, 96.499/95.732%; muon 98.236/97.481%, 98.519/97.988%, 98.504/97.972%; photon 53.017/38.621%, 46.772/36.134%, 54.447/42.046%; charged pion 94.703/93.553%, 94.530/93.964%, 97.100/96.703%. Photon reconstructability is therefore the most definition-dependent.

## 11. Relation to PID audit

A uniquely associated truth particle counts as reconstructed regardless of its PFO PID. Photon performance loss is overwhelmingly association/reconstruction: W has 33,729 unmatched and 137 ambiguous photons versus 2,656 associated-but-misidentified photons; P8O has 55,836, 162, and 2,600. For pions, electrons, and muons, the conditional PID-misidentified counts exceed the strict association losses. The two audits therefore locate different bottlenecks without redefining either denominator.

## 12. Conclusions

The W/P8O agreement is close for charged species (within 0.77 percentage points) and for tau-origin photons (0.86 pp), while the 12.40 pp global photon difference is caused by the much larger, very soft non-tau photon population in P8O. For every species the dominant outcome is associated_unique; among losses, unmatched dominates electrons and photons, ambiguity narrowly dominates W muons and charged pions, and unmatched dominates P8O muons.

The numerical conclusions, exact outcome partitions, binned denominators, and uncertainties are available in the CSV tables. The figure set contains 25 primary plots (13 pages) and 30 backup plots (15 pages). All 28 pages were rendered and visually inspected with no blank, clipped, distorted, or unreadable page found.

No ISR/FSR label is assigned to non-tau photons.
