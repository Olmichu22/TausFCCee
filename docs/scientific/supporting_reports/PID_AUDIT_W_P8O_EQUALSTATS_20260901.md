# Conditional PID audit: W versus P8O (2026-09-01)

## 1. Scope

Read-only PID audit of the exact frozen equal-statistics samples (18,000 events per sample). It reuses historical G and frozen L_direct/L_ancestor assignments. No REC file was reread.

## 2. PID denominator definition

This audit conditions on the existence of a unique associated representative PFO. It does **not** measure PFO reconstruction efficiency. `associated_unique` is the frozen inversion result with one representative PFO; unmatched and ambiguous outcomes are reported separately and excluded from the PID denominator.

`correct_PID` means reconstructed PFO PID category equals the selected truth species. `misidentified_PID` means a unique representative exists but its reconstructed PID differs. These two categories sum to 100% of associated-unique truth particles.

## 3. Association versus PID distinction

Unmatched is association/reconstruction loss, not misidentification. `ambiguous_multiple_pfo` is also not misidentification. The frozen PID mapping is |PDG| 11 electron, 13 muon, 22 photon, 211 charged pion, 310 K0S, 2112 neutron, 3122 Lambda; sentinel 999 is not a PID.

## 4. L_ancestor conditional PID matrices

| Truth | W associated | W correct | W accuracy | P8O associated | P8O correct | P8O accuracy | P8O-W |
|---|---:|---:|---:|---:|---:|---:|---:|
| electron | 6,725 | 6,042 | 89.844% | 6,662 | 5,958 | 89.433% | -0.411 pp |
| muon | 6,254 | 5,544 | 88.647% | 6,185 | 5,489 | 88.747% | +0.100 pp |
| photon | 40,478 | 37,822 | 93.438% | 40,627 | 38,027 | 93.600% | +0.162 pp |
| charged pion | 32,447 | 29,938 | 92.267% | 32,908 | 30,264 | 91.965% | -0.302 pp |

## 5. Electron PID

Dominant wrong destinations: W — photon 52.6%, charged pion 42.9%, neutron 3.7%; P8O — photon 60.4%, charged pion 33.4%, neutron 4.7%.
W: selected 6,969, associated 6,725, unmatched 234, ambiguous 10, misidentified 683. P8O: selected 6,959, associated 6,662, unmatched 275, ambiguous 22, misidentified 704.
Failures are localized at low scale and extreme angles: truth pT medians are 10.26->1.90 GeV W and 10.20->2.37 GeV P8O; the mis-ID truth |cos(theta)|>0.98 fractions are 18.6%/19.3%.

## 6. Muon PID

Dominant wrong destinations: W — charged pion 72.4%, neutron 26.9%, electron 0.3%; P8O — charged pion 67.2%, neutron 31.6%, photon 0.9%.
W: selected 6,349, associated 6,254, unmatched 44, ambiguous 51, misidentified 710. P8O: selected 6,313, associated 6,185, unmatched 71, ambiguous 57, misidentified 696.
Failures are strongly low-momentum and partly forward/backward: truth pT medians are 11.25->1.27 GeV W and 11.68->1.38 GeV P8O; mis-ID truth |cos(theta)|>0.98 is 16.3%/16.1%.

## 7. Photon PID

Dominant wrong destinations: W — electron 59.6%, neutron 20.1%, charged pion 14.1%; P8O — electron 60.2%, neutron 20.9%, charged pion 12.4%.
W: selected 74,344, associated 40,478, unmatched 33,729, ambiguous 137, misidentified 2,656. P8O: selected 96,625, associated 40,627, unmatched 55,836, ambiguous 162, misidentified 2,600.
Truth-energy medians are similar for correct and mis-ID (3.76->3.53 GeV W; 3.51->3.42 GeV P8O), while truth pT is modestly lower. The representative-PFO energy median drops to 1.54/1.54 GeV and pT to 0.95/0.92 GeV. No enhanced extreme-forward/backward concentration is seen in the mis-ID subset.

## 8. Charged-pion PID

Dominant wrong destinations: W — neutron 36.1%, electron 33.2%, muon 12.0%; P8O — neutron 38.0%, electron 30.4%, photon 11.9%.
W: selected 33,416, associated 32,447, unmatched 427, ambiguous 542, misidentified 2,509. P8O: selected 34,030, associated 32,908, unmatched 501, ambiguous 621, misidentified 2,644.
Truth pT medians fall from 7.74 to 2.99 GeV in W and from 7.40 to 3.15 GeV in P8O. The truth |cos(theta)|>0.98 fraction rises from about 0.1% to 27.5%/25.3% (W/P8O). Representative-PFO pT medians are 2.64/2.69 GeV for mis-ID, showing the same low-pT and forward/backward localization.

## 9. W versus P8O comparison

All four L_ancestor conditional PID accuracies agree within 0.42 percentage points. Destination patterns are likewise close; no meaningful generator-dependent PID shift is established.

## 10. G/L_direct/L_ancestor robustness

Association semantics visibly affect apparent electron and photon PID. Electron accuracy spans 82.045/85.495% (W/P8O) under G, 91.225/91.247% under L_direct, and 89.844/89.433% under L_ancestor. Photon accuracy spans 96.042/95.402%, 98.893/98.860%, and 93.438/93.600%. Muon and pion conditional accuracies are comparatively stable, within about 0.44 and 1.05 percentage points across definitions respectively. See `pid_truth_definition_robustness.csv` and the robustness figure.

## 11. Conclusions

Photon truth-to-reco loss is dominated by unmatched association, not wrong PID: W 33,729 unmatched versus 2,656 associated-but-wrong; P8O 55,836 versus 2,600. Ambiguous counts are 137 and 162 and remain separate.
For charged pions, associated-but-wrong PID is larger than unmatched or ambiguous separately: W 2,509 misidentified, 427 unmatched, 542 ambiguous; P8O 2,644, 501, 621.
Kinematic localization is reported descriptively in the correct/mis-ID truth- and reco-side summaries and problem-region table; no ISR/FSR or detailed interaction-topology label is inferred.

Figures: 15 primary (8 PDF pages) and 41 backup (21 PDF pages).
