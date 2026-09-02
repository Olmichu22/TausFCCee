# Photon-origin topology, reconstruction and PID audit (2026-09-01)

## 1. Scope

Read-only join of frozen photon_origin_v1, truth-to-PFO association outcomes, and conditional PID. P8O `explicit_ISR` retains positive genealogy evidence. W `parentless_unresolved` is not relabelled ISR. Direct tau daughters are not labelled FSR.

## 2. Frozen photon_origin_v1

All 170,969 selected-photon keys join one-to-one. The exact 18,000-event manifests and frozen categories are unchanged.

## 3. Origin populations and yields

| Sample/category | N truth | photons/event | Lancestor association | conditional photon PID |
|---|---:|---:|---:|---:|
| P8O electron_parent_non_tau_unresolved | 1 | 0.000056 | 0.000% | N/A |
| P8O explicit_ISR | 49,120 | 2.728889 | 3.530% | 94.348% |
| P8O tau_decay_pi0 | 39,727 | 2.207056 | 90.757% | 93.491% |
| P8O tau_direct_daughter_unresolved | 7,697 | 0.427611 | 36.001% | 94.551% |
| P8O tau_other_descendant | 80 | 0.004444 | 83.750% | 94.030% |
| W parentless_unresolved | 25,879 | 1.437722 | 4.695% | 94.239% |
| W tau_decay_pi0 | 39,585 | 2.199167 | 91.350% | 93.393% |
| W tau_direct_daughter_unresolved | 8,640 | 0.480000 | 33.854% | 93.573% |
| W tau_other_descendant | 240 | 0.013333 | 73.750% | 94.915% |

## 4. W parentless versus P8O explicit ISR truth topology

The two dominant non-tau populations are **partially similar**, never physically identified. Both are overwhelmingly soft and forward/backward, but P8O explicit ISR is substantially more extreme. W/P8O multiplicities are 1.438/2.729 photons per event. Median E is 3.17e-4/1.41e-6 GeV and median pT is 2.86e-5/9.88e-9 GeV. Fractions with E<0.1 GeV are 82.77%/90.46%; pT<0.1 GeV 96.14%/97.27%; |cos(theta)|>0.98 50.19%/83.75%. KS effect-size distances are 0.298 for E, 0.444 for pT and 0.367 for cos(theta). Thus W parentless photons may be described as ISR-like kinematically, while retaining their unresolved label.

## 5. Tau-decay pi0 reference

The frozen pi0-descendant population provides the well-defined tau-decay control: 39,585 W and 39,727 P8O photons, or 2.199/2.207 per event. Association efficiencies are 91.350%/90.757%, and conditional photon-PID accuracies are 93.393%/93.491%. Reconstructed-PFO median E is 3.814/3.628 GeV and median pT is 2.613/2.529 GeV. This confirms close W/P8O agreement.

## 6. Direct tau-daughter population

These are direct tau daughters in the stored event record; radiative interpretation remains unresolved. Counts are 8,640 W and 7,697 P8O (0.480/0.428 per event). Association is 33.854%/36.001%; conditional photon PID is 93.573%/94.551%. The strong loss is concentrated below 0.1 GeV, where association is only 1.321%/1.360%, while it reaches about 99% above 1 GeV. Associated-PFO median E and pT are 1.524/1.603 GeV and 1.091/1.176 GeV. Frozen E_gamma/E_tau and photon-tau angular variables remain available from the discovery artifact.

## 7. PFO association efficiency by origin

Association uses frozen L_ancestor only for primary results. P8O explicit ISR has 1,734/49,120 unique associations (3.530%), 47,374 unmatched and 12 ambiguous. W parentless has 1,215/25,879 (4.695%), 24,654 unmatched and 10 ambiguous. Pi0-descendant association is about 91%, while direct tau daughters are 34–36%. Unmatched and ambiguous remain outside the numerator and every origin closes exactly.

## 8. Efficiency versus E/pT/angle

Physical variables, logarithmic x-axes and geometric bins are used for E/p/pT. Efficiency uncertainties are Wald binomial one-sigma; plotted bins require N_truth>=25. For W parentless/P8O ISR, efficiency is 0.313%/0.230% below 0.1 GeV, 25.225%/34.264% over 0.1–1 GeV and 27.765%/37.083% above 1 GeV. It is 1.632%/1.243% below pT=0.1 GeV, rising above 94%/96% for pT>=1 GeV. At |cos(theta)|>0.999 it is only 0.110%/0.047%; in |cos(theta)|<0.8 it is 6.115%/14.022%. Their different integrated efficiencies are therefore largely phase-space-composition effects.

## 9. Conditional PID by origin

PID is conditioned only on associated_unique. P8O ISR photon-PID accuracy is 94.348%; its dominant wrong categories are neutron 3.114%, electron 1.326% and charged pion 0.865%. W parentless accuracy is 94.239%, with neutron 2.963%, electron 1.728% and charged pion 0.905%. Pi0 and direct-tau categories also lie at 93.4–94.6%. Origin dependence is therefore much larger in association than in conditional PID.

## 10. Associated-PFO kinematics

All reconstructed curves are explicitly titled ASSOCIATED RECONSTRUCTED PFO and use the frozen representative PFO rows. For W parentless/P8O ISR, reconstructed median E is 0.387/0.355 GeV and median pT is 0.162/0.152 GeV: the rare associated subsets look much more alike than the full truth populations.

## 11. Full origin -> association -> PID outcomes

The stacked outcome figure normalizes no-unique-PFO, ambiguity, associated photon PID, associated electron PID and associated other PID to all truth photons in each origin. Only 4.424% of all W parentless photons and 3.330% of all P8O ISR photons reach an associated PFO with photon PID. Corresponding values are 85.316%/84.845% for pi0 descendants and 31.678%/34.040% for direct tau daughters.

## 12. Association-definition robustness

G/Ldirect/Lancestor differences are association-definition dependence, not detector-efficiency variation. W parentless efficiencies are 4.475/3.938/4.695%; P8O ISR 3.361/3.013/3.530%; W/P8O pi0 descendants 88.768/78.388/91.350% and 82.483/77.811/90.757%; direct tau daughters 34.028/29.722/33.854% and 36.846/32.090/36.001%.

## 13. Limitations

Kinematic similarity cannot supply missing W genealogy. The factorization diagnoses performance but does not revise origin identity.

## 14. Conclusions

P8O explicit ISR: association 3.530%, conditional photon PID 94.348%. W parentless unresolved: association 4.695%, conditional photon PID 94.239%. Their very low end-to-end photon outcomes are driven by the overwhelmingly soft/beam-collinear truth populations and association loss, not conditional PID. The earlier W/P8 migration-matrix discrepancy is therefore understandable as an origin-mixture and association effect: P8O contains 1.898 times more dominant non-tau photons, concentrated in even lower-efficiency phase space, while tau-decay controls and conditional PID agree closely.
Figures: 16 primary (8 pages) and 27 backup (14 pages).
