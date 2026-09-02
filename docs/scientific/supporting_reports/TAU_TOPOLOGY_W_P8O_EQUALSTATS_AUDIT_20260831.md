# W/P8O equal-statistics tau-topology audit (2026-08-31)

## 1. Scope and sample definition

Read-only comparison of validated WHIZARD plus our reconstruction (W) and PYTHIA8 plus our reconstruction (P8O). No P8C input is used. No simulation, reconstruction, linker, or HitAnalysis step was rerun.

## 2. Equal-statistics event selection

- `N_W_available = 1,992,000`
- `N_P8O_available = 18,000`
- `N_AUDIT = 18,000` events per sample
- Ordering: numeric source suffix, then zero-based `event_in_file`; the first exact `N_AUDIT` events are used.

W ranges:
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_000242385/events_000242385_REC.edm4hep.root` events [0, 2000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_000358005/events_000358005_REC.edm4hep.root` events [0, 2000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_000497029/events_000497029_REC.edm4hep.root` events [0, 2000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_000879444/events_000879444_REC.edm4hep.root` events [0, 2000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_000898343/events_000898343_REC.edm4hep.root` events [0, 2000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_000966520/events_000966520_REC.edm4hep.root` events [0, 2000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_001218437/events_001218437_REC.edm4hep.root` events [0, 2000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_001434184/events_001434184_REC.edm4hep.root` events [0, 2000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/ILD20260821_2k/outputs/events_001714006/events_001714006_REC.edm4hep.root` events [0, 2000)

P8O ranges:
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_1_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_32_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_118_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_246_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_256_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_301_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_366_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_448_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_485_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_518_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_545_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_590_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_600_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_701_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_766_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_791_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_867_REC.edm4hep.root` events [0, 1000)
- `/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/samples_pythia_20260821/our_reco/outputs/out_sim_edm4hep_984_REC.edm4hep.root` events [0, 1000)

The composite event identifier is `source-file-id:event_in_file`, because `EventHeader.eventNumber` restarts between files.

## 3. Validation

All blocking checks pass: 18,000 unique composite events per sample, valid stored parent/daughter indices, no traversal cycles, present MCParticle/PandoraPFO collections, finite four-momentum-derived quantities, and exact truth-origin/PFO-charge/photon-subclass reconciliations. Fixed bins are recorded in `tables/histogram_binning.csv`.

Selected truth reuses the frozen implementations: W `extract_G_Ldirect_Lancestor.selected_indices` (including the stored SiTracks truth-link signal restriction when active) and P8O `extract_pythia_Lancestor.selected`. Both require generatorStatus 1, exclude neutrinos, and require nonzero momentum. Tau origin is recursively stored-parent ancestry to `abs(PDG)==15`, with deterministic nearest-tau depth and no inferred ancestry.

## 4. Tau topology and decay channels

Exactly two stored tau objects occur in W 61.394% and P8O 65.311% of events. The terminal decay-tau pair occurs in W 100.000% and P8O 99.994%. P8O's extra stored tau multiplicity is retained as a record-level topology difference, while decay channels use only terminal tau representatives (no stored tau daughter) to avoid counting successive copies twice.

- W: hadronic_1prong_with_pi0 36.797%, electron_leptonic 18.336%, muon_leptonic 17.583%, hadronic_1prong_no_pi0 12.256%, hadronic_3prong_no_pi0 9.744%, hadronic_3prong_with_pi0 5.153%, other_hadronic 0.131%
- P8O: hadronic_1prong_with_pi0 36.892%, electron_leptonic 17.688%, muon_leptonic 17.493%, hadronic_1prong_no_pi0 12.144%, hadronic_3prong_no_pi0 9.991%, hadronic_3prong_with_pi0 5.664%, other_hadronic 0.128%

Visible signatures use stable visible descendants, exclude neutrinos, and preserve each explicit pi0 as one object (its photon daughters are not separately counted in that signature).

- W leading signatures: pi0 pi+ (12.742%), pi- pi0 (12.478%), e+ (8.969%), e- (8.917%), mu- (8.847%), mu+ (8.733%), pi- (5.492%), pi+ (5.417%)
- P8O leading signatures: pi0 pi+ (12.666%), pi- pi0 (12.491%), e- (8.902%), mu- (8.863%), e+ (8.711%), mu+ (8.630%), pi+ (5.555%), pi- (5.422%)

## 5. Tau kinematics

The terminal-tau distributions are very close. Mean p is 44.960 GeV (W) versus 45.007 GeV (P8O); mean pT is 33.142 versus 33.279 GeV; mean theta is 89.961 versus 89.983 degrees. The p medians are 45.446 and 45.476 GeV and the p q05-q95 intervals are [43.349, 46.129] and [43.608, 45.816] GeV. Revised common-bin unit-area and event-normalized overlays are in `plot_revision_20260901/plots/primary/03-06_*`; full mean/std/median/q05/q50/q95 values remain in `tables/kinematic_summary.csv`.

## 6. Tau-origin truth population

- W: 96,297 particles (5.350/event); photon 50.329%, charged_pion 34.701%, electron 7.237%, muon 6.593%, other_charged_hadron 0.791%, neutral_hadron 0.349%.
- P8O: 96,011 particles (5.334/event); photon 49.478%, charged_pion 35.444%, electron 7.242%, muon 6.575%, other_charged_hadron 0.883%, neutral_hadron 0.378%.

## 7. Tau-origin photon subclasses

The mutually exclusive precedence is explicit-pi0 ancestor, then immediate-tau parent, then other stored tau-descendant chain.

- W: {'pi0_descendant': 39585, 'direct_tau_daughter': 8732, 'other_tau_descendant': 148}
- P8O: {'pi0_descendant': 39727, 'direct_tau_daughter': 7697, 'other_tau_descendant': 80}

Thus pi0 descendants constitute 81.677% (W) and 83.629% (P8O); direct tau daughters 18.017% and 16.203%; other chains 0.305% and 0.168%. Their inclusive photon kinematics are also close: mean energy 5.100 versus 5.018 GeV, mean pT 3.756 versus 3.727 GeV, and mean theta 90.539 versus 90.334 degrees.

Equal-subset tau-origin photon fractions are W 65.190% and P8O 49.163%. The full-sample 65.481% and 49.163% values are retained only as context, not imposed on this subset.

## 8. Non-tau-origin truth population

- W: 25,879 particles (1.438/event); photon 100.000%.
- P8O: 49,127 particles (2.729/event); photon 99.988%, electron 0.012%.

Here “non-tau-origin” means only “no recursively stored tau ancestor”; it is not labelled ISR/FSR or assigned another production mechanism.

## 9. Non-tau photon diagnostic

- W: 25,879 photons, 1.438/event, 34.810% of selected photons; stored-parent multiplicity {'zero': 25879}.
- P8O: 49,121 photons, 2.729/event, 50.837% of selected photons; stored-parent multiplicity {'one': 49121}.

Energy/p/pT/theta/cos(theta), threshold fractions, immediate-parent PDGs, nearest non-photon ancestors, and matched-scale 2D maps are tabulated/plotted without interpreting parentless photons as a specific physical mechanism.

P8O has 1.898 times the non-tau-photon yield. Its photons are softer and more forward/backward: mean energy 0.100 versus 0.217 GeV in W; mean pT 0.0199 versus 0.0312 GeV; `|cos(theta)|>0.98` fraction 83.748% versus 50.191%. All W photons have zero stored parents. All P8O photons have one: 47.027% have immediate parent PDG -11, 46.884% PDG +11, and 6.089% PDG 22; nearest non-photon ancestors divide 50.137%/-11 and 49.863%/+11. These are stored-record facts, not production-mechanism labels.

## 10. Reconstructed charged/neutral PFO population

- W: all 6.035/event, charged 2.729/event, neutral 3.306/event.
- P8O: all 6.114/event, charged 2.741/event, neutral 3.373/event.

The primary PFO audit uses all reconstructed PandoraPFOs and the strict stored reconstructed charge split, with no truth matching or fiducial cut. Kinematic and track/cluster summaries are in the machine-readable tables and backup plots.

Charged-PFO mean energy is 12.321 versus 12.146 GeV, mean pT 9.393 versus 9.279 GeV, and mean theta 90.123 versus 90.060 degrees. Neutral-PFO mean energy is 4.872 versus 4.691 GeV, mean pT 3.398 versus 3.309 GeV, and mean theta 90.558 versus 90.493 degrees. Multiplicity and spectra are therefore close, with small W-hardening and P8O-yield shifts.

### 10.1 Plotting and reconstructed-PID extension (2026-09-01)

This is a plotting/presentation extension of the frozen 18,000-event-per-sample audit, not a new event selection or truth analysis. Existing cached truth products were reused. The exact 27 REC ranges in `tables/sample_manifest_equalstats.csv` were reread read-only only for the `PandoraPFOs` fields PDG, charge, energy, momentum, tracks, and clusters. No truth association, reconstruction, linker, HitAnalysis, fiducial cut, or efficiency calculation was run.

The reconstructed identity mapping is the frozen migration-matrix destination mapping: absolute PDG 11 = electron, 13 = muon, 22 = photon, 211 = charged pion, 310 = K0S, 2112 = neutron, and 3122 = Lambda. Destination code 999 remains the migration-matrix `unmatched` sentinel and is not a reconstructed identity. No PFO had an unmapped reconstructed PDG, so no `other_reco_pid` entries were needed.

| PID | W count | W / event | P8O count | P8O / event |
|---|---:|---:|---:|---:|
| electron | 9,726 | 0.5403 | 9,599 | 0.5333 |
| muon | 5,877 | 0.3265 | 5,811 | 0.3228 |
| photon | 50,977 | 2.8321 | 51,876 | 2.8820 |
| charged pion | 33,518 | 1.8621 | 33,925 | 1.8847 |
| K0S | 273 | 0.0152 | 266 | 0.0148 |
| neutron | 7,936 | 0.4409 | 8,231 | 0.4573 |
| Lambda | 327 | 0.0182 | 348 | 0.0193 |

The PID sums are exactly 108,634 W and 110,056 P8O PFOs. The independent charge partitions remain exactly 49,121 charged + 59,513 neutral in W and 49,335 + 60,721 in P8O. Reconstructed photon-PFO yield is only 1.76% higher in P8O (an absolute +0.0499/event), despite the 1.898-fold P8O/W ratio in non-tau-origin truth photons. This is a reco-only population comparison and is not an efficiency statement. Photon-PFO spectra are close, with P8O slightly softer: median energy 2.339 versus 2.512 GeV and median pT 1.570 versus 1.674 GeV. Electron, muon, and charged-pion PFO counts, spectra, and angular distributions are likewise close. The largest relative rare-category shifts are neutron +3.7% and Lambda +6.4% in P8O, but the Lambda absolute shift is only 21 objects (0.00117/event); these are not treated as unexpectedly large physics differences.

The revised human-facing terminology is “selected visible truth particles from tau decays.” Both tau decay chains are included and neutrinos are excluded; the composition is not a branching fraction, and event-level multiplicity is summed over both taus. The stored/terminal tau-copy multiplicity plot was removed from the primary scientific set but retained under `plot_revision_20260901/plots/validation/`. Tau-photon genealogy remains the mutually exclusive frozen pi0-descendant/direct-tau-daughter/other-tau-descendant classification; direct tau daughters are not labelled FSR.

Every displayed angular population now has both theta and cos(theta). W and P8O always share bins and axis ranges. Linear momentum ranges use common fixed ranges derived from the combined validated distributions, with underflow/overflow counts recorded in the plot manifest and overflow included in the edge bin. Non-tau photons retain the physical E, p, and pT variables in GeV, using shared logarithmic axes and logarithmically spaced bins from the finite positive support; there were zero non-positive selected values. In each log-axis 1D plot, the left shape panel sums to unity independently for W and P8O and is labelled `Fraction of photons / log bin`, while the right yield panel preserves multiplicity and is labelled `Photons / event / log bin`. The primary 2D diagnostics use theta versus physical E/pT on logarithmic kinematic axes, with shared bins and shared logarithmic color normalization; absolute-cosine variants are backup diagnostics. No displayed variable is transformed to log10(X).

The revised set contains 38 primary and 45 backup plots. The contact sheets contain 19 and 23 pages respectively, and every rendered page was visually inspected for content, clipping, legibility, and aspect-ratio preservation. Outputs are:

- `artifacts/tau_topology_audit_W_P8O_equalstats_20260831/plot_revision_20260901/`
- `docs/TAU_TOPOLOGY_W_P8O_EQUALSTATS_PRIMARY_PLOTS_20260831.pdf`
- `docs/TAU_TOPOLOGY_W_P8O_EQUALSTATS_BACKUP_PLOTS_20260831.pdf`
- pre-revision PDFs and manifest: `artifacts/tau_topology_audit_W_P8O_equalstats_20260831/contact_sheets/archive_pre_plot_revision_20260901/`

## 11. Relation to migration matrices

The equal-statistics result independently tests the full-sample observation that W has a larger selected-photon fraction with stored tau ancestry. Differences in the complementary non-tau-origin photon multiplicity and phase space are therefore qualitatively relevant to the global photon migration-matrix discrepancy, but this descriptive audit alone does not establish causality.

## 12. Conclusions and explicit answers

Status is **PASS_WITH_WARNINGS**: all numerical/integrity validations pass, but P8O event `out_sim_edm4hep_1:838` has 11 stored tau objects and four terminal tau representatives. The other 17,999 P8O events and all 18,000 W events have exactly two terminal representatives.

1. **Basic tau-pair topology:** yes at terminal-decay level (100.000% W, 99.994% P8O), with the single P8O four-terminal-tau event above. Stored-copy multiplicities differ modestly: 2.772/event W and 2.859/event P8O.
2. **Tau p, pT, theta:** comparable; mean differences are 0.047 GeV, 0.137 GeV, and 0.022 degrees respectively.
3. **Decay-channel fractions:** comparable; the largest compact-category difference is electron leptonic, +0.648 percentage points in W, followed by 3-prong with pi0, -0.511 points.
4. **Stable visible contents:** comparable. The same eight signatures lead both samples and each leading-signature fraction differs by at most 0.26 percentage points.
5. **Tau-origin truth spectra:** comparable at this descriptive resolution. Total multiplicities are 5.350/event W and 5.334/event P8O and the principal species fractions differ by less than one percentage point.
6. **Tau-origin photon subclasses:** broadly comparable, with a 1.95-point larger pi0-descendant fraction in P8O and a 1.81-point larger direct-tau fraction in W.
7. **Largest truth-population difference:** non-tau-origin photons, not the tau-origin yield: 25,879 W versus 49,121 P8O, while tau-origin photons are 48,465 versus 47,504.
8. **Non-tau photon size:** 1.438/event and 34.810% of W selected photons versus 2.729/event and 50.837% of P8O selected photons.
9. **Non-tau photon kinematics:** P8O has the higher multiplicity but softer and more beam-aligned support; W has higher mean energy/p and pT and a substantially smaller extreme-forward/backward fraction.
10. **Stored parents:** W has no stored parent for any of these photons; every P8O photon has one stored parent, dominated by e+ or e-, with 6.089% having a photon immediate parent. No mechanism is inferred.
11. **Charged PFOs:** comparable: 2.729 versus 2.741/event, with W mean energy/pT only 1.4%/1.2% higher and essentially identical angular means.
12. **Neutral PFOs:** comparable: 3.306 versus 3.373/event, with W mean energy/pT about 3.9%/2.7% higher and essentially identical angular means.
13. **Migration-matrix connection:** qualitatively consistent. The nearly equal tau-origin photon counts but very different complementary non-tau population naturally shift the global ancestry fraction and can contribute to the observed global photon-composition/migration discrepancy; this audit does not prove causality.
14. **Are tau-origin products more similar than the global population?** yes descriptively: their total yield, species composition, photon subclasses, and main spectra are much closer than the global selected-photon population, whose difference is dominated by non-tau-origin photons.
