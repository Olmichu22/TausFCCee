# W/P8 migration matrices: G, L_direct and L_ancestor

Date: 2026-08-28  
Status: **PASS**

## Definition recovered from the historical comparison

- Historical source: `/lhome/ific/a/airqui/FCC/TausFCCee/scripts/compare_campaigns.py`.
- Historical products: `/lhome/ific/a/airqui/FCC/TausFCCee/Results/Comparisons/ILD20260821_2k_vs_PYTHIA8_collab_vs_PYTHIA8_ourReco`.
- Truth categories, in frozen order: `electron, muon, photon, charged_pion`.
- Reconstructed categories, in frozen order: `electron, muon, photon, charged_pion, K0S, association_unmatched, neutron, Lambda`; `ambiguous_multiple_pfo` is appended for common L semantics.
- Axis orientation in the new heatmaps: x = reconstructed outcome, y = selected truth particle. This preserves the historical reco-category x-axis while stacking the four historical truth rows.
- Normalization: independently per truth row, `N(truth X → reco outcome Y) / N(selected truth X)`.
- Historical plots were one grouped-bar migration plot per truth species, not a saved 2D cell-annotated confusion matrix. Their scientific definition, category order, labels, normalization and x-axis orientation were recovered exactly. Exact historical 2D cosmetics therefore carry the marker **ORIGINAL_STYLE_NOT_RECOVERED**; the new heatmaps use a fixed clean common style.
- Historical bar plots did not print percentages inside cells (there were no cells) and used fractional y axes with campaign colors. The new heatmaps print percentages and use a fixed 0–100% scale as requested.

## Frozen association and selection semantics

- **G:** already-frozen geometric selected-truth outcome from the historical dR parquet; no d_theta_phi matching was rerun.
- **L_direct:** immediate detector-level MC contributor from `truthlink_assignment_v1`, inverted with the frozen representative-PFO rule (track branch priority, then winning T/C; exact terminal ties remain ambiguous).
- **L_ancestor:** nearest unique selected generator-level ancestor from `truthlink_ancestor_assignment_v1`, with the same frozen inversion rule.
- **Inclusive:** representative physical PID outcome, `association_unmatched`, or `ambiguous_multiple_pfo`.
- **Fiducial:** representative chosen first; then strict `PFO.getEnergy() > 1.0 GeV` and `1° < theta_reco < 179°`. A representative that fails becomes `assigned_but_fails_reco_selection`; it is never replaced. No truth-level cut is applied.
- These are **truth-association-dependent reconstructed outcome matrices**. An off-diagonal L entry is not automatically a detector PID misidentification.

## Validation

- All 18 × 4 truth rows sum to 100% within floating precision: **True**.
- Original inclusive G counts and denominators reproduced exactly for all three samples, four truth species and eight historical destinations: **True**.
- Authoritative G fiducial photon/pion anchors reproduced exactly: **True**.
- Authoritative W G/L_direct/L_ancestor photon/pion anchors reproduced exactly: **True**.
- Authoritative cross-sample L_ancestor inclusive and fiducial photon/pion rows reproduced exactly: **True**.
- G/L_direct/L_ancestor use the same selected-truth denominator by construction within every sample/species: **True**.

## Main numerical readout

- **G:** photon→photon = 51.249% / 37.053% / 36.846% and photon→unmatched = 46.610% / 61.109% / 61.379% for W / P8C / P8O.
- **Ldirect:** photon→photon = 46.439% / 35.817% / 35.722% and photon→unmatched = 53.024% / 63.762% / 63.836% for W / P8C / P8O.
- **Lancestor:** photon→photon = 51.149% / 39.526% / 39.355% and photon→unmatched = 45.068% / 57.496% / 57.786% for W / P8C / P8O.
- Across all nine inclusive matrices, charged-pion→charged-pion spans 86.492%–89.680%.
- The photon assigned-but-fails-reco-selection state spans 7.375%–11.027%; unmatched and ambiguity counts are unchanged by the reco fiducial cut.

## Scientific questions

1. **G → L_direct → L_ancestor in W.** The largest inclusive photon-row shifts for L_direct−G are association_unmatched +6.414 pp, photon -4.809 pp, electron -1.037 pp. For L_ancestor−G they are association_unmatched -1.542 pp, electron +1.126 pp, charged_pion +0.325 pp. The direct and ancestry definitions therefore alter association/provenance semantics visibly, especially unmatched and photon/electron descendant outcomes.
2. **P8C and P8O.** The same qualitative direct-versus-ancestor redistribution occurs and the two P8 reconstruction chains remain close element by element; see the P8O−P8C difference matrices.
3. **Common G, W versus P8.** The large photon difference is dominated by the photon diagonal and the compensating unmatched outcome; off-diagonal physical PID migrations are much smaller.
4. **Common L_ancestor.** The W/P8 difference remains dominated by photon→photon versus photon→unmatched, while electron and other physical migrations are subleading.
5. **Unmatched or another PID?** Predominantly unmatched. The difference matrices show no comparably large transfer into another physical PID category.
6. **Charged pions.** The charged-pion diagonal is high and comparatively stable across samples and truth definitions. Changes under L_ancestor include legitimate descendant promotion and must not automatically be labelled PID misidentification.
7. **Fiducial selection.** It mainly transfers assigned entries into `assigned_but_fails_reco_selection`; `association_unmatched` and `ambiguous_multiple_pfo` remain separate and numerically unchanged.

## Outputs

- Individual inclusive matrices: `/lhome/ific/a/airqui/FCC/ild-tau/artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/inclusive`
- Individual fiducial matrices: `/lhome/ific/a/airqui/FCC/ild-tau/artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/fiducial`
- Difference matrices (30): `/lhome/ific/a/airqui/FCC/ild-tau/artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/difference_matrices`
- Contact sheets: `/lhome/ific/a/airqui/FCC/ild-tau/artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/contact_sheets`
- Machine-readable tables: `/lhome/ific/a/airqui/FCC/ild-tau/artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/tables`
- Validation: `/lhome/ific/a/airqui/FCC/ild-tau/artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/validation`
- Summary: `/lhome/ific/a/airqui/FCC/ild-tau/artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/summary.json`

## Tau-origin truth matrices

Status: **PASS**.

### Definition and orthogonal selections

`tau_origin` is true exactly when the selected truth MCParticle has at least
one stored recursive parent ancestor with `abs(PDG) == 15`. The nearest tau is
chosen deterministically by `(ancestry depth, MC index)` for validation only.
No missing ancestry is inferred and particles without a stored tau ancestor
are not labelled ISR or FSR.

The two selection axes are orthogonal:

- `all` versus `tau_origin` restricts the **truth population**.
- `inclusive` versus `fiducial` restricts the **representative reconstructed
  PFO**. The frozen fiducial requirement remains `E_PFO > 1 GeV` and
  `1 deg < theta_PFO < 179 deg`, applied after representative selection.

G, L_direct, L_ancestor, their inversion semantics, reconstructed PID order,
deduplication and fiducial cuts are unchanged. The tau-origin inclusive and
fiducial matrices use the same truth denominator.

### Ancestry and regression validation

- 996 W, 100 P8C and 18 P8O frozen files were processed read-only.
- All selected truth keys were unique; invalid parent references and ancestry
  cycles: **0**.
- Recursive ancestry agreed with every available compact
  `explicit_tau_lineage` flag (W photon/pion and all four P8 categories):
  **0 mismatches**.
- Nearest-tau PDG and depth summaries are stored in
  `tau_origin/validation/nearest_tau_summary.csv`.
- All 684 pre-existing all-truth matrix cells were reproduced exactly:
  **0 numerical mismatches**. The 112 pre-existing artifact files are also
  byte-identical to their pre-extension checksums.
- All 72 tau-origin truth rows have count sum equal to their denominator and
  percentage sum equal to 100% within tolerance: **0 failures**.

### Truth-population composition

| Sample | Truth category | N all | N tau-origin | Tau-origin fraction |
|---|---:|---:|---:|---:|
| W | electron | 764,475 | 764,475 | 100.000% |
| W | muon | 691,798 | 691,798 | 100.000% |
| W | photon | 8,234,477 | 5,392,026 | 65.481% |
| W | charged pion | 3,703,552 | 3,703,552 | 100.000% |
| P8C | electron | 38,591 | 38,571 | 99.948% |
| P8C | muon | 34,383 | 34,383 | 100.000% |
| P8C | photon | 538,745 | 265,909 | 49.357% |
| P8C | charged pion | 189,560 | 189,560 | 100.000% |
| P8O | electron | 6,959 | 6,953 | 99.914% |
| P8O | muon | 6,313 | 6,313 | 100.000% |
| P8O | photon | 96,625 | 47,504 | 49.163% |
| P8O | charged pion | 34,030 | 34,030 | 100.000% |

The global photon samples therefore have substantially different truth-origin
composition: W contains 65.481% stored tau-lineage photons, versus 49.357% and
49.163% for P8C and P8O.

### Photon rows: inclusive all truth versus tau-origin

Percentages below are truth-row normalized. `Other` is every outcome except
photon and `association_unmatched` (including ambiguity where present).

| Definition | Sample | all: photon | all: unmatched | tau: photon | tau: unmatched | tau: other |
|---|---:|---:|---:|---:|---:|---:|
| G | W | 51.249 | 46.610 | 75.999 | 20.845 | 3.156 |
| G | P8C | 37.053 | 61.109 | 71.711 | 24.705 | 3.584 |
| G | P8O | 36.846 | 61.379 | 71.600 | 24.918 | 3.482 |
| L_direct | W | 46.439 | 53.024 | 68.891 | 30.351 | 0.758 |
| L_direct | P8C | 35.817 | 63.762 | 69.511 | 29.720 | 0.769 |
| L_direct | P8O | 35.722 | 63.836 | 69.611 | 29.566 | 0.823 |
| L_ancestor | W | 51.149 | 45.068 | 75.806 | 18.596 | 5.598 |
| L_ancestor | P8C | 39.526 | 57.496 | 76.636 | 17.564 | 5.800 |
| L_ancestor | P8O | 39.355 | 57.786 | 76.606 | 17.811 | 5.583 |

For inclusive tau-origin photons, the W-minus-P8C / W-minus-P8O changes in
`photon` and `association_unmatched` are respectively:

- G: `(+4.288, -3.861)` / `(+4.399, -4.073)` pp.
- L_direct: `(-0.620, +0.631)` / `(-0.720, +0.785)` pp.
- L_ancestor: `(-0.830, +1.032)` / `(-0.800, +0.785)` pp.

The largest W/P8 difference among physical off-diagonal photon outcomes is
0.246 pp for G, 0.053 pp for L_direct and 0.127 pp for L_ancestor. Thus the
remaining G difference is still almost entirely photon versus unmatched, not
a migration into another physical PID.

### Photon rows: fiducial tau-origin

| Definition | Sample | photon | unmatched | fails reco selection | other retained outcomes | selected-reco coverage |
|---|---:|---:|---:|---:|---:|---:|
| G | W | 62.943 | 20.845 | 13.660 | 2.552 | 65.496 |
| G | P8C | 58.622 | 24.705 | 13.738 | 2.936 | 61.557 |
| G | P8O | 58.587 | 24.918 | 13.645 | 2.850 | 61.437 |
| L_direct | W | 56.905 | 30.351 | 12.302 | 0.443 | 57.306 |
| L_direct | P8C | 57.308 | 29.720 | 12.514 | 0.458 | 57.723 |
| L_direct | P8O | 57.517 | 29.566 | 12.431 | 0.486 | 57.951 |
| L_ancestor | W | 62.650 | 18.596 | 14.872 | 3.882 | 66.263 |
| L_ancestor | P8C | 63.229 | 17.564 | 15.166 | 4.041 | 66.979 |
| L_ancestor | P8O | 63.350 | 17.811 | 14.984 | 3.854 | 66.889 |

The L_ancestor selected-reco coverages reproduce the consolidated anchors
66.263%, 66.979% and 66.889% for W, P8C and P8O. Fiducial selection changes
assigned outcomes into `assigned_but_fails_reco_selection`; unmatched and
ambiguity remain separate.

### Charged-pion rows

Every selected charged pion in W, P8C and P8O has a stored tau ancestor.
Consequently each tau-origin charged-pion row is exactly identical to its
all-truth counterpart; the selection cannot further reduce the sample
differences. Full inclusive and fiducial rows, including every physical PID,
unmatched, ambiguity and failed-reco outcome, are retained in
`tau_origin/tables/charged_pion_migration_rows_with_tau_origin.csv`.

The largest W/P8 cell differences are 1.450/1.442 pp for inclusive G,
0.766/0.716 pp for inclusive L_direct, and 0.771/0.747 pp for inclusive
L_ancestor (W-P8C/W-P8O). In fiducial matrices the corresponding maxima are
1.277/1.202, 0.582/0.488 and 0.584/0.466 pp. P8O-P8C differs by at most
0.180 pp over all pion outcomes and definitions.

### Interpretation

1. Restricting to stored tau ancestry removes most of the global W/P8 photon
   discrepancy. The inclusive photon-diagonal difference drops from about
   14.2-14.4 pp to 4.3-4.4 pp for G, from 10.6-10.7 pp to 0.6-0.7 pp in
   magnitude for L_direct, and from 11.6-11.8 pp to 0.8-0.83 pp in magnitude
   for L_ancestor.
2. The conclusion is robust across G and L_ancestor, although G retains a
   visible approximately 4.3 pp photon-versus-unmatched residual while
   L_ancestor is sample-independent at the approximately 1 pp level.
3. L_direct behaves differently by making W very slightly less photon-like
   and more unmatched than P8 after the tau-origin restriction; the sign is
   reversed but the magnitude is below 0.8 pp.
4. P8C and P8O remain extremely close: the largest photon-outcome difference
   is 0.248 pp. Charged-pion matrices are already nearly sample-independent,
   but tau-origin does not change them because their tau-origin fraction is
   exactly 100% in all samples.

### Tau-origin outputs

- Root: `artifacts/migration_matrices_W_P8_G_Ldirect_Lancestor_20260828/tau_origin/`
- Matrices: 9 inclusive and 9 fiducial, each in PNG and PDF.
- Contact sheets: 2, each in PNG and PDF.
- Difference matrices: 30, each in PNG and PDF, using the frozen 11 pp and
  15 pp symmetric scales for truth-definition and sample effects.
- Machine-readable tables and validations are under `tau_origin/tables/` and
  `tau_origin/validation/`.

No linker, HitAnalysis, reconstruction, simulation or generator step was run.
No frozen association, PID, deduplication or fiducial definition was changed.

## Storage migration validation

- The complete 110-file, 12,157,594-byte derived-output tree was copied from its former Lustre location to the artifact root above.
- Relative paths, byte sizes and SHA-256 checksums were compared deterministically before removal: **0 mismatches**.
- No plot, table or scientific result was regenerated or changed; only active storage-path metadata was updated after the exact-copy check.
- The former derived-output tree on Lustre was removed only after the exact-copy validation passed. Frozen Lustre inputs and shared sample/REC/SIM paths were not changed.

No linker, HitAnalysis, reconstruction, simulation or generator step was run. No frozen scientific definition was changed.
