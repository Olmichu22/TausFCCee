# Maintained W-versus-MC comparison pipeline

`scripts/analysis/build_mc_comparison.py` is the maintained parameterized port
of the frozen Talk-2 builders.  Both configured comparisons use the same
scientific contract in `configs/analysis/fcc_mc_comparisons_v1.yaml`; only the
sample identity, paths, scope, presentation label and output root differ.

## Authoritative builder inventory

| Output family | Archived executable authority | Frozen input authority |
|---|---|---|
| Part 1/2 tau modes and terminal-tau kinematics | `audit_tau_topology_equalstats.py` | `tau_records.csv`, `tau_decay_channels.csv` |
| Part 1/2 selected truth and photon kinematics | `audit_tau_topology_equalstats.py`; `audit_pfo_reconstruction_efficiency_equalstats.py` | `truth_species_composition.csv`, `selected_truth_association_rows.parquet` |
| Part 3 G/L_direct/L_ancestor efficiency and binned inefficiency | `audit_pfo_reconstruction_efficiency_equalstats.py` | `selected_truth_association_rows.parquet` |
| Part 3 PFO coverage | `audit_pid_equalstats.py`; frozen assignment rows | G rows and L_direct/L_ancestor PFO assignments |
| Part 3b residuals and photon-origin split | `audit_pid_equalstats.py` | `pid_lancestor_associated_rows.parquet` joined to frozen tau-origin flags |
| Part 4 conditional PID | `audit_pid_equalstats.py` | uniquely associated representative-PFO rows; unmatched and ambiguous rows are excluded |

The authoritative snapshots are under
`archive/ild-tau_pre_consolidation_active_20260902/`.  The frozen-output mapping
is independently recorded by `talk2_material/talk2_plot_manifest.csv`.  The
active implementation ports the executable rules; it does not import runtime
code from the archive.

## Frozen contract

The 2026-09-01 definitions remain unchanged: `selected_truth_v1`, stored-parent
tau ancestry, theta/phi G with strict `d < 0.1` and reco-side deduplication,
`truthlink_assignment_v1`, `truthlink_ancestor_assignment_v1`, and the frozen
track-then-cluster representative-PFO ranking.  PID is conditional only on a
unique association.  Residuals, bins, visible ranges, under/overflow
normalization and omission of zero-missed points from logarithmic curves are
shared constants, not sample configuration.

## W/P8O regression

Generate into a disposable directory; `talk2_material` is read-only:

```bash
python scripts/analysis/build_mc_comparison.py \
  --comparison whizard_p8o_18k \
  --output-root /tmp/fcc_tau_w_p8o_regression \
  --validation-mode
```

Machine-readable products are compared to frozen Talk-2 CSVs.  PNG pixel
identity is deliberately not required.

## Final W/KKMCee 2k production

The W side is the canonical source `000242385` (2,000 events).  The KKMCee side
is the final 2,000-event TruthlinkV1 REC, never test10.  The comparison validates
the already-produced workflow preflight receipt:

`/lustre/ific.uv.es/prj/gl/abehep.flc/FCC/KKMCee/ILD20260908_2k/validation/kkmcee_2k_truthlink_preflight.json`

After the workflow has produced the configured L_direct and L_ancestor parquet
products under `kkmcee_material/derived/`, run:

```bash
cd /lhome/ific/a/airqui/FCC/TausFCCee
export FCC_TAU_KKMCEE_MATERIAL=/lhome/ific/a/airqui/FCC/kkmcee_material
python scripts/analysis/build_mc_comparison.py \
  --comparison whizard_kkmcee_2k \
  --output-root /lhome/ific/a/airqui/FCC/kkmcee_material/comparison_W2k_KKMCee2k
```

The command refuses a missing/empty input, a failed or mismatched 2,000-event
preflight receipt, or missing assignment products.  It calls maintained G,
inverts existing assignments, and never reruns reconstruction, the linker,
L_direct, or L_ancestor.

The four-sample bookkeeping audit remains separate:

```bash
export FCC_TAU_KKMCEE_MATERIAL=/lhome/ific/a/airqui/FCC/kkmcee_material
python scripts/analysis/audit_fcc_samples.py \
  --catalog configs/analysis/fcc_sample_catalog_v1.yaml \
  --output-root /lhome/ific/a/airqui/FCC/kkmcee_material/audit
```

## Stable P8O/P8H 10k comparison

`p8o_p8h_stable10k` compares the controlled 10,000-event samples with an
inclusive selected-truth denominator.  It consumes the existing integrated
TruthLinkV1, L_direct, and L_ancestor products; it does not rebuild them.  To
avoid running the unrequested PID family:

```bash
export FCC_TAU_PYTHIA_20260909_MATERIAL=/lhome/ific/a/airqui/FCC/pythia_20260909_material
python scripts/analysis/build_mc_comparison.py \
  --comparison p8o_p8h_stable10k \
  --output-root "$FCC_TAU_PYTHIA_20260909_MATERIAL/comparison_P8O_P8H_stable10k" \
  --families part3 part3b photon_diagnostic
```

P8O starts from historical SIM, whereas P8H simulation and both reconstruction
chains use stable Key4hep 2026-04-08.  Therefore the comparison is symmetric
from reconstruction onward, not at simulation level.

Current stable diagnostic anchors are approximately `12.25/12.20 mrad` for
the P8O/P8H selected-photon theta central-68 half-widths and `0.84 mrad` for
the W stable control. Moving Key4hep release did not remove the broad P8
component, changing from the historical P8 source to P8H did not remove it, and
the W stable control remained essentially consistent with its historical
nightly result. These are observable-level comparisons, not a claim that the
difference is intrinsic detector or ECAL angular resolution; the event-record
diagnostics show contributions from geometry, association, and topology.
