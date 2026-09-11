# Maintained FCC analysis recipes

Definitions referenced here are frozen in
[Scientific definitions](SCIENTIFIC_DEFINITIONS.md). Workflow products are
data inputs governed by `fcc_tau_association_v1`; TausFCCee does not implement
L_direct or L_ancestor.

The frozen-definition parameterized W/P8O and W/KKMCee comparison suite,
including exact preflight and regression commands, is documented in
[`MC_COMPARISON_PIPELINE.md`](MC_COMPARISON_PIPELINE.md).

## Six-species performance report

**Purpose.** Report maintained L_ancestor association efficiency and residuals
for electron, muon, photon, charged pion, charged kaon, and K0L, while retaining
all other selected PDGs in `other_selected_truth` inventory rows.

**Command.** Select a configured comparison and a new output root:

```bash
python scripts/analysis/build_mc_comparison.py \
  --comparison whizard_p8h_stable_performance \
  --output-root /path/to/new/performance-report \
  --families performance
```

The configured `truth_p_min`, `truth_theta_min_deg`, and
`truth_theta_max_deg` defaults are all `null`, so no extra acceptance is applied.
Optional CLI values are inclusive truth-side cuts and never use reconstructed
kinematics. Outputs include integrated efficiency, differential efficiency
versus truth p/theta/phi, the five maintained residuals, binning/configuration
provenance, and supported reconstructed-PID accounting. Zero-denominator and
underflow/overflow bins remain machine-readable without fake plotted points.

Charged-kaon association, efficiency, and residuals are reported as truth
performance. Kaon PID performance is not evaluated because the maintained
reconstructed-PID categorization has no dedicated kaon category. Photon
residuals describe the maintained representative PFO and are not
unconditionally intrinsic ECAL resolution.

## A. Basic G / HitAnalysis

**Purpose.** Produce geometric selected-truth association (G), efficiency,
confusion, and resolution outputs from REC files.

**Inputs.** A `fcc_hit_analysis_input_manifest_v1` YAML containing non-empty
REC paths in deterministic order; writable, non-existing output target.

**Command.**

```bash
export INPUT_MANIFEST=/path/to/hit_inputs.yaml
export OUTPUT_ROOT=/path/to/new/writable/analysis-output
scripts/analysis/run_hit_analysis.sh \
  --input-manifest "$INPUT_MANIFEST" --output-root "$OUTPUT_ROOT" \
  --workers 8 --prefix MY_SAMPLE_ --dedup-mode reco \
  --assoc-max-dr 0.1 --dry-run
```

After review, remove `--dry-run` to process data. The resulting command performs
real HitAnalysis processing and writes outputs; verify the manifest and ensure
the target directory does not exist first.

**Outputs.** Standard HitAnalysis Parquet/plots plus `run_metadata.json` under
the explicit output root.

**Validation.** Dry-run input count; no worker exceptions; unique composite
event keys; metadata manifest hash and G contract; output files readable.

**Common issues.** Missing REC, excessive worker count, pre-existing output,
confusing `AssocMaxDR=0.1` with tau `dRMax=0.4`, or changing dedup from the
validated `reco` mode. See [Troubleshooting](TROUBLESHOOTING.md).

## B. G versus L_direct versus L_ancestor

**Purpose.** Compare geometric association, immediate MC contribution, and
selected generator-level ancestry without treating them as interchangeable.

**Inputs.** G Parquet output, workflow L_direct/L_ancestor products, source REC
provenance, and a `fcc_tau_workflow_product_manifest_v1` manifest.

**Command.** First validate the interface and manifest:

```bash
export WORKFLOW_MANIFEST=/path/to/workflow_products.yaml
python scripts/validation/validate_fcc_contracts.py
python -c 'import sys; from modules.fcc_workflow_interface import load_product_manifest; print(load_product_manifest(sys.argv[1]))' \
  "$WORKFLOW_MANIFEST"
```

Then provide the campaign-level product paths through the migration config in
Recipe C. There is deliberately no `--workflow-manifest` option on the matrix
script; its actual CLI is `--config`.

**Outputs.** Joined, truth-normalized comparison matrices and validations from
Recipe C.

**Validation.** Contract versions, product columns, source/event/PFO joins,
common denominators, and frozen anchor tables must pass.

**Common issues.** Using derived event integers as identity, joining the wrong
reconstruction provenance, or mistaking a candidate relation for a final
L_direct assignment.

## C. Reproduce the frozen W/P8 migration-matrix study

**Purpose.** Reproduce the frozen inclusive and fiducial
G/L_direct/L_ancestor migration-matrix study for W, P8C, and P8O, including its
numerical anchors. This maintained recipe is not a generic arbitrary-campaign
API; a future generic tool would require a separate interface and validation.

**Inputs.** A `fcc_migration_matrix_inputs_v1` YAML with G products, workflow
campaign CSV manifests, fiducial and L_ancestor validation tables, historical
provenance sources, a new output root, and report path. Expected inventory is
W=996, P8C=100, P8O=18 files.

**Command.**

```bash
export MATRIX_CONFIG=/path/to/frozen_W_P8_migration_inputs.yaml
python scripts/analysis/produce_migration_matrices.py --help
python scripts/analysis/produce_migration_matrices.py --config "$MATRIX_CONFIG"
```

The processing command writes matrices, tables, plots, and a report. It refuses
an existing output/report; verify every frozen input and target path first.

**Outputs.** Absolute inclusive/fiducial matrices, difference matrices,
contact sheets, CSV tables, validation tables, `summary.json`, and a report.

**Validation.** Every truth row sums to 100%; original inclusive G matrices and
frozen fiducial/L_ancestor anchors reproduce exactly; all three definitions use
the same selected-truth denominators.

**Common issues.** Incomplete campaign counts, wrong manifest columns,
non-contiguous G source groups, missing representative PFO joins, or existing
output roots.

## D. Tau-origin migration extension

**Purpose.** Extend an already validated matrix set to selected truth whose
recursive stored ancestry contains a tau.

**Inputs.** The same matrix config as Recipe C, its completed base output, REC
genealogy, and compact frozen truth flags.

**Command.**

```bash
python scripts/analysis/extend_migration_matrices_tau_origin.py --help
python scripts/analysis/extend_migration_matrices_tau_origin.py \
  --config "$MATRIX_CONFIG" --workers 8
```

This command performs real postprocessing and writes an extension below the
configured base output. Verify the completed Recipe-C output and all REC paths
before running it.

**Outputs.** `tau_origin/` matrices, difference tables, ancestry validation,
and summary below the configured base output.

**Validation.** Stored and recursively derived tau flags agree; ancestry has
no cycles or invalid references; all-truth matrices reproduce the base tables;
tau-origin rows sum to 100%.

**Common issues.** Running before Recipe C, missing REC parent branches,
incompatible compact-truth columns, excessive workers, or a residual partial
output directory.

## E. Four-sample truth/bookkeeping audit

**Purpose.** Apply one unweighted audit implementation to W, P8C, P8O, and
KKMCee without collapsing the two PYTHIA8 reconstruction chains.  The catalog
is `configs/analysis/fcc_sample_catalog_v1.yaml`.  Its W row deliberately uses
the frozen Talk-2 W18k nine-file scope; P8C and P8O use their complete 100k and
18k scopes, respectively, and KKMCee expects the final 2k linked REC plus
workflow-owned assignments.

```bash
export FCC_TAU_KKMCEE_MATERIAL=/path/to/kkmcee_material
python scripts/analysis/audit_fcc_samples.py \
  --catalog configs/analysis/fcc_sample_catalog_v1.yaml \
  --audit W P8C P8O KKMCee --workers 8 \
  --output-root "$FCC_TAU_KKMCEE_MATERIAL/audit"
```

The target must not exist. The command reads REC truth/PFO collections and
already-produced L_direct/L_ancestor tables; it does not rebuild an
association. It aborts before scanning when any configured prerequisite is
missing or empty. Outputs include the wide CSV/JSON/Markdown audit plus
machine-readable generator-status, assignment-status, ancestor-depth,
parentless-species, and input-provenance tables.

The frozen W/P8O PFO-coverage regression is explicit and read-only:

```bash
python scripts/validation/validate_fcc_audit_regression.py \
  --audit-csv "$FCC_TAU_KKMCEE_MATERIAL/audit/truth_bookkeeping_audit.csv" \
  --coverage-csv /path/to/frozen/W_vs_PYTHIA8_truth_association_coverage.csv \
  --p8c-outcomes-csv /path/to/frozen/P8C_association_outcome_summary.csv \
  --output-csv "$FCC_TAU_KKMCEE_MATERIAL/validation/audit_regression.csv"
```

For the planned equal-statistics WHIZARD/KKMCee comparison, the catalog pins
W source `000242385`, the first deterministic source in the frozen Talk-2
manifest, and the final KKMCee source ID `700000001`. Both represent exactly
2,000 events. This catalog establishes inputs and provenance only; it does not
claim that the archived Talk-2 plot builders are a maintained generic
comparison API.

## F. Inspect one stored MC event record

`dump_mc_event_record.py` renders the complete stored `MCParticles` genealogy
and `PandoraPFOs` together with existing L_direct and L_ancestor assignment
rows. It never regenerates an association. Supply one REC, its two assignment
Parquets, the stable source ID, and a new output directory. For the stable
2026-04-08 products, use the matching environment:

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
python scripts/analysis/dump_mc_event_record.py \
  --rec "$REC" --ldirect "$LDIRECT" --lancestor "$LANCESTOR" \
  --source-file-id "$SOURCE_ID" --first-events 10 \
  --output-dir /path/to/new/event-records
```

Select one exact event with `--event-key "$SOURCE_ID:7"` (or use
`--events 3,7,14` and ranges such as `--events 100-109`). A diagnostic search
for photon pairs with `10 <= |delta theta| <= 20 mrad` is:

```bash
python scripts/analysis/dump_mc_event_record.py \
  --rec "$REC" --ldirect "$LDIRECT" --lancestor "$LANCESTOR" \
  --source-file-id "$SOURCE_ID" --find-pdg 22 \
  --find-method L_ancestor --abs-theta-residual 10:20 --max-events 20 \
  --output-dir /path/to/new/photon-event-records
```

Each event produces `.detail.txt`, `.tree.txt`, and JSON. `MC#37` and `PFO#12`
mean the actual `MCParticles[37]` and `PandoraPFOs[12]` collection indices;
the JSON also records the fully qualified source/event/object identities. Every
persisted PFO association remains visible. `[REP]` marks the unique PFO chosen
by the maintained underlying-direct T/C representative rule; a terminal tie is
shown as `ambiguous_multiple_pfo` and is not resolved arbitrarily. Residual
searches use ANY associated PFO by default. Add `--representative-only` to use
only that unique representative, as in the maintained Part3b-style reduction.
`[REP]` resolves MC-to-one-PFO reduction for one-truth/one-PFO observables; it
does not change PFO ownership. Several PFOs may each be uniquely assigned to
the same MC, and the inspector continues to display all of them. `[MATCH]`,
when present, is display metadata only.

MC production vertices/endpoints and neutral single-cluster directions from the
nominal IP and MC production vertex are diagnostic metadata only. They do not
replace the official PFO-minus-truth residual or alter search/association. A
parent endpoint/daughter vertex is marked as a local coordinate match only when
its stored coordinates agree within `1e-6 mm`, a numerical equality tolerance,
not a physics distance cut. Multiple clusters are listed and never silently
reduced. The identical inspector implementation applies to P8H and W. The
frozen `selected_truth_v1` requirement is `generatorStatus == 1`,
`p >= 1e-10 GeV`, and `abs(PDG) not in {12,14,16}`.

The v2 inspector supports persisted L_direct and L_ancestor only: G and
presentation filters such as hadronization collapsing are intentionally
unsupported. The genealogy
is shown exactly as stored, including shared nodes, disconnected components,
and explicitly marked cycles.

### Stable photon-theta diagnostic examples

The following are event-level facts from existing stable REC and assignment
products, not new physics classifications. For P8H source `910000000`, event
1, photon `MC#15` has three L_ancestor-associated PFOs: `PFO#1` has
`T=1000`, `C=995`, and `dtheta=+1.213303 mrad` and is `[REP]`; `PFO#3` has
`T=0`, `C=946`, and `+10.037690 mrad`; `PFO#4` has `T=0`, `C=1000`, and
`+23.178000 mrad`. Thus this truth photon passes an ANY-PFO `>10 mrad` test
but does not contribute to the representative-PFO `>10 mrad` tail.

In P8H event 5, photons `MC#27` and `MC#28` share the stored production vertex
`(-0.272001,-3.617490,-34.312023) mm`. Their representative pairs have:

| Pair | official dtheta [mrad] | cluster-IP [mrad] | cluster-from-MC-vertex [mrad] | PFO-vs-cluster-IP [mrad] |
|---|---:|---:|---:|---:|
| `MC#27/PFO#3` | +17.122258 | +16.849890 | -0.516398 | +0.272369 |
| `MC#28/PFO#2` | +17.265029 | +16.945914 | -1.460408 | +0.319115 |

For this event, the large official residual is dominated by comparing the truth
momentum direction at a displaced production vertex with a neutral
reconstructed direction approximately pointing from the nominal IP to the
cluster. This one-event observation must not be generalized into an intrinsic
ECAL angular-resolution statement. The maintained photon-theta observable is
the official residual of the maintained representative PFO associated with a
selected-truth photon; it can contain geometry, association, and topology
effects.

A capped W stable search found representative photon residuals in the same
10--20 mrad interval in events `262, 288, 371, 394, 448` and then stopped; this
is not a population count. Event 394 is largely corrected by pointing from its
displaced MC vertex, whereas event 371 retains about 16 mrad after that
diagnostic correction. Other inspected cases include mixed topology, ancestry,
and representative-PFO effects. These observations are not automatic classifier
labels and do not assign physics-process names to stored genealogy edges.

### P8H primary-vertex provenance

The standalone `out_0.hepmc` record already contains longitudinally displaced
primary vertices before detector simulation; examples include approximately
`z = -20.01, -49.15, -11.30, +11.09, -33.96, +50.78, -64.37, +83.47 mm`. REC
MC primary vertices preserve this scale up to later transformations and decays,
so Geant4/Pandora did not create the large primary-z displacement.

The Pythia card enables `Beams:allowVertexSpread` with
`Beams:sigmaVertexX=5.96e-3`, `Beams:sigmaVertexY=23.8e-6`,
`Beams:sigmaVertexZ=0.397`, and `Beams:sigmaTime=10.89`; Pythia spatial values
are in mm and time is in mm/c. Separately, the k4Gen example configuration uses
`GaussSmearVertex` with `xVertexSigma=yVertexSigma=0.5 mm`,
`zVertexSigma=40.0 mm`, and `tVertexSigma=180 ps`. The observed O(40 mm) z
spread is strong evidence for an additional large k4Gen vertex smearing on top
of Pythia's internal spread, rather than `sigmaVertexZ=0.397 mm` alone. This is
a provenance/diagnostic observation: the generation configuration is unchanged,
and no k4Gen time-unit bug is claimed without validation of the exact production
version and configuration.
