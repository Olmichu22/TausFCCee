# Maintained FCC analysis recipes

Definitions referenced here are frozen in
[Scientific definitions](SCIENTIFIC_DEFINITIONS.md). Workflow products are
data inputs governed by `fcc_tau_association_v1`; TausFCCee does not implement
L_direct or L_ancestor.

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
