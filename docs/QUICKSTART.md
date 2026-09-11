# FCC analysis quick start

This path starts from reconstructed EDM4hep REC files. TausFCCee produces G,
the analysis-side geometric selected-truth association, through HitAnalysis.
It can optionally consume L_direct and L_ancestor products exported by
FCC-tau-workflow; it does not consume a pre-produced G product in the standard
golden path. Commands that process REC data perform real analysis and write
outputs, so verify the input manifest and output root before running them.

## 1. Clone and prepare the environment

```bash
git clone YOUR_TAUSFCCEE_URL TausFCCee
cd TausFCCee
```

Pip installs the local helper package and Python dependencies only when the
appropriate extras are requested. ROOT, podio, EDM4hep, and Key4hep are an
external stack and cannot be installed by pip.

The frozen FCC study used Key4hep `2026-08-21`. IFIC example:

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/releases/2026-08-21/x86_64-almalinux9-gcc14.2.0-opt/key4hep-stack/2026-08-21-5qmpe6/setup.sh
export TAUSFCCEE_DEPENDENCY_ROOT="$PWD/../tausfccee-dependencies"
python -m venv --system-site-packages "$TAUSFCCEE_DEPENDENCY_ROOT/venv"
source "$TAUSFCCEE_DEPENDENCY_ROOT/venv/bin/activate"
python -m pip install --no-build-isolation -e .
```

At another site, source an explicitly chosen compatible stack. That is not a
claim of regression validation against the frozen IFIC environment. The
virtual environment inherits Key4hep's external packages without attempting
to write to the read-only CVMFS Python prefix.

## 2. Create an input manifest

Use a stable source ID and an explicit REC path:

```yaml
schema_version: fcc_hit_analysis_input_manifest_v1
sample: W
inputs:
  - source_file_id: "033851393"
    path: /data/events_033851393_REC.edm4hep.root
```

Save it as `$INPUT_MANIFEST`. The maintained example is
`configs/analysis/example_hit_analysis_manifest.yaml`.

## 3. Validate and dry-run HitAnalysis

```bash
export INPUT_MANIFEST=/path/to/hit_inputs.yaml
export OUTPUT_ROOT=/path/to/writable/analysis-output
scripts/analysis/run_hit_analysis.sh \
  --input-manifest "$INPUT_MANIFEST" \
  --output-root "$OUTPUT_ROOT" \
  --workers 8 \
  --prefix W_ \
  --dedup-mode reco \
  --assoc-max-dr 0.1 \
  --dry-run
```

The dry-run resolves and checks every REC, prints the exact command, and does
not run HitAnalysis. G uses a strict `d<0.1`; the separate tau cone remains
`dRMax=0.4`.

## 4. Run the maintained analysis

Remove only `--dry-run` after reviewing the command and confirming that its
output directory does not exist:

```bash
scripts/analysis/run_hit_analysis.sh \
  --input-manifest "$INPUT_MANIFEST" \
  --output-root "$OUTPUT_ROOT" \
  --workers 8 \
  --prefix W_ \
  --dedup-mode reco \
  --assoc-max-dr 0.1
```

This command performs real HitAnalysis processing and writes outputs. The
wrapper runs the hardened options `-v --min-energy-cuts 10 --all-plot 11 13 22
211`, refuses the known target directory if it already exists, and writes
`run_metadata.json` with manifest hash, repository commit, contract versions,
worker count, and matching setup.

## 5. Inspect outputs

```bash
find "$OUTPUT_ROOT" -maxdepth 3 -type f -print | sort
python -m json.tool \
  "$OUTPUT_ROOT/W_results0.4_tph0.0_tpi0.0_n0.0_g0.0/run_metadata.json"
```

Require worker completion without hidden exceptions, readable Parquet/plots,
unique `(source_file_id,event_in_file)` identities, and metadata matching the
input manifest and requested G settings.

In raw HitAnalysis Parquet, `source_file_id` is the zero-based input-manifest
ordinal. Map it through the manifest to the declared stable source ID before
joining workflow L_direct/L_ancestor products; see
[`SCIENTIFIC_DEFINITIONS.md`](SCIENTIFIC_DEFINITIONS.md#event-identity).

## 6. Load workflow L associations optionally

TausFCCee produces G locally through HitAnalysis. L_direct and L_ancestor are
instead exported by FCC-tau-workflow and consumed here through a data-only
product manifest. Validate that manifest without importing the workflow
repository:

```bash
export WORKFLOW_MANIFEST=/path/to/products.yaml
python -c 'import sys; from modules.fcc_workflow_interface import load_product_manifest; print(load_product_manifest(sys.argv[1]))' \
  "$WORKFLOW_MANIFEST"
```

The authoritative `fcc_tau_workflow_product_manifest_v1` schema and
`fcc_tau_association_v1` contract are owned by FCC-tau-workflow. A consumer-side
example is mirrored in `configs/interfaces/example_workflow_product_manifest.yaml`.
Direct and ancestor Parquet schemas can be checked with
`modules.fcc_workflow_interface.read_association_table`.

## 7. Validate a one-event G/L smoke result

For the maintained software-validation fixture, validate the G output and its
data-only join with L_direct/L_ancestor:

```bash
python scripts/validation/validate_hit_analysis_smoke.py \
  --hit-manifest /path/to/hit_analysis_smoke_input.yaml \
  --result-dir /path/to/W_results0.4_tph0.0_tpi0.0_n0.0_g0.0 \
  --workflow-manifest /path/to/smoke_products.yaml \
  --golden-g-pfo /path/to/independent_historical_pfo_comparison.parquet \
  --output /path/to/new/g_hit_analysis_smoke_summary.json
```

The output path must not already exist. The validator maps HitAnalysis's raw
zero-based input ordinal through the input manifest before joining the stable
workflow key, rejects duplicate keys, and checks the independent historical
PFO-to-MC mapping. It is a validator only and does not run HitAnalysis.

## 8. Reproduce the frozen W/P8 migration-matrix study

Copy `configs/analysis/example_migration_matrix_inputs.yaml`, replace every
placeholder with the frozen G products, workflow campaign manifests,
validation anchors, provenance paths, and a new output root, then run:

```bash
export MATRIX_CONFIG=/path/to/frozen_W_P8_migration_inputs.yaml
python scripts/analysis/produce_migration_matrices.py --config "$MATRIX_CONFIG"
```

This command performs real postprocessing and writes matrices, tables, plots,
and a report. It requires the complete frozen W, P8C, and P8O inputs, their
validated file counts, and the frozen regression anchors; it is not a generic
arbitrary-campaign API or a one-file smoke command. A future generic migration
tool would require a separate interface and validation. See
[Analysis recipes](ANALYSIS_RECIPES.md).
