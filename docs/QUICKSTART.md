# FCC analysis quick start

This path starts from reconstructed EDM4hep REC files and optionally consumes
versioned G/L_direct/L_ancestor products. G is geometric selected-truth
association; L_direct selects an immediate MC contributor; L_ancestor maps it
to nearest unique selected generator-level ancestry. Commands that would process REC data
are marked **NOT EXECUTED IN STAGE 4**.

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
python -m pip install -e .
```

At another site, source an explicitly chosen compatible stack. That is not a
claim of regression validation against the frozen IFIC environment.

## 2. Create an input manifest

Use a stable source ID and an explicit REC path:

```yaml
schema_version: fcc_hit_analysis_input_manifest_v1
sample: MY_SAMPLE
inputs:
  - source_file_id: "000000001"
    path: /data/my_sample_REC.edm4hep.root
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
  --prefix MY_SAMPLE_ \
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
  --prefix MY_SAMPLE_ \
  --dedup-mode reco \
  --assoc-max-dr 0.1
```

**NOT EXECUTED IN STAGE 4.** The wrapper runs the hardened options
`-v --min-energy-cuts 10 --all-plot 11 13 22 211`, refuses the known target
directory if it already exists, and writes `run_metadata.json` with manifest
hash, repository commit, contract versions, worker count, and matching setup.

## 5. Inspect outputs

```bash
find "$OUTPUT_ROOT" -maxdepth 3 -type f -print | sort
python -m json.tool \
  "$OUTPUT_ROOT/MY_SAMPLE_results0.4_tph0.0_tpi0.0_n0.0_g0.0/run_metadata.json"
```

Require worker completion without hidden exceptions, readable Parquet/plots,
unique `(source_file_id,event_in_file)` identities, and metadata matching the
input manifest and requested G settings.

In raw HitAnalysis Parquet, `source_file_id` is the zero-based input-manifest
ordinal. Map it through the manifest to the declared stable source ID before
joining workflow L_direct/L_ancestor products; see
[`SCIENTIFIC_DEFINITIONS.md`](SCIENTIFIC_DEFINITIONS.md#event-identity).

## 6. Load workflow associations optionally

The repository interface is data only. Validate a product manifest without
importing the workflow repository:

```bash
export WORKFLOW_MANIFEST=/path/to/products.yaml
python -c 'import sys; from modules.fcc_workflow_interface import load_product_manifest; print(load_product_manifest(sys.argv[1]))' \
  "$WORKFLOW_MANIFEST"
```

The accepted schema is shown in
`configs/interfaces/example_workflow_product_manifest.yaml`. Direct and
ancestor Parquet schemas can be checked with
`modules.fcc_workflow_interface.read_association_table`.

## 7. Validate a one-event G/L smoke result

For the maintained software-validation fixture, validate the G output and its
data-only join with L_direct/L_ancestor:

```bash
python scripts/validation/validate_hit_analysis_smoke.py \
  --hit-manifest /path/to/hit_analysis_smoke_input.yaml \
  --result-dir /path/to/MY_SAMPLE_results0.4_tph0.0_tpi0.0_n0.0_g0.0 \
  --workflow-manifest /path/to/smoke_products.yaml \
  --golden-g-pfo /path/to/independent_historical_pfo_comparison.parquet \
  --output /path/to/new/g_hit_analysis_smoke_summary.json
```

The output path must not already exist. The validator maps HitAnalysis's raw
zero-based input ordinal through the input manifest before joining the stable
workflow key, rejects duplicate keys, and checks the independent historical
PFO-to-MC mapping. It is a validator only and does not run HitAnalysis.

## 8. Migration matrices

Copy `configs/analysis/example_migration_matrix_inputs.yaml`, replace every
placeholder with the frozen G products, workflow campaign manifests,
validation anchors, provenance paths, and a new output root, then run:

```bash
python scripts/analysis/produce_migration_matrices.py --config "$MATRIX_CONFIG"
```

**NOT EXECUTED IN STAGE 4.** This maintained analysis requires complete W,
P8C, and P8O inputs with their validated file counts; it is not a one-file
smoke command. See [Analysis recipes](ANALYSIS_RECIPES.md).
