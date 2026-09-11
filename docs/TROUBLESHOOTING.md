# FCC analysis troubleshooting

## REC path or collection problem

- **Symptom:** manifest resolution fails, `podio` cannot open a file, or
  `MCParticles`/`PandoraPFOs` is missing.
- **Cause:** wrong path, empty/corrupt REC, or incompatible reconstruction.
- **Diagnostic:** inspect the manifest, run `test -s "$REC"` and
  `podio-dump -e 0 "$REC"`.
- **Resolution:** correct the manifest or use a compatible validated REC; do
  not rename collections in analysis code.

## Event identity collision or failed join

- **Symptom:** duplicate events, cross-file collisions, or G/L product joins
  fail.
- **Cause:** a fixed events-per-file arithmetic ID replaced the composite key,
  or sample/source provenance differs.
- **Diagnostic:** inspect `(sample,source_file_id,event_in_file)` and object
  indices on both sides.
- **Resolution:** restore the authoritative composite identity and deterministic
  manifest order. Cantor/packed integers are conveniences only.

## Worker count or worker failure

- **Symptom:** memory pressure, slow startup, or incomplete output.
- **Cause:** too many workers for files/memory, or a worker raised an exception.
- **Diagnostic:** compare `--workers` with input count/resources and read the
  first worker traceback.
- **Resolution:** reduce workers or fix the underlying event/schema problem.
  Hardened HitAnalysis propagates worker failures; do not ignore them.

## Output-root issue

- **Symptom:** wrapper refuses an output or writes cannot start.
- **Cause:** path exists, parent is unwritable, or a legacy `Results/` path was
  assumed.
- **Diagnostic:** inspect the printed dry-run command and `test -w` on the
  chosen parent.
- **Resolution:** pass a new explicit `--output-root`; preserve existing output
  for review rather than overwriting it.

## G matching or dedup mismatch

- **Symptom:** counts differ from the frozen setup.
- **Cause:** changed `--assoc-max-dr`, `--dedup-mode`, truth selection, or
  confusion between G and tau-cone radii.
- **Diagnostic:** inspect `run_metadata.json` and
  `configs/analysis/geometric_association_v1.yaml`.
- **Resolution:** use strict threshold 0.1 and `dedup-mode=reco`; keep tau
  `dRMax=0.4` separate.

## Workflow manifest/schema mismatch

- **Symptom:** unsupported schema/contract, duplicate source ID, or missing
  required entry field.
- **Cause:** wrong manifest type or an ad-hoc product list.
- **Diagnostic:** run `python scripts/validation/validate_fcc_contracts.py` and
  call `load_product_manifest` on the file.
- **Resolution:** use `fcc_tau_workflow_product_manifest_v1` with
  `fcc_tau_association_v1`; regenerate from authoritative product metadata.

## L product join failure

- **Symptom:** missing source, event, PFO, representative, or denominator join.
- **Cause:** mixed samples/reconstruction versions, incomplete manifests, or
  mismatched direct/ancestor products.
- **Diagnostic:** compare source basenames/IDs, contract columns, event counts,
  and `(event_in_file,pfo_index)` coverage.
- **Resolution:** pair products from the same source REC and provenance. Do not
  fill missing rows or recompute L inside TausFCCee.

## PID mapping issue

- **Symptom:** unsupported reconstructed PDG or sentinel treated as a particle.
- **Cause:** data outside the frozen mapping or duplicated local mapping logic.
- **Diagnostic:** call `reconstructed_pid_category` and inspect the PFO PDG.
- **Resolution:** use `modules/fcc_truth_definitions.py`; `999` remains an
  unmatched sentinel, never a PID.

## Missing optional scientific artifacts

- **Symptom:** a frozen report references a large table/plot tree not in Git.
- **Cause:** repository policy keeps large machine-readable artifacts external.
- **Diagnostic:** consult `docs/scientific/freezes/20260901/` manifests and
  hashes.
- **Resolution:** obtain the frozen external artifact from its custodian and
  verify its checksum. Do not regenerate it as part of an operational run.

## Matrix configuration or output failure

- **Symptom:** count assertion, anchor mismatch, missing provenance path, or
  existing output root.
- **Cause:** incomplete W/P8C/P8O inputs, wrong config schema, incompatible
  frozen products, or attempted overwrite.
- **Diagnostic:** run `--help`, validate all config paths and inspect the first
  assertion.
- **Resolution:** correct the explicit `fcc_migration_matrix_inputs_v1` config
  or choose a new output root. Never change expected scientific anchors to
  force a pass.
