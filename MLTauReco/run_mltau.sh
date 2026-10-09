#!/bin/bash
# Runs a command in the ml-tau container with the ParTauDETR code on the path:
# the ml-tau-model submodule (external/ml-tau-model), its ml-tau-data submodule
# (the ntupelizer that built the training data) and this directory.
#
#   MLTauReco/run_mltau.sh python3 MLTauReco/process_edm4hep.py <files> --out-dir <dir>
#
# The model needs torch, lightning, omegaconf and fastjet's Python bindings,
# which the key4hep stack does not provide, hence the container. Set the image
# with MLTAU_CONTAINER and any extra bind mounts (e.g. the data directories)
# with MLTAU_BINDS, a comma-separated apptainer -B list.
#
# Needs the submodules: git submodule update --init --recursive
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
model_repo="$(cd -- "$here/../external/ml-tau-model" && pwd)"
if [[ ! -d "$model_repo/mltau" || ! -d "$model_repo/ml-tau-data/ntupelizer" ]]; then
	echo "run_mltau.sh: submodules missing, run 'git submodule update --init --recursive'" >&2
	exit 1
fi
if [[ -z "${MLTAU_CONTAINER:-}" ]]; then
	echo "run_mltau.sh: set MLTAU_CONTAINER to an apptainer image with torch, lightning, omegaconf and fastjet" >&2
	exit 1
fi
binds=()
[[ -n "${MLTAU_BINDS:-}" ]] && binds=(-B "$MLTAU_BINDS")
pythonpath="$here:$model_repo:$model_repo/mltau:$model_repo/ml-tau-data"
# keras is not used, but something imports it and it crashes without a backend
exec apptainer exec "${binds[@]}" --env PYTHONPATH="$pythonpath" --env KERAS_BACKEND=torch "$MLTAU_CONTAINER" "$@"
