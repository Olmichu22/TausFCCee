#!/bin/bash
# CLD vs ILD con vertex smearing corregido — etapas 1+2 (árbol reco + histogramas MDecs).
# Cada detector pasa por runTreeHistPipeline_MDecs.py (árbol + 14 runs de histogramas:
# 7 modos × sin²θ_eff ∈ {0.2312, 0.2315}). Los detectores van EN SERIE para no
# superar 20 núcleos simultáneos (n_workers = 20 en los YAML).
cd /nfs/cms/arqolmo/TausFCCee
set +u; source setupKey4Hep.sh; set -u

LOGDIR=logs/CLD_ILD_smearingChanged; mkdir -p "$LOGDIR"
for det in ${DETS:-CLD ILD}; do
  echo "===== [$det pipeline] START $(date) ====="
  python RhoAnalysis/runTreeHistPipeline_MDecs.py \
      --pipeline-config "config/pipeline/CLD_ILD_smearingChanged/${det}_reco.yaml" "$@" \
      > "$LOGDIR/pipeline_${det}.log" 2>&1
  echo "===== [$det pipeline] END rc=$? $(date) ====="
done
echo "TREES+HISTS DONE $(date)"
