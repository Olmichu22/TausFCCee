#!/bin/bash
# Stage 2 (no cuts) of the CLD vs ILD study: MDecs histograms for the six
# hadron+lepton channels, for the reco and gen trees of both detectors.
cd /nfs/cms/arqolmo/TausFCCee
set +u; source setupKey4Hep.sh; set -u
LOGDIR=logs/CLD_ILD; mkdir -p "$LOGDIR"

for cfg in CLD_reco ILD_reco CLD_gen ILD_gen; do
  echo "===== [hist $cfg] START $(date) ====="
  python RhoAnalysis/runTreeHistPipeline_MDecs.py \
      --pipeline-config "config/pipeline/CLD_ILD/${cfg}_nocuts.yaml" \
      > "$LOGDIR/hist_${cfg}_nocuts.log" 2>&1
  echo "===== [hist $cfg] END rc=$? $(date) ====="
done
echo "HISTOGRAMS (nocuts) DONE $(date)"
