#!/bin/bash
# Stage 1 of the CLD vs ILD polarization study: build the MDecs trees.
#   * reco trees (analysisRHOTree_MDecs_parallel.py) for both detectors
#   * gen-level trees (genOnlyRHOTree_MDecs_parallel.py) for both detectors
# Channels kept: pi (0), rho (2 reco, 1 = rho con un solo fotón; 1 gen), a1 (10),
#                e (-11), mu (-13).
# Base config: config/default/taurecolong_optimal.yaml (tph 0.35, NeutronCut 3).
cd /nfs/cms/arqolmo/TausFCCee
set +u; source setupKey4Hep.sh; set -u

NW=16
CFG=config/default/taurecolong_optimal.yaml
MODES="0 1 2 10 -11 -13"
LOGDIR=logs/CLD_ILD
mkdir -p "$LOGDIR"

run () {  # run <label> <script> <sample> <prefix> [extra...]
  local label=$1 script=$2 sample=$3 prefix=$4; shift 4
  echo "===== [$label] START $(date) ====="
  python "RhoAnalysis/$script" \
      --sample "$sample" -c "$CFG" --decay-modes $MODES \
      --prefix "$prefix" --n-workers $NW -v "$@" \
      > "$LOGDIR/${label}.log" 2>&1
  local rc=$?
  echo "===== [$label] END rc=$rc $(date) ====="
  return $rc
}

run RECO_CLD analysisRHOTree_MDecs_parallel.py ztt_2M  PolAnalysis_RECO_CLD_ztt2M_
run RECO_ILD analysisRHOTree_MDecs_parallel.py ild_fcc PolAnalysis_RECO_ILD_fcc_
run GEN_CLD  genOnlyRHOTree_MDecs_parallel.py  ztt_2M  PolAnalysis_GEN_CLD_ztt2M_
run GEN_ILD  genOnlyRHOTree_MDecs_parallel.py  ild_fcc PolAnalysis_GEN_ILD_fcc_
echo "ALL DONE $(date)"
