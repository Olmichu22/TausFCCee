#!/bin/bash
# Stage 3 of the CLD vs ILD study: PSO optimization of the selection cuts
# (dR, meson P, lepton P, Z mass) for each detector and each optimized channel
# (pion, rho and the electron; a1 and the muon are out of the analysis for now).
# The optimizer itself is channel-agnostic: it splits signal/background with
# genTauID == --selectGEN, so the electron only needs its own `--other` set.
#
# The MDecs tree is first projected onto the legacy one-hemisphere schema the
# optimizer expects (RhoAnalysis/mdecsTreeToLegacy.py). Signal = target channel
# gen-correct; background = migrations inside the same sample (no external
# Bhabha/Zqq samples exist for ILD, so none are used for either detector).
cd /nfs/cms/arqolmo/TausFCCee
set +u; source setupKey4Hep.sh; set -u

OUTBASE=Results/CutOptimization_CLD_ILD
LOGDIR=logs/CLD_ILD
mkdir -p "$OUTBASE" "$LOGDIR"

CLD_DIR=Results/RhoAnalysis/PolAnalysis_RECO_CLD_ztt2M_tau_trained0.4_tph0.35_tpi0_n3_g0.0
ILD_DIR=Results/RhoAnalysis/PolAnalysis_RECO_ILD_fcc_tau_trained0.4_tph0.35_tpi0_n3_g0.0
TREENAME=TTree_MDecs_0_1_2_10_-11_-13_tau_traineddecayAll_0.4_tph0.35_tpi0_n3_g0.0.root

# channel: name:reco_id:gen_id:other_ids — id reco del hemisferio objetivo, id
# gen que la PSO toma como señal (rho reco 2 → gen 1, remapeado dentro de
# optimize_cuts) y los ids reco admitidos en el otro hemisferio.
# Para los canales hadrónicos el otro lado es leptónico (la muestra de pares
# clásica); para el electrón, que se analiza como canal inclusivo, el otro lado
# es cualquier especie menos el propio electrón, de modo que la asignación de
# slots (objetivo = dec0) sea inequívoca.
CHANNELS="pion:0:0:-11,-13 rho:2:2:-11,-13 ele:-11:-11:0,2,10,-13"

# Argumentos opcionales: nombres de canal a optimizar (por defecto, todos).
#   bash 03_optimize_cuts.sh ele    → relanza solo el electron
SELECT="$*"

for det in CLD ILD; do
  case $det in
    CLD) TREE="$CLD_DIR/$TREENAME";;
    ILD) TREE="$ILD_DIR/$TREENAME";;
  esac
  for ch in $CHANNELS; do
    IFS=: read -r name reco_id gen_id other_ids <<< "$ch"
    if [ -n "$SELECT" ] && [[ " $SELECT " != *" $name "* ]]; then continue; fi
    other_ids=${other_ids//,/ }
    legacy="$OUTBASE/legacy_${det}_${name}.root"
    outdir="$OUTBASE/${det}_${name}"
    echo "===== [$det $name] adapter $(date) ====="
    python RhoAnalysis/mdecsTreeToLegacy.py --tree-file "$TREE" \
        --target "$reco_id" --other $other_ids -o "$legacy" > "$LOGDIR/opt_${det}_${name}.log" 2>&1 || continue
    echo "===== [$det $name] PSO $(date) ====="
    python RhoAnalysis/optimize_cuts.py \
        --signal-root "$legacy" --selectGEN "$gen_id" \
        --tauPcut 0 --eff-target 0.90 \
        --mesonP-bounds 0 50 --lepP-bounds 0 50 --Zmass-bounds 0 100 \
        --particles 200 --iters 300 --outdir "$outdir" -v \
        >> "$LOGDIR/opt_${det}_${name}.log" 2>&1
    echo "===== [$det $name] rc=$? $(date) ====="
    tail -3 "$outdir/optimization_results.csv" 2>/dev/null
  done
done
echo "CUT OPTIMIZATION DONE $(date)"
