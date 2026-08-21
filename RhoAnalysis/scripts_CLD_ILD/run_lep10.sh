#!/bin/bash
# Leptonic inclusive channels (e + X, mu + X) with P_target > 10 GeV.
# Stage 2 (histograms) -> makeCosBins -> fitPolAssym, for CLD and ILD.
# Baseline to compare against: the *_nocuts runs (same setup, no momentum cut).
cd /nfs/cms/arqolmo/TausFCCee
set +u; source setupKey4Hep.sh; set -u

LOGDIR=logs/CLD_ILD; mkdir -p "$LOGDIR"
OUTBASE=Binned_histograms_MDecs/CLD_ILD
STEM='zTaum_Dec0Pgt10.0_lt100.0_sineff0.2312_vismOff_tau_traineddecayAll_0.4_tph0.35_tpi0_n3_g0.0'
LUMI_COMMON=1.354481      # 2e6 / 1476.58 pb / 1000  -> 2M taus
LUMI_FCC=7000             # 7 ab^-1, one FCC-ee Z-pole year

# ---------- paso 2: histogramas ----------
for det in CLD ILD; do
  echo "===== [hist $det lep10] START $(date) ====="
  python RhoAnalysis/runTreeHistPipeline_MDecs.py \
      --pipeline-config "config/pipeline/CLD_ILD/${det}_reco_lep10.yaml" \
      > "$LOGDIR/hist_${det}_lep10.log" 2>&1
  echo "===== [hist $det lep10] END rc=$? $(date) ====="
done

# ---------- pasos 4a/4b: binado + fits ----------
# det:sample_dir:ngen:lumi_pb:lumi_fb
DETS="CLD:Results/RhoAnalysis/PolAnalysis_RECO_CLD_ztt2M_tau_trained0.4_tph0.35_tpi0_n3_g0.0:2000000:1354.481301:1.354481
ILD:Results/RhoAnalysis/PolAnalysis_RECO_ILD_fcc_tau_trained0.4_tph0.35_tpi0_n3_g0.0:1918000:1298.947568:1.298948"

for row in $DETS; do
  IFS=: read -r det dir ngen lumi_pb lumi_fb <<< "$row"
  for ch in ele muon; do
    OUT="$OUTBASE/${det}_${ch}_lep10"
    LOG="$LOGDIR/fit_${det}_${ch}_lep10.log"
    : > "$LOG"
    echo "===== [$det $ch lep10] binning $(date) ====="
    python RhoAnalysis/makeCosBins_MDecs.py \
        --sample-dir "$dir" --stem "$STEM" \
        --target-decay "$ch" --other-decay all \
        --weights corr --bg-def rho_only \
        --signal-type Ztt --signal-ngen "$ngen" --signal-lumi-pb "$lumi_pb" \
        -o "$OUT" -v >> "$LOG" 2>&1 || { echo "  [FAIL binning] ver $LOG"; continue; }

    BINNED="$OUT/BINED_MDecs_${ch}_all_corr_rho_only.root"
    echo "===== [$det $ch lep10] fits $(date) ====="
    python RhoAnalysis/fitPolAssym.py -i "$BINNED" -o "$OUT/fit_MCstat" \
        --bg-mode total -v --extra-legend "$det $ch  ${lumi_fb} fb^{-1} (MC stat)" >> "$LOG" 2>&1
    python RhoAnalysis/fitPolAssym.py -i "$BINNED" -o "$OUT/fit_2Mtau" \
        --bg-mode total -v --lumi-base "$lumi_fb" --lumi-target "$LUMI_COMMON" \
        --extra-legend "$det $ch  2M #tau#tau" >> "$LOG" 2>&1
    python RhoAnalysis/fitPolAssym.py -i "$BINNED" -o "$OUT/fit_FCCyear" \
        --bg-mode total -v --lumi-base "$lumi_fb" --lumi-target "$LUMI_FCC" \
        --extra-legend "$det $ch  7 ab^{-1} (1 FCC-ee Z year)" >> "$LOG" 2>&1
    echo "  -> $OUT"
  done
done
echo "LEP10 DONE $(date)"
