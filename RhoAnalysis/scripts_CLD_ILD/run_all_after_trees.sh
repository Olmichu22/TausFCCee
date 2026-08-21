#!/bin/bash
# Master driver for the CLD vs ILD polarization study, stages 2 → 6.
# Assumes 01_build_trees.sh has already produced the four MDecs trees.
#
#   2  histograms without cuts (reco + gen, both detectors)
#   3  PSO cut optimization per detector and channel
#   3b optimized-cut pipeline configs
#   2b histograms with the optimized cuts (reco only)
#   4  cos(theta) binning + A_tau fit (correlated weights), both variants
#   5  CLD/ILD/gen comparison figures
#   6  summary tables
cd /nfs/cms/arqolmo/TausFCCee
set +u; source setupKey4Hep.sh; set -u
LOGDIR=logs/CLD_ILD; mkdir -p "$LOGDIR"
step () { echo; echo "############ $* — $(date) ############"; }

step "2 — histograms (no cuts)"
bash RhoAnalysis/scripts_CLD_ILD/02_histograms_nocuts.sh

step "4 — binning + fit (no cuts)"
python RhoAnalysis/scripts_CLD_ILD/04_binning_and_fit.py --variant nocuts \
    > "$LOGDIR/04_fit_nocuts.log" 2>&1
tail -5 "$LOGDIR/04_fit_nocuts.log"

step "5 — comparison figures (no cuts)"
python RhoAnalysis/scripts_CLD_ILD/05_compare_plots.py --variant nocuts \
    > "$LOGDIR/05_plots_nocuts.log" 2>&1
tail -5 "$LOGDIR/05_plots_nocuts.log"

step "3 — PSO cut optimization"
bash RhoAnalysis/scripts_CLD_ILD/03_optimize_cuts.sh

step "3b — optimized-cut configs"
python RhoAnalysis/scripts_CLD_ILD/03b_make_optcut_configs.py

step "2b — histograms (optimized cuts)"
for det in CLD ILD; do
  cfg="config/pipeline/CLD_ILD/${det}_reco_optcuts.yaml"
  [ -f "$cfg" ] || { echo "[skip] $cfg not found"; continue; }
  python RhoAnalysis/runTreeHistPipeline_MDecs.py --pipeline-config "$cfg" \
      > "$LOGDIR/hist_${det}_optcuts.log" 2>&1
  echo "hist $det optcuts rc=$?"
done

step "4 — binning + fit (optimized cuts)"
python RhoAnalysis/scripts_CLD_ILD/04_binning_and_fit.py --variant optcuts \
    --detectors CLD ILD > "$LOGDIR/04_fit_optcuts.log" 2>&1
tail -5 "$LOGDIR/04_fit_optcuts.log"

step "5 — comparison figures (optimized cuts)"
python RhoAnalysis/scripts_CLD_ILD/05_compare_plots.py --variant optcuts \
    > "$LOGDIR/05_plots_optcuts.log" 2>&1
tail -5 "$LOGDIR/05_plots_optcuts.log"

step "6 — summary"
python RhoAnalysis/scripts_CLD_ILD/06_summary.py --variants nocuts optcuts

echo; echo "############ ALL STAGES DONE — $(date) ############"
