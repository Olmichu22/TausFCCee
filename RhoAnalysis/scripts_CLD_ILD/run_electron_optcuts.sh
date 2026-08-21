#!/bin/bash
# Incremental driver: adds the ELECTRON channel to the optimized-cut variant of
# the CLD vs ILD study, reusing the pion/rho PSO results already on disk.
#
#   3   PSO cut optimization, electron only
#   3b  regenerate the optimized-cut pipeline configs (pion + rho + ele)
#   2b  histograms with the optimized cuts, run ele_incl only
#   4   cos(theta) binning + A_tau fit for the ele channel (optcuts)
#   5   comparison figures (all channels, optcuts)
#   6   summary tables (nocuts + optcuts)
cd /nfs/cms/arqolmo/TausFCCee
set +u; source setupKey4Hep.sh; set -u
LOGDIR=logs/CLD_ILD; mkdir -p "$LOGDIR"
step () { echo; echo "############ $* — $(date) ############"; }

step "3 — PSO cut optimization (electron)"
bash RhoAnalysis/scripts_CLD_ILD/03_optimize_cuts.sh ele

step "3b — optimized-cut configs"
python RhoAnalysis/scripts_CLD_ILD/03b_make_optcut_configs.py

step "2b — histograms (optimized cuts, ele_incl)"
for det in CLD ILD; do
  cfg="config/pipeline/CLD_ILD/${det}_reco_optcuts.yaml"
  [ -f "$cfg" ] || { echo "[skip] $cfg not found"; continue; }
  python RhoAnalysis/runTreeHistPipeline_MDecs.py --pipeline-config "$cfg" \
      --runs ele_incl > "$LOGDIR/hist_${det}_optcuts_ele.log" 2>&1
  echo "hist $det optcuts ele rc=$?"
done

step "4 — binning + fit (optimized cuts, ele)"
python RhoAnalysis/scripts_CLD_ILD/04_binning_and_fit.py --variant optcuts \
    --detectors CLD ILD --channels ele > "$LOGDIR/04_fit_optcuts_ele.log" 2>&1
tail -5 "$LOGDIR/04_fit_optcuts_ele.log"

step "5 — comparison figures (optimized cuts)"
python RhoAnalysis/scripts_CLD_ILD/05_compare_plots.py --variant optcuts \
    > "$LOGDIR/05_plots_optcuts.log" 2>&1
tail -5 "$LOGDIR/05_plots_optcuts.log"

step "6 — summary"
python RhoAnalysis/scripts_CLD_ILD/06_summary.py --variants nocuts optcuts

echo; echo "############ ELECTRON OPTCUTS DONE — $(date) ############"
