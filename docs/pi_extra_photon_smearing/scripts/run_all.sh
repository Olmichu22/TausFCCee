#!/bin/bash
# Rehace todo docs/pi_extra_photon para ztt_2M_smearing (Ztt_SM_2M_SigmaZ).
# Mismo analisis; solo cambian el tree de entrada (F=) y los symlinks EDM4hep.
cd "$(dirname "$(readlink -f "$0")")"
source /cvmfs/sw.hsf.org/key4hep/setup.sh >/dev/null 2>&1
python3 -c "import uproot, podio" || { echo "key4hep env no disponible"; exit 1; }

echo "== fase 1: extraccion $(date)"
python3 ptau_origin.py                > ptau_origin.log 2>&1 &
python3 baseline_full.py              > baseline_full.log 2>&1 &
python3 origen_extract.py tables_origen.npz > origen_extract.log 2>&1 &
python3 discrim_extract.py            > discrim_extract.log 2>&1 &
python3 discrim_hard_extract.py       > discrim_hard_extract.log 2>&1 &
python3 cutjust_extract.py            > cutjust_extract.log 2>&1 &
python3 cutjust_extract_mgg.py        > cutjust_extract_mgg.log 2>&1 &
python3 cutjust_extract_allpairs.py   > cutjust_extract_allpairs.log 2>&1 &
python3 cutjust_extract_pi0match.py   > cutjust_extract_pi0match.log 2>&1 &
python3 photon_links.py               > photon_links.log 2>&1 &
( cd /nfs/cms/arqolmo/TausFCCee && python3 "$OLDPWD/validate_corr.py" ) > validate_corr.log 2>&1 &
wait
echo "== fase 2: analisis $(date)"
python3 sec11_tables.py > sec11_tables.log 2>&1
python3 visP_no_separa.py > visP_no_separa.log 2>&1
python3 origen_analyze.py tables_origen.npz > origen_analyze.log 2>&1; mv -f RESULTS.md RESULTS_origen.md
python3 discrim_analysis.py      > discrim_analysis.log 2>&1;      mv -f RESULTS.md RESULTS_discrim.md
python3 discrim_hard_analysis.py > discrim_hard_analysis.log 2>&1
python3 summarize.py             > summarize.log 2>&1;             mv -f RESULTS.md RESULTS_raw.md
for s in cutjust_plot cutjust_plot_mgg cutjust_plot_allpairs cutjust_plot_pi0match cutjust_plot_pi0match_dist cutjust_plot_pi0match_mass; do
  python3 $s.py > $s.log 2>&1
done

echo "== fase 3: figuras $(date)"
F=../figs
cp -f fig_migracion_vs_GenTauP.png fig_visP_no_separa.png $F/
for f in fig_migracion_P_cos fig_cinematica_1d fig_P_vs_dR_2d fig_Pgamma_por_Ptau fig_ratio_por_Ptau fig_balance_momento fig_fsr_gen_vs_reco; do cp -f $f.png $F/origen_$f.png; done
for f in fig_link_weight fig_cluster_distance fig_energy fig_fsr_gen; do cp -f $f.png $F/raw_$f.png; done
for f in fig1_distribuciones fig2_malo_por_origen fig3_roc_foton fig4_roc_tau fig5_eff_vs_P fig6_migracion fig7_masa_tau; do cp -f $f.png $F/discrim_$f.png; done
for f in dist_duro roc_duro eff_vs_GenTauP eff_vs_VisP; do cp -f hard_fig_$f.png $F/discrim_hard_fig_$f.png; done
cp -f cutjust_fig*.png $F/
grep -l -E "Traceback|Error" *.log
echo "== fin $(date)"
