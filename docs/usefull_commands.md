```shell
python TauAnalysis/runShardedPlots.py --sample ild_fcc -j 60 -c '/nfs/cms/arqolmo/TausFCCee/config/default/taurecolong_optimal.yaml' --prefix ILD_FCC_2M_OptCutsCLD_ --hist-config '/nfs/cms/arqolmo/TausFCCee/config/histograms/tau_long_results.yml' -v --clean --overwrite
```

```shell
DIR=Results/RhoAnalysis/PolAnalysis_RECO_ILD_fcc_tau_trained0.4_tph0.35_tpi0_n3_g0.0
STEM=zTaum_dRgt3.06_5.0_Dec0Pgt2.0_lt45.57_Dec1Pgt10.0_lt41.16_cosAcc0.95_vismOff_tau_traineddecayAll_0.4_tph0.35_tpi0_n3_g0.0
OUT=Binned_histograms_MDecs/CLD_ILD/ILD_rho_lep_legacycuts

python RhoAnalysis/makeCosBins_MDecs.py \
    --sample-dir "$DIR" --stem "$STEM" \
    --target-decay rho --other-decay lep \
    --weights corr --bg-def ss \
    --signal-type Ztt \
    --signal-ngen 2000000 --signal-lumi-pb 1354.481301 \
    -o "$OUT" -v
```

```shell
BINNED=$OUT/BINED_MDecs_rho_lep_corr_ss.root

# (a) estadística MC de los 2M eventos simulados
python RhoAnalysis/fitPolAssym.py -i "$BINNED" \
    -o "$OUT/fit_MCstat" --bg-mode total --rebin 2 -v \
    --extra-legend "CLD  1.35 fb^{-1} (MC stat)"

# (b) 7 ab^-1 = 1 año de Z-pole FCC-ee
python RhoAnalysis/fitPolAssym.py -i "$BINNED" \
    -o "$OUT/fit_FCCyear" --bg-mode total --rebin 2 -v \
    --lumi-base 1.354481 --lumi-target 7000 \
    --extra-legend "CLD  7 ab^{-1} (1 FCC-ee Z year)"
```