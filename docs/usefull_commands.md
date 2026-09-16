```shell
python TauAnalysis/runShardedPlots.py --sample ild_fcc -j 60 -c '/nfs/cms/arqolmo/TausFCCee/config/default/taurecolong_optimal.yaml' --prefix ILD_FCC_2M_OptCutsCLD_ --hist-config '/nfs/cms/arqolmo/TausFCCee/config/histograms/tau_long_results.yml' -v --clean --overwrite
```

```shell
python HitAnalysis/particle_level_analisis_parallel.py --sample ztt_2M_smearing --prefix CLD_2M_ztt_smearingChanged_ --n-workers 60 --dedup-mode reco --fake-bin-by-reco --min-energy-cuts 10 --all-plot 22 211 13 11
```
```shell
python HitAnalysis/replot_from_parquet.py /nfs/cms/arqolmo/TausFCCee/Results/TauReco/CLD_2M_ztt_smearingChanged_results0.4_tph0.0_tpi0.0_n0.0_g0.0 --fake-bin-by-reco --all-plot 22 211 13 11 --min-energy-cuts 10
```
```shell
python TauAnalysis/TausCompletePlot.py -i /nfs/cms/arqolmo/TausFCCee/Results/TauReco/merged_CLD_FCC_2M_OptCutsCLD_smearingChanged_tau_trained0.4_tph0.35_tpi0_n3_g0.0 -p config/plots/taulong_plotconfig_results.yaml
```


```shell
# Grafos "Not reconstructed" (hEffiGen*ToLost = Unmatched + NoTau + Other).
# No los produce la reco: hay que regenerarlos tras cada produccion, y las
# entradas "* Not reconstructed" de config/plots/compare_CLD_ILD.yaml
# (anclas *cldmerged / *ildmerged) apuntan a este migration_merged.root.
RUN=merged_CLD_FCC_2M_OptCutsCLD_PionPhotonFSR_tau_trained0.4_tph0.35_tpi0_n3_g0.0
python TauAnalysis/mergeMigrationCategories.py \
-i "Results/TauReco/$RUN/tau_traineddecayAll_0.4_tph0.35_tpi0_n3_g0.0.root"
# -> Results/TauReco/$RUN/migration_merged.root  (-o para otra ruta)
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