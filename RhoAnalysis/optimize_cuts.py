# optimize_cuts.py
from __future__ import annotations

import argparse
import os
from dataclasses import asdict
import numpy as np

from modules.optimize_pso.data import BranchMap, RootSample, load_samples, load_dataset_weights
from modules.optimize_pso.metrics import LossConfig
from modules.optimize_pso.objective import CutObjective, CutParams
from modules.optimize_pso.pso import PSOConfig, pso_minimize
from modules.optimize_pso.plots import plot_variable_with_cuts, plot_loss_history
import matplotlib.pyplot as plt
import pandas as pd

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signal-root",
                    default="Results/RhoAnalysis/tau_trained0.4_tph0.35_tpi0_n3_g0.0/tau_traineddecay2_0.4_tph0.35_tpi0_n3_g0.0.root"
                    , help="ROOT con la señal (contiene señal+migraciones)")
    ap.add_argument("--bg-root", default=[], nargs="*",
                    help="ROOT de fondo externo (0 o más). Sin ellos el fondo son "
                         "solo las migraciones dentro del fichero de señal.")
    ap.add_argument("--bg-names", default=None, nargs="*",
                    help="Etiquetas para --bg-root (default: BG0, BG1, …)")
    ap.add_argument("--tree", default="outtree_original")
    ap.add_argument("--tauPcut", type=float, default=5, help="Corte fijo mínimo de recoMesonP (tauPCut)")
    ap.add_argument("--selectGEN", type=int, default=2, help="genTauID considerado señal (selectGEN)")

    # bounds (ajústalos a tus rangos físicos)
    ap.add_argument("--dR-bounds", type=float, nargs=2, default=[0.0, 6.0])
    ap.add_argument("--mesonP-bounds", type=float, nargs=2, default=[0.0, 100.0])
    ap.add_argument("--lepP-bounds", type=float, nargs=2, default=[0.0, 100.0])
    ap.add_argument("--Zmass-bounds", type=float, nargs=2, default=[0.0, 100.0])
    # Loss config
    ap.add_argument("--eff-target", type=float, default=0.90)
    ap.add_argument("--eff-lambda", type=float, default=5.0)
    ap.add_argument("--use-s-over-b", action="store_true", help="Si se activa: score = S/(B+eps). Si no: S/sqrt(S+B)")
    ap.add_argument(
    "--compare-yaml",
    default=None,
    help="YAML con luminosidad, xsec y n_events por dataset (sin él, peso 1)"
    )
    ap.add_argument("--outdir", default="validation_plots_d0/",
                    help="Directorio de salida de los plots de validación y del CSV")
    # PSO
    ap.add_argument("--particles", type=int, default=500)
    ap.add_argument("--iters", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument(
    "-v","--verbose",
    action="store_true",
    help="Print PSO progress"
)
    ap.add_argument("--use-styler", action="store_true", help="Usar styler de pandas para output de dataframes")

    return ap.parse_args()


def main():
    args = parse_args()

    # Un fichero de señal (señal + migraciones) + N ficheros de fondo externo.
    # Los nombres legacy Zee/Zqq/Bhabha se mantienen si se pasan 3 fondos sin
    # --bg-names, para no romper los comandos antiguos.
    bg_names = args.bg_names
    if bg_names is None:
        bg_names = (["Zee", "Zqq", "Bhabha"] if len(args.bg_root) == 3
                    else [f"BG{i}" for i in range(len(args.bg_root))])
    if len(bg_names) != len(args.bg_root):
        raise ValueError("--bg-names debe tener tantos elementos como --bg-root")
    samples = [RootSample(name="signal", path=args.signal_root, tree=args.tree,
                          is_signal_file=True)]
    samples += [RootSample(name=n, path=p, tree=args.tree, is_signal_file=False)
                for n, p in zip(bg_names, args.bg_root)]

    if args.use_styler:
        plt.style.use("/nfs/cms/arqolmo/SDHCAL_Energy/utils/newams.mplstyle")
    branches = BranchMap(
        # si tienes una rama de pesos (p.ej. "weight" o "genWeight") ponla aquí:
        weight=None
    )
    if args.compare_yaml:
        dataset_weights = load_dataset_weights(args.compare_yaml, args.selectGEN)
    else:
        # Sin YAML de luminosidad todos los sucesos pesan 1. OJO: el score
        # S/sqrt(S+B) escala como sqrt(w), así que el peso relativo de la
        # penalización de eficiencia (eff_lambda) cambia con la normalización;
        # para un resultado escalado a luminosidad, pasa --compare-yaml.
        dataset_weights = {s.path: 1.0 for s in samples}
    loaded = load_samples(samples, branches=branches, tauP_min_fixed=args.tauPcut, dataset_weights=dataset_weights)

    loss_cfg = LossConfig(
        eff_target=float(args.eff_target),
        eff_lambda=float(args.eff_lambda),
        use_s_over_b=bool(args.use_s_over_b),
    )
    selectGEN = args.selectGEN
    if selectGEN == 2:
        selectGEN = 1
        
    obj = CutObjective(
        loaded=loaded,
        branches=branches,
        selectGEN=int(selectGEN),
        loss_cfg=loss_cfg,
    )

    # Parámetros optimizados: [dR_min, dR_max, mesonP_min, mesonP_max, lepP_min, lepP_max]
    # Ojo: bounds de min/max por separado. Si el PSO propone min>max, objective lo penaliza.
    dR_lo, dR_hi = map(float, args.dR_bounds)
    m_lo, m_hi = map(float, args.mesonP_bounds)
    l_lo, l_hi = map(float, args.lepP_bounds)
    zm_lo, zm_hi = map(float, args.Zmass_bounds)

    bounds = [
        (dR_lo, dR_hi),  # dR_min
        (dR_lo, dR_hi),  # dR_max
        (m_lo, m_hi),    # mesonP_min
        (m_lo, m_hi),    # mesonP_max
        (l_lo, l_hi),    # lepP_min
        (l_lo, l_hi),    # lepP_max
        (zm_lo, zm_hi),  # Zmass_min
        (zm_lo, zm_hi),  # Zmass_max
    ]

    def f(x: np.ndarray) -> float:
        p = CutParams(
            dR_min=float(x[0]),
            dR_max=float(x[1]),
            mesonP_min=float(x[2]),
            mesonP_max=float(x[3]),
            lepP_min=float(x[4]),
            lepP_max=float(x[5]),
            zmass_min=float(x[6]),
            zmass_max=float(x[7]),
        )
        return obj.evaluate(p, modify_mask_inplace=False).loss

    pso_cfg = PSOConfig(
        n_particles=int(args.particles),
        n_iters=int(args.iters),
        seed=int(args.seed),
        verbose=bool(args.verbose),
        patience=30,
        cognitive=2,
        social=1,
        # velocity_clamp=0.2,  # opcional
    )

    res = pso_minimize(f, bounds=bounds, cfg=pso_cfg)

    best = CutParams(
        dR_min=float(res.best_x[0]),
        dR_max=float(res.best_x[1]),
        mesonP_min=float(res.best_x[2]),
        mesonP_max=float(res.best_x[3]),
        lepP_min=float(res.best_x[4]),
        lepP_max=float(res.best_x[5]),
        zmass_min=float(res.best_x[6]),
        zmass_max=float(res.best_x[7]),
        
    )
    best_eval = obj.evaluate(best)

    print("\n=== BEST CUTS ===")
    print(asdict(best))
    print("\n=== METRICS ===")
    print({
        "loss": best_eval.loss,
        "S": best_eval.S,
        "B": best_eval.B,
        "effS": best_eval.effS,
        **best_eval.details
    })
    print("\n=== PSO ===")
    print({"best_f": res.best_f, "n_steps": len(res.history_best)})
    print("\nPlotting validation histograms...")
    outdir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"
    os.makedirs(outdir, exist_ok=True)
    plot_variable_with_cuts(
        loaded, branches, selectGEN, best,
        var="dR", bins=60, outpath=outdir
    )

    plot_variable_with_cuts(
        loaded, branches, selectGEN, best,
        var="mesonP", bins=60, outpath=outdir
    )

    plot_variable_with_cuts(
        loaded, branches, selectGEN, best,
        var="lepP", bins=60, outpath=outdir
    )
    plot_variable_with_cuts(
        loaded, branches, selectGEN, best,
        var="ZMass", bins=60, outpath=outdir, range=(0, 100)
    )
    plot_loss_history(
    res.history_best,
    outpath=outdir,
)
    results_df = pd.DataFrame({
        "dR_min": [best.dR_min],
        "dR_max": [best.dR_max],
        "mesonP_min": [best.mesonP_min],
        "mesonP_max": [best.mesonP_max],
        "lepP_min": [best.lepP_min],
        "lepP_max": [best.lepP_max],
        "Zmass_min": [best.zmass_min],
        "Zmass_max": [best.zmass_max],
        "loss": [best_eval.loss],
        "S": [best_eval.S],
        "B": [best_eval.B],
        "effS": [best_eval.effS],
        **{k: [v] for k, v in best_eval.details.items()},
        "n_steps": [len(res.history_best)],
    })
    results_df.to_csv(outdir+"optimization_results.csv", index=False)
    print(f"Results saved to {outdir}optimization_results.csv")
if __name__ == "__main__":
    main()
