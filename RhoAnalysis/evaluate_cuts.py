# optimize_cuts.py
from __future__ import annotations

import argparse
from dataclasses import asdict
import numpy as np

from modules.optimize_pso.data import BranchMap, RootSample, load_samples, load_dataset_weights
from modules.optimize_pso.metrics import LossConfig
from modules.optimize_pso.objective import CutObjective, CutParams
from modules.optimize_pso.pso import PSOConfig, pso_minimize
from modules.optimize_pso.plots import plot_variable_with_cuts, plot_loss_history, plot_omega_shape_comparison_from_arrays
import pandas as pd
import os
import matplotlib.pyplot as plt
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signal-root",
                    default="Results/RhoAnalysis/tau_trained0.4_tph0.35_tpi0_n3_g0.0/tau_traineddecay2_0.4_tph0.35_tpi0_n3_g0.0.root"
                    , help="ROOT con la señal (contiene señal+migraciones)")
    ap.add_argument("--bg-root",default=["Results/RhoAnalysis/Zee_sample_tau_trained0.4_tph0.35_tpi0_n3_g0.0/tau_traineddecay2_0.4_tph0.35_tpi0_n3_g0.0.root",
                                         "Results/RhoAnalysis/Zqq_sampletau_trained0.4_tph0.35_tpi0_n3_g0.0/tau_traineddecay2_0.4_tph0.35_tpi0_n3_g0.0.root",
                                         "Results/RhoAnalysis/bhabha_sample_tau_trained0.4_tph0.35_tpi0_n3_g0.0/tau_traineddecay2_0.4_tph0.35_tpi0_n3_g0.0.root"],
                     nargs=3, help="3 ROOT de fondo externo")
    ap.add_argument("--tree", default="outtree_original")
    ap.add_argument("--tauPcut", type=float, default=5, help="Corte fijo mínimo de recoMesonP (tauPCut)")
    ap.add_argument("--selectGEN", type=int, default=2, help="genTauID considerado señal (selectGEN)")

    # bounds (ajústalos a tus rangos físicos)
    ap.add_argument("--dR-bounds", type=float, nargs=2, default=[0.0, 6.0])
    ap.add_argument("--mesonP-bounds", type=float, nargs=2, default=[0.0, 100.0])
    ap.add_argument("--lepP-bounds", type=float, nargs=2, default=[0.0, 100.0])
    ap.add_argument("--cut-list", type=str, default="", help="Fichero .csv con lista de cortes a evaluar")

    # Loss config
    # ap.add_argument("--eff-target", type=float, default=0.90)
    # ap.add_argument("--eff-lambda", type=float, default=5.0)
    # ap.add_argument("--use-s-over-b", action="store_true", help="Si se activa: score = S/(B+eps). Si no: S/sqrt(S+B)")
    ap.add_argument(
    "--compare-yaml",
    default="config/plots/Optimal_Variable/Rho Decay/OptimalVariableBK_vals.yaml",
    help="YAML con luminosidad, xsec y n_events"
    )
    # PSO
    ap.add_argument("--particles", type=int, default=500)
    ap.add_argument("--iters", type=int, default=10000)
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

    if args.selectGEN == 2:
        samples = [
            RootSample(name="signal", path=args.signal_root, tree=args.tree, is_signal_file=True),
            RootSample(name="Zee", path=args.bg_root[0], tree=args.tree, is_signal_file=False),
            RootSample(name="Zqq", path=args.bg_root[1], tree=args.tree, is_signal_file=False),
            RootSample(name="Bhabha", path=args.bg_root[2], tree=args.tree, is_signal_file=False),
        ]
    elif args.selectGEN == 0:
        samples = [
            RootSample(name="signal", path=args.signal_root.replace("decay2", "decay0"), tree=args.tree, is_signal_file=True),
            RootSample(name="Zee", path=args.bg_root[0].replace("decay2", "decay0"), tree=args.tree, is_signal_file=False),
            RootSample(name="Zqq", path=args.bg_root[1].replace("decay2", "decay0"), tree=args.tree, is_signal_file=False),
            RootSample(name="Bhabha", path=args.bg_root[2].replace("decay2", "decay0"), tree=args.tree, is_signal_file=False),
        ]
    branches = BranchMap(
        # si tienes una rama de pesos (p.ej. "weight" o "genWeight") ponla aquí:
        weight=None,
        omega_gen="genOmega",
        # omega_gen_p1="genOmega",
        # omega_gen_m1="Omega_GEN_SIGNAL_M1",
        # omega_signal="Omega_SIGNAL",
        omega_reco="omega",
        # recomeson_X="RecoMeson_X"
    )
    dataset_weights = load_dataset_weights(args.compare_yaml, args.selectGEN)
    loaded = load_samples(samples, branches=branches, tauP_min_fixed=args.tauPcut, dataset_weights=dataset_weights)
    for l in loaded:
        print(f"Loaded sample: {l.name}, events: {len(l.weights)}")
    loss_cfg = LossConfig(
        eff_target=float(0.9),
        eff_lambda=float(1),
        use_s_over_b=bool(False),
    )
    selectGEN = args.selectGEN
    if selectGEN == 2:
        selectGEN = 1
        optimal_variable = "omega"
    elif selectGEN == 0:
        optimal_variable = "meson_x"
    obj = CutObjective(
        loaded=loaded,
        branches=branches,
        selectGEN=selectGEN,
        loss_cfg=loss_cfg,
    )
    if args.use_styler:
        plt.style.use("/nfs/cms/arqolmo/SDHCAL_Energy/utils/newams.mplstyle")
        
    if args.cut_list:
        cut_df = pd.read_csv(args.cut_list)
        optimal_variable_shapes = []
        cut_param_list = []
        for index, row in cut_df.iterrows():
            p = CutParams(
                dR_min=row['dR_min'],
                dR_max=row['dR_max'],
                mesonP_min=row['mesonP_min'],
                mesonP_max=row['mesonP_max'],
                lepP_min=row['lepP_min'],
                lepP_max=row['lepP_max'],
                zmass_min=row['Zmass_min'],
                zmass_max=row['Zmass_max'],
            )
            cut_param_list.append(p)
    else:
        p = CutParams(
            dR_min=0,
            dR_max=5.,
            mesonP_min=0.0,
            mesonP_max=100.0,
            lepP_min=0.0,
            lepP_max=100.0,
            zmass_min=0.0,
            zmass_max=100.0,
        )
        cut_param_list = [p]
        optimal_variable_shapes = []
    for _, p in enumerate(cut_param_list):
        
        objective_result = obj.evaluate(p)
        sig = [s for s in loaded if s.is_signal_file][0]

        mask = sig.base_mask & (sig.genTauID == selectGEN)

        if optimal_variable == "omega":
            omega_vals = sig.arrays[branches.omega_reco][mask]
            omega_w = sig.weights[mask]

            optimal_variable_shapes.append({
                "cuts": p,
                "omega": omega_vals.copy(),
                "weights": omega_w.copy(),
            })
        elif optimal_variable == "meson_x":
            meson_x_vals = sig.mesonX[mask]
            meson_x_w = sig.weights[mask]
            
            optimal_variable_shapes.append({
                "cuts": p,
                "meson_x": meson_x_vals.copy(),
                "weights": meson_x_w.copy(),
            })
            # print("Signal: ",objective_result.S)
        # print(objective_result.details)
        result_values = {}
        total_events = 0
        for key in objective_result.details:
            if key != "S_baseline":
                result_values[key] = objective_result.details[key]
                total_events += objective_result.details[key]
        total_events += objective_result.S
        result_values["Signal"] = objective_result.S
        result_values["Total events"] = total_events
        frac_result_values = {}
        for key in result_values:
            if key != "Total events":
                frac_result_values[key+"_frac"] = result_values[key]/total_events
        
        result_values.update(frac_result_values)
        cut_values = asdict(p)
        result_values.update(cut_values)
        df = pd.DataFrame([result_values])
        # Ordenar por columna
        df = df.reindex(sorted(df.columns), axis=1)
        print("\n=== EVALUATION RESULTS ===")
        print(df.to_string(index=False))
        
        outdir = f"validation_plots_single_cut_{args.selectGEN}/"
        if os.path.exists(outdir + "evaluation_results.csv"):
            df_existing = pd.read_csv(outdir + "evaluation_results.csv")
            row = df.iloc[0]

            exists = df_existing.eq(row).all(axis=1).any()

            if not exists:
                df_existing = pd.concat([df_existing, df], ignore_index=True)
                df_existing.to_csv(outdir + "evaluation_results.csv", index=False)
            else:
                print("Results already exist in evaluation_results.csv.")
            # Check if any row matches new data
            # if not ((df_existing == df).all(axis=1)).any():
            #     df_existing = pd.concat([df_existing, df], ignore_index=True)
            #     df_existing.to_csv(outdir + "evaluation_results.csv", index=False)
            # else:
            #     print("Results already exist in evaluation_results.csv.")
        else:
            os.makedirs(outdir, exist_ok=True)
            df.to_csv(outdir + "evaluation_results.csv", index=False)
        
                # Parámetros optimizados: [dR_min, dR_max, mesonP_min, mesonP_max, lepP_min, lepP_max]
    #     # Ojo: bounds de min/max por separado. Si el PSO propone min>max, objective lo penaliza.
    #     dR_lo, dR_hi = map(float, args.dR_bounds)
    #     m_lo, m_hi = map(float, args.mesonP_bounds)
    #     l_lo, l_hi = map(float, args.lepP_bounds)

    #     bounds = [
    #         (dR_lo, dR_hi),  # dR_min
    #         (dR_lo, dR_hi),  # dR_max
    #         (m_lo, m_hi),    # mesonP_min
    #         (m_lo, m_hi),    # mesonP_max
    #         (l_lo, l_hi),    # lepP_min
    #         (l_lo, l_hi),    # lepP_max
    #     ]

    #     def f(x: np.ndarray) -> float:
    #         p = CutParams(
    #             dR_min=float(x[0]),
    #             dR_max=float(x[1]),
    #             mesonP_min=float(x[2]),
    #             mesonP_max=float(x[3]),
    #             lepP_min=float(x[4]),
    #             lepP_max=float(x[5]),
    #         )
    #         return obj.evaluate(p).loss

    #     pso_cfg = PSOConfig(
    #         n_particles=int(args.particles),
    #         n_iters=int(args.iters),
    #         seed=int(args.seed),
    #         verbose=bool(args.verbose),
    #         patience=200
    #         # velocity_clamp=0.2,  # opcional
    #     )

    #     res = pso_minimize(f, bounds=bounds, cfg=pso_cfg)

    #     best = CutParams(
    #         dR_min=float(res.best_x[0]),
    #         dR_max=float(res.best_x[1]),
    #         mesonP_min=float(res.best_x[2]),
    #         mesonP_max=float(res.best_x[3]),
    #         lepP_min=float(res.best_x[4]),
    #         lepP_max=float(res.best_x[5]),
    #     )
    #     best_eval = obj.evaluate(best)

    #     print("\n=== BEST CUTS ===")
    #     print(asdict(best))
    #     print("\n=== METRICS ===")
    #     print({
    #         "loss": best_eval.loss,
    #         "S": best_eval.S,
    #         "B": best_eval.B,
    #         "effS": best_eval.effS,
    #         **best_eval.details
    #     })
    #     print("\n=== PSO ===")
    #     print({"best_f": res.best_f, "n_steps": len(res.history_best)})
    #     print("\nPlotting validation histograms...")
        if len(cut_param_list) == 1:
            plot_variable_with_cuts(
                loaded, branches, selectGEN, p,
                var="dR", bins=60, outpath=outdir, range=(0,5)
            )

            plot_variable_with_cuts(
                loaded, branches, selectGEN, p,
                var="mesonP", bins=60, outpath=outdir, range=(0,50)
            )

            plot_variable_with_cuts(
                loaded, branches, selectGEN, p,
                var="lepP", bins=60, outpath=outdir, range=(0,50)
            )
            
            plot_variable_with_cuts(
                loaded, branches, selectGEN, p,
                var="ZMass", bins=60, outpath=outdir, range=(0,100)
            )
            if optimal_variable == "omega":
                plot_variable_with_cuts(
                    loaded,
                    branches,
                    selectGEN,
                    p,
                    var="omega",
                    bins=60,
                    range=(-1.0, 1.0),
                    outpath=outdir,
                )
            elif optimal_variable == "meson_x":
                 plot_variable_with_cuts(
                    loaded,
                    branches,
                    selectGEN,
                    p,
                    var="meson_x",
                    bins=60,
                    range=(0.0, 1.0),
                    outpath=outdir,
                )
            # plot_variable_with_cuts(
            # loaded,
            # branches,
            # selectGEN,
            # p,
            # var="omega",
            # bins=60,
            # range=(-1.0, 1.0),
            # outpath=outdir,
            # )


#     plot_loss_history(
#     res.history_best,
#     outpath="validation_plots/",
# )
    if len(optimal_variable_shapes) > 1:
        range_lim = (-1.0, 1.0) if optimal_variable == "omega" else (0.0, 1.0)
        plot_omega_shape_comparison_from_arrays(
            optimal_shapes=optimal_variable_shapes,
            bins=60,
            range=range_lim,
            outpath=outdir,
            optimal_var=optimal_variable,
            use_vars_as_labels=False
        )
if __name__ == "__main__":
    main()
