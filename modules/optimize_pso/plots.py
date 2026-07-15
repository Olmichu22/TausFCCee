# cutopt/plots.py
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from .data import LoadedSample, BranchMap
from .objective import CutParams


def plot_loss_history(
    loss_history,
    outpath: str,
    fname: str = "loss_evolution.png",
):
    os.makedirs(outpath, exist_ok=True)

    it = np.arange(1, len(loss_history) + 1)

    plt.figure(figsize=(7, 4))
    plt.plot(it, loss_history, lw=2)
    plt.xlabel("Iteración PSO")
    plt.ylabel("Best loss")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig(os.path.join(outpath, fname))
    plt.close()

def plot_variable_with_cuts(
    loaded: list[LoadedSample],
    branches: BranchMap,
    selectGEN: int,
    cuts: CutParams,
    var: str,
    bins: int = 60,
    range: tuple[float, float] | None = None,
    outpath: str = "",
):
    """
    var ∈ {"dR", "mesonP", "lepP"}
    Dibuja histogramas ponderados por dataset:
      - Señal
      - Migraciones (BG interno)
      - BG externos
    """

    if outpath:
        os.makedirs(outpath, exist_ok=True)

    # ------------------ separar samples ------------------
    sig = [s for s in loaded if s.is_signal_file][0]
    bgs = [s for s in loaded if not s.is_signal_file]

    # ------------------ selector de variable ------------------
    def get_array(s: LoadedSample):
        if var == "dR":
            return s.dR
        if var == "mesonP":
            return s.arrays[branches.recoMesonP]
        if var == "lepP":
            return s.arrays[branches.lepP]
        if var == "omega":
            # Señal y migraciones usan Omega_SIGNAL
            # if s.is_signal_file:
            #     return s.arrays[branches.omega_signal]
            # # Fondos externos usan OMEGA_BG
            return s.arrays[branches.omega_reco]
        if var == "ZMass":
            return s.ZMass

        if var == "meson_x":
            return s.mesonX
        raise ValueError(f"Variable desconocida: {var}")

    # ==========================================================
    # SIGNAL & MIGRATIONS (mismo fichero, distinto genTauID)
    # ==========================================================
    base = sig.base_mask

    is_signal = sig.genTauID == selectGEN
    is_migration = ~is_signal

    sig_mask = base & is_signal
    mig_mask = base & is_migration

    sig_vals = get_array(sig)[sig_mask]
    sig_w = sig.weights[sig_mask]

    mig_vals = get_array(sig)[mig_mask]
    mig_w = sig.weights[mig_mask]

    # ==========================================================
    # EXTERNAL BACKGROUND (otros ficheros)
    # ==========================================================
    bg_vals = []
    bg_w = []
    plt.figure(figsize=(7, 5))

    for bg in bgs:
        m = bg.base_mask
        # Ignoring Zee for the moment
        if bg.name == "Zee":
            continue
        # if bg.name!= "Zqq":
            # continue
        # bg_vals.append(get_array(bg)[m])
        # bg_w.append(bg.weights[m])

    # if len(bg_vals) > 0:
    #     bg_vals = np.concatenate(bg_vals)
    #     bg_w = np.concatenate(bg_w)
    # else:
    #     bg_vals = np.array([])
    #     bg_w = np.array([])

    # ==========================================================
    # PLOT
    # ==========================================================

        # BG externo
        plt.hist(
            get_array(bg)[m],
            bins=bins,
            range=range,
            weights=bg.weights[m],
            histtype="step",
            # alpha=0.35,
            label=bg.name,
            # color="tab:gray",
        )

    # Migraciones
    plt.hist(
        mig_vals,
        bins=bins,
        range=range,
        weights=mig_w,
        histtype="step",
        linewidth=1.,
        label="Migrations",
        color="tab:red",
    )

    # Señal
    plt.hist(
        sig_vals,
        bins=bins,
        range=range,
        weights=sig_w,
        histtype="step",
        linewidth=2.2,
        label="Signal",
        color="tab:green",
        linestyle="--",
    )

    # ------------------ líneas de corte ------------------
    if var == "dR":
        if cuts.dR_min != 0.0:
            plt.axvline(cuts.dR_min, color="k", linestyle="--")
        if cuts.dR_max < 5.0:
            plt.axvline(cuts.dR_max, color="k", linestyle="--")
        plt.xlabel("ΔR (rads)")
        plt.title("ΔR between meson and lepton")
    elif var == "mesonP":
        if cuts.mesonP_min != 0.0:
            plt.axvline(cuts.mesonP_min, color="k", linestyle="--")
        if cuts.mesonP_max < 50.0:
            plt.axvline(cuts.mesonP_max, color="k", linestyle="--")
        plt.xlabel("P (GeV/c)")
        plt.title("Meson momentum")
    elif var == "lepP":
        if cuts.lepP_min != 0.0:
            plt.axvline(cuts.lepP_min, color="k", linestyle="--")
        if cuts.lepP_max < 50.0:
            plt.axvline(cuts.lepP_max, color="k", linestyle="--")
        plt.xlabel("P (GeV/c)")
        plt.title("Lepton momentum")
    elif  var == "omega":
        plt.title(r"Optimal variable for $\rho$ decay")
        plt.xlabel(r"$\omega$")
    elif var == "ZMass":
        if cuts.zmass_min != 0.0:
            plt.axvline(cuts.zmass_min, color="k", linestyle="--")
        if cuts.zmass_max < 100.0:
            plt.axvline(cuts.zmass_max, color="k", linestyle="--")
        plt.xlabel("Z Mass (GeV/c²)")
        plt.title("Invariant mass of the Z boson")
    elif var == "meson_x":
        plt.xlabel(r"$\frac{E_{Meson}}{E_{Beam}}$")
        plt.title(r"Optimal variable for $\pi$ decay")

    plt.ylabel("Events (weighted)")
    # plt.yscale("log")
    plt.legend()
    plt.tight_layout()

    # if var == "dR":
        # plt.yscale("log")
    fname = f"histograma_cortes_{var}.png"
    plt.savefig(os.path.join(outpath, fname))
    plt.close()

from dataclasses import asdict

def plot_omega_shape_comparison_from_arrays(
    optimal_shapes: list[dict],
    bins: int = 60,
    range: tuple[float, float] = (-1.0, 1.0),
    outpath: str = "",
    fname: str = "omega_shape_comparison.png",
    optimal_var: str = "omega",
    use_vars_as_labels: bool = True,   
):
    """
    Dibuja la comparación de formas de omega reco a partir de arrays
    ya filtrados (post-corte).
    """

    if outpath:
        os.makedirs(outpath, exist_ok=True)

    if len(optimal_shapes) < 2:
        return

    # ------------------ detectar qué cortes cambian ------------------
    cut_dicts = [asdict(d["cuts"]) for d in optimal_shapes]
    keys = cut_dicts[0].keys()

    varying_keys = [
        k for k in keys
        if len({d[k] for d in cut_dicts}) > 1
    ]
    plt.figure(figsize=(7, 5))

    for i, d in enumerate(optimal_shapes):
        optimal = d[optimal_var]
        weights = d["weights"]
        cuts = d["cuts"]

        if optimal.size == 0:
            continue

        label_parts = [
            f"{k}={getattr(cuts, k):.2f}"
            for k in varying_keys
        ]
        if use_vars_as_labels:
            label = ", ".join(label_parts) if label_parts else f"Corte {i+1}"
        else:
            label = f"Corte {i+1}"

        plt.hist(
            optimal,
            bins=bins,
            range=range,
            weights=weights,
            # density=True,
            histtype="step",
            linewidth=2,
            label=label,
            # Línea discontinua
            # linestyle="--"
        )
    if optimal_var == "omega":
        plt.title(r"$\omega$ reco comparison after cuts")
        plt.xlabel(r"$\omega$")
    elif optimal_var == "meson_x":
        plt.xlabel(r"$\frac{E_{Meson}}{E_{Beam}}$")
        plt.title(r"Meson $\frac{E_{Meson}}{E_{Beam}}$ comparison after cuts")
    plt.ylabel("Events (weighted)")
    plt.legend(fontsize=9)
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.savefig(os.path.join(outpath, fname))
    plt.close()

