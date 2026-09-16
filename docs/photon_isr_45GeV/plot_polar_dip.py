"""Polar-angle distribution of generated photons around the 45 GeV efficiency dip.

Draws, for the energy bin where the photon efficiency drops (43.3-45.0 GeV) and
its two neighbours, the polar-angle spectrum of generated photons split by
reconstruction outcome.  The dip is a sample-composition effect: half of the
photons in that bin are beam-collinear ISR photons that escape down the beam
pipe.

Input is the particle-level association dataframe (dR matching) produced by
HitAnalysis/particle_level_analisis_parallel.py.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from matplotlib.ticker import AutoMinorLocator

# Bordes del binning de energia usado en plot_efficiency_vs_momentum (30 bins, 0-50)
ENERGY_BINS = [(41.667, 43.333), (43.333, 45.0), (45.0, 46.667)]

COLORS = {
    "gamma": "#2166AC",   # reconstruido como foton
    "other": "#F4A582",   # reconstruido, PID equivocado
    "lost": "#B2182B",    # sin contrapartida reco
}


def load_photons(path, p_min=40.0):
    """Read generated photons above ``p_min`` from the association parquet."""
    cols = ["Gen_pid", "Reco_pid", "Gen_Px", "Gen_Py", "Gen_Pz"]
    pf = pq.ParquetFile(path)
    chunks = []
    for i in range(pf.metadata.num_row_groups):
        df = pf.read_row_group(i, columns=cols).to_pandas()
        df = df[df["Gen_pid"].abs() == 22]
        p = np.sqrt(df.Gen_Px ** 2 + df.Gen_Py ** 2 + df.Gen_Pz ** 2)
        df = df.assign(P=p)[p > p_min]
        if len(df):
            chunks.append(df)
    out = pd.concat(chunks, ignore_index=True)
    cos_theta = out.Gen_Pz / out.P
    # Angulo polar en radianes, sin plegar: los dos picos ISR (+z y -z) quedan
    # en los extremos del rango
    out["theta"] = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    out["theta_beam"] = np.minimum(out["theta"], np.pi - out["theta"])
    out["outcome"] = np.where(
        out.Reco_pid.abs() == 999, "lost",
        np.where(out.Reco_pid.abs() == 22, "gamma", "other"),
    )
    return out


def make_figure(df, detector, outdir, dpi=200):
    edges = np.linspace(0.0, np.pi, 61)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    fig, axes = plt.subplots(
        1, 3, figsize=(13.5, 5.0), sharey=True,
        gridspec_kw=dict(wspace=0.06, left=0.065, right=0.985, top=0.80, bottom=0.145),
    )
    y_max = 0.0

    for ax, (lo, hi) in zip(axes, ENERGY_BINS):
        sub = df[(df.P >= lo) & (df.P < hi)]
        bottom = np.zeros(len(centers))
        for key, label in (
            ("gamma", r"reconstructed as $\gamma$"),
            ("other", "reconstructed, wrong PID"),
            ("lost", "not reconstructed"),
        ):
            counts, _ = np.histogram(sub.theta[sub.outcome == key], bins=edges)
            ax.bar(centers, counts, width=widths, bottom=bottom, align="center",
                   color=COLORS[key], edgecolor="white", linewidth=0.3,
                   label=label if ax is axes[0] else None, zorder=2)
            bottom += counts
        y_max = max(y_max, bottom.max())

        eff = float((sub.outcome == "gamma").mean())
        forward = float((sub.theta_beam < 0.035).mean())

        ax.set_xlim(0.0, np.pi)
        ax.set_xticks(np.arange(0, 5) * np.pi / 4.0)
        ax.set_xticklabels(["0", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
        ax.set_xlabel(r"$\theta_{\gamma}^{\mathrm{gen}}$  [rad]")
        ax.set_title(rf"$p_{{\gamma}}^{{\mathrm{{gen}}}} \in [{lo:.1f},\ {hi:.1f}]$ GeV",
                     fontsize=12, pad=8)
        ax.text(0.035, 0.965,
                f"N = {len(sub):,}\n"
                rf"$\varepsilon_{{\gamma}}$ = {eff:.3f}" + "\n"
                rf"$\theta_{{\mathrm{{beam}}}} < 0.035$ rad: {forward:.0%}",
                transform=ax.transAxes, va="top", ha="left", fontsize=10.5, zorder=5,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.95))

        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.tick_params(which="major", length=6, width=1.1)
        ax.tick_params(which="minor", length=3, width=0.8)
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

    axes[0].set_ylabel("Generated photons / bin")
    axes[0].set_ylim(0.0, 1.55 * y_max)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.985, 0.985),
               ncol=3, fontsize=11, frameon=False)

    fig.text(0.065, 0.925, f"{detector}", fontsize=15, fontweight="bold", va="bottom")
    fig.text(0.115, 0.929, r"$Z \to \tau^{+}\tau^{-}$, full simulation",
             fontsize=11.5, va="bottom", color="0.3")

    os.makedirs(outdir, exist_ok=True)
    stem = os.path.join(outdir, f"photon_polar_angle_dip_{detector}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", dpi=dpi)
        print(f"Saved -> {stem}.{ext}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True,
                    help="association_results_full_dR.parquet")
    ap.add_argument("-d", "--detector", default="ILD")
    ap.add_argument("-o", "--outdir", default="docs/photon_isr_45GeV")
    args = ap.parse_args()

    df = load_photons(args.input)
    make_figure(df, args.detector, args.outdir)


if __name__ == "__main__":
    main()
