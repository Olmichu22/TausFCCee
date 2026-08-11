#!/usr/bin/env python
"""
Regenera todos los plots de particle_level_analisis_parallel.py a partir de los
parquets de asociaciones ya escritos, sin releer los ficheros ROOT.

Uso típico:
    python HitAnalysis/replot_from_parquet.py <outputpath> --fake-bin-by-reco \
        --all-plot 22 211 13 11

<outputpath> es el directorio que contiene association_results_full_dR.parquet y
association_results_full_truthlink.parquet (los .pkl.gz también valen).  Los
plots se escriben en el mismo árbol de subdirectorios que produce el análisis;
usa --outdir para volcarlos en otro sitio y no pisar los originales.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.ConfusionMatrixParticleLevel import (plot_confusion_matrices,
                                                  plot_energy_distributions,
                                                  plot_efficiency_vs_momentum,
                                                  plot_fake_rate_vs_momentum)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from particle_level_analisis_parallel import build_association_structures

# Mismos bins que el análisis (particle_level_analisis_parallel.py:820-830)
BINS = [0, 1, 5, 10, 20, 30, 45, 100, np.inf]
E_BINS = np.linspace(0.0, 40.0, 40)


def load_df(outputpath, basename):
    for ext, loader in ((".parquet", pd.read_parquet), (".pkl.gz", pd.read_pickle)):
        path = os.path.join(outputpath, basename + ext)
        if os.path.exists(path):
            print(f"Leyendo {path} ...")
            return loader(path), path
    return None, None


def replot(full_df, tag, outdir, args, selection_note):
    """Reproduce las etapas de plot del análisis para una de las dos ramas."""
    if full_df is None or full_df.empty:
        print(f"[{tag}] DataFrame vacío o ausente. Saltado.")
        return

    n_events = int(full_df["event_id"].nunique())
    print(f"[{tag}] {len(full_df):,} filas, {n_events:,} eventos")

    assoc, edist = build_association_structures(
        full_df, BINS, E_BINS, fake_bin_by_reco=args.fake_bin_by_reco
    )

    def d(*parts):
        return os.path.join(outdir, *parts)

    plot_confusion_matrices(assoc, output_dir=d("confusion_matrices_particle_level", tag))
    plot_energy_distributions(edist, output_dir=d("energy_distributions", tag),
                              all_plot_pdgs=args.all_plot)

    plot_efficiency_vs_momentum(full_df, output_dir=d("efficiency_plots", tag))
    plot_efficiency_vs_momentum(full_df, output_dir=d("efficiency_plots_theta", tag),
                                plot_type="theta")

    plot_fake_rate_vs_momentum(full_df, output_dir=d("fake_rate_plots", tag),
                               selection_note=selection_note, n_events=n_events)
    plot_fake_rate_vs_momentum(full_df, output_dir=d("fake_rate_plots_theta", tag),
                               plot_type="theta",
                               selection_note=selection_note, n_events=n_events)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outputpath", help="Directorio con los parquets de asociaciones")
    p.add_argument("--outdir", default=None,
                   help="Directorio de salida (por defecto: el mismo outputpath)")
    p.add_argument("--fake-bin-by-reco", action="store_true", default=False,
                   help="Binear las filas fake por Reco_energy en las matrices")
    p.add_argument("--all-plot", type=int, nargs="*", default=[22], metavar="PDG",
                   help="PDGs con plots combinados de resolución")
    p.add_argument("--branch", choices=["dR", "truthlink", "both"], default="both",
                   help="Qué rama de matching replotear")
    p.add_argument("--selection-note", default="",
                   help="Nota de selección para los plots de fakes; por defecto se "
                        "reconstruye desde el config.yaml del directorio si existe")
    args = p.parse_args()

    outdir = args.outdir or args.outputpath
    os.makedirs(outdir, exist_ok=True)

    note = args.selection_note
    if not note:
        cfg_path = os.path.join(args.outputpath, "config.yaml")
        note = (f"replot desde parquet ({os.path.basename(os.path.normpath(args.outputpath))})"
                + (f"  |  cfg: {cfg_path}" if os.path.exists(cfg_path) else ""))

    branches = ["dR", "truthlink"] if args.branch == "both" else [args.branch]
    for tag in branches:
        df, path = load_df(args.outputpath, f"association_results_full_{tag}")
        if df is None:
            print(f"[{tag}] No se encontró association_results_full_{tag}.parquet/.pkl.gz")
            continue
        replot(df, tag, outdir, args, f"matching {tag}  |  {note}")

    print(f"Plots escritos en {outdir}")


if __name__ == "__main__":
    main()
