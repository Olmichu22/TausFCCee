#!/usr/bin/env python
"""
Regenera todos los plots de particle_level_analisis_parallel.py a partir de los
parquets de asociaciones ya escritos, sin releer los ficheros ROOT.

Uso típico:
    python HitAnalysis/replot_from_parquet.py <outputpath> --fake-bin-by-reco \
        --all-plot 22 211 13 11 --min-energy-cuts 10 20

Comparación de dos (o más) directorios, p.ej. dos detectores:
    python HitAnalysis/replot_from_parquet.py <outCLD> <outILD> --labels CLD ILD \
        --min-energy-cuts 10

Con más de un <outputpath> el script pasa a modo comparación: en vez de rehacer
todos los plots, superpone las distribuciones de resolución de momento y de
theta de cada directorio en un mismo PNG (momentum_resolution/ y
theta_resolution/, resolution_matched_<pid>.png y
resolution_pred_<pid>.png) junto con los perfiles de resolución vs energía
(energy_distributions/, residual_resolution_<metrica>_<serie>.png) y el fake
rate y el fake yield (fake_rate_plots/ y fake_rate_plots_theta/,
fake_rate_<pid>.png y fake_yield_<pid>.png).

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
                                                  plot_fake_rate_vs_momentum,
                                                  plot_momentum_resolution,
                                                  plot_theta_resolution,
                                                  compare_momentum_resolution,
                                                  compare_theta_resolution,
                                                  compare_energy_resolution,
                                                  compare_fake_rate_vs_momentum)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from particle_level_analisis_parallel import (build_association_structures,
                                              apply_min_energy_cut,
                                              energy_cut_label)

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

    plot_momentum_resolution(full_df, output_dir=d("momentum_resolution", tag))
    plot_theta_resolution(full_df, output_dir=d("theta_resolution", tag))

    # ── Variantes con corte mínimo en energía (--min-energy-cuts) ────────────
    # Mismo tratamiento que el análisis: matriz integrada (un único bin
    # E > umbral) y resoluciones de momento y theta sobre el mismo subconjunto.
    for threshold in (args.min_energy_cuts or []):
        cut_label = energy_cut_label(threshold)
        df_cut = apply_min_energy_cut(full_df, threshold, mode=args.min_energy_var)
        print(f"[{tag}] corte E > {threshold:g} GeV (var={args.min_energy_var}): "
              f"{len(df_cut):,} de {len(full_df):,} filas")
        if df_cut.empty:
            print(f"[{tag}] corte E > {threshold:g} GeV: sin filas, se omiten los plots")
            continue

        # fake_bin_by_reco=True: con un único bin no reparte los fakes, solo
        # evita que caigan en la etiqueta "nan" y generen una matriz espuria.
        assoc_cut, _ = build_association_structures(
            df_cut, [threshold, np.inf], E_BINS, fake_bin_by_reco=True
        )
        plot_confusion_matrices(
            assoc_cut,
            output_dir=d("confusion_matrices_particle_level", tag, cut_label),
        )
        plot_momentum_resolution(
            df_cut,
            output_dir=d("momentum_resolution", tag, cut_label),
            title_suffix=f" | E > {threshold:g} GeV",
        )
        plot_theta_resolution(
            df_cut,
            output_dir=d("theta_resolution", tag, cut_label),
            title_suffix=f" | E > {threshold:g} GeV",
        )


def compare(dfs_by_label, tag, outdir, args):
    """Superpone la resolución de momento y theta de varios directorios para una rama."""
    for label, df in dfs_by_label.items():
        print(f"[{tag}] {label}: {len(df):,} filas, "
              f"{int(df['event_id'].nunique()):,} eventos")

    def d(*parts):
        return os.path.join(outdir, *parts)

    def compare_all(dfs, subdir_parts, title_suffix=""):
        compare_momentum_resolution(
            dfs, output_dir=d("momentum_resolution", *subdir_parts),
            normalize=not args.no_normalize, pids=args.compare_pids,
            legend_fontsize=args.legend_fontsize, title_suffix=title_suffix,
        )
        compare_theta_resolution(
            dfs, output_dir=d("theta_resolution", *subdir_parts),
            normalize=not args.no_normalize, pids=args.compare_pids,
            legend_fontsize=args.legend_fontsize, title_suffix=title_suffix,
        )
        # Perfiles de resolución vs energía: hay que rehacer las estructuras
        # por dataset, son ellas (no el DataFrame) lo que consume el plot.
        edists = {}
        for label, df in dfs.items():
            _, edist = build_association_structures(
                df, BINS, E_BINS, fake_bin_by_reco=args.fake_bin_by_reco
            )
            edists[label] = edist
        compare_energy_resolution(
            edists, output_dir=d("energy_distributions", *subdir_parts),
            all_plot_pdgs=args.all_plot, title_suffix=title_suffix,
            legend_fontsize=args.legend_fontsize,
        )
        n_events = {label: int(df["event_id"].nunique()) for label, df in dfs.items()}
        compare_fake_rate_vs_momentum(
            dfs, output_dir=d("fake_rate_plots", *subdir_parts),
            n_events_by_label=n_events, pids=args.compare_pids,
            legend_fontsize=args.legend_fontsize, title_suffix=title_suffix,
        )
        compare_fake_rate_vs_momentum(
            dfs, output_dir=d("fake_rate_plots_theta", *subdir_parts),
            plot_type="theta", n_events_by_label=n_events, pids=args.compare_pids,
            legend_fontsize=args.legend_fontsize, title_suffix=title_suffix,
        )

    compare_all(dfs_by_label, [tag])

    for threshold in (args.min_energy_cuts or []):
        cut_label = energy_cut_label(threshold)
        cut_dfs = {}
        for label, df in dfs_by_label.items():
            df_cut = apply_min_energy_cut(df, threshold, mode=args.min_energy_var)
            print(f"[{tag}] {label}: corte E > {threshold:g} GeV "
                  f"(var={args.min_energy_var}): {len(df_cut):,} de {len(df):,} filas")
            if not df_cut.empty:
                cut_dfs[label] = df_cut
        if not cut_dfs:
            print(f"[{tag}] corte E > {threshold:g} GeV: sin filas, se omiten los plots")
            continue
        compare_all(cut_dfs, [tag, cut_label], title_suffix=f" | E > {threshold:g} GeV")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outputpath", nargs="+",
                   help="Directorio(s) con los parquets de asociaciones. Con más de "
                        "uno se activa el modo comparación")
    p.add_argument("--outdir", default=None,
                   help="Directorio de salida (por defecto: el mismo outputpath; en "
                        "modo comparación, <primer outputpath>/comparison)")
    p.add_argument("--fake-bin-by-reco", action="store_true", default=False,
                   help="Binear las filas fake por Reco_energy en las matrices")
    p.add_argument("--all-plot", type=int, nargs="*", default=[22], metavar="PDG",
                   help="PDGs con plots combinados de resolución")
    p.add_argument("--branch", choices=["dR", "truthlink", "both"], default="both",
                   help="Qué rama de matching replotear")
    p.add_argument("--min-energy-cuts", type=float, nargs="*", default=[], metavar="GEV",
                   help="Umbrales de energía (GeV). Para cada valor se rehacen la "
                        "matriz de confusión integrada (un único bin E > umbral) y "
                        "las resoluciones de momento y theta en un subdirectorio "
                        "Emin_<valor>GeV")
    p.add_argument("--min-energy-var", choices=["auto", "gen", "reco"], default="auto",
                   help="Energía usada por --min-energy-cuts: 'auto' (Gen_energy, o "
                        "Reco_energy para los fakes), 'gen' o 'reco'")
    p.add_argument("--labels", nargs="*", default=None, metavar="LABEL",
                   help="Modo comparación: nombre de cada directorio en la leyenda "
                        "(por defecto, el basename del directorio)")
    p.add_argument("--compare-pids", type=int, nargs="*", default=None, metavar="PDG",
                   help="Modo comparación: limitar los plots a estos |PDG|")
    p.add_argument("--legend-fontsize", type=float, default=11,
                   help="Modo comparación: tamaño de fuente de la leyenda (11 por defecto)")
    p.add_argument("--no-normalize", action="store_true", default=False,
                   help="Modo comparación: superponer cuentas absolutas en vez de "
                        "densidades (por defecto se normaliza al área)")
    p.add_argument("--selection-note", default="",
                   help="Nota de selección para los plots de fakes; por defecto se "
                        "reconstruye desde el config.yaml del directorio si existe")
    args = p.parse_args()

    branches = ["dR", "truthlink"] if args.branch == "both" else [args.branch]

    # ── Modo comparación: varios directorios superpuestos ────────────────────
    if len(args.outputpath) > 1:
        labels = args.labels or [os.path.basename(os.path.normpath(o))
                                 for o in args.outputpath]
        if len(labels) != len(args.outputpath):
            p.error(f"--labels tiene {len(labels)} entradas y se han dado "
                    f"{len(args.outputpath)} directorios")

        outdir = args.outdir or os.path.join(args.outputpath[0], "comparison")
        os.makedirs(outdir, exist_ok=True)

        for tag in branches:
            dfs_by_label = {}
            for label, path in zip(labels, args.outputpath):
                df, _ = load_df(path, f"association_results_full_{tag}")
                if df is None:
                    print(f"[{tag}] {label}: no se encontró "
                          f"association_results_full_{tag}.parquet/.pkl.gz en {path}")
                    continue
                dfs_by_label[label] = df
            if len(dfs_by_label) < 2:
                print(f"[{tag}] menos de dos directorios con datos, se omite")
                continue
            compare(dfs_by_label, tag, outdir, args)

        print(f"Plots de comparación escritos en {outdir}")
        return

    # ── Modo normal: un único directorio ─────────────────────────────────────
    outputpath = args.outputpath[0]
    outdir = args.outdir or outputpath
    os.makedirs(outdir, exist_ok=True)

    note = args.selection_note
    if not note:
        cfg_path = os.path.join(outputpath, "config.yaml")
        note = (f"replot desde parquet ({os.path.basename(os.path.normpath(outputpath))})"
                + (f"  |  cfg: {cfg_path}" if os.path.exists(cfg_path) else ""))

    for tag in branches:
        df, path = load_df(outputpath, f"association_results_full_{tag}")
        if df is None:
            print(f"[{tag}] No se encontró association_results_full_{tag}.parquet/.pkl.gz")
            continue
        replot(df, tag, outdir, args, f"matching {tag}  |  {note}")

    print(f"Plots escritos en {outdir}")


if __name__ == "__main__":
    main()
