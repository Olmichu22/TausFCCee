"""
Detector-performance parametrisation of the full simulation, per particle type.

For the correctly identified matched pairs (|Gen_pid| == |Reco_pid|) of one
species it measures, in bins of the true energy/momentum and of |cos θ|:

  - the response  E_reco / E_true  (scale)  and the bias of every residual;
  - the resolution of (E_reco − E_true)/E_true, (|p_reco| − |p_gen|)/|p_gen|,
    θ_reco − θ_gen and φ_reco − φ_gen;

and fits the energy/momentum dependence with the usual FCC-ee baseline forms
(each term added in quadrature, ⊕):

  energy    σ_E/E = a / √E ⊕ b           (calorimetric, neutral particles)
  momentum  σ_p/p = a · p_T ⊕ b          (tracker, charged particles)
  θ, φ      σ     = a / p ⊕ b   [mrad]   (multiple scattering / cluster position)

The result is written as a summary table (YAML, CSV and PNG) next to one plot
per observable with the measured points and the fitted curve.

The track impact parameter of the reference table is not available here: the
association DataFrames only carry four-momenta.
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml

from modules.ConfusionMatrixParticleLevel import (_residual_frame, _resolution_value,
                                                  _resolution_error, pid_label)


# Nombres admitidos en la línea de comandos además del PDG
PARTICLE_ALIASES = {
    "photon": 22, "gamma": 22,
    "electron": 11, "e": 11,
    "muon": 13, "mu": 13,
    "pion": 211, "pi": 211,
    "kaon": 321, "k": 321,
    "proton": 2212, "p": 2212,
    "neutron": 2112, "n": 2112,
    "k0l": 130, "klong": 130,
}

CHARGED_PDGS = {11, 13, 211, 321, 2212}

# Bordes por defecto (GeV): finos a baja energía, donde varía la resolución
DEFAULT_FIT_E_BINS = [0.5, 1, 1.5, 2, 3, 4, 5, 7, 10, 13, 16, 20, 25, 30, 35, 40, 46]
# |cos θ|: barril, transición y endcap
DEFAULT_FIT_COS_BINS = [0.0, 0.7, 0.9, 1.0]

# Modelos σ(x) = (a·x^n) ⊕ b.  "xvar" es la variable gen frente a la que se
# binea y se ajusta; el exponente n fija la forma.
FIT_MODELS = {
    "stochastic": dict(n=-0.5, formula=r"\frac{%s}{\sqrt{E}} \oplus %s", xlabel="E$_{gen}$ [GeV]"),
    "linear":     dict(n=1.0,  formula=r"%s \cdot p_T \oplus %s",        xlabel="p$_{T,gen}$ [GeV]"),
    "inverse":    dict(n=-1.0, formula=r"\frac{%s}{p} \oplus %s",        xlabel="|p$_{gen}$| [GeV]"),
}


def parse_particle(value):
    """argparse type: accept a PDG code (sign ignored) or a name such as 'photon'."""
    key = str(value).strip().lower()
    if key in PARTICLE_ALIASES:
        return PARTICLE_ALIASES[key]
    try:
        return abs(int(key))
    except ValueError:
        import argparse
        raise argparse.ArgumentTypeError(
            f"Unknown particle {value!r}: use a PDG code or one of {sorted(PARTICLE_ALIASES)}"
        )


def _observables(pid):
    """Observables to fit for a species: (name, residual column, x column, model, unit, label)."""
    if pid in CHARGED_PDGS:
        main = ("momentum", "p_res", "Gen_PT", "linear", "", r"$\sigma_p / p$")
    else:
        main = ("energy", "E_res", "Gen_energy", "stochastic", "", r"$\sigma_E / E$")
    return [
        main,
        ("theta", "theta_res", "Gen_P", "inverse", "mrad", r"$\sigma_\theta$"),
        ("phi",   "phi_res",   "Gen_P", "inverse", "mrad", r"$\sigma_\phi$"),
    ]


def matched_frame(full_df, pid):
    """
    Rows correctly identified as `pid` with every residual used by the fits.

    Adds to the matched frame of _residual_frame the relative energy residual,
    the response E_reco/E_true, the transverse momentum and the φ residual
    (wrapped to (−π, π], in mrad).
    """
    _, df_res = _residual_frame(full_df, ctx="resolution_fits")
    if df_res is None:
        return None
    df = df_res.loc[(df_res["Gen_pid"] == pid) & (df_res["Reco_pid"] == pid)].copy()
    if df.empty:
        return df

    for col in ("Gen_energy", "Reco_energy"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["E_res"] = (df["Reco_energy"] - df["Gen_energy"]) / df["Gen_energy"]
    df["E_resp"] = df["Reco_energy"] / df["Gen_energy"]
    df["Gen_PT"] = np.hypot(df["Gen_Px"], df["Gen_Py"])
    df["Gen_abscos"] = np.abs(np.cos(df["Gen_theta"]))

    dphi = np.arctan2(df["Reco_Py"], df["Reco_Px"]) - np.arctan2(df["Gen_Py"], df["Gen_Px"])
    df["phi_res"] = (np.mod(dphi + np.pi, 2 * np.pi) - np.pi) * 1e3
    return df


def binned_profile(df, res_col, x_col, edges, metric="std90", min_entries=50):
    """
    Resolution, bias and response of `res_col` in bins of `x_col`.

    Returns a DataFrame with one row per populated bin: x (mean of x_col in the
    bin), x_lo, x_hi, n, sigma, sigma_err, bias (median residual), the residual
    quantiles q05, q25, q75 and q95 (box and whiskers of the box plot) and
    response (mean E_reco/E_true, only meaningful for the energy).
    """
    rows = []
    x = df[x_col].to_numpy(dtype=float)
    r = df[res_col].to_numpy(dtype=float)
    resp = df["E_resp"].to_numpy(dtype=float)
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (x >= lo) & (x < hi) & np.isfinite(r)
        n = int(sel.sum())
        if n < min_entries:
            continue
        sigma = _resolution_value(r[sel], metric)
        err = _resolution_error(r[sel], metric)
        if sigma is None or err is None:
            continue
        resp_sel = resp[sel][np.isfinite(resp[sel])]
        q05, q25, q75, q95 = np.percentile(r[sel], [5, 25, 75, 95])
        rows.append(dict(
            x=float(np.mean(x[sel])), x_lo=float(lo), x_hi=float(hi), n=n,
            sigma=float(sigma), sigma_err=float(err),
            bias=float(np.median(r[sel])),
            q05=float(q05), q25=float(q25), q75=float(q75), q95=float(q95),
            response=float(np.mean(resp_sel)) if resp_sel.size else np.nan,
            response_err=(float(np.std(resp_sel) / np.sqrt(resp_sel.size))
                          if resp_sel.size > 1 else np.nan),
        ))
    return pd.DataFrame(rows)


def fit_quadrature(prof, model):
    """
    Fit σ(x) = (a·x^n) ⊕ b to a binned profile.

    Returns dict(a, a_err, b, b_err, chi2, ndf) or None when there are fewer
    than three points or the fit fails.  a and b are constrained to be ≥ 0.
    """
    if prof is None or len(prof) < 3:
        return None
    from scipy.optimize import curve_fit

    n = FIT_MODELS[model]["n"]

    def f(x, a, b):
        return np.sqrt((a * np.power(x, n)) ** 2 + b ** 2)

    x = prof["x"].to_numpy()
    y = prof["sigma"].to_numpy()
    ey = prof["sigma_err"].to_numpy()
    ey = np.where(ey > 0, ey, np.max(ey[ey > 0]) if np.any(ey > 0) else 1.0)

    # Arranque: término dominante en el extremo donde manda cada uno
    i_a = np.argmin(x) if n < 0 else np.argmax(x)
    i_b = np.argmax(x) if n < 0 else np.argmin(x)
    p0 = [max(y[i_a] / np.power(x[i_a], n), 1e-6), max(y[i_b] * 0.5, 1e-6)]
    try:
        popt, pcov = curve_fit(f, x, y, p0=p0, sigma=ey, absolute_sigma=True,
                               bounds=([0.0, 0.0], [np.inf, np.inf]), maxfev=20000)
    except Exception as exc:
        print(f"[resolution_fits] fit did not converge ({model}): {exc}")
        return None
    perr = np.sqrt(np.clip(np.diag(pcov), 0, None))
    chi2 = float(np.sum(((y - f(x, *popt)) / ey) ** 2))
    return dict(a=float(popt[0]), a_err=float(perr[0]),
                b=float(popt[1]), b_err=float(perr[1]),
                chi2=chi2, ndf=int(len(x) - 2))


def _model_curve(model, a, b, x):
    return np.sqrt((a * np.power(x, FIT_MODELS[model]["n"])) ** 2 + b ** 2)


def _fmt_term(value, err, unit):
    """Relative terms in %, angular ones in mrad; a term compatible with zero prints as 0."""
    if not np.isfinite(err) or value < err:
        return "0"
    return f"{100 * value:.3g}\\%" if unit == "" else f"{value:.3g}\\,\\mathrm{{mrad}}"


def _formula(model, fit, unit):
    if fit is None:
        return "not fitted (< 3 populated bins)"
    return "$" + FIT_MODELS[model]["formula"] % (_fmt_term(fit["a"], fit["a_err"], unit),
                                                 _fmt_term(fit["b"], fit["b_err"], unit)) + "$"


def _plot_profiles(profiles, fits, model, unit, ylabel, title, fname, log_x=True):
    """Measured σ per |cos θ| region with the fitted curve of each region."""
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    colors = ["black", "#4363d8", "#e6194b", "#3cb44b", "#f58231", "#911eb4"]
    scale = 100.0 if unit == "" else 1.0
    for color, (region, prof) in zip(colors, profiles.items()):
        if prof is None or prof.empty:
            continue
        ax.errorbar(prof["x"], scale * prof["sigma"], yerr=scale * prof["sigma_err"],
                    xerr=[prof["x"] - prof["x_lo"], prof["x_hi"] - prof["x"]],
                    fmt="o", ms=4, color=color, capsize=2, label=region)
        fit = fits.get(region)
        if fit is not None:
            xx = np.geomspace(prof["x_lo"].min(), prof["x_hi"].max(), 200)
            ax.plot(xx, scale * _model_curve(model, fit["a"], fit["b"], xx), color=color,
                    lw=1.5, label=f"{_formula(model, fit, unit)}  "
                                  f"($\\chi^2$/ndf = {fit['chi2']:.1f}/{fit['ndf']})")
    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel(FIT_MODELS[model]["xlabel"])
    ax.set_ylabel(ylabel + (" [%]" if unit == "" else f" [{unit}]"))
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=11, ncol=2, framealpha=0.85)
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved resolution plot → {fname}")


def _plot_bias(profiles, model, unit, ylabel, title, fname, response=False):
    """Scale (E_reco/E_true) or median residual per bin and |cos θ| region."""
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    colors = ["black", "#4363d8", "#e6194b", "#3cb44b", "#f58231", "#911eb4"]
    scale = 100.0 if (unit == "" and not response) else 1.0
    for color, (region, prof) in zip(colors, profiles.items()):
        if prof is None or prof.empty:
            continue
        if response:
            ax.errorbar(prof["x"], prof["response"], yerr=prof["response_err"],
                        fmt="o", ms=4, color=color, capsize=2, label=region)
        else:
            ax.plot(prof["x"], scale * prof["bias"], "o-", ms=4, color=color, label=region)
    ax.axhline(1.0 if response else 0.0, color="grey", ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(FIT_MODELS[model]["xlabel"])
    if response:
        ax.set_ylabel(r"$\langle E_{reco} / E_{true} \rangle$")
    else:
        # Lo que se dibuja es la mediana del residuo, no la resolución σ
        residual = _residual_label(ylabel)
        ax.set_ylabel(f"median {residual}" + (" [%]" if unit == "" else f" [{unit}]"))
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=11, framealpha=0.85)
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _residual_label(ylabel):
    """Residual whose distribution a σ label refers to."""
    return {
        r"$\sigma_E / E$": r"$(E_{reco} - E_{gen}) / E_{gen}$",
        r"$\sigma_p / p$": r"$(|p_{reco}| - |p_{gen}|) / |p_{gen}|$",
        r"$\sigma_\theta$": r"$\theta_{reco} - \theta_{gen}$",
        r"$\sigma_\phi$": r"$\phi_{reco} - \phi_{gen}$",
    }.get(ylabel, ylabel)


def _region_slug(region):
    """File-name tag of a |cos θ| region: 'all' or e.g. 'cos0.7-0.9'."""
    if region == "all":
        return "all"
    parts = region.split()
    return f"cos{parts[0]}-{parts[-1]}"


def _plot_box(prof, edges, model, unit, ylabel, title, fname, color="black"):
    """
    Box plot of the residual per bin for one |cos θ| region, with the bias below.

    Upper panel: the line of each box is the median, the box spans the
    interquartile range and the whiskers the central 90 % (5th–95th
    percentile).  Outliers are not drawn.  The bins sit on a categorical axis,
    equally spaced.  The y range follows the boxes: it extends at most twice
    the span covered by the boxes beyond them, so that a long tail does not
    flatten everything else.  A whisker that leaves the axis is cut at the edge
    and labelled with its actual value.

    Lower panel: the same median on its own scale, i.e. its distance to zero.
    The error bar is the Gaussian estimate of the uncertainty of a median,
    1.2533 · σ / √n with σ = IQR / 1.349.
    """
    if prof is None or prof.empty:
        return
    edges = [float(e) for e in edges]
    n_bins = len(edges) - 1
    scale = 100.0 if unit == "" else 1.0
    units = " [%]" if unit == "" else f" [{unit}]"

    pos = np.array([edges.index(lo) for lo in prof["x_lo"]], dtype=float)
    med = scale * prof["bias"].to_numpy()
    q05, q25, q75, q95 = (scale * prof[c].to_numpy() for c in ("q05", "q25", "q75", "q95"))

    # Rango en y: el de las cajas más, como mucho, dos veces su extensión
    span = max(q75.max() - q25.min(), 1e-12)
    y_lo = max(q05.min(), q25.min() - 2 * span)
    y_hi = min(q95.max(), q75.max() + 2 * span)
    pad = 0.05 * (y_hi - y_lo)
    clipped_lo, clipped_hi = q05.min() < y_lo, q95.max() > y_hi
    # Hueco extra para las etiquetas de los bigotes cortados
    y_lo -= pad * (2.5 if clipped_lo else 1.0)
    y_hi += pad * (2.5 if clipped_hi else 1.0)

    fig, (ax, ax_med) = plt.subplots(
        2, 1, sharex=True, figsize=(max(8.0, 0.55 * n_bins + 2.0), 7.0),
        gridspec_kw=dict(height_ratios=[2.2, 1.0]), constrained_layout=True)

    stats = [dict(med=m, q1=a, q3=b, whislo=lo, whishi=hi, fliers=[])
             for m, a, b, lo, hi in zip(med, q25, q75, q05, q95)]
    ax.bxp(stats, positions=pos, widths=0.55, showfliers=False,
           patch_artist=True, manage_ticks=False,
           boxprops=dict(facecolor=color, alpha=0.35, edgecolor=color, lw=1.0),
           medianprops=dict(color=color, lw=2.0),
           whiskerprops=dict(color=color, lw=1.0),
           capprops=dict(color=color, lw=1.0))
    # Bigote fuera del eje: se anota su valor real en el borde
    for x, lo, hi in zip(pos, q05, q95):
        for value, edge, va, out in ((hi, y_hi, "top", hi > y_hi),
                                     (lo, y_lo, "bottom", lo < y_lo)):
            if out:
                ax.annotate(f"{value:.3g}", (x, edge),
                            xytext=(0, -2 if va == "top" else 2),
                            textcoords="offset points", fontsize=8, color=color,
                            ha="center", va=va,
                            bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    ax.axhline(0.0, color="grey", ls=":")
    ax.set_ylim(y_lo, y_hi)
    ax.set_ylabel(_residual_label(ylabel) + units)
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)
    ax.text(0.99, 0.97,
            "line: median  |  box: 25–75 %  |  whiskers: 5–95 %"
            + ("\nnumbers: whiskers beyond the axis" if clipped_lo or clipped_hi else ""),
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))

    # Panel inferior: la mediana respecto a cero, en su propia escala
    med_err = 1.2533 * (q75 - q25) / 1.349 / np.sqrt(prof["n"].to_numpy())
    ax_med.errorbar(pos, med, yerr=med_err, fmt="o-", ms=4, color=color, capsize=2)
    ax_med.axhline(0.0, color="grey", ls=":")
    ax_med.set_ylabel("median − 0" + units)
    ax_med.grid(True, alpha=0.3)
    ax_med.set_xlim(-0.5, n_bins - 0.5)
    ax_med.set_xticks(range(n_bins))
    ax_med.set_xticklabels([f"{lo:g}–{hi:g}" for lo, hi in zip(edges[:-1], edges[1:])],
                           rotation=45, ha="right")
    ax_med.set_xlabel(FIT_MODELS[model]["xlabel"])
    fig.align_ylabels([ax, ax_med])
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_dist_grid(df, res_col, x_col, edges, model, unit, ylabel, title, fname,
                    color="black", min_entries=50, n_hist_bins=60):
    """
    Grid with the residual distribution of every bin of `x_col`, one panel each.

    Every panel is drawn around its own core, median ± 4 σ with σ = IQR / 1.349
    (widened to contain zero), because the width changes by more than an order
    of magnitude along the energy range; the fraction of entries left outside
    is quoted.  The solid vertical line marks zero, the exact response, and the
    dashed one the median.  Bins with fewer than `min_entries` entries are left
    empty.
    """
    edges = [float(e) for e in edges]
    n_bins = len(edges) - 1
    scale = 100.0 if unit == "" else 1.0
    units = " [%]" if unit == "" else f" [{unit}]"
    x = df[x_col].to_numpy(dtype=float)
    r = scale * df[res_col].to_numpy(dtype=float)
    var = FIT_MODELS[model]["xlabel"].replace(" [GeV]", "")

    n_cols = min(4, n_bins)
    n_rows = int(np.ceil(n_bins / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.0 * n_cols, 2.8 * n_rows),
                             constrained_layout=True, squeeze=False)
    drawn = False
    for ax, lo, hi in zip(axes.flat, edges[:-1], edges[1:]):
        vals = r[(x >= lo) & (x < hi) & np.isfinite(r)]
        ax.set_title(f"{lo:g} ≤ {var} < {hi:g} GeV", fontsize=10)
        ax.tick_params(labelsize=8)
        if vals.size < min_entries:
            ax.text(0.5, 0.5, f"n = {vals.size:,} < {min_entries}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color="grey")
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        q25, med, q75 = np.percentile(vals, [25, 50, 75])
        half = 4.0 * max((q75 - q25) / 1.349, 1e-12)
        # El cero siempre dentro del panel, con algo de margen
        x_min = min(med - half, -0.1 * half)
        x_max = max(med + half, 0.1 * half)
        outside = float(np.mean((vals < x_min) | (vals > x_max)))
        ax.hist(vals, bins=n_hist_bins, range=(x_min, x_max), histtype="stepfilled",
                facecolor=color, alpha=0.35, edgecolor=color, lw=1.0)
        ax.axvline(0.0, color="black", lw=1.2)
        ax.axvline(med, color=color, ls="--", lw=1.5)
        ax.set_xlim(x_min, x_max)
        ax.set_yticks([])
        ax.text(0.98, 0.96,
                f"n = {vals.size:,}\nmedian = {med:+.3g}\noutside: {100 * outside:.1f} %",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5))
        drawn = True
    for ax in axes.flat[n_bins:]:
        ax.axis("off")
    if not drawn:
        plt.close(fig)
        return
    fig.supxlabel(_residual_label(ylabel) + units, fontsize=12)
    fig.supylabel("entries", fontsize=12)
    fig.suptitle(title + "\nsolid line: zero  |  dashed line: median  |  "
                 "range: median ± 4σ (σ = IQR/1.349)", fontsize=12)
    fig.savefig(fname, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _plot_vs_costheta(df, observables, cos_edges, metric, min_entries, title, fname):
    """Resolution integrated over energy as a function of |cos θ_gen|."""
    edges = np.linspace(0.0, 1.0, 21) if cos_edges is None else cos_edges
    fig, axes = plt.subplots(1, len(observables), figsize=(5 * len(observables), 4.5),
                             constrained_layout=True)
    for ax, (name, res_col, _, _, unit, ylabel) in zip(np.atleast_1d(axes), observables):
        prof = binned_profile(df, res_col, "Gen_abscos", edges, metric, min_entries)
        if prof.empty:
            continue
        scale = 100.0 if unit == "" else 1.0
        ax.errorbar(prof["x"], scale * prof["sigma"], yerr=scale * prof["sigma_err"],
                    xerr=[prof["x"] - prof["x_lo"], prof["x_hi"] - prof["x"]],
                    fmt="o", ms=4, color="black", capsize=2)
        ax.set_xlabel(r"|cos $\theta_{gen}$|")
        ax.set_ylabel(ylabel + (" [%]" if unit == "" else f" [{unit}]"))
        ax.grid(True, alpha=0.3)
    fig.suptitle(title)
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_summary_table(rows, title, fname):
    """Render the fitted parametrisations as a table, in the style of the FCC-ee baseline."""
    fig, ax = plt.subplots(figsize=(10, 0.6 + 0.55 * len(rows)))
    ax.axis("off")
    ax.text(0.0, 1.0, title, fontsize=12, weight="bold", va="top", transform=ax.transAxes)
    for i, (label, region, formula) in enumerate(rows):
        y = 1.0 - (i + 1.3) / (len(rows) + 1.3)
        ax.text(0.0, y, label, fontsize=11, va="center", transform=ax.transAxes)
        ax.text(0.33, y, region, fontsize=9, color="grey", va="center", transform=ax.transAxes)
        ax.text(0.55, y, formula, fontsize=12, va="center", transform=ax.transAxes)
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fit_detector_performance(full_df, pids, output_dir, e_bins=None, cos_bins=None,
                             metric="std90", min_entries=50, title_suffix="",
                             points_only=False):
    """
    Fit the detector-performance parametrisation for every PDG in `pids`.

    Only the correctly identified matched pairs (|Gen_pid| == |Reco_pid|) enter.
    Writes, under output_dir/<pid>/:
      - fit_<observable>.png        σ vs E/p per |cos θ| region with the fit;
      - bias_<observable>.png       median residual per bin;
      - boxplot_<observable>_<region>.png
                                    one per |cos θ| region: box plot of the
                                    residual per bin (median, 25–75 % box,
                                    5–95 % whiskers) with the median below;
      - response_energy.png         ⟨E_reco/E_true⟩ per bin;
      - resolution_vs_costheta.png  σ integrated in energy vs |cos θ|;
      - dist_<observable>_<region>.png
                                    one per |cos θ| region: grid with the
                                    residual distribution of every bin, with
                                    lines at zero and at the median;
      - profiles.csv                every binned point;
    and output_dir/fit_summary.{yaml,csv,png} with the fitted terms.

    With `points_only` no fit is performed: the σ vs E/p plot only shows the
    measured points and is written as points_<observable>.png, and no
    fit_summary.* is produced, so the outputs of a previous fit are kept.
    """
    e_bins = list(DEFAULT_FIT_E_BINS if e_bins is None else e_bins)
    cos_bins = list(DEFAULT_FIT_COS_BINS if cos_bins is None else cos_bins)
    os.makedirs(output_dir, exist_ok=True)

    regions = {"all": (0.0, 1.0)}
    for lo, hi in zip(cos_bins[:-1], cos_bins[1:]):
        regions[f"{lo:g} ≤ |cos θ| < {hi:g}"] = (lo, hi)

    summary, table_rows = [], []
    for pid in pids:
        df = matched_frame(full_df, pid)
        if df is None or df.empty:
            print(f"[resolution_fits] no correctly identified {pid_label(pid)} rows. Skipping.")
            continue
        lbl = pid_label(pid)
        pid_dir = os.path.join(output_dir, str(pid))
        os.makedirs(pid_dir, exist_ok=True)
        print(f"[resolution_fits] {lbl}: {len(df):,} matched gen=reco pairs")

        observables = _observables(pid)
        all_profiles = []
        for name, res_col, x_col, model, unit, ylabel in observables:
            profiles, fits = {}, {}
            # Mismo color por región que en el plot de bias
            box_colors = ["black", "#4363d8", "#e6194b", "#3cb44b", "#f58231", "#911eb4"]
            for color, (region, (lo, hi)) in zip(box_colors, regions.items()):
                # El último borde es inclusivo para no perder |cos θ| = 1
                sub = df.loc[(df["Gen_abscos"] >= lo)
                             & ((df["Gen_abscos"] < hi) | (hi >= 1.0))]
                prof = binned_profile(sub, res_col, x_col, e_bins, metric, min_entries)
                profiles[region] = prof
                fits[region] = None if points_only else fit_quadrature(prof, model)
                _plot_dist_grid(sub, res_col, x_col, e_bins, model, unit, ylabel,
                                f"{lbl} → {lbl}: residual distribution, "
                                f"{region}{title_suffix}",
                                os.path.join(pid_dir,
                                             f"dist_{name}_{_region_slug(region)}.png"),
                                color=color, min_entries=min_entries)
                if not prof.empty:
                    all_profiles.append(prof.assign(observable=name, region=region))

                fit = fits[region]
                summary.append(dict(
                    pid=int(pid), particle=lbl, observable=name, region=region,
                    model=model, metric=metric, n_pairs=int(len(sub)),
                    n_bins=int(len(prof)),
                    a=None if fit is None else fit["a"],
                    a_err=None if fit is None else fit["a_err"],
                    b=None if fit is None else fit["b"],
                    b_err=None if fit is None else fit["b_err"],
                    chi2=None if fit is None else fit["chi2"],
                    ndf=None if fit is None else fit["ndf"],
                    unit="relative" if unit == "" else unit,
                ))
                table_rows.append((f"{lbl}  {ylabel}", region, _formula(model, fit, unit)))

            title = f"{lbl} → {lbl}: {ylabel} ({metric}){title_suffix}"
            _plot_profiles(profiles, fits, model, unit, ylabel, title,
                           os.path.join(pid_dir, f"{'points' if points_only else 'fit'}_{name}.png"))
            _plot_bias(profiles, model, unit, ylabel, f"{lbl} → {lbl}: bias{title_suffix}",
                       os.path.join(pid_dir, f"bias_{name}.png"))
            # Un box plot por región
            for color, (region, prof) in zip(box_colors, profiles.items()):
                _plot_box(prof, e_bins, model, unit, ylabel,
                          f"{lbl} → {lbl}: residual distribution, {region}{title_suffix}",
                          os.path.join(pid_dir,
                                       f"boxplot_{name}_{_region_slug(region)}.png"),
                          color=color)
            if name in ("energy", "momentum"):
                _plot_bias(profiles, model, unit, ylabel,
                           f"{lbl} → {lbl}: energy scale{title_suffix}",
                           os.path.join(pid_dir, "response_energy.png"), response=True)

        _plot_vs_costheta(df, observables, None, metric, min_entries,
                          f"{lbl} → {lbl}: resolution vs |cos θ| ({metric}){title_suffix}",
                          os.path.join(pid_dir, "resolution_vs_costheta.png"))
        if all_profiles:
            pd.concat(all_profiles, ignore_index=True).to_csv(
                os.path.join(pid_dir, "profiles.csv"), index=False)

    # Sin ajuste no hay términos que resumir
    if not summary or points_only:
        return None

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(output_dir, "fit_summary.csv"), index=False)
    with open(os.path.join(output_dir, "fit_summary.yaml"), "w") as f:
        yaml.safe_dump(dict(e_bins=[float(e) for e in e_bins],
                            cos_bins=[float(c) for c in cos_bins],
                            metric=metric, min_entries=int(min_entries),
                            fits=summary), f, sort_keys=False, allow_unicode=True)
    _plot_summary_table(table_rows, f"full-sim detector performance ({metric}){title_suffix}",
                        os.path.join(output_dir, "fit_summary.png"))
    print(f"[resolution_fits] summary → {os.path.join(output_dir, 'fit_summary.yaml')}")
    return summary_df
