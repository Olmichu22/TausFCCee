"""Machine-readable and presentation products for ``fcc_mc_comparison_v1``."""
from __future__ import annotations

from collections import Counter
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from modules.fcc_mc_comparison import (
    ASSOCIATIONS, OUTCOMES, PID_CATEGORIES, TRUTH_SPECIES, SampleData,
    hist_with_flow, quantile_summary, write_csv,
)

COLORS = ("#2369bd", "#d85b2a")
METHOD_COLORS = {"G": "#2369bd", "Ldirect": "#d85b2a", "Lancestor": "#2a9d55",
                 "L_direct": "#d85b2a", "L_ancestor": "#2a9d55"}
SPECIES_LABEL = {"electron": "Electron", "muon": "Muon", "photon": "Photon",
                 "charged_pion": "Charged pion"}
SHORT = {"electron": "e", "muon": "mu", "photon": "gamma", "charged_pion": "pi"}
PRESENTATION_RANGES = {
    "momentum": {"electron": (-.2, .1), "muon": (-.05, .05), "photon": (-.4, .4), "charged_pion": (-.05, .05)},
    "theta_mrad": {"electron": (-3., 3.), "muon": (-2., 2.), "photon": (-4., 4.), "charged_pion": (-3., 3.)},
}


def _truth_rows(data: SampleData, species: str, tau: bool | None = None):
    rows = [row for row in data.truth if row["truth_species"] == species]
    return rows if tau is None else [row for row in rows if bool(row["tau_ancestor"]) is tau]


def _outcomes(data: SampleData, method: str, species: str, tau: bool | None = True):
    rows = [row for row in data.outcomes if row["truth_definition"] == method and row["truth_species"] == species]
    return rows if tau is None else [row for row in rows if bool(row["tau_ancestor"]) is tau]


def _pairs(data: SampleData, method: str, species: str, tau: bool | None = None):
    rows = [row for row in data.pairs if row["association_method"] == method and row["truth_species"] == species]
    return rows if tau is None else [row for row in rows if bool(row["tau_ancestor"]) is tau]


def _hist_table(samples, selector, variables, bins_by_variable):
    rows = []
    for variable in variables:
        for data in samples:
            population = selector(data)
            values = [float(row[variable]) for row in population]
            for item in hist_with_flow(values, bins_by_variable[variable], len(values)):
                rows.append({"sample": data.presentation_label, "variable": variable, **item})
    return rows


def _step_plot(path: Path, rows: list[dict], variables: tuple[str, ...], title: str, labels: dict):
    fig, axes = plt.subplots(1, len(variables), figsize=(4.6 * len(variables), 4.2), constrained_layout=True)
    if len(variables) == 1:
        axes = [axes]
    samples = list(dict.fromkeys(row["sample"] for row in rows))
    for axis, variable in zip(axes, variables):
        for color, sample in zip(COLORS, samples):
            selected = [row for row in rows if row["sample"] == sample and row["variable"] == variable and row["visible"]]
            edges = [float(selected[0]["bin_low"])] + [float(row["bin_high"]) for row in selected]
            values = [float(row["fraction"]) for row in selected]
            axis.stairs(values, edges, label=sample, color=color, linewidth=1.7)
        if variable in {"p", "pt", "truth_p"} and min(edges) > 0 and max(edges) / min(edges) > 1e4:
            axis.set_xscale("log")
        axis.set_xlabel(labels[variable]); axis.set_ylabel("Fraction / bin"); axis.grid(alpha=.2); axis.legend()
    fig.suptitle(title)
    path.parent.mkdir(parents=True, exist_ok=True); fig.savefig(path, dpi=170); plt.close(fig)


def build_part12(samples: list[SampleData], root: Path, second_token: str) -> list[Path]:
    out = root / "part12"; out.mkdir(parents=True, exist_ok=True); written = []
    suffix = f"W_vs_{second_token}"
    # Decay modes.
    decay_rows = []
    categories = sorted(set().union(*(set(data.decay_modes) for data in samples)),
                        key=lambda key: -sum(data.decay_modes[key] for data in samples))
    for data in samples:
        total = sum(data.decay_modes.values())
        for category in categories:
            count = data.decay_modes[category]
            decay_rows.append({"sample": data.presentation_label, "category": category,
                               "count": count, "fraction_of_taus": count / total})
    path = out / f"tau_decay_modes_{suffix}.csv"; write_csv(path, decay_rows); written.append(path)
    _bar_categories(out / f"tau_decay_modes_{suffix}.png", decay_rows, "category", "fraction_of_taus", "Tau decay modes")

    # Selected visible multiplicities.
    names = (("tau_origin", "electron", "e from tau decays"),
             ("tau_origin", "muon", "μ from tau decays"),
             ("tau_origin", "photon", "γ from tau decays"),
             ("tau_origin", "charged_pion", "π± from tau decays"))
    multiplicity = []
    for data in samples:
        used = sum(data.truth_counts[origin, species] for origin, species, _ in names)
        tau_total = sum(value for (origin, _), value in data.truth_counts.items() if origin == "tau_origin")
        for origin, species, label in names:
            count = data.truth_counts[origin, species]
            multiplicity.append({"sample": data.presentation_label, "category": label, "count": count,
                                 "multiplicity_per_event": count / data.expected_events})
        other = tau_total - used
        multiplicity.append({"sample": data.presentation_label, "category": "Other visible tau-decay particles",
                             "count": other, "multiplicity_per_event": other / data.expected_events})
        non_tau = data.truth_counts["non_tau_origin", "photon"]
        multiplicity.append({"sample": data.presentation_label, "category": "Non-tau γ", "count": non_tau,
                             "multiplicity_per_event": non_tau / data.expected_events})
    path = out / f"selected_truth_multiplicity_{suffix}.csv"; write_csv(path, multiplicity); written.append(path)
    _bar_categories(out / f"selected_truth_multiplicity_{suffix}.png", multiplicity, "category", "multiplicity_per_event",
                    "Selected visible truth multiplicity")

    bins_linear = {"p": np.linspace(0, 100, 101), "pt": np.linspace(0, 100, 101),
                   "theta": np.linspace(0, 180, 91)}
    labels = {"p": "p [GeV]", "pt": "pT [GeV]", "theta": "theta [deg]"}
    kine_specs = [("terminal_tau_kinematics", lambda d: d.terminal_taus, "Terminal tau kinematics")]
    for species in TRUTH_SPECIES:
        kine_specs.append((f"tau_{species}_truth_kinematics", lambda d, s=species: _truth_rows(d, s, True),
                           f"Tau-origin {SPECIES_LABEL[species].lower()} truth kinematics"))
    for stem, selector, title in kine_specs:
        rows = _hist_table(samples, selector, ("p", "pt", "theta"), bins_linear)
        path = out / f"{stem}_{suffix}.csv"; write_csv(path, rows); written.append(path)
        _step_plot(out / f"{stem}_{suffix}.png", rows, ("p", "pt", "theta"), title, labels)

    standard = {"p": np.geomspace(1e-4, 1e2, 61), "pt": np.geomspace(1e-4, 1e2, 61),
                "theta": np.linspace(0, 180, 91)}
    rows = _hist_table(samples, lambda d: _truth_rows(d, "photon", False), ("p", "pt", "theta"), standard)
    path = out / f"non_tau_photon_kinematics_{suffix}.csv"; write_csv(path, rows); written.append(path)
    _step_plot(out / f"non_tau_photon_kinematics_{suffix}.png", rows, ("p", "pt", "theta"), "Non-tau photon kinematics", labels)
    extended = {"p": np.geomspace(1e-10, 1e2, 121), "pt": np.geomspace(1e-12, 1e2, 121),
                "theta": np.linspace(0, 180, 73)}
    rows = _hist_table(samples, lambda d: _truth_rows(d, "photon", False), ("p", "pt", "theta"), extended)
    path = out / f"non_tau_photon_kinematics_extended_{suffix}.csv"; write_csv(path, rows); written.append(path)
    _step_plot(out / f"non_tau_photon_kinematics_extended_{suffix}.png", rows, ("p", "pt", "theta"),
               "Extended non-tau photon kinematics", labels)
    return written


def _bar_categories(path, rows, category, value, title):
    samples = list(dict.fromkeys(row["sample"] for row in rows)); cats = list(dict.fromkeys(row[category] for row in rows))
    x = np.arange(len(cats)); width = .8 / len(samples)
    fig, axis = plt.subplots(figsize=(max(8, .8 * len(cats)), 4.8), constrained_layout=True)
    for index, (sample, color) in enumerate(zip(samples, COLORS)):
        lookup = {row[category]: float(row[value]) for row in rows if row["sample"] == sample}
        axis.bar(x + (index - (len(samples)-1)/2) * width, [lookup.get(cat, 0) for cat in cats], width,
                 label=sample, color=color)
    axis.set_xticks(x, cats, rotation=30, ha="right"); axis.set_ylabel(value.replace("_", " "))
    axis.set_title(title); axis.grid(axis="y", alpha=.2); axis.legend(); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170); plt.close(fig)


def _integrated_rows(data: SampleData, method: str, tau: bool | None = True):
    rows = []
    for species in TRUTH_SPECIES:
        population = _outcomes(data, method, species, tau)
        counts = Counter(row["outcome"] for row in population); n = len(population)
        associated, unmatched, ambiguous = (counts[name] for name in OUTCOMES)
        missed = unmatched + ambiguous
        rows.append({"sample": data.presentation_label, "association_method": method.replace("Ldirect", "L_direct").replace("Lancestor", "L_ancestor"),
                     "species": species, "N_truth": n, "N_associated_unique": associated,
                     "N_unmatched": unmatched, "N_ambiguous": ambiguous, "N_missed": missed,
                     "efficiency_percent": "" if not n else 100 * associated / n,
                     "inefficiency_percent": "" if not n else 100 * missed / n})
    return rows


def _binned_inefficiency(data: SampleData, method: str, variable: str, tau: bool | None = True):
    bins = np.geomspace(1e-4, 1e2, 61) if variable == "p" else np.linspace(0, 180, 73)
    rows = []
    for species in TRUTH_SPECIES:
        population = _outcomes(data, method, species, tau)
        epsilon = np.finfo(float).eps * max(1, abs(float(bins[-1])))
        clipped = lambda values: np.clip(np.asarray(values, float), bins[0] + epsilon, bins[-1] - epsilon)
        all_counts, _ = np.histogram(clipped([row[variable] for row in population]), bins=bins)
        by_outcome = {name: np.histogram(clipped([row[variable] for row in population if row["outcome"] == name]), bins=bins)[0]
                      for name in OUTCOMES}
        for low, high, n, associated, unmatched, ambiguous in zip(
                bins[:-1], bins[1:], all_counts, by_outcome["associated_unique"],
                by_outcome["association_unmatched"], by_outcome["ambiguous_multiple_pfo"]):
            missed = int(unmatched + ambiguous)
            rows.append({"sample": data.presentation_label,
                "association_method": method.replace("Ldirect", "L_direct").replace("Lancestor", "L_ancestor"),
                "species": species, "variable": variable, "bin_low": low, "bin_high": high,
                "N_truth": int(n), "N_associated_unique": int(associated), "N_unmatched": int(unmatched),
                "N_ambiguous": int(ambiguous), "N_missed": missed,
                "inefficiency_percent": "" if n == 0 else 100 * missed / n,
                "drawn": bool(n and missed > 0)})
    return rows


def build_part3(samples: list[SampleData], root: Path, second_token: str,
                first_token: str = "W", truth_tau: bool | None = True) -> list[Path]:
    out = root / "part3"; (out / "method_comparison").mkdir(parents=True, exist_ok=True)
    (out / "non_tau_backup").mkdir(parents=True, exist_ok=True); written = []
    tokens = (first_token, second_token)
    for data, token in zip(samples, tokens):
        all_integrated = []
        for method in ASSOCIATIONS:
            integrated = _integrated_rows(data, method, truth_tau)
            all_integrated.extend(integrated)
            path = out / f"{token}_{method}_reco_efficiency.csv"; write_csv(path, integrated); written.append(path)
            _bar_categories(out / f"{token}_{method}_reco_efficiency.png", integrated, "species", "efficiency_percent",
                            f"{data.presentation_label}: {method} association efficiency")
            for variable in ("p", "theta"):
                binned = _binned_inefficiency(data, method, variable, truth_tau)
                path = out / f"{token}_{method}_reco_inefficiency_vs_{variable}.csv"; write_csv(path, binned); written.append(path)
                _line_inefficiency(out / f"{token}_{method}_reco_inefficiency_vs_{variable}.png", binned, variable,
                                   f"{data.presentation_label}: {method} inefficiency")
        label_token = "WHIZARD" if token == "W" else token
        method_integrated = [{**row, "association_method": row["association_method"].replace("L_direct", "Ldirect").replace("L_ancestor", "Lancestor")}
                             for row in all_integrated]
        path = out / "method_comparison" / f"{label_token}_reco_efficiency_G_Ldirect_Lancestor.csv"
        write_csv(path, method_integrated); written.append(path)
        _method_efficiency(out / "method_comparison" / f"{label_token}_reco_efficiency_G_Ldirect_Lancestor.png",
                           all_integrated, data.presentation_label)
        for variable in ("p", "theta"):
            combined = [{**row, "association_method": row["association_method"].replace("L_direct", "Ldirect").replace("L_ancestor", "Lancestor")}
                        for method in ASSOCIATIONS for row in _binned_inefficiency(data, method, variable, truth_tau)]
            for row in combined:
                row["efficiency_percent"] = "" if not row["N_truth"] else 100 * row["N_associated_unique"] / row["N_truth"]
                row["drawn_on_log_plot"] = row["drawn"]
            path = out / "method_comparison" / f"{label_token}_reco_inefficiency_vs_{variable}_G_Ldirect_Lancestor.csv"
            write_csv(path, combined); written.append(path)
            _line_inefficiency(out / "method_comparison" / f"{label_token}_reco_inefficiency_vs_{variable}_G_Ldirect_Lancestor.png",
                               combined, variable, f"{data.presentation_label}: G/L_direct/L_ancestor")
    # Two-sample Lancestor and non-tau photon comparisons.
    combined = [row for data in samples for row in _integrated_rows(data, "Lancestor", truth_tau)]
    stem = f"{first_token}_vs_{second_token}_Lancestor_reco_efficiency"
    path = out / f"{stem}.csv"; write_csv(path, combined); written.append(path)
    _bar_categories(out / f"{stem}.png", combined, "species", "efficiency_percent", "L_ancestor efficiency")
    for variable in ("p", "theta"):
        rows = [row for data in samples for row in _binned_inefficiency(data, "Lancestor", variable, truth_tau)]
        path = out / f"{first_token}_vs_{second_token}_Lancestor_reco_inefficiency_vs_{variable}.csv"; write_csv(path, rows); written.append(path)
        _line_inefficiency(path.with_suffix(".png"), rows, variable, "L_ancestor inefficiency")
    for method in ASSOCIATIONS:
        rows = [row for data in samples for row in _integrated_rows(data, method, False) if row["species"] == "photon"]
        directory = out if method == "Lancestor" else out / "non_tau_backup"
        stem = f"{first_token}_vs_{second_token}_{method}_non_tau_photon_reco_efficiency"
        path = directory / f"{stem}.csv"; write_csv(path, rows); written.append(path)
        _bar_categories(path.with_suffix(".png"), rows, "species", "efficiency_percent", "Non-tau photon efficiency")
        for variable in ("p", "theta"):
            binned = [row for data in samples for row in _binned_inefficiency(data, method, variable, False)
                      if row["species"] == "photon"]
            path = directory / f"{first_token}_vs_{second_token}_{method}_non_tau_photon_reco_inefficiency_vs_{variable}.csv"
            write_csv(path, binned); written.append(path); _line_inefficiency(path.with_suffix(".png"), binned, variable,
                                                                              "Non-tau photon inefficiency")
    return written


def _line_inefficiency(path, rows, variable, title):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), constrained_layout=True)
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    methods = list(dict.fromkeys(r["association_method"] for r in rows))
    by_method = len(samples) == 1 and len(methods) > 1
    series = methods if by_method else samples
    for axis, species in zip(axes.flat, TRUTH_SPECIES):
        for index, name in enumerate(series):
            selected = [r for r in rows if r["species"] == species and str(r["drawn"]).lower() == "true" and
                        r["association_method" if by_method else "sample"] == name]
            axis.plot([math.sqrt(float(r["bin_low"])*float(r["bin_high"])) if variable == "p" else (float(r["bin_low"])+float(r["bin_high"]))/2 for r in selected],
                      [float(r["inefficiency_percent"]) for r in selected], marker="o", ms=2, lw=1,
                      color=METHOD_COLORS.get(name, COLORS[index % len(COLORS)]), label=name)
        if variable == "p": axis.set_xscale("log")
        if any(r["species"] == species and str(r["drawn"]).lower() == "true" for r in rows):
            axis.set_yscale("log")
        axis.set_title(SPECIES_LABEL[species]); axis.grid(alpha=.2)
        if any(r["species"] == species and str(r["drawn"]).lower() == "true" for r in rows): axis.legend()
    fig.suptitle(title); path.parent.mkdir(parents=True, exist_ok=True); fig.savefig(path, dpi=170); plt.close(fig)


def _method_efficiency(path, rows, label):
    x=np.arange(len(TRUTH_SPECIES)); width=.24
    fig,axis=plt.subplots(figsize=(8.5,4.8),constrained_layout=True)
    for index,method in enumerate(("G","Ldirect","Lancestor")):
        lookup={row["species"]:float(row["efficiency_percent"]) for row in rows
                if row["association_method"].replace("_", "") == method}
        axis.bar(x+(index-1)*width,[lookup[s] for s in TRUTH_SPECIES],width,label=method,color=METHOD_COLORS[method])
    axis.set_xticks(x,[SPECIES_LABEL[s] for s in TRUTH_SPECIES]);axis.set_ylabel("Efficiency [%]")
    axis.set_title(f"{label}: association methods");axis.grid(axis="y",alpha=.2);axis.legend()
    path.parent.mkdir(parents=True,exist_ok=True);fig.savefig(path,dpi=170);plt.close(fig)


def _residual(pair, residual):
    if residual == "momentum": return (float(pair["reco_p"]) - float(pair["truth_p"])) / float(pair["truth_p"])
    return (float(pair["reco_theta"]) - float(pair["truth_theta"])) * math.pi / 180 * 1000


def build_part3b(samples: list[SampleData], root: Path, second_token: str,
                 first_token: str = "W") -> list[Path]:
    out = root / "part3b"; out.mkdir(parents=True, exist_ok=True); written = []; suffix = f"{first_token}_vs_{second_token}"
    summaries_by_key = {}
    for residual in ("momentum", "theta_mrad"):
        rows = []
        for species in TRUTH_SPECIES:
            for data in samples:
                population = [row for row in _pairs(data, "Lancestor", species) if not row.get("aggregate_pid_only")]
                values = [_residual(row, residual) for row in population]
                low, high = PRESENTATION_RANGES[residual][species]; bins = np.linspace(low, high, 121)
                for item in hist_with_flow(values, bins, len(values)):
                    rows.append({"sample": data.presentation_label, "species": species, "residual": residual, **item})
                summary = quantile_summary(values)
                outside = np.mean((np.asarray(values) < low) | (np.asarray(values) > high))
                summaries_by_key[data.presentation_label, species, residual] = {"sample": data.presentation_label, "species": species, "residual": residual,
                                  **summary, "visible_xmin": low, "visible_xmax": high,
                                  "fraction_outside_visible_range": float(outside)}
        stem = f"{suffix}_pfo_{'momentum' if residual == 'momentum' else 'theta'}_residuals"
        path = out / f"{stem}.csv"; write_csv(path, rows); written.append(path)
        _residual_plot(path.with_suffix(".png"), rows, residual, samples)
    summaries = [summaries_by_key[data.presentation_label, species, residual]
                 for data in samples for species in TRUTH_SPECIES for residual in ("momentum", "theta_mrad")]
    path = out / f"{suffix}_pfo_residual_summary.csv"; write_csv(path, summaries); written.append(path)
    # Photon truth-origin presentation product for both samples.
    for data, token in zip(samples, (first_token, second_token)):
        rows = []
        for residual, bins in (("momentum", np.linspace(-1, 1.1, 121)),
                               ("theta_mrad", np.linspace(-40, 40, 121))):
            for tau, origin in ((True, "tau-origin photons"), (False, "non-tau photons")):
                population = [r for r in _pairs(data, "Lancestor", "photon", tau) if not r.get("aggregate_pid_only")]
                values = [_residual(row, residual) for row in population]
                for item in hist_with_flow(values, bins, len(values)):
                    visible = item.pop("visible")
                    rows.append({"sample": data.internal_name, "species": "photon", "origin": origin,
                                 "residual": residual, **item, "visible_in_main_plot": visible})
        stem = f"{token}_photon_pfo_residuals_by_truth_origin"
        path = out / f"{stem}.csv"; write_csv(path, rows); written.append(path)
        _photon_origin_plot(path.with_suffix(".png"), rows, data.presentation_label)
    return written


def build_photon_diagnostic(samples: list[SampleData], root: Path,
                            first_token: str, second_token: str) -> list[Path]:
    """Kinematic diagnostics for already-associated L_ancestor photons."""
    out = root / "part3b/photon_diagnostic"; out.mkdir(parents=True, exist_ok=True)
    written = []; pairs = {data.presentation_label: _pairs(data, "Lancestor", "photon", None) for data in samples}
    all_p = np.asarray([float(row["truth_p"]) for rows in pairs.values() for row in rows])
    p_bins = np.geomspace(float(all_p.min()), float(all_p.max()), 61)
    theta_bins = np.linspace(0, 180, 73)
    hist_rows = []
    for data in samples:
        rows = pairs[data.presentation_label]
        for variable, bins in (("truth_p", p_bins), ("truth_theta", theta_bins)):
            values = [float(row[variable]) for row in rows]
            for item in hist_with_flow(values, bins, len(values)):
                hist_rows.append({"sample": data.presentation_label, "variable": variable, **item})
    stem = f"{first_token}_vs_{second_token}_associated_photon_truth_kinematics"
    path = out / f"{stem}.csv"; write_csv(path, hist_rows); written.append(path)
    _step_plot(path.with_suffix(".png"), hist_rows, ("truth_p", "truth_theta"),
               "L_ancestor associated photon truth kinematics",
               {"truth_p": "p_truth [GeV]", "truth_theta": "theta_truth [deg]"})

    y_bins = np.linspace(-100, 100, 201); matrices = {}; count_rows = []; summary = []
    for data in samples:
        rows = pairs[data.presentation_label]
        p_values = np.asarray([float(row["truth_p"]) for row in rows])
        residuals = np.asarray([_residual(row, "theta_mrad") for row in rows])
        matrix, _, _ = np.histogram2d(p_values, residuals, bins=(p_bins, y_bins)); matrices[data.presentation_label] = matrix
        for ix in range(len(p_bins)-1):
            for iy in range(len(y_bins)-1):
                count_rows.append({"sample": data.presentation_label, "p_low": p_bins[ix], "p_high": p_bins[ix+1],
                                   "delta_theta_low_mrad": y_bins[iy], "delta_theta_high_mrad": y_bins[iy+1],
                                   "count": int(matrix[ix, iy])})
        for origin, selected in (("inclusive", rows),
                                 ("tau-origin", [r for r in rows if r["tau_ancestor"]]),
                                 ("non-tau", [r for r in rows if not r["tau_ancestor"]]),
                                 ("parentless", [r for r in rows if r.get("parentless")])):
            values = [_residual(row, "theta_mrad") for row in selected]
            item = quantile_summary(values) if values else {"N": 0, "median": "", "q16": "", "q84": "", "central68_halfwidth": ""}
            summary.append({"sample": data.presentation_label, "origin": origin, **item})
    path2d = out / f"{first_token}_vs_{second_token}_photon_delta_theta_vs_truth_p_2d.csv"
    write_csv(path2d, count_rows); written.append(path2d)
    positive = np.concatenate([m[m > 0] for m in matrices.values()]); vmax = float(positive.max())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    image = None
    for axis, data in zip(axes, samples):
        image = axis.pcolormesh(p_bins, y_bins, matrices[data.presentation_label].T,
                                norm=LogNorm(vmin=1, vmax=vmax), cmap="viridis")
        axis.set_xscale("log"); axis.set_xlabel("p_truth [GeV]"); axis.set_ylabel("theta_PFO - theta_truth [mrad]")
        axis.set_title(data.presentation_label)
    fig.colorbar(image, ax=axes, label="2D bin population")
    fig.suptitle("Associated photon polar-angle residual versus truth momentum")
    fig.savefig(path2d.with_suffix(".png"), dpi=170); plt.close(fig)

    assoc_rows = []
    for data in samples:
        values = [float(row["truth_p"]) for row in pairs[data.presentation_label]]
        counts, _ = np.histogram(values, p_bins)
        for low, high, count in zip(p_bins[:-1], p_bins[1:], counts):
            assoc_rows.append({"sample": data.presentation_label, "p_low": low, "p_high": high,
                               "associated_count": int(count), "associated_per_event": count / data.expected_events})
    path_counts = out / f"{first_token}_vs_{second_token}_associated_photon_counts_vs_truth_p.csv"
    write_csv(path_counts, assoc_rows); written.append(path_counts)
    fig, axis = plt.subplots(figsize=(7.5, 4.8), constrained_layout=True)
    for color, data in zip(COLORS, samples):
        selected = [r for r in assoc_rows if r["sample"] == data.presentation_label]
        edges = [selected[0]["p_low"]] + [r["p_high"] for r in selected]
        axis.stairs([r["associated_per_event"] for r in selected], edges, label=data.presentation_label, color=color)
    axis.set_xscale("log"); axis.set_xlabel("p_truth [GeV]"); axis.set_ylabel("Associated photons / event / log bin")
    axis.grid(alpha=.2); axis.legend(); axis.set_title("L_ancestor associated photons versus truth momentum")
    fig.savefig(path_counts.with_suffix(".png"), dpi=170); plt.close(fig)
    path_summary = out / f"{first_token}_vs_{second_token}_photon_diagnostic_summary.csv"
    write_csv(path_summary, summary); written.append(path_summary)
    return written


def _residual_plot(path, rows, residual, samples):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), constrained_layout=True)
    for axis, species in zip(axes.flat, TRUTH_SPECIES):
        for color, data in zip(COLORS, samples):
            selected = [r for r in rows if r["sample"] == data.presentation_label and r["species"] == species and r["visible"]]
            edges = [selected[0]["bin_low"]] + [r["bin_high"] for r in selected]
            axis.stairs([r["fraction"] for r in selected], edges, color=color, label=data.presentation_label)
        axis.axvline(0, color="0.4", lw=.7); axis.set_title(SPECIES_LABEL[species]); axis.grid(alpha=.2)
    fig.legend(loc="upper center", ncol=2); fig.suptitle("PFO momentum residuals" if residual == "momentum" else "PFO polar-angle residuals")
    path.parent.mkdir(parents=True, exist_ok=True); fig.savefig(path, dpi=170); plt.close(fig)


def _photon_origin_plot(path, rows, label):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), constrained_layout=True)
    for axis, residual in zip(axes, ("momentum", "theta_mrad")):
        for color, origin in zip(COLORS, ("tau-origin photons", "non-tau photons")):
            selected = [r for r in rows if r["origin"] == origin and r["residual"] == residual and r["visible_in_main_plot"]]
            edges = [selected[0]["bin_low"]] + [r["bin_high"] for r in selected]
            axis.stairs([r["fraction"] for r in selected], edges, label=f"{origin} (N={sum(r['count'] for r in rows if r['origin']==origin and r['residual']==residual)})", color=color)
        axis.axvline(0, color="0.4", lw=.7); axis.grid(alpha=.2); axis.legend()
    fig.suptitle(f"Photon PFO residuals by truth origin\n({label})"); fig.savefig(path, dpi=170); plt.close(fig)


def _pid_rows(data: SampleData, method: str, tau: bool | None = True):
    rows = []
    for species in TRUTH_SPECIES:
        population = _pairs(data, method, species, tau)
        counts = Counter(row["reco_category"] for row in population); denominator = len(population)
        for reco in PID_CATEGORIES:
            rows.append({"sample": data.presentation_label, "association_method": method,
                         "truth_species": species, "reco_pid": reco, "count": counts[reco],
                         "conditional_fraction_percent": 100 * counts[reco] / denominator,
                         "conditional_denominator": denominator})
    return rows


def build_part4(samples: list[SampleData], root: Path, second_token: str) -> list[Path]:
    out = root / "part4"; method_out = out / "method_comparison"; method_out.mkdir(parents=True, exist_ok=True); written=[]
    tokens = ("W", second_token)
    correct = []
    for data, token in zip(samples, tokens):
        all_methods = []
        for method in ASSOCIATIONS:
            rows = _pid_rows(data, method)
            all_methods.extend(rows)
            path = method_out / f"{token}_{method}_conditional_pid_matrix.csv"; write_csv(path, rows); written.append(path)
            _pid_matrix(path.with_suffix(".png"), rows, f"{data.presentation_label}: {method}")
            if method == "Lancestor":
                for species in TRUTH_SPECIES:
                    diagonal = next(r for r in rows if r["truth_species"] == species and r["reco_pid"] == species)
                    correct.append({"sample": data.presentation_label, "truth_species": species,
                                    "correct_PID": diagonal["count"],
                                    "associated_unique_denominator": diagonal["conditional_denominator"],
                                    "conditional_PID_efficiency_percent": diagonal["conditional_fraction_percent"]})
        label_token = "WHIZARD" if token == "W" else token
        path = method_out / f"{label_token}_conditional_pid_matrices_G_Ldirect_Lancestor.csv"
        write_csv(path, all_methods); written.append(path)
        _pid_composite(path.with_suffix(".png"), all_methods, data.presentation_label)
        lrows = [r for r in all_methods if r["association_method"] == "Lancestor"]
        simple = [{"sample": r["sample"], "truth_species": r["truth_species"], "reco_pid": r["reco_pid"],
                   "count": r["count"], "associated_unique_denominator": r["conditional_denominator"],
                   "fraction_percent": r["conditional_fraction_percent"]} for r in lrows]
        path = out / f"{token}_conditional_pid_matrix.csv"; write_csv(path, simple); written.append(path)
        _pid_matrix(path.with_suffix(".png"), lrows, f"{data.presentation_label}: L_ancestor conditional PID")
        non_tau = []
        for method in ASSOCIATIONS:
            population = _pairs(data, method, "photon", False); counts=Counter(r["reco_category"] for r in population); n=len(population)
            for reco in PID_CATEGORIES:
                non_tau.append({"sample": data.presentation_label, "association_method": method,
                                "truth_origin": "non_tau_photon", "reco_pid": reco, "count": counts[reco],
                                "conditional_fraction_percent": 100*counts[reco]/n, "conditional_denominator": n})
        path = method_out / f"{label_token}_non_tau_photon_pid_G_Ldirect_Lancestor.csv"
        write_csv(path, non_tau); written.append(path); _pid_non_tau(path.with_suffix(".png"), non_tau, data.presentation_label)
    path = out / f"W_vs_{second_token}_conditional_pid_efficiency.csv"; write_csv(path, correct); written.append(path)
    _bar_categories(path.with_suffix(".png"), correct, "truth_species", "conditional_PID_efficiency_percent", "Conditional correct-PID efficiency")
    coverage = [row for data in samples for row in data.pfo_coverage]
    path = method_out / f"W_vs_{second_token}_truth_association_coverage.csv"; write_csv(path, coverage); written.append(path)
    _bar_categories(path.with_suffix(".png"), coverage, "association_method", "coverage_percent", "PFO-to-selected-truth coverage")
    return written


def _pid_matrix(path, rows, title):
    matrix=np.zeros((len(TRUTH_SPECIES),len(PID_CATEGORIES)))
    for row in rows: matrix[TRUTH_SPECIES.index(row["truth_species"]),PID_CATEGORIES.index(row["reco_pid"])]=float(row["conditional_fraction_percent"])
    fig,axis=plt.subplots(figsize=(8,4.5),constrained_layout=True); image=axis.imshow(matrix,vmin=0,vmax=100,cmap="Blues",aspect="auto")
    axis.set_xticks(range(len(PID_CATEGORIES)),PID_CATEGORIES,rotation=35,ha="right");axis.set_yticks(range(len(TRUTH_SPECIES)),TRUTH_SPECIES)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]): axis.text(j,i,f"{matrix[i,j]:.1f}",ha="center",va="center",fontsize=8)
    axis.set_title(title);fig.colorbar(image,ax=axis,label="Conditional fraction [%]");fig.savefig(path,dpi=170);plt.close(fig)


def _pid_composite(path, rows, label):
    fig,axes=plt.subplots(1,3,figsize=(17,4.5),constrained_layout=True)
    image=None
    for axis,method in zip(axes,ASSOCIATIONS):
        matrix=np.zeros((len(TRUTH_SPECIES),len(PID_CATEGORIES)))
        for row in rows:
            if row["association_method"]==method: matrix[TRUTH_SPECIES.index(row["truth_species"]),PID_CATEGORIES.index(row["reco_pid"])]=float(row["conditional_fraction_percent"])
        image=axis.imshow(matrix,vmin=0,vmax=100,cmap="Blues",aspect="auto");axis.set_title(method)
        axis.set_xticks(range(len(PID_CATEGORIES)),PID_CATEGORIES,rotation=35,ha="right");axis.set_yticks(range(len(TRUTH_SPECIES)),TRUTH_SPECIES)
    fig.colorbar(image,ax=axes,label="Conditional fraction [%]");fig.suptitle(f"{label}: conditional PID");fig.savefig(path,dpi=170);plt.close(fig)


def _pid_non_tau(path, rows, label):
    _bar_categories(path, rows, "reco_pid", "conditional_fraction_percent", f"{label}: non-tau photon conditional PID")


def build_all(samples: list[SampleData], root: Path, second_token: str) -> list[Path]:
    paths=[]
    paths.extend(build_part12(samples,root,second_token)); paths.extend(build_part3(samples,root,second_token))
    paths.extend(build_part3b(samples,root,second_token)); paths.extend(build_part4(samples,root,second_token))
    return paths
