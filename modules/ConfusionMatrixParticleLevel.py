import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
import itertools
from collections import defaultdict


# ─────────────────────────────────────────────
#  Helpers
# ─────────────────────────────────────────────

def parse_association_results(association_results_df):
    """
    Parse the flat dict  {  "<|gen_pid|>_<|reco_pid|>_<interval>": count, ... }
    into a nested structure:
        { energy_bin_str: { (gen_pid, reco_pid): count } }
    """
    parsed = defaultdict(dict)
    for key, count in association_results_df.items():
        parts = key.split("_", 2)          # split on the first two underscores only
        if len(parts) != 3:
            continue
        gen_pid, reco_pid, energy_bin = parts
        parsed[energy_bin][(int(gen_pid), int(reco_pid))] = count
    return parsed


def build_matrix(bin_data, gen_ids, reco_ids):
    """
    Build a 2-D numpy array (len(gen_ids) x len(reco_ids)) from the dict
    { (gen_pid, reco_pid): count }.
    """
    mat = np.zeros((len(gen_ids), len(reco_ids)), dtype=float)
    for (gid, rid), count in bin_data.items():
        if gid in gen_ids and rid in reco_ids:
            i = gen_ids.index(gid)
            j = reco_ids.index(rid)
            mat[i, j] = count
    return mat


def safe_normalise_rows(mat):
    """Normalise each ROW (gen) → efficiency."""
    row_sums = mat.sum(axis=1, keepdims=True)
    return np.divide(mat, row_sums, where=row_sums != 0, out=np.zeros_like(mat))


def safe_normalise_cols(mat):
    """Normalise each COLUMN (reco) → purity."""
    col_sums = mat.sum(axis=0, keepdims=True)
    return np.divide(mat, col_sums, where=col_sums != 0, out=np.zeros_like(mat))


def pid_label(pid):
    """Human-readable label for a PDG particle ID (absolute value expected)."""
    pid_map = {
        11:   "e±",
        13:   "μ±",
        15:   "τ±",
        22:   "γ",
        111:  "π⁰",
        211:  "π±",
        130:  "K⁰L",
        310:  "K⁰S",
        321:  "K±",
        2112: "n",
        2212: "p",
        3122: "Λ",
        3112: "Σ⁻",
        3222: "Σ⁺",
        3312: "Ξ⁻",
        3322: "Ξ⁰",
    }
    return pid_map.get(pid, str(pid))


def _energy_bin_center(energy_bin):
    parts = energy_bin.strip("()[]").split(",")
    if len(parts) == 2:
        try:
            lo = float(parts[0])
            hi = float(parts[1])
            if not np.isfinite(hi):
                return lo * 1.5   # open-ended bin: use 1.5× the lower edge
            return (lo + hi) / 2
        except ValueError:
            return 0.5
    return float(parts[0])


def _extract_true_reco_pairs(samples):
    true_energies = []
    residuals = []

    for sample in samples:
        if isinstance(sample, dict):
            true_energy = sample.get("Gen_energy", sample.get("true_energy"))
            reco_energy = sample.get("Reco_energy", sample.get("reco_energy"))
        elif isinstance(sample, (tuple, list)) and len(sample) >= 2:
            true_energy, reco_energy = sample[0], sample[1]
        else:
            continue

        if true_energy is None or reco_energy is None:
            continue

        true_energy = float(true_energy)
        reco_energy = float(reco_energy)
        if not np.isfinite(true_energy) or not np.isfinite(reco_energy):
            continue
        if true_energy == 0:
            continue

        true_energies.append(true_energy)
        residuals.append((reco_energy - true_energy) / true_energy)

    return np.asarray(true_energies, dtype=float), np.asarray(residuals, dtype=float)


# def _std90(values):
#     values = np.asarray(values, dtype=float)
#     if values.size == 0:
#         return None
#     low, high = np.percentile(values, [5, 95])
#     central = values[(values >= low) & (values <= high)]
#     if central.size == 0:
#         return None
#     return float(np.std(central))

def _std90(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return None
    
    x = np.sort(values)
    n_low = int(len(x) * 0.1)
    n_quant = len(x) - n_low
    
    if n_quant == len(x):
        return float(np.std(x))
    
    distances = x[n_quant:] - x[:n_low]
    i_start = np.argmin(distances)
    central = x[i_start:i_start + n_quant]
    
    return float(np.std(central))
def _resolution_value(residuals, mode):
    residuals = np.asarray(residuals, dtype=float)
    residuals = residuals[np.isfinite(residuals)]
    if residuals.size == 0:
        return None

    if mode == "std":
        return float(np.std(residuals))
    if mode == "iqr":
        p16, p84 = np.percentile(residuals, [16, 84])
        return float((p84 - p16) / 2.0)
    if mode == "std90":
        return _std90(residuals)
    raise ValueError(f"Unknown resolution mode: {mode}")


def _resolution_error(residuals, mode):
    """Statistical uncertainty on the resolution metric.

    std   → σ / √(2(n−1))   exact for Gaussian, good approximation otherwise
    iqr   → value / √n       rough approximation
    std90 → σ₉₀ / √(2(n₉₀−1)) same formula applied to the central-90 % subset
    """
    residuals = np.asarray(residuals, dtype=float)
    residuals = residuals[np.isfinite(residuals)]
    n = residuals.size
    if n < 2:
        return None

    if mode == "std":
        sigma = float(np.std(residuals))
        return sigma / np.sqrt(2 * (n - 1))
    if mode == "iqr":
        p16, p84 = np.percentile(residuals, [16, 84])
        value = (p84 - p16) / 2.0
        return float(value / np.sqrt(n))
    if mode == "std90":
        low, high = np.percentile(residuals, [5, 95])
        central = residuals[(residuals >= low) & (residuals <= high)]
        n90 = central.size
        if n90 < 2:
            return None
        sigma90 = float(np.std(central))
        return sigma90 / np.sqrt(2 * (n90 - 1))
    raise ValueError(f"Unknown resolution mode: {mode}")


def _plot_resolution_curve(series_by_bin, energy_bins, output_dir, filename, title, ylabel, label, metric_mode, dpi, log_x=True):
    data_to_plot = []
    for energy_bin in energy_bins:
        residuals = series_by_bin.get(energy_bin, [])
        value = _resolution_value(residuals, metric_mode)
        if value is None:
            continue
        err = _resolution_error(residuals, metric_mode)
        data_to_plot.append((energy_bin, value, err))

    if not data_to_plot:
        return False

    bins, resolution_values, errors = zip(*data_to_plot)
    bin_centers = [_energy_bin_center(b) for b in bins]
    yerr = [e if e is not None else 0.0 for e in errors]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(
        bin_centers,
        resolution_values,
        yerr=yerr,
        fmt="o-",
        capsize=5,
        linewidth=1.0,
        markersize=5,
        label=label,
    )
    if log_x:
        ax.set_xscale("log")
    else:
        ax.set_xticks(bin_centers)
        ax.set_xticklabels([str(round(b, 1)) for b in bin_centers], rotation=45, ha="right")
    ax.set_xlabel("True energy bin (GeV)", loc="left")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both" if log_x else "major")
    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    # Ticks on both y-axes, with minor ticks
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="y", which="both", left=True, right=True)
    ax.tick_params(axis="y", which="major", length=6, width=1.2)
    ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    plt.tight_layout()
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=dpi)
    print(f"Saved energy resolution plot → {output_path}")
    plt.close()
    return True


def _plot_combined_resolution_curves(series_by_bin, energy_bins, output_dir, filename, title, dpi, log_x=True):
    metrics = [
        ("std",   "std",        "#0055CC", "σ((E_reco - E_true) / E_true)"),
        ("iqr",   "IQR84-16/2", "#E84800", "IQR84-16 / 2 of ((E_reco - E_true) / E_true)"),
        ("std90", "std90",      "#008A00", "std90 of ((E_reco - E_true) / E_true)"),
    ]

    fig, ax = plt.subplots(figsize=(8, 6))
    plotted_any = False
    for metric_mode, label, color, _ in metrics:
        data_to_plot = []
        for energy_bin in energy_bins:
            residuals = series_by_bin.get(energy_bin, [])
            value = _resolution_value(residuals, metric_mode)
            if value is None:
                continue
            err = _resolution_error(residuals, metric_mode)
            data_to_plot.append((energy_bin, value, err))

        if not data_to_plot:
            continue

        bins, resolution_values, errors = zip(*data_to_plot)
        bin_centers = [_energy_bin_center(b) for b in bins]
        yerr = [e if e is not None else 0.0 for e in errors]

        ax.errorbar(
            bin_centers,
            resolution_values,
            yerr=yerr,
            fmt="o-",
            capsize=4,
            linewidth=1.5,
            markersize=4,
            color=color,
            label=label,
        )
        plotted_any = True

    if not plotted_any:
        plt.close()
        return False

    if log_x:
        ax.set_xscale("log")
    else:
        tick_bins = sorted(energy_bins, key=_energy_bin_center)
        tick_positions = [_energy_bin_center(b) for b in tick_bins]
        tick_labels = [str(round(pos, 1)) for pos in tick_positions]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=45, ha="right")
    ax.set_ylim(bottom=0, top=0.5)
    ax.set_xlabel("True energy bin (GeV)", loc="left")
    ax.set_ylabel("Resolution of ((E_reco - E_true) / E_true)")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both" if log_x else "major")
    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    # Ticks on both y-axes, with minor ticks
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="y", which="both", left=True, right=True)
    ax.tick_params(axis="y", which="major", length=6, width=1.2)
    ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    plt.tight_layout()
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=dpi)
    print(f"Saved combined energy resolution plot → {output_path}")
    plt.close()
    return True


# ─────────────────────────────────────────────
#  Core plotting function
# ─────────────────────────────────────────────

def plot_energy_distributions(energy_distribution_results, output_dir=".", dpi=150, all_plot_pdgs=None):
    """
    Plot resolution versus true energy for each migration type.
    The dict is expected to have the structure:
    {"migrationtype_energybin": [(E_true, E_reco), ...], ...}
    where "migrationtype" is a string like "key = str(abs(gen_pid)) + "_" + str(abs(reco_pid))"
    and "energybin" is a string from pd.cut binning on E_true (e.g. "(0, 1]", "(1, 5]", etc.).
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    list_keys = sorted(energy_distribution_results.keys())
    migration_types = sorted(set(k.rsplit("_", 1)[0] for k in list_keys))
    energy_bins = sorted(
        set(k.rsplit("_", 1)[1] for k in list_keys),
        key=lambda x: float(x.strip("()[]").split(",")[0]),
    )
    if all_plot_pdgs is None:
        all_plot_pdgs = [22]
    _all_pdgs_str = [str(int(p)) for p in all_plot_pdgs]
    aggregated_results      = {p: defaultdict(list) for p in _all_pdgs_str}
    aggregated_true_results = {p: defaultdict(list) for p in _all_pdgs_str}
    for key, energies in energy_distribution_results.items():
        migration, energy_bin = key.rsplit("_", 1)
        migration_parts = migration.split("_", 1)
        if len(migration_parts) != 2:
            continue
        gen_pid_str, reco_pid = migration_parts
        for pdg_str in _all_pdgs_str:
            if gen_pid_str != "999" and reco_pid == pdg_str:
                aggregated_results[pdg_str][energy_bin].extend(energies)
            if reco_pid != "999" and gen_pid_str == pdg_str:
                aggregated_true_results[pdg_str][energy_bin].extend(energies)

    for migration in migration_types:
        gen_pid_str, reco_pid_str = migration.split("_", 1)
        if "999" in (gen_pid_str, reco_pid_str):
            continue
        if gen_pid_str != reco_pid_str and reco_pid_str != "22" and gen_pid_str != "22":
            continue

        series_by_bin = {}
        for egen_bin in energy_bins:
            key = f"{migration}_{egen_bin}"
            pairs = energy_distribution_results.get(key, [])
            _, residuals = _extract_true_reco_pairs(pairs)
            if residuals.size == 0:
                continue
            series_by_bin[egen_bin] = residuals

        if not series_by_bin:
            continue

        _plot_resolution_curve(
            series_by_bin,
            energy_bins,
            output_dir,
            f"residual_resolution_std_{migration}.png",
            f"Energy resolution (std) for migration: {migration}",
            "σ((E_reco - E_true) / E_true)",
            migration,
            "std",
            dpi,
        )
        _plot_resolution_curve(
            series_by_bin,
            energy_bins,
            output_dir,
            f"residual_resolution_iqr84_16_{migration}.png",
            f"Energy resolution (IQR84-16/2) for migration: {migration}",
            "IQR84-16 / 2 of ((E_reco - E_true) / E_true)",
            migration,
            "iqr",
            dpi,
        )
        _plot_resolution_curve(
            series_by_bin,
            energy_bins,
            output_dir,
            f"residual_resolution_std90_{migration}.png",
            f"Energy resolution (std90) for migration: {migration}",
            "std90 of ((E_reco - E_true) / E_true)",
            migration,
            "std90",
            dpi,
        )

        _plot_combined_resolution_curves(
            series_by_bin,
            energy_bins,
            output_dir,
            f"residual_resolution_combined_{migration}.png",
            f"Residual resolution comparison for migration: {migration}",
            dpi,
        )

    for pdg_str in _all_pdgs_str:
        pdg_lbl = pid_label(int(pdg_str))

        if aggregated_results[pdg_str]:
            series_by_bin = {}
            for energy_bin, pairs in aggregated_results[pdg_str].items():
                _, residuals = _extract_true_reco_pairs(pairs)
                if residuals.size == 0:
                    continue
                series_by_bin[energy_bin] = residuals

            if series_by_bin:
                _plot_resolution_curve(
                    series_by_bin, energy_bins, output_dir,
                    f"residual_resolution_std_all_{pdg_str}.png",
                    f"Energy resolution (std) for all reconstructed {pdg_lbl}",
                    "σ((E_reco - E_true) / E_true)",
                    f"all_{pdg_str}", "std", dpi,
                )
                _plot_resolution_curve(
                    series_by_bin, energy_bins, output_dir,
                    f"residual_resolution_iqr84_16_all_{pdg_str}.png",
                    f"Energy resolution (IQR84-16/2) for all reconstructed {pdg_lbl}",
                    "IQR84-16 / 2 of ((E_reco - E_true) / E_true)",
                    f"all_{pdg_str}", "iqr", dpi,
                )
                _plot_resolution_curve(
                    series_by_bin, energy_bins, output_dir,
                    f"residual_resolution_std90_all_{pdg_str}.png",
                    f"Energy resolution (std90) for all reconstructed {pdg_lbl}",
                    "std90 of ((E_reco - E_true) / E_true)",
                    f"all_{pdg_str}", "std90", dpi,
                )
                _plot_combined_resolution_curves(
                    series_by_bin, energy_bins, output_dir,
                    f"residual_resolution_combined_all_{pdg_str}.png",
                    f"Residual resolution comparison for all reconstructed {pdg_lbl}",
                    dpi,
                )

        if aggregated_true_results[pdg_str]:
            true_series_by_bin = {}
            for energy_bin, pairs in aggregated_true_results[pdg_str].items():
                _, residuals = _extract_true_reco_pairs(pairs)
                if residuals.size == 0:
                    continue
                true_series_by_bin[energy_bin] = residuals

            if true_series_by_bin:
                _plot_resolution_curve(
                    true_series_by_bin, energy_bins, output_dir,
                    f"residual_resolution_std_all_true_{pdg_str}.png",
                    f"Energy resolution (std) for all true {pdg_lbl} (any reco)",
                    "σ((E_reco - E_true) / E_true)",
                    f"all_true_{pdg_str}", "std", dpi,
                )
                _plot_resolution_curve(
                    true_series_by_bin, energy_bins, output_dir,
                    f"residual_resolution_iqr84_16_all_true_{pdg_str}.png",
                    f"Energy resolution (IQR84-16/2) for all true {pdg_lbl} (any reco)",
                    "IQR84-16 / 2 of ((E_reco - E_true) / E_true)",
                    f"all_true_{pdg_str}", "iqr", dpi,
                )
                _plot_resolution_curve(
                    true_series_by_bin, energy_bins, output_dir,
                    f"residual_resolution_std90_all_true_{pdg_str}.png",
                    f"Energy resolution (std90) for all true {pdg_lbl} (any reco)",
                    "std90 of ((E_reco - E_true) / E_true)",
                    f"all_true_{pdg_str}", "std90", dpi,
                )
                _plot_combined_resolution_curves(
                    true_series_by_bin, energy_bins, output_dir,
                    f"residual_resolution_combined_all_true_{pdg_str}.png",
                    f"Residual resolution comparison for all true {pdg_lbl} (any reco)",
                    dpi,
                )


def plot_efficiency_vs_momentum(
    full_df,
    output_dir=".",
    dpi=150,
    n_bins=30,
    p_min=0.0,
    p_max=50.0,
    plot_type="default",
):
    """
    Plot n(gen→reco)/n_gen vs |p_gen| for every (Gen_pid, Reco_pid) pair found in
    full_df.  Produces:
      - One PNG per (Gen_pid, Reco_pid): efficiency_{gen_pid}_{reco_pid}.png
      - One global PNG per Gen_pid with all reco destinations overlaid plus a
        dashed total line: efficiency_global_{gen_pid}.png
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if plot_type == "default":
        bin_name = "p_bin"
        target_col = "Gen_P"
    elif plot_type == "theta":
        bin_name = "theta_bin"
        target_col = "Gen_theta"

    required_cols = {"Gen_pid", "Reco_pid", "Gen_Px", "Gen_Py", "Gen_Pz"}
    missing = required_cols - set(full_df.columns)
    if missing:
        print(f"[plot_efficiency_vs_momentum] Missing columns: {missing}. Skipping.")
        return

    if full_df.empty:
        print("[plot_efficiency_vs_momentum] Empty DataFrame. Skipping.")
        return

    # Forzar dtype numérico: los tipos cppyy de edm4hep pueden quedar como object en pandas
    for _col in ["Gen_Px", "Gen_Py", "Gen_Pz"]:
        full_df = full_df.copy()
        full_df[_col] = pd.to_numeric(full_df[_col], errors="coerce")

    df = full_df.loc[
        (full_df["Gen_pid"] != -999)
        & np.isfinite(full_df["Gen_Px"])
        & np.isfinite(full_df["Gen_Py"])
        & np.isfinite(full_df["Gen_Pz"])
    ].copy()

    if df.empty:
        print("[plot_efficiency_vs_momentum] No valid gen rows. Skipping.")
        return

    df["Gen_pid"] = df["Gen_pid"].abs()
    df["Reco_pid"] = df["Reco_pid"].abs()

    df["Gen_P"] = np.sqrt(df["Gen_Px"] ** 2 + df["Gen_Py"] ** 2 + df["Gen_Pz"] ** 2)
    if plot_type == "theta":
        df["Gen_theta"] = np.arccos(
            np.clip(df["Gen_Pz"] / df["Gen_P"], -1.0, 1.0)
        )
        p_min, p_max = 0.0, np.pi
    

    edges = np.linspace(p_min, p_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    df[bin_name] = pd.cut(
        df[target_col], bins=edges, labels=False, right=True, include_lowest=True
    )
    df = df.dropna(subset=[bin_name])
    df[bin_name] = df[bin_name].astype(int)

    n_gen_series = df.groupby(["Gen_pid", bin_name]).size().rename("n_gen")
    n_pair_series = (
        df.groupby(["Gen_pid", "Reco_pid", bin_name]).size().rename("n_pair")
    )

    counts_df = n_pair_series.reset_index().merge(
        n_gen_series.reset_index(), on=["Gen_pid", bin_name], how="left"
    )
    counts_df["eff"] = counts_df["n_pair"] / counts_df["n_gen"]
    counts_df["eff_err"] = np.sqrt(
        np.clip(
            counts_df["eff"] * (1.0 - counts_df["eff"]) / counts_df["n_gen"],
            0.0,
            None,
        )
    )

    def _eff_arrays(sub):
        eff_vals = np.full(n_bins, np.nan)
        eff_errs = np.full(n_bins, np.nan)
        valid = sub[[bin_name, "eff", "eff_err"]].copy()
        valid[bin_name] = valid[bin_name].astype(int)
        valid = valid[(valid[bin_name] >= 0) & (valid[bin_name] < n_bins)]
        eff_vals[valid[bin_name].values] = valid["eff"].values
        eff_errs[valid[bin_name].values] = np.where(
            np.isfinite(valid["eff_err"].values), valid["eff_err"].values, 0.0
        )
        return eff_vals, eff_errs

    def _style_ax(ax):
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis="y", which="both", left=True, right=True)
        ax.tick_params(axis="y", which="major", length=6, width=1.2)
        ax.tick_params(axis="y", which="minor", length=3, width=0.8)
        ax.grid(True)

    # ── Individual plots ──────────────────────────────────────────────────────
    for (gen_pid, reco_pid), sub in counts_df.groupby(["Gen_pid", "Reco_pid"]):
        eff_vals, eff_errs = _eff_arrays(sub)
        if not np.any(np.isfinite(eff_vals)):
            continue

        gen_lbl = pid_label(int(gen_pid))
        reco_lbl = pid_label(int(reco_pid)) if reco_pid != -999 else "unmatched"
        kind = "Efficiency" if gen_pid == reco_pid else "Migration"

        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
        ax.errorbar(
            centers,
            np.where(np.isfinite(eff_vals), eff_vals, np.nan),
            yerr=np.where(np.isfinite(eff_errs), eff_errs, 0.0),
            fmt="o-",
            capsize=4,
            linewidth=1.5,
            markersize=4,
        )
        ax.set_xlim(p_min, p_max)
        y_top = min(1.15, max(0.1, float(np.nanmax(eff_vals)) * 1.4))
        ax.set_ylim(0, y_top)
        if plot_type == "theta":
            ax.set_xlabel("θ_gen [rad]")
        else:
            ax.set_xlabel("|p_gen| [GeV]")
        ax.set_ylabel("n(gen→reco) / n_gen")
        ax.set_title(f"{kind}: {gen_lbl} → {reco_lbl}")
        _style_ax(ax)

        fname = os.path.join(output_dir, f"efficiency_{int(gen_pid)}_{int(reco_pid)}.png")
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved efficiency plot → {fname}")
        plt.close()

    # ── Global plots (one per Gen_pid) ────────────────────────────────────────
    prop_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for gen_pid, gen_group in counts_df.groupby("Gen_pid"):
        gen_lbl = pid_label(int(gen_pid))
        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

        total_eff = np.zeros(n_bins)
        color_idx = 0

        for reco_pid, sub in gen_group.groupby("Reco_pid"):
            eff_vals, eff_errs = _eff_arrays(sub)
            if not np.any(np.isfinite(eff_vals)):
                continue

            reco_lbl = pid_label(int(reco_pid)) if reco_pid != -999 else "unmatched"
            color = prop_cycle[color_idx % len(prop_cycle)]
            color_idx += 1

            ax.errorbar(
                centers,
                np.where(np.isfinite(eff_vals), eff_vals, np.nan),
                yerr=np.where(np.isfinite(eff_errs), eff_errs, 0.0),
                fmt="o-",
                capsize=3,
                linewidth=1.2,
                markersize=3,
                color=color,
                label=f"{gen_lbl} → {reco_lbl}",
            )
            total_eff += np.where(np.isfinite(eff_vals), eff_vals, 0.0)

        ax.plot(
            centers,
            total_eff,
            linestyle="--",
            linewidth=2.0,
            color="black",
            label="Total",
        )
        ax.set_xlim(p_min, p_max)
        ax.set_ylim(0, 1.15)
        if plot_type == "theta":
            ax.set_xlabel("θ_gen [rad]")
        else:
            ax.set_xlabel("|p_gen| [GeV]")
        ax.set_ylabel("n(gen→reco) / n_gen")
        ax.set_title(f"Efficiency & migrations: {gen_lbl}")
        ax.legend(fontsize=8, loc="best")
        _style_ax(ax)

        fname = os.path.join(output_dir, f"efficiency_global_{int(gen_pid)}.png")
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved global efficiency plot → {fname}")
        plt.close()


def plot_confusion_matrices(
    association_results_df,
    output_dir=".",
    figsize_per_cell=0.3,
    min_cell_size=4,
    cmap_abs="Blues",
    cmap_eff="Greens",
    cmap_pur="Oranges",
    annotate=True,
    save_individual=True,
    save_combined=True,
    dpi=150,
):
    """
    Parameters
    ----------
    association_results_df : dict
        Keys: "<|gen_pid|>_<|reco_pid|>_<energy_interval>"
        Values: integer counts
    output_dir : str
        Path to the folder where all output files will be saved.
        The folder is created automatically if it does not exist.
        Individual files are named  confusion_matrix_bin_<N>.pdf  (one per
        energy bin) and the combined file is  confusion_matrix_all_bins.pdf.
    figsize_per_cell : float
        Approximate inches per matrix cell (for auto-sizing).
    min_cell_size : int
        Minimum number of cells to use for figsize calculation.
    cmap_abs / cmap_eff / cmap_pur : str
        Matplotlib colormap names for the three matrix types.
    annotate : bool
        Write numeric values inside cells.
    save_individual : bool
        Save one PDF per energy bin containing the three matrices.
    save_combined : bool
        Save a single large PDF with all bins stacked.
    dpi : int
        Resolution for raster output.
    """
    os.makedirs(output_dir, exist_ok=True)

    parsed = parse_association_results(association_results_df)
    if not parsed:
        print("[plot_confusion_matrices] association_results_df is empty – nothing to plot.")
        return

    # ── Collect ALL pids across every bin ─────────────────────────────────────
    all_gen_pids  = sorted({gid for bd in parsed.values() for (gid, _) in bd})
    all_reco_pids = sorted({rid for bd in parsed.values() for (_, rid) in bd})

    energy_bins = sorted(parsed.keys())

    # ── Colour / style constants ───────────────────────────────────────────────
    TITLE_FONT  = dict(fontsize=11, fontweight="bold", color="#1a1a2e")
    LABEL_FONT  = dict(fontsize=9,  color="#2d2d2d")
    ANNOT_FONT_ABS = dict(fontsize=7.5, ha="center", va="center")
    SPINE_COLOR = "#cccccc"

    def _make_axes(n_rows, n_cols):
        """Compute figure width/height for a given matrix shape."""
        cells   = max(max(n_rows, n_cols), min_cell_size)
        side    = cells * figsize_per_cell
        w       = side + 2.5   # extra for labels
        h       = side + 1.5
        return w, h

    def _draw_matrix(ax, matrix, gen_ids, reco_ids, cmap, title,
                     fmt_abs=True, vmin=None, vmax=None):
        """Draw one confusion matrix on *ax*."""
        n_gen  = len(gen_ids)
        n_reco = len(reco_ids)

        im = ax.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=vmin if vmin is not None else 0,
            vmax=vmax if vmax is not None else (matrix.max() if matrix.max() > 0 else 1),
        )

        # Colourbar
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=7)
        if not fmt_abs:
            cb.set_label("fraction", fontsize=7, color="#555")

        # Axes ticks
        ax.set_xticks(range(n_reco))
        ax.set_yticks(range(n_gen))
        ax.set_xticklabels([pid_label(r) for r in reco_ids],
                           rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels([pid_label(g) for g in gen_ids], fontsize=8)

        ax.set_xlabel("Reco particle", **LABEL_FONT)
        ax.set_ylabel("Gen particle",  **LABEL_FONT)
        ax.set_title(title, **TITLE_FONT, pad=8)

        # Spines
        for spine in ax.spines.values():
            spine.set_edgecolor(SPINE_COLOR)

        # Grid lines between cells
        ax.set_xticks(np.arange(-0.5, n_reco, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n_gen,  1), minor=True)
        ax.grid(which="minor", color=SPINE_COLOR, linewidth=0.5)
        ax.tick_params(which="minor", bottom=False, left=False)

        # Annotations
        if annotate:
            thresh = matrix.max() / 2.0 if matrix.max() > 0 else 0.5
            for i in range(n_gen):
                for j in range(n_reco):
                    val = matrix[i, j]
                    if val == 0:
                        continue
                    color = "white" if val > thresh else "#333333"
                    txt = (f"{int(val):,}" if fmt_abs else f"{val:.2f}")
                    ax.text(j, i, txt, color=color, **ANNOT_FONT_ABS)

        return im

    # ── Build figures ──────────────────────────────────────────────────────────
    combined_figs = []

    # ── General matrix (all energy bins combined) ─────────────────────────────
    bd_general = defaultdict(int)
    for bd in parsed.values():
        for pair, count in bd.items():
            bd_general[pair] += count

    mat_abs_general = build_matrix(bd_general, all_gen_pids, all_reco_pids)
    mat_eff_general = safe_normalise_rows(mat_abs_general)
    mat_pur_general = safe_normalise_cols(mat_abs_general)

    n_g_general = len(all_gen_pids)
    n_r_general = len(all_reco_pids)
    w_general, h_general = _make_axes(n_g_general, n_r_general)

    general_matrices = [
        (mat_abs_general, cmap_abs, "Absolute counts",          True,  None, None, "absolute"),
        (mat_eff_general, cmap_eff, "Efficiency  (norm. by Gen)", False, 0,    1,   "efficiency"),
        (mat_pur_general, cmap_pur, "Purity  (norm. by Reco)",   False, 0,    1,   "purity"),
    ]

    general_figs = []
    for mat, cmap, title, fmt_abs, vmin, vmax, suffix in general_matrices:
        fig, ax = plt.subplots(figsize=(w_general + 2.5, h_general + 1.5), constrained_layout=True)
        fig.suptitle(
            f"All energy bins combined  —  {title}",
            fontsize=13, fontweight="bold", color="#0d0d1a",
        )
        _draw_matrix(ax, mat, all_gen_pids, all_reco_pids,
                    cmap, title, fmt_abs=fmt_abs, vmin=vmin, vmax=vmax)

        if save_individual:
            fname = os.path.join(output_dir, f"confusion_matrix_general_{suffix}.png")
            fig.savefig(fname, bbox_inches="tight", dpi=dpi)
            print(f"  Saved → {fname}")

        general_figs.append(fig)

    combined_figs.extend(general_figs)

    for bin_idx, ebin in enumerate(energy_bins):
            bd = parsed[ebin]

            # Only keep pids actually present in this bin
            gen_ids_bin  = sorted({gid for (gid, _) in bd})
            reco_ids_bin = sorted({rid for (_, rid) in bd})

            mat_abs = build_matrix(bd, gen_ids_bin, reco_ids_bin)
            mat_eff = safe_normalise_rows(mat_abs)   # norm by gen  → efficiency
            mat_pur = safe_normalise_cols(mat_abs)   # norm by reco → purity

            n_g = len(gen_ids_bin)
            n_r = len(reco_ids_bin)
            w, h = _make_axes(n_g, n_r)

            matrices = [
                (mat_abs, cmap_abs, "Absolute counts",          True,  None, None, "absolute"),
                (mat_eff, cmap_eff, "Efficiency  (norm. by Gen)", False, 0,    1,   "efficiency"),
                (mat_pur, cmap_pur, "Purity  (norm. by Reco)",   False, 0,    1,   "purity"),
            ]

            bin_figs = []
            for mat, cmap, title, fmt_abs, vmin, vmax, suffix in matrices:
                fig, ax = plt.subplots(figsize=(w + 2.5, h + 1.5), constrained_layout=True)
                fig.suptitle(
                    f"Energy bin:  {ebin} GeV  —  {title}",
                    fontsize=13, fontweight="bold", color="#0d0d1a",
                )
                _draw_matrix(ax, mat, gen_ids_bin, reco_ids_bin,
                            cmap, title, fmt_abs=fmt_abs, vmin=vmin, vmax=vmax)

                if save_individual:
                    fname = os.path.join(output_dir, f"confusion_matrix_bin_{bin_idx:02d}_{suffix}.png")
                    fig.savefig(fname, bbox_inches="tight", dpi=dpi)
                    print(f"  Saved → {fname}")

                bin_figs.append(fig)

            combined_figs.extend(bin_figs)

    # ── Combined figure (all bins stacked) ────────────────────────────────────
    if save_combined and combined_figs:
        from matplotlib.backends.backend_pdf import PdfPages
        combined_fname = os.path.join(output_dir, "confusion_matrix_all_bins.pdf")
        with PdfPages(combined_fname) as pdf:
            for fig in combined_figs:
                pdf.savefig(fig, bbox_inches="tight", dpi=dpi)
        print(f"\n  Combined PDF saved → {combined_fname}")

    plt.close("all")
    return combined_figs


# ─────────────────────────────────────────────
#  Quick self-test with synthetic data
# ─────────────────────────────────────────────

if __name__ == "__main__":

    import numpy as np

    rng = np.random.default_rng(42)

    # Particle IDs (absolute values, as stored in the dict)
    pids = [11, 13, 211, 22, 2212]

    bins_labels = [
        "(0, 1]", "(1, 5]", "(5, 10]", "(10, 20]",
        "(20, 30]", "(30, 45]", "(45, 100]", "(100, inf]",
    ]

    synthetic = {}
    for ebin in bins_labels:
        # Not every bin has every particle
        n_gen  = rng.integers(2, len(pids) + 1)
        n_reco = rng.integers(2, len(pids) + 1)
        gen_sub  = rng.choice(pids, size=n_gen,  replace=False).tolist()
        reco_sub = rng.choice(pids, size=n_reco, replace=False).tolist()

        for gid in gen_sub:
            for rid in reco_sub:
                count = int(rng.integers(0, 200))
                if count == 0:
                    continue
                key = f"{gid}_{rid}_{ebin}"
                synthetic[key] = count

    plot_confusion_matrices(
        synthetic,
        output_dir="confusion_matrices_output",
        save_individual=False,
        save_combined=True,
        annotate=True,
    )
