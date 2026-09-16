import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from matplotlib.lines import Line2D
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


# Color fijo por PDG (valor absoluto) para que una misma especie tenga siempre
# el mismo color en todos los plots globales: e→e es del mismo color en
# efficiency_global_11.png, en efficiency_global_22.png (γ→e) y en
# fake_rate_global.png.  Sin esto el color dependía del orden de aparición del
# groupby y cambiaba de fichero a fichero.
PID_COLORS = {
    11:   "#e6194b",   # e±      rojo
    13:   "#3cb44b",   # μ±      verde
    15:   "#911eb4",   # τ±      morado
    22:   "#4363d8",   # γ       azul
    111:  "#f58231",   # π⁰      naranja
    211:  "#0f766e",   # π±      verde azulado oscuro (nada de cian: se lee mal)
    130:  "#f032e6",   # K⁰L     magenta
    310:  "#a3c700",   # K⁰S     lima oscuro
    321:  "#c2185b",   # K±      rosa oscuro
    2112: "#b8860b",   # n       dorado oscuro (el gris chocaba con "unmatched")
    2212: "#9a6324",   # p       marrón
    3122: "#808000",   # Λ       oliva
    3112: "#000075",   # Σ⁻      azul marino
    3222: "#800000",   # Σ⁺      granate
    3312: "#1b7a4a",   # Ξ⁻      verde oscuro
    3322: "#d2691e",   # Ξ⁰      teja
    -999: "#808080",   # unmatched / fake   gris
    999:  "#808080",   # idem tras el abs() de Reco_pid
}

# Paleta de reserva para PDGs no listados; se asigna de forma determinista a
# partir del propio PDG, así que tampoco depende del orden de aparición.
# Todos son tonos oscuros/saturados y sin cian: los pasteles no se leen sobre
# fondo blanco cuando la curva es fina.
_PID_FALLBACK_COLORS = [
    "#6b6b6b", "#4d4d4d", "#8c564b", "#8c8c00", "#3b6ea5",
    "#b03a2e", "#6a3d9a", "#2e7d32", "#c1660a", "#1f4e79",
]


def pid_color(pid):
    """Color estable para un PDG (se espera el valor absoluto, salvo -999)."""
    pid = int(pid)
    if pid in PID_COLORS:
        return PID_COLORS[pid]
    return _PID_FALLBACK_COLORS[abs(pid) % len(_PID_FALLBACK_COLORS)]


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


def _pairs_slow_path(samples):
    """Recorrido fila a fila, para dicts o secuencias con valores ausentes."""
    rows = []
    for sample in samples:
        if isinstance(sample, dict):
            true_energy = sample.get("Gen_energy", sample.get("true_energy"))
            reco_energy = sample.get("Reco_energy", sample.get("reco_energy"))
        elif isinstance(sample, (tuple, list, np.ndarray)) and len(sample) >= 2:
            true_energy, reco_energy = sample[0], sample[1]
        else:
            continue

        if true_energy is None or reco_energy is None:
            continue
        rows.append((float(true_energy), float(reco_energy)))

    return np.asarray(rows, dtype=float).reshape(-1, 2)


def _as_pair_array(samples):
    """
    Normaliza cualquiera de los formatos aceptados a un array (N, 2) de
    [E_true, E_reco]: array (N, 2) ya construido (camino rápido, el que produce
    build_association_structures), lista de tuplas/listas, o lista de dicts.
    """
    if isinstance(samples, np.ndarray):
        return samples.astype(float, copy=False).reshape(-1, 2)
    if len(samples) == 0:
        return np.empty((0, 2), dtype=float)
    try:
        return np.asarray(samples, dtype=float).reshape(-1, 2)
    except (TypeError, ValueError):
        return _pairs_slow_path(samples)


def _extract_true_reco_pairs(samples):
    pairs = _as_pair_array(samples)
    if pairs.shape[0] == 0:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)

    true_energies = pairs[:, 0]
    reco_energies = pairs[:, 1]
    keep = (
        np.isfinite(true_energies)
        & np.isfinite(reco_energies)
        & (true_energies != 0)
    )
    true_energies = true_energies[keep]
    reco_energies = reco_energies[keep]

    return true_energies, (reco_energies - true_energies) / true_energies


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
ENERGY_RESOLUTION_METRICS = [
    ("std",   "std",          "σ((E_reco - E_true) / E_true)",                  "std"),
    ("iqr",   "iqr84_16",     "IQR84-16 / 2 of ((E_reco - E_true) / E_true)",   "IQR84-16/2"),
    ("std90", "std90",        "std90 of ((E_reco - E_true) / E_true)",          "std90"),
]


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


def _energy_resolution_series(energy_distribution_results, all_plot_pdgs=None):
    """
    Build the (series_by_bin) sets behind the resolution-vs-energy profiles.

    Returns (energy_bins, specs) where specs is a list of
    (stub, description, series_by_bin):

      - one entry per migration gen_pid -> reco_pid that is kept
        (same species, or anything involving photons),
      - one "all_<pdg>"      entry per requested PDG: everything reconstructed
        as that species,
      - one "all_true_<pdg>" entry per requested PDG: every true particle of
        that species, whatever it was reconstructed as.

    stub feeds the file name and description the plot title, so the single and
    the comparison plots stay in step.
    """
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
        pair_array = None
        for pdg_str in _all_pdgs_str:
            wants_reco = gen_pid_str != "999" and reco_pid == pdg_str
            wants_true = reco_pid != "999" and gen_pid_str == pdg_str
            if not (wants_reco or wants_true):
                continue
            # se normaliza una sola vez por clave, no una vez por PDG
            if pair_array is None:
                pair_array = _as_pair_array(energies)
            if wants_reco:
                aggregated_results[pdg_str][energy_bin].append(pair_array)
            if wants_true:
                aggregated_true_results[pdg_str][energy_bin].append(pair_array)

    specs = []

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

        if series_by_bin:
            specs.append((migration, f"migration: {migration}", series_by_bin))

    for pdg_str in _all_pdgs_str:
        pdg_lbl = pid_label(int(pdg_str))
        for aggregated, stub, description in (
            (aggregated_results[pdg_str], f"all_{pdg_str}",
             f"all reconstructed {pdg_lbl}"),
            (aggregated_true_results[pdg_str], f"all_true_{pdg_str}",
             f"all true {pdg_lbl} (any reco)"),
        ):
            if not aggregated:
                continue
            series_by_bin = {}
            for energy_bin, chunks in aggregated.items():
                _, residuals = _extract_true_reco_pairs(np.concatenate(chunks))
                if residuals.size == 0:
                    continue
                series_by_bin[energy_bin] = residuals
            if series_by_bin:
                specs.append((stub, description, series_by_bin))

    return energy_bins, specs


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

    energy_bins, specs = _energy_resolution_series(energy_distribution_results, all_plot_pdgs)

    for stub, description, series_by_bin in specs:
        for metric_mode, fname_tag, ylabel, title_tag in ENERGY_RESOLUTION_METRICS:
            _plot_resolution_curve(
                series_by_bin,
                energy_bins,
                output_dir,
                f"residual_resolution_{fname_tag}_{stub}.png",
                f"Energy resolution ({title_tag}) for {description}",
                ylabel,
                stub,
                metric_mode,
                dpi,
            )

        _plot_combined_resolution_curves(
            series_by_bin,
            energy_bins,
            output_dir,
            f"residual_resolution_combined_{stub}.png",
            f"Residual resolution comparison for {description}",
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
    legend_fontsize=13,
):
    """
    Plot n(gen→reco)/n_gen vs |p_gen| for every (Gen_pid, Reco_pid) pair found in
    full_df.  Produces:
      - One PNG per (Gen_pid, Reco_pid): efficiency_{gen_pid}_{reco_pid}.png
      - One global PNG per Gen_pid with all reco destinations overlaid plus a
        dashed total line: efficiency_global_{gen_pid}.png
    The global plots use a two-column legend with font size *legend_fontsize*.
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
        reco_lbl = pid_label(int(reco_pid)) if abs(int(reco_pid)) != 999 else "unmatched"
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
    # El color lo fija el Reco_pid (pid_color), no el orden del groupby: así la
    # curva "→ e±" es del mismo color en todos los ficheros.
    for gen_pid, gen_group in counts_df.groupby("Gen_pid"):
        gen_lbl = pid_label(int(gen_pid))
        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

        total_eff = np.zeros(n_bins)
        handles = []

        # La curva diagonal (gen == reco) va siempre la primera de la leyenda,
        # el resto de migraciones detrás por PDG creciente: así la eficiencia
        # propiamente dicha se lee arriba del todo en todos los ficheros.
        reco_groups = sorted(
            gen_group.groupby("Reco_pid"),
            key=lambda item: (0 if int(item[0]) == int(gen_pid) else 1, int(item[0])),
        )

        for reco_pid, sub in reco_groups:
            eff_vals, eff_errs = _eff_arrays(sub)
            if not np.any(np.isfinite(eff_vals)):
                continue

            reco_lbl = pid_label(int(reco_pid)) if abs(int(reco_pid)) != 999 else "unmatched"
            color = pid_color(reco_pid)

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
            handles.append(
                Line2D([], [], color=color, marker="o", markersize=4, linewidth=1.2,
                       label=f"{gen_lbl} → {reco_lbl}")
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
        handles.append(Line2D([], [], color="black", linestyle="--", linewidth=2.0,
                              label="Total"))
        ax.set_xlim(p_min, p_max)
        ax.set_ylim(0, 1.15)
        if plot_type == "theta":
            ax.set_xlabel("θ_gen [rad]")
        else:
            ax.set_xlabel("|p_gen| [GeV]")
        ax.set_ylabel("n(gen→reco) / n_gen")
        ax.set_title(f"Efficiency & migrations: {gen_lbl}")
        # Handles explícitos: matplotlib pondría el "Total" (Line2D) por delante
        # de los errorbars y la diagonal dejaría de ser la primera entrada.
        ax.legend(handles=handles, fontsize=legend_fontsize, ncol=2, loc="best",
                  columnspacing=1.2, handlelength=1.8)
        _style_ax(ax)

        fname = os.path.join(output_dir, f"efficiency_global_{int(gen_pid)}.png")
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved global efficiency plot → {fname}")
        plt.close()


def _fake_rate_counts(full_df, plot_type="default", n_bins=30, p_min=0.0,
                      p_max=50.0, ctx="plot_fake_rate_vs_momentum"):
    """
    Count PFOs and fakes per (Reco_pid, momentum/theta bin).

    Returns None if full_df cannot provide them, otherwise a dict with the
    counts_df (columns Reco_pid, <bin>, n_reco, n_fake, fake, fake_err), the
    bin centers and width, and the axis range/label: everything both the
    single-dataset plots and the comparison overlays need.
    """
    if plot_type == "theta":
        bin_name = "theta_bin"
        target_col = "Reco_theta"
    else:
        bin_name = "p_bin"
        target_col = "Reco_P"

    required_cols = {"Gen_pid", "Reco_pid", "Reco_Px", "Reco_Py", "Reco_Pz"}
    missing = required_cols - set(full_df.columns)
    if missing:
        print(f"[{ctx}] Missing columns: {missing}. Skipping.")
        return None

    if full_df.empty:
        print(f"[{ctx}] Empty DataFrame. Skipping.")
        return None

    df = full_df[full_df["Reco_pid"] != -999].copy()

    # Forzar dtype numérico: los tipos cppyy de edm4hep pueden quedar como object
    for _col in ["Reco_Px", "Reco_Py", "Reco_Pz"]:
        df[_col] = pd.to_numeric(df[_col], errors="coerce")

    df = df[
        np.isfinite(df["Reco_Px"])
        & np.isfinite(df["Reco_Py"])
        & np.isfinite(df["Reco_Pz"])
    ]

    if df.empty:
        print(f"[{ctx}] No valid reco rows. Skipping.")
        return None

    # ── Deduplicación: un PFO = una fila ──────────────────────────────────────
    # Con --dedup-mode gen (y con el matching por dR) un mismo PFO puede ser el
    # mejor candidato de varios gen y aparecer repetido; contar filas inflaría
    # el denominador.  Las filas CON match gen van primero para que un PFO
    # emparejado nunca acabe contado como fake.
    df["_matched"] = (df["Gen_pid"] != -999).astype(int)
    if {"event_id", "reco"}.issubset(df.columns):
        df = df.sort_values("_matched", ascending=False).drop_duplicates(
            subset=["event_id", "reco"], keep="first"
        )

    df["Reco_pid"] = df["Reco_pid"].abs()

    # |p_reco| a partir del momento, NO de Reco_energy: en la rama dR
    # Reco_energy es P(masa 0) para los matched pero la energía real del PFO
    # para los fakes, y mezclarlas sesgaría numerador contra denominador.
    df["Reco_P"] = np.sqrt(df["Reco_Px"] ** 2 + df["Reco_Py"] ** 2 + df["Reco_Pz"] ** 2)
    if plot_type == "theta":
        df["Reco_theta"] = np.arccos(
            np.clip(df["Reco_Pz"] / df["Reco_P"], -1.0, 1.0)
        )
        p_min, p_max = 0.0, np.pi

    edges = np.linspace(p_min, p_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]

    df[bin_name] = pd.cut(
        df[target_col], bins=edges, labels=False, right=True, include_lowest=True
    )
    df = df.dropna(subset=[bin_name])
    if df.empty:
        print(f"[{ctx}] No reco rows inside the binning range. Skipping.")
        return None
    df[bin_name] = df[bin_name].astype(int)

    n_reco_series = df.groupby(["Reco_pid", bin_name]).size().rename("n_reco")
    n_fake_series = (
        df[df["_matched"] == 0].groupby(["Reco_pid", bin_name]).size().rename("n_fake")
    )

    counts_df = n_reco_series.reset_index().merge(
        n_fake_series.reset_index(), on=["Reco_pid", bin_name], how="left"
    )
    counts_df["n_fake"] = counts_df["n_fake"].fillna(0.0)
    counts_df["fake"] = counts_df["n_fake"] / counts_df["n_reco"]
    # Error binomial sobre el denominador de PFOs
    counts_df["fake_err"] = np.sqrt(
        np.clip(counts_df["fake"] * (1.0 - counts_df["fake"]) / counts_df["n_reco"], 0.0, None)
    )


    return {
        "counts_df": counts_df,
        "centers": centers,
        "edges": edges,
        "bin_width": bin_width,
        "bin_name": bin_name,
        "n_bins": n_bins,
        "p_min": p_min,
        "p_max": p_max,
        "xlabel": "θ_reco [rad]" if plot_type == "theta" else "|p_reco| [GeV]",
    }


def plot_fake_rate_vs_momentum(
    full_df,
    output_dir=".",
    dpi=150,
    n_bins=30,
    p_min=0.0,
    p_max=50.0,
    plot_type="default",
    selection_note="",
    n_events=None,
):
    """
    Plot the fake rate n(reco sin match gen) / n_reco vs |p_reco| (o theta_reco)
    para cada Reco_pid presente en full_df.  Es el complemento reco-side de
    plot_efficiency_vs_momentum: aquí el denominador son PFOs, no MCParticles.

    IMPORTANTE: "fake" significa "PFO sin contrapartida gen DENTRO del conjunto
    gen seleccionado" (filtro de generatorStatus, max_gen_pdg, neutrinos, y en
    la rama dR también el cono dR<0.1 y el filtro de señal en detector).  No es
    un fake rate absoluto; usa `selection_note` para dejar constancia de la
    configuración con la que se produjo el plot.

    Produce:
      - fake_rate_{reco_pid}.png   : tasa de fakes por PDG reco
      - fake_rate_global.png       : todos los PDG superpuestos
      - fake_yield_{reco_pid}.png  : fakes por evento (absoluto), si n_events

    Parameters
    ----------
    full_df : pandas.DataFrame
        DataFrame de asociaciones con las columnas Gen_pid, Reco_pid, reco,
        event_id y Reco_Px/Py/Pz.
    n_events : int or None
        Número total de eventos procesados.  Si se pasa, se generan además los
        plots de rendimiento absoluto (fakes por evento).
    """
    os.makedirs(output_dir, exist_ok=True)

    counts = _fake_rate_counts(full_df, plot_type=plot_type, n_bins=n_bins,
                               p_min=p_min, p_max=p_max)
    if counts is None:
        return
    counts_df = counts["counts_df"]
    centers   = counts["centers"]
    bin_width = counts["bin_width"]
    bin_name  = counts["bin_name"]
    p_min, p_max = counts["p_min"], counts["p_max"]

    def _bin_arrays(sub, value_col, err_col):
        vals = np.full(n_bins, np.nan)
        errs = np.full(n_bins, np.nan)
        valid = sub[[bin_name, value_col, err_col]].copy()
        valid[bin_name] = valid[bin_name].astype(int)
        valid = valid[(valid[bin_name] >= 0) & (valid[bin_name] < n_bins)]
        vals[valid[bin_name].values] = valid[value_col].values
        errs[valid[bin_name].values] = np.where(
            np.isfinite(valid[err_col].values), valid[err_col].values, 0.0
        )
        return vals, errs

    def _style_ax(ax):
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis="y", which="both", left=True, right=True)
        ax.tick_params(axis="y", which="major", length=6, width=1.2)
        ax.tick_params(axis="y", which="minor", length=3, width=0.8)
        ax.grid(True)

    def _annotate_selection(fig):
        if selection_note:
            fig.suptitle(selection_note, fontsize=7, color="#555555", y=1.02)

    xlabel = counts["xlabel"]

    # ── Individual plots ──────────────────────────────────────────────────────
    for reco_pid, sub in counts_df.groupby("Reco_pid"):
        fake_vals, fake_errs = _bin_arrays(sub, "fake", "fake_err")
        if not np.any(np.isfinite(fake_vals)):
            continue

        reco_lbl = pid_label(int(reco_pid))

        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
        ax.errorbar(
            centers,
            fake_vals,
            yerr=np.where(np.isfinite(fake_errs), fake_errs, 0.0),
            fmt="o-",
            capsize=4,
            linewidth=1.5,
            markersize=4,
            color="#c0392b",
        )
        ax.set_xlim(p_min, p_max)
        y_top = min(1.15, max(0.1, float(np.nanmax(fake_vals)) * 1.4))
        ax.set_ylim(0, y_top)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("n(sin match gen) / n_reco")
        ax.set_title(f"Fake rate: {reco_lbl}")
        _style_ax(ax)
        _annotate_selection(fig)

        fname = os.path.join(output_dir, f"fake_rate_{int(reco_pid)}.png")
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved fake rate plot → {fname}")
        plt.close()

    # ── Global plot (todos los Reco_pid superpuestos) ─────────────────────────
    # Sin línea "Total": cada curva tiene su propio denominador.
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    drawn = False

    for reco_pid, sub in counts_df.groupby("Reco_pid"):
        fake_vals, fake_errs = _bin_arrays(sub, "fake", "fake_err")
        if not np.any(np.isfinite(fake_vals)):
            continue
        color = pid_color(reco_pid)
        drawn = True
        ax.errorbar(
            centers,
            fake_vals,
            yerr=np.where(np.isfinite(fake_errs), fake_errs, 0.0),
            fmt="o-",
            capsize=3,
            linewidth=1.2,
            markersize=3,
            color=color,
            label=pid_label(int(reco_pid)),
        )

    if drawn:
        ax.set_xlim(p_min, p_max)
        ax.set_ylim(0, 1.15)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("n(sin match gen) / n_reco")
        ax.set_title("Fake rate por tipo de PFO")
        ax.legend(fontsize=8, loc="best")
        _style_ax(ax)
        _annotate_selection(fig)
        fname = os.path.join(output_dir, "fake_rate_global.png")
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved global fake rate plot → {fname}")
    plt.close()

    # ── Rendimiento absoluto: fakes por evento ────────────────────────────────
    # El ratio puede ser engañoso a bajo p, donde el denominador también cae;
    # para el tau reco lo que importa es cuántos PFOs espurios hay por evento.
    if not n_events:
        return

    counts_df["yield"] = counts_df["n_fake"] / (float(n_events) * bin_width)
    counts_df["yield_err"] = np.sqrt(counts_df["n_fake"]) / (float(n_events) * bin_width)

    x_unit = "rad" if plot_type == "theta" else "GeV"
    for reco_pid, sub in counts_df.groupby("Reco_pid"):
        yield_vals, yield_errs = _bin_arrays(sub, "yield", "yield_err")
        if not np.any(np.isfinite(yield_vals)) or np.nanmax(yield_vals) <= 0:
            continue

        reco_lbl = pid_label(int(reco_pid))
        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
        ax.errorbar(
            centers,
            yield_vals,
            yerr=np.where(np.isfinite(yield_errs), yield_errs, 0.0),
            fmt="o-",
            capsize=4,
            linewidth=1.5,
            markersize=4,
            color="#8e44ad",
        )
        ax.set_xlim(p_min, p_max)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"fakes / evento / {x_unit}")
        ax.set_title(f"Fake yield: {reco_lbl}  ({n_events} eventos)")
        _style_ax(ax)
        _annotate_selection(fig)

        fname = os.path.join(output_dir, f"fake_yield_{int(reco_pid)}.png")
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        print(f"Saved fake yield plot → {fname}")
        plt.close()


def _fit_text_to_cells(ax, texts, n_cols, n_rows, fill_w=0.7, fill_h=0.5):
    """
    Agranda por igual las anotaciones *texts* de un imshow hasta el mayor tamaño
    de fuente con el que la más ancha sigue cabiendo en su celda (fill_w/fill_h
    son la fracción de anchura/altura de celda que puede ocupar el texto).
    """
    if not texts:
        return
    fig = ax.figure
    # Dibujar antes de medir: constrained_layout fija aquí el tamaño final del eje
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    box = ax.get_window_extent(renderer)
    cell_w = box.width / n_cols
    cell_h = box.height / n_rows
    scale = min(
        min(fill_w * cell_w / bb.width, fill_h * cell_h / bb.height)
        for bb in (t.get_window_extent(renderer) for t in texts)
    )
    size = texts[0].get_fontsize() * scale
    for t in texts:
        t.set_fontsize(size)


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
    LABEL_FONT  = dict(fontsize=14, color="#2d2d2d")
    TICK_FONTSIZE = 13
    CBAR_FONTSIZE = 12
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
        cb.ax.tick_params(labelsize=CBAR_FONTSIZE)
        if not fmt_abs:
            cb.set_label("fraction", fontsize=CBAR_FONTSIZE, color="#555")

        # Axes ticks
        ax.set_xticks(range(n_reco))
        ax.set_yticks(range(n_gen))
        ax.set_xticklabels([pid_label(r) for r in reco_ids],
                           rotation=45, ha="right", fontsize=TICK_FONTSIZE)
        ax.set_yticklabels([pid_label(g) for g in gen_ids], fontsize=TICK_FONTSIZE)

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
            texts = []
            for i in range(n_gen):
                for j in range(n_reco):
                    val = matrix[i, j]
                    if val == 0:
                        continue
                    color = "white" if val > thresh else "#333333"
                    txt = (f"{int(val):,}" if fmt_abs else f"{val:.2f}")
                    texts.append(ax.text(j, i, txt, color=color, **ANNOT_FONT_ABS))
            _fit_text_to_cells(ax, texts, n_cols=n_reco, n_rows=n_gen)

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


# ─────────────────────────────────────────────
#  Resolution  (momentum and polar angle)
# ─────────────────────────────────────────────

# Cada observable se describe una sola vez: la columna de residuo que añade
# _residual_frame, cómo rotular el eje y el título, con cuántos decimales citar
# las estadísticas y qué límites usar en el rango automático.  Todo lo demás
# (histogramas, comparación, ficheros de salida) es común a los dos.
RESOLUTION_SPECS = {
    "p": dict(
        column="p_res",
        title_kind="Momentum resolution",
        xlabel="(|p$_{reco}$| − |p$_{gen}$|) / |p$_{gen}$|",
        xlabel_ascii="(|p$_{reco}$| - |p$_{gen}$|) / |p$_{gen}$|",
        stat_fmt=".4f",
        range_kwargs=dict(min_half_width=0.02, hard_lo=-1.0, hard_hi=3.0),
    ),
    "theta": dict(
        column="theta_res",
        title_kind="Theta resolution",
        xlabel=r"$\theta_{reco}$ − $\theta_{gen}$ [mrad]",
        xlabel_ascii=r"$\theta_{reco}$ - $\theta_{gen}$ [mrad]",
        stat_fmt=".3f",
        # El residuo en θ es absoluto y no tiene tope físico como el relativo en
        # |p| (que no puede bajar de −1), así que no se recorta a mano.
        range_kwargs=dict(min_half_width=1.0, hard_lo=-np.inf, hard_hi=np.inf),
    ),
}


def _resolution_spec(observable):
    try:
        return RESOLUTION_SPECS[observable]
    except KeyError:
        raise ValueError(
            f"Unknown resolution observable: {observable!r}. "
            f"Expected one of {sorted(RESOLUTION_SPECS)}."
        ) from None


def _auto_residual_range(residuals, n_sigma=5.0, min_half_width=0.02,
                         hard_lo=-1.0, hard_hi=3.0):
    """Pick an x-range that shows the bulk of the distribution.

    The width comes from a robust core estimate — median ± n_sigma × σ_rob with
    σ_rob = (p84 − p16)/2 — instead of the plain min/max or a high percentile,
    so a per-mille tail of catastrophic mismeasurements cannot squash the peak
    into a single bin.  The result is clipped to the data extent and to
    physically sensible limits (a residual cannot go below −1).
    """
    med = float(np.median(residuals))
    p16, p84 = np.percentile(residuals, [16, 84])
    sigma_rob = float(p84 - p16) / 2.0

    if not np.isfinite(sigma_rob) or sigma_rob <= 0:
        sigma_rob = float(np.std(residuals))
    if not np.isfinite(sigma_rob) or sigma_rob <= 0:
        sigma_rob = min_half_width

    lo = med - n_sigma * sigma_rob
    hi = med + n_sigma * sigma_rob

    # No dejar marco vacío a los lados si los datos no llegan hasta ahí
    lo = max(lo, float(np.min(residuals)))
    hi = min(hi, float(np.max(residuals)))

    # Mantener el 0 dentro del marco: es la referencia visual del plot
    lo = min(lo, -min_half_width)
    hi = max(hi, min_half_width)

    return max(lo, hard_lo), min(hi, hard_hi)


def _plot_residual_hist(residuals, title, fname, dpi, n_bins, r_range=None,
                        spec=None, ctx="plot_momentum_resolution"):
    """Draw one residual histogram with its statistics box.

    spec selects the observable (see RESOLUTION_SPECS): it fixes the x-label,
    the number of decimals of the statistics box and the automatic range.
    """
    spec = spec or RESOLUTION_SPECS["p"]
    fmt = spec["stat_fmt"]
    residuals = np.asarray(residuals, dtype=float)
    residuals = residuals[np.isfinite(residuals)]
    if residuals.size < 2:
        print(f"[{ctx}] Too few entries for {fname}. Skipping.")
        return

    if r_range is None:
        lo, hi = _auto_residual_range(residuals, **spec["range_kwargs"])
    else:
        lo, hi = r_range

    inside = residuals[(residuals >= lo) & (residuals <= hi)]
    if inside.size < 2:
        print(f"[{ctx}] No entries inside ({lo}, {hi}) for {fname}. Skipping.")
        return

    counts, edges = np.histogram(inside, bins=n_bins, range=(lo, hi))
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    ax.hist(inside, bins=edges, histtype="stepfilled", color="#4363d8",
            alpha=0.35, edgecolor="#4363d8", linewidth=1.4)
    ax.errorbar(centers, counts, yerr=np.sqrt(counts), fmt="none",
                ecolor="#4363d8", elinewidth=1.0, capsize=2)

    # Las estadísticas se citan sobre TODAS las entradas, no solo las visibles;
    # el número de las que se van fuera del marco se indica aparte.
    n_out = residuals.size - inside.size
    stats = [
        f"entries = {residuals.size}",
        f"mean = {np.mean(residuals):+{fmt}}",
        f"std = {np.std(residuals):{fmt}}",
    ]
    std90 = _std90(residuals)
    if std90 is not None:
        stats.append(f"std90 = {std90:{fmt}}")
    stats.append(f"median = {np.median(residuals):+{fmt}}")
    if n_out:
        stats.append(f"out of range = {n_out}")

    ax.axvline(0.0, color="black", linestyle=":", linewidth=1.2)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(spec["xlabel"])
    ax.set_ylabel("Entries")
    ax.set_title(title)
    ax.text(0.02, 0.98, "\n".join(stats), transform=ax.transAxes,
            va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="both", which="both", top=True, right=True)
    ax.grid(True, alpha=0.3)

    plt.savefig(fname, dpi=dpi, bbox_inches="tight")
    print(f"Saved {spec['title_kind'].lower()} plot → {fname}")
    plt.close()


def _residual_frame(full_df, ctx="plot_momentum_resolution"):
    """
    Common preprocessing of the resolution plots: returns (df, df_res) with the
    absolute PDGs, |p| and theta of both partners and the residuals p_res and
    theta_res, or (None, None) if full_df cannot provide it.

    df      : every row, with Gen_pid/Reco_pid in absolute value and Gen_P/Reco_P
    df_res  : only the gen-reco matched rows, with the extra columns
              Gen_theta/Reco_theta (rad), p_res (relative) and theta_res (mrad)
    """
    required_cols = {"Gen_pid", "Reco_pid",
                     "Gen_Px", "Gen_Py", "Gen_Pz",
                     "Reco_Px", "Reco_Py", "Reco_Pz"}
    missing = required_cols - set(full_df.columns)
    if missing:
        print(f"[{ctx}] Missing columns: {missing}. Skipping.")
        return None, None

    if full_df.empty:
        print(f"[{ctx}] Empty DataFrame. Skipping.")
        return None, None

    df = full_df.copy()
    # Forzar dtype numérico: los tipos cppyy de edm4hep pueden quedar como object
    for _col in ["Gen_Px", "Gen_Py", "Gen_Pz", "Reco_Px", "Reco_Py", "Reco_Pz"]:
        df[_col] = pd.to_numeric(df[_col], errors="coerce")

    df["Gen_pid"] = df["Gen_pid"].abs()
    df["Reco_pid"] = df["Reco_pid"].abs()

    df["Gen_P"] = np.sqrt(df["Gen_Px"] ** 2 + df["Gen_Py"] ** 2 + df["Gen_Pz"] ** 2)
    df["Reco_P"] = np.sqrt(df["Reco_Px"] ** 2 + df["Reco_Py"] ** 2 + df["Reco_Pz"] ** 2)

    valid = (
        np.isfinite(df["Gen_P"])
        & np.isfinite(df["Reco_P"])
        & (df["Gen_P"] > 0)
        & (df["Gen_pid"] != 999)
        & (df["Reco_pid"] != 999)
    )
    df_res = df.loc[valid].copy()
    if df_res.empty:
        print(f"[{ctx}] No gen-reco matched rows. Skipping.")
        return df, None

    df_res["p_res"] = (df_res["Reco_P"] - df_res["Gen_P"]) / df_res["Gen_P"]

    # θ = arccos(pz/|p|).  El residuo va en absoluto (mrad) y no en relativo:
    # θ pasa por cero y no tiene una escala propia como |p|, así que (θ_reco −
    # θ_gen)/θ_gen explotaría en el forward.  Gen_P > 0 lo garantiza la máscara
    # de arriba; Reco_P == 0 dejaría θ_reco indefinido, y esas filas salen NaN
    # y las descartan los propios plots.
    gen_cos = np.clip(df_res["Gen_Pz"].to_numpy() / df_res["Gen_P"].to_numpy(), -1.0, 1.0)
    reco_p = df_res["Reco_P"].to_numpy()
    reco_cos = np.full(reco_p.shape, np.nan)
    ok = np.isfinite(reco_p) & (reco_p > 0)
    reco_cos[ok] = np.clip(df_res["Reco_Pz"].to_numpy()[ok] / reco_p[ok], -1.0, 1.0)

    df_res["Gen_theta"] = np.arccos(gen_cos)
    df_res["Reco_theta"] = np.arccos(reco_cos)
    df_res["theta_res"] = (df_res["Reco_theta"] - df_res["Gen_theta"]) * 1e3
    return df, df_res


def _plot_resolution_by_species(
    full_df,
    observable,
    output_dir=".",
    dpi=150,
    n_bins=100,
    r_min=None,
    r_max=None,
    min_entries=20,
    title_suffix="",
    ctx="plot_momentum_resolution",
):
    """Shared body of plot_momentum_resolution / plot_theta_resolution.

    observable picks the residual (see RESOLUTION_SPECS); everything else --
    the species loop, the file names and the "matched"/"pred" split -- is the
    same for both.
    """
    spec = _resolution_spec(observable)
    col = spec["column"]
    kind = spec["title_kind"]

    os.makedirs(output_dir, exist_ok=True)

    # r_min/r_max fijan el rango a mano; si no se dan, cada especie escoge el
    # suyo a partir de sus percentiles (ver _auto_residual_range).
    r_range = None if (r_min is None or r_max is None) else (r_min, r_max)

    df, df_res = _residual_frame(full_df, ctx=ctx)
    if df_res is None:
        return

    # ── a) Partículas correctamente identificadas ────────────────────────────
    correct = df_res.loc[df_res["Gen_pid"] == df_res["Reco_pid"]]
    for pid, sub in correct.groupby("Gen_pid"):
        if len(sub) < min_entries:
            continue
        lbl = pid_label(int(pid))
        _plot_residual_hist(
            sub[col].to_numpy(),
            title=f"{kind}, correctly identified: {lbl} → {lbl}{title_suffix}",
            fname=os.path.join(output_dir, f"resolution_matched_{int(pid)}.png"),
            dpi=dpi, n_bins=n_bins, r_range=r_range, spec=spec, ctx=ctx,
        )

    # ── b) Todo lo reconstruido como esa especie (acierte o no) ──────────────
    # Denominador informativo: PFOs de esa especie sin gen asociado (fakes), que
    # no pueden entrar en el histograma por no tener |p_gen| de referencia.
    n_no_gen = (
        df.loc[(df["Gen_pid"] == 999) | ~np.isfinite(df["Gen_P"]) | (df["Gen_P"] <= 0)]
        .groupby("Reco_pid").size()
    )
    for pid, sub in df_res.groupby("Reco_pid"):
        if len(sub) < min_entries:
            continue
        lbl = pid_label(int(pid))
        n_fake = int(n_no_gen.get(pid, 0))
        purity = (sub["Gen_pid"] == pid).mean()
        _plot_residual_hist(
            sub[col].to_numpy(),
            title=(f"{kind}, reconstructed as {lbl} "
                   f"(purity {purity:.1%}, {n_fake} without gen match){title_suffix}"),
            fname=os.path.join(output_dir, f"resolution_pred_{int(pid)}.png"),
            dpi=dpi, n_bins=n_bins, r_range=r_range, spec=spec, ctx=ctx,
        )


def plot_momentum_resolution(
    full_df,
    output_dir=".",
    dpi=150,
    n_bins=100,
    r_min=None,
    r_max=None,
    min_entries=20,
    title_suffix="",
):
    """
    Plot the relative momentum resolution (|p_reco| − |p_gen|) / |p_gen| for
    every particle species found in full_df.  Two PNGs per species:

      - resolution_matched_{pid}.png : only correctly identified particles,
        i.e. rows with |Gen_pid| == |Reco_pid| == pid.
      - resolution_pred_{pid}.png    : every PFO reconstructed as that species
        (|Reco_pid| == pid) that has a gen match, regardless of whether the
        PID is right.

    Rows without a gen or reco partner (unmatched / fakes) carry no residual
    and are dropped; the number of such reco-as-pid rows is reported in the
    "pred" plot title as they cannot enter the distribution.

    title_suffix is appended to every plot title, e.g. to state an extra
    selection applied upstream (" | E > 10 GeV").
    """
    _plot_resolution_by_species(
        full_df, "p", output_dir=output_dir, dpi=dpi, n_bins=n_bins,
        r_min=r_min, r_max=r_max, min_entries=min_entries,
        title_suffix=title_suffix, ctx="plot_momentum_resolution",
    )


def plot_theta_resolution(
    full_df,
    output_dir=".",
    dpi=150,
    n_bins=100,
    r_min=None,
    r_max=None,
    min_entries=20,
    title_suffix="",
):
    """
    Same as plot_momentum_resolution, but for the polar angle: the absolute
    residual θ_reco − θ_gen in mrad, with θ = arccos(p_z / |p|).  The residual
    is absolute rather than relative because θ crosses zero and has no scale of
    its own, so a relative one would diverge in the forward region.

    Written to its own directory (one PNG per species, same file names), so it
    sits next to the momentum plots instead of overwriting them.
    """
    _plot_resolution_by_species(
        full_df, "theta", output_dir=output_dir, dpi=dpi, n_bins=n_bins,
        r_min=r_min, r_max=r_max, min_entries=min_entries,
        title_suffix=title_suffix, ctx="plot_theta_resolution",
    )

COMPARISON_COLORS = ["#4363d8", "#e6194b", "#3cb44b", "#f58231",
                     "#911eb4", "#8c564b", "#0f766e", "#808000"]


def _comparison_residual_range(residual_sets, r_min=None, r_max=None,
                               range_kwargs=None):
    """Common x-range for an overlay: the union of the per-dataset auto ranges."""
    if r_min is not None and r_max is not None:
        return r_min, r_max
    range_kwargs = range_kwargs or {}
    los, his = [], []
    for res in residual_sets:
        lo, hi = _auto_residual_range(res, **range_kwargs)
        los.append(lo)
        his.append(hi)
    return min(los), max(his)


def _plot_residual_comparison(residuals_by_label, title, fname, dpi, n_bins,
                              r_range=None, normalize=True, extra_by_label=None,
                              legend_fontsize=11, spec=None,
                              ctx="compare_momentum_resolution"):
    """Overlay several residual distributions of the same observable.

    spec selects the observable (see RESOLUTION_SPECS): x-label, statistics
    format and automatic range.
    """
    spec = spec or RESOLUTION_SPECS["p"]
    fmt = spec["stat_fmt"]
    clean = {}
    for label, res in residuals_by_label.items():
        res = np.asarray(res, dtype=float)
        res = res[np.isfinite(res)]
        if res.size >= 2:
            clean[label] = res
    if not clean:
        print(f"[{ctx}] Too few entries for {fname}. Skipping.")
        return

    lo, hi = _comparison_residual_range(clean.values(), *(r_range or (None, None)),
                                        range_kwargs=spec["range_kwargs"])
    edges = np.linspace(lo, hi, n_bins + 1)

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    legend_entries = []
    for i, (label, res) in enumerate(clean.items()):
        color = COMPARISON_COLORS[i % len(COMPARISON_COLORS)]
        inside = res[(res >= lo) & (res <= hi)]
        if inside.size < 2:
            print(f"[{ctx}] {label}: no entries inside "
                  f"({lo:.3f}, {hi:.3f}). Skipping this dataset.")
            continue
        # density=True para que dos muestras con estadística muy distinta sigan
        # siendo comparables de forma; las estadísticas de la leyenda se citan
        # sobre TODAS las entradas, no solo las visibles.
        ax.hist(inside, bins=edges, histtype="step", color=color, linewidth=1.8,
                density=normalize)
        std90 = _std90(res)
        # Una línea por bloque de estadísticas: en una sola línea la leyenda se
        # sale del marco y hay que encogerla hasta hacerla ilegible.
        lines = [str(label)]
        if extra_by_label and label in extra_by_label:
            lines.append(extra_by_label[label])
        lines.append(f"N = {res.size}, mean = {np.mean(res):+{fmt}}")
        stat_line = f"std = {np.std(res):{fmt}}"
        if std90 is not None:
            stat_line += f", std90 = {std90:{fmt}}"
        lines.append(stat_line)
        legend_entries.append((color, "\n".join(lines)))

    if not legend_entries:
        plt.close()
        return

    handles = [plt.Line2D([], [], color=c, linewidth=2.2, label=t)
               for c, t in legend_entries]
    ax.legend(handles=handles, loc="upper left", fontsize=legend_fontsize,
              framealpha=0.9, labelspacing=0.9, handlelength=1.6,
              borderpad=0.6, fancybox=True)
    # Sitio para el recuadro: si no, tapa el pico de la distribución.
    ax.set_ylim(top=ax.get_ylim()[1] * (1.25 + 0.12 * len(legend_entries)))

    ax.axvline(0.0, color="black", linestyle=":", linewidth=1.2)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(spec["xlabel_ascii"])
    ax.set_ylabel("Normalized entries" if normalize else "Entries")
    ax.set_title(title)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="both", which="both", top=True, right=True)
    ax.grid(True, alpha=0.3)

    plt.savefig(fname, dpi=dpi, bbox_inches="tight")
    print(f"Saved {spec['title_kind'].lower()} comparison -> {fname}")
    plt.close()


def _compare_resolution_by_species(
    dfs_by_label,
    observable,
    output_dir=".",
    dpi=150,
    n_bins=100,
    r_min=None,
    r_max=None,
    min_entries=20,
    title_suffix="",
    normalize=True,
    pids=None,
    legend_fontsize=11,
    ctx="compare_momentum_resolution",
):
    """Shared body of compare_momentum_resolution / compare_theta_resolution."""
    spec = _resolution_spec(observable)
    col = spec["column"]
    kind = spec["title_kind"]

    os.makedirs(output_dir, exist_ok=True)
    r_range = None if (r_min is None or r_max is None) else (r_min, r_max)

    matched = {}   # {pid: {label: residuals}}
    pred = {}
    pred_info = {}  # {pid: {label: (purity, n_fake)}}

    for label, full_df in dfs_by_label.items():
        df, df_res = _residual_frame(full_df, ctx=f"{ctx}[{label}]")
        if df_res is None:
            continue

        correct = df_res.loc[df_res["Gen_pid"] == df_res["Reco_pid"]]
        for pid, sub in correct.groupby("Gen_pid"):
            if len(sub) < min_entries:
                continue
            matched.setdefault(int(pid), {})[label] = sub[col].to_numpy()

        n_no_gen = (
            df.loc[(df["Gen_pid"] == 999) | ~np.isfinite(df["Gen_P"]) | (df["Gen_P"] <= 0)]
            .groupby("Reco_pid").size()
        )
        for pid, sub in df_res.groupby("Reco_pid"):
            if len(sub) < min_entries:
                continue
            pid_i = int(pid)
            pred.setdefault(pid_i, {})[label] = sub[col].to_numpy()
            pred_info.setdefault(pid_i, {})[label] = (
                float((sub["Gen_pid"] == pid).mean()), int(n_no_gen.get(pid, 0))
            )

    if pids is not None:
        keep = {int(p) for p in pids}
        matched = {k: v for k, v in matched.items() if k in keep}
        pred = {k: v for k, v in pred.items() if k in keep}

    for pid, per_label in sorted(matched.items()):
        lbl = pid_label(pid)
        _plot_residual_comparison(
            per_label,
            title=f"{kind}, correctly identified: {lbl} -> {lbl}{title_suffix}",
            fname=os.path.join(output_dir, f"resolution_matched_{pid}.png"),
            dpi=dpi, n_bins=n_bins, r_range=r_range, normalize=normalize,
            legend_fontsize=legend_fontsize, spec=spec, ctx=ctx,
        )

    for pid, per_label in sorted(pred.items()):
        lbl = pid_label(pid)
        # La pureza y los fakes son por dataset, así que van en la leyenda de
        # cada curva en lugar del título.
        extra = {
            label: (f"purity {pred_info[pid][label][0]:.1%}, "
                    f"{pred_info[pid][label][1]} without gen match")
            for label in per_label
        }
        _plot_residual_comparison(
            per_label, extra_by_label=extra,
            title=f"{kind}, reconstructed as {lbl}{title_suffix}",
            fname=os.path.join(output_dir, f"resolution_pred_{pid}.png"),
            dpi=dpi, n_bins=n_bins, r_range=r_range, normalize=normalize,
            legend_fontsize=legend_fontsize, spec=spec, ctx=ctx,
        )


def compare_momentum_resolution(
    dfs_by_label,
    output_dir=".",
    dpi=150,
    n_bins=100,
    r_min=None,
    r_max=None,
    min_entries=20,
    title_suffix="",
    normalize=True,
    pids=None,
    legend_fontsize=11,
):
    """
    Overlay the relative momentum resolution of several datasets (e.g. two
    detectors) in the same plots, one PNG per species:

      - resolution_matched_{pid}.png : rows with |Gen_pid| == |Reco_pid| == pid.
      - resolution_pred_{pid}.png    : rows with |Reco_pid| == pid and a gen
        match, regardless of whether the PID is right.

    dfs_by_label
        Ordered mapping {label: full_df}; the label names the curve in the
        legend, which also carries entries / mean / std / std90 per dataset.
    normalize
        Draw densities instead of raw counts, so samples of different size stay
        comparable in shape.  Set to False to compare absolute yields.
    pids
        Restrict the output to these |PDG| codes; by default every species that
        reaches min_entries in at least one dataset is plotted.

    Species missing from a dataset (or below min_entries there) simply do not
    contribute a curve, so an asymmetric comparison is still drawn.
    """
    _compare_resolution_by_species(
        dfs_by_label, "p", output_dir=output_dir, dpi=dpi, n_bins=n_bins,
        r_min=r_min, r_max=r_max, min_entries=min_entries,
        title_suffix=title_suffix, normalize=normalize, pids=pids,
        legend_fontsize=legend_fontsize, ctx="compare_momentum_resolution",
    )


def compare_theta_resolution(
    dfs_by_label,
    output_dir=".",
    dpi=150,
    n_bins=100,
    r_min=None,
    r_max=None,
    min_entries=20,
    title_suffix="",
    normalize=True,
    pids=None,
    legend_fontsize=11,
):
    """
    Same as compare_momentum_resolution, but overlaying the absolute polar-angle
    residual θ_reco − θ_gen (mrad) of several datasets.  Same file names, in a
    directory of its own.
    """
    _compare_resolution_by_species(
        dfs_by_label, "theta", output_dir=output_dir, dpi=dpi, n_bins=n_bins,
        r_min=r_min, r_max=r_max, min_entries=min_entries,
        title_suffix=title_suffix, normalize=normalize, pids=pids,
        legend_fontsize=legend_fontsize, ctx="compare_theta_resolution",
    )


def _plot_resolution_curve_comparison(series_by_label, energy_bins, output_dir,
                                      filename, title, ylabel, metric_mode, dpi,
                                      log_x=True, legend_fontsize=11):
    """Overlay the resolution-vs-energy profile of several datasets."""
    fig, ax = plt.subplots(figsize=(8, 6))
    plotted_any = False
    all_centers = []

    for i, (label, series_by_bin) in enumerate(series_by_label.items()):
        data_to_plot = []
        for energy_bin in energy_bins:
            residuals = series_by_bin.get(energy_bin, [])
            value = _resolution_value(residuals, metric_mode)
            if value is None:
                continue
            data_to_plot.append((energy_bin, value,
                                 _resolution_error(residuals, metric_mode),
                                 np.size(residuals)))
        if not data_to_plot:
            print(f"[compare_energy_resolution] {label}: no bins for {filename}. "
                  "Skipping this dataset.")
            continue

        bins, values, errors, counts = zip(*data_to_plot)
        centers = [_energy_bin_center(b) for b in bins]
        all_centers.extend(centers)
        ax.errorbar(
            centers, values,
            yerr=[e if e is not None else 0.0 for e in errors],
            fmt="o-", capsize=4, linewidth=1.5, markersize=5,
            color=COMPARISON_COLORS[i % len(COMPARISON_COLORS)],
            label=f"{label}\n{len(bins)} bins, N = {int(sum(counts))}",
        )
        plotted_any = True

    if not plotted_any:
        plt.close()
        return False

    if log_x:
        ax.set_xscale("log")
    else:
        ticks = sorted(set(all_centers))
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(round(t, 1)) for t in ticks], rotation=45, ha="right")
    ax.set_xlabel("True energy bin (GeV)", loc="left")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=legend_fontsize, framealpha=0.9, labelspacing=0.9,
              handlelength=1.6, borderpad=0.6, fancybox=True)
    ax.grid(True, which="both" if log_x else "major")
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="y", which="both", left=True, right=True)
    ax.tick_params(axis="y", which="major", length=6, width=1.2)
    ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    plt.tight_layout()
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=dpi)
    print(f"Saved energy resolution comparison -> {output_path}")
    plt.close()
    return True


def compare_energy_resolution(
    edists_by_label,
    output_dir=".",
    dpi=150,
    all_plot_pdgs=None,
    title_suffix="",
    log_x=True,
    legend_fontsize=11,
):
    """
    Overlay the resolution-vs-true-energy profiles of several datasets (e.g. two
    detectors), one PNG per series and metric:

        residual_resolution_{std,iqr84_16,std90}_{stub}.png

    with the same stubs as plot_energy_distributions (a migration such as
    "22_22", or the aggregates "all_<pdg>" / "all_true_<pdg>"), so the
    comparison directory mirrors the single-dataset one.

    edists_by_label
        Ordered mapping {label: energy_distribution_results}, i.e. the second
        output of build_association_structures for each dataset.

    A series present in only one dataset is still drawn, with the single curve.
    The combined std/IQR/std90 plot has no comparison counterpart: three metrics
    times N datasets in one frame is unreadable.
    """
    os.makedirs(output_dir, exist_ok=True)

    # {stub: (description, {label: series_by_bin})}, en el orden en que aparecen
    series_by_stub = {}
    energy_bins = []
    for label, edist in edists_by_label.items():
        bins, specs = _energy_resolution_series(edist, all_plot_pdgs)
        for b in bins:
            if b not in energy_bins:
                energy_bins.append(b)
        for stub, description, series_by_bin in specs:
            entry = series_by_stub.setdefault(stub, (description, {}))
            entry[1][label] = series_by_bin

    energy_bins.sort(key=lambda x: float(x.strip("()[]").split(",")[0]))

    for stub, (description, per_label) in series_by_stub.items():
        for metric_mode, fname_tag, ylabel, title_tag in ENERGY_RESOLUTION_METRICS:
            _plot_resolution_curve_comparison(
                per_label, energy_bins, output_dir,
                f"residual_resolution_{fname_tag}_{stub}.png",
                f"Energy resolution ({title_tag}) for {description}{title_suffix}",
                ylabel, metric_mode, dpi, log_x=log_x,
                legend_fontsize=legend_fontsize,
            )


def compare_fake_rate_vs_momentum(
    dfs_by_label,
    output_dir=".",
    dpi=150,
    n_bins=30,
    p_min=0.0,
    p_max=50.0,
    plot_type="default",
    n_events_by_label=None,
    pids=None,
    legend_fontsize=11,
    title_suffix="",
):
    """
    Overlay the fake rate (and the absolute fake yield) of several datasets,
    one PNG per reconstructed species:

      - fake_rate_{reco_pid}.png  : n(reco without gen match) / n_reco
      - fake_yield_{reco_pid}.png : fakes per event and per GeV (or rad), only
        for the datasets whose event count is given in n_events_by_label

    Same definition and binning as plot_fake_rate_vs_momentum, so the curves are
    directly comparable with the single-dataset plots; "fake" still means "PFO
    without a gen partner INSIDE the selected gen set", which makes the
    comparison meaningful only between runs made with the same gen selection.

    dfs_by_label
        Ordered mapping {label: full_df}.
    n_events_by_label
        Mapping {label: n_events}.  The yield plot needs every dataset to have
        one; datasets without it are dropped from that plot (not from the rate).
    """
    os.makedirs(output_dir, exist_ok=True)
    n_events_by_label = n_events_by_label or {}

    per_label = {}
    for label, full_df in dfs_by_label.items():
        counts = _fake_rate_counts(full_df, plot_type=plot_type, n_bins=n_bins,
                                   p_min=p_min, p_max=p_max,
                                   ctx=f"compare_fake_rate_vs_momentum[{label}]")
        if counts is None:
            continue
        df = counts["counts_df"]
        n_events = n_events_by_label.get(label)
        if n_events:
            df["yield"] = df["n_fake"] / (float(n_events) * counts["bin_width"])
            df["yield_err"] = np.sqrt(df["n_fake"]) / (float(n_events) * counts["bin_width"])
        per_label[label] = counts

    if not per_label:
        print("[compare_fake_rate_vs_momentum] No dataset with valid reco rows. Skipping.")
        return

    ref = next(iter(per_label.values()))
    centers, bin_name = ref["centers"], ref["bin_name"]
    x_lo, x_hi, xlabel = ref["p_min"], ref["p_max"], ref["xlabel"]
    x_unit = "rad" if plot_type == "theta" else "GeV"

    reco_pids = sorted({int(pid)
                        for c in per_label.values()
                        for pid in c["counts_df"]["Reco_pid"].unique()})
    if pids is not None:
        keep = {int(p) for p in pids}
        reco_pids = [pid for pid in reco_pids if pid in keep]

    def _series(counts, reco_pid, value_col, err_col):
        """(values, errors) padded to the binning, or None if not available."""
        df = counts["counts_df"]
        if value_col not in df.columns:
            return None
        sub = df[df["Reco_pid"] == reco_pid]
        if sub.empty:
            return None
        vals = np.full(counts["n_bins"], np.nan)
        errs = np.full(counts["n_bins"], np.nan)
        idx = sub[bin_name].astype(int).values
        inside = (idx >= 0) & (idx < counts["n_bins"])
        vals[idx[inside]] = sub[value_col].values[inside]
        errs[idx[inside]] = np.where(np.isfinite(sub[err_col].values[inside]),
                                     sub[err_col].values[inside], 0.0)
        if not np.any(np.isfinite(vals)):
            return None
        return vals, errs

    def _overlay(value_col, err_col, reco_pid, ylabel, title, fname,
                 y_top=None, extra_by_label=None):
        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
        drawn, maxima = False, []
        for i, (label, counts) in enumerate(per_label.items()):
            series = _series(counts, reco_pid, value_col, err_col)
            if series is None:
                continue
            vals, errs = series
            legend = str(label)
            if extra_by_label and label in extra_by_label:
                legend += f"\n{extra_by_label[label]}"
            n_fake = int(counts["counts_df"].loc[
                counts["counts_df"]["Reco_pid"] == reco_pid, "n_fake"].sum())
            n_reco = int(counts["counts_df"].loc[
                counts["counts_df"]["Reco_pid"] == reco_pid, "n_reco"].sum())
            legend += f"\nfakes = {n_fake} / {n_reco} PFOs"
            ax.errorbar(
                counts["centers"], vals,
                yerr=np.where(np.isfinite(errs), errs, 0.0),
                fmt="o-", capsize=4, linewidth=1.5, markersize=4,
                color=COMPARISON_COLORS[i % len(COMPARISON_COLORS)],
                label=legend,
            )
            maxima.append(float(np.nanmax(vals)))
            drawn = True

        if not drawn:
            plt.close()
            return

        ax.set_xlim(x_lo, x_hi)
        if y_top is None:
            ax.set_ylim(bottom=0)
        else:
            ax.set_ylim(0, min(y_top, max(0.1, max(maxima) * 1.6)))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=legend_fontsize, framealpha=0.9, labelspacing=0.9,
                  handlelength=1.6, borderpad=0.6, fancybox=True, loc="best")
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis="y", which="both", left=True, right=True)
        ax.tick_params(axis="y", which="major", length=6, width=1.2)
        ax.tick_params(axis="y", which="minor", length=3, width=0.8)
        ax.grid(True)
        path = os.path.join(output_dir, fname)
        plt.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"Saved fake rate comparison -> {path}")
        plt.close()

    for reco_pid in reco_pids:
        reco_lbl = pid_label(reco_pid)
        _overlay("fake", "fake_err", reco_pid,
                 "n(sin match gen) / n_reco",
                 f"Fake rate: {reco_lbl}{title_suffix}",
                 f"fake_rate_{reco_pid}.png", y_top=1.15)
        # El yield es absoluto, así que el número de eventos de cada dataset es
        # parte de la definición de la curva y va en la leyenda.
        _overlay("yield", "yield_err", reco_pid,
                 f"fakes / evento / {x_unit}",
                 f"Fake yield: {reco_lbl}{title_suffix}",
                 f"fake_yield_{reco_pid}.png",
                 extra_by_label={label: f"{n} eventos"
                                 for label, n in n_events_by_label.items() if n})


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
