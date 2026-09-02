#!/usr/bin/env python3
"""Extend frozen migration matrices with a stored-tau-ancestor truth subset.

This reads frozen G/Ldirect/Lancestor outcomes and MCParticle parent relations.
It does not run or alter association, linking, reconstruction, simulation, or
generation.  Existing all-truth artifacts are used as immutable regressions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pyarrow.parquet as pq
import uproot

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts.analysis import produce_migration_matrices as migration  # noqa: E402
from scripts.analysis.produce_migration_matrices import (  # noqa: E402
    DEFINITIONS,
    OUTCOME_LABELS,
    SAMPLES,
    TRUTH_BY_PDG,
    TRUTH_LABELS,
    TRUTH_ORDER,
    annotate,
    category_from_pdg,
    fiducial_outcome,
    iter_g_groups,
    l_outcome,
    load_frozen_l,
    manifest_index,
    outcomes,
    read_csv,
    style_axis,
    write_csv,
)

FINAL = TAU_ROOT = None
G_PATHS = {}
MANIFESTS = {}


def configure(path: Path) -> None:
    global FINAL, TAU_ROOT, G_PATHS, MANIFESTS
    migration.configure(path)
    FINAL = migration.FINAL
    TAU_ROOT = FINAL / "tau_origin"
    G_PATHS = migration.G_PATHS
    MANIFESTS = migration.MANIFESTS

POPULATIONS = ("inclusive_tau_origin", "fiducial_tau_origin")
ROOT_BRANCHES = (
    "MCParticles/MCParticles.PDG",
    "MCParticles/MCParticles.parents_begin",
    "MCParticles/MCParticles.parents_end",
    "_MCParticles_parents/_MCParticles_parents.index",
)


def reco_selection(population: str) -> str:
    if population.startswith("inclusive_"):
        return "inclusive"
    if population.startswith("fiducial_"):
        return "fiducial"
    raise ValueError(population)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compact_truth_path(sample: str, meta: dict[str, str]) -> Path:
    if sample == "W":
        return Path(meta["truth_out"])
    ancestor = Path(meta["ancestor_output"])
    return ancestor.parents[1] / "selected_truth" / f"{sample}_suffix{meta['suffix']}_selected_truth.parquet"


def load_compact_tau_flags(sample: str, meta: dict[str, str]):
    path = compact_truth_path(sample, meta)
    if sample == "W":
        columns = ["event_in_file", "mc_index", "truth_pdg", "explicit_tau_lineage"]
        index_col, pdg_col = "mc_index", "truth_pdg"
    else:
        columns = ["event_in_file", "truth_index", "pdg", "explicit_tau_lineage"]
        index_col, pdg_col = "truth_index", "pdg"
    rows = pq.read_table(path, columns=columns).to_pylist()
    flags = {}
    pdgs = {}
    for row in rows:
        pdg = int(row[pdg_col])
        if abs(pdg) not in TRUTH_BY_PDG:
            continue
        key = (int(row["event_in_file"]), int(row[index_col]))
        if key in flags:
            raise AssertionError(f"duplicate compact truth key: {sample} {path} {key}")
        if row["explicit_tau_lineage"] is None:
            raise AssertionError(f"null explicit_tau_lineage: {sample} {path} {key}")
        flags[key] = bool(row["explicit_tau_lineage"])
        pdgs[key] = pdg
    return path, flags, pdgs


def nearest_tau_for_particle(pdgs, parents, start: int):
    """Return nearest stored tau with cycle protection and deterministic ties."""
    queue = deque((parent, 1, frozenset({start})) for parent in parents[start])
    best_depth = {}
    taus = []
    while queue:
        index, depth, path = queue.popleft()
        if index in path:
            raise ValueError(f"cycle in MC ancestry reachable from index {start}")
        if index < 0 or index >= len(pdgs):
            raise IndexError(f"invalid parent index {index} for collection of size {len(pdgs)}")
        if index in best_depth and best_depth[index] <= depth:
            continue
        best_depth[index] = depth
        if abs(int(pdgs[index])) == 15:
            taus.append((depth, index, int(pdgs[index])))
        next_path = path | {index}
        queue.extend((parent, depth + 1, next_path) for parent in parents[index])
    return min(taus) if taus else None


def derive_tau_records(source_rec: Path, truth_rows: list[tuple]):
    selected = {}
    by_event = defaultdict(list)
    for event, gen, _reco, gen_pdg, _reco_pdg in truth_rows:
        if abs(int(gen_pdg)) not in TRUTH_BY_PDG:
            continue
        key = (int(event), int(gen))
        if key in selected:
            raise AssertionError(f"duplicate selected truth key in G rows: {source_rec} {key}")
        selected[key] = int(gen_pdg)
        by_event[int(event)].append(int(gen))

    with uproot.open(source_rec) as root_file:
        tree = root_file["events"]
        arrays = tree.arrays(list(ROOT_BRANCHES), library="ak")
    pdg_events = arrays[ROOT_BRANCHES[0]]
    begin_events = arrays[ROOT_BRANCHES[1]]
    end_events = arrays[ROOT_BRANCHES[2]]
    relation_events = arrays[ROOT_BRANCHES[3]]
    if by_event and max(by_event) >= len(pdg_events):
        raise AssertionError(f"G event outside REC range: {source_rec}")

    result = {}
    for event, indices in by_event.items():
        pdgs = pdg_events[event].to_list()
        begins = begin_events[event].to_list()
        ends = end_events[event].to_list()
        relations = relation_events[event].to_list()
        if not (len(pdgs) == len(begins) == len(ends)):
            raise AssertionError(f"MCParticle branch-length mismatch: {source_rec} event {event}")
        parents = []
        for begin, end in zip(begins, ends):
            if not (0 <= int(begin) <= int(end) <= len(relations)):
                raise AssertionError(f"invalid MC parent relation range: {source_rec} event {event}")
            parents.append([int(value) for value in relations[int(begin):int(end)]])
        for index in indices:
            if index < 0 or index >= len(pdgs):
                raise AssertionError(f"selected MC index out of range: {source_rec} {(event, index)}")
            key = (event, index)
            if int(pdgs[index]) != selected[key]:
                raise AssertionError(f"G/REC truth PDG mismatch: {source_rec} {key}")
            nearest = nearest_tau_for_particle(pdgs, parents, index)
            result[key] = {
                "has_tau_ancestor": nearest is not None,
                "nearest_tau_index": None if nearest is None else nearest[1],
                "nearest_tau_pdg": None if nearest is None else nearest[2],
                "tau_ancestry_depth": None if nearest is None else nearest[0],
            }
    if len(result) != len(selected):
        raise AssertionError(f"ancestry coverage mismatch: {source_rec}")
    return result, selected


def process_file(task):
    sample, source_name, truth_rows, meta = task
    source_rec = Path(meta["source_rec"] if sample == "W" else meta["input_REC"])
    tau_records, selected_pdgs = derive_tau_records(source_rec, truth_rows)
    compact_path, compact_flags, compact_pdgs = load_compact_tau_flags(sample, meta)
    selected_keys = set(selected_pdgs)
    expected_compact = {
        key for key, pdg in selected_pdgs.items()
        if sample != "W" or abs(pdg) in {22, 211}
    }
    if set(compact_flags) != expected_compact:
        raise AssertionError(
            f"compact/G selected-truth key mismatch: {sample} {source_name} "
            f"missing={len(expected_compact-set(compact_flags))} extra={len(set(compact_flags)-expected_compact)}"
        )
    compact_mismatches = 0
    for key, flag in compact_flags.items():
        if compact_pdgs[key] != selected_pdgs[key]:
            raise AssertionError(f"compact/G truth PDG mismatch: {sample} {source_name} {key}")
        compact_mismatches += int(flag != tau_records[key]["has_tau_ancestor"])
    if compact_mismatches:
        raise AssertionError(f"stored/derived tau-lineage mismatch: {sample} {source_name} {compact_mismatches}")

    pfo_map, direct_by_truth, ancestor_by_truth, _observed = load_frozen_l(sample, meta)
    all_counts = Counter()
    tau_counts = Counter()
    all_denominators = Counter()
    tau_denominators = Counter()
    depth_counts = Counter()
    for event, gen, reco, gen_pdg, reco_pdg in truth_rows:
        truth_pdg = abs(int(gen_pdg))
        if truth_pdg not in TRUTH_BY_PDG:
            continue
        truth = TRUTH_BY_PDG[truth_pdg]
        truth_key = (int(event), int(gen))
        ancestry = tau_records[truth_key]
        all_denominators[truth] += 1
        if ancestry["has_tau_ancestor"]:
            tau_denominators[truth] += 1
            depth_counts[(truth, ancestry["tau_ancestry_depth"], ancestry["nearest_tau_pdg"])] += 1

        if abs(int(reco_pdg)) == 999:
            g_inc, g_pfo = "association_unmatched", None
        else:
            pfo_key = (int(event), int(reco))
            g_pfo = pfo_map.get(pfo_key)
            if g_pfo is None:
                raise AssertionError(f"G representative missing from frozen PFOs: {sample} {source_name} {pfo_key}")
            g_inc = category_from_pdg(reco_pdg)
            if abs(int(g_pfo["pfo_type"])) != abs(int(reco_pdg)):
                raise AssertionError(f"G reco/PFO type mismatch: {sample} {source_name} {pfo_key}")
        g_fid = fiducial_outcome(g_inc, g_pfo)
        direct_inc, direct_pfo = l_outcome(direct_by_truth.get(truth_key, []), pfo_map, int(event))
        direct_fid = fiducial_outcome(direct_inc, direct_pfo)
        ancestor_inc, ancestor_pfo = l_outcome(ancestor_by_truth.get(truth_key, []), pfo_map, int(event))
        ancestor_fid = fiducial_outcome(ancestor_inc, ancestor_pfo)
        for definition, inclusive, fiducial in (
            ("G", g_inc, g_fid),
            ("Ldirect", direct_inc, direct_fid),
            ("Lancestor", ancestor_inc, ancestor_fid),
        ):
            all_counts[(definition, "inclusive", truth, inclusive)] += 1
            all_counts[(definition, "fiducial", truth, fiducial)] += 1
            if ancestry["has_tau_ancestor"]:
                tau_counts[(definition, "inclusive", truth, inclusive)] += 1
                tau_counts[(definition, "fiducial", truth, fiducial)] += 1

    return {
        "sample": sample,
        "source_name": source_name,
        "source_rec": str(source_rec),
        "compact_truth": str(compact_path),
        "selected_keys": len(selected_keys),
        "compact_keys": len(compact_flags),
        "tau_keys": sum(tau_denominators.values()),
        "compact_mismatches": compact_mismatches,
        "all_counts": dict(all_counts),
        "tau_counts": dict(tau_counts),
        "all_denominators": dict(all_denominators),
        "tau_denominators": dict(tau_denominators),
        "depth_counts": dict(depth_counts),
    }


def task_stream():
    manifest_maps = {sample: manifest_index(sample) for sample in SAMPLES}
    for sample in SAMPLES:
        expected = {"W": 996, "P8C": 100, "P8O": 18}[sample]
        count = 0
        for source_name, truth_rows in iter_g_groups(G_PATHS[sample]):
            meta = manifest_maps[sample].get(Path(source_name).name)
            if meta is None:
                raise AssertionError(f"{sample} G source absent from frozen manifest: {source_name}")
            count += 1
            yield sample, source_name, truth_rows, meta
        if count != expected:
            raise AssertionError(f"{sample} G file groups {count} != {expected}")


def merge_counter(target: Counter, sample: str, values: dict):
    for key, value in values.items():
        target[(sample, *key)] += int(value)


def run_parallel(n_workers: int):
    all_counts = Counter()
    tau_counts = Counter()
    all_denominators = Counter()
    tau_denominators = Counter()
    depth_counts = Counter()
    validation_rows = []
    completed = Counter()
    tasks = iter(task_stream())
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        pending = set()
        for _ in range(2 * n_workers):
            try:
                pending.add(pool.submit(process_file, next(tasks)))
            except StopIteration:
                break
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                result = future.result()
                sample = result["sample"]
                merge_counter(all_counts, sample, result["all_counts"])
                merge_counter(tau_counts, sample, result["tau_counts"])
                for truth, value in result["all_denominators"].items():
                    all_denominators[(sample, truth)] += int(value)
                for truth, value in result["tau_denominators"].items():
                    tau_denominators[(sample, truth)] += int(value)
                for key, value in result["depth_counts"].items():
                    depth_counts[(sample, *key)] += int(value)
                validation_rows.append({
                    "sample": sample,
                    "source_name": result["source_name"],
                    "source_rec": result["source_rec"],
                    "compact_truth": result["compact_truth"],
                    "selected_truth_keys": result["selected_keys"],
                    "compact_truth_keys": result["compact_keys"],
                    "tau_origin_keys": result["tau_keys"],
                    "stored_vs_recursive_mismatches": result["compact_mismatches"],
                    "invalid_parent_references": 0,
                    "ancestry_cycles": 0,
                    "pass": True,
                })
                completed[sample] += 1
                sample_expected = {"W": 996, "P8C": 100, "P8O": 18}[sample]
                if completed[sample] % 25 == 0 or completed[sample] == sample_expected:
                    print(f"ANCESTRY {sample}: {completed[sample]}/{sample_expected}", flush=True)
                try:
                    pending.add(pool.submit(process_file, next(tasks)))
                except StopIteration:
                    pass
    expected = Counter({"W": 996, "P8C": 100, "P8O": 18})
    if completed != expected:
        raise AssertionError(f"processed files mismatch: {completed} != {expected}")
    validation_rows.sort(key=lambda row: (SAMPLES.index(row["sample"]), row["source_name"]))
    return all_counts, tau_counts, all_denominators, tau_denominators, depth_counts, validation_rows


def tau_matrix(tau_counts, tau_denominators, sample, definition, population):
    selection = reco_selection(population)
    order = outcomes(selection)
    raw = np.zeros((len(TRUTH_ORDER), len(order)), dtype=np.int64)
    values = np.zeros_like(raw, dtype=float)
    for i, truth in enumerate(TRUTH_ORDER):
        denominator = tau_denominators[(sample, truth)]
        if denominator <= 0:
            raise AssertionError(f"zero tau-origin denominator: {sample} {truth}")
        for j, outcome in enumerate(order):
            count = tau_counts[(sample, definition, selection, truth, outcome)]
            raw[i, j] = count
            values[i, j] = 100.0 * count / denominator
    return raw, values


def plot_absolute_tau(path_stem: Path, values, sample: str, definition: str, population: str):
    selection = reco_selection(population)
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    image = ax.imshow(values, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    annotate(ax, values)
    style_axis(ax, selection)
    label = definition.replace("Ldirect", r"$L_{direct}$").replace("Lancestor", r"$L_{ancestor}$")
    ax.set_title(f"{sample} — {label} — {population}\ntruth-normalized migration fraction [%]", fontsize=13)
    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label("migration fraction [%]")
    fig.text(0.5, 0.012, "Tau-origin selected truth subset; each truth row sums to 100%.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_contact_tau(path_stem: Path, matrices: dict, population: str):
    selection = reco_selection(population)
    fig, axes = plt.subplots(3, 3, figsize=(25, 13.5), sharex=True, sharey=True)
    image = None
    for i, sample in enumerate(SAMPLES):
        for j, definition in enumerate(DEFINITIONS):
            ax = axes[i, j]
            values = matrices[(sample, definition, population)]
            image = ax.imshow(values, cmap="Blues", vmin=0, vmax=100, aspect="auto")
            annotate(ax, values)
            style_axis(ax, selection)
            if i != 2:
                ax.set_xlabel("")
            if j != 0:
                ax.set_ylabel("")
            label = definition.replace("Ldirect", r"$L_{direct}$").replace("Lancestor", r"$L_{ancestor}$")
            ax.set_title(f"{sample} — {label}", fontsize=14, fontweight="bold")
    fig.suptitle(f"W / P8 truth-normalized migration matrices — {population}", fontsize=20, fontweight="bold")
    color_axis = fig.add_axes([0.925, 0.14, 0.012, 0.70])
    colorbar = fig.colorbar(image, cax=color_axis)
    colorbar.set_label("migration fraction [%]")
    fig.text(0.5, 0.008, "Rows: tau-origin selected truth particle. Columns: reconstructed outcome. Fixed scale 0–100%; each row sums to 100%.", ha="center", fontsize=11)
    fig.subplots_adjust(left=0.055, right=0.90, bottom=0.12, top=0.91, wspace=0.12, hspace=0.28)
    fig.savefig(path_stem.with_suffix(".png"), dpi=240, bbox_inches="tight")
    fig.savefig(path_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_difference_tau(path_stem: Path, values, title: str, population: str, limit: float):
    selection = reco_selection(population)
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    image = ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto")
    annotate(ax, values, difference=True, limit=limit)
    style_axis(ax, selection)
    ax.set_title(f"{title} — {population}\ntruth-normalized migration difference [percentage points]", fontsize=13)
    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label("difference [percentage points]")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-workers", type=int, default=8)
    args = parser.parse_args()
    if args.n_workers < 1:
        raise ValueError("--n-workers must be positive")
    if not FINAL.is_dir():
        raise FileNotFoundError(FINAL)
    if TAU_ROOT.exists():
        raise FileExistsError(f"refusing to overwrite {TAU_ROOT}")
    existing_table = FINAL / "tables/migration_matrices_all.csv"
    existing_summary = FINAL / "summary.json"
    if not existing_table.is_file() or not existing_summary.is_file():
        raise FileNotFoundError("validated all-truth artifacts are incomplete")

    print(f"PROCESS frozen outcomes and recursive stored ancestry with {args.n_workers} workers", flush=True)
    all_counts, tau_counts, all_denominators, tau_denominators, depth_counts, ancestry_rows = run_parallel(args.n_workers)

    regression_rows = []
    for row in read_csv(existing_table):
        sample = row["sample"]
        definition = row["truth_definition"]
        selection = row["selection"]
        truth = row["truth_category"]
        outcome = row["reco_category"]
        observed_count = all_counts[(sample, definition, selection, truth, outcome)]
        observed_denominator = all_denominators[(sample, truth)]
        passed = observed_count == int(row["count"]) and observed_denominator == int(row["truth_denominator"])
        regression_rows.append({
            "sample": sample,
            "truth_definition": definition,
            "selection": selection,
            "truth_category": truth,
            "reco_category": outcome,
            "expected_count": row["count"],
            "observed_count": observed_count,
            "expected_denominator": row["truth_denominator"],
            "observed_denominator": observed_denominator,
            "pass": passed,
        })
    if not all(row["pass"] for row in regression_rows):
        raise AssertionError("existing all-truth numerical regression failed")

    staging = Path(tempfile.mkdtemp(prefix=".tau_origin.partial.", dir=FINAL))
    try:
        for relative in ("inclusive", "fiducial", "contact_sheets", "difference_matrices", "tables", "validation", "reports"):
            (staging / relative).mkdir()

        denominator_rows = []
        for sample in SAMPLES:
            for truth in TRUTH_ORDER:
                n_all = all_denominators[(sample, truth)]
                n_tau = tau_denominators[(sample, truth)]
                denominator_rows.append({
                    "sample": sample,
                    "truth_category": truth,
                    "N_all": n_all,
                    "N_tau_origin": n_tau,
                    "fraction_tau_origin": f"{n_tau / n_all:.17g}",
                })
        write_csv(staging / "tables/tau_origin_denominators.csv", denominator_rows)

        matrices = {}
        matrix_rows = []
        matrix_validation = []
        for population in POPULATIONS:
            selection = reco_selection(population)
            for sample in SAMPLES:
                for definition in DEFINITIONS:
                    raw, values = tau_matrix(tau_counts, tau_denominators, sample, definition, population)
                    matrices[(sample, definition, population)] = values
                    for i, truth in enumerate(TRUTH_ORDER):
                        denominator = tau_denominators[(sample, truth)]
                        count_sum = int(raw[i].sum())
                        percent_sum = float(values[i].sum())
                        passed = count_sum == denominator and abs(percent_sum - 100.0) < 1e-10
                        matrix_validation.append({
                            "sample": sample,
                            "truth_definition": definition,
                            "population": population,
                            "truth_category": truth,
                            "denominator": denominator,
                            "count_sum": count_sum,
                            "percent_sum": f"{percent_sum:.15g}",
                            "difference_from_100": f"{percent_sum - 100.0:.15g}",
                            "pass": passed,
                        })
                        for j, outcome in enumerate(outcomes(selection)):
                            matrix_rows.append({
                                "sample": sample,
                                "truth_definition": definition,
                                "population": population,
                                "truth_category": truth,
                                "reco_category": outcome,
                                "count": int(raw[i, j]),
                                "truth_denominator": denominator,
                                "fraction": f"{values[i, j] / 100.0:.17g}",
                                "percent": f"{values[i, j]:.15g}",
                            })
        if not all(row["pass"] for row in matrix_validation):
            raise AssertionError("tau-origin matrix row-sum validation failed")

        combined_rows = []
        for row in read_csv(existing_table):
            combined_rows.append({
                "sample": row["sample"],
                "truth_definition": row["truth_definition"],
                "population": f"{row['selection']}_all",
                "truth_category": row["truth_category"],
                "reco_category": row["reco_category"],
                "count": row["count"],
                "truth_denominator": row["truth_denominator"],
                "fraction": row["fraction"],
                "percent": row["percent"],
            })
        combined_rows.extend(matrix_rows)
        write_csv(staging / "tables/migration_matrices_with_tau_origin.csv", combined_rows)
        write_csv(staging / "tables/photon_migration_rows_with_tau_origin.csv", [row for row in combined_rows if row["truth_category"] == "photon"])
        write_csv(staging / "tables/charged_pion_migration_rows_with_tau_origin.csv", [row for row in combined_rows if row["truth_category"] == "charged_pion"])
        write_csv(staging / "validation/all_truth_numerical_regression.csv", regression_rows)
        write_csv(staging / "validation/tau_origin_matrix_validation.csv", matrix_validation)
        write_csv(staging / "validation/ancestry_validation_by_file.csv", ancestry_rows)
        depth_rows = [
            {
                "sample": sample,
                "truth_category": truth,
                "tau_ancestry_depth": depth,
                "nearest_tau_pdg": pdg,
                "count": count,
            }
            for (sample, truth, depth, pdg), count in sorted(
                depth_counts.items(), key=lambda item: (SAMPLES.index(item[0][0]), TRUTH_ORDER.index(item[0][1]), item[0][2], item[0][3])
            )
        ]
        write_csv(staging / "validation/nearest_tau_summary.csv", depth_rows)

        print("PLOT 18 tau-origin absolute matrices", flush=True)
        for population in POPULATIONS:
            target = staging / reco_selection(population)
            for sample in SAMPLES:
                for definition in DEFINITIONS:
                    stem = target / f"migration_matrix_{sample}_{definition}_{population}"
                    plot_absolute_tau(stem, matrices[(sample, definition, population)], sample, definition, population)
        print("PLOT 2 tau-origin contact sheets", flush=True)
        for population in POPULATIONS:
            stem = staging / "contact_sheets" / f"migration_matrix_grid_W_P8_G_Ldirect_Lancestor_{population}"
            plot_contact_tau(stem, matrices, population)

        root_summary = json.loads(existing_summary.read_text())
        limits = root_summary["difference_color_limits_pp"]
        difference_specs = []
        for population in POPULATIONS:
            for sample in SAMPLES:
                for left in ("Ldirect", "Lancestor"):
                    values = matrices[(sample, left, population)] - matrices[(sample, "G", population)]
                    difference_specs.append(("truth_definition_effect", sample, left, "G", population, values))
            for definition in DEFINITIONS:
                for left, right in (("P8C", "W"), ("P8O", "W"), ("P8O", "P8C")):
                    values = matrices[(left, definition, population)] - matrices[(right, definition, population)]
                    difference_specs.append(("sample_effect", definition, left, right, population, values))
        difference_rows = []
        print(f"PLOT {len(difference_specs)} tau-origin difference matrices", flush=True)
        for family, scope, left, right, population, values in difference_specs:
            selection = reco_selection(population)
            limit = float(limits[family])
            if family == "truth_definition_effect":
                title = f"{scope}: {left} − {right}"
                filename = f"difference_{scope}_{left}_minus_{right}_{population}"
                left_sample = right_sample = scope
                left_definition, right_definition = left, right
            else:
                title = f"{left} − {right}: {scope}"
                filename = f"difference_{left}_minus_{right}_{scope}_{population}"
                left_sample, right_sample = left, right
                left_definition = right_definition = scope
            plot_difference_tau(staging / "difference_matrices" / filename, values, title, population, limit)
            for i, truth in enumerate(TRUTH_ORDER):
                for j, outcome in enumerate(outcomes(selection)):
                    difference_rows.append({
                        "comparison_family": family,
                        "population": population,
                        "left_sample": left_sample,
                        "right_sample": right_sample,
                        "left_truth_definition": left_definition,
                        "right_truth_definition": right_definition,
                        "truth_category": truth,
                        "reco_category": outcome,
                        "difference_fraction": f"{values[i, j] / 100.0:.17g}",
                        "difference_percentage_points": f"{values[i, j]:.15g}",
                        "family_color_limit_pp": limit,
                    })
        write_csv(staging / "tables/difference_matrices_tau_origin.csv", difference_rows)

        validation = {
            "status": "PASS",
            "processed_frozen_files": {sample: sum(row["sample"] == sample for row in ancestry_rows) for sample in SAMPLES},
            "selected_truth_duplicates": 0,
            "invalid_parent_references": 0,
            "ancestry_cycles": 0,
            "stored_vs_recursive_tau_flag_mismatches": sum(int(row["stored_vs_recursive_mismatches"]) for row in ancestry_rows),
            "all_truth_regression_cells": len(regression_rows),
            "all_truth_regression_mismatches": sum(not row["pass"] for row in regression_rows),
            "tau_origin_matrix_rows": len(matrix_validation),
            "tau_origin_matrix_sum_failures": sum(not row["pass"] for row in matrix_validation),
            "common_denominator_across_definitions": True,
            "common_denominator_inclusive_fiducial": True,
            "association_rerun": False,
        }
        summary = {
            "status": "PASS",
            "definition": "selected truth MCParticle has at least one stored recursive parent ancestor with abs(PDG)==15",
            "output_root": str(TAU_ROOT),
            "populations": ["inclusive_all", "fiducial_all", *POPULATIONS],
            "tau_origin_denominators": {
                sample: {truth: tau_denominators[(sample, truth)] for truth in TRUTH_ORDER}
                for sample in SAMPLES
            },
            "matrix_counts": {"inclusive_tau_origin": 9, "fiducial_tau_origin": 9, "difference": len(difference_specs)},
            "contact_sheets": [
                str(TAU_ROOT / "contact_sheets/migration_matrix_grid_W_P8_G_Ldirect_Lancestor_inclusive_tau_origin.png"),
                str(TAU_ROOT / "contact_sheets/migration_matrix_grid_W_P8_G_Ldirect_Lancestor_fiducial_tau_origin.png"),
            ],
            "difference_color_limits_pp": limits,
            "validation": validation,
            "source_products": {
                "G": {sample: str(path) for sample, path in G_PATHS.items()},
                "L_manifests": {sample: str(path) for sample, path in MANIFESTS.items()},
                "all_truth_matrix_table": str(existing_table),
            },
            "forbidden_steps": {
                "linker_rerun": False,
                "HitAnalysis": False,
                "reconstruction": False,
                "simulation": False,
                "generator": False,
                "association_definition_changed": False,
            },
        }
        (staging / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        (staging / "validation/validation_summary.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
        os.replace(staging, TAU_ROOT)
        print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    except Exception:
        print(f"FAILED staging retained for audit: {staging}", file=sys.stderr, flush=True)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    cli = parser.parse_args()
    configure(cli.config)
    sys.argv = [sys.argv[0], "--n-workers", str(cli.workers)]
    main()
