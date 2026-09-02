#!/usr/bin/env python3
"""Build truth-normalized migration matrices from frozen FCC-tau products.

This is read-only scientific postprocessing.  It does not run matching,
linking, HitAnalysis, reconstruction, simulation, or generation.  G outcomes
come from the historical dR parquets.  L_direct and L_ancestor outcomes are
the frozen PFO assignments inverted with the frozen representative-PFO rule.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pyarrow.parquet as pq
import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from modules.fcc_truth_definitions import RECO_PID_BY_ABS_PDG  # noqa: E402
from modules.fcc_workflow_interface import truthlink_representative  # noqa: E402


FINAL = REPORT = FID = LANC = None
HIST_SCRIPT = HIST_COMPARISON = None
G_PATHS = {}
MANIFESTS = {}
THETA_MIN = math.radians(1.0)
THETA_MAX = math.radians(179.0)


def configure(path: Path) -> None:
    global FINAL, REPORT, FID, LANC, HIST_SCRIPT, HIST_COMPARISON, G_PATHS, MANIFESTS
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict) or data.get("schema_version") != "fcc_migration_matrix_inputs_v1":
        raise ValueError("unsupported migration-matrix configuration")
    FINAL = Path(data["output_root"]).expanduser().resolve()
    REPORT = Path(data["report_path"]).expanduser().resolve()
    FID = Path(data["validation_inputs"]["fiducial_tables"]).expanduser().resolve()
    LANC = Path(data["validation_inputs"]["lancestor_tables"]).expanduser().resolve()
    HIST_SCRIPT = Path(data["provenance"]["historical_definition_source"]).expanduser().resolve()
    HIST_COMPARISON = Path(data["provenance"]["historical_comparison_reference"]).expanduser().resolve()
    G_PATHS = {sample: Path(value).expanduser().resolve() for sample, value in data["g_products"].items()}
    MANIFESTS = {sample: Path(value).expanduser().resolve() for sample, value in data["workflow_manifests"].items()}
    if set(G_PATHS) != {"W", "P8C", "P8O"} or set(MANIFESTS) != set(G_PATHS):
        raise ValueError("G products and workflow manifests must define W, P8C and P8O")

SAMPLES = ("W", "P8C", "P8O")
DEFINITIONS = ("G", "Ldirect", "Lancestor")
SELECTIONS = ("inclusive", "fiducial")
TRUTH_LABELS = {
    "electron": r"$e$",
    "muon": r"$\mu$",
    "photon": r"$\gamma$",
    "charged_pion": r"$\pi^{\pm}$",
}
OUTCOME_LABELS = {
    "electron": r"$e$",
    "muon": r"$\mu$",
    "photon": r"$\gamma$",
    "charged_pion": r"$\pi^{\pm}$",
    "K0S": r"$K^0_S$",
    "association_unmatched": "unmatched",
    "neutron": r"$n$",
    "Lambda": r"$\Lambda$",
    "ambiguous_multiple_pfo": "ambiguous\nmultiple PFO",
    "assigned_but_fails_reco_selection": "assigned, fails\nreco selection",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and fields is None:
        raise ValueError(f"cannot infer fields for empty table: {path}")
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


TRUTH_BY_PDG = {pdg: RECO_PID_BY_ABS_PDG[pdg] for pdg in (11, 13, 22, 211)}
HIST_DEST_PDGS = [11, 13, 22, 211, 310, 999, 2112, 3122]
HIST_DEST_LABELS = {
    11: "e", 13: "mu", 22: "gamma", 211: "pi", 310: "K0S",
    999: "unmatched", 2112: "n", 3122: "Lambda",
}
DEST_NAME = dict(RECO_PID_BY_ABS_PDG)
DEST_NAME[999] = "association_unmatched"

BASE_OUTCOMES = tuple(DEST_NAME[pdg] for pdg in HIST_DEST_PDGS) + ("ambiguous_multiple_pfo",)
TRUTH_ORDER = tuple(TRUTH_BY_PDG[pdg] for pdg in TRUTH_BY_PDG)


def outcomes(selection: str) -> tuple[str, ...]:
    if selection == "inclusive":
        return BASE_OUTCOMES
    return BASE_OUTCOMES + ("assigned_but_fails_reco_selection",)


def category_from_pdg(pdg: int) -> str:
    value = abs(int(pdg))
    if value not in DEST_NAME or value == 999:
        raise AssertionError(f"reconstructed PDG outside historical categories: {pdg}")
    return DEST_NAME[value]


def fiducial_outcome(inclusive: str, pfo: dict | None) -> str:
    if inclusive in {"association_unmatched", "ambiguous_multiple_pfo"}:
        return inclusive
    if pfo is None:
        raise AssertionError("assigned outcome without representative PFO")
    energy = float(pfo["pfo_energy"])
    theta = float(pfo["pfo_theta"])
    if energy > 1.0 and THETA_MIN < theta < THETA_MAX:
        return inclusive
    return "assigned_but_fails_reco_selection"


def manifest_index(sample: str) -> dict[str, dict[str, str]]:
    rows = read_csv(MANIFESTS[sample])
    expected = {"W": 996, "P8C": 100, "P8O": 18}[sample]
    if len(rows) != expected:
        raise AssertionError(f"{sample} manifest rows {len(rows)} != {expected}")
    key_col = "source_rec" if sample == "W" else "input_REC"
    result = {Path(row[key_col]).name: row for row in rows}
    if len(result) != len(rows):
        raise AssertionError(f"duplicate source basename in {sample} manifest")
    return result


def normalize_weight(value) -> int:
    return 0 if value is None else int(value)


def load_frozen_l(sample: str, meta: dict[str, str]):
    """Load one frozen PFO assignment block; no association is recomputed."""
    direct_by_truth = defaultdict(list)
    ancestor_by_truth = defaultdict(list)
    pfo_by_key = {}
    observed = set()
    if sample == "W":
        direct_cols = [
            "event_in_file", "pfo_index", "pfo_type", "pfo_energy", "pfo_theta",
            "truthlink_status", "assigned_mc_index", "track_permille", "cluster_permille",
        ]
        ancestor_cols = ["event_in_file", "pfo_index", "ancestor_mc_index"]
        direct_rows = pq.read_table(meta["direct"], columns=direct_cols).to_pylist()
        ancestor_rows = pq.read_table(meta["ancestor_out"], columns=ancestor_cols).to_pylist()
        ancestor_map = {
            (int(row["event_in_file"]), int(row["pfo_index"])): row["ancestor_mc_index"]
            for row in ancestor_rows
        }
        if len(ancestor_map) != len(ancestor_rows) or len(direct_rows) != len(ancestor_rows):
            raise AssertionError("W direct/ancestor frozen PFO join mismatch")
        rows = []
        for row in direct_rows:
            key = (int(row["event_in_file"]), int(row["pfo_index"]))
            if key not in ancestor_map:
                raise AssertionError("missing W ancestor PFO key")
            item = dict(row)
            item["direct_status"] = row["truthlink_status"]
            item["direct_mc_index"] = row["assigned_mc_index"]
            item["ancestor_mc_index"] = ancestor_map[key]
            rows.append(item)
    else:
        columns = [
            "event_in_file", "pfo_index", "pfo_type", "pfo_energy", "pfo_theta",
            "direct_status", "direct_mc_index", "ancestor_mc_index",
            "track_permille", "cluster_permille",
        ]
        rows = pq.read_table(meta["ancestor_output"], columns=columns).to_pylist()

    for raw in rows:
        event = int(raw["event_in_file"])
        pfo_index = int(raw["pfo_index"])
        pdg = abs(int(raw["pfo_type"]))
        observed.add(pdg)
        if pdg not in HIST_DEST_PDGS or pdg == 999:
            raise AssertionError(f"{sample} frozen PFO type {pdg} outside historical reco inventory")
        row = {
            "pfo_index": pfo_index,
            "pfo_type": int(raw["pfo_type"]),
            "pfo_energy": float(raw["pfo_energy"]),
            "pfo_theta": float(raw["pfo_theta"]),
            "track_permille": normalize_weight(raw["track_permille"]),
            "cluster_permille": normalize_weight(raw["cluster_permille"]),
        }
        key = (event, pfo_index)
        if key in pfo_by_key:
            raise AssertionError(f"duplicate frozen PFO key {sample} {key}")
        pfo_by_key[key] = row
        if raw["direct_status"] == "assigned" and raw["direct_mc_index"] is not None:
            direct_by_truth[(event, int(raw["direct_mc_index"]))].append(row)
        if raw["ancestor_mc_index"] is not None:
            ancestor_by_truth[(event, int(raw["ancestor_mc_index"]))].append(row)
    return pfo_by_key, direct_by_truth, ancestor_by_truth, observed


def l_outcome(rows: list[dict], pfo_by_key: dict, event: int) -> tuple[str, dict | None]:
    representative = truthlink_representative(rows)
    if representative["status"] == "unmatched":
        return "association_unmatched", None
    if representative["status"] == "ambiguous_multiple_pfo":
        return "ambiguous_multiple_pfo", None
    key = (event, int(representative["pfo_index"]))
    pfo = pfo_by_key.get(key)
    if pfo is None:
        raise AssertionError(f"representative PFO missing: {key}")
    return category_from_pdg(pfo["pfo_type"]), pfo


def iter_g_groups(path: Path) -> Iterable[tuple[str, list[tuple]]]:
    columns = ["source_file", "event_in_file", "gen", "reco", "Gen_pid", "Reco_pid"]
    current_source = None
    current_rows = []
    seen_sources = set()
    for batch in pq.ParquetFile(path).iter_batches(columns=columns, batch_size=131_072):
        data = {name: batch.column(i).to_pylist() for i, name in enumerate(columns)}
        for values in zip(*(data[name] for name in columns)):
            source, event, gen, reco, gen_pdg, reco_pdg = values
            source = str(source)
            if current_source is None:
                current_source = source
            if source != current_source:
                if current_source in seen_sources:
                    raise AssertionError(f"non-contiguous source group in {path}: {current_source}")
                seen_sources.add(current_source)
                yield current_source, current_rows
                current_source = source
                current_rows = []
            current_rows.append((int(event), int(gen), int(reco), int(gen_pdg), int(reco_pdg)))
    if current_source is not None:
        if current_source in seen_sources:
            raise AssertionError(f"repeated final source group in {path}: {current_source}")
        yield current_source, current_rows


def process_all():
    counts = Counter()
    denominators = Counter()
    file_counts = Counter()
    observed_pfo_pdgs = defaultdict(set)
    manifest_maps = {sample: manifest_index(sample) for sample in SAMPLES}

    for sample in SAMPLES:
        print(f"PROCESS {sample}: {G_PATHS[sample]}", flush=True)
        expected_files = {"W": 996, "P8C": 100, "P8O": 18}[sample]
        for file_number, (source_name, truth_rows) in enumerate(iter_g_groups(G_PATHS[sample]), start=1):
            meta = manifest_maps[sample].get(Path(source_name).name)
            if meta is None:
                raise AssertionError(f"{sample} G source absent from frozen manifest: {source_name}")
            pfo_map, direct_by_truth, ancestor_by_truth, observed = load_frozen_l(sample, meta)
            observed_pfo_pdgs[sample].update(observed)
            seen_truth = set()
            for event, gen, reco, gen_pdg, reco_pdg in truth_rows:
                truth_pdg = abs(gen_pdg)
                if truth_pdg not in TRUTH_BY_PDG:
                    continue
                truth = TRUTH_BY_PDG[truth_pdg]
                truth_key = (event, gen)
                if truth_key in seen_truth:
                    raise AssertionError(f"duplicate selected truth key in {sample}/{source_name}: {truth_key}")
                seen_truth.add(truth_key)
                denominators[(sample, truth)] += 1

                if abs(reco_pdg) == 999:
                    g_inc, g_pfo = "association_unmatched", None
                else:
                    pfo_key = (event, reco)
                    g_pfo = pfo_map.get(pfo_key)
                    if g_pfo is None:
                        raise AssertionError(f"G representative missing from frozen PFOs: {sample} {source_name} {pfo_key}")
                    g_inc = category_from_pdg(reco_pdg)
                    if abs(int(g_pfo["pfo_type"])) != abs(reco_pdg):
                        raise AssertionError(f"G reco category/PFO type mismatch: {sample} {pfo_key}")
                g_fid = fiducial_outcome(g_inc, g_pfo)

                direct_inc, direct_pfo = l_outcome(direct_by_truth.get(truth_key, []), pfo_map, event)
                direct_fid = fiducial_outcome(direct_inc, direct_pfo)
                ancestor_inc, ancestor_pfo = l_outcome(ancestor_by_truth.get(truth_key, []), pfo_map, event)
                ancestor_fid = fiducial_outcome(ancestor_inc, ancestor_pfo)

                for definition, inc, fid in (
                    ("G", g_inc, g_fid),
                    ("Ldirect", direct_inc, direct_fid),
                    ("Lancestor", ancestor_inc, ancestor_fid),
                ):
                    counts[(sample, definition, "inclusive", truth, inc)] += 1
                    counts[(sample, definition, "fiducial", truth, fid)] += 1
            file_counts[sample] += 1
            if file_number % 50 == 0 or file_number == expected_files:
                print(f"  {sample}: {file_number}/{expected_files} frozen files", flush=True)
        if file_counts[sample] != expected_files:
            raise AssertionError(f"{sample} G file groups {file_counts[sample]} != {expected_files}")
    return counts, denominators, file_counts, observed_pfo_pdgs


def collapsed_count(counts, sample, definition, selection, truth, outcome):
    if outcome == "gamma":
        outcome = "photon"
    if outcome == "other":
        if truth == "photon":
            excluded = {"photon", "electron", "charged_pion", "association_unmatched", "ambiguous_multiple_pfo", "assigned_but_fails_reco_selection"}
        elif truth == "charged_pion":
            excluded = {"charged_pion", "electron", "muon", "photon", "K0S", "neutron", "association_unmatched", "ambiguous_multiple_pfo", "assigned_but_fails_reco_selection"}
        else:
            raise ValueError(truth)
        return sum(
            counts[(sample, definition, selection, truth, item)]
            for item in outcomes(selection) if item not in excluded
        )
    return counts[(sample, definition, selection, truth, outcome)]


def validate_historical_g(counts, denominators) -> list[dict]:
    """Require exact reproduction of all original inclusive G migration rows."""
    campaign = {"W": "WHIZARD", "P8C": "PYTHIA8_COLLAB", "P8O": "PYTHIA8_OURRECO"}
    rows = []
    for pdg, truth in TRUTH_BY_PDG.items():
        source = HIST_COMPARISON / f"tables/migration_{truth}.csv"
        for expected in read_csv(source):
            reco_pdg = int(expected["reco_pdg"])
            outcome = DEST_NAME[reco_pdg]
            for sample in SAMPLES:
                prefix = campaign[sample]
                expected_count = int(expected[f"{prefix}_count"])
                expected_den = int(expected[f"{prefix}_N_truth"])
                actual_count = counts[(sample, "G", "inclusive", truth, outcome)]
                actual_den = denominators[(sample, truth)]
                passed = actual_count == expected_count and actual_den == expected_den
                rows.append({
                    "sample": sample,
                    "truth_category": truth,
                    "reco_category": outcome,
                    "expected_count": expected_count,
                    "actual_count": actual_count,
                    "expected_denominator": expected_den,
                    "actual_denominator": actual_den,
                    "pass": passed,
                })
                if not passed:
                    raise AssertionError(f"historical G reproduction failed: {rows[-1]}")
    return rows


def validate_headline_table(counts, denominators, path: Path, sample_name_map: dict[str, str], definition_map: dict[str, str], truths: set[str]):
    rows = []
    useful = set(BASE_OUTCOMES) | {"gamma", "other", "assigned_but_fails_reco_selection"}
    for expected in read_csv(path):
        truth = expected["species"]
        if truth not in truths:
            continue
        metric = expected["metric"]
        if metric not in useful:
            continue
        sample_method = expected["sample_or_method"]
        sample = sample_name_map.get(sample_method)
        definition = definition_map.get(sample_method)
        if sample is None or definition is None:
            continue
        selection = expected["selection"]
        actual = collapsed_count(counts, sample, definition, selection, truth, metric)
        expected_count = int(expected["count"])
        expected_den = int(expected["truth_denominator"])
        passed = actual == expected_count and denominators[(sample, truth)] == expected_den
        record = {
            "source": str(path), "sample": sample, "truth_definition": definition,
            "selection": selection, "truth_category": truth, "outcome": metric,
            "expected_count": expected_count, "actual_count": actual,
            "expected_denominator": expected_den, "actual_denominator": denominators[(sample, truth)],
            "pass": passed,
        }
        rows.append(record)
        if not passed:
            raise AssertionError(f"headline anchor failed: {record}")
    return rows


def validate_lancestor_tables(counts, denominators):
    rows = []
    for truth, prefix in (("photon", "photon"), ("charged_pion", "pion")):
        for selection in SELECTIONS:
            path = LANC / f"{prefix}_migration_{selection}_Lancestor.csv"
            for expected in read_csv(path):
                sample = expected["sample"]
                outcome = expected["outcome"]
                actual = collapsed_count(counts, sample, "Lancestor", selection, truth, outcome)
                expected_count = int(expected["count"])
                expected_den = int(expected["denominator"])
                passed = actual == expected_count and denominators[(sample, truth)] == expected_den
                record = {
                    "source": str(path), "sample": sample, "truth_definition": "Lancestor",
                    "selection": selection, "truth_category": truth, "outcome": outcome,
                    "expected_count": expected_count, "actual_count": actual,
                    "expected_denominator": expected_den, "actual_denominator": denominators[(sample, truth)],
                    "pass": passed,
                }
                rows.append(record)
                if not passed:
                    raise AssertionError(f"Lancestor anchor failed: {record}")
    return rows


def validate_authoritative_anchors(counts, denominators):
    rows = []
    rows += validate_headline_table(
        counts, denominators, FID / "headline_photon_W_P8.csv",
        {"W_G": "W", "P8C_G": "P8C", "P8O_G": "P8O"},
        {"W_G": "G", "P8C_G": "G", "P8O_G": "G"}, {"photon"},
    )
    rows += validate_headline_table(
        counts, denominators, FID / "headline_pion_W_P8.csv",
        {"W_G": "W", "P8C_G": "P8C", "P8O_G": "P8O"},
        {"W_G": "G", "P8C_G": "G", "P8O_G": "G"}, {"charged_pion"},
    )
    rows += validate_headline_table(
        counts, denominators, FID / "headline_photon_G_Ldirect_Lancestor.csv",
        {"G": "W", "Ldirect": "W", "Lancestor": "W"},
        {"G": "G", "Ldirect": "Ldirect", "Lancestor": "Lancestor"}, {"photon"},
    )
    rows += validate_headline_table(
        counts, denominators, FID / "headline_pion_G_Ldirect_Lancestor.csv",
        {"G": "W", "Ldirect": "W", "Lancestor": "W"},
        {"G": "G", "Ldirect": "Ldirect", "Lancestor": "Lancestor"}, {"charged_pion"},
    )
    rows += validate_lancestor_tables(counts, denominators)
    return rows


def matrix_array(counts, denominators, sample, definition, selection):
    order = outcomes(selection)
    values = np.zeros((len(TRUTH_ORDER), len(order)), dtype=float)
    raw = np.zeros_like(values, dtype=np.int64)
    for i, truth in enumerate(TRUTH_ORDER):
        den = denominators[(sample, truth)]
        if den <= 0:
            raise AssertionError(f"zero truth denominator: {sample} {truth}")
        for j, outcome in enumerate(order):
            value = counts[(sample, definition, selection, truth, outcome)]
            raw[i, j] = value
            values[i, j] = 100.0 * value / den
    return raw, values


def cell_text(value: float, difference: bool = False) -> str:
    if abs(value) < 1e-14:
        return "0"
    if abs(value) < 0.01:
        return ("−" if value < 0 else "+" if difference else "") + "<0.01"
    if difference:
        return f"{value:+.2f}".replace("-", "−")
    return f"{value:.2f}"


def annotate(ax, values, difference=False, limit=100.0):
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = float(values[i, j])
            if difference:
                color = "white" if abs(value) > 0.55 * limit else "black"
            else:
                color = "white" if value > 52 else "black"
            ax.text(j, i, cell_text(value, difference), ha="center", va="center", fontsize=7.2, color=color)


def style_axis(ax, selection):
    order = outcomes(selection)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([OUTCOME_LABELS[x] for x in order], rotation=42, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(TRUTH_ORDER)))
    ax.set_yticklabels([TRUTH_LABELS[x] for x in TRUTH_ORDER], fontsize=10)
    ax.set_xlabel("Reconstructed outcome", fontsize=10)
    ax.set_ylabel("Selected truth particle", fontsize=10)
    ax.set_xticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(TRUTH_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)


def plot_absolute(path_stem: Path, values, sample, definition, selection):
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    image = ax.imshow(values, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    annotate(ax, values)
    style_axis(ax, selection)
    ax.set_title(f"{sample} — {definition.replace('Ldirect', r'$L_{direct}$').replace('Lancestor', r'$L_{ancestor}$')} — {selection}\ntruth-normalized migration fraction [%]", fontsize=13)
    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label("migration fraction [%]")
    fig.text(0.5, 0.012, "Truth-association-dependent reconstructed outcome matrix; each truth row sums to 100%.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_contact(path_stem: Path, matrices: dict, selection: str):
    fig, axes = plt.subplots(3, 3, figsize=(25, 13.5), sharex=True, sharey=True)
    image = None
    for i, sample in enumerate(SAMPLES):
        for j, definition in enumerate(DEFINITIONS):
            ax = axes[i, j]
            values = matrices[(sample, definition, selection)]
            image = ax.imshow(values, cmap="Blues", vmin=0, vmax=100, aspect="auto")
            annotate(ax, values)
            style_axis(ax, selection)
            if i != 2:
                ax.set_xlabel("")
            if j != 0:
                ax.set_ylabel("")
            label = definition.replace("Ldirect", r"$L_{direct}$").replace("Lancestor", r"$L_{ancestor}$")
            ax.set_title(f"{sample} — {label}", fontsize=14, fontweight="bold")
    fig.suptitle(f"W / P8 truth-normalized migration matrices — {selection}", fontsize=20, fontweight="bold")
    color_axis = fig.add_axes([0.925, 0.14, 0.012, 0.70])
    colorbar = fig.colorbar(image, cax=color_axis)
    colorbar.set_label("migration fraction [%]")
    fig.text(0.5, 0.008, "Rows: selected truth particle. Columns: reconstructed outcome. Fixed scale 0–100%; each row sums to 100%.", ha="center", fontsize=11)
    fig.subplots_adjust(left=0.055, right=0.90, bottom=0.12, top=0.91, wspace=0.12, hspace=0.28)
    fig.savefig(path_stem.with_suffix(".png"), dpi=240, bbox_inches="tight")
    fig.savefig(path_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_difference(path_stem: Path, values, title: str, selection: str, limit: float):
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    image = ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto")
    annotate(ax, values, difference=True, limit=limit)
    style_axis(ax, selection)
    ax.set_title(f"{title} — {selection}\ntruth-normalized migration difference [percentage points]", fontsize=13)
    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label("difference [percentage points]")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def nice_limit(value: float) -> float:
    return max(1.0, math.ceil(value * 2.0) / 2.0)


def top_changes(matrices, sample, left, right, selection, truth, n=3):
    i = TRUTH_ORDER.index(truth)
    delta = matrices[(sample, left, selection)][i] - matrices[(sample, right, selection)][i]
    pairs = sorted(zip(outcomes(selection), delta), key=lambda item: abs(item[1]), reverse=True)
    return pairs[:n]


def pct(matrices, sample, definition, selection, truth, outcome):
    return float(matrices[(sample, definition, selection)][TRUTH_ORDER.index(truth), outcomes(selection).index(outcome)])


def findings_markdown(matrices):
    lines = []
    for definition in DEFINITIONS:
        w_gamma = pct(matrices, "W", definition, "inclusive", "photon", "photon")
        c_gamma = pct(matrices, "P8C", definition, "inclusive", "photon", "photon")
        o_gamma = pct(matrices, "P8O", definition, "inclusive", "photon", "photon")
        w_un = pct(matrices, "W", definition, "inclusive", "photon", "association_unmatched")
        c_un = pct(matrices, "P8C", definition, "inclusive", "photon", "association_unmatched")
        o_un = pct(matrices, "P8O", definition, "inclusive", "photon", "association_unmatched")
        lines.append(
            f"- **{definition}:** photon→photon = {w_gamma:.3f}% / {c_gamma:.3f}% / {o_gamma:.3f}% and "
            f"photon→unmatched = {w_un:.3f}% / {c_un:.3f}% / {o_un:.3f}% for W / P8C / P8O."
        )
    pion_values = [pct(matrices, s, d, "inclusive", "charged_pion", "charged_pion") for s in SAMPLES for d in DEFINITIONS]
    lines.append(f"- Across all nine inclusive matrices, charged-pion→charged-pion spans {min(pion_values):.3f}%–{max(pion_values):.3f}%.")
    fail_values = [pct(matrices, s, d, "fiducial", "photon", "assigned_but_fails_reco_selection") for s in SAMPLES for d in DEFINITIONS]
    lines.append(f"- The photon assigned-but-fails-reco-selection state spans {min(fail_values):.3f}%–{max(fail_values):.3f}%; unmatched and ambiguity counts are unchanged by the reco fiducial cut.")
    return "\n".join(lines)


def make_report(final: Path, matrices, denominators, validation, diff_count):
    w_direct = top_changes(matrices, "W", "Ldirect", "G", "inclusive", "photon")
    w_ancestor = top_changes(matrices, "W", "Lancestor", "G", "inclusive", "photon")
    def change_text(items):
        return ", ".join(f"{name} {value:+.3f} pp" for name, value in items)
    text = f"""# W/P8 migration matrices: G, L_direct and L_ancestor

Date: 2026-08-28
Status: **PASS**

## Definition recovered from the historical comparison

- Historical source: `{HIST_SCRIPT}`.
- Historical products: `{HIST_COMPARISON}`.
- Truth categories, in frozen order: `{', '.join(TRUTH_ORDER)}`.
- Reconstructed categories, in frozen order: `{', '.join(BASE_OUTCOMES[:-1])}`; `ambiguous_multiple_pfo` is appended for common L semantics.
- Axis orientation in the new heatmaps: x = reconstructed outcome, y = selected truth particle. This preserves the historical reco-category x-axis while stacking the four historical truth rows.
- Normalization: independently per truth row, `N(truth X → reco outcome Y) / N(selected truth X)`.
- Historical plots were one grouped-bar migration plot per truth species, not a saved 2D cell-annotated confusion matrix. Their scientific definition, category order, labels, normalization and x-axis orientation were recovered exactly. Exact historical 2D cosmetics therefore carry the marker **ORIGINAL_STYLE_NOT_RECOVERED**; the new heatmaps use a fixed clean common style.
- Historical bar plots did not print percentages inside cells (there were no cells) and used fractional y axes with campaign colors. The new heatmaps print percentages and use a fixed 0–100% scale as requested.

## Frozen association and selection semantics

- **G:** already-frozen geometric selected-truth outcome from the historical dR parquet; no d_theta_phi matching was rerun.
- **L_direct:** immediate detector-level MC contributor from `truthlink_assignment_v1`, inverted with the frozen representative-PFO rule (track branch priority, then winning T/C; exact terminal ties remain ambiguous).
- **L_ancestor:** nearest unique selected generator-level ancestor from `truthlink_ancestor_assignment_v1`, with the same frozen inversion rule.
- **Inclusive:** representative physical PID outcome, `association_unmatched`, or `ambiguous_multiple_pfo`.
- **Fiducial:** representative chosen first; then strict `PFO.getEnergy() > 1.0 GeV` and `1° < theta_reco < 179°`. A representative that fails becomes `assigned_but_fails_reco_selection`; it is never replaced. No truth-level cut is applied.
- These are **truth-association-dependent reconstructed outcome matrices**. An off-diagonal L entry is not automatically a detector PID misidentification.

## Validation

- All 18 × 4 truth rows sum to 100% within floating precision: **{validation['all_matrix_sums_100']}**.
- Original inclusive G counts and denominators reproduced exactly for all three samples, four truth species and eight historical destinations: **{validation['historical_g_exact']}**.
- Authoritative G fiducial photon/pion anchors reproduced exactly: **{validation['g_fiducial_anchors']}**.
- Authoritative W G/L_direct/L_ancestor photon/pion anchors reproduced exactly: **{validation['w_truth_definition_anchors']}**.
- Authoritative cross-sample L_ancestor inclusive and fiducial photon/pion rows reproduced exactly: **{validation['lancestor_anchors']}**.
- G/L_direct/L_ancestor use the same selected-truth denominator by construction within every sample/species: **{validation['common_denominators']}**.

## Main numerical readout

{findings_markdown(matrices)}

## Scientific questions

1. **G → L_direct → L_ancestor in W.** The largest inclusive photon-row shifts for L_direct−G are {change_text(w_direct)}. For L_ancestor−G they are {change_text(w_ancestor)}. The direct and ancestry definitions therefore alter association/provenance semantics visibly, especially unmatched and photon/electron descendant outcomes.
2. **P8C and P8O.** The same qualitative direct-versus-ancestor redistribution occurs and the two P8 reconstruction chains remain close element by element; see the P8O−P8C difference matrices.
3. **Common G, W versus P8.** The large photon difference is dominated by the photon diagonal and the compensating unmatched outcome; off-diagonal physical PID migrations are much smaller.
4. **Common L_ancestor.** The W/P8 difference remains dominated by photon→photon versus photon→unmatched, while electron and other physical migrations are subleading.
5. **Unmatched or another PID?** Predominantly unmatched. The difference matrices show no comparably large transfer into another physical PID category.
6. **Charged pions.** The charged-pion diagonal is high and comparatively stable across samples and truth definitions. Changes under L_ancestor include legitimate descendant promotion and must not automatically be labelled PID misidentification.
7. **Fiducial selection.** It mainly transfers assigned entries into `assigned_but_fails_reco_selection`; `association_unmatched` and `ambiguous_multiple_pfo` remain separate and numerically unchanged.

## Outputs

- Individual inclusive matrices: `{final / 'inclusive'}`
- Individual fiducial matrices: `{final / 'fiducial'}`
- Difference matrices ({diff_count}): `{final / 'difference_matrices'}`
- Contact sheets: `{final / 'contact_sheets'}`
- Machine-readable tables: `{final / 'tables'}`
- Validation: `{final / 'validation'}`
- Summary: `{final / 'summary.json'}`

No linker, HitAnalysis, reconstruction, simulation or generator step was run. No frozen scientific definition was changed.
"""
    tmp = REPORT.with_name(f".{REPORT.name}.partial.{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, REPORT)


def main():
    if FINAL.exists():
        raise FileExistsError(f"refusing to overwrite {FINAL}")
    if REPORT.exists():
        raise FileExistsError(f"refusing to overwrite {REPORT}")
    for path in [*G_PATHS.values(), *MANIFESTS.values()]:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    counts, denominators, file_counts, observed = process_all()
    print("VALIDATE original inclusive G matrices", flush=True)
    historical_g = validate_historical_g(counts, denominators)
    print("VALIDATE authoritative fiducial and L_ancestor anchors", flush=True)
    anchors = validate_authoritative_anchors(counts, denominators)

    parent = FINAL.parent
    staging = Path(tempfile.mkdtemp(prefix=f".{FINAL.name}.partial.", dir=parent))
    try:
        for name in ("inclusive", "fiducial", "difference_matrices", "contact_sheets", "tables", "validation", "reports"):
            (staging / name).mkdir()

        all_rows = []
        validation_rows = []
        matrices = {}
        for selection in SELECTIONS:
            for sample in SAMPLES:
                for definition in DEFINITIONS:
                    raw, values = matrix_array(counts, denominators, sample, definition, selection)
                    matrices[(sample, definition, selection)] = values
                    for i, truth in enumerate(TRUTH_ORDER):
                        denominator = denominators[(sample, truth)]
                        matrix_sum = float(values[i].sum())
                        validation_rows.append({
                            "sample": sample, "truth_definition": definition, "selection": selection,
                            "truth_category": truth, "denominator": denominator,
                            "matrix_sum": f"{matrix_sum:.15g}",
                            "difference_from_100": f"{matrix_sum - 100.0:.15g}",
                            "count_sum": int(raw[i].sum()), "pass": int(raw[i].sum()) == denominator and abs(matrix_sum - 100.0) < 1e-10,
                        })
                        for j, outcome in enumerate(outcomes(selection)):
                            all_rows.append({
                                "sample": sample, "truth_definition": definition, "selection": selection,
                                "truth_category": truth, "reco_category": outcome,
                                "count": int(raw[i, j]), "truth_denominator": denominator,
                                "fraction": f"{values[i, j] / 100.0:.17g}", "percent": f"{values[i, j]:.15g}",
                            })
        if not all(row["pass"] for row in validation_rows):
            raise AssertionError("matrix row sum validation failed")

        write_csv(staging / "tables/migration_matrices_all.csv", all_rows)
        write_csv(staging / "tables/photon_migration_rows.csv", [row for row in all_rows if row["truth_category"] == "photon"])
        write_csv(staging / "tables/charged_pion_migration_rows.csv", [row for row in all_rows if row["truth_category"] == "charged_pion"])
        write_csv(staging / "validation/matrix_validation.csv", validation_rows)
        write_csv(staging / "validation/historical_G_reproduction.csv", historical_g)
        write_csv(staging / "validation/authoritative_anchor_reproduction.csv", anchors)

        print("PLOT 18 absolute matrices", flush=True)
        for selection in SELECTIONS:
            target = staging / selection
            for sample in SAMPLES:
                for definition in DEFINITIONS:
                    stem = target / f"migration_matrix_{sample}_{definition}_{selection}"
                    plot_absolute(stem, matrices[(sample, definition, selection)], sample, definition, selection)

        print("PLOT contact sheets", flush=True)
        for selection in SELECTIONS:
            stem = staging / "contact_sheets" / f"migration_matrix_grid_W_P8_G_Ldirect_Lancestor_{selection}"
            plot_contact(stem, matrices, selection)

        difference_specs = []
        within_values = []
        for selection in SELECTIONS:
            for sample in SAMPLES:
                for left in ("Ldirect", "Lancestor"):
                    value = matrices[(sample, left, selection)] - matrices[(sample, "G", selection)]
                    difference_specs.append(("truth_definition_effect", sample, left, "G", selection, value))
                    within_values.append(value)
        sample_values = []
        for selection in SELECTIONS:
            for definition in DEFINITIONS:
                for left, right in (("P8C", "W"), ("P8O", "W"), ("P8O", "P8C")):
                    value = matrices[(left, definition, selection)] - matrices[(right, definition, selection)]
                    difference_specs.append(("sample_effect", definition, left, right, selection, value))
                    sample_values.append(value)
        within_limit = nice_limit(max(float(np.max(np.abs(value))) for value in within_values))
        sample_limit = nice_limit(max(float(np.max(np.abs(value))) for value in sample_values))
        difference_rows = []
        print(f"PLOT {len(difference_specs)} difference matrices", flush=True)
        for family, scope, left, right, selection, values in difference_specs:
            limit = within_limit if family == "truth_definition_effect" else sample_limit
            if family == "truth_definition_effect":
                title = f"{scope}: {left} − {right}"
                filename = f"difference_{scope}_{left}_minus_{right}_{selection}"
                left_sample = right_sample = scope
                left_definition, right_definition = left, right
            else:
                title = f"{left} − {right}: {scope}"
                filename = f"difference_{left}_minus_{right}_{scope}_{selection}"
                left_sample, right_sample = left, right
                left_definition = right_definition = scope
            plot_difference(staging / "difference_matrices" / filename, values, title, selection, limit)
            for i, truth in enumerate(TRUTH_ORDER):
                for j, outcome in enumerate(outcomes(selection)):
                    difference_rows.append({
                        "comparison_family": family, "selection": selection,
                        "left_sample": left_sample, "right_sample": right_sample,
                        "left_truth_definition": left_definition, "right_truth_definition": right_definition,
                        "truth_category": truth, "reco_category": outcome,
                        "difference_fraction": f"{values[i, j] / 100.0:.17g}",
                        "difference_percentage_points": f"{values[i, j]:.15g}",
                        "family_color_limit_pp": limit,
                    })
        write_csv(staging / "tables/difference_matrices.csv", difference_rows)

        original_audit = {
            "status": "PASS",
            "historical_source": str(HIST_SCRIPT),
            "historical_source_sha256": sha256(HIST_SCRIPT),
            "historical_products": str(HIST_COMPARISON),
            "truth_categories": list(TRUTH_ORDER),
            "reco_category_pdgs": HIST_DEST_PDGS,
            "reco_category_labels": [HIST_DEST_LABELS[pdg] for pdg in HIST_DEST_PDGS],
            "axis_orientation": {"x": "reconstructed outcome", "y": "selected truth particle"},
            "normalization": "N(truth X -> reco outcome Y) / N(selected truth X), independently per truth row",
            "historical_rendering": "one grouped-bar plot per truth species; no 2D cell-annotated matrix found",
            "exact_historical_2d_style": "ORIGINAL_STYLE_NOT_RECOVERED",
            "new_cell_format": "0 for zero; <0.01 for nonzero below 0.01%; otherwise 2 decimals",
            "absolute_color_scale_percent": [0, 100],
        }
        (staging / "reports/original_format_audit.json").write_text(json.dumps(original_audit, indent=2, sort_keys=True) + "\n")

        validation = {
            "all_matrix_sums_100": True,
            "historical_g_exact": all(row["pass"] for row in historical_g),
            "g_fiducial_anchors": all(row["pass"] for row in anchors if "headline_" in row["source"] and row["truth_definition"] == "G"),
            "w_truth_definition_anchors": all(row["pass"] for row in anchors if "G_Ldirect_Lancestor" in row["source"]),
            "lancestor_anchors": all(row["pass"] for row in anchors if "migration_" in Path(row["source"]).name),
            "common_denominators": True,
        }
        summary = {
            "status": "PASS",
            "version": "migration_matrices_W_P8_G_Ldirect_Lancestor_20260828",
            "output_root": str(FINAL),
            "matrix_definition": "truth-association-dependent reconstructed outcome; truth-row normalized",
            "truth_categories": list(TRUTH_ORDER),
            "inclusive_reco_categories": list(outcomes("inclusive")),
            "fiducial_reco_categories": list(outcomes("fiducial")),
            "axis_orientation": {"x": "reconstructed outcome", "y": "selected truth particle"},
            "selection": {"energy": "PFO.getEnergy() > 1.0 GeV", "theta": "1 deg < theta_reco < 179 deg", "truth_level_cut": False, "representative_before_cut": True},
            "matrix_counts": {"inclusive": 9, "fiducial": 9, "difference": len(difference_specs)},
            "contact_sheets": [
                str(FINAL / "contact_sheets/migration_matrix_grid_W_P8_G_Ldirect_Lancestor_inclusive.png"),
                str(FINAL / "contact_sheets/migration_matrix_grid_W_P8_G_Ldirect_Lancestor_fiducial.png"),
            ],
            "selected_truth_denominators": {sample: {truth: denominators[(sample, truth)] for truth in TRUTH_ORDER} for sample in SAMPLES},
            "processed_frozen_files": dict(file_counts),
            "observed_frozen_pfo_pdgs": {sample: sorted(values) for sample, values in observed.items()},
            "difference_color_limits_pp": {"truth_definition_effect": within_limit, "sample_effect": sample_limit},
            "original_format": original_audit,
            "validation": validation,
            "source_products": {
                "G": {sample: str(path) for sample, path in G_PATHS.items()},
                "L_manifests": {sample: str(path) for sample, path in MANIFESTS.items()},
                "authoritative_fiducial_tables": str(FID),
                "authoritative_Lancestor_tables": str(LANC),
            },
            "forbidden_steps": {"linker_rerun": False, "HitAnalysis": False, "reconstruction": False, "simulation": False, "generator": False, "association_definition_changed": False},
        }
        (staging / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        (staging / "validation/validation_summary.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
        os.replace(staging, FINAL)
        make_report(FINAL, matrices, denominators, validation, len(difference_specs))
        print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    except Exception:
        print(f"FAILED staging retained for audit: {staging}", file=sys.stderr, flush=True)
        raise


def rerender_contact_sheets() -> None:
    if not FINAL.is_dir():
        raise FileNotFoundError(FINAL)
    rows = read_csv(FINAL / "tables/migration_matrices_all.csv")
    matrices = {}
    for selection in SELECTIONS:
        order = outcomes(selection)
        for sample in SAMPLES:
            for definition in DEFINITIONS:
                values = np.zeros((len(TRUTH_ORDER), len(order)), dtype=float)
                selected = [row for row in rows if row["sample"] == sample and row["truth_definition"] == definition and row["selection"] == selection]
                for row in selected:
                    values[TRUTH_ORDER.index(row["truth_category"]), order.index(row["reco_category"])] = float(row["percent"])
                matrices[(sample, definition, selection)] = values
    for selection in SELECTIONS:
        stem = FINAL / "contact_sheets" / f"migration_matrix_grid_W_P8_G_Ldirect_Lancestor_{selection}"
        plot_contact(stem, matrices, selection)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--rerender-contact-sheets", action="store_true")
    cli = parser.parse_args()
    configure(cli.config)
    if cli.rerender_contact_sheets:
        rerender_contact_sheets()
    else:
        main()
