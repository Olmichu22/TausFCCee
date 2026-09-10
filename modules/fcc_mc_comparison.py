"""Frozen 2026-09-01 W-versus-MC comparison primitives.

The functions in this module are a parameterized port of the executable
builders preserved in ``archive/ild-tau_pre_consolidation_active_20260902``:
``audit_tau_topology_equalstats.py``, ``audit_pfo_reconstruction_efficiency_equalstats.py``
and ``audit_pid_equalstats.py``.  Sample names and paths are data; the selection,
ancestry, association-outcome, representative-PFO and PID rules are not.
"""
from __future__ import annotations

from collections import Counter, defaultdict, deque
import csv
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pyarrow.parquet as pq
import yaml

from modules.fcc_truth_definitions import (
    has_tau_origin,
    reconstructed_pid_category,
    selected_truth_particle,
)

COMPARISON_SCHEMA = "fcc_mc_comparison_v1"
TRUTH_SPECIES = ("electron", "muon", "photon", "charged_pion")
PID_CATEGORIES = ("electron", "muon", "photon", "charged_pion", "K0S", "neutron", "Lambda")
ASSOCIATIONS = ("G", "Ldirect", "Lancestor")
OUTCOMES = ("associated_unique", "association_unmatched", "ambiguous_multiple_pfo")
PDG_TO_TRUTH = {11: "electron", 13: "muon", 22: "photon", 211: "charged_pion"}


def read_csv(path: Path) -> list[dict]:
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: Iterable[dict], fields: Iterable[str] | None = None) -> None:
    materialized = list(rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fields) if fields is not None else list(materialized[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(materialized)


def expand_path(value: str) -> Path:
    expanded = os.path.expandvars(str(value))
    if "$" in expanded:
        raise ValueError(f"unresolved environment variable: {value}")
    return Path(expanded)


def load_comparison(path: Path, name: str) -> dict:
    payload = yaml.safe_load(Path(path).read_text())
    if payload.get("schema_version") != COMPARISON_SCHEMA:
        raise ValueError("unsupported comparison configuration")
    try:
        comparison = payload["comparisons"][name]
    except KeyError as error:
        raise ValueError(f"unknown comparison: {name}") from error
    if len(comparison.get("samples", [])) != 2:
        raise ValueError("a comparison must contain exactly two samples")
    return {"name": name, "contract": payload["scientific_contract"], **comparison}


def kine_values(energy: float, px: float, py: float, pz: float) -> dict:
    p = math.sqrt(float(px) ** 2 + float(py) ** 2 + float(pz) ** 2)
    pt = math.hypot(float(px), float(py))
    return {
        "energy": float(energy), "p": p, "pt": pt,
        "theta": math.degrees(math.atan2(pt, float(pz))) if p else 0.0,
        "costheta": float(pz) / p if p else 1.0,
    }


def kine_particle(particle) -> dict:
    momentum = particle.getMomentum()
    return kine_values(particle.getEnergy(), momentum.x, momentum.y, momentum.z)


def object_index(obj) -> int:
    return int(obj.getObjectID().index)


def truth_species(pdg: int, charge: float) -> str:
    known = PDG_TO_TRUTH.get(abs(int(pdg)))
    if known:
        return known
    if abs(int(pdg)) > 100 and float(charge) != 0:
        return "other_charged_hadron"
    if abs(int(pdg)) > 100:
        return "neutral_hadron"
    return "other"


def truthlink_representative(rows: list[dict]) -> dict:
    """Exact port of the frozen W/P8O representative-PFO inversion."""
    if not rows:
        return {"status": "unmatched", "pfo_index": None, "multiplicity": 0}
    if len(rows) == 1:
        return {"status": "assigned", "pfo_index": int(rows[0]["pfo_index"]), "multiplicity": 1}
    track = [row for row in rows if int(row["track_permille"]) > 0]
    if track:
        best_t = max(int(row["track_permille"]) for row in track)
        stage = [row for row in track if int(row["track_permille"]) == best_t]
        best_c = max(int(row["cluster_permille"]) for row in stage)
        winners = [row for row in stage if int(row["cluster_permille"]) == best_c]
    else:
        best_c = max(int(row["cluster_permille"]) for row in rows)
        winners = [row for row in rows if int(row["cluster_permille"]) == best_c]
    if len(winners) != 1:
        return {"status": "ambiguous_multiple_pfo", "pfo_index": None,
                "multiplicity": len(rows), "tied_pfo_count": len(winners)}
    return {"status": "assigned", "pfo_index": int(winners[0]["pfo_index"]),
            "multiplicity": len(rows)}


def strict_outcome(rows: list[dict]) -> tuple[str, int | None]:
    representative = truthlink_representative(rows)
    if representative["status"] == "unmatched":
        return "association_unmatched", None
    if representative["status"] == "ambiguous_multiple_pfo":
        return "ambiguous_multiple_pfo", None
    return "associated_unique", int(representative["pfo_index"])


def _visible_signature(mc, daughters: list[list[int]], tau_index: int) -> tuple[list[int], int]:
    leaves, neutrinos, seen = [], [], set()
    pending = list(daughters[tau_index])
    while pending:
        index = pending.pop()
        if index in seen:
            continue
        seen.add(index)
        particle = mc[index]
        pdg = int(particle.getPDG())
        if abs(pdg) in {12, 14, 16}:
            neutrinos.append(index)
        elif pdg == 111:
            leaves.append(index)
        elif int(particle.getGeneratorStatus()) == 1 or not daughters[index]:
            leaves.append(index)
        else:
            pending.extend(daughters[index])
    return leaves, len(neutrinos)


def classify_terminal_tau(mc, daughters: list[list[int]], tau_index: int) -> str:
    """Exact frozen visible-terminal classification (explicit pi0 is one object)."""
    leaves, _ = _visible_signature(mc, daughters, tau_index)
    pdgs = [int(mc[index].getPDG()) for index in leaves]
    if any(abs(pdg) == 11 for pdg in pdgs):
        return "electron_leptonic"
    if any(abs(pdg) == 13 for pdg in pdgs):
        return "muon_leptonic"
    prongs = sum(abs(float(mc[index].getCharge())) > 1e-12
                 for index in leaves if int(mc[index].getPDG()) != 111)
    n_pi0 = sum(int(mc[index].getPDG()) == 111 for index in leaves)
    if prongs == 1:
        return "hadronic_1prong_with_pi0" if n_pi0 else "hadronic_1prong_no_pi0"
    if prongs == 3:
        return "hadronic_3prong_with_pi0" if n_pi0 else "hadronic_3prong_no_pi0"
    if prongs > 0:
        return "other_hadronic"
    return "other_or_unresolved"


@dataclass
class SampleData:
    internal_name: str
    presentation_label: str
    expected_events: int
    truth: list[dict]
    outcomes: list[dict]
    pairs: list[dict]
    terminal_taus: list[dict]
    decay_modes: Counter
    truth_counts: Counter
    pfo_coverage: list[dict]
    provenance: dict


def _normalize_frozen_rows(rows: list[dict], sample: str) -> list[dict]:
    normalized = []
    for row in rows:
        if row["sample"] != sample:
            continue
        item = dict(row)
        for key in ("event_in_file", "truth_index", "truth_pdg"):
            if key in item:
                item[key] = int(item[key])
        for key in ("energy", "p", "pt", "theta", "costheta"):
            if key in item:
                item[key] = float(item[key])
            for prefix in ("truth_", "reco_"):
                compound = prefix + key
                if compound in item:
                    item[compound] = float(item[compound])
        if "tau_ancestor" in item:
            item["tau_ancestor"] = str(item["tau_ancestor"]).lower() == "true"
        normalized.append(item)
    return normalized


def load_frozen_sample(spec: dict) -> SampleData:
    """Load authoritative frozen intermediate products without recomputation."""
    internal = spec["internal_name"]
    expected = int(spec["expected_events"])
    sources = spec["frozen_products"]
    outcome_rows = _normalize_frozen_rows(
        pq.read_table(expand_path(sources["selected_truth_outcomes"])).to_pylist(), internal)
    truth = [row for row in outcome_rows if row["truth_definition"] == "Lancestor"]
    pair_rows = _normalize_frozen_rows(
        pq.read_table(expand_path(sources["lancestor_pairs"])).to_pylist(), internal)
    truth_lookup = {(row["source_file"], row["event_in_file"], row["truth_index"]): row for row in truth}
    pairs = []
    for row in pair_rows:
        key = (row["source_file"], row["event_in_file"], row["truth_index"])
        anchor = truth_lookup[key]
        row["tau_ancestor"] = anchor["tau_ancestor"]
        row["association_method"] = "Lancestor"
        pairs.append(row)
    terminal = []
    for row in read_csv(expand_path(sources["tau_records"])):
        if row["sample"] != internal or str(row["terminal_decay_tau"]).lower() != "true":
            continue
        terminal.append({"p": float(row["p"]), "pt": float(row["pt"]),
                         "theta": float(row["theta_deg"])})
    decay = Counter()
    for row in read_csv(expand_path(sources["tau_decay_channels"])):
        if row["sample"] == internal:
            decay[row["category"]] = int(row["count"])
    truth_counts = Counter()
    for row in read_csv(expand_path(sources["truth_composition"])):
        if row["sample"] == internal:
            truth_counts[row["origin"], row["category"]] = int(row["count"])
    coverage = []
    if sources.get("coverage"):
        for row in read_csv(expand_path(sources["coverage"])):
            if row["sample"] == spec["presentation_label"]:
                coverage.append(dict(row))
    provenance = {"adapter": "frozen_products", "inputs": sources}
    if sources.get("presentation_pid_root"):
        root = expand_path(sources["presentation_pid_root"])
        token = "W" if internal == "W" else spec["presentation_label"]
        for method in ("G", "Ldirect"):
            for row in read_csv(root / f"{token}_{method}_conditional_pid_matrix.csv"):
                for _ in range(int(row["count"])):
                    pairs.append({"sample": internal, "association_method": method,
                                  "truth_species": row["truth_species"], "reco_category": row["reco_pid"],
                                  "tau_ancestor": True, "aggregate_pid_only": True})
        non_tau_token = "WHIZARD" if internal == "W" else spec["presentation_label"]
        for row in read_csv(root / f"{non_tau_token}_non_tau_photon_pid_G_Ldirect_Lancestor.csv"):
            if row["association_method"] == "Lancestor":
                continue
            for _ in range(int(row["count"])):
                pairs.append({"sample": internal, "association_method": row["association_method"],
                              "truth_species": "photon", "reco_category": row["reco_pid"],
                              "tau_ancestor": False, "aggregate_pid_only": True})
    return SampleData(internal, spec["presentation_label"], expected, truth, outcome_rows,
                      pairs, terminal, decay, truth_counts, coverage, provenance)


def _assignment_maps(direct_path: Path, ancestor_path: Path):
    direct_rows = pq.read_table(direct_path).to_pylist()
    ancestor_rows = pq.read_table(ancestor_path).to_pylist()
    ancestor_by_pfo = {(int(row["event_in_file"]), int(row["pfo_index"])): row.get("ancestor_mc_index")
                       for row in ancestor_rows}
    if len(ancestor_by_pfo) != len(ancestor_rows) or len(direct_rows) != len(ancestor_rows):
        raise AssertionError("direct/ancestor PFO-row join is not one-to-one")
    pfo_map, direct, ancestor, by_event, ancestor_by_event = {}, defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
    for row in direct_rows:
        event, pfo = int(row["event_in_file"]), int(row["pfo_index"])
        values = kine_values(row["pfo_energy"], row["pfo_px"], row["pfo_py"], row["pfo_pz"])
        item = {"pfo_index": pfo, "pfo_type": int(row["pfo_type"]),
                "reco_category": reconstructed_pid_category(int(row["pfo_type"])),
                "track_permille": int(row.get("track_permille") or 0),
                "cluster_permille": int(row.get("cluster_permille") or 0), **values}
        pfo_map[event, pfo] = item
        by_event[event].append(row)
        if row.get("truthlink_status") == "assigned" and row.get("assigned_mc_index") is not None:
            direct[event, int(row["assigned_mc_index"])].append(item)
        promoted = ancestor_by_pfo[event, pfo]
        ancestor_by_event[event].append(promoted)
        if promoted is not None:
            ancestor[event, int(promoted)].append(item)
    return pfo_map, direct, ancestor, by_event, ancestor_by_event


def extract_workflow_sample(spec: dict) -> SampleData:
    """Build normalized rows from a final REC and frozen workflow assignments.

    This calls the maintained HitAnalysis G implementation and only inverts
    already-produced L_direct/L_ancestor PFO rows.  It never rebuilds either L.
    """
    import podio.root_io as root_io
    from modules.NeutralRecover import get_reco_mc_links_by_dR

    internal, label = spec["internal_name"], spec["presentation_label"]
    expected = int(spec["expected_events"])
    rec = expand_path(spec["rec"])
    direct_path, ancestor_path = expand_path(spec["direct"]), expand_path(spec["ancestor"])
    for path in (rec, direct_path, ancestor_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    receipt = spec.get("preflight_receipt")
    if receipt:
        payload = json.loads(expand_path(receipt).read_text())
        if payload.get("status") != "PASS" or int(payload.get("events", -1)) != expected:
            raise RuntimeError("TruthlinkV1 preflight receipt does not cover the configured scope")
        if Path(payload["input"]).resolve() != rec.resolve():
            raise RuntimeError("preflight receipt/input mismatch")

    pfo_map, direct, ancestor, direct_by_event, ancestor_by_event = _assignment_maps(direct_path, ancestor_path)
    truth, outcomes, pairs, terminal, decay = [], [], [], [], Counter()
    truth_counts = Counter()
    coverage = {method: Counter() for method in ASSOCIATIONS}
    n_events = 0
    reader = root_io.Reader(str(rec))
    for event_index, event in enumerate(reader.get("events")):
        n_events += 1
        mc = list(event.get("MCParticles")); pfos = list(event.get("PandoraPFOs"))
        pdgs = [int(particle.getPDG()) for particle in mc]
        parents = [[object_index(parent) for parent in particle.getParents()] for particle in mc]
        daughters = [[object_index(child) for child in particle.getDaughters()] for particle in mc]
        selected = [index for index, particle in enumerate(mc) if selected_truth_particle(particle)]
        selected_set = set(selected)
        gframe = get_reco_mc_links_by_dR(event, {}, {}, max_dR=0.1, dedup_mode="reco")
        grows = {int(row.gen): row for row in gframe.itertuples() if int(row.gen) >= 0}
        matched_g_pfos = {int(row.reco) for row in gframe.itertuples() if int(row.gen) >= 0 and int(row.reco) >= 0}
        coverage["G"]["usable"] += len(matched_g_pfos)
        coverage["G"]["all"] += len(pfos)

        for row in direct_by_event.get(event_index, []):
            coverage["Ldirect"]["all"] += 1
            status, assigned = row.get("truthlink_status"), row.get("assigned_mc_index")
            if status == "assigned" and assigned is not None and int(assigned) in selected_set:
                coverage["Ldirect"]["usable"] += 1
            elif "ambiguous" in str(status):
                coverage["Ldirect"]["ambiguous"] += 1
            elif status == "assigned":
                coverage["Ldirect"]["non_analysis"] += 1
        for promoted in ancestor_by_event.get(event_index, []):
            coverage["Lancestor"]["all"] += 1
            coverage["Lancestor"]["usable"] += int(promoted is not None)

        terminal_indices = [index for index, pdg in enumerate(pdgs)
                            if abs(pdg) == 15 and not any(abs(pdgs[child]) == 15 for child in daughters[index])]
        for index in terminal_indices:
            terminal.append(kine_particle(mc[index]))
            decay[classify_terminal_tau(mc, daughters, index)] += 1

        event_pfos = {object_index(pfo): pfo for pfo in pfos}
        for index in selected:
            particle = mc[index]
            species = truth_species(pdgs[index], particle.getCharge())
            tau = has_tau_origin(pdgs, parents, index)
            kine = kine_particle(particle)
            base = {"sample": internal, "source_file": rec.name,
                    "source_file_id": str(spec.get("source_file_id", "")), "event_in_file": event_index,
                    "truth_index": index, "truth_species": species, "truth_pdg": pdgs[index],
                    "tau_ancestor": tau, "parentless": not parents[index], **kine}
            truth.append(base)
            truth_counts[("tau_origin" if tau else "non_tau_origin", species)] += 1
            if species not in TRUTH_SPECIES:
                continue
            grow = grows.get(index)
            g_pfo = None if grow is None or int(grow.reco) < 0 else int(grow.reco)
            method_rows = {
                "G": ("association_unmatched", None) if g_pfo is None else ("associated_unique", g_pfo),
                "Ldirect": strict_outcome(direct.get((event_index, index), [])),
                "Lancestor": strict_outcome(ancestor.get((event_index, index), [])),
            }
            for method, (outcome, pfo_index) in method_rows.items():
                outcomes.append({**base, "truth_definition": method, "outcome": outcome})
                if outcome != "associated_unique":
                    continue
                if method == "G":
                    pfo = event_pfos[pfo_index]
                    reco = kine_particle(pfo)
                    category = reconstructed_pid_category(int(pfo.getPDG()))
                else:
                    item = pfo_map[event_index, pfo_index]
                    reco = {key: item[key] for key in ("energy", "p", "pt", "theta", "costheta")}
                    category = item["reco_category"]
                pairs.append({**base, "association_method": method,
                              "representative_pfo_index": pfo_index, "reco_category": category,
                              **{f"truth_{key}": kine[key] for key in ("energy", "p", "pt", "theta", "costheta")},
                              **{f"reco_{key}": reco[key] for key in ("energy", "p", "pt", "theta", "costheta")}})
    if n_events != expected:
        raise AssertionError(f"{internal}: expected {expected} events, found {n_events}")
    coverage_rows = []
    for method in ASSOCIATIONS:
        count = coverage[method]
        coverage_rows.append({"sample": label, "association_method": method,
            "N_all_PFO": count["all"], "N_usable_truth_link": count["usable"],
            "N_not_usable": count["all"] - count["usable"],
            "coverage_percent": 100 * count["usable"] / count["all"],
            "N_ambiguous": count["ambiguous"], "N_linked_to_non_analysis_MC": count["non_analysis"],
            "definition": "unique PFO link to selected stable analysis-level truth particle",
            "source": "maintained G rows / Ldirect assignments / Lancestor assignments"})
    return SampleData(internal, label, expected, truth, outcomes, pairs, terminal, decay,
                      truth_counts, coverage_rows, {"adapter": "workflow_products",
                                      "source_file_id": str(spec.get("source_file_id", "")), "rec": str(rec),
                                      "direct": str(direct_path), "ancestor": str(ancestor_path)})


def load_sample(spec: dict) -> SampleData:
    adapter = spec.get("adapter")
    if adapter == "frozen_products":
        return load_frozen_sample(spec)
    if adapter == "workflow_products":
        return extract_workflow_sample(spec)
    if adapter == "workflow_manifest":
        manifest = read_csv(expand_path(spec["manifest"]))
        columns = spec["columns"]
        include = set(map(str, spec.get("include_source_ids", [])))
        parts = []
        for row in manifest:
            source_id = str(row[columns["source_id"]])
            if include and source_id not in include:
                continue
            part = {**spec, "adapter": "workflow_products", "expected_events": int(spec["events_per_file"]),
                    "source_file_id": source_id,
                    "rec": row[columns["rec"]], "direct": row[columns["direct"]],
                    "ancestor": row[columns["ancestor"]]}
            part.pop("manifest", None); part.pop("columns", None); part.pop("include_source_ids", None)
            parts.append(extract_workflow_sample(part))
        if include and len(parts) != len(include):
            raise AssertionError("workflow manifest does not contain every configured source")
        expected = int(spec["expected_events"])
        if sum(part.expected_events for part in parts) != expected:
            raise AssertionError("workflow manifest event scope mismatch")
        decay = Counter(); truth_counts = Counter()
        for part in parts:
            decay.update(part.decay_modes); truth_counts.update(part.truth_counts)
        return SampleData(spec["internal_name"], spec["presentation_label"], expected,
            [row for part in parts for row in part.truth],
            [row for part in parts for row in part.outcomes],
            [row for part in parts for row in part.pairs],
            [row for part in parts for row in part.terminal_taus],
            decay, truth_counts,
            _sum_coverage(parts, spec["presentation_label"]),
            {"adapter": "workflow_manifest", "manifest": spec["manifest"], "parts": len(parts),
             "source_file_ids": sorted(include),
             "resolved_parts": [part.provenance for part in parts]})
    raise ValueError(f"unsupported sample adapter: {adapter}")


def _sum_coverage(parts: list[SampleData], label: str) -> list[dict]:
    combined = defaultdict(Counter)
    for part in parts:
        for row in part.pfo_coverage:
            method = row["association_method"]
            for key in ("N_all_PFO", "N_usable_truth_link", "N_not_usable", "N_ambiguous",
                        "N_linked_to_non_analysis_MC"):
                combined[method][key] += int(row[key])
    rows = []
    for method in ASSOCIATIONS:
        count = combined[method]
        rows.append({"sample": label, "association_method": method, **count,
            "coverage_percent": 100 * count["N_usable_truth_link"] / count["N_all_PFO"],
            "definition": "unique PFO link to selected stable analysis-level truth particle",
            "source": "maintained G rows / Ldirect assignments / Lancestor assignments"})
    return rows


def hist_with_flow(values: Iterable[float], bins: np.ndarray, denominator: int) -> list[dict]:
    values = np.asarray(list(values), dtype=float)
    counts, _ = np.histogram(values, bins=bins)
    result = [{"bin_low": -math.inf, "bin_high": float(bins[0]),
               "count": int(np.count_nonzero(values < bins[0]))}]
    result.extend({"bin_low": float(low), "bin_high": float(high), "count": int(count)}
                  for low, high, count in zip(bins[:-1], bins[1:], counts))
    result.append({"bin_low": float(bins[-1]), "bin_high": math.inf,
                   "count": int(np.count_nonzero(values > bins[-1]))})
    for index, row in enumerate(result):
        row["fraction"] = row["count"] / denominator if denominator else 0.0
        row["visible"] = 0 < index < len(result) - 1
    if sum(row["count"] for row in result) != denominator:
        raise AssertionError("histogram flow accounting failed")
    return result


def quantile_summary(values: Iterable[float]) -> dict:
    array = np.asarray(list(values), dtype=float)
    q16, median, q84 = np.quantile(array, (.16, .5, .84))
    return {"N": len(array), "median": float(median), "q16": float(q16), "q84": float(q84),
            "central68_halfwidth": float((q84 - q16) / 2)}
