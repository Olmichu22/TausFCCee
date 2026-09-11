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
from modules.fcc_workflow_interface import truthlink_representative

COMPARISON_SCHEMA = "fcc_mc_comparison_v1"
# Legacy four-panel products remain unchanged; new machine-readable performance
# products use the complete explicit nominal inventory below.
TRUTH_SPECIES = ("electron", "muon", "photon", "charged_pion")
NOMINAL_TRUTH_SPECIES = (
    "electron", "muon", "photon", "charged_pion", "charged_kaon", "K0L",
)
PID_CATEGORIES = ("electron", "muon", "photon", "charged_pion", "K0S", "neutron", "Lambda")
ASSOCIATIONS = ("G", "Ldirect", "Lancestor")
OUTCOMES = ("associated_unique", "association_unmatched", "ambiguous_multiple_pfo")
PDG_TO_TRUTH = {
    11: "electron", 13: "muon", 22: "photon", 211: "charged_pion",
    321: "charged_kaon", 130: "K0L",
}
FIDUCIAL_FIELDS = ("truth_p_min", "truth_theta_min_deg", "truth_theta_max_deg")
DEFAULT_EFFICIENCY_BINS = {
    "p": (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0),
    "theta": tuple(float(value) for value in range(0, 181, 10)),
    "phi": tuple(float(-math.pi + index * math.pi / 6) for index in range(13)),
}


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
    configured = dict(payload.get("performance", {}))
    override = dict(comparison.get("performance", {}))
    configured.update({key: value for key, value in override.items() if key != "efficiency_bins"})
    configured["efficiency_bins"] = {
        **configured.get("efficiency_bins", {}), **override.get("efficiency_bins", {}),
    }
    return {"name": name, "contract": payload["scientific_contract"], **comparison,
            "performance": normalize_performance_config(configured)}


def kine_values(energy: float, px: float, py: float, pz: float) -> dict:
    px, py, pz = float(px), float(py), float(pz)
    p = math.sqrt(px ** 2 + py ** 2 + pz ** 2)
    pt = math.hypot(px, py)
    return {
        "energy": float(energy), "p": p, "pt": pt,
        "theta": math.degrees(math.atan2(pt, pz)) if p else 0.0,
        "phi": math.atan2(py, px),
        "costheta": pz / p if p else 1.0,
        "px": px, "py": py, "pz": pz,
    }


def kine_particle(particle) -> dict:
    momentum = particle.getMomentum()
    return kine_values(particle.getEnergy(), momentum.x, momentum.y, momentum.z)


def object_index(obj) -> int:
    return int(obj.getObjectID().index)


def nominal_truth_species(pdg: int) -> str | None:
    """Return the explicit nominal performance category, never an inferred one."""
    pdg = int(pdg)
    if abs(pdg) in (11, 13, 211, 321):
        return PDG_TO_TRUTH[abs(pdg)]
    if pdg in (22, 130):
        return PDG_TO_TRUTH[pdg]
    return None


def truth_species(pdg: int, charge: float) -> str:
    known = nominal_truth_species(pdg)
    if known:
        return known
    if abs(int(pdg)) > 100 and float(charge) != 0:
        return "other_charged_hadron"
    if abs(int(pdg)) > 100:
        return "neutral_hadron"
    return "other"


def wrap_delta_phi(phi_reco: float, phi_truth: float) -> float:
    """Return phi_reco-phi_truth in the canonical interval (-pi, pi]."""
    wrapped = (float(phi_reco) - float(phi_truth) + math.pi) % (2 * math.pi) - math.pi
    return math.pi if wrapped <= -math.pi else wrapped


def clamped_acos_mrad(cosine: float) -> float:
    return math.acos(max(-1.0, min(1.0, float(cosine)))) * 1000.0


def angle3d_mrad(truth_xyz: Iterable[float], reco_xyz: Iterable[float]) -> float | None:
    truth = tuple(float(value) for value in truth_xyz)
    reco = tuple(float(value) for value in reco_xyz)
    if len(truth) != 3 or len(reco) != 3:
        raise ValueError("3D opening angle requires exactly three components per vector")
    truth_norm = math.sqrt(sum(value * value for value in truth))
    reco_norm = math.sqrt(sum(value * value for value in reco))
    if truth_norm == 0.0 or reco_norm == 0.0:
        return None
    cosine = sum(left * right for left, right in zip(truth, reco)) / (truth_norm * reco_norm)
    return clamped_acos_mrad(cosine)


def relative_residual(reco: float, truth: float) -> float | None:
    return None if float(truth) == 0.0 else (float(reco) - float(truth)) / float(truth)


def performance_residuals(pair: dict) -> dict:
    """Additive residual columns for one maintained truth/representative pair."""
    result = {
        "dp_over_p": relative_residual(pair["reco_p"], pair["truth_p"]),
        "dtheta_mrad": (float(pair["reco_theta"]) - float(pair["truth_theta"])) * math.pi / 180 * 1000,
        "de_over_e": relative_residual(pair["reco_energy"], pair["truth_energy"]),
        "dphi_mrad": None,
        "angle3d_mrad": None,
    }
    if "reco_phi" in pair and "truth_phi" in pair:
        result["dphi_mrad"] = wrap_delta_phi(pair["reco_phi"], pair["truth_phi"]) * 1000
    vector_keys = tuple(f"{side}_{axis}" for side in ("truth", "reco") for axis in ("px", "py", "pz"))
    if all(key in pair for key in vector_keys):
        result["angle3d_mrad"] = angle3d_mrad(
            (pair["truth_px"], pair["truth_py"], pair["truth_pz"]),
            (pair["reco_px"], pair["reco_py"], pair["reco_pz"]),
        )
    return result


def selected_truth_inventory(rows: Iterable[dict]) -> list[dict]:
    """Count every selected-truth row as nominal species or explicit other PDG."""
    nominal = Counter()
    other = Counter()
    for row in rows:
        pdg = int(row["truth_pdg"])
        species = nominal_truth_species(pdg)
        if species is None:
            other[pdg] += 1
        else:
            nominal[species] += 1
    result = [
        {"category": "nominal_species", "species": species, "pdg": "", "count": nominal[species]}
        for species in NOMINAL_TRUTH_SPECIES
    ]
    result.extend(
        {"category": "other_selected_truth", "species": "other_selected_truth", "pdg": pdg, "count": count}
        for pdg, count in sorted(other.items())
    )
    if sum(int(row["count"]) for row in result) != sum(nominal.values()) + sum(other.values()):
        raise AssertionError("selected-truth inventory accounting failed")
    return result


def normalize_performance_config(config: dict | None = None) -> dict:
    """Return validated optional truth-fiducial and differential-bin settings."""
    configured = dict(config or {})
    result = {field: configured.get(field) for field in FIDUCIAL_FIELDS}
    for field in FIDUCIAL_FIELDS:
        if result[field] is not None:
            result[field] = float(result[field])
    if result["truth_p_min"] is not None and result["truth_p_min"] < 0:
        raise ValueError("truth_p_min must be non-negative or null")
    low, high = result["truth_theta_min_deg"], result["truth_theta_max_deg"]
    if low is not None and high is not None and low > high:
        raise ValueError("truth theta minimum exceeds maximum")
    result["association_method"] = str(configured.get("association_method", "Lancestor"))
    if result["association_method"] not in ASSOCIATIONS:
        raise ValueError("performance association_method is not maintained")
    supplied_bins = configured.get("efficiency_bins", {})
    result["efficiency_bins"] = {}
    for variable, defaults in DEFAULT_EFFICIENCY_BINS.items():
        edges = tuple(float(value) for value in supplied_bins.get(variable, defaults))
        if len(edges) < 2 or any(not left < right for left, right in zip(edges[:-1], edges[1:])):
            raise ValueError(f"{variable} efficiency bins must be strictly increasing")
        result["efficiency_bins"][variable] = edges
    p_edges, theta_edges, phi_edges = (result["efficiency_bins"][name]
                                        for name in ("p", "theta", "phi"))
    if p_edges[0] > 0:
        raise ValueError("truth-p bins must start at or below zero")
    if theta_edges[0] > 0 or theta_edges[-1] < 180:
        raise ValueError("truth-theta bins must cover [0,180] degrees")
    if not (math.isclose(phi_edges[0], -math.pi) and math.isclose(phi_edges[-1], math.pi)):
        raise ValueError("truth-phi bins must span exactly [-pi,+pi]")
    return result


def _truth_fiducial_accepts_normalized(row: dict, configured: dict) -> bool:
    if configured["truth_p_min"] is not None:
        truth_p = float(row["p"] if "p" in row else row["truth_p"])
        if truth_p < configured["truth_p_min"]:
            return False
    if configured["truth_theta_min_deg"] is not None or configured["truth_theta_max_deg"] is not None:
        truth_theta = float(row["theta"] if "theta" in row else row["truth_theta"])
        if (configured["truth_theta_min_deg"] is not None
                and truth_theta < configured["truth_theta_min_deg"]):
            return False
        if (configured["truth_theta_max_deg"] is not None
                and truth_theta > configured["truth_theta_max_deg"]):
            return False
    return True


def truth_fiducial_accepts(row: dict, config: dict | None = None) -> bool:
    """Apply optional inclusive cuts to truth p/theta only."""
    return _truth_fiducial_accepts_normalized(row, normalize_performance_config(config))


def fiducial_outcomes(rows: Iterable[dict], config: dict | None = None) -> list[dict]:
    """Filter already-selected truth outcome rows without consulting reco values."""
    configured = normalize_performance_config(config)
    return [row for row in rows if _truth_fiducial_accepts_normalized(row, configured)]


def fiducial_provenance(config: dict | None = None) -> dict:
    configured = normalize_performance_config(config)
    return {field: "none" if configured[field] is None else configured[field]
            for field in FIDUCIAL_FIELDS}


def canonical_truth_phi(phi: float) -> float:
    """Return one truth phi in the same canonical interval (-pi,pi] as dphi."""
    return wrap_delta_phi(float(phi), 0.0)


def binned_efficiency(rows: Iterable[dict], variable: str,
                      bins: Iterable[float]) -> list[dict]:
    """Count maintained associated-unique successes in truth-variable bins."""
    population = list(rows)
    edges = np.asarray(tuple(float(value) for value in bins), dtype=float)
    if len(edges) < 2 or np.any(edges[1:] <= edges[:-1]):
        raise ValueError("efficiency bins must be strictly increasing")
    values = np.asarray([
        canonical_truth_phi(row["phi"]) if variable == "phi" else float(row[variable])
        for row in population
    ], dtype=float)
    successes = np.asarray([row["outcome"] == "associated_unique" for row in population], dtype=bool)
    denominator, _ = np.histogram(values, bins=edges)
    numerator, _ = np.histogram(values[successes], bins=edges)
    under = values < edges[0]
    over = values > edges[-1]
    records = [{"bin_kind": "underflow", "bin_low": -math.inf, "bin_high": float(edges[0]),
                "N_truth": int(under.sum()), "N_associated_unique": int((under & successes).sum())}]
    records.extend(
        {"bin_kind": "regular", "bin_low": float(low), "bin_high": float(high),
         "N_truth": int(total), "N_associated_unique": int(passed)}
        for low, high, total, passed in zip(edges[:-1], edges[1:], denominator, numerator)
    )
    records.append({"bin_kind": "overflow", "bin_low": float(edges[-1]), "bin_high": math.inf,
                    "N_truth": int(over.sum()), "N_associated_unique": int((over & successes).sum())})
    if sum(row["N_truth"] for row in records) != len(population):
        raise AssertionError(f"truth-{variable} binning does not account for the denominator")
    if (any(row["N_associated_unique"] > row["N_truth"] for row in records)
            or sum(row["N_associated_unique"] for row in records) != int(successes.sum())):
        raise AssertionError("efficiency numerator is not a subset of the denominator")
    for row in records:
        total, passed = row["N_truth"], row["N_associated_unique"]
        row.update({"variable": variable,
                    "efficiency_percent": "" if total == 0 else 100 * passed / total})
    return records


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
        for key in ("energy", "p", "pt", "theta", "phi", "costheta", "px", "py", "pz"):
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
        row.update(performance_residuals(row))
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

    It only loads the configured maintained association methods. G is computed
    only when explicitly requested; L_direct/L_ancestor are always read from
    existing assignment products and are never rebuilt here.
    """
    import podio.root_io as root_io

    internal, label = spec["internal_name"], spec["presentation_label"]
    methods = tuple(spec.get("association_methods", ASSOCIATIONS))
    if not methods or any(method not in ASSOCIATIONS for method in methods):
        raise ValueError("association_methods contains an unsupported method")
    collect_topology = bool(spec.get("collect_tau_topology", True))
    if "G" in methods:
        from modules.NeutralRecover import get_reco_mc_links_by_dR
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
    coverage = {method: Counter() for method in methods}
    n_events = 0
    reader = root_io.Reader(str(rec))
    for event_index, event in enumerate(reader.get("events")):
        n_events += 1
        mc = list(event.get("MCParticles"))
        pfos = list(event.get("PandoraPFOs")) if "G" in methods else []
        pdgs = [int(particle.getPDG()) for particle in mc]
        parents = [[object_index(parent) for parent in particle.getParents()] for particle in mc]
        daughters = [[object_index(child) for child in particle.getDaughters()] for particle in mc]
        selected = [index for index, particle in enumerate(mc) if selected_truth_particle(particle)]
        selected_set = set(selected)
        grows = {}
        if "G" in methods:
            gframe = get_reco_mc_links_by_dR(event, {}, {}, max_dR=0.1, dedup_mode="reco")
            grows = {int(row.gen): row for row in gframe.itertuples() if int(row.gen) >= 0}
            matched_g_pfos = {int(row.reco) for row in gframe.itertuples()
                              if int(row.gen) >= 0 and int(row.reco) >= 0}
            coverage["G"]["usable"] += len(matched_g_pfos)
            coverage["G"]["all"] += len(pfos)

        for row in (direct_by_event.get(event_index, []) if "Ldirect" in methods else ()):
            coverage["Ldirect"]["all"] += 1
            status, assigned = row.get("truthlink_status"), row.get("assigned_mc_index")
            if status == "assigned" and assigned is not None and int(assigned) in selected_set:
                coverage["Ldirect"]["usable"] += 1
            elif "ambiguous" in str(status):
                coverage["Ldirect"]["ambiguous"] += 1
            elif status == "assigned":
                coverage["Ldirect"]["non_analysis"] += 1
        for promoted in (ancestor_by_event.get(event_index, []) if "Lancestor" in methods else ()):
            coverage["Lancestor"]["all"] += 1
            coverage["Lancestor"]["usable"] += int(promoted is not None)

        if collect_topology:
            terminal_indices = [index for index, pdg in enumerate(pdgs)
                                if abs(pdg) == 15
                                and not any(abs(pdgs[child]) == 15 for child in daughters[index])]
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
            if species not in NOMINAL_TRUTH_SPECIES:
                continue
            grow = grows.get(index)
            g_pfo = None if grow is None or int(grow.reco) < 0 else int(grow.reco)
            method_rows = {}
            if "G" in methods:
                method_rows["G"] = (("association_unmatched", None) if g_pfo is None
                                    else ("associated_unique", g_pfo))
            if "Ldirect" in methods:
                method_rows["Ldirect"] = strict_outcome(direct.get((event_index, index), []))
            if "Lancestor" in methods:
                method_rows["Lancestor"] = strict_outcome(ancestor.get((event_index, index), []))
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
                    reco = {key: item[key] for key in (
                        "energy", "p", "pt", "theta", "phi", "costheta", "px", "py", "pz",
                    )}
                    category = item["reco_category"]
                pair = {**base, "association_method": method,
                        "representative_pfo_index": pfo_index, "reco_category": category,
                        **{f"truth_{key}": kine[key] for key in (
                            "energy", "p", "pt", "theta", "phi", "costheta", "px", "py", "pz",
                        )},
                        **{f"reco_{key}": reco[key] for key in (
                            "energy", "p", "pt", "theta", "phi", "costheta", "px", "py", "pz",
                        )}}
                pair.update(performance_residuals(pair))
                pairs.append(pair)
    if n_events != expected:
        raise AssertionError(f"{internal}: expected {expected} events, found {n_events}")
    coverage_rows = []
    for method in methods:
        count = coverage[method]
        all_pfos = count["all"]
        coverage_rows.append({"sample": label, "association_method": method,
            "N_all_PFO": all_pfos, "N_usable_truth_link": count["usable"],
            "N_not_usable": all_pfos - count["usable"],
            "coverage_percent": "" if not all_pfos else 100 * count["usable"] / all_pfos,
            "N_ambiguous": count["ambiguous"], "N_linked_to_non_analysis_MC": count["non_analysis"],
            "definition": "unique PFO link to selected stable analysis-level truth particle",
            "source": "maintained G rows / Ldirect assignments / Lancestor assignments"})
    return SampleData(internal, label, expected, truth, outcomes, pairs, terminal, decay,
                      truth_counts, coverage_rows, {"adapter": "workflow_products",
                                      "source_file_id": str(spec.get("source_file_id", "")), "rec": str(rec),
                                      "direct": str(direct_path), "ancestor": str(ancestor_path),
                                      "association_methods": list(methods),
                                      "collect_tau_topology": collect_topology})


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
    for method in combined:
        count = combined[method]
        all_pfos = count["N_all_PFO"]
        rows.append({"sample": label, "association_method": method, **count,
            "coverage_percent": "" if not all_pfos else 100 * count["N_usable_truth_link"] / all_pfos,
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
