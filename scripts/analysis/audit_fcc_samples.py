#!/usr/bin/env python3
"""Read-only, configuration-driven truth/PFO bookkeeping audit.

This deliberately derives truth bookkeeping from the REC MC graph and consumes
the already-produced L_direct/L_ancestor tables.  It never rebuilds either
assignment.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from modules.fcc_sample_catalog import load_catalog, resolve_records  # noqa: E402
from modules.fcc_truth_definitions import has_tau_origin, selected_truth_particle  # noqa: E402

SPECIES = {11: "electron", 13: "muon", 22: "photon", 211: "charged_pion"}


def _index(obj) -> int:
    return int(obj.getObjectID().index)


def scan_rec(path: Path) -> dict:
    from podio import root_io

    counts = Counter()
    status = Counter()
    parentless_species = Counter()
    selected_species = Counter()
    selected_per_event = []
    reader = root_io.Reader(str(path))
    for event in reader.get("events"):
        event_selected = 0
        counts["events"] += 1
        mc = event.get("MCParticles")
        pfos = event.get("PandoraPFOs")
        counts["mc"] += len(mc)
        counts["pfos"] += len(pfos)
        by_index = {_index(p): p for p in mc}
        parents = {i: [_index(q) for q in p.getParents()] for i, p in by_index.items()}
        pdgs = {i: int(p.getPDG()) for i, p in by_index.items()}
        for p in mc:
            status[str(int(p.getGeneratorStatus()))] += 1
            if not selected_truth_particle(p):
                continue
            event_selected += 1
            counts["selected"] += 1
            pdg = int(p.getPDG())
            species = SPECIES.get(abs(pdg), "other")
            selected_species[species] += 1
            if not p.getParents():
                counts["selected_parentless"] += 1
                parentless_species[species] += 1
            if abs(pdg) == 22:
                counts["photons"] += 1
                ordered = sorted(by_index)
                if ordered != list(range(len(ordered))):
                    raise ValueError(f"non-contiguous MCParticle indices in {path}")
                tau = has_tau_origin([pdgs[i] for i in ordered], [parents[i] for i in ordered], _index(p))
                counts["photons_tau" if tau else "photons_non_tau"] += 1
                if not p.getParents():
                    counts["photons_parentless"] += 1
        selected_per_event.append(event_selected)
    return {
        "counts": dict(counts),
        "generator_status": dict(status),
        "parentless_species": dict(parentless_species),
        "selected_species": dict(selected_species),
        "selected_per_event": selected_per_event,
    }


def selected_path(record: dict, sample: str) -> Path | None:
    if record.get("selected_truth"):
        return Path(record["selected_truth"])
    ancestor = Path(record["ancestor"])
    name = ancestor.name.replace("_ancestor_assignment.parquet", "_selected_truth.parquet")
    candidate = ancestor.parents[1] / "selected_truth" / name
    return candidate if candidate.exists() else None


def audit_sample(name: str, config: dict, workers: int = 1):
    records = resolve_records(config)
    rec_paths = [Path(record["source_rec"]) for record in records]
    prerequisites = [Path(record[key]) for record in records for key in ("source_rec", "direct", "ancestor")]
    missing = [str(path) for path in prerequisites if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(f"{name}: missing/non-empty prerequisite(s): {missing}")
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            scanned_records = list(pool.map(scan_rec, rec_paths))
    else:
        scanned_records = [scan_rec(path) for path in rec_paths]
    raw = Counter()
    gen_status = Counter()
    parentless_species = Counter()
    selected_species = Counter()
    selected_per_event = []
    direct_status = Counter()
    ancestor_status = Counter()
    depths = Counter()
    provenance = []
    compact_selected = 0
    for record, scanned in zip(records, scanned_records):
        required = [Path(record[k]) for k in ("source_rec", "direct", "ancestor")]
        raw.update(scanned["counts"])
        gen_status.update(scanned["generator_status"])
        parentless_species.update(scanned["parentless_species"])
        selected_species.update(scanned["selected_species"])
        selected_per_event.extend(scanned["selected_per_event"])
        direct = pq.read_table(required[1])
        ancestor = pq.read_table(required[2])
        if len(direct) != scanned["counts"].get("pfos", 0) or len(ancestor) != len(direct):
            raise ValueError(f"{name}/{record['source_id']}: PFO/direct/ancestor row mismatch")
        direct_status.update("<null>" if v is None else str(v) for v in direct["truthlink_status"].to_pylist())
        ancestor_status.update("<null>" if v is None else str(v) for v in ancestor["ancestor_status"].to_pylist())
        depths.update("<null>" if v is None else str(int(v)) for v in ancestor["ancestor_depth"].to_pylist())
        sp = selected_path(record, name)
        if sp:
            compact_selected += pq.read_metadata(sp).num_rows
        provenance.append({
            "sample_internal_name": name,
            "source_file_id": record["source_id"],
            "source_rec": str(required[0]),
            "direct_assignment": str(required[1]),
            "ancestor_assignment": str(required[2]),
            "selected_truth": str(sp or ""),
            "source_rec_bytes": required[0].stat().st_size,
        })

    expected = int(config["expected_events"])
    if raw["events"] != expected:
        raise ValueError(f"{name}: expected {expected} events, found {raw['events']}")
    assigned = direct_status["assigned"]
    same = ancestor_status["same_direct_selected"]
    promoted = ancestor_status["promoted_unique_ancestor"]
    no_ancestor = ancestor_status["ancestor_no_selected_ancestor"] + ancestor_status["no_selected_ancestor"]
    valid = same + promoted
    unresolved = raw["pfos"] - valid - no_ancestor
    errors = sum(v for k, v in ancestor_status.items() if "error" in k or "cycle" in k)
    def frac(n, d):
        return n / d if d else None
    compact_scope = config.get("selected_truth_product_scope", "all_selected_truth")
    compact_expected = (selected_species["photon"] + selected_species["charged_pion"]
                        if compact_scope == "photon_and_charged_pion" else raw["selected"])
    row = {
        "sample_internal_name": name,
        "presentation_label": config["presentation_label"],
        "generator": config["generator"],
        "reconstruction_chain": config["reconstruction_chain"],
        "audit_scope": config["audit_scope"],
        "N_input_files": len(records),
        "N_events": raw["events"],
        "total_mc_particles": raw["mc"],
        "total_selected_truth": raw["selected"],
        "selected_electrons": selected_species["electron"],
        "selected_muons": selected_species["muon"],
        "selected_charged_pions": selected_species["charged_pion"],
        "selected_other": selected_species["other"],
        "selected_truth_per_event": raw["selected"] / raw["events"],
        "selected_truth_per_event_median": float(median(selected_per_event)),
        "selected_truth_parentless": raw["selected_parentless"],
        "selected_truth_parentless_fraction": frac(raw["selected_parentless"], raw["selected"]),
        "selected_photons": raw["photons"],
        "selected_tau_origin_photons": raw["photons_tau"],
        "selected_non_tau_photons": raw["photons_non_tau"],
        "tau_origin_photon_fraction": frac(raw["photons_tau"], raw["photons"]),
        "non_tau_photon_fraction": frac(raw["photons_non_tau"], raw["photons"]),
        "selected_parentless_photons": raw["photons_parentless"],
        "selected_parentless_photon_fraction": frac(raw["photons_parentless"], raw["photons"]),
        "total_pfos": raw["pfos"],
        "ldirect_assigned_any_mc": assigned,
        "ldirect_assigned_fraction_all_pfos": frac(assigned, raw["pfos"]),
        "ldirect_ambiguous_unresolved": raw["pfos"] - assigned,
        "ldirect_direct_mc_selected": same,
        "ldirect_direct_mc_not_selected": assigned - same,
        "lancestor_depth_0": depths["0"],
        "lancestor_depth_1": depths["1"],
        "lancestor_depth_2": depths["2"],
        "lancestor_depth_ge_3": sum(v for k, v in depths.items() if k not in {"<null>", "0", "1", "2"} and int(k) >= 3),
        "lancestor_promoted": promoted,
        "lancestor_promoted_fraction_assigned": frac(promoted, assigned),
        "lancestor_no_selected_ancestor": no_ancestor,
        "lancestor_unresolved_ambiguous": unresolved,
        "lancestor_cycle_error": errors,
        "lancestor_final_unique_valid_selected": valid,
        "lancestor_final_unique_valid_fraction_all_pfos": frac(valid, raw["pfos"]),
        "compact_selected_truth_rows": compact_selected,
        "compact_selected_truth_scope": compact_scope,
        "compact_selected_matches_declared_scope": compact_selected == compact_expected,
    }
    status_rows = ([{"sample_internal_name": name, "domain": "L_direct", "status": k, "count": v} for k, v in sorted(direct_status.items())]
                   + [{"sample_internal_name": name, "domain": "L_ancestor", "status": k, "count": v} for k, v in sorted(ancestor_status.items())])
    generator_rows = [{"sample_internal_name": name, "generatorStatus": k, "count": v} for k, v in sorted(gen_status.items(), key=lambda x: int(x[0]))]
    depth_rows = [{"sample_internal_name": name, "ancestor_depth": k, "count": v} for k, v in sorted(depths.items())]
    species_rows = [{"sample_internal_name": name, "species": k,
                     "selected_count": selected_species[k],
                     "selected_parentless_count": parentless_species[k]}
                    for k in sorted(set(selected_species) | set(parentless_species))]
    return row, generator_rows, status_rows, depth_rows, provenance, species_rows


def write_csv(path: Path, rows: list[dict]):
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def render_markdown(rows: list[dict]) -> str:
    out = ["# FCC sample truth/bookkeeping audit", "", "Generated from REC MC/PFO collections and frozen L_direct/L_ancestor assignments.", "Non-tau photons are selected photons not stored as tau descendants; no ISR interpretation is made.", "", "| Sample | Scope | Events | Selected truth | Photons (tau/non-tau) | PFOs | L_direct assigned | L_ancestor valid |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        out.append(f"| {r['sample_internal_name']} | {r['audit_scope']} | {r['N_events']:,} | {r['total_selected_truth']:,} | {r['selected_photons']:,} ({r['selected_tau_origin_photons']:,}/{r['selected_non_tau_photons']:,}) | {r['total_pfos']:,} | {r['ldirect_assigned_any_mc']:,} | {r['lancestor_final_unique_valid_selected']:,} |")
    out += ["", "Scientific identifiers: `selected_truth_v1`, `truthlink_assignment_v1`, `truthlink_ancestor_assignment_v1`.", "All counts are unweighted."]
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, default=REPO / "configs/analysis/fcc_sample_catalog_v1.yaml")
    parser.add_argument("--audit", nargs="+", default=["W", "P8C", "P8O", "KKMCee"])
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args(argv)
    if args.output_root.exists():
        raise FileExistsError(f"refusing to overwrite existing audit directory: {args.output_root}")
    catalog = load_catalog(args.catalog)
    unknown = set(args.audit) - set(catalog["samples"])
    if unknown:
        raise ValueError(f"unknown sample(s): {sorted(unknown)}")
    stage = Path(tempfile.mkdtemp(prefix="fcc-audit-", dir=str(args.output_root.parent)))
    try:
        summaries=[]; generators=[]; statuses=[]; depths=[]; provenance=[]; species=[]
        for name in args.audit:
            result = audit_sample(name, catalog["samples"][name], args.workers)
            summaries.append(result[0]); generators += result[1]; statuses += result[2]; depths += result[3]; provenance += result[4]; species += result[5]
        write_csv(stage / "truth_bookkeeping_audit.csv", summaries)
        (stage / "truth_bookkeeping_audit.json").write_text(json.dumps({"schema_version":"fcc_truth_bookkeeping_audit_v1", "generated_utc":datetime.now(timezone.utc).isoformat(), "samples":summaries}, indent=2) + "\n")
        (stage / "truth_bookkeeping_audit.md").write_text(render_markdown(summaries))
        write_csv(stage / "generator_status_distribution.csv", generators)
        write_csv(stage / "assignment_status_distribution.csv", statuses)
        write_csv(stage / "ancestor_depth_distribution.csv", depths)
        write_csv(stage / "input_provenance.csv", provenance)
        write_csv(stage / "parentless_selected_by_species.csv", species)
        shutil.move(str(stage), args.output_root)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
