"""Plain-text inspection model for stored FCC MC genealogy and PFO assignments.

This module consumes persisted L_direct and L_ancestor rows. It does not
recalculate an assignment; representative annotations reuse the maintained
analysis helper over the persisted direct T/C weights.
"""
from __future__ import annotations

from collections import defaultdict
import math
from typing import Iterable

from modules.fcc_workflow_interface import truthlink_representative


LOCAL_VERTEX_TOLERANCE_MM = 1.0e-6


def momentum_theta(energy: float, momentum) -> dict[str, float]:
    px, py, pz = float(momentum.x), float(momentum.y), float(momentum.z)
    p = math.sqrt(px * px + py * py + pz * pz)
    theta = math.degrees(math.atan2(math.hypot(px, py), pz)) if p else 0.0
    return {"energy_GeV": float(energy), "p_GeV": p, "theta_deg": theta}


def theta_residual_mrad(reco_theta_deg: float, truth_theta_deg: float) -> float:
    """Frozen FCC-tau polar-angle residual; theta is deliberately not wrapped."""
    return (float(reco_theta_deg) - float(truth_theta_deg)) * math.pi / 180.0 * 1000.0


def momentum_residual(reco_p: float, truth_p: float) -> float | None:
    return None if float(truth_p) == 0.0 else (float(reco_p) - float(truth_p)) / float(truth_p)


def _point(position) -> dict[str, float]:
    x, y, z = float(position.x), float(position.y), float(position.z)
    return {"x_mm": x, "y_mm": y, "z_mm": z, "rxy_mm": math.hypot(x, y)}


def _point_theta_deg(point: dict[str, float], origin: dict[str, float] | None = None) -> float:
    origin = origin or {"x_mm": 0.0, "y_mm": 0.0, "z_mm": 0.0}
    dx = point["x_mm"] - origin["x_mm"]
    dy = point["y_mm"] - origin["y_mm"]
    dz = point["z_mm"] - origin["z_mm"]
    return math.degrees(math.atan2(math.hypot(dx, dy), dz))


def cluster_geometry_diagnostics(pfo: dict, truth: dict) -> dict | None:
    """Return diagnostic pointing quantities without changing the official residual."""
    if float(pfo["charge"]) != 0.0:
        return None
    clusters = pfo.get("clusters", [])
    if not clusters:
        return {"status": "no_cluster", "clusters": []}
    if len(clusters) != 1:
        return {"status": "multiple_clusters", "clusters": clusters}
    cluster = clusters[0]
    position = cluster["position"]
    theta_ip = _point_theta_deg(position)
    theta_vertex = _point_theta_deg(position, truth["vertex"])
    return {
        "status": "single_cluster",
        "cluster_index": cluster["index"],
        "cluster_position": position,
        "theta_cluster_IP_deg": theta_ip,
        "theta_cluster_from_MC_vertex_deg": theta_vertex,
        "dtheta_cluster_IP_mrad": theta_residual_mrad(theta_ip, truth["theta_deg"]),
        "dtheta_cluster_from_MC_vertex_mrad": theta_residual_mrad(theta_vertex, truth["theta_deg"]),
        "dtheta_PFO_vs_cluster_IP_mrad": theta_residual_mrad(pfo["theta_deg"], theta_ip),
    }


def edge_geometry(parent: dict, daughter: dict) -> dict:
    """Compare stored endpoint/vertex coordinates using numerical equality only."""
    endpoint, vertex = parent["endpoint"], daughter["vertex"]
    separation = math.sqrt(sum(
        (endpoint[axis] - vertex[axis]) ** 2 for axis in ("x_mm", "y_mm", "z_mm")
    ))
    return {
        "daughter_index": int(daughter["index"]),
        "endpoint_vertex_separation_mm": separation,
        "local_coordinate_match": separation <= LOCAL_VERTEX_TOLERANCE_MM,
        "local_coordinate_tolerance_mm": LOCAL_VERTEX_TOLERANCE_MM,
    }


def invert_assignments(pfos: Iterable[dict], method: str) -> dict[int, list[int]]:
    """Invert persisted unique PFO->MC indices without representative reduction."""
    key = "ldirect" if method == "L_direct" else "lancestor" if method == "L_ancestor" else None
    if key is None:
        raise ValueError(f"unsupported association method: {method}")
    inverse: dict[int, list[int]] = defaultdict(list)
    for pfo in pfos:
        mc_index = pfo[key].get("mc_index")
        if mc_index is not None:
            inverse[int(mc_index)].append(int(pfo["index"]))
    return {mc: sorted(indices) for mc, indices in sorted(inverse.items())}


def residual_match(
    pfo: dict,
    mc_particles: list[dict],
    *,
    pdg: int,
    method: str,
    minimum_abs_mrad: float,
    maximum_abs_mrad: float,
    representative_only: bool = False,
) -> dict | None:
    key = "ldirect" if method == "L_direct" else "lancestor" if method == "L_ancestor" else None
    if key is None:
        raise ValueError(f"unsupported association method: {method}")
    assignment = pfo[key]
    if representative_only and not assignment.get("representative", False):
        return None
    mc_index = assignment.get("mc_index")
    if mc_index is None or not 0 <= int(mc_index) < len(mc_particles):
        return None
    mc = mc_particles[int(mc_index)]
    if int(mc["pdg"]) != int(pdg):
        return None
    residual = theta_residual_mrad(pfo["theta_deg"], mc["theta_deg"])
    if not minimum_abs_mrad <= abs(residual) <= maximum_abs_mrad:
        return None
    return {
        "mc_index": int(mc_index),
        "pfo_index": int(pfo["index"]),
        "pdg": int(mc["pdg"]),
        "method": method,
        "truth_theta_deg": float(mc["theta_deg"]),
        "reco_theta_deg": float(pfo["theta_deg"]),
        "delta_theta_mrad": residual,
    }


def _object_index(obj) -> int:
    return int(obj.getObjectID().index)


def build_event_record(
    event,
    *,
    source_file_id: str,
    event_in_file: int,
    direct_rows: Iterable[dict],
    ancestor_rows: Iterable[dict],
) -> dict:
    mc_objects = list(event.get("MCParticles"))
    pfo_objects = list(event.get("PandoraPFOs"))

    direct = _rows_by_pfo(direct_rows, "L_direct")
    ancestor = _rows_by_pfo(ancestor_rows, "L_ancestor")
    expected_pfos = set(range(len(pfo_objects)))
    if set(direct) != expected_pfos:
        raise ValueError("L_direct rows do not cover PandoraPFO collection indices exactly")
    if set(ancestor) != expected_pfos:
        raise ValueError("L_ancestor rows do not cover PandoraPFO collection indices exactly")
    mc_particles = []
    for index, particle in enumerate(mc_objects):
        kine = momentum_theta(particle.getEnergy(), particle.getMomentum())
        mc_particles.append({
            "index": index,
            "identity": f"{source_file_id}:{event_in_file}:MC:{index}",
            "pdg": int(particle.getPDG()),
            "generator_status": int(particle.getGeneratorStatus()),
            "parents": [_object_index(parent) for parent in particle.getParents()],
            "daughters": [_object_index(child) for child in particle.getDaughters()],
            "vertex": _point(particle.getVertex()),
            "endpoint": _point(particle.getEndpoint()),
            **kine,
        })

    pfos = []
    for index, particle in enumerate(pfo_objects):
        kine = momentum_theta(particle.getEnergy(), particle.getMomentum())
        drow = direct.get(index)
        arow = ancestor.get(index)
        ldirect = _direct_assignment(drow)
        lancestor = _ancestor_assignment(arow)
        item = {
            "index": index,
            "identity": f"{source_file_id}:{event_in_file}:PFO:{index}",
            "type": int(particle.getPDG()),
            "charge": float(particle.getCharge()),
            "clusters": [
                {"index": _object_index(cluster), "position": _point(cluster.getPosition())}
                for cluster in particle.getClusters()
            ],
            **kine,
            "ldirect": ldirect,
            "lancestor": lancestor,
        }
        item["ldirect"]["diagnostics"] = _diagnostics(item, ldirect, mc_particles)
        item["lancestor"]["diagnostics"] = _diagnostics(item, lancestor, mc_particles)
        pfos.append(item)

    direct_inverse = invert_assignments(pfos, "L_direct")
    ancestor_inverse = invert_assignments(pfos, "L_ancestor")
    _mark_representatives(pfos, "L_direct")
    _mark_representatives(pfos, "L_ancestor")
    for mc in mc_particles:
        mc["ldirect_pfos"] = direct_inverse.get(mc["index"], [])
        mc["lancestor_pfos"] = ancestor_inverse.get(mc["index"], [])
        mc["daughter_edges"] = [
            edge_geometry(mc, mc_particles[child])
            for child in dict.fromkeys(mc["daughters"])
            if 0 <= child < len(mc_particles)
        ]

    return {
        "source_file_id": str(source_file_id),
        "event_in_file": int(event_in_file),
        "event_identity": f"{source_file_id}:{event_in_file}",
        "mc_particles": mc_particles,
        "pandora_pfos": pfos,
        "selection_matches": [],
    }


def _mark_representatives(pfos: list[dict], method: str) -> None:
    """Annotate every inversion group via the maintained analysis helper."""
    key = "ldirect" if method == "L_direct" else "lancestor" if method == "L_ancestor" else None
    if key is None:
        raise ValueError(f"unsupported association method: {method}")
    grouped: dict[int, list[dict]] = defaultdict(list)
    for pfo in pfos:
        assignment = pfo[key]
        assignment["representative"] = False
        assignment["representative_status"] = "unmatched"
        mc_index = assignment.get("mc_index")
        if mc_index is not None:
            grouped[int(mc_index)].append({
                "pfo_index": int(pfo["index"]),
                "track_permille": int(pfo["ldirect"].get("track_permille", 0)),
                "cluster_permille": int(pfo["ldirect"].get("cluster_permille", 0)),
            })
    by_index = {int(pfo["index"]): pfo for pfo in pfos}
    for rows in grouped.values():
        result = truthlink_representative(rows)
        for row in rows:
            by_index[row["pfo_index"]][key]["representative_status"] = result["status"]
        winner = result.get("pfo_index")
        if winner is not None:
            by_index[int(winner)][key]["representative"] = True


def _rows_by_pfo(rows: Iterable[dict], label: str) -> dict[int, dict]:
    result = {}
    for row in rows:
        pfo = int(row["pfo_index"])
        if pfo in result:
            raise ValueError(f"duplicate {label} row for PFO#{pfo}")
        result[pfo] = row
    return result


def _direct_assignment(row: dict | None) -> dict:
    if row is None:
        raise ValueError("missing L_direct assignment row")
    index = row.get("assigned_mc_index")
    return {
        "status": str(row["truthlink_status"]),
        "mc_index": None if index is None else int(index),
        "track_permille": int(row.get("track_permille") or 0),
        "cluster_permille": int(row.get("cluster_permille") or 0),
    }


def _ancestor_assignment(row: dict | None) -> dict:
    if row is None:
        raise ValueError("missing L_ancestor assignment row")
    index = row.get("ancestor_mc_index")
    result = {
        "status": str(row["ancestor_status"]),
        "mc_index": None if index is None else int(index),
    }
    if "ancestor_depth" in row:
        depth = row.get("ancestor_depth")
        result["ancestor_depth"] = None if depth is None else int(depth)
    return result


def _diagnostics(pfo: dict, assignment: dict, mc_particles: list[dict]) -> dict | None:
    index = assignment.get("mc_index")
    if index is None or not 0 <= int(index) < len(mc_particles):
        return None
    truth = mc_particles[int(index)]
    result = {
        "delta_theta_mrad": theta_residual_mrad(pfo["theta_deg"], truth["theta_deg"]),
        "momentum_residual": momentum_residual(pfo["p_GeV"], truth["p_GeV"]),
    }
    cluster = cluster_geometry_diagnostics(pfo, truth)
    if cluster is not None:
        result["neutral_cluster_geometry"] = cluster
    return result


def render_detail(record: dict) -> str:
    lines = [
        f"EVENT RECORD source_file_id={record['source_file_id']} event_in_file={record['event_in_file']}",
        f"EVENT IDENTITY {record['event_identity']}",
    ]
    for match in record.get("selection_matches", []):
        lines.extend([
            "", "SELECTION MATCH",
            f"  MC#{match['mc_index']} -> PFO#{match['pfo_index']}",
            f"  PDG: {match['pdg']}",
            f"  method: {match['method']}",
            f"  representative: {match.get('representative', False)}",
            f"  truth theta: {match['truth_theta_deg']:.6f} deg",
            f"  reco theta: {match['reco_theta_deg']:.6f} deg",
            f"  delta theta: {match['delta_theta_mrad']:+.6f} mrad",
        ])

    pfo_lookup = {pfo["index"]: pfo for pfo in record["pandora_pfos"]}
    lines.extend(["", "A. MC PARTICLES"])
    for mc in record["mc_particles"]:
        lines.append(
            f"MC#{mc['index']} pdg={mc['pdg']} gs={mc['generator_status']} "
            f"E={mc['energy_GeV']:.6f} GeV theta={mc['theta_deg']:.6f} deg "
            f"vtx={_format_point(mc['vertex'])} end={_format_point(mc['endpoint'])} "
            f"parents={mc['parents']} daughters={mc['daughters']} "
            f"D={_pfo_summaries(mc['ldirect_pfos'], pfo_lookup, 'ldirect')} "
            f"A={_pfo_summaries(mc['lancestor_pfos'], pfo_lookup, 'lancestor')}"
        )

    lines.extend(["", "B. PANDORA PFOs"])
    for pfo in record["pandora_pfos"]:
        lines.extend([
            f"PFO#{pfo['index']} type={pfo['type']} q={pfo['charge']:.6g} "
            f"E={pfo['energy_GeV']:.6f} GeV theta={pfo['theta_deg']:.6f} deg",
            f"  D: {_assignment_summary(pfo['ldirect'], record['mc_particles'])}",
            f"  A: {_assignment_summary(pfo['lancestor'], record['mc_particles'])}",
        ])
    return "\n".join(lines) + "\n"


def _format_point(point: dict) -> str:
    return (f"({point['x_mm']:.6f},{point['y_mm']:.6f},{point['z_mm']:.6f};"
            f"Rxy={point['rxy_mm']:.6f}) mm")


def _association_markers(assignment: dict) -> str:
    markers = "[REP]" if assignment.get("representative", False) else ""
    if assignment.get("search_match", False):
        markers += "[MATCH]"
    if assignment.get("representative_status") == "ambiguous_multiple_pfo":
        markers += "[REP:ambiguous_multiple_pfo]"
    return markers


def _pfo_summaries(indices: list[int], lookup: dict[int, dict], key: str) -> str:
    values = [
        f"PFO#{index}(E={lookup[index]['energy_GeV']:.6f},theta={lookup[index]['theta_deg']:.6f})"
        f"{_association_markers(lookup[index][key])}"
        for index in indices
    ]
    return "[" + ", ".join(values) + "]"


def _assignment_summary(assignment: dict, mc_particles: list[dict]) -> str:
    status, index = assignment["status"], assignment.get("mc_index")
    markers = _association_markers(assignment)
    if index is None:
        return f"status={status} MC=None {markers}".rstrip()
    if not 0 <= int(index) < len(mc_particles):
        return f"status={status} MC#{index} [invalid reference] {markers}".rstrip()
    mc = mc_particles[int(index)]
    text = (f"status={status} MC#{index} pdg={mc['pdg']} gs={mc['generator_status']} "
            f"E={mc['energy_GeV']:.6f} GeV theta={mc['theta_deg']:.6f} deg")
    diagnostics = assignment.get("diagnostics")
    if diagnostics:
        text += " " + _diagnostic_summary(diagnostics)
        text += _cluster_detail_suffix(diagnostics.get("neutral_cluster_geometry"))
    if "ancestor_depth" in assignment:
        text += f" depth={assignment['ancestor_depth']}"
    return (text + (" " + markers if markers else "")).rstrip()


def _diagnostic_summary(diagnostics: dict) -> str:
    momentum = diagnostics["momentum_residual"]
    momentum_text = "undefined" if momentum is None else f"{momentum:+.6f}"
    return f"dtheta={diagnostics['delta_theta_mrad']:+.6f} mrad dp/p={momentum_text}"


def _cluster_detail_suffix(cluster: dict | None) -> str:
    if cluster is None:
        return ""
    if cluster["status"] == "no_cluster":
        return " cluster_geometry=no_cluster"
    if cluster["status"] == "multiple_clusters":
        indices = [item["index"] for item in cluster["clusters"]]
        return f" cluster_geometry=multiple_clusters(indices={indices})"
    return (
        f" cluster#{cluster['cluster_index']} pos={_format_point(cluster['cluster_position'])}"
        f" theta_cluster_IP={cluster['theta_cluster_IP_deg']:.6f} deg"
        f" theta_cluster_from_vtx={cluster['theta_cluster_from_MC_vertex_deg']:.6f} deg"
        f" dtheta_cluster_IP={cluster['dtheta_cluster_IP_mrad']:+.6f} mrad"
        f" dtheta_cluster_from_vtx={cluster['dtheta_cluster_from_MC_vertex_mrad']:+.6f} mrad"
        f" dtheta_PFO_vs_cluster_IP={cluster['dtheta_PFO_vs_cluster_IP_mrad']:+.6f} mrad"
    )


def render_genealogy(record: dict) -> str:
    particles = record["mc_particles"]
    pfo_lookup = {pfo["index"]: pfo for pfo in record.get("pandora_pfos", [])}
    valid = set(range(len(particles)))
    roots = [mc["index"] for mc in particles if not any(parent in valid for parent in mc["parents"])]
    shown: set[int] = set()
    lines = [
        f"MC GENEALOGY source_file_id={record['source_file_id']} event_in_file={record['event_in_file']}",
    ]

    def visit(index: int, prefix: str, connector: str, active: set[int], edge: dict | None = None) -> None:
        position = prefix + connector
        if index not in valid:
            lines.append(f"{position}MC#{index} [invalid reference]")
            return
        if index in active:
            lines.append(f"{position}MC#{index} [cycle]{_edge_suffix(edge)}")
            return
        if index in shown:
            lines.append(f"{position}MC#{index} [already shown]{_edge_suffix(edge)}")
            return
        mc = particles[index]
        lines.append(position + _tree_node(mc) + _edge_suffix(edge))
        shown.add(index)
        next_active = set(active)
        next_active.add(index)
        child_prefix = prefix + ("│   " if connector == "├── " else "    " if connector == "└── " else "")
        annotation_prefix = child_prefix if connector else prefix + "    "
        if mc.get("ldirect_pfos") or mc.get("lancestor_pfos"):
            lines.append(annotation_prefix + f"vtx={_format_point(mc['vertex'])}")
            lines.append(annotation_prefix + f"end={_format_point(mc['endpoint'])}")
        for annotation in _tree_annotations(mc, pfo_lookup):
            lines.append(annotation_prefix + annotation)
        edges = {item["daughter_index"]: item for item in mc.get("daughter_edges", [])}
        children = list(dict.fromkeys(int(child) for child in mc["daughters"]))
        for offset, child in enumerate(children):
            child_connector = "└── " if offset == len(children) - 1 else "├── "
            visit(child, child_prefix, child_connector, next_active, edges.get(child))

    for offset, root in enumerate(roots):
        if offset:
            lines.append("")
        visit(root, "", "", set())
    for index in range(len(particles)):
        if index not in shown:
            lines.append("\nDISCONNECTED/CYCLIC COMPONENT")
            visit(index, "", "", set())
    return "\n".join(lines) + "\n"


def _edge_suffix(edge: dict | None) -> str:
    if edge is None:
        return ""
    separation = edge["endpoint_vertex_separation_mm"]
    label = "local-coordinate-match" if edge["local_coordinate_match"] else "non-local-coordinate-match"
    return f" [edge {label} separation={separation:.6g} mm]"


def _tree_node(mc: dict) -> str:
    return (f"MC#{mc['index']} pdg={mc['pdg']} gs={mc['generator_status']} "
            f"E={mc['energy_GeV']:.6f} GeV theta={mc['theta_deg']:.6f} deg")


def _tree_annotations(mc: dict, lookup: dict[int, dict]) -> list[str]:
    annotations = []
    for label, indices, key in (
        ("D", mc.get("ldirect_pfos", []), "ldirect"),
        ("A", mc.get("lancestor_pfos", []), "lancestor"),
    ):
        for index in indices:
            pfo = lookup[index]
            diagnostics = pfo[key].get("diagnostics")
            diagnostic = "" if diagnostics is None else " " + _diagnostic_summary(diagnostics)
            cluster = "" if diagnostics is None else _tree_cluster_summary(
                diagnostics.get("neutral_cluster_geometry")
            )
            annotations.append(
                f"[{label}:PFO#{index} E={pfo['energy_GeV']:.6f} GeV "
                f"theta={pfo['theta_deg']:.6f} deg{diagnostic}{cluster}]"
                f"{_association_markers(pfo[key])}"
            )
    return annotations


def _tree_cluster_summary(cluster: dict | None) -> str:
    if cluster is None:
        return ""
    if cluster["status"] == "no_cluster":
        return " cluster_geometry=no_cluster"
    if cluster["status"] == "multiple_clusters":
        indices = [item["index"] for item in cluster["clusters"]]
        return f" cluster_geometry=multiple_clusters(indices={indices})"
    return (
        f" dtheta_cluster_IP={cluster['dtheta_cluster_IP_mrad']:+.6f} mrad"
        f" dtheta_cluster_from_vtx={cluster['dtheta_cluster_from_MC_vertex_mrad']:+.6f} mrad"
        f" dtheta_PFO_vs_cluster_IP={cluster['dtheta_PFO_vs_cluster_IP_mrad']:+.6f} mrad"
    )
