import math

import pytest

from modules.fcc_event_record import (
    _mark_representatives,
    cluster_geometry_diagnostics,
    edge_geometry,
    invert_assignments,
    render_detail,
    render_genealogy,
    residual_match,
    theta_residual_mrad,
)


def point(x=0.0, y=0.0, z=0.0):
    return {"x_mm": x, "y_mm": y, "z_mm": z, "rxy_mm": math.hypot(x, y)}


def mc(index, *, pdg=22, parents=None, daughters=None, direct=None, ancestor=None,
       vertex=None, endpoint=None, theta=30.0):
    return {
        "index": index,
        "identity": f"000000001:7:MC:{index}",
        "pdg": pdg,
        "generator_status": 1,
        "energy_GeV": 5.0,
        "p_GeV": 5.0,
        "theta_deg": theta,
        "parents": list(parents or []),
        "daughters": list(daughters or []),
        "vertex": vertex or point(),
        "endpoint": endpoint or point(),
        "daughter_edges": [],
        "ldirect_pfos": list(direct or []),
        "lancestor_pfos": list(ancestor or []),
    }


def pfo(index, *, direct_mc=None, ancestor_mc=None, track=0, cluster=100,
        energy=4.9, theta=30.744845, charge=0.0, clusters=None, pfo_type=22):
    diagnostics = {"delta_theta_mrad": 13.0, "momentum_residual": -0.0125}
    return {
        "index": index,
        "identity": f"000000001:7:PFO:{index}",
        "type": pfo_type,
        "charge": charge,
        "energy_GeV": energy,
        "p_GeV": energy,
        "theta_deg": theta,
        "clusters": list(clusters or []),
        "ldirect": {
            "status": "assigned" if direct_mc is not None else "orphan_no_relation",
            "mc_index": direct_mc,
            "track_permille": track,
            "cluster_permille": cluster,
            "diagnostics": diagnostics if direct_mc is not None else None,
        },
        "lancestor": {
            "status": "same_direct_selected" if ancestor_mc is not None else "ancestor_no_selected_ancestor",
            "mc_index": ancestor_mc,
            "ancestor_depth": 0 if ancestor_mc is not None else None,
            "diagnostics": diagnostics if ancestor_mc is not None else None,
        },
    }


def record(particles, pfos=None):
    result = {
        "source_file_id": "000000001",
        "event_in_file": 7,
        "event_identity": "000000001:7",
        "mc_particles": particles,
        "pandora_pfos": list(pfos or []),
        "selection_matches": [],
    }
    _mark_representatives(result["pandora_pfos"], "L_direct")
    _mark_representatives(result["pandora_pfos"], "L_ancestor")
    return result


def test_inversion_keeps_all_pfos_for_one_mc():
    pfos = [
        {"index": 8, "ldirect": {"mc_index": 3}, "lancestor": {"mc_index": 1}},
        {"index": 2, "ldirect": {"mc_index": 3}, "lancestor": {"mc_index": 1}},
    ]
    assert invert_assignments(pfos, "L_direct") == {3: [2, 8]}


def test_direct_and_ancestor_may_point_to_different_mc_indices():
    item = {"index": 4, "ldirect": {"mc_index": 9}, "lancestor": {"mc_index": 2}}
    assert invert_assignments([item], "L_direct") == {9: [4]}
    assert invert_assignments([item], "L_ancestor") == {2: [4]}


def test_mc_detail_is_one_logical_line_per_particle():
    particles = [mc(0, daughters=[1], direct=[0]), mc(1, parents=[0], ancestor=[0])]
    text = render_detail(record(particles, [pfo(0, direct_mc=0, ancestor_mc=1)]))
    rows = [line for line in text.splitlines() if line.startswith("MC#")]
    assert len(rows) == 2
    assert all(" pdg=" in line and " parents=" in line and " daughters=" in line for line in rows)
    assert "D=[PFO#0(E=4.900000,theta=30.744845)[REP]]" in rows[0]
    assert "A=[PFO#0(E=4.900000,theta=30.744845)[REP]]" in rows[1]


def test_pfo_detail_has_compact_direct_and_ancestor_records():
    text = render_detail(record([mc(0, direct=[0], ancestor=[0])], [pfo(0, direct_mc=0, ancestor_mc=0)]))
    assert "PFO#0 type=22 q=0 E=4.900000 GeV theta=30.744845 deg" in text
    assert "D: status=assigned MC#0 pdg=22 gs=1" in text
    assert "A: status=same_direct_selected MC#0 pdg=22 gs=1" in text
    assert "dtheta=+13.000000 mrad dp/p=-0.012500" in text


def test_genealogy_uses_conventional_tree_indentation():
    particles = [
        mc(0, daughters=[1, 2]),
        mc(1, parents=[0], daughters=[3]),
        mc(2, parents=[0]),
        mc(3, parents=[1]),
    ]
    text = render_genealogy(record(particles))
    assert "├── MC#1 pdg=22" in text
    assert "│   └── MC#3 pdg=22" in text
    assert "└── MC#2 pdg=22" in text
    assert "+-->" not in text


def test_tree_annotations_show_all_pfo_kinematics_and_residuals():
    particles = [mc(0, direct=[0, 1], ancestor=[1])]
    pfos = [pfo(0, direct_mc=0), pfo(1, direct_mc=0, ancestor_mc=0)]
    text = render_genealogy(record(particles, pfos))
    assert "[D:PFO#0 E=4.900000 GeV theta=30.744845 deg dtheta=+13.000000 mrad dp/p=-0.012500]" in text
    assert "[D:PFO#1 E=4.900000 GeV theta=30.744845 deg dtheta=+13.000000 mrad dp/p=-0.012500]" in text
    assert "[A:PFO#1 E=4.900000 GeV theta=30.744845 deg dtheta=+13.000000 mrad dp/p=-0.012500]" in text


def test_shared_node_is_not_expanded_twice():
    particles = [
        mc(0, daughters=[1, 2]),
        mc(1, parents=[0], daughters=[3]),
        mc(2, parents=[0], daughters=[3]),
        mc(3, parents=[1, 2]),
    ]
    text = render_genealogy(record(particles))
    assert text.count("MC#3 pdg=22") == 1
    assert "MC#3 [already shown]" in text


def test_cycle_is_marked_and_terminates():
    particles = [mc(0, parents=[1], daughters=[1]), mc(1, parents=[0], daughters=[0])]
    text = render_genealogy(record(particles))
    assert "DISCONNECTED/CYCLIC COMPONENT" in text
    assert "MC#0 [cycle]" in text


def test_theta_residual_uses_degree_to_mrad_conversion():
    assert theta_residual_mrad(31.0, 30.0) == pytest.approx(math.pi / 180.0 * 1000.0)


def test_residual_predicate_selects_ten_to_twenty_mrad():
    particles = [mc(0)]
    item = {
        "index": 5,
        "theta_deg": 30.0 + math.degrees(0.013),
        "ldirect": {"mc_index": None},
        "lancestor": {"mc_index": 0},
    }
    match = residual_match(
        item,
        particles,
        pdg=22,
        method="L_ancestor",
        minimum_abs_mrad=10.0,
        maximum_abs_mrad=20.0,
    )
    assert match is not None
    assert match["delta_theta_mrad"] == pytest.approx(13.0)
    assert residual_match(
        item,
        particles,
        pdg=13,
        method="L_ancestor",
        minimum_abs_mrad=10.0,
        maximum_abs_mrad=20.0,
    ) is None


def test_maintained_representative_is_marked_without_hiding_other_pfos():
    pfos = [
        pfo(9, ancestor_mc=0, track=400, cluster=800, energy=99.0, theta=80.0, pfo_type=22),
        pfo(2, ancestor_mc=0, track=500, cluster=10, energy=0.1, theta=30.1, pfo_type=13),
        pfo(1, ancestor_mc=0, track=0, cluster=999, energy=999.0, theta=120.0, pfo_type=11),
    ]
    _mark_representatives(pfos, "L_ancestor")
    assert [item["lancestor"]["representative"] for item in pfos] == [False, True, False]
    text = render_detail(record([mc(0, ancestor=[1, 2, 9])], pfos))
    assert text.count("PFO#") >= 6
    mc_row = next(line for line in text.splitlines() if line.startswith("MC#0 "))
    assert "PFO#2(E=0.100000,theta=30.100000)[REP]" in mc_row
    assert "PFO#1(E=999.000000,theta=120.000000)[REP]" not in mc_row


def test_terminal_weight_tie_is_ambiguous_not_index_resolved():
    pfos = [pfo(8, ancestor_mc=0, track=3, cluster=7), pfo(2, ancestor_mc=0, track=3, cluster=7)]
    _mark_representatives(pfos, "L_ancestor")
    assert not any(item["lancestor"]["representative"] for item in pfos)
    assert {item["lancestor"]["representative_status"] for item in pfos} == {"ambiguous_multiple_pfo"}
    text = render_genealogy(record([mc(0, ancestor=[2, 8])], pfos))
    assert text.count("[REP:ambiguous_multiple_pfo]") == 2


def test_residual_search_default_any_pfo_and_representative_only():
    truth = [mc(0)]
    nonrep = pfo(4, ancestor_mc=0, theta=30.0 + math.degrees(0.013), track=1)
    rep = pfo(5, ancestor_mc=0, theta=30.0 + math.degrees(0.003), track=2)
    _mark_representatives([nonrep, rep], "L_ancestor")
    kwargs = dict(pdg=22, method="L_ancestor", minimum_abs_mrad=10.0, maximum_abs_mrad=20.0)
    assert residual_match(nonrep, truth, **kwargs) is not None
    assert residual_match(nonrep, truth, representative_only=True, **kwargs) is None
    assert residual_match(rep, truth, representative_only=True, **kwargs) is None
    rep["theta_deg"] = 30.0 + math.degrees(0.013)
    assert residual_match(rep, truth, representative_only=True, **kwargs) is not None


def test_mc_vertex_and_endpoint_are_rendered():
    particle = mc(0, vertex=point(1, 2, 3), endpoint=point(4, 5, 6))
    text = render_detail(record([particle]))
    assert "vtx=(1.000000,2.000000,3.000000;Rxy=2.236068) mm" in text
    assert "end=(4.000000,5.000000,6.000000;Rxy=6.403124) mm" in text


def test_neutral_single_cluster_geometry_is_numerically_correct():
    truth = mc(0, vertex=point(0, 0, -10), theta=90.0)
    item = pfo(0, ancestor_mc=0, theta=91.0, clusters=[{"index": 7, "position": point(100, 0, 0)}])
    result = cluster_geometry_diagnostics(item, truth)
    assert result["status"] == "single_cluster"
    assert result["dtheta_cluster_IP_mrad"] == pytest.approx(0.0)
    assert result["dtheta_cluster_from_MC_vertex_mrad"] == pytest.approx(
        (math.degrees(math.atan2(100, 10)) - 90.0) * math.pi / 180 * 1000
    )
    assert result["dtheta_PFO_vs_cluster_IP_mrad"] == pytest.approx(math.pi / 180 * 1000)


def test_cluster_geometry_charged_none_and_multiple_are_transparent():
    truth = mc(0)
    assert cluster_geometry_diagnostics(pfo(0, charge=1.0), truth) is None
    assert cluster_geometry_diagnostics(pfo(0), truth)["status"] == "no_cluster"
    clusters = [{"index": i, "position": point(i + 1, 0, 0)} for i in (3, 4)]
    result = cluster_geometry_diagnostics(pfo(0, clusters=clusters), truth)
    assert result == {"status": "multiple_clusters", "clusters": clusters}


def test_edge_geometry_distinguishes_coordinate_match_from_separation():
    parent = mc(0, endpoint=point(1, 2, 3))
    local = edge_geometry(parent, mc(1, vertex=point(1, 2, 3)))
    separated = edge_geometry(parent, mc(2, vertex=point(1, 2, 4)))
    assert local["local_coordinate_match"] and local["endpoint_vertex_separation_mm"] == 0
    assert not separated["local_coordinate_match"]
    assert separated["endpoint_vertex_separation_mm"] == pytest.approx(1.0)
