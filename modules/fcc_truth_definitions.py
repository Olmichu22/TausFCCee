"""Frozen FCC analysis truth and reconstructed-PID definitions.

These helpers contain no workflow implementation and consume only values or
EDM-like particle objects.  Versions correspond to the 2026-09-01 scientific
freeze.
"""
from __future__ import annotations

from collections import deque
import math


SELECTED_TRUTH_VERSION = "selected_truth_v1"
PHOTON_ORIGIN_VERSION = "photon_origin_v1"
PID_MAPPING_VERSION = "fcc_reco_pid_v1"
TAU_ORIGIN_VERSION = "stored_parent_tau_origin_v1"

NEUTRINO_ABS_PDGS = frozenset({12, 14, 16})
MIN_SELECTED_MOMENTUM_GEV = 1.0e-10
RECO_PID_BY_ABS_PDG = {
    11: "electron",
    13: "muon",
    22: "photon",
    211: "charged_pion",
    310: "K0S",
    2112: "neutron",
    3122: "Lambda",
}


def selected_truth_values(generator_status: int, pdg: int, momentum) -> bool:
    """Return the frozen G selected-truth decision from primitive values."""
    norm = math.sqrt(sum(float(value) ** 2 for value in momentum))
    return (
        int(generator_status) == 1
        and abs(int(pdg)) not in NEUTRINO_ABS_PDGS
        and norm >= MIN_SELECTED_MOMENTUM_GEV
    )


def selected_truth_particle(particle) -> bool:
    """Apply :func:`selected_truth_values` to an EDM4hep-like MCParticle."""
    momentum = particle.getMomentum()
    return selected_truth_values(
        particle.getGeneratorStatus(),
        particle.getPDG(),
        (momentum.x, momentum.y, momentum.z),
    )


def reconstructed_pid_category(pdg: int) -> str:
    """Map a physical reconstructed PDG to the frozen analysis category.

    Sentinel 999 represents an association outcome and is deliberately not a
    reconstructed PID.
    """
    value = abs(int(pdg))
    if value == 999:
        raise ValueError("999 is an association sentinel, not a reconstructed PID")
    try:
        return RECO_PID_BY_ABS_PDG[value]
    except KeyError as error:
        raise ValueError(f"unsupported reconstructed PDG: {pdg}") from error


def _ancestor_indices(parents: list[list[int]], start: int) -> set[int]:
    """Return stored recursive ancestors, rejecting cycles and bad indices."""
    size = len(parents)
    if not 0 <= int(start) < size:
        raise IndexError(start)
    found: set[int] = set()
    pending = deque((int(parent), frozenset({int(start)})) for parent in parents[int(start)])
    while pending:
        index, path = pending.popleft()
        if not 0 <= index < size:
            raise IndexError(index)
        if index in path:
            raise ValueError("cycle in MC ancestry")
        if index in found:
            continue
        found.add(index)
        pending.extend((int(parent), path | {index}) for parent in parents[index])
    return found


def has_tau_origin(pdgs: list[int], parents: list[list[int]], start: int) -> bool:
    """Whether stored recursive ancestry contains a particle with |PDG|=15."""
    if len(pdgs) != len(parents):
        raise ValueError("PDG and parent collections have different lengths")
    return any(abs(int(pdgs[index])) == 15 for index in _ancestor_indices(parents, start))


def descriptive_photon_origin(
    *, tau_origin: bool, pi0_ancestor: bool, direct_tau_daughter: bool,
    immediate_parent_pdgs: list[int],
) -> str:
    """Frozen descriptive category used before photon_origin_v1."""
    if tau_origin:
        if pi0_ancestor:
            return "tau_decay_pi0"
        if direct_tau_daughter:
            return "tau_direct_daughter"
        return "tau_other_descendant"
    if not immediate_parent_pdgs:
        return "parentless_non_tau"
    if len(immediate_parent_pdgs) > 1:
        return "multiple_parent_non_tau"
    return {
        11: "electron_parent_non_tau",
        -11: "positron_parent_non_tau",
        22: "photon_chain_non_tau",
    }.get(int(immediate_parent_pdgs[0]), "other_resolved_non_tau")


def photon_origin_v1(descriptive_class: str, beam_chain_positive: bool) -> str:
    """Return the frozen photon_origin_v1 label without ISR/FSR reinterpretation."""
    if descriptive_class == "tau_decay_pi0":
        return "tau_decay_pi0"
    if descriptive_class == "tau_direct_daughter":
        return "tau_direct_daughter_unresolved"
    if descriptive_class == "tau_other_descendant":
        return "tau_other_descendant"
    if beam_chain_positive:
        return "explicit_ISR"
    if descriptive_class == "parentless_non_tau":
        return "parentless_unresolved"
    if descriptive_class == "electron_parent_non_tau":
        return "electron_parent_non_tau_unresolved"
    if descriptive_class == "positron_parent_non_tau":
        return "positron_parent_non_tau_unresolved"
    if descriptive_class == "photon_chain_non_tau":
        return "photon_chain_non_tau_unresolved"
    return "resolved_non_tau_other"
