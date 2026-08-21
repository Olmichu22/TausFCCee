"""
analysisRHOTree_MDecs_parallel.py

Versión reco-level del árbol MDecs. Lee ficheros EDM4hep, reconstruye ambos
hemisferios tau (hadrónicos via extractTauDecays + leptónicos via muonReco/
electronReco), los empareja con gen-taus por dR, y escribe un TTree con la
misma estructura tau1_*/tau2_* que genOnlyRHOTree_MDecs_parallel.py.

Compatible directamente con RhoHistFromTree_MDecs_parallel.py.

USO:
    python analysisRHOTree_MDecs_parallel.py \\
        -c config/default/taurecolong.yaml \\
        --decay-modes 0 2 10 -11 -13 \\
        --n-workers 4

--decay-modes acepta varias desintegraciones RECO (ρ=2, π=0, a1=10, e=-11, µ=-13)
y escribe un único árbol con todos los pares cuyos dos hemisferios estén en la
lista; el emparejamiento por par lo hace después RhoHistFromTree_MDecs (mismo
pipeline que el gen-only). Sin la flag se aceptan todas las desintegraciones.
"""

import ctypes
import logging
import math
import multiprocessing
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import ROOT
import yaml
from podio import root_io
from ROOT import TFile, TTree

from modules.TauDecays import extractTauDecays
from modules import myutils, tauReco
from modules import optimalVariabRho, weightsPol
from modules.rhoParallelUtils import (
    split_filenames, _make_p4_from_const, _setup_worker_logging, merge_partial_root_files,
)


# ── Config ────────────────────────────────────────────────────────────────────

_DEFAULT_CONFIG = "config/default/taurecolong.yaml"
_OUTPUT_BASE    = "Results/RhoAnalysis/"
# Subcarpeta dentro de <outputpath>/logs/ donde se guardan los logs de este script
_LOG_SOURCE     = "analysisRHOTree_MDecs"

_SHARED_VARIABS = ["GenZMass", "GenZVisMass", "ZMass", "beamE"]

_TAU_SCALAR_SUFFIXES = [
    "P", "E", "M", "Theta", "Phi",
    "visP", "visE", "visM", "visTheta", "visPhi",
    "pionP", "pionE", "pionM", "pionTheta", "pionPhi",
    "lepP", "lepE", "lepTheta", "lepPhi", "lepPDG",
    "decayID", "tauPDG", "cos_theta", "cos_psi", "cos_beta",
    "omega", "weight_P1", "weight_M1",
    "cos_theta_tau", "optimalVar", "isElectron", "nPhotons",
    "recoVisP", "recoVisE", "recoVisM", "recoVisTheta", "recoVisPhi",
    "recoPionP", "recoPionE", "recoPionM", "recoPionTheta", "recoPionPhi",
    "recoTauID", "recoCharge",
    "recoLepP", "recoLepE", "recoLepTheta", "recoLepPhi", "recoLepPDG",
    "reco_weight_P1", "reco_weight_M1",
]

_TAU_VECTOR_SUFFIXES = [
    "photons_E", "photons_theta", "photons_phi",
    "reco_photons_E", "reco_photons_theta", "reco_photons_phi",
]

_VARIABS = (
    _SHARED_VARIABS
    + [f"tau1_{s}" for s in _TAU_SCALAR_SUFFIXES]
    + [f"tau2_{s}" for s in _TAU_SCALAR_SUFFIXES]
)
_VECTOR_VARIABS = (
    [f"tau1_{s}" for s in _TAU_VECTOR_SUFFIXES]
    + [f"tau2_{s}" for s in _TAU_VECTOR_SUFFIXES]
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _remap_to_gen(decay_pair):
    """Convierte IDs reco a IDs gen: 2→1. Devuelve frozenset."""
    return frozenset(1 if d == 2 else d for d in decay_pair)


def split_mlpf(mlpf_results, file_chunks):
    mlpf_chunks = []
    file_offset = 0
    for chunk in file_chunks:
        lo = file_offset * 1000
        hi = (file_offset + len(chunk)) * 1000
        sub = {k - lo: v for k, v in mlpf_results.items() if lo <= k < hi}
        mlpf_chunks.append(sub)
        file_offset += len(chunk)
    return mlpf_chunks


def _extract_pion_p4(daughters):
    """Extrae el pion cargado (|PDG|==211) de un dict de daughters."""
    pion = ROOT.TLorentzVector()
    pion.SetXYZM(0, 0, 0, 0)
    if daughters is None:
        return pion
    for key in daughters:
        if abs(daughters[key].getPDG()) in (211, 321, 323): #SOLVED BUG
            pion = _make_p4_from_const(daughters[key])
            break
    return pion


def _tau_pdg_from_charge(charge):
    """PDG (15 / -15) of the tau in a reco hemisphere, from its reconstructed charge.

    The polarization P(z) is defined with z = cos(theta_tau-), so every weight
    needs to know which hemisphere it is looking at. |q| != 1 is not expected
    (see docs/plan_signo_costheta_taum.md, section 4.3); q >= 0 falls back to tau+
    and is counted separately in the worker log.
    """
    return 15 if charge < 0 else -15


def _build_reco_candidates(recoTaus, recoElectrons, recoMuons,
                            minPTauElectron, minPTauMuon):
    """
    Combina hadrónicos, muones y electrones en una lista uniforme.
    Cada entrada: {tau_id, vis_p4, pion_p4, lep_p4, is_electron, consts}
    """
    candidates = []

    for tau in recoTaus:
        consts = tau.getDaughters()
        candidates.append({
            "tau_id":      tau.getID(),
            "vis_p4":      tau.getMomentum(),
            "pion_p4":     _extract_pion_p4(consts),
            "lep_p4":      None,
            "is_electron": False,
            "consts":      consts,
            "charge":      tau.getCharge(),
        })

    for mu_key in recoMuons:
        mu_p4 = recoMuons[mu_key].getMomentum()
        if mu_p4.P() > minPTauMuon:
            candidates.append({
                "tau_id": -13, "vis_p4": mu_p4,
                "pion_p4": None, "lep_p4": mu_p4,
                "is_electron": False, "consts": {},
                "charge": recoMuons[mu_key].getCharge(),
            })

    for e_key in recoElectrons:
        e_p4 = recoElectrons[e_key].getMomentum()
        if e_p4.P() > minPTauElectron:
            candidates.append({
                "tau_id": -11, "vis_p4": e_p4,
                "pion_p4": None, "lep_p4": e_p4,
                "is_electron": True, "consts": {},
                "charge": recoElectrons[e_key].getCharge(),
            })

    return candidates


def _select_decay_pair(candidates, decay_pair):
    """
    Devuelve (cand0, cand1) cuyas tau_ids coinciden con decay_pair[0] y decay_pair[1].
    Para pares simétricos requiere 2 candidatos del mismo tipo.
    Ambos seleccionados por mayor P. Devuelve (None, None) si no hay par válido.
    """
    id0, id1 = decay_pair
    pool0 = [c for c in candidates if c["tau_id"] == id0]
    pool1 = [c for c in candidates if c["tau_id"] == id1]

    if id0 == id1:
        if len(pool0) < 2:
            return None, None
        sorted_pool = sorted(pool0, key=lambda c: c["vis_p4"].P(), reverse=True)
        return sorted_pool[0], sorted_pool[1]

    if not pool0 or not pool1:
        return None, None

    cand0 = max(pool0, key=lambda c: c["vis_p4"].P())
    cand1 = max(pool1, key=lambda c: c["vis_p4"].P())
    return cand0, cand1


def _select_decay_modes(candidates, genTaus, reco_filter):
    """
    Selecciona los dos hemisferios del evento entre los candidatos reco cuyos
    tau_id están en reco_filter (o todos si reco_filter es None).

    Importante: la ρ reco tiene tau_id == 2 (la gen usa 1). El filtro compara
    contra los ids RECO tal cual, así que para aceptar ρ hay que incluir 2.

    Para asegurar que los dos candidatos vienen de hemisferios distintos se usa
    el matching a gen-taus: se queda el candidato líder en P por cada gen-tau
    (nearest sin exclusión, descartando duplicados del mismo hemisferio).
    Sin gen-taus disponibles cae a los dos candidatos líderes en P.

    Devuelve (cand0, gen_idx0, cand1, gen_idx1) o None.
    """
    pool = [c for c in candidates
            if reco_filter is None or c["tau_id"] in reco_filter]
    if len(pool) < 2:
        return None
    pool.sort(key=lambda c: c["vis_p4"].P(), reverse=True)

    if len(genTaus) >= 2:
        chosen = {}  # gen_idx -> candidato líder en P de ese hemisferio
        for c in pool:
            gidx = _match_gen_tau(c["vis_p4"], genTaus, set())
            if gidx == -1 or gidx in chosen:
                continue
            chosen[gidx] = c
            if len(chosen) == 2:
                break
        if len(chosen) < 2:
            return None
        (g0, c0), (g1, c1) = list(chosen.items())
        return c0, g0, c1, g1

    return pool[0], -1, pool[1], -1


def _match_gen_tau(vis_p4, genTaus, used_indices):
    """Devuelve el índice gen más cercano en dR, excluyendo used_indices. -1 si vacío."""
    best_dr, best_idx = 10.0, -1
    for g, tau in genTaus.items():
        if g in used_indices:
            continue
        dr = myutils.dRAngle(tau.getMomentum(), vis_p4)
        if dr < best_dr:
            best_dr, best_idx = dr, g
    return best_idx


def _fill_tau_branches_mdecs(branches, prefix, reco_cand, gen_tau_obj, beamE, sin_eff):
    """
    Fill the {prefix}_* branches of one hemisphere.
    Gen branches  <- gen_tau_obj (if not None).
    Reco branches <- reco_cand.

    Returns 1 if the reconstructed charge of this hemisphere has |q| != 1, else 0,
    so the caller can keep a running count for the log (see the docs plan, 4.3).
    """
    # ── Reco branches ────────────────────────────────────────────────────────
    vis_p4  = reco_cand["vis_p4"]
    pion_p4 = reco_cand["pion_p4"] if reco_cand["pion_p4"] is not None else ROOT.TLorentzVector()
    lep_p4  = reco_cand["lep_p4"]  if reco_cand["lep_p4"]  is not None else ROOT.TLorentzVector()
    reco_id = reco_cand["tau_id"]

    reco_charge = reco_cand.get("charge", -1.0)
    reco_tau_pdg = _tau_pdg_from_charge(reco_charge)
    anomalous_charge = 0 if abs(reco_charge) == 1 else 1

    branches[f"{prefix}_recoTauID"].value     = float(reco_id)
    branches[f"{prefix}_recoCharge"].value    = float(reco_charge)
    branches[f"{prefix}_recoVisP"].value      = vis_p4.P()
    branches[f"{prefix}_recoVisE"].value      = vis_p4.E()
    branches[f"{prefix}_recoVisM"].value      = vis_p4.M()
    branches[f"{prefix}_recoVisTheta"].value  = vis_p4.Theta()
    branches[f"{prefix}_recoVisPhi"].value    = vis_p4.Phi()
    branches[f"{prefix}_recoPionP"].value     = pion_p4.P()
    branches[f"{prefix}_recoPionE"].value     = pion_p4.E()
    branches[f"{prefix}_recoPionM"].value     = pion_p4.M()
    branches[f"{prefix}_recoPionTheta"].value = pion_p4.Theta()
    branches[f"{prefix}_recoPionPhi"].value   = pion_p4.Phi()

    is_lep_reco = reco_id in (-11, -13)
    branches[f"{prefix}_recoLepP"].value     = lep_p4.P()     if is_lep_reco else 0.0
    branches[f"{prefix}_recoLepE"].value     = lep_p4.E()     if is_lep_reco else 0.0
    branches[f"{prefix}_recoLepTheta"].value = lep_p4.Theta() if is_lep_reco else 0.0
    branches[f"{prefix}_recoLepPhi"].value   = lep_p4.Phi()   if is_lep_reco else 0.0
    branches[f"{prefix}_recoLepPDG"].value   = float(abs(reco_id)) if is_lep_reco else 0.0

    # Fotones reco (hadrónicos)
    for const_key, const in reco_cand["consts"].items():
        if abs(const.getPDG()) == 22:
            ph = _make_p4_from_const(const)
            branches[f"{prefix}_reco_photons_E"].push_back(ph.E())
            branches[f"{prefix}_reco_photons_theta"].push_back(ph.Theta())
            branches[f"{prefix}_reco_photons_phi"].push_back(ph.Phi())

    # ── Reco-level polarization weights (collinear proxy) ────────────────────
    # The tau direction is not observable (neutrino escapes). Approximate it
    # by the visible system direction with E = E_beam (collinear approximation).
    # This gives x = E_vis_reco / E_beam, which is the correct energy fraction
    # for the Alcaraz (2026) weight formulas.
    rw_P1 = rw_M1 = 1.0
    if beamE > 0:
        tau_E_proxy = beamE
        tau_P_proxy = math.sqrt(max(0.0, tau_E_proxy**2 - weightsPol._M_TAU**2))
        tau_proxy_p4 = ROOT.TLorentzVector()
        tau_proxy_p4.SetPxPyPzE(
            tau_P_proxy * math.sin(vis_p4.Theta()) * math.cos(vis_p4.Phi()),
            tau_P_proxy * math.sin(vis_p4.Theta()) * math.sin(vis_p4.Phi()),
            tau_P_proxy * math.cos(vis_p4.Theta()),
            tau_E_proxy,
        )
        if reco_id in (1, 2):
            # ρ reco completo (2) y ρ sin un fotón (1) son ambos ρ: ESTÁNDAR = variable
            # óptima ω reco (wVariabRECO con cinemática reco), no la H simplificada.
            _, _, _, omega_reco = optimalVariabRho.wVariabRECO(vis_p4, pion_p4, beamE)
            rw_P1 = weightsPol.newAtauRhoOmega(tau_proxy_p4, omega_reco, +1,
                                               tau_pdg=reco_tau_pdg, sin_eff=sin_eff)
            rw_M1 = weightsPol.newAtauRhoOmega(tau_proxy_p4, omega_reco, -1,
                                               tau_pdg=reco_tau_pdg, sin_eff=sin_eff)
        elif reco_id in (0, 10):
            # π/a1: su observable óptimo ya es H (z para el π); no hay ω no trivial.
            rw_P1 = weightsPol.newAtau(tau_proxy_p4, vis_p4, reco_id, +1,
                                       tau_pdg=reco_tau_pdg, sin_eff=sin_eff)
            rw_M1 = weightsPol.newAtau(tau_proxy_p4, vis_p4, reco_id, -1,
                                       tau_pdg=reco_tau_pdg, sin_eff=sin_eff)
        elif reco_id in (-11, -13):
            rw_P1 = weightsPol.newAtauLep(lep_p4, tau_proxy_p4, beamE, +1,
                                          tau_pdg=reco_tau_pdg, sin_eff=sin_eff)
            rw_M1 = weightsPol.newAtauLep(lep_p4, tau_proxy_p4, beamE, -1,
                                          tau_pdg=reco_tau_pdg, sin_eff=sin_eff)
    branches[f"{prefix}_reco_weight_P1"].value = rw_P1
    branches[f"{prefix}_reco_weight_M1"].value = rw_M1

    # ── Gen branches ─────────────────────────────────────────────────────────
    if gen_tau_obj is None:
        # tauPDG y decayID llevan centinela, NO 0.0. Pasa de verdad y en volumen:
        # con menos de 2 gen-taus en el evento _select_decay_modes devuelve
        # gen_idx = -1 y ambos hemisferios se quedan sin gen — el caso de las
        # muestras con has_gen_taus: false (bhabha), donde es TODO el sample.
        #   - tauPDG: es la fuente del signo de z = cos(θ_τ⁻); un 0 se leería como
        #     τ⁺ en silencio en el hist stage, el modo de fallo que este convenio
        #     elimina. _gen_tau_pdg lanza ante cualquier cosa que no sea ±15.
        #   - decayID: un 0 es "pión gen", así que estos hemisferios se colaban
        #     como SIGNAL en el canal del pión (_tau_category compara decayID con
        #     el canal). Con -999 caen a BGOther, que es lo que son. -999 ya es el
        #     valor "desconocido" en todo el hist stage (los .get(..., -999)).
        #   - weight_P1/M1: 1.0, el neutro MULTIPLICATIVO, no 0.0. Sin gen-tau no
        #     hay repesado posible, y lo correcto es dejar el evento inalterado, no
        #     borrarlo del histograma. Además el hist stage cae a este producto para
        #     el peso joint cuando falta el gen (_compute_joint_weights), así que un
        #     0.0 aquí vaciaría en silencio las variantes corr_P1/corr_M1.
        _no_gen_defaults = {
            "omega": -999.0, "optimalVar": -999.0,
            "tauPDG": -999.0, "decayID": -999.0,
            "weight_P1": 1.0, "weight_M1": 1.0,
        }
        for sfx in _TAU_SCALAR_SUFFIXES:
            if sfx.startswith("reco"):
                continue
            val = _no_gen_defaults.get(sfx, 0.0)
            branches[f"{prefix}_{sfx}"].value = val
        return anomalous_charge

    decayID = gen_tau_obj.getID()
    tauP4   = gen_tau_obj.getMomentum()
    gen_vis = gen_tau_obj.getvisMomentum()

    branches[f"{prefix}_P"].value           = tauP4.P()
    branches[f"{prefix}_E"].value           = tauP4.E()
    branches[f"{prefix}_M"].value           = tauP4.M()
    branches[f"{prefix}_Theta"].value       = tauP4.Theta()
    branches[f"{prefix}_Phi"].value         = tauP4.Phi()
    branches[f"{prefix}_cos_theta_tau"].value = math.cos(tauP4.Theta())
    branches[f"{prefix}_visP"].value        = gen_vis.P()
    branches[f"{prefix}_visE"].value        = gen_vis.E()
    branches[f"{prefix}_visM"].value        = gen_vis.M()
    branches[f"{prefix}_visTheta"].value    = gen_vis.Theta()
    branches[f"{prefix}_visPhi"].value      = gen_vis.Phi()
    branches[f"{prefix}_decayID"].value     = float(decayID)
    gen_tau_pdg = int(gen_tau_obj.getPDG())
    branches[f"{prefix}_tauPDG"].value      = float(gen_tau_pdg)

    # Daughters: pion cargado + fotones gen
    daughters  = gen_tau_obj.getDaughters()
    gen_pion   = ROOT.TLorentzVector()
    gen_pion.SetXYZM(0, 0, 0, 0)
    for key in daughters:
        if abs(daughters[key].getPDG()) in (211, 321, 323): #SOLVED BUG
            gen_pion = _make_p4_from_const(daughters[key])
            break

    branches[f"{prefix}_pionP"].value     = gen_pion.P()
    branches[f"{prefix}_pionE"].value     = gen_pion.E()
    branches[f"{prefix}_pionM"].value     = gen_pion.M()
    branches[f"{prefix}_pionTheta"].value = gen_pion.Theta()
    branches[f"{prefix}_pionPhi"].value   = gen_pion.Phi()

    n_gen_photons = 0.0
    for key in daughters:
        if abs(daughters[key].getPDG()) == 111:
            for ph in daughters[key].getDaughters():
                if ph.getGeneratorStatus() != 1:
                    continue
                n_gen_photons += 1
                phP4 = _make_p4_from_const(ph)
                branches[f"{prefix}_photons_E"].push_back(phP4.E())
                branches[f"{prefix}_photons_theta"].push_back(phP4.Theta())
                branches[f"{prefix}_photons_phi"].push_back(phP4.Phi())
    branches[f"{prefix}_nPhotons"].value = n_gen_photons

    # Polarización (dispatch por decayID)
    if decayID == 1:
        (cos_theta, cos_psi, cos_beta, omega,
         w_P1, w_M1) = optimalVariabRho.wVariab(
            tauP4, gen_vis, gen_pion, beamE, sin_eff=sin_eff, tau_pdg=gen_tau_pdg)
    elif decayID == 0:  # π gen: peso con cos(θ*) geométrico (boost exacto, α=1)
        cts = weightsPol.cosThetaStar(tauP4, gen_vis)
        w_P1 = weightsPol.newAtauFromH(tauP4, cts, +1, tau_pdg=gen_tau_pdg, sin_eff=sin_eff)
        w_M1 = weightsPol.newAtauFromH(tauP4, cts, -1, tau_pdg=gen_tau_pdg, sin_eff=sin_eff)
        cos_theta = cts
        cos_psi = cos_beta = 0.0
        omega = -999.0
    elif decayID == 10:  # a1 gen: z_R vía newAtau
        w_P1 = weightsPol.newAtau(tauP4, gen_vis, decayID, +1, tau_pdg=gen_tau_pdg, sin_eff=sin_eff)
        w_M1 = weightsPol.newAtau(tauP4, gen_vis, decayID, -1, tau_pdg=gen_tau_pdg, sin_eff=sin_eff)
        cos_theta = cos_psi = cos_beta = 0.0
        omega = -999.0
    elif decayID in (-11, -13):
        w_P1 = weightsPol.newAtauLep(gen_vis, tauP4, beamE, +1, tau_pdg=gen_tau_pdg, sin_eff=sin_eff)
        w_M1 = weightsPol.newAtauLep(gen_vis, tauP4, beamE, -1, tau_pdg=gen_tau_pdg, sin_eff=sin_eff)
        cos_theta = cos_psi = cos_beta = 0.0
        omega = -999.0
    else:
        # Unsupported gen modes must not populate the physical omega spectrum.
        # Keep them as a sentinel so they stay out of the central histogram bin.
        cos_theta = cos_psi = cos_beta = 0.0
        omega = -999.0
        w_P1 = w_M1 = 1.0

    # Variable óptima (observable): definición ÚNICA vía helper compartido.
    opt_var = optimalVariabRho.optimal_var(decayID, tauP4, gen_vis, gen_pion, beamE,
                                           tau_pdg=gen_tau_pdg)

    branches[f"{prefix}_cos_theta"].value  = cos_theta
    branches[f"{prefix}_cos_psi"].value    = cos_psi
    branches[f"{prefix}_cos_beta"].value   = cos_beta
    branches[f"{prefix}_omega"].value      = omega
    branches[f"{prefix}_weight_P1"].value  = w_P1
    branches[f"{prefix}_weight_M1"].value  = w_M1
    branches[f"{prefix}_optimalVar"].value = opt_var

    # Leptón gen
    is_lep_gen = decayID in (-11, -13)
    branches[f"{prefix}_lepP"].value      = gen_vis.P()     if is_lep_gen else 0.0
    branches[f"{prefix}_lepE"].value      = gen_vis.E()     if is_lep_gen else 0.0
    branches[f"{prefix}_lepTheta"].value  = gen_vis.Theta() if is_lep_gen else 0.0
    branches[f"{prefix}_lepPhi"].value    = gen_vis.Phi()   if is_lep_gen else 0.0
    branches[f"{prefix}_lepPDG"].value    = float(abs(decayID)) if is_lep_gen else 0.0
    branches[f"{prefix}_isElectron"].value = float(decayID == -11)

    return anomalous_charge


# ── Worker ────────────────────────────────────────────────────────────────────

def process_chunk_mdecs(filenames_chunk, mlpf_chunk, global_event_offset,
                         config_bundle, worker_id):
    outputpath       = config_bundle["outputpath"]
    reco_filter      = config_bundle["reco_filter"]
    dRMax            = config_bundle["dRMax"]
    tauPCut          = config_bundle["tauPCut"]
    minPTauPhoton    = config_bundle["minPTauPhoton"]
    minPTauPion      = config_bundle["minPTauPion"]
    PNeutron         = config_bundle["PNeutron"]
    generalPCut      = config_bundle["generalPCut"]
    photon_config    = config_bundle["photon_config"]
    sin_eff          = config_bundle["sin_eff"]
    minPTauElectron  = config_bundle["minPTauElectron"]
    minPTauMuon      = config_bundle["minPTauMuon"]
    lepton_xor_p     = config_bundle.get("lepton_xor_p", 0.0)
    gen_taus_sample  = config_bundle.get("gen_taus_sample", True)
    gatr_results_path = config_bundle["gatr_results_path"]
    test_pfo         = config_bundle["test_pfo"]
    fileOutName_base  = config_bundle["fileOutName_base"]

    loggers        = _setup_worker_logging(outputpath, worker_id, _LOG_SOURCE)
    logger_process = loggers["processing"]
    logger_io      = loggers["io"]

    # Salida parcial
    partial_path = os.path.join(outputpath,
                                f"partial_worker{worker_id}_{fileOutName_base}.root")
    partial_outfile = TFile(partial_path, "RECREATE")
    tree = TTree("outtree_original", "MDecs reco+gen variables")

    branches = {}
    for var in _VARIABS:
        branches[var] = ctypes.c_double(0.0)
        tree.Branch(var, ctypes.addressof(branches[var]), f"{var}/D")
    for var in _VECTOR_VARIABS:
        branches[var] = ROOT.std.vector("double")()
        tree.Branch(var, branches[var])

    totalEvents    = 0
    selectedEvents = 0
    eventid        = global_event_offset
    skipped_files  = []
    # Hemisferios con |q_reco| != 1: no se espera ninguno (el signo de z sale de la
    # carga), pero se cuenta para que el caso no pase inadvertido si algún día
    # aparece con estadística no despreciable. Ver docs/plan_signo_costheta_taum.md.
    anomalousCharge = 0

    for filename in filenames_chunk:
        try:
            file_reader = root_io.Reader([filename])
        except Exception as exc:
            logger_io.warning("Worker %d: skipping %s: %s", worker_id, filename, exc)
            skipped_files.append(filename)
            continue

        for event in file_reader.get("events"):
            if totalEvents % 500 == 0:
                logger_process.info("Worker %d: processed %d events", worker_id, totalEvents)
            totalEvents += 1

            # Clear vectors
            for var in _VECTOR_VARIABS:
                branches[var].clear()

            mc_particles = event.get("MCParticles")
            beamE = mc_particles[0].getEnergy()

            # Gen taus
            genTaus = tauReco.findAllGenTaus(mc_particles) if gen_taus_sample else {}
            nGenTaus = len(genTaus)

            GenZMass = GenZVisMass = 0.0
            if nGenTaus == 2:
                gen_list = list(genTaus.values())
                ZGenP4  = gen_list[0].getMomentum() + gen_list[1].getMomentum()
                ZVisP4  = gen_list[0].getvisMomentum() + gen_list[1].getvisMomentum()
                GenZMass    = ZGenP4.M()
                GenZVisMass = ZVisP4.M()

            # Reco reconstruction
            pfos = event.get("PandoraPFOs")
            recoTau_raw, recoElectrons, recoMuons, _ = extractTauDecays(
                gatr_results_path, mlpf_chunk, eventid, pfos,
                dRMax, minPTauPhoton, minPTauPion, PNeutron,
                generalPCut, photon_config, False, test_pfo, logger_process,
            )
            recoTaus = myutils.sort_by_P(recoTau_raw)

            # Single-lepton XOR cut (replicates legacy: exactamente 1 electrón XOR 1 muón con P>threshold)
            if lepton_xor_p > 0.0:
                n_e  = sum(1 for e_k in recoElectrons  if recoElectrons[e_k].getMomentum().P()  > lepton_xor_p)
                n_mu = sum(1 for mu_k in recoMuons     if recoMuons[mu_k].getMomentum().P()     > lepton_xor_p)
                if not ((n_e == 1) ^ (n_mu == 1)):
                    eventid += 1
                    continue

            # Candidate list + selección por decay-modes (dos hemisferios distintos)
            candidates = _build_reco_candidates(
                recoTaus, recoElectrons, recoMuons, minPTauElectron, minPTauMuon)
            sel = _select_decay_modes(candidates, genTaus, reco_filter)
            if sel is None:
                eventid += 1
                continue
            cand0, gen_idx_0, cand1, gen_idx_1 = sel

            # Optional momentum cut
            if (cand0["vis_p4"].P() < tauPCut and cand1["vis_p4"].P() < tauPCut):
                eventid += 1
                continue

            # ZMass reco
            ZMass = (cand0["vis_p4"] + cand1["vis_p4"]).M()

            # Order tau1/tau2 by energy (tau1 = higher E); conserva su gen-match
            if cand0["vis_p4"].E() >= cand1["vis_p4"].E():
                tau1_cand, gen_tau1 = cand0, genTaus.get(gen_idx_0)
                tau2_cand, gen_tau2 = cand1, genTaus.get(gen_idx_1)
            else:
                tau1_cand, gen_tau1 = cand1, genTaus.get(gen_idx_1)
                tau2_cand, gen_tau2 = cand0, genTaus.get(gen_idx_0)

            # Fill shared branches
            branches["GenZMass"].value    = GenZMass
            branches["GenZVisMass"].value = GenZVisMass
            branches["ZMass"].value       = ZMass
            branches["beamE"].value       = beamE

            # Fill per-tau branches
            anomalousCharge += _fill_tau_branches_mdecs(
                branches, "tau1", tau1_cand, gen_tau1, beamE, sin_eff)
            anomalousCharge += _fill_tau_branches_mdecs(
                branches, "tau2", tau2_cand, gen_tau2, beamE, sin_eff)

            selectedEvents += 1
            tree.Fill()
            eventid += 1

    if skipped_files:
        logger_io.warning("Worker %d: %d file(s) skipped", worker_id, len(skipped_files))

    partial_outfile.cd()
    tree.Write()
    partial_outfile.Close()

    if anomalousCharge:
        logger_io.warning("Worker %d: %d reco hemisphere(s) with |q| != 1 "
                          "(sign of z taken as tau+ for those)",
                          worker_id, anomalousCharge)
    logger_io.info("Worker %d: done. Events=%d Selected=%d",
                   worker_id, totalEvents, selectedEvents)
    return partial_path, totalEvents, selectedEvents


# ── Merge ─────────────────────────────────────────────────────────────────────

# ── Main ──────────────────────────────────────────────────────────────────────

def my_hook(parser):
    parser.add_argument("--sin-eff", type=float, default=None,
                        help="Effective sin^2 theta_W for polarization weights")
    parser.add_argument("--n-workers", type=int, default=None,
                        help="Parallel workers (default: min(n_files, n_cpus))")
    parser.add_argument("--decay-modes", type=int, nargs="+", default=None,
                        metavar="ID",
                        help="Lista de decayIDs RECO permitidos (ej: 0 2 -11 -13 10). "
                             "La ρ reco es 2 (la gen es 1). Se acepta cualquier par "
                             "cuyos DOS hemisferios estén en la lista; el emparejamiento "
                             "por par lo hace después RhoHistFromTree_MDecs. "
                             "Default: todas las desintegraciones.")
    parser.add_argument("--decay-pair", type=int, nargs=2, default=None,
                        metavar=("ID0", "ID1"),
                        help="[DEPRECATED] Par reco único (ej: 2 -13). Si se da y no hay "
                             "--decay-modes, se trata como lista de dos modos.")
    parser.add_argument("-e", "--electron-cut", type=float, default=0.0,
                        help="Minimum electron P [GeV] (default: 0)")
    parser.add_argument("-u", "--muon-cut", type=float, default=0.0,
                        help="Minimum muon P [GeV] (default: 0)")
    parser.add_argument("--lepton-xor-p", type=float, default=0.0,
                        help="Si >0, aplica corte legacy: exactamente 1 electrón XOR 1 muón "
                             "con P > este valor [GeV]. Replica analysisRHOTree.py (default: 0 = desactivado)")
    parser.add_argument("--sys-err", type=str,
                        default="config/systematics/err_sys.yml")
    parser.add_argument("--test-extremes", action="store_true",
                        help="Test photon energy resolution extremes (unused in MDecs)")


def main():
    general_configs = myutils.setup_analysis_config(_DEFAULT_CONFIG, _OUTPUT_BASE,
                                                    parser_hook=my_hook,
                                                    log_subdir=_LOG_SOURCE)
    loggers    = general_configs["loggers"]
    run_config = general_configs["config"]
    args       = general_configs["args"]

    # decay_modes (ids RECO; ρ=2): CLI > YAML. Compat: --decay-pair → lista de 2.
    decay_modes = args.decay_modes or run_config.get("general", {}).get("decay_modes")
    if decay_modes is None and args.decay_pair:
        decay_modes = list(args.decay_pair)
    reco_filter = set(decay_modes) if decay_modes else None
    if reco_filter is not None:
        loggers["config"].info("decay_modes (ids reco, ρ=2): %s", sorted(reco_filter))
    else:
        loggers["config"].info("No decay_modes — se aceptan todas las desintegraciones")

    dRMax         = run_config["cuts"]["dRMax"]
    tauPCut       = run_config["cuts"]["tauCut"]
    minPTauPhoton = run_config["cuts"]["TauPhotonPCut"]
    minPTauPion   = run_config["cuts"]["TauPionPCut"]
    PNeutron      = run_config["cuts"]["NeutronCut"]
    generalPCut   = run_config["cuts"]["generalPCut"]
    outputpath    = general_configs["outputpath"]

    sys_errors    = run_config.get("systematics_errors", {})
    photon_config = sys_errors.get("photon_config", {})
    sample        = run_config["general"]["sample"]
    test_arg      = general_configs["flags"]["test"]
    gen_taus_sample = general_configs.get("has_gen_taus", True)

    gatr_results_path = args.gatr_result

    decay_tag        = "All" if not decay_modes else "_".join(str(d) for d in decay_modes)
    raw_fileOutName  = os.path.join(outputpath, general_configs["fileOutName"])
    fileOutName_base = f"TTree_MDecs_{decay_tag}_{Path(raw_fileOutName).stem}"
    fileOutName      = os.path.join(outputpath, f"{fileOutName_base}.root")

    filenames, mlpf_results = myutils.get_root_trees_path(
        sample, gatr_results_path, loggers, test_arg, args,
        skip_root_validation=True,
    )
    loggers["io"].info("Found %d input files", len(filenames))
    if not filenames:
        loggers["io"].error("No input files found. Aborting.")
        sys.exit(1)

    if args.input_list and len(args.input_list) == 1:
        n_workers = 1
    else:
        n_workers = args.n_workers or min(len(filenames), os.cpu_count() or 1)
    loggers["io"].info("Using %d workers for %d files", n_workers, len(filenames))

    config_bundle = {
        "reco_filter":      reco_filter,
        "dRMax":            dRMax,
        "tauPCut":          tauPCut,
        "minPTauPhoton":    minPTauPhoton,
        "minPTauPion":      minPTauPion,
        "PNeutron":         PNeutron,
        "generalPCut":      generalPCut,
        "photon_config":    photon_config,
        "sin_eff":          args.sin_eff,
        "minPTauElectron":  args.electron_cut,
        "minPTauMuon":      args.muon_cut,
        "lepton_xor_p":     args.lepton_xor_p,
        "gen_taus_sample":  gen_taus_sample,
        "gatr_results_path": gatr_results_path,
        "test_pfo":         args.test_pfo,
        "outputpath":       outputpath,
        "fileOutName_base": fileOutName_base,
    }

    totals = {"events": 0, "selected": 0}

    if n_workers == 1:
        loggers["io"].info("Sequential mode.")
        mlpf_empty = {}
        partial_path, tot, sel = process_chunk_mdecs(
            filenames, mlpf_empty, 0, config_bundle, 0)
        import shutil
        shutil.move(partial_path, fileOutName)
        totals["events"]   = tot
        totals["selected"] = sel
        loggers["io"].info("Sequential done: events=%d selected=%d", tot, sel)
    else:
        file_chunks   = split_filenames(filenames, n_workers)
        mlpf_chunks   = split_mlpf(mlpf_results, file_chunks)
        n_chunks      = len(file_chunks)
        event_offsets = [sum(len(file_chunks[j]) for j in range(i)) * 1000
                         for i in range(n_chunks)]

        ctx           = multiprocessing.get_context("fork")
        partial_files = []
        t_start       = time.time()

        loggers["io"].info("Launching %d workers...", n_chunks)
        with ProcessPoolExecutor(max_workers=n_chunks, mp_context=ctx) as executor:
            futures = {
                executor.submit(process_chunk_mdecs,
                                file_chunks[i], mlpf_chunks[i],
                                event_offsets[i], config_bundle, i): i
                for i in range(n_chunks)
            }
            for n_done, future in enumerate(as_completed(futures), start=1):
                wid = futures[future]
                try:
                    partial_path, tot, sel = future.result()
                    partial_files.append(partial_path)
                    totals["events"]   += tot
                    totals["selected"] += sel
                    elapsed = time.time() - t_start
                    print(f"  [{n_done}/{n_chunks}] worker {wid} done "
                          f"({tot} events, {elapsed:.1f}s)", flush=True)
                except Exception as exc:
                    loggers["io"].error("Worker %d failed: %s", wid, exc)
                    raise

        loggers["io"].info("Merging %d partial files → %s",
                           len(partial_files), fileOutName)
        merge_partial_root_files(partial_files, fileOutName)

        for p in partial_files:
            try:
                os.remove(p)
            except OSError:
                pass

        loggers["io"].info("Total events=%d selected=%d (%.1fs)",
                           totals["events"], totals["selected"],
                           time.time() - t_start)

    # Resumen CSV + config
    results_df = pd.DataFrame({
        "TotalEvents":    [totals["events"]],
        "SelectedEvents": [totals["selected"]],
    })
    results_df.to_csv(
        os.path.join(outputpath, f"results_summary_{decay_tag}.csv"), index=False)

    output_config_file = os.path.join(outputpath, "config.yaml")
    run_config["args"]       = vars(args)
    run_config["run_Type"]   = "analysisRHOTree_MDecs"
    run_config["general"]["decay_modes"] = decay_modes
    with open(output_config_file, "w") as f:
        yaml.dump(run_config, f)

    loggers["io"].info("Config saved to %s", output_config_file)
    loggers["io"].info("Output file: %s", fileOutName)


if __name__ == "__main__":
    main()
