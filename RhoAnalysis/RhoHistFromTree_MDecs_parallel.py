#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RhoHistFromTree_MDecs_parallel.py

Versión MDecs de rhoHistFromTree_parallel.py.

Lee árboles con variables tau1_*/tau2_* (generados por genOnlyRHOTree_parallel.py)
y genera histogramas para pares de desintegraciones configurables (decay_pair).

Las categorías usan doble sufijo {cat_este}_{cat_otro} por hemisferio:
  - Primer sufijo: estado del hemisferio al que pertenece la variable
  - Segundo sufijo: estado del hemisferio opuesto
La clasificación usa decayID para histogramas Gen/Matched, recoTauID para Reco.

Los histogramas se nombran con _dec0 / _dec1 según el hemisferio.
Los histogramas compartidos (ZMass, cross-hemisphere) no llevan sufijo de hemisferio.

Requiere en el YAML de config:
  general:
    decay_pair: [decID_0, decID_1]   # e.g. [2, -13] para rho + muon

USO:
    python RhoAnalysis/RhoHistFromTree_MDecs_parallel.py \\
        --tree-file Results/RhoAnalysis/.../tau_trainedAll_*.root \\
        -c config/default/taurecolong_mdecs.yaml \\
        --hist-config-mdecs config/histograms/rho_analysis_config_mdecs.yml \\
        -d -777 -v --n-workers 4
"""
import logging
import math
import multiprocessing
import numpy as np
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import yaml
import ROOT
from ROOT import TFile

from modules import myutils, weightsPol, optimalVariabRho
from modules.rhoTreeUtils import extract_scalars_optional, make_p4

# ── Constantes ────────────────────────────────────────────────────────────────

_DEFAULT_CONFIG = "config/default/taurecolong.yaml"
_OUTPUT_BASE    = "Results/RhoAnalysis/"
# Subcarpeta dentro de <outputpath>/logs/ donde se guardan los logs de este script
_LOG_SOURCE     = "RhoHistFromTree_MDecs"

_SHARED_MDECS = [
    ("GenZMass",    "GenZMass",    float),
    ("GenZVisMass", "GenZVisMass", float),
    ("ZRecoMass",   "ZMass",       float),
    ("beamE",       "beamE",       float),
]

_TAU_KEYS_MDECS = [
    "P", "E", "M", "Theta", "Phi",
    "visP", "visE", "visM", "visTheta", "visPhi",
    "pionP", "pionE", "pionM", "pionTheta", "pionPhi",
    "lepP", "lepE", "lepTheta", "lepPhi", "lepPDG",
    "decayID", "tauPDG", "genHelicity", "cos_theta", "cos_psi", "cos_beta",
    "omega", "weight_P1", "weight_M1",
    "cos_theta_tau", "optimalVar", "isElectron", "nPhotons",
    "recoVisP", "recoVisE", "recoVisM", "recoVisTheta", "recoVisPhi",
    "recoPionP", "recoPionE", "recoPionM", "recoPionTheta", "recoPionPhi",
    "recoTauID", "recoCharge",
    "recoLepP", "recoLepE", "recoLepTheta", "recoLepPhi", "recoLepPDG",
]

# Mapa de decayID a nombre de categoría BG
_BG_MAP = {-13: "BGMuon", -11: "BGEle", 0: "BGPion", 1: "BGRho", 10: "BGA1"}

# ── Pesos por hemisferio ───────────────────────────────────────────────────────
# Firma: (dec_vars, other_vars, base_weight) → effective_weight
# dec_vars: vars del hemisferio cuya variable estamos llenando
# other_vars: vars del otro hemisferio

WEIGHT_VALUES_MDECS = {
    "nominal":      lambda dv, ov, w: w,
    # Gen-level per-tau weights (exact tau kinematics)
    "P1":           lambda dv, ov, w: w * dv.get("weight_P1",           1.0),
    "M1":           lambda dv, ov, w: w * dv.get("weight_M1",           1.0),
    # Gen-level joint two-tau weights — Alcaraz (2026) eqs. (9)/(13), includes cross-term
    "corr_P1":      lambda dv, ov, w: w * dv.get("weight_corr_P1",      1.0),
    "corr_M1":      lambda dv, ov, w: w * dv.get("weight_corr_M1",      1.0),
    "other_P1":     lambda dv, ov, w: w * ov.get("weight_P1",           1.0),
    "other_M1":     lambda dv, ov, w: w * ov.get("weight_M1",           1.0),
    # Reco-level per-tau weights (visible-system proxy for tau direction)
    "reco_P1":      lambda dv, ov, w: w * dv.get("reco_weight_P1",      1.0),
    "reco_M1":      lambda dv, ov, w: w * dv.get("reco_weight_M1",      1.0),
    # Reco-level joint two-tau weights
    "reco_corr_P1": lambda dv, ov, w: w * dv.get("reco_weight_corr_P1", 1.0),
    "reco_corr_M1": lambda dv, ov, w: w * dv.get("reco_weight_corr_M1", 1.0),
}

# ── Reglas de llenado por hemisferio ──────────────────────────────────────────
# Formato: (level, base_name, x_fn(v, sh), y_fn(v, sh) | None, [cond_fn(v, sh)])
# v = vars del hemisferio, sh = shared_vars
# El nombre real en el YAML es base_name + "_dec0" o "_dec1".

# Predicados de carga del tau de ESTE hemisferio (PDG: 15 = τ⁻, −15 = τ⁺).
def _has_gen_tau(v):
    """True if this hemisphere has a matched gen tau (tauPDG == +-15).

    Hemispheres without one carry the -999 sentinel written by
    analysisRHOTree_MDecs_parallel.py (events with fewer than 2 gen taus, i.e.
    non-tautau backgrounds). Nothing gen-level can be computed for them.
    """
    return abs(int(v.get("tauPDG", 0))) == 15


def _gen_tau_pdg(v):
    """PDG (15 / -15) of the gen tau of this hemisphere, for the sign of z.

    P(z) is always defined with z = cos(theta_tau-), so every weight needs the
    hemisphere charge; see weightsPol._z_taum and docs/plan_signo_costheta_taum.md.

    Raises on anything that is not +-15 instead of defaulting to a charge: an
    unmatched hemisphere silently read as a tau+ is the exact failure mode this
    convention removes. Guard the call with _has_gen_tau.
    """
    pdg = int(v.get("tauPDG", 0))
    if abs(pdg) != 15:
        raise ValueError(
            f"tauPDG={pdg} is not +-15: this hemisphere has no matched gen tau, "
            "so the sign of z = cos(theta_tau-) is undefined. Guard with "
            "_has_gen_tau before asking for it."
        )
    return pdg


def _reco_tau_sign(v):
    """+1 if this hemisphere is a tau- (charge < 0), -1 if it is a tau+.

    Multiplies cos(theta_vis_reco) to build z = cos(theta_tau-) at reco level.
    |q| != 1 is assumed not to happen (docs plan, 4.3); q >= 0 reads as tau+.
    """
    return 1.0 if v.get("recoCharge", -1.0) < 0 else -1.0


def _reco_tau_pdg(v):
    """PDG (15 / -15) of the tau of this reco hemisphere, from recoCharge."""
    return 15 if v.get("recoCharge", -1.0) < 0 else -15


def _is_taup(v):  # τ⁺  (carga +1)
    return int(v.get("tauPDG", 0)) == -15
def _is_taum(v):  # τ⁻  (carga −1)
    return int(v.get("tauPDG", 0)) == 15

# Helicidad COMBINADA en carga, definida según la helicidad del τ⁻ (convención Cepeda):
# agrupa el τ⁻ con su pareja τ⁺ (de helicidad cruda opuesta por la correlación de espín).
# HelNeg = τ⁻ con hel<0 (∪ τ⁺ con hel>0); HelPos = τ⁻ con hel>0 (∪ τ⁺ con hel<0).
def _is_helneg(v):  # helicidad del τ⁻ negativa
    pdg = int(v.get("tauPDG", 0)); hel = v.get("genHelicity", 0.0)
    return (pdg == 15 and hel < 0) or (pdg == -15 and hel > 0)
def _is_helpos(v):  # helicidad del τ⁻ positiva
    pdg = int(v.get("tauPDG", 0)); hel = v.get("genHelicity", 0.0)
    return (pdg == 15 and hel > 0) or (pdg == -15 and hel < 0)

FILL_RULES_PER_DEC = [
    # ── Gen ───────────────────────────────────────────────────────────────────
    ("Gen",     "Omega",              lambda v, sh: v["omega"],                    None),
    ("Gen",     "CosTheta_GEN",       lambda v, sh: v["cos_theta"],               None),
    ("Gen",     "CosPsi_GEN",         lambda v, sh: v["cos_psi"],                 None),
    ("Gen",     "CosThetaTau_GEN",    lambda v, sh: v["cos_theta_tau"],           None),
    ("Gen",     "CosThetaMeson_GEN",  lambda v, sh: math.cos(v["visTheta"]),      None),
    ("Gen",     "CosThetaMeson_Reco", lambda v, sh: math.cos(v["recoVisTheta"]),  None),
    ("Gen",     "DecayType",          lambda v, sh: v["decayID"],                 None),
    ("Gen",     "OptimalVar",         lambda v, sh: v["optimalVar"],              None),
    ("Gen",     "VisP",               lambda v, sh: v["visP"],                    None),
    ("Gen",     "TauP",               lambda v, sh: v["P"],                       None),
    ("Gen",     "Omega_ZGenMass",     lambda v, sh: v["omega"],                   lambda v, sh: sh["GenZMass"]),
    ("Gen",     "Omega_ZVisMass",     lambda v, sh: v["omega"],                   lambda v, sh: sh["GenZVisMass"]),
    # Condicionados: solo cuando omega < 0
    ("Gen",     "VisP_omega_neg",     lambda v, sh: v["visP"],                    None,
                                      lambda v, sh: v["omega"] < 0),
    ("Gen",     "TauP_omega_neg",     lambda v, sh: v["P"],                       None,
                                      lambda v, sh: v["omega"] < 0),
    ("Gen",     "TauTheta_omega_neg", lambda v, sh: v["Theta"],                   None,
                                      lambda v, sh: v["omega"] < 0),
    ("Gen",     "TauPhi_omega_neg",   lambda v, sh: v["Phi"],                     None,
                                      lambda v, sh: v["omega"] < 0),
    # Omega en bins de masa resonante del rho (visM gen) — closure por masa (test 6.1.1)
    ("Gen",     "Omega_RhoMass_Bin1", lambda v, sh: v["omega"],                   None,
                                      lambda v, sh: 0.40 <= v["visM"] < 0.60),
    ("Gen",     "Omega_RhoMass_Bin2", lambda v, sh: v["omega"],                   None,
                                      lambda v, sh: 0.60 <= v["visM"] < 0.70),
    ("Gen",     "Omega_RhoMass_Bin3", lambda v, sh: v["omega"],                   None,
                                      lambda v, sh: 0.70 <= v["visM"] < 0.80),
    ("Gen",     "Omega_RhoMass_Bin4", lambda v, sh: v["omega"],                   None,
                                      lambda v, sh: 0.80 <= v["visM"] < 0.90),
    ("Gen",     "Omega_RhoMass_Bin5", lambda v, sh: v["omega"],                   None,
                                      lambda v, sh: 0.90 <= v["visM"] < 1.05),
    ("Gen",     "Omega_RhoMass_Bin6", lambda v, sh: v["omega"],                   None,
                                      lambda v, sh: 1.05 <= v["visM"] < 1.50),
    # ── Reco ──────────────────────────────────────────────────────────────────
    ("Reco",    "Omega_Reco",         lambda v, sh: v.get("_omega_reco", -999.0),  None),
    ("Reco",    "CosTheta",           lambda v, sh: v["cos_theta"],               None),
    ("Reco",    "CosPsi",             lambda v, sh: v["cos_psi"],                 None),
    ("Reco",    "RecoVisEOverBeamE",  lambda v, sh: v["recoVisE"] / sh["beamE"]
                                                    if sh["beamE"] else 0.0,      None),
    ("Reco",    "RecoVisCosTheta",    lambda v, sh: math.cos(v["recoVisTheta"]),  None),
    ("Reco",    "RecoVisP",           lambda v, sh: v["recoVisP"],                None),
    ("Reco",    "RecoDecayType",      lambda v, sh: float(v["recoTauID"]),        None),
    ("Reco",    "Optimal_X",          lambda v, sh: v.get("_optimal_x", 0.0),    None),
    # 2D var-óptima-reco × cosθ-visible-reco para el fit de polarización (A_τ/A_e).
    # X = variable óptima reco del canal (ω reco para ρ; x reco para π/lep); Y = cosθ visible reco.
    ("Reco",    "OptimalReco_vs_CosThetaVis",
                                       lambda v, sh: v.get("_optimal_unified_reco", -999.0),
                                       lambda v, sh: _reco_tau_sign(v) * math.cos(v["recoVisTheta"])),
    # Variable óptima del ρ vía red neuronal (modelo v2), salida ∈[0,1]. Solo se
    # llena en el hemisferio ρ de pares ρ-leptón (con --use-nn-optimal); en el
    # resto vale -999. NO se usa para repesado; los histogramas de ω no cambian.
    ("Reco",    "OptimalNN_Reco",     lambda v, sh: v.get("_optimal_nn", -999.0),  None),
    ("Reco",    "OptimalNN_vs_CosThetaVis",
                                       lambda v, sh: v.get("_optimal_nn", -999.0),
                                       lambda v, sh: _reco_tau_sign(v) * math.cos(v["recoVisTheta"])),
    # ── Matched ───────────────────────────────────────────────────────────────
    ("Matched", "VisEOverBeamE",      lambda v, sh: v["visE"] / sh["beamE"]
                                                    if sh["beamE"] else 0.0,      None),
    ("Matched", "VisCosTheta",        lambda v, sh: math.cos(v["visTheta"]),      None),

    # ── Separación por carga τ⁺/τ⁻ de los observables de polarización ──────────
    # Mismo x que la variable inclusiva, con cond_fn extra sobre tauPDG.
    # Para los Omega_RhoMass_Bin* se mantiene además el corte de masa visM.
    ("Gen",  "Omega_plus",       lambda v, sh: v["omega"],      None, lambda v, sh: _is_taup(v)),
    ("Gen",  "Omega_minus",      lambda v, sh: v["omega"],      None, lambda v, sh: _is_taum(v)),
    ("Gen",  "OptimalVar_plus",  lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taup(v)),
    ("Gen",  "OptimalVar_minus", lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taum(v)),
    ("Gen",  "Omega_RhoMass_Bin1_plus",  lambda v, sh: v["omega"], None, lambda v, sh: 0.40 <= v["visM"] < 0.60 and _is_taup(v)),
    ("Gen",  "Omega_RhoMass_Bin1_minus", lambda v, sh: v["omega"], None, lambda v, sh: 0.40 <= v["visM"] < 0.60 and _is_taum(v)),
    ("Gen",  "Omega_RhoMass_Bin2_plus",  lambda v, sh: v["omega"], None, lambda v, sh: 0.60 <= v["visM"] < 0.70 and _is_taup(v)),
    ("Gen",  "Omega_RhoMass_Bin2_minus", lambda v, sh: v["omega"], None, lambda v, sh: 0.60 <= v["visM"] < 0.70 and _is_taum(v)),
    ("Gen",  "Omega_RhoMass_Bin3_plus",  lambda v, sh: v["omega"], None, lambda v, sh: 0.70 <= v["visM"] < 0.80 and _is_taup(v)),
    ("Gen",  "Omega_RhoMass_Bin3_minus", lambda v, sh: v["omega"], None, lambda v, sh: 0.70 <= v["visM"] < 0.80 and _is_taum(v)),
    ("Gen",  "Omega_RhoMass_Bin4_plus",  lambda v, sh: v["omega"], None, lambda v, sh: 0.80 <= v["visM"] < 0.90 and _is_taup(v)),
    ("Gen",  "Omega_RhoMass_Bin4_minus", lambda v, sh: v["omega"], None, lambda v, sh: 0.80 <= v["visM"] < 0.90 and _is_taum(v)),
    ("Gen",  "Omega_RhoMass_Bin5_plus",  lambda v, sh: v["omega"], None, lambda v, sh: 0.90 <= v["visM"] < 1.05 and _is_taup(v)),
    ("Gen",  "Omega_RhoMass_Bin5_minus", lambda v, sh: v["omega"], None, lambda v, sh: 0.90 <= v["visM"] < 1.05 and _is_taum(v)),
    ("Gen",  "Omega_RhoMass_Bin6_plus",  lambda v, sh: v["omega"], None, lambda v, sh: 1.05 <= v["visM"] < 1.50 and _is_taup(v)),
    ("Gen",  "Omega_RhoMass_Bin6_minus", lambda v, sh: v["omega"], None, lambda v, sh: 1.05 <= v["visM"] < 1.50 and _is_taum(v)),
    # Helicidad por carga: τ⁻ hel=-1 (tipo M1), τ⁻ hel=+1 (tipo P1); τ⁺ invertido.
    ("Gen",  "Omega_minus_NegHel", lambda v, sh: v["omega"], None, lambda v, sh: _is_taum(v) and v.get("genHelicity", 0.0) < 0),
    ("Gen",  "Omega_minus_PosHel", lambda v, sh: v["omega"], None, lambda v, sh: _is_taum(v) and v.get("genHelicity", 0.0) > 0),
    ("Gen",  "Omega_plus_NegHel",  lambda v, sh: v["omega"], None, lambda v, sh: _is_taup(v) and v.get("genHelicity", 0.0) < 0),
    ("Gen",  "Omega_plus_PosHel",  lambda v, sh: v["omega"], None, lambda v, sh: _is_taup(v) and v.get("genHelicity", 0.0) > 0),
    # ("Gen", "OptimalVar_PionMCut", lambda v, sh: v["optimalVar"], None, lambda v, sh: v.get("pionM", 0.0) < 0.2),
    # ("Gen", "OptimalVar_minus_NegHel_pionMCut", lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taum(v) and v.get("genHelicity", 0.0) < 0 and v.get("pionM", 0.0) < 0.2),
    # ("Gen", "OptimalVar_minus_PosHel_pionMCut", lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taum(v) and v.get("genHelicity", 0.0) > 0 and v.get("pionM", 0.0) < 0.2),    
    # ("Gen", "OptimalVar_plus_NegHel_pionMCut",  lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taup(v) and v.get("genHelicity", 0.0) < 0 and v.get("pionM", 0.0) < 0.2),
    # ("Gen", "OptimalVar_plus_PosHel_pionMCut",  lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taup(v) and v.get("genHelicity", 0.0) > 0 and v.get("pionM", 0.0) < 0.2),
    # ("Gen", "OptimalVar_PosHel_pionMCut",  lambda v, sh: v["optimalVar"], None, lambda v, sh: ((_is_taup(v) and v.get("genHelicity", 0.0) > 0) or (_is_taum(v) and v.get("genHelicity", 0.0) < 0)) and v.get("pionM", 0.0) < 0.2),
    # ("Gen", "OptimalVar_NegHel_pionMCut",  lambda v, sh: v["optimalVar"], None, lambda v, sh: ((_is_taup(v) and v.get("genHelicity", 0.0) < 0) or (_is_taum(v) and v.get("genHelicity", 0.0) > 0)) and v.get("pionM", 0.0) < 0.2),
    ("Gen", "OptimalVar_PosHel",  lambda v, sh: v["optimalVar"], None, lambda v, sh: ((_is_taup(v) and v.get("genHelicity", 0.0) > 0) or (_is_taum(v) and v.get("genHelicity", 0.0) < 0))),
    ("Gen", "OptimalVar_NegHel",  lambda v, sh: v["optimalVar"], None, lambda v, sh: ((_is_taup(v) and v.get("genHelicity", 0.0) < 0) or (_is_taum(v) and v.get("genHelicity", 0.0) > 0))),
    ("Gen",  "OptimalVar_minus_NegHel", lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taum(v) and v.get("genHelicity", 0.0) < 0),
    ("Gen",  "OptimalVar_minus_PosHel", lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taum(v) and v.get("genHelicity", 0.0) > 0),
    ("Gen",  "OptimalVar_plus_NegHel",  lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taup(v) and v.get("genHelicity", 0.0) < 0),
    ("Gen",  "OptimalVar_plus_PosHel",  lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_taup(v) and v.get("genHelicity", 0.0) > 0),
    # Helicidad COMBINADA en carga (convención Cepeda, basada en la helicidad del τ⁻).
    ("Gen",  "Omega_HelNeg",      lambda v, sh: v["omega"],      None, lambda v, sh: _is_helneg(v)),
    ("Gen",  "Omega_HelPos",      lambda v, sh: v["omega"],      None, lambda v, sh: _is_helpos(v)),
    ("Gen",  "OptimalVar_HelNeg", lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_helneg(v)),
    ("Gen",  "OptimalVar_HelPos", lambda v, sh: v["optimalVar"], None, lambda v, sh: _is_helpos(v)),
    ("Reco", "Optimal_X_plus",   lambda v, sh: v.get("_optimal_x", 0.0), None, lambda v, sh: _is_taup(v)),
    ("Reco", "Optimal_X_minus",  lambda v, sh: v.get("_optimal_x", 0.0), None, lambda v, sh: _is_taum(v)),
    # Omega_Reco en 6 bins de masa reco del rho (recoVisM) — mismo binado que gen, pero reco
    ("Reco", "Omega_Reco_RhoMass_Bin1", lambda v, sh: v.get("_omega_reco", -999.0), None,
                                         lambda v, sh: 0.40 <= v["recoVisM"] < 0.60),
    ("Reco", "Omega_Reco_RhoMass_Bin2", lambda v, sh: v.get("_omega_reco", -999.0), None,
                                         lambda v, sh: 0.60 <= v["recoVisM"] < 0.70),
    ("Reco", "Omega_Reco_RhoMass_Bin3", lambda v, sh: v.get("_omega_reco", -999.0), None,
                                         lambda v, sh: 0.70 <= v["recoVisM"] < 0.80),
    ("Reco", "Omega_Reco_RhoMass_Bin4", lambda v, sh: v.get("_omega_reco", -999.0), None,
                                         lambda v, sh: 0.80 <= v["recoVisM"] < 0.90),
    ("Reco", "Omega_Reco_RhoMass_Bin5", lambda v, sh: v.get("_omega_reco", -999.0), None,
                                         lambda v, sh: 0.90 <= v["recoVisM"] < 1.05),
    ("Reco", "Omega_Reco_RhoMass_Bin6", lambda v, sh: v.get("_omega_reco", -999.0), None,
                                         lambda v, sh: 1.05 <= v["recoVisM"] < 1.50),
]

# Reglas para histogramas compartidos (no por hemisferio).
# Formato: (level, hist_name, x_fn(v0, v1, sh), y_fn(v0, v1, sh) | None)
# v0 = vars_dec0, v1 = vars_dec1, sh = shared_vars
# Categoría se calcula desde la perspectiva de dec0.

FILL_RULES_SHARED = [
    ("Gen",  "GenZMass",               lambda v0, v1, sh: sh["GenZMass"],    None),
    ("Gen",  "GenVisZMass",            lambda v0, v1, sh: sh["GenZVisMass"], None),
    ("Reco", "RecoZMass",              lambda v0, v1, sh: sh["ZRecoMass"],   None),
    ("Reco", "DeltaR_dec0_dec1",       None,                                 None),  # calculado inline
    # Cross-hemisphere 2D
    ("Gen",  "VisP_dec0_vs_dec1",      lambda v0, v1, sh: v0["visP"],        lambda v0, v1, sh: v1["visP"]),
    ("Gen",  "Omega_dec0_vs_VisP_dec1",lambda v0, v1, sh: v0["omega"],       lambda v0, v1, sh: v1["visP"]),
    ("Gen",  "Omega_dec0_vs_Omega_dec1",  lambda v0, v1, sh: v0["omega"],        lambda v0, v1, sh: v1["omega"]),
    ("Gen",  "WeightP1_dec0_vs_dec1",    lambda v0, v1, sh: v0["weight_P1"],    lambda v0, v1, sh: v1["weight_P1"]),
    ("Gen",  "WeightM1_dec0_vs_dec1",    lambda v0, v1, sh: v0["weight_M1"],    lambda v0, v1, sh: v1["weight_M1"]),
]


# ── Funciones auxiliares ───────────────────────────────────────────────────────

def write_histograms_recursive(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            write_histograms_recursive(v)
    else:
        try:
            obj.Write()
        except AttributeError:
            pass


def split_entry_ranges(n_entries, n_workers):
    k, rem = divmod(n_entries, n_workers)
    ranges, start = [], 0
    for i in range(n_workers):
        end = start + k + (1 if i < rem else 0)
        if start < end:
            ranges.append((start, end))
        start = end
    return ranges


def _reattach_histograms(nested, tfile):
    if isinstance(nested, dict):
        for v in nested.values():
            _reattach_histograms(v, tfile)
    elif isinstance(nested, ROOT.TH1):
        nested.SetDirectory(tfile)


def _flatten_histograms(nested, result=None):
    if result is None:
        result = {}
    if isinstance(nested, dict):
        for v in nested.values():
            _flatten_histograms(v, result)
    else:
        try:
            result[nested.GetName()] = nested
        except AttributeError:
            pass
    return result


def _assign_single_decay(tau1_vars, tau2_vars, id_gen, id_reco=None, only_gen=False):
    """Asigna un único tau a 'este' y el otro a 'otro' según el canal objetivo.

    only_gen=False (reco real): selecciona por id RECO exacto (recoTauID == id_reco;
        1 y 2 son canales distintos, sin remap), como analysisRHOTree.py
        (`recoTauID==selectDecay`). La señal/fondo se decide luego por gen.
    only_gen=True (árbol gen-only): selecciona por decayID (id gen, id_gen ya con
        remap 2→1).

    Si tau1 coincide: devuelve (tau1, tau2).
    Si solo tau2 coincide: devuelve (tau2, tau1).
    Si ambos coinciden: devuelve (tau1, tau2) — una sola vez, por convención.
    Si ninguno coincide: devuelve (None, None).
    """
    if only_gen:
        target, key = id_gen, "decayID"
    else:
        target, key = id_reco, "recoTauID"
    t1 = int(tau1_vars.get(key, -999))
    t2 = int(tau2_vars.get(key, -999))
    if t1 == target:
        return tau1_vars, tau2_vars
    if t2 == target:
        return tau2_vars, tau1_vars
    return None, None


def _recompute_weights(tau_vars, beamE, sin_eff, use_omega=False, use_costheta_pion=True):
    """
    Recalcula weight_P1/weight_M1 con un sin²θ_eff dado.
    Usa cinemática gen exacta (P, Theta, Phi, E del tau gen).
    Para rho con use_omega=True usa newAtauFromH con la variable óptima ω gen almacenada.
    Para pión con use_costheta_pion=True usa cos(θ*) geométrico recalculado al vuelo
    (boost exacto desde las cinemáticas gen — funciona también con árboles antiguos);
    en caso contrario usa newAtau con H_V/z_R (Alcaraz 2026 eq. 4-5).

    With no matched gen tau (tauPDG = -999) there is nothing to recompute: the
    weights are left exactly as the producer wrote them, rather than inventing a
    sign for z.
    """
    if not _has_gen_tau(tau_vars):
        return
    decay_id = int(tau_vars.get("decayID", -999))
    tau_pdg  = _gen_tau_pdg(tau_vars)
    if decay_id in (0, 1, 10):  # hadrónico (pion, rho, a1)
        tauP4 = make_p4(tau_vars["P"],    tau_vars["Theta"],    tau_vars["Phi"],    tau_vars["E"])
        visP4 = make_p4(tau_vars["visP"], tau_vars["visTheta"], tau_vars["visPhi"], tau_vars["visE"])
        if use_omega and decay_id == 1:
            omega = tau_vars.get("omega", -999.0)
            w_P1  = weightsPol.newAtauFromH(tauP4, omega, +1, tau_pdg=tau_pdg, sin_eff=sin_eff)
            w_M1  = weightsPol.newAtauFromH(tauP4, omega, -1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        elif use_costheta_pion and decay_id == 0:
            cts   = weightsPol.cosThetaStar(tauP4, visP4)
            w_P1  = weightsPol.newAtauFromH(tauP4, cts, +1, tau_pdg=tau_pdg, sin_eff=sin_eff)
            w_M1  = weightsPol.newAtauFromH(tauP4, cts, -1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        else:
            w_P1  = weightsPol.newAtau(tauP4, visP4, decay_id, +1, tau_pdg=tau_pdg, sin_eff=sin_eff)
            w_M1  = weightsPol.newAtau(tauP4, visP4, decay_id, -1, tau_pdg=tau_pdg, sin_eff=sin_eff)
    elif decay_id in (-11, -13):  # leptónico
        tauP4 = make_p4(tau_vars["P"],    tau_vars["Theta"],    tau_vars["Phi"],    tau_vars["E"])
        visP4 = make_p4(tau_vars["visP"], tau_vars["visTheta"], tau_vars["visPhi"], tau_vars["visE"])
        w_P1  = weightsPol.newAtauLep(visP4, tauP4, beamE, +1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        w_M1  = weightsPol.newAtauLep(visP4, tauP4, beamE, -1, tau_pdg=tau_pdg, sin_eff=sin_eff)
    else:
        w_P1 = w_M1 = 1.0
    tau_vars["weight_P1"] = w_P1
    tau_vars["weight_M1"] = w_M1


def _recompute_optimal_var(tau_vars, beamE):
    """Recalcula optimalVar (gen) al vuelo con la definición ÚNICA del helper compartido
    (optimalVariabRho.optimal_var): π/a1/lep = E_vis/E_τ, ρ = ω.

    Permite reutilizar árboles antiguos con la definición actual del observable (p.ej.
    el cambio del pión de E_π a E_vis). Puramente cinemático → no depende de sin_eff;
    para ρ/lep coincide con lo almacenado, para π corrige árboles previos.

    With no matched gen tau (tauPDG = -999) the gen observable does not exist: -999.
    """
    if not _has_gen_tau(tau_vars):
        tau_vars["optimalVar"] = -999.0
        return
    decay_id = int(tau_vars.get("decayID", -999))
    tauP4  = make_p4(tau_vars["P"],     tau_vars["Theta"],     tau_vars["Phi"],     tau_vars["E"])
    visP4  = make_p4(tau_vars["visP"],  tau_vars["visTheta"],  tau_vars["visPhi"],  tau_vars["visE"])
    pionP4 = make_p4(tau_vars["pionP"], tau_vars["pionTheta"], tau_vars["pionPhi"], tau_vars["pionE"])
    tau_vars["optimalVar"] = optimalVariabRho.optimal_var(
        decay_id, tauP4, visP4, pionP4, beamE, tau_pdg=_gen_tau_pdg(tau_vars))


def _reco_htype(reco_id):
    """Tipo hadrónico para weightsPol (espera 0/1/10): ρ reco (2) → 1.

    A nivel RECO la ρ tiene recoTauID==2 (la gen usa decayID==1). weightsPol
    (_compute_H, newAtau) solo entiende 0/1/10, así que hay que remapear 2→1.
    """
    return 1 if reco_id == 2 else reco_id


def _recompute_reco_weights(tau_vars, beamE, sin_eff, use_omega=False):
    """
    Calcula reco_weight_P1/reco_weight_M1 con cinemática reco.

    El tau no es observable directamente (falta el neutrino). Se usa la
    aproximación collinear: tau proxy con E=beamE y dirección del visible reco.
    Así x = E_vis_reco / E_beam (fracción de energía correcta para las fórmulas).

    Para leptónico se usa el leptón reco como visible y el proxy como tau.
    Con use_omega=True y rho reco: usa wVariabRECO para obtener ω y lo pasa a
    newAtauRhoOmega en lugar de H_V = alpha_V * z_R.
    """
    reco_id  = int(tau_vars.get("recoTauID", -999))
    tau_pdg  = _reco_tau_pdg(tau_vars)
    if beamE <= 0:
        tau_vars["reco_weight_P1"] = 1.0
        tau_vars["reco_weight_M1"] = 1.0
        return

    # Tau proxy: misma dirección que el visible reco, E = beamE (collinear approx)
    tau_E = beamE
    tau_P = math.sqrt(max(0.0, tau_E**2 - weightsPol._M_TAU**2))
    tau_proxy_P4 = make_p4(tau_P, tau_vars["recoVisTheta"],
                            tau_vars["recoVisPhi"], tau_E)

    if reco_id in (0, 1, 2, 10):  # hadrónico reco (ρ reco = 2)
        vis_P4 = make_p4(tau_vars["recoVisP"], tau_vars["recoVisTheta"],
                          tau_vars["recoVisPhi"], tau_vars["recoVisE"])
        if use_omega and reco_id in (1, 2):
            pion_P4 = make_p4(tau_vars["recoPionP"], tau_vars["recoPionTheta"],
                               tau_vars["recoPionPhi"], tau_vars["recoPionE"])
            _, _, _, omega_reco = optimalVariabRho.wVariabRECO(vis_P4, pion_P4, beamE)
            w_P1 = weightsPol.newAtauRhoOmega(tau_proxy_P4, omega_reco, +1,
                                              tau_pdg=tau_pdg, sin_eff=sin_eff)
            w_M1 = weightsPol.newAtauRhoOmega(tau_proxy_P4, omega_reco, -1,
                                              tau_pdg=tau_pdg, sin_eff=sin_eff)
        else:
            htype = _reco_htype(reco_id)
            w_P1 = weightsPol.newAtau(tau_proxy_P4, vis_P4, htype, +1,
                                      tau_pdg=tau_pdg, sin_eff=sin_eff)
            w_M1 = weightsPol.newAtau(tau_proxy_P4, vis_P4, htype, -1,
                                      tau_pdg=tau_pdg, sin_eff=sin_eff)
    elif reco_id in (-11, -13):  # leptónico reco
        lep_P4 = make_p4(tau_vars["recoLepP"], tau_vars["recoLepTheta"],
                          tau_vars["recoLepPhi"], tau_vars["recoLepE"])
        w_P1 = weightsPol.newAtauLep(lep_P4, tau_proxy_P4, beamE, +1,
                                     tau_pdg=tau_pdg, sin_eff=sin_eff)
        w_M1 = weightsPol.newAtauLep(lep_P4, tau_proxy_P4, beamE, -1,
                                     tau_pdg=tau_pdg, sin_eff=sin_eff)
    else:
        w_P1 = w_M1 = 1.0
    tau_vars["reco_weight_P1"] = w_P1
    tau_vars["reco_weight_M1"] = w_M1


def _get_H_for_joint(tau_vars, beamE, use_omega, use_costheta_pion=True):
    """
    Devuelve el observable de spin H a usar en la fórmula joint.
    - Rho con use_omega=True: usa omega almacenado en el árbol (variable óptima completa).
    - Pión con use_costheta_pion=True: usa cos(θ*) geométrico (boost exacto al vuelo).
    - Rho con use_omega=False, pión con flag off, a1: usa H_V/z_R de _compute_H.
    - Leptónico: usa H_ell (no depende de los flags).
    """
    decay_id = int(tau_vars.get("decayID", -999))
    if use_omega and decay_id == 1:
        return tau_vars.get("omega", 0.0)
    if decay_id in (0, 1, 10):
        tauP4 = make_p4(tau_vars["P"],    tau_vars["Theta"], tau_vars["Phi"],    tau_vars["E"])
        visP4 = make_p4(tau_vars["visP"], tau_vars["visTheta"], tau_vars["visPhi"], tau_vars["visE"])
        if use_costheta_pion and decay_id == 0:
            return weightsPol.cosThetaStar(tauP4, visP4)
        H = weightsPol._compute_H(visP4, tauP4, decay_id)
        return H if H is not None else 0.0
    if decay_id in (-11, -13):
        visP4 = make_p4(tau_vars["visP"], tau_vars["visTheta"], tau_vars["visPhi"], tau_vars["visE"])
        return weightsPol._compute_H_lep(visP4, beamE)
    return 0.0


def _compute_joint_weights(vars_dec0, vars_dec1, beamE, sin_eff, use_omega=False,
                           use_costheta_pion=True):
    """
    Peso conjunto dos-tau (Alcaraz 2026 eqs. 9, 13 y 16) con cinemática gen.
    Incluye el término cruzado H·H' ausente en el producto independiente.
    Almacena weight_corr_P1/M1 (idéntico en ambos hemisferios).

    Se aplica a TODAS las combinaciones con observable definido: had-had,
    had-lep y lep-lep (eq. 16 general). El producto solo se usa como último
    recurso para decays no soportados.

    Con use_omega=True: usa omega almacenado en el árbol para taus rho,
    H_V/H_ell para los demás. La fórmula siempre es:
      W = [1 + P_new*(H + H') + H*H'] / [1 + P_sm*(H + H') + H*H']
    """
    id0 = int(vars_dec0.get("decayID", -999))
    id1 = int(vars_dec1.get("decayID", -999))
    is_had = lambda d: d in (0, 1, 10)
    is_lep = lambda d: d in (-11, -13)
    # El peso joint necesita z = cos(θ_τ⁻), es decir el gen-tau de al menos un
    # hemisferio. Sin él (fondos no-ττ: tauPDG = -999) se cae al producto de los
    # pesos por-tau, igual que con los decays no soportados.
    has_gen = _has_gen_tau(vars_dec0) and _has_gen_tau(vars_dec1)

    def p4(v):
        return make_p4(v["P"], v["Theta"], v["Phi"], v["E"])

    for New_Atau, suffix in [(+1.0, "P1"), (-1.0, "M1")]:
        if not has_gen:
            w = (vars_dec0.get(f"weight_{suffix}", 1.0) *
                 vars_dec1.get(f"weight_{suffix}", 1.0))
        elif is_had(id0) and (is_had(id1) or is_lep(id1)):
            H  = _get_H_for_joint(vars_dec0, beamE, use_omega, use_costheta_pion)
            Hp = _get_H_for_joint(vars_dec1, beamE, use_omega, use_costheta_pion)
            w  = weightsPol.newAtauJoint(p4(vars_dec0), H, Hp, New_Atau,
                                         tau_pdg=_gen_tau_pdg(vars_dec0), sin_eff=sin_eff)
        elif is_lep(id0) and is_had(id1):
            H  = _get_H_for_joint(vars_dec0, beamE, use_omega, use_costheta_pion)
            Hp = _get_H_for_joint(vars_dec1, beamE, use_omega, use_costheta_pion)
            w  = weightsPol.newAtauJoint(p4(vars_dec1), Hp, H, New_Atau,
                                         tau_pdg=_gen_tau_pdg(vars_dec1), sin_eff=sin_eff)
        elif is_lep(id0) and is_lep(id1):
            # lep-lep: fórmula joint general (Alcaraz eq. 16) con H = H_ell para ambos
            # τ. El término cruzado H·H' (correlación de espín) es físico también aquí;
            # el producto independiente sería incorrecto.
            H  = _get_H_for_joint(vars_dec0, beamE, use_omega, use_costheta_pion)
            Hp = _get_H_for_joint(vars_dec1, beamE, use_omega, use_costheta_pion)
            w  = weightsPol.newAtauJoint(p4(vars_dec0), H, Hp, New_Atau,
                                         tau_pdg=_gen_tau_pdg(vars_dec0), sin_eff=sin_eff)
        else:
            # Solo decays no soportados (p.ej. id -2): producto como último recurso.
            w = (vars_dec0.get(f"weight_{suffix}", 1.0) *
                 vars_dec1.get(f"weight_{suffix}", 1.0))
        vars_dec0[f"weight_corr_{suffix}"] = w
        vars_dec1[f"weight_corr_{suffix}"] = w


def _reco_tau_proxy(tau_vars, beamE):
    """Tau proxy reco: dirección del visible reco, E=beamE (collinear approx)."""
    tau_E = beamE
    tau_P = math.sqrt(max(0.0, tau_E**2 - weightsPol._M_TAU**2))
    return make_p4(tau_P, tau_vars["recoVisTheta"], tau_vars["recoVisPhi"], tau_E)


def _compute_omega_reco(tau_vars, beamE):
    """
    ω reco = variable óptima reconstruida del canal ρ vía wVariabRECO
    (misma fórmula angular que gen pero con θ_ρ analítico y 4-vectores reco).
    Solo definida para ρ reco (recoTauID 1 ó 2: ρ con un fotón perdido vs ρ
    completo — ambos vienen de un ρ, se calcula la var óptima igual); para los
    demás canales devuelve -999.0, igual que la convención de la rama gen `omega`.
    """
    if beamE <= 0 or int(tau_vars.get("recoTauID", -999)) not in (1, 2):
        return -999.0
    vis_P4  = make_p4(tau_vars["recoVisP"],  tau_vars["recoVisTheta"],
                      tau_vars["recoVisPhi"], tau_vars["recoVisE"])
    pion_P4 = make_p4(tau_vars["recoPionP"], tau_vars["recoPionTheta"],
                      tau_vars["recoPionPhi"], tau_vars["recoPionE"])
    ct_reco, cp_reco, _, omega_reco = optimalVariabRho.wVariabRECO(vis_P4, pion_P4, beamE)
    # Guardar cos_theta/cos_psi reco para el corte de bordes de ω (legacy: abs(==1) =
    # caso clampeado, cinemática fuera del rango τ→ρ permitido → enriquecido en fakes).
    tau_vars["_ct_reco"] = ct_reco
    tau_vars["_cp_reco"] = cp_reco
    return omega_reco


def _read_reco_photons(entry, prefix):
    """Lee los vectores de fotones reco de un hemisferio como lista de (E, theta, phi).

    Los vectores no están en _TAU_KEYS_MDECS (que son escalares); se leen aparte
    solo cuando hace falta (path NN). Devuelve [] si las ramas no existen.
    """
    try:
        Es     = getattr(entry, f"{prefix}_reco_photons_E")
        thetas = getattr(entry, f"{prefix}_reco_photons_theta")
        phis   = getattr(entry, f"{prefix}_reco_photons_phi")
    except AttributeError:
        return []
    return [(float(Es[j]), float(thetas[j]), float(phis[j])) for j in range(len(Es))]


def _compute_nn_optimal(vars_rho, vars_lep, nn_model):
    """Variable óptima reco del ρ vía la red (modelo v2): salida sigmoide ∈[0,1].

    Combina el sistema visible (meson) y los 2 fotones de mayor E del hemisferio ρ
    con el leptón del hemisferio opuesto. Devuelve -999.0 si no hay ≥2 fotones reco.
    Solo se usa como variable óptima (NO para repesado).
    """
    from modules import mlpPolInference
    photons = vars_rho.get("_reco_photons", [])
    if len(photons) < 2:
        return -999.0
    # Modelo por defecto (optimal_lepton): primeras 4 features = PIÓN cargado reco.
    pion = (vars_rho["recoPionE"], vars_rho["recoPionTheta"],
            vars_rho["recoPionPhi"], vars_rho["recoPionP"])
    lep = (vars_lep["recoLepE"], vars_lep["recoLepTheta"],
           vars_lep["recoLepPhi"], vars_lep["recoLepP"])
    feats = mlpPolInference.build_features(pion, photons, lep)
    return float(nn_model.predict(feats))


def _maybe_set_nn_optimal(vars_a, vars_b, nn_model):
    """Asigna `_optimal_nn` al hemisferio ρ del par (ρ,leptón); -999 al leptón.

    Acepta el par en cualquier orden. Si no es un par ρ-leptón no hace nada.
    """
    def _is_rho(v):  return int(v.get("recoTauID", -999)) in (1, 2)
    def _is_lep(v):  return int(v.get("recoTauID", -999)) in (-11, -13)

    vars_a["_optimal_nn"] = -999.0
    vars_b["_optimal_nn"] = -999.0
    if _is_rho(vars_a) and _is_lep(vars_b):
        vars_a["_optimal_nn"] = _compute_nn_optimal(vars_a, vars_b, nn_model)
    elif _is_lep(vars_a) and _is_rho(vars_b):
        vars_b["_optimal_nn"] = _compute_nn_optimal(vars_b, vars_a, nn_model)


def _get_H_for_joint_reco(tau_vars, beamE, use_omega):
    """
    Observable de spin H a usar en la fórmula joint con cinemática reco.
    Análogo a _get_H_for_joint pero con ramas reco y el tau proxy collinear.
    - Rho con use_omega=True: ω reco vía wVariabRECO (variable óptima completa).
    - Rho con use_omega=False, pion, a1: H_V = alpha_V*z_R de _compute_H.
    - Leptónico: H_ell = variable óptima leptónica (no hay otra; x_ell es toda la info).
    """
    reco_id = int(tau_vars.get("recoTauID", -999))
    if use_omega and reco_id in (1, 2):
        vis_P4  = make_p4(tau_vars["recoVisP"],  tau_vars["recoVisTheta"],
                          tau_vars["recoVisPhi"], tau_vars["recoVisE"])
        pion_P4 = make_p4(tau_vars["recoPionP"], tau_vars["recoPionTheta"],
                          tau_vars["recoPionPhi"], tau_vars["recoPionE"])
        _, _, _, omega_reco = optimalVariabRho.wVariabRECO(vis_P4, pion_P4, beamE)
        return omega_reco
    if reco_id in (0, 1, 2, 10):
        vis_P4 = make_p4(tau_vars["recoVisP"], tau_vars["recoVisTheta"],
                         tau_vars["recoVisPhi"], tau_vars["recoVisE"])
        H = weightsPol._compute_H(vis_P4, _reco_tau_proxy(tau_vars, beamE), _reco_htype(reco_id))
        return H if H is not None else 0.0
    if reco_id in (-11, -13):
        lep_P4 = make_p4(tau_vars["recoLepP"], tau_vars["recoLepTheta"],
                         tau_vars["recoLepPhi"], tau_vars["recoLepE"])
        return weightsPol._compute_H_lep(lep_P4, beamE)
    return 0.0


def _compute_reco_joint_weights(vars_dec0, vars_dec1, beamE, sin_eff, use_omega=False):
    """
    Peso conjunto dos-tau con cinemática reco (proxy collinear).
    Misma lógica que _compute_joint_weights pero usando el tau proxy
    (E=beamE, dirección=recoVisTheta/Phi) y observables reco.
    Con use_omega=True usa ω reco (wVariabRECO) para taus rho, consistente
    con _recompute_reco_weights. Almacena reco_weight_corr_P1/M1.
    """
    reco_id0 = int(vars_dec0.get("recoTauID", -999))
    reco_id1 = int(vars_dec1.get("recoTauID", -999))
    is_had = lambda d: d in (0, 1, 2, 10)  # ρ reco = 2
    is_lep = lambda d: d in (-11, -13)

    if beamE <= 0:
        for suffix in ("P1", "M1"):
            vars_dec0[f"reco_weight_corr_{suffix}"] = 1.0
            vars_dec1[f"reco_weight_corr_{suffix}"] = 1.0
        return

    for New_Atau, suffix in [(+1.0, "P1"), (-1.0, "M1")]:
        if is_had(reco_id0) and (is_had(reco_id1) or is_lep(reco_id1)):
            H  = _get_H_for_joint_reco(vars_dec0, beamE, use_omega)
            Hp = _get_H_for_joint_reco(vars_dec1, beamE, use_omega)
            w  = weightsPol.newAtauJoint(_reco_tau_proxy(vars_dec0, beamE), H, Hp,
                                         New_Atau, tau_pdg=_reco_tau_pdg(vars_dec0),
                                         sin_eff=sin_eff)
        elif is_lep(reco_id0) and is_had(reco_id1):
            H  = _get_H_for_joint_reco(vars_dec0, beamE, use_omega)
            Hp = _get_H_for_joint_reco(vars_dec1, beamE, use_omega)
            w  = weightsPol.newAtauJoint(_reco_tau_proxy(vars_dec1, beamE), Hp, H,
                                         New_Atau, tau_pdg=_reco_tau_pdg(vars_dec1),
                                         sin_eff=sin_eff)
        elif is_lep(reco_id0) and is_lep(reco_id1):
            # lep-lep reco: fórmula joint general (eq. 16), H_ell para ambos.
            H  = _get_H_for_joint_reco(vars_dec0, beamE, use_omega)
            Hp = _get_H_for_joint_reco(vars_dec1, beamE, use_omega)
            w  = weightsPol.newAtauJoint(_reco_tau_proxy(vars_dec0, beamE), H, Hp,
                                         New_Atau, tau_pdg=_reco_tau_pdg(vars_dec0),
                                         sin_eff=sin_eff)
        else:
            # Solo decays no soportados: producto como último recurso.
            w = (vars_dec0.get(f"reco_weight_{suffix}", 1.0) *
                 vars_dec1.get(f"reco_weight_{suffix}", 1.0))
        vars_dec0[f"reco_weight_corr_{suffix}"] = w
        vars_dec1[f"reco_weight_corr_{suffix}"] = w


def _assign_hemispheres(tau1_vars, tau2_vars, id0_gen, id1_gen,
                        reco_id0=None, reco_id1=None, only_gen=False):
    """Asigna tau1/tau2 a dec0/dec1 (order-independent).

    only_gen=False (reco real): selecciona/empareja por id RECO exacto
        (recoTauID == reco_id0/reco_id1; 1 y 2 son canales distintos, sin remap),
        como analysisRHOTree.py (`recoTauID==selectDecay`). La señal/fondo se
        decide después por verdad gen.
    only_gen=True (árbol gen-only): empareja por decayID (ids gen, con remap 2→1
        ya aplicado en id0_gen/id1_gen).

    Devuelve (vars_dec0, vars_dec1) o (None, None) si el evento no encaja con
    el par esperado. Para par simétrico, asigna tau1→dec0, tau2→dec1.
    """
    if only_gen:
        a0, a1, key = id0_gen, id1_gen, "decayID"
    else:
        a0, a1, key = reco_id0, reco_id1, "recoTauID"
    t1_id = int(tau1_vars.get(key, -999))
    t2_id = int(tau2_vars.get(key, -999))

    if t1_id == a0 and t2_id == a1:
        return tau1_vars, tau2_vars
    if t2_id == a0 and t1_id == a1:
        return tau2_vars, tau1_vars
    return None, None


def _classify_hemisphere(dec_vars, expected_gen_id, use_reco=False,
                          expected_reco_id=None, only_gen=False):
    """Devuelve la categoría de un hemisferio.

    Returns: una de 'SIGNAL', 'BGMuon', 'BGEle', 'BGPion', 'BGRho', 'BGA1', 'BGOther'
    Si expected_gen_id es None, siempre devuelve la categoría BG real (modo single-decay).

    only_gen=True  → el árbol es gen-only (recoTauID contiene IDs gen); se remap 2→1
                     tanto en la rama gen como reco.
    only_gen=False → árbol reco real. El evento YA está seleccionado por id reco
                     aguas arriba (_assign_hemispheres / _assign_single_decay), donde
                     1 y 2 son canales distintos. Aquí la SEÑAL/FONDO se decide por
                     la VERDAD GEN (decayID), igual que analysisRHOTree.py
                     (`genTauID==selectGEN`): SIGNAL si el gen coincide con el canal,
                     y el desglose BG{Muon,Ele,Pion,Rho,A1,Other} es por decayID gen.
                     Todos los fondos así definidos están reconstruidos como el canal
                     (p.ej. ρ) → tienen variable óptima/ω definida.
    """
    if use_reco and not only_gen:
        # Reco real: señal/fondo por VERDAD GEN (decayID), no por id reco.
        actual = int(dec_vars.get("decayID", -999))
        if expected_gen_id is not None and actual == expected_gen_id:
            return "SIGNAL"
        return _BG_MAP.get(actual, "BGOther")
    else:
        # Gen path, o reco de árbol gen-only (only_gen=True).
        # Aquí no se remapea el decayID gen real: si el árbol contiene 2,
        # ese modo se clasifica como fondo y no como rho señal.
        actual = int(dec_vars.get("recoTauID" if use_reco else "decayID", -999))
        if only_gen and use_reco and actual == 2:
            actual = 1
        if expected_gen_id is not None and actual == expected_gen_id:
            return "SIGNAL"
        return _BG_MAP.get(actual, "BGOther")


def _get_fill_categories(this_cat, other_cat):
    """Devuelve la lista de categorías (doble sufijo) a rellenar para este evento.

    Incluye siempre ALL_ALL.
    Incluye la combinación exacta.
    Si this_cat es un BG específico, incluye también BG_{other_cat}.
    Si other_cat es un BG específico, incluye también {this_cat}_BG.
    """
    cats = ["ALL_ALL"]
    exact = f"{this_cat}_{other_cat}"
    if exact != "ALL_ALL":
        cats.append(exact)

    this_is_bg = this_cat not in ("SIGNAL",)
    other_is_bg = other_cat not in ("SIGNAL",)

    if this_is_bg and exact != f"BG_{other_cat}":
        cats.append(f"BG_{other_cat}")
    if other_is_bg and exact != f"{this_cat}_BG":
        cats.append(f"{this_cat}_BG")
    if this_is_bg and other_is_bg:
        cats.append("BG_BG")

    return cats


def _fill_hemisphere(hists, dec_vars, other_vars, shared_vars, dec_idx,
                     this_cat_gen, other_cat_gen, this_cat_reco, other_cat_reco,
                     weight):
    """Rellena todos los histogramas de un hemisferio según FILL_RULES_PER_DEC.

    dec_idx: 0 o 1 (para construir el sufijo _dec0 / _dec1)
    La clasificación usa gen para niveles Gen/Matched y reco para Reco.
    """
    dec_suffix = f"dec{dec_idx}"

    for rule in FILL_RULES_PER_DEC:
        level, base_name, x_fn, y_fn = rule[:4]
        cond_fn = rule[4] if len(rule) > 4 else None

        if cond_fn is not None and not cond_fn(dec_vars, shared_vars):
            continue

        var_name = f"{base_name}_{dec_suffix}"

        if level not in hists or var_name not in hists[level]:
            continue

        # Clasificación por nivel
        if level == "Reco":
            this_cat, other_cat = this_cat_reco, other_cat_reco
        else:
            this_cat, other_cat = this_cat_gen, other_cat_gen

        cats_to_fill = _get_fill_categories(this_cat, other_cat)

        x_val = x_fn(dec_vars, shared_vars)
        y_val = y_fn(dec_vars, shared_vars) if y_fn is not None else None

        for cat in cats_to_fill:
            if cat not in hists[level][var_name]:
                continue
            for w_name, w_hist in hists[level][var_name][cat].items():
                eff_w = WEIGHT_VALUES_MDECS[w_name](dec_vars, other_vars, weight)
                if y_val is None:
                    w_hist.Fill(x_val, eff_w)
                else:
                    w_hist.Fill(x_val, y_val, eff_w)


def _fill_shared(hists, vars_dec0, vars_dec1, shared_vars,
                 cat_dec0_gen, cat_dec1_gen, cat_dec0_reco, cat_dec1_reco,
                 weight, dR=None):
    """Rellena histogramas compartidos según FILL_RULES_SHARED.

    La categoría se calcula desde la perspectiva de dec0.
    """
    for rule in FILL_RULES_SHARED:
        level, hist_name, x_fn, y_fn = rule[:4]
        cond_fn = rule[4] if len(rule) > 4 else None

        if cond_fn is not None and not cond_fn(vars_dec0, vars_dec1, shared_vars):
            continue

        if level not in hists or hist_name not in hists[level]:
            continue

        if level == "Reco":
            this_cat, other_cat = cat_dec0_reco, cat_dec1_reco
        else:
            this_cat, other_cat = cat_dec0_gen, cat_dec1_gen

        cats_to_fill = _get_fill_categories(this_cat, other_cat)

        # DeltaR se maneja de forma especial
        if hist_name == "DeltaR_dec0_dec1":
            if dR is None:
                continue
            for cat in cats_to_fill:
                if cat not in hists[level][hist_name]:
                    continue
                for w_name, w_hist in hists[level][hist_name][cat].items():
                    eff_w = WEIGHT_VALUES_MDECS[w_name](vars_dec0, vars_dec1, weight)
                    w_hist.Fill(dR, eff_w)
            continue

        x_val = x_fn(vars_dec0, vars_dec1, shared_vars)
        y_val = y_fn(vars_dec0, vars_dec1, shared_vars) if y_fn is not None else None

        for cat in cats_to_fill:
            if cat not in hists[level][hist_name]:
                continue
            for w_name, w_hist in hists[level][hist_name][cat].items():
                eff_w = WEIGHT_VALUES_MDECS[w_name](vars_dec0, vars_dec1, weight)
                if y_val is None:
                    w_hist.Fill(x_val, eff_w)
                else:
                    w_hist.Fill(x_val, y_val, eff_w)


def _fill_zvismassbins(hists, vars_dec0, vars_dec1, shared_vars,
                       cat_dec0_gen, cat_dec1_gen, weight):
    """Rellena los histogramas de bins de ZVisMass (lógica condicional)."""
    z_vis = shared_vars.get("GenZVisMass", 0.0)
    bin_ranges = [(0, 40, 1), (40, 70, 2), (70, 100, 3)]

    for lo, hi, idx in bin_ranges:
        if not (lo <= z_vis < hi):
            continue
        for dec_idx, dec_vars, other_vars, this_cat, other_cat in [
            (0, vars_dec0, vars_dec1, cat_dec0_gen, cat_dec1_gen),
            (1, vars_dec1, vars_dec0, cat_dec1_gen, cat_dec0_gen),
        ]:
            var_name = f"Omega_ZVisMass_Bin{idx}_dec{dec_idx}"
            if "Gen" not in hists or var_name not in hists["Gen"]:
                continue
            cats_to_fill = _get_fill_categories(this_cat, other_cat)
            for cat in cats_to_fill:
                if cat not in hists["Gen"][var_name]:
                    continue
                for w_name, w_hist in hists["Gen"][var_name][cat].items():
                    eff_w = WEIGHT_VALUES_MDECS[w_name](dec_vars, other_vars, weight)
                    w_hist.Fill(dec_vars["omega"], eff_w)
        break


# ── Función principal de llenado por rango ────────────────────────────────────

_REQUIRED_CHARGE_BRANCHES = ["tau1_tauPDG", "tau2_tauPDG", "tau1_recoCharge", "tau2_recoCharge"]


def _require_charge_branches(trees):
    """Abort if the input trees predate the z = cos(theta_tau-) convention.

    The per-hemisphere charge fixes the sign of z in every polarization weight and
    in the fit axis. Branch reading falls back to 0.0 for missing branches, which
    would silently tag every hemisphere as a tau+ — the exact silent-failure mode
    this convention was introduced to remove. Fail loudly instead: such trees have
    to be regenerated (see docs/plan_signo_costheta_taum.md, phase 2).
    """
    for tree_key, tree in trees.items():
        present = {b.GetName() for b in tree.GetListOfBranches()}
        missing = [b for b in _REQUIRED_CHARGE_BRANCHES if b not in present]
        if missing:
            raise RuntimeError(
                f"Tree '{tree_key}' is missing the charge branches {missing}. "
                "It was produced before the z = cos(theta_tau-) convention and "
                "must be regenerated with analysisRHOTree_MDecs_parallel.py / "
                "genOnlyRHOTree_MDecs_parallel.py."
            )


def process_tree_range_mdecs(trees, root_histograms_super,
                              weight, decay_pair,
                              cuts_cfg, logger_process, other_BG_id,
                              start_entry, end_entry,
                              single_decay_id=None, only_gen=False,
                              sin_eff=None, compute_weights=False,
                              use_omega=False, use_costheta_pion=True,
                              use_nn_optimal=False, nn_model_path=None):
    """Procesa entradas [start_entry, end_entry) de los árboles MDecs.

    Devuelve dict de contadores: totalEvents, selectedEvents, sumWeights, ...
    """
    # Carga (cacheada por proceso) del MLP en numpy para la variable óptima NN.
    nn_model = None
    if use_nn_optimal:
        from modules import mlpPolInference
        nn_model = mlpPolInference.load_mlp(nn_model_path)
    id0_gen = 1 if decay_pair[0] == 2 else decay_pair[0]
    id1_gen = 1 if decay_pair[1] == 2 else decay_pair[1]
    single_decay_id_gen = (1 if single_decay_id == 2 else single_decay_id) if single_decay_id is not None else None

    tauPCut    = cuts_cfg.get("tauPCut",   0.0)
    meson_cut  = cuts_cfg.get("meson_cut",  [0.0, np.inf])
    lepton_cut = cuts_cfg.get("lepton_cut", [0.0, np.inf])
    zmass_cut  = cuts_cfg.get("zmass_cut",  [0.0, np.inf])
    angle_sep  = cuts_cfg.get("angle_sep",  [0.0, np.inf])
    cos_acc    = cuts_cfg.get("cos_acceptance", 1.0)  # corte aceptancia |cos(θ_mesón)|<=cos_acc (legacy: 0.95)
    omega_border = cuts_cfg.get("omega_border_cut", False)  # legacy: descarta |cos_theta_reco|==1 o |cos_psi_reco|==1
    extra_cuts = cuts_cfg.get("extra_cuts", [])
    vism_cut   = cuts_cfg.get("vism_cut", 0.0)  # corte masa visible del pión (decayID 0): visM < umbral; 0 = desactivado

    def _vism_ok(vd):
        """True si el hemisferio pasa el corte de masa visible (solo aplica a piones gen)."""
        if vism_cut <= 0:
            return True
        if int(vd.get("decayID", -999)) != 0:
            return True
        return vd.get("visM", 0.0) < vism_cut

    _require_charge_branches(trees)

    totalEvents    = 0
    selectedEvents = 0
    sumWeights     = 0.0
    sumWeightsP1   = 0.0
    sumWeightsM1   = 0.0

    for tree_key, tree in trees.items():
        root_hists = root_histograms_super[tree_key]

        for i in range(start_entry, end_entry):
            tree.GetEntry(i)
            entry = tree

            if tree_key == "original":
                totalEvents += 1

            # Extraer ramas compartidas
            shared_vars = extract_scalars_optional(entry, _SHARED_MDECS, default=0.0)
            beamE = shared_vars["beamE"]

            # Extraer ramas por tau
            tau1_vars = {k: float(getattr(entry, f"tau1_{k}", 0.0)) for k in _TAU_KEYS_MDECS}
            tau2_vars = {k: float(getattr(entry, f"tau2_{k}", 0.0)) for k in _TAU_KEYS_MDECS}

            # Fotones reco (vectores, no escalares): solo si se usa la red.
            if use_nn_optimal:
                tau1_vars["_reco_photons"] = _read_reco_photons(entry, "tau1")
                tau2_vars["_reco_photons"] = _read_reco_photons(entry, "tau2")

            if compute_weights and sin_eff is not None:
                _recompute_weights(tau1_vars, beamE, sin_eff, use_omega=use_omega,
                                   use_costheta_pion=use_costheta_pion)
                _recompute_weights(tau2_vars, beamE, sin_eff, use_omega=use_omega,
                                   use_costheta_pion=use_costheta_pion)
                # Variable óptima (observable) también al vuelo: misma definición que los
                # generadores (helper compartido), reutilizable con árboles antiguos.
                _recompute_optimal_var(tau1_vars, beamE)
                _recompute_optimal_var(tau2_vars, beamE)

            # Reco-level weights: always computed (visible-system proxy for tau direction)
            _sin = sin_eff if sin_eff is not None else 0.2312
            _recompute_reco_weights(tau1_vars, beamE, _sin, use_omega=use_omega)
            _recompute_reco_weights(tau2_vars, beamE, _sin, use_omega=use_omega)

            if single_decay_id_gen is not None:
                # ── Modo single-decay ──────────────────────────────────────────
                vars_this, vars_other = _assign_single_decay(
                    tau1_vars, tau2_vars, single_decay_id_gen,
                    id_reco=single_decay_id, only_gen=only_gen)
                if vars_this is None:
                    continue

                # Joint two-tau weights (Alcaraz 2026 eqs. 9/13/16): SIEMPRE la fórmula
                # joint completa con término cruzado, en gen y reco. El otro hemisferio
                # está disponible (vars_other). El producto de pesos per-tau es incorrecto
                # (doble-cuenta la polarización de producción) y ya no se usa salvo, dentro
                # de _compute_joint_weights, para decays no soportados.
                _compute_joint_weights(vars_this, vars_other, beamE, _sin, use_omega=use_omega,
                                       use_costheta_pion=use_costheta_pion)
                _compute_reco_joint_weights(vars_this, vars_other, beamE, _sin, use_omega=use_omega)

                for vd in (vars_this, vars_other):
                    vd["_optimal_x"] = (2.0 * vd["recoVisE"] / beamE - 1.0) if beamE else 0.0
                    vd["_omega_reco"] = _compute_omega_reco(vd, beamE)
                    # Variable óptima reco unificada por canal: ω reco para ρ (id 1 ó 2),
                    # x reco para el resto.
                    vd["_optimal_unified_reco"] = (
                        vd["_omega_reco"] if int(vd.get("recoTauID", -999)) in (1, 2) else vd["_optimal_x"])
                # Variable óptima del ρ vía red (solo par ρ-leptón); no toca ω.
                if use_nn_optimal:
                    _maybe_set_nn_optimal(vars_this, vars_other, nn_model)

                # Corte masa visible del pión (anti-contaminación): visM < umbral
                if not (_vism_ok(vars_this) and _vism_ok(vars_other)):
                    continue
                if vars_this["recoVisP"] < tauPCut:
                    continue
                zmass = shared_vars["ZRecoMass"]
                if not (zmass_cut[0] <= zmass <= zmass_cut[1]):
                    continue
                if not (meson_cut[0] <= vars_this["recoVisP"] <= meson_cut[1]):
                    continue
                # Aceptancia legacy: |cos(θ_mesón)| <= cos_acc (analysisRHOTree.py, 0.95)
                if abs(math.cos(vars_this["recoVisTheta"])) > cos_acc:
                    continue
                if not (lepton_cut[0] <= vars_other["recoVisP"] <= lepton_cut[1]):
                    continue

                p4_this = make_p4(vars_this["recoVisP"], vars_this["recoVisTheta"],
                                  vars_this["recoVisPhi"], vars_this["recoVisE"])
                p4_other = make_p4(vars_other["recoVisP"], vars_other["recoVisTheta"],
                                   vars_other["recoVisPhi"], vars_other["recoVisE"])
                dR = myutils.dRAngle(p4_this, p4_other)
                if not (angle_sep[0] <= dR <= angle_sep[1]):
                    continue
                # Corte de bordes de ω (legacy analysisRHOTree.py): descarta el ρ si su
                # cos_theta/cos_psi reco quedan clampeados a ±1 (cinemática no física).
                if omega_border:
                    _ct = vars_this.get("_ct_reco"); _cp = vars_this.get("_cp_reco")
                    if _ct is not None and (abs(_ct) >= 1.0 or abs(_cp) >= 1.0):
                        continue

                if extra_cuts:
                    skip = False
                    flat = {**vars_this, **{f"d1_{k}": v for k, v in vars_other.items()},
                            **shared_vars}
                    for expr in extra_cuts:
                        try:
                            if not eval(expr, {"__builtins__": {}}, flat):
                                skip = True
                                break
                        except Exception:
                            skip = True
                            break
                    if skip:
                        continue

                cat_this_gen   = _classify_hemisphere(vars_this,  single_decay_id_gen, use_reco=False,
                                                     only_gen=only_gen)
                cat_other_gen  = _classify_hemisphere(vars_other, None,                use_reco=False,
                                                     only_gen=only_gen)
                cat_this_reco  = _classify_hemisphere(vars_this,  single_decay_id_gen, use_reco=True,
                                                     expected_reco_id=single_decay_id,
                                                     only_gen=only_gen)
                cat_other_reco = _classify_hemisphere(vars_other, None,                use_reco=True,
                                                     only_gen=only_gen)

                if tree_key == "original":
                    selectedEvents += 1
                    if cat_this_gen == "SIGNAL":
                        sumWeights   += weight
                        sumWeightsP1 += weight * vars_this.get("weight_P1", 1.0)
                        sumWeightsM1 += weight * vars_this.get("weight_M1", 1.0)
                    else:
                        other_BG_id[f"{cat_this_gen}_{cat_other_gen}"] = (
                            other_BG_id.get(f"{cat_this_gen}_{cat_other_gen}", 0) + 1)

                # Solo dec_idx=0; _dec1 queda vacío en este modo
                _fill_hemisphere(root_hists, vars_this, vars_other, shared_vars, 0,
                                 cat_this_gen, cat_other_gen,
                                 cat_this_reco, cat_other_reco, weight)

                _fill_shared(root_hists, vars_this, vars_other, shared_vars,
                             cat_this_gen, cat_other_gen,
                             cat_this_reco, cat_other_reco,
                             weight, dR=dR)

                _fill_zvismassbins(root_hists, vars_this, vars_other, shared_vars,
                                   cat_this_gen, cat_other_gen, weight)

            else:
                # ── Modo pair (comportamiento original) ───────────────────────
                vars_dec0, vars_dec1 = _assign_hemispheres(
                    tau1_vars, tau2_vars, id0_gen, id1_gen,
                    reco_id0=decay_pair[0], reco_id1=decay_pair[1], only_gen=only_gen)
                if vars_dec0 is None:
                    continue

                # Joint two-tau weights (Alcaraz 2026 eqs. 9/13/16): SIEMPRE joint completo
                # (gen y reco), independiente de --compute-weights.
                _compute_joint_weights(vars_dec0, vars_dec1, beamE, _sin, use_omega=use_omega,
                                       use_costheta_pion=use_costheta_pion)
                _compute_reco_joint_weights(vars_dec0, vars_dec1, beamE, _sin, use_omega=use_omega)

                for vd in (vars_dec0, vars_dec1):
                    vd["_optimal_x"] = (2.0 * vd["recoVisE"] / beamE - 1.0) if beamE else 0.0
                    vd["_omega_reco"] = _compute_omega_reco(vd, beamE)
                    # Variable óptima reco unificada por canal: ω reco para ρ (id 1 ó 2),
                    # x reco para el resto.
                    vd["_optimal_unified_reco"] = (
                        vd["_omega_reco"] if int(vd.get("recoTauID", -999)) in (1, 2) else vd["_optimal_x"])
                # Variable óptima del ρ vía red (solo par ρ-leptón); no toca ω.
                if use_nn_optimal:
                    _maybe_set_nn_optimal(vars_dec0, vars_dec1, nn_model)

                # Corte masa visible del pión (anti-contaminación): visM < umbral
                if not (_vism_ok(vars_dec0) and _vism_ok(vars_dec1)):
                    continue
                if vars_dec0["recoVisP"] < tauPCut:
                    continue
                zmass = shared_vars["ZRecoMass"]
                if not (zmass_cut[0] <= zmass <= zmass_cut[1]):
                    continue
                if not (meson_cut[0] <= vars_dec0["recoVisP"] <= meson_cut[1]):
                    continue
                # Aceptancia legacy: |cos(θ_mesón)| <= cos_acc (analysisRHOTree.py, 0.95)
                if abs(math.cos(vars_dec0["recoVisTheta"])) > cos_acc:
                    continue
                if not (lepton_cut[0] <= vars_dec1["recoVisP"] <= lepton_cut[1]):
                    continue

                p4_dec0 = make_p4(vars_dec0["recoVisP"], vars_dec0["recoVisTheta"],
                                   vars_dec0["recoVisPhi"], vars_dec0["recoVisE"])
                p4_dec1 = make_p4(vars_dec1["recoVisP"], vars_dec1["recoVisTheta"],
                                   vars_dec1["recoVisPhi"], vars_dec1["recoVisE"])
                dR = myutils.dRAngle(p4_dec0, p4_dec1)
                if not (angle_sep[0] <= dR <= angle_sep[1]):
                    continue
                # Corte de bordes de ω (legacy analysisRHOTree.py): descarta el ρ si su
                # cos_theta/cos_psi reco quedan clampeados a ±1 (cinemática no física).
                if omega_border:
                    _ct = vars_dec0.get("_ct_reco"); _cp = vars_dec0.get("_cp_reco")
                    if _ct is not None and (abs(_ct) >= 1.0 or abs(_cp) >= 1.0):
                        continue

                if extra_cuts:
                    skip = False
                    flat = {**vars_dec0, **{f"d1_{k}": v for k, v in vars_dec1.items()},
                            **shared_vars}
                    for expr in extra_cuts:
                        try:
                            if not eval(expr, {"__builtins__": {}}, flat):
                                skip = True
                                break
                        except Exception:
                            skip = True
                            break
                    if skip:
                        continue

                cat_dec0_gen  = _classify_hemisphere(vars_dec0, id0_gen, use_reco=False,
                                                     only_gen=only_gen)
                cat_dec1_gen  = _classify_hemisphere(vars_dec1, id1_gen, use_reco=False,
                                                     only_gen=only_gen)
                cat_dec0_reco = _classify_hemisphere(vars_dec0, id0_gen, use_reco=True,
                                                     expected_reco_id=decay_pair[0],
                                                     only_gen=only_gen)
                cat_dec1_reco = _classify_hemisphere(vars_dec1, id1_gen, use_reco=True,
                                                     expected_reco_id=decay_pair[1],
                                                     only_gen=only_gen)

                if tree_key == "original":
                    selectedEvents += 1
                    if cat_dec0_gen == "SIGNAL" and cat_dec1_gen == "SIGNAL":
                        sumWeights   += weight
                        sumWeightsP1 += weight * vars_dec0.get("weight_P1", 1.0)
                        sumWeightsM1 += weight * vars_dec0.get("weight_M1", 1.0)
                    else:
                        other_BG_id[f"{cat_dec0_gen}_{cat_dec1_gen}"] = (
                            other_BG_id.get(f"{cat_dec0_gen}_{cat_dec1_gen}", 0) + 1)

                _fill_hemisphere(root_hists, vars_dec0, vars_dec1, shared_vars, 0,
                                 cat_dec0_gen, cat_dec1_gen, cat_dec0_reco, cat_dec1_reco, weight)
                _fill_hemisphere(root_hists, vars_dec1, vars_dec0, shared_vars, 1,
                                 cat_dec1_gen, cat_dec0_gen, cat_dec1_reco, cat_dec0_reco, weight)

                _fill_shared(root_hists, vars_dec0, vars_dec1, shared_vars,
                             cat_dec0_gen, cat_dec1_gen, cat_dec0_reco, cat_dec1_reco,
                             weight, dR=dR)

                _fill_zvismassbins(root_hists, vars_dec0, vars_dec1, shared_vars,
                                   cat_dec0_gen, cat_dec1_gen, weight)

    return {
        "totalEvents":    totalEvents,
        "selectedEvents": selectedEvents,
        "sumWeights":     sumWeights,
        "sumWeightsP1":   sumWeightsP1,
        "sumWeightsM1":   sumWeightsM1,
    }


# ── Worker ────────────────────────────────────────────────────────────────────

def process_chunk_stage2_mdecs(input_root, tree_keys, entry_range, config_bundle, worker_id):
    """Worker: procesa un rango de entradas y escribe histogramas parciales."""
    outputpath = config_bundle["outputpath"]

    root_logger = logging.getLogger()
    for h in root_logger.handlers[:]:
        root_logger.removeHandler(h)
        h.close()
    log_dir = os.path.join(outputpath, "logs", _LOG_SOURCE)
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"worker_{worker_id}.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        force=True,
    )
    logger = logging.getLogger(f"worker_{worker_id}")
    logger.info("Worker %d: entradas [%d, %d)", worker_id, entry_range[0], entry_range[1])

    hist_config      = dict(config_bundle["hist_config"])
    hist_config.pop("_all_cats", None)
    decay_pair       = config_bundle["decay_pair"]
    weight           = config_bundle["weight"]
    cuts_cfg         = config_bundle["cuts_cfg"]
    fileOutName_base = config_bundle["fileOutName_base"]
    single_decay_id  = config_bundle.get("single_decay_id", None)
    only_gen         = config_bundle.get("only_gen", False)
    sin_eff          = config_bundle.get("sin_eff", None)
    compute_weights  = config_bundle.get("compute_weights", False)
    use_omega        = config_bundle.get("use_omega", False)
    use_costheta_pion = config_bundle.get("use_costheta_pion", True)
    use_nn_optimal   = config_bundle.get("use_nn_optimal", False)
    nn_model_path    = config_bundle.get("nn_model_path", None)

    infile = TFile.Open(input_root, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Worker {worker_id}: no se puede abrir {input_root}")

    trees = {}
    for key in tree_keys:
        t = infile.Get(f"outtree_{key}")
        if isinstance(t, ROOT.TTree):
            trees[key] = t
        else:
            logger.warning("Worker %d: árbol 'outtree_%s' no encontrado", worker_id, key)

    partial_path = os.path.join(
        outputpath, f"partial_histos_worker{worker_id}_{fileOutName_base}.root")
    pfile = TFile(partial_path, "RECREATE")
    pfile.cd()

    base_hists = myutils.build_histogram_registry(hist_config)
    root_histograms_super = {"original": base_hists}

    if "min_err" in trees:
        root_histograms_super["min_err"] = myutils.clone_histograms_with_suffix(base_hists, "_min")
        _reattach_histograms(root_histograms_super["min_err"], pfile)
    if "max_err" in trees:
        root_histograms_super["max_err"] = myutils.clone_histograms_with_suffix(base_hists, "_max")
        _reattach_histograms(root_histograms_super["max_err"], pfile)

    other_BG_id = {}
    start_entry, end_entry = entry_range
    counters = process_tree_range_mdecs(
        trees=trees,
        root_histograms_super=root_histograms_super,
        weight=weight,
        decay_pair=decay_pair,
        cuts_cfg=cuts_cfg,
        logger_process=logger,
        other_BG_id=other_BG_id,
        start_entry=start_entry,
        end_entry=end_entry,
        single_decay_id=single_decay_id,
        only_gen=only_gen,
        sin_eff=sin_eff,
        compute_weights=compute_weights,
        use_omega=use_omega,
        use_costheta_pion=use_costheta_pion,
        use_nn_optimal=use_nn_optimal,
        nn_model_path=nn_model_path,
    )

    infile.Close()

    pfile.cd()
    for tree_key in root_histograms_super:
        write_histograms_recursive(root_histograms_super[tree_key])
    pfile.Close()

    logger.info("Worker %d: terminado. Events=%d Selected=%d",
                worker_id, counters["totalEvents"], counters["selectedEvents"])
    return partial_path, counters, other_BG_id


# ── Merge ─────────────────────────────────────────────────────────────────────

def merge_histogram_dicts_from_files(partial_files, root_histograms_super):
    all_main = {}
    for subtree in root_histograms_super.values():
        _flatten_histograms(subtree, all_main)

    for partial_path in partial_files:
        pfile = TFile.Open(partial_path, "READ")
        if not pfile or pfile.IsZombie():
            print(f"[WARN] No se puede abrir parcial {partial_path}")
            continue
        for name, main_hist in all_main.items():
            partial_hist = pfile.Get(name)
            if partial_hist and isinstance(partial_hist, ROOT.TH1):
                main_hist.Add(partial_hist)
        pfile.Close()


# ── Hook y main ───────────────────────────────────────────────────────────────

def my_hook(parser):
    parser.add_argument("--tree-file", type=str, required=True,
        help="Input ROOT file con outtree_original (y opcionalmente outtree_min_err/max_err)")
    parser.add_argument("--ang", type=float, default=[0.0, np.inf], nargs="+",
        help="Separación angular entre decays (por defecto: 0.0 a infinito)")
    parser.add_argument("--meson-cut", type=float, default=[0.0, np.inf], nargs="+",
        help="Rango de corte en P del dec0 (por defecto: sin corte)")
    parser.add_argument("--lepton-cut", type=float, default=[0.0, np.inf], nargs="+",
        help="Rango de corte en P del dec1 (por defecto: sin corte)")
    parser.add_argument("--zmass-cut", type=float, default=[0.0, np.inf], nargs="+",
        help="Rango de masa Z reco (por defecto: sin corte)")
    parser.add_argument("--cos-acceptance", type=float, default=1.0,
        help="Corte de aceptancia legacy: rechaza el evento si |cos(θ_mesón_vis)| > este "
             "valor (analysisRHOTree.py usa 0.95). Default 1.0 = desactivado.")
    parser.add_argument("--omega-border-cut", action="store_true", default=False,
        help="Corte de bordes de ω legacy: descarta el ρ si su cos_theta o cos_psi reco "
             "(wVariabRECO) quedan clampeados a ±1 (cinemática fuera del rango τ→ρ, "
             "enriquecida en fakes). Añade 'omegaBorder_' al nombre. Default: desactivado.")
    parser.add_argument("--hist-config-mdecs", type=str,
        default="config/histograms/rho_analysis_config_mdecs.yml",
        help="Config YAML de histogramas MDecs")
    parser.add_argument("--n-workers", type=int, default=None,
        help="Número de workers paralelos (default: min(n_entries, n_cpus))")
    parser.add_argument("--cut", type=str, nargs="+", default=[], metavar="EXPR",
        help="Expresiones Python extra sobre vars_dict (vars de dec0 + d1_* para dec1)")
    parser.add_argument("--decay-pair", type=int, nargs=2, default=None,
        metavar=("DECID_0", "DECID_1"),
        help="Par de decayIDs a analizar (sobrescribe general.decay_pair del YAML). "
             "Ejemplo: --decay-pair 2 -13")
    parser.add_argument("--single-decay", type=int, default=None,
        metavar="DECID",
        help="Modo single-decay: selecciona todos los eventos donde al menos un tau "
             "tenga este decayID (ej. 2=rho, -13=muón). El otro lado se detecta "
             "automáticamente y sus pesos se usan en corr_P1/corr_M1. "
             "Incompatible con --decay-pair.")
    parser.add_argument("--only-gen", action="store_true", default=False,
        help="Indica que el árbol es gen-only (generado con genOnlyRHOTree): "
             "las ramas reco contienen IDs gen, por lo que se aplica el remap "
             "2→1 también en la clasificación reco. "
             "Por defecto (árbol reco real) el remap solo se aplica en gen.")
    parser.add_argument("--sin-eff", type=float, default=None, metavar="SIN2THETA",
        help="Valor de sin²θ_eff para recalcular weight_P1/weight_M1 on-the-fly "
             "(requiere --compute-weights). Por defecto usa los pesos del árbol.")
    parser.add_argument("--compute-weights", action="store_true", default=False,
        help="Recalcula weight_P1/weight_M1 desde las cinemáticas almacenadas "
             "en lugar de leer los pesos del árbol. Necesario para usar --sin-eff.")
    # Por defecto se usa la variable óptima ω (fórmula completa, con cosPsi/cosBeta) como
    # H para el canal rho, tanto a nivel reco (newAtauRhoOmega) como gen (al recalcular con
    # --compute-weights). Para desactivarlo y usar la H simplificada H_V = alpha_V*z_R: --no-omega-weights.
    parser.add_argument("--no-omega-weights", dest="omega_weights", action="store_false",
        help="Desactiva ω: para el canal rho usa la H simplificada H_V = alpha_V * z_R "
             "(newAtau) en lugar de la variable óptima ω. Añade 'noOmega_' al nombre de salida. "
             "Por defecto ω está ACTIVADO (reco y gen).")
    parser.set_defaults(omega_weights=True)
    # Por defecto, para el canal pión gen se usa cos(θ*) geométrico (boost exacto, α=1) como H,
    # recalculado al vuelo con --compute-weights (funciona también con árboles antiguos). El reco
    # del pión SIEMPRE usa z_R (aproximación). Para usar z_R también en gen: --no-costheta-pion.
    parser.add_argument("--no-costheta-pion", dest="costheta_pion", action="store_false",
        help="Desactiva cos(θ*) geométrico para el pión gen: usa la H analítica z_R "
             "(newAtau) en lugar del ángulo del boost. Añade 'zRpion_' al nombre de salida. "
             "Por defecto cos(θ*) está ACTIVADO (solo gen; reco siempre z_R).")
    parser.set_defaults(costheta_pion=True)
    parser.add_argument("--vism-cut", type=float, default=0.2,
        help="Corte anti-contaminación del pión: descarta el hemisferio pión (decayID 0) "
             "si su masa visible visM >= este valor (un π limpio tiene visM≈m_π=0.14; los "
             "contaminados con neutros extra suben). El observable sigue siendo E_vis. "
             "Por defecto 0.2 (ACTIVADO). Pon 0 para desactivarlo.")
    parser.add_argument("--merge-hemispheres", action="store_true", default=False,
        help="Para muestras con la MISMA partícula en los dos hemisferios: funde cada par "
             "de histogramas {base}_dec0 + {base}_dec1 en uno solo {base} (suma de ambos "
             "hemisferios) y guarda SOLO los fusionados (sin los dec0/dec1 por separado). "
             "Los cross-hemisferio (*_dec0_vs_dec1) se conservan. Añade 'merged_' al nombre.")
    parser.add_argument("--use-nn-optimal", action="store_true", default=False,
        help="Solo para pares ρ-leptón: calcula la variable óptima reco del ρ con la red "
             "neuronal (modelo v2) en lugar de ω. Su salida (∈[0,1]) llena los histogramas "
             "OptimalNN_* (NO se usa para repesado ni toca los histogramas de ω).")
    parser.add_argument("--nn-model", type=str,
        default="MLPolResults/train_pol_results_optimal_lepton/mlp_numpy.npz",
        help="Ruta al .npz del MLP exportado a numpy (ver RhoAnalysis/MLP/exportMLPToNumpy.py). "
             "Por defecto el modelo optimal_lepton (features del pión).")


def main():
    general_configs = myutils.setup_analysis_config(
        _DEFAULT_CONFIG, _OUTPUT_BASE, parser_hook=my_hook,
        log_subdir=_LOG_SOURCE)
    loggers    = general_configs["loggers"]
    run_config = general_configs["config"]
    args       = general_configs["args"]
    logger_config  = loggers["config"]
    logger_io      = loggers["io"]
    logger_process = loggers["processing"]

    single_decay_id = args.single_decay
    only_gen        = args.only_gen
    sin_eff         = args.sin_eff
    compute_weights = args.compute_weights
    use_omega       = args.omega_weights
    use_costheta_pion = args.costheta_pion
    merge_hemispheres = args.merge_hemispheres
    use_nn_optimal  = args.use_nn_optimal
    nn_model_path   = args.nn_model

    if single_decay_id is not None and args.decay_pair is not None:
        logger_io.error("--single-decay y --decay-pair son incompatibles; usa solo uno.")
        sys.exit(1)

    if single_decay_id is not None:
        # Modo single-decay: decay_pair solo se usa internamente como dummy en process_tree
        decay_pair = [single_decay_id, single_decay_id]
        logger_config.info("single_decay_id: %d", single_decay_id)
    else:
        # Leer decay_pair: CLI tiene prioridad sobre el YAML
        if args.decay_pair is not None:
            decay_pair = args.decay_pair
        else:
            decay_pair = run_config.get("general", {}).get("decay_pair")
        if decay_pair is None or len(decay_pair) != 2:
            logger_io.error(
                "Especifica el par de desintegraciones con --decay-pair DECID_0 DECID_1, "
                "con 'general.decay_pair: [decID_0, decID_1]' en el YAML de config, "
                "o usa --single-decay DECID para modo single-decay."
            )
            sys.exit(1)
        decay_pair = [int(d) for d in decay_pair]
        logger_config.info("decay_pair: %s", decay_pair)

    if merge_hemispheres and single_decay_id is None and decay_pair[0] != decay_pair[1]:
        logger_io.warning(
            "--merge-hemispheres con par asimétrico %s: se fundirán hemisferios de distinta "
            "especie; normalmente solo tiene sentido con la misma partícula en ambos.", decay_pair)

    tauPCut = run_config["cuts"]["tauCut"]

    angle_sep  = args.ang;        angle_sep  = [angle_sep[0],  np.inf] if len(angle_sep)  == 1 else angle_sep
    meson_cut  = args.meson_cut;  meson_cut  = [meson_cut[0],  np.inf] if len(meson_cut)  == 1 else meson_cut
    lepton_cut = args.lepton_cut; lepton_cut = [lepton_cut[0], np.inf] if len(lepton_cut) == 1 else lepton_cut
    zmass_cut  = args.zmass_cut;  zmass_cut  = [zmass_cut[0],  np.inf] if len(zmass_cut)  == 1 else zmass_cut

    input_root = args.tree_file
    if not os.path.isfile(input_root):
        logger_io.error("Input ROOT file %s not found", input_root)
        sys.exit(1)

    # Cargar config de histogramas
    with open(args.hist_config_mdecs, "r") as f:
        hist_config = yaml.safe_load(f)
    # Eliminar anclas YAML auxiliares que no son secciones de histogramas
    hist_config.pop("_all_cats", None)

    # Detectar árboles y número de entradas
    infile = TFile.Open(input_root, "READ")
    if not infile or infile.IsZombie():
        logger_io.error("Could not open %s", input_root)
        sys.exit(1)
    tree_keys = []
    for key in ["original", "min_err", "max_err"]:
        t = infile.Get(f"outtree_{key}")
        if isinstance(t, ROOT.TTree) and t.GetEntries() > 0:
            tree_keys.append(key)
            logger_io.info("Found tree outtree_%s with %d entries", key, t.GetEntries())
    if "original" not in tree_keys:
        logger_io.error("outtree_original not found or empty in %s", input_root)
        sys.exit(1)
    n_entries = infile.Get("outtree_original").GetEntries()
    infile.Close()

    weight = 1.0
    cuts_cfg = {
        "tauPCut":   tauPCut,
        "meson_cut": meson_cut,
        "lepton_cut": lepton_cut,
        "zmass_cut": zmass_cut,
        "angle_sep": angle_sep,
        "cos_acceptance": args.cos_acceptance,
        "omega_border_cut": args.omega_border_cut,
        "extra_cuts": args.cut,
        "vism_cut": args.vism_cut,
    }

    outputpath = os.path.dirname(input_root)
    if single_decay_id is not None:
        out_prefix = f"HistosMDecs_single{single_decay_id}_"
    else:
        out_prefix = f"HistosMDecs_{decay_pair[0]}_{decay_pair[1]}_"
    # Marca de versión del convenio de signo de z. Con zTaum (z = cos θ_τ⁻ referido
    # SIEMPRE al τ⁻) los resultados dejan de ser comparables con los anteriores,
    # que usaban el cos θ del hemisferio: la marca evita pisarlos. Va DESPUÉS de los
    # ids del par para no romper _parse_pair_ids de makeCosBins_MDecs.py.
    out_prefix += "zTaum_"
    if angle_sep[0] > 0:
        out_prefix += f"dRgt{angle_sep[0]}_{angle_sep[1]}_"
    if meson_cut[0] > 0 or meson_cut[1] < 100:
        out_prefix += f"Dec0Pgt{meson_cut[0]}_lt{meson_cut[1]}_"
    if lepton_cut[0] > 0 or lepton_cut[1] < 100:
        out_prefix += f"Dec1Pgt{lepton_cut[0]}_lt{lepton_cut[1]}_"
    if zmass_cut[0] > 0 or zmass_cut[1] < 200:
        out_prefix += f"Zmassgt{zmass_cut[0]}_lt{zmass_cut[1]}_"
    if args.cos_acceptance < 1.0:
        out_prefix += f"cosAcc{args.cos_acceptance}_"
    if args.omega_border_cut:
        out_prefix += "omegaBorder_"
    if args.cut:
        safe = "_".join(e.replace(" ", "").replace("==", "eq").replace(">", "gt").replace("<", "lt")
                        for e in args.cut)
        out_prefix += f"cut_{safe}_"
    if sin_eff is not None:
        out_prefix += f"sineff{sin_eff}_"
    if not use_omega:
        out_prefix += "noOmega_"
    if not use_costheta_pion:
        out_prefix += "zRpion_"
    if use_nn_optimal:
        out_prefix += "NN_"
    if merge_hemispheres:
        out_prefix += "merged_"
    if args.vism_cut <= 0:
        out_prefix += "vismOff_"
    else:
        out_prefix += f"vism{args.vism_cut:g}_"
    fileOutName = os.path.join(outputpath, out_prefix + general_configs["fileOutName"])
    fileOutName_base = Path(fileOutName).stem
    os.makedirs(outputpath, exist_ok=True)

    n_workers = args.n_workers or min(n_entries, os.cpu_count() or 1)
    logger_io.info("n_entries=%d, n_workers=%d", n_entries, n_workers)

    # Modo secuencial
    if n_workers == 1:
        logger_io.info("Sequential mode.")
        base_hists = myutils.build_histogram_registry(hist_config)
        root_histograms_super = {"original": base_hists}

        infile = TFile.Open(input_root, "READ")
        trees = {}
        for key in tree_keys:
            t = infile.Get(f"outtree_{key}")
            if isinstance(t, ROOT.TTree):
                trees[key] = t
        if "min_err" in trees:
            root_histograms_super["min_err"] = myutils.clone_histograms_with_suffix(base_hists, "_min")
        if "max_err" in trees:
            root_histograms_super["max_err"] = myutils.clone_histograms_with_suffix(base_hists, "_max")

        other_BG_id = {}
        counters = process_tree_range_mdecs(
            trees=trees,
            root_histograms_super=root_histograms_super,
            weight=weight,
            decay_pair=decay_pair,
            cuts_cfg=cuts_cfg,
            logger_process=logger_process,
            other_BG_id=other_BG_id,
            start_entry=0,
            end_entry=n_entries,
            single_decay_id=single_decay_id,
            only_gen=only_gen,
            sin_eff=sin_eff,
            compute_weights=compute_weights,
            use_omega=use_omega,
            use_costheta_pion=use_costheta_pion,
        )
        infile.Close()
        _write_output_mdecs(fileOutName, root_histograms_super, counters, other_BG_id, logger_io,
                            merge_hemispheres=merge_hemispheres)
        return

    # Modo paralelo
    entry_ranges = split_entry_ranges(n_entries, n_workers)
    n_chunks     = len(entry_ranges)

    config_bundle = {
        "hist_config":      hist_config,
        "decay_pair":       decay_pair,
        "weight":           weight,
        "cuts_cfg":         cuts_cfg,
        "outputpath":       outputpath,
        "fileOutName_base": fileOutName_base,
        "single_decay_id":  single_decay_id,
        "only_gen":         only_gen,
        "sin_eff":          sin_eff,
        "compute_weights":  compute_weights,
        "use_omega":        use_omega,
        "use_costheta_pion": use_costheta_pion,
        "use_nn_optimal":   use_nn_optimal,
        "nn_model_path":    nn_model_path,
    }

    ctx = multiprocessing.get_context("fork")
    partial_files = []
    all_counters  = []
    all_bg_ids    = []
    t_start = time.time()

    logger_io.info("Launching %d workers...", n_chunks)
    with ProcessPoolExecutor(max_workers=n_chunks, mp_context=ctx) as executor:
        futures = {
            executor.submit(
                process_chunk_stage2_mdecs,
                input_root, tree_keys, entry_ranges[i], config_bundle, i,
            ): i
            for i in range(n_chunks)
        }
        for n_done, future in enumerate(as_completed(futures), start=1):
            wid = futures[future]
            try:
                partial_path, counters, bg_ids = future.result()
                partial_files.append(partial_path)
                all_counters.append(counters)
                all_bg_ids.append(bg_ids)
                elapsed = time.time() - t_start
                print(f"  [{n_done}/{n_chunks}] worker {wid} terminado "
                      f"({counters['totalEvents']} eventos, {elapsed:.1f}s)", flush=True)
            except Exception as exc:
                logger_io.error("Worker %d falló: %s", wid, exc)
                raise

    # Reconstruir histogramas vacíos y fusionar
    base_hists = myutils.build_histogram_registry(hist_config)
    root_histograms_super = {"original": base_hists}
    if "min_err" in tree_keys:
        root_histograms_super["min_err"] = myutils.clone_histograms_with_suffix(base_hists, "_min")
    if "max_err" in tree_keys:
        root_histograms_super["max_err"] = myutils.clone_histograms_with_suffix(base_hists, "_max")

    logger_io.info("Merging %d partial files...", len(partial_files))
    merge_histogram_dicts_from_files(partial_files, root_histograms_super)

    for p in partial_files:
        try:
            os.remove(p)
        except OSError:
            pass

    merged_counters = {
        "totalEvents":    sum(c["totalEvents"]    for c in all_counters),
        "selectedEvents": sum(c["selectedEvents"]  for c in all_counters),
        "sumWeights":     sum(c["sumWeights"]      for c in all_counters),
        "sumWeightsP1":   sum(c["sumWeightsP1"]    for c in all_counters),
        "sumWeightsM1":   sum(c["sumWeightsM1"]    for c in all_counters),
    }
    merged_bg_ids = {}
    for d in all_bg_ids:
        for k, v in d.items():
            merged_bg_ids[k] = merged_bg_ids.get(k, 0) + v

    logger_io.info("Total: events=%d selected=%d (%.1fs)",
                   merged_counters["totalEvents"], merged_counters["selectedEvents"],
                   time.time() - t_start)
    
    _write_output_mdecs(fileOutName, root_histograms_super, merged_counters, merged_bg_ids, logger_io,
                        merge_hemispheres=merge_hemispheres)


def _write_merged_hemispheres(nested, logger_io):
    """Funde los histogramas por-hemisferio y escribe SOLO los fusionados.

    Para cada par {base}_dec0{rest} + {base}_dec1{rest} crea {base}{rest} = dec0 + dec1
    (suma de ambos hemisferios) y lo escribe sin el sufijo _dec0/_dec1. Pensado para
    muestras con la MISMA partícula en los dos hemisferios (p.ej. --decay-pair 0 0).

    Los histogramas sin pareja por-hemisferio se escriben tal cual: los compartidos
    (ZMass…) y los cross-hemisferio (*_dec0_vs_dec1, DeltaR_dec0_dec1, *_dec0_vs_*_dec1),
    cuya pareja con _dec0→_dec1 no existe, así que no se funden.

    Debe llamarse con el directorio de salida ya activo (outfile.cd()); los Clone se
    crean ahí y se escriben.
    """
    flat = _flatten_histograms(nested)
    consumed = set()
    n_merged = 0
    for name, h in flat.items():
        if "_dec0" not in name or name in consumed:
            continue
        name1 = name.replace("_dec0", "_dec1", 1)
        if name1 == name or name1 not in flat:
            continue  # no es par por-hemisferio (p.ej. *_dec0_vs_dec1)
        merged = h.Clone(name.replace("_dec0", "", 1))
        merged.Add(flat[name1])
        merged.Write()
        consumed.add(name)
        consumed.add(name1)
        n_merged += 1
    for name, h in flat.items():
        if name not in consumed:
            h.Write()
    logger_io.info("Hemisferios fusionados: %d pares dec0+dec1 -> 1 histograma (sin sufijo dec)",
                   n_merged)


def _write_output_mdecs(fileOutName, root_histograms_super, counters, other_BG_id, logger_io,
                        merge_hemispheres=False):
    import pandas as pd
    logger_io.info("Writing output ROOT file %s", fileOutName)
    logger_io.info("Events=%d Selected=%d", counters["totalEvents"], counters["selectedEvents"])

    outfile = TFile(fileOutName, "RECREATE")
    outfile.cd()
    for tree_key in root_histograms_super:
        if merge_hemispheres:
            _write_merged_hemispheres(root_histograms_super[tree_key], logger_io)
        else:
            write_histograms_recursive(root_histograms_super[tree_key])
                
            # hist_to_merge={
            #     "OptimalVar_NegHel_dec0_ALL_ALL":
            #     ["OptimalVar_minus_NegHel_dec0_ALL_ALL",
            #     "OptimalVar_plus_PosHel_dec0_ALL_ALL"],
            #     "OptimalVar_NegHel_dec1_ALL_ALL":
            #     ["OptimalVar_minus_NegHel_dec1_ALL_ALL",
            #     "OptimalVar_plus_PosHel_dec1_ALL_ALL"],
            #     "OptimalVar_PosHel_dec0_ALL_ALL":
            #     ["OptimalVar_minus_PosHel_dec0_ALL_ALL",
            #     "OptimalVar_plus_NegHel_dec0_ALL_ALL"],
            #     "OptimalVar_PosHel_dec1_ALL_ALL":
            #     ["OptimalVar_minus_PosHel_dec1_ALL_ALL",
            #     "OptimalVar_plus_NegHel_dec1_ALL_ALL"],
            #     "Omega_NegHel_dec0_ALL_ALL"
            #     : ["Omega_minus_NegHel_dec0_ALL_ALL",
            #        "Omega_plus_PosHel_dec0_ALL_ALL"],
            #     "Omega_NegHel_dec1_ALL_ALL"
            #     : ["Omega_minus_NegHel_dec1_ALL_ALL",
            #        "Omega_plus_PosHel_dec1_ALL_ALL"],
            #     "Omega_PosHel_dec0_ALL_ALL"
            #     : ["Omega_minus_PosHel_dec0_ALL_ALL",
            #        "Omega_plus_NegHel_dec0_ALL_ALL"],
            #     "Omega_PosHel_dec1_ALL_ALL"
            #     : ["Omega_minus_PosHel_dec1_ALL_ALL",
            #        "Omega_plus_NegHel_dec1_ALL_ALL"],
            # })
    outfile.Close()

    csv_name = fileOutName.replace(".root", "_otherBGid.csv")
    df = pd.DataFrame(sorted(other_BG_id.items()), columns=["category", "count"])
    df.to_csv(csv_name, index=False)
    logger_io.info("BG category counts saved to %s", csv_name)
    logger_io.info("Done. Results in %s", fileOutName)


if __name__ == "__main__":
    main()
