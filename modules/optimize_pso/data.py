# cutopt/data.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import numpy as np
import ROOT
import uproot

import yaml

def load_dataset_weights(compare_yaml_path: str, selectGEN: int) -> dict:
    with open(compare_yaml_path, "r") as f:
        cfg = yaml.safe_load(f)

    lumi = float(cfg["global"]["luminosity"])

    weights = {}
    for ds in cfg["datasets"]:
        if selectGEN == 0:
            path = ds["path"].replace("decay2", "decay0")
        else:
            path = ds["path"]
        xsec = float(ds["xsec"])
        n_events = float(ds["n_events"])
        weights[path] = lumi * xsec / n_events

    return weights

def _eta_from_theta(theta: np.ndarray) -> np.ndarray:
    # eta = -ln(tan(theta/2))
    # clamp para evitar infs numéricos extremos si theta≈0 o theta≈pi
    eps = 1e-12
    t = np.clip(theta, eps, np.pi - eps)
    return -np.log(np.tan(t / 2.0))


def _delta_phi(phi1: np.ndarray, phi2: np.ndarray) -> np.ndarray:
    dphi = phi1 - phi2
    # wrap a [-pi, pi]
    return (dphi + np.pi) % (2.0 * np.pi) - np.pi


def delta_r_from_spherical(
    theta1: np.ndarray, phi1: np.ndarray,
    theta2: np.ndarray, phi2: np.ndarray
) -> np.ndarray:
    dtheta = theta1-theta2
    dphi = _delta_phi(phi1, phi2)
    return np.sqrt(dtheta * dtheta + dphi * dphi)

def invariant_mass_from_spherical(
    p1: np.ndarray, theta1: np.ndarray, phi1: np.ndarray, E1: np.ndarray,
    p2: np.ndarray, theta2: np.ndarray, phi2: np.ndarray, E2: np.ndarray,
) -> np.ndarray:
    """
    Calcula de forma vectorizada la masa invariante del sistema formado
    por dos partículas dadas en coordenadas esféricas.

    Parámetros
    ----------
    p1, theta1, phi1, E1 : np.ndarray
        Momento, ángulo polar, ángulo azimutal y energía de la partícula 1
    p2, theta2, phi2, E2 : np.ndarray
        Momento, ángulo polar, ángulo azimutal y energía de la partícula 2

    Retorna
    -------
    np.ndarray
        Masa invariante del sistema (>= 0)
    """

    # Componentes cartesianas partícula 1
    px1 = p1 * np.sin(theta1) * np.cos(phi1)
    py1 = p1 * np.sin(theta1) * np.sin(phi1)
    pz1 = p1 * np.cos(theta1)

    # Componentes cartesianas partícula 2
    px2 = p2 * np.sin(theta2) * np.cos(phi2)
    py2 = p2 * np.sin(theta2) * np.sin(phi2)
    pz2 = p2 * np.cos(theta2)

    # Suma de cuadrivectores
    E  = E1 + E2
    px = px1 + px2
    py = py1 + py2
    pz = pz1 + pz2

    # Masa invariante
    m2 = E*E - (px*px + py*py + pz*pz)

    # Protección numérica
    return np.sqrt(np.maximum(m2, 0.0))


@dataclass(frozen=True)
class BranchMap:
    # Nombres de ramas. Ajusta aquí si en tus trees difieren.
    beamE: str = "beamE"
    recoMesonE: str = "recoMesonE"
    recoMesonP: str = "recoMesonP"
    recoMesonTheta: str = "recoMesonTheta"
    recoMesonPhi: str = "recoMesonPhi"
    lepP: str = "lepP"
    lepTheta: str = "lepTheta"
    lepPhi: str = "lepPhi"
    lepE: str = "lepE"
    genTauID: str = "genTauID"
    # opcional
    weight: Optional[str] = None  # por defecto weight=1
    omega_gen: Optional[str] = None
    omega_gen_p1: Optional[str] = None
    omega_gen_m1: Optional[str] = None
    omega_reco: Optional[str] = None
    
    # omega_bg: Optional[str] = None
    
    


@dataclass
class RootSample:
    name: str
    path: str
    tree: str = "outtree_original"
    is_signal_file: bool = False  # True solo para el fichero “signal”


# @dataclass
# class LoadedSample:
#     name: str
#     is_signal_file: bool
#     arrays: Dict[str, np.ndarray]
#     # precomputados
#     dR: np.ndarray
#     weights: np.ndarray
#     genTauID: np.ndarray

#     # baseline mask opcional (cortes fijos)
#     base_mask: np.ndarray

@dataclass
class LoadedSample:
    name: str
    is_signal_file: bool
    arrays: Dict[str, np.ndarray]
    dR: np.ndarray
    ZMass: np.ndarray        # ⬅️ variable añadida
    mesonX: np.ndarray     # ⬅️ variable añadida
    weights: np.ndarray          # ⬅️ peso FINAL por evento
    genTauID: np.ndarray
    base_mask: np.ndarray

def load_samples(
    samples: Sequence[RootSample],
    branches: BranchMap,
    tauP_min_fixed: float,
    dataset_weights: dict,
) -> List[LoadedSample]:
    """
    Carga los ROOT a arrays una sola vez y precomputa dR.
    Aplica un baseline mask fijo: recoMesonP >= tauP_min_fixed
    """
    out: List[LoadedSample] = []

    needed = [
        branches.beamE,
        branches.recoMesonE,
        branches.recoMesonP,
        branches.recoMesonTheta,
        branches.recoMesonPhi,
        branches.lepP,
        branches.lepTheta,
        branches.lepPhi,
        branches.lepE,
        branches.genTauID,
    ]
    if branches.weight:
        needed.append(branches.weight)
    if branches.omega_gen:
        needed.append(branches.omega_gen)
    if branches.omega_gen_p1:
        needed.append(branches.omega_gen_p1)
        needed.append(branches.omega_gen_m1)
    if branches.omega_reco:
        needed.append(branches.omega_reco)
    # if branches.recomeson_X:
    #     needed.append(branches.recomeson_X)
    # if branches.meson_X_p1:
    #     needed.append(branches.meson_X_p1)  
    #     needed.append(branches.meson_X_m1)
    #     needed.append(branches.omega_signal)
    # if branches.omega_bg:
    #     needed.append(branches.omega_bg)

    for s in samples:
        with uproot.open(s.path) as f:
          t = f[s.tree]
          arr = t.arrays(needed, library="np")

        # Convertir a float64/int64 donde convenga
        beamE = np.asarray(arr[branches.beamE], dtype=np.float64)
        recoMesonP = np.asarray(arr[branches.recoMesonP], dtype=np.float64)
        recoMesonTheta = np.asarray(arr[branches.recoMesonTheta], dtype=np.float64)
        recoMesonPhi = np.asarray(arr[branches.recoMesonPhi], dtype=np.float64)
        recoMesonE = np.asarray(arr[branches.recoMesonE], dtype=np.float64)
        lepTheta = np.asarray(arr[branches.lepTheta], dtype=np.float64)
        lepPhi = np.asarray(arr[branches.lepPhi], dtype=np.float64)
        lepP = np.asarray(arr[branches.lepP], dtype=np.float64)
        lepE = np.asarray(arr[branches.lepE], dtype=np.float64)
        
        ZMass = invariant_mass_from_spherical(
            recoMesonP, recoMesonTheta, recoMesonPhi, recoMesonE,
            lepP, lepTheta, lepPhi, lepE,
        )
        
        mesonX = recoMesonE / beamE
        dR = delta_r_from_spherical(recoMesonTheta, recoMesonPhi, lepTheta, lepPhi)

        genTauID = np.asarray(arr[branches.genTauID], dtype=np.int64)

        if branches.weight:
            w_evt = np.asarray(arr[branches.weight], dtype=np.float64)
        else:
            w_evt = np.ones_like(recoMesonP, dtype=np.float64)

        w_dataset = dataset_weights[s.path]   # ⬅️ CLAVE
        weights = w_evt * w_dataset

        base_mask = recoMesonP >= float(tauP_min_fixed)

        out.append(
            LoadedSample(
                name=s.name,
                is_signal_file=s.is_signal_file,
                arrays={k: np.asarray(v) for k, v in arr.items()},
                dR=dR,
                ZMass=ZMass,
                mesonX=mesonX,
                weights=weights,
                genTauID=genTauID,
                base_mask=base_mask,
            )
        )

    return out
