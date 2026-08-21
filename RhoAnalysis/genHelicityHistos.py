#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genHelicityHistos.py

Script sencillo y autónomo. Partiendo del TTree gen-only generado desde la
muestra KKMC 10M (ramas tau1_*/tau2_* por evento), produce y guarda en un
rootfile:

  1) cos(theta) del tau (z = cos de la theta polar del tau), para las
     variaciones SM, Hel+1, Hel-1, P1, M1, corr_P1, corr_M1, separado por
     carga del tau (tau- = minus, tau+ = plus).

  2) fz = |z| / (1 + z^2), con las mismas variaciones y cargas, y además
     separado en forward (z>0) y backward (z<0).

Solo se rellenan los hemisferios cuyo decayID coincide con --decay.

Variaciones (peso aplicado a cada tau del decay elegido):
    SM       -> peso 1 (eventos KKMC no pesados)
    HelP1    -> peso 1, solo taus con genHelicity > 0
    HelM1    -> peso 1, solo taus con genHelicity < 0
    P1       -> x weight_P1   (single-tau, A_tau=+1, leido del arbol)
    M1       -> x weight_M1   (single-tau, A_tau=-1, leido del arbol)
    corr_P1  -> x peso joint dos-tau (Alcaraz eq.9/13/16), calculado al vuelo
    corr_M1  -> x peso joint dos-tau

USO:
    python RhoAnalysis/genHelicityHistos.py \\
        --tree-file Results/RhoAnalysis/.../tau_trainedAll_*.root \\
        -o Results/RhoAnalysis/genHelicityHistos_rho.root \\
        --decay 1
"""
import argparse
import math
import sys
from pathlib import Path

import ROOT

# Permite ejecutar como script desde la raiz del repo (modules.*)
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from modules import weightsPol
from modules.rhoTreeUtils import make_p4

# ── Configuracion de histogramas ────────────────────────────────────────────

VARIATIONS = ["SM", "HelP1", "HelM1", "P1", "M1", "corrP1", "corrM1"]
CHARGES    = ["minus", "plus"]   # minus = tau- (PDG 15), plus = tau+ (PDG -15)

# Etiqueta de decay para titulo/nombre de fichero
_DECAY_NAME = {0: "pion", 1: "rho", 10: "a1", -11: "ele", -13: "muon"}


def book_histograms(nbins):
    """Crea todos los TH1 y los devuelve en un dict {name: TH1D}."""
    h = {}
    for var in VARIATIONS:
        for ch in CHARGES:
            name = f"CosThetaTau_{var}_{ch}"
            h[name] = ROOT.TH1D(name, f"cos#theta_{{#tau}} {var} {ch};cos#theta_{{#tau}};Eventos",
                                nbins, -1.0, 1.0)
            for pre in ("GENFZ", "GENFZ_FW", "GENFZ_BW"):
                fname = f"{pre}_{var}_{ch}"
                h[fname] = ROOT.TH1D(fname, f"fz {pre} {var} {ch};f_{{z}};Eventos",
                                     nbins, 0.0, 0.5)
    return h


# ── Lectura del arbol ───────────────────────────────────────────────────────

_TAU_SCALARS = [
    "cos_theta_tau", "tauPDG", "genHelicity", "weight_P1", "weight_M1",
    "decayID", "omega", "P", "E", "Theta", "Phi",
    "visP", "visE", "visTheta", "visPhi",
]


def read_tau(entry, prefix):
    """Lee las ramas {prefix}_* en un dict."""
    return {k: float(getattr(entry, f"{prefix}_{k}")) for k in _TAU_SCALARS}


# ── Observable de spin H para el peso joint ─────────────────────────────────

def _H_for_joint(tv, beamE):
    """H de la formula joint (mismo criterio que el pipeline MDecs, use_omega=True,
    use_costheta_pion=True): rho->omega del arbol, pion->cos(theta*) exacto,
    a1->H_V, leptonico->H_ell. Otros -> 0."""
    decay_id = int(tv["decayID"])
    if decay_id == 1:                      # rho: variable optima omega (del arbol)
        return tv["omega"]
    tauP4 = make_p4(tv["P"],   tv["Theta"],   tv["Phi"],   tv["E"])
    visP4 = make_p4(tv["visP"], tv["visTheta"], tv["visPhi"], tv["visE"])
    if decay_id == 0:                      # pion: cos(theta*) geometrico exacto
        return weightsPol.cosThetaStar(tauP4, visP4)
    if decay_id == 10:                     # a1: H_V = alpha_V * z_R
        H = weightsPol._compute_H(visP4, tauP4, 10)
        return H if H is not None else 0.0
    if decay_id in (-11, -13):             # leptonico: H_ell
        return weightsPol._compute_H_lep(visP4, beamE)
    return 0.0


def compute_corr_weights(tau_m, tau_p, beamE, sin_eff):
    """Pesos joint dos-tau (corr_P1, corr_M1). El tau- es la referencia para
    z=cosTheta (Alcaraz eq.9). Devuelve (w_corr_P1, w_corr_M1)."""
    H_m = _H_for_joint(tau_m, beamE)   # tau-
    H_p = _H_for_joint(tau_p, beamE)   # tau+
    tauMinusP4 = make_p4(tau_m["P"], tau_m["Theta"], tau_m["Phi"], tau_m["E"])
    w_p1 = weightsPol.newAtauJoint(tauMinusP4, H_m, H_p, +1.0, tau_pdg=15, sin_eff=sin_eff)
    w_m1 = weightsPol.newAtauJoint(tauMinusP4, H_m, H_p, -1.0, tau_pdg=15, sin_eff=sin_eff)
    return w_p1, w_m1


# ── Llenado ─────────────────────────────────────────────────────────────────

def fill_tau(h, tv, w_corr_p1, w_corr_m1):
    """Rellena las 7 variaciones x (cosTheta + fz/FW/BW) para un tau."""
    charge = "minus" if int(tv["tauPDG"]) == 15 else "plus"
    z   = tv["cos_theta_tau"]
    absz = abs(z)
    fz  = absz / (1.0 + absz * absz)
    hel = tv["genHelicity"]

    # (variacion, peso, condicion_extra)
    entries = [
        ("SM",     1.0,               True),
        ("HelP1",  1.0,               hel > 0),
        ("HelM1",  1.0,               hel < 0),
        ("P1",     tv["weight_P1"],   True),
        ("M1",     tv["weight_M1"],   True),
        ("corrP1", w_corr_p1,         True),
        ("corrM1", w_corr_m1,         True),
    ]

    for var, w, cond in entries:
        if not cond:
            continue
        h[f"CosThetaTau_{var}_{charge}"].Fill(z, w)
        h[f"GENFZ_{var}_{charge}"].Fill(fz, w)
        if z > 0:
            h[f"GENFZ_FW_{var}_{charge}"].Fill(fz, w)
        else:
            h[f"GENFZ_BW_{var}_{charge}"].Fill(fz, w)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree-file", required=True, help="Rootfile con el TTree gen-only")
    ap.add_argument("-o", "--output", required=True, help="Rootfile de salida con los histogramas")
    ap.add_argument("--tree-name", default="outtree_original", help="Nombre del TTree (default: outtree_original)")
    ap.add_argument("--decay", type=int, default=1,
                    help="decayID a rellenar: 0=pion, 1=rho, 10=a1, -11=ele, -13=muon (default: 1)")
    ap.add_argument("--sin-eff", type=float, default=0.2312, help="sin^2(theta_eff) para los pesos (default: 0.2312)")
    ap.add_argument("--nbins", type=int, default=50, help="Numero de bins (default: 50)")
    ap.add_argument("--max-entries", type=int, default=-1, help="Limita entradas procesadas (debug)")
    args = ap.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.SetDefaultSumw2(True)

    infile = ROOT.TFile.Open(args.tree_file)
    if not infile or infile.IsZombie():
        sys.exit(f"ERROR: no se pudo abrir {args.tree_file}")
    tree = infile.Get(args.tree_name)
    if not isinstance(tree, ROOT.TTree):
        sys.exit(f"ERROR: no se encontro el TTree '{args.tree_name}' en {args.tree_file}")

    n_entries = tree.GetEntries()
    if args.max_entries > 0:
        n_entries = min(n_entries, args.max_entries)
    dname = _DECAY_NAME.get(args.decay, str(args.decay))
    print(f"[genHelicityHistos] {n_entries} entradas | decay={args.decay} ({dname}) | sin_eff={args.sin_eff}")

    h = book_histograms(args.nbins)

    n_filled = 0
    for i in range(n_entries):
        tree.GetEntry(i)
        tau1 = read_tau(tree, "tau1")
        tau2 = read_tau(tree, "tau2")
        beamE = float(getattr(tree, "beamE", 0.0))

        # Identifica tau-/tau+ para el peso joint (referencia z = tau-)
        if int(tau1["tauPDG"]) == 15:
            tau_m, tau_p = tau1, tau2
        else:
            tau_m, tau_p = tau2, tau1
        w_corr_p1, w_corr_m1 = compute_corr_weights(tau_m, tau_p, beamE, args.sin_eff)

        # Rellena solo los hemisferios cuyo decay coincide con --decay
        for tv in (tau1, tau2):
            if int(tv["decayID"]) == args.decay:
                fill_tau(h, tv, w_corr_p1, w_corr_m1)
                n_filled += 1

        if (i + 1) % 200000 == 0:
            print(f"  ... {i + 1}/{n_entries} entradas")

    print(f"[genHelicityHistos] hemisferios rellenados (decay={args.decay}): {n_filled}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    outfile = ROOT.TFile(args.output, "RECREATE")
    for hist in h.values():
        hist.Write()
    outfile.Close()
    infile.Close()
    print(f"[genHelicityHistos] {len(h)} histogramas guardados en {args.output}")


if __name__ == "__main__":
    main()
