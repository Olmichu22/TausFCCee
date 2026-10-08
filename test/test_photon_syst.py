"""Comprobaciones de las variaciones sistemáticas de fotones de TTreesTausLong.

Ejecutar (con key4hep cargado):  python test/test_photon_syst.py
"""
import math
import os
import sys

import numpy as np
import ROOT

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules import tauReco
from TauAnalysis.TTreesTausLong import (
    VariedParticle, photon_variations, vary_photons,
)


class FakePFO:
    """PFO minimo con la interfaz que usa tauReco (getMomentum/getMass/getPDG/getCharge)."""

    def __init__(self, pdg, p, theta, phi):
        self.pdg = pdg
        self.mass = 0.13957 if abs(pdg) == 211 else 0.0
        self.charge = math.copysign(1, pdg) if abs(pdg) == 211 else 0
        self.p4 = ROOT.TLorentzVector()
        self.p4.SetXYZM(p * math.sin(theta) * math.cos(phi),
                        p * math.sin(theta) * math.sin(phi),
                        p * math.cos(theta), self.mass)

    def getMomentum(self):
        return self.p4

    def getMass(self):
        return self.mass

    def getPDG(self):
        return self.pdg

    def getCharge(self):
        return self.charge


def rho_event():
    """pi- + dos fotones de un pi0 dentro del cono, y un pi+ en el hemisferio opuesto."""
    return [
        FakePFO(-211, 20.0, 1.2, 0.3),
        FakePFO(22, 8.0, 1.25, 0.35),
        FakePFO(22, 5.0, 1.15, 0.25),
        FakePFO(211, 30.0, math.pi - 1.2, 0.3 + math.pi),
    ]


def find(parts):
    return tauReco.findAllTaus(parts, 0.4, 0.0, 0.0, 0.0, 0.0)


def check_validation():
    assert photon_variations(None) == []
    assert photon_variations({"energy": None}) == []
    vs = photon_variations({"energy": {"factor": 0.1, "const": 0.01},
                            "direction": {"sigma_theta": 0.02}})
    assert [v[0] for v in vs] == ["PhEnUp", "PhEnDown", "PhDir"], vs
    for bad in ({"energy": {"const": 0.01}},                      # falta factor
                {"energy": {"factor": 0, "const": 0.01, "x": 1}},  # clave desconocida
                {"direction": {"const_theta": 0.01}},              # claves antiguas
                {"direction": {"sigma_theta": -1}},                # negativo
                {"energy": {"factor": "a", "const": 0}},           # no numérico
                {"pion": {}},                                      # sección desconocida
                {"energy": {"factor": 0, "const": 0.01, "name": "a"}},  # name solo en listas
                {"direction": [{"sigma_theta": 0.01}, {"sigma_theta": 0.01}]},  # sufijo repetido
                {"direction": [{"sigma_theta": 0.01, "name": "a"},
                               {"sigma_theta": 0.02, "name": "a"}]},  # name repetido
                {"direction": [{"sigma_theta": 0.01, "name": "a-b"}]},  # name no válido
                {"direction": [0.01]}):                            # punto que no es un dict
        try:
            photon_variations(bad)
        except ValueError:
            continue
        raise AssertionError(f"{bad} should have been rejected")
    print("validation OK")


def check_scan():
    assert photon_variations({"energy": []}) == []
    vs = photon_variations({
        "energy": [{"factor": 0, "const": 0.005},
                   {"factor": 0.1, "const": 0.01},
                   {"factor": 0, "const": 0.02, "name": "c2pct"}],
        "direction": [{"sigma_theta": 0.001}, {"sigma_theta": 1e-5}],
    })
    assert [v[0] for v in vs] == [
        "PhEnUp_c0p005", "PhEnDown_c0p005",
        "PhEnUp_f0p1_c0p01", "PhEnDown_f0p1_c0p01",
        "PhEnUp_c2pct", "PhEnDown_c2pct",
        "PhDir_s0p001", "PhDir_s1em05",
    ], vs
    # El name no llega al cfg que reciben las funciones de tauReco
    assert all("name" not in v[2] for v in vs)
    assert vs[4][2] == {"factor": 0.0, "const": 0.02}
    # Una sola semilla: sufijos sin tocar; la energía no lleva semilla
    assert [v[3] for v in vs] == [None] * 6 + [12345, 12345]

    # Varias semillas: solo se repite la dirección
    cfg = {"energy": [{"factor": 0, "const": 0.01}],
           "direction": [{"sigma_theta": 0.001}, {"sigma_theta": 0.01}]}
    vs = photon_variations(cfg, [1, 2])
    assert [(v[0], v[3]) for v in vs] == [
        ("PhEnUp_c0p01", None), ("PhEnDown_c0p01", None),
        ("PhDir_s0p001_seed1", 1), ("PhDir_s0p001_seed2", 2),
        ("PhDir_s0p01_seed1", 1), ("PhDir_s0p01_seed2", 2),
    ], vs
    vs = photon_variations({"direction": {"sigma_theta": 0.01}}, [1, 2])
    assert [v[0] for v in vs] == ["PhDir_seed1", "PhDir_seed2"], vs
    for bad_seeds in ([], [1, 1], [-1], [1.5]):
        try:
            photon_variations(cfg, bad_seeds)
        except ValueError:
            continue
        raise AssertionError(f"seeds {bad_seeds} should have been rejected")

    # Misma semilla para dos sigma_theta: mismo eje por fotón y ángulo ∝ sigma
    parts = rho_event()
    small = vary_photons(parts, "direction", {"sigma_theta": 0.001}, np.random.default_rng([7, 42]))
    large = vary_photons(parts, "direction", {"sigma_theta": 0.002}, np.random.default_rng([7, 42]))
    n_checked = 0
    for o, s, l in zip(parts, small, large):
        if abs(o.getPDG()) != 22:
            continue
        v0 = o.getMomentum().Vect()
        vs_, vl = s.getMomentum().Vect(), l.getMomentum().Vect()
        assert abs(vl.Angle(v0) / vs_.Angle(v0) - 2.0) < 1e-3
        # Desplazamientos transversos paralelos: mismo lado del cono
        ds, dl = vs_ - v0, vl - v0
        assert ds.Angle(dl) < 1e-2, ds.Angle(dl)
        n_checked += 1
    assert n_checked == 2
    # Semillas distintas: realizaciones distintas del mismo sigma_theta
    other = vary_photons(parts, "direction", {"sigma_theta": 0.001}, np.random.default_rng([8, 42]))
    assert any(s.getMomentum() != o.getMomentum()
               for p, s, o in zip(parts, small, other) if abs(p.getPDG()) == 22)
    print("scan OK")


def check_varied_particle():
    orig = FakePFO(22, 10.0, 1.0, 0.0)
    new_p4 = ROOT.TLorentzVector(1, 2, 3, math.sqrt(14))
    vp = VariedParticle(orig, new_p4)
    assert vp.getMomentum() is new_p4
    assert vp.getMass() == 0.0 and vp.getPDG() == 22 and vp.getCharge() == 0
    assert vp == vp and not (vp == orig) and vp != orig
    print("VariedParticle OK")


def check_energy_scale():
    parts = rho_event()
    cfg = {"factor": 0.0, "const": 0.01}
    up = vary_photons(parts, "energy_up", cfg, None)
    down = vary_photons(parts, "energy_down", cfg, None)
    for o, u, d in zip(parts, up, down):
        if abs(o.getPDG()) != 22:
            assert u is o and d is o
            continue
        p = o.getMomentum().P()
        assert abs(u.getMomentum().P() - 1.01 * p) < 1e-9 * p
        assert abs(d.getMomentum().P() - 0.99 * p) < 1e-9 * p
        assert abs(u.getMomentum().Angle(o.getMomentum().Vect())) < 1e-12

    nom, t_up, t_down = find(parts), find(up), find(down)
    photons_E = 13.0
    for t_nom, t_var, sign in ((nom, t_up, +1), (nom, t_down, -1)):
        assert len(t_nom) == len(t_var)
        assert t_nom[0].getID() == t_var[0].getID()
        dE = t_var[0].getMomentum().E() - t_nom[0].getMomentum().E()
        assert abs(dE - sign * 0.01 * photons_E) < 1e-6, dE
    print("energy scale OK")


def check_closure_and_direction():
    parts = rho_event()
    nom = find(parts)
    for kind, cfg in (("energy_up", {"factor": 0.0, "const": 0.0}),
                      ("direction", {"sigma_theta": 0.0})):
        var = find(vary_photons(parts, kind, cfg, np.random.default_rng(1)))
        assert len(var) == len(nom)
        for k in nom:
            assert var[k].getID() == nom[k].getID()
            assert (var[k].getMomentum() - nom[k].getMomentum()).P() < 1e-6
    # Smearing reproducible con la misma semilla y |p| conservado
    a = vary_photons(parts, "direction", {"sigma_theta": 0.01}, np.random.default_rng([7, 42]))
    b = vary_photons(parts, "direction", {"sigma_theta": 0.01}, np.random.default_rng([7, 42]))
    for o, x, y in zip(parts, a, b):
        assert x.getMomentum() == y.getMomentum()
        assert abs(x.getMomentum().P() - o.getMomentum().P()) < 1e-6
    print("closure + direction OK")


def check_cone_axis():
    # pion blando; fotón duro a dR=0.3 (arrastra el eje "running" hacia phi>0);
    # fotón blando a dR=0.35 en el lado opuesto: dentro con eje fijo en el pión,
    # fuera con el eje acumulado.
    parts = [
        FakePFO(-211, 5.0, 1.2, 0.0),
        FakePFO(22, 20.0, 1.2, 0.30),
        FakePFO(22, 1.0, 1.2, -0.35),
    ]
    running = tauReco.findAllTaus(parts, 0.4, 0.0, 0.0, 0.0, 0.0)
    default = tauReco.findAllTaus(parts, 0.4, 0.0, 0.0, 0.0, 0.0, cone_axis="running")
    lead = tauReco.findAllTaus(parts, 0.4, 0.0, 0.0, 0.0, 0.0, cone_axis="lead")
    assert running[0].getnConst() == default[0].getnConst() == 2, running[0].getnConst()
    assert lead[0].getnConst() == 3, lead[0].getnConst()
    # El 4-momento del pión no se toca en ningún modo
    assert abs(parts[0].getMomentum().P() - 5.0) < 1e-9
    # "lead": el maxCone es la distancia al pión del constituyente más lejano
    assert abs(lead[0].getMaxCone() - 0.35) < 1e-9, lead[0].getMaxCone()
    try:
        tauReco.findAllTaus(parts, 0.4, 0.0, 0.0, 0.0, 0.0, cone_axis="pion")
    except ValueError:
        pass
    else:
        raise AssertionError("invalid cone_axis should raise")
    print("cone axis OK")


if __name__ == "__main__":
    check_validation()
    check_scan()
    check_varied_particle()
    check_energy_scale()
    check_closure_and_direction()
    check_cone_axis()
    print("ALL OK")
