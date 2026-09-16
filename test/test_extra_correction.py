"""Comprobaciones de tauReco.extraTauRecoCorrection con eventos sinteticos.

Ejecutar (con key4hep cargado):  python test/test_extra_correction.py
"""
import math
import os
import sys

import ROOT

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules import tauReco


class FakePFO:
    """PFO minimo con la interfaz que usa tauReco (getMomentum/getMass/getPDG/getCharge)."""

    def __init__(self, pdg, p, theta, phi, mass=None, charge=None):
        if mass is None:
            mass = 0.13957 if abs(pdg) == 211 else 0.0
        if charge is None:
            charge = math.copysign(1, pdg) if abs(pdg) == 211 else 0
        self.pdg = pdg
        self.mass = mass
        self.charge = charge
        self.p4 = ROOT.TLorentzVector()
        self.p4.SetXYZM(p * math.sin(theta) * math.cos(phi),
                        p * math.sin(theta) * math.sin(phi),
                        p * math.cos(theta), mass)

    def getMomentum(self):
        return self.p4

    def getMass(self):
        return self.mass

    def getPDG(self):
        return self.pdg

    def getCharge(self):
        return self.charge


CFG = {"enable": True, "modes": ["pion_photon_fsr"]}

THETA = math.pi / 2


def reco_id(pfos, extra_correction=None):
    taus = tauReco.findAllTaus(pfos, 0.4, 0., 0., 0., 0.,
                               extra_correction=extra_correction)
    assert len(taus) == 1, f"esperado 1 tau, obtenidos {len(taus)}"
    return taus[0].getID(), taus[0].getMomentum().P()


def case(name, pfos, id_base, id_corr):
    got_base, p_base = reco_id(pfos)
    got_corr, p_corr = reco_id(pfos, CFG)
    ok = (got_base == id_base) and (got_corr == id_corr)
    print(f"[{'OK ' if ok else 'FAIL'}] {name}: base ID {got_base} (P {p_base:.2f}) "
          f"-> corregido ID {got_corr} (P {p_corr:.2f}); esperado {id_base} -> {id_corr}")
    return ok


def main():
    ok = True

    # 1) FSR duro de la linea del tau: pion 20 GeV + foton 15 GeV a dR 0.2.
    #    m(pi+gamma) ~ 3.4 GeV, muy por encima de m_rho -> se quita.
    ok &= case("FSR duro",
               [FakePFO(211, 20., THETA, 0.), FakePFO(22, 15., THETA + 0.2, 0.)],
               1, 0)

    # 2) Fragmento del shower del pion: foton de 0.5 GeV pegado al pion.
    ok &= case("shower blando",
               [FakePFO(211, 20., THETA, 0.), FakePFO(22, 0.5, THETA + 0.04, 0.)],
               1, 0)

    # 3) rho: los dos fotones dan la masa del pi0 -> intocables.
    #    2 fotones de 5 GeV separados 0.027 rad dan m_gg = 0.135 GeV.
    dgamma = 2 * math.asin(0.1349768 / (2 * math.sqrt(5. * 5.)))
    ok &= case("rho (pi0 real)",
               [FakePFO(211, 20., THETA, 0.),
                FakePFO(22, 5., THETA + 0.1, 0.),
                FakePFO(22, 5., THETA + 0.1 + dgamma, 0.)],
               2, 2)

    # 4) pi0 asimetrico con un foton perdido: el superviviente es duro pero
    #    m(pi+gamma) se queda por debajo de m_rho -> no se toca.
    ok &= case("pi0 con foton perdido",
               [FakePFO(211, 20., THETA, 0.), FakePFO(22, 6., THETA + 0.035, 0.)],
               1, 1)

    # 5) Sin configuracion no cambia nada (regresion de la reconstruccion base).
    pfos = [FakePFO(211, 20., THETA, 0.), FakePFO(22, 15., THETA + 0.2, 0.)]
    id_none, _ = reco_id(pfos, None)
    id_disabled, _ = reco_id(pfos, {"enable": False, "modes": ["pion_photon_fsr"]})
    ok_cfg = id_none == 1 and id_disabled == 1
    print(f"[{'OK ' if ok_cfg else 'FAIL'}] sin correccion / desactivada: {id_none}, {id_disabled} (esperado 1, 1)")
    ok &= ok_cfg

    # 6) Modo desconocido -> error explicito.
    try:
        reco_id(pfos, {"modes": ["no_existe"]})
        print("[FAIL] modo desconocido: no lanzo excepcion")
        ok = False
    except KeyError as exc:
        print(f"[OK ] modo desconocido lanza KeyError: {str(exc)[:60]}...")

    print("\nRESULTADO:", "todo OK" if ok else "HAY FALLOS")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
