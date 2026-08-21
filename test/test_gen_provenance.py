"""Unit tests for the gen-level provenance helpers of ``modules.tauReco``.

These run on stub MCParticles, so they need ROOT but not podio/edm4hep: they can
be executed outside the key4hep stack with ``python test/test_gen_provenance.py``.
They cover the three things the tagging has to get right:

  * neutral products that no decay-mode counter registers (K0_L, n, Lambda…),
  * where a gen photon comes from (pi0 vs radiation),
  * primary taus vs radiative ``tau -> gamma -> tau tau`` chains, including the
    ``e+e- -> gamma* -> tau tau`` case, which also hangs from a photon and must
    NOT be flagged as secondary.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules import tauReco
from modules.tauReco import (
    PHOTON_ORIGIN_CHARGED_RAD,
    PHOTON_ORIGIN_NOT_A_PHOTON,
    PHOTON_ORIGIN_OTHER,
    PHOTON_ORIGIN_PI0,
    PHOTON_ORIGIN_SIMULATION,
    PHOTON_ORIGIN_TAU_FSR,
)


# ── MCParticle stubs ──────────────────────────────────────────────────────────

class _ObjectID:
    def __init__(self, index):
        self.index = index


class _Momentum:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class FakeMC:
    """Minimal stand-in for an edm4hep MCParticle."""

    _next_index = [0]

    def __init__(self, pdg, status=1, charge=0.0, mass=0.0, momentum=(0.0, 0.0, 1.0)):
        self.pdg = pdg
        self.status = status
        self.charge = charge
        self.mass = mass
        self.momentum = _Momentum(*momentum)
        self.parents = []
        self.daughters = []
        self.index = FakeMC._next_index[0]
        FakeMC._next_index[0] += 1

    def child(self, other):
        """Attach *other* as a daughter of self and return it."""
        self.daughters.append(other)
        other.parents.append(self)
        return other

    def getPDG(self):
        return self.pdg

    def getGeneratorStatus(self):
        return self.status

    def getCharge(self):
        return self.charge

    def getMass(self):
        return self.mass

    def getMomentum(self):
        return self.momentum

    def getParents(self):
        return self.parents

    def getDaughters(self):
        return self.daughters

    def getObjectID(self):
        return _ObjectID(self.index)

    def isCreatedInSimulation(self):
        return self.status == 0


def _tau_decay(tau, with_extra_neutral=None):
    """Give *tau* a simple 1-prong decay, optionally plus an extra neutral."""
    tau.child(FakeMC(-211 if tau.pdg > 0 else 211, status=1, charge=-1.0 if tau.pdg > 0 else 1.0,
                     mass=0.13957, momentum=(0.0, 0.0, 10.0)))
    tau.child(FakeMC(16 if tau.pdg > 0 else -16, status=1, momentum=(0.0, 0.0, 2.0)))
    if with_extra_neutral is not None:
        tau.child(FakeMC(with_extra_neutral, status=1, charge=0.0,
                         mass=0.497, momentum=(0.0, 1.0, 5.0)))
    return tau


# ── Scenario builders ─────────────────────────────────────────────────────────

def build_primary_z_event():
    """e+e- -> Z -> tau tau, with one copy of each tau before the decay."""
    beam = FakeMC(11, status=4, charge=-1.0)
    z = beam.child(FakeMC(23, status=3))
    taus = []
    for pdg in (15, -15):
        copy = z.child(FakeMC(pdg, status=3, charge=-1.0 if pdg > 0 else 1.0))
        final = copy.child(FakeMC(pdg, status=2, charge=-1.0 if pdg > 0 else 1.0,
                                  mass=1.777, momentum=(0.0, 0.0, 40.0)))
        taus.append(_tau_decay(final))
    return beam, taus


def build_gamma_star_event():
    """e+e- -> gamma* -> tau tau: primary taus that hang from a photon."""
    beam = FakeMC(11, status=4, charge=-1.0)
    gamma = beam.child(FakeMC(22, status=3))
    taus = []
    for pdg in (15, -15):
        final = gamma.child(FakeMC(pdg, status=2, charge=-1.0 if pdg > 0 else 1.0,
                                   mass=1.777, momentum=(0.0, 0.0, 40.0)))
        taus.append(_tau_decay(final))
    return beam, taus


def build_secondary_tau_event():
    """tau -> gamma -> tau tau on top of a primary Z pair."""
    beam, primary_taus = build_primary_z_event()
    radiator = primary_taus[0]
    gamma = radiator.child(FakeMC(22, status=3, momentum=(0.0, 0.0, 20.0)))
    secondaries = []
    for pdg in (15, -15):
        final = gamma.child(FakeMC(pdg, status=2, charge=-1.0 if pdg > 0 else 1.0,
                                   mass=1.777, momentum=(0.0, 0.0, 8.0)))
        secondaries.append(_tau_decay(final))
    return beam, primary_taus, secondaries, gamma


def flatten(root):
    """Depth-first list of the whole decay tree, parents before daughters."""
    out, stack, seen = [], [root], set()
    while stack:
        cur = stack.pop(0)
        if cur.index in seen:
            continue
        seen.add(cur.index)
        out.append(cur)
        stack.extend(cur.getDaughters())
    return out


# ── Tests ─────────────────────────────────────────────────────────────────────

def test_is_extra_neutral():
    assert tauReco.is_extra_neutral(FakeMC(130, charge=0.0)), "K0_L must be tagged"
    assert tauReco.is_extra_neutral(FakeMC(2112, charge=0.0)), "neutron must be tagged"
    assert tauReco.is_extra_neutral(FakeMC(3122, charge=0.0)), "Lambda must be tagged"
    assert not tauReco.is_extra_neutral(FakeMC(22, charge=0.0)), "photons are tagged apart"
    assert not tauReco.is_extra_neutral(FakeMC(111, charge=0.0)), "pi0 is counted"
    assert not tauReco.is_extra_neutral(FakeMC(16, charge=0.0)), "neutrinos are skipped"
    assert not tauReco.is_extra_neutral(FakeMC(211, charge=1.0)), "charged never tagged"
    print("ok  is_extra_neutral")


def test_photon_origin():
    tau = FakeMC(15, status=2, charge=-1.0)
    pi0 = tau.child(FakeMC(111, status=2))
    g_pi0 = pi0.child(FakeMC(22, status=1))
    g_fsr = tau.child(FakeMC(22, status=1))
    pion = tau.child(FakeMC(-211, status=1, charge=-1.0))
    g_pion = pion.child(FakeMC(22, status=1))
    g_sim = pion.child(FakeMC(22, status=0))

    assert tauReco.classify_photon_origin(g_pi0) == PHOTON_ORIGIN_PI0
    assert tauReco.classify_photon_origin(g_fsr) == PHOTON_ORIGIN_TAU_FSR
    assert tauReco.classify_photon_origin(g_pion) == PHOTON_ORIGIN_CHARGED_RAD
    assert tauReco.classify_photon_origin(g_sim) == PHOTON_ORIGIN_SIMULATION
    assert tauReco.classify_photon_origin(pion) == PHOTON_ORIGIN_NOT_A_PHOTON

    # El ancestro agrupa los dos gammas del mismo pi0.
    g_pi0_b = pi0.child(FakeMC(22, status=1))
    _, parent_pdg, ancestor_idx = tauReco.photon_origin_info(g_pi0)
    _, parent_pdg_b, ancestor_idx_b = tauReco.photon_origin_info(g_pi0_b)
    assert parent_pdg == parent_pdg_b == 111
    assert ancestor_idx == ancestor_idx_b == pi0.index

    # Un fotón sin ancestros conocidos no revienta.
    assert tauReco.classify_photon_origin(FakeMC(22, status=1)) == PHOTON_ORIGIN_OTHER
    print("ok  photon origin")


def test_tau_origin_primary_z():
    _, taus = build_primary_z_event()
    for tau in taus:
        origin_pdg, is_secondary, mother_idx, photon_idx = tauReco.classify_tau_origin(tau)
        assert origin_pdg == 23, origin_pdg
        assert not is_secondary
        assert mother_idx == -1 and photon_idx == -1
    print("ok  primary tau (Z)")


def test_tau_origin_gamma_star():
    """The crux: a primary tau from gamma* also has a photon parent."""
    _, taus = build_gamma_star_event()
    for tau in taus:
        origin_pdg, is_secondary, mother_idx, photon_idx = tauReco.classify_tau_origin(tau)
        assert origin_pdg == 22, origin_pdg
        assert not is_secondary, "gamma* taus must NOT be flagged as secondary"
        assert mother_idx == -1 and photon_idx == -1
    print("ok  primary tau (gamma*)")


def test_tau_origin_secondary():
    _, primary_taus, secondaries, gamma = build_secondary_tau_event()
    mother = primary_taus[0]
    for tau in secondaries:
        origin_pdg, is_secondary, mother_idx, photon_idx = tauReco.classify_tau_origin(tau)
        assert origin_pdg == 22, origin_pdg
        assert is_secondary, "tau -> gamma -> tau tau must be flagged"
        assert mother_idx == mother.index, (mother_idx, mother.index)
        assert photon_idx == gamma.index, (photon_idx, gamma.index)
    print("ok  secondary tau")


def test_find_all_gen_taus_links():
    """findAllGenTaus must resolve the mother through the copy chain."""
    beam, primary_taus, secondaries, gamma = build_secondary_tau_event()
    gen_taus = tauReco.findAllGenTaus(flatten(beam))

    assert len(gen_taus) == 4, f"expected 4 status-2 taus, got {len(gen_taus)}"

    by_mc_idx = {gen_taus[k].getMCIdx(): k for k in gen_taus}
    mother_key = by_mc_idx[primary_taus[0].index]

    for tau in secondaries:
        key = by_mc_idx[tau.index]
        gen = gen_taus[key]
        assert gen.getIsSecondary()
        assert gen.getOriginPDG() == 22
        assert gen.getMotherTauKey() == mother_key, (gen.getMotherTauKey(), mother_key)
        assert gen.getMotherTauMCIdx() == primary_taus[0].index
        assert gen.getRadPhotonMCIdx() == gamma.index

    for tau in primary_taus:
        gen = gen_taus[by_mc_idx[tau.index]]
        assert not gen.getIsSecondary()
        assert gen.getMotherTauKey() == -1
    print("ok  findAllGenTaus mother links")


def test_extra_neutrals_do_not_change_the_id():
    """A tau -> pi K0_L nu keeps ID 0 but comes out flagged."""
    plain = _tau_decay(FakeMC(15, status=2, charge=-1.0, mass=1.777,
                              momentum=(0.0, 0.0, 40.0)))
    with_k0l = _tau_decay(FakeMC(15, status=2, charge=-1.0, mass=1.777,
                                 momentum=(0.0, 0.0, 40.0)), with_extra_neutral=130)

    plain_data = tauReco.visTauGen(plain)
    k0l_data = tauReco.visTauGen(with_k0l)

    assert plain_data["ID"] == 0, plain_data["ID"]
    assert k0l_data["ID"] == 0, "the decay-mode encoding must not change"
    assert not plain_data["hasExtraNeutrals"]
    assert k0l_data["hasExtraNeutrals"], "K0_L must raise the flag"
    assert [p.getPDG() for p in k0l_data["extraNeutrals"].values()] == [130]
    # Sigue formando parte de const y del momento visible, como antes.
    assert 130 in [p.getPDG() for p in k0l_data["const"].values()]
    assert k0l_data["visP4"].P() > plain_data["visP4"].P()
    print("ok  extra neutrals keep the decay ID")


def test_const_origin_is_filled():
    tau = FakeMC(15, status=2, charge=-1.0, mass=1.777, momentum=(0.0, 0.0, 40.0))
    tau.child(FakeMC(-211, status=1, charge=-1.0, mass=0.13957, momentum=(0.0, 0.0, 10.0)))
    tau.child(FakeMC(16, status=1, momentum=(0.0, 0.0, 2.0)))
    tau.child(FakeMC(22, status=1, momentum=(0.0, 0.5, 3.0)))

    data = tauReco.visTauGen(tau)
    origins = data["constOrigin"]
    consts = data["const"]
    assert len(origins) == len(consts)
    for key, part in consts.items():
        if abs(part.getPDG()) == 22:
            assert origins[key] == PHOTON_ORIGIN_TAU_FSR, origins[key]
        else:
            assert origins[key] == PHOTON_ORIGIN_NOT_A_PHOTON
    print("ok  constOrigin")


if __name__ == "__main__":
    test_is_extra_neutral()
    test_photon_origin()
    test_tau_origin_primary_z()
    test_tau_origin_gamma_star()
    test_tau_origin_secondary()
    test_find_all_gen_taus_links()
    test_extra_neutrals_do_not_change_the_id()
    test_const_origin_is_filled()
    print("\nAll gen-provenance tests passed.")
