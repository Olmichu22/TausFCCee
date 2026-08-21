"""Unit tests for ``modules.genDecayModes``, the true gen-level decay label.

Like ``test_gen_provenance.py`` these run on stub MCParticles, so they need
ROOT but not podio/edm4hep and can be run with
``python test/test_gen_decay_modes.py``.

They pin down the four canonicalisation rules that make the label usable as a
grouping key: neutrinos and tau FSR dropped, tau+ folded onto tau-, stable
ordering, and the fact that the label never disturbs the visible-topology ID.
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)

from modules import genDecayModes, tauReco
from test_gen_provenance import FakeMC


def _tau(charge=-1.0):
    return FakeMC(15 if charge < 0 else -15, status=2, charge=charge, mass=1.77682)


def test_neutrinos_and_fsr_dropped():
    """nu_tau and photons hanging off the tau must not reach the label."""
    tau = _tau()
    tau.child(FakeMC(16, status=1))                      # nu_tau
    tau.child(FakeMC(22, status=1))                      # FSR del tau
    tau.child(FakeMC(-211, status=1, charge=-1.0))
    tau.child(FakeMC(111, status=2))

    pdgs = genDecayModes.canonical_daughters(tau)
    assert pdgs == [-211, 111], pdgs
    assert genDecayModes.decay_label(pdgs) == "-211,111"
    assert genDecayModes.true_mode(pdgs) == genDecayModes.MODE_PI_PI0
    print("ok  neutrinos and FSR dropped")


def test_leptonic_collapses_to_single_label():
    """Dropping every neutrino makes tau -> l nu nu a one-entry label."""
    tau = _tau()
    tau.child(FakeMC(16, status=1))
    tau.child(FakeMC(-12, status=1))
    tau.child(FakeMC(11, status=1, charge=-1.0))

    pdgs = genDecayModes.canonical_daughters(tau)
    assert pdgs == [11], pdgs
    assert genDecayModes.true_mode(pdgs) == genDecayModes.MODE_E
    print("ok  leptonic label")


def test_charge_folding_keeps_self_conjugates():
    """tau+ folds onto the tau- convention; the omega is left alone."""
    minus = _tau(charge=-1.0)
    minus.child(FakeMC(16, status=1))
    minus.child(FakeMC(223, status=2))                   # omega
    minus.child(FakeMC(-321, status=1, charge=-1.0))     # K-

    plus = _tau(charge=+1.0)
    plus.child(FakeMC(-16, status=1))
    plus.child(FakeMC(223, status=2))
    plus.child(FakeMC(321, status=1, charge=+1.0))       # K+

    label_minus = genDecayModes.decay_label(genDecayModes.canonical_daughters(minus))
    label_plus = genDecayModes.decay_label(genDecayModes.canonical_daughters(plus))

    assert label_minus == "-321,223", label_minus
    assert label_plus == label_minus, (label_plus, label_minus)
    assert genDecayModes.true_mode(label_plus) == genDecayModes.MODE_OMEGA
    print("ok  charge folding")


def test_ordering_is_canonical():
    """The same multiset in a different order gives the same label."""
    first, second = _tau(), _tau()
    for pdg in (16, 111, -211, 111):
        first.child(FakeMC(pdg, status=1 if pdg != 111 else 2))
    for pdg in (111, 16, 111, -211):
        second.child(FakeMC(pdg, status=1 if pdg != 111 else 2))

    label = genDecayModes.decay_label(genDecayModes.canonical_daughters(first))
    assert label == genDecayModes.decay_label(genDecayModes.canonical_daughters(second))
    assert label == "-211,111,111", label
    print("ok  canonical ordering")


def test_simulation_secondaries_ignored():
    """status 0 products come from the detector sim, not from the decay."""
    tau = _tau()
    tau.child(FakeMC(16, status=1))
    tau.child(FakeMC(-211, status=1, charge=-1.0))
    tau.child(FakeMC(2112, status=0))                    # neutrón de simulación

    assert genDecayModes.canonical_daughters(tau) == [-211]
    print("ok  simulation secondaries ignored")


def test_unknown_label_is_flagged():
    """Anything outside MODE_TABLE must fall back to MODE_UNKNOWN, not crash."""
    tau = _tau()
    tau.child(FakeMC(16, status=1))
    tau.child(FakeMC(-4122, status=2, charge=-1.0))      # inventado

    mode = genDecayModes.true_mode(genDecayModes.canonical_daughters(tau))
    assert mode == genDecayModes.MODE_UNKNOWN, mode
    assert genDecayModes.mode_name(mode) == "unknown"
    print("ok  unknown label")


def test_visTauGen_exposes_mode_without_touching_id():
    """The K0_L case: ID stays 0 (1-prong topology) while the mode says K0 pi."""
    tau = _tau()
    tau.child(FakeMC(16, status=1))
    tau.child(FakeMC(-211, status=1, charge=-1.0, mass=0.13957))
    tau.child(FakeMC(-311, status=2, mass=0.49761)).child(
        FakeMC(130, status=1, mass=0.49761))

    data = tauReco.visTauGen(tau)

    assert data["ID"] == 0, data["ID"]
    assert data["hasExtraNeutrals"] is True
    assert data["decayDaughterPDG"] == [-311, -211], data["decayDaughterPDG"]
    assert data["trueMode"] == genDecayModes.MODE_K0_PI, data["trueMode"]
    print("ok  visTauGen exposes the mode without touching the ID")


def test_mode_table_labels_are_canonical():
    """Every key of MODE_TABLE must already be in canonical order."""
    for label in genDecayModes.MODE_TABLE:
        pdgs = [int(p) for p in label.split(",")]
        ordered = sorted(pdgs, key=lambda p: (abs(p), p), reverse=True)
        assert pdgs == ordered, label
    print(f"ok  {len(genDecayModes.MODE_TABLE)} MODE_TABLE labels are canonical")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("\nall gen decay mode tests passed")
