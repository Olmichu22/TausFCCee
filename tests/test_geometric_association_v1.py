import math
import unittest

import ROOT

from modules.NeutralRecover import get_reco_mc_links_by_dR
from modules.myutils import dRAngle


class _Vector:
    def __init__(self, theta, phi=0.0):
        self.x = math.sin(theta) * math.cos(phi)
        self.y = math.sin(theta) * math.sin(phi)
        self.z = math.cos(theta)


class _ObjectID:
    def __init__(self, index): self.index = index


class _Particle:
    def __init__(self, index, pdg, theta, *, status=1):
        self._id = _ObjectID(index); self._pdg = pdg
        self._momentum = _Vector(theta); self._status = status
    def getObjectID(self): return self._id
    def getPDG(self): return self._pdg
    def getMomentum(self): return self._momentum
    def getGeneratorStatus(self): return self._status
    def getEnergy(self): return 1.0


class _Event:
    def __init__(self, truth, pfos): self._data = {"MCParticles": truth, "PandoraPFOs": pfos}
    def get(self, name):
        if name in self._data: return self._data[name]
        if name in {"CalohitMCTruthLink", "SiTracksMCTruthLink"}: return []
        raise KeyError(name)


def _p4(theta, phi):
    vector = _Vector(theta, phi); result = ROOT.TLorentzVector()
    result.SetXYZM(vector.x, vector.y, vector.z, 0.0); return result


class GeometricAssociationContractTest(unittest.TestCase):
    def test_angle_wrapping(self):
        self.assertAlmostEqual(dRAngle(_p4(math.pi / 2, math.pi - .01), _p4(math.pi / 2, -math.pi + .01)), .02)

    def test_threshold_is_strict(self):
        truth = _Particle(1, 22, math.pi / 2)
        pfo = _Particle(2, 22, math.pi / 2 + .1)
        exact = dRAngle(_p4(math.pi / 2, 0), _p4(math.pi / 2 + .1, 0))
        rows = get_reco_mc_links_by_dR(_Event([truth], [pfo]), {}, {}, max_dR=exact, dedup_mode="reco")
        self.assertEqual(set(rows["gen"]), {1, -999})
        self.assertEqual(set(rows["reco"]), {-999, 2})

    def test_reco_dedup_keeps_nearest_and_loser_is_unmatched(self):
        truth = [_Particle(1, 22, math.pi / 2), _Particle(2, 22, math.pi / 2 + .02)]
        pfo = _Particle(3, 22, math.pi / 2 + .018)
        rows = get_reco_mc_links_by_dR(_Event(truth, [pfo]), {}, {}, dedup_mode="reco")
        matched = rows[rows["reco"] == 3]
        self.assertEqual(matched.iloc[0]["gen"], 2)
        self.assertTrue(((rows["gen"] == 1) & (rows["reco"] == -999)).any())

    def test_unmatched_truth_and_fake_pfo(self):
        rows = get_reco_mc_links_by_dR(
            _Event([_Particle(1, 22, math.pi / 2)], [_Particle(2, 22, 0.0)]),
            {}, {}, max_dR=.1, dedup_mode="reco",
        )
        self.assertTrue(((rows["gen"] == 1) & (rows["reco"] == -999)).any())
        self.assertTrue(((rows["gen"] == -999) & (rows["reco"] == 2)).any())


if __name__ == "__main__":
    unittest.main()
