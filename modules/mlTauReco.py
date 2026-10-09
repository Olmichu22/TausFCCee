"""ParTauDETR taus, as a drop-in for tauReco.findAllTaus.

The model runs beforehand, outside the key4hep stack (MLTauReco/, see its
README): MLTauReco/process_edm4hep.py writes one `<input name>_mltaus.parquet`
per EDM4hep file, one row per jet. This module reads those predictions back
inside the event loop and returns, per event, the {index: RecoParticle} dict
findAllTaus builds from its pion cones:

    from modules import mlTauReco
    mltaus = mlTauReco.MLTauReader("/path/to/predictions")
    for filename in filenames:
        for local_event, event in enumerate(root_io.Reader([filename]).get("events")):
            recoTau_raw = mltaus.findAllTaus(filename, local_event)

How the taus compare with tauReco.buildTauFromPion:
  - getMomentum(): the visible tau p4, the sum of the predicted daughters.
  - getID(): in buildTauFromPion's photon-counting convention, with each
    predicted pi0 counted as two photons: 2 * n_pi0 for 1 prong,
    10 + 2 * n_pi0 for 3 prongs, -1 otherwise. The ceil(ID / 2) of
    TTreesTausLong (RecoTauDM) and tauReco.MatchRecoGenTau then recover the
    pi0 count exactly. `tau.n_pi0` holds it directly.
  - getCharge(), getPDG(): as in findAllTaus (15 for negative charge, else -15).
  - getDaughters(): {i: MLTauDaughter}, the predicted pi+- / pi0 (PDG +-211 /
    111), not PFOs. Hence no photons among them, and the extra corrections
    (tauReco.extraTauRecoCorrection) do not apply. The jet's PFOs are in
    `tau.cand_pfo_idx`, as indices into the event's PandoraPFOs.
  - getMaxCone(): the largest angle between a daughter and the tau.
  - extra attributes: tau_score (p(tau)), n_pi0, jet_idx, jet_p4, cand_pfo_idx.

Jets whose tau score is below `tau_score_cut` are not returned. The default
cut is the 90%-efficiency working point of the model, measured on hadronic
taus in Z -> tautau against jets in Z -> qq; the model bundle's bundle.json
lists the others. Jets containing a reconstructed electron or muon are not
returned either, unless veto_lepton_jets=False: the model never saw such jets
in training, and electronReco / muonReco reconstruct the leptonic taus.
"""

import json
import os

import awkward as ak
import ROOT

from modules.ParticleObjects import RecoParticle

# The bundled model in MLTauReco/model, which MLTauReco/process_edm4hep.py
# also uses by default (DEFAULT_MODEL_DIR there).
DEFAULT_MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "MLTauReco", "model")
DEFAULT_WORKING_POINT = "0.90"


class MLTauDaughter:
    """A predicted tau daughter: charged pion (meson class 0) or pi0 (class 1)."""

    def __init__(self, px, py, pz, energy, charge, meson_class, objectness):
        self.p4 = ROOT.TLorentzVector(px, py, pz, energy)
        self.charge = int(charge)
        self.meson_class = int(meson_class)
        self.objectness = float(objectness)
        self.pdg = 111 if meson_class == 1 else (211 if charge >= 0 else -211)

    def getMomentum(self):
        return self.p4

    def getMass(self):
        return self.p4.M()

    def getPDG(self):
        return self.pdg

    def getCharge(self):
        return self.charge


def photon_equivalent_id(n_charged, n_pi0):
    """buildTauFromPion's ID (assignTauID), with every pi0 as two photons."""
    if n_charged == 1:
        return 2 * n_pi0
    if n_charged == 3:
        return 10 + 2 * n_pi0
    return -1


def working_points(model_dir=None):
    """{signal efficiency: {"tau_score_threshold", "qq_jet_misid"}} of the model."""
    with open(os.path.join(model_dir or DEFAULT_MODEL_DIR, "bundle.json")) as f:
        return json.load(f)["tau_id_working_points"]


class MLTauReader:
    """
    Per-event taus from the `<input name>_mltaus.parquet` files in
    `prediction_dir`, loaded lazily, one input file at a time.
    """

    def __init__(self, prediction_dir, tau_score_cut=None, veto_lepton_jets=True, model_dir=None):
        self.prediction_dir = prediction_dir
        self.veto_lepton_jets = veto_lepton_jets
        if tau_score_cut is None:
            tau_score_cut = working_points(model_dir)[DEFAULT_WORKING_POINT]["tau_score_threshold"]
        self.tau_score_cut = float(tau_score_cut)
        self._file = None
        self._rows = []
        self._by_event = {}

    def _load(self, root_filename):
        name = os.path.splitext(os.path.basename(root_filename))[0]
        if name == self._file:
            return
        path = os.path.join(self.prediction_dir, f"{name}_mltaus.parquet")
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"No ParTauDETR predictions for {root_filename}: expected {path} "
                "(run MLTauReco/process_edm4hep.py on the file first)"
            )
        rows = ak.from_parquet(path)
        self._rows = ak.to_list(rows)
        self._by_event = {}
        for i, ev in enumerate(ak.to_numpy(rows.event_idx)):
            self._by_event.setdefault(int(ev), []).append(i)
        self._file = name

    def jets(self, root_filename, local_event):
        """All jets of the event as dicts, every column, no selection."""
        self._load(root_filename)
        return [self._rows[i] for i in self._by_event.get(int(local_event), [])]

    def findAllTaus(self, root_filename, local_event):
        taus = {}
        for row in self.jets(root_filename, local_event):
            if row["tau_score"] < self.tau_score_cut:
                continue
            if self.veto_lepton_jets and row["n_reco_leptons"] > 0:
                continue
            if not row["daughter_px"]:
                continue
            const = {
                i: MLTauDaughter(*vals)
                for i, vals in enumerate(
                    zip(
                        row["daughter_px"], row["daughter_py"], row["daughter_pz"],
                        row["daughter_energy"], row["daughter_charge"],
                        row["daughter_meson_class"], row["daughter_objectness"],
                    )
                )
            }
            p4 = ROOT.TLorentzVector(row["tau_px"], row["tau_py"], row["tau_pz"], row["tau_energy"])
            max_cone = max(d.getMomentum().Angle(p4.Vect()) for d in const.values())
            tau = RecoParticle(
                p4, photon_equivalent_id(row["tau_n_charged"], row["tau_n_neutral"]),
                row["tau_charge"], max_cone, len(const), const,
            )
            tau.setPDG(15 if row["tau_charge"] < 0 else -15)
            tau.tau_score = row["tau_score"]
            tau.n_pi0 = row["tau_n_neutral"]
            tau.jet_idx = row["jet_idx"]
            tau.jet_p4 = ROOT.TLorentzVector()
            tau.jet_p4.SetPtEtaPhiE(row["jet_pt"], row["jet_eta"], row["jet_phi"], row["jet_energy"])
            tau.cand_pfo_idx = list(row["cand_pfo_idx"])
            taus[len(taus)] = tau
        return taus
