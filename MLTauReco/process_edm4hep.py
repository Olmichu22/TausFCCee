"""Reconstruct taus with ParTauDETR directly from CLD EDM4hep files.

The ML counterpart of `modules/tauReco.findAllTaus`. Instead of
seeding a cone on every charged pion, every event is clustered into jets
(ee_genkt, R = 0.4, all PandoraPFOs) and each jet goes through the model, which
predicts per jet:
  - p(tau), the tau-ID score (cut on it to select taus);
  - the visible tau daughters: 4-momentum, charge, and class (charged pion /
    neutral pion);
  - from those, the visible tau 4-momentum, charge, and decay mode.

The jet inputs are built exactly as for the training ntuples: the reco side of
ml-tau-data's ntupelizer (external/ml-tau-model/ml-tau-data, which
run_mltau.sh puts on the path) clusters the jets and computes the
jet-axis-signed impact parameters.
The model's own dataloader then turns them into features, applying the same
input scaler. Nothing reads the generator record, so the inputs do not need
MC truth and jets are not required to match a gen jet. In the training ntuples
every jet had a matching gen jet, and jets containing a reconstructed e or mu
were removed. Here they are kept and flagged (`n_reco_leptons`), and the
tau-ID score decides.

One parquet file per input file, one row per jet. `event_idx` is the entry
number in the input file, i.e. the position podio's root_io.Reader yields the
event at, so the rows join back onto an event loop over the same file (see
modules/mlTauReco.py).

Runs on CPU, in the ml-tau container:
  MLTauReco/run_mltau.sh python3 MLTauReco/process_edm4hep.py \
      <EDM4hep .root files> --out-dir <output dir>
"""

import argparse
import glob
import json
import os
import time

import awkward as ak
import numpy as np
import torch
import uproot
import vector
from omegaconf import open_dict

from ntupelizer.tools import clustering as cl
from ntupelizer.tools import general as g
from ntupelizer.tools import lifetime as lt
from ntupelizer.tools import particle_filters as pfl
from ntupelizer.tools.ntupelizing import EDM4HEPNtupelizer

from mltau.models.ParTauDETR_module import ParTauDETRModule
from mltau.tools.evaluation.decode_ParTauDETR import _predicted_components, tau_scores
from mltau.tools.io.ParTauDETR_dataloader import ParticleTransformerDETRDataset

# The bundled U1-log model: checkpoint, input scaler and bundle.json.
DEFAULT_MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "model")

PFO = "PandoraPFOs"
TRACKS = "SiTracks_Refitted"
VERTICES = "PrimaryVertices"
TRACK_STATES = f"_{TRACKS}_trackStates"
BRANCHES = [
    "EventHeader.eventNumber",
    *(
        f"{PFO}.{f}"
        for f in (
            "PDG", "energy", "mass", "charge",
            "momentum.x", "momentum.y", "momentum.z",
            "tracks_begin", "tracks_end",
        )
    ),
    f"_{PFO}_tracks.index",
    *(
        f"{TRACK_STATES}.{f}"
        for f in (
            "D0", "phi", "omega", "Z0", "tanLambda",
            "referencePoint.x", "referencePoint.y", "referencePoint.z",
            "covMatrix.values[21]",
        )
    ),
    *(f"{VERTICES}.position.{c}" for c in "xyz"),
]
# What the signed-IP computation needs; the model only reads signed dz/dxy and
# their errors (see ParTauDETR_dataloader._NEEDED_COLUMNS).
LIFETIME_VARS = [
    "dz", "dxy", "dz_error", "dxy_error",
    "pca_x", "pca_y", "pca_z", "vertex_x", "vertex_y", "vertex_z",
]
SIGNED_LIFETIME_VARS = ["dz", "dxy"]
# Unbound on purpose: neither helper reads the instance.
_constituent_property = EDM4HEPNtupelizer.get_jet_constituent_property
_candid_from_pdg = EDM4HEPNtupelizer.get_candid_from_pdg

PION_MASS = {0: 0.13957, 1: 0.13498}  # meson class -> decoded daughter mass


def load_model(model_dir: str):
    """
    The bundled checkpoint, with its scaler resolved inside `model_dir`.

    The checkpoint stores the training config; pointing output_dir and the
    scaler path at `model_dir` makes the bundle relocatable. The dataset refuses a scaler fitted on
    differently defined features, so a mismatched bundle fails loudly here.
    """
    model = ParTauDETRModule.load_from_checkpoint(
        os.path.join(model_dir, "model.ckpt"), map_location="cpu"
    )
    model.eval()
    cfg = model.cfg
    with open_dict(cfg):
        cfg.output_dir = os.path.abspath(model_dir)
        cfg.training.input_scaling.scaler_path = os.path.join(
            cfg.output_dir, "scaler", "cand_feature_scaler.npz"
        )
    ds = ParticleTransformerDETRDataset.for_arrays(cfg)
    with open(os.path.join(model_dir, "bundle.json")) as f:
        info = json.load(f)
    return model, ds, info


def read_events(path: str, entry_start=None, entry_stop=None) -> ak.Array:
    with uproot.open(path) as f:
        return f["events"].arrays(BRANCHES, entry_start=entry_start, entry_stop=entry_stop)


def _dummy_gen_columns(n_jets: int) -> dict:
    """
    Generator-level columns with background-jet values (no tau, no daughters).

    The dataloader builds training targets next to the inputs and so needs
    these columns. They only feed the targets and the is_tau label; the model
    inputs and its outputs do not depend on them.
    """
    zero_p4 = ak.zip({k: np.zeros(n_jets) for k in ("pt", "eta", "phi", "energy")})
    empty = ak.Array([[] for _ in range(n_jets)])
    empty_p4 = ak.zip(
        {k: ak.values_astype(empty, np.float64) for k in ("pt", "eta", "phi", "energy")}
    )
    return {
        "gen_jet_p4": zero_p4,
        "gen_jet_tau_p4": zero_p4,
        "gen_jet_tau_vis_daughter_p4s": empty_p4,
        "gen_jet_tau_vis_daughter_pdgs": ak.values_astype(empty, np.int64),
        "gen_jet_tau_vis_daughter_charges": ak.values_astype(empty, np.float64),
        "gen_jet_tau_decaymode": np.full(n_jets, -1, dtype=np.int64),
        "gen_jet_tau_charge": np.full(n_jets, -999.0),
    }


def build_jets(events: ak.Array, entry_offset: int = 0) -> ak.Array:
    """
    Per-jet model inputs, as the training ntupelizer (ml-tau-data
    EDM4HEPNtupelizer.ntupelize) builds them, minus everything generator-level.
    """
    reco_particles, reco_p4 = pfl.RecoParticleFilter(arrays=events, p_type=PFO).results
    reco_jets, constituent_idx = cl.RecoJetClusterer(
        particles=reco_particles, particles_p4=reco_p4
    ).results
    n_per_jet = ak.num(constituent_idx, axis=-1)

    lifetime_all = lt.find_all_track_pcas(
        events=events,
        reco_particle_collection=PFO,
        track_collection=TRACKS,
        vertex_collection=VERTICES,
        valid_particle_mask=reco_particles["PDG"] != 0,
    )
    lifetime = lt.assign_lifetime_vars_to_jets(
        all_particle_lifetime_info=lifetime_all,
        reco_jet_constituent_indices=constituent_idx,
        reco_jets=reco_jets,
        lifetime_vars=LIFETIME_VARS,
        signed_lifetime_vars=SIGNED_LIFETIME_VARS,
    )

    def per_cand(prop):
        return _constituent_property(None, prop, constituent_idx, n_per_jet)

    event_idx = ak.local_index(reco_jets.px, axis=0) + entry_offset
    jets = {
        "event_idx": ak.broadcast_arrays(event_idx, reco_jets.px)[0],
        "event_number": ak.broadcast_arrays(
            ak.firsts(events["EventHeader.eventNumber"]), reco_jets.px
        )[0],
        "jet_idx": ak.local_index(reco_jets.px, axis=1),
        "reco_jet_p4": g.reinitialize_p4(reco_jets),
        "reco_cand_p4s": g.reinitialize_p4(per_cand(reco_p4)),
        "reco_cand_pdgs": per_cand(_candid_from_pdg(None, reco_particles)),
        "reco_cand_charges": per_cand(reco_particles.charge),
        # Index into the event's PandoraPFOs, to link back to the PFO objects.
        "reco_cand_pfo_idx": constituent_idx,
        "reco_cand_signed_dz": lifetime["signed_dz"],
        "reco_cand_dz_error": lifetime["dz_error"],
        "reco_cand_signed_dxy": lifetime["signed_dxy"],
        "reco_cand_dxy_error": lifetime["dxy_error"],
    }
    # vector's internal {rho, phi, eta, t} -> the ntuples' {pt, eta, phi, energy}
    for key in ("reco_jet_p4", "reco_cand_p4s"):
        jets[key] = ak.zip(
            {k: getattr(jets[key], k) for k in ("pt", "eta", "phi", "energy")}
        )
    data = ak.Array({k: ak.flatten(v, axis=1) for k, v in jets.items()})
    for k, v in _dummy_gen_columns(len(data)).items():
        data[k] = v
    return data


def predict(model, ds, jets: ak.Array, threshold: float, mass_from_class: bool) -> dict:
    """
    Model outputs per jet: p(tau), and the predicted daughters (queries whose
    objectness passes `threshold`), sorted by pT.

    Daughters are written for every jet, also those the tagger rejects:
    objectness is trained on tau jets only, so on non-tau jets the daughter set
    means nothing, and the tau-ID score is what tells the two apart.
    """
    with torch.no_grad():
        batch = ds.build_tensors(jets)
        outputs = model.forward(batch)[0]
        reco_jet_p4s = batch[6]
        score = tau_scores(outputs)
        obj, p4, charge, meson_class = _predicted_components(
            outputs, reco_jet_p4s, mass_from_class=mass_from_class
        )
    keep = obj >= threshold
    px, py, pz, en = (ak.to_numpy(getattr(p4, c)) for c in ("px", "py", "pz", "energy"))
    keep_np = keep.numpy()
    counts = keep_np.sum(axis=1)

    def jag(x):
        return ak.unflatten(x[keep_np], counts)

    dau = {
        "px": jag(px), "py": jag(py), "pz": jag(pz), "energy": jag(en),
        "charge": jag(charge.numpy().astype(np.int8)),
        "meson_class": jag(meson_class.numpy().astype(np.int8)),
        "objectness": jag(obj.numpy().astype(np.float32)),
    }
    dau_pt = np.hypot(dau["px"], dau["py"])
    order = ak.argsort(dau_pt, axis=1, ascending=False)
    dau = {k: v[order] for k, v in dau.items()}

    tau = vector.awk(
        ak.zip({c: ak.sum(dau[c], axis=1) for c in ("px", "py", "pz", "energy")})
    )
    n_charged = ak.to_numpy(ak.sum(dau["meson_class"] == 0, axis=1))
    n_neutral = ak.to_numpy(ak.sum(dau["meson_class"] == 1, axis=1))
    # Same convention as training/evaluation: 5 * (n_charged - 1) + n_neutral,
    # i.e. 0 = h, 1 = h pi0, 2 = h 2pi0, 10 = 3h, 11 = 3h pi0; -1 = no charged
    # daughter. Even n_charged gives a valid index but no physical tau decay.
    dm = np.where(n_charged > 0, 5 * (n_charged - 1) + n_neutral, -1)
    has_dau = counts > 0
    return {
        "tau_score": score.numpy().astype(np.float32),
        "tau_px": ak.to_numpy(tau.px), "tau_py": ak.to_numpy(tau.py),
        "tau_pz": ak.to_numpy(tau.pz), "tau_energy": ak.to_numpy(tau.energy),
        "tau_pt": ak.to_numpy(tau.pt),
        "tau_eta": np.where(has_dau, ak.to_numpy(tau.eta), 0.0),
        "tau_phi": ak.to_numpy(tau.phi),
        "tau_mass": np.where(has_dau, ak.to_numpy(tau.mass), 0.0),
        "tau_charge": ak.to_numpy(ak.sum(dau["charge"], axis=1)).astype(np.int8),
        "tau_decaymode": dm.astype(np.int16),
        "tau_n_charged": n_charged.astype(np.int8),
        "tau_n_neutral": n_neutral.astype(np.int8),
        **{f"daughter_{k}": v for k, v in dau.items()},
    }


def process_file(path, out_path, model, ds, info, events_per_chunk, threshold):
    t0 = time.time()
    with uproot.open(path) as f:
        n_events = f["events"].num_entries
    parts = []
    for start in range(0, n_events, events_per_chunk):
        stop = min(start + events_per_chunk, n_events)
        jets = build_jets(read_events(path, start, stop), entry_offset=start)
        if len(jets) == 0:
            continue
        pred = predict(model, ds, jets, threshold, info["mass_from_class"])
        out = {
            "event_idx": jets.event_idx,
            "event_number": jets.event_number,
            "jet_idx": jets.jet_idx,
            **{f"jet_{k}": jets.reco_jet_p4[k] for k in ("pt", "eta", "phi", "energy")},
            "jet_n_cands": ak.num(jets.reco_cand_pdgs),
            "n_reco_leptons": ak.sum(
                (jets.reco_cand_pdgs == 11) | (jets.reco_cand_pdgs == 13), axis=1
            ),
            "cand_pfo_idx": jets.reco_cand_pfo_idx,
            **pred,
        }
        parts.append(ak.Array(out))
    data = ak.concatenate(parts) if parts else ak.Array([])
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    ak.to_parquet(data, out_path, compression="zstd")
    print(
        f"{os.path.basename(path)}: {n_events} events, {len(data)} jets "
        f"-> {out_path} ({time.time() - t0:.0f} s)",
        flush=True,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("inputs", nargs="+", help="EDM4hep ROOT files or glob patterns")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--model-dir", default=DEFAULT_MODEL_DIR)
    ap.add_argument("--events-per-chunk", type=int, default=1000)
    ap.add_argument(
        "--threshold", type=float, default=None,
        help="daughter objectness threshold; default: the checkpoint's calibrated value",
    )
    ap.add_argument("--threads", type=int, default=0, help="torch CPU threads; 0 = all")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    torch.set_num_threads(args.threads or max(1, os.cpu_count() or 1))
    model, ds, info = load_model(args.model_dir)
    threshold = (
        args.threshold
        if args.threshold is not None
        else float(getattr(model, "score_threshold_calibrated", 0.5))
    )
    paths = sorted({p for pat in args.inputs for p in (glob.glob(pat) or [pat])})
    print(f"model {info['name']}, {len(paths)} input file(s), objectness threshold {threshold}")
    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        out_path = os.path.join(args.out_dir, f"{name}_mltaus.parquet")
        if os.path.exists(out_path) and not args.overwrite:
            print(f"{out_path} exists, skipped")
            continue
        process_file(path, out_path, model, ds, info, args.events_per_chunk, threshold)


if __name__ == "__main__":
    main()
