import sys
import os
import math
import time
import subprocess
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import ROOT
from ROOT import TFile, TTree
import numpy as np
from podio import root_io
import edm4hep
import pprint
import yaml
import logging

from modules import tauReco, electronReco, muonReco, myutils


_NEUTRINO_PDGS = {12, 14, 16}


# ── File splitting helpers ─────────────────────────────────────────────────────

def split_filenames(filenames, n_workers):
    k, rem = divmod(len(filenames), n_workers)
    chunks, start = [], 0
    for i in range(n_workers):
        end = start + k + (1 if i < rem else 0)
        if start < end:
            chunks.append(filenames[start:end])
        start = end
    return chunks


def split_mlpf(mlpf_results, file_chunks):
    """Split mlpf_results into per-worker sub-dicts with renormalized keys."""
    mlpf_chunks = []
    file_offset = 0
    for chunk in file_chunks:
        lo = file_offset * 1000
        hi = (file_offset + len(chunk)) * 1000
        sub = {k - lo: v for k, v in mlpf_results.items() if lo <= k < hi}
        mlpf_chunks.append(sub)
        file_offset += len(chunk)
    return mlpf_chunks


# ── Gen-level decay-tree helper ────────────────────────────────────────────────

def get_final_state_constituents(gen_tau_consts):
    """Recursively traverse the decay tree from a gen tau's direct daughters.

    Takes the const dict from GenParticle.getDaughters() (raw MCParticle objects)
    and returns the final-state (generatorStatus==1, non-neutrino) leaves,
    each tagged with the index of their pi0 ancestor within this tau (-1 if none).

    Returns: list of (MCParticle, pi0_ancestor_idx)
    """
    results = []
    pi0_counter = [0]

    def _recurse(particle, pi0_idx):
        pdg = abs(particle.getPDG())
        if pdg in _NEUTRINO_PDGS:
            return
        daughters = list(particle.getDaughters())
        status = particle.getGeneratorStatus()
        if not daughters or status == 1:
            results.append((particle, pi0_idx))
            return
        for d in daughters:
            if abs(d.getPDG()) == 111:
                idx = pi0_counter[0]
                pi0_counter[0] += 1
                _recurse(d, idx)
            else:
                _recurse(d, pi0_idx)

    for key in sorted(gen_tau_consts.keys()):
        part = gen_tau_consts[key]
        pdg = abs(part.getPDG())
        if pdg in _NEUTRINO_PDGS:
            continue
        if pdg == 111:
            idx = pi0_counter[0]
            pi0_counter[0] += 1
            _recurse(part, idx)
        else:
            _recurse(part, -1)
    return results


# ── RecoMCTruthLink photon matching helper ─────────────────────────────────────

def build_truth_links(event, filter_gen_status=True, max_gen_pdg=10000,
                      weight_mode="decoded", dedup_mode="reco"):
    """Build bidirectional gen<->reco index maps from RecoMCTruthLink.

    Returns dict with:
      'gen_to_reco': {mc_obj_index: pfo_obj_index}
      'reco_to_gen': {pfo_obj_index: mc_obj_index}
    Uses getObjectID().index for both sides.
    """
    try:
        links = list(event.get("RecoMCTruthLink"))
    except Exception:
        return {"gen_to_reco": {}, "reco_to_gen": {}}

    if not links:
        return {"gen_to_reco": {}, "reco_to_gen": {}}

    rows = []
    for link in links:
        try:
            if hasattr(link, "getSim") and hasattr(link, "getRec"):
                sim_obj = link.getSim()
                rec_obj = link.getRec()
            else:
                sim_obj = link.getTo()
                rec_obj = link.getFrom()

            gen_status = sim_obj.getGeneratorStatus()
            gen_pdg = sim_obj.getPDG()
            if filter_gen_status and gen_status != 1:
                continue
            if abs(int(gen_pdg)) in _NEUTRINO_PDGS:
                continue
            if max_gen_pdg is not None and abs(int(gen_pdg)) > max_gen_pdg:
                continue

            gen_idx = sim_obj.getObjectID().index
            reco_idx = rec_obj.getObjectID().index
            weight = getattr(link, "getWeight", lambda: 0.0)()
            rows.append((gen_idx, reco_idx, float(weight)))
        except Exception:
            continue

    if not rows:
        return {"gen_to_reco": {}, "reco_to_gen": {}}

    if weight_mode == "decoded":
        def _eff_w(w):
            enc = int(w)
            tw = (enc % 10000) / 1000.0
            cw = (enc // 10000) / 1000.0
            return tw if tw > 0.0 else cw
        rows = [(g, r, _eff_w(w)) for g, r, w in rows]

    if dedup_mode == "reco":
        best = {}
        for g, r, w in rows:
            if r not in best or w > best[r][1]:
                best[r] = (g, w)
        pairs = [(g, r) for r, (g, _) in best.items()]
    else:
        best = {}
        for g, r, w in rows:
            if g not in best or w > best[g][1]:
                best[g] = (r, w)
        pairs = [(g, r) for g, (r, _) in best.items()]

    return {
        "gen_to_reco": {g: r for g, r in pairs},
        "reco_to_gen": {r: g for g, r in pairs},
    }


# ── Worker ─────────────────────────────────────────────────────────────────────

def process_chunk(filenames_chunk, mlpf_chunk, global_event_offset,
                  config_bundle, worker_id):
    """Process a chunk of ROOT files and write results to a temporary TFile.

    Returns the path of the written temp file.
    """
    # Worker-local logging (avoid concurrent writes to the parent's log files)
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()
    log_file = os.path.join(config_bundle["outputpath"], f"worker_{worker_id}.log")
    logging.basicConfig(
        filename=log_file, level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s", force=True,
    )
    logger = logging.getLogger(f"worker_{worker_id}")
    logger.info("Worker %d started, processing %d files", worker_id, len(filenames_chunk))

    # Extract config
    cuts              = config_bundle["cuts"]
    dRMax             = cuts["dRMax"]
    minPTauPhoton     = cuts["TauPhotonPCut"]
    minPTauPion       = cuts["TauPionPCut"]
    PNeutron          = cuts["NeutronCut"]
    generalPCut       = cuts["generalPCut"]
    dRMatch           = cuts["dRMatch"]
    selectDecay       = config_bundle["selectDecay"]
    pfobjects         = config_bundle["pfobjects"]
    genparts          = config_bundle["genparts"]
    gatr_results_path = config_bundle["gatr_results_path"]
    test_pfo          = config_bundle["test_pfo"]
    outputpath        = config_bundle["outputpath"]
    cut_string        = config_bundle["cut_string"]
    sample            = config_bundle["sample"]
    weight_mode       = config_bundle.get("weight_mode", "decoded")
    dedup_mode        = config_bundle.get("dedup_mode", "reco")
    filter_gen_status = config_bundle.get("filter_gen_status", True)
    max_gen_pdg       = config_bundle.get("max_gen_pdg", 10000)

    # Temporary output file — created inside the worker after fork
    tmp_path = os.path.join(outputpath, f"tmp_chunk_{worker_id}.root")
    outfile_tmp = TFile(tmp_path, "RECREATE")
    tree = TTree("Tau_tree", f"Tree {cut_string}_Sample_{sample}")

    # ── Branch declarations ──────────────────────────────────────────────────
    numGenTaus    = np.array([0], dtype=int)
    numRecoTaus   = np.array([0], dtype=int)
    numGenPhotons = np.array([0], dtype=int)
    numRecoPhotons= np.array([0], dtype=int)
    beamE         = np.array([0.0], dtype=np.float64)

    # Gen tau
    GenEventId       = ROOT.std.vector("int")()
    GenTauPt         = ROOT.std.vector("float")()
    GenVisTauPt      = ROOT.std.vector("float")()
    GenTauP          = ROOT.std.vector("float")()
    GenVisTauP       = ROOT.std.vector("float")()
    GenTauType       = ROOT.std.vector("int")()
    GenVisTauMass    = ROOT.std.vector("float")()
    GenTauQ          = ROOT.std.vector("float")()
    GenTauEta        = ROOT.std.vector("float")()
    GenTauTheta      = ROOT.std.vector("float")()
    GenTauPhi        = ROOT.std.vector("float")()
    GenTauMass       = ROOT.std.vector("float")()
    GenVisTauEta     = ROOT.std.vector("float")()
    GenVisTauTheta   = ROOT.std.vector("float")()
    GenVisTauPhi     = ROOT.std.vector("float")()
    GenTauDR         = ROOT.std.vector("float")()
    GenTauNConsts    = ROOT.std.vector("int")()
    GenTauNConstKey  = ROOT.std.vector("int")()
    GenTauConstKey   = ROOT.std.vector("int")()
    GenMatchedKey    = ROOT.std.vector("int")()
    GenConstP        = ROOT.std.vector("float")()
    GenConstTheta    = ROOT.std.vector("float")()
    GenConstEta      = ROOT.std.vector("float")()
    GenConstPDG      = ROOT.std.vector("int")()
    GenConstPi0Key   = ROOT.std.vector("int")()
    GenConstPhi      = ROOT.std.vector("float")()

    # Reco tau
    RecoTauPt        = ROOT.std.vector("float")()
    RecoTauP         = ROOT.std.vector("float")()
    RecoTauMass      = ROOT.std.vector("float")()
    RecoTauType      = ROOT.std.vector("int")()
    RecoTauDM        = ROOT.std.vector("int")()
    RecoTauQ         = ROOT.std.vector("float")()
    RecoTauEta       = ROOT.std.vector("float")()
    RecoTauTheta     = ROOT.std.vector("float")()
    RecoTauPhi       = ROOT.std.vector("float")()
    RecoTauDR        = ROOT.std.vector("float")()
    RecoTauNConsts   = ROOT.std.vector("int")()
    RecoTauNConstKey = ROOT.std.vector("int")()
    RecoTauConstKey  = ROOT.std.vector("int")()
    RecoMatchedKey   = ROOT.std.vector("int")()
    RecoConstP       = ROOT.std.vector("float")()
    RecoConstTheta   = ROOT.std.vector("float")()
    RecoConstEta     = ROOT.std.vector("float")()
    RecoConstPDG     = ROOT.std.vector("int")()
    RecoConstPhi     = ROOT.std.vector("float")()
    # Gen photons (event-level)
    GenPhotonP       = ROOT.std.vector("float")()
    GenPhotonPt      = ROOT.std.vector("float")()
    GenPhotonEta     = ROOT.std.vector("float")()
    GenPhotonTheta   = ROOT.std.vector("float")()
    GenPhotonPhi     = ROOT.std.vector("float")()
    GenPhotonMCIdx   = ROOT.std.vector("int")()
    GenPhotonTauKey  = ROOT.std.vector("int")()

    # Reco photons (event-level, from PandoraPFOs)
    RecoPhotonP            = ROOT.std.vector("float")()
    RecoPhotonPt           = ROOT.std.vector("float")()
    RecoPhotonEta          = ROOT.std.vector("float")()
    RecoPhotonTheta        = ROOT.std.vector("float")()
    RecoPhotonPhi          = ROOT.std.vector("float")()
    RecoPhotonPFOIdx       = ROOT.std.vector("int")()
    RecoPhotonTauKey       = ROOT.std.vector("int")()
    RecoPhotonGenMatchIdx  = ROOT.std.vector("int")()

    tree.Branch("numGenTaus",     numGenTaus,     "numGenTaus/I")
    tree.Branch("numRecoTaus",    numRecoTaus,    "numRecoTaus/I")
    tree.Branch("numGenPhotons",  numGenPhotons,  "numGenPhotons/I")
    tree.Branch("numRecoPhotons", numRecoPhotons, "numRecoPhotons/I")
    tree.Branch("beamE",          beamE,          "beamE/D")
    # Gen tau
    tree.Branch("GenEventId",      GenEventId)
    tree.Branch("GenTauPt",        GenTauPt)
    tree.Branch("GenVisTauPt",     GenVisTauPt)
    tree.Branch("GenTauP",         GenTauP)
    tree.Branch("GenVisTauP",      GenVisTauP)
    tree.Branch("GenTauType",      GenTauType)
    tree.Branch("GenVisTauMass",   GenVisTauMass)
    tree.Branch("GenTauQ",         GenTauQ)
    tree.Branch("GenTauEta",       GenTauEta)
    tree.Branch("GenTauTheta",     GenTauTheta)
    tree.Branch("GenTauPhi",       GenTauPhi)
    tree.Branch("GenTauMass",      GenTauMass)
    tree.Branch("GenVisTauEta",    GenVisTauEta)
    tree.Branch("GenVisTauTheta",  GenVisTauTheta)
    tree.Branch("GenVisTauPhi",    GenVisTauPhi)
    tree.Branch("GenTauDR",        GenTauDR)
    tree.Branch("GenTauNConsts",   GenTauNConsts)
    tree.Branch("GenTauNConstKey", GenTauNConstKey)
    tree.Branch("GenTauConstKey",  GenTauConstKey)
    tree.Branch("GenMatchedKey",   GenMatchedKey)
    tree.Branch("GenConstP",       GenConstP)
    tree.Branch("GenConstTheta",   GenConstTheta)
    tree.Branch("GenConstEta",     GenConstEta)
    tree.Branch("GenConstPDG",     GenConstPDG)
    tree.Branch("GenConstPi0Key",  GenConstPi0Key)
    tree.Branch("GenConstPhi",     GenConstPhi)
    # Reco tau
    tree.Branch("RecoTauPt",        RecoTauPt)
    tree.Branch("RecoTauP",         RecoTauP)
    tree.Branch("RecoTauMass",      RecoTauMass)
    tree.Branch("RecoTauType",      RecoTauType)
    tree.Branch("RecoTauDM",        RecoTauDM)
    tree.Branch("RecoTauQ",         RecoTauQ)
    tree.Branch("RecoTauEta",       RecoTauEta)
    tree.Branch("RecoTauTheta",     RecoTauTheta)
    tree.Branch("RecoTauPhi",       RecoTauPhi)
    tree.Branch("RecoTauDR",        RecoTauDR)
    tree.Branch("RecoTauNConsts",   RecoTauNConsts)
    tree.Branch("RecoTauNConstKey", RecoTauNConstKey)
    tree.Branch("RecoTauConstKey",  RecoTauConstKey)
    tree.Branch("RecoMatchedKey",   RecoMatchedKey)
    tree.Branch("RecoConstP",       RecoConstP)
    tree.Branch("RecoConstTheta",   RecoConstTheta)
    tree.Branch("RecoConstEta",     RecoConstEta)
    tree.Branch("RecoConstPDG",     RecoConstPDG)
    tree.Branch("RecoConstPhi",     RecoConstPhi)
    # Gen photons
    tree.Branch("GenPhotonP",      GenPhotonP)
    tree.Branch("GenPhotonPt",     GenPhotonPt)
    tree.Branch("GenPhotonEta",    GenPhotonEta)
    tree.Branch("GenPhotonTheta",  GenPhotonTheta)
    tree.Branch("GenPhotonPhi",    GenPhotonPhi)
    tree.Branch("GenPhotonMCIdx",  GenPhotonMCIdx)
    tree.Branch("GenPhotonTauKey", GenPhotonTauKey)
    # Reco photons
    tree.Branch("RecoPhotonP",           RecoPhotonP)
    tree.Branch("RecoPhotonPt",          RecoPhotonPt)
    tree.Branch("RecoPhotonEta",         RecoPhotonEta)
    tree.Branch("RecoPhotonTheta",       RecoPhotonTheta)
    tree.Branch("RecoPhotonPhi",         RecoPhotonPhi)
    tree.Branch("RecoPhotonPFOIdx",      RecoPhotonPFOIdx)
    tree.Branch("RecoPhotonTauKey",      RecoPhotonTauKey)
    tree.Branch("RecoPhotonGenMatchIdx", RecoPhotonGenMatchIdx)

    all_vectors = [
        GenEventId, GenTauPt, GenVisTauPt, GenTauP, GenVisTauP, GenTauType, GenVisTauMass,
        GenTauQ, GenTauEta, GenTauTheta, GenTauPhi, GenTauMass,
        GenVisTauEta, GenVisTauTheta, GenVisTauPhi,
        GenTauDR, GenTauNConsts, GenTauNConstKey,
        GenTauConstKey, GenMatchedKey, GenConstP, GenConstTheta, GenConstEta, GenConstPDG, GenConstPhi,
        GenConstPi0Key,
        RecoTauPt, RecoTauP, RecoTauMass, RecoTauType, RecoTauDM, RecoTauQ, RecoTauEta,
        RecoTauTheta, RecoTauPhi, RecoTauDR, RecoTauNConsts, RecoTauNConstKey, RecoTauConstKey,
        RecoMatchedKey, RecoConstP, RecoConstTheta, RecoConstEta, RecoConstPDG, RecoConstPhi,
        GenPhotonP, GenPhotonPt, GenPhotonEta, GenPhotonTheta, GenPhotonPhi,
        GenPhotonMCIdx, GenPhotonTauKey,
        RecoPhotonP, RecoPhotonPt, RecoPhotonEta, RecoPhotonTheta, RecoPhotonPhi,
        RecoPhotonPFOIdx, RecoPhotonTauKey, RecoPhotonGenMatchIdx,
    ]

    # ── Event loop ──────────────────────────────────────────────────────────
    cumulative_local_eventid = 0
    for filename in filenames_chunk:
        file_reader = root_io.Reader([filename])
        file_local_eventid = -1

        for file_local_eventid, event in enumerate(file_reader.get("events")):
            local_eventid = cumulative_local_eventid + file_local_eventid
            event_id_global = global_event_offset + local_eventid
            if local_eventid % 500 == 0:
                logger.info("Worker %d: local event %d (global %d)",
                            worker_id, local_eventid, event_id_global)

            mc_particles = event.get(genparts)
            pfos = event.get(pfobjects)  # PandoraPFOs — always used for reco photons
            beamE[0] = mc_particles[0].getEnergy()

            # Reco tau reconstruction (MLPF or PFO)
            if gatr_results_path is not None and not test_pfo:
                particles = mlpf_chunk.get(local_eventid, {})
                recoTau_raw = tauReco.findAllTaus(
                    particles, dRMax, minPTauPhoton, minPTauPion,
                    PNeutron, generalPCut, charge_condition=False
                )
                recoElectrons = electronReco.findAllElectrons(particles, generalPCut)
                recoMuons = muonReco.findAllMuons(particles, generalPCut)
            else:
                recoTau_raw = tauReco.findAllTaus(
                    pfos, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut
                )
                recoElectrons = electronReco.findAllElectrons(pfos, generalPCut)
                recoMuons = muonReco.findAllMuons(pfos, generalPCut)

            genTaus = tauReco.findAllGenTaus(mc_particles)
            nGenTaus = len(genTaus)

            recoTaus = {}
            pidx = 0
            for t in range(len(recoTau_raw)):
                recoTaus[pidx] = recoTau_raw[t]; pidx += 1
            for e in range(len(recoElectrons)):
                recoTaus[pidx] = recoElectrons[e]; pidx += 1
            for m in range(len(recoMuons)):
                recoTaus[pidx] = recoMuons[m]; pidx += 1
            nRecoTaus = len(recoTaus)

            for v in all_vectors:
                v.clear()
            numGenTaus[0]  = nGenTaus
            numRecoTaus[0] = nRecoTaus

            # ── Pre-compute final-state constituents for all gen taus ──────
            # Cached to avoid double traversal (used both for branch filling
            # and for building the photon tau-membership index).
            tau_final_consts = {
                i: get_final_state_constituents(genTaus[i].getDaughters())
                for i in range(nGenTaus)
            }

            # ── Gen photon → gen tau membership index ─────────────────────
            # Maps MCParticle object-index → gen tau array index
            tau_photon_mc_idx = {}
            for i, consts in tau_final_consts.items():
                for (part, _) in consts:
                    if abs(part.getPDG()) == 22:
                        try:
                            tau_photon_mc_idx[part.getObjectID().index] = i
                        except Exception:
                            pass

            # ── Reco photon → reco tau membership index ────────────────────
            # Maps rounded momentum tuple → reco tau array index.
            # Works for both PFO and MLPF constituents; when MLPF is used
            # and reco photons come from PandoraPFOs, this will yield no
            # matches (RecoPhotonTauKey = -1 for all), which is correct.
            reco_tau_photon_set = {}
            for j in range(nRecoTaus):
                for _, const_c in recoTaus[j].getDaughters().items():
                    try:
                        if abs(const_c.getPDG()) != 22:
                            continue
                        p = const_c.getMomentum()
                        try:
                            key = (round(p.x, 5), round(p.y, 5), round(p.z, 5))
                        except AttributeError:
                            key = (round(p.X(), 5), round(p.Y(), 5), round(p.Z(), 5))
                        reco_tau_photon_set[key] = j
                    except Exception:
                        continue

            # ── Fill gen tau branches ──────────────────────────────────────
            for i in range(nGenTaus):
                genVisTauP4 = genTaus[i].getvisMomentum()
                genTauP4    = genTaus[i].getMomentum()

                GenEventId.push_back(event_id_global)
                GenTauPt.push_back(genTauP4.Pt())
                GenVisTauPt.push_back(genVisTauP4.Pt())
                GenTauP.push_back(genTauP4.P())
                GenVisTauP.push_back(genVisTauP4.P())
                GenVisTauMass.push_back(genVisTauP4.M())
                GenTauType.push_back(genTaus[i].getID())
                GenTauQ.push_back(genTaus[i].getCharge())
                GenTauEta.push_back(genTauP4.Eta())
                GenTauTheta.push_back(genTauP4.Theta())
                GenTauPhi.push_back(genTauP4.Phi())
                GenTauMass.push_back(genTauP4.M())
                GenVisTauEta.push_back(genVisTauP4.Eta())
                GenVisTauTheta.push_back(genVisTauP4.Theta())
                GenVisTauPhi.push_back(genVisTauP4.Phi())
                GenTauDR.push_back(genTaus[i].getMaxAngle())
                GenTauNConstKey.push_back(i)

                final_consts = tau_final_consts[i]
                GenTauNConsts.push_back(len(final_consts))
                for (part, pi0_idx) in final_consts:
                    constP4 = ROOT.TLorentzVector()
                    constP4.SetXYZM(
                        part.getMomentum().x, part.getMomentum().y,
                        part.getMomentum().z, part.getMass(),
                    )
                    GenTauConstKey.push_back(i)
                    GenConstPDG.push_back(part.getPDG())
                    GenConstP.push_back(constP4.P())
                    GenConstTheta.push_back(constP4.Theta())
                    GenConstEta.push_back(constP4.Eta())
                    GenConstPhi.push_back(constP4.Phi())
                    GenConstPi0Key.push_back(pi0_idx)

            # ── Fill reco tau branches ─────────────────────────────────────
            for i in range(nRecoTaus):
                recoTauP4      = recoTaus[i].getMomentum()
                recoTauId      = recoTaus[i].getID()
                recoTauNConsts = recoTaus[i].getnConst()
                recoTauConsts  = recoTaus[i].getDaughters()

                recoDM = recoTauId
                if 0 <= recoTauId < 10:
                    recoDM = math.ceil(recoTauId / 2)
                elif recoTauId >= 10:
                    recoDM = 10 + math.ceil((recoTauId - 10) / 2)

                RecoTauPt.push_back(recoTauP4.Pt())
                RecoTauP.push_back(recoTauP4.P())
                RecoTauMass.push_back(recoTauP4.M())
                RecoTauType.push_back(recoTauId)
                RecoTauDM.push_back(recoDM)
                RecoTauQ.push_back(recoTaus[i].getCharge())
                RecoTauEta.push_back(recoTauP4.Eta())
                RecoTauTheta.push_back(recoTauP4.Theta())
                RecoTauPhi.push_back(recoTauP4.Phi())
                RecoTauDR.push_back(recoTaus[i].getMaxCone())
                RecoTauNConstKey.push_back(i)
                RecoTauNConsts.push_back(recoTauNConsts)

                for c in range(recoTauNConsts):
                    const_c = recoTauConsts[c]
                    constP4 = ROOT.TLorentzVector()
                    try:
                        constP4.SetXYZM(
                            const_c.getMomentum().x, const_c.getMomentum().y,
                            const_c.getMomentum().z, const_c.getMass(),
                        )
                    except AttributeError:
                        constP4.SetXYZM(
                            const_c.getMomentum().X(), const_c.getMomentum().Y(),
                            const_c.getMomentum().Z(), const_c.getMass(),
                        )
                    RecoTauConstKey.push_back(i)
                    RecoConstPDG.push_back(const_c.getPDG())
                    RecoConstP.push_back(constP4.P())
                    RecoConstTheta.push_back(constP4.Theta())
                    RecoConstEta.push_back(constP4.Eta())
                    RecoConstPhi.push_back(constP4.Phi())
            # ── Gen–reco tau matching ──────────────────────────────────────
            nTausType = 0
            for i in range(nGenTaus):
                findMatch, nTausType = tauReco.MatchRecoGenTau(
                    genTaus[i], recoTaus, nTausType,
                    maxDRMatch=dRMatch, selectDecay=selectDecay,
                )
                GenMatchedKey.push_back(i)
                RecoMatchedKey.push_back(findMatch)

            # ── RecoMCTruthLink photon matching ────────────────────────────
            truth_links = build_truth_links(
                event,
                filter_gen_status=filter_gen_status,
                max_gen_pdg=max_gen_pdg,
                weight_mode=weight_mode,
                dedup_mode=dedup_mode,
            )
            reco_to_gen = truth_links["reco_to_gen"]

            # ── Gen photons (all event-level, generatorStatus==1) ──────────
            gen_photon_list = []
            for mc_part in mc_particles:
                try:
                    if abs(mc_part.getPDG()) != 22 or mc_part.getGeneratorStatus() != 1:
                        continue
                    mc_idx  = mc_part.getObjectID().index
                    tau_key = tau_photon_mc_idx.get(mc_idx, -1)
                    p4 = ROOT.TLorentzVector()
                    p4.SetXYZM(
                        mc_part.getMomentum().x, mc_part.getMomentum().y,
                        mc_part.getMomentum().z, mc_part.getMass(),
                    )
                    gen_photon_list.append((mc_idx, tau_key, p4))
                except Exception:
                    continue

            # mc_obj_index → position in GenPhoton* arrays (for reco→gen lookup)
            gen_mc_to_pos = {mc_idx: pos for pos, (mc_idx, _, _) in enumerate(gen_photon_list)}

            for (mc_idx, tau_key, p4) in gen_photon_list:
                GenPhotonP.push_back(p4.P())
                GenPhotonPt.push_back(p4.Pt())
                GenPhotonEta.push_back(p4.Eta())
                GenPhotonTheta.push_back(p4.Theta())
                GenPhotonPhi.push_back(p4.Phi())
                GenPhotonMCIdx.push_back(mc_idx)
                GenPhotonTauKey.push_back(tau_key)
            numGenPhotons[0] = len(gen_photon_list)

            # ── Reco photons (from PandoraPFOs — consistent with TruthLink) ─
            n_reco_photons = 0
            for pfo in pfos:
                try:
                    if abs(pfo.getPDG()) != 22:
                        continue

                    pfo_idx = pfo.getObjectID().index

                    # Tau membership via momentum signature
                    p = pfo.getMomentum()
                    try:
                        mom_key = (round(p.x, 5), round(p.y, 5), round(p.z, 5))
                    except AttributeError:
                        mom_key = (round(p.X(), 5), round(p.Y(), 5), round(p.Z(), 5))
                    tau_key = reco_tau_photon_set.get(mom_key, -1)

                    # Gen match via RecoMCTruthLink
                    gen_mc_idx = reco_to_gen.get(pfo_idx, -1)
                    gen_pos    = gen_mc_to_pos.get(gen_mc_idx, -1)

                    p4 = ROOT.TLorentzVector()
                    try:
                        p4.SetXYZM(p.x, p.y, p.z, pfo.getMass())
                    except AttributeError:
                        p4.SetXYZM(p.X(), p.Y(), p.Z(), pfo.getMass())

                    RecoPhotonP.push_back(p4.P())
                    RecoPhotonPt.push_back(p4.Pt())
                    RecoPhotonEta.push_back(p4.Eta())
                    RecoPhotonTheta.push_back(p4.Theta())
                    RecoPhotonPhi.push_back(p4.Phi())
                    RecoPhotonPFOIdx.push_back(pfo_idx)
                    RecoPhotonTauKey.push_back(tau_key)
                    RecoPhotonGenMatchIdx.push_back(gen_pos)
                    n_reco_photons += 1
                except Exception:
                    continue
            numRecoPhotons[0] = n_reco_photons

            tree.Fill()

        if file_local_eventid >= 0:
            cumulative_local_eventid += file_local_eventid + 1

    outfile_tmp.cd()
    tree.Write()
    outfile_tmp.Close()
    logger.info("Worker %d finished → %s", worker_id, tmp_path)
    return tmp_path


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    default_config = "config/default/taurecolong.yaml"
    outputbasepath = "Results/TauReco/"

    def _parser_hook(parser):
        parser.add_argument(
            "--n-workers", type=int, default=None,
            help="Number of parallel workers (default: available CPUs)",
        )
        parser.add_argument(
            "--dedup-mode", choices=["gen", "reco"], default="reco",
            help="Deduplication side for RecoMCTruthLink photon matching",
        )
        parser.add_argument(
            "--skip-gen-status-filter", action="store_true", default=False,
            help="Skip generatorStatus==1 filter on MCParticles",
        )
        parser.add_argument(
            "--max-gen-pdg", type=int, default=10000,
            help="Ignore gen particles with |PDG| > this value",
        )
        parser.add_argument(
            "--weight-mode", choices=["raw", "decoded"], default="decoded",
            help="How to interpret RecoMCTruthLink weights",
        )

    general_configs = myutils.setup_analysis_config(
        default_config, outputbasepath, parser_hook=_parser_hook
    )

    loggers        = general_configs["loggers"]
    run_config     = general_configs["config"]
    logger_config  = loggers["config"]
    logger_io      = loggers["io"]
    logger_process = loggers["processing"]

    args = general_configs["args"]

    dRMax         = run_config["cuts"]["dRMax"]
    minPTauPhoton = run_config["cuts"]["TauPhotonPCut"]
    minPTauPion   = run_config["cuts"]["TauPionPCut"]
    PNeutron      = run_config["cuts"]["NeutronCut"]
    dRMatch       = run_config["cuts"]["MatchedGenMinDR"]
    generalPCut   = run_config["cuts"]["generalPCut"]
    selectDecay   = general_configs["decay"]
    cut_string    = general_configs["decay_str"]
    sample        = run_config["general"]["sample"]
    outputpath    = general_configs["outputpath"]
    fileOutName   = os.path.join(outputpath, "Tree_" + general_configs["fileOutName"])
    test_arg      = general_configs["flags"]["test"]
    gatr_path     = general_configs["args"].gatr_result

    logger_config.info("Configuration loaded!")
    logger_config.info("Configuration:\n%s", pprint.pformat(general_configs, indent=4))

    filenames, mlpf_results = myutils.get_root_trees_path(
        sample, gatr_path, loggers, test_arg, args
    )
    if test_arg:
        logger_io.info("Test mode: limiting to 10 files.")
        filenames = filenames[:10]
    if not filenames:
        logger_io.error("No ROOT files found. Aborting.")
        sys.exit(1)
    logger_io.info("Read %d files", len(filenames))

    config_bundle = {
        "cuts": {
            "dRMax":         dRMax,
            "TauPhotonPCut": minPTauPhoton,
            "TauPionPCut":   minPTauPion,
            "NeutronCut":    PNeutron,
            "generalPCut":   generalPCut,
            "dRMatch":       dRMatch,
        },
        "selectDecay":       selectDecay,
        "cut_string":        cut_string,
        "sample":            sample,
        "genparts":          "MCParticles",
        "pfobjects":         "PandoraPFOs",
        "gatr_results_path": gatr_path,
        "test_pfo":          getattr(args, "test_pfo", False),
        "outputpath":        outputpath,
        "dedup_mode":        args.dedup_mode,
        "filter_gen_status": not args.skip_gen_status_filter,
        "max_gen_pdg":       args.max_gen_pdg,
        "weight_mode":       args.weight_mode,
    }

    n_workers   = args.n_workers or min(len(filenames), os.cpu_count() or 1)
    n_workers   = max(1, n_workers)
    file_chunks = split_filenames(filenames, n_workers)
    mlpf_chunks = split_mlpf(mlpf_results, file_chunks)

    event_offsets = []
    acc = 0
    for chunk in file_chunks:
        event_offsets.append(acc)
        acc += len(chunk) * 1000

    logger_io.info("Launching %d workers over %d files", n_workers, len(filenames))
    for i, chunk in enumerate(file_chunks):
        logger_io.info("  Worker %d: %d files, event offset %d", i, len(chunk), event_offsets[i])

    # Fork before any TFile is created in main to avoid ROOT state issues
    ctx = multiprocessing.get_context("fork")
    tmp_paths = []
    t_start = time.time()
    n_chunks = len(file_chunks)

    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
        futures = {
            executor.submit(
                process_chunk,
                file_chunks[i], mlpf_chunks[i],
                event_offsets[i], config_bundle, i,
            ): i
            for i in range(n_chunks)
        }
        for n_done, future in enumerate(as_completed(futures), start=1):
            wid = futures[future]
            elapsed = time.time() - t_start
            try:
                tmp_path = future.result()
                tmp_paths.append(tmp_path)
                logger_process.info(
                    "Worker %d done (%d/%d) | %.0fs elapsed", wid, n_done, n_chunks, elapsed
                )
            except Exception as exc:
                logger_process.error("Worker %d raised: %s", wid, exc)
            bar_len = 30
            filled  = int(bar_len * n_done / n_chunks)
            bar     = "█" * filled + "░" * (bar_len - filled)
            rem     = elapsed / n_done * (n_chunks - n_done)
            print(
                f"\r[{bar}] {n_done}/{n_chunks} | elapsed: {elapsed:.0f}s | remaining ~{rem:.0f}s",
                end="", flush=True,
            )
    print()
    total_elapsed = time.time() - t_start
    print(f"All workers done in {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")

    # Merge temporary files into the final output using hadd (fork-safe).
    # TChain::Merge from Python crashes after fork due to ROOT's Cling state
    # being corrupted in the parent process — hadd runs in a clean subprocess.
    logger_io.info("Merging %d temp files → %s", len(tmp_paths), fileOutName)
    hadd_cmd = ["hadd", "-f", fileOutName] + sorted(tmp_paths)
    result = subprocess.run(hadd_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger_io.error("hadd failed (exit %d):\n%s", result.returncode, result.stderr)
        sys.exit(1)
    logger_io.info("Output written: %s", fileOutName)

    for p in tmp_paths:
        try:
            os.remove(p)
        except Exception:
            pass

    # Save config snapshot
    out_cfg = os.path.join(outputpath, "config.yaml")
    lbl = run_config.get("output", {}).get("outputlabels")
    if not isinstance(lbl, list):
        run_config.setdefault("output", {})["outputlabels"] = [] if lbl is None else [lbl]
    with open(out_cfg, "w") as f:
        yaml.dump(run_config, f)
    logger_io.info("Config saved → %s", out_cfg)
    logger_io.info("End of job")


if __name__ == "__main__":
    main()
