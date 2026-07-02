from collections import defaultdict
import math
import pandas as pd
import numpy as np
import ROOT
from modules.myutils import dRAngle
from itertools import combinations
import os
import edm4hep
# setCharge(...)
#  |      void edm4hep::MutableReconstructedParticle::setCharge(float value)
#  |  
#  |  setCovMatrix(...)
#  |      void edm4hep::MutableReconstructedParticle::setCovMatrix(edm4hep::CovMatrix4f value)
#  |      void edm4hep::MutableReconstructedParticle::setCovMatrix(float value, edm4hep::FourMomCoords dimI, edm4hep::FourMomCoords dimJ)
#  |  
#  |  setDecayVertex(...)
#  |      void edm4hep::MutableReconstructedParticle::setDecayVertex(const edm4hep::Vertex& value)
#  |  
#  |  setEnergy(...)
#  |      void edm4hep::MutableReconstructedParticle::setEnergy(float value)
#  |  setMass(...)
#  |      void edm4hep::MutableReconstructedParticle::setMass(float value)
#  |  
#  |  setMomentum(...)
#  |      void edm4hep::MutableReconstructedParticle::setMomentum(edm4hep::Vector3f value)
#  |  
#  |  setPDG(...)
#  |      void edm4hep::MutableReconstructedParticle::setPDG(int32_t value)
DETECTOR_TYPES = {
    'INNER_TRACKER': 0,
    'ECAL': 1,
    'HCAL': 2,
    'MUON_TRACKER': 3
}
DETECTOR_ID_TO_TYPE = {
    0: 'INNER_TRACKER',
    1: 'ECAL',  
    2: 'HCAL',
    3: 'MUON_TRACKER'
}

# Estados de generador válidos (solo partículas estables finales)
VALID_GEN_STATUS = {1}

# PDG de neutrinos (ignorar en análisis gen)
NEUTRINO_PDGS = {12, 14, 16}

# Códigos PDG → nombre legible (sin distinguir partícula/antipartícula)
PDG_NAMES = {
    11: "epm",
    13: "μpm",
    22: "gamma",
    111: "π0",
    211: "πpm",
    321: "Kpm",
    2212: "p_bar_p",
    2112: "n_bar_n",
    130: "K_L0",
    310: "K_S0",
    15: "τpm",
}

# Colecciones de calorímetro con su tipo de detector
ECAL_COLLECTIONS = [
    ("ECALBarrel", DETECTOR_TYPES['ECAL']),
    ("ECALEndcap", DETECTOR_TYPES['ECAL']),
    ("ECALOther", DETECTOR_TYPES['ECAL']),
]
HCAL_COLLECTIONS = [
    ("HCALBarrel", DETECTOR_TYPES['HCAL']),
    ("HCALEndcap", DETECTOR_TYPES['HCAL']),
    ("HCALOther", DETECTOR_TYPES['HCAL']),
]
MUON_COLLECTIONS = [
    ("MUON", DETECTOR_TYPES['MUON_TRACKER']),
]
TRACK_COLLECTION = "SiTracks_Refitted"

ALL_CALO_COLLECTIONS = ECAL_COLLECTIONS + HCAL_COLLECTIONS + MUON_COLLECTIONS

c_light = 2.99792458e8
Bz_clic = 4.0
Bz_cld = 2.0
mchp = 0.139570

def omega_to_pt(omega, isclic):
    if isclic:
        Bz = Bz_clic
    else:
        Bz = Bz_cld
    a = c_light * 1e3 * 1e-15
    return a * Bz / abs(omega)

def get_track_momentum(trackstate, isclic=True):
    pt = omega_to_pt(trackstate.omega, isclic)
    phi = trackstate.phi
    pz = trackstate.tanLambda * pt
    px = pt * math.cos(phi)
    py = pt * math.sin(phi)
    p = math.sqrt(px * px + py * py + pz * pz)
    energy = math.sqrt(p * p + mchp * mchp)
    theta = math.acos(pz / p)
    return p, theta, phi, energy, px, py, pz

def pdg_name(pid):
    """Nombre legible para un código PDG."""
    pid_int = int(pid)
    return PDG_NAMES.get(pid_int, str(pid_int))


class ParticleWHits:
    
    def __init__(self, pfo, associated_hits):
        self.pfo = pfo
        self.associated_hits = associated_hits
        self.hits_tracker = None
        self.hits_ecal = None
        self.hits_hcal = None
        self.hits_muon = None
        self.separate_hits_by_detector()
    
    def separate_hits_by_detector(self):
        self.hits_tracker = []
        self.hits_ecal = []
        self.hits_hcal = []
        self.hits_muon = []
        for type_hit, hit_coords in self.associated_hits.items():
            if type_hit == str(DETECTOR_TYPES['INNER_TRACKER']):
                self.hits_tracker = hit_coords
            elif type_hit == str(DETECTOR_TYPES['ECAL']):
                self.hits_ecal = hit_coords
            elif type_hit == str(DETECTOR_TYPES['HCAL']):
                self.hits_hcal = hit_coords
            elif type_hit == str(DETECTOR_TYPES['MUON_TRACKER']):
                self.hits_muon = hit_coords
    
    def get_numpy_hits(self):
        
        tracker = np.array(self.hits_tracker) if self.hits_tracker else np.array([])
        ecal = np.array(self.hits_ecal) if self.hits_ecal else np.array([])
        hcal = np.array(self.hits_hcal) if self.hits_hcal else np.array([])
        muon = np.array(self.hits_muon) if self.hits_muon else np.array([])
        return {"0": tracker, "1": ecal, "2": hcal, "3": muon}

def debug_reco_tau(reco_tau, hit_type_map):
    constituents = reco_tau.getDaughters()
    particle_hit_collection = []
    for i, const in constituents.items():
        hit_collection = {"0": [], "1":[], "2": [], "3": []}
        clusters = const.getClusters()
        for cluster in clusters:
            hits = cluster.getHits()
            for hit in hits:
                hit_obj_id = hit.getObjectID()
                key = (hit_obj_id.collectionID, hit_obj_id.index)
                if key in hit_type_map:
                    detector_type = hit_type_map[key]
                    hit_x, hit_y, hit_z = hit.getPosition().x, hit.getPosition().y, hit.getPosition().z
                    hit_collection[str(detector_type)].append((hit_x, hit_y, hit_z))
        # print(hit_collection)
        particle_hits = ParticleWHits(const, hit_collection)
        particle_hit_collection.append(particle_hits)
    
    return particle_hit_collection
                    
def plot_debug_reco_tau(particle_hit_collection, non_assoc_tracks_df, output_dir, event_idx):
    # Plot hits and tracks in 3D
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    colors = ['r', 'g', 'b', 'y']
    for i, particle_hits in enumerate(particle_hit_collection):
        # Individual plot
        # big figure
        fig = plt.figure()
        fig.set_size_inches(10, 8)
        ax = fig.add_subplot(111, projection='3d')
        hits = particle_hits.get_numpy_hits()
        # print(hits)
        for det_type, hit_array in hits.items():
            if hit_array.size > 0:
                # color = colors[int(DETECTOR)]
                ax.scatter(hit_array[:, 0], hit_array[:, 1], hit_array[:, 2], c=colors[int(det_type)], label=f"Particle {i} - {DETECTOR_ID_TO_TYPE[int(det_type)]}")                                
        assoc_pfo = particle_hits.pfo
        if assoc_pfo.getPDG() == 2112:
            # Show tracks as vectors from origin with direction given by momentum
            track_info_list = []
            for track_row in non_assoc_tracks_df.itertuples():
                px = track_row.px
                py = track_row.py
                pz = track_row.pz
                # energy = track_row.energy
                # theta = track_row.theta
                # phi = track_row.phi
                # charge = track_row.charge
                # Scale momentum for visualization
                scale = 10  # Adjust as needed for better visualization
                ax.quiver(0, 0, 0, px*scale, py*scale, pz*scale, color='k', label=f"Track {track_row.track_idx}")
                track_info_list.append(
                    f"Track {track_row.track_idx} |\n "
                    f"p={track_row.p:.2f} GeV,\n "
                    f"θ={track_row.theta:.2f}, φ={track_row.phi:.2f},\n "
                    f"E={track_row.energy:.2f} GeV, \n"
                    f"q={track_row.charge}\n"
                )
            track_info = "\n".join(track_info_list)                            
            ax.text2D(0.65, 0.75, track_info,
                    transform=ax.transAxes,
                    fontsize=9,
                    bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.6))
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        output_file = f"{output_dir}/event_{event_idx}_particle_{i}.png"
        # Add info about the particle and the tracks
        p = math.sqrt(assoc_pfo.getMomentum().x**2 + assoc_pfo.getMomentum().y**2 + assoc_pfo.getMomentum().z**2)
        theta = math.acos(assoc_pfo.getMomentum().z / p) if p > 0 else 0
        phi = math.atan2(assoc_pfo.getMomentum().y, assoc_pfo.getMomentum().x)
        n_hits = sum(len(hit_array) for hit_array in hits.values())
        info_text = (
            f"Particle {i}\n"
            f"PDG: {assoc_pfo.getPDG()}\n"
            f"E: {assoc_pfo.getEnergy():.2f} GeV\n"
            f"|p|: {p:.2f} GeV\n"
            f"θ: {theta:.2f}\n"
            f"φ: {phi:.2f}\n"
            f"N hits: {n_hits}"
        )

        ax.text2D(0.02, 0.98, info_text,
                transform=ax.transAxes,
                fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))        
        plt.savefig(output_file)
        plt.close(fig)
    # Whole event plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    summary_info = "Whole event\n"
    for i, particle_hits in enumerate(particle_hit_collection):
        hits = particle_hits.get_numpy_hits()
        for det_type, hit_array in hits.items():
            if hit_array.size > 0:
                ax.scatter(hit_array[:, 0], hit_array[:, 1], hit_array[:, 2], c=colors[int(det_type)], label=f"Particle {i} - {DETECTOR_ID_TO_TYPE[int(det_type)]}")
        summary_info += f"Particle {i}: PDG {particle_hits.pfo.getPDG()}, E {particle_hits.pfo.getEnergy():.2f} GeV, N hits {sum(len(hit_array) for hit_array in hits.values())}\n"
    for track_row in non_assoc_tracks_df.itertuples():
        px = track_row.px
        py = track_row.py
        pz = track_row.pz
        scale = 10  # Adjust as needed for better visualization
        ax.quiver(0, 0, 0, px*scale, py*scale, pz*scale, color='k', label=f"Track {track_row.track_idx}")
        summary_info += f"Track {track_row.track_idx}: p {track_row.p:.2f} GeV, theta {track_row.theta:.2f}, phi {track_row.phi:.2f}, energy {track_row.energy:.2f} GeV, charge {track_row.charge}\n"
    ax.text2D(0.02, 0.98, summary_info,
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment='top',
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    output_file = f"{output_dir}/event_{event_idx}_whole_event.png"
    plt.savefig(output_file)
    plt.close(fig)
    
        
        
        
def process_one_neutron(pfos, idx_pfo_neutron,
                        hit_type_map,
                        non_associated_tracks_df,
                        event_idx,
                        df_reco_mc_links,
                        genParticles,
                        sitracks,
                        neutral_recover,
                        path_to_save,
                        save_extra_data=False,
                        cut_in_DR=0.1,
                        cut_in_energy=0.5,
                        apply_filter_in_energy=False,
                        logger_process=None):
    """Procesa un PFO candidato a neutrón, buscando tracks no asociados compatibles y comparando con partículas generadas.
    Guarda información detallada de la asociación para análisis posterior.
    """
  
    track_centroid_assoc = {"Pfo_idx": [],
                            "Track_idx": [],
                            "Ang_dist_hcal": [],
                            "Compatible_hcal": [],
                            "Ang_dist_ecal": [],
                            "Compatible_ecal": [],
                            "Ang_dist_total": [],
                            "Compatible_total": [],
                            "Energy_pfo": [],
                            "Energy_track": [],
                            "Compatible_energy": []}
    
    # Get pfo with idx_pfo_neutron
    for pfo in pfos:

        pfo_idx = pfo.getObjectID().index
        if pfo_idx == idx_pfo_neutron:
            break
        
    clusters = pfo.getClusters()
    for cluster in clusters:
        cluster_hits = cluster.getHits()
        hits_hcal = []
        hits_ecal = []
        for hit in cluster_hits:
            hit_obj_id = hit.getObjectID()
            key = (hit_obj_id.collectionID, hit_obj_id.index)
            if key in hit_type_map:
                detector_type = hit_type_map[key]
                hit_x = hit.getPosition().x
                hit_y = hit.getPosition().y
                hit_z = hit.getPosition().z
                if detector_type == DETECTOR_TYPES['ECAL']:
                    hits_ecal.append([hit_x, hit_y, hit_z])
                elif detector_type == DETECTOR_TYPES['HCAL']:
                    hits_hcal.append([hit_x, hit_y, hit_z])

    logger_process.debug(f"N hits HCAL: {len(hits_hcal)}, N hits ECAL: {len(hits_ecal)}")
    logger_process.debug(f"Mean position HCAL hits: {np.mean(hits_hcal, axis=0) if hits_hcal else 'N/A'}, Mean position ECAL hits: {np.mean(hits_ecal, axis=0) if hits_ecal else 'N/A'}")

    if hits_hcal:
        hits_hcal = np.array(hits_hcal)


    if hits_ecal:
        hits_ecal = np.array(hits_ecal)

    for track_row in non_associated_tracks_df.itertuples():
        si_track_idx = track_row.track_idx
        hcal_signal = len(hits_hcal) > 0
        ecal_signal = len(hits_ecal) > 0
        if ecal_signal and hcal_signal:
            total_hits = np.concatenate((hits_hcal, hits_ecal), axis=0)
            total_centroid = np.mean(total_hits, axis=0)
            ecal_centroid = np.mean(hits_ecal, axis=0)
            hcal_centroid = np.mean(hits_hcal, axis=0)
        elif hcal_signal:
            ecal_centroid = np.array([0, 0, 0])
            hcal_centroid = np.mean(hits_hcal, axis=0)
            total_centroid = np.mean(hits_hcal, axis=0)
        elif ecal_signal:
            hcal_centroid = np.array([0, 0, 0])
            ecal_centroid = np.mean(hits_ecal, axis=0)
            total_centroid = np.mean(hits_ecal, axis=0)
        
        
        track_p4 = ROOT.TLorentzVector()
        px = track_row.px
        py = track_row.py
        pz = track_row.pz
        track_p4.SetXYZM(px, py, pz, 0.0)
        track_centroid_assoc["Pfo_idx"].append(pfo.getObjectID().index)
        track_centroid_assoc["Track_idx"].append(si_track_idx)
        
        track_centroid_assoc["Energy_track"].append(track_row.energy)
        track_centroid_assoc["Energy_pfo"].append(pfo.getEnergy())
        track_centroid_assoc["Compatible_energy"].append(abs(track_row.energy - pfo.getEnergy())/pfo.getEnergy() < cut_in_energy)
        
        total_p4 = ROOT.TLorentzVector()
        total_p4.SetXYZM(total_centroid[0], total_centroid[1], total_centroid[2], 0.0)
        
        dR_total = dRAngle(track_p4, total_p4)                    
        track_centroid_assoc["Compatible_total"].append(dR_total < cut_in_DR)
        track_centroid_assoc["Ang_dist_total"].append(dR_total)
        
        if hcal_signal:
            hcal_p4 = ROOT.TLorentzVector()
            hcal_p4.SetXYZM(hcal_centroid[0], hcal_centroid[1], hcal_centroid[2], 0.0)
            dR = dRAngle(track_p4, hcal_p4)
            track_centroid_assoc["Ang_dist_hcal"].append(dR)
            track_centroid_assoc["Compatible_hcal"].append(dR < cut_in_DR)
        else:
            track_centroid_assoc["Ang_dist_hcal"].append(np.nan)
            track_centroid_assoc["Compatible_hcal"].append(False)
            
        if ecal_signal:
            ecal_p4 = ROOT.TLorentzVector()
            ecal_p4.SetXYZM(ecal_centroid[0], ecal_centroid[1], ecal_centroid[2], 0.0)
            dR_ecal = dRAngle(track_p4, ecal_p4)
            track_centroid_assoc["Ang_dist_ecal"].append(dR_ecal)
            track_centroid_assoc["Compatible_ecal"].append(dR_ecal < cut_in_DR)
        else:
            track_centroid_assoc["Ang_dist_ecal"].append(np.nan)
            track_centroid_assoc["Compatible_ecal"].append(False)

                
    trk_assoc_df = pd.DataFrame(track_centroid_assoc)
    if save_extra_data:
        trk_assoc_df.to_csv(f"{path_to_save}event_{event_idx}_pfo_{pfo.getObjectID().index}_track_centroid_assoc.csv", index=False)
    
    if apply_filter_in_energy:
      condition = trk_assoc_df["Compatible_energy"] & trk_assoc_df["Compatible_total"]
    else:
      condition = trk_assoc_df["Compatible_total"]
      
    if condition.any():
        # Check that only one track is compatible
        compatible_tracks = trk_assoc_df[condition]
        if len(compatible_tracks) == 1:
            # Obtain Track momentum
            logger_process.debug(f"Compatible track found for PFO {pfo.getObjectID().index}")
            track_idx = compatible_tracks.iloc[0]["Track_idx"]
            # Get track info from unassociated tracks dataframe
            track_info = non_associated_tracks_df[non_associated_tracks_df["track_idx"] == track_idx].iloc[0]
            p = track_info["p"]
            px = track_info["px"]
            py = track_info["py"]
            pz = track_info["pz"]
            energy = track_info["energy"]
            theta = track_info["theta"]
            phi = track_info["phi"]
            charge = track_info["charge"]
            # Look for the pion associated to the neutron in the gen particles
            pfo_idx = pfo.getObjectID().index
            
            # Get Generator Particle info
            gen_idx = df_reco_mc_links[df_reco_mc_links["reco"] == pfo_idx]["gen"].values
            if len(gen_idx) > 0:
                # Check if any is 211
                gen_pids = df_reco_mc_links[df_reco_mc_links["reco"] == pfo_idx]["Gen_pid"].values
                if 211 in gen_pids:
                    # Find the index of the pion (211)
                    pion_idx = np.where(gen_pids == 211)[0][0]
                    gen_idx = gen_idx[pion_idx]
                else:
                    # If no pion is found, set gen_idx to -1 or handle appropriately
                    gen_idx = gen_idx[0]
            gen_pid = df_reco_mc_links[df_reco_mc_links["reco"] == pfo_idx]["Gen_pid"].values[0]
            for gen_part in genParticles:
                if gen_part.getObjectID().index == gen_idx:
                    break
                    
            gen_momentum = gen_part.getMomentum()
            gen_px, gen_py, gen_pz = gen_momentum.x, gen_momentum.y, gen_momentum.z
            mass = gen_part.getMass()
            gen_p4 = ROOT.TLorentzVector()
            gen_p4.SetXYZM(gen_px, gen_py, gen_pz, mass)
            # Compare gen pion momentum with track momentum
            track_p4 = ROOT.TLorentzVector()
            track_p4.SetPxPyPzE(px, py, pz, energy)
            
            # Compare resolution in P
            res_p = (track_p4.P() - gen_p4.P()) / gen_p4.P()
            # Compare resolution in theta
            res_theta = (track_p4.Theta() - gen_p4.Theta()) / gen_p4.Theta()
            # Compare resolution in phi
            res_phi = (track_p4.Phi() - gen_p4.Phi()) / gen_p4.Phi()
            neutral_recover["event_idx"].append(event_idx)
            neutral_recover["gen_idx"].append(gen_idx)
            neutral_recover["reco_idx"].append(pfo_idx)
            neutral_recover["Track_idx"].append(track_idx)
            neutral_recover["gen_pid"].append(gen_pid)
            neutral_recover["gen_energy"].append(gen_part.getEnergy())
            neutral_recover["reco_energy"].append(pfo.getEnergy())
            neutral_recover["res_p"].append(res_p)
            neutral_recover["res_theta"].append(res_theta)
            neutral_recover["res_phi"].append(res_phi)
            associated_track_info = {
                "event_idx": event_idx,
                "pfo_idx": pfo_idx,
                "track_idx": track_idx,
                "gen_idx": gen_idx,
                "gen_pid": gen_pid,
                "track_charge": charge,
                "track_p": track_p4.P(),
                "track_theta": track_p4.Theta(),
                "track_phi": track_p4.Phi(),
                "track_energy": track_p4.E(),
                "gen_p": gen_p4.P(),
                "gen_theta": gen_p4.Theta(),
                "gen_phi": gen_p4.Phi(),
                "gen_energy": gen_p4.E(),
                "res_p": res_p,
                "res_theta": res_theta,
                "res_phi": res_phi}
            assoc_df = pd.DataFrame([associated_track_info])
            if save_extra_data:
                assoc_df.to_csv(f"{path_to_save}event_{event_idx}_pfo_{pfo.getObjectID().index}_track_assoc_detailed.csv", index=False)

        elif len(compatible_tracks) > 1:
            logger_process.debug(f"Multiple compatible tracks found for PFO {pfo.getObjectID().index}: {compatible_tracks} event idx: {event_idx}")
            # Obtain energy of the pfo
            # Try to combine two tracks, sum energy and compare with the pfo energy
            tracks_energy = []
            for track_idx in compatible_tracks["Track_idx"]:
                trackstate = sitracks[int(track_idx)].getTrackStates()[0]
                p, theta, phi, energy, px, py, pz = get_track_momentum(trackstate, isclic=False)
                tracks_energy.append(energy)
            # Try 2 by 2 combinations
            track_info_df = pd.DataFrame({
                "Track_idx": compatible_tracks["Track_idx"].values,
                "Energy": tracks_energy
            })
            if save_extra_data:
                track_info_df.to_csv(f"{path_to_save}event_{event_idx}_pfo_{pfo.getObjectID().index}_compatible_tracks_energy.csv", index=False)
            best_combination = None
            best_diff = float("inf")
            combs_results = {"Combination": [], "Combined_energy": [], "Energy_diff": []}
            for comb in combinations(range(len(tracks_energy)), 2):
                comb_energy = sum(tracks_energy[i] for i in comb)
                diff = abs(comb_energy - pfo.getEnergy())
                combs_results["Combination"].append(comb)
                combs_results["Combined_energy"].append(comb_energy)
                combs_results["Energy_diff"].append(diff)
                if diff < best_diff:
                    best_diff = diff
                    best_combination = comb
            if best_combination is not None:
                logger_process.debug(f"Best combination for PFO {pfo.getObjectID().index}: {best_combination}, energy difference: {best_diff}")
            else:
                logger_process.debug(f"No valid combination found for PFO {pfo.getObjectID().index}")

            # gen assocs
            associated_track_info = {
                "event_idx": [event_idx]*2,
                "pfo_idx": [pfo.getObjectID().index]*2,
                "track_idxs": best_combination,
                "track_energies": [tracks_energy[i] for i in best_combination],
                "pfo_energy": [pfo.getEnergy()]*2,}
            assoc_df = pd.DataFrame(associated_track_info)
            combs_results_df = pd.DataFrame(combs_results)
            if save_extra_data:
                assoc_df.to_csv(f"{path_to_save}event_{event_idx}_pfo_{pfo.getObjectID().index}_track_assoc_combination_detailed.csv", index=False)
                combs_results_df.to_csv(f"{path_to_save}event_{event_idx}_pfo_{pfo.getObjectID().index}_track_assoc_combination_results.csv", index=False)
    else:
        logger_process.debug(f"No compatible tracks found for PFO {pfo.getObjectID().index}")

def build_hit_type_map(event):
    """
    Construye un diccionario que mapea (collectionID, index) -> detector_type.
    Esto permite identificar el tipo de cada hit cuando se accede desde PFO clusters.
    """
    hit_type_map = {}  # (collectionID, index) -> detector_type
    hit_energy_map = {}  # (collectionID, index) -> energy
    
    for coll_name, detector_type in ALL_CALO_COLLECTIONS:
        try:
            coll = event.get(coll_name)
            for hit in coll:
                obj_id = hit.getObjectID()
                key = (obj_id.collectionID, obj_id.index)
                hit_type_map[key] = detector_type
                hit_energy_map[key] = hit.getEnergy()
        except Exception:
            pass
    
    return hit_type_map, hit_energy_map

# def analyze_pfos(event,df_reco_mc_links, hit_type_map, hit_energy_map, verbose=False, logger_process=None):
#     """
#     Analiza los PandoraPFOs (partículas reconstruidas).
    
#     Para cada PFO cuenta los hits en ECAL, HCAL, tracks y muon usando el mapeo previo.
#     """
#     try:
#         pfos = event.get("PandoraPFOs")
#     except Exception:
#         logger_process.error("    No se encontró colección PandoraPFOs")
#         return []
    
#     particles = []
    
#     for pfo_idx, pfo in enumerate(pfos):
#         try:
#             # Información básica del PFO
#             momentum = pfo.getMomentum()
#             px, py, pz = momentum.x, momentum.y, momentum.z
#             p = math.sqrt(px**2 + py**2 + pz**2)
#             energy = pfo.getEnergy()
#             pid = pfo.getPDG()
#             charge = pfo.getCharge()
            
#             theta = math.acos(pz / p) if p > 0 else 0
#             phi = math.atan2(py, px)
            
#             # Usar abs(pid) para no distinguir partícula/antipartícula
#             pid = abs(pid)
            
#             # Contar tracks asociados
#             tracks = pfo.getTracks()
#             n_tracks = len(tracks) if tracks else 0
            
#             # Contar hits de calorímetro por tipo usando el mapeo
#             clusters = pfo.getClusters()

#             n_ecal = 0
#             n_hcal = 0
#             n_muon = 0
#             E_ecal = 0.0
#             E_hcal = 0.0
#             E_muon = 0.0
            
#             for cluster in clusters:
#                 cluster_hits = cluster.getHits()
                
#                 for hit in cluster_hits:
#                     hit_obj_id = hit.getObjectID()
#                     key = (hit_obj_id.collectionID, hit_obj_id.index)
                    
#                     # Buscar el tipo de detector en el mapeo
#                     if key in hit_type_map:
#                         detector_type = hit_type_map[key]
#                         hit_energy = hit_energy_map.get(key, hit.getEnergy())
                        
#                         if detector_type == DETECTOR_TYPES['ECAL']:
#                             n_ecal += 1
#                             E_ecal += hit_energy
#                         elif detector_type == DETECTOR_TYPES['HCAL']:
#                             n_hcal += 1
#                             E_hcal += hit_energy
#                         elif detector_type == DETECTOR_TYPES['MUON_TRACKER']:
#                             n_muon += 1
#                             E_muon += hit_energy
            
#             # Calcular ratio
#             ratio_hcal_ecal = E_hcal / E_ecal if E_ecal > 0 else np.nan
            
#             # Obtener información de matching desde df_reco_mc_links
#             match_rows = df_reco_mc_links[df_reco_mc_links["reco"] == pfo_idx]
#             if not match_rows.empty:
#                 matched_gen_idx = match_rows.iloc[0]["gen"]
#                 is_matched = matched_gen_idx != -999
#             else:
#                 matched_gen_idx = -999
#                 is_matched = False
            
#             particles.append({
#                 "pfo_id": pfo_idx,
#                 "pid": pid,
#                 "pid_name": pdg_name(pid),
#                 "charge": charge,
#                 "E_reco": energy,
#                 "p": p,
#                 "theta": theta,
#                 "phi": phi,
#                 "n_track": n_tracks,
#                 "n_ecal": n_ecal,
#                 "n_hcal": n_hcal,
#                 "n_muon": n_muon,
#                 "n_total": n_tracks + n_ecal + n_hcal + n_muon,
#                 "E_ecal": E_ecal,
#                 "E_hcal": E_hcal,
#                 "E_muon": E_muon,
#                 "ratio_hcal_ecal": ratio_hcal_ecal,
#                 "is_matched": is_matched,
#                 "matched_gen_idx": matched_gen_idx
#             })
            
#         except Exception as e:
#               logger_process.error(f"    Error procesando PFO {pfo_idx}: {e}")
#     # exit(0)
#     return particles


def get_reco_mc_links_by_dR(event, hit_type_map, hit_energy_map, logger_process=None):
    """
    Asocia partículas gen con partículas reco usando distancia dR en el plano (theta, phi).
    
    Para cada partícula gen (status 1, sin neutrinos, con señal en detector):
    - Encuentra la partícula reco más cercana usando dRAngle(p1, p2)
    - Ignora la carga
    
    Devuelve un DataFrame similar a get_reco_mc_links pero con asociaciones basadas en dR.
    """
    try:
        mc_particles = event.get("MCParticles")
    except Exception:
        return pd.DataFrame()
    
    try:
        pfos = event.get("PandoraPFOs")
    except Exception:
        logger_process.error("    No se encontró colección PandoraPFOs")
        return pd.DataFrame()
    
    # Construir mc_stats para filtrar partículas sin señal
    try:
        calo_truth_links = list(event.get("CalohitMCTruthLink"))
    except Exception:
        calo_truth_links = []
    
    try:
        track_truth_links = list(event.get("SiTracksMCTruthLink"))
    except Exception:
        track_truth_links = []
    
    mc_stats = defaultdict(lambda: {'n_track': 0, 'n_ecal': 0, 'n_hcal': 0, 'n_muon': 0})
    
    for link in calo_truth_links:
        try:
            hit_obj = link.getRec()
            mc_obj = link.getSim()
            mc_idx = mc_obj.getObjectID().index
            hit_obj_id = hit_obj.getObjectID()
            key = (hit_obj_id.collectionID, hit_obj_id.index)
            if key in hit_type_map:
                detector_type = hit_type_map[key]
                if detector_type == DETECTOR_TYPES['ECAL']:
                    mc_stats[mc_idx]['n_ecal'] += 1
                elif detector_type == DETECTOR_TYPES['HCAL']:
                    mc_stats[mc_idx]['n_hcal'] += 1
                elif detector_type == DETECTOR_TYPES['MUON_TRACKER']:
                    mc_stats[mc_idx]['n_muon'] += 1
        except Exception:
            pass
    
    for link in track_truth_links:
        try:
            mc_obj = link.getSim()
            mc_idx = mc_obj.getObjectID().index
            mc_stats[mc_idx]['n_track'] += 1
        except Exception:
            pass
    all_0 = True    
    for mc_id in mc_stats.keys():
        for ke in mc_stats[mc_id]:
            if mc_stats[mc_idx][ke]!=0:
                all_0 = False
            
    # Construir lista de partículas gen válidas (status 1, sin neutrinos, con señal)
    valid_gen_particles = []
    # print("Obtaining gen particles")
    for idex, part in enumerate(mc_particles):
        # print(idex)
        try:
            gen_status = part.getGeneratorStatus()
            # print(gen_status)
            if gen_status not in VALID_GEN_STATUS:
                continue
            
            pid_raw = part.getPDG()
            if abs(pid_raw) in NEUTRINO_PDGS:
                continue
            # print(pid_raw)
            
            momentum = part.getMomentum()
            p = math.sqrt(momentum.x**2 + momentum.y**2 + momentum.z**2)
            if p < 1e-10:
                continue
            
            idx = part.getObjectID().index
            
            # Filtrar partículas sin señal en el detector
            stats = mc_stats.get(idx, {'n_track': 0, 'n_ecal': 0, 'n_hcal': 0, 'n_muon': 0})
            n_total = stats['n_track'] + stats['n_ecal'] + stats['n_hcal'] + stats['n_muon']
            if n_total == 0 and not all_0:
                continue  # No interacciona con el detector
            
            valid_gen_particles.append((idx, part, momentum, pid_raw))
        except Exception:
            pass
    # print("Valid gen parts", valid_gen_particles)
    # Construir lista de PFOs
    valid_pfos = []
    for pfo in pfos:
        try:
            momentum = pfo.getMomentum()
            px, py, pz = momentum.x, momentum.y, momentum.z
            p = math.sqrt(px**2 + py**2 + pz**2)
            if p < 1e-10:
                continue
            
            pfo_idx = pfo.getObjectID().index
            pid = pfo.getPDG()
            energy = pfo.getEnergy()
            valid_pfos.append((pfo_idx, momentum, pid, energy))
        except Exception:
            pass
    # print("Valid pfos", valid_pfos)
    
    # Para cada partícula gen, encontrar la reco más cercana
    reco_gen_link = {"gen": [], "reco": [], "dR": [], "Gen_pid": [], "Reco_pid": [], "Gen_energy": [], "Reco_energy": [],
                     "Gen_Px": [], "Gen_Py": [], "Gen_Pz": [], "Reco_Px": [], "Reco_Py": [], "Reco_Pz": []}
    
    for gen_idx, gen_part, gen_momentum, gen_pid in valid_gen_particles:
        # Crear TLorentzVector para gen
        gen_p4 = ROOT.TLorentzVector()
        gen_p4.SetXYZT(gen_momentum.x, gen_momentum.y, gen_momentum.z, gen_part.getEnergy())
        
        min_dR = 0.1
        best_reco_idx = -999
        best_reco_pid = -999
        best_reco_p4 = None
        # Buscar el reco más cercano
        for reco_idx, reco_momentum, reco_pid, reco_energy in valid_pfos:
            # Crear TLorentzVector para reco
            # Nota: no tenemos energía directa del momentum, así que usamos E = sqrt(p^2 + m^2)
            # Para simplificar, asumimos masa 0 (fotón) o calculamos E ~ p para relativistas
            # p_reco = math.sqrt(reco_momentum.x**2 + reco_momentum.y**2 + reco_momentum.z**2)
            reco_p4 = ROOT.TLorentzVector()
            reco_p4.SetXYZM(reco_momentum.x, reco_momentum.y, reco_momentum.z, 0.0)
            
            # Calcular dR
            dR = dRAngle(gen_p4, reco_p4)
            
            if dR < min_dR:
                min_dR = dR
                best_reco_idx = reco_idx
                best_reco_pid = reco_pid
                best_reco_p4 = reco_p4
        
        # Agregar asociación
        reco_gen_link["gen"].append(gen_idx)
        reco_gen_link["reco"].append(best_reco_idx)
        reco_gen_link["dR"].append(min_dR)
        reco_gen_link["Gen_pid"].append(gen_pid)
        reco_gen_link["Reco_pid"].append(best_reco_pid)

        reco_gen_link["Gen_energy"].append(gen_p4.P())
        reco_gen_link["Gen_Px"].append(gen_momentum.x)
        reco_gen_link["Gen_Py"].append(gen_momentum.y)
        reco_gen_link["Gen_Pz"].append(gen_momentum.z)

        reco_gen_link["Reco_energy"].append(best_reco_p4.P() if best_reco_p4 else np.nan)
        reco_gen_link["Reco_Px"].append(best_reco_p4.Px() if best_reco_p4 else np.nan)
        reco_gen_link["Reco_Py"].append(best_reco_p4.Py() if best_reco_p4 else np.nan)
        reco_gen_link["Reco_Pz"].append(best_reco_p4.Pz() if best_reco_p4 else np.nan)

    df_reco_mc_links = pd.DataFrame(reco_gen_link)
    # print("df_reco_mc_links", df_reco_mc_links) 
    if df_reco_mc_links.empty:
        return df_reco_mc_links
    
    df_reco_mc_links = df_reco_mc_links.sort_values("gen").reset_index(drop=True)
    
    # Para cada gen, mantener solo el match con mínimo dR
    df_reco_mc_links = df_reco_mc_links.loc[df_reco_mc_links.groupby("gen")["dR"].idxmin()].reset_index(drop=True)
    
    # Añadir gen sin match (los que no están en la lista)
    gen_indices = set(df_reco_mc_links["gen"])
    new_rows = []
    for gen_idx, gen_part, gen_momentum, gen_pid in valid_gen_particles:
        if gen_idx not in gen_indices:
            new_rows.append({"gen": gen_idx, "reco": -999, "dR": np.nan, "Gen_pid": gen_pid, "Reco_pid": -999,
                              "Gen_energy": gen_part.getEnergy(),
                              "Gen_Px": gen_momentum.x, "Gen_Py": gen_momentum.y, "Gen_Pz": gen_momentum.z,
                              "Reco_Px": np.nan, "Reco_Py": np.nan, "Reco_Pz": np.nan})
    
    # Añadir reco sin match (los que no están en la lista)
    reco_indices = set(df_reco_mc_links["reco"])
    for reco_idx, reco_momentum, reco_pid, reco_energy in valid_pfos:
        if reco_idx not in reco_indices:
            new_rows.append({"gen": -999, "reco": reco_idx, "dR": np.nan, "Gen_pid": -999, "Reco_pid": reco_pid,
                              "Reco_energy": reco_energy,
                              "Gen_Px": np.nan, "Gen_Py": np.nan, "Gen_Pz": np.nan,
                              "Reco_Px": reco_momentum.x, "Reco_Py": reco_momentum.y, "Reco_Pz": reco_momentum.z})
    
    if new_rows:
        df_reco_mc_links = pd.concat([df_reco_mc_links, pd.DataFrame(new_rows)], ignore_index=True)
    
    # Eliminar columna dR para mantener compatibilidad con el resto del código
    df_reco_mc_links = df_reco_mc_links.drop(columns=["dR"])
    
    return df_reco_mc_links

def recover_pion_from_neutrals(particles_pfos, event, event_idx, logger_process, neutral_recover_cfg):
    
    """
    Processes a ROOT file and returns DataFrames of gen and reco particles.
    Uses associations by dR distance instead of RecoMCTruthLink.
    
    arguments:
- particles_pfos: list of reconstructed PFOs in the event
- event: ROOT event to process
- event_idx: event index (for logging and saving)
- logger_process: logger for debug/info messages
- neutral_recover_cfg: configuration for neutral recovery, with options such as:
    - save_extra_data: bool, if True saves intermediate DataFrames for further analysis
    - cut_in_DR: float, dR threshold to consider an association compatible
    - cut_in_energy: float, energy difference threshold to consider compatible
    - apply_filter_in_energy: bool, if True requires energy compatibility in addition to dR to consider an association valid
returns:
- pfos: modified pfos with information of recovered tracks (if compatible one is found)
- neutral_recover_df: DataFrame with information of associations found between neutron PFOs and non-associated tracks, including resolutions in p, theta, phi with respect to particles

    """
    
    neutral_recover = {"event_idx":[], "gen_idx": [], "reco_idx": [], "Track_idx": [], "gen_pid": [], "gen_energy": [], "reco_energy":[], "res_p": [], "res_theta": [], "res_phi": []}
        
    # Construir mapeo de hits una vez por evento
    hit_type_map, hit_energy_map = build_hit_type_map(event)
    
    # Usar asociaciones por dR
    df_reco_mc_links = get_reco_mc_links_by_dR(event, hit_type_map, hit_energy_map, logger_process=logger_process)
    
    pfos = particles_pfos
    genParticles = event.get("MCParticles")
    non_associated_tracks_df = pd.DataFrame()
    try:
        if (df_reco_mc_links["Reco_pid"] == 2112).any(): # neutron case
        # idx where reco_id is 2112
            # idx_reco_neutron = df_reco_mc_links[df_reco_mc_links["Reco_pid"] == 2112].index
            associated_tracks = []
            for pfo in pfos:
                # Get track
                pfo_track = pfo.getTracks()
                n_tracks = len(pfo_track) if pfo_track else 0

                for i in range(n_tracks):
                    associated_tracks.append(pfo_track[i].getObjectID().index)

            sitracks = event.get("SiTracks_Refitted")
            non_associated_tracks = {"track_idx":[], "p":[], "px":[], "py":[], "pz":[], "theta":[], "phi":[], "energy":[], "charge":[]}
            for i in range(len(sitracks)):
                si_track_idx = sitracks[i].getObjectID().index
                if si_track_idx not in associated_tracks:
                    trackstate = sitracks[i].getTrackStates()[0]
                    p, theta, phi, energy, px, py, pz = get_track_momentum(trackstate, isclic=False)
                    non_associated_tracks["track_idx"].append(si_track_idx)
                    non_associated_tracks["p"].append(p)
                    non_associated_tracks["px"].append(px)
                    non_associated_tracks["py"].append(py)
                    non_associated_tracks["pz"].append(pz)
                    non_associated_tracks["theta"].append(theta)
                    non_associated_tracks["phi"].append(phi)
                    non_associated_tracks["energy"].append(energy)
                    non_associated_tracks["charge"].append(-1 if trackstate.omega < 0 else 1)
            non_associated_tracks_df = pd.DataFrame(non_associated_tracks)
            # n unique pfo id with pid 2112
            idx_pfo_neutron = df_reco_mc_links[df_reco_mc_links["Reco_pid"] == 2112]["reco"].unique()
            
            
            
            save_extra_data = neutral_recover_cfg.get("save_extra_data", False)
            cut_in_DR = neutral_recover_cfg.get("cut_in_DR", 0.1)
            cut_in_energy = neutral_recover_cfg.get("cut_in_energy", 0.5)
            apply_filter_in_energy = neutral_recover_cfg.get("apply_filter_in_energy", False)
            
            if save_extra_data:
                base_save_dir = neutral_recover_cfg.get("save_dir", "neutral_recover_extra_data/")
                save_dir = f"{base_save_dir}event_{event_idx}/"
                if not os.path.exists(save_dir):
                    os.makedirs(save_dir)
            else:
                save_dir = None
                
            logger_process.debug(f"IDX PFO neutron: {idx_pfo_neutron}, Event idx: {event_idx}")
            if len(idx_pfo_neutron) == 1 and not non_associated_tracks_df.empty:

                if save_extra_data:    
                    df_reco_mc_links.to_csv(f"{save_dir}event_{event_idx}_reco_mc_links.csv", index=False)
                
                # One reco neutron case
                process_one_neutron(pfos,
                                    idx_pfo_neutron[0],
                                    hit_type_map,
                                    non_associated_tracks_df,
                                    event_idx,
                                    df_reco_mc_links,
                                    genParticles,
                                    sitracks,
                                    neutral_recover,
                                    save_dir,
                                    save_extra_data=save_extra_data,
                                    cut_in_DR=cut_in_DR,
                                    cut_in_energy=cut_in_energy,
                                    apply_filter_in_energy=apply_filter_in_energy,
                                    logger_process=logger_process)
                
            elif len(idx_pfo_neutron) > 1 and not non_associated_tracks_df.empty:
                if save_extra_data:
                    df_reco_mc_links.to_csv(f"{save_dir}event_{event_idx}_reco_mc_links.csv", index=False)
                for idx in idx_pfo_neutron:
                    process_one_neutron(pfos, idx,
                                        hit_type_map,
                                        non_associated_tracks_df,
                                        event_idx,
                                        df_reco_mc_links,
                                        genParticles,
                                        sitracks,
                                        neutral_recover,
                                        save_dir,
                                        save_extra_data=save_extra_data,
                                        cut_in_DR=cut_in_DR,
                                        cut_in_energy=cut_in_energy,
                                        apply_filter_in_energy=apply_filter_in_energy,
                                        logger_process=logger_process)
                
        neutral_recover_df = pd.DataFrame(neutral_recover)
        
        
        # add file name to every row of neutral_recover_df
        neutral_recover_df["event_idx"] = event_idx

        
        # Change pfo level information
        changes_info = {"event": [], "reco_idx":[], "track_idx":[]}
        reco_unique_reco_pion_neutron_idx = neutral_recover_df["reco_idx"].unique()
        for pfo in pfos:
            pfo_idx = pfo.getObjectID().index
            if pfo_idx in reco_unique_reco_pion_neutron_idx:
                # Track momentum
                assoc_rows = neutral_recover_df[neutral_recover_df["reco_idx"] == pfo_idx]
                if assoc_rows.empty:
                    continue
                assoc_track = assoc_rows["Track_idx"].iloc[0]
                # Track energy and momentum
                track_matches = non_associated_tracks_df[non_associated_tracks_df["track_idx"] == assoc_track]
                if track_matches.empty:
                    continue
                track_info = track_matches.iloc[0]
                pfo.setEnergy(track_info["energy"])
                px = track_info["px"]
                py = track_info["py"]
                pz = track_info["pz"]
                pfo.setMomentum(edm4hep.Vector3f(px, py, pz))
                pfo.setEnergy(track_info["energy"])
                pfo.setMass(mchp)   # masa del pion cargado
                charge = track_info["charge"]
                pfo.setPDG(int(charge*211))
                pfo.setCharge(charge)
                changes_info["reco_idx"].append(pfo_idx)
                changes_info["track_idx"].append(assoc_track)
                changes_info["event"].append(event_idx)
    except Exception as e:
        logger_process.debug(f"Error en asociación de neutrones: {e}")
        changes_info = {"event": [], "reco_idx":[], "track_idx":[]}
        neutral_recover_df = pd.DataFrame()
    
    extra_info_dict = {"non_associated_tracks_df": non_associated_tracks_df, "hit_type_map": hit_type_map, "df_reco_mc_links": df_reco_mc_links}
    return pfos, changes_info, extra_info_dict



def _run_main_test(root_path, max_events=None):
    try:
        from podio import root_io
    except ImportError as exc:
        raise RuntimeError(
            "No se pudo importar podio.root_io. Ejecuta este test en un entorno Key4HEP."
        ) from exc
    import logging
    
    logging.basicConfig(level=logging.DEBUG, format="[%(levelname)s] %(message)s")
    logger_process = logging.getLogger("NeutralRecoverTest")
    
    reader = root_io.Reader(root_path)
    events = list(reader.get("events"))
    total_events = len(events)
    if total_events == 0:
        raise RuntimeError(f"El archivo no contiene eventos en la colección 'events': {root_path}")

    neutral_recover_cfg = {
        "save_extra_data": True,
        "cut_in_DR": 0.1,
        "cut_in_energy": 0.5,
        "apply_filter_in_energy": False,
        "save_dir": "neutral_recover_images_test/"
    }

    print(f"[NeutralRecover test] Archivo: {root_path}")
    print(f"[NeutralRecover test] Nº eventos totales: {total_events}")
    print(f"[NeutralRecover test] Eventos a probar: {max_events}")
    print(f"[NeutralRecover test] Config: {neutral_recover_cfg}")

    failures = []
    for i, event in enumerate(events):
        if i >= max_events:
            break
        event_idx = i
        try:

            pfos = event.get("PandoraPFOs")
            new_pfos, changes_info, extra_info_dict = recover_pion_from_neutrals(
                pfos,
                event,
                event_idx,
                logger_process,
                neutral_recover_cfg,
            )

            if not isinstance(changes_info, dict):
                raise RuntimeError("changes_info no es un dict")
            for key in ("event", "reco_idx", "track_idx"):
                if key not in changes_info:
                    raise RuntimeError(f"Falta la clave '{key}' en changes_info")

            print(
                f"  - Evento {event_idx}: "
                f"PFOs={len(pfos)} "
                f"cambios={changes_info}"
            )
            for idx in changes_info["reco_idx"]:
              for pfo in new_pfos:
                if pfo.getObjectID().index == idx:
                    print(f"    - PFO {idx} nuevo PDG: {pfo.getPDG()}")
                    break
        except Exception as exc:
            failures.append((event_idx, str(exc)))
            print(f"  - Evento {event_idx}: ERROR -> {exc}")

    if failures:
        summary = "; ".join([f"evento {idx}: {msg}" for idx, msg in failures])
        raise RuntimeError(f"Fallaron {len(failures)} evento(s): {summary}")

    print(
        "[NeutralRecover test] OK: "
        f"Procesados {len(events)} eventos, fallaron {len(failures)} eventos."
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Test de recover_pion_from_neutrals en eventos límite de un archivo ROOT EDM4HEP."
    )
    parser.add_argument(
        "--root-file",
        required=True,
        help="Ruta al archivo .root con eventos de prueba.",
    )
    parser.add_argument(
        "--max-events",
        default=999999,
        type=int,
        help="Lista opcional de eventos separados por coma (ej: 0,4,8). Si no se indica, usa eventos límite.",
    )

    args = parser.parse_args()

    _run_main_test(args.root_file, max_events=args.max_events)