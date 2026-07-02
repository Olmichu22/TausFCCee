from modules import tauReco, electronReco, muonReco, NeutralRecover
import ROOT
import edm4hep
from typing import Optional

def extractTauDecays(gatr_results_path,
                     mlpf_results,
                     eventid,
                     pfos,
                     dRMax,
                     minPTauPhoton,
                     minPTauPion,
                     PNeutron,
                     generalPCut,
                     foton_config,
                     test_extremes,
                     test_pfo,
                     logger_process,
                     neutral_recover_cfg: dict=dict(),
                     event=None,
                     only_association=False):
  if gatr_results_path is not None and not test_pfo:
    logger_process.debug("Using GATr results for event %d", eventid)
    particles = mlpf_results.get(eventid, {})
    charge_condition = False
  else:
    logger_process.debug("Using PandoraPFO results for event %d", eventid)
    particles = pfos
    charge_condition = True
  
  # print(neutral_recover_cfg)
  # exit(0)
  if neutral_recover_cfg.get("enable", False):
      logger_process.debug("Applying neutron recovery for event %d", eventid)
      if not only_association:
        particles, recover_extra_info, extra_info_dict = NeutralRecover.recover_pion_from_neutrals(particles, event, eventid, logger_process, neutral_recover_cfg)
      else:
        _, recover_extra_info, extra_info_dict = NeutralRecover.recover_pion_from_neutrals(particles, event, eventid, logger_process, neutral_recover_cfg)

        
  recoTau = tauReco.findAllTaus(particles,
                                dRMax,
                                minPTauPhoton,
                                minPTauPion,
                                PNeutron,
                                generalPCut,
                                charge_condition=charge_condition)
  recoElectrons = electronReco.findAllElectrons(particles, generalPCut)
  recoMuons = muonReco.findAllMuons(particles, generalPCut)
  
  
  def test_extremes_for_photons(particles, test_pfo, gatr_results_path, foton_config, syst_type="energy"):
    if test_pfo or gatr_results_path is None:
      particles_w_extremes_max = [p.clone() for p in particles]
      particles_w_extremes_min = [p.clone() for p in particles]
    else:
      particles_w_extremes_max = [p.copy() for p in particles]
      particles_w_extremes_min = [p.copy() for p in particles]
    for ind, part in enumerate(particles):
        pdg = abs(part.getPDG())
    
        if pdg == 22:
          if test_pfo or gatr_results_path is None:
            # If using PandoraPFOs, need to build TLorentzVector
            p4 = ROOT.TLorentzVector()
            p4.SetXYZM(part.getMomentum().x, part.getMomentum().y, part.getMomentum().z, part.getMass())
          else:
            # Using GATr results, can use directly because part is a RecoParticle
            p4 = part.getMomentum()
          if syst_type == "energy":
            p4_max, p4_min = tauReco.electromagnetic_energy_error_p4_extremes(p4, foton_config)
            if test_pfo or gatr_results_path is None:
              # Update momentum in the cloned PandoraPFOs
              new_p_max = edm4hep.Vector3f()
              new_p_max.x = p4_max.X()
              new_p_max.y = p4_max.Y()
              new_p_max.z = p4_max.Z()
              particles_w_extremes_max[ind].setMomentum(new_p_max)
              
              new_p_min = edm4hep.Vector3f()
              new_p_min.x = p4_min.X()
              new_p_min.y = p4_min.Y()
              new_p_min.z = p4_min.Z()
              particles_w_extremes_min[ind].setMomentum(new_p_min)
              
            else:
              particles_w_extremes_max[ind].setMomentum(p4_max)
              particles_w_extremes_min[ind].setMomentum(p4_min)
              
          elif syst_type == "direction":
            p4_smeared = tauReco.electromagnetic_direction_error_p4_extremes(p4, foton_config)
            if test_pfo or gatr_results_path is None:
              # Update momentum in the cloned PandoraPFOs
              new_p_smeared = edm4hep.Vector3f()
              new_p_smeared.x = p4_smeared.X()
              new_p_smeared.y = p4_smeared.Y()
              new_p_smeared.z = p4_smeared.Z()
              particles_w_extremes_max[ind].setMomentum(new_p_smeared)
              # particles_w_extremes_min[ind].setMomentum(new_p_smeared)
            else:
              particles_w_extremes_max[ind].setMomentum(p4_smeared)
              # particles_w_extremes_min[ind].setMomentum(p4_smeared)
          else:
            raise ValueError(f"Invalid syst_type {syst_type} for testing extremes of photons")
          if syst_type == "energy":
            logger_process.debug(f"Momentum changes max: {p4_max.P()} min: {p4_min.P()}")
          elif syst_type == "direction":
            logger_process.debug(f"Direction changes: {p4_smeared.Theta()}, {p4_smeared.Phi()}")
    
    if syst_type == "energy":
      recoTau_max = tauReco.findAllTaus(
          particles_w_extremes_max, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut, charge_condition=charge_condition
      )
      recoTau_min = tauReco.findAllTaus(
          particles_w_extremes_min, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut, charge_condition=charge_condition
      )
    elif syst_type == "direction":
      recoTau_max = tauReco.findAllTaus(
          particles_w_extremes_max, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut, charge_condition=charge_condition
      )
      recoTau_min = None
    return recoTau_max, recoTau_min
  
  # # Check extremes of photon energy
  # if test_extremes:
  #     if test_pfo or gatr_results_path is None:
  #       particles_w_extremes_max = [p.clone() for p in particles]
  #       particles_w_extremes_min = [p.clone() for p in particles]
  #     else:
  #       particles_w_extremes_max = [p.copy() for p in particles]
  #       particles_w_extremes_min = [p.copy() for p in particles]
  #     for ind, part in enumerate(particles):
  #         pdg = abs(part.getPDG())
      
  #         if pdg == 22:
  #           if test_pfo or gatr_results_path is None:
  #             # If using PandoraPFOs, need to build TLorentzVector
  #             p4 = ROOT.TLorentzVector()
  #             p4.SetXYZM(part.getMomentum().x, part.getMomentum().y, part.getMomentum().z, part.getMass())
  #           else:
  #             # Using GATr results, can use directly because part is a RecoParticle
  #             p4 = part.getMomentum()
  #           if foton_config.get("energy", {}):
  #             p4_max, p4_min = tauReco.electromagnetic_energy_error_p4_extremes(p4, foton_config.get("energy", {}))
  #           if test_pfo or gatr_results_path is None:
  #             # Update momentum in the cloned PandoraPFOs
  #             new_p_max = edm4hep.Vector3f()
  #             new_p_max.x = p4_max.X()
  #             new_p_max.y = p4_max.Y()
  #             new_p_max.z = p4_max.Z()
  #             particles_w_extremes_max[ind].setMomentum(new_p_max)
              
  #             new_p_min = edm4hep.Vector3f()
  #             new_p_min.x = p4_min.X()
  #             new_p_min.y = p4_min.Y()
  #             new_p_min.z = p4_min.Z()
  #             particles_w_extremes_min[ind].setMomentum(new_p_min)
              
  #           else:
  #             particles_w_extremes_max[ind].setMomentum(p4_max)
  #             particles_w_extremes_min[ind].setMomentum(p4_min)
      
  #           logger_process.debug(f"Momentum changes max: {p4_max.P()} min: {p4_min.P()}")
  #     recoTau_max = tauReco.findAllTaus(
  #         particles_w_extremes_max, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut, charge_condition=charge_condition
  #     )
  #     recoTau_min = tauReco.findAllTaus(
  #         particles_w_extremes_min, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut, charge_condition=charge_condition
  #     )
  # else: 
  #     recoTau_max = None
  #     recoTau_min = None
  tau_extremes = {"energy": {"max": None, "min": None}, "direction": {"max": None, "min": None}}
  if foton_config.get("energy", {}) and test_extremes:
      recoTau_max, recoTau_min = test_extremes_for_photons(particles, test_pfo, gatr_results_path, foton_config["energy"], syst_type="energy")
      tau_extremes["energy"]["max"] = recoTau_max
      tau_extremes["energy"]["min"] = recoTau_min
  if foton_config.get("direction", {}) and test_extremes:
      recoTau_max, recoTau_min = test_extremes_for_photons(particles, test_pfo, gatr_results_path, foton_config["direction"], syst_type="direction")
      tau_extremes["direction"]["max"] = recoTau_max
      tau_extremes["direction"]["min"] = recoTau_min
    
  if neutral_recover_cfg.get("return_hit_type_map", False):
      return recoTau, recoElectrons, recoMuons, tau_extremes, recover_extra_info, extra_info_dict
  return recoTau, recoElectrons, recoMuons, tau_extremes