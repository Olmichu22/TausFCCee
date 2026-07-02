import argparse
import copy
import logging
import math
import os
import pickle
import pprint
import sys
from array import array
from pathlib import Path
from modules.TauDecays import extractTauDecays
from modules.NeutralRecover import debug_reco_tau, plot_debug_reco_tau
from modules.ConfusionMatrixParticleLevel import plot_confusion_matrices, plot_energy_distributions
import edm4hep
import numpy as np
import pandas as pd
import ROOT
import yaml
from podio import root_io
from ROOT import TH1F, TH2F, TFile, TTree

from modules import (ParticleObjects, electronReco, muonReco, myutils, pi0Reco,
                     tauReco, particleMatch)
from modules.ParticleObjects import RecoParticle

def write_histograms_recursive(obj):
    """
    Recorre un diccionario anidado y ejecuta `.Write()` en cada objeto tipo ROOT histogram.
    """
    if isinstance(obj, dict):
        for value in obj.values():
            write_histograms_recursive(value)
    else:
        # Si no es diccionario, asumimos que es un histograma de ROOT
        try:
            obj.Write()
        except AttributeError:
            print(f"Objeto {obj} no tiene método .Write(). Ignorado.")

# ----------------------------------------------------------------------------
def my_hook(parser):
    parser.add_argument("--sys-err", type=str, default="config/systematics/err_sys.yml", help="YAML file with systematics errors to apply")
    parser.add_argument("--test-extremes", action="store_true", help="Test the extremes of the photon energy resolution")

# Load config (necessary for set up the logger)
default_config = "config/default/taurecolong.yaml"
# Output Configuration
outputbasepath = "Results/TauReco/"

general_configs = myutils.setup_analysis_config(default_config, outputbasepath, parser_hook = my_hook)


loggers = general_configs["loggers"]

run_config = general_configs["config"]
neutral_recover_cfg = run_config.get("neutral_recover", {})
# config = myutils.load_yaml_config(args.config, default_config)
logger_config = loggers["config"]
logger_io = loggers["io"]
logger_process = loggers["processing"]
logger_pi0mass = loggers["pi0mass"]

args =  general_configs["args"]
test_pfo = args.test_pfo

# Cut Configuration
dRMax=run_config["cuts"]["dRMax"]
minPTauPhoton =run_config["cuts"]["TauPhotonPCut"]
minPTauPion = run_config["cuts"]["TauPionPCut"]
PNeutron = run_config["cuts"]["NeutronCut"]
dRMatch = run_config["cuts"]["MatchedGenMinDR"]
generalPCut = run_config["cuts"]["generalPCut"]

selectDecay=general_configs["decay"]

outputpath = general_configs["outputpath"]
fileOutName = os.path.join(general_configs["outputpath"], general_configs["fileOutName"])



# Continue with the rest of configs
sys_errors = general_configs["config"].get("systematics_errors", {})
photon_config = sys_errors.get("photon_config", {})
test_extremes = args.test_extremes

# ------------------------------------------------------------------------
# General Configuration
sample=run_config["general"]["sample"]
matched_cm_arg = general_configs["flags"]["matched_cm"]
test_arg = general_configs["flags"]["test"]

logger_config.info("Configuration loaded!")
logger_config.info("Configuration:\n%s", pprint.pformat(general_configs, indent=4))


# ------------------------------------------------------------------------
gatr_results_path = general_configs["args"].gatr_result

filenames, mlpf_results = myutils.get_root_trees_path(sample,
                                                      gatr_results_path,
                                                      loggers,
                                                      test_arg,
                                                      args)
reader = root_io.Reader(filenames)
logger_io.info("Read %d files", len(filenames))
logger_io.info("First %s files.", filenames[:10]) 

# Configs and reading finished
# ----------------------------------------------------------------------


# collections to use
genparts = "MCParticles"
pfobjects = "PandoraPFOs"
# pfobjects ="TightSelectedPandoraPFOs"



histogram_config = general_configs.get("histograms_config", {})
root_histograms = myutils.set_up_root_histograms(histogram_config)


if test_extremes:
    logger_process.info("Testing extremes is enabled.")

    root_histograms_super = {
        "original": root_histograms,
        "min_err": myutils.clone_histograms_with_suffix(root_histograms, "_min"),
        "max_err": myutils.clone_histograms_with_suffix(root_histograms, "_max")
    }

else:
    root_histograms_super = {"original": root_histograms}
    

if test_extremes:
    logger_process.info("Testing extremes is enabled.")

    root_histograms_super = {
        "original": root_histograms,
        "min_err": myutils.clone_histograms_with_suffix(root_histograms, "_min"),
        "max_err": myutils.clone_histograms_with_suffix(root_histograms, "_max")
    }

else:
    root_histograms_super = {"original": root_histograms}
    
result_labels = {}
result_labels["tau1"] = []
result_labels["tau2"] = []
result_labels["id-tau1"] = []
result_labels["id-tau2"] = []

if test_extremes:
    results_labels_super = {
        "original": result_labels,
        "min_err": copy.deepcopy(result_labels),
        "max_err": copy.deepcopy(result_labels)
    }
else:
    results_labels_super = {"original": result_labels}

true_predicted_label = {"GenID": [], "True": [], "Predicted": [], "PhotonPredicted": [],
                        "Countpions":[],"Countphotons": [], "Countneutrons": []}
# unmatched_true_label = {}

if test_extremes:
    true_predicted_label_super = {
        "original": true_predicted_label,
        "min_err": copy.deepcopy(true_predicted_label),
        "max_err": copy.deepcopy(true_predicted_label)
    }
else:
    true_predicted_label_super = {"original": true_predicted_label}


countEvents = 0
# unmatched_reco_pions_match_list = []
# unmatched_reco_pions_P_per_miss = {"Non_matched":root_histograms["Reco"]["Events"]["UnmatchedPionsP"]}
# unmatched_reco_pions_theta_per_miss = {}

# unmatched_gen_pions_match_list = []
# unmatched_gen_pions_P_per_miss = dict()
# unmatched_gen_pions_theta_per_miss = dict()
association_results_df = dict()
energy_distribution_results = dict()
bins = [0, 1, 5, 10, 20, 30, 45, 100, np.inf]
energy_max = 60
energy_min = 0 
# 50 bins
e_bins = np.linspace(energy_min, energy_max, 51)

association_df_cumulated = dict()
for eventid, event in enumerate(reader.get("events")):
    # if gatr_results_path is not None and eventid > len(gatr_results) - 1:
    #     logger_process.info("Reached the end of GATr results, stopping processing.")
    #     break
    logger_process.debug("Processing event %d", eventid)
    if countEvents % 1000 == 0:
        logger_process.info("Processing event %d", countEvents)
    countEvents += 1

    mc_particles = event.get(genparts)
    pfos = event.get(pfobjects)

    genTaus = tauReco.findAllGenTaus(mc_particles)
    nGenTaus = len(genTaus)

    logger_process.debug(
        "Found %d gen taus. Details:\n%s",
        nGenTaus,
        "\n".join("GenTau %d: %s" % (i, tau) for i, tau in genTaus.items()),
    )
    if neutral_recover_cfg.get("return_hit_type_map", False):
        recoTau, recoElectrons, recoMuons, tau_extremes, recover_extra_info, extra_info_dict = extractTauDecays(gatr_results_path,
                                                                                    mlpf_results,
                                                                                    eventid,
                                                                                    pfos,
                                                                                    dRMax,
                                                                                    minPTauPhoton,
                                                                                    minPTauPion,
                                                                                    PNeutron,
                                                                                    generalPCut,
                                                                                    photon_config,
                                                                                    test_extremes,
                                                                                    test_pfo,
                                                                                    logger_process,
                                                                                    neutral_recover_cfg=neutral_recover_cfg,
                                                                                    event=event,
                                                                                    only_association=True)
    else:
        recoTau, recoElectrons, recoMuons, tau_extremes = extractTauDecays(gatr_results_path,
                                                                                    mlpf_results,
                                                                                    eventid,
                                                                                    pfos,
                                                                                    dRMax,
                                                                                    minPTauPhoton,
                                                                                    minPTauPion,
                                                                                    PNeutron,
                                                                                    generalPCut,
                                                                                    photon_config,
                                                                                    test_extremes,
                                                                                    test_pfo,
                                                                                    logger_process,
                                                                                    neutral_recover_cfg=neutral_recover_cfg,
                                                                                    event=event)

    if neutral_recover_cfg.get("return_hit_type_map", False):
        hit_type_map = extra_info_dict.get("hit_type_map", {})
        non_assoc_tracks_df = extra_info_dict.get("non_associated_tracks_df", pd.DataFrame())
        df_reco_mc_links = extra_info_dict.get("df_reco_mc_links", pd.DataFrame())
    else:
        hit_type_map = None
        non_assoc_tracks_df = pd.DataFrame()
        df_reco_mc_links = pd.DataFrame()
    
    # tau_extremes is a dict {"energy": {...}, "direction": {...}} populated only
    # when test_extremes is enabled; the per-extreme reco taus are consumed solely
    # by the (currently disabled) test-extremes block below.
    recoTaus_extremes = {
        "original": recoTau,
        "min_err": None,
        "max_err": None
    }
    # print(df_reco_mc_links.head())
    for row in df_reco_mc_links.itertuples():
        gen_pid = row.Gen_pid
        reco_pid = row.Reco_pid
        base_key = str(abs(gen_pid)) + "_" + str(abs(reco_pid))
        # print(key)
        gen_energy = row.Gen_energy
        reco_energy = row.Reco_energy
        energy_bin = pd.cut([gen_energy], bins=bins)[0]
        energy_bin_for_dist = pd.cut([gen_energy], bins=e_bins)[0]
        key = base_key + "_" + str(energy_bin)
        key_dist = base_key + "_" + str(energy_bin_for_dist)
        if key not in association_results_df:
            association_results_df[key] = 0
        association_results_df[key] += 1
        if not np.isfinite(gen_energy) or not np.isfinite(reco_energy) or gen_energy == 0:
            continue
        if key_dist not in energy_distribution_results:
            energy_distribution_results[key_dist] = []
        energy_distribution_results[key_dist].append((gen_energy, reco_energy))
    # exit(0)
    association_df_cumulated[eventid] = df_reco_mc_links
    # Photon resolution
    # gen_comp = particleMatch.GetGenTauDecayProducts(mc_particles)
                
    # for key in root_histograms_super.keys():
        
        
    #     # # Select the type of reconstructed taus to use
    #     # logger_process.debug("Processing extremes type: %s", key)
    #     # recoTau = recoTaus_extremes[key]
    #     root_histograms = root_histograms_super[key]
    #     root_histograms["Gen"][""]
        # true_predicted_label = true_predicted_label_super[key]
        # result_labels = results_labels_super[key]
        
        # nRecoTaus = len(recoTau)
        # nRecoElectrons = len(recoElectrons)
        # nRecoMuons = len(recoMuons)

        # recoTaus = {}
        # pidx = 0
        # for taui in range(nRecoTaus):
        #     recoTaus[pidx] = recoTau[taui]
        #     pidx += 1
        # for elei in range(nRecoElectrons):
        #     recoTaus[pidx] = recoElectrons[elei]
        #     pidx += 1
        # for mui in range(nRecoMuons):
        #     recoTaus[pidx] = recoMuons[mui]
        #     pidx += 1
        # nRecoTaus = len(recoTaus)

        # logger_process.debug(
        #     "Found %d reconstructed taus. Details:\n%s",
        #     nRecoTaus,
        #     "\n".join("RecoTau %d: %s" % (i, tau) for i, tau in recoTaus.items()),
        # )

        # for tau_idx, reco_tau in recoTaus.items():
        #   decay_id = reco_tau.getID()
        #   if decay_id == -1:
        #     neutral_hit_collection = debug_reco_tau(reco_tau, hit_type_map)
        #     # Generate plots
        #     out_dir = os.path.join(outputpath, f"event_{eventid}_tau_{tau_idx}_decay_{decay_id}")
        #     os.makedirs(out_dir, exist_ok=True)
        #     plot_debug_reco_tau(neutral_hit_collection, non_assoc_tracks_df=non_assoc_tracks_df, output_dir=out_dir, event_idx=eventid)
        #     # for const_id, const in reco_tau.getDaughters().items():
        #     #   if abs(const.getPDG()) == 2112:
        #     #     print("Found tau with ID -1 and a neutron among its daughters:")
        #     #     print("Constituents:")
        #     #     for S_const_id, S_const in reco_tau.getDaughters().items():
        #     #       print(S_const.getPDG())
        #     exit(0)


logger_io.info("Processed %d events", countEvents)

# outfile = ROOT.TFile(fileOutName, "RECREATE")
# for key in root_histograms_super:
#     root_histograms = root_histograms_super[key]
#     true_predicted_label = true_predicted_label_super[key]
#     result_labels = results_labels_super[key]
#     suffix = "" if key == "original" else f"_{key}"
#     myutils.write_plot_config(root_histograms, outputpath, suffix)
# # exit(0)
#     root_histograms = myutils.calc_efficiency(root_histograms, histogram_config, suffix)
#     write_histograms_recursive(root_histograms)




#     decaystr = "decayAll" if selectDecay == -777 else "decay{}".format(selectDecay)
#     true_predicted_label_output_file = outputpath + f"true_predicted_label_{decaystr}{suffix}.csv"
#     true_predicted_label_df = pd.DataFrame(true_predicted_label)
#     true_predicted_label_df.to_csv(true_predicted_label_output_file, index=False)

    # result_labels = pd.DataFrame(result_labels)
    # result_labels.to_csv(outputpath + f"result_labels_pfo{suffix}.csv", index=False)

# unmatched_reco_pions_match_list_df = pd.DataFrame({"PDGID": unmatched_reco_pions_match_list})
# unmatched_reco_pions_match_count = unmatched_reco_pions_match_list_df["PDGID"].value_counts()
# Guardamos
# unmatched_reco_pions_match_count.to_csv(outputpath + "unmatched_reco_pions_match_count.csv", index=True) 

# unmatched_gen_pions_match_list_df = pd.DataFrame({"PDGID": unmatched_gen_pions_match_list})
# unmatched_gen_pions_match_count = unmatched_gen_pions_match_list_df["PDGID"].value_counts()
# Guardamos
# unmatched_gen_pions_match_count.to_csv(outputpath + "unmatched_gen_pions_match_count.csv", index=True)
# Check if config["output"]["outputlabels"] is a list

# Associations whole df
# add event id column
for eventid, df in association_df_cumulated.items():
    df["event_id"] = eventid
association_results_full_df = pd.concat(association_df_cumulated.values(), ignore_index=True)
association_results_full_df.to_csv(outputpath + "association_results_full.csv", index=False)


images_dir = os.path.join(outputpath, "confusion_matrices_particle_level")
plot_confusion_matrices(association_results_df, output_dir=images_dir)

energy_distribution_results_dir = os.path.join(outputpath, "energy_distributions")
plot_energy_distributions(energy_distribution_results, output_dir=energy_distribution_results_dir)
if type(run_config["output"]["outputlabels"]) is not list:
    if run_config["output"]["outputlabels"] is None:
        run_config["output"]["outputlabels"] = []
    else:
        run_config["output"]["outputlabels"] = [run_config["output"]["outputlabels"]]
# if true_predicted_label_output_file not in run_config["output"]["outputlabels"]:
    # run_config["output"]["outputlabels"].append(outputpath + f"true_predicted_label_{decaystr}.csv")

output_config_file = outputpath + "config.yaml"
with open(output_config_file, "w") as file:
    yaml.dump(run_config, file)
    logger_io.info("Configuration file saved to %s", output_config_file)



# outfile = ROOT.TFile(fileOutName, "RECREATE")
# for key in root_histograms_super:
#     root_histograms = root_histograms_super[key]
#     true_predicted_label = true_predicted_label_super[key]
#     result_labels = results_labels_super[key]
#     suffix = "" if key == "original" else f"_{key}"
#     myutils.write_plot_config(root_histograms, outputpath, suffix)
# # exit(0)
#     root_histograms = myutils.calc_efficiency(root_histograms, histogram_config, suffix)
#     write_histograms_recursive(root_histograms)

logger_io.info("Output file %s", fileOutName)
logger_io.info("End of job")
# outfile.Close()
