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


# Modos de desintegración seguidos en los plots de migración, con la etiqueta
# usada en el nombre del histograma. Los dos modos leptónicos llevan palabra en
# vez del id con signo para no meter un menos en los nombres de ROOT.
MIGRATION_MODE_TAG = {0: "0", 1: "1", 2: "2", 3: "3", 4: "4",
                      10: "10", 11: "11", 12: "12", -11: "Elec", -13: "Muon"}

# Modos gen que ya tienen su propio bloque de llenado en el bucle de eventos.
MIGRATION_MODES_ALREADY_FILLED = (0, 1, 2, 10)


# A partir de este numero de fotones se agrupan todos juntos: entre todos no
# llegan al 2% de la muestra.
MIGRATION_MAX_PHOTONS = 5


def migration_label(recoTauId):
    """
    Map a raw reconstructed tau id to its migration-histogram label.

    The reco side of the migration plots counts PHOTONS, not pi0s, so this uses
    ``recoTauId`` rather than the ``recoDM`` that feeds the confusion matrix:
    one and two photons land in "h1g" and "h2g" instead of being merged.

    Ids are ``nPhotons`` for one-prong taus and ``10 + nPhotons`` for three-prong
    ones. Everything from MIGRATION_MAX_PHOTONS photons up is pooled into a
    single bucket per prong count. ``-21`` (three prongs plus a neutral hadron
    absorbed into the tau) gets its own "3hNeutron" bucket and anything else
    ends up in "Other"; the "Unmatched" and "NoTau" labels are filled directly
    at their call sites.
    """
    if recoTauId == -11:
        return "Elec"
    if recoTauId == -13:
        return "Muon"
    if recoTauId == -21:
        return "3hNeutron"
    if recoTauId < 0:
        return "Other"
    prong = "h" if recoTauId < 10 else "3h"
    n_photons = recoTauId if recoTauId < 10 else recoTauId - 10
    if n_photons == 0:
        return prong
    if n_photons >= MIGRATION_MAX_PHOTONS:
        return f"{prong}{MIGRATION_MAX_PHOTONS}gPlus"
    return f"{prong}{n_photons}g"


def fill_decay_migration(root_histograms, genTauId, label, genTauP, genTauTheta):
    """
    Fill the "where does this gen decay mode end up" histograms for one gen tau.

    Both the momentum and the polar angle version are filled, so the migration
    can be read against either variable.

    Args:
        root_histograms (dict): Nested dict of ROOT histograms.
        genTauId (int): Generator level decay mode of the tau.
        label (str): Reco outcome, i.e. a decay mode or "Unmatched"/"NoTau"/"Other".
        genTauP (float): Momentum of the generator level tau.
        genTauTheta (float): Polar angle of the generator level tau.
    """
    gen_tag = MIGRATION_MODE_TAG.get(genTauId)
    if gen_tag is None:
        return
    for variable, value in (("TauP", genTauP), ("TauTheta", genTauTheta)):
        hist = root_histograms["Matched"]["Events"].get(f"{variable}{gen_tag}To{label}")
        if hist is not None:
            hist.Fill(value)


def fill_migration_denominators(root_histograms, genTauId, genTauP, genTauTheta):
    """
    Fill the Gen denominators of the migration plots for one gen tau.

    Only the modes without a dedicated block in the event loop are handled here;
    the rest are already filled alongside their other gen distributions.
    """
    gen_tag = MIGRATION_MODE_TAG.get(genTauId)
    if gen_tag is None or genTauId in MIGRATION_MODES_ALREADY_FILLED:
        return
    root_histograms["Gen"]["Events"][f"TauP{gen_tag}"].Fill(genTauP)
    root_histograms["Gen"]["Events"][f"TauTheta{gen_tag}"].Fill(genTauTheta)

# ----------------------------------------------------------------------------
def my_hook(parser):
    parser.add_argument("--sys-err", type=str, default="config/systematics/err_sys.yml", help="YAML file with systematics errors to apply")
    parser.add_argument("--test-extremes", action="store_true", help="Test the extremes of the photon energy resolution")


# filenames = ['/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_0.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_1.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_10.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_11.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_12.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_13.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_14.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_15.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_16.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_17.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_18.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_19.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_2.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_20.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_21.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_22.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_23.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_24.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_26.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_27.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_28.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_29.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_3.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_30.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_31.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_32.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_33.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_35.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_36.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_37.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_38.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_39.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_4.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_40.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_41.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_44.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_45.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_46.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_47.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_48.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_49.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_5.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_50.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_52.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_53.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_54.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_55.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_56.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_58.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_59.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_6.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_60.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_61.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_62.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_64.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_65.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_66.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_67.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_68.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_69.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_7.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_70.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_72.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_73.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_74.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_75.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_76.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_77.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_78.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_79.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_8.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_80.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_81.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_82.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_83.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_85.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_86.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_87.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_88.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_89.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_9.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_90.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_91.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_92.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_93.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_94.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_95.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_96.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_97.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_98.root', '/pnfs/ciemat.es/data/calice/local/arqolmo/ILDFCCFullSim/Ztautau/out_reco_edm4hep_99.root']
# reader = root_io.Reader(filenames)
# for i, event in enumerate(reader.get("events")):
    # print(i)
# exit(0)
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
dRMatch = run_config["cuts"]["MatchedGenMaxDR"]
generalPCut = run_config["cuts"]["generalPCut"]

# Correcciones extra sobre el tau reconstruido (tauReco.extraTauRecoCorrection).
extra_correction = run_config.get("extra_reco_correction") or {}
if extra_correction.get("enable") and extra_correction.get("modes"):
    unknown = [m for m in extra_correction["modes"]
               if m not in tauReco.availableExtraCorrections()]
    if unknown:
        logger_config.error("Modo(s) de correccion desconocido(s): %s. Disponibles: %s",
                            unknown, tauReco.availableExtraCorrections())
        sys.exit(1)
    logger_config.info("Correcciones extra activas: %s (params: %s)",
                       extra_correction["modes"], extra_correction.get("params", {}))
else:
    extra_correction = None

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
                                                      args,
                                                      skip_root_validation=False)
# print(filenames)
# for file in filenames:
    # check_existence = Path(file).is_file()
    # print(f"Checking file: {file}, exists: {check_existence}")
reader = root_io.Reader(filenames)
# for i, event in enumerate(reader.get("events")):
    # print("Event %d: %s", i, event)
# print("Fuera del bucle")
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
                                                                                   event=event,
                                                                                   extra_correction=extra_correction)

    recoTaus_extremes = {
        "original": recoTau,
        # "min_err": recoTau_min,
        # "max_err": recoTau_max
    }
    # Photon resolution
    gen_comp = particleMatch.GetGenTauDecayProducts(mc_particles)
    gen_photons = {k: v for k, v in gen_comp.items() if abs(v.getPDG()) == 22}
    for reco_particle in pfos:
        pdg = abs(reco_particle.getPDG())
        if pdg == 22:
            if isinstance(reco_particle.getMomentum(), ROOT.TLorentzVector):            
                recoPhotonP4 = reco_particle.getMomentum()
            else:
                recoPhotonP4 = ROOT.TLorentzVector()
                recoPhotonP4.SetXYZM(
                    reco_particle.getMomentum().x,
                    reco_particle.getMomentum().y,
                    reco_particle.getMomentum().z,
                    0
                )
            min_DR = 9999
            max_dr = 0.4
            matched_gen_photon = -1
            for gen_idx, gen_photon in gen_photons.items():
                genPhotonP4 = gen_photon.getMomentum()
                dR = myutils.dRAngle(recoPhotonP4, genPhotonP4)
                if dR < min_DR and dR < max_dr:
                    min_DR = dR
                    matched_gen_photon = gen_idx
            if matched_gen_photon >= 0:
                genPhotonP4 = gen_photons[matched_gen_photon].getMomentum()
                photonPRes = (recoPhotonP4.P() - genPhotonP4.P()) / (genPhotonP4.P() + 1e-8)
                
                root_histograms = root_histograms_super["original"]
                root_histograms["Matched"]["Resolution"]["PhotonPRes"].Fill(photonPRes)
            
    # print(root_histograms["Matched"]["Resolution"]["PhotonPRes"].GetEntries())        
    # Neutral hadrons (Not interested in this way)
    # if gatr_results_path is not None and not general_configs["args"].test_pfo:
    #     particles = mlpf_results.get("particles", {})
    # else:
    #     particles = pfos

    
    # for particle in particles:
    #     recoParticlecharge = particle.getCharge()
    #     if recoParticlecharge == 0:
    #         recoParticleP4 = particle.getMomentum()
    #         if isinstance(recoParticleP4, ROOT.TLorentzVector):
    #             root_histograms["Reco"]["Events"]["NeutralHadronP"].Fill(recoParticleP4.P())
    #         else:
    #             recoParticleP4 = ROOT.TLorentzVector()
    #             recoParticleP4.SetXYZM(
    #                 particle.getMomentum().x,
    #                 particle.getMomentum().y,
    #                 particle.getMomentum().z,
    #                 particle.getMass()
    #             )
    #             root_histograms["Reco"]["Events"]["NeutralHadronP"].Fill(recoParticleP4.P())
                
    for key in root_histograms_super.keys():
        # Select the type of reconstructed taus to use
        logger_process.debug("Processing extremes type: %s", key)
        recoTau = recoTaus_extremes[key]
        root_histograms = root_histograms_super[key]
        true_predicted_label = true_predicted_label_super[key]
        result_labels = results_labels_super[key]
        
        nRecoTaus = len(recoTau)
        nRecoElectrons = len(recoElectrons)
        nRecoMuons = len(recoMuons)

        recoTaus = {}
        pidx = 0
        for taui in range(nRecoTaus):
            recoTaus[pidx] = recoTau[taui]
            pidx += 1
        for elei in range(nRecoElectrons):
            recoTaus[pidx] = recoElectrons[elei]
            pidx += 1
        for mui in range(nRecoMuons):
            recoTaus[pidx] = recoMuons[mui]
            pidx += 1
        nRecoTaus = len(recoTaus)

        logger_process.debug(
            "Found %d reconstructed taus. Details:\n%s",
            nRecoTaus,
            "\n".join("RecoTau %d: %s" % (i, tau) for i, tau in recoTaus.items()),
        )

        foundGen = False

        nGenTausType = 0
        nTausType = 0
        nGenTausHad = 0
        
        if nGenTaus == 0:
            result_labels["tau1"].append(-999)
            result_labels["tau2"].append(-999)
            result_labels["id-tau1"].append(-999)
            result_labels["id-tau2"].append(-999)
        
        GenTausP4 = []
        GenTausVisP4 = []
        RecoTausP4 = []
        for i in range(0, nGenTaus):
            genVisTauP4 = genTaus[
                i
            ].getvisMomentum()  
            genTauId = genTaus[i].getID()
            genTauQ = genTaus[i].getCharge()
            genTauP4 = genTaus[i].getMomentum()
            genTauDR = genTaus[
                i
            ].getMaxAngle()  # Maximum angle between the tau and its constituents
            genTauNConsts = genTaus[i].getnConst()
            genTauConsts = genTaus[i].getDaughters()
            
            GenTausP4.append(genTauP4)
            GenTausVisP4.append(genVisTauP4)
            
            if nGenTaus > 2 and i <=1:
                result_labels[f"tau{i+1}"].append(-999)
            elif nGenTaus <= 2:
                result_labels[f"tau{i+1}"].append(genTauId)
            if nGenTaus == 1:
                result_labels["tau2"].append(-999)
            
            # remove leptonic decays
            if genTauId >= 0:
                nGenTausHad += 1

            # if genTauId < 0:
            #    continue

            # nGenTausHad+=1

            # pick only a decay mode in particular if you want
            if selectDecay != -777 and selectDecay != genTauId:
                continue

            nGenTausType += 1
            foundGen = True
            
            # Get Photon Information at generation level
            cum_p4 = ROOT.TLorentzVector()
            cum_p4.SetXYZM(
                0, 0, 0, 0
            )
            for c in range(0, genTauNConsts):
                const = genTauConsts[c]
                dauP4 = ROOT.TLorentzVector()
                try:
                    dauP4.SetXYZM(
                        const.getMomentum().x,
                        const.getMomentum().y,
                        const.getMomentum().z,
                        const.getMass(),
                    )
                except AttributeError:
                    # If the dau does not have a mass, we assume it is a photon
                    dauP4.SetXYZM(
                        const.getMomentum().X(),
                        const.getMomentum().Y(),
                        const.getMomentum().Z(),
                        const.getMass(),
                    )
                if genTauId == 0:
                    cum_p4 += dauP4
                # Si se trata de un pion neutro, lo añadimos a la suma de fotones
                if const.getPDG() == 111:
                    photons = const.getDaughters()
                    for photon in photons:
                        photonP4 = ROOT.TLorentzVector()
                        photonP4.SetXYZM(
                            photon.getMomentum().x,
                            photon.getMomentum().y,
                            photon.getMomentum().z,
                            photon.getMass(),
                        )
                        root_histograms["Gen"]["Events"]["PhotonP"].Fill(photonP4.P())
                        root_histograms["Gen"]["Events"]["PhotonTheta"].Fill(photonP4.Theta())
                elif abs(const.getPDG()) == 211:  # Pion
                    root_histograms["Gen"]["Events"]["PionPid"].Fill(dauP4.P(), genTauId)
                    root_histograms["Gen"]["Events"]["PionP"].Fill(dauP4.P())
                    if genTauId == 0:
                        root_histograms["Gen"]["Events"]["PionP0"].Fill(dauP4.P())
                    elif genTauId == 1:
                        root_histograms["Gen"]["Events"]["PionP1"].Fill(dauP4.P())
                    elif genTauId == 2:
                        root_histograms["Gen"]["Events"]["PionP2"].Fill(dauP4.P())
                    elif genTauId == 10:
                        root_histograms["Gen"]["Events"]["PionP10"].Fill(dauP4.P())
                    if dauP4.P() > 5:
                        root_histograms["Gen"]["Events"]["PionTheta"].Fill(dauP4.Theta())
                    # If the photon is matched, fill the histogram

            # # # P4 Tau filters
            # if genVisTauP4.P() < 5:
            #     continue
            
            ###############################################################
            # IMPORTANT!!! CUT IN 0.95
            if abs(math.cos(genVisTauP4.Theta())) > 0.95:
                continue
            ################################################################
            
            
            
            # print ("Gen",genTauP4.P(),genVisTauP4.P(),genVisTauP4.Theta(),genVisTauP4.Phi(),genTauId,genTauQ,genTauDR,genTauNConsts)

            # Fill histograms
            root_histograms["Gen"]["Events"]["TauPt"].Fill(genTauP4.Pt())  # Transverse momentum
            root_histograms["Gen"]["Events"]["TauVisPt"].Fill(genVisTauP4.Pt())  # Visible transverse momentum
            root_histograms["Gen"]["Events"]["TauP"].Fill(genTauP4.P())  # Momentum
            root_histograms["Gen"]["Events"]["TauVisP"].Fill(genVisTauP4.P())  # Visible momentum
            root_histograms["Gen"]["Events"]["TauVisMass"].Fill(genVisTauP4.M())  # Visible mass
            root_histograms["Gen"]["Events"]["TauType"].Fill(genTauId)  # Tau decay type
            root_histograms["Gen"]["Events"]["TauQ"].Fill(genTauQ)  # Tau charge
            root_histograms["Gen"]["Events"]["TauEta"].Fill(genTauP4.Eta())  # Pseudo-rapidity
            root_histograms["Gen"]["Events"]["TauTheta"].Fill(genTauP4.Theta())  # Theta angle

            if genTauId >= 0:
                root_histograms["Gen"]["Events"]["HadronTauP"].Fill(genTauP4.P())  # Hadronic tau momentum
            if genTauId == 0:
                root_histograms["Gen"]["Events"]["TauTheta0"].Fill(genTauP4.Theta())
                root_histograms["Gen"]["Events"]["TauP0"].Fill(genTauP4.P())
                root_histograms["Gen"]["Events"]["TauVisP0"].Fill(genVisTauP4.P())
            elif genTauId == 1:
                root_histograms["Gen"]["Events"]["TauTheta1"].Fill(genTauP4.Theta())
                root_histograms["Gen"]["Events"]["TauP1"].Fill(genTauP4.P())
                root_histograms["Gen"]["Events"]["TauVisP1"].Fill(genVisTauP4.P())
            elif genTauId == 2:
                root_histograms["Gen"]["Events"]["TauTheta2"].Fill(genTauP4.Theta())
                root_histograms["Gen"]["Events"]["TauP2"].Fill(genTauP4.P())
                root_histograms["Gen"]["Events"]["TauVisP2"].Fill(genVisTauP4.P())
            elif genTauId == 10:
                root_histograms["Gen"]["Events"]["TauTheta10"].Fill(genTauP4.Theta())
                root_histograms["Gen"]["Events"]["TauP10"].Fill(genTauP4.P())
                root_histograms["Gen"]["Events"]["TauVisP10"].Fill(genVisTauP4.P())

            # Denominadores de los plots de migracion para los modos sin bloque propio
            fill_migration_denominators(
                root_histograms, genTauId, genTauP4.P(), genTauP4.Theta()
            )

            root_histograms["Gen"]["Events"]["TauDR"].Fill(genTauDR)  # Angle of Tau Constituents
            countPionsRun = 0

            # print ("all GEN")
            # Look inside the generator level tau: check the constituents (decay products)

            # Compare with reconstructed taus using angle matching
            findMatch, nTausType = tauReco.MatchRecoGenTau(
                genTaus[i], recoTaus, nTausType, maxDRMatch=dRMatch, selectDecay=selectDecay
            )
            # For each generator level tau, find the reconstructed tau that is closest:
            # if not matched_cm:
            true_predicted_label["GenID"].append(str(eventid) + str(i))
            true_predicted_label["True"].append(genTauId)

            # If you have not found it, continue: this is a efficiency loss
            if findMatch == -1:
                logger_process.debug("No match found for gen tau %s", genTaus[i])

                if nGenTaus > 2 and i <=1:
                    result_labels[f"id-tau{i+1}"].append(-999)
                elif nGenTaus <= 2:
                    result_labels[f"id-tau{i+1}"].append(-2)
                
                if nGenTaus == 1:
                    result_labels["id-tau2"].append(-999)
                
                # if not matched_cm:
                    # true_predicted_label["Predicted"].append(-1)
                true_predicted_label["Predicted"].append(-2)
                true_predicted_label["PhotonPredicted"].append(-2)
                true_predicted_label["Countpions"].append(-999)
                true_predicted_label["Countphotons"].append(-999)
                true_predicted_label["Countneutrons"].append(-999)
                fill_decay_migration(
                    root_histograms, genTauId, "Unmatched", genTauP4.P(), genTauP4.Theta()
                )
                continue
            

            logger_process.debug("Found matched tau. Details:\n%s", recoTaus[findMatch])
            # if matched_cm:
            #     true_predicted_label["GenID"].append(str(eventid) + str(i))
            #     true_predicted_label["True"].append(genTauId)
            # now, get the kinematics of the matched reco tau
            recoTauP4 = recoTaus[findMatch].getMomentum()
            recoTauId = recoTaus[findMatch].getID()
            recoTauQ = recoTaus[findMatch].getCharge()
            recoTauDR = recoTaus[findMatch].getMaxCone()
            recoTauConsts = recoTaus[findMatch].getDaughters()
            recoTauNConsts = recoTaus[findMatch].getnConst()

            if recoTauId == -1:
                logger_process.debug("Match is not tau for gen tau %s", genTaus[i])

                if nGenTaus > 2 and i <=1:
                    result_labels[f"id-tau{i+1}"].append(-1)
                elif nGenTaus <= 2:
                    result_labels[f"id-tau{i+1}"].append(-1)
                
                if nGenTaus == 1:
                    result_labels["id-tau2"].append(-1)
                
                # if not matched_cm:
                    # true_predicted_label["Predicted"].append(-1)
                true_predicted_label["Predicted"].append(-1)
                true_predicted_label["PhotonPredicted"].append(-1)
                n_photons = 0
                n_pions = 0
                n_neutrons = 0
                for comp in range(0, recoTauNConsts):
                    const = recoTauConsts[comp]
                    if abs(const.getPDG()) == 211:
                        n_pions += 1
                    elif const.getPDG() == 22:
                        n_photons += 1
                    elif const.getPDG() == 2112:
                        n_neutrons += 1

                true_predicted_label["Countpions"].append(n_pions)
                true_predicted_label["Countphotons"].append(n_photons)
                true_predicted_label["Countneutrons"].append(n_neutrons)
                fill_decay_migration(
                    root_histograms, genTauId, "NoTau", genTauP4.P(), genTauP4.Theta()
                )
                continue

            RecoTausP4.append(recoTauP4)
            # Reassign the recoTauId to the recoDM
            # This need to be checked as there exist other ID at Gen Level
            recoDM = recoTauId
            n_pi0s = 0
            
            if recoTauId < 10 and recoTauId >= 0:
                nPhotons = recoTauId
                n_pi0s = math.ceil(nPhotons / 2)
                recoDM = n_pi0s
            elif recoTauId >= 10:
                nPhotons = recoTauId - 10
                n_pi0s = math.ceil(nPhotons / 2)
                recoDM = 10 + n_pi0s
            
            if nGenTaus > 2 and i <=1:
                result_labels[f"id-tau{i+1}"].append(-999)
            elif nGenTaus <= 2:
                result_labels[f"id-tau{i+1}"].append(recoDM)

            
            true_predicted_label["Predicted"].append(recoDM)
            true_predicted_label["PhotonPredicted"].append(recoTauId)

            fill_decay_migration(
                root_histograms, genTauId, migration_label(recoTauId),
                genTauP4.P(), genTauP4.Theta()
            )

            root_histograms["Matched"]["Resolution"]["TauP"].Fill((recoTauP4.P() - genTauP4.P()) / genTauP4.P())
            root_histograms["Matched"]["Resolution"]["TauPt"].Fill((recoTauP4.Pt() - genTauP4.Pt()) / genTauP4.Pt())
            # hMatchedTausChargeRes.Fill(abs(recoTauQ) - abs(genTauQ) / abs(genTauQ))
            root_histograms["Matched"]["Resolution"]["TauMaxAngle"].Fill((recoTauDR - genTauDR) / genTauDR)
            root_histograms["Matched"]["Resolution"]["TauNComp"].Fill((recoTauNConsts - genTauNConsts) / genTauNConsts)

            #          print ("Reco?",recoTauP4.P(),recoTauId,recoTauQ,recoTauDR,recoTauNConsts)

            # Now that we have a matched (gen,reco) pair, more checks for efficiency and resolution

            countPionsRun = 0
            # print ("Matched GEN!")
            # GEN: Look inside the tau, constituents:
            for c in range(0, genTauNConsts):
                const = genTauConsts[c]
                
                dauP4 = ROOT.TLorentzVector()
                try:
                    dauP4.SetXYZM(
                        const.getMomentum().x,
                        const.getMomentum().y,
                        const.getMomentum().z,
                        const.getMass(),
                    )
                except AttributeError:
                    # If the dau does not have a mass, we assume it is a photon
                    dauP4.SetXYZM(
                        const.getMomentum().X(),
                        const.getMomentum().Y(),
                        const.getMomentum().Z(),
                        const.getMass(),
                    )
                # Si se trata de un pion neutro, lo añadimos a la suma de fotones
                if const.getPDG() == 111:
                    photons = const.getDaughters()
                    for photon in photons:
                        photonP4 = ROOT.TLorentzVector()
                        photonP4.SetXYZM(
                            photon.getMomentum().x,
                            photon.getMomentum().y,
                            photon.getMomentum().z,
                            photon.getMass(),
                        )
                        root_histograms["Matched"]["Events"]["PhotonP"].Fill(photonP4.P())
                        root_histograms["Matched"]["Events"]["PhotonTheta"].Fill(photonP4.Theta())
                elif abs(const.getPDG()) == 211:  # Pion
                    root_histograms["Matched"]["Events"]["PionP"].Fill(dauP4.P())
                    if dauP4.P() > 5:
                        root_histograms["Matched"]["Events"]["PionTheta"].Fill(dauP4.Theta())
                    # If the photon is matched, fill the histogram


            countPionsRun = 0
            # Init empyu TLorentzVector for the photon momentum
            photonCumulativeP4 = ROOT.TLorentzVector()
            photonCumulativeP4.SetXYZM(0, 0, 0, 0)
            # RECO:  Look inside the tau, constituents:
            # Filling histograms for the matched tau with reco level information
            recoTaus_photons = {}
            n_photons = 0
            n_pions = 0
            n_neutrons = 0
            for c in range(0, recoTauNConsts):
                const = recoTauConsts[c]
                constP4 = ROOT.TLorentzVector()
                try:
                    constP4.SetXYZM(
                        const.getMomentum().x,
                        const.getMomentum().y,
                        const.getMomentum().z,
                        const.getMass(),
                    )
                except AttributeError:
                    # If the const does not have a mass, we assume it is a photon
                    constP4.SetXYZM(
                        const.getMomentum().X(),
                        const.getMomentum().Y(),
                        const.getMomentum().Z(),
                        const.getMass(),
                    )

                # This may be an error as we are confusing photons with neutrons

                logger_pi0mass.debug(f"Evaluating const to get photon with PDGID {const.getPDG()} and charge {const.getCharge()}")
                if const.getCharge() == 0 and const.getPDG() != 22:
                    root_histograms["Reco"]["Events"]["NeutralHadronTauP"].Fill(constP4.P())
                if const.getPDG() == 22:
                    photonCumulativeP4 += constP4
                    recoTaus_photons[n_photons] = const
                    n_photons += 1
                    root_histograms["Reco"]["Events"]["PhotonP"].Fill(constP4.P())
                    root_histograms["Reco"]["Events"]["PhotonTheta"].Fill(constP4.Theta())
                elif abs(const.getPDG()) == 211:  # Pion
                    root_histograms["Reco"]["Events"]["PionP"].Fill(constP4.P())
                    if constP4.P() < minPTauPion:
                        logger_process.warning(f"Encontrado pion con momento menor {constP4.P()}")
                        logger_process.warning(const)
                    # if constP4.P() < 15:
                    #     logger_process.warning(f"Encontrado pion con momento mayor")
                    #     logger_process.warning(const)
                    
                        
                    n_pions += 1
                    if constP4.P() > 5: # No se pq está este corte
                        root_histograms["Reco"]["Events"]["PionTheta"].Fill(constP4.Theta())

                elif abs(const.getPDG()) == 2112:  # Neutron
                    n_neutrons += 1
            
            true_predicted_label["Countpions"].append(n_pions)
            true_predicted_label["Countphotons"].append(n_photons)
            true_predicted_label["Countneutrons"].append(n_neutrons)
            

            # Exploratorio de las migraciones 0 -> h1g y 0 -> h2g: un pi+- gen
            # que acaba reconstruido con uno o dos fotones de mas. Se rellena
            # una entrada por foton, asi que los casos de 2 fotones contribuyen
            # dos veces a los histogramas de foton y de angulo.
            if genTauId == 0 and recoTauId in [1]:
                migPionP4 = None
                migPhotonsP4 = []
                for c in range(0, recoTauNConsts):
                    const = recoTauConsts[c]
                    constP4 = ROOT.TLorentzVector()
                    try:
                        constP4.SetXYZM(
                            const.getMomentum().x,
                            const.getMomentum().y,
                            const.getMomentum().z,
                            const.getMass(),
                        )
                    except AttributeError:
                        constP4.SetXYZM(
                            const.getMomentum().X(),
                            const.getMomentum().Y(),
                            const.getMomentum().Z(),
                            const.getMass(),
                        )
                    if const.getPDG() == 22:
                        migPhotonsP4.append(ROOT.TLorentzVector(constP4))
                    elif abs(const.getPDG()) == 211 and (
                        migPionP4 is None or constP4.P() > migPionP4.P()
                    ):
                        # El pion de la desintegracion es el mas energetico
                        migPionP4 = ROOT.TLorentzVector(constP4)

                migTag = migration_label(recoTauId)
                if migPionP4 is not None:
                    root_histograms["Matched"]["Events"][f"Mig0To{migTag}PionP"].Fill(
                        migPionP4.P()
                    )
                    root_histograms["Matched"]["Events"][f"Mig0To{migTag}PionPPhotonP"].Fill(
                        migPionP4.P(), migPhotonsP4[0].P()
                    )
                    root_histograms["Matched"]["Events"][f"Mig0To{migTag}PionPPhotonPRelative"].Fill(
                                            migPhotonsP4[0].P() / (migPionP4.P() + 1e-6)
                                        )
                for migPhotonP4 in migPhotonsP4:
                    root_histograms["Matched"]["Events"][f"Mig0To{migTag}PhotonP"].Fill(
                        migPhotonP4.P()
                    )
                    if migPionP4 is not None:
                        root_histograms["Matched"]["Events"][f"Mig0To{migTag}PionPhotonAngle"].Fill(
                            migPionP4.Angle(migPhotonP4.Vect())
                        )
            
            if genTauId == 1 and recoTauId in [1]:
                            migPionP4 = None
                            migPhotonsP4 = []
                            for c in range(0, recoTauNConsts):
                                const = recoTauConsts[c]
                                constP4 = ROOT.TLorentzVector()
                                try:
                                    constP4.SetXYZM(
                                        const.getMomentum().x,
                                        const.getMomentum().y,
                                        const.getMomentum().z,
                                        const.getMass(),
                                    )
                                except AttributeError:
                                    constP4.SetXYZM(
                                        const.getMomentum().X(),
                                        const.getMomentum().Y(),
                                        const.getMomentum().Z(),
                                        const.getMass(),
                                    )
                                if const.getPDG() == 22:
                                    migPhotonsP4.append(ROOT.TLorentzVector(constP4))
                                elif abs(const.getPDG()) == 211 and (
                                    migPionP4 is None or constP4.P() > migPionP4.P()
                                ):
                                    # El pion de la desintegracion es el mas energetico
                                    migPionP4 = ROOT.TLorentzVector(constP4)
            
                            migTag = migration_label(recoTauId)
                            if migPionP4 is not None:
                                root_histograms["Matched"]["Events"][f"Mig1To{migTag}PionP"].Fill(
                                    migPionP4.P()
                                )
                                root_histograms["Matched"]["Events"][f"Mig1To{migTag}PionPPhotonP"].Fill(
                                    migPionP4.P(), migPhotonsP4[0].P()
                                )
                                root_histograms["Matched"]["Events"][f"Mig1To{migTag}PionPPhotonPRelative"].Fill(
                                                        migPhotonsP4[0].P() / (migPionP4.P() + 1e-6)
                                                    )
                            for migPhotonP4 in migPhotonsP4:
                                root_histograms["Matched"]["Events"][f"Mig1To{migTag}PhotonP"].Fill(
                                    migPhotonP4.P()
                                )
                                if migPionP4 is not None:
                                    root_histograms["Matched"]["Events"][f"Mig1To{migTag}PionPhotonAngle"].Fill(
                                        migPionP4.Angle(migPhotonP4.Vect())
                                    )

            if n_pi0s > 0:
                logger_pi0mass.debug(
                    f"Found {n_pi0s} pi0s ({n_photons} photons) in the matched reco with Id {recoTauId} tau with real Id {genTauId}"
                )

                #  Rho decay Pi P
                if  genTauId == 1 and recoTauId == 2:
                    photonCumP4 = ROOT.TLorentzVector()
                    photonCumP4.SetXYZM(0, 0, 0, 0)
                    
                    for c in range(0, recoTauNConsts):
                        const = recoTauConsts[c]
                        constP4 = ROOT.TLorentzVector()
                        try:
                            constP4.SetXYZM(
                                const.getMomentum().x,
                                const.getMomentum().y,
                                const.getMomentum().z,
                                const.getMass(),
                            )
                        except AttributeError:
                            # If the const does not have a mass, we assume it is a photon
                            constP4.SetXYZM(
                                const.getMomentum().X(),
                                const.getMomentum().Y(),
                                const.getMomentum().Z(),
                                const.getMass(),
                            )
                        if const.getCharge() == 0:
                            root_histograms["Reco"]["Events"]["RhoTwoPhotonDecayP"].Fill(constP4.P())
                            photonCumP4 += constP4
                        else:
                            root_histograms["Reco"]["Events"]["RhoPiDecayP"].Fill(constP4.P())
                            pionP4 = constP4
                            
                    root_histograms["Reco"]["Events"]["RhoTwoPhotonDecaySumP"].Fill(photonCumP4.P())
                    root_histograms["Reco"]["Events"]["RhoTwoPhotonDecayPiPhotonSumP"].Fill(
                        pionP4.P(), photonCumP4.P()
                    )
                    # Gen level
                    photonCumP4 = ROOT.TLorentzVector()
                    photonCumP4.SetXYZM(0, 0, 0, 0)
                    for c in range(0, genTauNConsts):
                        const = genTauConsts[c]
                        constP4 = ROOT.TLorentzVector()
                        try:
                            constP4.SetXYZM(
                                const.getMomentum().x,
                                const.getMomentum().y,
                                const.getMomentum().z,
                                const.getMass(),
                            )
                        except AttributeError:
                            # If the const does not have a mass, we assume it is a photon
                            constP4.SetXYZM(
                                const.getMomentum().X(),
                                const.getMomentum().Y(),
                                const.getMomentum().Z(),
                                const.getMass(),
                            )
                        if const.getCharge() == 0:
                            dau = const.getDaughters()
                            for dau in const.getDaughters():
                                photonP = ROOT.TLorentzVector()
                                try:
                                    photonP.SetXYZM(
                                        dau.getMomentum().x,
                                        dau.getMomentum().y,
                                        dau.getMomentum().z,
                                        dau.getMass(),
                                    )
                                except AttributeError:
                                    # If the dau does not have a mass, we assume it is a photon
                                    photonP.SetXYZM(
                                        dau.getMomentum().X(),
                                        dau.getMomentum().Y(),
                                        dau.getMomentum().Z(),
                                        dau.getMass(),
                                    )
                                root_histograms["Gen"]["Events"]["RhoTwoPhotonDecayP"].Fill(photonP.P())
                                photonCumP4 += photonP
                        else:
                            root_histograms["Gen"]["Events"]["RhoPiDecayP"].Fill(constP4.P())
                            pionP4 = constP4
                    root_histograms["Gen"]["Events"]["RhoTwoPhotonDecaySumP"].Fill(photonCumP4.P())
                    root_histograms["Gen"]["Events"]["RhoTwoPhotonDecayPiPhotonSumP"].Fill(
                        pionP4.P(), photonCumP4.P()
                    )
                if  recoTauId == 1 and genTauId == 1:
                    for c in range(0, recoTauNConsts):
                        const = recoTauConsts[c]
                        constP4 = ROOT.TLorentzVector()
                        try:
                            constP4.SetXYZM(
                                const.getMomentum().x,
                                const.getMomentum().y,
                                const.getMomentum().z,
                                const.getMass(),
                            )
                        except AttributeError:
                            constP4.SetXYZM(
                                const.getMomentum().X(),
                                const.getMomentum().Y(),
                                const.getMomentum().Z(),
                                const.getMass(),
                            )
                        if const.getCharge() == 0:
                            photonP = constP4
                            root_histograms["Reco"]["Events"]["RhoOnePhotonDecayPhotonP"].Fill(constP4.P())
                        else:
                            root_histograms["Reco"]["Events"]["RhoOnePhotonDecayPiP"].Fill(constP4.P())
                            pionP4 = constP4
                    root_histograms["Reco"]["Events"]["RhoOnePhotonDecayPiPhotonP"].Fill(pionP4.P(), photonP.P())
                    # Gen level
                    photonCumP4 = ROOT.TLorentzVector()
                    photonCumP4.SetXYZM(0, 0, 0, 0)
                    photonsP = []
                    for c in range(0, genTauNConsts):
                        const = genTauConsts[c]
                        constP4 = ROOT.TLorentzVector()
                        try:
                            constP4.SetXYZM(
                                const.getMomentum().x,
                                const.getMomentum().y,
                                const.getMomentum().z,
                                const.getMass(),
                            )
                        except AttributeError:
                            constP4.SetXYZM(
                                const.getMomentum().X(),
                                const.getMomentum().Y(),
                                const.getMomentum().Z(),
                                const.getMass(),
                            )
                        if const.getPDG() == 111:
                            pi0daughters = const.getDaughters()
                            for dau in pi0daughters:
                                photonP = ROOT.TLorentzVector()
                                try:
                                    photonP.SetXYZM(
                                        dau.getMomentum().x,
                                        dau.getMomentum().y,
                                        dau.getMomentum().z,
                                        dau.getMass(),
                                    )
                                except AttributeError:
                                    photonP.SetXYZM(
                                        dau.getMomentum().X(),
                                        dau.getMomentum().Y(),
                                        dau.getMomentum().Z(),
                                        dau.getMass(),
                                    )
                                root_histograms["Gen"]["Events"]["RhoOnePhotonDecayPhotonP"].Fill(photonP.P())
                                photonCumP4 += photonP
                                photonsP.append(photonP)
                        else:
                            root_histograms["Gen"]["Events"]["RhoOnePhotonDecayPiP"].Fill(constP4.P())
                            pionP4 = constP4
                            
                    root_histograms["Gen"]["Events"]["RhoOnePhotonDecayPhotonSumP"].Fill(photonCumP4.P())
                    root_histograms["Gen"]["Events"]["RhoOnePhotonDecayPiPhotonSumP"].Fill(pionP4.P(), photonCumP4.P())
                    ang = myutils.dRAngle(photonsP[0], photonsP[1])
                    root_histograms["Gen"]["Events"]["RhoOnePhotonDecayPhotonAng"].Fill(ang)

                # print(recoTaus_photons.keys())
                if n_photons == 3 and (genTauId == 2 or genTauId == 1):
                    logger_pi0mass.debug(
                        f"Found 3 photons in the matched reco tau with real Id {genTauId}"
                    )
                    pi0Mass_strmass, noMatchedPhotons = pi0Reco.getPi0Mass(recoTaus_photons, strategy = {"mass":-1})
                    if pi0Mass_strmass:
                        root_histograms["Reco"]["Events"]["ConstPi0MassFromPhotonMass"].Fill(pi0Mass_strmass)
                    
                    if noMatchedPhotons:
                        for ide, photon in noMatchedPhotons.items(): 
                            PhotonP4 = ROOT.TLorentzVector()
                            try:
                                PhotonP4.SetXYZM(
                                    photon.getMomentum().x,
                                    photon.getMomentum().y,
                                    photon.getMomentum().z,
                                    photon.getMass(),
                                )
                            except AttributeError:
                                # If the photon does not have a mass, we assume it is a photon
                                PhotonP4.SetXYZM(
                                    photon.getMomentum().X(),
                                    photon.getMomentum().Y(),
                                    photon.getMomentum().Z(),
                                    photon.getMass(),
                                )
                            if genTauId == 2:
                                root_histograms["Reco"]["Events"]["ConstlessPhotonPa1strMass"].Fill(PhotonP4.P())
                                root_histograms["Reco"]["Events"]["ConstlessPhotonPa1strMassZoom"].Fill(PhotonP4.P())

                            elif genTauId == 1:
                                root_histograms["Reco"]["Events"]["ConstxtraPhotonPrhostrMass"].Fill(PhotonP4.P())
                                root_histograms["Reco"]["Events"]["ConstxtraPhotonPrhostrMassZoom"].Fill(PhotonP4.P())
                                
                        matched_keys = [key for key in range(3) if key not in noMatchedPhotons.keys()]
                        first_matched_P = ROOT.TLorentzVector()
                        if matched_keys:
                            try:
                                first_matched_P.SetXYZM(
                                    recoTaus_photons[matched_keys[0]].getMomentum().x,
                                    recoTaus_photons[matched_keys[0]].getMomentum().y,
                                    recoTaus_photons[matched_keys[0]].getMomentum().z,
                                    recoTaus_photons[matched_keys[0]].getMass(),
                                )
                                second_matched_P = ROOT.TLorentzVector()
                                second_matched_P.SetXYZM(
                                    recoTaus_photons[matched_keys[1]].getMomentum().x,
                                    recoTaus_photons[matched_keys[1]].getMomentum().y,
                                    recoTaus_photons[matched_keys[1]].getMomentum().z,
                                    recoTaus_photons[matched_keys[1]].getMass(),
                                )
                            except AttributeError:
                                # If the photon does not have a mass, we assume it is a photon
                                first_matched_P.SetXYZM(
                                    recoTaus_photons[matched_keys[0]].getMomentum().X(),
                                    recoTaus_photons[matched_keys[0]].getMomentum().Y(),
                                    recoTaus_photons[matched_keys[0]].getMomentum().Z(),
                                    recoTaus_photons[matched_keys[0]].getMass(),
                                )
                                second_matched_P = ROOT.TLorentzVector()
                                second_matched_P.SetXYZM(
                                    recoTaus_photons[matched_keys[1]].getMomentum().X(),
                                    recoTaus_photons[matched_keys[1]].getMomentum().Y(),
                                    recoTaus_photons[matched_keys[1]].getMomentum().Z(),
                                    recoTaus_photons[matched_keys[1]].getMass(),
                                )
                            non_matched_P = PhotonP4
                            # Moment of the photons
                            root_histograms["Reco"]["Events"]["ThreePhotonMatchOnestrMassP"].Fill(first_matched_P.P())
                            root_histograms["Reco"]["Events"]["ThreePhotonMatchTwostrMassP"].Fill(second_matched_P.P())
                            root_histograms["Reco"]["Events"]["ThreePhotonNoMatchstrMassP"].Fill(non_matched_P.P())

                    pi0Mass_strdist, noMatchedPhotons = pi0Reco.getPi0Mass(recoTaus_photons, strategy = {"distance":-1})
                    # print(pi0Mass_strdist, noMatchedPhotons)
                    if pi0Mass_strdist:
                        root_histograms["Reco"]["Events"]["ConstPi0MassFromPhotonDist"].Fill(pi0Mass_strdist)
                    if noMatchedPhotons:
                        for ide, photon in noMatchedPhotons.items(): 
                            PhotonP4 = ROOT.TLorentzVector()
                            try:
                                PhotonP4.SetXYZM(
                                    photon.getMomentum().x,
                                    photon.getMomentum().y,
                                    photon.getMomentum().z,
                                    photon.getMass(),
                                )
                            except AttributeError:
                                # If the photon does not have a mass, we assume it is a photon
                                PhotonP4.SetXYZM(
                                    photon.getMomentum().X(),
                                    photon.getMomentum().Y(),
                                    photon.getMomentum().Z(),
                                    photon.getMass(),
                                )
                            if genTauId == 2:
                                root_histograms["Reco"]["Events"]["ConstlessPhotonPa1strDist"].Fill(PhotonP4.P())
                            elif genTauId == 1:
                                root_histograms["Reco"]["Events"]["ConstxtraPhotonPrhostrDist"].Fill(PhotonP4.P())
                        matched_keys = [key for key in range(3) if key not in noMatchedPhotons.keys()]
                        first_matched_P = ROOT.TLorentzVector()
                        # print("no matched photons", noMatchedPhotons.keys())
                        # print("matched Keys", matched_keys)
                        try:
                            first_matched_P.SetXYZM(
                                recoTaus_photons[matched_keys[0]].getMomentum().x,
                                recoTaus_photons[matched_keys[0]].getMomentum().y,
                                recoTaus_photons[matched_keys[0]].getMomentum().z,
                                recoTaus_photons[matched_keys[0]].getMass(),
                            )
                            second_matched_P = ROOT.TLorentzVector()
                            second_matched_P.SetXYZM(
                                recoTaus_photons[matched_keys[1]].getMomentum().x,
                                recoTaus_photons[matched_keys[1]].getMomentum().y,
                                recoTaus_photons[matched_keys[1]].getMomentum().z,
                                recoTaus_photons[matched_keys[1]].getMass(),
                            )
                        except AttributeError:
                            # If the photon does not have a mass, we assume it is a photon
                            first_matched_P.SetXYZM(
                                recoTaus_photons[matched_keys[0]].getMomentum().X(),
                                recoTaus_photons[matched_keys[0]].getMomentum().Y(),
                                recoTaus_photons[matched_keys[0]].getMomentum().Z(),
                                recoTaus_photons[matched_keys[0]].getMass(),
                            )
                            second_matched_P = ROOT.TLorentzVector()
                            second_matched_P.SetXYZM(
                                recoTaus_photons[matched_keys[1]].getMomentum().X(),
                                recoTaus_photons[matched_keys[1]].getMomentum().Y(),
                                recoTaus_photons[matched_keys[1]].getMomentum().Z(),
                                recoTaus_photons[matched_keys[1]].getMass(),
                            )
                        non_matched_P = PhotonP4
                        # Moment of the photons
                        root_histograms["Reco"]["Events"]["ThreePhotonMatchOnestrDistP"].Fill(first_matched_P.P())
                        root_histograms["Reco"]["Events"]["ThreePhotonMatchTwostrDistP"].Fill(second_matched_P.P())
                        root_histograms["Reco"]["Events"]["ThreePhotonNoMatchstrDistP"].Fill(non_matched_P.P())
                                
                elif n_photons == 1 and (genTauId == 0 or genTauId == 1):
                    logger_pi0mass.debug(
                        f"Found 1 photons in the matched reco tau with real Id {genTauId}"
                    )
                    if genTauId == 0:
                        root_histograms["Reco"]["Events"]["ConstxtraPhotonPi"].Fill(photonCumulativeP4.P())
                        root_histograms["Reco"]["Events"]["ConstxtraPhotonPiZoom"].Fill(photonCumulativeP4.P())
                    elif genTauId == 1:
                        root_histograms["Reco"]["Events"]["ConstlessPhotonPrho"].Fill(photonCumulativeP4.P())
                        root_histograms["Reco"]["Events"]["ConstlessPhotonPrhoZoom"].Fill(photonCumulativeP4.P())
                        
                
                elif n_photons == 2:
                    if genTauId == 1:
                        logger_pi0mass.debug(
                            f"Found 2 photons in the matched reco tau with real Id {genTauId}"
                        )
                        root_histograms["Reco"]["Events"]["ConstPi0MassRhoGen"].Fill(photonCumulativeP4.M())
                        photon1_P4 = ROOT.TLorentzVector()
                        try:
                            photon1_P4.SetXYZM(
                                recoTaus_photons[0].getMomentum().x,
                                recoTaus_photons[0].getMomentum().y,
                                recoTaus_photons[0].getMomentum().z,
                                recoTaus_photons[0].getMass(),
                            )
                            photon2_P4 = ROOT.TLorentzVector()
                            photon2_P4.SetXYZM(
                                recoTaus_photons[1].getMomentum().x,
                                recoTaus_photons[1].getMomentum().y,
                                recoTaus_photons[1].getMomentum().z,
                                recoTaus_photons[1].getMass(),
                            )
                        except AttributeError:
                            # If the photon does not have a mass, we assume it is a photon
                            photon1_P4.SetXYZM(
                                recoTaus_photons[0].getMomentum().X(),
                                recoTaus_photons[0].getMomentum().Y(),
                                recoTaus_photons[0].getMomentum().Z(),
                                recoTaus_photons[0].getMass(),
                            )
                            photon2_P4 = ROOT.TLorentzVector()
                            photon2_P4.SetXYZM(
                                recoTaus_photons[1].getMomentum().X(),
                                recoTaus_photons[1].getMomentum().Y(),
                                recoTaus_photons[1].getMomentum().Z(),
                                recoTaus_photons[1].getMass(),
                            )
                        ang_dist = myutils.dRAngle(
                            photon1_P4, photon2_P4
                        )
                        root_histograms["Reco"]["Events"]["ConstTwoPhotonAngDist"].Fill(ang_dist)
                    root_histograms["Reco"]["Events"]["ConstPi0MassRhoReco"].Fill(photonCumulativeP4.M())
                


            root_histograms["Reco"]["Events"]["TauType"].Fill(recoTauId)
            root_histograms["Reco"]["Events"]["TauPt"].Fill(recoTauP4.Pt())
            root_histograms["Reco"]["Events"]["TauP"].Fill(recoTauP4.P())

            root_histograms["Reco"]["Events"]["TauMass"].Fill(recoTauP4.M())
            root_histograms["Reco"]["Events"]["TauQ"].Fill(recoTauQ)
            root_histograms["Reco"]["Events"]["TauEta"].Fill(recoTauP4.Eta())
            root_histograms["Reco"]["Events"]["TauTheta"].Fill(recoTauP4.Theta())

            root_histograms["Reco"]["Events"]["TauDR"].Fill(recoTauDR)

            root_histograms["Matched"]["Events"]["TauPt"].Fill(genTauP4.Pt())
            root_histograms["Matched"]["Events"]["TauVisPt"].Fill(genVisTauP4.Pt())
            root_histograms["Matched"]["Events"]["TauP"].Fill(genTauP4.P())
            root_histograms["Matched"]["Events"]["TauVisP"].Fill(genVisTauP4.P())

            root_histograms["Matched"]["Events"]["TauVisMass"].Fill(genVisTauP4.M())
            root_histograms["Matched"]["Events"]["TauType"].Fill(genTauId)
            root_histograms["Matched"]["Events"]["TauQ"].Fill(genTauQ)
            root_histograms["Matched"]["Events"]["TauEta"].Fill(genTauP4.Eta())
            root_histograms["Matched"]["Events"]["TauTheta"].Fill(genTauP4.Theta())

            root_histograms["Matched"]["Events"]["TauDR"].Fill(genTauDR)

            root_histograms["GenVSReco"]["Events"]["TauPt"].Fill(recoTauP4.Pt(), genVisTauP4.Pt())
            root_histograms["GenVSReco"]["Events"]["TauMass"].Fill(recoTauP4.M(), genVisTauP4.M())
            root_histograms["GenVSReco"]["Events"]["TauType"].Fill(recoTauId, genTauId)
            root_histograms["GenVSReco"]["Events"]["TauP"].Fill(recoTauP4.P(), genVisTauP4.P())

            root_histograms["GenVSReco"]["Events"]["TauDR"].Fill(recoTauDR, genTauDR)

            if genTauId >=0:
                root_histograms["Matched"]["Events"]["HadronTauP"].Fill(genTauP4.P())
            
            if genTauId == 0:
                root_histograms["Matched"]["Events"]["TauP0"].Fill(genTauP4.P())
                root_histograms["Matched"]["Events"]["TauVisP0"].Fill(genVisTauP4.P())
                root_histograms["Matched"]["Events"]["TauTheta0"].Fill(genTauP4.Theta())
                root_histograms["Reco"]["Resolution"]["TauTheta0"].Fill(
                    (recoTauP4.Theta() - genVisTauP4.Theta())
                    / (genVisTauP4.Theta()+1e-10)
                )
                root_histograms["Reco"]["Resolution"]["TauP0"].Fill(
                    (recoTauP4.P() - genVisTauP4.P()) / (genVisTauP4.P()+1e-10)
                )
                if recoDM == 0:
                    root_histograms["Matched"]["Events"]["TauP0Correct"].Fill(genTauP4.P())
                    
                    root_histograms["Matched"]["Resolution"]["TauP0"].Fill(
                        (recoTauP4.P() - genVisTauP4.P()) / (genVisTauP4.P()+1e-10)
                    )
                    root_histograms["Matched"]["Resolution"]["TauTheta0"].Fill(
                        (recoTauP4.Theta() - genVisTauP4.Theta())
                        / (genVisTauP4.Theta()+1e-10)
                    )
                    
            elif genTauId == 1:
                root_histograms["Matched"]["Events"]["TauP1"].Fill(genTauP4.P())
                root_histograms["Matched"]["Events"]["TauVisP1"].Fill(genVisTauP4.P())
                root_histograms["Matched"]["Events"]["TauTheta1"].Fill(genTauP4.Theta())
                root_histograms["Reco"]["Resolution"]["TauTheta1"].Fill(
                    (recoTauP4.Theta() - genVisTauP4.Theta())
                    / (genVisTauP4.Theta()+1e-10)
                )
                root_histograms["Reco"]["Resolution"]["TauP1"].Fill(
                    (recoTauP4.P() - genVisTauP4.P()) / (genVisTauP4.P()+1e-10)
                )
                if recoDM == 1:
                    root_histograms["Matched"]["Events"]["TauP1Correct"].Fill(genTauP4.P())
                    
                    root_histograms["Matched"]["Resolution"]["TauP1"].Fill(
                        (recoTauP4.P() - genVisTauP4.P()) / (genVisTauP4.P()+1e-10)
                    )
                    root_histograms["Matched"]["Resolution"]["TauTheta1"].Fill(
                        (recoTauP4.Theta() - genVisTauP4.Theta())
                        / (genVisTauP4.Theta()+1e-10)
                    )
            elif genTauId == 2:
                root_histograms["Matched"]["Events"]["TauP2"].Fill(genTauP4.P())
                root_histograms["Matched"]["Events"]["TauVisP2"].Fill(genVisTauP4.P())
                root_histograms["Matched"]["Events"]["TauTheta2"].Fill(genTauP4.Theta())
                root_histograms["Reco"]["Resolution"]["TauTheta2"].Fill(
                    (recoTauP4.Theta() - genVisTauP4.Theta())
                    / (genVisTauP4.Theta()+1e-10)
                )
                root_histograms["Reco"]["Resolution"]["TauP2"].Fill(
                    (recoTauP4.P() - genVisTauP4.P()) / (genVisTauP4.P()+1e-10)
                )
                if recoDM == 2:
                    root_histograms["Matched"]["Events"]["TauP2Correct"].Fill(genTauP4.P())
                    
                    root_histograms["Matched"]["Resolution"]["TauP2"].Fill(
                        (recoTauP4.P() - genVisTauP4.P()) / (genVisTauP4.P()+1e-10)
                    )
                    root_histograms["Matched"]["Resolution"]["TauTheta2"].Fill(
                        (recoTauP4.Theta() - genVisTauP4.Theta())
                        / (genVisTauP4.Theta()+1e-10)
                    )
            elif genTauId == 10:
                root_histograms["Matched"]["Events"]["TauP10"].Fill(genTauP4.P())
                root_histograms["Matched"]["Events"]["TauVisP10"].Fill(genVisTauP4.P())
                root_histograms["Matched"]["Events"]["TauTheta10"].Fill(genTauP4.Theta())
                root_histograms["Reco"]["Resolution"]["TauTheta10"].Fill(
                    (recoTauP4.Theta() - genVisTauP4.Theta())
                    / (genVisTauP4.Theta()+1e-10)
                )
                root_histograms["Reco"]["Resolution"]["TauP10"].Fill(
                    (recoTauP4.P() - genVisTauP4.P()) / (genVisTauP4.P()+1e-10)
                )
                if recoDM == 10:
                    root_histograms["Matched"]["Events"]["TauP10Correct"].Fill(genTauP4.P())
                    
                    root_histograms["Matched"]["Resolution"]["TauP10"].Fill(
                        (recoTauP4.P() - genVisTauP4.P()) / (genVisTauP4.P()+1e-10)
                    )
                    root_histograms["Matched"]["Resolution"]["TauTheta10"].Fill(
                        (recoTauP4.Theta() - genVisTauP4.Theta())
                        / (genVisTauP4.Theta()+1e-10)
                    )
            
            # Fill the histograms for the matched reco tau with gen level information
            if recoDM == 0:
                root_histograms["Reco"]["Events"]["TauTheta0"].Fill(recoTauP4.Theta())
                root_histograms["Reco"]["Events"]["TauP0"].Fill(recoTauP4.P())
            elif recoDM == 1:
                root_histograms["Reco"]["Events"]["TauTheta1"].Fill(recoTauP4.Theta())
                root_histograms["Reco"]["Events"]["TauP1"].Fill(recoTauP4.P())
                root_histograms["Reco"]["Events"]["RhoMass"].Fill(recoTauP4.M())
            elif recoDM == 2:
                root_histograms["Reco"]["Events"]["TauTheta2"].Fill(recoTauP4.Theta())
                root_histograms["Reco"]["Events"]["TauP2"].Fill(recoTauP4.P())
                root_histograms["Reco"]["Events"]["a1OnePionMass"].Fill(recoTauP4.M())
            elif recoDM == 10:
                root_histograms["Reco"]["Events"]["TauTheta10"].Fill(recoTauP4.Theta())
                root_histograms["Reco"]["Events"]["TauP10"].Fill(recoTauP4.P())
                root_histograms["Reco"]["Events"]["a1ThreePionMass"].Fill(recoTauP4.M())
                
            
            
            # Resolution plots:
            if genVisTauP4.P() != 0:
                root_histograms["Reco"]["Resolution"]["TauPt"].Fill((recoTauP4.Pt() - genVisTauP4.Pt()) / genVisTauP4.Pt())
                root_histograms["Reco"]["Resolution"]["TauMass"].Fill((recoTauP4.M() - genVisTauP4.M()) / genVisTauP4.M())
                root_histograms["Reco"]["Resolution"]["TauP"].Fill((recoTauP4.P() - genVisTauP4.P()) / genVisTauP4.P())

        
        if len(GenTausP4)==2:
            cumP4 = GenTausP4[0] + GenTausP4[1]
            cumVisP4 = GenTausVisP4[0] + GenTausVisP4[1]

            root_histograms["Gen"]["Events"]["ZMass"].Fill(cumP4.M())
            root_histograms["Gen"]["Events"]["ZVisMass"].Fill(cumVisP4.M())
            gen_vis_angle_between_taus = myutils.dRAngle(GenTausVisP4[0], GenTausVisP4[1])
            gen_angle_between_taus = myutils.dRAngle(GenTausP4[0], GenTausP4[1])

            root_histograms["Gen"]["Events"]["TauPairVisAngle"].Fill(gen_vis_angle_between_taus)
            root_histograms["Gen"]["Events"]["TauPairAngle"].Fill(gen_angle_between_taus)

        if len(RecoTausP4)==2:
            cumRecoP4 = RecoTausP4[0] + RecoTausP4[1]
            root_histograms["Reco"]["Events"]["ZMass"].Fill(cumRecoP4.M())
            angle_between_taus = myutils.dRAngle(RecoTausP4[0], RecoTausP4[1])
            root_histograms["Reco"]["Events"]["TauPairAngle"].Fill(angle_between_taus)
        # print ("Taus???",nGenTaus,nTaus)
        root_histograms["Reco"]["Events"]["NTaus"].Fill(nRecoTaus)
        root_histograms["Reco"]["Events"]["NTausType"].Fill(nTausType)
        root_histograms["Gen"]["Events"]["NTausType"].Fill(nGenTausType)
        root_histograms["Gen"]["Events"]["NTaus"].Fill(nGenTausHad)


# Do efficiencies (divide matched gen by all gen)
# hEffiGenPi0Mass = root_histograms["Matched"]["Events"]["ConstPi0Mass"].Clone()
# hEffiGenPi0Mass.SetName("hEffiGenPi0Mass")
# hEffiGenPi0Mass.Divide(root_histograms["Gen"]["Events"]["ConstPi0Mass"])
# root_histograms["Matched"]["Effi"]["Pi0Mass"] = hEffiGenPi0Mass

# hEffiGenTauPt = root_histograms["Matched"]["Events"]["TauPt"].Clone()
# hEffiGenTauPt.SetName("hEffiGenTauPt")
# hEffiGenTauPt.Divide(root_histograms["Gen"]["Events"]["TauPt"])
# root_histograms["Matched"]["Effi"]["TauPt"] = hEffiGenTauPt

# hEffiGenVisTauPt = root_histograms["Matched"]["Events"]["TauVisPt"].Clone()
# hEffiGenVisTauPt.SetName("hEffiGenVisTauPt")
# hEffiGenVisTauPt.Divide(root_histograms["Gen"]["Events"]["TauVisPt"])
# root_histograms["Matched"]["Effi"]["VisTauPt"] = hEffiGenVisTauPt

# hEffiGenTauP = root_histograms["Matched"]["Events"]["TauP"].Clone()
# hEffiGenTauP.SetName("hEffiGenTauP")
# hEffiGenTauP.Divide(root_histograms["Gen"]["Events"]["TauP"])
# root_histograms["Matched"]["Effi"]["TauP"] = hEffiGenTauP


# hEffiGenVisTauP = root_histograms["Matched"]["Events"]["TauVisP"].Clone()
# hEffiGenVisTauP.SetName("hEffiGenVisTauP")
# hEffiGenVisTauP.Divide(root_histograms["Gen"]["Events"]["TauVisP"])
# root_histograms["Matched"]["Effi"]["VisTauP"] = hEffiGenVisTauP

# hEffiGenVisTauMass = root_histograms["Matched"]["Events"]["TauVisMass"].Clone()
# hEffiGenVisTauMass.SetName("hEffiGenVisTauMass")
# hEffiGenVisTauMass.Divide(root_histograms["Gen"]["Events"]["TauVisMass"])
# root_histograms["Matched"]["Effi"]["VisTauMass"] = hEffiGenVisTauMass

# hEffiGenTauEta = root_histograms["Matched"]["Events"]["TauEta"].Clone()
# hEffiGenTauEta.SetName("hEffiGenTauEta")
# hEffiGenTauEta.Divide(root_histograms["Gen"]["Events"]["TauEta"])
# root_histograms["Matched"]["Effi"]["TauEta"] = hEffiGenTauEta

# hEffiGenTauTheta = root_histograms["Matched"]["Events"]["TauTheta"].Clone()
# hEffiGenTauTheta.SetName("hEffiGenTauTheta")
# hEffiGenTauTheta.Divide(root_histograms["Gen"]["Events"]["TauTheta"])
# root_histograms["Matched"]["Effi"]["TauTheta"] = hEffiGenTauTheta

# hEffiGenTauType = root_histograms["Matched"]["Events"]["TauType"].Clone()
# hEffiGenTauType.SetName("hEffiGenTauType")
# hEffiGenTauType.Divide(root_histograms["Gen"]["Events"]["TauType"])
# root_histograms["Matched"]["Effi"]["TauType"] = hEffiGenTauType

# hEffiGenPhotonP = ROOT.TGraphAsymmErrors()
# hEffiGenPhotonP.SetName("hEffiGenPhotonP")
# hEffiGenPhotonP.Divide(root_histograms["Matched"]["Events"]["PhotonP"], root_histograms["Gen"]["Events"]["PhotonP"], "cl=0.683 b(1,1) mode")
# root_histograms["Matched"]["Effi"]["PhotonP"] = hEffiGenPhotonP

# hEffiGenPhotonTheta = ROOT.TGraphAsymmErrors()
# hEffiGenPhotonTheta.SetName("hEffiGenPhotonTheta")
# hEffiGenPhotonTheta.Divide(root_histograms["Matched"]["Events"]["PhotonTheta"], root_histograms["Gen"]["Events"]["PhotonTheta"], "cl=0.683 b(1,1) mode")
# root_histograms["Matched"]["Effi"]["PhotonTheta"] = hEffiGenPhotonTheta

# hEffiGenPionP = ROOT.TGraphAsymmErrors()
# hEffiGenPionP.Divide(root_histograms["Matched"]["Events"]["PionP"], root_histograms["Gen"]["Events"]["PionP"], "cl=0.683 b(1,1) mode")
# hEffiGenPionP.SetName("hEffiGenTauPionP")
# root_histograms["Matched"]["Effi"]["PionP"] = hEffiGenPionP

# hEffiGenPionTheta = ROOT.TGraphAsymmErrors()
# hEffiGenPionTheta.SetName("hEffiGenTauPionTheta")
# hEffiGenPionTheta.Divide(root_histograms["Matched"]["Events"]["PionTheta"], root_histograms["Gen"]["Events"]["PionTheta"], "cl=0.683 b(1,1) mode")
# root_histograms["Matched"]["Effi"]["PionTheta"] = hEffiGenPionTheta

# hEffiAllPionsP = ROOT.TGraphAsymmErrors()
# hEffiAllPionsP.Divide(root_histograms["Matched"]["Events"]["AllPionsP"], root_histograms["Gen"]["Events"]["AllPionsP"], "cl=0.683 b(1,1) mode")
# hEffiAllPionsP.SetName("hEffiGenAllPionsP")
# root_histograms["Matched"]["Effi"]["AllPionsP"] = hEffiAllPionsP

# root_histograms["Gen"]["Events"]["AllPionsP"].Sumw2() # COMPARAR CON MARÍA
# root_histograms["Matched"]["Events"]["AllPionsP"].Sumw2()

# eff_PionP_Pion=root_histograms["Matched"]["Events"]["AllPionsP"].Clone()
# eff_PionP_Pion.SetName("eff_PionP_Pion")
# eff_PionP_Pion.Divide(root_histograms["Gen"]["Events"]["AllPionsP"])
# root_histograms["Matched"]["Effi"]["AllPionsP_check"] = eff_PionP_Pion

# hEffiAllPionsTheta = ROOT.TGraphAsymmErrors()
# hEffiAllPionsTheta.Divide(root_histograms["Matched"]["Events"]["AllPionsTheta"], root_histograms["Gen"]["Events"]["AllPionsTheta"], "cl=0.683 b(1,1) mode")
# hEffiAllPionsTheta.SetName("hEffiGenAllPionsTheta")
# root_histograms["Matched"]["Effi"]["AllPionsTheta"] = hEffiAllPionsTheta

# hGenUnmatchedPionsPNorm = ROOT.TGraphAsymmErrors()
# hGenUnmatchedPionsPNorm.SetName("hGenUnmatchedPionsPNorm")
# hGenUnmatchedPionsPNorm.Divide(root_histograms["Gen"]["Events"]["UnmatchedPionsP"], root_histograms["Gen"]["Events"]["AllPionsP"], "cl=0.683 b(1,1) mode")
# root_histograms["Gen"]["Effi"]["UnmatchedPionsP"] = hGenUnmatchedPionsPNorm

# hGenUnmatchedPionsThetaNorm = ROOT.TGraphAsymmErrors()
# hGenUnmatchedPionsThetaNorm.SetName("hGenUnmatchedPionsThetaNorm")
# hGenUnmatchedPionsThetaNorm.Divide(root_histograms["Gen"]["Events"]["UnmatchedPionsTheta"], root_histograms["Gen"]["Events"]["AllPionsTheta"], "cl=0.683 b(1,1) mode")
# root_histograms["Gen"]["Effi"]["UnmatchedPionsTheta"] = hGenUnmatchedPionsThetaNorm

# print("ERROR POR BINS EN UNMATCHED")
# hGenMissmatchedPionsPNorm = ROOT.TGraphAsymmErrors()
# hGenMissmatchedPionsPNorm.SetName("hGenMissmatchedPionsPNorm")
# hGenMissmatchedPionsPNorm.Divide(root_histograms["Gen"]["Events"]["MissmatchedPionsP"], root_histograms["Gen"]["Events"]["AllPionsP"], "cl=0.683 b(1,1) mode")
# root_histograms["Gen"]["Effi"]["MissmatchedPionsP"] = hGenMissmatchedPionsPNorm

# hGenMissmatchedPionsThetaNorm = ROOT.TGraphAsymmErrors()
# hGenMissmatchedPionsThetaNorm.SetName("hGenMissmatchedPionsThetaNorm")
# hGenMissmatchedPionsThetaNorm.Divide(root_histograms["Gen"]["Events"]["MissmatchedPionsTheta"], root_histograms["Gen"]["Events"]["AllPionsTheta"], "cl=0.683 b(1,1) mode")
# root_histograms["Gen"]["Effi"]["MissmatchedPionsTheta"] = hGenMissmatchedPionsThetaNorm

# hEffiAllPhotonsP = ROOT.TGraphAsymmErrors()
# hEffiAllPhotonsP.SetName("hEffiGenAllPhotonsP")
# hEffiAllPhotonsP.Divide(root_histograms["Matched"]["Events"]["AllPhotonsP"], root_histograms["Gen"]["Events"]["AllPhotonsP"], "cl=0.683 b(1,1) mode")
# root_histograms["Matched"]["Effi"]["AllPhotonsP"] = hEffiAllPhotonsP

# hEffiAllPhotonsTheta = ROOT.TGraphAsymmErrors()
# hEffiAllPhotonsTheta.SetName("hEffiGenAllPhotonsTheta")
# hEffiAllPhotonsTheta.Divide(root_histograms["Matched"]["Events"]["AllPhotonsTheta"], root_histograms["Gen"]["Events"]["AllPhotonsTheta"], "cl=0.683 b(1,1) mode")
# root_histograms["Matched"]["Effi"]["AllPhotonsTheta"] = hEffiAllPhotonsTheta

# hGenUnmatchedPhotonsPNorm = ROOT.TGraphAsymmErrors()
# hGenUnmatchedPhotonsPNorm.SetName("hGenUnmatchedPhotonsPNorm")
# hGenUnmatchedPhotonsPNorm.Divide(root_histograms["Gen"]["Events"]["UnmatchedPhotonsP"], root_histograms["Gen"]["Events"]["AllPhotonsP"], "cl=0.683 b(1,1) mode")
# root_histograms["Gen"]["Effi"]["UnmatchedPhotonsP"] = hGenUnmatchedPhotonsPNorm

# hGenUnmatchedPhotonsThetaNorm = ROOT.TGraphAsymmErrors()
# hGenUnmatchedPhotonsThetaNorm.SetName("hGenUnmatchedPhotonsThetaNorm")
# hGenUnmatchedPhotonsThetaNorm.Divide(root_histograms["Gen"]["Events"]["UnmatchedPhotonsTheta"], root_histograms["Gen"]["Events"]["AllPhotonsTheta"], "cl=0.683 b(1,1) mode")
# root_histograms["Gen"]["Effi"]["UnmatchedPhotonsTheta"] = hGenUnmatchedPhotonsThetaNorm

# hGenMissmatchedPhotonsPNorm = ROOT.TGraphAsymmErrors()
# hGenMissmatchedPhotonsPNorm.SetName("hGenMissmatchedPhotonsPNorm")
# hGenMissmatchedPhotonsPNorm.Divide(root_histograms["Gen"]["Events"]["MissmatchedPhotonsP"], root_histograms["Gen"]["Events"]["AllPhotonsP"], "cl=0.683 b(1,1) mode")
# root_histograms["Gen"]["Effi"]["MissmatchedPhotonsP"] = hGenMissmatchedPhotonsPNorm

# hGenMissmatchedPhotonsThetaNorm = ROOT.TGraphAsymmErrors()
# hGenMissmatchedPhotonsThetaNorm.SetName("hGenMissmatchedPhotonsThetaNorm")
# hGenMissmatchedPhotonsThetaNorm.Divide(root_histograms["Gen"]["Events"]["MissmatchedPhotonsTheta"], root_histograms["Gen"]["Events"]["AllPhotonsTheta"], "cl=0.683 b(1,1) mode")
# root_histograms["Gen"]["Effi"]["MissmatchedPhotonsTheta"] = hGenMissmatchedPhotonsThetaNorm


# root_histograms["Matched"]["Events"]["HadronTauP"].Sumw2()
# root_histograms["Gen"]["Events"]["HadronTauP"].Sumw2()

# hEffiGenHadronTauP = ROOT.TGraphAsymmErrors()
# hEffiGenHadronTauP.Divide(root_histograms["Matched"]["Events"]["HadronTauP"], root_histograms["Gen"]["Events"]["HadronTauP"], "cl=0.683 b(1,1) mode")
# hEffiGenHadronTauP.SetName("hEffiGenHadronTauP")
# root_histograms["Matched"]["Effi"]["HadronTauP"] = hEffiGenHadronTauP

# hEffiGenAllPionsPCut = root_histograms["Matched"]["Events"]["AllPionsPCut"].Clone()
# hEffiGenAllPionsPCut.SetName("hEffiGenAllPionsPCut")
# hEffiGenAllPionsPCut.Divide(root_histograms["Gen"]["Events"]["AllPionsPCut"])
# root_histograms["Matched"]["Effi"]["AllPionsPCut"] = hEffiGenAllPionsPCut

# hEffiGenAllPionsPionThetaCut = root_histograms["Matched"]["Events"]["AllPionsThetaCut"].Clone()
# hEffiGenAllPionsPionThetaCut.SetName("hEffiGenAllPionsPionThetaCut")
# hEffiGenAllPionsPionThetaCut.Divide(root_histograms["Gen"]["Events"]["AllPionsThetaCut"])
# root_histograms["Matched"]["Effi"]["AllPionsThetaCut"] = hEffiGenAllPionsPionThetaCut

# # Theta angle per decay
# hEffiGenTauTheta0 = root_histograms["Matched"]["Events"]["TauTheta0"].Clone()
# hEffiGenTauTheta0.SetName("hEffiGenTauTheta0")
# hEffiGenTauTheta0.Divide(root_histograms["Gen"]["Events"]["TauTheta0"])
# root_histograms["Matched"]["Effi"]["TauTheta0"] = hEffiGenTauTheta0

# hEffiGenTauTheta1 = root_histograms["Matched"]["Events"]["TauTheta1"].Clone()
# hEffiGenTauTheta1.SetName("hEffiGenTauTheta1")
# hEffiGenTauTheta1.Divide(root_histograms["Gen"]["Events"]["TauTheta1"])
# root_histograms["Matched"]["Effi"]["TauTheta1"] = hEffiGenTauTheta1

# hEffiGenTauTheta2 = root_histograms["Matched"]["Events"]["TauTheta2"].Clone()
# hEffiGenTauTheta2.SetName("hEffiGenTauTheta2")
# hEffiGenTauTheta2.Divide(root_histograms["Gen"]["Events"]["TauTheta2"])
# root_histograms["Matched"]["Effi"]["TauTheta2"] = hEffiGenTauTheta2

# hEffiGenTauTheta10 = root_histograms["Matched"]["Events"]["TauTheta10"].Clone()
# hEffiGenTauTheta10.SetName("hEffiGenTauTheta10")
# hEffiGenTauTheta10.Divide(root_histograms["Gen"]["Events"]["TauTheta10"])
# root_histograms["Matched"]["Effi"]["TauTheta10"] = hEffiGenTauTheta10

# root_histograms["Matched"]["Events"]["TauP0"].Sumw2()
# root_histograms["Matched"]["Events"]["TauP1"].Sumw2()
# root_histograms["Matched"]["Events"]["TauP2"].Sumw2()
# root_histograms["Matched"]["Events"]["TauP10"].Sumw2()
# root_histograms["Gen"]["Events"]["TauP0"].Sumw2()
# root_histograms["Gen"]["Events"]["TauP1"].Sumw2()
# root_histograms["Gen"]["Events"]["TauP2"].Sumw2()
# root_histograms["Gen"]["Events"]["TauP10"].Sumw2()

# hEffiGenTauP0 = ROOT.TGraphAsymmErrors()
# hEffiGenTauP0.Divide(root_histograms["Matched"]["Events"]["TauP0"], root_histograms["Gen"]["Events"]["TauP0"], "cl=0.683 b(1,1) mode")
# hEffiGenTauP0.SetName("hEffiGenTauP0")
# root_histograms["Matched"]["Effi"]["TauP0"] = hEffiGenTauP0

# hEffiGenTauP1 = ROOT.TGraphAsymmErrors()
# hEffiGenTauP1.Divide(root_histograms["Matched"]["Events"]["TauP1"], root_histograms["Gen"]["Events"]["TauP1"], "cl=0.683 b(1,1) mode")
# hEffiGenTauP1.SetName("hEffiGenTauP1")
# root_histograms["Matched"]["Effi"]["TauP1"] = hEffiGenTauP1

# hEffiGenTauP1Correct = ROOT.TGraphAsymmErrors()
# hEffiGenTauP1Correct.Divide(root_histograms["Matched"]["Events"]["TauP1Correct"], root_histograms["Gen"]["Events"]["TauP1"], "cl=0.683 b(1,1) mode")
# hEffiGenTauP1Correct.SetName("hEffiGenTauP1Correct")
# root_histograms["Matched"]["Effi"]["TauP1Correct"] = hEffiGenTauP1Correct

# hEffiGenTauP2 = ROOT.TGraphAsymmErrors()
# hEffiGenTauP2.Divide(root_histograms["Matched"]["Events"]["TauP2"], root_histograms["Gen"]["Events"]["TauP2"], "cl=0.683 b(1,1) mode")
# hEffiGenTauP2.SetName("hEffiGenTauP2")
# root_histograms["Matched"]["Effi"]["TauP2"] = hEffiGenTauP2

# hEffiGenTauP10 = ROOT.TGraphAsymmErrors()
# hEffiGenTauP10.Divide(root_histograms["Matched"]["Events"]["TauP10"], root_histograms["Gen"]["Events"]["TauP10"], "cl=0.683 b(1,1) mode")
# hEffiGenTauP10.SetName("hEffiGenTauP10")
# root_histograms["Matched"]["Effi"]["TauP10"] = hEffiGenTauP10

# # Vis Momentum per decay
# hEffiGenTauVisP0 = root_histograms["Matched"]["Events"]["TauVisP0"].Clone()
# hEffiGenTauVisP0.SetName("hEffiGenTauVisP0")
# hEffiGenTauVisP0.Divide(root_histograms["Gen"]["Events"]["TauVisP0"])
# root_histograms["Matched"]["Effi"]["TauVisP0"] = hEffiGenTauVisP0
# hEffiGenTauVisP1 = root_histograms["Matched"]["Events"]["TauVisP1"].Clone()
# hEffiGenTauVisP1.SetName("hEffiGenTauVisP1")
# hEffiGenTauVisP1.Divide(root_histograms["Gen"]["Events"]["TauVisP1"])
# root_histograms["Matched"]["Effi"]["TauVisP1"] = hEffiGenTauVisP1
# hEffiGenTauVisP2 = root_histograms["Matched"]["Events"]["TauVisP2"].Clone()
# hEffiGenTauVisP2.SetName("hEffiGenTauVisP2")
# hEffiGenTauVisP2.Divide(root_histograms["Gen"]["Events"]["TauVisP2"])
# root_histograms["Matched"]["Effi"]["TauVisP2"] = hEffiGenTauVisP2
# hEffiGenTauVisP10 = root_histograms["Matched"]["Events"]["TauVisP10"].Clone()
# hEffiGenTauVisP10.SetName("hEffiGenTauVisP10")
# hEffiGenTauVisP10.Divide(root_histograms["Gen"]["Events"]["TauVisP10"])
# root_histograms["Matched"]["Effi"]["TauVisP10"] = hEffiGenTauVisP10


logger_io.info("Processed %d events", countEvents)

outfile = ROOT.TFile(fileOutName, "RECREATE")
for key in root_histograms_super:
    root_histograms = root_histograms_super[key]
    true_predicted_label = true_predicted_label_super[key]
    result_labels = results_labels_super[key]
    suffix = "" if key == "original" else f"_{key}"
    myutils.write_plot_config(root_histograms, outputpath, suffix)
# exit(0)
    root_histograms = myutils.calc_efficiency(root_histograms, histogram_config, suffix)
    write_histograms_recursive(root_histograms)




    decaystr = "decayAll" if selectDecay == -777 else "decay{}".format(selectDecay)
    true_predicted_label_output_file = outputpath + f"true_predicted_label_{decaystr}{suffix}.csv"
    true_predicted_label_df = pd.DataFrame(true_predicted_label)
    true_predicted_label_df.to_csv(true_predicted_label_output_file, index=False)

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



if type(run_config["output"]["outputlabels"]) is not list:
    if run_config["output"]["outputlabels"] is None:
        run_config["output"]["outputlabels"] = []
    else:
        run_config["output"]["outputlabels"] = [run_config["output"]["outputlabels"]]
if true_predicted_label_output_file not in run_config["output"]["outputlabels"]:
    run_config["output"]["outputlabels"].append(outputpath + f"true_predicted_label_{decaystr}.csv")

output_config_file = outputpath + "config.yaml"
with open(output_config_file, "w") as file:
    yaml.dump(run_config, file)
    logger_io.info("Configuration file saved to %s", output_config_file)



logger_io.info("Output file %s", outputpath + fileOutName)
logger_io.info("End of job")
outfile.Close()
