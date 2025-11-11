import sys
import math
import ROOT
import yaml
import os
import argparse
import copy
import logging
import pandas as pd
import pickle
import numpy as np
import pprint
from pathlib import Path
from modules import myutils




logger_io = logging.getLogger('io')

PLOT_1D_TEMPLATE = {"title": "", "x": "", "y": "", "fit": False, "fitrange": [0,0]}

PLOT_2D_TEMPLATE = {"title": "", "x": "", "y": ""}

def write_plot_config(root_histograms, outputpath):
    """Write a the plot_config_file based on histogram info present in root_histograms

    Args:
        root_histograms (dict): Dictionary containing ROOT histograms.
        outputpath (str): Path to the output directory.
    """
    # Get all 1D histograms
    histograms_1d = []
    histograms_2d = []
    def extract_histograms(obj):
        if isinstance(obj, dict):
            for key, value in obj.items():
                extract_histograms(value)
        else:
            if not isinstance(obj, ROOT.TH1):
                histograms_2d.append((obj.GetName()))
            else:
                histograms_1d.append((obj.GetName()))
    extract_histograms(root_histograms)
    plot_config_name = os.path.join(outputpath, "plot_config.yaml")
    if os.path.exists(plot_config_name):
        logger_io.warning(f"{plot_config_name} exists and will be updated.")
        with open(plot_config_name, "r") as file:
            plot_config = yaml.safe_load(file)
    else:
        plot_config = {"variabs_hist": [], "plot_titles_config_hist":{}, "variabs_2d": [], "plot_titles_config_2d":{}}
    for var in histograms_1d:
        if var not in plot_config["variabs_hist"]:
            plot_config["variabs_hist"].append(var)
            plot_config["plot_titles_config_hist"][var] = copy.deepcopy(PLOT_1D_TEMPLATE)
    for var in histograms_2d:
        if var not in plot_config["variabs_2d"]:
            plot_config["variabs_2d"].append(var)
            plot_config["plot_titles_config_2d"][var] = copy.deepcopy(PLOT_2D_TEMPLATE)
    with open(plot_config_name, "w") as file:
        yaml.dump(plot_config, file, sort_keys=False)

def write_histograms_recursive(obj):
    """
    Traverse a nested dictionary and call `.Write()` on each object of type ROOT histogram.
    """
    if isinstance(obj, dict):
        for value in obj.values():
            write_histograms_recursive(value)
    else:
        # If it is not a dict, assume it is a root hist
        try:
            obj.Write()
        except AttributeError:
            print(f"Object {obj} has not .Write() method. Ignored.")


def set_up_root_histograms(histograms_config):
    """
    Configure root histograms based on the input config.
    Args:
        histograms_config (dict): Dict with hist configuration.
    Returns:
        dict: Nested dict with root histograms.
    """
    
    root_histograms = {}     
    for hist_level in histograms_config.keys():
        root_histograms[hist_level] = {}
        for category in histograms_config[hist_level].keys():
            root_histograms[hist_level][category] = {}
    
    for hist_level, category in histograms_config.items():
        for hist_class, hist_configs in category.items():
            if hist_configs is None:
                continue
            if hist_class == "Effi":
                continue
            for hist_name, hist_params in hist_configs.items():
                name = hist_params.get("name", hist_name)
                hist_type = hist_params.get("type", "1D")
                if hist_type == "1D":
                    bins = hist_params.get("bins", 10)
                    range_0, range_1 = hist_params.get("range", (0, 1))
                    root_histograms[hist_level][hist_class][hist_name] = ROOT.TH1F(
                        name, "", bins, range_0, range_1
                    )
                elif hist_type == "2D":
                    bins_x = hist_params.get("x_bins", 10)
                    bins_y = hist_params.get("y_bins", 10)
                    range_0x, range_1x = hist_params.get("x_range", [0, 1])
                    range_0y, range_1y = hist_params.get("y_range", [0, 1])
                    root_histograms[hist_level][hist_class][hist_name] = ROOT.TH2F(
                        name, "", bins_x, range_0x, range_1x, bins_y, range_0y, range_1y
                    )
    return root_histograms

def calc_efficiency(root_histograms, histograms_config):

    for hist_level, category in histograms_config.items():
        for hist_class, hist_configs in category.items():
            if hist_configs is None:
                continue
            if hist_class != "Effi":
                continue
            for hist_name, hist_params in hist_configs.items():
                name = hist_params.get("name", hist_name)
                num_hist_name = hist_params.get("numerator")
                num_hist_level, num_name = num_hist_name.split("/")
                denom_hist_name = hist_params.get("denominator")
                denom_hist_level, denom_name = denom_hist_name.split("/")
                try:
                    logger_io.debug(f"Calculating efficenciy for {name}: {num_hist_level}/{num_name} over {denom_hist_level}/{denom_name}")
                    numerator = root_histograms[num_hist_level]["Events"][num_name]
                    denominator = root_histograms[denom_hist_level]["Events"][denom_name]
                    efficiency_hist = ROOT.TGraphAsymmErrors(numerator, denominator, "cl=0.683 b(1,1) mode")
                    efficiency_hist.SetName(name)
                    root_histograms[hist_level][hist_class][hist_name] = efficiency_hist
                except KeyError as e:
                    logger_io.error(f"Error when calculating efficenciy for {name}: {e}")
                    continue
    return root_histograms

def associate_reco_with_gen_taus(gen_taus,
                                 reco_tau) -> tuple[int, float]:
    """

    Args:
        gen_taus (dict[GenParticle]): Dict of gen particles.
        reco_tau (dict[RecoParticle]): Dict of reco particles.

    Returns:
        (gen_tau_key, cos_gen_reco) (tupple[int, float]): Key to the nearest gen tau, and angle with respect to reco tau.
    """
    
    # Obtener dirección de cada tau
    tau_directions = []
    for key, tau in gen_taus.items():
        px = tau.getMomentum().X()
        py = tau.getMomentum().Y()
        pz = tau.getMomentum().Z()
        tau_directions.append((px, py, pz))
    
    reco_tau_direction = [reco_tau.getMomentum().X(),
                          reco_tau.getMomentum().Y(),
                          reco_tau.getMomentum().Z()]
    
    # Calcular cosenos de ángulos entre direcciones
    cos_r_tau1 = np.dot(reco_tau_direction, tau_directions[0]) / (np.linalg.norm(reco_tau_direction) * np.linalg.norm(tau_directions[0]))
    cos_r_tau2 = np.dot(reco_tau_direction, tau_directions[1]) / (np.linalg.norm(reco_tau_direction) * np.linalg.norm(tau_directions[1]))
    
    # El hemisferio 1 corresponde al tau 1 si el coseno es mayor
    if cos_r_tau1 > cos_r_tau2:
        return list(gen_taus.keys())[0], cos_r_tau1
    else:
        return list(gen_taus.keys())[1], cos_r_tau2

# I'm sure this exists already
def dRAngle(p1,p2):
    """
    Calculate the angle between two particles in the eta-phi plane
    Args:
        p1 (TLorentzVector): 4-momentum vector of the first particle
        p2 (TLorentzVector): 4-momentum vector of the second particle
    Returns:
        float: angle between the two particles in the theta-phi plane
    """
    dphi=p1.Phi()-p2.Phi()
    if (dphi>math.pi) : dphi=2*math.pi-dphi
    if (dphi<-math.pi) : dphi=2*math.pi+dphi
    dtheta=p1.Theta()-p2.Theta()
    dR=math.sqrt(dtheta*dtheta+dphi*dphi)
    return dR

# trick to prevent broken files (should not be a problem at CIEMAT)
def open_root_file(file_path):
    try:
        # Suppress ROOT's default error messages to the terminal
        ROOT.gErrorIgnoreLevel = ROOT.kError

        # Attempt to open the ROOT file in "READ" mode without auto-recovery
        root_file = ROOT.TFile.Open(file_path, "READ")

        # Check if the file is a zombie
        if not root_file or root_file.IsZombie():
            logger_io.error(f"Error: '{file_path}' is a zombie or could not be opened.")
            raise IOError(f"Error: '{file_path}' is a zombie or could not be opened.")
        
        # Check if file is recoverable (potentially corrupted)
        if root_file.TestBit(ROOT.TFile.kRecovered):
            logger_io.error(f"Warning: '{file_path}' is corrupted and has been recovered.")
            raise IOError(f"Error: '{file_path}' is corrupted and has been recovered.")
        
        #print(f"'{file_path}' opened successfully.")
        logger_io.debug("File '%s' opened successfully.", file_path)
        return root_file

    except Exception as e:
        logger_io.error("Error opening file '%s': %s", file_path, e)
        # print(f"Error: {e}")
        return None

# Fuction to sort by tau P
def sort_by_P(Tau):
    tau_with_P = []

    for i in range(0,len(Tau)):
        tau_with_P.append((Tau[i], Tau[i].getMomentum().P()))
    
    # Sort the list based on the P() value in descending order
    sorted_tau_with_P = sorted(tau_with_P, key=lambda x: x[1], reverse=True)
    
    # Extract only the sorted Tau[i] objects from the tuples
    sortedTau = [tau for tau, _ in sorted_tau_with_P]
   
    return sortedTau

def load_yaml_config(config_file, default_config):
    """
    Load the YAML configuration file if it exists.
    Args:
            args (argparse.Namespace): command line arguments
            config_file (str): path to the YAML configuration file
    Returns:
            dict: configuration parameters
    """
    if config_file is not None and os.path.exists(config_file):
        with open(config_file, "r") as file:
            config = yaml.safe_load(file)
            # print(f"Loaded configuration parameters from '{config_file}'.")
    elif default_config:
        if not os.path.exists(default_config):
            raise FileNotFoundError(f"Error: '{default_config}' does not exist. A valid default configuration file is required.")
        with open(default_config, "r") as file:
            config = yaml.safe_load(file)
            # print(f"Loaded default configuration parameters from '{default_config}'.")
    else:
        raise FileNotFoundError(f"Error: A valid default configuration file is required.")
    return config

def load_yaml_config(config_file, default_config):
    """
    Load the YAML configuration file if it exists.
    Args:
        args (argparse.Namespace): command line arguments
        config_file (str): path to the YAML configuration file
    Returns:
        config (dict): configuration parameters
    """
    if config_file is not None and os.path.exists(config_file):
        with open(config_file, "r") as file:
            config = yaml.safe_load(file)
            # print(f"Loaded configuration parameters from '{config_file}'.")
    elif default_config:
        if not os.path.exists(default_config):
            raise FileNotFoundError(f"Error: '{default_config}' does not exist. A valid default configuration file is required.")
        with open(default_config, "r") as file:
            config = yaml.safe_load(file)
            # print(f"Loaded default configuration parameters from '{default_config}'.")
    else:
        raise FileNotFoundError(f"Error: A valid default configuration file is required.")
    return config


def setup_analysis_config(
    default_config: str = "config/default/taurecolong.yaml",
    output_base: str = "Results/TauReco/",
    parser_hook=None,
    exp = False,
    particle_analysis = False
):
    """
    Encapsulates argument configuration, config loading, cut application,
    output path setup, and logging initialization.

    Args:
        default_config (str, optional): Path to the default YAML configuration file. 
            Defaults to "config/default/taurecolong.yaml".
        output_base (str, optional): Base directory for output results. 
            Defaults to "Results/TauReco/".
        parser_hook (_type_, optional): Optional hook to modify the argument parser. 
            Defaults to None.
        exp (bool, optional): Whether to enable experimental analysis mode. 
            Defaults to False.
        particle_analysis (bool, optional): Whether to enable particle-level analysis. 
            Defaults to False.

    Raises:
        ValueError: Raised if `--test-pfo` is used without `--gatr-result`.

    Returns:
        dict: Dictionary containing:
            - **args** (*argparse.Namespace*): Parsed command-line arguments.
            - **config** (*dict*): Updated configuration dictionary.
            - **outputpath** (*str*): Created output path.
            - **fileOutName** (*str*): Output file name.
            - **loggers** (*dict*): Dictionary of loggers (config, io, processing, pi0mass).
    """
    # Argument parser setup
    parser = argparse.ArgumentParser(
        description="Configure the analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    
    
    parser.add_argument("-f", "--sample")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-d", "--decay", type=int)
    parser.add_argument("-p", "--TauPhotonPCut", type=float)
    parser.add_argument("-i", "--TauPionPCut", type=float)
    parser.add_argument("-t","--tauCut",default=2,type=float)
    parser.add_argument("-R", "--dRMax", type=float)
    parser.add_argument("-n", "--NeutronCut", type=float)
    parser.add_argument("-g", "--generalPCut", type=float)
    parser.add_argument("-r", "--MatchedGenMaxDR", type=float)
    parser.add_argument(
        "-m", "--matchedCM",
        default="True",
        type=str,
        help="Use only matched taus to compute confusion matrix.",
    )
    parser.add_argument(
        "--test",
        type=str,
        help="Run in test mode with limited number of files",
    )
    parser.add_argument(
        "-c", "--config", type=str, help="Configuration file"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level: -v for INFO, -vv for DEBUG",
    )
    parser.add_argument(
        "--gatr-result",
        type=str,
        help="Path to GATR result for the analysis.",
    )
    parser.add_argument(
        "--test-pfo",
        action="store_true",
        help="Use this flag to test the PFOs in same files as GATr.",
    )
    parser.add_argument(
        "--hist-config",
        type=str,
        default="config/histograms/particles_config.yml",
        help="Path to the histogram configuration file.",
    )
    parser.add_argument(
        "--prefix", required = False
    )

    if parser_hook is not None:
        parser_hook(parser)

    args = parser.parse_args()

    # Load config
    config = load_yaml_config(args.config, default_config)
    histograms_config = load_yaml_config(args.hist_config, None)

    # Cut Configuration
    cuts = config.get("cuts", {})
    for key in ["tauCut", "dRMax", "TauPhotonPCut", "TauPionPCut", "NeutronCut", "MatchedGenMaxDR", "generalPCut"]:
        val = getattr(args, key) if getattr(args, key) is not None else cuts.get(key)
        cuts[key] = val
    config["cuts"] = cuts

    # Decay selection
    decay_list = config.setdefault("general", {}).setdefault("decay", [])
    select_decay = args.decay if args.decay is not None else decay_list[0]
    if args.decay is not None and args.decay not in decay_list:
        decay_list.append(args.decay)
    config["general"]["decay"] = decay_list

    # Output filename
    outfile = args.outfile or config["general"].get("outfile")
    config["general"]["outfile"] = outfile

    # Build strings
    def _first(val):
        return val[0] if isinstance(val, list) else val

    dr = _first(cuts["dRMax"])
    tph = _first(cuts["TauPhotonPCut"])
    tpi = _first(cuts["TauPionPCut"])
    npe = _first(cuts["NeutronCut"])
    gpc = _first(cuts["generalPCut"])
    mdr = _first(cuts["MatchedGenMaxDR"])

    suffix = f"_{dr}_tph{tph}_tpi{tpi}_n{npe}_g{gpc}"
    decay_str = f"decay{select_decay}" + suffix
    if select_decay == -777:
        decay_str = "decayAll" + suffix
    file_out = f"{outfile}{decay_str}.root"

    # Output path logic
    base = output_base + outfile + suffix[1:] + "/"
    if args.gatr_result and args.test_pfo and not args.prefix:
        path = output_base + "PFO_" + outfile + suffix[1:] + "/"
    if args.gatr_result and args.prefix:
        path = output_base + args.prefix + outfile + suffix[1:] + "/"
    elif args.test_pfo:
        raise ValueError("Cannot use --test-pfo without --gatr-result.")
    else:
        path = base
    if args.gatr_result:
        path = "GATr_" + path
    if particle_analysis:
        path = "ParticleEval_" + path
    os.makedirs(path, exist_ok=True)

    config.setdefault("output", {}).setdefault("outputfile", [])
    if not config["output"].get("outputfile") is None:
        if file_out not in config["output"]["outputfile"]:
            config["output"]["outputfile"].append(file_out)
    else:
        config["output"]["outputfile"] = [file_out]
    config["output"]["outputpath"] = path

    if not exp:
        # Logging
        lvl = logging.WARNING if args.verbose == 0 else logging.INFO if args.verbose == 1 else logging.DEBUG
        handlers = []
        if args.verbose < 2:
            handlers = [logging.StreamHandler(sys.stdout), logging.FileHandler(os.path.join(path, "app.log"), mode="w")]
        elif args.verbose == 2:
            sh = logging.StreamHandler(sys.stdout); sh.setLevel(logging.DEBUG)
            fh = logging.FileHandler(os.path.join(path, "app.log"), mode="w"); fh.setLevel(logging.DEBUG)
            handlers = [sh, fh]
        else:
            handlers = [logging.FileHandler(os.path.join(path, "app.log"), mode="w")]

        logging.basicConfig(
            level=lvl,
            format="%(asctime)s, %(levelname)s, [%(name)s] - %(message)s",
            handlers=handlers
        )

        loggers = {
            "config": logging.getLogger("config"),
            "io": logging.getLogger("io"),
            "processing": logging.getLogger("processing"),
            "pi0mass": logging.getLogger("pi0mass")
        }
        loggers["config"].info("Configuration loaded!")
        loggers["config"].info("Configuration:\n%s", pprint.pformat(config, indent=4))

    else:
        loggers = {}

    # General args to config
    for key in ["sample", "matchedCM", "test"]:
        config["general"][key] = getattr(args, key) if getattr(args, key) is not None else config["general"].get(key)

    # Convert flags
    matched_cm = True if config["general"]["matchedCM"] == "True" else False
    test_mode = True if config["general"]["test"] == "True" else False

    
    return {
        "args": args,
        "config": config,
        "histograms_config": histograms_config,
        "outputpath": path,
        "fileOutName": file_out,
        "loggers": loggers,
        "decay": select_decay,
        "flags": {
            "matched_cm": matched_cm,
            "test": test_mode
        },
        "decay_str": decay_str
    }

def get_single_root_file(file_path, loggers):
    """
    Reads a single ROOT file from the specified path and returns it in the same 
    format as `get_root_trees_path` (i.e., as a list of file paths).

    Args:
        file_path (str): Full path to the .root file to be read.
        loggers (dict): Dictionary of loggers with at least the "io" key used for 
                        logging information, warnings, and errors.

    Raises:
        SystemExit: If the file does not exist or cannot be opened.

    Returns:
        list[str]: A list containing the valid ROOT file path.
    """
    filenames = []

    if not os.path.exists(file_path):
        loggers["io"].error("Specified ROOT file does not exist: %s", file_path)
        sys.exit(1)

    loggers["io"].info("Reading ROOT file: %s", file_path)
    my_file = Path(file_path)

    if my_file.is_file():
        root_file = myutils.open_root_file(file_path)
        if not root_file or root_file.IsZombie():
            loggers["io"].warning("File %s is a zombie or could not be opened.", file_path)
            sys.exit(1)
        filenames.append(file_path)
    else:
        loggers["io"].error("Specified path is not a file: %s", file_path)
        sys.exit(1)

    return filenames


def get_root_trees_path(sample, gatr_results_path, loggers, test, path=None, file_prefix=None):
    """
    Loads ROOT file paths and associated GATr (Graph Analysis Training results) predictions 
    for a given sample. Handles both local GATr result files and simulation-only workflows.

    Depending on whether `gatr_results_path` is provided, it either:
      - Loads GATr prediction files and corresponding ROOT simulation files listed in a CSV file.
      - Or, if `gatr_results_path` is None, reads simulation ROOT files directly from a predefined path.

    Args:
        sample (str): Name of the dataset or sample to process (used when `gatr_results_path` is None).
        gatr_results_path (str or None): Path to a CSV file containing columns `prediction_file` and 
            `simulation_file`. If None, simulation files are loaded from the default path.
        loggers (dict): Dictionary of loggers with at least the `"io"` key used for logging 
            information, warnings, and errors.
        test (bool): If True, limits the processing to only one file for quick testing.
        path (str): If provided, uses this path to look for root_files.
        file_prefix (str): If provided, uses this prefix to look for root_files.

    Raises:
        SystemExit: If `gatr_results_path` is provided but the path does not exist.

    Returns:
        tuple:
            - **filenames** (*list[str]*): List of valid ROOT file paths to be processed.
            - **mlpf_results** (*dict*): Dictionary mapping unique event IDs to MLPF/GATr predictions 
              (empty if no `gatr_results_path` is provided).
    """
    mlpf_results = {}
    
    if gatr_results_path is not None:
        if not os.path.exists(gatr_results_path):
            loggers["io"].error("GATr results path %s does not exist.", gatr_results_path)
            sys.exit(1)
        else:
            loggers["io"].info("Using GATr results from %s", gatr_results_path)
        # abrimos archivo configuracion yml
        mlpf_config = pd.read_csv(gatr_results_path)
        filenames = []
        n_files = 0
        n_preds = 1
        for i, row in enumerate(mlpf_config.iterrows()):
            if test == True and i > 0:
                break
            
            mlpf_predictions_path = row[1]["prediction_file"]
            simulation_path = row[1]["simulation_file"]
            my_file = Path(simulation_path)
            loggers["io"].debug("Reading file %s", simulation_path)
            if my_file.is_file():
                root_file = myutils.open_root_file(simulation_path)
                if not root_file or root_file.IsZombie():
                    loggers["io"].warning("File %s is a zombie or could not be opened.", simulation_path)
                    continue
                filenames.append(simulation_path)
            
            with open(mlpf_predictions_path, "rb") as f:
                mlpf_preds_i = pickle.load(f)
            
            loggers["io"].debug("Read %d GATr results", len(mlpf_results))
                
            for key, value in mlpf_preds_i.items():
                key_id = n_files*1000 + key - 1
                mlpf_results[key_id] = value
                n_preds += 1
            n_files += 1

        loggers["io"].info("Total predictions loaded: %d", n_preds)
            
    else:
        # Simulation files
        if not path:
            path = "/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/"
            dir_path = path + "/" + sample
            file = "out_reco_edm4hep_edm4hep"
        else:
            dir_path = path
            file = file_prefix
        filenames = []

        nfiles = len(os.listdir(dir_path))

        nfiles = 1000
        if test == True:
            nfiles = 1

        loggers["io"].info("Reading files from %s", dir_path)
        for i in range(1, nfiles + 1):
            filename = dir_path + "/" + file + "_{}.root".format(i)
            loggers["io"].debug("Reading file %s", filename)
            my_file = Path(filename)
            if my_file.is_file():
                root_file = myutils.open_root_file(filename)
                if not root_file or root_file.IsZombie():
                    logger_io.warning("File %s is a zombie or could not be opened.", filename)
                    continue
                filenames.append(filename)
    return filenames, mlpf_results
