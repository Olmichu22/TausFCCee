import sys
import math
import ROOT
import yaml
import os
import shutil
import argparse
import glob
import copy
import logging
import pandas as pd
import pickle
import numpy as np
import pprint
from pathlib import Path
from modules.ParticleObjects import GenParticle, RecoParticle



def fit_sigma_energy(sigma_results):
    """Fit sigma(E) vs E using sigma(E) = a/sqrt(E) ⊕ b

    Args:
        sigma_results (dict): Dict containing 'energy', 'sigma' and 'sigma_err' lists.
    Returns:
        popt (array): Best-fit parameters [a, b]
        pcov (2D array): Covariance matrix
    """

    # Convert inputs to numpy arrays
    E = np.asarray(sigma_results['E_center'])
    sigma = np.asarray(sigma_results['sigma'])
    sigma_err = np.asarray(sigma_results['sigma_err'])

    # Model: sigma(E) = sqrt( (a/sqrt(E))^2 + b^2 )
    def sigma_model(E, a, b):
        return np.sqrt((a / np.sqrt(E))**2 + b**2)

    # Initial guess
    p0 = [1.0, 0.1]

    # Fit
    popt, pcov = curve_fit(
        sigma_model,
        E,
        sigma,
        sigma=sigma_err,
        absolute_sigma=True,
        p0=p0
    )

    return popt, pcov

def compute_sigma_from_hist(hist, edges, values):
    """
    Compute FWHM by fitting a Gaussian to the histogram produced by numpy,
    assuming the distribution is approximately normal. The fit range is fixed
    to [-0.2, 0.2].

    Parameters
    ----------
    hist : np.ndarray
        Array of bin counts from numpy.histogram.
    edges : np.ndarray
        Bin edges returned by numpy.histogram (same length = len(hist)+1).

    Returns
    -------
    fwhm : float
        Full width at half maximum (FWHM) = 2.355 * sigma_fit.
    fwhm_err : float
        Uncertainty in FWHM from the sigma fit uncertainty.
    """

    # ---- Crear histograma ROOT a partir del histograma numpy ----
    nbins = len(hist)

    # Crear TH1F
    h = ROOT.TH1F("h_tmp", "resolution", nbins, float(edges[0]), float(edges[-1]))

    # Rellenar con los contenidos de numpy
    for i in range(nbins):
        h.SetBinContent(i+1, float(hist[i]))

    # Si no hay estadística suficiente, devolver NaN
    if h.GetEntries() < 5:
        return np.nan, np.nan

    # ---- Definir función gaussiana ----
    gaus = ROOT.TF1("gaus_tmp", "gaus", -0.2, 0.2)

    # Parámetros iniciales razonables (evita fallos de ajuste)
    peak = h.GetMaximum()
    mean_est = h.GetMean()
    sigma_est = h.GetRMS()

    gaus.SetParameters(peak, mean_est, sigma_est)

    # ---- Ajuste gaussiano en rango [-0.2, 0.2] ----
    fit_result = h.Fit(gaus, "SQN", "", -0.2, 0.2)
    # S = silent, Q = quiet, N = no drawing

    if int(fit_result) != 0:
        # Ajuste fallido
        return np.nan, np.nan
    # Extraer media y su error
    mean_val = np.mean(values)
    mean_err = np.std(values) / math.sqrt(len(values))
    # ---- Extraer sigma y error ----
    sigma = gaus.GetParameter(2)
    sigma_err = gaus.GetParError(2)

    # ---- Calcular FWHM ----
    # fwhm = 2.355 * sigma
    # fwhm_err = 2.355 * sigma_err
    fwhm = sigma
    fwhm_err = sigma_err
    return fwhm, fwhm_err

logger_io = logging.getLogger('io')

PLOT_1D_TEMPLATE = {"title": "", "x": "", "y": "", "fit": False, "fitrange": [0,0]}

PLOT_2D_TEMPLATE = {"title": "", "x": "", "y": ""}



def clone_histograms_with_suffix(hist_dict, suffix):
    """
    Clone all ROOT histograms inside an arbitrarily nested dictionary
    and append a suffix ('_min' or '_max') to their internal ROOT name.

    The dictionary structure and keys remain identical.
    """
    new_dict = {}

    for key, value in hist_dict.items():
        # Caso 1: El valor es otro subdiccionario -> recursión
        if isinstance(value, dict):
            new_dict[key] = clone_histograms_with_suffix(value, suffix)

        # Caso 2: Es un histograma ROOT (TH1, TH2…)
        elif isinstance(value, ROOT.TH1) or isinstance(value, ROOT.TH2):
            # Clon seguro del histograma
            cloned = value.Clone()
            cloned.SetDirectory(0)

            # Cambiar nombre interno del histograma ROOT
            cloned.SetName(value.GetName() + suffix)

            new_dict[key] = cloned

        else:
            # Por si apareciera algo inesperado
            new_dict[key] = value

    return new_dict


def write_plot_config(root_histograms, outputpath, suffix=""):
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
    plot_config_name = os.path.join(outputpath, f"plot_config{suffix}.yaml")
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

def _make_hist(root_name, var_cfg):
    """Create a single ROOT TH1F or TH2F from v2-config parameters."""
    hist_type = var_cfg.get("type", "1D")
    if hist_type == "1D":
        bins = var_cfg.get("bins", 10)
        r0, r1 = var_cfg.get("range", [0, 1])
        return ROOT.TH1F(root_name, "", bins, r0, r1)
    elif hist_type == "2D":
        bx = var_cfg.get("x_bins", 10)
        by = var_cfg.get("y_bins", 10)
        rx0, rx1 = var_cfg.get("x_range", [0, 1])
        ry0, ry1 = var_cfg.get("y_range", [0, 1])
        return ROOT.TH2F(root_name, "", bx, rx0, rx1, by, ry0, ry1)
    raise ValueError(f"Unknown histogram type '{hist_type}' for '{root_name}'")


def build_histogram_registry(config: dict) -> dict:
    """Build a 4-level histogram registry from a v2 YAML config.

    Reads the ``histograms`` section of *config* and expands each variable
    entry into one ROOT histogram per (variable, category, weight) triple.

    Root histogram naming convention::

        {variable}_{category}              # weight == nominal
        {variable}_{category}_{weight}     # weight in [P1, M1, …]

    Args:
        config: dict loaded from rho_analysis_config_v2.yml.

    Returns:
        hists[level][variable][category][weight] → ROOT.TH1F / ROOT.TH2F
    """
    defaults = config["defaults"]["categories"]
    hists = {}

    for level, variables in config["histograms"].items():
        hists[level] = {}
        for var_name, var_cfg in variables.items():
            hists[level][var_name] = {}
            cats_raw = var_cfg["categories"]

            # ``categories`` may be a list (use global defaults) or a dict
            # (per-category weight overrides).
            if isinstance(cats_raw, list):
                cat_configs = {cat: defaults[cat] for cat in cats_raw}
            else:
                cat_configs = {
                    cat: (override if override is not None else defaults[cat])
                    for cat, override in cats_raw.items()
                }

            for cat, cat_cfg in cat_configs.items():
                hists[level][var_name][cat] = {}
                weight_variants = ["nominal"] + list(cat_cfg.get("weights", []))
                for w in weight_variants:
                    suffix = "" if w == "nominal" else f"_{w}"
                    root_name = f"{var_name}_{cat}{suffix}"
                    hists[level][var_name][cat][w] = _make_hist(root_name, var_cfg)

    return hists


def calc_efficiency(root_histograms, histograms_config, suffix = ""):

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
                    print(name + suffix)
                    efficiency_hist.SetName(name +  suffix)
                    root_histograms[hist_level][hist_class][hist_name] = efficiency_hist
                except KeyError as e:
                    logger_io.error(f"Error when calculating efficenciy for {name}: {e}")
                    continue
    return root_histograms

def associate_reco_with_gen_taus(gen_taus: dict[GenParticle],
                                 reco_tau: dict[RecoParticle]) -> tuple[int, float]:
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
    particle_analysis = False,
    log_subdir = None
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
    parser.add_argument("-t","--tauCut",default=0,type=float)
    parser.add_argument("-R", "--dRMax", type=float)
    parser.add_argument("-n", "--NeutronCut", type=float)
    parser.add_argument("-g", "--generalPCut", type=float)
    parser.add_argument("-r", "--MatchedGenMinDR", type=float)
    parser.add_argument(
        "-m", "--matchedCM",
        default="True",
        type=str,
        help="Use only matched taus to compute confusion matrix.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
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
        "--prefix", required=False
    )
    parser.add_argument(
        "--samples-config",
        type=str,
        default="config/samples/samples.yaml",
        help="YAML con el mapeo sample → ruta en disco.",
    )
    parser.add_argument(
        "--input-list",
        nargs="+",
        metavar="FILE",
        default=None,
        help="Uno o varios ficheros ROOT (rutas absolutas). Omite el escaneo de directorio.",
    )
    parser.add_argument(
        "--shard",
        type=str,
        metavar="I/N",
        default=None,
        help="Procesa solo el shard I de N (0-indexado). Reparte los ficheros de "
             "entrada en N grupos estriados; pensado para lanzar N procesos en "
             "paralelo y fusionar la salida después (ver runShardedPlots.py).",
    )

    if parser_hook is not None:
        parser_hook(parser)

    args = parser.parse_args()

    # Load config
    config = load_yaml_config(args.config, default_config)
    histograms_config = load_yaml_config(args.hist_config, None)

    # systematics error (if exists)
    if hasattr(args, "sys_err"):
        sys_err_file = args.sys_err
        if os.path.exists(sys_err_file):
            with open(sys_err_file, "r") as file:
                sys_err_config = yaml.safe_load(file)
            config["systematics_errors"] = sys_err_config
        else:
            raise FileNotFoundError(f"Error: Systematics error file '{sys_err_file}' does not exist.")
    
    
    # Cut Configuration
    cuts = config.get("cuts", {})
    for key in ["tauCut", "dRMax", "TauPhotonPCut", "TauPionPCut", "NeutronCut", "MatchedGenMinDR", "generalPCut"]:
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
    mdr = _first(cuts["MatchedGenMinDR"])

    suffix = f"_{dr}_tph{tph}_tpi{tpi}_n{npe}_g{gpc}"
    decay_str = f"decay{select_decay}" + suffix
    if select_decay == -777:
        decay_str = "decayAll" + suffix
    file_out = f"{outfile}{decay_str}.root"

    # Output path logic
    base = output_base + outfile + suffix[1:] + "/"
    if args.gatr_result and args.test_pfo and not args.prefix:
        path = output_base + "PFO_" + outfile + suffix[1:] + "/"
    elif args.gatr_result and args.prefix:
        path = output_base + args.prefix + outfile + suffix[1:] + "/"
    elif args.test_pfo:
        raise ValueError("Cannot use --test-pfo without --gatr-result.")
    if args.input_list and args.gatr_result:
        raise ValueError("Cannot use --input-list together with --gatr-result.")
    elif args.prefix:
        path = output_base + args.prefix + outfile + suffix[1:] + "/"
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

    # Carpeta de logs: si se indica log_subdir → path/logs/<log_subdir>/.
    # Se vacía en cada run para que solo contenga los logs de la última ejecución.
    if log_subdir:
        log_dir = os.path.join(path, "logs", log_subdir)
        if os.path.isdir(log_dir):
            # best-effort: en NFS, borrar un log aún abierto por otro proceso deja
            # un fichero ".nfs*" ocupado (Errno 16) que abortaría toda la corrida.
            # ignore_errors lo tolera; los logs nuevos se reescriben con mode="w".
            shutil.rmtree(log_dir, ignore_errors=True)
        os.makedirs(log_dir, exist_ok=True)
    else:
        log_dir = path

    if not exp:
        # Logging
        lvl = logging.WARNING if args.verbose == 0 else logging.INFO if args.verbose == 1 else logging.DEBUG
        app_log = os.path.join(log_dir, f"app_{decay_str}.log")
        handlers = []
        if args.verbose < 2:
            handlers = [logging.StreamHandler(sys.stdout), logging.FileHandler(app_log, mode="w")]
        elif args.verbose == 2:
            sh = logging.StreamHandler(sys.stdout); sh.setLevel(logging.DEBUG)
            fh = logging.FileHandler(app_log, mode="w"); fh.setLevel(logging.DEBUG)
            handlers = [sh, fh]
        else:
            handlers = [logging.FileHandler(app_log, mode="w")]

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
    test_mode = args.test

    # Leer has_gen_taus desde samples.yaml para el sample activo
    has_gen_taus = False
    _samples_cfg = getattr(args, "samples_config", "config/samples/samples.yaml")
    if os.path.exists(_samples_cfg):
        with open(_samples_cfg) as _f:
            _sdb = yaml.safe_load(_f)
        _current_sample = config["general"]["sample"]
        _entry = _sdb.get("samples", {}).get(_current_sample, {})
        has_gen_taus = bool(_entry.get("has_gen_taus", False))
        if not has_gen_taus:
            for _e in _sdb.get("samples", {}).values():
                if _current_sample and _current_sample.lower() in [a.lower() for a in _e.get("aliases", [])]:
                    has_gen_taus = bool(_e.get("has_gen_taus", False))
                    break

    return {
        "args": args,
        "config": config,
        "histograms_config": histograms_config,
        "outputpath": path,
        "logdir": log_dir,
        "fileOutName": file_out,
        "loggers": loggers,
        "decay": select_decay,
        "flags": {
            "matched_cm": matched_cm,
            "test": test_mode
        },
        "decay_str": decay_str,
        "has_gen_taus": has_gen_taus,
    }


def parse_shard_spec(shard):
    """Parse a ``--shard I/N`` string into the tuple ``(I, N)``.

    Returns ``None`` when *shard* is None/empty (sharding disabled).
    """
    if not shard:
        return None
    try:
        index_str, total_str = str(shard).split("/")
        shard_index, n_shards = int(index_str), int(total_str)
    except ValueError:
        raise ValueError(f"--shard debe tener el formato I/N (recibido: {shard!r})")
    if n_shards < 1:
        raise ValueError(f"--shard N debe ser >= 1 (recibido: {shard!r})")
    if not 0 <= shard_index < n_shards:
        raise ValueError(f"--shard I debe cumplir 0 <= I < N (recibido: {shard!r})")
    return shard_index, n_shards


def apply_shard(filenames, mlpf_results, shard, loggers):
    """Keep only the files belonging to one shard, remapping the GATr predictions.

    Files are split with a stride (``filenames[I::N]``) so the shards stay
    balanced even when file sizes vary.

    ``mlpf_results`` is keyed as ``file_position * 1000 + local_event`` — the
    same convention the event loop relies on when it looks up predictions by a
    running event id over ``root_io.Reader(filenames)``. Dropping files shifts
    every position, so the keys are rebuilt against the new positions.
    """
    spec = parse_shard_spec(shard)
    if spec is None:
        return filenames, mlpf_results

    shard_index, n_shards = spec
    if n_shards == 1:
        return filenames, mlpf_results

    kept = [(pos, name) for pos, name in enumerate(filenames) if pos % n_shards == shard_index]

    if mlpf_results:
        new_pos_of = {old_pos: new_pos for new_pos, (old_pos, _) in enumerate(kept)}
        remapped = {}
        for key, value in mlpf_results.items():
            old_pos, local = divmod(key, 1000)
            new_pos = new_pos_of.get(old_pos)
            if new_pos is not None:
                remapped[new_pos * 1000 + local] = value
        mlpf_results = remapped

    filenames = [name for _, name in kept]
    loggers["io"].info(
        "Shard %d/%d: %d file(s) of the original list, %d prediction(s).",
        shard_index, n_shards, len(filenames), len(mlpf_results),
    )
    return filenames, mlpf_results


def get_root_trees_path(sample, gatr_results_path, loggers, test, args=None, skip_root_validation: bool = False):
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
        args (argparse.Namespace or None): Optional argument parser namespace. Used to read
            ``args.input_list`` (explicit file list) and ``args.samples_config`` (YAML path).

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
            if test == True and i > 5:
                break
            
            mlpf_predictions_path = row[1]["prediction_file"]
            simulation_path = row[1]["simulation_file"]
            my_file = Path(simulation_path)
            loggers["io"].debug("Reading file %s", simulation_path)
            if my_file.is_file():
                if not skip_root_validation:
                    root_file = open_root_file(simulation_path)
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
            
    elif args is not None and args.input_list:
        # Explicit file list — skip directory scan entirely
        if test:
            loggers["io"].warning("--test has no effect when --input-list is used.")
        filenames = []
        for f in args.input_list:
            my_file = Path(f)
            if not my_file.is_file():
                loggers["io"].warning("File %s not found, skipping.", f)
                continue
            if not skip_root_validation:
                root_file = open_root_file(f)
                if not root_file or root_file.IsZombie():
                    loggers["io"].warning("File %s is a zombie or could not be opened.", f)
                    continue
            filenames.append(f)
        loggers["io"].info("Input-list mode: %d file(s) to process.", len(filenames))
    else:
        # Resolve sample → directory via samples YAML
        samples_config_path = args.samples_config if args is not None else "config/samples/samples.yaml"
        if not os.path.exists(samples_config_path):
            loggers["io"].error("Samples config file %s not found.", samples_config_path)
            sys.exit(1)
        with open(samples_config_path, "r") as f:
            samples_db = yaml.safe_load(f)

        default_base = samples_db.get("default_base", "")
        default_prefix = samples_db.get("default_file_prefix", "out_reco_edm4hep_edm4hep")

        # Build flat lookup: alias_lower → entry
        lookup = {}
        for name, entry in samples_db.get("samples", {}).items():
            lookup[name.lower()] = entry
            for alias in entry.get("aliases", []):
                lookup[alias.lower()] = entry

        entry = lookup.get(sample.lower())
        if entry is None:
            loggers["io"].error("Sample '%s' not found in %s.", sample, samples_config_path)
            sys.exit(1)

        entry_prefix = entry.get("file_prefix", default_prefix)
        entry_glob = entry.get("file_glob", samples_db.get("default_file_glob"))

        # Resolve one or more source directories for this sample.
        # Supported keys (singular = 1 dir, plural = list of dirs; can be combined):
        #   path  / paths   -> absolute directory path(s)
        #   folder/ folders -> name(s) relative to default_base
        # List elements may be plain strings or dicts {path|folder, file_prefix,
        # file_glob} to give a given directory its own prefix/glob.
        #
        # Two ways of naming the input files:
        #   * file_prefix -> files are numbered sequentially: "<prefix>_<i>.root"
        #   * file_glob   -> shell pattern matched inside the directory, for
        #                    samples whose file names carry arbitrary run ids
        #                    (e.g. "events_*_REC.edm4hep.root"). Takes priority
        #                    over file_prefix when both are given.
        dir_specs = []  # list of (dir_path, file_prefix, file_glob)

        def _add_dir(value, is_folder):
            if isinstance(value, (list, tuple)):
                for v in value:
                    _add_dir(v, is_folder)
            elif isinstance(value, dict):
                sub_prefix = value.get("file_prefix", entry_prefix)
                sub_glob = value.get("file_glob", entry_glob)
                if "path" in value:
                    dir_specs.append((value["path"], sub_prefix, sub_glob))
                elif "folder" in value:
                    dir_specs.append((os.path.join(default_base, value["folder"]), sub_prefix, sub_glob))
                else:
                    loggers["io"].warning("Ignoring dir spec without 'path'/'folder': %r", value)
            else:
                dp = os.path.join(default_base, value) if is_folder else value
                dir_specs.append((dp, entry_prefix, entry_glob))

        for v in entry.get("paths", []):
            _add_dir(v, is_folder=False)
        for v in entry.get("folders", []):
            _add_dir(v, is_folder=True)
        if "path" in entry:
            _add_dir(entry["path"], is_folder=False)
        if "folder" in entry:
            _add_dir(entry["folder"], is_folder=True)

        if not dir_specs:
            loggers["io"].error(
                "Sample '%s' has no 'path(s)'/'folder(s)' in %s.", sample, samples_config_path
            )
            sys.exit(1)

        bad_indices = set(entry.get("bad_file_indices", []))
        if bad_indices:
            loggers["io"].info("Skipping %d bad file indices for sample '%s': %s",
                               len(bad_indices), sample, sorted(bad_indices))

        filenames = []
        remaining = 100 if test else None  # --test caps the total number of files
        for dir_path, file_prefix, file_glob in dir_specs:
            if remaining is not None and remaining <= 0:
                break
            if not os.path.isdir(dir_path):
                loggers["io"].warning("Directory %s not found, skipping.", dir_path)
                continue

            if file_glob:
                # Glob mode: candidate names come from the directory listing.
                candidates = sorted(
                    f for f in glob.glob(os.path.join(dir_path, file_glob))
                    if os.path.isfile(f)
                )
                if not candidates:
                    loggers["io"].warning(
                        "Pattern '%s' matched no file in %s.", file_glob, dir_path
                    )
            else:
                # Sequential mode: names are rebuilt as "<prefix>_<i>.root".
                nfiles = sum(
                    1
                    for fname in os.listdir(dir_path)
                    if fname.endswith(".root") and os.path.isfile(os.path.join(dir_path, fname))
                )
                candidates = [
                    os.path.join(dir_path, f"{file_prefix}_{i}.root")
                    for i in range(1, nfiles + 1)
                ]

            if remaining is not None:
                candidates = candidates[:remaining]

            loggers["io"].info("Reading files from %s (%d files)", dir_path, len(candidates))
            # bad_file_indices son 1-based sobre esta lista de candidatos
            for i, filename in enumerate(candidates, start=1):
                if i in bad_indices:
                    loggers["io"].debug("Skipping bad file index %d (%s)", i, filename)
                    continue
                loggers["io"].debug("Reading file %s", filename)
                my_file = Path(filename)
                if my_file.is_file():
                    if not skip_root_validation:
                        root_file = open_root_file(filename)
                        if not root_file or root_file.IsZombie():
                            loggers["io"].warning("File %s is a zombie or could not be opened.", filename)
                            continue
                    filenames.append(filename)

            if remaining is not None:
                remaining = 100 - len(filenames)

        loggers["io"].info("Total files to process for sample '%s': %d", sample, len(filenames))

    if args is not None:
        filenames, mlpf_results = apply_shard(
            filenames, mlpf_results, getattr(args, "shard", None), loggers
        )
    return filenames, mlpf_results


def compute_photon_resolution_two_by_two(reco_photons, gen_photons):
    """
    Calcula la resolución (E_reco - E_gen) / E_gen para dos fotones RECO y dos GEN.

    Parámetros:
        reco_photons: lista de 2 TLorentzVector
        gen_photons:  lista de 2 TLorentzVector

    Retorna:
        (res1, res2): resoluciones asociadas tras minimizar ΔR total
    """

    if len(reco_photons) != 2 or len(gen_photons) != 2:
        raise ValueError("Se requieren exactamente 2 fotones reco y 2 gen")

    r1, r2 = reco_photons
    g1, g2 = gen_photons

    # ------------------------------
    # Combinación A: r1↔g1 , r2↔g2
    # ------------------------------
    dR11 = dRAngle(r1, g1)
    dR22 = dRAngle(r2, g2)
    costA = dR11 + dR22

    # ------------------------------
    # Combinación B: r1↔g2 , r2↔g1
    # ------------------------------
    dR12 = dRAngle(r1, g2)
    dR21 = dRAngle(r2, g1)
    costB = dR12 + dR21

    # Determinar emparejamiento óptimo
    if costA <= costB:
        pairs = [(r1, g1), (r2, g2)]
    else:
        pairs = [(r1, g2), (r2, g1)]

    # Calcular resoluciones
    resolutions = []
    angles = []
    gen_energies = []
    for reco_p4, gen_p4 in pairs:
        E_reco = reco_p4.P()
        E_gen  = gen_p4.P()
        res = (E_reco - E_gen) / E_gen if E_gen != 0 else 999.
        resolutions.append(res)
        angles.append(gen_p4.Theta())
        gen_energies.append(E_gen)

    return tuple(resolutions), tuple(angles), tuple(gen_energies)



def update_resolution_hist(resolution, E_gen, theta,
                           theta_bins_rad,
                           energy_bins,
                           hist_dict,
                           bin_edges):
    """
    Update the resolution histogram dictionary based on angle and true energy.

    Parameters
    ----------
    resolution : float
        Value of (E_reco - E_gen) / E_gen.
    E_gen : float
        True photon energy.
    theta : float
        Photon incident angle in radians.
    theta_bins_rad : dict
        Dictionary defining detector regions and their angular ranges in radians.
        Example: {"barrel": (theta_min, theta_max), ...}
    energy_bins : list of [float, float]
        List of energy bin intervals.
    hist_dict : dict
        Nested dictionary of histograms, structured as hist_dict[region][energy_bin].
        Each entry is a numpy array containing histogram counts.
    bin_edges : np.ndarray
        Array of histogram bin edges.

    Returns
    -------
    dict
        The updated histogram dictionary.
    """

    # ---- Determinar la región del detector según theta ----
    region = None
    for reg, (tmin, tmax) in theta_bins_rad.items():
        if tmin <= theta < tmax:
            region = reg
            break

    # Si no pertenece a ninguna región conocida, no actualizar nada
    if region is None:
        return hist_dict

    # ---- Determinar el bin de energía ----
    energy_bin = None
    for i, (emin, emax) in enumerate(energy_bins):
        if emin <= E_gen < emax:
            energy_bin = i
            break

    # Si la energía no cae en un bin válido, no actualizar nada
    if energy_bin is None:
        return hist_dict
    
    hist_dict[region][energy_bin].append(resolution)

    return hist_dict