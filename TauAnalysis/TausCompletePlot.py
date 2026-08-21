import ROOT 
import argparse
import yaml
import pandas as pd
import os
import warnings
from modules.plotting import plot_1D_hist, plot_2D_hist, plot_hist_together, plot_cm, plot_hist_zoom

parser = argparse.ArgumentParser(description="Configure the plot",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", default="Results/TauReco/effis0.4_tph0.35_tpi0_n3_g0.0", help="Input path where data to process is stored")
parser.add_argument("-d","--decay", default=-777, type=str) #,"decay0","decay1","decay10"]
parser.add_argument("-p","--plotconfig",default="config/plots/taulong_plotconfig_results.yaml", type=str)
parser.add_argument(
    "-t", "--typeplot", 
    nargs="+",  # Permite múltiples valores separados por espacios
    default=["All"],  # Valor predeterminado
    type=str, 
    help="List of types of plots to generate. Choose between All, 1D, 2D, CM (Confusion Matrix), CMP (Confusion Matrix with Photons)"
)
parser.add_argument(
    "-s",
    "--same",
    default="True",
    type=str,
    help="If True, Reco, Gen and MatchedGen are plotted in the same plot",
)

type_plot = ["1D", "2D", "CM", "CMP"]

args = parser.parse_args()
inputBasePath = args.input
decay = args.decay
plot_config_path = args.plotconfig
typeplot = args.typeplot
if "All" in typeplot:
  typeplot = type_plot

if set(typeplot) - set(type_plot) != set():
  raise Exception(f"One type plot not valid. Choose between {type_plot}")

# Load Plot Config File
with open(plot_config_path, "r") as file:
  plot_config = yaml.safe_load(file)

# Load config file of the analysis
with open(inputBasePath+"/config.yaml", "r") as yaml_file:
  config = yaml.safe_load(yaml_file)

# print(config)
# Do not show statistics
ROOT.gStyle.SetOptStat(0)


decay_str = "decayAll" if decay == -777 else f"decay{decay}"

outputpath = inputBasePath+"/Images_"+decay_str

if not os.path.exists(outputpath):
    os.makedirs(outputpath)

input_file = None
labels_file = None
# Look for the output file in the outputfile list
for outputfile in config["output"]["outputfile"]:

  if decay_str in outputfile:
    input_file = str(outputfile)  # Ensure file is a string
    break
if config["output"].get("outputlabels") is not None:
  for outputfile in config["output"]["outputlabels"]:
    if decay_str in outputfile:
      labels_file = str(outputfile)  # Ensure labels_file is a string
      break
else:
  labels_file = None

if input_file == None:
  raise Exception(f"Output file not found for decay {decay_str}")
if labels_file == None:
  # Warning if labels file is not found
  warnings.warn(f"\n\n==============\n Labels file not found for decay {decay_str}. CM plot will not be generated.\n==============\n")

file_path = config["output"]["outputpath"]+"Histos_"+input_file
print("\n=====\n Trying to open file: ", file_path, "\n=====\n")
if not os.path.exists(file_path):
  file_path = config["output"]["outputpath"]+input_file
  print(file_path)
  if not os.path.exists(file_path):
    raise Exception(f"Input file {input_file} not found in path {config['output']['outputpath']} or {config['output']['outputpath']}Histos_")
true_predict_labels = labels_file

print(f"Input file: {file_path}")
# exit(0)
# Open the file
file = ROOT.TFile(file_path)

# def formatHisto(file,variab,rename,titleX,color=ROOT.kBlack, linestyle=1, fillstyle=3001):
#   histo = file.Get(variab)
#   histo.SetName(rename)
#   if titleX != None:
#     histo.SetXTitle(titleX)
#   histo.SetLineColor(color)
#   # histo.SetLineWidth(2)
#   histo.SetLineStyle(linestyle)
#   if fillstyle != None:
#     histo.SetFillStyle(fillstyle)
#     histo.SetFillColor(color)
#   #histo.SetMarkerColor(color)
#   #histo.SetMarkerStyle(20)
#   # histo.Sumw2()
#   return histo

# Select case for the type of plot
variabs1D = plot_config.get("variabs_hist", [])
labels1D = plot_config.get("plot_titles_config_hist", {})
variabs2D = plot_config.get("variabs_2d", [])
labels2D = plot_config.get("plot_titles_config_2d", {})

def check_variab_in_file(file, variabs):
  missing_variabs = []
  for variab in variabs:
    histo = file.Get(variab)
    if histo is None:
      missing_variabs.append(variab)
  return missing_variabs

if labels_file != None:
  results_df = pd.read_csv(true_predict_labels)

for typeplt in typeplot:
  if typeplt == "1D":
    variabs = plot_config.get("variabs_hist", [])
    labels = plot_config.get("plot_titles_config_hist", {})
    if args.same == "True":
      variabs_and_config = plot_config.get("plot_together", dict())
      plot_hist_together(file, variabs_and_config, outputpath)
    normalize = plot_config.get("norm", False)
    plot_1D_hist(file, variabs, labels, outputpath, normalize)

    zoom_config = plot_config.get("Zoom", None)
    if zoom_config:
      plot_hist_zoom(file, zoom_config, outputpath)
  elif typeplt == "2D" and variabs2D:
    variabs = plot_config.get("variabs_2d", [])
    labels = plot_config.get("plot_titles_config_2d", {})
    plot_2D_hist(file, variabs, labels, outputpath)
  elif typeplt == "CM" and labels_file != None:
    results_df = pd.read_csv(true_predict_labels)
    plot_cm_configs = plot_config.get("cm_config",{})
    plot_cm(results_df, outputpath, plot_config=plot_cm_configs)
  elif typeplt == "CMP" and labels_file != None:
    results_df = pd.read_csv(true_predict_labels)
    plot_cm_configs = plot_config.get("cm_config",{})
    plot_cm(results_df, outputpath, plotphotons=True, plot_config=plot_cm_configs)