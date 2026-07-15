import numpy as np
import os
import sys
import math
import pprint
import yaml
import copy
import ROOT
# import matplotlib.pyplot as plt
from ROOT import TFile
# No necesitamos TTree, TH1F, TH2F explícitamente, solo ROOT.TTree & co.

from modules import myutils  # mismo sistema de config que el código original
from modules import pi0Reco
from modules.plotting import plot_sigma_vs_energy_root
from modules import optimalVariabRho
import argparse
# ---------------------------------------------------------------------
# Función auxiliar: escribir recursivamente todos los histogramas ROOT
# contenidos en diccionarios anidados.
# ---------------------------------------------------------------------


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "-i", "--tree-file", type=str, required=True, help="Input ROOT file with outtree_*"
    )
    argparser.add_argument(
        "-o", "--output-dir", type=str, default="MLPolResults/datasets", help="Directory to save output datasets"
    )
    argparser.add_argument("-p","--pol", default=1, type=int, help="Polarization to use (1, -1 or 0 for SM)")
    argparser.add_argument("--gen", action="store_true", help="Use generator-level variables instead of reconstructed")
    argparser.add_argument("--use-pion", action="store_true", help="Use pion variables instead of meson variables")
    argparser.add_argument("--suffix", type=str, default="", help="Suffix to add to output file names")
    args = argparser.parse_args()
    input_root = args.tree_file
    output_dir = args.output_dir
    pol = args.pol
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # input_root_path = os.path.dirname(input_root)
    if not os.path.isfile(input_root):
        print("Input ROOT file %s not found", input_root)
        sys.exit(1)

    print("Using input tree file: %s", input_root)

    # -----------------------------------------------------------------
    # Configuración de histogramas (misma que el análisis original)
    # -----------------------------------------------------------------

    # Abrir fichero ROOT de entrada
    infile = TFile.Open(input_root, "READ")
    if not infile or infile.IsZombie():
        print("Could not open input ROOT file %s", input_root)
        sys.exit(1)


    tree_name = "outtree_original"
    tree = infile.Get(tree_name)
    if not tree:
      print(f"Tree '{tree_name}' not found in file {input_root}")
      sys.exit(1)
    n_entries = tree.GetEntries()
    # 14 variables: rho(4) + photon1(3) + photon2(3) + lepton(4)
    train_data_array = np.zeros((n_entries, 14))
    n_entries_decay = 0
    for i in range(n_entries):
        tree.GetEntry(i)
        entry = tree

        try:
            if not args.gen:
                recoTauID = int(entry.recoTauID)
                if recoTauID != 2:
                    continue
                n_entries_decay += 1
                if args.use_pion:
                    recoMesonE     = float(entry.recoPionE)
                    recoMesonTheta = float(entry.recoPionTheta)
                    recoMesonPhi   = float(entry.recoPionPhi)
                    recoMesonP     = float(entry.recoPionP)
                else:
                    recoMesonE     = float(entry.recoMesonE)
                    recoMesonTheta = float(entry.recoMesonTheta)
                    recoMesonPhi   = float(entry.recoMesonPhi)
                    recoMesonP     = float(entry.recoMesonP)

                try:
                    reco_photons_phi   = [float(entry.reco_photons_phi[j])   for j in range(2)]
                    reco_photons_theta = [float(entry.reco_photons_theta[j]) for j in range(2)]
                    reco_photons_E     = [float(entry.reco_photons_E[j])     for j in range(2)]
                except Exception as e:
                    print(f"Error extracting photon variables for entry {i}: {e}")
                    sys.exit(1)

                # sort photons by E (highest first)
                idx_sort = np.argsort(reco_photons_E)[::-1]

                lepE     = float(entry.lepE)
                lepTheta = float(entry.lepTheta)
                lepPhi   = float(entry.lepPhi)
                lepP     = float(entry.lepP)

                train_data_array[n_entries_decay - 1,  0] = recoMesonE
                train_data_array[n_entries_decay - 1,  1] = recoMesonTheta
                train_data_array[n_entries_decay - 1,  2] = recoMesonPhi
                train_data_array[n_entries_decay - 1,  3] = recoMesonP
                train_data_array[n_entries_decay - 1,  4] = reco_photons_E[idx_sort[0]]
                train_data_array[n_entries_decay - 1,  5] = reco_photons_theta[idx_sort[0]]
                train_data_array[n_entries_decay - 1,  6] = reco_photons_phi[idx_sort[0]]
                train_data_array[n_entries_decay - 1,  7] = reco_photons_E[idx_sort[1]]
                train_data_array[n_entries_decay - 1,  8] = reco_photons_theta[idx_sort[1]]
                train_data_array[n_entries_decay - 1,  9] = reco_photons_phi[idx_sort[1]]
                train_data_array[n_entries_decay - 1, 10] = lepE
                train_data_array[n_entries_decay - 1, 11] = lepTheta
                train_data_array[n_entries_decay - 1, 12] = lepPhi
                train_data_array[n_entries_decay - 1, 13] = lepP
            else:
                gen_ID = int(entry.genTauID)
                if gen_ID != 1:
                    continue
                n_entries_decay += 1
                if args.use_pion:
                    genMesonE     = float(entry.genPionE)
                    genMesonTheta = float(entry.genPionTheta)
                    genMesonPhi   = float(entry.genPionPhi)
                    genMesonP     = float(entry.genPionP)
                else:
                    genMesonE     = float(entry.genMesonE)
                    genMesonTheta = float(entry.genMesonTheta)
                    genMesonPhi   = float(entry.genMesonPhi)
                    genMesonP     = float(entry.genMesonP)

                try:
                    gen_photons_phi   = [float(entry.gen_photons_phi[j])   for j in range(2)]
                    gen_photons_theta = [float(entry.gen_photons_theta[j]) for j in range(2)]
                    gen_photons_E     = [float(entry.gen_photons_E[j])     for j in range(2)]
                except Exception as e:
                    print(f"Error extracting photon variables for entry {i}: {e}")
                    sys.exit(1)

                # sort photons by E (highest first)
                idx_sort = np.argsort(gen_photons_E)[::-1]

                lepE     = float(entry.genLepE)
                lepTheta = float(entry.genLepTheta)
                lepPhi   = float(entry.genLepPhi)
                lepP     = float(entry.genLepP)

                train_data_array[i,  0] = genMesonE
                train_data_array[i,  1] = genMesonTheta
                train_data_array[i,  2] = genMesonPhi
                train_data_array[i,  3] = genMesonP
                train_data_array[i,  4] = gen_photons_E[idx_sort[0]]
                train_data_array[i,  5] = gen_photons_theta[idx_sort[0]]
                train_data_array[i,  6] = gen_photons_phi[idx_sort[0]]
                train_data_array[i,  7] = gen_photons_E[idx_sort[1]]
                train_data_array[i,  8] = gen_photons_theta[idx_sort[1]]
                train_data_array[i,  9] = gen_photons_phi[idx_sort[1]]
                train_data_array[i, 10] = lepE
                train_data_array[i, 11] = lepTheta
                train_data_array[i, 12] = lepPhi
                train_data_array[i, 13] = lepP
        except Exception as e:
            print(f"Error processing entry {i}: {e}")
            sys.exit(1)
    print(f"Total entries in tree: {n_entries}, entries with recoTauID=2: {n_entries_decay}")
    train_data_array = train_data_array[:n_entries_decay]

    prefix = "gen" if args.gen else "reco"
    column_names = np.array([
        f"{prefix}MesonE", f"{prefix}MesonTheta", f"{prefix}MesonPhi", f"{prefix}MesonP",
        "photon1_E", "photon1_theta", "photon1_phi",
        "photon2_E", "photon2_theta", "photon2_phi",
        "lepE", "lepTheta", "lepPhi", "lepP",
    ])

    str_out = f"train_data_pol{pol}_{prefix}_{args.suffix}.npz"
    output_file = os.path.join(output_dir, str_out)
    np.savez(output_file, data=train_data_array, columns=column_names)
    print(f"Saved training data array to {output_file}")
    
if __name__ == "__main__":   main()