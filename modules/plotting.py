import ROOT
ROOT.gROOT.SetBatch(True)
import numpy as np
import ctypes
import os
os.environ["XDG_CACHE_HOME"] = "/nfs/cms/arqolmo/ExamplesFCCFullSim/tmp/mplconfig"
os.makedirs(os.environ["XDG_CACHE_HOME"], exist_ok=True)

os.environ["MPLCONFIGDIR"] = "/nfs/cms/arqolmo/ExamplesFCCFullSim/tmp/mplconfig"
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import pandas as pd
import os
from array import array
# from modules.myutils import compute_sigma_from_hist

PI0 = "π⁰"
PI = "π"
MU = "μ"
E = "e"
N = "n"
NEUTRINO = "ν"
TAU = "τ"
GAMMA = "γ"

RPI0 = "\\pi ^{0}"
RPI = "\\pi"
RMU = "\\mu"
RE = "e"
RN = "n"
RNEUTRINO = "\\nu"
RTAU = "\\tau"
RGAMMA = "\\gamma"
RRHO = "\\rho"  
RA1 = "a_{1}"
import ROOT
import numpy as np




# ──────────────────────────────────────────────────────────────────────────────
# Funciones de plotting
# ──────────────────────────────────────────────────────────────────────────────

def plot_energy_ratio(df, level_label, output_dir, pdg_name):
    """Histogramas de fracción de energía ECAL y HCAL por tipo de partícula."""
    if df.empty or "E_ecal" not in df.columns or "E_hcal" not in df.columns:
        return

    pids = sorted(df["pid"].unique())
    n = len(pids)
    if n == 0:
        return

    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)

    for idx, pid in enumerate(pids):
        ax = axes[idx // ncols][idx % ncols]
        subset = df[df["pid"] == pid].copy()
        
        # Calcular energía total de calorímetros y fracciones
        E_total = subset["E_ecal"] + subset["E_hcal"]
        mask = E_total > 0  # Evitar división por cero
        
        frac_ecal = (subset.loc[mask, "E_ecal"] / E_total[mask]).dropna()
        frac_hcal = (subset.loc[mask, "E_hcal"] / E_total[mask]).dropna()
        
        if frac_ecal.empty and frac_hcal.empty:
            ax.set_title(f"{pdg_name(pid)} (no data)")
            ax.set_xlim(0, 1)
            continue
        
        # Histogramas superpuestos
        bins = np.linspace(0, 1, 51)
        if not frac_ecal.empty:
            ax.hist(frac_ecal, bins=bins, alpha=0.6, label="ECAL", 
                    color="tab:blue", edgecolor="black", linewidth=0.5)
        if not frac_hcal.empty:
            ax.hist(frac_hcal, bins=bins, alpha=0.6, label="HCAL", 
                    color="tab:orange", edgecolor="black", linewidth=0.5)
        
        ax.set_xlim(0, 1)
        ax.set_xlabel("E / E_total")
        ax.set_ylabel("Particles")
        ax.set_title(f"{pdg_name(pid)} (N={mask.sum()})")
        ax.legend(fontsize=8)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(f"ECAL / HCAL fraction — {level_label}", fontsize=14, y=1.02)
    fig.tight_layout()
    fname = os.path.join(output_dir, f"energy_fraction_{level_label.lower().replace(' ', '_')}.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {fname}")


def plot_hit_distributions(df, level_label, output_dir, pdg_name):
    """Histogramas del nº de hits por subdetector para cada tipo de partícula."""
    if df.empty:
        return
    colors = {"n_track": "tab:green", "n_ecal": "tab:blue", "n_hcal": "tab:orange", "n_muon": "tab:purple"}
    hit_cols = [c for c in df.columns if c.startswith("n_") and c not in ["n_total"]]
    if not hit_cols:
        return

    pids = sorted(df["pid"].unique())
    n = len(pids)
    if n == 0:
        return

    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

    for idx, pid in enumerate(pids):
        ax = axes[idx // ncols][idx % ncols]
        subset = df[df["pid"] == pid]
        name = pdg_name(pid)

        # -------- NUEVO: binning común por PID --------
        max_global = 0
        for col in hit_cols:
            if col in subset.columns:
                vals = subset[col].dropna()
                if not vals.empty:
                    max_global = max(max_global, int(vals.max()))

        if max_global == 0:
            continue

        # bins = np.arange(0, max_global + 2, 20)
        bins = np.linspace(0, max_global + 1, 50)
        # ----------------------------------------------

        # -------- CAMBIO: acumular y apilar --------
        data = []
        labels = []
        cols = []

        for col in hit_cols:
            if col in subset.columns:
                vals = subset[col].dropna()
                if vals.empty or vals.max() == 0:
                    continue
                data.append(vals)
                labels.append(col.replace("n_", "").upper())
                cols.append(colors.get(col, "gray"))

        if data:
            ax.hist(
                data,
                bins=bins,
                stacked=True,
                alpha=0.6,
                edgecolor="black",
                linewidth=0.5,
                color=cols,
                align="left",
                label=labels,
            )
        # -------------------------------------------

        ax.set_xlabel("Nº hits")
        ax.set_ylabel("Particles")
        ax.set_yscale("log")
        ax.set_title(f"{name} (N={len(subset)})")
        ax.legend(fontsize=7)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(f"Hits per subdetector — {level_label}", fontsize=14, y=1.02)
    fig.tight_layout()
    fname = os.path.join(output_dir, f"hits_per_subdet_{level_label.lower().replace(' ', '_')}.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {fname}")



def plot_energy_by_subdet(df, level_label, output_dir, pdg_name):
    """Histograms of E_ECAL and E_HCAL by particle type."""
    if df.empty:
        return
    
    energy_cols = [c for c in ["E_ecal", "E_hcal"] if c in df.columns]
    if not energy_cols:
        return

    pids = sorted(df["pid"].unique())
    n = len(pids)
    if n == 0:
        return

    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

    colors = {"E_ecal": "tab:blue", "E_hcal": "tab:orange"}

    for idx, pid in enumerate(pids):
        ax = axes[idx // ncols][idx % ncols]
        subset = df[df["pid"] == pid]
        name = pdg_name(pid)

        # -------- NUEVO: bins comunes por PID --------
        max_global = 0.0
        for col in energy_cols:
            vals = subset[col].dropna()
            if not vals.empty:
                max_global = max(max_global, vals.max())

        if max_global == 0:
            continue

        bins = np.linspace(0, max_global, 50)
        # --------------------------------------------

        for col in energy_cols:
            vals = subset[col].dropna()
            if vals.empty:
                continue

            ax.hist(
                vals,
                bins=bins,
                alpha=0.6,
                label=col.replace("E_", "").upper(),
                color=colors.get(col, "gray"),
                edgecolor="black",
                linewidth=0.5,
            )

        ax.set_xlabel("Energy [GeV]")
        ax.set_ylabel("Particles")
        ax.set_title(f"{name} (N={len(subset)})")
        ax.legend(fontsize=8)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(f"ECAL / HCAL energy by particle — {level_label}", fontsize=14, y=1.02)
    fig.tight_layout()
    fname = os.path.join(output_dir, f"energy_ecal_hcal_{level_label.lower().replace(' ', '_')}.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {fname}")


def plot_ecal_hcal_ratio(df, level_label, output_dir, pdg_name):
    """Histograma del ratio E_ECAL / E_HCAL por tipo de partícula."""
    if df.empty or "E_ecal" not in df.columns or "E_hcal" not in df.columns:
        return

    pids = sorted(df["pid"].unique())
    n = len(pids)
    if n == 0:
        return

    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)

    for idx, pid in enumerate(pids):
        ax = axes[idx // ncols][idx % ncols]
        subset = df[df["pid"] == pid].copy()
        
        # Calcular ratio ECAL/HCAL (solo donde HCAL > 0)
        mask = subset["E_hcal"] > 0
        ratio = (subset.loc[mask, "E_ecal"] / subset.loc[mask, "E_hcal"]).dropna()
        
        if ratio.empty:
            ax.set_title(f"{pdg_name(pid)} (no data)")
            continue
        
        # Limitar ratio para visualización (valores muy altos -> ECAL dominante)
        ratio_clipped = ratio.clip(upper=10)
        
        ax.hist(ratio_clipped, bins=50, alpha=0.75, color="tab:green", 
                edgecolor="black", linewidth=0.5)
        ax.set_xlabel("E_ECAL / E_HCAL")
        ax.set_ylabel("Particles")
        ax.set_title(f"{pdg_name(pid)} (N={len(ratio)})")
        ax.axvline(x=1, color="red", linestyle="--", linewidth=1, alpha=0.7)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(f"E_ECAL / E_HCAL ratio — {level_label}", fontsize=14, y=1.02)
    fig.tight_layout()
    fname = os.path.join(output_dir, f"ratio_ecal_hcal_{level_label.lower().replace(' ', '_')}.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {fname}")


def plot_hits_ecal_vs_hcal(df, level_label, output_dir, pdg_name):
    """Histograma comparativo de número de hits ECAL vs HCAL por tipo de partícula."""
    if df.empty or "n_ecal" not in df.columns or "n_hcal" not in df.columns:
        return

    pids = sorted(df["pid"].unique())
    n = len(pids)
    if n == 0:
        return

    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

    for idx, pid in enumerate(pids):
        ax = axes[idx // ncols][idx % ncols]
        subset = df[df["pid"] == pid]
        name = pdg_name(pid)
        
        n_ecal = subset["n_ecal"].dropna()
        n_hcal = subset["n_hcal"].dropna()
        
        if n_ecal.empty and n_hcal.empty:
            ax.set_title(f"{name} (no data)")
            continue
        
        # Bins comunes
        max_val = max(n_ecal.max() if not n_ecal.empty else 0, 
                      n_hcal.max() if not n_hcal.empty else 0)
        if max_val == 0:
            continue
        bins = np.linspace(0, max_val + 1, 50)
        
        if not n_ecal.empty:
            ax.hist(n_ecal, bins=bins, alpha=0.6, label="ECAL", 
                    color="tab:blue", edgecolor="black", linewidth=0.5)
        if not n_hcal.empty:
            ax.hist(n_hcal, bins=bins, alpha=0.6, label="HCAL", 
                    color="tab:orange", edgecolor="black", linewidth=0.5)
        
        ax.set_xlabel("Nº hits")
        ax.set_ylabel("Particles")
        ax.set_yscale("log")
        ax.set_title(f"{name} (N={len(subset)})")
        ax.legend(fontsize=8)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(f"ECAL vs HCAL hits — {level_label}", fontsize=14, y=1.02)
    fig.tight_layout()
    fname = os.path.join(output_dir, f"hits_ecal_vs_hcal_{level_label.lower().replace(' ', '_')}.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {fname}")


def plot_sigma_vs_energy_root(hist_dict, bin_edges, energy_bins, region_order=None):
    """
    Compute sigma for each (region, energy bin) histogram and
    create a ROOT TCanvas showing sigma vs energy with individual TGraphs.

    Returns the canvas, a list of graphs, and the sigma_results.
    """

    if region_order is None:
        region_order = list(hist_dict.keys())

    sigma_results = {region: {} for region in region_order}

    canvas = ROOT.TCanvas("c_sigma", "Sigma vs Energy", 900, 700)

    # Lista de gráficos individuales
    graphs = []

    colors = {
        "barrel": ROOT.kRed+1,
        "endcap": ROOT.kBlue+1,
        "transition": ROOT.kGreen+2
    }

    first_graph_drawn = False  # para la primera llamada a Draw()

    for region in region_order:
        E_centers = []
        E_errors = []
        sigmas = []
        sigma_errors = []

        last_bin = len(energy_bins) - 1

        for i, (E_min, E_max) in enumerate(energy_bins):

            if (i != last_bin) and (not np.isfinite(E_max)):
                continue

            E_max_finite = E_max if np.isfinite(E_max) else 25.0
            E_center = 0.5 * (E_min + E_max_finite)
            E_err = 0.5 * (E_max_finite - E_min)

            if region not in hist_dict or i not in hist_dict[region]:
                continue

            values_list = hist_dict[region][i]
            hist, edges = np.histogram(values_list, bins=50, range=(-0.2, 0.2))

            # print(f"Region: {region}, Energy bin: {E_min}-{E_max_finite}, Number of entries: {len(values_list)}")
            plt.hist(values_list, bins=50)
            plt.title(f"Region: {region}, Energy bin: {E_min}-{E_max_finite}")
            plt.savefig(f"hist_region_{region}_energy_{E_min}_{E_max_finite}.png")
            plt.close()

            sigma, sigma_err = compute_sigma_from_hist(hist, edges, values_list)
            print(f"Computed sigma: {sigma} ± {sigma_err}")
            if np.isnan(sigma):
                continue

            E_centers.append(E_center)
            E_errors.append(E_err)
            sigmas.append(sigma)
            sigma_errors.append(sigma_err)

            sigma_results[region][i] = {
                "E_center": E_center,
                "E_err": E_err,
                "sigma": sigma,
                "sigma_err": sigma_err
            }

        if len(E_centers) == 0 or np.all(np.isnan(sigmas)):
            continue

        n = len(E_centers)
        x = np.array(E_centers, dtype='float64')
        y = np.array(sigmas, dtype='float64')
        ex = np.array(E_errors, dtype='float64')
        ey = np.array(sigma_errors, dtype='float64')

        graph = ROOT.TGraphErrors(n, x, y, ex, ey)
        graph.SetLineColor(colors.get(region, ROOT.kBlack))
        graph.SetMarkerColor(colors.get(region, ROOT.kBlack))
        graph.SetMarkerStyle(20)
        graph.SetLineWidth(2)
        graph.SetTitle(region)

        graphs.append((graph, region))

        # Dibujar cada gráfico individualmente
        if not first_graph_drawn:
            graph.Draw("AP")
            first_graph_drawn = True
        else:
            graph.Draw("P SAME")

    if not graphs:
        return None, [], {}

    first_graph = graphs[0][0]
    first_graph.GetXaxis().SetTitle("Energy [units]")
    first_graph.GetYaxis().SetTitle("Sigma of resolution")
    first_graph.GetXaxis().SetLimits(0, 26)
    first_graph.SetMinimum(0.0)
    first_graph.SetMaximum(0.5)

    # Crear leyenda
    legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    for g, region in graphs:
        legend.AddEntry(g, region, "p")

    # IMPORTANTE: Dibujar ANTES de Modified/Update y añadir al canvas
    legend.Draw()
    canvas.GetListOfPrimitives().Add(legend)

    canvas.Modified()
    canvas.Update()

    return canvas, graphs, sigma_results

def id_to_key_root(event_id, photons=False):
  if photons:
    if event_id < 0:
      if event_id == -13:
        key = f"{RMU}"
      elif event_id == -11:
        key = f"{RE}"
      elif event_id <= -20:
        key = f"h{RN}"
      elif event_id == -1:
        key = "Unknown"
      elif event_id == -2:
        key = "Unmatched"
      else:
        key = "Unknown ID"      
    elif event_id == 0:
      key = f"h"
    elif event_id < 10:
      key = f"h{event_id}{RGAMMA}"
    elif event_id == 10:
      key = f"3h"
    else:
      key = f"3h{event_id-10}{RGAMMA}"

  else:
    if event_id < 0:
      if event_id == -13:
        key = f"{RTAU} \\rightarrow {RMU}2{RNEUTRINO}"
      elif event_id == -11:
        key = f"{RTAU} \\rightarrow {RE}2{RNEUTRINO}"
      elif event_id <= -20:
        key = f"{RPI}{RN}"
      elif event_id == -1:
        key = "Unknown"
      elif event_id == -2:
        key = "Unmatched"
      else:
        key = "Unknown ID"  
    elif event_id == 0:
      key = f"{RTAU} \\rightarrow {RPI}{RNEUTRINO}"
    elif event_id < 10:
      if event_id == 1:
        key = f"{RRHO} \\rightarrow {RPI}{RPI0}{RNEUTRINO}"
      elif event_id == 2:
        key = f"{RA1} \\rightarrow {RPI}{event_id}{RPI0}{RNEUTRINO}"
      else :
        key = f"{RTAU} \\rightarrow {RPI}{event_id}{RPI0}{RNEUTRINO}"
    elif event_id == 10:
      key = f"{RA1} \\rightarrow 3{RPI}{RNEUTRINO}"
    else:
      key = f"{RA1} \\rightarrow {3}{RPI}{event_id-10}{RPI0}{RNEUTRINO}"
  return key


def id_to_key(event_id, photons=False):
  if photons:
    if event_id < 0:
      if event_id == -13:
        key = f"{MU}"
      elif event_id == -11:
        key = f"{E}"
      elif event_id <= -20:
        key = f"h{N}"
      elif event_id == -1:
        key = "Unknown"
      elif event_id == -2:
        key = "Unmatched"
      else:
        key = "Unknown ID"  
    elif event_id == 0:
      key = f"h"
    elif event_id == 1:
      key = f"h{GAMMA}"
    elif event_id < 10:
      key = f"h{event_id}{GAMMA}"
    elif event_id == 10:
      key = f"3h"
    else:
      key = f"3h{event_id-10}{GAMMA}"

  else:
    if event_id < 0:
      if event_id == -13:
        key = f"{TAU} → {MU}2{NEUTRINO}"
      elif event_id == -11:
        key = f"{TAU} → {E}2{NEUTRINO}"
      elif event_id <= -20:
        key = f"{PI}{N}"
      elif event_id == -1:
        key = "Unknown"
      elif event_id == -2:
        key = "Unmatched"
      else:
        key = "Unknown ID"  
    elif event_id == 0:
      key = f"{TAU} → {PI}{NEUTRINO}"
    elif event_id < 10:
      key = f"{TAU} → {PI}{event_id}{PI0}{NEUTRINO}"
    elif event_id == 10:
      key = f"{TAU} → 3{PI}{NEUTRINO}"
    else:
      key = f"{TAU} → {3}{PI}{event_id-10}{PI0}{NEUTRINO}"
  return key

def plot_1D_hist(file, variabs, labels, outputpath, normalize):
    """
    Plots individual 1D histograms.
    Each histogram is retrieved from the ROOT file, formatted using the provided labels,
    and then saved as a PNG file in the "1D" subfolder of outputpath.
    """
    # Create output folder for 1D plots
    out_dir = os.path.join(outputpath, "1D")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # Create a canvas for drawing
    c = ROOT.TCanvas("c_1D", "1D Histograms", 900, 700)
    # Loop over each variable
    for var in variabs:
        histo = file.Get(var)
        if not histo:
            print(f"Warning: Histogram '{var}' not found in the ROOT file.")
            continue
        # If label config exists, set axis titles and overall title
        if var in labels:
            cfg = labels[var]
            # Si hist es un ROOT.TGraphAsymmErrors
            if isinstance(histo, ROOT.TGraphAsymmErrors):
              histo.GetXaxis().SetTitle(cfg.get("x", ""))
              histo.GetYaxis().SetTitle(cfg.get("y", ""))
            else:
              histo.SetXTitle(cfg.get("x", ""))
              histo.SetYTitle(cfg.get("y", ""))
            histo.GetYaxis().SetMaxDigits(2)
            title = cfg.get("title", "")
            if title:
              # print(title)
              if "Decay" in title and "(" in title:
                # Buscamos decay y el id entre ()
                title_text = title.split("Decay")[0].split("(")[0]
                decay = title.split("(")[1].split(")")[0]
                decay_id = decay.split(" ")[1]
                decay_str = id_to_key_root(int(decay_id))
                title = "\\text{" + title_text + "}(" + decay_str+")"
            histo.SetTitle(title)
        else:
          cfg = {}
        c.Clear()
        if cfg.get("fit", False):
          ROOT.gStyle.SetOptFit(1111)
          

        if "effi" in var.lower():
          # Configuración para gráficos de eficiencia (que pueden ser TGraphAsymmErrors)
          histo.SetMarkerStyle(20)   # círculos sólidos
          histo.SetMarkerSize(1.1)
          histo.SetMarkerColor(9)
        if normalize:
          if "effi" in var.lower():
            if isinstance(histo, ROOT.TGraphAsymmErrors):
              histo.DrawNormalized("AP")  # A: crea los ejes, P: muestra puntos, E: muestra barras de error
            else:
              histo.DrawNormalized("P")  # Solo puntos para histogramas normales de eficiencia
          else:
            histo.DrawNormalized("HIST")  # Para histogramas regulares, usar HIST
        else:
          if "effi" in var.lower():
            if isinstance(histo, ROOT.TGraphAsymmErrors):
              histo.Draw("AP")  # A: crea los ejes, P: muestra puntos, E: muestra barras de error
            else:
              histo.Draw("P")  # Solo puntos para histogramas normales de eficiencia
            histo.GetYaxis().SetRangeUser(0.0, 1.1)
            histo.GetYaxis().SetNoExponent(ROOT.kTRUE)
          else:
            histo.Draw("HIST")
        if cfg.get("fit", False):
          histo.SetLineColor(ROOT.kBlack)
          # Definir el ajuste Crystal Ball en el rango deseado
          # cb = ROOT.TF1("cb", "crystalball", 0.05, 0.25)

          # Parámetros: mean, sigma, alpha, n
          # cb.SetParameters(0.135, 0.01, 1.5, 5)

          # Ajustar histograma
          # fitres = histo.Fit(cb, "R")
          
          if cfg.get("fitrange", None):
            fit_range = cfg["fitrange"]
            fitres = histo.Fit("gaus", "Q","", fit_range[0], fit_range[1])  # devuelve TFitResultPtr
          else:
            fitres = histo.Fit("gaus", "Q")  # devuelve TFitResultPtr
            
          # f = histo.GetFunction("cb")
          f = histo.GetFunction("gaus")
          
          if f:
            f.SetLineWidth(2)
            # aseguramos que se vea sobre lo ya dibujado
            f.Draw("same")
            c.Update()
        out_file = os.path.join(out_dir, f"{var}.png")
        if cfg.get("print_mean", False):
          # Get absolute mean of the histogram (not from the fit)
          values = []
          for bin in range(1, histo.GetNbinsX() + 1):
              bin_center = abs(histo.GetBinCenter(bin))
              bin_content = histo.GetBinContent(bin)
              values.extend([bin_center] * int(bin_content))
          if values:
            mean = np.mean(values)
          # mean = histo.GetMean()
          print("===============================")
          print(f"Mean of histogram '{var}': {mean}")
          print("===============================")
        if cfg.get("logy", False):
          c.SetLogy()
        c.SaveAs(out_file)
        print(f"Saved histogram '{var}' as '{out_file}'")
    ROOT.gStyle.SetOptFit()
    c.Close()


def plot_hist_zoom(file, zoom_config, outputpath):
    """
    Plots groups of 1D histograms together on a single canvas.
    For each group defined in the plot_together configuration, all histograms are drawn
    on the same canvas with different colors and a legend.
    The resulting canvas is saved in a "together" subfolder inside the "1D" folder.
    """
    # Create output folder for together plots under 1D/together
    out_dir = os.path.join(outputpath, "1D", "zoom")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # Iterate over each group in the together configuration
    for group, cfg in zoom_config.items():
        c = ROOT.TCanvas(f"c_zoom_{group}", group, 900, 700)
        position = cfg.get("position", None)
        if position:
          legend = ROOT.TLegend(position[0], position[1], position[2], position[3])
        else:
          legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
        legend.SetTextSize(0.03)
        legend.SetBorderSize(1)
        legend.SetFillStyle(0)
        first = True
        normalize = cfg.get("norm", False)
        global_max_value = 0
        for histo_name in cfg["main"]:
          histo = file.Get(histo_name)
          if not histo:
            continue
          max_histo = histo.GetMaximum()
          if max_histo > global_max_value:
            global_max_value = max_histo
        
        main = cfg.get("main", None)
        # For each histogram name in the group configuration
        for i, var in enumerate(main):
            histo = file.Get(var)
            if not histo:
                print(f"Warning: Histogram '{var}' not found in the ROOT file.")
                continue
            # Set the common X and Y axis titles from the group config

            
            # Assign a distinct line color for each histogram (simple scheme)
            if i == 0:
                histo.SetLineColor(ROOT.kRed)
            elif i == 1:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.kBlue)
            elif i == 2:
                histo.SetLineStyle(3)
                histo.SetLineColor(ROOT.kGreen+2)
            elif i == 3:
                histo.SetLineColor(ROOT.kMagenta)
                histo.SetLineStyle(4)
            elif i == 4:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.kYellow)
            elif i == 5:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.Spring)
            elif i == 6:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.kAzure)
            elif i == 7:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.kViolet)
                
            histo.SetLineWidth(2)
            
            # Draw the first histogram normally; then draw others with "same"
            if first:
              histo.SetXTitle(cfg.get("x", ""))
              histo.SetYTitle(cfg.get("y", ""))
              # Optionally, set the title of the histogram (or leave it empty)
              title = cfg.get("title", "")
              if title:
                if "Decay" in title and "(" in title:
                  # Buscamos decay y el id entre ()
                  title_text = title.split("Decay")[0].split("(")[0]
                  decay = title.split("(")[1].split(")")[0]
                  decay_id = decay.split(" ")[1]
                  decay_str = id_to_key_root(int(decay_id))
                  title = "\\text{" + title_text + "}(" + decay_str+")"
              histo.SetTitle(title)
              max_val = histo.GetMaximum()
              if normalize:
                print(f"Max value before scaling: {max_val}")
                if max_val != 0:
                # Escalamos para que el máximo sea 1.
                  histo.Scale(1.0 / max_val)
                histo.GetYaxis().SetRangeUser(0., 1.1)
              else:
                histo.GetYaxis().SetRangeUser(0, 1.1 * global_max_value)
              histo.Draw("HIST")
              first = False
            else:
              max_val = histo.GetMaximum()
              if normalize:
                if max_val != 0:
                # Escalamos para que el máximo sea 1.
                  histo.Scale(1.0 / max_val)
                # histo_norm.GetYaxis().SetRangeUser(0, 1)
                

              histo.Draw("HIST same")

                # Set axis from 0 to 1
                
            # Add legend entry using the corresponding label from config, if provided
            label = cfg["labels"][i] if i < len(cfg["labels"]) else var
            legend.AddEntry(histo, label, "l")
        
        legend.Draw()
        first = True
      
        zoom_pos = cfg.get("zoom_position", [0.5, 0.6, 0.9, 0.9])
        zoom = ROOT.TPad("zoom", "zoom", zoom_pos[0], zoom_pos[1], zoom_pos[2], zoom_pos[3])
        zoom.Draw()
        zoom.cd()
        zoom_plots = cfg.get("zoom", None)
        for i, var in enumerate(zoom_plots):
            histo = file.Get(var)
            if not histo:
                print(f"Warning: Histogram '{var}' not found in the ROOT file.")
                continue
            # Set the common X and Y axis titles from the group config

            
            # Assign a distinct line color for each histogram (simple scheme)
            if i == 0:
                histo.SetLineColor(ROOT.kRed)
            elif i == 1:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.kBlue)
            elif i == 2:
                histo.SetLineStyle(3)
                histo.SetLineColor(ROOT.kGreen+2)
            elif i == 3:
                histo.SetLineColor(ROOT.kMagenta)
                histo.SetLineStyle(4)
            elif i == 4:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.kYellow)
            elif i == 5:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.Spring)
            elif i == 6:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.kAzure)
            elif i == 7:
                histo.SetLineStyle(2)
                histo.SetLineColor(ROOT.kViolet)
                
            histo.SetLineWidth(2)
            
            # Draw the first histogram normally; then draw others with "same"
            if first:
              max_val = histo.GetMaximum()
              if normalize:
                print(f"Max value before scaling: {max_val}")
                if max_val != 0:
                # Escalamos para que el máximo sea 1.
                  histo.Scale(1.0 / max_val)
                histo.GetYaxis().SetRangeUser(0., 1.1)
              else:
                histo.GetYaxis().SetRangeUser(0, 1.1 * global_max_value)
              histo.SetTitle("")
              histo.Draw("HIST")
              first = False
            else:
              max_val = histo.GetMaximum()
              if normalize:
                if max_val != 0:
                # Escalamos para que el máximo sea 1.
                  histo.Scale(1.0 / max_val)
                # histo_norm.GetYaxis().SetRangeUser(0, 1)
                

              histo.Draw("HIST same")
          
        
        
        out_file = os.path.join(out_dir, f"{group}.png")
        c.SaveAs(out_file)
        print(f"Saved group '{group}' as '{out_file}'")
        c.Close()

def get_graph_max(graph):
    """
    Returns Y maxium in a TGraphAsymmErrors.
    If input is not TGraphAsymmErrors, the function returns GetMaximum().
    """
    if isinstance(graph, ROOT.TGraphAsymmErrors):
        n_points = graph.GetN()
        if n_points == 0:
            return -1111.0
        y_vals = [graph.GetPointY(i) for i in range(n_points)]
        return max(y_vals)
    else:
        return graph.GetMaximum()

def plot_hist_together(file, together_config, outputpath):
    """
    Plots groups of 1D histograms together on a single canvas.
    For each group defined in the plot_together configuration, all histograms are drawn
    on the same canvas with different colors and a legend.
    The resulting canvas is saved in a "together" subfolder inside the "1D" folder.
    """
    # Create output folder for together plots under 1D/together
    out_dir = os.path.join(outputpath, "1D", "together")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # Iterate over each group in the together configuration
    for group, cfg in together_config.items():
        if "effi" in group.lower():
          scatter_plot = True
        else:
          scatter_plot = False
        c = ROOT.TCanvas(f"c_together_{group}", group, 900, 700)
        c.SetLeftMargin(0.12)   # margen izquierdo (0.1-0.2 suele ser bueno)
        c.SetRightMargin(0.05)  # margen derecho
        c.SetBottomMargin(0.12) # margen inferior
        c.SetTopMargin(0.08)    # margen superior
        position = cfg.get("position", None)
        if position:
          legend = ROOT.TLegend(position[0], position[1], position[2], position[3])
        else:
          legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
        
        ncolumns = cfg.get("ncolumns", 1)
        legend.SetNColumns(ncolumns)
        
        legend_text_size = cfg.get("legend_txt_size", 0.03)
        legend.SetTextSize(legend_text_size)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        first = True
        
        draw_mode = cfg.get("draw_mode", "hist")
        normalize = cfg.get("norm", False)
        global_max_value = 0
        rebin = cfg.get("rebin", None)
        # print(global_max_value)
        # print("===============================")
        for histo_name in cfg["variabs"]:
          # original_histo = file.Get(histo_name)
          histo = file.Get(histo_name)
          if not histo:
            continue
          # histo = original_histo.Clone()
          
          # Apply rebin if specified
          if rebin:
            histo.Rebin(rebin)
          # print(histo)
          max_histo = get_graph_max(histo)
          # print(max_histo)
          # print("===============================")
          if max_histo > global_max_value:
            global_max_value = max_histo
        
        # For each histogram name in the group configuration
        linewidth = cfg.get("linewidth", 2)
        markersize = cfg.get("markersize", 1.2)
        drawn_histos = {}
        fill_graphs = []
        if cfg.get("gradient", False):
            # seet gStylePalette
            ROOT.gStyle.SetPalette(ROOT.kSolar)
            
            # n_vars = len(cfg["variabs"])
            # gradient_colors = [ROOT.TColor.GetColorGradient(i / max(1, n_vars - 1), 3) for i in range(n_vars)]
        for i, var in enumerate(cfg["variabs"]):
            histo = file.Get(var)
            if not histo:
                print(f"Warning: Histogram '{var}' not found in the ROOT file.")
                continue
            # Set the common X and Y axis titles from the group config

            if not cfg.get("gradient", False):
            # Assign a distinct line color for each histogram (simple scheme)
                if i == 0:
                  # histo.SetLineWidth(linewidth)
                  if scatter_plot:
                    histo.SetMarkerStyle(20)
                  if cfg.get("fill", False):
                    histo.SetFillColor(ROOT.kBlue)
                  histo.SetLineColor(ROOT.kBlue)
                  histo.SetMarkerColor(ROOT.kBlue)
                elif i == 1:
                  if not scatter_plot:
                    histo.SetLineStyle(1)
                  else:
                    histo.SetMarkerStyle(20)   # círculos sólidos
                  if cfg.get("fill", False):
                    histo.SetFillColor(ROOT.kRed)
                  histo.SetLineColor(ROOT.kRed)
                  histo.SetMarkerColor(ROOT.kRed)
                  # histo.SetLineWidth(linewidth)
                elif i == 2:
                  if not scatter_plot:
                    histo.SetLineStyle(1)
                  else:
                    histo.SetMarkerStyle(20)   # círculos sólidos
                  if cfg.get("fill", False):
                    histo.SetFillColor(ROOT.kGreen+2)
                  histo.SetLineColor(ROOT.kGreen+2)
                  histo.SetMarkerColor(ROOT.kGreen+2)
                  # histo.SetLineWidth(linewidth)
                elif i == 3:
                    if not scatter_plot:
                        histo.SetLineStyle(2)
                    else:
                        histo.SetMarkerStyle(21)
                    # histo.SetLineColor(ROOT.kOrange+7)
                    histo.SetLineColor(ROOT.kOrange+7)
                    histo.SetMarkerColor(ROOT.kOrange+7)
                    if cfg.get("fill", False):
                        histo.SetFillColor(ROOT.kOrange-3)

                elif i == 4:
                  if not scatter_plot:
                      histo.SetLineStyle(2)
                  else:
                      histo.SetMarkerStyle(22)
                  histo.SetLineColor(ROOT.kViolet+1)
                  # histo.SetLineColor(ROOT.kViolet+1)
                  histo.SetMarkerColor(ROOT.kViolet+1)
                  if cfg.get("fill", False):
                      histo.SetFillColor(ROOT.kViolet-5)

                elif i == 5:
                  if not scatter_plot:
                      histo.SetLineStyle(2)
                  else:
                      histo.SetMarkerStyle(23)
                  histo.SetLineColor(ROOT.kTeal+2) 
                  # histo.SetLineColor(ROOT.kTeal+2)
                  histo.SetMarkerColor(ROOT.kTeal+2)
                  if cfg.get("fill", False):
                      histo.SetFillColor(ROOT.kTeal-5)

                elif i == 6:
                  if not scatter_plot:
                      histo.SetLineStyle(9)
                  else:
                      histo.SetMarkerStyle(24)
                  histo.SetLineColor(ROOT.kPink+9)   # rojo vino oscuro
                  histo.SetMarkerColor(ROOT.kPink+9)
                  if cfg.get("fill", False):
                      histo.SetFillColor(ROOT.kPink-2)
                elif i == 7:
                    if not scatter_plot:
                        histo.SetLineStyle(3)
                    else:
                        histo.SetMarkerStyle(25)
                    histo.SetLineColor(ROOT.kBlue+3)
                    histo.SetMarkerColor(ROOT.kBlue+3)
                    if cfg.get("fill", False):
                        histo.SetFillColor(ROOT.kBlue-7)
            # else:
                # histo.SetLineColor(gradient_colors[i])
                # histo.SetMarkerColor(gradient_colors[i])
                # if cfg.get("fill", False):
                
                    # histo.SetFillColor(gradient_colors[i])
            histo.SetLineWidth(linewidth)
            histo.SetMarkerSize(markersize)
            if cfg.get("fill", False):
              histo.SetFillStyle(3002)
            # Draw the first histogram normally; then draw others with "same"
            drawn_histos[var] = histo
            if first:
              if isinstance(histo, ROOT.TGraphAsymmErrors):
                histo.GetXaxis().SetTitle(cfg.get("x", ""))
                histo.GetYaxis().SetTitle(cfg.get("y", ""))
              else:
                histo.SetXTitle(cfg.get("x", ""))
                histo.SetYTitle(cfg.get("y", ""))
              # Optionally, set the title of the histogram (or leave it empty)
              title = cfg.get("title", "")
              if title:
                if "Decay" in title and "(" in title:
                  # Buscamos decay y el id entre ()
                  title_text = title.split("Decay")[0].split("(")[0]
                  decay = title.split("(")[1].split(")")[0]
                  decay_id = decay.split(" ")[1]
                  decay_str = id_to_key_root(int(decay_id))
                  title = "\\text{" + title_text + "}(" + decay_str+")"
              histo.SetTitle(title)
              max_val = histo.GetMaximum()
              if normalize:
                print(f"Max value before scaling: {max_val}")
                if max_val != 0:
                # Escalamos para que el máximo sea 1.
                  histo.Scale(1.0 / max_val)
                if cfg.get("logy", False):
                  histo.GetYaxis().SetRangeUser(0.001, 1.1)
                else:
                  histo.GetYaxis().SetRangeUser(0., 1.1)
              else:
                if cfg.get("logy", False):
                  histo.GetYaxis().SetRangeUser(0.001, 1.1 * global_max_value)
                else:
                  histo.GetYaxis().SetRangeUser(0, 1.1 * global_max_value)
              
              if scatter_plot:
                # histo.Sumw2()
                if isinstance(histo, ROOT.TGraphAsymmErrors):
                    draw_opt = "AP PMC PLC" if cfg.get("gradient", False) else "AP"
                    if cfg.get("fill", False) and cfg.get("gradient", False):
                        draw_opt += "PFC"
                #   histo.Draw("AP")
                else:
                    draw_opt = "P PMC" if cfg.get("gradient", False) else "P"
                    # histo.Draw("P")
                # histo.Draw("L SAME")
              else:
                if draw_mode == "points":
                    histo.SetMarkerStyle(20)
                    draw_opt = "P PMC" if cfg.get("gradient", False) else "P"

                elif draw_mode == "line":
                    draw_opt = "L PLC" if cfg.get("gradient", False) else "L"

                else:  # hist por defecto
                    draw_opt = "HIST PLC" if cfg.get("gradient", False) else "HIST"
                    if cfg.get("fill", False) and cfg.get("gradient", False):
                        draw_opt += " PFC"
              histo.Draw(draw_opt)
              first = False
            else:
              max_val = histo.GetMaximum()
              if normalize:
                if max_val != 0:
                # Escalamos para que el máximo sea 1.
                  histo.Scale(1.0 / max_val)
                # histo_norm.GetYaxis().SetRangeUser(0, 1)
                
              if scatter_plot:
                # histo.Sumw2()
                draw_opt = "P SAME PMC PLC" if cfg.get("gradient", False) else "P SAME"
                if cfg.get("fill", False) and cfg.get("gradient", False):
                    draw_opt += "PFC"
                # histo.Draw("P same")
                # histo.Draw("L same")
              else:
                if draw_mode == "points":
                    histo.SetMarkerStyle(20)
                    draw_opt = "P SAME PMC" if cfg.get("gradient", False) else "P SAME"
                elif draw_mode == "line":
                    draw_opt = "L SAME PLC" if cfg.get("gradient", False) else "L SAME"
                else:
                    draw_opt = "HIST SAME PLC" if cfg.get("gradient", False) else "HIST SAME"
                    if cfg.get("fill", False) and cfg.get("gradient", False):
                        draw_opt += " PFC"
              histo.Draw(draw_opt)

                # Set axis from 0 to 1
                
            # Add legend entry using the corresponding label from config, if provided
            label = cfg["labels"][i] if i < len(cfg["labels"]) else var
            legend.AddEntry(histo, label, "l")
        
        # ===== Fill-between functionality =====
        if "fill_between" in cfg:
            fb = cfg["fill_between"] 

            reference_name = fb["reference"]
            targets = fb["histos"]
            fill_style = fb.get("fillstyle", 3004)
            fill_colors = fb.get("colors", ROOT.kGray+1)

            if reference_name not in drawn_histos:
                print(f"[WARNING] reference '{reference_name}' not found for fill_between")
            else:
                href = drawn_histos[reference_name]

                # Comprobamos que referencia y targets sean TH1 (solo tiene sentido para histos)
                if not isinstance(href, ROOT.TH1):
                    print(f"[WARNING] reference '{reference_name}' is not a TH1, skip polygon fill_between")
                else:
                    nb = href.GetNbinsX()

                    for k, hname in enumerate(targets):
                        if hname not in drawn_histos:
                            print(f"[WARNING] target '{hname}' not found for fill_between")
                            continue
                        htar = drawn_histos[hname]
                        if not isinstance(htar, ROOT.TH1):
                            print(f"[WARNING] target '{hname}' is not a TH1, skip in polygon fill_between")
                            continue

                        fill_color = (
                            fill_colors[k]
                            if isinstance(fill_colors, list)
                            else fill_colors
                        )

                        # --- Modo HIST: polígono siguiendo bordes de bin ---
                        if draw_mode.lower() == "hist":
                            # 4 puntos por bin: 2 para borde superior, 2 para borde inferior
                            g = ROOT.TGraph(4 * nb)
                            g.SetName(f"fill_between_{reference_name}_{hname}")
                            idx = 0

                            # Parte superior: x_low -> x_high manteniendo y_up constante
                            for b in range(1, nb + 1):
                                x_low = href.GetBinLowEdge(b)
                                x_high = x_low + href.GetBinWidth(b)
                                yref = href.GetBinContent(b)
                                ytar = htar.GetBinContent(b)
                                y_up = max(yref, ytar)
                                # punto en x_low
                                g.SetPoint(idx, x_low, y_up)
                                idx += 1
                                # punto en x_high
                                g.SetPoint(idx, x_high, y_up)
                                idx += 1

                            # Parte inferior: x_high -> x_low (orden inverso)
                            for b in range(nb, 0, -1):
                                x_low = href.GetBinLowEdge(b)
                                x_high = x_low + href.GetBinWidth(b)
                                yref = href.GetBinContent(b)
                                ytar = htar.GetBinContent(b)
                                y_low = min(yref, ytar)
                                # punto en x_high
                                g.SetPoint(idx, x_high, y_low)
                                idx += 1
                                # punto en x_low
                                g.SetPoint(idx, x_low, y_low)
                                idx += 1

                        # --- Otros modos (line/points): comportamiento antiguo por centros ---
                        else:
                            g = ROOT.TGraph(2 * nb)
                            g.SetName(f"fill_between_{reference_name}_{hname}")
                            idx = 0

                            # Upper boundary (por centros)
                            for b in range(1, nb + 1):
                                x = href.GetBinCenter(b)
                                yref = href.GetBinContent(b)
                                ytar = htar.GetBinContent(b)
                                g.SetPoint(idx, x, max(yref, ytar))
                                idx += 1

                            # Lower boundary (por centros, en orden inverso)
                            for b in range(nb, 0, -1):
                                x = href.GetBinCenter(b)
                                yref = href.GetBinContent(b)
                                ytar = htar.GetBinContent(b)
                                g.SetPoint(idx, x, min(yref, ytar))
                                idx += 1

                        g.SetFillColor(fill_color)
                        g.SetFillStyle(fill_style)
                        g.SetLineColor(fill_color)
                        g.Draw("F SAME")
                        fill_graphs.append(g)
        legend.Draw()
        if cfg.get("gridy", False):
          c.SetGridy()
        if cfg.get("logy", False):
          c.SetLogy()
        out_file = os.path.join(out_dir, f"{group}.png")
        c.SaveAs(out_file)
        print(f"Saved group '{group}' as '{out_file}'")
        c.Close()


def plot_2D_hist(file, variabs, labels, outputpath):
    """
    Plots 2D histograms.
    Each 2D histogram is drawn with the "COLZ" option and saved as a PNG file
    in the "2D" subfolder of outputpath.
    """
    # Create output folder for 2D plots
    out_dir = os.path.join(outputpath, "2D")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    c = ROOT.TCanvas("c_2D", "2D Histograms", 900, 700)
    for var in variabs:
        histo = file.Get(var)
        if not histo:
            print(f"Warning: 2D Histogram '{var}' not found in the ROOT file.")
            continue
        if var in labels:
            cfg = labels[var]
            histo.SetXTitle(cfg.get("x", ""))
            histo.SetYTitle(cfg.get("y", ""))
            histo.SetTitle(cfg.get("title", ""))
            if cfg.get("rebinx", None):
                histo.RebinX(cfg["rebinx"])
            if cfg.get("rebiny", None):
                histo.RebinY(cfg["rebiny"])
        c.Clear()
        histo.Draw("COLZ")
        if cfg.get("name", None):
          out_name = cfg["name"]
        else:
          out_name = var
        out_file = os.path.join(out_dir, f"{out_name}.png")
        c.SaveAs(out_file)
        print(f"Saved 2D histogram '{var}' as '{out_file}'")
    c.Close()

def plot_cm(results_df, outputpath, plotphotons=False, plot_config={}):
    """
    Generates and saves two confusion matrix plots:
      1. With absolute values.
      2. With normalized values (per actual class) expressed as percentages.
    
    If plotphotons is True, uses the 'PhotonPredicted' column instead of 'Predicted'.
    The resulting confusion matrix may be rectangular if the true and predicted classes differ.
    
    Matrices are saved in the same folder as before, with an added suffix.
    """
    # Extract true labels and predicted labels based on the flag
    y_true = results_df['True']
    if plotphotons:
        y_pred = results_df['PhotonPredicted']
        suffix = "_photons"
    else:
        y_pred = results_df['Predicted']
        suffix = ""
    
    # Compute the confusion matrix:
    if plotphotons:
        # Use pandas crosstab to allow a rectangular matrix
        cm_df = pd.crosstab(y_true, y_pred)
        cm = cm_df.values
        classes_true = cm_df.index.values
        classes_pred = cm_df.columns.values
        mapped_classes_true = [id_to_key(cls, photons=False) for cls in classes_true]
        # print(classes_pred[:10])
        mapped_classes_pred = [id_to_key(cls, photons=True) for cls in classes_pred]
        if "decays" in plot_config:
          # Select only the decays in the config file
          decays = plot_config["decays"]
          photondecays = plot_config["photonDecays"]
          cm = cm_df.loc[decays, photondecays].copy()
          classes_true = cm.index.values
          classes_pred = cm.columns.values
          cm = cm.values
          mapped_classes_true = [id_to_key(cls, photons=False) for cls in classes_true]
          mapped_classes_pred = [id_to_key(cls, photons=True) for cls in classes_pred]
        # print(mapped_classes_pred[:10])
    else:
        # Use sklearn's confusion_matrix for a square matrix
        classes = np.unique(np.concatenate((y_true.values, y_pred.values)))
        cm = confusion_matrix(y_true, y_pred, labels=classes)
        cm_df = pd.DataFrame(cm, index=classes, columns=classes)
        
        if "decays" in plot_config:
          # Select only the decays in the config file
          decays = plot_config["decays"]
          cm = cm_df.loc[decays, decays].copy()
          cm = cm.values
          classes = decays
        mapped_classes_true = [id_to_key(cls, photons=False) for cls in classes]
        mapped_classes_pred = mapped_classes_true  # Same for both axes
    
    # Create output directory if it doesn't exist
    cm_dir = os.path.join(outputpath, "CM")
    if not os.path.exists(cm_dir):
        os.makedirs(cm_dir)

    # --- Absolute values plot ---
    if "decays" in plot_config:
      plt.figure(figsize=(8, 6))
      fontsize = 10
    else:
      plt.figure(figsize=(12, 8))
      fontsize = 8
    plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
    plt.title("Confusion Matrix (Absolute Values)")
    if plot_config.get("colorbar", True):
      plt.colorbar()
    if plotphotons:
        xtick_marks = np.arange(len(mapped_classes_pred))
        ytick_marks = np.arange(len(mapped_classes_true))
        plt.xticks(xtick_marks, mapped_classes_pred, rotation=45)
        plt.yticks(ytick_marks, mapped_classes_true)
    else:
        tick_marks = np.arange(len(mapped_classes_true))
        plt.xticks(tick_marks, mapped_classes_true, rotation=45)
        plt.yticks(tick_marks, mapped_classes_true)
    
    thresh = cm.max() / 2.
    # Annotate each cell with the absolute value
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, format(cm[i, j], 'd'),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black", fontsize=fontsize)
    
    plt.ylabel('Actual Label')
    plt.xlabel('Predicted Label')
    plt.tight_layout()
    plt.savefig(os.path.join(cm_dir, "confusion_matrix_absolute" + suffix + ".png"))
    plt.close()
    
    # --- Normalized values plot (per actual label) ---    
    cm_normalized = cm_df.to_numpy().astype('float') / cm_df.to_numpy().sum(axis=1)[:, np.newaxis]
    cm_normalized = np.nan_to_num(cm_normalized)  # Replace NaN with 0 for rows with zero sum
    cm_normalized = pd.DataFrame(cm_normalized, index=cm_df.index, columns=cm_df.columns)
    if "decays" in plot_config:
      # Select only the decays in the config file
      if plotphotons:
        cm_normalized = cm_normalized.loc[decays, photondecays].copy()
        classes_true = cm_normalized.index.values
        classes_pred = cm_normalized.columns.values
        cm_normalized = cm_normalized.values

        mapped_classes_true = [id_to_key(cls, photons=False) for cls in classes_true]
        mapped_classes_pred = [id_to_key(cls, photons=True) for cls in classes_pred]
      else:
        cm_normalized = cm_normalized.loc[decays, decays].copy()
        classes_true = cm_normalized.index.values
        classes_pred = cm_normalized.columns.values
        cm_normalized = cm_normalized.values

        mapped_classes_true = [id_to_key(cls, photons=False) for cls in classes_true]
        mapped_classes_pred = [id_to_key(cls, photons=False) for cls in classes_pred]
      plt.figure(figsize=(8, 6))
      fontsize = 10
    else:
      plt.figure(figsize=(12, 8))
      fontsize = 8
    plt.imshow(cm_normalized, interpolation='nearest', cmap=plt.cm.Blues)
    plt.title("Confusion Matrix (Normalized)")
    if plot_config.get("colorbar", True):
      plt.colorbar()
    if plotphotons:
        plt.xticks(np.arange(len(mapped_classes_pred)), mapped_classes_pred, rotation=45)
        plt.yticks(np.arange(len(mapped_classes_true)), mapped_classes_true)
    else:
        tick_marks = np.arange(len(mapped_classes_true))
        plt.xticks(tick_marks, mapped_classes_true, rotation=45)
        plt.yticks(tick_marks, mapped_classes_true)
    # cm_normalized = cm_normalized.to_numpy()
    thresh_norm = cm_normalized.max() / 2.
    # Annotate each cell with the percentage
    for i in range(cm_normalized.shape[0]):
        for j in range(cm_normalized.shape[1]):
            percentage = cm_normalized[i, j] * 100
            plt.text(j, i, f"{percentage:.1f}%",
                     horizontalalignment="center",
                     color="white" if cm_normalized[i, j] > thresh_norm else "black", fontsize=fontsize)
    
    plt.ylabel('Actual Label')
    plt.xlabel('Predicted Label')
    plt.tight_layout()
    plt.savefig(os.path.join(cm_dir, "confusion_matrix_normalized" + suffix + ".png"))
    plt.close()
    
    print(f"Saved confusion matrices to '{cm_dir}'")




def plot_absolute(canvas, graphs, absolute_keys, colors, xaxis,
                           min_absolute_value, max_absolute_value,
                           title, outputpath, mig_str):
    canvas.cd()
    # Dibujar cada gráfico
    for i, key in enumerate(absolute_keys):
        color = colors[i % len(colors)]
        graph = graphs[key]
        
        graph.SetTitle(key)
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(20)
        graph.SetLineWidth(2)
        
        if i == 0:
            graph.GetXaxis().SetTitle(xaxis)
            graph.GetYaxis().SetTitle("Counts")
            graph.GetYaxis().SetRangeUser(0.9 * min_absolute_value, 1.1 * max_absolute_value)
            graph.Draw("alp")  # Dibuja ejes, línea y puntos
        else:
            graph.Draw("lp")   # Dibuja línea y puntos sin reconfigurar ejes

    canvas.BuildLegend()
    canvas.Update()

    # Agrega un recuadro de texto (TPaveText) centrado en la parte superior
    t = ROOT.TPaveText(0.2, 0.95, 0.8, 0.99, "NDC")
    t.SetTextAlign(22)  # Centrado horizontal y vertical
    t.SetTextSize(0.04)
    t.SetFillStyle(0)   # Fondo transparente
    t.SetBorderSize(0)  # Sin borde
    t.AddText(title)
    t.Draw()

    # Guarda el canvas como imagen
    canvas.SaveAs(outputpath + f"graphs_plot_{mig_str}.png")


def plot_metric(canvas, graphs, metric_keys, colors, xaxis,
                       title, outputpath, metric, mig_str):
    canvas.cd()
    for i, key in enumerate(metric_keys):
        color = colors[i % len(colors)]
        graph = graphs[key]
        
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(20)
        graph.SetLineWidth(2)
        graph.SetTitle(key)
        graph.SetName(key)
        
        if i == 0:
            graph.GetXaxis().SetTitle(xaxis)
            graph.GetYaxis().SetTitle("Normalized Counts")
            graph.GetYaxis().SetRangeUser(0, 1.2)
            graph.Draw("alp")
        else:
            graph.Draw("lp")
    
    canvas.BuildLegend()
    canvas.Update()
    
    t = ROOT.TPaveText(0.2, 0.95, 0.8, 0.99, "NDC")
    t.SetTextAlign(22)
    t.SetTextSize(0.04)
    t.SetFillStyle(0)
    t.SetBorderSize(0)
    t.AddText(title)
    t.Draw()

    canvas.Update()
    canvas.SaveAs(outputpath + f"graphs_plot_{metric}_{mig_str}.png")




def _draw_one_object_1d(obj, first, draw_as_scatter):
    """
    Dibuja un objeto ROOT 1D (TH1 o TGraph) con la opción adecuada.
    
    Si el histograma tiene fill configurado (FillStyle != 0), se fuerza
    el uso de 'HIST' para que el relleno sea visible.
    """
    if isinstance(obj, ROOT.TGraphAsymmErrors):
        opt = "AP" if first else "P same"
    else:
        # Detectar si el histograma tiene fill configurado
        has_fill = isinstance(obj, ROOT.TH1) and obj.GetFillStyle() != 0
        
        if has_fill:
            # Con fill, siempre usar HIST para que se vea el relleno
            opt = "HIST" if first else "HIST same"
        elif draw_as_scatter:
            opt = "P" if first else "P same"
        else:
            opt = "HIST" if first else "HIST same"
    obj.Draw(opt)
    return opt
def get_total_entries(obj, raw=False):
    """
    Returns total entries for TH1 or TGraph.

    raw=False (default): suma de pesos (TH1 -> Integral; TGraph -> GetN).
    raw=True:            conteo crudo SIN ponderar (TH1 -> GetEntries; TGraph -> GetN).
    Útil cuando el mismo conjunto de eventos se rellena con pesos distintos
    (p.ej. variantes de polarización): Integral difiere, GetEntries no.
    """
    if isinstance(obj, ROOT.TH1):
        return obj.GetEntries() if raw else obj.Integral()
    elif isinstance(obj, ROOT.TGraph):
        return obj.GetN()
    return 0

def save_chi2_csv(chi2_results, outdir, fig_name):
    """
    Saves Chi2Test results (chi2/ndf and p-value) to a CSV.
    chi2_results: list of dicts with keys Reference, Dataset, Chi2_ndf, pvalue, ndf
    """
    if not chi2_results:
        return
    df = pd.DataFrame(chi2_results)
    csv_path = os.path.join(outdir, f"{fig_name}_chi2.csv")
    df.to_csv(csv_path, index=False)
    print(f"[OK] Saved chi2 CSV: {csv_path}")


def save_entries_csv(entries_per_dataset, percentages, total_entries, outdir, fig_name):
    """
    Saves histogram entries data (absolute and percentage) to a CSV file.
    
    Args:
        entries_per_dataset: dict with label -> number of entries
        percentages: dict with label -> percentage of total
        total_entries: total number of entries across all datasets
        outdir: output directory path
        fig_name: name of the figure (used for CSV filename)
    """
    if not entries_per_dataset:
        return
    
    csv_data = []
    for label in entries_per_dataset.keys():
        n_entries = entries_per_dataset[label]
        pct = percentages.get(label, 0.0)
        csv_data.append({
            "Label": label,
            "Entries": int(n_entries),
            "Percentage": round(pct, 2)
        })
    
    # Add total row
    csv_data.append({
        "Label": "TOTAL",
        "Entries": int(total_entries),
        "Percentage": 100.0
    })
    
    df = pd.DataFrame(csv_data)
    csv_path = os.path.join(outdir, f"{fig_name}_entries.csv")
    df.to_csv(csv_path, index=False)
    print(f"[OK] Saved entries CSV: {csv_path}")

def _apply_style_1d(obj, color, linestyle, markerstyle, markersize, linewidth,
                    fill=False, fillstyle=3004, fillalpha=1.0):
    """
    Aplica estilo a un objeto ROOT (TH1 o TGraph).
    
    Args:
        obj: Objeto ROOT (TH1 o TGraph)
        color: Color de línea y marcador (string evaluable, e.g. 'ROOT.kRed')
        linestyle: Estilo de línea (int)
        markerstyle: Estilo de marcador (int)
        markersize: Tamaño del marcador (float)
        linewidth: Grosor de línea (int)
        fill: Si True, aplica relleno al histograma (bool)
        fillstyle: Estilo del relleno (int, default 3004 = rayado diagonal)
        fillalpha: Transparencia del relleno 0.0 (transparente) a 1.0 (opaco)
    """
    line_color = eval(color)
    obj.SetLineColor(line_color)
    obj.SetMarkerColor(line_color)
    obj.SetLineStyle(linestyle)
    obj.SetMarkerStyle(markerstyle)
    obj.SetMarkerSize(markersize)
    obj.SetLineWidth(linewidth)
    
    # Aplicar fill si está habilitado y el objeto es TH1
    if fill and isinstance(obj, ROOT.TH1):
        if fillalpha < 1.0:
            # Usar transparencia (requiere ROOT >= 6)
            # SetFillColorAlpha toma color y alpha (0-1)
            obj.SetFillColorAlpha(line_color, fillalpha)
        else:
            obj.SetFillColor(line_color)
        obj.SetFillStyle(fillstyle)

def plot_compare_1D_across_files(files_info, plots, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    for fig_name, cfg in plots.items():
        # Tamaño del título del plot (TPaveText). Se aplica vía gStyle antes de dibujar
        # y se restaura al final para no afectar a otros plots.
        _orig_title_fs = ROOT.gStyle.GetTitleFontSize()
        if "title_size" in cfg:
            ROOT.gStyle.SetTitleFontSize(float(cfg["title_size"]))
        c = ROOT.TCanvas(f"c_compare_{fig_name}", fig_name, 900, 700)
        if cfg.get("gridx", False):
            c.SetGridx()
        if cfg.get("gridy", False):
            c.SetGridy()
        c.SetLeftMargin(0.14)
        c.SetRightMargin(0.05)
        c.SetBottomMargin(0.12)
        c.SetTopMargin(0.08)
        leg_pos = cfg.get("legend", [0.65, 0.70, 0.88, 0.88])
        legend = ROOT.TLegend(*leg_pos)
        legend_text_size = cfg.get("legend_txt_size", 0.03)
        legend_margin = cfg.get("legend_margin", 0.2)
        legend.SetMargin(legend_margin)
        legend.SetNColumns(cfg.get("legend_ncolumns", 1))
        legend.SetTextSize(legend_text_size)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        # ¿Gráfico de eficiencia? -> eje 0..1.1 y sin notación científica
        # Heurística por nombre del plot
        is_eff = ("effi" in fig_name.lower()) or ("eff" in fig_name.lower())

        first = True
        global_max = 0.0
        global_min = None
        draw_as_scatter = False  # si alguna serie es TGraphAsymmErrors
        objs = []                # [(obj, label, kind, force_draw)], kind in {"main", "common"}

        # --- Normalización por defecto para las series principales
        norm_main = cfg.get("normalize", "none").lower()
        if norm_main not in ("none", "max", "integral"):
            norm_main = "none"
        rebin = cfg.get("rebin", None)
        # Modo de dibujo global para las series principales (ej. "HIST E", "E1", "P").
        # Vacío = comportamiento automático de _draw_one_object_1d.
        main_draw_mode = str(cfg.get("draw_mode", "")).upper()

        # --- 1) Recuperar, estilizar y (si procede) normalizar las series principales
        list_graph_values = dict()
        entries_per_dataset = {}

        ignore_entries_for  = cfg.get("ignore_entries_for", [])
        # --- Procesar entradas con "sum_histos" (histogramas sumados de uno o varios archivos) ---
        sum_entries = cfg.get("sum_histos", {})
        # Crear diccionario de datasets por label para acceso rápido
        ds_by_label = {ds["label"]: ds for ds in files_info}
        
        for sum_label, sum_cfg in sum_entries.items():
            # sum_cfg puede tener:
            #   histos: lista de nombres de histogramas (todos del mismo archivo "from")
            #   O
            #   sources: lista de {from: "dataset_label", histo: "histo_name", weight: opcional}
            
            sum_histo = None
            
            if "sources" in sum_cfg:
                # Modo multi-archivo: cada source especifica de qué dataset y qué histograma
                for src_entry in sum_cfg["sources"]:
                    src_label = src_entry.get("from")
                    hname = src_entry.get("histo")
                    
                    if src_label not in ds_by_label:
                        print(f"[WARN] Plot '{fig_name}': sum_histos source from='{src_label}' no coincide con ningún dataset.")
                        continue
                    
                    src_ds = ds_by_label[src_label]
                    f_src = src_ds["file"]
                    h = f_src.Get(hname)
                    
                    src_weight = src_entry.get("weight", src_ds.get("weight", 1.0))
                    if not h:
                        print(f"[WARN] '{hname}' not found in '{src_label}' for sum '{sum_label}'.")
                        continue
                    
                    if sum_histo is None:
                        sum_histo = h.Clone(f"sum_{sum_label}_clone")
                        if src_weight != 1.0:
                            sum_histo.Scale(src_weight)
                    else:
                        h_clone = h.Clone(f"{hname}_{src_label}_temp")
                        if src_weight != 1.0:
                            h_clone.Scale(src_weight)
                        sum_histo.Add(h_clone)
            
            elif "histos" in sum_cfg:
                # Modo archivo único: todos los histogramas del mismo dataset
                histo_names = sum_cfg["histos"]
                src_label = sum_cfg.get("from", files_info[0]["label"] if files_info else None)
                
                if src_label not in ds_by_label:
                    print(f"[WARN] Plot '{fig_name}': sum_histos.from='{src_label}' no coincide con ningún dataset.")
                    continue
                
                src_ds = ds_by_label[src_label]
                f_src = src_ds["file"]
                
                for i, hname in enumerate(histo_names):
                    h = f_src.Get(hname)
                    if not h:
                        print(f"[WARN] '{hname}' not found in '{src_label}' for sum '{sum_label}'.")
                        continue
                    
                    if sum_histo is None:
                        sum_histo = h.Clone(f"sum_{sum_label}_clone")
                    else:
                        sum_histo.Add(h)
            
            if sum_histo is None:
                print(f"[WARN] No se pudo crear histograma sumado '{sum_label}'.")
                continue
            
            # Aplicar peso global si está definido
            global_weight = sum_cfg.get("weight", 1.0)
            if global_weight != 1.0 and isinstance(sum_histo, ROOT.TH1):
                sum_histo.Scale(global_weight)
            
            # Aplicar estilo (usar valores por defecto si no se especifican)
            s_color = sum_cfg.get("color", "ROOT.kBlack")
            s_linestyle = sum_cfg.get("linestyle", 1)
            s_markerstyle = sum_cfg.get("markerstyle", 20)
            s_markersize = sum_cfg.get("markersize", 1.5)
            s_linewidth = sum_cfg.get("linewidth", 2)
            s_fill = sum_cfg.get("fill", False)
            s_fillstyle = sum_cfg.get("fillstyle", 3004)
            s_fillalpha = sum_cfg.get("fillalpha", 1.0)
            _apply_style_1d(
                sum_histo, s_color, s_linestyle, s_markerstyle, s_markersize, s_linewidth,
                fill=s_fill, fillstyle=s_fillstyle, fillalpha=s_fillalpha
            )
            
            if sum_label not in ignore_entries_for:
                entries_per_dataset[sum_label] = get_total_entries(
                    sum_histo, raw=str(cfg.get("show_entries", False)).lower() == "raw")
            
            # Normalización
            if isinstance(sum_histo, ROOT.TH1):
                if rebin and isinstance(rebin, int) and rebin > 1:
                    sum_histo.Rebin(rebin)
                if norm_main == "max":
                    m = sum_histo.GetMaximum()
                    if m > 0:
                        sum_histo.Scale(1.0 / m)
                elif norm_main == "integral":
                    integ = sum_histo.Integral()
                    if integ > 0:
                        sum_histo.Scale(1.0 / integ)
                global_max = max(global_max, sum_histo.GetMaximum())
                ymin = sum_histo.GetMinimum()
                global_min = ymin if global_min is None else min(global_min, ymin)
            
            if sum_cfg.get("label"):
              sum_label_name = sum_cfg["label"]
            else:
              sum_label_name = sum_label
            
            # Determinar modo de dibujo
            # draw_with_errors: dibuja solo puntos con barras de error estadísticas (sin histograma en bins)
            sum_draw_mode = sum_cfg.get("draw", "")
            draw_errors_overlay = False
            if not sum_draw_mode:
                if sum_cfg.get("draw_with_errors", False):
                    # Solo puntos con errores en X e Y, sin histograma en bins
                    sum_draw_mode = "E1 P"
                    draw_errors_overlay = False
                elif sum_cfg.get("markerstyle"):
                    sum_draw_mode = "HIST P"  # Mostrar línea + puntos
            
            objs.append((sum_histo, sum_label_name, "sum", sum_draw_mode))
            
            # Ya no necesitamos el overlay porque dibujamos directamente con errores
            if draw_errors_overlay:
                sum_histo_errors = sum_histo.Clone(f"sum_{sum_label}_errors_clone")
                # El clon hereda el estilo, solo cambiamos el modo de dibujo
                objs.append((sum_histo_errors, None, "sum_errors_overlay", "E1 P"))
        for ds in files_info:
            label = ds["label"]
            f = ds["file"]

            # nombre del objeto por dataset o común
            if "per_dataset" in cfg:
                hname = cfg["per_dataset"].get(label, None)
                if hname is None:
                    print(f"[WARN] Plot '{fig_name}': dataset '{label}' has not 'per_dataset'.")
                    continue
            else:
                hname = cfg["name"]

            obj = f.Get(hname)
            if not obj:
                print(f"[WARN] '{hname}' not found in data '{label}'.")
                continue

            # Detectamos si es TGraphAsymmErrors
            if isinstance(obj, ROOT.TGraphAsymmErrors):
                draw_as_scatter = True

            # Clonamos para no alterar el original
            o = obj.Clone(f"{hname}_{label}_clone")
            weight = ds.get("weight", 1.0)
            if weight != 1.0:
              if isinstance(o, ROOT.TH1):
                  o.Scale(weight)
            if isinstance(o, ROOT.TGraphAsymmErrors):
              n_points = o.GetN()
              x_vals = [float(o.GetPointX(i)) for i in range(n_points)]
              y_vals = [float(o.GetPointY(i)) for i in range(n_points)]
              y_err_lo = [float(o.GetErrorYlow(i)) for i in range(n_points)]
              y_err_hi = [float(o.GetErrorYhigh(i)) for i in range(n_points)]

    # Guarda también errores (low/high)
              list_graph_values[label] = (x_vals, y_vals, y_err_lo, y_err_hi)
              # list_graph_values[label] = (x_vals, y_vals)
            # Flag global de sombreado: con cfg['shade_edges'] activo, los histogramas que no
            # piden su propio 'fill' se rellenan con un patrón rayado (borde sombreado).
            # OJO: CompareAlgs inyecta fill=False por defecto en cada dataset, así que aquí se
            # combina con OR (no se puede usar ds.get('fill', _shade)). shade_style/shade_alpha
            # ajustan patrón y transparencia (alpha<1 puede no renderizar en PNG; hatch a 1.0 sí).
            _shade = cfg.get("shade_edges", False)
            _ds_fill = ds.get("fill", False)
            if _shade and not _ds_fill:
                _fill_on, _fill_style, _fill_alpha = True, cfg.get("shade_style", 3354), cfg.get("shade_alpha", 1.0)
            else:
                _fill_on, _fill_style, _fill_alpha = _ds_fill, ds.get("fillstyle", 3004), ds.get("fillalpha", 1.0)
            _apply_style_1d(
                o,
                ds["color"],
                ds["linestyle"],
                ds["markerstyle"],
                ds["markersize"],
                ds.get("linewidth", 2),
                fill=_fill_on,
                fillstyle=_fill_style,
                fillalpha=_fill_alpha,
            )
            if label not in ignore_entries_for:
              entries_per_dataset[label] = get_total_entries(
                  o, raw=str(cfg.get("show_entries", False)).lower() == "raw")

            # Normalización (series principales)
            if isinstance(o, ROOT.TH1):
                if rebin and isinstance(rebin, int) and rebin > 1:
                    o.Rebin(rebin)
                if norm_main == "max":
                    m = o.GetMaximum()
                    if m > 0:
                        o.Scale(1.0 / m)
                elif norm_main == "integral":
                    integ = o.Integral()
                    if integ > 0:
                        o.Scale(1.0 / integ)
                global_max = max(global_max, o.GetMaximum())
                ymin = o.GetMinimum()
                global_min = ymin if global_min is None else min(global_min, ymin)
                # Populate list_graph_values for diff panel (after normalization)
                nb = o.GetNbinsX()
                x_vals  = [o.GetBinCenter(b)  for b in range(1, nb + 1)]
                y_vals  = [o.GetBinContent(b) for b in range(1, nb + 1)]
                y_errs  = [o.GetBinError(b)   for b in range(1, nb + 1)]
                list_graph_values[label] = (x_vals, y_vals, y_errs, y_errs)
            else:
                # TGraph: estimamos extremos Y a partir de puntos
                n = o.GetN()
                if n > 0:
                    ys = [o.GetPointY(i) for i in range(n)]
                    if ys:
                        global_max = max(global_max, max(ys))
                        mn = min(ys)
                        global_min = mn if global_min is None else min(global_min, mn)

            objs.append((o, label, "main", main_draw_mode))
        total_entries = sum(entries_per_dataset.values()) if entries_per_dataset else 0

        percentages = {}
        if total_entries > 0:
            for k, v in entries_per_dataset.items():
                percentages[k] = 100.0 * v / total_entries
        # --- 2) Recuperar la SERIE COMÚN (opcional) y añadirla UNA sola vez
        common_cfg_list = cfg.get("common_list", None)
        if common_cfg_list:
          for common_cfg_key in common_cfg_list:  
              common_cfg = common_cfg_list.get(common_cfg_key, None)
              # dataset fuente (o primero)
              src_label = common_cfg.get("from", files_info[0]["label"] if files_info else None)
              src = next((ds for ds in files_info if ds["label"] == src_label), None)
              if src is None:
                  print(f"[WARN] Plot '{fig_name}': common.from='{src_label}' no coincide con ningún dataset.")
              else:
                  # nombre del histo común (por dataset o único)
                  if "per_dataset" in common_cfg:
                      hname_c = common_cfg["per_dataset"].get(src_label, None)
                  else:
                      hname_c = common_cfg.get("name", None)

                  if not hname_c:
                      print(f"[WARN] Plot '{fig_name}': 'common' sin 'name' ni 'per_dataset'.")
                  else:
                      obj_c = src["file"].Get(hname_c)
                      if not obj_c:
                          print(f"[WARN] Plot '{fig_name}': no se encontró común '{hname_c}' en '{src_label}'.")
                      else:
                          oc = obj_c.Clone(f"{hname_c}_{src_label}_common_clone")

                          # Estilo propio del común (por defecto: negro, línea discontinua, marcador hueco)
                          c_color       = common_cfg.get("color", ROOT.kBlack)
                          c_linestyle   = common_cfg.get("linestyle", 2)
                          c_markerstyle = common_cfg.get("markerstyle", 24)
                          c_markersize  = common_cfg.get("markersize", 1.5)
                          c_linewidth   = common_cfg.get("linewidth", 2)
                          c_fill        = common_cfg.get("fill", False)
                          c_fillstyle   = common_cfg.get("fillstyle", 3004)
                          c_fillalpha   = common_cfg.get("fillalpha", 1.0)

                          _apply_style_1d(
                              oc, c_color, c_linestyle, c_markerstyle, c_markersize, c_linewidth,
                              fill=c_fill, fillstyle=c_fillstyle, fillalpha=c_fillalpha
                          )

                          # Normalización específica del común (si falta, hereda la principal)
                          norm_c = common_cfg.get("normalize", norm_main).lower()
                          if norm_c not in ("none", "max", "integral"):
                              norm_c = "none"

                          if isinstance(oc, ROOT.TH1):
                              # aplicar mismo rebin que a las series principales
                              if rebin and isinstance(rebin, int) and rebin > 1:
                                  oc.Rebin(rebin)
                              if norm_c == "max":
                                  m = oc.GetMaximum()
                                  if m > 0:
                                      oc.Scale(1.0 / m)
                              elif norm_c == "integral":
                                  integ = oc.Integral()
                                  if integ > 0:
                                      oc.Scale(1.0 / integ)

                              global_max = max(global_max, oc.GetMaximum())
                              ymin = oc.GetMinimum()
                              global_min = ymin if global_min is None else min(global_min, ymin)
                          else:
                              # Si fuera TGraph
                              n = oc.GetN()
                              if n > 0:
                                  ys = [oc.GetPointY(i) for i in range(n)]
                                  if ys:
                                      global_max = max(global_max, max(ys))
                                      mn = min(ys)
                                      global_min = mn if global_min is None else min(global_min, mn)

                          # Etiqueta y posible modo de dibujo forzado
                          common_label = common_cfg.get("label", "Common")
                          force_draw   = common_cfg.get("draw", "").upper()  # "", "HIST", "P"

                          # Añade al final para que el común no tape las curvas principales
                          objs.append((oc, common_label, "common", force_draw))
        # dict label -> TH1 (series principales, post-normalización) para chi2_test
        main_histos = {label: o for (o, label, kind, _) in objs if kind == "main" and isinstance(o, ROOT.TH1)}

        diff_graph = None
        diff_graphs = []          # lista de TGraphAsymmErrors (uno por par)
        diff_ymin = None
        diff_ymax = None
        difference_line = cfg.get("diff_line", None)
        # ¿El panel inferior es cociente (valor a comparar / referencia) en vez de resta?
        # Se activa con diff_line: {mode: ratio} (alias: div/divide/cociente).
        _dcfg_mode = difference_line if isinstance(difference_line, dict) else {}
        diff_is_ratio = str(_dcfg_mode.get("mode", "diff")).lower() in ("ratio", "div", "divide", "cociente")
        # --- Diferencia punto a punto (se dibujará en un PAD inferior) ---
        if difference_line is not None and len(list_graph_values) >= 2:
            def _get_yerrs(values_tuple, n):
              # values_tuple puede ser (x, y) o (x, y, yerr_lo, yerr_hi)
              if len(values_tuple) >= 4:
                  lo = np.asarray(values_tuple[2], dtype=float)[:n]
                  hi = np.asarray(values_tuple[3], dtype=float)[:n]
              else:
                  lo = np.zeros(n, dtype=float)
                  hi = np.zeros(n, dtype=float)
              return lo, hi

            # --- Determinar la lista de pares a representar ---
            dcfg = difference_line if isinstance(difference_line, dict) else {}
            _all_labels = list(list_graph_values.keys())
            pairs = dcfg.get("pairs", None)
            if pairs:
                # Lista explícita de pares [ref, other]
                diff_pairs = [tuple(p) for p in pairs]
            else:
                # Compatibilidad hacia atrás: un único par con los dos primeros labels
                diff_pairs = [(_all_labels[0], _all_labels[1])]

            # Para el rango Y automático (sobre todos los pares)
            global_dmin = None
            global_dmax = None

            for pair_idx, pair in enumerate(diff_pairs):
                lbl_a, lbl_b = pair[0], pair[1]
                if lbl_a not in list_graph_values or lbl_b not in list_graph_values:
                    print(f"[WARN] Plot '{fig_name}': diff_line pair "
                          f"({lbl_a}, {lbl_b}) labels not found in datasets; skipping.")
                    continue

                x0 = np.asarray(list_graph_values[lbl_a][0], dtype=float)
                y0 = np.asarray(list_graph_values[lbl_a][1], dtype=float)
                x1 = np.asarray(list_graph_values[lbl_b][0], dtype=float)
                y1 = np.asarray(list_graph_values[lbl_b][1], dtype=float)

                min_len = min(len(x0), len(x1), len(y0), len(y1))
                x0, y0 = x0[:min_len], y0[:min_len]
                x1, y1 = x1[:min_len], y1[:min_len]

                # Si x difiere entre métodos, aquí podrías interpolar (np.interp); por ahora, por índice.
                # mode='ratio' → cociente y1/y0 = (valor a comparar lbl_b) / (referencia lbl_a);
                # por defecto resta y0 - y1.
                if diff_is_ratio:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        diff_vals = np.where(y0 != 0.0, y1 / y0, np.nan)
                else:
                    diff_vals = y0 - y1

                e0_lo, e0_hi = _get_yerrs(list_graph_values[lbl_a], min_len)
                e1_lo, e1_hi = _get_yerrs(list_graph_values[lbl_b], min_len)

                if diff_is_ratio:
                    # Error del cociente R=y1/y0: sigma_R = |R|·sqrt((s1/y1)^2 + (s0/y0)^2)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        rel0_lo = np.where(y0 != 0.0, e0_lo / np.abs(y0), 0.0)
                        rel0_hi = np.where(y0 != 0.0, e0_hi / np.abs(y0), 0.0)
                        rel1_lo = np.where(y1 != 0.0, e1_lo / np.abs(y1), 0.0)
                        rel1_hi = np.where(y1 != 0.0, e1_hi / np.abs(y1), 0.0)
                        diff_err_lo = np.nan_to_num(np.abs(diff_vals) * np.sqrt(rel1_lo**2 + rel0_lo**2))
                        diff_err_hi = np.nan_to_num(np.abs(diff_vals) * np.sqrt(rel1_hi**2 + rel0_hi**2))
                else:
                    diff_err_lo = np.sqrt(e0_lo**2 + e1_lo**2)
                    diff_err_hi = np.sqrt(e0_hi**2 + e1_hi**2)
                g = ROOT.TGraphAsymmErrors(min_len)
                for i in range(min_len):
                    xi  = float(x0[i])
                    yi  = float(diff_vals[i])
                    if not np.isfinite(yi):            # cociente con referencia 0 → punto en la línea base
                        yi = 1.0 if diff_is_ratio else 0.0
                    elo = float(diff_err_lo[i])
                    ehi = float(diff_err_hi[i])
                    g.SetPoint(i, xi, yi)
                    # errores en X = 0; en Y asimétricos (low/high)
                    g.SetPointError(i, 0.0, 0.0, elo, ehi)

                # Estilo: color/linestyle/marker del dataset que se RESTA (lbl_b),
                # no el de referencia (la diferencia es y(lbl_a) - y(lbl_b)).
                g.SetName(f"gdiff_{fig_name}_{pair_idx}")
                g.SetTitle("")
                ref_ds = ds_by_label.get(lbl_b, None)
                if ref_ds is not None:
                    _apply_style_1d(
                        g,
                        ref_ds["color"],
                        ref_ds["linestyle"],
                        ref_ds["markerstyle"],
                        ref_ds["markersize"],
                        ref_ds.get("linewidth", 2),
                    )
                else:
                    # Fallback al estilo gris original
                    g.SetMarkerStyle(20)
                    g.SetMarkerSize(1.0)
                    g.SetMarkerColor(ROOT.kGray+2)
                    g.SetLineColor(ROOT.kGray+2)
                    g.SetLineStyle(1)

                diff_graphs.append(g)

                # Acumular extremos para el rango Y automático (nan-safe para el cociente)
                if diff_vals.size:
                    _lo_arr = diff_vals - diff_err_lo
                    _hi_arr = diff_vals + diff_err_hi
                    if np.any(np.isfinite(_lo_arr)) and np.any(np.isfinite(_hi_arr)):
                        dmin_p = float(np.nanmin(_lo_arr))
                        dmax_p = float(np.nanmax(_hi_arr))
                        global_dmin = dmin_p if global_dmin is None else min(global_dmin, dmin_p)
                        global_dmax = dmax_p if global_dmax is None else max(global_dmax, dmax_p)

            # Conservar 'diff_graph' como primer elemento para no romper
            # las referencias del pad superior (rango X) y compatibilidad.
            diff_graph = diff_graphs[0] if diff_graphs else None

            # --- Rango Y del panel inferior ---
            if "y_range" in dcfg:
                diff_ymin, diff_ymax = float(dcfg["y_range"][0]), float(dcfg["y_range"][1])
            else:
                if global_dmin is not None and global_dmax is not None:
                    dmin, dmax = global_dmin, global_dmax
                    pad = 0.1 * max(1e-9, (dmax - dmin) if (dmax != dmin) else abs(dmax) + 1.0)
                    diff_ymin, diff_ymax = dmin - pad, dmax + pad
                else:
                    diff_ymin, diff_ymax = (0.9, 1.1) if diff_is_ratio else (-0.1, 0.1)
            # gx = array('d', x0.tolist())
            # gy = array('d', diff_vals.tolist())
            # diff_graph = ROOT.TGraph(len(gx), gx, gy)
            # diff_graph.SetName(f"gdiff_{fig_name}")
            # diff_graph.SetTitle("")
            # diff_graph.SetMarkerStyle(20)
            # diff_graph.SetMarkerSize(1.0)
            # diff_graph.SetMarkerColor(ROOT.kGray+2)
            # diff_graph.SetLineColor(ROOT.kGray+2)

            # # Rango Y del panel inferior
            # dcfg = difference_line if isinstance(difference_line, dict) else {}
            # if "y_range" in dcfg:
            #     diff_ymin, diff_ymax = float(dcfg["y_range"][0]), float(dcfg["y_range"][1])
            # else:
            #     # auto (incluye margen)
            #     if diff_vals.size:
            #         dmin, dmax = float(np.min(diff_vals)), float(np.max(diff_vals))
            #         pad = 0.1 * max(1e-9, (dmax - dmin) if (dmax != dmin) else abs(dmax) + 1.0)
            #         diff_ymin, diff_ymax = dmin - pad, dmax + pad
            #     else:
            #         diff_ymin, diff_ymax = -0.1, 0.1



        # --- 3) Dibujar
        # --- 3) Dibujar en dos pads si hay diferencia ---
        # Layout: top (70%) para principales, bottom (30%) para diferencia
        have_diff_panel = diff_graph is not None
        diff_panel_xrange = None  # rango X compartido entre pad superior e inferior

        if have_diff_panel:
            pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.30, 1.0, 1.0)
            pad_bot = ROOT.TPad("pad_bot", "pad_bot", 0.0, 0.00, 1.0, 0.30)
            pad_top.SetTickx(1)
            pad_top.SetBottomMargin(0.02)   # que no moleste al pad de abajo
            pad_bot.SetTopMargin(0.05)
            pad_bot.SetBottomMargin(0.35)   # deja sitio a las etiquetas del eje x
            if cfg.get("gridx", False): pad_top.SetGridx()
            if cfg.get("gridy", False): pad_top.SetGridy()
            if cfg.get("logy", False):
                pad_top.SetLogy()

            pad_top.Draw()
            pad_bot.Draw()
            pad_top.cd()
        else:
            # dibujo normal en el canvas completo
            if cfg.get("logy", False):
                c.SetLogy()
                # --- Añadir entradas extra a la leyenda (solo texto, sin marca ni línea) ---
        extra_legend = cfg.get("extra_legend", [])
        if isinstance(extra_legend, str):
            extra_legend = [extra_legend]  # permite pasar un único string
        for extra_text in extra_legend:
            # Creamos un objeto "nulo" invisible para añadir solo texto
            null_obj = ROOT.TObject()
            legend.AddEntry(null_obj, extra_text, "")
        # Ejes/títulos se toman del primer objeto que dibujamos
        for i, (o, label, kind, force_draw) in enumerate(objs):
            if i == 0:
                # Títulos ejes
                if isinstance(o, ROOT.TGraphAsymmErrors):
                    o.GetXaxis().SetTitle(cfg.get("x", ""))
                    o.GetYaxis().SetTitle(cfg.get("y", ""))
                    # Rango Y
                    if "y_range" in cfg:
                        ymin, ymax = cfg["y_range"]
                        o.GetYaxis().SetRangeUser(ymin, ymax)
                    else:
                        if is_eff:
                          if "diff_line" in cfg.keys():
                            o.GetYaxis().SetRangeUser(0.0, 1.1)
                          else:
                            o.GetYaxis().SetRangeUser(0.0, 1.1)
                              
                          o.GetYaxis().SetNoExponent(ROOT.kTRUE)
                        else:
                            if cfg.get("logy", False):
                                # evitar rangos inválidos
                                lo = max(1e-6, 0.001 if global_min is None else max(1e-6, 0.5*max(1e-6, global_min)))
                                hi = 1.2*max(1e-9, global_max)
                                o.GetYaxis().SetRangeUser(lo, hi)
                            else:
                                o.GetYaxis().SetRangeUser(0.0, 1.2*max(1e-9, global_max))
                else:
                    o.SetXTitle(cfg.get("x", ""))
                    o.SetYTitle(cfg.get("y", ""))
                    if "y_range" in cfg:
                        ymin, ymax = cfg["y_range"]
                        o.GetYaxis().SetRangeUser(ymin, ymax)
                    else:
                        if is_eff:
                            o.GetYaxis().SetRangeUser(0.0, 1.1)
                            o.GetYaxis().SetNoExponent(ROOT.kTRUE)
                        else:
                            if cfg.get("logy", False):
                                o.GetYaxis().SetRangeUser(0.001, 1.2*max(1e-9, global_max))
                            else:
                                o.GetYaxis().SetRangeUser(0.0, 1.2*max(1e-9, global_max))

                    if is_eff or cfg.get("no_exponent_y", False):
                        o.GetYaxis().SetNoExponent(ROOT.kTRUE)

                # Tamaños de título de ejes (pad principal). El título X solo se aplica aquí
                # si NO hay panel inferior; si lo hay, el eje X visible es el de abajo.
                if "ytitle_size" in cfg:
                    o.GetYaxis().SetTitleSize(float(cfg["ytitle_size"]))
                if "xtitle_size" in cfg and not have_diff_panel:
                    o.GetXaxis().SetTitleSize(float(cfg["xtitle_size"]))

                # Divisiones de ticks y longitud (pad principal).
                # yndiv/xndiv: número de divisiones (ej. 510 = 10 primarias, 5 secundarias).
                # ytick_length/xtick_length: longitud del tick como fracción del pad.
                # ylabel_size/xlabel_size: tamaño de las etiquetas numéricas.
                if "yndiv" in cfg:
                    o.GetYaxis().SetNdivisions(int(cfg["yndiv"]))
                if "ytick_length" in cfg:
                    o.GetYaxis().SetTickLength(float(cfg["ytick_length"]))
                if "ylabel_size" in cfg:
                    o.GetYaxis().SetLabelSize(float(cfg["ylabel_size"]))
                if not have_diff_panel:
                    if "xndiv" in cfg:
                        o.GetXaxis().SetNdivisions(int(cfg["xndiv"]))
                    if "xtick_length" in cfg:
                        o.GetXaxis().SetTickLength(float(cfg["xtick_length"]))
                    if "xlabel_size" in cfg:
                        o.GetXaxis().SetLabelSize(float(cfg["xlabel_size"]))

                # Título del canvas/objeto
                o.SetTitle(cfg.get("title", ""))
                # Ocultar etiquetas y título del eje X en el pad superior si hay panel de diferencia
                if have_diff_panel:
                  # Alinear el eje X de ambos pads usando el rango REAL de los datos.
                  # Para un TH1 hay que usar SetRangeUser (zoom): SetLimits remapea los
                  # bordes de bin y deforma/desplaza el histograma (lo sacaba fuera de +-1
                  # al copiar el rango auto-margen del TGraph de diferencia).
                  if isinstance(o, ROOT.TH1):
                      x_min = o.GetXaxis().GetXmin()
                      x_max = o.GetXaxis().GetXmax()
                      o.GetXaxis().SetRangeUser(x_min, x_max)
                  else:
                      x_min = diff_graph.GetXaxis().GetXmin()
                      x_max = diff_graph.GetXaxis().GetXmax()
                      o.GetXaxis().SetLimits(x_min, x_max)
                  diff_panel_xrange = (x_min, x_max)
                  xax = o.GetXaxis()
                  xax.SetNdivisions(int(cfg.get("xndiv", 510)))
                  xax.SetTickLength(float(cfg.get("xtick_length", 0.03)))

                  # Ocultar números/título en el pad superior (manteniendo ticks)
                  xax.SetLabelSize(0)
                  xax.SetTitleSize(0)


                # Dibujo primero
                if force_draw:
                    o.Draw(force_draw)
                else:
                    _draw_one_object_1d(o, True, draw_as_scatter)
            else:
                # Resto: dibujar con SAME (o forzar modo si se indicó)
                if isinstance(o, ROOT.TLine):
                  # Asegúrate de que ya haya un marco (algún objeto 2D) dibujado antes.
                  # Si es el primero (raro), tendrías que crear un frame, pero aquí asumimos que NO es el primero.
                  o.Draw("same")
                elif force_draw:
                    o.Draw(f"{force_draw} same")
                else:
                    _draw_one_object_1d(o, False, draw_as_scatter)

            # Leyenda (si es TGraph, usamos 'lp'; si forzó "P", 'p'; si histo, 'l')
            # Saltar objetos de overlay (sin etiqueta) para la leyenda
            if label is None:
                continue
                
            if isinstance(o, ROOT.TGraphAsymmErrors):
                legopt = "lp"
            else:
                # Detectar si el histograma tiene fill configurado
                has_fill = isinstance(o, ROOT.TH1) and o.GetFillStyle() != 0 and o.GetFillColor() != 0
                
                # Si force_draw contiene "P", mostrar punto en leyenda
                if force_draw and "P" in force_draw.upper():
                    legopt = "lp" if ("HIST" in force_draw.upper() or "L" in force_draw.upper()) else "p"
                    if has_fill:
                        legopt = "f" + legopt  # Añadir fill a la leyenda
                elif kind == "sum" and any(obj[2] == "sum_errors_overlay" and obj[0].GetName().startswith(f"sum_{label}") for obj in objs):
                    # Este histograma tiene un overlay de errores, mostrar como lp
                    legopt = "lep"
                    if has_fill:
                        legopt = "f" + legopt
                else:
                    legopt = "l" if not draw_as_scatter else "lp"
                    if has_fill:
                        legopt = "f" + legopt  # Añadir fill a la leyenda
            if isinstance(o, ROOT.TLine):
              legend.AddEntry(o, label, "l")
            else:
              # legend.AddEntry(o, label, legopt)
              show_entries = cfg.get("show_entries", False)
              # ignore_entries_for = cfg.get("ignore_entries_for", [])
              if show_entries and label in entries_per_dataset and label not in ignore_entries_for:
                  n = entries_per_dataset[label]
                  p = percentages.get(label, 0.0)
                  if show_entries == "absolute" or str(show_entries).lower() == "raw":
                    # 'raw' = conteo crudo (GetEntries) ya calculado en entries_per_dataset
                    label_ext = f"{label}  (N={int(n)})"
                  elif show_entries == "percent":
                    label_ext = f"{label}  ({p:.1f}%)"
                  else:
                    label_ext = f"{label}  (N={int(n)}, {p:.1f}%)"
              else:
                  label_ext = label

              legend.AddEntry(o, label_ext, legopt)

        # --- Añadir entradas extra a la leyenda (solo texto, sin marca ni línea) ---
        # extra_legend = cfg.get("extra_legend", [])
        # if isinstance(extra_legend, str):
        #     extra_legend = [extra_legend]  # permite pasar un único string
        # for extra_text in extra_legend:
        #     # Creamos un objeto "nulo" invisible para añadir solo texto
        #     null_obj = ROOT.TObject()
        #     legend.AddEntry(null_obj, extra_text, "")

        legend.Draw()

        # Etiqueta de canal / anotación libre (campo `label` del plot).
        # Soporta LaTeX de ROOT (TLatex). Posición configurable con `label_pos: [x, y]`
        # en coordenadas NDC; por defecto esquina superior izquierda (0.15, 0.85).
        _plot_label = cfg.get("label", None)
        if _plot_label is not None:
            _lx, _ly = cfg.get("label_pos", [0.15, 0.85])
            _lt = ROOT.TLatex(_lx, _ly, str(_plot_label))
            _lt.SetNDC(True)
            _lt.SetTextSize(float(cfg.get("label_size", 0.045)))
            _lt.SetTextFont(int(cfg.get("label_font", 42)))
            _lt.Draw()

        if have_diff_panel:
          pad_bot.cd()

          # Si quieres que el eje X solo se rotule aquí, reduce el tamaño de labels en el pad superior:
          # (ya hicimos SetBottomMargin en pad_top; esto suele bastar)
          dcfg2 = difference_line if isinstance(difference_line, dict) else {}
          title = dcfg2.get("label", "ratio" if diff_is_ratio else "#Delta")
          # Dibuja el primer scatter con ejes propios; el resto con "P SAME"
          diff_graph.Draw("AP")
          # Forzar el MISMO rango X que el pad superior para que los ejes queden alineados
          if diff_panel_xrange is not None:
              x_min, x_max = diff_panel_xrange
          else:
              x_min = diff_graph.GetXaxis().GetXmin()
              x_max = diff_graph.GetXaxis().GetXmax()
          diff_graph.GetXaxis().SetLimits(x_min, x_max)
          diff_graph.GetXaxis().SetNdivisions(int(cfg.get("diff_xndiv", cfg.get("xndiv", 510))))
          diff_graph.GetXaxis().SetTickLength(float(cfg.get("diff_xtick_length", cfg.get("xtick_length", 0.03))))
          diff_graph.GetYaxis().SetNdivisions(int(cfg.get("diff_yndiv", 505)))
          diff_graph.GetYaxis().SetTickLength(float(cfg.get("diff_ytick_length", 0.03)))
          diff_graph.GetXaxis().SetTitle(cfg.get("x", ""))
          diff_graph.GetYaxis().SetTitle(title)
          diff_graph.GetYaxis().SetRangeUser(diff_ymin, diff_ymax)

          # Tamaños de título y etiquetas (pad inferior). Configurables desde YAML.
          diff_graph.GetXaxis().SetTitleSize(float(cfg.get("xtitle_size", 0.10)))
          diff_graph.GetYaxis().SetTitleSize(float(cfg.get("diff_ytitle_size", 0.10)))
          diff_graph.GetXaxis().SetLabelSize(float(cfg.get("diff_xlabel_size", 0.08)))
          diff_graph.GetYaxis().SetLabelSize(float(cfg.get("diff_ylabel_size", 0.08)))
          diff_graph.GetYaxis().SetTitleOffset(0.4)

          # Dibujar las líneas de diferencia restantes (un par adicional cada una),
          # cada una con el color/estilo del dataset de referencia del par.
          for g in diff_graphs[1:]:
              g.Draw("P SAME")

          # Línea horizontal de referencia: y=1 para cociente, y=0 para resta
          x_min = diff_graph.GetXaxis().GetXmin()
          x_max = diff_graph.GetXaxis().GetXmax()
          _yref = 1.0 if diff_is_ratio else 0.0
          y0line = ROOT.TLine(x_min, _yref, x_max, _yref)
          y0line.SetLineStyle(2)
          y0line.SetLineColor(ROOT.kGray+1)
          y0line.Draw("same")

          pad_bot.Update()
        c.Update()

        # Guardar CSV con datos de entradas si show_entries está activo
        show_entries = cfg.get("show_entries", False)
        if show_entries and entries_per_dataset:
            save_entries_csv(entries_per_dataset, percentages, total_entries, outdir, fig_name)

        # Chi2Test entre la referencia (primer dataset) y el resto
        # Opción: chi2_test: true  (referencia = primer dataset de per_dataset / datasets)
        # o chi2_test: "label de referencia"  para elegir explícitamente.
        chi2_cfg = cfg.get("chi2_test", False)
        if chi2_cfg and len(main_histos) >= 2:
            labels_ordered = list(main_histos.keys())
            if isinstance(chi2_cfg, str) and chi2_cfg in main_histos:
                ref_label = chi2_cfg
            else:
                ref_label = labels_ordered[0]
            h_ref = main_histos[ref_label]
            chi2_results = []
            for lbl, h in main_histos.items():
                if lbl == ref_label:
                    continue
                # Chi2Test con WW (ambos histogramas pesados/normalizados)
                chi2_val = ctypes.c_double(0)
                ndf_val  = ctypes.c_int(0)
                igood    = ctypes.c_int(0)
                pval = h_ref.Chi2TestX(h, chi2_val, ndf_val, igood, "WW P")
                ndf = int(ndf_val.value)
                chi2_ndf = chi2_val.value / ndf if ndf > 0 else float("nan")
                chi2_results.append({
                    "Plot":      fig_name,
                    "Reference": ref_label,
                    "Dataset":   lbl,
                    "Chi2":      round(chi2_val.value, 4),
                    "NDF":       ndf,
                    "Chi2_NDF":  round(chi2_ndf, 4),
                    "pvalue":    round(float(pval), 6),
                })
                print(f"  chi2/ndf({ref_label} vs {lbl}) = {chi2_ndf:.3f}  (p={pval:.4f})")
            save_chi2_csv(chi2_results, outdir, fig_name)

        # Guardar
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, f"{fig_name}.png")
        c.SaveAs(outpath)
        print(f"[OK] Saved: {outpath}")
        c.Close()
        if "title_size" in cfg:
            ROOT.gStyle.SetTitleFontSize(_orig_title_fs)

