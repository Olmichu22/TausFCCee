import ROOT
ROOT.gROOT.SetBatch(True)
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import pandas as pd
import os
from array import array

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
              print(title)
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
        legend.SetTextSize(0.03)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        first = True
        normalize = cfg.get("norm", False)
        global_max_value = 0
        # print(global_max_value)
        # print("===============================")
        for histo_name in cfg["variabs"]:
          histo = file.Get(histo_name)
          if not histo:
            continue
          # print(histo)
          max_histo = get_graph_max(histo)
          # print(max_histo)
          # print("===============================")
          if max_histo > global_max_value:
            global_max_value = max_histo
        
        # For each histogram name in the group configuration
        linewidth = cfg.get("linewidth", 2)
        markersize = cfg.get("markersize", 1.2)
        for i, var in enumerate(cfg["variabs"]):
            histo = file.Get(var)
            if not histo:
                print(f"Warning: Histogram '{var}' not found in the ROOT file.")
                continue
            # Set the common X and Y axis titles from the group config

            
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
                histo.SetLineStyle(4)
              else:
                histo.SetMarkerStyle(20)   # círculos sólidos
                 
              histo.SetLineColor(ROOT.kMagenta)
              histo.SetMarkerColor(ROOT.kMagenta)
              histo.SetMarkerStyle(20)   # círculos sólidos
            elif i == 4:
              if not scatter_plot:
                histo.SetLineStyle(2)
              else:
                histo.SetMarkerStyle(20)   # círculos sólidos
              histo.SetLineColor(ROOT.kCyan)
              histo.SetMarkerColor(ROOT.kCyan)
              # histo.SetMarkerStyle(34)   # círculos sólidos
            elif i == 5:
              if not scatter_plot:
                histo.SetLineStyle(2)
              else:
                histo.SetMarkerStyle(20)   # círculos sólidos
              histo.SetLineColor(ROOT.kSpring)
              histo.SetMarkerColor(ROOT.kSpring)
            elif i == 6:
              if not scatter_plot:
                histo.SetLineStyle(2)
              else:
                histo.SetMarkerStyle(20)   # círculos sólidos
              histo.SetLineColor(ROOT.kYellow)
              histo.SetMarkerColor(ROOT.kYellow)
            elif i == 7:
              if not scatter_plot:
                histo.SetLineStyle(2)
              else:
                histo.SetMarkerStyle(20)   # círculos sólidos
              histo.SetLineColor(ROOT.kAzure+3)
              histo.SetMarkerColor(ROOT.kAzure+3)
                
            histo.SetLineWidth(linewidth)
            histo.SetMarkerSize(markersize)
            if cfg.get("fill", False):
              histo.SetFillStyle(3002)
            # Draw the first histogram normally; then draw others with "same"
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
                  histo.Draw("AP")
                else:
                  histo.Draw("P")
                # histo.Draw("L SAME")
              else:
                histo.Draw("HIST")
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
                histo.Draw("P same")
                # histo.Draw("L same")
              else:
                histo.Draw("HIST same")

                # Set axis from 0 to 1
                
            # Add legend entry using the corresponding label from config, if provided
            label = cfg["labels"][i] if i < len(cfg["labels"]) else var
            legend.AddEntry(histo, label, "l")
        
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
        c.Clear()
        histo.Draw("COLZ")
        out_file = os.path.join(out_dir, f"{var}.png")
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
    if isinstance(obj, ROOT.TGraphAsymmErrors):
        opt = "AP" if first else "P same"
    else:
        if draw_as_scatter:
            opt = "P" if first else "P same"
        else:
            opt = "HIST" if first else "HIST same"
    obj.Draw(opt)
    return opt

def _apply_style_1d(obj, color, linestyle, markerstyle, markersize, linewidth):
    obj.SetLineColor(color)
    obj.SetMarkerColor(color)
    obj.SetLineStyle(linestyle)
    obj.SetMarkerStyle(markerstyle)
    obj.SetMarkerSize(markersize)
    obj.SetLineWidth(linewidth)

def plot_compare_1D_across_files(files_info, plots, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for fig_name, cfg in plots.items():
        c = ROOT.TCanvas(f"c_compare_{fig_name}", fig_name, 900, 700)
        if cfg.get("gridx", False):
            c.SetGridx()
        if cfg.get("gridy", False):
            c.SetGridy()
        leg_pos = cfg.get("legend", [0.65, 0.70, 0.88, 0.88])
        legend = ROOT.TLegend(*leg_pos)
        legend.SetTextSize(0.03)
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

        # --- 1) Recuperar, estilizar y (si procede) normalizar las series principales
        list_graph_values = dict()
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
            if isinstance(o, ROOT.TGraphAsymmErrors):
              n_points = o.GetN()
              x_vals = [float(o.GetPointX(i)) for i in range(n_points)]
              y_vals = [float(o.GetPointY(i)) for i in range(n_points)]
              y_err_lo = [float(o.GetErrorYlow(i)) for i in range(n_points)]
              y_err_hi = [float(o.GetErrorYhigh(i)) for i in range(n_points)]

    # Guarda también errores (low/high)
              list_graph_values[label] = (x_vals, y_vals, y_err_lo, y_err_hi)
              # list_graph_values[label] = (x_vals, y_vals)
            _apply_style_1d(o, ds["color"], ds["linestyle"], ds["markerstyle"], ds["markersize"], ds.get("linewidth", 2))

            # Normalización (series principales)
            if isinstance(o, ROOT.TH1):
                if norm_main == "max":
                    m = o.GetMaximum()
                    if m > 0:
                        o.Scale(1.0 / m)
                elif norm_main == "integral":
                    integ = o.Integral()
                    if integ > 0:
                        o.Scale(1.0 / integ)
                if rebin and isinstance(rebin, int) and rebin > 1:
                    o.Rebin(rebin)
                global_max = max(global_max, o.GetMaximum())
                ymin = o.GetMinimum()
                global_min = ymin if global_min is None else min(global_min, ymin)
            else:
                # TGraph: estimamos extremos Y a partir de puntos
                n = o.GetN()
                if n > 0:
                    ys = [o.GetPointY(i) for i in range(n)]
                    if ys:
                        global_max = max(global_max, max(ys))
                        mn = min(ys)
                        global_min = mn if global_min is None else min(global_min, mn)

            objs.append((o, label, "main", ""))  # "" => sin modo de dibujo forzado

        # --- 2) Recuperar la SERIE COMÚN (opcional) y añadirla UNA sola vez
        common_cfg = cfg.get("common", None)
        if common_cfg:
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

                        _apply_style_1d(oc, c_color, c_linestyle, c_markerstyle, c_markersize, c_linewidth)

                        # Normalización específica del común (si falta, hereda la principal)
                        norm_c = common_cfg.get("normalize", norm_main).lower()
                        if norm_c not in ("none", "max", "integral"):
                            norm_c = "none"

                        if isinstance(oc, ROOT.TH1):
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
        diff_graph = None
        diff_ymin = None
        diff_ymax = None
        difference_line = cfg.get("diff_line", None)
        # --- Diferencia punto a punto (se dibujará en un PAD inferior) ---
        if difference_line is not None and len(list_graph_values) >= 2:
            _labels = list(list_graph_values.keys())
            x0 = np.asarray(list_graph_values[_labels[0]][0], dtype=float)
            y0 = np.asarray(list_graph_values[_labels[0]][1], dtype=float)
            x1 = np.asarray(list_graph_values[_labels[1]][0], dtype=float)
            y1 = np.asarray(list_graph_values[_labels[1]][1], dtype=float)
            
            min_len = min(len(x0), len(x1), len(y0), len(y1))
            x0, y0 = x0[:min_len], y0[:min_len]
            x1, y1 = x1[:min_len], y1[:min_len]
            
            # Si x difiere entre métodos, aquí podrías interpolar (np.interp); por ahora, resta por índice
            diff_vals = y0 - y1
            def _get_yerrs(values_tuple, n):
              # values_tuple puede ser (x, y) o (x, y, yerr_lo, yerr_hi)
              if len(values_tuple) >= 4:
                  lo = np.asarray(values_tuple[2], dtype=float)[:n]
                  hi = np.asarray(values_tuple[3], dtype=float)[:n]
              else:
                  lo = np.zeros(n, dtype=float)
                  hi = np.zeros(n, dtype=float)
              return lo, hi

            e0_lo, e0_hi = _get_yerrs(list_graph_values[_labels[0]], min_len)
            e1_lo, e1_hi = _get_yerrs(list_graph_values[_labels[1]], min_len)

            diff_err_lo = np.sqrt(e0_lo**2 + e1_lo**2)
            diff_err_hi = np.sqrt(e0_hi**2 + e1_hi**2)
            diff_graph = ROOT.TGraphAsymmErrors(min_len)
            for i in range(min_len):
                xi  = float(x0[i])
                yi  = float(diff_vals[i])
                elo = float(diff_err_lo[i])
                ehi = float(diff_err_hi[i])
                diff_graph.SetPoint(i, xi, yi)
                # errores en X = 0; en Y asimétricos (low/high)
                diff_graph.SetPointError(i, 0.0, 0.0, elo, ehi)

            # Estilo
            diff_graph.SetName(f"gdiff_{fig_name}")
            diff_graph.SetTitle("")
            diff_graph.SetMarkerStyle(20)
            diff_graph.SetMarkerSize(1.0)
            diff_graph.SetMarkerColor(ROOT.kGray+2)
            diff_graph.SetLineColor(ROOT.kGray+2)
            diff_graph.SetLineStyle(1)
            # Opcional: banda semitransparente alrededor de los puntos
            # diff_graph.SetFillColorAlpha(ROOT.kGray+1, 0.25)

            # --- Rango Y del panel inferior ---
            dcfg = difference_line if isinstance(difference_line, dict) else {}
            if "y_range" in dcfg:
                diff_ymin, diff_ymax = float(dcfg["y_range"][0]), float(dcfg["y_range"][1])
            else:
                if diff_vals.size:
                    dmin, dmax = float(np.min(diff_vals - diff_err_lo)), float(np.max(diff_vals + diff_err_hi))
                    pad = 0.1 * max(1e-9, (dmax - dmin) if (dmax != dmin) else abs(dmax) + 1.0)
                    diff_ymin, diff_ymax = dmin - pad, dmax + pad
                else:
                    diff_ymin, diff_ymax = -0.1, 0.1
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

        if have_diff_panel:
            pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.30, 1.0, 1.0)
            pad_bot = ROOT.TPad("pad_bot", "pad_bot", 0.0, 0.00, 1.0, 0.30)
            pad_top.SetTickx(1)
            pad_top.SetBottomMargin(0.02)   # que no moleste al pad de abajo
            pad_bot.SetTopMargin(0.05)
            pad_bot.SetBottomMargin(0.35)   # deja sitio a las etiquetas del eje x
            if cfg.get("gridx", False): pad_top.SetGridx()
            if cfg.get("gridy", False): pad_top.SetGridy()

            pad_top.Draw()
            pad_bot.Draw()
            pad_top.cd()
        else:
            # dibujo normal en el canvas completo
            if cfg.get("logy", False):
                c.SetLogy()

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

                # Título del canvas/objeto
                o.SetTitle(cfg.get("title", ""))
                # Ocultar etiquetas y título del eje X en el pad superior si hay panel de diferencia
                if have_diff_panel:
                                  # Usa el rango del gráfico de diferencia para forzar el mismo X en ambos pads
                  x_min = diff_graph.GetXaxis().GetXmin()
                  x_max = diff_graph.GetXaxis().GetXmax()

                  # Forzar mismo rango y mismas divisiones/ticklength en el pad superior
                  # (SetLimits para TGraph; para TH1 también funciona al haber frame creado)
                  o.GetXaxis().SetLimits(x_min, x_max)
                  xax = o.GetXaxis()
                  xax.SetNdivisions(510)   # mismas divisiones que abajo
                  xax.SetTickLength(0.03)  # misma longitud de tick

                  # Ocultar números/título en el pad superior (manteniendo ticks)
                  xax.SetLabelSize(0)
                  xax.SetTitleSize(0)


                # Dibujo primero
                if force_draw in ("HIST", "P"):
                    o.Draw(force_draw)
                else:
                    _draw_one_object_1d(o, True, draw_as_scatter)
            else:
                # Resto: dibujar con SAME (o forzar modo si se indicó)
                if isinstance(o, ROOT.TLine):
                  # Asegúrate de que ya haya un marco (algún objeto 2D) dibujado antes.
                  # Si es el primero (raro), tendrías que crear un frame, pero aquí asumimos que NO es el primero.
                  o.Draw("same")
                elif force_draw in ("HIST", "P"):
                    o.Draw(f"{force_draw} same")
                else:
                    _draw_one_object_1d(o, False, draw_as_scatter)

            # Leyenda (si es TGraph, usamos 'lp'; si forzó "P", 'p'; si histo, 'l')
            if isinstance(o, ROOT.TGraphAsymmErrors):
                legopt = "lp"
            else:
                legopt = "p" if force_draw == "P" else ("l" if not draw_as_scatter else "lp")
            if isinstance(o, ROOT.TLine):
              legend.AddEntry(o, label, "l")
            else:
              legend.AddEntry(o, label, legopt)

        legend.Draw()
        if have_diff_panel:
          pad_bot.cd()

          # Si quieres que el eje X solo se rotule aquí, reduce el tamaño de labels en el pad superior:
          # (ya hicimos SetBottomMargin en pad_top; esto suele bastar)
          title = difference_line["label"]
          # Dibuja el scatter con ejes propios
          diff_graph.Draw("AP")
          x_min = diff_graph.GetXaxis().GetXmin()
          x_max = diff_graph.GetXaxis().GetXmax()
          diff_graph.GetXaxis().SetLimits(x_min, x_max)
          diff_graph.GetXaxis().SetNdivisions(510)
          diff_graph.GetXaxis().SetTickLength(0.03)
          diff_graph.GetXaxis().SetTitle(cfg.get("x", ""))
          diff_graph.GetYaxis().SetTitle(title)
          diff_graph.GetYaxis().SetRangeUser(diff_ymin, diff_ymax)

          # Fuente un poco más pequeña para encajar
          diff_graph.GetXaxis().SetTitleSize(0.10)
          diff_graph.GetYaxis().SetTitleSize(0.10)
          diff_graph.GetXaxis().SetLabelSize(0.08)
          diff_graph.GetYaxis().SetLabelSize(0.08)
          diff_graph.GetYaxis().SetTitleOffset(0.4)

          # Línea horizontal y=0
          x_min = diff_graph.GetXaxis().GetXmin()
          x_max = diff_graph.GetXaxis().GetXmax()
          y0line = ROOT.TLine(x_min, 0.0, x_max, 0.0)
          y0line.SetLineStyle(2)
          y0line.SetLineColor(ROOT.kGray+1)
          y0line.Draw("same")

          pad_bot.Update()
        c.Update()

        # Guardar
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, f"{fig_name}.png")
        c.SaveAs(outpath)
        print(f"[OK] Saved: {outpath}")
        c.Close()

