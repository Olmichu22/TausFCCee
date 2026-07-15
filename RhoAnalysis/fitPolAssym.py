 #!/usr/bin/env python3
import ROOT
from ROOT import TFile, TCanvas, TMinuit, TLegend, gStyle, gPad
import ctypes
import math
import numpy as np
import argparse
import os

ROOT.gROOT.SetBatch(True)
gStyle.SetOptStat(0)

# =====================================================
# MAIN
# =====================================================
def main(
    input_file,
    nBins=20,
    rebin=2,
    force_perfect_agreement=False,
    use_likelihood=True,
    outdir=".",
    verbose=False,
    no_bg=False,
    bg_mode="total",  # "total", "mig", "ext", "split"
    extra_legend=None,  # Lista de textos extra para la leyenda
    lumi_base=None,     # Luminosidad base de los datos (fb^-1)
    lumi_target=None,   # Luminosidad objetivo para reescalar (fb^-1)
):
    fitTable = {
    "bin": [],
    "cos_min": [],
    "cos_max": [],
    "scale_signal": [],
    "scale_signal_err": [],
    "pol": [],
    "pol_err": [],
    "scale_bg": [],
    "scale_bg_err": [],
    "scale_bg_mig": [],
    "scale_bg_mig_err": [],
    "scale_bg_ext": [],
    "scale_bg_ext_err": [],
    }
    os.makedirs(outdir, exist_ok=True)
    images_outdir = os.path.join(outdir, "fit_images")
    os.makedirs(images_outdir, exist_ok=True)
    tag = "_perfect" if force_perfect_agreement else "_def"
    if no_bg:
        tag += "_nobg"
    elif bg_mode == "mig":
        tag += "_migonly"
    elif bg_mode == "ext":
        tag += "_extonly"
    elif bg_mode == "split":
        tag += "_split"
    out_root = os.path.join(outdir, f"fitPol_results{tag}.root")
    fout = ROOT.TFile(out_root, "RECREATE")

    
    
    file = TFile.Open(input_file)
    if not file or file.IsZombie():
        raise RuntimeError(f"No se pudo abrir {input_file}")

    # Calcular factor de escala de luminosidad
    lumi_scale = 1.0
    if lumi_base is not None and lumi_target is not None:
        if lumi_base <= 0:
            raise ValueError("lumi_base debe ser > 0")
        lumi_scale = lumi_target / lumi_base
        if verbose:
            print(f"[INFO] Reescalando luminosidad: {lumi_base:.2f} fb^-1 → {lumi_target:.2f} fb^-1 (factor: {lumi_scale:.3f})")

    binSize = 100.0 / nBins
    var = "histo"
    title = "#omega_{#rho}"

    vectorPol = {}
    vectorPolError = {}
    vectorCosThetaMin = {}
    vectorCosThetaMax = {}
    def _max_with_err(h):
        m = 0.0
        for i in range(1, h.GetNbinsX() + 1):
            m = max(m, h.GetBinContent(i) + h.GetBinError(i))
        return m
    # =================================================
    # FIT POR BIN
    # =================================================
    def dofit(bin, fullRange=False):

        if not fullRange:
            binName = f"_{bin}"
            vectorCosThetaMin[bin] = bin * binSize / 100 * 2 - 1
            vectorCosThetaMax[bin] = (bin + 1) * binSize / 100 * 2 - 1
        else:
            binName = "_full"
            vectorCosThetaMin[bin] = -1
            vectorCosThetaMax[bin] = 1

        hist_bg = file.Get(var + "_BG" + binName)
        hist_bg_mig = file.Get(var + "_BG_migrations" + binName)
        hist_bg_zqq = file.Get(var + "_BG_Zqq" + binName)
        hist_bg_bhabha = file.Get(var + "_BG_Bhabha" + binName)
        hist_m1 = file.Get(var + "_SIGNAL_M1" + binName)
        hist_p1 = file.Get(var + "_SIGNAL_P1" + binName)
        hist_sm = file.Get(var + "_SIGNAL" + binName)

        if not all([hist_m1, hist_p1, hist_sm]):
            raise RuntimeError(f"Histogramas de señal faltantes en bin {bin}")
        
        # Aplicar reescalado de luminosidad si está definido
        if lumi_scale != 1.0:
            hist_m1.Scale(lumi_scale)
            hist_p1.Scale(lumi_scale)
            hist_sm.Scale(lumi_scale)
            if hist_bg:
                hist_bg.Scale(lumi_scale)
            if hist_bg_mig:
                hist_bg_mig.Scale(lumi_scale)
            if hist_bg_zqq:
                hist_bg_zqq.Scale(lumi_scale)
            if hist_bg_bhabha:
                hist_bg_bhabha.Scale(lumi_scale)
        
        # Construir BG externo (Zqq + Bhabha) si existen
        hist_bg_ext = None
        if hist_bg_zqq and hist_bg_bhabha:
            hist_bg_ext = hist_bg_zqq.Clone(f"hist_bg_ext_{bin}")
            hist_bg_ext.Add(hist_bg_bhabha)
        elif hist_bg_zqq:
            hist_bg_ext = hist_bg_zqq.Clone(f"hist_bg_ext_{bin}")
        elif hist_bg_bhabha:
            hist_bg_ext = hist_bg_bhabha.Clone(f"hist_bg_ext_{bin}")

        # Verificar disponibilidad según modo
        if not no_bg:
            if bg_mode == "total" and not hist_bg:
                raise RuntimeError(f"histo_BG faltante en bin {bin}")
            if bg_mode == "mig" and not hist_bg_mig:
                raise RuntimeError(f"histo_BG_migrations faltante en bin {bin}")
            if bg_mode == "ext" and not hist_bg_ext:
                raise RuntimeError(f"Fondos externos (Zqq/Bhabha) faltantes en bin {bin}")
            if bg_mode == "split" and (not hist_bg_mig or not hist_bg_ext):
                raise RuntimeError(f"Fondos mig/ext faltantes en bin {bin}")

        # Datos
        if force_perfect_agreement:
            # hist_data = hist_bg.Clone("hist_data")
            ctheta = (vectorCosThetaMax[bin] + vectorCosThetaMin[bin]) / 2
            Atautheo = 0.14955426
            Aetheo = 0.14955426
            perfectPtau = -(Atautheo * (1 + ctheta**2) + 2 * Aetheo * ctheta) / (
                1 + ctheta**2 + 2 * Aetheo * Atautheo * ctheta
            )
            # hist_data.Add(hist_p1, (1 + perfectPtau) / 2)
            # hist_data.Add(hist_m1, (1 - perfectPtau) / 2)
            # if force_perfect_agreement:
            if no_bg:
                hist_data = hist_m1.Clone("hist_data")
                hist_data.Scale((1. - perfectPtau)/2)
                hist_data.Add(hist_p1, (1. + perfectPtau)/2)
            else:
                hist_data = hist_bg.Clone("hist_data")
                hist_data.Add(hist_p1, (1. + perfectPtau)/2)
                hist_data.Add(hist_m1, (1. - perfectPtau)/2)
        else:
            hist_data = hist_sm.Clone("hist_data")
            if not no_bg:
                if bg_mode == "total":
                    hist_data.Add(hist_bg)
                elif bg_mode == "mig":
                    hist_data.Add(hist_bg_mig)
                elif bg_mode == "ext":
                    hist_data.Add(hist_bg_ext)
                elif bg_mode == "split":
                    hist_data.Add(hist_bg_mig)
                    hist_data.Add(hist_bg_ext)

        hists_to_rebin = [hist_data, hist_m1, hist_p1]
        if hist_bg:
            hists_to_rebin.append(hist_bg)
        if hist_bg_mig:
            hists_to_rebin.append(hist_bg_mig)
        if hist_bg_ext:
            hists_to_rebin.append(hist_bg_ext)
        for h in hists_to_rebin:
            h.Rebin(rebin)

        # =================================================
        # FUNCIÓN DE COSTE
        # =================================================
        def my_minimization(npar, gin, f, par, iflag):
            val = 0.0
            for i in range(
                hist_data.GetXaxis().FindBin(-1),
                hist_data.GetXaxis().FindBin(1.4) + 1,
            ):
                observed = hist_data.GetBinContent(i)
                
                # Calcular contribución de BG según modo
                bg_contrib = 0.0
                if not no_bg:
                    if bg_mode == "total":
                        bg_contrib = par[2] * hist_bg.GetBinContent(i)
                    elif bg_mode == "mig":
                        bg_contrib = par[2] * hist_bg_mig.GetBinContent(i)
                    elif bg_mode == "ext":
                        bg_contrib = par[2] * hist_bg_ext.GetBinContent(i)
                    elif bg_mode == "split":
                        bg_contrib = par[2] * hist_bg_mig.GetBinContent(i) + par[3] * hist_bg_ext.GetBinContent(i)
                
                Nm1 = (1 - par[1]) / 2 * hist_m1.GetBinContent(i)
                Np1 = (1 + par[1]) / 2 * hist_p1.GetBinContent(i)

                expected = bg_contrib + par[0] * (Nm1 + Np1)
                if expected <= 0:
                    continue

                if use_likelihood:
                    val += 2 * (expected - observed * math.log(expected))
                else:
                    val += (observed - expected) ** 2 / expected

            f.value = val

        # =================================================
        # MINUIT
        # =================================================
        if no_bg:
            npar = 2
        elif bg_mode == "split":
            npar = 4  # scale_signal, pol, scale_bg_mig, scale_bg_ext
        else:
            npar = 3  # scale_signal, pol, scale_bg
        
        minuit = TMinuit(npar)
        minuit.SetFCN(my_minimization)
        minuit.SetPrintLevel(1)  # Set to 1 for standard output

        minuit.DefineParameter(0, "scale_signal", 1.0, 0.01, 0, 10)
        minuit.DefineParameter(1, "pol", 0.15, 0.0001, -1, 1)
        if not no_bg:
            if bg_mode == "split":
                minuit.DefineParameter(2, "scale_bg_mig", 1.0, 0.01, 0, 10)
                minuit.DefineParameter(3, "scale_bg_ext", 1.0, 0.01, 0, 10)
            else:
                minuit.DefineParameter(2, "scale_bg", 1.0, 0.01, 0, 10)

        migrad_result = minuit.Command("MIGRAD")
        if migrad_result != 0:
            print("Minimization did not converge.")
        # Resultados
        pval = ctypes.c_double(0)
        perr = ctypes.c_double(0)

        minuit.GetParameter(1, pval, perr)
        vectorPol[bin] = pval.value
        vectorPolError[bin] = perr.value

        if verbose:
            print(f"[BIN {bin}] P_tau = {pval.value:.5f} ± {perr.value:.5f}")
        # =============================
        # Construcción del histograma ajustado
        # =============================
        p_scale = ctypes.c_double(0)
        e_scale = ctypes.c_double(0)
        p_pol   = ctypes.c_double(0)
        e_pol   = ctypes.c_double(0)
        p_bg    = ctypes.c_double(0)
        e_bg = ctypes.c_double(0)

        minuit.GetParameter(0, p_scale, e_scale)
        minuit.GetParameter(1, p_pol,   e_pol)
        
        p_bg_mig = ctypes.c_double(0)
        e_bg_mig = ctypes.c_double(0)
        p_bg_ext = ctypes.c_double(0)
        e_bg_ext = ctypes.c_double(0)
        
        if not no_bg:
            if bg_mode == "split":
                minuit.GetParameter(2, p_bg, e_bg)  # reuse for mig
                p_bg_mig.value, e_bg_mig.value = p_bg.value, e_bg.value
                minuit.GetParameter(3, p_bg_ext, e_bg_ext)
                p_bg.value = 0.0  # no single scale_bg in split mode
                e_bg.value = 0.0
            else:
                minuit.GetParameter(2, p_bg, e_bg)
                if bg_mode == "mig":
                    p_bg_mig.value, e_bg_mig.value = p_bg.value, e_bg.value
                elif bg_mode == "ext":
                    p_bg_ext.value, e_bg_ext.value = p_bg.value, e_bg.value
        else:
            p_bg.value = 0.0
            e_bg.value = 0.0

        
        fitTable["bin"].append(bin)
        fitTable["cos_min"].append(vectorCosThetaMin[bin])
        fitTable["cos_max"].append(vectorCosThetaMax[bin])

        fitTable["scale_signal"].append(p_scale.value)
        fitTable["scale_signal_err"].append(e_scale.value)

        fitTable["pol"].append(p_pol.value)
        fitTable["pol_err"].append(e_pol.value)

        fitTable["scale_bg"].append(p_bg.value)
        fitTable["scale_bg_err"].append(e_bg.value)
        fitTable["scale_bg_mig"].append(p_bg_mig.value)
        fitTable["scale_bg_mig_err"].append(e_bg_mig.value)
        fitTable["scale_bg_ext"].append(p_bg_ext.value)
        fitTable["scale_bg_ext_err"].append(e_bg_ext.value)

        hist_m1_scaled = hist_m1.Clone()
        hist_p1_scaled = hist_p1.Clone()

        hist_m1_scaled.Scale(p_scale.value * (1 - p_pol.value) / 2)
        hist_p1_scaled.Scale(p_scale.value * (1 + p_pol.value) / 2)

        # Construir histograma de fit total según modo
        hist_bg_mig_plot = None
        hist_bg_ext_plot = None
        hist_bg_plot = None
        
        if no_bg:
            hist_fit_total = hist_m1_scaled.Clone("hist_fit_total")
            hist_fit_total.Add(hist_p1_scaled)
        else:
            if bg_mode == "total":
                hist_fit_total = hist_bg.Clone("hist_fit_total")
                hist_fit_total.Scale(p_bg.value)
                hist_bg_plot = hist_fit_total.Clone("hist_bg_plot")
            elif bg_mode == "mig":
                hist_fit_total = hist_bg_mig.Clone("hist_fit_total")
                hist_fit_total.Scale(p_bg.value)
                hist_bg_mig_plot = hist_fit_total.Clone("hist_bg_mig_plot")
            elif bg_mode == "ext":
                hist_fit_total = hist_bg_ext.Clone("hist_fit_total")
                hist_fit_total.Scale(p_bg.value)
                hist_bg_ext_plot = hist_fit_total.Clone("hist_bg_ext_plot")
            elif bg_mode == "split":
                # Crear histogramas escalados separados para dibujar
                hist_bg_mig_plot = hist_bg_mig.Clone("hist_bg_mig_plot")
                hist_bg_mig_plot.Scale(p_bg_mig.value)
                hist_bg_ext_plot = hist_bg_ext.Clone("hist_bg_ext_plot")
                hist_bg_ext_plot.Scale(p_bg_ext.value)
                # Fit total = mig + ext + señal
                hist_fit_total = hist_bg_mig_plot.Clone("hist_fit_total")
                hist_fit_total.Add(hist_bg_ext_plot)
            
            hist_fit_total.Add(hist_m1_scaled)
            hist_fit_total.Add(hist_p1_scaled)

        # =============================
        # Plot por bin
        # =============================
        suffix = "full" if fullRange else str(bin)
        canvas = TCanvas(f"canvas_{suffix}", "", 800, 600)
        canvas.SetLeftMargin(0.12)   # Aumentar margen izquierdo
        canvas.SetRightMargin(0.05)
        canvas.SetBottomMargin(0.12)
        # Quitar relleno explícitamente antes de dibujar
        hists_to_style = [hist_data, hist_fit_total, hist_m1_scaled, hist_p1_scaled]
        if hist_bg_plot:
            hists_to_style.append(hist_bg_plot)
        if hist_bg_mig_plot:
            hists_to_style.append(hist_bg_mig_plot)
        if hist_bg_ext_plot:
            hists_to_style.append(hist_bg_ext_plot)
        for h in hists_to_style:
            h.SetFillStyle(0)
            h.SetFillColor(0)
        
        # Calcular y_max considerando todos los histogramas de fondo
        y_max_candidates = [
            _max_with_err(hist_data),
            hist_fit_total.GetMaximum(),
            hist_m1_scaled.GetMaximum(),
            hist_p1_scaled.GetMaximum(),
        ]
        if hist_bg_plot:
            y_max_candidates.append(hist_bg_plot.GetMaximum())
        if hist_bg_mig_plot:
            y_max_candidates.append(hist_bg_mig_plot.GetMaximum())
        if hist_bg_ext_plot:
            y_max_candidates.append(hist_bg_ext_plot.GetMaximum())
        y_max = max(y_max_candidates)

        hist_data.SetMinimum(0.0)
        hist_data.SetMaximum(1.25 * y_max)
        
        hist_data.SetMarkerStyle(20)
        hist_data.SetLineColor(ROOT.kBlack)
        hist_data.Draw("E")

        hist_fit_total.SetLineColor(ROOT.kRed)
        hist_fit_total.SetLineWidth(2)
        hist_fit_total.Draw("HIST SAME")

        # Dibujar fondos según el modo
        if hist_bg_plot:
            # Modo "total"
            hist_bg_plot.SetLineColor(ROOT.kBlue)
            hist_bg_plot.SetFillStyle(1001)
            hist_bg_plot.SetFillColorAlpha(ROOT.kBlue, 0.5)
            hist_bg_plot.Draw("HIST SAME")
        if hist_bg_mig_plot:
            # Modo "mig" o "split"
            hist_bg_mig_plot.SetLineColor(ROOT.kMagenta)
            hist_bg_mig_plot.SetLineStyle(1)
            hist_bg_mig_plot.SetFillStyle(1001)
            hist_bg_mig_plot.SetFillColorAlpha(ROOT.kMagenta, 0.5)
            hist_bg_mig_plot.Draw("HIST SAME")
        if hist_bg_ext_plot:
            # Modo "ext" o "split"
            hist_bg_ext_plot.SetLineColor(ROOT.kBlue+2)
            hist_bg_ext_plot.SetFillStyle(1001)
            hist_bg_ext_plot.SetFillColorAlpha(ROOT.kBlue+2, 0.5)  
            # hist_bg_ext_plot.SetLineStyle(2)
            hist_bg_ext_plot.Draw("HIST SAME")

        hist_m1_scaled.SetLineColor(ROOT.kGreen+4)
        hist_m1_scaled.SetLineStyle(7)
        
        hist_p1_scaled.SetLineColor(ROOT.kGreen+4)
        hist_p1_scaled.SetLineStyle(10)
        hist_m1_scaled.Draw("HIST SAME")
        hist_p1_scaled.Draw("HIST SAME")

        leg = TLegend(0.38,0.6,0.88,0.88)
        #Ajustamos tamaño de letra
        leg.SetTextSize(0.045)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.AddEntry(hist_data,"Pseudo-data","lep")
        leg.AddEntry(hist_fit_total,"Total fit","l")
        # Añadir entradas de leyenda según el modo de fondo
        if hist_bg_plot:
            leg.AddEntry(hist_bg_plot, "Background", "l")
        if hist_bg_mig_plot:
            leg.AddEntry(hist_bg_mig_plot, "BG (migrations)", "l")
        if hist_bg_ext_plot:
            leg.AddEntry(hist_bg_ext_plot, "BG (external)", "l")
        leg.AddEntry(hist_m1_scaled,"A_{#tau}=+1","l")
        leg.AddEntry(hist_p1_scaled,"A_{#tau}=-1","l")
        # Añadir texto extra informativo a la leyenda
        if extra_legend:
            for extra_text in extra_legend:
                leg.AddEntry(ROOT.nullptr, extra_text, "")
        leg.Draw()
        # Títulos ejes
        hist_data.GetXaxis().SetTitle("\\omega_{\\rho}")
        hist_data.GetYaxis().SetTitle("Events")

        canvas.SaveAs(os.path.join(images_outdir, f"fit_result_minuit_{suffix}{tag}.png"))
    # =================================================
    # LOOP BINS
    # =================================================
    for b in range(nBins):
        dofit(b, False)
    dofit(nBins, True)

    # =================================================
    # TABLA DE CONTEOS DE EVENTOS (usando histos _full)
    # =================================================
    def _integral(hname):
        h = file.Get(hname)
        if not h:
            return 0.0
        return h.Integral() * lumi_scale

    n_sig  = _integral(f"{var}_SIGNAL_full")
    n_mig  = _integral(f"{var}_BG_migrations_full")
    n_zqq  = _integral(f"{var}_BG_Zqq_full")
    n_bhab = _integral(f"{var}_BG_Bhabha_full")
    n_ext  = n_zqq + n_bhab
    n_tot  = n_sig + n_mig + n_ext

    counts_txt = os.path.join(outdir, f"event_counts{tag}.txt")
    with open(counts_txt, "w") as fc:
        fc.write(f"{'Source':<22} {'Events':>12} {'%':>8}\n")
        fc.write("-" * 44 + "\n")
        fc.write(f"{'Signal (SS)':<22} {n_sig:>12.1f} {100*n_sig/n_tot if n_tot else 0:>7.2f}%\n")
        fc.write(f"{'BG migrations':<22} {n_mig:>12.1f} {100*n_mig/n_tot if n_tot else 0:>7.2f}%\n")
        if n_zqq > 0:
            fc.write(f"{'BG Zqq':<22} {n_zqq:>12.1f} {100*n_zqq/n_tot if n_tot else 0:>7.2f}%\n")
        if n_bhab > 0:
            fc.write(f"{'BG Bhabha':<22} {n_bhab:>12.1f} {100*n_bhab/n_tot if n_tot else 0:>7.2f}%\n")
        if n_ext > 0:
            fc.write(f"{'BG external (total)':<22} {n_ext:>12.1f} {100*n_ext/n_tot if n_tot else 0:>7.2f}%\n")
        fc.write("-" * 44 + "\n")
        fc.write(f"{'Total':<22} {n_tot:>12.1f} {'100.00%':>8}\n")
    if verbose:
        with open(counts_txt) as fc:
            print("[INFO] Tabla de conteos:\n" + fc.read())

    def my_minimization2(npar, gin, f, par, iflag):
        val = 0.0
        for b in range(nBins):
            ctheta = 0.5 * (vectorCosThetaMin[b] + vectorCosThetaMax[b])
            obs = vectorPol[b]
            err = vectorPolError[b]
            if err == 0:
                continue
            Atau, Ae = par[0], par[1]
            exp = -(Atau*(1+ctheta**2)+2*Ae*ctheta)/(1+ctheta**2+2*Ae*Atau*ctheta)
            val += (obs-exp)**2 / err**2
        f.value = val

    minuit2 = ROOT.TMinuit(2)
    minuit2.SetFCN(my_minimization2)
    minuit2.DefineParameter(0,"Atau",0.14,1e-5,0,0.2)
    minuit2.DefineParameter(1,"Ae",0.16,1e-5,0,0.2)
    minuit2.Command("MIGRAD")

    pval = ctypes.c_double(0)
    perr = ctypes.c_double(0)

    minuit2.GetParameter(0,pval,perr)
    Atau, Atau_err = pval.value, perr.value

    minuit2.GetParameter(1,pval,perr)
    Ae, Ae_err = pval.value, perr.value

    gv_ga = (1 - math.sqrt(1 - Ae*Ae)) / Ae
    term1 = (Ae*Ae) / math.sqrt(1 - Ae*Ae)
    term2 = 1 - math.sqrt(1 - Ae*Ae)
    dgv_dAe = (term1 - term2) / (Ae*Ae)
    gv_ga_err = abs(dgv_dAe) * Ae_err

    sin2theta_eff = (1 - gv_ga) / 4
    sin2theta_eff_err = gv_ga_err / 4
    
    summary_txt = os.path.join(outdir, f"asymmetry_results{tag}.txt")

    with open(summary_txt, "w") as f:
        f.write("# Asymmetry fit results\n")
        f.write(f"Atau        {Atau:.6f}  {Atau_err:.6f}\n")
        f.write(f"Ae          {Ae:.6f}  {Ae_err:.6f}\n")
        f.write(f"gv/ga       {gv_ga:.6f}  {gv_ga_err:.6f}\n")
        f.write(f"sin2theta_eff  {sin2theta_eff:.6f}  {sin2theta_eff_err:.6f}\n")

    table_txt = os.path.join(outdir, f"fit_parameters_per_bin{tag}.txt")

    with open(table_txt, "w") as f:
        f.write(
            "# bin  cos_min  cos_max  "
            "scale_sig  err  "
            "pol  err  "
            "scale_bg  err  "
            "chi2\n"
        )
        for i in range(len(fitTable["bin"])):
            f.write(
                f"{fitTable['bin'][i]:2d}  "
                f"{fitTable['cos_min'][i]: .4f}  "
                f"{fitTable['cos_max'][i]: .4f}  "
                f"{fitTable['scale_signal'][i]: .5f}  "
                f"{fitTable['scale_signal_err'][i]: .5f}  "
                f"{fitTable['pol'][i]: .6f}  "
                f"{fitTable['pol_err'][i]: .6f}  "
                f"{fitTable['scale_bg'][i]: .5f}  "
                f"{fitTable['scale_bg_err'][i]: .5f}  "
                # f"{fitTable['chi2'][i]: .2f}\n"
            )

    
    graph = ROOT.TGraphErrors()

    for bin in range(nBins):
        x = 0.5 * (vectorCosThetaMin[bin] + vectorCosThetaMax[bin])
        xerr = 0.5 * (vectorCosThetaMax[bin] - vectorCosThetaMin[bin])
        y = vectorPol[bin]
        yerr = vectorPolError[bin]

        graph.SetPoint(bin, x, y)
        graph.SetPointError(bin, xerr, yerr)

    graph.SetName("P_tau_vs_costheta")
    graph.GetXaxis().SetTitle("cos #theta_{#tau}")
    graph.GetYaxis().SetTitle("P_{#tau}")
    # title axis size
    graph.GetXaxis().SetTitleSize(0.05)
    graph.GetYaxis().SetTitleSize(0.05)
    
    fout.cd()
    graph.Write()
    c = ROOT.TCanvas("c_pol", "", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)
    graph.SetMarkerStyle(20)
    graph.Draw("AP")
    
    fit_func = ROOT.TF1(
    "fit_func",
    "-([0]*(1+x*x)+2*[1]*x)/(1+x*x+2*[1]*[0]*x)",
    -1, 1
)
    fit_func.SetParameters(Atau, Ae)
    fit_func.SetLineColor(ROOT.kRed)
    fit_func.SetLineWidth(2)
    fit_func.Draw("SAME")
    
    leg = ROOT.TLegend(0.5,0.89,0.85,0.65)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetLineWidth(0)
    leg.AddEntry(graph,"Pseudo Data","pl")
    leg.AddEntry(fit_func,"Fit for A_{#tau} and A_{e}","l")
    # Añadimos extra_leg
    if extra_legend:
        for extra_text in extra_legend:
            leg.AddEntry(ROOT.nullptr, extra_text, "")
    leg.Draw()

    png = os.path.join(images_outdir, f"polarization_vs_costheta{tag}.png")
    pdf = os.path.join(images_outdir, f"polarization_vs_costheta{tag}.pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)
                #     f"{fitTable['bin'][i]:2d}  "
                # f"{fitTable['cos_min'][i]: .4f}  "
                # f"{fitTable['cos_max'][i]: .4f}  "
                # f"{fitTable['scale_signal'][i]: .5f}  "
                # f"{fitTable['scale_signal_err'][i]: .5f}  "
                # f"{fitTable['pol'][i]: .6f}  "
                # f"{fitTable['pol_err'][i]: .6f}  "
                # f"{fitTable['scale_bg'][i]: .5f}  "

        
    import pandas as pd 
    pol_results_df = pd.DataFrame(fitTable)
    # Redondear todas las columnas a 3 decimales
    pol_results_df = pol_results_df.round(3) 
    pol_results_df.to_csv(os.path.join(outdir, f"polarization_results{tag}.csv"), index=False)
    # txt = os.path.join(outdir, f"polarization_results{tag}.txt")
    # with open(txt, "w") as f:
    #     f.write("# bin  cosTheta_min  cosTheta_max  P_tau  P_tau_err\n")
    #     for b in range(nBins):
    #         f.write(
    #             f"{b:2d}  "
    #             f"{vectorCosThetaMin[b]: .4f}  "
    #             f"{vectorCosThetaMax[b]: .4f}  "
    #             f"{vectorPol[b]: .6f}  "
    #             f"{vectorPolError[b]: .6f}\n"
    #         )

    file.Close()
    fout.Close()
    if verbose:
        print("[OK] Polarization fit results saved.")



# =====================================================
# CLI
# =====================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser("Polarization fit with scalable BG")

    parser.add_argument("-i","--input",default="Binned_histograms/BINED_template_Histos_dRgt3.06_5.0_MesonPgt0.0_lt45.57_LeptonPgt0.0_lt41.16_tau_traineddecay2_0.4_tph0.35_tpi0_n3_g0.0.root",
                        help="Input BINED ROOT file")
    parser.add_argument("--nBins", type=int, default=20)
    parser.add_argument("--rebin", type=int, default=2)
    parser.add_argument("--perfect", action="store_true")
    parser.add_argument("--chi2", action="store_true")
    parser.add_argument("-o", "--outdir", default="./Binned_histograms/")
    parser.add_argument("--no-bg", action="store_true", help="Fit without background component")
    parser.add_argument("--bg-mode", choices=["total", "mig", "ext", "split"], default="total",
                        help="Background mode: total (BG completo), mig (solo migraciones), ext (solo externo), split (mig+ext con scales separados)")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--extra-legend", nargs="+", default=["7 fb^{-1}\\text{  }\\sqrt{s}=91GeV"],
                        help="Texto(s) extra para mostrar en la leyenda (informativo, sin gráfico asociado)")
    parser.add_argument("--lumi-base", type=float, default=None,
                        help="Luminosidad base a la que están escalados los datos (fb^-1)")
    parser.add_argument("--lumi-target", type=float, default=None,
                        help="Luminosidad objetivo para reescalar señal y fondos (fb^-1)")

    args = parser.parse_args()

    main(
        input_file=args.input,
        nBins=args.nBins,
        rebin=args.rebin,
        force_perfect_agreement=args.perfect,
        use_likelihood=not args.chi2,
        outdir=args.outdir,
        verbose=args.verbose,
        no_bg=args.no_bg,
        bg_mode=args.bg_mode,
        extra_legend=args.extra_legend,
        lumi_base=args.lumi_base,
        lumi_target=args.lumi_target,
    )