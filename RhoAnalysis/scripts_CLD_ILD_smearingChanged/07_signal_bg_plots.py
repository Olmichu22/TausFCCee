#!/usr/bin/env python3
"""07_signal_bg_plots.py — señal y fondo por canal, CLD vs ILD (smearing corregido).

Dos histogramas de barras sobre los cuatro canales de medida de la polarización:

  pi + X (cortes optimos), rho + lepton (cortes legacy), e + X y mu + X (sin cortes)

  1. Numero bruto de eventos de senal (SS)
  2. Porcentaje de fondo (migraciones) sobre el total

Los conteos se leen de las carpetas `fit_MCstat` (estadistica MC real) y se escalan
aqui a 2M eventos Z->tautau generados con el factor 2.0e6/N_gen (ILD: 1.59; CLD: 1),
de modo que ambas muestras son directamente comparables. Los conteos no dependen del
sin^2(theta)_eff, asi que se leen los de 0.2312.

OJO: no se usan los event_counts_def.txt de `fit_2Mtau`/`fit_FCCyear` porque
fitPolAssym.py aplica alli el factor de luminosidad dos veces (escala los
histogramas y vuelve a multiplicar la integral por lumi_scale), de modo que su
columna de eventos sale multiplicada por lumi_scale^2. Los porcentajes y los
resultados del fit no estan afectados.

  python RhoAnalysis/scripts_CLD_ILD_smearingChanged/07_signal_bg_plots.py
"""
import argparse
import csv
import glob
import os

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

REPO = "/nfs/cms/arqolmo/TausFCCee"
N_TAU_COMMON = 2.0e6
DET_DIR = {
    "CLD": "Results/RhoAnalysis/PolAnalysis_RECO_CLD_ztt2M_smearingChanged_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    "ILD": "Results/RhoAnalysis/PolAnalysis_RECO_ILD_fcc_smearingChanged_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
}
# Azul/naranja de la paleta Okabe-Ito (segura para daltonismo)
COLORS = {"CLD": "#0072B2", "ILD": "#E69F00"}
DETECTORS = ["CLD", "ILD"]
# modo -> etiqueta del eje X
MODES = [
    ("pion_optcuts",       "#pi + X (opt cuts)"),
    ("rho_lep_legacycuts", "#rho + lep (legacy)"),
    ("ele_nocuts",         "e + X (no cuts)"),
    ("muon_nocuts",        "#mu + X (no cuts)"),
]


def read_counts(path):
    """Devuelve (senal, fondo, total) de un event_counts_def.txt."""
    sig = bg = tot = None
    for line in open(path):
        parts = line.split()
        if line.startswith("Signal"):
            sig = float(parts[2])
        elif line.startswith("BG migrations"):
            bg = float(parts[2])
        elif line.startswith("Total"):
            tot = float(parts[1])
    if sig is None or bg is None or tot is None:
        raise SystemExit(f"No se pudo leer {path}")
    return sig, bg, tot


def scale_to_2M(det):
    """Factor N_gen -> 2M eventos generados, leido del results_summary del arbol."""
    hits = glob.glob(os.path.join(REPO, DET_DIR[det], "results_summary_*.csv"))
    if not hits:
        raise SystemExit(f"No results_summary_*.csv para {det}")
    with open(hits[0]) as f:
        ngen = float(next(csv.DictReader(f))["TotalEvents"])
    return N_TAU_COMMON / ngen


def make_hist(name, values, color):
    h = ROOT.TH1F(name, "", len(MODES), 0, len(MODES))
    for i, (_, label) in enumerate(MODES, start=1):
        h.SetBinContent(i, values[i - 1])
        h.SetBinError(i, 0.0)
        h.GetXaxis().SetBinLabel(i, label)
    col = ROOT.TColor.GetColor(color)
    h.SetFillColor(col)
    h.SetLineColor(col)
    h.SetLineWidth(1)
    h.SetBarWidth(0.34)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleOffset(1.35)
    return h


def draw(values, ytitle, fmt, outstem, title, subtitle):
    c = ROOT.TCanvas("c_" + outstem, "", 900, 600)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.09)
    c.SetBottomMargin(0.12)

    ymax = max(max(values[d]) for d in DETECTORS)
    hists = []
    for j, det in enumerate(DETECTORS):
        h = make_hist(f"h_{outstem}_{det}", values[det], COLORS[det])
        h.SetBarOffset(0.16 + 0.34 * j)
        h.SetMaximum(1.28 * ymax)
        h.SetMinimum(0.0)
        h.GetYaxis().SetTitle(ytitle)
        h.Draw("bar" if j == 0 else "bar same")
        hists.append(h)

    # Etiqueta directa sobre cada barra (evita leer el valor del eje)
    lat = ROOT.TLatex()
    lat.SetTextSize(0.030)
    lat.SetTextAlign(21)
    for j, det in enumerate(DETECTORS):
        lat.SetTextColor(ROOT.TColor.GetColor(COLORS[det]))
        for i, v in enumerate(values[det]):
            lat.DrawLatex(i + 0.33 + 0.34 * j, v + 0.025 * ymax, fmt(v))

    # Arriba a la izquierda: en ambos plots las barras crecen hacia la derecha
    leg = ROOT.TLegend(0.17, 0.74, 0.37, 0.89)
    leg.SetTextSize(0.034)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    for det, h in zip(DETECTORS, hists):
        leg.AddEntry(h, det, "f")
    leg.Draw()

    txt = ROOT.TLatex()
    txt.SetNDC()
    txt.SetTextSize(0.038)
    txt.SetTextFont(62)
    txt.DrawLatex(0.13, 0.935, title)
    txt.SetTextFont(42)
    txt.SetTextSize(0.030)
    txt.SetTextAlign(31)
    txt.SetTextColor(ROOT.kGray + 2)
    txt.DrawLatex(0.96, 0.935, subtitle)

    outdir = os.path.join(REPO, "Results/CLD_ILD_smearingChanged_Summary/plots")
    os.makedirs(outdir, exist_ok=True)
    for ext in ("png", "pdf"):
        c.SaveAs(os.path.join(outdir, f"{outstem}.{ext}"))
    return os.path.join(outdir, outstem + ".png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sin-eff", default="0.2312",
                    help="Carpeta de la que se leen los conteos (son iguales para ambos)")
    ap.add_argument("--binbase",
                    default="Binned_histograms_MDecs/CLD_ILD_smearingChanged_sineff{sin_eff}")
    args = ap.parse_args()
    base = os.path.join(REPO, args.binbase.format(sin_eff=args.sin_eff))

    signal = {d: [] for d in DETECTORS}
    bgfrac = {d: [] for d in DETECTORS}
    factor = {d: scale_to_2M(d) for d in DETECTORS}
    for det in DETECTORS:
        print(f"[INFO] {det}: scale factor to 2M = {factor[det]:.4f}")
    print(f"{'channel':22s} {'det':4s} {'signal(MC)':>12s} {'signal(2M)':>12s} {'BG %':>7s}")
    for mode, _ in MODES:
        for det in DETECTORS:
            sig, bg, tot = read_counts(
                os.path.join(base, f"{det}_{mode}", "fit_MCstat", "event_counts_def.txt"))
            signal[det].append(sig * factor[det])
            bgfrac[det].append(100.0 * bg / tot)
            print(f"{mode:22s} {det:4s} {sig:12.1f} {sig * factor[det]:12.1f} "
                  f"{100.0 * bg / tot:6.2f}%")

    signal_k = {d: [v / 1000.0 for v in signal[d]] for d in DETECTORS}
    sub = "Scaled to 2M generated Z#rightarrow#tau#tau"
    p1 = draw(signal_k, "Signal events (SS)  [#times10^{3}]", lambda v: f"{v:.0f}k",
              "signal_events_2Mtau", "Selected signal events per channel", sub)
    p2 = draw(bgfrac, "Background (migrations)  [%]", lambda v: f"{v:.1f}%",
              "bg_fraction_2Mtau", "Background fraction per channel", sub)
    print(f"\n[OK] {p1}\n[OK] {p2}")


if __name__ == "__main__":
    main()
