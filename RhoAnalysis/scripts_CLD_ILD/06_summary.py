#!/usr/bin/env python3
"""06_summary.py — collect the CLD vs ILD polarization results into one table.

Scans the outputs of stages 3 and 4 and writes a CSV + a markdown table with,
per detector and channel:
  * A_tau and A_e from `asymmetry_results_def.txt`, at the MC-equivalent
    luminosity, at the common 2M-tau normalization, and extrapolated to one
    FCC-ee Z-pole year
  * the optimized cuts and their S, B, efficiency (from `optimization_results.csv`)
"""
import argparse
import csv
import glob
import os

REPO = "/nfs/cms/arqolmo/TausFCCee"
DETECTORS = ["CLD", "ILD", "CLDgen", "ILDgen"]
# Canales de la etapa 4 (ver 04_binning_and_fit.py): rho + lepton explicito y
# los inclusivos rho/pi/e/mu. a1 y el id 1 quedan fuera por ahora.
CHANNELS = ["rho_lep", "rho", "pion", "ele", "muon"]
# Canales con cortes PSO (el muon aun no se optimiza).
CUT_CHANNELS = ["pion", "rho", "ele"]


def read_asym(path):
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 3:
            out[parts[0]] = (float(parts[1]), float(parts[2]))
    return out


def read_cuts(det, ch):
    path = os.path.join(REPO, "Results/CutOptimization_CLD_ILD",
                        f"{det}_{ch}", "optimization_results.csv")
    if not os.path.isfile(path):
        return {}
    with open(path) as f:
        return next(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="+", default=["nocuts"])
    ap.add_argument("--binbase", default="Binned_histograms_MDecs/CLD_ILD")
    ap.add_argument("-o", "--outdir", default="Results/CLD_ILD_Summary")
    args = ap.parse_args()

    outdir = os.path.join(REPO, args.outdir)
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for variant in args.variants:
        for det in DETECTORS:
            for ch in CHANNELS:
                base = os.path.join(REPO, args.binbase, f"{det}_{ch}_{variant}")
                for scope, sub in (("MC stat", "fit_MCstat"),
                                   ("2M generated tautau", "fit_2Mtau"),
                                   ("7 ab^-1 (1 FCC-ee Z year)", "fit_FCCyear")):
                    a = read_asym(os.path.join(base, sub, "asymmetry_results_def.txt"))
                    if not a:
                        continue
                    rows.append({
                        "variant": variant, "detector": det, "channel": ch,
                        "lumi_scope": scope,
                        "Atau": a.get("Atau", ("", ""))[0],
                        "Atau_err": a.get("Atau", ("", ""))[1],
                        "Ae": a.get("Ae", ("", ""))[0],
                        "Ae_err": a.get("Ae", ("", ""))[1],
                        "sin2theta_eff": a.get("sin2theta_eff", ("", ""))[0],
                        "sin2theta_eff_err": a.get("sin2theta_eff", ("", ""))[1],
                    })

    csv_path = os.path.join(outdir, "asymmetry_summary.csv")
    if rows:
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    cut_rows = []
    for det in DETECTORS:
        for ch in CUT_CHANNELS:
            c = read_cuts(det, ch)
            if c:
                cut_rows.append({"detector": det, "channel": ch, **c})
    cuts_csv = os.path.join(outdir, "optimized_cuts.csv")
    if cut_rows:
        keys = sorted({k for r in cut_rows for k in r})
        keys = ["detector", "channel"] + [k for k in keys if k not in ("detector", "channel")]
        with open(cuts_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(cut_rows)

    md = ["# CLD vs ILD — tau polarization summary", "",
          "## Asymmetry fits", "",
          "| variant | detector | channel | luminosity | A_tau | A_e |",
          "|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['variant']} | {r['detector']} | {r['channel']} | {r['lumi_scope']} | "
                  f"{r['Atau']:.4f} ± {r['Atau_err']:.4f} | {r['Ae']:.4f} ± {r['Ae_err']:.4f} |")
    if cut_rows:
        md += ["", "## Optimized cuts (PSO)", "",
               "| detector | channel | dR | meson P [GeV] | lepton P [GeV] | m_Z [GeV] | eff_S | S | B |",
               "|---|---|---|---|---|---|---|---|---|"]
        for c in cut_rows:
            md.append(
                f"| {c['detector']} | {c['channel']} | "
                f"{float(c['dR_min']):.2f}–{float(c['dR_max']):.2f} | "
                f"{float(c['mesonP_min']):.2f}–{float(c['mesonP_max']):.2f} | "
                f"{float(c['lepP_min']):.2f}–{float(c['lepP_max']):.2f} | "
                f"{float(c['Zmass_min']):.2f}–{float(c['Zmass_max']):.2f} | "
                f"{float(c['effS']):.3f} | {float(c['S']):.0f} | {float(c['B']):.0f} |")
    md_path = os.path.join(outdir, "summary.md")
    open(md_path, "w").write("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\n[OK] {csv_path}\n[OK] {cuts_csv}\n[OK] {md_path}")


if __name__ == "__main__":
    main()
