#!/usr/bin/env python3
"""04_binning_and_fit.py — cos(theta) binning + A_tau fit, CLD vs ILD con smearing corregido.

Replica la etapa 4 de scripts_CLD_ILD (04_binning_and_fit.py, run_lep10.sh y la
receta legacycuts de docs/usefull_commands.md) sobre los árboles
PolAnalysis_RECO_{CLD_ztt2M,ILD_fcc}_smearingChanged_*, para los modos

  rho_lep_legacycuts, ele_nocuts, muon_nocuts, ele_lep10, muon_lep10, pion_optcuts

y para cada sin²θ_eff producido en la etapa 2 (0.2312 y 0.2315). Salida:

  Binned_histograms_MDecs/CLD_ILD_smearingChanged_sineff<S>/<DET>_<modo>/

Por carpeta: makeCosBins_MDecs.py (pesos SIEMPRE correlados, lumi-weight = 1) y
fitPolAssym.py a luminosidad MC, a 2M taus generados y a 1 año FCC-ee (7 ab^-1).

Las carpetas son independientes y se procesan en paralelo (--jobs, ≤ 20).
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

REPO = "/nfs/cms/arqolmo/TausFCCee"
XSEC_PB = 1476.58
LUMI_FCC_YEAR_FB = 7000.
N_TAU_COMMON = 2.0e6
TAG = "tau_traineddecayAll_0.4_tph0.35_tpi0_n3_g0.0"
# Leyenda pegada a la esquina superior derecha, con letra menor que la de por
# defecto (0.038), para que no caiga sobre los histogramas.
LEGEND = ["--legend-fit", "0.62", "0.58", "0.98", "0.92",
          "--legend-pol", "0.62", "0.68", "0.98", "0.92",
          "--legend-size", "0.030"]

DET_DIR = {
    "CLD": "Results/RhoAnalysis/PolAnalysis_RECO_CLD_ztt2M_smearingChanged_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    "ILD": "Results/RhoAnalysis/PolAnalysis_RECO_ILD_fcc_smearingChanged_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
}
# Cortes de la PSO original por detector (config/pipeline/CLD_ILD/<DET>_reco_optcuts.yaml)
PION_CUTS = {
    "CLD": "dRgt0.0135_5.9882_Dec0Pgt2.496_lt50.0_Zmassgt1.9033_lt100.0_",
    "ILD": "dRgt0.0308_6.0_Dec0Pgt2.4962_lt50.0_Zmassgt0.2667_lt100.0_",
}
LEGACY_CUTS = "dRgt3.06_5.0_Dec0Pgt2.0_lt45.57_Dec1Pgt10.0_lt41.16_cosAcc0.95_"

# modo -> (target, other, prefijo de fichero, tokens de corte (str o dict por det))
MODES = {
    "rho_lep_legacycuts": ("rho",  "lep", "HistosMDecs_2_-11_",     LEGACY_CUTS),
    "ele_nocuts":         ("ele",  "all", "HistosMDecs_single-11_", ""),
    "muon_nocuts":        ("muon", "all", "HistosMDecs_single-13_", ""),
    "ele_lep10":          ("ele",  "all", "HistosMDecs_single-11_", "Dec0Pgt10.0_lt100.0_"),
    "muon_lep10":         ("muon", "all", "HistosMDecs_single-13_", "Dec0Pgt10.0_lt100.0_"),
    "pion_optcuts":       ("pion", "all", "HistosMDecs_single0_",   PION_CUTS),
}


def n_generated(sample_dir):
    hits = glob.glob(os.path.join(REPO, sample_dir, "results_summary_*.csv"))
    if not hits:
        raise SystemExit(f"No results_summary_*.csv in {sample_dir}")
    with open(hits[0]) as f:
        return float(next(csv.DictReader(f))["TotalEvents"])


def run(cmd, log):
    with open(log, "a") as fh:
        fh.write("\n$ " + " ".join(cmd) + "\n")
        fh.flush()
        return subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT).returncode


def process(det, mode, sin_eff, outbase, logdir, fits_only=False):
    target, other, fprefix, cuts = MODES[mode]
    cuts = cuts[det] if isinstance(cuts, dict) else cuts
    sample_dir = DET_DIR[det]
    stem = f"zTaum_{cuts}sineff{sin_eff}_vismOff_{TAG}"
    hist = os.path.join(REPO, sample_dir, fprefix + stem + ".root")
    if not os.path.isfile(hist):
        return det, mode, sin_eff, f"MISSING {hist}"

    ngen = n_generated(sample_dir)
    lumi_pb = ngen / XSEC_PB
    lumi_fb = lumi_pb / 1000.0
    lumi_common_fb = N_TAU_COMMON / XSEC_PB / 1000.0

    outdir = os.path.join(outbase.format(sin_eff=sin_eff), f"{det}_{mode}")
    os.makedirs(os.path.join(REPO, outdir), exist_ok=True)
    log = os.path.join(logdir, f"fit_{det}_{mode}_sineff{sin_eff}.log")
    open(log, "a" if fits_only else "w").close()
    py = sys.executable

    # Pesos siempre correlados (joint two-tau, término cruzado).
    bg_def = "ss" if other == "lep" else "rho_only"
    binned = os.path.join(outdir, f"BINED_MDecs_{target}_{other}_corr_{bg_def}.root")
    if fits_only:
        if not os.path.isfile(os.path.join(REPO, binned)):
            return det, mode, sin_eff, f"MISSING {binned} (--fits-only)"
    else:
        rc = run([py, "RhoAnalysis/makeCosBins_MDecs.py",
                  "--sample-dir", sample_dir, "--stem", stem,
                  "--target-decay", target, "--other-decay", other,
                  "--weights", "corr", "--bg-def", bg_def,
                  "--signal-type", "Ztt",
                  "--signal-ngen", f"{ngen:.0f}", "--signal-lumi-pb", f"{lumi_pb:.6f}",
                  "-o", outdir, "-v"], log)
        if rc:
            return det, mode, sin_eff, f"FAIL makeCosBins rc={rc} (ver {log})"
    # Leyenda corta: detector + estadística en fb^-1 (origen entre paréntesis)
    rcs = [
        run([py, "RhoAnalysis/fitPolAssym.py", "-i", binned,
             "-o", os.path.join(outdir, "fit_MCstat"),
             "--bg-mode", "total", "--rebin", "2", "-v", *LEGEND,
             "--extra-legend", f"{det}  {lumi_fb:.2f} fb^{{-1}} (MC stat)"], log),
        run([py, "RhoAnalysis/fitPolAssym.py", "-i", binned,
             "-o", os.path.join(outdir, "fit_2Mtau"),
             "--bg-mode", "total", "--rebin", "2", "-v", *LEGEND,
             "--lumi-base", f"{lumi_fb:.6f}", "--lumi-target", f"{lumi_common_fb:.6f}",
             "--extra-legend", f"{det}  {lumi_common_fb:.2f} fb^{{-1}} (2M #tau#tau)"], log),
        run([py, "RhoAnalysis/fitPolAssym.py", "-i", binned,
             "-o", os.path.join(outdir, "fit_FCCyear"),
             "--bg-mode", "total", "--rebin", "2", "-v", *LEGEND,
             "--lumi-base", f"{lumi_fb:.6f}", "--lumi-target", f"{LUMI_FCC_YEAR_FB}",
             "--extra-legend", f"{det}  {LUMI_FCC_YEAR_FB:.0f} fb^{{-1}} (FCC year)"], log),
    ]
    status = "OK" if not any(rcs) else f"FAIL rcs={rcs} (ver {log})"
    return det, mode, sin_eff, status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sin-eff", nargs="+", default=["0.2312", "0.2315"])
    ap.add_argument("--detectors", nargs="+", default=["CLD", "ILD"])
    ap.add_argument("--modes", nargs="+", default=list(MODES))
    ap.add_argument("--outbase", default="Binned_histograms_MDecs/CLD_ILD_smearingChanged_sineff{sin_eff}")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--fits-only", action="store_true",
                    help="Reutiliza los BINED_*.root existentes y solo rehace los fits/plots")
    args = ap.parse_args()
    if args.jobs > 20:
        raise SystemExit("--jobs > 20 no permitido")

    logdir = os.path.join(REPO, "logs/CLD_ILD_smearingChanged")
    os.makedirs(logdir, exist_ok=True)
    tasks = [(d, m, s) for s in args.sin_eff for d in args.detectors for m in args.modes]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(process, d, m, s, args.outbase, logdir, args.fits_only)
                for d, m, s in tasks]
        for fut in as_completed(futs):
            det, mode, s, status = fut.result()
            print(f"[{status[:4]}] sineff={s} {det}_{mode}: {status}", flush=True)


if __name__ == "__main__":
    main()
