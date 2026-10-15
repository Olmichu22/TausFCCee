#!/usr/bin/env python3
"""04_binning_and_fit.py — cos(theta) binning + A_tau fit, TauTool vs ParTauDETR (CLD).

Same recipe as scripts_CLD_ILD_smearingChanged/04_binning_and_fit.py, on the
rho-lepton channel of the two CLD reconstructions of ztt_2M_smearing:

  TauTool : cone taus (PolAnalysis_RECO_CLD_ztt2M_smearingChanged_*, 2M events)
  MLTau   : ParTauDETR taus (PolAnalysis_RECO_CLD_ztt2M_smearingChanged_MLTau_*, 2M events)

with and without the legacy cuts. Output:

  Binned_histograms_MDecs/CLD_TauTool_vs_MLTau_sineff<S>/<ALG>_rho_lep_<cuts>/

With --bg-def other than the default "ss" (signal = both hemispheres gen-correct)
the folder gets the definition as a suffix, e.g. <ALG>_rho_lep_<cuts>_rho_only/
(signal = rho hemisphere gen-correct, whatever the lepton), so nothing is
overwritten.

Per folder: makeCosBins_MDecs.py (correlated weights) and fitPolAssym.py at MC
luminosity, at 2M generated taus and at one FCC-ee year (7 ab^-1). Both
algorithms run on the same 2M events (files 1-2000).
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
LEGEND = ["--legend-fit", "0.62", "0.58", "0.98", "0.92",
          "--legend-pol", "0.62", "0.68", "0.98", "0.92",
          "--legend-size", "0.030"]

ALG_DIR = {
    "TauTool": "Results/RhoAnalysis/PolAnalysis_RECO_CLD_ztt2M_smearingChanged_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    "MLTau": "Results/RhoAnalysis/PolAnalysis_RECO_CLD_ztt2M_smearingChanged_MLTau_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
}
LEGACY_CUTS = "dRgt3.06_5.0_Dec0Pgt2.0_lt45.57_Dec1Pgt10.0_lt41.16_cosAcc0.95_"
# modo -> tokens de corte en el nombre del fichero de histogramas
MODES = {"rho_lep_legacycuts": LEGACY_CUTS, "rho_lep_nocuts": ""}


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


def process(alg, mode, sin_eff, outbase, logdir, bg_def="ss"):
    sample_dir = ALG_DIR[alg]
    stem = f"zTaum_{MODES[mode]}sineff{sin_eff}_vismOff_{TAG}"
    for lep in ("-11", "-13"):
        hist = os.path.join(REPO, sample_dir, f"HistosMDecs_2_{lep}_{stem}.root")
        if not os.path.isfile(hist):
            return alg, mode, f"MISSING {hist}"

    ngen = n_generated(sample_dir)
    lumi_pb = ngen / XSEC_PB
    lumi_fb = lumi_pb / 1000.0
    lumi_common_fb = N_TAU_COMMON / XSEC_PB / 1000.0

    # "ss" conserva el nombre de carpeta original; el resto lleva sufijo
    tag = f"{alg}_{mode}" + ("" if bg_def == "ss" else f"_{bg_def}")
    outdir = os.path.join(outbase.format(sin_eff=sin_eff), tag)
    os.makedirs(os.path.join(REPO, outdir), exist_ok=True)
    log = os.path.join(logdir, f"fit_{tag}_sineff{sin_eff}.log")
    open(log, "w").close()
    py = sys.executable

    binned = os.path.join(outdir, f"BINED_MDecs_rho_lep_corr_{bg_def}.root")
    rc = run([py, "RhoAnalysis/makeCosBins_MDecs.py",
              "--sample-dir", sample_dir, "--stem", stem,
              "--target-decay", "rho", "--other-decay", "lep",
              "--weights", "corr", "--bg-def", bg_def,
              "--signal-type", "Ztt",
              "--signal-ngen", f"{ngen:.0f}", "--signal-lumi-pb", f"{lumi_pb:.6f}",
              "-o", outdir, "-v"], log)
    if rc:
        return alg, mode, f"FAIL makeCosBins rc={rc} (ver {log})"
    fit = [py, "RhoAnalysis/fitPolAssym.py", "-i", binned,
           "--bg-mode", "total", "--rebin", "2", "-v", *LEGEND]
    rcs = [
        run(fit + ["-o", os.path.join(outdir, "fit_MCstat"),
                   "--extra-legend", f"CLD {alg}  {lumi_fb:.2f} fb^{{-1}} (MC stat)"], log),
        run(fit + ["-o", os.path.join(outdir, "fit_2Mtau"),
                   "--lumi-base", f"{lumi_fb:.6f}", "--lumi-target", f"{lumi_common_fb:.6f}",
                   "--extra-legend", f"CLD {alg}  {lumi_common_fb:.2f} fb^{{-1}} (2M #tau#tau)"], log),
        run(fit + ["-o", os.path.join(outdir, "fit_FCCyear"),
                   "--lumi-base", f"{lumi_fb:.6f}", "--lumi-target", f"{LUMI_FCC_YEAR_FB}",
                   "--extra-legend", f"CLD {alg}  {LUMI_FCC_YEAR_FB:.0f} fb^{{-1}} (FCC year)"], log),
    ]
    return alg, mode, "OK" if not any(rcs) else f"FAIL rcs={rcs} (ver {log})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sin-eff", default="0.2315")
    ap.add_argument("--algs", nargs="+", default=list(ALG_DIR))
    ap.add_argument("--modes", nargs="+", default=list(MODES))
    ap.add_argument("--outbase", default="Binned_histograms_MDecs/CLD_TauTool_vs_MLTau_sineff{sin_eff}")
    ap.add_argument("--bg-def", default="ss", choices=["ss", "sa", "rho_only"],
                    help="Background definition passed to makeCosBins_MDecs.py (default: ss)")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    logdir = os.path.join(REPO, "logs/CLD_TauTool_vs_MLTau")
    os.makedirs(logdir, exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(process, a, m, args.sin_eff, args.outbase, logdir, args.bg_def)
                for a in args.algs for m in args.modes]
        for fut in as_completed(futs):
            alg, mode, status = fut.result()
            print(f"[{status[:4]}] {alg}_{mode}: {status}", flush=True)


if __name__ == "__main__":
    main()
