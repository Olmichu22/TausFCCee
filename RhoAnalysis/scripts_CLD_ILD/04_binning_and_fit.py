#!/usr/bin/env python3
"""04_binning_and_fit.py — cos(theta) binning + A_tau fit for the CLD vs ILD study.

For every detector (CLD, ILD) and every analysed channel — rho against the
leptonic hemispheres (rho_lep) and the inclusive rho, pion, electron and muon
channels, whose other hemisphere is unrestricted — it runs:

  1. `makeCosBins_MDecs.py` with the CORRELATED weight templates
     (`--weights corr`, i.e. the two-hemisphere formula of eq. 9/14/16) and a
     luminosity weight of exactly 1 per event (`--signal-lumi-pb ngen/xsec`),
     so the binned templates carry the raw MC statistics.
  2. `fitPolAssym.py` three times (`--lumi-base` → `--lumi-target`):
     at the MC-equivalent luminosity, at a common normalization of 2M generated
     Z->tautau events, and rescaled to one FCC-ee Z-pole year.

The MC-equivalent luminosity is L = N_gen / sigma(Z->tautau), with N_gen read
from the `results_summary_*.csv` written by the tree stage.

The 2M-tau pass puts both detectors on the same generated statistics without
extrapolating: CLD is simulated with exactly 2,000,000 events (factor 1.000) and
ILD with 1,918,000 (factor 1.043), so the errors it quotes are the simulated
ones, not a sqrt(L) projection.

Usage:
    python RhoAnalysis/scripts_CLD_ILD/04_binning_and_fit.py [--variant nocuts]
    python RhoAnalysis/scripts_CLD_ILD/04_binning_and_fit.py --fits-only
"""
import argparse
import csv
import glob
import os
import subprocess
import sys

REPO = "/nfs/cms/arqolmo/TausFCCee"
XSEC_PB = 1476.58        # sigma(e+e- -> Z -> tau tau) at the Z pole, as in makeCosBins
LUMI_FCC_YEAR_FB = 7000.  # 7 ab^-1: target for one FCC-ee Z-pole year
N_TAU_COMMON = 2.0e6      # common normalization: generated Z->tautau events

DET_DIR = {
    "CLD": "Results/RhoAnalysis/PolAnalysis_RECO_CLD_ztt2M_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    "ILD": "Results/RhoAnalysis/PolAnalysis_RECO_ILD_fcc_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    # Gen-level references (perfect detector): same machinery on the gen trees.
    "CLDgen": "Results/RhoAnalysis/PolAnalysis_GEN_CLD_ztt2M_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    "ILDgen": "Results/RhoAnalysis/PolAnalysis_GEN_ILD_fcc_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
}
# Canales analizados, como (target-decay, other-decay) de makeCosBins_MDecs.py.
# Los inclusivos ("all") se leen de los ficheros single-decay
# HistosMDecs_single{id}_*, y rho_lep de los pares rho x e / rho x mu.
# a1 y el id 1 (rho con un solo foton) quedan fuera por ahora.
CHANNELS = {
    "rho_lep": ("rho",  "lep"),
    "rho":     ("rho",  "all"),
    "pion":    ("pion", "all"),
    "ele":     ("ele",  "all"),
    "muon":    ("muon", "all"),
}
DECAY_ID = {"rho": 2, "pion": 0, "ele": -11, "muon": -13}


def n_generated(sample_dir):
    """Total generated events of the sample, from the tree-stage summary CSV."""
    hits = glob.glob(os.path.join(REPO, sample_dir, "results_summary_*.csv"))
    if not hits:
        raise SystemExit(f"No results_summary_*.csv in {sample_dir}")
    with open(hits[0]) as f:
        row = next(csv.DictReader(f))
    return float(row["TotalEvents"])


def is_cut_stem(stem):
    """True if the stem encodes kinematic cuts (see the out_prefix rules of
    RhoHistFromTree_MDecs_parallel.py)."""
    return any(tok in stem for tok in ("dRgt", "Dec0Pgt", "Dec1Pgt", "Zmassgt"))


def find_stem(sample_dir, target, other, variant, lepton_id=-11):
    """Stem shared by all HistosMDecs files of a channel (cuts are encoded in it).

    `variant` selects between the no-cut and the optimized-cut production, which
    live side by side in the same directory. An inclusive channel (`other`
    == "all") is served by the single-decay file of the target; otherwise the
    pair file target x lepton is used.
    """
    if other == "all":
        prefix = f"HistosMDecs_single{DECAY_ID[target]}_"
    else:
        prefix = f"HistosMDecs_{DECAY_ID[target]}_{lepton_id}_"
    pat = os.path.join(REPO, sample_dir, prefix + "*.root")
    want_cuts = variant != "nocuts"
    hits = [p for p in sorted(glob.glob(pat), key=os.path.getmtime)
            if is_cut_stem(os.path.basename(p)) == want_cuts]
    if not hits:
        print(f"  [WARN] no {variant} histogram file matching {pat}")
        return None
    base = os.path.basename(hits[-1])[:-len(".root")]
    return base[len(prefix):]


def run(cmd, log):
    print("  $ " + " ".join(cmd), flush=True)
    with open(log, "a") as fh:
        fh.write("\n$ " + " ".join(cmd) + "\n")
        fh.flush()
        rc = subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc:
        print(f"  [FAIL rc={rc}] see {log}", flush=True)
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="nocuts",
                    help="Label of the histogram variant (nocuts / optcuts)")
    ap.add_argument("--outbase", default="Binned_histograms_MDecs/CLD_ILD")
    ap.add_argument("--detectors", nargs="+", default=["CLD", "ILD"])
    ap.add_argument("--channels", nargs="+", default=list(CHANNELS))
    ap.add_argument("--fits-only", action="store_true",
                    help="Reuse the existing BINED_*.root and rerun only the fits")
    args = ap.parse_args()

    lumi_common_fb = N_TAU_COMMON / XSEC_PB / 1000.0

    logdir = os.path.join(REPO, "logs/CLD_ILD")
    os.makedirs(logdir, exist_ok=True)

    for det in args.detectors:
        sample_dir = DET_DIR[det]
        ngen = n_generated(sample_dir)
        lumi_pb = ngen / XSEC_PB              # weight = lumi*xsec/ngen = 1
        lumi_fb = lumi_pb / 1000.0
        print(f"\n===== {det}: N_gen={ngen:.0f}, L_equiv={lumi_fb:.3f} fb^-1 =====")

        for ch in args.channels:
            target, other = CHANNELS[ch]
            outdir = os.path.join(args.outbase, f"{det}_{ch}_{args.variant}")
            log = os.path.join(logdir, f"fit_{det}_{ch}_{args.variant}.log")
            open(log, "w" if not args.fits_only else "a").close()
            # En modo single-decay el otro hemisferio se clasifica siempre como
            # fondo (_classify_hemisphere con expected_gen_id=None), asi que
            # SIGNAL_SIGNAL esta vacio y la senal inclusiva vive en SIGNAL_BG:
            # es justo lo que suma bg_def='rho_only'. Los pares siguen con 'ss'.
            bg_def = "ss" if other == "lep" else "rho_only"
            # Nombre que escribe makeCosBins_MDecs.py: objetivo + tag del otro
            # hemisferio, NO el nombre del canal (rho_lep → rho + lep).
            binned = os.path.join(outdir,
                                  f"BINED_MDecs_{target}_{other}_corr_{bg_def}.root")

            if not args.fits_only:
                stem = find_stem(sample_dir, target, other, args.variant)
                if stem is None:
                    print(f"-- {det} {ch}: no {args.variant} histograms, skipped")
                    continue
                print(f"-- {det} {ch}: stem={stem}")
                if run([sys.executable, "RhoAnalysis/makeCosBins_MDecs.py",
                        "--sample-dir", sample_dir, "--stem", stem,
                        "--target-decay", target, "--other-decay", other,
                        "--weights", "corr", "--bg-def", bg_def,
                        "--signal-type", "Ztt",
                        "--signal-ngen", f"{ngen:.0f}",
                        "--signal-lumi-pb", f"{lumi_pb:.6f}",
                        "-o", outdir, "-v"], log):
                    continue
            elif not os.path.isfile(os.path.join(REPO, binned)):
                print(f"-- {det} {ch}: no {binned}, skipped")
                continue

            # (a) MC-equivalent luminosity
            run([sys.executable, "RhoAnalysis/fitPolAssym.py", "-i", binned,
                 "-o", os.path.join(outdir, "fit_MCstat"), "--bg-mode", "total", "-v",
                 "--extra-legend", f"{det}  {lumi_fb:.2f} fb^{{-1}} (MC stat)"], log)
            # (b) both detectors at the same 2M generated taus, no extrapolation
            run([sys.executable, "RhoAnalysis/fitPolAssym.py", "-i", binned,
                 "-o", os.path.join(outdir, "fit_2Mtau"), "--bg-mode", "total", "-v",
                 "--lumi-base", f"{lumi_fb:.6f}", "--lumi-target", f"{lumi_common_fb:.6f}",
                 "--extra-legend", f"{det}  2M #tau#tau (x{N_TAU_COMMON / ngen:.3f})"], log)
            # (c) one FCC-ee Z-pole year
            run([sys.executable, "RhoAnalysis/fitPolAssym.py", "-i", binned,
                 "-o", os.path.join(outdir, "fit_FCCyear"), "--bg-mode", "total", "-v",
                 "--lumi-base", f"{lumi_fb:.6f}", "--lumi-target", f"{LUMI_FCC_YEAR_FB}",
                 "--extra-legend", f"{det}  7 ab^{{-1}} (1 FCC-ee Z year)"], log)


if __name__ == "__main__":
    main()
