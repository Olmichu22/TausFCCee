#!/usr/bin/env python3
"""05_compare_plots.py — CLD vs ILD comparison figures for the polarization study.

For each analysed channel — rho against the leptonic hemisphere (rho_lep) and
the inclusive rho, pion, electron and muon channels — it:

  1. collects the histogram files of the channel (`hadd` of the e and mu pair
     files for rho_lep; the single-decay file for the inclusive ones);
  2. writes CompareAlgs YAML configs into `config/plots/CLD_ILD/`;
  3. runs `TauAnalysis/CompareAlgs.py` on each of them.

Figures produced per channel:
  * `optvar_all` — optimal variable, the four datasets together
    (CLD reco, ILD reco, CLD gen, ILD gen)
  * `optvar_reco` / `optvar_gen` — the same, reco-only and gen-only
  * `kinematics`  — visible momentum, reco Z mass, dR between hemispheres
  * `templates_CLD` / `templates_ILD` — nominal vs correlated P1/M1 reweighted
    templates, i.e. the input of the A_tau fit

Usage:
    python RhoAnalysis/scripts_CLD_ILD/05_compare_plots.py [--variant nocuts]
"""
import argparse
import glob
import os
import subprocess
import sys

REPO = "/nfs/cms/arqolmo/TausFCCee"
DIRS = {
    ("CLD", "reco"): "Results/RhoAnalysis/PolAnalysis_RECO_CLD_ztt2M_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    ("ILD", "reco"): "Results/RhoAnalysis/PolAnalysis_RECO_ILD_fcc_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    ("CLD", "gen"):  "Results/RhoAnalysis/PolAnalysis_GEN_CLD_ztt2M_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
    ("ILD", "gen"):  "Results/RhoAnalysis/PolAnalysis_GEN_ILD_fcc_tau_trained0.4_tph0.35_tpi0_n3_g0.0",
}
# Canales analizados, como (target-decay, other-decay). Los inclusivos ("all")
# se leen del fichero single-decay del objetivo; rho_lep, de los pares rho x l.
# a1 y el id 1 (rho con un solo foton) quedan fuera por ahora.
CHANNELS = {
    "rho_lep": ("rho",  "lep"),
    "rho":     ("rho",  "all"),
    "pion":    ("pion", "all"),
    "ele":     ("ele",  "all"),
    "muon":    ("muon", "all"),
}
DECAY_ID = {"rho": 2, "pion": 0, "ele": -11, "muon": -13}
TITLES = {"rho_lep": r"\rho^{\pm}", "rho": r"\rho^{\pm}", "pion": r"\pi^{\pm}",
          "ele": r"e^{\pm}", "muon": r"\mu^{\pm}"}
# Etiqueta del hemisferio acompanante segun el --other-decay del canal.
PARTNER = {"lep": r"\ell", "all": "X"}
# Name of the optimal observable per channel: the full omega for the rho and
# the analytic z_R (energy fraction) for the pion and the leptons.
OBS_LABEL = {"rho_lep": r"\omega_{\rho}", "rho": r"\omega_{\rho}",
             "pion": r"z_{\pi}", "ele": r"z_{e}", "muon": r"z_{\mu}"}
COLORS = {("CLD", "reco"): "ROOT.kRed", ("ILD", "reco"): "ROOT.kBlue",
          ("CLD", "gen"): "ROOT.kOrange+7", ("ILD", "gen"): "ROOT.kAzure+7"}
MARKERS = {("CLD", "reco"): 20, ("ILD", "reco"): 21,
           ("CLD", "gen"): 24, ("ILD", "gen"): 25}
# Optimal variable of the analysed hemisphere, PER CHANNEL — the same choice the
# fit makes in `_optimal_unified_reco` (X axis of OptimalReco_vs_CosThetaVis):
# the full omega for the rho (`Omega_Reco`, filled only for recoTauID 1/2) and
# x = 2E_vis/E_beam - 1 (`Optimal_X`) for the pion and the leptons. Using
# `Optimal_X` for the rho would show a flat, pion-like distribution instead of
# the true omega. In a gen-only tree the reco branches mirror the truth, so the
# same names give the generator-level distributions. The target hemisphere is
# always dec0, both in the pair runs and in the single-decay ones.
OPTVAR = {"rho_lep": "Omega_Reco_dec0", "rho": "Omega_Reco_dec0",
          "pion": "Optimal_X_dec0", "ele": "Optimal_X_dec0",
          "muon": "Optimal_X_dec0"}


def is_cut_stem(name):
    """True if the file name encodes kinematic cuts (out_prefix of the hist stage)."""
    return any(tok in name for tok in ("dRgt", "Dec0Pgt", "Dec1Pgt", "Zmassgt"))


def merge_channel(det, level, ch, variant, outdir):
    """hadd the files of one channel; return the merged path.

    An inclusive channel is a single file (the single-decay run); rho_lep is the
    hadd of its e and mu pair files.
    """
    src = os.path.join(REPO, DIRS[(det, level)])
    target, other = CHANNELS[ch]
    tid = DECAY_ID[target]
    if other == "all":
        patterns = [f"HistosMDecs_single{tid}_*.root"]
    else:
        patterns = [f"HistosMDecs_{tid}_{lep}_*.root" for lep in (-11, -13)]
    want_cuts = variant != "nocuts"
    files = []
    for pat in patterns:
        hits = [p for p in sorted(glob.glob(os.path.join(src, pat)),
                                  key=os.path.getmtime)
                if is_cut_stem(os.path.basename(p)) == want_cuts]
        if not hits:
            print(f"  [WARN] no file for {det}/{level} {ch} ({pat})")
            continue
        files.append(hits[-1])
    if not files:
        return None
    out = os.path.join(outdir, f"{det}_{level}_{ch}_{variant}.root")
    subprocess.run(["hadd", "-f", out] + files, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return out


def ds_block(label, path, det, level):
    return (f"  - path: {path}\n"
            f"    label: {label}\n"
            f"    color: {COLORS[(det, level)]}\n"
            f"    linestyle: 1\n"
            f"    markerstyle: {MARKERS[(det, level)]}\n"
            f"    markersize: 1.0\n"
            f"    linewidth: 3\n")


def _q(text):
    """Single-quoted YAML scalar: ROOT LaTeX is full of backslashes (illegal in
    double-quoted YAML) and of ': ' inside \\text{...} (illegal unquoted)."""
    text = str(text).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        text = text[1:-1]
    return "'" + text.replace("'", "''") + "'"


def plot_block(key, per_dataset, title, xlab, ylab="Norm. events",
               normalize="integral", extra=""):
    lines = [f"  {key}:"]
    lines.append("    per_dataset:")
    for lbl, name in per_dataset.items():
        lines.append(f"      {_q(lbl)}: {name}")
    lines.append(f"    title: {_q(title)}")
    lines.append(f"    x: {_q(xlab)}")
    lines.append(f"    y: {_q(ylab)}")
    lines.append(f'    normalize: "{normalize}"')
    lines.append("    legend: [0.60, 0.72, 0.90, 0.90]")
    if extra:
        lines.append(extra)
    return "\n".join(lines) + "\n"


def write_and_run(cfg_text, cfg_path, outdir, log):
    with open(cfg_path, "w") as f:
        f.write(cfg_text)
    cmd = [sys.executable, "TauAnalysis/CompareAlgs.py", "-c", cfg_path, "-o", outdir]
    with open(log, "a") as fh:
        fh.write("\n$ " + " ".join(cmd) + "\n")
        fh.flush()
        rc = subprocess.run(cmd, cwd=REPO, stdout=fh,
                            stderr=subprocess.STDOUT).returncode
    print(f"  {'OK ' if rc == 0 else 'FAIL'} {os.path.basename(cfg_path)} -> {outdir}")
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="nocuts")
    ap.add_argument("--outbase", default="TauPolOutputs/CLD_ILD")
    args = ap.parse_args()

    cfgdir = os.path.join(REPO, "config/plots/CLD_ILD")
    merged_dir = os.path.join(REPO, args.outbase, "merged_histos")
    logdir = os.path.join(REPO, "logs/CLD_ILD")
    for d in (cfgdir, merged_dir, logdir):
        os.makedirs(d, exist_ok=True)
    log = os.path.join(logdir, f"compare_plots_{args.variant}.log")
    open(log, "w").close()

    for ch in CHANNELS:
        print(f"\n===== channel {ch} ({args.variant}) =====")
        merged = {}
        for det, level in DIRS:
            # gen histograms have no cut variant: always take the nocuts ones
            v = "nocuts" if level == "gen" else args.variant
            p = merge_channel(det, level, ch, v, merged_dir)
            if p:
                merged[(det, level)] = p
        if not merged:
            continue

        outdir = os.path.join(REPO, args.outbase, f"{ch}_{args.variant}")
        lab = {(d, l): f"{d} {'reco' if l == 'reco' else 'gen'}" for d, l in merged}
        had = TITLES[ch]
        partner = PARTNER[CHANNELS[ch][1]]
        # En los ficheros single-decay el otro hemisferio siempre cae en BG, asi
        # que la senal del canal inclusivo es SIGNAL_BG, no SIGNAL_SIGNAL.
        sigcat = "SIGNAL_SIGNAL" if CHANNELS[ch][1] == "lep" else "SIGNAL_BG"

        def datasets(keys):
            return "datasets:\n" + "".join(
                ds_block(lab[k], merged[k], k[0], k[1]) for k in keys)

        # ── optimal variable: all four, then reco-only and gen-only ──────────
        groups = {
            "optvar_all":  list(merged),
            "optvar_reco": [k for k in merged if k[1] == "reco"],
            "optvar_gen":  [k for k in merged if k[1] == "gen"],
        }
        for gname, keys in groups.items():
            if not keys:
                continue
            txt = datasets(keys) + "\nplots:\n"
            txt += plot_block(
                f"{gname}_signal",
                {lab[k]: f"{OPTVAR[ch]}_{sigcat}" for k in keys},
                rf"\text{{Optimal variable, }} {had} + {partner} \text{{ (signal)}}",
                OBS_LABEL[ch])
            txt += plot_block(
                f"{gname}_allall",
                {lab[k]: f"{OPTVAR[ch]}_ALL_ALL" for k in keys},
                rf"\text{{Optimal variable, }} {had} + {partner} \text{{ (all selected)}}",
                OBS_LABEL[ch])
            write_and_run(txt, os.path.join(cfgdir, f"{ch}_{gname}_{args.variant}.yaml"),
                          outdir, log)

        # ── kinematics (reco only) ───────────────────────────────────────────
        reco_keys = [k for k in merged if k[1] == "reco"]
        if reco_keys:
            txt = datasets(reco_keys) + "\nplots:\n"
            txt += plot_block("visP_hadron",
                              {lab[k]: f"RecoVisP_dec0_{sigcat}" for k in reco_keys},
                              rf"\text{{Visible P of the }} {had}", '"P [GeV/c]"')
            txt += plot_block("visP_lepton",
                              {lab[k]: f"RecoVisP_dec1_{sigcat}" for k in reco_keys},
                              rf"\text{{Visible P of the other hemisphere ({partner})}}", '"P [GeV/c]"')
            txt += plot_block("reco_zmass",
                              {lab[k]: "RecoZMass_ALL_ALL" for k in reco_keys},
                              r"\text{Reconstructed } m_{Z}",
                              r'"m_{Z} [GeV/c^{2}]"')
            txt += plot_block("deltaR",
                              {lab[k]: "DeltaR_dec0_dec1_ALL_ALL" for k in reco_keys},
                              r"\Delta R \text{ between hemispheres}", r'"\Delta R"')
            txt += plot_block("cos_theta_vis",
                              {lab[k]: f"RecoVisCosTheta_dec0_{sigcat}" for k in reco_keys},
                              rf"\cos\theta \text{{ of the }} {had}", r'"\cos\theta_{vis}"')
            write_and_run(txt, os.path.join(cfgdir, f"{ch}_kinematics_{args.variant}.yaml"),
                          outdir, log)

        # ── fit templates: nominal vs correlated P1/M1, per detector ─────────
        for det, level in reco_keys:
            path = merged[(det, level)]
            txt = ("datasets:\n"
                   f"  - path: {path}\n    label: SM\n    color: ROOT.kBlack\n"
                   "    linestyle: 1\n    markerstyle: 20\n    markersize: 1.0\n    linewidth: 3\n"
                   f"  - path: {path}\n    label: A_tau=+1\n    color: ROOT.kRed\n"
                   "    linestyle: 2\n    markerstyle: 21\n    markersize: 1.0\n    linewidth: 3\n"
                   f"  - path: {path}\n    label: A_tau=-1\n    color: ROOT.kBlue\n"
                   "    linestyle: 2\n    markerstyle: 22\n    markersize: 1.0\n    linewidth: 3\n"
                   "\nplots:\n")
            # Same weights the fit uses: gen-level correlated morphing (`_corr_*`)
            # applied to the reco observable. The suffix names the ORIGIN
            # scenario, so `_corr_M1` is the template that mimics A_tau = +1.
            txt += plot_block(
                f"templates_{det}",
                {"SM": f"{OPTVAR[ch]}_{sigcat}",
                 "A_tau=+1": f"{OPTVAR[ch]}_{sigcat}_corr_M1",
                 "A_tau=-1": f"{OPTVAR[ch]}_{sigcat}_corr_P1"},
                rf"\text{{{det}: reweighted templates, }} {had} + {partner}",
                OBS_LABEL[ch], ylab="Norm. events")
            write_and_run(txt, os.path.join(cfgdir, f"{ch}_templates_{det}_{args.variant}.yaml"),
                          outdir, log)


if __name__ == "__main__":
    main()
