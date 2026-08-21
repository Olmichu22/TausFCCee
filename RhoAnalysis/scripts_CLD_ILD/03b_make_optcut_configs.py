#!/usr/bin/env python3
"""03b_make_optcut_configs.py — turn the PSO results into hist-stage configs.

Reads `optimization_results.csv` written by `optimize_cuts.py` for every
detector and channel and writes `config/pipeline/CLD_ILD/<DET>_reco_optcuts.yaml`,
a copy of the no-cuts pipeline where each channel run carries its own optimized
`ang` (dR), `meson_cut`, `lepton_cut` and `zmass_cut`.

Runs written, mirroring the no-cuts configs: the inclusive single-decay channels
of pion, rho and electron, plus the rho x e / rho x mu pairs. The inclusive runs
drop `lepton_cut`, which constrains dec1: in single-decay mode the other
hemisphere is any species, not necessarily the partner the PSO optimized
against. For the electron channel the PSO's `mesonP` bounds are the cut on the
electron momentum itself (it occupies the dec0 slot of the adapter), so
`meson_cut` carries them. The inclusive muon channel is not written: the PSO
does not cover it yet, so its optimized production would repeat the no-cuts one.
"""
import csv
import os

REPO = "/nfs/cms/arqolmo/TausFCCee"
OPT = "Results/CutOptimization_CLD_ILD"
CFG = "config/pipeline/CLD_ILD"
# Canales con cortes optimizados por la PSO (a1 y muon fuera por ahora).
# Ojo: la clave es el nombre del run en el config (<clave>_incl), y para el
# electron el nombre del canal aguas abajo (04/05/06) es "ele".
CHANNELS = {"pion": 0, "rho": 2, "ele": -11}
LEPTONS = {"el": -11, "mu": -13}

HEADER = """# ═══════════════════════════════════════════════════════════════════════════
# Estudio CLD vs ILD — etapa 2 con los CORTES OPTIMIZADOS por PSO
# (RhoAnalysis/optimize_cuts.py, ver Results/CutOptimization_CLD_ILD/).
# Generado por RhoAnalysis/scripts_CLD_ILD/03b_make_optcut_configs.py: no editar
# a mano, se regenera al relanzar la optimización.
# ═══════════════════════════════════════════════════════════════════════════
"""


def read_cuts(det, ch):
    path = os.path.join(REPO, OPT, f"{det}_{ch}", "optimization_results.csv")
    if not os.path.isfile(path):
        print(f"  [WARN] missing {path}")
        return None
    with open(path) as f:
        return next(csv.DictReader(f))


def main():
    for det, prefix in (("CLD", "PolAnalysis_RECO_CLD_ztt2M_"),
                        ("ILD", "PolAnalysis_RECO_ILD_fcc_")):
        outdir = (f"Results/RhoAnalysis/{prefix}"
                  "tau_trained0.4_tph0.35_tpi0_n3_g0.0")
        tree = (f"{outdir}/TTree_MDecs_0_1_2_10_-11_-13_"
                "tau_traineddecayAll_0.4_tph0.35_tpi0_n3_g0.0.root")
        txt = HEADER + f"""
pipeline:

  tree:
    enabled: false
    mode: reco
    prefix: {prefix}
    config: config/default/taurecolong_optimal.yaml

  tree_file: {tree}

  hist_common:
    hist_config_mdecs: config/histograms/rho_analysis_config_mdecs.yml
    config: config/default/taurecolong_optimal.yaml
    n_workers: 16
    only_gen: false
    compute_weights: false
    no_omega_weights: false
    sin_eff: 0.2312
    cos_acceptance: null
    omega_border_cut: false
    cut: []
    verbose: 1

  hist_runs:
"""
        n_ok = 0
        for ch, cid in CHANNELS.items():
            c = read_cuts(det, ch)
            if c is None:
                continue
            n_ok += 1
            ang = f"[{float(c['dR_min']):.4f}, {float(c['dR_max']):.4f}]"
            meson = f"[{float(c['mesonP_min']):.4f}, {float(c['mesonP_max']):.4f}]"
            lepton = f"[{float(c['lepP_min']):.4f}, {float(c['lepP_max']):.4f}]"
            zmass = f"[{float(c['Zmass_min']):.4f}, {float(c['Zmass_max']):.4f}]"

            # Canal inclusivo: sin lepton_cut (dec1 no es necesariamente un leptón)
            txt += f"""    - name: {ch}_incl
      single_decay: {cid}
      vism_cut: 0.
      ang: {ang}
      meson_cut: {meson}
      zmass_cut: {zmass}

"""
            # Pares canal x leptón: solo el rho los usa en la etapa 4 (rho_lep)
            if ch != "rho":
                continue
            for lname, lid in LEPTONS.items():
                txt += f"""    - name: {ch}_{lname}
      decay_pair: [{cid}, {lid}]
      vism_cut: 0.
      ang: {ang}
      meson_cut: {meson}
      lepton_cut: {lepton}
      zmass_cut: {zmass}

"""
        path = os.path.join(REPO, CFG, f"{det}_reco_optcuts.yaml")
        if n_ok == 0:
            print(f"[SKIP] {det}: no PSO results, config not written")
            continue
        open(path, "w").write(txt)
        print(f"[OK] {path} ({n_ok} channels)")


if __name__ == "__main__":
    main()
