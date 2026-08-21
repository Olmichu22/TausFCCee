"""
genOnlyRHOTree_parallel.py

Versión gen-only de analysisRHOTree_parallel.py.
No reconstruye taus a nivel reco. Selecciona a nivel generador eventos con
exactamente 2 gen-taus cuyos decayIDs estén en la lista decay_modes.

Nomenclatura de ramas: tau1_* (mayor energía) y tau2_* (menor energía).
Las ramas "reco" se rellenan con los valores gen para compatibilidad directa
con rhoHistFromTree_parallel.py.

USO:
    python genOnlyRHOTree_parallel.py [mismos args que analysisRHOTree_parallel.py]
        [--decay-modes ID [ID ...]] [--n-workers N]

NOTA: --decay-modes acepta cualquier lista de decayIDs a nivel reco (2 se
      remapea a 1 gen). Sin --decay-modes se aceptan todos los pares.
"""

import ctypes
import logging
import math
import multiprocessing
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import ROOT
import yaml
from podio import root_io
from ROOT import TFile, TTree

from modules import myutils, tauReco
from modules import optimalVariabRho
from modules import weightsPol
from modules.rhoParallelUtils import (
    split_filenames, _make_p4_from_const, _setup_worker_logging, merge_partial_root_files,
)


# ── Config ────────────────────────────────────────────────────────────────────

_DEFAULT_CONFIG = "config/default/taurecolong.yaml"
_OUTPUT_BASE    = "Results/RhoAnalysis/"
# Subcarpeta dentro de <outputpath>/logs/ donde se guardan los logs de este script
_LOG_SOURCE     = "genOnlyRHOTree_MDecs"

# Sufijos de ramas escalares por tau (se prefijan con "tau1_" / "tau2_")
_TAU_SCALAR_SUFFIXES = [
    "P", "E", "M", "Theta", "Phi",
    "visP", "visE", "visM", "visTheta", "visPhi",
    "pionP", "pionE", "pionM", "pionTheta", "pionPhi",
    "pionPDG",
    "lepP", "lepE", "lepTheta", "lepPhi", "lepPDG",
    "decayID",
    "tauPDG",
    "genHelicity",
    "cos_theta", "cos_psi", "cos_beta",
    "omega",
    "weight_P1", "weight_M1",
    "cos_theta_tau",
    "optimalVar",
    "isElectron",
    "nPhotons",
    "recoVisP", "recoVisE", "recoVisM", "recoVisTheta", "recoVisPhi",
    "recoPionP", "recoPionE", "recoPionM", "recoPionTheta", "recoPionPhi",
    "recoTauID",
    "recoCharge",
    "recoLepP", "recoLepE", "recoLepTheta", "recoLepPhi", "recoLepPDG",
]

# Sufijos de ramas vectoriales por tau
_TAU_VECTOR_SUFFIXES = [
    "photons_E", "photons_theta", "photons_phi",
    # Espejo reco de los fotones gen (gen-only): igual que el resto de ramas reco
    # son copia del gen. Necesario para que el path NN de RhoHistFromTree
    # (--use-nn-optimal) pueda correr también sobre árboles gen (closure).
    "reco_photons_E", "reco_photons_theta", "reco_photons_phi",
]

# Ramas compartidas (sin prefijo tau)
_SHARED_VARIABS = ["GenZMass", "GenZVisMass", "ZMass", "beamE"]

# Listas completas construidas con prefijos
_VARIABS = (
    _SHARED_VARIABS
    + [f"tau1_{s}" for s in _TAU_SCALAR_SUFFIXES]
    + [f"tau2_{s}" for s in _TAU_SCALAR_SUFFIXES]
)
_VECTOR_VARIABS = (
    [f"tau1_{s}" for s in _TAU_VECTOR_SUFFIXES]
    + [f"tau2_{s}" for s in _TAU_VECTOR_SUFFIXES]
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _remap_to_gen(decay_modes):
    """Convierte IDs reco a IDs gen: 2 → 1 (rho reco → rho gen). Devuelve set o None."""
    if decay_modes is None:
        return None
    return {1 if d == 2 else d for d in decay_modes}


def _fill_tau_branches(branches, prefix, tauObj, beamE, sin_eff):
    """
    Rellena todas las ramas {prefix}_* para un gen-tau dado.
    Las ramas "reco" se copian de los valores gen (árbol gen-only).
    Fallbacks: 0 para variables no aplicables al decay, -999 para variables
    cuya fórmula no está definida para ese decay.
    """
    decayID = tauObj.getID()
    tauP4   = tauObj.getMomentum()    # TLorentzVector
    visP4   = tauObj.getvisMomentum() # TLorentzVector

    # ── Tau 4-momentum ────────────────────────────────────────────────────────
    branches[f"{prefix}_P"].value         = tauP4.P()
    branches[f"{prefix}_E"].value         = tauP4.E()
    branches[f"{prefix}_M"].value         = tauP4.M()
    branches[f"{prefix}_Theta"].value     = tauP4.Theta()
    branches[f"{prefix}_Phi"].value       = tauP4.Phi()
    branches[f"{prefix}_cos_theta_tau"].value = math.cos(tauP4.Theta())

    # ── Visible P4 ────────────────────────────────────────────────────────────
    branches[f"{prefix}_visP"].value      = visP4.P()
    branches[f"{prefix}_visE"].value      = visP4.E()
    branches[f"{prefix}_visM"].value      = visP4.M()
    branches[f"{prefix}_visTheta"].value  = visP4.Theta()
    branches[f"{prefix}_visPhi"].value    = visP4.Phi()
    # reco mirror
    branches[f"{prefix}_recoVisP"].value      = visP4.P()
    branches[f"{prefix}_recoVisE"].value      = visP4.E()
    branches[f"{prefix}_recoVisM"].value      = visP4.M()
    branches[f"{prefix}_recoVisTheta"].value  = visP4.Theta()
    branches[f"{prefix}_recoVisPhi"].value    = visP4.Phi()

    # ── Decay ID ──────────────────────────────────────────────────────────────
    branches[f"{prefix}_decayID"].value   = float(decayID)
    tau_pdg = int(tauObj.getPDG())
    branches[f"{prefix}_tauPDG"].value    = float(tau_pdg)
    branches[f"{prefix}_recoTauID"].value = float(decayID)
    # Espejo reco de la carga (árbol gen-only): tau⁻ (PDG 15) tiene carga −1.
    branches[f"{prefix}_recoCharge"].value = -1.0 if tau_pdg == 15 else 1.0
    hel = tauObj.getHelicity()
    branches[f"{prefix}_genHelicity"].value = float(hel) if hel is not None else -999.0

    # ── Daughters: charged pion y fotones de π⁰ ───────────────────────────────
    daughters = tauObj.getDaughters()
    pionP4 = ROOT.TLorentzVector()
    pionP4.SetXYZM(0.0, 0.0, 0.0, 0.0)
    pion_found = False
    pion_pdg = -999
    for key in daughters:
        # if abs(daughters[key].getPDG()) == 211:
        if abs(daughters[key].getPDG()) in (211, 321, 323):  # kaones y piones tratados como rho #SOLVED BUG
            pionP4 = _make_p4_from_const(daughters[key])
            pion_found = True
            pion_pdg = daughters[key].getPDG()
            break

    # if decayID == 1:
    #     daughter_pdgs = [int(daughters[key].getPDG()) for key in daughters]
    #     print(
    #         f"[rho-debug] {prefix} tauPDG={int(tauObj.getPDG())} "
    #         f"nDaughters={len(daughters)} foundPion={pion_found} "
    #         f"pionPDG={pion_pdg} pionM={pionP4.M():.6f} "
    #         f"daughtersPDG={daughter_pdgs}"
    #     )

    branches[f"{prefix}_pionP"].value      = pionP4.P()
    branches[f"{prefix}_pionE"].value      = pionP4.E()
    branches[f"{prefix}_pionM"].value      = pionP4.M()
    branches[f"{prefix}_pionTheta"].value  = pionP4.Theta()
    branches[f"{prefix}_pionPhi"].value    = pionP4.Phi()
    branches[f"{prefix}_pionPDG"].value    = float(pion_pdg)
    # reco mirror
    branches[f"{prefix}_recoPionP"].value      = pionP4.P()
    branches[f"{prefix}_recoPionE"].value      = pionP4.E()
    branches[f"{prefix}_recoPionM"].value      = pionP4.M()
    branches[f"{prefix}_recoPionTheta"].value  = pionP4.Theta()
    branches[f"{prefix}_recoPionPhi"].value    = pionP4.Phi()

    # Fotones gen desde π⁰
    nPhotons = 0.0
    for key in daughters:
        if abs(daughters[key].getPDG()) == 111:
            pi0_daughters = daughters[key].getDaughters()
            for photon_const in pi0_daughters:
                if photon_const.getGeneratorStatus() != 1:
                    continue
                nPhotons += 1
                ph = _make_p4_from_const(photon_const)
                branches[f"{prefix}_photons_E"].push_back(ph.E())
                branches[f"{prefix}_photons_theta"].push_back(ph.Theta())
                branches[f"{prefix}_photons_phi"].push_back(ph.Phi())
                # Espejo reco (gen-only): copia de los fotones gen.
                branches[f"{prefix}_reco_photons_E"].push_back(ph.E())
                branches[f"{prefix}_reco_photons_theta"].push_back(ph.Theta())
                branches[f"{prefix}_reco_photons_phi"].push_back(ph.Phi())
    branches[f"{prefix}_nPhotons"].value = nPhotons

    # ── Variables de polarización (dispatch por decayID) ──────────────────────
    if decayID == 1:  # ρ: fórmula completa con ángulos + ω
        (cos_theta, cos_psi, cos_beta, gen_w,
         weight_P1, weight_M1) = optimalVariabRho.wVariab(
            tauP4, visP4, pionP4, beamE, sin_eff=sin_eff, tau_pdg=tau_pdg)
        branches[f"{prefix}_cos_theta"].value  = cos_theta
        branches[f"{prefix}_cos_psi"].value    = cos_psi
        branches[f"{prefix}_cos_beta"].value   = cos_beta
        branches[f"{prefix}_omega"].value      = gen_w
        branches[f"{prefix}_weight_P1"].value  = weight_P1
        branches[f"{prefix}_weight_M1"].value  = weight_M1

    elif decayID == 0:  # π: peso gen con cos(θ*) geométrico (boost exacto, α=1)
        cts = weightsPol.cosThetaStar(tauP4, visP4)
        weight_P1 = weightsPol.newAtauFromH(tauP4, cts, +1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        weight_M1 = weightsPol.newAtauFromH(tauP4, cts, -1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        branches[f"{prefix}_cos_theta"].value  = cts
        branches[f"{prefix}_cos_psi"].value    = 0.0
        branches[f"{prefix}_cos_beta"].value   = 0.0
        branches[f"{prefix}_omega"].value      = -999.0
        branches[f"{prefix}_weight_P1"].value  = weight_P1
        branches[f"{prefix}_weight_M1"].value  = weight_M1

    elif decayID == 10:  # a1: sólo pesos (z_R vía newAtau), sin ángulos ρ
        weight_P1 = weightsPol.newAtau(tauP4, visP4, decayID, +1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        weight_M1 = weightsPol.newAtau(tauP4, visP4, decayID, -1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        branches[f"{prefix}_cos_theta"].value  = 0.0
        branches[f"{prefix}_cos_psi"].value    = 0.0
        branches[f"{prefix}_cos_beta"].value   = 0.0
        branches[f"{prefix}_omega"].value      = -999.0
        branches[f"{prefix}_weight_P1"].value  = weight_P1
        branches[f"{prefix}_weight_M1"].value  = weight_M1

    elif decayID in (-11, -13):  # leptónico
        weight_P1 = weightsPol.newAtauLep(visP4, tauP4, beamE, +1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        weight_M1 = weightsPol.newAtauLep(visP4, tauP4, beamE, -1, tau_pdg=tau_pdg, sin_eff=sin_eff)
        branches[f"{prefix}_cos_theta"].value  = 0.0
        branches[f"{prefix}_cos_psi"].value    = 0.0
        branches[f"{prefix}_cos_beta"].value   = 0.0
        branches[f"{prefix}_omega"].value      = -999.0
        branches[f"{prefix}_weight_P1"].value  = weight_P1
        branches[f"{prefix}_weight_M1"].value  = weight_M1

    else:  # decay no contemplado
        branches[f"{prefix}_cos_theta"].value  = 0.0
        branches[f"{prefix}_cos_psi"].value    = 0.0
        branches[f"{prefix}_cos_beta"].value   = 0.0
        branches[f"{prefix}_omega"].value      = -999.0
        branches[f"{prefix}_weight_P1"].value  = 1.0
        branches[f"{prefix}_weight_M1"].value  = 1.0

    # Variable óptima (observable): definición ÚNICA vía helper compartido
    # (π/a1/lep = E_vis/E_τ, ρ = ω). Evita drift con el hist stage.
    branches[f"{prefix}_optimalVar"].value = optimalVariabRho.optimal_var(
        decayID, tauP4, visP4, pionP4, beamE, tau_pdg=tau_pdg)

    # ── Leptón visible (solo para decays leptónicos) ──────────────────────────
    if decayID in (-11, -13):
        branches[f"{prefix}_lepP"].value      = visP4.P()
        branches[f"{prefix}_lepE"].value      = visP4.E()
        branches[f"{prefix}_lepTheta"].value  = visP4.Theta()
        branches[f"{prefix}_lepPhi"].value    = visP4.Phi()
        branches[f"{prefix}_lepPDG"].value    = float(abs(decayID))
        branches[f"{prefix}_isElectron"].value = float(decayID == -11)
        # reco mirrors
        branches[f"{prefix}_recoLepP"].value      = visP4.P()
        branches[f"{prefix}_recoLepE"].value      = visP4.E()
        branches[f"{prefix}_recoLepTheta"].value  = visP4.Theta()
        branches[f"{prefix}_recoLepPhi"].value    = visP4.Phi()
        branches[f"{prefix}_recoLepPDG"].value    = float(abs(decayID))
    else:
        branches[f"{prefix}_lepP"].value      = 0.0
        branches[f"{prefix}_lepE"].value      = 0.0
        branches[f"{prefix}_lepTheta"].value  = 0.0
        branches[f"{prefix}_lepPhi"].value    = 0.0
        branches[f"{prefix}_lepPDG"].value    = 0.0
        branches[f"{prefix}_isElectron"].value = 0.0
        branches[f"{prefix}_recoLepP"].value      = 0.0
        branches[f"{prefix}_recoLepE"].value      = 0.0
        branches[f"{prefix}_recoLepTheta"].value  = 0.0
        branches[f"{prefix}_recoLepPhi"].value    = 0.0
        branches[f"{prefix}_recoLepPDG"].value    = 0.0


# ── Worker ────────────────────────────────────────────────────────────────────

def process_chunk_gen_only(filenames_chunk, config_bundle, worker_id):
    """
    Procesa un chunk de ficheros ROOT a nivel gen.
    Selecciona pares de taus donde ambos decayIDs están en gen_filter.
    Las ramas reco se rellenan con los valores gen.
    """
    outputpath       = config_bundle["outputpath"]
    gen_filter       = config_bundle["gen_filter"]   # set o None (acepta todo)
    sin_eff          = config_bundle["sin_eff"]
    helicity        = config_bundle.get("helicity", False)
    fileOutName_base = config_bundle["fileOutName_base"]

    loggers        = _setup_worker_logging(outputpath, worker_id, _LOG_SOURCE)
    logger_process = loggers["processing"]
    logger_io      = loggers["io"]

    
    genparts = "MCParticles"

    # ── Crear TTree de salida parcial ─────────────────────────────────────────
    partial_path = os.path.join(outputpath,
                                f"partial_worker{worker_id}_{fileOutName_base}.root")
    partial_outfile = TFile(partial_path, "RECREATE")

    tree = TTree("outtree_original", "gen-only processed variables")
    branches = {}
    for var in _VARIABS:
        branches[var] = ctypes.c_double(0.0)
        tree.Branch(var, ctypes.addressof(branches[var]), f"{var}/D")
    for var in _VECTOR_VARIABS:
        branches[var] = ROOT.std.vector("double")()
        tree.Branch(var, branches[var])

    totalEvents    = 0
    selectedEvents = 0
    skipped_files  = []

    # ── Event loop ────────────────────────────────────────────────────────────
    for filename in filenames_chunk:
        try:
            file_reader = root_io.Reader([filename])
        except Exception as exc:
            logger_io.warning("Worker %d: skipping unreadable file %s: %s",
                              worker_id, filename, exc)
            skipped_files.append(filename)
            continue

        for event in file_reader.get("events"):
            if totalEvents % 500 == 0:
                logger_process.info("Worker %d: processed %d events",
                                    worker_id, totalEvents)
            totalEvents += 1

            # Limpiar vectores al inicio de cada evento
            for var in _VECTOR_VARIABS:
                branches[var].clear()

            mc_particles = event.get(genparts)
            beamE = mc_particles[0].getEnergy()

            genTaus = tauReco.findAllGenTaus(mc_particles, getHelicity=helicity)
            n_taus = len(genTaus)

            if n_taus < 2:
                continue
            if n_taus > 2:
                logger_process.warning(
                    "Worker %d: event %d has %d gen taus (expected 2) — skipping",
                    worker_id, totalEvents, n_taus)
                continue

            tau_list = list(genTaus.values())  # exactamente 2
            ids = [t.getID() for t in tau_list]

            # Ambos taus deben tener su decayID en gen_filter
            if gen_filter is not None and not all(i in gen_filter for i in ids):
                continue

            # Ordenar por energía descendente: tau1 = mayor energía
            tau_list.sort(key=lambda t: t.getMomentum().E(), reverse=True)
            tau1_obj, tau2_obj = tau_list

            # ── Cinemática compartida ─────────────────────────────────────────
            tau1P4  = tau1_obj.getMomentum()
            tau2P4  = tau2_obj.getMomentum()
            tau1vis = tau1_obj.getvisMomentum()
            tau2vis = tau2_obj.getvisMomentum()

            ZGenP4      = tau1P4 + tau2P4
            ZVisP4      = tau1vis + tau2vis
            GenZMass    = ZGenP4.M()
            GenZVisMass = ZVisP4.M()

            branches["GenZMass"].value    = GenZMass
            branches["GenZVisMass"].value = GenZVisMass
            branches["ZMass"].value       = GenZVisMass
            branches["beamE"].value       = beamE

            # ── Ramas por tau ─────────────────────────────────────────────────
            _fill_tau_branches(branches, "tau1", tau1_obj, beamE, sin_eff)
            _fill_tau_branches(branches, "tau2", tau2_obj, beamE, sin_eff)

            selectedEvents += 1
            tree.Fill()

    if skipped_files:
        logger_io.warning("Worker %d: %d file(s) skipped due to read errors",
                          worker_id, len(skipped_files))

    partial_outfile.cd()
    tree.Write()
    partial_outfile.Close()

    logger_io.info("Worker %d: done. Events=%d Selected=%d",
                   worker_id, totalEvents, selectedEvents)
    return partial_path, totalEvents, selectedEvents


# ── Merge ─────────────────────────────────────────────────────────────────────

# ── Main ──────────────────────────────────────────────────────────────────────

def my_hook(parser):
    parser.add_argument("--sin-eff", type=float, default=None,
                        help="Effective sin^2 theta_W for weight calculation")
    parser.add_argument("--n-workers", type=int, default=None,
                        help="Número de procesos paralelos (default: min(n_files, n_cpus))")
    parser.add_argument("--decay-modes", type=int, nargs="+", default=None,
                        metavar="ID",
                        help="Lista de decayIDs permitidos para AMBOS taus del par "
                             "(ej: 0 2 -11 -13). 2 se remapea a 1 (rho gen). "
                             "Default: acepta todos los pares.")
    parser.add_argument("--helicity", action="store_true",
                        help="Si se activa, se rellenan ramas de helicidad y (si la muestra las tiene)")


def main():
    general_configs = myutils.setup_analysis_config(_DEFAULT_CONFIG, _OUTPUT_BASE,
                                                    parser_hook=my_hook,
                                                    log_subdir=_LOG_SOURCE)
    loggers    = general_configs["loggers"]
    run_config = general_configs["config"]
    args       = general_configs["args"]

    outputpath  = general_configs["outputpath"]
    fileOutName = os.path.join(outputpath, general_configs["fileOutName"])
    sample      = run_config["general"]["sample"]
    test_arg    = general_configs["flags"]["test"]

    fileOutName_base = Path(fileOutName).stem

    # decay_modes: CLI tiene prioridad sobre config YAML
    decay_modes = args.decay_modes or run_config.get("general", {}).get("decay_modes")
    gen_filter  = _remap_to_gen(decay_modes)
    if gen_filter is not None:
        loggers["config"].info("decay_modes (reco): %s  →  gen_filter: %s",
                               decay_modes, gen_filter)
    else:
        loggers["config"].info("No decay_modes specified — all pairs accepted")

    filenames, _ = myutils.get_root_trees_path(
        sample, None, loggers, test_arg, args,
        skip_root_validation=True,
    )
    loggers["io"].info("Found %d input files", len(filenames))
    if not filenames:
        loggers["io"].error("No input files found. Aborting.")
        sys.exit(1)

    if args.input_list and len(args.input_list) == 1:
        n_workers = 1
    else:
        n_workers = args.n_workers or min(len(filenames), os.cpu_count() or 1)
    loggers["io"].info("Using %d workers for %d files", n_workers, len(filenames))

    config_bundle = {
        "selectDecay":      -777,       # siempre -777; código downstream lo detecta
        "gen_filter":       gen_filter,
        "sin_eff":          args.sin_eff,
        "outputpath":       outputpath,
        "fileOutName_base": fileOutName_base,
        "helicity":        args.helicity,
    }

    if n_workers == 1:
        loggers["io"].info("Sequential mode (n_workers=1).")
        partial_path, tot, sel = process_chunk_gen_only(filenames, config_bundle, 0)
        import shutil
        shutil.move(partial_path, fileOutName)
        loggers["io"].info("Sequential done: events=%d selected=%d", tot, sel)
        totals = {"events": tot, "selected": sel}
    else:
        file_chunks = split_filenames(filenames, n_workers)
        n_chunks    = len(file_chunks)

        ctx = multiprocessing.get_context("fork")
        partial_files = []
        totals = {"events": 0, "selected": 0}
        t_start = time.time()

        loggers["io"].info("Launching %d workers...", n_chunks)
        with ProcessPoolExecutor(max_workers=n_chunks, mp_context=ctx) as executor:
            futures = {
                executor.submit(process_chunk_gen_only,
                                file_chunks[i], config_bundle, i): i
                for i in range(n_chunks)
            }
            for n_done, future in enumerate(as_completed(futures), start=1):
                wid = futures[future]
                try:
                    partial_path, tot, sel = future.result()
                    partial_files.append(partial_path)
                    totals["events"]   += tot
                    totals["selected"] += sel
                    elapsed = time.time() - t_start
                    print(f"  [{n_done}/{n_chunks}] worker {wid} terminado "
                          f"({tot} eventos, {elapsed:.1f}s acumulado)", flush=True)
                except Exception as exc:
                    loggers["io"].error("Worker %d falló: %s", wid, exc)
                    raise

        loggers["io"].info("Merging %d partial files → %s",
                           len(partial_files), fileOutName)
        merge_partial_root_files(partial_files, fileOutName)

        for p in partial_files:
            try:
                os.remove(p)
            except OSError:
                pass

        loggers["io"].info("Total events=%d selected=%d (%.1fs)",
                           totals["events"], totals["selected"],
                           time.time() - t_start)

    results_df = pd.DataFrame({
        "TotalEvents":    [totals["events"]],
        "SelectedEvents": [totals["selected"]],
    })
    decay_tag = ("_".join(str(d) for d in sorted(decay_modes))
                 if decay_modes else "all")
    results_df.to_csv(
        os.path.join(outputpath, f"results_summary_{decay_tag}.csv"), index=False)

    output_config_file = os.path.join(outputpath, "config.yaml")
    run_config["args"]         = vars(args)
    run_config["run_Type"]     = "genOnlyRHOTree"
    run_config["general"]["decay_modes"] = decay_modes
    with open(output_config_file, "w") as f:
        yaml.dump(run_config, f)
    loggers["io"].info("Config guardado en %s", output_config_file)
    loggers["io"].info("Fichero de salida: %s", fileOutName)


if __name__ == "__main__":
    main()
