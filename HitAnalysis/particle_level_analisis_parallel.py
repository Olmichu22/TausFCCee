import argparse
import logging
import math
import multiprocessing
import os
import pickle
import pprint
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import ROOT
import yaml
from podio import root_io
from ROOT import TH1F, TH2F, TFile, TTree
import edm4hep

from modules.TauDecays import extractTauDecays
from modules.NeutralRecover import debug_reco_tau, plot_debug_reco_tau, get_reco_mc_links_by_dR
from modules.ConfusionMatrixParticleLevel import plot_confusion_matrices, plot_energy_distributions, plot_efficiency_vs_momentum
from modules import (ParticleObjects, electronReco, muonReco, myutils, pi0Reco,
                     tauReco, particleMatch)
from modules.ParticleObjects import RecoParticle


# ── Helpers (module-level, usados tanto en main como en workers) ──────────────

def write_histograms_recursive(obj):
    """Recorre un dict anidado y llama .Write() en cada histograma ROOT."""
    if isinstance(obj, dict):
        for value in obj.values():
            write_histograms_recursive(value)
    else:
        try:
            obj.Write()
        except AttributeError:
            print(f"Objeto {obj} no tiene método .Write(). Ignorado.")


def split_filenames(filenames, n_workers):
    """Divide la lista de ficheros en n_workers chunks lo más iguales posible."""
    k, rem = divmod(len(filenames), n_workers)
    chunks, start = [], 0
    for i in range(n_workers):
        end = start + k + (1 if i < rem else 0)
        if start < end:
            chunks.append(filenames[start:end])
        start = end
    return chunks


def split_mlpf(mlpf_results, file_chunks):
    """
    Divide mlpf_results en sub-dicts con claves renormalizadas a índice local.

    En myutils.get_root_trees_path los eventos se indexan como:
        key_id = n_files * 1000 + local_key - 1
    Por tanto el chunk que empieza en el fichero global f_offset contiene claves
    en el rango [f_offset*1000, (f_offset + len(chunk))*1000).
    Las renormalizamos a 0..len(chunk)*1000 para que el eventid local del worker
    encaje directamente con mlpf_chunk.get(local_eventid).
    """
    mlpf_chunks = []
    file_offset = 0
    for chunk in file_chunks:
        lo = file_offset * 1000
        hi = (file_offset + len(chunk)) * 1000
        sub = {k - lo: v for k, v in mlpf_results.items() if lo <= k < hi}
        mlpf_chunks.append(sub)
        file_offset += len(chunk)
    return mlpf_chunks


# ── Función worker (debe ser picklable → nivel de módulo) ─────────────────────

def process_chunk(filenames_chunk, mlpf_chunk, global_event_offset,
                  config_bundle, worker_id):
    """
    Lee un subconjunto de ficheros ROOT y devuelve un DataFrame con las
    asociaciones gen-reco por evento.

    Parámetros
    ----------
    filenames_chunk : list[str]
        Paths de los ficheros ROOT que procesa este worker.
    mlpf_chunk : dict
        Subconjunto de mlpf_results con claves renormalizadas a índice local.
    global_event_offset : int
        Desplazamiento para calcular event_id globalmente único.
    config_bundle : dict
        Parámetros de configuración serializables (sin loggers ni objetos ROOT).
    worker_id : int
        Índice del worker (para logging y trazabilidad).

    Devuelve
    --------
    tuple[pd.DataFrame, pd.DataFrame]
        (df_links_dr, df_links_truth)
        DataFrames de asociaciones para matching por dR y por RecoMCTruthLink.
    """
    # ── Logging propio del worker ─────────────────────────────────────────────
    # Con fork el worker hereda los handlers del padre; los eliminamos para
    # evitar escrituras concurrentes en el mismo fichero de log.
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()

    log_file = os.path.join(config_bundle["outputpath"], f"worker_{worker_id}.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        force=True,
    )
    logger = logging.getLogger(f"worker_{worker_id}")
    logger.info("Worker %d arrancado, procesando %d ficheros", worker_id, len(filenames_chunk))

    # ── Extraer parámetros del bundle ─────────────────────────────────────────
    neutral_recover_cfg = config_bundle["neutral_recover_cfg"]
    cuts                = config_bundle["cuts"]
    photon_config       = config_bundle["photon_config"]
    test_extremes       = config_bundle["test_extremes"]
    test_pfo            = config_bundle["test_pfo"]
    gatr_results_path   = config_bundle["gatr_results_path"]
    pfobjects           = config_bundle["pfobjects"]
    dedup_mode          = config_bundle.get("dedup_mode", "reco")
    filter_gen_status   = config_bundle.get("filter_gen_status", True)
    max_gen_pdg         = config_bundle.get("max_gen_pdg", None)
    weight_mode         = config_bundle.get("weight_mode", "decoded")

    df_list_dr = []
    df_list_truth = []

    # PDGs de neutrinos (no dejan señal en detector)
    _NEUTRINO_PDGS = {12, 14, 16}
    # Solo partículas estables finales
    _VALID_GEN_STATUS = {1}

    expected_cols = [
        "gen", "reco", "Gen_pid", "Reco_pid", "Gen_energy", "Reco_energy", "event_id",
        "Gen_Px", "Gen_Py", "Gen_Pz", "Reco_Px", "Reco_Py", "Reco_Pz",
        "source_file", "event_in_file",
    ]

    # Contadores de diagnóstico para el resumen final del worker
    _truth_diag = {
        "n_events_total": 0,
        "n_events_with_links": 0,
        "total_links_raw": 0,
        "total_links_parsed": 0,
        "first_events_with_links": [],   # primeros event_ids globales con links
        "collection_used": None,
        "collection_found_at_event": None,
    }

    def _normalize_links_df(df_links, event_id_global, source_file="", event_in_file=-1):
        if df_links is None or df_links.empty:
            return pd.DataFrame(columns=expected_cols)

        df_out = df_links.copy()
        for col in expected_cols:
            if col not in df_out.columns:
                df_out[col] = np.nan

        df_out = df_out[expected_cols].copy()
        df_out["event_id"] = event_id_global
        df_out["source_file"] = source_file
        df_out["event_in_file"] = event_in_file
        return df_out

    def _safe_get_value(obj, method_name, default=np.nan):
        try:
            method = getattr(obj, method_name, None)
            if method is None:
                return default
            return method()
        except Exception:
            return default

    def _safe_get_index(obj, default=-999):
        try:
            obj_id = obj.getObjectID()
            return obj_id.index
        except Exception:
            return default

    def _get_reco_mc_links_truth(event, local_eventid, event_id_global):
        """
        Construye df_reco_mc_links usando RecoMCTruthLink (asociación por peso
        de depósito de energía, igual que analyze_particles.py).

        Filtra:
          - Solo MCParticles con generatorStatus == 1 (estables finales)
          - Excluye neutrinos (PDG 12, 14, 16)

        Nota: requiere podio compatible con el formato del fichero. Los ficheros
        ILD escritos con podio > 1.1.x usan podio::LinkData y no son legibles
        con la versión instalada (1.1.0); para ellos el reader devuelve 0 eventos.
        """
        _truth_diag["n_events_total"] += 1

        try:
            mc_particles = event.get("MCParticles")
        except Exception:
            mc_particles = []

        try:
            pfos = event.get("PandoraPFOs")
        except Exception:
            pfos = []

        reco_mc_truth_links = []
        used_collection = None

        # Intentar leer la colección (probamos el nombre estándar y cualquier
        # variante que aparezca en las colecciones disponibles)
        if _truth_diag["collection_used"] is not None:
            # Ya sabemos qué nombre funciona: usarlo directamente
            candidate_collections = [_truth_diag["collection_used"]]
        else:
            candidate_collections = ["RecoMCTruthLink"]
            try:
                available = list(event.getAvailableCollections())
                candidate_collections += [
                    n for n in available if "RecoMCTruthLink" in n and n != "RecoMCTruthLink"
                ]
            except Exception:
                pass
            candidate_collections = list(dict.fromkeys(candidate_collections))

        for coll_name in candidate_collections:
            try:
                reco_mc_truth_links = list(event.get(coll_name))
                used_collection = coll_name
                # Memorizar para eventos futuros
                if _truth_diag["collection_used"] is None:
                    _truth_diag["collection_used"] = coll_name
                    _truth_diag["collection_found_at_event"] = local_eventid
                break
            except Exception:
                continue

        n_links_raw = len(reco_mc_truth_links)
        _truth_diag["total_links_raw"] += n_links_raw

        reco_gen_link = {
            "gen": [], "reco": [], "weight": [],
            "Gen_pid": [], "Reco_pid": [],
            "Gen_energy": [], "Reco_energy": [],
            "Gen_Px": [], "Gen_Py": [], "Gen_Pz": [],
            "Reco_Px": [], "Reco_Py": [], "Reco_Pz": []
        }

        parsed_links = 0
        for link in reco_mc_truth_links:
            try:
                # check if link has getSim and getRec methods
                if not (hasattr(link, "getSim") and hasattr(link, "getRec")):
                    sim_obj = link.getTo()
                    rec_obj = link.getFrom()
                else:
                    sim_obj = link.getSim()
                    rec_obj = link.getRec()
                
                gen_status = _safe_get_value(sim_obj, "getGeneratorStatus", default=-1)
                gen_pdg    = _safe_get_value(sim_obj, "getPDG", default=0)
                if filter_gen_status and gen_status not in _VALID_GEN_STATUS:
                    continue
                abs_pdg = abs(int(gen_pdg))
                if abs_pdg in _NEUTRINO_PDGS:
                    continue
                if max_gen_pdg is not None and abs_pdg > max_gen_pdg:
                    continue
                reco_gen_link["gen"].append(_safe_get_index(sim_obj, default=-999))
                reco_gen_link["reco"].append(_safe_get_index(rec_obj, default=-999))
                reco_gen_link["weight"].append(_safe_get_value(link, "getWeight", default=0.0))
                reco_gen_link["Gen_pid"].append(gen_pdg)
                reco_gen_link["Reco_pid"].append(_safe_get_value(rec_obj, "getPDG", default=np.nan))
                reco_gen_link["Gen_energy"].append(_safe_get_value(sim_obj, "getEnergy", default=np.nan))
                reco_gen_link["Gen_Px"].append(_safe_get_value(sim_obj, "getMomentum", default=edm4hep.Vector3f()).x)
                reco_gen_link["Gen_Py"].append(_safe_get_value(sim_obj, "getMomentum", default=edm4hep.Vector3f()).y)
                reco_gen_link["Gen_Pz"].append(_safe_get_value(sim_obj, "getMomentum", default=edm4hep.Vector3f()).z)
                reco_gen_link["Reco_energy"].append(_safe_get_value(rec_obj, "getEnergy", default=np.nan))
                reco_gen_link["Reco_Px"].append(_safe_get_value(rec_obj, "getMomentum", default=edm4hep.Vector3f()).x)
                reco_gen_link["Reco_Py"].append(_safe_get_value(rec_obj, "getMomentum", default=edm4hep.Vector3f()).y)
                reco_gen_link["Reco_Pz"].append(_safe_get_value(rec_obj, "getMomentum", default=edm4hep.Vector3f()).z)
                parsed_links += 1
            except Exception:
                continue

        _truth_diag["total_links_parsed"] += parsed_links
        if parsed_links > 0:
            _truth_diag["n_events_with_links"] += 1
            if len(_truth_diag["first_events_with_links"]) < 5:
                _truth_diag["first_events_with_links"].append(event_id_global)

        df_reco_mc_links = pd.DataFrame(reco_gen_link)
        if df_reco_mc_links.empty:
            return pd.DataFrame(columns=list(reco_gen_link.keys()))

        # Columna de peso efectivo para ranking (--weight-mode):
        #   'raw'     → usar el valor codificado directamente
        #   'decoded' → decodificar track/cluster weight (Bohdan Dudar encoding):
        #               encodedWeight = int(clusterW*1000)*10000 + int(trackW*1000)
        #               preferir trackWeight si > 0, si no usar clusterWeight
        if weight_mode == "decoded":
            encoded = df_reco_mc_links["weight"].astype(int)
            track_w   = (encoded % 10000) / 1000.0
            cluster_w = (encoded // 10000) / 1000.0
            df_reco_mc_links["_eff_weight"] = np.where(track_w > 0.0, track_w, cluster_w)
            weight_col = "_eff_weight"
        else:
            weight_col = "weight"

        # Deduplicación configurable (--dedup-mode):
        #   'gen'  → un reco por gen, max weight  (legacy; permite reco duplicados)
        #   'reco' → un gen por reco, max weight  (resuelve fusión de fotones)
        if dedup_mode == "reco":
            df_reco_mc_links = df_reco_mc_links.sort_values("reco").reset_index(drop=True)
            df_reco_mc_links = df_reco_mc_links.loc[
                df_reco_mc_links.groupby("reco")[weight_col].idxmax()
            ].reset_index(drop=True)
        else:
            df_reco_mc_links = df_reco_mc_links.sort_values("gen").reset_index(drop=True)
            df_reco_mc_links = df_reco_mc_links.loc[
                df_reco_mc_links.groupby("gen")[weight_col].idxmax()
            ].reset_index(drop=True)
        df_reco_mc_links.drop(
            columns=[c for c in ["weight", "_eff_weight"] if c in df_reco_mc_links.columns],
            inplace=True,
        )

        gen_indices  = set(df_reco_mc_links["gen"])
        reco_indices = set(df_reco_mc_links["reco"])
        new_rows = []

        # Gen sin match en reco
        for part in mc_particles:
            try:
                gen_status = part.getGeneratorStatus()
                if filter_gen_status and gen_status not in _VALID_GEN_STATUS:
                    continue
                gen_pdg = part.getPDG()
                abs_pdg = abs(int(gen_pdg))
                if abs_pdg in _NEUTRINO_PDGS:
                    continue
                if max_gen_pdg is not None and abs_pdg > max_gen_pdg:
                    continue
                idx = part.getObjectID().index
                if idx not in gen_indices:
                    new_rows.append(
                        {
                            "gen": idx,
                            "reco": -999,
                            "Gen_pid": gen_pdg,
                            "Reco_pid": -999,
                            "Gen_energy": part.getEnergy(),
                            "Reco_energy": np.nan,
                            "Gen_Px": part.getMomentum().x,
                            "Gen_Py": part.getMomentum().y,
                            "Gen_Pz": part.getMomentum().z,
                            "Reco_Px": np.nan,
                            "Reco_Py": np.nan,
                            "Reco_Pz": np.nan,
                        }
                    )
            except Exception:
                continue

        for pfo in pfos:
            try:
                reco_idx = pfo.getObjectID().index
                if reco_idx not in reco_indices:
                    new_rows.append(
                        {
                            "gen": -999,
                            "reco": reco_idx,
                            "Gen_pid": -999,
                            "Reco_pid": pfo.getPDG(),
                            "Gen_energy": np.nan,
                            "Reco_energy": pfo.getEnergy(),
                            "Gen_Px": np.nan,
                            "Gen_Py": np.nan,
                            "Gen_Pz": np.nan,
                            "Reco_Px": pfo.getMomentum().x,
                            "Reco_Py": pfo.getMomentum().y,
                            "Reco_Pz": pfo.getMomentum().z,

                        }
                    )
            except Exception:
                continue

        if new_rows:
            df_reco_mc_links = pd.concat([df_reco_mc_links, pd.DataFrame(new_rows)], ignore_index=True)

        return df_reco_mc_links

    # ── Reader propio (se crea dentro del worker, después del fork) ───────────
    # Iteramos fichero a fichero para poder registrar source_file y event_in_file
    cumulative_local_eventid = 0
    for filename in filenames_chunk:
        source_file = os.path.basename(filename)
        file_reader = root_io.Reader([filename])
        file_local_eventid = -1
        for file_local_eventid, event in enumerate(file_reader.get("events")):
            local_eventid = cumulative_local_eventid + file_local_eventid
            if local_eventid % 1000 == 0:
                logger.info("Worker %d: evento local %d", worker_id, local_eventid)

            pfos = event.get(pfobjects)

            try:
                # Ejecutar matching MLPF por dR siempre que haya resultados MLPF disponibles,
                # independientemente del flag return_hit_type_map (que controla visualización de hits).
                run_mlpf_match = (gatr_results_path is not None) and (mlpf_chunk.get(local_eventid) is not None)
                if run_mlpf_match:
                    # return_hit_type_map=True garantiza que extractTauDecays devuelva 6 valores
                    # (incluyendo extra_info_dict con df_reco_mc_links); sin él solo devuelve 4.
                    _neutral_cfg_assoc = {**neutral_recover_cfg, "return_hit_type_map": True}
                    (_, _, _, _, _,
                     extra_info_dict) = extractTauDecays(
                        gatr_results_path, mlpf_chunk, local_eventid,
                        pfos,
                        cuts["dRMax"], cuts["minPTauPhoton"], cuts["minPTauPion"],
                        cuts["PNeutron"], cuts["generalPCut"],
                        photon_config, test_extremes, test_pfo, logger,
                        neutral_recover_cfg=_neutral_cfg_assoc,
                        event=event,
                        only_association=True,
                    )
                    df_reco_mc_links_dr = extra_info_dict.get("df_reco_mc_links", pd.DataFrame())
                else:
                    # Sin MLPF: matching por dR entre MCParticles y PandoraPFOs.
                    # Se pasan mapas de hits vacíos → no hay filtrado por señal en detector.
                    df_reco_mc_links_dr = get_reco_mc_links_by_dR(
                        event, {}, {}, logger_process=logger
                    )

                df_reco_mc_links_truth = _get_reco_mc_links_truth(
                    event, local_eventid, global_event_offset + local_eventid
                )

            except Exception as exc:
                logger.error("Worker %d: error en evento local %d: %s",
                             worker_id, local_eventid, exc)
                continue

            event_id_global = global_event_offset + local_eventid
            df_out_dr = _normalize_links_df(df_reco_mc_links_dr, event_id_global,
                                            source_file=source_file,
                                            event_in_file=file_local_eventid)
            if not df_out_dr.empty:
                df_list_dr.append(df_out_dr)

            df_out_truth = _normalize_links_df(df_reco_mc_links_truth, event_id_global,
                                               source_file=source_file,
                                               event_in_file=file_local_eventid)
            if not df_out_truth.empty:
                df_list_truth.append(df_out_truth)

        cumulative_local_eventid += file_local_eventid + 1

    # ── Resumen de diagnóstico de truth-links para este worker ───────────────
    logger.info(
        "Worker %d finalizado: %d eventos con datos dR, %d con datos truth",
        worker_id,
        len(df_list_dr),
        len(df_list_truth),
    )
    logger.info(
        "Worker %d truth-link stats: colección='%s' (encontrada en evento local %s) | "
        "%d/%d eventos con links | %d links brutos | %d links tras filtros | "
        "primeros event_ids con links: %s",
        worker_id,
        _truth_diag["collection_used"],
        _truth_diag["collection_found_at_event"],
        _truth_diag["n_events_with_links"],
        _truth_diag["n_events_total"],
        _truth_diag["total_links_raw"],
        _truth_diag["total_links_parsed"],
        _truth_diag["first_events_with_links"],
    )
    if _truth_diag["n_events_with_links"] == 0 and _truth_diag["n_events_total"] > 0:
        logger.warning(
            "Worker %d: ningún evento tuvo RecoMCTruthLink con datos. "
            "Posibles causas: (1) ficheros ILD incompatibles con podio 1.1.0 "
            "(usar podio::LinkData requiere versión más reciente); "
            "(2) la colección existe pero está vacía en todos los eventos.",
            worker_id,
        )

    df_dr = pd.concat(df_list_dr, ignore_index=True) if df_list_dr else pd.DataFrame(columns=expected_cols)
    df_truth = pd.concat(df_list_truth, ignore_index=True) if df_list_truth else pd.DataFrame(columns=expected_cols)
    return df_dr, df_truth


# ── Merge de DataFrames parciales ─────────────────────────────────────────────

def merge_results(partial_dfs):
    """Concatena los DataFrames parciales de todos los workers."""
    valid = [df for df in partial_dfs if not df.empty]
    if valid:
        return pd.concat(valid, ignore_index=True)
    return pd.DataFrame(columns=["gen", "reco", "Gen_pid", "Reco_pid",
                                  "Gen_energy", "Reco_energy", "event_id",
                                  "Gen_Px", "Gen_Py", "Gen_Pz", "Reco_Px", "Reco_Py", "Reco_Pz",
                                  "source_file", "event_in_file"])


# ── Estructuras de asociación (vectorizado) ───────────────────────────────────

def build_association_structures(full_df, bins, e_bins):
    """
    Calcula association_results_df y energy_distribution_results a partir del
    DataFrame global usando operaciones vectorizadas sobre pandas.

    Construye claves independientes para asociación y resolución de energía:
    - key       = "{|gen_pid|}_{|reco_pid|}_{ebin_coarse}"
    - key_dist  = "{|gen_pid|}_{|reco_pid|}_{ebin_fine}"

    Devuelve
    --------
    association_results_df : dict  {key: count}
    energy_distribution_results : dict  {key: [(E_true, E_reco), ...]}
    """
    if full_df.empty:
        return {}, {}

    df = full_df.copy()
    df["_gen_pid_abs"]  = df["Gen_pid"].abs().astype(int).astype(str)
    df["_reco_pid_abs"] = df["Reco_pid"].abs().astype(int).astype(str)
    df["_pid_key"]      = df["_gen_pid_abs"] + "_" + df["_reco_pid_abs"]

    df["_ebin"]      = pd.cut(df["Gen_energy"], bins=bins).astype(str)
    df["_ebin_dist"] = pd.cut(df["Gen_energy"], bins=e_bins).astype(str)

    # key replica: "{pid}_{ebin_coarse}"
    df["_assoc_key"] = df["_pid_key"] + "_" + df["_ebin"]

    association_results_df = df.groupby("_assoc_key").size().to_dict()

    df_dist = df[np.isfinite(df["Gen_energy"]) & np.isfinite(df["Reco_energy"]) & (df["Gen_energy"] != 0)].copy()
    df_dist["_energy_pair"] = list(zip(df_dist["Gen_energy"], df_dist["Reco_energy"]))
    df_dist["_dist_key"] = df_dist["_pid_key"] + "_" + df_dist["_ebin_dist"]

    energy_distribution_results = (
        df_dist.groupby("_dist_key")["_energy_pair"]
        .apply(list)
        .to_dict()
    )

    return association_results_df, energy_distribution_results


# ── Placeholder: relleno de histogramas ROOT (implementación futura) ──────────

def fill_particle_level_histograms(full_df, root_histograms, histogram_config):
    """
    Rellena los histogramas ROOT definidos en histogram_config con los valores
    del DataFrame global full_df.

    La firma sigue el mismo patrón que analysisRHOTree.py:
        root_histograms[level][class][name].Fill(value)

    La correspondencia (Gen_pid, Reco_pid) → histograma se definirá en el YAML
    bajo la sección histograms_config, igual que en analysisRHOTree.py.
    El DataFrame full_df tiene columnas:
        event_id, Gen_pid, Reco_pid, Gen_energy, Reco_energy

    TODO: implementar el mapeo cuando se defina histograms_config en el YAML.
    """
    pass


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    default_config = "config/default/taurecolong.yaml"
    outputbasepath = "Results/TauReco/"

    def my_hook(parser):
        parser.add_argument(
            "--sys-err", type=str, default="config/systematics/err_sys.yml",
            help="YAML con errores sistemáticos",
        )
        parser.add_argument(
            "--test-extremes", action="store_true",
            help="Test de los extremos de la resolución en energía del fotón",
        )
        parser.add_argument(
            "--n-workers", type=int, default=None,
            help="Número de workers paralelos (por defecto: núcleos disponibles)",
        )
        parser.add_argument(
            "--dedup-mode",
            choices=["gen", "reco"],
            default="reco",
            help=(
                "Deduplication side for RecoMCTruthLink matching. "
                "'gen' (default): one reco per gen, max weight — legacy behaviour. "
                "'reco': one gen per reco, max weight — avoids duplicate PFOs "
                "in photon-fusion cases."
            ),
        )
        parser.add_argument(
            "--skip-gen-status-filter",
            action="store_true",
            default=False,
            help=(
                "If set, skip the generatorStatus==1 filter on MCParticles. "
                "By default only stable final-state particles are considered. "
                "Disable when intermediate particles leave calorimeter hits "
                "and their links should be preserved."
            ),
        )
        parser.add_argument(
            "--max-gen-pdg",
            type=int,
            default=10000,
            help=(
                "If set, ignore gen particles with |PDG| > max-gen-pdg. "
                "Useful to exclude nuclear fragments and other exotic codes "
                "(e.g. 1000020030). Standard hadrons have |PDG| <= ~3500; "
                "a value of 10000 excludes all exotic/nuclear codes."
            ),
        )
        parser.add_argument(
            "--weight-mode",
            choices=["raw", "decoded"],
            default="decoded",
            help=(
                "How to interpret RecoMCTruthLink weights. "
                "'raw': use the encoded value as-is for ranking. "
                "'decoded'(default): decode track/cluster weights "
                "(encodedWeight = int(clusterW*1000)*10000 + int(trackW*1000)); "
                "prefer track weight when available, fall back to cluster weight."
            ),
        )
        parser.add_argument(
            "--all-plot",
            type=int,
            nargs="*",
            default=[22],
            metavar="PDG",
            help=(
                "PDG codes for which combined_all and combined_all_true resolution "
                "plots are generated (e.g. --all-plot 22 211 13). Default: 22 (photons)."
            ),
        )

    general_configs = myutils.setup_analysis_config(
        default_config, outputbasepath, parser_hook=my_hook
    )

    loggers             = general_configs["loggers"]
    run_config          = general_configs["config"]
    neutral_recover_cfg = run_config.get("neutral_recover", {})

    logger_config  = loggers["config"]
    logger_io      = loggers["io"]
    logger_process = loggers["processing"]

    args          = general_configs["args"]
    test_pfo      = args.test_pfo
    test_extremes = args.test_extremes

    # Cortes
    dRMax         = run_config["cuts"]["dRMax"]
    minPTauPhoton = run_config["cuts"]["TauPhotonPCut"]
    minPTauPion   = run_config["cuts"]["TauPionPCut"]
    PNeutron      = run_config["cuts"]["NeutronCut"]
    generalPCut   = run_config["cuts"]["generalPCut"]

    sys_errors    = run_config.get("systematics_errors", {})
    photon_config = sys_errors.get("photon_config", {})

    outputpath        = general_configs["outputpath"]
    fileOutName       = os.path.join(outputpath, general_configs["fileOutName"])
    sample            = run_config["general"]["sample"]
    test_arg          = general_configs["flags"]["test"]
    gatr_results_path = args.gatr_result

    logger_config.info("Configuración cargada!")
    logger_config.info("Configuración:\n%s", pprint.pformat(general_configs, indent=4))

    # ── Carga de ficheros (solo en main, antes del fork) ─────────────────────
    filenames, mlpf_results = myutils.get_root_trees_path(
        sample, gatr_results_path, loggers, test_arg, args
    )
    logger_io.info("Leídos %d ficheros", len(filenames))
    logger_io.info("Primeros ficheros: %s", filenames[:10])

    if not filenames:
        logger_io.error("No se encontraron ficheros ROOT. Abortando.")
        sys.exit(1)

    # Bins para estructuras de asociación
    bins   = [0, 1, 5, 10, 20, 30, 45, 100, np.inf]

    # Binning no lineal: pocos bines a baja energía y muchos hacia 50 GeV.
    # alpha < 1 concentra los bordes en la parte alta del rango.
    n_edges = 40  # 50 bines
    e_min = 0.
    e_max = 40.0
    alpha = 0.8
    # t = np.linspace(0.0, 1.0, n_edges)
    # e_bins = e_min + (e_max - e_min) * np.power(t, alpha)
    e_bins = np.linspace(e_min, e_max, n_edges)
    # ── Bundle de configuración serializable (sin loggers, sin ROOT) ──────────
    config_bundle = {
        "neutral_recover_cfg": neutral_recover_cfg,
        "cuts": {
            "dRMax":         dRMax,
            "minPTauPhoton": minPTauPhoton,
            "minPTauPion":   minPTauPion,
            "PNeutron":      PNeutron,
            "generalPCut":   generalPCut,
        },
        "photon_config":     photon_config,
        "test_extremes":     test_extremes,
        "test_pfo":          test_pfo,
        "genparts":          "MCParticles",
        "pfobjects":         "PandoraPFOs",
        "gatr_results_path": gatr_results_path,
        "outputpath":        outputpath,
        "dedup_mode":        args.dedup_mode,
        "filter_gen_status": not args.skip_gen_status_filter,
        "max_gen_pdg":       args.max_gen_pdg,
        "weight_mode":       args.weight_mode,
    }

    # ── División del trabajo ──────────────────────────────────────────────────
    n_workers   = args.n_workers or min(len(filenames), os.cpu_count() or 1)
    n_workers   = max(1, n_workers)
    file_chunks = split_filenames(filenames, n_workers)
    mlpf_chunks = split_mlpf(mlpf_results, file_chunks)

    # Offset global de event_id por worker, consistente con el esquema de
    # indexación de mlpf_results (1000 eventos por fichero en myutils).
    event_offsets = []
    acc = 0
    for chunk in file_chunks:
        event_offsets.append(acc)
        acc += len(chunk) * 1000

    logger_io.info("Lanzando %d workers sobre %d ficheros", n_workers, len(filenames))
    for i, chunk in enumerate(file_chunks):
        logger_io.info("  Worker %d: %d ficheros, offset eventos %d",
                       i, len(chunk), event_offsets[i])

    # ── Ejecución paralela ────────────────────────────────────────────────────
    # Usamos fork (default en Linux): los workers se forkan antes de que se
    # cree cualquier TFile ROOT en el proceso principal, por lo que no hay
    # estado ROOT compartido que pueda corromperse.
    ctx = multiprocessing.get_context("fork")
    partial_dfs_dr = []
    partial_dfs_truth = []
    total_rows_dr = 0
    total_rows_truth = 0
    n_chunks = len(file_chunks)

    t_start = time.time()

    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
        futures = {
            executor.submit(
                process_chunk,
                file_chunks[i],
                mlpf_chunks[i],
                event_offsets[i],
                config_bundle,
                i,
            ): i
            for i in range(n_chunks)
        }

        for n_done, future in enumerate(as_completed(futures), start=1):
            wid = futures[future]
            try:
                df_dr, df_truth = future.result()
                total_rows_dr += len(df_dr)
                total_rows_truth += len(df_truth)
                partial_dfs_dr.append(df_dr)
                partial_dfs_truth.append(df_truth)
            except Exception as exc:
                logger_process.error("Worker %d lanzó excepción: %s", wid, exc)
                df_dr = pd.DataFrame()
                df_truth = pd.DataFrame()

            # ── Progreso en stdout ────────────────────────────────────────────
            elapsed   = time.time() - t_start
            avg_per_w = elapsed / n_done
            remaining = avg_per_w * (n_chunks - n_done)
            bar_len   = 30
            filled    = int(bar_len * n_done / n_chunks)
            bar       = "█" * filled + "░" * (bar_len - filled)
            print(
                f"\r[{bar}] {n_done}/{n_chunks} workers  "
                f"| filas dR: {total_rows_dr:>8,}  "
                f"| filas truth: {total_rows_truth:>8,}  "
                f"| transcurrido: {elapsed:>6.0f}s  "
                f"| restante estimado: {remaining:>6.0f}s  "
                f"| último worker: {wid} (dR={len(df_dr):,}, truth={len(df_truth):,})",
                end="", flush=True,
            )
            logger_process.info(
                "Worker %d completado (%d/%d) | dR=%d filas | truth=%d filas | %.0fs transcurridos",
                wid, n_done, n_chunks, len(df_dr), len(df_truth), elapsed,
            )

    print()  # salto de línea tras la barra de progreso
    total_elapsed = time.time() - t_start
    print(f"Todos los workers completados en {total_elapsed:.1f}s "
          f"({total_elapsed/60:.1f} min)")

    # ── Merge de resultados ───────────────────────────────────────────────────
    full_df_dr = merge_results(partial_dfs_dr)
    full_df_truth = merge_results(partial_dfs_truth)
    logger_io.info("Total filas tras merge (dR): %d", len(full_df_dr))
    logger_io.info("Total filas tras merge (RecoMCTruthLink): %d", len(full_df_truth))

    # ── Estructuras de asociación (vectorizado) ───────────────────────────────
    association_results_df_dr, energy_distribution_results_dr = build_association_structures(
        full_df_dr, bins, e_bins
    )
    association_results_df_truth, energy_distribution_results_truth = build_association_structures(
        full_df_truth, bins, e_bins
    )

    # ── Histogramas ROOT (creados DESPUÉS del fork, solo en main) ─────────────
    histogram_config = general_configs.get("histograms_config", {})
    root_histograms  = myutils.set_up_root_histograms(histogram_config)

    fill_particle_level_histograms(full_df_dr, root_histograms, histogram_config)

    # ── Escritura de salidas ──────────────────────────────────────────────────

    # DataFrame completo de asociaciones
    out_csv_dr = os.path.join(outputpath, "association_results_full_dR.csv")
    full_df_dr.to_csv(out_csv_dr, index=False)
    logger_io.info("DataFrame de asociaciones (dR) guardado en %s", out_csv_dr)

    out_csv_truth = os.path.join(outputpath, "association_results_full_truthlink.csv")
    full_df_truth.to_csv(out_csv_truth, index=False)
    logger_io.info("DataFrame de asociaciones (RecoMCTruthLink) guardado en %s", out_csv_truth)

    # Matrices de confusión + resolución en energía
    images_dir_dr = os.path.join(outputpath, "confusion_matrices_particle_level", "dR")
    plot_confusion_matrices(association_results_df_dr, output_dir=images_dir_dr)
    images_dir_truth = os.path.join(outputpath, "confusion_matrices_particle_level", "truthlink")
    plot_confusion_matrices(association_results_df_truth, output_dir=images_dir_truth)

    energy_dist_dir_dr = os.path.join(outputpath, "energy_distributions", "dR")
    plot_energy_distributions(energy_distribution_results_dr, output_dir=energy_dist_dir_dr, all_plot_pdgs=args.all_plot)
    energy_dist_dir_truth = os.path.join(outputpath, "energy_distributions", "truthlink")
    plot_energy_distributions(energy_distribution_results_truth, output_dir=energy_dist_dir_truth, all_plot_pdgs=args.all_plot)

    eff_dir_dr = os.path.join(outputpath, "efficiency_plots", "dR")
    plot_efficiency_vs_momentum(full_df_dr, output_dir=eff_dir_dr)
    eff_dir_truth = os.path.join(outputpath, "efficiency_plots", "truthlink")
    plot_efficiency_vs_momentum(full_df_truth, output_dir=eff_dir_truth)
    
    eff_dir_dr_theta = os.path.join(outputpath, "efficiency_plots_theta", "dR")
    plot_efficiency_vs_momentum(full_df_dr, output_dir=eff_dir_dr_theta, plot_type="theta")
    eff_dir_truth_theta = os.path.join(outputpath, "efficiency_plots_theta", "truthlink", "theta")
    plot_efficiency_vs_momentum(full_df_truth, output_dir=eff_dir_truth_theta, plot_type="theta")

    # Fichero ROOT con histogramas (descomentar cuando fill_particle_level_histograms
    # esté implementado con las secciones histograms_config del YAML)
    # outfile = ROOT.TFile(fileOutName, "RECREATE")
    # root_histograms = myutils.calc_efficiency(root_histograms, histogram_config)
    # write_histograms_recursive(root_histograms)
    # myutils.write_plot_config(root_histograms, outputpath)
    # outfile.Close()

    # Snapshot de la configuración usada
    if not isinstance(run_config["output"]["outputlabels"], list):
        if run_config["output"]["outputlabels"] is None:
            run_config["output"]["outputlabels"] = []
        else:
            run_config["output"]["outputlabels"] = [run_config["output"]["outputlabels"]]

    output_config_file = os.path.join(outputpath, "config.yaml")
    with open(output_config_file, "w") as f:
        yaml.dump(run_config, f)
    logger_io.info("Configuración guardada en %s", output_config_file)
    logger_io.info("Fin del job")


if __name__ == "__main__":
    main()
