"""Helpers de infraestructura compartidos por los drivers paralelos stage-1 del
análisis del rho (``analysisRHOTree_MDecs_parallel`` y ``genOnlyRHOTree_MDecs_parallel``).

Solo contiene utilidades genéricas (reparto de ficheros entre workers, logging por
worker, construcción de TLorentzVector desde constituyentes EDM4hep y merge de ROOT
files parciales). La lógica específica de cada script (CLI ``my_hook``, remapeo de
decayIDs ``_remap_to_gen``, llenado de ramas, etc.) permanece en cada script porque
difiere entre gen-only y reco.
"""

import logging
import os
import subprocess

import ROOT


def split_filenames(filenames, n_workers):
    """Reparte ``filenames`` en ``n_workers`` chunks lo más equilibrados posible."""
    k, rem = divmod(len(filenames), n_workers)
    chunks, start = [], 0
    for i in range(n_workers):
        end = start + k + (1 if i < rem else 0)
        if start < end:
            chunks.append(filenames[start:end])
        start = end
    return chunks


def _make_p4_from_const(const):
    """Construye un TLorentzVector desde una partícula EDM4hep."""
    p4 = ROOT.TLorentzVector()
    try:
        p4.SetXYZM(const.getMomentum().x, const.getMomentum().y,
                   const.getMomentum().z, const.getMass())
    except AttributeError:
        p4.SetXYZM(const.getMomentum().X(), const.getMomentum().Y(),
                   const.getMomentum().Z(), const.getMass())
    return p4


def _setup_worker_logging(outputpath, worker_id, log_source):
    """Configura logging por worker en ``outputpath/logs/<log_source>/worker_<id>.log``.

    ``log_source`` es la constante ``_LOG_SOURCE`` de cada script (distinta entre
    gen-only y reco), de modo que cada driver escribe en su propio subdirectorio.
    """
    root_logger = logging.getLogger()
    for h in root_logger.handlers[:]:
        root_logger.removeHandler(h)
        h.close()
    log_dir = os.path.join(outputpath, "logs", log_source)
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"worker_{worker_id}.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        force=True,
    )
    lg = logging.getLogger(f"worker_{worker_id}")
    return {"processing": lg, "io": lg, "config": lg}


def merge_partial_root_files(partial_files, final_output):
    """Fusiona ROOT files parciales con ``hadd``; si falla, usa ``TFileMerger``."""
    cmd = ["hadd", "-f", final_output] + partial_files
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        return
    print(f"[WARN] hadd falló (rc={result.returncode}), usando TFileMerger:\n{result.stderr}")
    merger = ROOT.TFileMerger(False)
    merger.OutputFile(final_output, "RECREATE")
    for f in partial_files:
        merger.AddFile(f)
    merger.Merge()
