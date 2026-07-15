#!/usr/bin/env python3
"""
Pipeline MDecs: 1 generación de árbol + N ejecuciones de histogramas.

Orquesta los scripts existentes (no los modifica):
  - genOnlyRHOTree_MDecs_parallel.py   (tree.mode: gen)
  - analysisRHOTree_MDecs_parallel.py  (tree.mode: reco)
  - RhoHistFromTree_MDecs_parallel.py  (paso de histogramas, 1..N runs)

Toda la configuración viene de un YAML (ver
config/pipeline/tree_hist_pipeline_example.yaml). El paso de árbol es
opcional: con `tree.enabled: false` (o --hist-only) solo se generan
histogramas a partir de `pipeline.tree_file`.

Uso:
  python RhoAnalysis/runTreeHistPipeline_MDecs.py \
      --pipeline-config config/pipeline/tree_hist_pipeline_example.yaml [--dry-run]
"""

import argparse
import glob
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

import yaml

_REPO_ROOT   = Path(__file__).resolve().parent.parent
_OUTPUT_BASE = "Results/RhoAnalysis/"

_TREE_SCRIPTS = {
    "gen":  "RhoAnalysis/genOnlyRHOTree_MDecs_parallel.py",
    "reco": "RhoAnalysis/analysisRHOTree_MDecs_parallel.py",
}
_HIST_SCRIPT = "RhoAnalysis/RhoHistFromTree_MDecs_parallel.py"

# Claves de cortes estándar que admite setup_analysis_config como flag CLI.
_CUT_KEYS = ["tauCut", "dRMax", "TauPhotonPCut", "TauPionPCut",
             "NeutronCut", "generalPCut", "MatchedGenMinDR"]

# ── Mapeo clave YAML → (flag CLI, tipo) ──────────────────────────────────────
# Tipos: "scalar" (un valor), "list" (nargs+), "flag" (booleano store_true).
# Solo se exponen flags con efecto real; los flags muertos de
# setup_analysis_config (--decay, --gatr-result, --test-pfo, --matchedCM,
# --test, --hist-config) y los deprecated/unused de los scripts
# (--decay-pair de analysisRHOTree, --test-extremes) quedan fuera.

_TREE_COMMON_SPEC = {
    "sample":         ("--sample",         "scalar"),
    "prefix":         ("--prefix",         "scalar"),
    "config":         ("--config",         "scalar"),
    "outfile":        ("--outfile",        "scalar"),
    "n_workers":      ("--n-workers",      "scalar"),
    "decay_modes":    ("--decay-modes",    "list"),
    "sin_eff":        ("--sin-eff",        "scalar"),
    "samples_config": ("--samples-config", "scalar"),
    "input_list":     ("--input-list",     "list"),
}

_TREE_GEN_SPEC = {
    "helicity": ("--helicity", "flag"),
}

_TREE_RECO_SPEC = {
    "electron_cut": ("--electron-cut", "scalar"),
    "muon_cut":     ("--muon-cut",     "scalar"),
    "lepton_xor_p": ("--lepton-xor-p", "scalar"),
    "sys_err":      ("--sys-err",      "scalar"),
}

_HIST_SPEC = {
    "hist_config_mdecs": ("--hist-config-mdecs", "scalar"),
    "config":            ("--config",            "scalar"),
    "n_workers":         ("--n-workers",         "scalar"),
    "decay_pair":        ("--decay-pair",        "list"),
    "single_decay":      ("--single-decay",      "scalar"),
    "sin_eff":           ("--sin-eff",           "scalar"),
    "compute_weights":   ("--compute-weights",   "flag"),
    "no_omega_weights":  ("--no-omega-weights",  "flag"),
    "cos_acceptance":    ("--cos-acceptance",    "scalar"),
    "omega_border_cut":  ("--omega-border-cut",  "flag"),
    "ang":               ("--ang",               "list"),
    "meson_cut":         ("--meson-cut",         "list"),
    "lepton_cut":        ("--lepton-cut",        "list"),
    "zmass_cut":         ("--zmass-cut",         "list"),
    "cut":               ("--cut",               "list"),
    "merge_hemispheres": ("--merge-hemispheres", "flag"),
    "vism_cut":          ("--vism-cut",          "scalar"),
    "use_nn_optimal":    ("--use-nn-optimal",    "flag"),
    "nn_model":          ("--nn-model",          "scalar"),
}

# Claves extra (gestionadas a mano, no son flags directos)
_TREE_EXTRA_KEYS = {"enabled", "mode", "cuts", "verbose"}
_HIST_EXTRA_KEYS = {"name", "only_gen", "verbose"}


# ── Carga y validación del YAML ──────────────────────────────────────────────

def _check_keys(section, allowed, where):
    unknown = set(section) - allowed
    if unknown:
        sys.exit(f"[pipeline] Claves desconocidas en {where}: {sorted(unknown)}\n"
                 f"           Permitidas: {sorted(allowed)}")


def load_pipeline_config(path):
    if not os.path.isfile(path):
        sys.exit(f"[pipeline] No existe el config: {path}")
    with open(path) as f:
        raw = yaml.safe_load(f)
    cfg = (raw or {}).get("pipeline")
    if cfg is None:
        sys.exit(f"[pipeline] El YAML debe tener una sección raíz 'pipeline:' ({path})")

    _check_keys(cfg, {"tree", "tree_file", "hist_common", "hist_runs"}, "pipeline")

    tree = cfg.get("tree") or {}
    _check_keys(tree,
                set(_TREE_COMMON_SPEC) | set(_TREE_GEN_SPEC) | set(_TREE_RECO_SPEC)
                | _TREE_EXTRA_KEYS,
                "pipeline.tree")
    tree.setdefault("enabled", True)
    tree.setdefault("mode", "gen")
    if tree["mode"] not in _TREE_SCRIPTS:
        sys.exit(f"[pipeline] tree.mode debe ser 'gen' o 'reco' (recibido: {tree['mode']})")
    cuts = tree.get("cuts") or {}
    _check_keys(cuts, set(_CUT_KEYS), "pipeline.tree.cuts")

    hist_common = cfg.get("hist_common") or {}
    _check_keys(hist_common, (set(_HIST_SPEC) | _HIST_EXTRA_KEYS) - {"name"},
                "pipeline.hist_common")

    runs = cfg.get("hist_runs") or []
    if not isinstance(runs, list) or not runs:
        sys.exit("[pipeline] pipeline.hist_runs debe ser una lista con al menos un run")
    merged_runs = []
    for i, run in enumerate(runs):
        name = run.get("name", f"run{i}")
        _check_keys(run, set(_HIST_SPEC) | _HIST_EXTRA_KEYS, f"hist_runs[{name}]")
        merged = {**hist_common, **run}
        merged["name"] = name
        # Sin config explícito, heredar el del paso de árbol (el default
        # interno del script, taurecolong.yaml, ya no existe en el repo)
        if merged.get("config") is None:
            merged["config"] = tree.get("config")
        has_pair   = merged.get("decay_pair") is not None
        has_single = merged.get("single_decay") is not None
        if has_pair == has_single:
            sys.exit(f"[pipeline] hist_runs[{name}]: especifica exactamente uno de "
                     f"'decay_pair' o 'single_decay'")
        if has_pair and len(merged["decay_pair"]) != 2:
            sys.exit(f"[pipeline] hist_runs[{name}]: decay_pair debe tener 2 elementos")
        # Si el run define single_decay, anula un decay_pair heredado (y viceversa)
        merged["decay_pair" if has_single else "single_decay"] = None
        merged_runs.append(merged)
    cfg["hist_runs"] = merged_runs
    cfg["tree"] = tree
    return cfg


# ── Construcción de comandos ─────────────────────────────────────────────────

def _append_opts(cmd, section, spec):
    for key, (flag, kind) in spec.items():
        val = section.get(key)
        if val is None:
            continue
        if kind == "flag":
            if val:
                cmd.append(flag)
        elif kind == "list":
            if val:  # lista vacía → sin flag (nargs='+' no admite 0 valores)
                cmd.append(flag)
                cmd.extend(str(v) for v in val)
        else:
            cmd.extend([flag, str(val)])


def _append_verbose(cmd, section):
    v = int(section.get("verbose") or 0)
    if v:
        cmd.append("-" + "v" * min(v, 2))


def build_tree_cmd(tree):
    cmd = [sys.executable, _TREE_SCRIPTS[tree["mode"]]]
    _append_opts(cmd, tree, _TREE_COMMON_SPEC)
    _append_opts(cmd, tree, _TREE_GEN_SPEC if tree["mode"] == "gen" else _TREE_RECO_SPEC)
    for key, val in (tree.get("cuts") or {}).items():
        cmd.extend([f"--{key}", str(val)])
    _append_verbose(cmd, tree)
    return cmd


def build_hist_cmd(run, tree_file, only_gen):
    cmd = [sys.executable, _HIST_SCRIPT, "--tree-file", tree_file]
    if only_gen:
        cmd.append("--only-gen")
    _append_opts(cmd, run, _HIST_SPEC)
    _append_verbose(cmd, run)
    return cmd


def resolve_only_gen(run, tree):
    """only_gen: 'auto' → True si el árbol es gen-only (regla crítica:
    árboles de genOnlyRHOTree SIEMPRE con --only-gen; de analysisRHOTree nunca)."""
    val = run.get("only_gen", "auto")
    if val in (None, "auto"):
        return tree["mode"] == "gen"
    return bool(val)


# ── Localización del árbol generado ──────────────────────────────────────────
# Reproduce la lógica de nombres de myutils.setup_analysis_config
# (modules/myutils.py:625-668) sin ejecutarla.

def _first(val):
    return val[0] if isinstance(val, list) else val


def predict_tree_path(tree):
    base_cfg_path = tree.get("config") or "config/default/taurecolong.yaml"
    with open(_REPO_ROOT / base_cfg_path) as f:
        base_cfg = yaml.safe_load(f)

    cuts = dict(base_cfg.get("cuts") or {})
    # Los overrides CLI pasan por argparse type=float → 3 se vuelve 3.0 en el sufijo
    for key, val in (tree.get("cuts") or {}).items():
        cuts[key] = float(val)

    tph = _first(cuts.get("TauPhotonPCut"))
    tpi = _first(cuts.get("TauPionPCut"))
    npe = _first(cuts.get("NeutronCut"))
    gpc = _first(cuts.get("generalPCut"))
    dr  = _first(cuts.get("dRMax"))
    suffix = f"_{dr}_tph{tph}_tpi{tpi}_n{npe}_g{gpc}"

    outfile = tree.get("outfile") or base_cfg["general"].get("outfile")
    decay   = (base_cfg.get("general", {}).get("decay") or [-777])[0]
    decay_str = ("decayAll" if decay == -777 else f"decay{decay}") + suffix

    outdir = _REPO_ROOT / _OUTPUT_BASE / f"{tree.get('prefix') or ''}{outfile}{suffix[1:]}"

    if tree["mode"] == "gen":
        fname = f"{outfile}{decay_str}.root"
    else:
        decay_modes = tree.get("decay_modes")
        decay_tag = "All" if not decay_modes else "_".join(str(d) for d in decay_modes)
        fname = f"TTree_MDecs_{decay_tag}_{outfile}{decay_str}.root"
    return outdir, outdir / fname


def locate_tree_file(tree, t_start):
    outdir, predicted = predict_tree_path(tree)
    if predicted.is_file():
        return str(predicted)
    print(f"[pipeline] Árbol predicho no encontrado: {predicted}\n"
          f"           Buscando en {outdir} ...")
    pattern = "TTree_MDecs_*.root" if tree["mode"] == "reco" else "*decay*.root"
    candidates = [
        p for p in glob.glob(str(outdir / pattern))
        if not os.path.basename(p).startswith("HistosMDecs_")
        and "partial" not in os.path.basename(p)
        and os.path.getmtime(p) >= t_start
    ]
    if not candidates:
        sys.exit(f"[pipeline] No se encontró ningún árbol generado en {outdir} "
                 f"(patrón {pattern}, mtime posterior al inicio del paso 1). Abortando.")
    newest = max(candidates, key=os.path.getmtime)
    print(f"[pipeline] Usando el más reciente: {newest}")
    return newest


# ── Ejecución ────────────────────────────────────────────────────────────────

def run_cmd(cmd, label, dry_run):
    print(f"\n{'─' * 79}\n[pipeline] {label}:\n  {shlex.join(cmd)}\n{'─' * 79}",
          flush=True)
    if dry_run:
        return 0, 0.0
    t0 = time.time()
    result = subprocess.run(cmd, cwd=_REPO_ROOT)
    return result.returncode, time.time() - t0


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline MDecs: generación de árbol + N runs de histogramas",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pipeline-config", required=True,
                        help="YAML del pipeline (config/pipeline/*.yaml)")
    parser.add_argument("--hist-only", action="store_true",
                        help="Salta el paso de árbol aunque tree.enabled sea true")
    parser.add_argument("--tree-file", type=str, default=None,
                        help="Override de pipeline.tree_file")
    parser.add_argument("--runs", nargs="+", metavar="NAME", default=None,
                        help="Ejecuta solo los hist_runs con estos nombres")
    parser.add_argument("--dry-run", action="store_true",
                        help="Imprime los comandos sin ejecutarlos")
    args = parser.parse_args()

    cfg  = load_pipeline_config(args.pipeline_config)
    tree = cfg["tree"]
    runs = cfg["hist_runs"]

    if args.runs:
        known = {r["name"] for r in runs}
        missing = set(args.runs) - known
        if missing:
            sys.exit(f"[pipeline] --runs desconocidos: {sorted(missing)}. "
                     f"Disponibles: {sorted(known)}")
        runs = [r for r in runs if r["name"] in args.runs]

    do_tree   = tree["enabled"] and not args.hist_only
    tree_file = args.tree_file or cfg.get("tree_file")
    if not do_tree and not tree_file and not args.dry_run:
        sys.exit("[pipeline] Modo solo-histogramas sin árbol: indica "
                 "pipeline.tree_file en el YAML o --tree-file.")

    summary = []

    # Paso 1: árbol
    t_start = time.time()
    if do_tree:
        rc, dt = run_cmd(build_tree_cmd(tree), f"PASO 1 — árbol ({tree['mode']})",
                         args.dry_run)
        summary.append((f"tree ({tree['mode']})", rc, dt))
        if rc != 0:
            print(f"\n[pipeline] El paso de árbol falló (rc={rc}). Abortando.")
            _print_summary(summary, tree_file)
            sys.exit(rc)

    # Localizar el árbol para los histogramas
    if not tree_file:
        if args.dry_run:
            _, predicted = predict_tree_path(tree)
            tree_file = str(predicted)
            print(f"\n[pipeline] (dry-run) Árbol esperado: {tree_file}")
        else:
            tree_file = locate_tree_file(tree, t_start)
    if not args.dry_run and not os.path.isfile(tree_file):
        sys.exit(f"[pipeline] El árbol no existe: {tree_file}")

    # Paso 2: histogramas (un fallo no detiene los siguientes runs)
    for run in runs:
        only_gen = resolve_only_gen(run, tree)
        cmd = build_hist_cmd(run, tree_file, only_gen)
        rc, dt = run_cmd(cmd, f"PASO 2 — histogramas [{run['name']}]", args.dry_run)
        summary.append((f"hist [{run['name']}]", rc, dt))
        if rc != 0:
            print(f"[pipeline] hist [{run['name']}] falló (rc={rc}); "
                  f"continúo con los siguientes runs.")

    _print_summary(summary, tree_file)
    sys.exit(1 if any(rc != 0 for _, rc, _ in summary) else 0)


def _print_summary(summary, tree_file):
    print(f"\n{'═' * 79}\n[pipeline] Resumen  (árbol: {tree_file})")
    for label, rc, dt in summary:
        status = "OK  " if rc == 0 else f"FAIL(rc={rc})"
        print(f"  {status:12s} {label:40s} {dt:8.1f}s")
    print("═" * 79)


if __name__ == "__main__":
    main()
