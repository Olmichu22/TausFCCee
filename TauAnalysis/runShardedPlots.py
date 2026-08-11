#!/usr/bin/env python
"""Run plotTausLongResults.py sharded across N processes and merge the output.

The analysis script itself is not modified: each shard is a plain subprocess
that receives ``--shard I/N`` (handled in ``myutils.get_root_trees_path``) plus
a ``shardI_`` output prefix, so the N runs never touch each other's files.

Merging is not a plain ``hadd``. The TH1/TH2 histograms are additive and do
merge that way, but the ``Effi`` entries are ``TGraphAsymmErrors`` — ratios,
which cannot be summed. They are recomputed here from the merged numerator and
denominator histograms.

Each run stages its shards under ``Results/TauReco/_shards_<prefix>/``, so
concurrent runs are independent as long as they are given different
``--prefix`` values. Nothing outside that staging directory is ever deleted.

Usage:
    python TauAnalysis/runShardedPlots.py -j 16 \
        --hist-config config/histograms/tau_long_results.yml [...plot args]

Every argument this script does not recognise is forwarded verbatim to
plotTausLongResults.py.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import time

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import ROOT

ROOT.gROOT.SetBatch(True)

from modules import myutils

DEFAULT_SCRIPT = "TauAnalysis/plotTausLongResults.py"
DEFAULT_HIST_CONFIG = "config/histograms/tau_long_results.yml"
# setup_analysis_config builds "Results/TauReco/<prefix><outfile><cuts>/", and
# prepends "GATr_" to the whole string when --gatr-result is given.
OUTPUT_BASE = "Results/TauReco/"
GATR_OUTPUT_BASE = "GATr_Results/TauReco/"
# Per-shard output goes into a staging directory of its own, one per --prefix,
# so two launcher runs with different prefixes never see each other's files.
# --prefix accepts a path separator, and setup_analysis_config's makedirs
# creates the nested path, so this needs no change on the analysis side.
STAGING_FMT = "_shards_{scope}"
SHARD_PREFIX = "shard"
MERGED_PREFIX = "merged_"
# Histogram-name suffixes produced by --test-extremes (see
# myutils.clone_histograms_with_suffix). "" is the nominal set.
VARIATION_SUFFIXES = ["", "_min", "_max"]


# ── Launching ─────────────────────────────────────────────────────────────────

def take_option(passthrough, name):
    """Pop ``--name VALUE`` / ``--name=VALUE`` out of *passthrough*, or return None."""
    for i, token in enumerate(passthrough):
        if token == name:
            if i + 1 >= len(passthrough):
                return None
            value = passthrough[i + 1]
            del passthrough[i:i + 2]
            return value
        if token.startswith(name + "="):
            value = token.split("=", 1)[1]
            del passthrough[i]
            return value
    return None


def peek_option(passthrough, name):
    """Read ``--name VALUE`` / ``--name=VALUE`` without removing it."""
    return take_option(list(passthrough), name)


def staging_dir(user_prefix):
    """Directory holding this run's per-shard output, scoped by --prefix."""
    return os.path.join(OUTPUT_BASE, STAGING_FMT.format(scope=user_prefix or "default"))


def shard_dirs(staging, shard_index):
    """Existing output directories belonging to a given shard index.

    Restricted to directories so the sibling ``shardI_launcher.log`` files,
    which share the prefix, are not mistaken for output.
    """
    pattern = os.path.join(staging, f"{SHARD_PREFIX}{shard_index}_*")
    return sorted(p for p in glob.glob(pattern) if os.path.isdir(p))


def launch_shards(script, n_shards, passthrough, python_exe, user_prefix, staging):
    """Start one subprocess per shard and wait for all of them.

    Returns the list of per-shard output directories, in shard order.
    """
    staging_leaf = os.path.basename(staging)
    os.makedirs(staging, exist_ok=True)
    procs = []
    for i in range(n_shards):
        cmd = [python_exe, script, "--shard", f"{i}/{n_shards}"]
        # The prefix nests the run under its staging directory and tags the
        # shard; any --prefix the user passed is kept, appended after ours.
        cmd += ["--prefix", f"{staging_leaf}/{SHARD_PREFIX}{i}_{user_prefix}"]
        cmd += passthrough
        print(f"[shard {i}] {' '.join(cmd)}")
        log_path = os.path.join(staging, f"{SHARD_PREFIX}{i}_launcher.log")
        log_file = open(log_path, "w")
        procs.append((i, subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.STDOUT), log_file))

    t_start = time.time()
    failed = []
    done = 0
    for i, proc, log_file in procs:
        rc = proc.wait()
        log_file.close()
        done += 1
        elapsed = time.time() - t_start
        status = "ok" if rc == 0 else f"FAILED (exit {rc})"
        print(f"[shard {i}] {status} | {done}/{len(procs)} done | {elapsed:.0f}s elapsed")
        if rc != 0:
            failed.append(i)

    if failed:
        print(f"\n{len(failed)} shard(s) failed: {failed}")
        for i in failed:
            print(f"  log: {os.path.join(staging, f'{SHARD_PREFIX}{i}_launcher.log')}")
        sys.exit(1)

    dirs = []
    for i in range(n_shards):
        matches = shard_dirs(staging, i)
        if len(matches) != 1:
            print(f"Expected exactly one output directory for shard {i}, found {matches}")
            sys.exit(1)
        dirs.append(matches[0])
    return dirs


# ── Merging ───────────────────────────────────────────────────────────────────

def single_file(directory, pattern, what):
    matches = sorted(glob.glob(os.path.join(directory, pattern)))
    if len(matches) != 1:
        print(f"Expected exactly one {what} in {directory}, found {matches}")
        sys.exit(1)
    return matches[0]


def hadd_shards(shard_dirs_list, raw_path):
    """Sum the per-shard ROOT files. Correct for TH1/TH2, bogus for the TGraphs."""
    inputs = [single_file(d, "*.root", "ROOT output") for d in shard_dirs_list]
    cmd = ["hadd", "-f", raw_path] + inputs
    print(f"\nhadd {len(inputs)} file(s) → {raw_path}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"hadd failed (exit {result.returncode}):\n{result.stderr}")
        sys.exit(1)
    return inputs


def read_histograms(root_file, histogram_config, suffix):
    """Rebuild the nested root_histograms dict from a merged file.

    Keys mirror the YAML config (that is what calc_efficiency indexes by),
    while the objects are looked up by their ROOT name plus *suffix*.
    Returns None when this variation is absent from the file.
    """
    histograms = {}
    found = 0
    for level, categories in histogram_config.items():
        histograms[level] = {}
        for hist_class, entries in categories.items():
            histograms[level][hist_class] = {}
            if entries is None or hist_class == "Effi":
                continue
            for key, params in entries.items():
                obj = root_file.Get(params.get("name", key) + suffix)
                if not obj:
                    continue
                obj.SetDirectory(0)
                histograms[level][hist_class][key] = obj
                found += 1
    return histograms if found else None


def write_histograms_recursive(obj, seen):
    if isinstance(obj, dict):
        for value in obj.values():
            write_histograms_recursive(value, seen)
        return
    name = obj.GetName()
    if name in seen:
        return
    seen.add(name)
    obj.Write()


def merge_root(shard_dirs_list, histogram_config, out_dir, out_name):
    """hadd the shard files, then recompute the efficiencies on the totals."""
    raw_path = os.path.join(out_dir, "merged_raw.root")
    hadd_shards(shard_dirs_list, raw_path)

    raw = ROOT.TFile.Open(raw_path)
    variations = {}
    for suffix in VARIATION_SUFFIXES:
        histograms = read_histograms(raw, histogram_config, suffix)
        if histograms is not None:
            variations[suffix] = histograms
    raw.Close()

    if not variations:
        print(f"No histograms from {DEFAULT_HIST_CONFIG} were found in {raw_path}.")
        sys.exit(1)

    # A fresh file, so the meaningless hadd'ed TGraphs are dropped rather than
    # left behind as a second key with the same name.
    out_path = os.path.join(out_dir, out_name)
    outfile = ROOT.TFile(out_path, "RECREATE")
    seen = set()
    for suffix, histograms in variations.items():
        label = "nominal" if suffix == "" else suffix.lstrip("_")
        print(f"Recomputing efficiencies on the merged histograms ({label})")
        histograms = myutils.calc_efficiency(histograms, histogram_config, suffix)
        write_histograms_recursive(histograms, seen)
    outfile.Close()
    os.remove(raw_path)
    print(f"Merged histograms → {out_path}")
    return out_path


def merge_csvs(shard_dirs_list, out_dir):
    """Concatenate the per-shard label CSVs, one output per distinct filename.

    A ``shard`` column is added because GenID is built from a per-run event
    counter (``str(eventid) + str(tau_index)``), which restarts in every shard;
    the pair (shard, GenID) is what identifies a row after merging.
    """
    by_name = {}
    for shard_index, directory in enumerate(shard_dirs_list):
        for path in sorted(glob.glob(os.path.join(directory, "true_predicted_label_*.csv"))):
            by_name.setdefault(os.path.basename(path), []).append((shard_index, path))

    for name, entries in by_name.items():
        frames = []
        for shard_index, path in entries:
            frame = pd.read_csv(path)
            frame["shard"] = shard_index
            frames.append(frame)
        merged = pd.concat(frames, ignore_index=True)
        paths = [p for _, p in entries]
        out_path = os.path.join(out_dir, name)
        merged.to_csv(out_path, index=False)
        print(f"Merged {len(paths)} CSV(s), {len(merged)} rows → {out_path}")


def copy_side_files(reference_dir, out_dir):
    """Copy the config/plot-config written by every shard. They are identical."""
    for pattern in ("plot_config*.yaml", "config.yaml"):
        for path in sorted(glob.glob(os.path.join(reference_dir, pattern))):
            shutil.copy2(path, os.path.join(out_dir, os.path.basename(path)))


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Shard plotTausLongResults.py across processes and merge the results.",
        epilog="Unrecognised arguments are forwarded to the analysis script.",
    )
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 1,
                        help="Number of shards / parallel processes (default: all CPUs)")
    parser.add_argument("--script", default=DEFAULT_SCRIPT,
                        help=f"Analysis script to shard (default: {DEFAULT_SCRIPT})")
    parser.add_argument("--overwrite", action="store_true",
                        help="Delete this run's pre-existing staging directory before running")
    parser.add_argument("--clean", action="store_true",
                        help="Delete this run's staging directory after a successful merge")
    parser.add_argument("--merge-only", action="store_true",
                        help="Skip the run and merge the shard directories already on disk")
    parser.add_argument("--python", default=sys.executable,
                        help="Python interpreter used for the shard subprocesses")
    args, passthrough = parser.parse_known_args()

    if args.jobs < 1:
        parser.error("-j/--jobs must be >= 1")

    # The launcher needs the same histogram config the shards will use, to know
    # which names to read back and how to rebuild the efficiencies.
    hist_config_path = peek_option(passthrough, "--hist-config") or DEFAULT_HIST_CONFIG
    if not os.path.exists(hist_config_path):
        parser.error(f"Histogram config not found: {hist_config_path}")
    with open(hist_config_path) as f:
        histogram_config = yaml.safe_load(f)

    # --prefix is consumed here and re-emitted with the shard tag in front,
    # so the two do not fight over the same argparse destination.
    user_prefix = take_option(passthrough, "--prefix") or ""

    global OUTPUT_BASE
    if peek_option(passthrough, "--gatr-result"):
        OUTPUT_BASE = GATR_OUTPUT_BASE
    os.makedirs(OUTPUT_BASE, exist_ok=True)

    staging = staging_dir(user_prefix)
    if not user_prefix:
        print("Warning: no --prefix given, so this run stages under "
              f"{staging}. Concurrent runs without a distinct --prefix will "
              "overwrite each other.")

    if args.merge_only:
        dirs = []
        for i in range(args.jobs):
            matches = shard_dirs(staging, i)
            if len(matches) != 1:
                parser.error(f"--merge-only: expected one directory for shard {i} "
                             f"under {staging}, found {matches}")
            dirs.append(matches[0])
    else:
        # Only this run's staging directory is ever removed. The merged output
        # is left alone: it is rewritten in place, and deleting by glob here is
        # what previously let two concurrent runs destroy each other's results.
        if os.path.isdir(staging) and not args.overwrite:
            print(f"Staging directory already exists: {staging}")
            print("Re-run with --overwrite to replace it, or move it aside.")
            sys.exit(1)
        shutil.rmtree(staging, ignore_errors=True)

        print(f"Launching {args.jobs} shard(s) of {args.script}\n")
        t_start = time.time()
        dirs = launch_shards(args.script, args.jobs, passthrough, args.python,
                             user_prefix, staging)
        print(f"\nAll shards done in {time.time() - t_start:.0f}s")

    # The merged directory mirrors the shard-0 name, with the shard tag swapped.
    shard0_name = os.path.basename(dirs[0].rstrip("/"))
    out_dir = os.path.join(OUTPUT_BASE,
                           MERGED_PREFIX + shard0_name[len(f"{SHARD_PREFIX}0_"):])
    os.makedirs(out_dir, exist_ok=True)

    out_name = os.path.basename(single_file(dirs[0], "*.root", "ROOT output"))
    merge_root(dirs, histogram_config, out_dir, out_name)
    merge_csvs(dirs, out_dir)
    copy_side_files(dirs[0], out_dir)

    if args.clean:
        # The whole staging directory, logs included — nothing outside it.
        shutil.rmtree(staging, ignore_errors=True)
        print(f"Removed staging directory {staging}")

    print(f"\nOutput: {out_dir}")


if __name__ == "__main__":
    main()
