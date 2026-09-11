#!/usr/bin/env python3
"""Atomically write a small HitAnalysis provenance sidecar."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--input-count", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--dedup-mode", required=True)
    parser.add_argument("--assoc-max-dr", type=float, required=True)
    args = parser.parse_args()
    commit = subprocess.run(
        ["git", "-C", str(args.repo_root), "rev-parse", "HEAD"],
        check=True, text=True, stdout=subprocess.PIPE,
    ).stdout.strip()
    data = {
        "schema_version": "fcc_hit_analysis_run_metadata_v1",
        "sample": args.sample,
        "source_manifest": str(args.source_manifest.resolve()),
        "source_manifest_sha256": sha256(args.source_manifest),
        "input_count": args.input_count,
        "repo_commit": commit,
        "analysis_definition_version": "fcc_hit_analysis_hardened_v1",
        "geometric_association_contract": "geometric_association_v1",
        "workflow_association_contract": "fcc_tau_association_v1",
        "workers": args.workers,
        "dedup_mode": args.dedup_mode,
        "assoc_max_dr": args.assoc_max_dr,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.partial.{os.getpid()}")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, args.output)


if __name__ == "__main__":
    main()
