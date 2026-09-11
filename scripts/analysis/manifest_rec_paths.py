#!/usr/bin/env python3
"""Validate a HitAnalysis input manifest and print deterministic REC paths."""
from __future__ import annotations

import argparse
from pathlib import Path

import yaml


SCHEMA = "fcc_hit_analysis_input_manifest_v1"


def load(path: Path) -> dict:
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict) or data.get("schema_version") != SCHEMA:
        raise ValueError(f"expected schema_version={SCHEMA}")
    if not isinstance(data.get("sample"), str) or not data["sample"]:
        raise ValueError("manifest requires a sample")
    entries = data.get("inputs")
    if not isinstance(entries, list) or not entries:
        raise ValueError("manifest requires non-empty inputs")
    seen_ids = set(); seen_paths = set()
    for entry in entries:
        if not isinstance(entry, dict) or set(("source_file_id", "path")) - set(entry):
            raise ValueError(f"invalid input entry: {entry!r}")
        source_id = str(entry["source_file_id"]); rec = str(entry["path"])
        if source_id in seen_ids or rec in seen_paths or "\n" in rec:
            raise ValueError("duplicate/invalid manifest identity or path")
        seen_ids.add(source_id); seen_paths.add(rec)
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--sample-only", action="store_true")
    args = parser.parse_args()
    data = load(args.manifest)
    if args.sample_only:
        print(data["sample"])
    else:
        for entry in data["inputs"]:
            print(entry["path"])


if __name__ == "__main__":
    main()
