#!/usr/bin/env python3
"""Dump stored MC genealogy and persisted L_direct/L_ancestor assignments."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from modules.fcc_event_record import (  # noqa: E402
    build_event_record,
    render_detail,
    render_genealogy,
    residual_match,
)
from modules.fcc_workflow_interface import validate_association_columns  # noqa: E402


def parse_events(value: str) -> list[int]:
    result: set[int] = set()
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start_text, stop_text = token.split("-", 1)
            start, stop = int(start_text), int(stop_text)
            if start < 0 or stop < start:
                raise argparse.ArgumentTypeError(f"invalid event range: {token}")
            result.update(range(start, stop + 1))
        else:
            event = int(token)
            if event < 0:
                raise argparse.ArgumentTypeError("event numbers must be non-negative")
            result.add(event)
    if not result:
        raise argparse.ArgumentTypeError("empty event selection")
    return sorted(result)


def parse_abs_range(value: str) -> tuple[float, float]:
    try:
        low_text, high_text = value.split(":", 1)
        low, high = float(low_text), float(high_text)
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError("expected MIN:MAX") from error
    if low < 0 or high < low:
        raise argparse.ArgumentTypeError("absolute residual range must satisfy 0 <= MIN <= MAX")
    return low, high


def read_assignment_rows(path: Path, level: str, source_file_id: str, events: list[int] | None) -> list[dict]:
    import pyarrow.parquet as pq

    schema = pq.ParquetFile(path).schema_arrow
    validate_association_columns(schema.names, level)
    filters = [("source_file_id", "=", str(source_file_id))]
    if events is not None:
        filters.append(("event_in_file", "in", list(map(int, events))))
    return pq.read_table(path, filters=filters).to_pylist()


def rows_by_event(rows: list[dict]) -> dict[int, list[dict]]:
    result: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        result[int(row["event_in_file"])].append(row)
    return result


def write_record(record: dict, output_dir: Path) -> list[Path]:
    stem = f"event_{record['source_file_id']}_{record['event_in_file']}"
    paths = [output_dir / f"{stem}.detail.txt", output_dir / f"{stem}.tree.txt", output_dir / f"{stem}.json"]
    existing = [path for path in paths if path.exists()]
    if existing:
        raise FileExistsError(f"refusing to overwrite: {existing[0]}")
    paths[0].write_text(render_detail(record))
    paths[1].write_text(render_genealogy(record))
    paths[2].write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return paths


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Inspect stored MCParticles genealogy and existing L_direct/L_ancestor assignments."
    )
    result.add_argument("--rec", required=True, type=Path)
    result.add_argument("--ldirect", required=True, type=Path)
    result.add_argument("--lancestor", required=True, type=Path)
    result.add_argument("--source-file-id", required=True)
    result.add_argument("--output-dir", required=True, type=Path)
    selection = result.add_mutually_exclusive_group()
    selection.add_argument("--first-events", type=int)
    selection.add_argument("--events", type=parse_events)
    selection.add_argument("--event-key")
    result.add_argument("--find-pdg", type=int)
    result.add_argument("--find-method", choices=("L_direct", "L_ancestor"))
    result.add_argument("--abs-theta-residual", type=parse_abs_range)
    result.add_argument(
        "--representative-only", action="store_true",
        help="restrict a residual search to the maintained unique representative PFO",
    )
    result.add_argument("--max-events", type=int, default=20)
    return result


def main() -> int:
    args = parser().parse_args()
    for path in (args.rec, args.ldirect, args.lancestor):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if args.first_events is not None and args.first_events < 1:
        raise ValueError("--first-events must be positive")
    if args.max_events < 1:
        raise ValueError("--max-events must be positive")

    search_values = (args.find_pdg, args.find_method, args.abs_theta_residual)
    searching = any(value is not None for value in search_values)
    if searching and not all(value is not None for value in search_values):
        raise ValueError("residual search requires --find-pdg, --find-method and --abs-theta-residual")
    if args.representative_only and not searching:
        raise ValueError("--representative-only requires a residual search")
    if searching and (args.events is not None or args.event_key is not None or args.first_events is not None):
        raise ValueError("residual search cannot be combined with explicit/default event selectors")

    selected: list[int] | None
    if searching:
        selected = None
    elif args.event_key:
        try:
            key_source, event_text = args.event_key.split(":", 1)
            event = int(event_text)
        except ValueError as error:
            raise ValueError("--event-key must be SOURCE_FILE_ID:EVENT_IN_FILE") from error
        if key_source != str(args.source_file_id) or event < 0:
            raise ValueError("--event-key does not match --source-file-id or has a negative event")
        selected = [event]
    elif args.events is not None:
        selected = args.events
    else:
        count = 10 if args.first_events is None else args.first_events
        selected = list(range(count))

    direct = rows_by_event(read_assignment_rows(args.ldirect, "direct", args.source_file_id, selected))
    ancestor = rows_by_event(read_assignment_rows(args.lancestor, "ancestor", args.source_file_id, selected))
    args.output_dir.mkdir(parents=True, exist_ok=True)

    from podio import root_io

    wanted = None if selected is None else set(selected)
    found: set[int] = set()
    written: list[Path] = []
    low, high = args.abs_theta_residual if searching else (None, None)
    reader = root_io.Reader(str(args.rec))
    for event_index, event in enumerate(reader.get("events")):
        if wanted is not None and event_index not in wanted:
            if wanted and event_index > max(wanted):
                break
            continue
        record = build_event_record(
            event,
            source_file_id=str(args.source_file_id),
            event_in_file=event_index,
            direct_rows=direct.get(event_index, []),
            ancestor_rows=ancestor.get(event_index, []),
        )
        if searching:
            matches = []
            for pfo in record["pandora_pfos"]:
                match = residual_match(
                    pfo,
                    record["mc_particles"],
                    pdg=args.find_pdg,
                    method=args.find_method,
                    minimum_abs_mrad=low,
                    maximum_abs_mrad=high,
                    representative_only=args.representative_only,
                )
                if match is not None:
                    key = "ldirect" if args.find_method == "L_direct" else "lancestor"
                    pfo[key]["search_match"] = True
                    matches.append(match)
            if not matches:
                continue
            record["selection_matches"] = matches
        written.extend(write_record(record, args.output_dir))
        found.add(event_index)
        if searching and len(found) >= args.max_events:
            break
        if wanted is not None and found == wanted:
            break

    if wanted is not None and found != wanted:
        raise IndexError(f"requested events not found in REC: {sorted(wanted - found)}")
    print(f"wrote {len(found)} event record(s), {len(written)} file(s), to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
