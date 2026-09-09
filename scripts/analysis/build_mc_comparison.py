#!/usr/bin/env python3
"""Build the maintained frozen-definition W-versus-second-MC plot suite."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from modules.fcc_mc_comparison import load_comparison, load_sample, write_csv  # noqa: E402
from modules.fcc_mc_comparison_outputs import build_all  # noqa: E402


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--comparison", required=True,
                        choices=("whizard_p8o_18k", "whizard_kkmcee_2k"))
    result.add_argument("--config", type=Path, default=REPO / "configs/analysis/fcc_mc_comparisons_v1.yaml")
    result.add_argument("--output-root", type=Path, required=True)
    result.add_argument("--validation-mode", action="store_true",
                        help="compare generated CSV products against frozen Talk-2 references")
    result.add_argument("--overwrite", action="store_true")
    return result


def _cell_equal(left, right) -> bool:
    if left == right:
        return True
    try:
        a, b = float(left), float(right)
    except (TypeError, ValueError):
        return str(left) == str(right)
    if math.isnan(a) and math.isnan(b):
        return True
    return math.isclose(a, b, rel_tol=2e-12, abs_tol=2e-14)


def compare_csv(generated: Path, reference: Path) -> tuple[bool, int, str]:
    import csv
    with generated.open(newline="") as stream:
        new_rows = list(csv.DictReader(stream)); new_fields = list(new_rows[0]) if new_rows else []
    with reference.open(newline="") as stream:
        old_rows = list(csv.DictReader(stream)); old_fields = list(old_rows[0]) if old_rows else []
    if new_fields != old_fields:
        return False, 0, f"columns differ: {new_fields} != {old_fields}"
    if len(new_rows) != len(old_rows):
        return False, 0, f"row count differs: {len(new_rows)} != {len(old_rows)}"
    differences = 0
    first = ""
    for index, (new, old) in enumerate(zip(new_rows, old_rows), start=2):
        for field in new_fields:
            if not _cell_equal(new[field], old[field]):
                differences += 1
                if not first:
                    first = f"row {index} field {field}: {new[field]!r} != {old[field]!r}"
    return differences == 0, differences, first


def regression(generated: list[Path], output_root: Path, reference_root: Path) -> list[dict]:
    families = {"part12": 0, "part3": 0, "part3b": 0, "part4": 0}
    rows = []
    for path in generated:
        if path.suffix != ".csv":
            continue
        relative = path.relative_to(output_root)
        reference = reference_root / relative
        if not reference.is_file():
            continue
        passed, differences, detail = compare_csv(path, reference)
        family = relative.parts[0]; families[family] += 1
        rows.append({"family": family, "product": str(relative), "reference": str(reference),
                     "status": "PASS" if passed else "FAIL", "cell_discrepancies": differences,
                     "first_discrepancy": detail})
    missing = [family for family, count in families.items() if count == 0]
    if missing:
        raise AssertionError(f"no regression products checked for families: {missing}")
    return rows


def _read_rows(path: Path) -> list[dict]:
    import csv
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_summary(output_root: Path, comparison: dict, samples) -> None:
    token = comparison["output_token"]
    coverage = _read_rows(output_root / f"part4/method_comparison/W_vs_{token}_truth_association_coverage.csv")
    integrated = []
    for file_token in ("W", token):
        for method in ("G", "Ldirect", "Lancestor"):
            integrated.extend(_read_rows(output_root / f"part3/{file_token}_{method}_reco_efficiency.csv"))
    residuals = _read_rows(output_root / f"part3b/W_vs_{token}_pfo_residual_summary.csv")
    pid = _read_rows(output_root / f"part4/W_vs_{token}_conditional_pid_efficiency.csv")
    payload = {"comparison": comparison["name"],
        "event_scopes": {sample.internal_name: sample.expected_events for sample in samples},
        "pfo_coverage": coverage, "integrated_association": integrated,
        "residual_summary": residuals, "conditional_pid_diagonal": pid,
        "caveat": "Truth particle multiplicities and photon ancestry categories are generator-record dependent. No ISR/FSR interpretation is assigned at this stage."}
    destination = output_root / "summary/comparison_summary.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    lines = [f"# {samples[0].presentation_label} versus {samples[1].presentation_label}", "",
             f"Scope: {samples[0].expected_events:,} versus {samples[1].expected_events:,} events.", "",
             "The suite uses the frozen 2026-09-01 selected-truth, G, L_direct, L_ancestor, tau-origin, representative-PFO, PID and residual definitions.", "",
             "> Truth particle multiplicities and photon ancestry categories are generator-record dependent. No ISR/FSR interpretation is assigned at this stage.", "",
             "Machine-readable headline values are in `summary/comparison_summary.json`; exact inputs are in `manifest/provenance.json`.", ""]
    (output_root / "README.md").write_text("\n".join(lines))


def main() -> None:
    args = parser().parse_args()
    comparison = load_comparison(args.config, args.comparison)
    if args.validation_mode and "regression_reference" not in comparison:
        raise ValueError("validation mode is only configured for the frozen W/P8O comparison")
    if args.output_root.exists() and any(args.output_root.iterdir()) and not args.overwrite:
        raise FileExistsError(f"non-empty output root (use --overwrite): {args.output_root}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    samples = [load_sample(spec) for spec in comparison["samples"]]
    if [sample.expected_events for sample in samples] != [int(spec["expected_events"]) for spec in comparison["samples"]]:
        raise AssertionError("configured event-scope mismatch")
    generated = build_all(samples, args.output_root, comparison["output_token"])
    manifest = [{"path": str(path.relative_to(args.output_root)), "kind": path.suffix.lstrip(".")}
                for path in sorted(generated)]
    write_csv(args.output_root / "manifest/generated_products.csv", manifest)
    provenance = {"comparison": args.comparison, "scientific_contract": comparison["contract"],
                  "samples": [sample.provenance for sample in samples],
                  "event_scopes": {sample.internal_name: sample.expected_events for sample in samples},
                  "test10_scientific_use": False}
    (args.output_root / "manifest/provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    write_summary(args.output_root, comparison, samples)
    if args.validation_mode:
        validation = regression(generated, args.output_root, Path(comparison["regression_reference"]))
        write_csv(args.output_root / "validation/talk2_numerical_regression.csv", validation)
        failures = [row for row in validation if row["status"] != "PASS"]
        summary = {"status": "PASS" if not failures else "FAIL", "products_checked": len(validation),
                   "failures": len(failures), "pixel_identity_required": False}
        (args.output_root / "validation/regression_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        if failures:
            raise AssertionError(f"{len(failures)} numerical regression products failed")
    print(json.dumps({"status": "PASS", "comparison": args.comparison,
                      "generated_products": len(generated), "output_root": str(args.output_root)}, indent=2))


if __name__ == "__main__":
    main()
