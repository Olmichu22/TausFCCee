#!/usr/bin/env python3
"""Compare generic-audit W/P8O PFO bookkeeping with frozen Talk-2 anchors."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path


def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-csv", type=Path, required=True)
    parser.add_argument("--coverage-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--p8c-outcomes-csv", type=Path)
    args = parser.parse_args()
    if args.output_csv.exists():
        raise FileExistsError(args.output_csv)
    audit = {row["sample_internal_name"]: row for row in rows(args.audit_csv)}
    frozen = rows(args.coverage_csv)
    labels = {"W": "WHIZARD", "P8O": "PYTHIA8"}
    fields = {"G": None, "Ldirect": "ldirect_direct_mc_selected", "Lancestor": "lancestor_final_unique_valid_selected"}
    output = []
    for sample, label in labels.items():
        for method, field in fields.items():
            old = next(row for row in frozen if row["sample"] == label and row["association_method"] == method)
            checks = [("total_pfos", int(old["N_all_PFO"]), int(audit[sample]["total_pfos"]))]
            if field:
                checks.append((f"{method}_usable_selected_truth", int(old["N_usable_truth_link"]), int(audit[sample][field])))
            for quantity, historical, recomputed in checks:
                output.append({"sample_internal_name": sample, "association_method": method,
                               "quantity": quantity, "historical": historical,
                               "recomputed": recomputed, "exact": historical == recomputed,
                               "historical_source": str(args.coverage_csv)})
    if args.p8c_outcomes_csv:
        species_fields = {"electron":"selected_electrons", "muon":"selected_muons",
                          "photon":"selected_photons", "charged_pion":"selected_charged_pions"}
        for old in rows(args.p8c_outcomes_csv):
            field = species_fields[old["truth_species"]]
            historical = int(old["N_selected_truth"]); recomputed = int(audit["P8C"][field])
            output.append({"sample_internal_name":"P8C", "association_method":"Lancestor",
                           "quantity":f"selected_{old['truth_species']}", "historical":historical,
                           "recomputed":recomputed, "exact":historical == recomputed,
                           "historical_source":str(args.p8c_outcomes_csv)})
    with args.output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader(); writer.writerows(output)
    if not all(row["exact"] for row in output):
        raise SystemExit("audit regression mismatch")
    print(f"PASS: {len(output)} exact frozen audit comparisons")


if __name__ == "__main__":
    main()
