#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "usage: $0 --input-manifest FILE --output-root DIR [--workers N] [--prefix TAG] [--dedup-mode reco|gen] [--assoc-max-dr X] [--dry-run]" >&2
}
fail() { echo "ERROR: $*" >&2; exit 2; }

manifest=""; output_root=""; workers=40; prefix=""; dedup_mode=reco; assoc_max_dr=0.1; dry_run=0
while (( $# )); do
    case "$1" in
        --input-manifest) (( $# >= 2 )) || fail "--input-manifest requires a value"; manifest=$2; shift 2 ;;
        --output-root) (( $# >= 2 )) || fail "--output-root requires a value"; output_root=$2; shift 2 ;;
        --workers) (( $# >= 2 )) || fail "--workers requires a value"; workers=$2; shift 2 ;;
        --prefix) (( $# >= 2 )) || fail "--prefix requires a value"; prefix=$2; shift 2 ;;
        --dedup-mode) (( $# >= 2 )) || fail "--dedup-mode requires a value"; dedup_mode=$2; shift 2 ;;
        --assoc-max-dr) (( $# >= 2 )) || fail "--assoc-max-dr requires a value"; assoc_max_dr=$2; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) fail "unknown argument: $1" ;;
    esac
done
[[ -n "$manifest" && -n "$output_root" ]] || { usage; fail "manifest and output root are required"; }
[[ "$workers" =~ ^[1-9][0-9]*$ ]] || fail "workers must be positive"
[[ "$dedup_mode" == reco || "$dedup_mode" == gen ]] || fail "invalid dedup mode"
[[ -z "$prefix" || "$prefix" =~ ^[A-Za-z0-9_.-]+$ ]] || fail "unsafe prefix"

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
manifest=$(readlink -f "$manifest"); output_root=$(readlink -m "$output_root")
test -r "$manifest" || fail "manifest is not readable: $manifest"
mapfile -t input_files < <(python "$repo_root/scripts/analysis/manifest_rec_paths.py" "$manifest")
sample=$(python "$repo_root/scripts/analysis/manifest_rec_paths.py" "$manifest" --sample-only)
(( ${#input_files[@]} > 0 )) || fail "manifest contains no inputs"
for input_file in "${input_files[@]}"; do test -s "$input_file" || fail "missing or empty REC: $input_file"; done
prefix=${prefix:-${sample}_}

command=(python "$repo_root/HitAnalysis/particle_level_analisis_parallel.py"
    --input-list "${input_files[@]}" --output-root "$output_root" --prefix "$prefix"
    --n-workers "$workers" --dedup-mode "$dedup_mode" --assoc-max-dr "$assoc_max_dr"
    -v --min-energy-cuts 10 --all-plot 11 13 22 211)
printf 'Command:'; printf ' %q' "${command[@]}"; printf '\n'
printf 'Sample: %s\nInputs: %d\nOutput root: %s\n' "$sample" "${#input_files[@]}" "$output_root"
(( dry_run == 0 )) || exit 0
test ! -e "$output_root/${prefix}results0.4_tph0.0_tpi0.0_n0.0_g0.0" || fail "output exists; refusing overwrite"
mkdir -p "$output_root"
(cd "$repo_root" && "${command[@]}")
python "$repo_root/scripts/analysis/write_run_metadata.py" \
    --repo-root "$repo_root" --output "$output_root/${prefix}results0.4_tph0.0_tpi0.0_n0.0_g0.0/run_metadata.json" \
    --sample "$sample" --source-manifest "$manifest" --input-count "${#input_files[@]}" \
    --workers "$workers" --dedup-mode "$dedup_mode" --assoc-max-dr "$assoc_max_dr"
