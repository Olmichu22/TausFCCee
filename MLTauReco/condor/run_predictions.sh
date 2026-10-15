#!/bin/bash
# Predicciones ParTauDETR sobre una lista de ficheros EDM4hep, repartida entre
# NPROC procesos que comparten la GPU del job. Relanzable: los ficheros con
# salida ya escrita se saltan.
#
#   run_predictions.sh <lista de ficheros> <directorio de salida> [NPROC] [THREADS]
set -uo pipefail

LIST="$1"
OUT="$2"
NPROC="${3:-16}"
THREADS="${4:-2}"

# Ruta absoluta: condor copia el ejecutable al spool, así que no se puede
# deducir el repo a partir de la ubicación de este script.
REPO="${MLTAU_REPO:-/nfs/cms/arqolmo/TausFCCee}"
LOGS="$REPO/MLTauReco/condor/logs"
mkdir -p "$LOGS"
# stdout/stderr del job se quedan en el spool de condor: se duplican aquí
exec > >(tee -a "$LOGS/job_$(date +%Y%m%d_%H%M%S).log") 2>&1

export MLTAU_CONTAINER="${MLTAU_CONTAINER:-/nfs/cms/arqolmo/containers/mltau.sif}"
# /nfs: repo, modelo y contenedor; /pnfs: entradas y salida
export MLTAU_BINDS="/nfs,/pnfs"
export MLTAU_APPTAINER_ARGS="--nv"

echo "host $(hostname), CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv 2>&1 || true
echo "$(wc -l < "$LIST") files -> $OUT, $NPROC processes x $THREADS threads"

mkdir -p "$OUT" "$LOGS"
# dCache no deja sobrescribir: solo se copia la primera vez
[[ -e "$OUT/input_files.txt" ]] || cp "$LIST" "$OUT/input_files.txt"
cd "$REPO" || exit 1

WORK="$(mktemp -d "${_CONDOR_SCRATCH_DIR:-/tmp}/mltau_parts.XXXXXX")"
split -d -n "r/$NPROC" "$LIST" "$WORK/part_"

pids=()
for part in "$WORK"/part_*; do
    # shellcheck disable=SC2046  # una ruta por línea, sin espacios
    MLTauReco/run_mltau.sh python3 MLTauReco/process_edm4hep.py $(cat "$part") \
        --out-dir "$OUT" --threads "$THREADS" > "$LOGS/worker_$(basename "$part").log" 2>&1 &
    pids+=($!)
done

status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done
rm -rf "$WORK"
echo "done: $(ls "$OUT" | grep -c '_mltaus.parquet$') parquet files in $OUT, status $status"
exit $status
