#!/bin/bash

# Title: submit-adamixture-cv.sh
# Description: Submits the ADAMIXTURE cross-validation pipeline: a SLURM array
#              (one K per task, adamixture-cv-array.sh) followed by an aggregation
#              job (adamixture-cv-collect.sh) held with afterany so it summarizes
#              whatever finished. Arguments are passed through to both scripts.
# Usage: bash bash/submit-adamixture-cv.sh <data_path> [name]
#          Edit the K range via --array in adamixture-cv-array.sh (default 2-10).
# Developer: Pavel Salazar-Fernandez <pavel.salazar@galatea.bio>
# Version: 0.1, 2026-06-22

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "[ERROR] Usage: bash $0 <data_path> [name]" >&2
  exit 1
fi

SCRIPT_DIR=$(dirname "$(realpath "$0")")
mkdir -p logs

ARRAY_JOB=$(sbatch --parsable "${SCRIPT_DIR}/adamixture-cv-array.sh" "$@")
echo "[!] Submitted CV array job (K per task): ${ARRAY_JOB}"

COLLECT_JOB=$(sbatch --parsable --dependency=afterany:${ARRAY_JOB} \
  "${SCRIPT_DIR}/adamixture-cv-collect.sh" "$@")
echo "[!] Submitted CV collect job: ${COLLECT_JOB} (afterany ${ARRAY_JOB})"

echo ""
echo "[!] Monitor with: squeue -j ${ARRAY_JOB},${COLLECT_JOB}"
echo "[!] CV table will be written next to the input as <name>_adamixture_cv/<name>.cv_errors.tsv"
