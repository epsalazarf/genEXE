#!/bin/bash
#SBATCH --job-name=adamix_cv
#SBATCH --array=2-10
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --export=ALL
#SBATCH --output=logs/adamix_cv_%A_%a.out
#SBATCH --error=logs/adamix_cv_%A_%a.err

# Title: adamixture-cv-array.sh
# Description: Cross-validation for ADAMIXTURE, parallelized as a SLURM array —
#              one K per array task (K = array index). Each task runs a single-K
#              k-fold CV (`-k K --cv N`) so the whole K range finishes in roughly
#              one task's wall time instead of summing sequentially. CV retrains
#              every fold, so it is the slow part: keeping it off the sweep job
#              and fanning it out per K is the point of this script.
#              Results go to a *separate* "<name>_adamixture_cv" folder so the
#              sweep's .P/.Q files are not overwritten. Plotting is intentionally
#              omitted (do it afterwards with adamixture-plot on the .Q files).
# Usage:  mkdir -p logs
#         sbatch bash/adamixture-cv-array.sh <data_path> [name]
#           Edit --array above to change the K range (must start at >=2).
#           <data_path>  Path to genotypes (.bed | .pgen | .vcf | .vcf.gz)
#           [name]       Optional run name (default: input basename)
# Developer: Pavel Salazar-Fernandez <pavel.salazar@galatea.bio>
# Dependencies: ADAMIXTURE (conda env)
# Version: 0.1, 2026-06-22

set -euo pipefail

# ----------------------------------------------------------------------------
# Configuration (edit here)
# ----------------------------------------------------------------------------
ENV_NAME=ADAMIX          # conda env name holding the ADAMIXTURE install
ENV_PATH=${ENVDIR:+${ENVDIR%/}/${ENV_NAME}}          # built from $ENVDIR (export ENVDIR=/path/to/envs)
CV_FOLDS=5               # cross-validation folds (5 = ADAMIXTURE default)
SEED=42                  # random seed (matches the sweep)
OUT_SUFFIX=_adamixture_cv  # results subfolder (kept separate from the sweep)
USE_SCRATCH=0            # 1 = stage input + run on FENIX scratch, copy results back
SCRATCH_ROOT=${SCRATCH:-}  # FENIX scratch (from $SCRATCH; empty -> scratch is skipped)

K=${SLURM_ARRAY_TASK_ID:?This script must be submitted as a SLURM array (sbatch --array=...)}

# ----------------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------------
if [ $# -lt 1 ]; then
  echo "[ERROR] Missing input. Usage: sbatch $0 <data_path> [name]" >&2
  exit 1
fi

DATA_PATH=$(realpath "$1")
if [ ! -f "${DATA_PATH}" ]; then
  echo "[ERROR] Input file not found: ${DATA_PATH}" >&2
  exit 1
fi

INPUT_DIR=$(dirname "${DATA_PATH}")
INPUT_FILE=$(basename "${DATA_PATH}")

# Strip the genotype extension to obtain the dataset basename / prefix
case "${INPUT_FILE}" in
  *.vcf.gz) STEM="${INPUT_FILE%.vcf.gz}";   COMPANIONS=("${DATA_PATH}" "${DATA_PATH}.tbi" "${DATA_PATH}.csi") ;;
  *.vcf)    STEM="${INPUT_FILE%.vcf}";      COMPANIONS=("${DATA_PATH}") ;;
  *.bed)    STEM="${INPUT_FILE%.bed}";      COMPANIONS=("${INPUT_DIR}/${STEM}.bed" "${INPUT_DIR}/${STEM}.bim" "${INPUT_DIR}/${STEM}.fam") ;;
  *.pgen)   STEM="${INPUT_FILE%.pgen}";     COMPANIONS=("${INPUT_DIR}/${STEM}.pgen" "${INPUT_DIR}/${STEM}.pvar" "${INPUT_DIR}/${STEM}.psam") ;;
  *) echo "[ERROR] Unsupported input '${INPUT_FILE}'. Use .bed, .pgen, .vcf or .vcf.gz" >&2; exit 1 ;;
esac

NAME=${2:-${STEM}}
OUTDIR="${INPUT_DIR}/${NAME}${OUT_SUFFIX}"
mkdir -p "${OUTDIR}"

echo "[!] ADAMIXTURE cross-validation (array task K=${K})"
echo " >  input    : ${DATA_PATH}"
echo " >  name     : ${NAME}"
echo " >  K        : ${K}  (CV ${CV_FOLDS}-fold, seed ${SEED})"
echo " >  results  : ${OUTDIR}"

# ----------------------------------------------------------------------------
# Environment: activate the ADAMIXTURE env (CPU only — GPU is blocked on FENIX)
# ----------------------------------------------------------------------------
if [ -z "${ENV_PATH}" ]; then
  echo "[ERROR] \$ENVDIR is not set. Export it before submitting, e.g.:" >&2
  echo "          export ENVDIR=/mnt/data/fsanchezq/esalazarf/envs" >&2
  echo "        (the script then uses \$ENVDIR/${ENV_NAME})" >&2
  exit 1
fi

export PATH="${ENV_PATH}/bin:${PATH}"
if ! command -v adamixture >/dev/null 2>&1; then
  echo "[ERROR] 'adamixture' not found in ${ENV_PATH}/bin" >&2
  exit 1
fi

THREADS=${SLURM_CPUS_PER_TASK:-8}
echo "[!] Running on CPU with ${THREADS} thread(s)."

# ----------------------------------------------------------------------------
# Optional: stage input + run on scratch, then copy results back
# ----------------------------------------------------------------------------
RUN_DIR="${OUTDIR}"
RUN_DATA="${DATA_PATH}"
if [ "${USE_SCRATCH}" -eq 1 ] && [ -z "${SCRATCH_ROOT}" ]; then
  echo "[!] USE_SCRATCH=1 but \$SCRATCH is unset/empty — running in place (no scratch)."
fi

if [ "${USE_SCRATCH}" -eq 1 ] && [ -n "${SCRATCH_ROOT}" ]; then
  # $SCRATCH is already user-specific; per-task subdir keeps array tasks isolated
  SCRATCH_JOB="${SCRATCH_ROOT%/}/job_${SLURM_JOB_ID:-$$}"
  echo "[!] Scratch staging enabled: ${SCRATCH_JOB}"
  mkdir -p "${SCRATCH_JOB}/tmp" "${SCRATCH_JOB}/out"
  export TMPDIR="${SCRATCH_JOB}/tmp"          # torch/numpy temp spill -> fast local disk
  for f in "${COMPANIONS[@]}"; do
    [ -f "${f}" ] && cp -v "${f}" "${SCRATCH_JOB}/"
  done
  RUN_DATA="${SCRATCH_JOB}/${INPUT_FILE}"
  RUN_DIR="${SCRATCH_JOB}/out"
  # Copy this K's results back to the shared CV folder on exit (success or fail)
  trap 'echo "[!] Copying K='"${K}"' results from scratch back to ${OUTDIR}"; \
        cp -a "${RUN_DIR}/." "${OUTDIR}/" 2>/dev/null || true; \
        rm -rf "${SCRATCH_JOB}"' EXIT
fi

# ----------------------------------------------------------------------------
# Run single-K cross-validation
#   Plotting is omitted on purpose; ADAMIXTURE may still emit a default per-K
#   png, which is harmless — do the real (labeled) plot later via adamixture-plot.
# ----------------------------------------------------------------------------
echo "[!] Launching ADAMIXTURE CV for K=${K}..."
adamixture \
  -k "${K}" \
  --cv "${CV_FOLDS}" \
  --data_path "${RUN_DATA}" \
  --save_dir "${RUN_DIR}" \
  --name "${NAME}" \
  -t "${THREADS}" \
  -s "${SEED}" \
  2>&1 | tee "${RUN_DIR}/${NAME}.cv.${K}.log"

echo "[!] Completed >> ADAMIXTURE CV for K=${K} (results in ${OUTDIR})"
echo "[!] Collect CV errors across K once all tasks finish, e.g.:"
echo "      grep -iH 'cv error' ${OUTDIR}/${NAME}.cv.*.log"
