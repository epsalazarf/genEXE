#!/bin/bash
#SBATCH --job-name=adamix_train
#SBATCH --array=2-10
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --export=ALL
#SBATCH --output=logs/adamix_train_%A_%a.out
#SBATCH --error=logs/adamix_train_%A_%a.err

# Title: adamixture-train-array.sh
# Description: ADAMIXTURE multi-K training (the ancestry inference), parallelized
#              as a SLURM array — one K per task (K = array index). All K train at
#              once, so wall time is ~the slowest single K instead of the sum
#              (~minutes-to-hours per K vs ~10h sequential in adamixture-sweep.sh).
#              No cross-validation (training peaks ~8G, fits 16G comfortably).
#              ADAMIXTURE auto-emits a per-K .png preview, which is kept on purpose.
#              The aligned multi-K combined plot is done afterwards with
#              adamixture-plot on the .Q files. Results go to "<name>_adamixture/".
# Usage:  mkdir -p logs
#         sbatch bash/adamixture-train-array.sh <data_path> [name]
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
SEED=42                  # random seed (matches the sweep / ADMIXTURE -s 42)
OUT_SUFFIX=_adamixture   # results subfolder (same as the sweep — this is the inference)
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

echo "[!] ADAMIXTURE training (array task K=${K})"
echo " >  input    : ${DATA_PATH}"
echo " >  name     : ${NAME}"
echo " >  K        : ${K}  (seed ${SEED}, no CV)"
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
  # Copy this K's results back to the shared output folder on exit (success or fail)
  trap 'echo "[!] Copying K='"${K}"' results from scratch back to ${OUTDIR}"; \
        cp -a "${RUN_DIR}/." "${OUTDIR}/" 2>/dev/null || true; \
        rm -rf "${SCRATCH_JOB}"' EXIT
fi

# ----------------------------------------------------------------------------
# Train a single K. ADAMIXTURE auto-generates a per-K .png preview (kept).
# The aligned multi-K plot is produced later via adamixture-plot on the .Q files.
# ----------------------------------------------------------------------------
echo "[!] Launching ADAMIXTURE training for K=${K}..."
adamixture \
  -k "${K}" \
  --data_path "${RUN_DATA}" \
  --save_dir "${RUN_DIR}" \
  --name "${NAME}" \
  -t "${THREADS}" \
  -s "${SEED}" \
  2>&1 | tee "${RUN_DIR}/${NAME}.${K}.log"

echo "[!] Completed >> ADAMIXTURE training for K=${K} (results in ${OUTDIR})"
