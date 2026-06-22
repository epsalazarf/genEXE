#!/bin/bash
#SBATCH --job-name=adamix_cv_sweep
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=96:00:00
#SBATCH --export=ALL
#SBATCH --output=logs/adamix_cv_sweep_%j.out
#SBATCH --error=logs/adamix_cv_sweep_%j.err

# Title: adamixture-cv-sweep.sh
# Description: Sequential ADAMIXTURE multi-K sweep (K=2..10) WITH k-fold cross-
#              validation, as a single long job. Trades wall time for schedulability:
#              one 128G job queues far more easily than many concurrent high-mem
#              array tasks. Expect a long run (~2-4 days) — meant to be left over a
#              weekend. CV is memory-hungry (per-entry prediction over the whole
#              genotype matrix), hence 128G; training alone needs only ~8G.
#              Results go to a *separate* "<name>_adamixture_cv" folder so the
#              training-only sweep/array .P/.Q files are not overwritten. The CV
#              error per K is printed to the run log (grep hint at the end).
# Usage:  mkdir -p logs
#         sbatch bash/adamixture-cv-sweep.sh <data_path> [name]
#           <data_path>  Path to genotypes (.bed | .pgen | .vcf | .vcf.gz)
#           [name]       Optional run name (default: input basename)
#         Check the partition walltime cap first:  sinfo -o '%P %l'
# Developer: Pavel Salazar-Fernandez <pavel.salazar@galatea.bio>
# Dependencies: ADAMIXTURE (conda env)
# Version: 0.1, 2026-06-22

set -euo pipefail

# ----------------------------------------------------------------------------
# Configuration (edit here)
# ----------------------------------------------------------------------------
ENV_NAME=ADAMIX          # conda env name holding the ADAMIXTURE install
ENV_PATH=${ENVDIR:+${ENVDIR%/}/${ENV_NAME}}          # built from $ENVDIR (export ENVDIR=/path/to/envs)
MIN_K=2                  # sweep lower bound
MAX_K=10                 # sweep upper bound
CV_FOLDS=5               # cross-validation folds (5 = ADAMIXTURE default)
SEED=42                  # random seed (matches the training runs)
PLOT_FORMAT=png          # combined sweep plot format (png | pdf | jpg)
PLOT_DPI=300             # combined sweep plot resolution
USE_SCRATCH=0            # 1 = stage input + run on FENIX scratch, copy results back
SCRATCH_ROOT=${SCRATCH:-}  # FENIX scratch (from $SCRATCH; empty -> scratch is skipped)

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
OUTDIR="${INPUT_DIR}/${NAME}_adamixture_cv"
mkdir -p "${OUTDIR}"

echo "[!] ADAMIXTURE CV sweep (sequential)"
echo " >  input    : ${DATA_PATH}"
echo " >  name     : ${NAME}"
echo " >  K range  : ${MIN_K}..${MAX_K}  (CV ${CV_FOLDS}-fold, seed ${SEED})"
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
  # $SCRATCH is already user-specific; per-job subdir keeps concurrent runs isolated
  SCRATCH_JOB="${SCRATCH_ROOT%/}/job_${SLURM_JOB_ID:-$$}"
  echo "[!] Scratch staging enabled: ${SCRATCH_JOB}"
  mkdir -p "${SCRATCH_JOB}/tmp" "${SCRATCH_JOB}/out"
  export TMPDIR="${SCRATCH_JOB}/tmp"          # torch/numpy temp spill -> fast local disk
  for f in "${COMPANIONS[@]}"; do
    [ -f "${f}" ] && cp -v "${f}" "${SCRATCH_JOB}/"
  done
  RUN_DATA="${SCRATCH_JOB}/${INPUT_FILE}"
  RUN_DIR="${SCRATCH_JOB}/out"
  # Copy results back to the data-adjacent output folder on exit (success or fail)
  trap 'echo "[!] Copying results from scratch back to ${OUTDIR}"; \
        cp -a "${RUN_DIR}/." "${OUTDIR}/" 2>/dev/null || true; \
        rm -rf "${SCRATCH_JOB}"' EXIT
fi

# ----------------------------------------------------------------------------
# Run the sweep with cross-validation (sequential over all K)
# ----------------------------------------------------------------------------
echo "[!] Launching ADAMIXTURE CV sweep (this is a long run)..."
adamixture \
  --min_k "${MIN_K}" --max_k "${MAX_K}" \
  --cv "${CV_FOLDS}" \
  --data_path "${RUN_DATA}" \
  --save_dir "${RUN_DIR}" \
  --name "${NAME}" \
  --plot "${PLOT_FORMAT}" "${PLOT_DPI}" \
  -t "${THREADS}" \
  -s "${SEED}" \
  2>&1 | tee "${RUN_DIR}/${NAME}.${MIN_K}_${MAX_K}.cv.log"

echo "[!] Completed >> ADAMIXTURE CV sweep (results in ${OUTDIR})"
echo "[!] Inspect CV errors per K with:"
echo "      grep -iE 'cv|cross' ${OUTDIR}/${NAME}.${MIN_K}_${MAX_K}.cv.log"
