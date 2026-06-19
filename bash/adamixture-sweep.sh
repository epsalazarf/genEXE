#!/bin/bash
#SBATCH --job-name=adamix_sweep
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --export=ALL
#SBATCH --output=logs/adamix_sweep_%j.out
#SBATCH --error=logs/adamix_sweep_%j.err

# Title: adamixture-sweep.sh
# Description: Unsupervised ADAMIXTURE multi-K sweep (K=2..10) with k-fold cross-
#              validation and a combined sweep plot. Runs on CPU by default.
#              GPU is opt-in (USE_GPU=1): ADAMIXTURE JIT-compiles a CUDA extension
#              on first use, and the NVIDIA HPC SDK stores libcudart outside the
#              path PyTorch guesses, so the GPU branch also exposes that lib dir.
#              Results are written to a subfolder next to the input data.
#              An optional toggle stages I/O on the FENIX scratch volume.
# Usage:  mkdir -p logs
#         # CPU (default):
#         sbatch bash/adamixture-sweep.sh <data_path> [name]
#         # GPU (opt-in): request a GPU and enable USE_GPU on the same line:
#         sbatch --partition=avx512 --gres=gpu:1 --export=ALL,USE_GPU=1 \
#                bash/adamixture-sweep.sh <data_path> [name]
#           <data_path>  Path to genotypes (.bed | .pgen | .vcf | .vcf.gz)
#           [name]       Optional run name (default: input basename)
# Developer: Pavel Salazar-Fernandez <pavel.salazar@galatea.bio>
# Dependencies: ADAMIXTURE (conda env), nvidiahpc/25.1-0 (CUDA, for GPU runs)
# Version: 0.2, 2026-06-18

set -euo pipefail

# ----------------------------------------------------------------------------
# Configuration (edit here)
# ----------------------------------------------------------------------------
ENV_NAME=ADAMIX          # conda env name holding the ADAMIXTURE install
ENV_PATH=${ENVDIR:+${ENVDIR%/}/${ENV_NAME}}          # built from $ENVDIR (export ENVDIR=/path/to/envs)
USE_GPU=${USE_GPU:-0}    # 1 = run on GPU (also needs --gres=gpu:1; see header usage)
CUDA_HOME=${CUDA_HOME:-} # optional: HPC SDK CUDA dir holding lib64/libcudart.so (auto-detected if empty)
MIN_K=2                  # sweep lower bound
MAX_K=10                 # sweep upper bound
CV_FOLDS=5               # cross-validation folds (5 = ADAMIXTURE default)
SEED=42                  # random seed (matches ADMIXTURE -s 42)
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
OUTDIR="${INPUT_DIR}/${NAME}_adamixture"
mkdir -p "${OUTDIR}"

echo "[!] ADAMIXTURE multi-K sweep"
echo " >  input    : ${DATA_PATH}"
echo " >  name     : ${NAME}"
echo " >  K range  : ${MIN_K}..${MAX_K}  (CV ${CV_FOLDS}-fold, seed ${SEED})"
echo " >  results  : ${OUTDIR}"

# ----------------------------------------------------------------------------
# Environment: load CUDA toolkit and activate the ADAMIXTURE env
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

# ----------------------------------------------------------------------------
# Device selection
#   CPU (default): multi-threaded, always works.
#   GPU (USE_GPU=1): load the CUDA toolkit and expose libcudart so the JIT
#   CUDA-extension link step finds -lcudart (the HPC SDK keeps it outside the
#   path PyTorch guesses). Falls back to CPU if CUDA still isn't usable.
# ----------------------------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-4}
DEVICE_ARGS=()

if [ "${USE_GPU}" -eq 1 ]; then
  module load nvidiahpc/25.1-0 2>/dev/null || echo "[!] WARNING: could not load nvidiahpc module." >&2

  # Locate the HPC SDK CUDA dir (the one with lib64/libcudart.so) unless given
  if [ -z "${CUDA_HOME}" ]; then
    for c in /opt/nvidia/hpc_sdk/Linux_x86_64/*/cuda /opt/nvidia/hpc_sdk/Linux_x86_64/*/cuda/*/; do
      if [ -e "${c%/}/lib64/libcudart.so" ]; then CUDA_HOME="${c%/}"; break; fi
    done
  fi
  if [ -n "${CUDA_HOME}" ] && [ -e "${CUDA_HOME}/lib64/libcudart.so" ]; then
    export CUDA_HOME
    export LIBRARY_PATH="${CUDA_HOME}/lib64:${LIBRARY_PATH:-}"      # link-time (-lcudart)
    export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH:-}" # run-time
    echo "[!] CUDA_HOME=${CUDA_HOME} (libcudart exposed for the JIT build)."
  else
    echo "[!] WARNING: libcudart.so not found under /opt/nvidia/hpc_sdk." >&2
    echo "    Set CUDA_HOME to the dir whose lib64/ holds libcudart.so. Try:" >&2
    echo "      find /opt/nvidia/hpc_sdk -name libcudart.so 2>/dev/null" >&2
  fi

  if python -c "import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)" 2>/dev/null; then
    DEVICE_ARGS=(--device gpu)
    echo "[!] Running on GPU."
    nvidia-smi || true
  else
    echo "[!] WARNING: USE_GPU=1 but CUDA not visible to torch — using CPU (${THREADS} threads)." >&2
  fi
else
  echo "[!] Running on CPU with ${THREADS} thread(s). (Set USE_GPU=1 + --gres=gpu:1 for GPU.)"
fi

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
# Run the sweep
# ----------------------------------------------------------------------------
echo "[!] Launching ADAMIXTURE..."
adamixture \
  --min_k "${MIN_K}" --max_k "${MAX_K}" \
  --cv "${CV_FOLDS}" \
  --data_path "${RUN_DATA}" \
  --save_dir "${RUN_DIR}" \
  --name "${NAME}" \
  --plot "${PLOT_FORMAT}" "${PLOT_DPI}" \
  -t "${THREADS}" \
  -s "${SEED}" \
  "${DEVICE_ARGS[@]}" \
  2>&1 | tee "${RUN_DIR}/${NAME}.${MIN_K}_${MAX_K}.log"

echo "[!] Completed >> ADAMIXTURE sweep (results in ${OUTDIR})"
