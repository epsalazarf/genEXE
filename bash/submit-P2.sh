#!/bin/bash

# Title: submit-P2.sh
# Description: Submits the P2 pipeline from the data root directory (./).
#              P2a runs as a SLURM array (one job per chromosome, 1-22) for
#              VCF conversion, HWE filtering, and LD pruning. P2b is held until
#              all P2a tasks succeed, then concatenates chromosomes and runs PCA.
# Usage: bash scripts/submit-P2.sh   (run from the data root directory)
# Developer: Pavel Salazar-Fernandez <pavel.salazar@galatea.bio>
# Version: 0.2, 2026-06-03

SCRIPT_DIR=$(dirname "$(realpath "$0")")

mkdir -p logs pgen_files

JOB_A=$(sbatch --parsable ${SCRIPT_DIR}/P2a-vcf-hwe-ld.sh)
echo "[!] Submitted P2a (VCF + HWE + LD) array job: ${JOB_A}"

JOB_B=$(sbatch --parsable --dependency=afterok:${JOB_A} ${SCRIPT_DIR}/P2b-concat-pca.sh)
echo "[!] Submitted P2b (concat + PCA + cleanup) job: ${JOB_B} (depends on ${JOB_A})"

echo ""
echo "[!] Monitor with: squeue -j ${JOB_A},${JOB_B}"
echo "[!] Watch logs:   tail -f logs/P2a_HCmaf1_LD_${JOB_A}_*.out"
