#!/bin/bash
#SBATCH --job-name=H1K2_HCmaf1_LD
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=logs/P2a_HCmaf1_LD_%A_%a.out
#SBATCH --error=logs/P2a_HCmaf1_LD_%A_%a.err

# Title: P2a-vcf-hwe-ld.sh
# Description: Converts HCmaf1 VCF files (chr 1-22) to PLINK2 format, applies HWE filter,
#              and performs LD pruning. Parallelized per chromosome via SLURM array.
#              Run via submit-P2.sh from the data root directory (./).
# Developer: Pavel Salazar-Fernandez <pavel.salazar@galatea.bio>
# Dependencies: PLINK2
# Version: 0.2, 2026-06-03

CHR=${SLURM_ARRAY_TASK_ID}
MANA=4
WORKDIR=${SLURM_SUBMIT_DIR}
INPUT=${WORKDIR}/H1K2-v1.0.chr${CHR}.HCmaf1.vcf.gz
PREFIX=${WORKDIR}/pgen_files/H1K2-v1.0.chr${CHR}.HCmaf1

mkdir -p ${WORKDIR}/pgen_files

echo "[!] CHR${CHR}: Starting VCF conversion, HWE filtering, and LD pruning..."

# Step 1: VCF to PLINK2 format
echo "[&] CHR${CHR}: Converting VCF to PLINK2..."
plink2 --vcf ${INPUT} \
  --make-pgen --vcf-half-call m \
  --out ${PREFIX} \
  --threads ${MANA}

# Step 2: HWE filtering
echo "[&] CHR${CHR}: Applying HWE filter..."
plink2 --pfile ${PREFIX} \
  --hardy \
  --hwe 1e-25 keep-fewhet \
  --make-pgen \
  --out ${PREFIX}_HWE \
  --threads ${MANA}

# Step 3a: LD pruning scan (on HWE-filtered set)
echo "[&] CHR${CHR}: Scanning for LD..."
plink2 --pfile ${PREFIX}_HWE \
  --indep-pairwise 200kb 1 0.5 \
  --out ${PREFIX}_LD \
  --threads ${MANA}

# Step 3b: Extract LD-pruned SNPs
echo "[&] CHR${CHR}: Extracting LD-pruned SNPs..."
plink2 --pfile ${PREFIX}_HWE \
  --extract ${PREFIX}_LD.prune.in \
  --make-pgen \
  --out ${PREFIX}_LD \
  --threads ${MANA}

echo "[!] CHR${CHR}: Done."
