#!/bin/bash
#SBATCH --job-name=H1K2_HCmaf1_PCA
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=logs/P2b_HCmaf1_PCA_%j.out
#SBATCH --error=logs/P2b_HCmaf1_PCA_%j.err

# Title: P2b-concat-pca.sh
# Description: Concatenates LD-pruned PLINK2 files for chr 1-22, computes genome-wide PCA,
#              then cleans up intermediate pgen files (HWE and LD sets), keeping logs
#              and the direct VCF-to-pgen conversion files.
#              Depends on P2a-vcf-hwe-ld.sh completing successfully. Run via submit-P2.sh.
# Developer: Pavel Salazar-Fernandez <pavel.salazar@galatea.bio>
# Dependencies: PLINK2
# Version: 0.2, 2026-06-03

MANA=4
WORKDIR=${SLURM_SUBMIT_DIR}
PGEN_DIR=${WORKDIR}/pgen_files
MERGE_LIST=${PGEN_DIR}/H1K2-v1.0.HCmaf1.merge.txt

# Build chromosome file list
: > ${MERGE_LIST}
for i in {1..22}; do
  echo "${PGEN_DIR}/H1K2-v1.0.chr${i}.HCmaf1_LD" >> ${MERGE_LIST}
done
echo "[!] Merge list created: ${MERGE_LIST}"
echo " > Files listed: $(wc -l < ${MERGE_LIST})"

# Step 3: Concatenate chromosomes
echo "[!] STEP 3: Concatenating chromosomes..."
plink2 --pmerge-list ${MERGE_LIST} \
  --out ${PGEN_DIR}/H1K2-v1.0.achr.HCmaf1_LD \
  --threads ${MANA}

# Abort if merge failed
if [ $? -ne 0 ]; then
  echo "[ERROR] Merge failed — skipping PCA and cleanup." >&2
  exit 1
fi

# Step 4: Genome-wide PCA
echo "[!] STEP 4: Calculating genome-wide PCA..."
plink2 --pfile ${PGEN_DIR}/H1K2-v1.0.achr.HCmaf1_LD \
  --pca \
  --out ${WORKDIR}/H1K2-v1.0.achr.HCmaf1_LD.pca \
  --threads ${MANA}

# Step 5: Housekeeping — remove intermediate pgen sets, keep logs and VCF-derived pgens
echo "[!] STEP 5: Cleaning up intermediate pgen files..."
for i in {1..22}; do
  BASE=${PGEN_DIR}/H1K2-v1.0.chr${i}.HCmaf1

  # Remove HWE-filtered pgen set
  rm -f ${BASE}_HWE.pgen ${BASE}_HWE.pvar ${BASE}_HWE.psam
  echo "    [x] Removed: chr${i} HWE pgen set"

  # Remove LD-pruned pgen set (used for merge, no longer needed)
  rm -f ${BASE}_LD.pgen ${BASE}_LD.pvar ${BASE}_LD.psam
  echo "    [x] Removed: chr${i} LD pgen set"

  # Kept: ${BASE}.{pgen,pvar,psam}  (direct VCF-to-pgen)
  # Kept: ${BASE}_HWE.log, ${BASE}_LD.log, ${BASE}_LD.prune.{in,out}
done
echo "[!] STEP 5: Cleanup complete. VCF-derived pgens and all logs retained."

echo "[!] Completed >> PLINK2 CONCAT + PCA"
