#!/bin/bash
#SBATCH --job-name=H1K2_HCmaf1
#SBATCH --array=1-21
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=logs/H1K2_HCmaf1_%A_%a.out
#SBATCH --error=logs/H1K2_HCmaf1_%A_%a.err

CHR=${SLURM_ARRAY_TASK_ID}
RENAME_CHRS=/home/esalazarf/GitHub/alpha/Chi/misc/rename-chrs-rmv.txt
SITES=H1K2-v1.0.achr.HighConf_maf1.sites.vcf.gz
INPUT=../H1K2-v1.0.chr${CHR}.ARPnorm.vcf.gz
FILTERED=H1K2-v1.0.chr${CHR}.HCmaf1.vcf.gz
ANNOTATED=H1K2-v1.0.chr${CHR}.HCmaf1_varid.vcf.gz

module load bcftools

bcftools view \
    -T ${SITES} \
    ${INPUT} \
    -Oz -o ${FILTERED} \
    --threads 4 -W

bcftools annotate \
    --rename-chrs ${RENAME_CHRS} \
    --set-id '%CHROM-%POS-%REF-%ALT' \
    ${FILTERED} \
    -Oz -o ${ANNOTATED} \
    --threads 4 -W

mv ${ANNOTATED}     ${FILTERED}
mv ${ANNOTATED}.csi ${FILTERED}.csi
