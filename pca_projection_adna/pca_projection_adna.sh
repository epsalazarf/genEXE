#!/usr/bin/env bash
# =============================================================================
# pca_projection_adna.sh
# PCA projection of ancient DNA samples onto a modern reference panel.
#
# METHODOLOGY
#   1. QC and LD-prune the reference panel (skipped if outputs exist).
#   2. Compute PC allele weights on the reference panel exclusively
#      (--pca allele-wts), so ancient samples never influence the PC axes.
#   3. Intersect variants between reference and ancient dataset.
#   4. Project ancient samples onto reference PCs via --score with
#      no-mean-imputation (critical: missing genotypes must not be filled
#      from reference frequencies in aDNA contexts).
#   5. (Optional) Produce a quick diagnostic plot in R.
#
# REFERENCE
#   PLINK2 projection:  https://www.cog-genomics.org/plink/2.0/score#pca_project
#   Methodological basis: Barrie et al. 2024, Nature; Patterson et al. 2006.
#
# DEPENDENCIES
#   - plink2 >= 2.0.0-a.6 (allele-wts + acount syntax)
#   - R >= 4.0 with: ggplot2, data.table (optional plot)
#   - Standard coreutils (awk, sort, comm)
#
# USAGE
#   bash pca_projection_adna.sh [OPTIONS]
#
# OPTIONS (all have defaults; edit CONFIG section or pass as env vars)
#   REF_PFILE          Prefix of reference panel PLINK2 pfile (no extension)
#   ANCIENT_PFILE      Prefix of ancient dataset PLINK2 pfile (no extension)
#   OUTDIR             Output directory (will be created if absent)
#   THREADS            CPU threads for PLINK2
#   N_PCS              Number of principal components to compute/project
#   REF_GENO           Max per-variant missingness in reference panel (0-1)
#   REF_MIND           Max per-sample missingness in reference panel (0-1)
#   ANC_GENO           Max per-variant missingness in ancient dataset (0-1)
#                      (set to 1.0 to disable; ancient samples are sparse)
#   ANC_MIND           Max per-sample missingness in ancient dataset (0-1)
#                      (set to 1.0 to disable)
#   MAF_REF            Minor allele frequency filter for reference panel
#   LD_WIN             LD pruning window size (kb)
#   LD_STEP            LD pruning step size (variants)
#   LD_R2              LD pruning r² threshold
#   EXCLUDE_AMBIGUOUS  "yes" to drop A/T and C/G SNPs; "no" to retain them
#   KING_CUTOFF        Kinship coefficient cutoff for relatedness (reference)
#   MAKE_PLOT          "yes" to run R plotting after projection
#   REF_FA             Path to reference FASTA (optional; for --ref-from-fa)
#
# OUTPUT FILES (under OUTDIR/)
#   ref_qc.*           QC'd, LD-pruned reference panel (pgen/pvar/psam)
#   ref_pca.acount     Allele counts used for PCA (required for projection)
#   ref_pca.eigenvec            Sample-level PC coordinates (reference)
#   ref_pca.eigenvec.allele     Allele weights for projection
#   ref_pca.eigenval            Eigenvalues
#   intersect_variants.txt      Variant IDs shared by ref and ancient datasets
#   ancient_proj.sscore         Projected PC scores for ancient samples
#   pca_projection.pdf          Diagnostic plot (if MAKE_PLOT=yes)
#
# NOTES ON COLUMN INDICES (--score)
#   The .eigenvec.allele file produced by recent PLINK2 versions has:
#     col 1: #CHROM  col 2: ID  col 3: REF  col 4: ALT
#     col 5: A1 (the counted allele)  col 6+: PC weights (PC1..PCN)
#   Therefore: --score ref_pca.eigenvec.allele 2 5 ... --score-col-nums 6-$((5+N_PCS))
#   VERIFY this against your actual .eigenvec.allele header before running.
#
# =============================================================================
set -euo pipefail

# =============================================================================
# CONFIG — override via environment or edit here
# =============================================================================
REF_PFILE="${REF_PFILE:-/path/to/modern_reference}"
ANCIENT_PFILE="${ANCIENT_PFILE:-/path/to/ancient_dataset}"
OUTDIR="${OUTDIR:-./pca_projection_out}"
THREADS="${THREADS:-4}"
N_PCS="${N_PCS:-10}"

# Missingness — separate thresholds for reference vs ancient
REF_GENO="${REF_GENO:-0.05}"   # exclude variants >5% missing in reference
REF_MIND="${REF_MIND:-0.10}"   # exclude reference samples >10% missing
ANC_GENO="${ANC_GENO:-1.0}"    # no variant missingness filter on ancients
ANC_MIND="${ANC_MIND:-1.0}"    # no sample missingness filter on ancients

# MAF and LD pruning (reference panel only)
MAF_REF="${MAF_REF:-0.05}"
LD_WIN="${LD_WIN:-1000}"
LD_STEP="${LD_STEP:-50}"
LD_R2="${LD_R2:-0.1}"

# Strand-ambiguous SNP handling: "yes" = exclude A/T and C/G SNPs
EXCLUDE_AMBIGUOUS="${EXCLUDE_AMBIGUOUS:-yes}"

# Kinship cutoff for reference panel relatedness pruning
KING_CUTOFF="${KING_CUTOFF:-0.0884}"   # ~2nd degree (KING default)

# Optional quick R plot
MAKE_PLOT="${MAKE_PLOT:-yes}"

# Optional reference FASTA (leave empty to skip --ref-from-fa)
REF_FA="${REF_FA:-}"

# =============================================================================
# HELPERS
# =============================================================================
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }

check_pfile() {
    local prefix="$1"
    [[ -f "${prefix}.pgen" && -f "${prefix}.pvar" && -f "${prefix}.psam" ]]
}

check_dep() {
    command -v "$1" &>/dev/null || { log "ERROR: '$1' not found in PATH."; exit 1; }
}

# =============================================================================
# PREFLIGHT
# =============================================================================
check_dep plink2
check_dep awk
check_dep sort
check_dep comm
[[ "${MAKE_PLOT}" == "yes" ]] && check_dep Rscript

check_pfile "${REF_PFILE}"     || { log "ERROR: Reference pfile not found: ${REF_PFILE}.{pgen,pvar,psam}"; exit 1; }
check_pfile "${ANCIENT_PFILE}" || { log "ERROR: Ancient pfile not found: ${ANCIENT_PFILE}.{pgen,pvar,psam}"; exit 1; }

mkdir -p "${OUTDIR}"

# Build optional --ref-from-fa flag
FA_FLAG=""
if [[ -n "${REF_FA}" ]]; then
    [[ -f "${REF_FA}" ]] || { log "ERROR: REF_FA set but file not found: ${REF_FA}"; exit 1; }
    FA_FLAG="--ref-from-fa --fa ${REF_FA}"
    log "Reference FASTA provided; will use --ref-from-fa."
else
    log "No REF_FA set; skipping --ref-from-fa. REF/ALT orientation will follow input encoding."
fi

# Column indices for --score (verify against your .eigenvec.allele header)
SCORE_ID_COL=2
SCORE_A1_COL=5
SCORE_START_COL=6
SCORE_END_COL=$(( SCORE_START_COL + N_PCS - 1 ))

log "Configuration summary:"
log "  REF_PFILE        = ${REF_PFILE}"
log "  ANCIENT_PFILE    = ${ANCIENT_PFILE}"
log "  OUTDIR           = ${OUTDIR}"
log "  N_PCS            = ${N_PCS}"
log "  EXCLUDE_AMBIG    = ${EXCLUDE_AMBIGUOUS}"
log "  REF missingness  = geno ${REF_GENO} / mind ${REF_MIND}"
log "  ANC missingness  = geno ${ANC_GENO} / mind ${ANC_MIND}"
log "  Score cols       = ID:${SCORE_ID_COL} A1:${SCORE_A1_COL} PCs:${SCORE_START_COL}-${SCORE_END_COL}"

# =============================================================================
# STEP 1 — REFERENCE PANEL QC AND LD PRUNING
# (skipped if ref_qc.pgen already exists)
# =============================================================================
REF_QC="${OUTDIR}/ref_qc"
REF_KING="${OUTDIR}/ref_king"

if check_pfile "${REF_QC}"; then
    log "STEP 1 SKIPPED: QC'd reference pfile already exists (${REF_QC})."
else
    log "STEP 1: Reference panel QC and LD pruning..."

    # 1a. Basic QC: MAF, missingness
    log "  1a. Applying MAF and missingness filters..."
    plink2 \
        --pfile "${REF_PFILE}" \
        --human \
        --maf "${MAF_REF}" \
        --geno "${REF_GENO}" \
        --mind "${REF_MIND}" \
        --autosome \
        --snps-only just-acgt \
        ${FA_FLAG} \
        --make-pgen \
        --threads "${THREADS}" \
        --out "${OUTDIR}/ref_mafgeno"

    # 1b. Remove strand-ambiguous SNPs (A/T and C/G) if requested
    if [[ "${EXCLUDE_AMBIGUOUS}" == "yes" ]]; then
        log "  1b. Removing strand-ambiguous SNPs (A/T, C/G)..."
        # Extract variant IDs where REF/ALT are complementary ambiguous pairs.
        # Expects ID format CHR-POS-REF-ALT (col 3); parses REF/ALT from the ID
        # rather than cols 4/5 to be robust to pvar formatting differences.
        awk 'NR>1 {
            id=$3;
            n=split(id, parts, "-");
            if (n >= 4) {
                ref=parts[n-1]; alt=parts[n];
                if ((ref=="A" && alt=="T") || (ref=="T" && alt=="A") ||
                    (ref=="C" && alt=="G") || (ref=="G" && alt=="C"))
                    print id
            }
        }' "${OUTDIR}/ref_mafgeno.pvar" > "${OUTDIR}/ambiguous_snps.txt"
        N_AMB=$(wc -l < "${OUTDIR}/ambiguous_snps.txt")
        log "  Found ${N_AMB} strand-ambiguous variants to exclude."
        plink2 \
            --pfile "${OUTDIR}/ref_mafgeno" \
            --exclude "${OUTDIR}/ambiguous_snps.txt" \
            --make-pgen \
            --threads "${THREADS}" \
            --out "${OUTDIR}/ref_noamb"
        REF_STEP1="${OUTDIR}/ref_noamb"
    else
        log "  1b. Retaining strand-ambiguous SNPs (EXCLUDE_AMBIGUOUS=no)."
        REF_STEP1="${OUTDIR}/ref_mafgeno"
    fi

    # 1c. LD pruning
    log "  1c. LD pruning (window=${LD_WIN}kb, step=${LD_STEP}, r2=${LD_R2})..."
    plink2 \
        --pfile "${REF_STEP1}" \
        --indep-pairwise "${LD_WIN}kb" "${LD_STEP}" "${LD_R2}" \
        --threads "${THREADS}" \
        --out "${OUTDIR}/ref_ld"

    # 1d. Relatedness pruning (reference panel only)
    log "  1d. Relatedness pruning (KING cutoff=${KING_CUTOFF})..."
    plink2 \
        --pfile "${REF_STEP1}" \
        --extract "${OUTDIR}/ref_ld.prune.in" \
        --king-cutoff "${KING_CUTOFF}" \
        --threads "${THREADS}" \
        --out "${REF_KING}"

    # 1e. Build final QC'd reference pfile
    log "  1e. Building final QC'd reference pfile..."
    plink2 \
        --pfile "${REF_STEP1}" \
        --extract "${OUTDIR}/ref_ld.prune.in" \
        --keep "${REF_KING}.king.cutoff.in.id" \
        --make-pgen \
        --threads "${THREADS}" \
        --out "${REF_QC}"

    log "STEP 1 DONE: $(wc -l < "${REF_QC}.psam") samples, $(wc -l < "${REF_QC}.pvar") variants retained."
fi

# =============================================================================
# STEP 2 — PCA ON REFERENCE PANEL (allele weights)
# (skipped if ref_pca.eigenvec.allele already exists)
# =============================================================================
REF_PCA="${OUTDIR}/ref_pca"

if [[ -f "${REF_PCA}.eigenvec.allele" && -f "${REF_PCA}.acount" ]]; then
    log "STEP 2 SKIPPED: Reference PCA outputs already exist (${REF_PCA}.eigenvec.allele)."
else
    log "STEP 2: Computing PCA allele weights on reference panel (N_PCS=${N_PCS})..."
    # --freq counts and --pca allele-wts must be in the same call.
    # approx (randomized SVD) is recommended for large reference panels.
    plink2 \
        --pfile "${REF_QC}" \
        --freq counts \
        --pca approx allele-wts "${N_PCS}" \
        --threads "${THREADS}" \
        --out "${REF_PCA}"

    log "STEP 2 DONE: Allele weights written to ${REF_PCA}.eigenvec.allele"
    log "  Verify column layout of .eigenvec.allele before projection:"
    head -1 "${REF_PCA}.eigenvec.allele" | tr '\t' '\n' | nl
fi

# =============================================================================
# STEP 3 — VARIANT INTERSECTION (reference ∩ ancient)
# =============================================================================
INTERSECT_VARS="${OUTDIR}/intersect_variants.txt"

if [[ -f "${INTERSECT_VARS}" ]]; then
    log "STEP 3 SKIPPED: Variant intersection file already exists (${INTERSECT_VARS})."
else
    log "STEP 3: Computing variant intersection between reference and ancient dataset..."

    # Extract variant IDs from each dataset (skip header, col 3 = ID in .pvar)
    awk 'NR>1 {print $3}' "${REF_QC}.pvar"          | sort > "${OUTDIR}/_ref_ids.txt"
    awk 'NR>1 {print $3}' "${ANCIENT_PFILE}.pvar"   | sort > "${OUTDIR}/_anc_ids.txt"

    comm -12 "${OUTDIR}/_ref_ids.txt" "${OUTDIR}/_anc_ids.txt" > "${INTERSECT_VARS}"

    N_SHARED=$(wc -l < "${INTERSECT_VARS}")
    N_REF=$(wc -l < "${OUTDIR}/_ref_ids.txt")
    N_ANC=$(wc -l < "${OUTDIR}/_anc_ids.txt")
    log "  Reference variants:  ${N_REF}"
    log "  Ancient variants:    ${N_ANC}"
    log "  Shared variants:     ${N_SHARED}"

    if [[ "${N_SHARED}" -lt 10000 ]]; then
        log "WARNING: Fewer than 10,000 shared variants. PCA quality may be compromised."
        log "  Consider relaxing QC thresholds or checking variant ID format concordance."
    fi

    rm -f "${OUTDIR}/_ref_ids.txt" "${OUTDIR}/_anc_ids.txt"
fi

# =============================================================================
# STEP 4 — ANCIENT SAMPLE MISSINGNESS FILTER (if thresholds < 1.0)
# =============================================================================
ANC_FILTERED="${OUTDIR}/ancient_filtered"

if check_pfile "${ANC_FILTERED}"; then
    log "STEP 4 SKIPPED: Filtered ancient pfile already exists (${ANC_FILTERED})."
else
    log "STEP 4: Applying missingness filters to ancient dataset..."
    log "  geno=${ANC_GENO}, mind=${ANC_MIND}"

    ANC_FLAGS=""
    [[ $(awk "BEGIN{print (${ANC_GENO} < 1.0)}") -eq 1 ]] && ANC_FLAGS="${ANC_FLAGS} --geno ${ANC_GENO}"
    [[ $(awk "BEGIN{print (${ANC_MIND} < 1.0)}") -eq 1 ]] && ANC_FLAGS="${ANC_FLAGS} --mind ${ANC_MIND}"

    if [[ -z "${ANC_FLAGS}" ]]; then
        log "  Both ANC_GENO and ANC_MIND are 1.0; skipping filter, symlinking input."
        for ext in pgen pvar psam; do
            ln -sf "$(realpath "${ANCIENT_PFILE}.${ext}")" "${ANC_FILTERED}.${ext}"
        done
    else
        plink2 \
            --pfile "${ANCIENT_PFILE}" \
            --extract "${INTERSECT_VARS}" \
            --human \
            ${ANC_FLAGS} \
            --make-pgen \
            --threads "${THREADS}" \
            --out "${ANC_FILTERED}"
    fi

    log "STEP 4 DONE: $(wc -l < "${ANC_FILTERED}.psam") ancient samples retained."
fi

# =============================================================================
# STEP 5 — PROJECT ANCIENT SAMPLES ONTO REFERENCE PCs
# =============================================================================
ANC_PROJ="${OUTDIR}/ancient_proj"

if [[ -f "${ANC_PROJ}.sscore" ]]; then
    log "STEP 5 SKIPPED: Projection output already exists (${ANC_PROJ}.sscore)."
else
    log "STEP 5: Projecting ancient samples onto reference PCs..."
    log "  Score file columns: ID=${SCORE_ID_COL}, A1=${SCORE_A1_COL}, PCs=${SCORE_START_COL}-${SCORE_END_COL}"
    log "  Using no-mean-imputation (mandatory for aDNA; missing data not imputed)."

    plink2 \
        --pfile "${ANC_FILTERED}" \
        --extract "${INTERSECT_VARS}" \
        --read-freq "${REF_PCA}.acount" \
        --score "${REF_PCA}.eigenvec.allele" \
            "${SCORE_ID_COL}" \
            "${SCORE_A1_COL}" \
            header-read \
            no-mean-imputation \
            variance-standardize \
        --score-col-nums "${SCORE_START_COL}-${SCORE_END_COL}" \
        --human \
        --threads "${THREADS}" \
        --out "${ANC_PROJ}"

    log "STEP 5 DONE: Projected scores in ${ANC_PROJ}.sscore"
    log "  Column layout of .sscore:"
    head -1 "${ANC_PROJ}.sscore" | tr '\t' '\n' | nl
fi

# =============================================================================
# STEP 6 — OPTIONAL R DIAGNOSTIC PLOT
# =============================================================================
if [[ "${MAKE_PLOT}" != "yes" ]]; then
    log "STEP 6 SKIPPED: MAKE_PLOT=${MAKE_PLOT}."
else
    log "STEP 6: Generating PCA diagnostic plot in R..."

    Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
# Resolve paths from shell environment
outdir      <- Sys.getenv("OUTDIR",   "./pca_projection_out")
ref_eigvec  <- file.path(outdir, "ref_pca.eigenvec")
anc_scores  <- file.path(outdir, "ancient_proj.sscore")
plot_out    <- file.path(outdir, "pca_projection.pdf")

# Load reference PC coordinates
ref <- fread(ref_eigvec)
# Rename: first two cols are #FID IID, then PC1..PCN
old_names <- names(ref)
new_names <- c("FID", "IID", paste0("PC", seq_len(ncol(ref) - 2)))
setnames(ref, old_names, new_names)
ref[, source := "Reference (modern)"]

# Load ancient projected scores
anc <- fread(anc_scores)
# .sscore columns: #FID IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM PC1_AVG ... PCN_AVG
# Keep FID, IID and PC columns
pc_cols <- grep("_AVG$", names(anc), value = TRUE)
n_pcs   <- length(pc_cols)
anc_sub <- anc[, c("#FID", "IID", pc_cols), with = FALSE]
setnames(anc_sub, c("#FID", pc_cols), c("FID", paste0("PC", seq_len(n_pcs))))
anc_sub[, source := "Ancient (projected)"]

# Align column names for rbind
shared_cols <- intersect(names(ref), names(anc_sub))
combined <- rbind(ref[, ..shared_cols], anc_sub[, ..shared_cols])

# Plot PC1 vs PC2
p <- ggplot(combined, aes(x = PC1, y = PC2, colour = source)) +
    geom_point(data = combined[source == "Reference (modern)"],
               alpha = 0.4, size = 1.2) +
    geom_point(data = combined[source == "Ancient (projected)"],
               alpha = 0.9, size = 2, shape = 17) +
    scale_colour_manual(values = c("Reference (modern)" = "#7B9EA6",
                                   "Ancient (projected)" = "#C0392B")) +
    labs(title = "PCA projection: ancient onto modern reference",
         subtitle = sprintf("Reference n=%d  |  Ancient projected n=%d",
                            nrow(ref), nrow(anc_sub)),
         x = "PC1", y = "PC2", colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

ggsave(plot_out, p, width = 7, height = 6)
cat(sprintf("Plot saved: %s\n", plot_out))
RSCRIPT

    log "STEP 6 DONE."
fi

# =============================================================================
# SUMMARY
# =============================================================================
log "Pipeline complete. Output directory: ${OUTDIR}"
log ""
log "Key output files:"
log "  ${OUTDIR}/ref_pca.eigenvec          Reference PC coordinates"
log "  ${OUTDIR}/ref_pca.eigenvec.allele   Allele weights (projection basis)"
log "  ${OUTDIR}/ref_pca.eigenval          Eigenvalues"
log "  ${OUTDIR}/intersect_variants.txt    Shared variant IDs"
log "  ${OUTDIR}/ancient_proj.sscore       Ancient projected PC scores"
[[ "${MAKE_PLOT}" == "yes" ]] && log "  ${OUTDIR}/pca_projection.pdf        Diagnostic plot"
log ""
log "IMPORTANT: Before interpreting projections, verify:"
log "  1. .eigenvec.allele column indices match SCORE_ID_COL/SCORE_A1_COL"
log "     (header printed to log at STEP 2 and STEP 5 above)"
log "  2. Variant ID format is consistent between ref and ancient .pvar files"
log "     (expected: CHR-POS-REF-ALT with hyphens; mismatches silently reduce shared set)"
log "  3. Check ALLELE_CT column in .sscore to confirm sufficient coverage"
log "     per ancient sample (very low counts = unreliable projections)"
