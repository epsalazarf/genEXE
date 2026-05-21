# PCA Projection of Ancient DNA onto a Modern Reference Panel

**Script:** `pca_projection_adna.sh`  
**Stack:** PLINK2 + R (ggplot2, data.table)  
**Scope:** Whole-genome SNP data; human autosomes; PLINK2 binary input (`.pgen`/`.pvar`/`.psam`)

---

## Table of Contents

1. [Rationale](#1-rationale)
2. [Prerequisites](#2-prerequisites)
3. [Quick Start](#3-quick-start)
4. [Configuration Reference](#4-configuration-reference)
5. [Pipeline Steps](#5-pipeline-steps)
6. [Output Files](#6-output-files)
7. [Critical Notes](#7-critical-notes)
8. [Known Failure Modes](#8-known-failure-modes)
9. [References](#9-references)

---

## 1. Rationale

### Why project ancient samples rather than include them in PCA directly?

Standard PCA computes principal components by maximizing variance across all samples simultaneously. When ancient DNA (aDNA) samples are included in that computation, two problems arise:

**Sparse genotyping distorts the variance structure.** Ancient samples typically have genome-wide missingness of 50–95%, depending on preservation and sequencing depth. The covariance matrix is then estimated from an incomplete data structure where missingness is not random — it is shaped by post-mortem damage, coverage heterogeneity, and the genomic regions targeted by capture arrays. Including these samples in the eigendecomposition can shift PC axes away from the population structure signal present in the modern reference and toward artifacts of data sparsity.

**Modern population structure is the reference frame.** In temporal studies of allele frequency trajectories (e.g., tracking immune-related variants across the Holocene), the goal is to place ancient samples *relative to* modern continental populations. The PC axes should be defined by the modern panel and remain fixed; ancient samples are then scored against those axes. This is analogous to predicting on held-out data using a model trained on the training set — the axes are not refitted to the test set.

The correct approach is therefore:

1. Compute PC **allele weights** on the modern reference panel alone.
2. Apply those weights to ancient samples via a dot product (the `--score` mechanism in PLINK2), yielding projected coordinates on the pre-defined axes.

This is the method implemented in this script, following the approach used in Barrie et al. (2024) and consistent with the framework in Patterson et al. (2006).

### Why `--pca allele-wts` and not `--pca var-wts`?

`allele-wts` produces a `.eigenvec.allele` file with one row per allele and columns for each PC weight. This format is directly consumable by `--score` without any reformatting. `var-wts` (the older form) is superseded in current PLINK2 versions for this use case.

### Why not merge the datasets before running PCA?

Merging requires full allele reconciliation, strand resolution, and produces a combined genotype matrix where the ancient samples' sparse rows would again influence the covariance structure. The projection approach avoids this entirely: the ancient samples are only used at scoring time, never at decomposition time.

---

## 2. Prerequisites

### Software

| Tool | Minimum version | Notes |
|------|----------------|-------|
| `plink2` | 2.0.0-a.6 | `allele-wts` and `acount` syntax; `--score` column layout changed in recent versions — see [Critical Notes](#7-critical-notes) |
| `R` | 4.0 | Only required if `MAKE_PLOT=yes` |
| `ggplot2` | any current | R package |
| `data.table` | any current | R package |
| `awk`, `sort`, `comm` | POSIX | Standard coreutils |

### Input data requirements

- **Reference panel:** PLINK2 pfile (`.pgen`/`.pvar`/`.psam`). Should be a well-characterized modern population panel (e.g., 1000 Genomes, HGDP, or a curated Eurasian panel). Sample size should be sufficient for stable eigendecomposition — at minimum several hundred samples; >1,000 preferred.
- **Ancient dataset:** PLINK2 pfile. Samples should already have passed upstream QC (damage filtering, contamination estimation, minimum coverage filtering). This script does not perform aDNA-specific QC such as mapDamage correction or contamination estimation.
- **Variant ID convention:** Both datasets must use the same ID format. This script expects `CHR-POS-REF-ALT` (hyphen-delimited). Mismatched conventions are a silent failure mode — see [Known Failure Modes](#8-known-failure-modes).

---

## 3. Quick Start

```bash
# Minimal invocation — edit paths, everything else uses defaults
REF_PFILE=/path/to/modern_reference \
ANCIENT_PFILE=/path/to/ancient_dataset \
OUTDIR=./my_pca_run \
bash pca_projection_adna.sh
```

To override any parameter without editing the script, pass it as an environment variable:

```bash
REF_PFILE=/data/hgdp \
ANCIENT_PFILE=/data/bronze_age \
OUTDIR=./bronze_pca \
THREADS=8 \
EXCLUDE_AMBIGUOUS=no \
ANC_GENO=0.95 \
MAKE_PLOT=yes \
bash pca_projection_adna.sh
```

The script is **idempotent**: any step whose output files already exist is skipped. This means you can interrupt and re-run safely, and reuse a completed reference panel preparation (`ref_qc.*`, `ref_pca.*`) for multiple ancient datasets by changing only `ANCIENT_PFILE` and `OUTDIR`.

---

## 4. Configuration Reference

All parameters have defaults and can be overridden via environment variables.

### Paths

| Variable | Default | Description |
|----------|---------|-------------|
| `REF_PFILE` | `/path/to/modern_reference` | Prefix of reference PLINK2 pfile (no extension) |
| `ANCIENT_PFILE` | `/path/to/ancient_dataset` | Prefix of ancient PLINK2 pfile (no extension) |
| `OUTDIR` | `./pca_projection_out` | Output directory; created if absent |
| `REF_FA` | *(empty)* | Path to reference FASTA for `--ref-from-fa`; optional |

### Compute

| Variable | Default | Description |
|----------|---------|-------------|
| `THREADS` | `4` | CPU threads passed to PLINK2 |
| `N_PCS` | `10` | Number of principal components to compute and project |

### Reference panel QC

| Variable | Default | Description |
|----------|---------|-------------|
| `MAF_REF` | `0.05` | Minor allele frequency cutoff for reference panel |
| `REF_GENO` | `0.05` | Maximum per-variant missingness in reference (0–1) |
| `REF_MIND` | `0.10` | Maximum per-sample missingness in reference (0–1) |
| `LD_WIN` | `1000` | LD pruning window size in kb |
| `LD_STEP` | `50` | LD pruning step size in variants |
| `LD_R2` | `0.1` | LD pruning r² threshold |
| `KING_CUTOFF` | `0.0884` | Kinship coefficient cutoff for relatedness pruning (~2nd degree) |

### Ancient dataset filters

| Variable | Default | Description |
|----------|---------|-------------|
| `ANC_GENO` | `1.0` | Maximum per-variant missingness in ancient dataset; `1.0` disables |
| `ANC_MIND` | `1.0` | Maximum per-sample missingness in ancient dataset; `1.0` disables |

Setting both to `1.0` (the default) applies no missingness filter to ancient samples, which is appropriate when samples have already been filtered upstream and further exclusion would remove too many individuals. Use a threshold such as `ANC_MIND=0.98` if you want to exclude samples with extreme genome-wide missingness.

### Strand-ambiguous SNPs

| Variable | Default | Description |
|----------|---------|-------------|
| `EXCLUDE_AMBIGUOUS` | `yes` | `yes` = exclude A/T and C/G SNPs; `no` = retain them |

See [Critical Notes](#7-critical-notes) for the rationale and trade-offs of each choice.

### Output control

| Variable | Default | Description |
|----------|---------|-------------|
| `MAKE_PLOT` | `yes` | `yes` = run R diagnostic plot after projection |

---

## 5. Pipeline Steps

### Step 1 — Reference panel QC and LD pruning

*Skipped if `ref_qc.pgen` already exists.*

**1a. Basic variant and sample filters**

Applied to the reference panel only:

- `--maf`: removes low-frequency variants that contribute noise rather than population structure signal.
- `--geno` / `--mind`: removes variants and samples with high missingness.
- `--autosome`: restricts to autosomes 1–22, excluding X, Y, and MT which have different inheritance and ploidy properties.
- `--snps-only just-acgt`: excludes indels and any non-ACGT allele encodings.
- `--ref-from-fa` (optional): reconciles REF/ALT orientation against the reference FASTA. Important if the reference panel was processed independently and may have inconsistent REF assignment.

**1b. Strand-ambiguous SNP removal (conditional)**

A/T and C/G SNPs cannot be strand-resolved by allele identity alone because the complement of A is T and the complement of C is G. When merging or aligning datasets from different laboratories or capture designs, these variants may be silently flipped, introducing noise into the PCA. The conservative default (`EXCLUDE_AMBIGUOUS=yes`) removes them entirely. Setting `EXCLUDE_AMBIGUOUS=no` retains them, which is acceptable if you are confident both datasets were processed against the same strand convention. The script parses REF and ALT from the `CHR-POS-REF-ALT` variant ID rather than from `.pvar` columns 4/5 to ensure consistency.

**1c. LD pruning**

PCA assumes that input variants are approximately independent. Linkage disequilibrium (LD) between nearby variants inflates the contribution of densely-genotyped or high-LD genomic regions (notably the HLA region and pericentromeric blocks) to the first few PCs. The default parameters (`1000kb`, step 50, r²=0.1) are intentionally conservative — more aggressive than the commonly cited 500kb/50/0.2 — to better capture broad continental-scale structure rather than fine-scale within-population relatedness.

If your study requires fine-scale within-population structure (e.g., distinguishing Early European Farmer sub-clusters), you may want to relax LD pruning and increase N_PCS accordingly.

**1d. Relatedness pruning**

Closely related individuals in the reference panel inflate certain PCs by creating sample pairs with extreme covariance. The KING kinship cutoff of 0.0884 corresponds approximately to second-degree relatives (half-siblings, grandparent-grandchild). Adjust downward (e.g., `0.0442` for third-degree) if your reference panel is known to contain extended family structure.

### Step 2 — PCA allele weights on the reference panel

*Skipped if `ref_pca.eigenvec.allele` and `ref_pca.acount` already exist.*

```bash
plink2 \
    --pfile ref_qc \
    --freq counts \
    --pca approx allele-wts 10 \
    --out ref_pca
```

The `--freq counts` and `--pca allele-wts` flags must appear in the same PLINK2 call. This produces:

- `ref_pca.acount`: allele count file used to compute allele frequencies for standardization during projection.
- `ref_pca.eigenvec`: sample-level PC coordinates for the reference panel.
- `ref_pca.eigenvec.allele`: per-allele PC weights used in Step 5.
- `ref_pca.eigenval`: eigenvalues.

The `approx` modifier uses a randomized SVD algorithm (IRAM/Implicitly Restarted Arnoldi Method) rather than exact eigendecomposition. For reference panels of several thousand samples this is computationally necessary and introduces negligible approximation error. For panels under ~500 samples, omit `approx` for exact decomposition.

**The ancient samples play no role in this step.** The PC axes are defined entirely by the modern reference.

### Step 3 — Variant intersection

The intersection is computed by extracting variant IDs from `ref_qc.pvar` and `ancient.pvar`, sorting both lists, and taking the common set via `comm -12`. Only variants present in both datasets after reference QC and LD pruning are used for projection.

The script emits a warning if fewer than 10,000 shared variants are found. Below this threshold, PC resolution is typically insufficient to separate continental populations reliably, and the projections should be interpreted with caution.

A common reason for low intersection counts is variant ID format mismatch — see [Known Failure Modes](#8-known-failure-modes).

### Step 4 — Ancient dataset missingness filter

*Skipped if `ancient_filtered.pgen` already exists.*

Ancient samples are filtered separately from the reference panel, using the `ANC_GENO` and `ANC_MIND` thresholds, which default to `1.0` (no filtering). The intersection variant list from Step 3 is applied here so that the filtered ancient pfile already contains only shared variants, reducing the data volume passed to Step 5.

If both thresholds are `1.0`, the step creates symlinks from the input pfile rather than copying it, avoiding unnecessary I/O.

### Step 5 — Projection via `--score`

```bash
plink2 \
    --pfile ancient_filtered \
    --extract intersect_variants.txt \
    --read-freq ref_pca.acount \
    --score ref_pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
    --score-col-nums 6-15 \
    --out ancient_proj
```

Each ancient sample's projected PC coordinate is computed as a weighted sum of its allele dosages, where the weights are the allele-specific PC loadings from the reference eigendecomposition. The result is written to `ancient_proj.sscore`.

Key flags:

- `--read-freq ref_pca.acount`: supplies the reference allele frequencies used for centering and standardization. These must come from the reference panel, not from the ancient samples.
- `no-mean-imputation`: missing genotypes in ancient samples are **not** replaced by the reference mean dosage. This is critical — see [Critical Notes](#7-critical-notes).
- `variance-standardize`: standardizes the score by the expected variance under Hardy-Weinberg equilibrium given the reference allele frequencies.
- `--score-col-nums 6-15`: selects PC weight columns 6 through 15 (for N_PCS=10). The column range is computed automatically from `N_PCS`.

### Step 6 — Diagnostic plot (optional)

Produces a PDF scatterplot of PC1 vs PC2 with reference samples shown as semi-transparent circles and ancient projected samples as solid triangles. This is a quick visual sanity check, not a publication-ready figure. The R script reads paths from the `OUTDIR` environment variable passed by the shell.

---

## 6. Output Files

All outputs are written to `OUTDIR/`.

| File | Description |
|------|-------------|
| `ref_qc.pgen/pvar/psam` | QC'd, LD-pruned, relatedness-pruned reference panel |
| `ref_pca.eigenvec` | Reference sample PC coordinates |
| `ref_pca.eigenvec.allele` | Per-allele PC weights (projection basis) |
| `ref_pca.eigenval` | Eigenvalues |
| `ref_pca.acount` | Reference allele counts (required for projection standardization) |
| `intersect_variants.txt` | Variant IDs shared between reference and ancient dataset |
| `ancient_filtered.pgen/pvar/psam` | Ancient dataset restricted to shared variants (symlink if no missingness filter applied) |
| `ancient_proj.sscore` | Projected PC scores for ancient samples |
| `pca_projection.pdf` | Diagnostic plot (if `MAKE_PLOT=yes`) |
| `ambiguous_snps.txt` | Strand-ambiguous variant IDs excluded (if `EXCLUDE_AMBIGUOUS=yes`) |

Intermediate files (`ref_mafgeno.*`, `ref_noamb.*`, `ref_ld.*`, `ref_king.*`) are retained in `OUTDIR/` and serve as checkpoints for the idempotency logic. They may be deleted after a successful run to recover disk space.

---

## 7. Critical Notes

### Column indices in `--score` may vary by PLINK2 version

The `.eigenvec.allele` file layout changed in recent PLINK2 releases. The current format (as of v2.0.0-a.6+) is:

```
#CHROM  ID  REF  ALT  A1  PC1  PC2  ...  PCN
col 1   2   3    4    5   6    7    ...  5+N
```

The script sets `SCORE_ID_COL=2`, `SCORE_A1_COL=5`, and PC columns starting at 6. **The script prints the full column header at STEP 2 and STEP 5 of the log.** Verify these before interpreting output. If your PLINK2 version produces a different layout, adjust the top-level variables `SCORE_ID_COL` and `SCORE_A1_COL` in the CONFIG section.

### `no-mean-imputation` is mandatory for aDNA

The default PLINK2 `--score` behavior imputes missing genotypes using the mean dosage derived from reference allele frequencies: for a variant with reference allele frequency p, a missing genotype is treated as dosage 2p. For modern samples with low, approximately-random missingness, this is a reasonable approximation.

For ancient samples, this is incorrect for two reasons. First, missingness is non-random: it correlates with GC content, proximity to repetitive elements, and post-mortem fragmentation patterns. Imputing from the reference mean introduces a systematic bias that pulls ancient samples toward the reference centroid, artifactually compressing their spread along PC axes. Second, for samples with very high missingness (>80%), the imputed contribution can dominate the actual genotyped signal, effectively replacing the sample's genetic identity with reference frequencies.

`no-mean-imputation` ensures that missing sites contribute zero to the score and that the projected coordinate reflects only actually-observed genotypes. The trade-off is that samples with very low coverage (few observed sites) will have high-variance projections — this is the correct behavior, and should be interpreted alongside the `ALLELE_CT` column in the `.sscore` output.

### `variance-standardize` and the monomorphic variant edge case

`variance-standardize` divides each allele's contribution by the square root of its expected variance under HWE: `sqrt(2 * p * (1-p))`, where p is the reference allele frequency from `--read-freq`. If a variant has frequency exactly 0 or 1 in the reference panel (or if numerical precision produces a NaN), PLINK2 will throw an error even when `--read-freq` is supplied. This can occur when projecting a single sample at a time, or when the ancient dataset contains alleles absent from the reference.

If this error appears, the options are:

1. Filter the `.eigenvec.allele` score file to remove variants with reference frequency 0 or 1 before projection.
2. Drop `variance-standardize` entirely and use `no-mean-imputation` alone. The projected coordinates will be on a different scale than the reference `.eigenvec`, but can be rescaled post-hoc by dividing by the square root of the corresponding eigenvalue.

### Variant ID format concordance

The intersection in Step 3 is a string comparison on variant IDs. This script expects the `CHR-POS-REF-ALT` convention with hyphens throughout. If the reference panel uses colons (`CHR:POS:REF:ALT`), rsIDs, or any other scheme while the ancient dataset uses a different convention, `comm -12` will silently return zero or near-zero shared variants. The script warns if the shared count falls below 10,000, but does not attempt to reconcile IDs automatically.

To pre-check concordance before running:

```bash
# Sample 20 IDs from each dataset and compare format
awk 'NR>1 && NR<=21 {print $3}' reference.pvar
awk 'NR>1 && NR<=21 {print $3}' ancient.pvar
```

If formats differ, normalize both datasets to `CHR-POS-REF-ALT` using `bcftools annotate --set-id` or PLINK2's `--set-all-var-ids` before running this script.

### The `ALLELE_CT` column in `.sscore` is a quality indicator

The `.sscore` output contains an `ALLELE_CT` column giving the number of alleles scored for each sample (i.e., twice the number of non-missing variants). For a projection based on, say, 80,000 shared variants, a sample with `ALLELE_CT` of 160,000 has complete coverage; a sample with 4,000 has covered only ~2.5% of sites. Projections based on very low `ALLELE_CT` values are unreliable and should not be interpreted as representative population placements. A reasonable minimum threshold for interpretation is approximately 10,000–20,000 alleles (5,000–10,000 variants), though this depends on the PC being examined and the population contrast being resolved.

### PC scale differences between reference and projected samples

The reference `.eigenvec` coordinates and the projected `.sscore` coordinates are not guaranteed to be on the same scale, even with `variance-standardize`. The reference coordinates are normalized by the eigenvalue during the SVD; the projected scores are raw dot products divided by variance. To place both on the same scale for visualization, divide the projected scores by `sqrt(eigenvalue_i)` for PC_i, where eigenvalues are read from `ref_pca.eigenval`. The diagnostic R plot does not perform this rescaling — it is adequate for a sanity check but not for publication.

---

## 8. Known Failure Modes

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `intersect_variants.txt` nearly empty | Variant ID format mismatch | Normalize IDs in both datasets before running |
| Ancient samples cluster at origin in plot | `no-mean-imputation` missing, or very low `ALLELE_CT` | Verify flag is present; check `ALLELE_CT` per sample |
| `variance-standardize failure` error | Monomorphic variant in score file | Filter zero/one-frequency variants from `.eigenvec.allele` |
| `--sort-vars` conflict error | Attempted to combine with a reporting flag | Run conversion and reporting in separate PLINK2 calls |
| Ancient samples compress toward reference centroid | Mean imputation active | Ensure `no-mean-imputation` is passed to `--score` |
| PC axes dominated by HLA or centromeric regions | LD pruning too lenient | Tighten `LD_R2` or add an explicit exclude region file for HLA (chr6:25–35Mb) |
| Very few samples retained after KING pruning | Reference panel has family structure | Raise `KING_CUTOFF` or supply a pre-filtered reference |

---

## 9. References

- **Barrie, W. et al.** (2024). Elevated genetic risk for multiple sclerosis emerged in steppe pastoralist populations. *Nature*, 625, 321–328. — Primary methodological template for temporal aDNA PRS and PCA projection.
- **Patterson, N. et al.** (2006). Population structure and eigenanalysis. *PLOS Genetics*, 2(12), e190. — Foundational reference for PCA in population genetics; establishes the standardization approach used by PLINK2.
- **Chang, C.C. et al.** (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, 4, 7. — PLINK2 primary citation.
- **PLINK2 documentation — PCA projection:** https://www.cog-genomics.org/plink/2.0/score#pca_project
