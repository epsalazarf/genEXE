#!/usr/bin/env Rscript
# pvar_merge_summary.R — Summarize and visualize pvar_merge_check.sh TSV output.
#
# Usage:
#   Rscript pvar_merge_summary.R <merge_check.tsv> [output_prefix]
#
# Outputs:
#   <prefix>_summary.tsv    Counts and percentages by set category
#   <prefix>_barplot.png    Stacked 100% bar plot per chromosome (A4 landscape PNG)
#   <prefix>_barplot.pdf    Same plot as PDF
#   <prefix>_extract.txt    VARID list of perfect matches for PLINK2 --extract
#
# Match codes from pvar_merge_check.sh:
#    2  Perfect match  — same CHROM, POS, REF, ALT
#    1  Semi-match     — same CHROM, POS, REF; different ALT
#   -1  REF mismatch   — same CHROM, POS; different REF allele
#    0  Missing        — present in only one dataset
#
# Dependencies: data.table, ggplot2, scales

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript pvar_merge_summary.R <merge_check.tsv> [output_prefix]\n")
  quit(status = 1)
}

tsv_file <- args[1]
prefix   <- if (length(args) >= 2) args[2] else tools::file_path_sans_ext(tsv_file)

if (!file.exists(tsv_file)) stop("File not found: ", tsv_file)

cat("Input  :", tsv_file, "\n")
cat("Prefix :", prefix, "\n\n")

# ── Load & classify ────────────────────────────────────────────────────────────
dt <- fread(tsv_file, colClasses = list(character = "CHROM"))

cat_levels <- c("DS2 only", "DS1 only", "REF mismatch", "Semi-match", "Perfect match")
cat_colors <- c(
  "Perfect match" = "#43A047",
  "Semi-match"    = "#FB8C00",
  "REF mismatch"  = "#E53935",
  "DS1 only"      = "#64B5F6",
  "DS2 only"      = "#BA68C8"
)

dt[, category := fcase(
  MATCH ==  2,                        "Perfect match",
  MATCH ==  1,                        "Semi-match",
  MATCH == -1,                        "REF mismatch",
  MATCH ==  0 & DS1 == 1 & DS2 == 0, "DS1 only",
  MATCH ==  0 & DS1 == 0 & DS2 == 1, "DS2 only"
)]
dt[, category := factor(category, levels = cat_levels)]

# ── 1. Summary table ──────────────────────────────────────────────────────────
n_union  <- nrow(dt)
n_ds1    <- dt[DS1 == 1, .N]
n_ds2    <- dt[DS2 == 1, .N]
n_inter  <- dt[MATCH == 2, .N]   # intersection = perfect matches only
n_semi   <- dt[MATCH == 1, .N]
n_rmis   <- dt[MATCH == -1, .N]
n_ds1only <- dt[MATCH == 0 & DS1 == 1, .N]
n_ds2only <- dt[MATCH == 0 & DS2 == 1, .N]

pct <- function(num, den) round(num / den * 100, 2)
na_char <- "—"

summary_dt <- data.table(
  Category   = c(
    "Union (DS1 ∪ DS2)",
    "DS1 variants",
    "DS2 variants",
    "Intersection — Perfect match [2]",
    "Shared site, diff ALT — Semi-match [1]",
    "Shared site, diff REF — REF mismatch [-1]",
    "DS1 only — absent in DS2 [0]",
    "DS2 only — absent in DS1 [0]"
  ),
  Count      = c(n_union, n_ds1, n_ds2, n_inter, n_semi, n_rmis, n_ds1only, n_ds2only),
  `Pct_Union` = c(100, pct(n_ds1, n_union), pct(n_ds2, n_union),
                   pct(n_inter, n_union), pct(n_semi, n_union),
                   pct(n_rmis, n_union), pct(n_ds1only, n_union), pct(n_ds2only, n_union)),
  `Pct_DS1`  = c(pct(n_union, n_ds1), 100, NA,
                  pct(n_inter, n_ds1), pct(n_semi, n_ds1),
                  pct(n_rmis, n_ds1), pct(n_ds1only, n_ds1), NA),
  `Pct_DS2`  = c(pct(n_union, n_ds2), NA, 100,
                  pct(n_inter, n_ds2), pct(n_semi, n_ds2),
                  pct(n_rmis, n_ds2), NA, pct(n_ds2only, n_ds2))
)
setnames(summary_dt, c("Category", "Count", "% of Union", "% of DS1", "% of DS2"))

summ_path <- paste0(prefix, "_summary.tsv")
fwrite(summary_dt, summ_path, sep = "\t", na = na_char)
cat("Summary table:", summ_path, "\n")

# Print to console as well
cat("\n")
print(summary_dt, na.print = na_char)
cat("\n")

# ── 2. Stacked 100% bar plot ──────────────────────────────────────────────────
chrom_sort <- function(x) {
  u      <- unique(x)
  num    <- suppressWarnings(as.integer(u))
  extras <- c(X = 23, Y = 24, XY = 25, MT = 26, M = 26)
  rank   <- ifelse(!is.na(num), num, extras[u])
  rank   <- ifelse(is.na(rank), 99, rank)
  u[order(rank)]
}

chr_sorted <- chrom_sort(dt$CHROM)
chr_levels <- c("All", chr_sorted)

# Combine "All" + per-chromosome rows
plot_dt <- rbindlist(list(
  dt[, .(CHROM = "All", category)],
  dt[, .(CHROM, category)]
))
plot_dt[, CHROM := factor(CHROM, levels = rev(chr_levels))]  # rev → "All" at top after flip

# Compute proportions per CHROM group
prop_dt <- plot_dt[, .N, by = .(CHROM, category)]
prop_dt[, prop := N / sum(N), by = CHROM]

# Expand to all category × CHROM combinations (ensures full legend even for 0-count cats)
full_grid <- CJ(CHROM = levels(plot_dt$CHROM), category = factor(cat_levels, levels = cat_levels))
prop_dt   <- merge(full_grid, prop_dt, by = c("CHROM", "category"), all.x = TRUE)
prop_dt[is.na(prop), `:=`(N = 0L, prop = 0)]
prop_dt[, CHROM    := factor(CHROM,    levels = rev(chr_levels))]
prop_dt[, category := factor(category, levels = cat_levels)]

n_chrom <- length(chr_sorted)

p <- ggplot(prop_dt, aes(x = CHROM, y = prop, fill = category)) +
  geom_bar(stat = "identity", width = 0.72) +
  geom_text(
    data = prop_dt[prop >= 0.02],
    aes(label = paste0(round(prop * 100, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 2.6, color = "white", fontface = "bold"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = cat_colors,
    breaks = rev(cat_levels),   # legend: Perfect match on top
    labels = c(
      "Perfect match" = "Perfect match [2]",
      "Semi-match"    = "Semi-match [1]",
      "REF mismatch"  = "REF mismatch [-1]",
      "DS1 only"      = "DS1 only [0]",
      "DS2 only"      = "DS2 only [0]"
    )
  ) +
  coord_flip() +
  labs(
    title    = "Variant Merge Compatibility by Chromosome",
    subtitle = sprintf(
      "DS1: %s  |  DS2: %s  |  Union: %s  |  Intersection (perfect): %s (%.1f%% of union)",
      format(n_ds1, big.mark = ","), format(n_ds2, big.mark = ","),
      format(n_union, big.mark = ","), format(n_inter, big.mark = ","),
      pct(n_inter, n_union)
    ),
    x    = NULL,
    y    = "Proportion of variants",
    fill = "Match category"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position       = "bottom",
    legend.title          = element_text(face = "bold", size = 10),
    legend.key.size       = unit(0.45, "cm"),
    legend.text           = element_text(size = 9),
    legend.margin         = margin(t = 4),
    plot.title            = element_text(face = "bold", size = 13),
    plot.subtitle         = element_text(size = 9, color = "grey40"),
    axis.text.y           = element_text(size = 9),
    axis.text.x           = element_text(size = 8),
    panel.grid.major.y    = element_blank(),
    panel.grid.minor      = element_blank(),
    plot.margin           = margin(10, 15, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = FALSE))

# A4 landscape: 297 × 210 mm = 11.69 × 8.27 in
png_path <- paste0(prefix, "_barplot.png")
pdf_path <- paste0(prefix, "_barplot.pdf")

ggsave(png_path, p, width = 11.69, height = 8.27, units = "in", dpi = 150)
ggsave(pdf_path, p, width = 11.69, height = 8.27, units = "in", device = "pdf")
cat("Bar plot PNG:", png_path, "\n")
cat("Bar plot PDF:", pdf_path, "\n")

# ── 3. Extract list (perfect matches for PLINK2 --extract) ────────────────────
extract_path <- paste0(prefix, "_extract.txt")
fwrite(dt[MATCH == 2, .(VARID)], extract_path, col.names = FALSE)
cat("Extract list :", extract_path, "\n")
cat("  ->", format(n_inter, big.mark = ","), "perfect-match variants\n")
