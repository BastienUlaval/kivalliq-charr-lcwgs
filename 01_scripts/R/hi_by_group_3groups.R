# =============================================================================
# hi_by_group_3groups.R — Hybrid index of the three Rankin Inlet PCA groups
#   (Main group, PC1 cluster, PC2 tail) defined by 45_plot_rankin_chrom_cor.R
#
# Usage (from angsd_pipeline2/):
#   $RSCRIPT 01_scripts/R/hi_by_group_3groups.R <hybrid_index_table>
#   <hybrid_index_table> = the per-individual table used for the triangle plot
#   (one row per fish, with a 'bam' column and a hybrid-index column)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Give the per-individual hybrid index table as argument")
OUT <- "04_pca/rankin_groups_check"

g  <- fread(file.path(OUT, "rankin_pca_3groups.tsv"))
hi <- fread(args[1])
cat("Columns of the hybrid index table:", paste(names(hi), collapse = ", "), "\n")

# sample key: BAM basename without extension (same rule as the PCA scripts)
idc <- if ("bam" %in% names(hi)) "bam" else names(hi)[1]
hi[, key := sub("\\..*$", "", basename(get(idc)))]
# hybrid index column: first column whose name looks like a hybrid index
hc <- grep("^(hi|hybrid|hybrid_?index|h_index)$|hybrid", names(hi), ignore.case = TRUE, value = TRUE)[1]
if (is.na(hc)) stop("No hybrid index column found; edit 'hc' in the script")
cat("Using", idc, "as sample ID and", hc, "as hybrid index\n")

d <- merge(g, hi[, .(key, HI = get(hc))], by.x = "sample", by.y = "key")
cat("Matched", nrow(d), "of", nrow(g), "Rankin fish\n")
d[, group := factor(group, levels = c("Main group", "PC1 cluster", "PC2 tail"))]

cat("\n--- Hybrid index by group ---\n")
print(d[, .(n = .N, mean = round(mean(HI), 3), median = round(median(HI), 3),
            sd = round(sd(HI), 3)), by = group][order(group)])

cat("\n--- Each small group vs the main group (Wilcoxon) ---\n")
for (k in c("PC1 cluster", "PC2 tail")) {
    w <- wilcox.test(d[group == k, HI], d[group == "Main group", HI])
    cat(sprintf("%-12s vs Main group: W = %.0f, p = %.3g\n", k, w$statistic, w$p.value))
}
cat("\n--- Kruskal-Wallis, three groups ---\n")
print(kruskal.test(HI ~ group, data = d))

fwrite(d, file.path(OUT, "rankin_3groups_hybrid_index.tsv"), sep = "\t")
cat("Written", file.path(OUT, "rankin_3groups_hybrid_index.tsv"), "\n")
