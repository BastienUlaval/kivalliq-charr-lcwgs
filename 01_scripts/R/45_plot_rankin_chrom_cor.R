# =============================================================================
# 45_plot_rankin_chrom_cor.R — Visualise the per-chromosome check of the
# second Rankin Inlet PCA group (output of 44_rankin_pca_groups.sh)
#
# Panel a: |r| between genome-wide PC1 and each chromosome's PCs (one bar per
#          chromosome). Inversion = one bar near 1, the others near 0.
#          Genome-wide signal = most bars high.
# Panel b: individual scores, genome-wide PC1 vs PC1 of one chromosome
#          (highest and lowest |r|), coloured by PCA group.
#
# Run from angsd_pipeline2/:
#   $RSCRIPT 01_scripts/R/45_plot_rankin_chrom_cor.R
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

OUT  <- "04_pca/rankin_groups_check"
CDIR <- file.path(OUT, "by_chrom")
COV  <- list.files("04_pca", pattern = "^rankin_.*_prunednosex\\.cov$", full.names = TRUE)[1]
BAM  <- "02_info/bamlists/rankin.bamlist"

pcs <- function(f, n = 2) {
    e <- eigen(as.matrix(read.table(f)), symmetric = TRUE)
    pc <- e$vectors[, 1:n, drop = FALSE]; colnames(pc) <- paste0("PC", 1:n); pc
}
G   <- pcs(COV)
# orient the axes so that the outlying fish have negative scores (the sign of a PC is arbitrary)
orient <- function(x) if (max(x) - median(x) > median(x) - min(x)) -x else x
G[, 1] <- orient(G[, 1]); G[, 2] <- orient(G[, 2])
ids <- sub("\\..*$", "", basename(readLines(BAM)))     # same order as the .cov rows
gt  <- fread(file.path(OUT, "rankin_pca_groups.tsv"))
m   <- match(ids, gt$sample)
# Three groups defined by position in the Rankin PCA (thresholds read on the figure)
TH_PC1 <- -0.10   # separate cluster on PC1
TH_PC2 <- -0.075  # tail of the main cloud on PC2
grp <- ifelse(G[, 1] < TH_PC1, "PC1 cluster", ifelse(G[, 2] < TH_PC2, "PC2 tail", "Main group"))
grp <- factor(grp, levels = c("Main group", "PC1 cluster", "PC2 tail"))
# population code: the 3-letter code found anywhere in the sample ID or in the pop column
CODES <- c("CRB", "DIA", "MEL", "AUL", "AKL")                       # geographic order
find_code <- function(x) { x <- toupper(x); hit <- sapply(CODES, grepl, x = x)
                           if (is.matrix(hit)) apply(hit, 1, function(h) if (any(h)) CODES[which(h)[1]] else NA)
                           else if (any(hit)) CODES[which(hit)[1]] else NA }
pop <- find_code(paste(ids, gt$pop[m]))
if (anyNA(pop)) { cat("Sample IDs without a population code (first 5):\n"); print(head(ids[is.na(pop)], 5)) }
pop <- factor(pop, levels = CODES)
# Population colours: taken from fig_common.R, the file used by
# 38_Fig3_regional_structure.R, so every fish has the colour it has in Fig. 3c
fc <- "01_scripts/R/fig_common.R"
POP_COLORS <- NULL
if (file.exists(fc)) {
    POP_COLORS <- tryCatch(local({ source(fc, local = TRUE); POP_COLORS }), error = function(e) NULL)
}
if (is.null(POP_COLORS) || !exists("POP_COLORS")) {   # same values as in fig_common.R
    POP_COLORS <- c("AUL" = "#E69F00", "MEL" = "#D55E00", "DIA" = "#CC6677",
                    "AKL" = "#882255", "CRB" = "#AA4499")
}
POP_COLORS <- POP_COLORS[levels(pop)]
SHAPES <- c("Main group" = 16, "PC1 cluster" = 17, "PC2 tail" = 15)

# ---- panel a: |r| per chromosome --------------------------------------------
res <- fread(file.path(OUT, "rankin_pca_by_chrom.tsv"))
res[, chrom := factor(chrom, levels = chrom[order(r_PC1)])]
med <- median(res$r_PC1)
pa <- ggplot(res, aes(r_PC1, chrom)) +
    geom_col(fill = "grey60", width = 0.7) +
    geom_vline(xintercept = med, linetype = 2, colour = "firebrick") +
    annotate("text", x = med, y = nrow(res) + 0.6, label = sprintf("median = %.2f", med),
             hjust = 1.05, vjust = 0, size = 3, colour = "firebrick") +
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "|r| with genome-wide PC1", y = NULL, tag = "a") +
    theme_bw(base_size = 9) + theme(panel.grid.major.y = element_blank())

# ---- panel c: scatter for the highest and lowest |r| chromosomes -------------
pick <- res[c(which.max(r_PC1), which.min(r_PC1))]
sc <- rbindlist(lapply(seq_len(nrow(pick)), function(i) {
    L <- pcs(file.path(CDIR, paste0(pick$chrom[i], ".cov")))
    k <- which.max(abs(cor(G[, 1], L)))                 # chromosome PC that matches best
    x <- L[, k] * sign(cor(G[, 1], L[, k]))             # flip sign so r is positive
    data.table(chrom = sprintf("%s (|r| = %.2f)", pick$chrom[i], abs(cor(G[, 1], x))),
               x = x, y = G[, 1], group = grp, pop = pop)
}))
pb <- ggplot(sc, aes(x, y, colour = pop, shape = group)) +
    geom_point(size = 1.6, alpha = 0.85) +
    facet_wrap(~chrom, ncol = 1, scales = "free_x") +
    scale_colour_manual(values = POP_COLORS, name = "Population") +
    scale_shape_manual(values = SHAPES, name = NULL) +
    labs(x = "Chromosome PC score", y = "Genome-wide PC1 score", tag = "c") +
    guides(colour = guide_legend(override.aes = list(size = 2.5, shape = 16), nrow = 1),
           shape  = guide_legend(override.aes = list(size = 2.5))) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.box = "vertical")

# ---- panel b: Rankin regional PCA (as Fig. 3c), shape = PCA group ----------
pcd <- data.table(PC1 = G[, 1], PC2 = G[, 2], pop = pop, group = grp)
pc <- ggplot(pcd, aes(PC1, PC2, colour = pop, shape = group)) +
    geom_hline(yintercept = TH_PC2, linetype = 3, colour = "grey60") +
    geom_vline(xintercept = TH_PC1, linetype = 3, colour = "grey60") +
    geom_point(size = 1.6, alpha = 0.85) +
    scale_colour_manual(values = POP_COLORS, name = "Population") +
    scale_shape_manual(values = SHAPES, name = NULL) +
    labs(x = "PC1", y = "PC2", tag = "b", title = "Rankin Inlet PCA") +
    guides(colour = "none", shape = "none") +
    theme_bw(base_size = 9) + theme(plot.title = element_text(size = 8, hjust = 0.5))
cat("\n--- PCA group by population ---\n"); print(addmargins(table(pcd$group, pcd$pop)))
cat("\n--- New groups vs k-means groups of script 44 ---\n"); print(table(pcd$group, gt$group[m]))
fwrite(data.table(sample = ids, pop = pop, PC1 = G[, 1], PC2 = G[, 2], group = grp),
       file.path(OUT, "rankin_pca_3groups.tsv"), sep = "\t")

if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    ggsave(file.path(OUT, "Fig_rankin_chrom_cor.png"), ((pa | (pc / pb + plot_layout(heights = c(1, 2)))) + plot_layout(widths = c(1, 1), guides = "collect")) &
               theme(legend.position = "bottom", legend.box = "vertical",
                     legend.text = element_text(size = 7), legend.title = element_text(size = 7)),
           width = 18, height = 20, units = "cm", dpi = 300)
} else {
    ggsave(file.path(OUT, "Fig_rankin_chrom_cor_a.png"), pa, width = 9, height = 14, units = "cm", dpi = 300)
    ggsave(file.path(OUT, "Fig_rankin_chrom_cor_b.png"), pb, width = 9, height = 14, units = "cm", dpi = 300)
    ggsave(file.path(OUT, "Fig_rankin_chrom_cor_c.png"), pc, width = 9, height = 8, units = "cm", dpi = 300)
}
cat("Written in", OUT, "\n")
