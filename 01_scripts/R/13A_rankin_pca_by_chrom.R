# =============================================================================
# 13A_rankin_pca_by_chrom.R — Who are the two groups in the Rankin Inlet PCA?
#
# Mode "groups"  : PC scores from the Rankin covariance matrix, k-means groups,
#                  composition by population, depth (mapped reads) and F by group.
# Mode "chroms"  : compares the genome-wide PC1/PC2 with the PC1/PC2 of each
#                  per-chromosome PCAngsd run. One chromosome with |r| close to 1
#                  and the others near 0  -> inversion / haploblock.
#                  Moderate |r| on every chromosome -> genome-wide signal
#                  (relatives, another stock, or a technical batch).
#
# Usage:
#   Rscript 13A_rankin_pca_by_chrom.R groups <cov> <bamlist> <out_dir> [mapped_reads.tsv] [F_individual.tsv]
#   Rscript 13A_rankin_pca_by_chrom.R chroms <cov> <chrom_dir> <out_dir>
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

pcs_from_cov <- function(f, n = 4) {
    C <- as.matrix(read.table(f))
    e <- eigen(C, symmetric = TRUE)
    pc <- e$vectors[, 1:n, drop = FALSE]
    colnames(pc) <- paste0("PC", 1:n)
    list(pc = pc, var = round(100 * e$values[1:n] / sum(e$values), 2))
}
sample_id <- function(bam) sub("\\..*$", "", basename(bam))
pop_of    <- function(id)  sub("s_.*$", "", id)

if (mode == "groups") {
    cov_f <- args[2]; bam_f <- args[3]; out <- args[4]
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    bams <- readLines(bam_f)
    g <- pcs_from_cov(cov_f)
    d <- data.frame(sample = sample_id(bams), pop = pop_of(sample_id(bams)), g$pc)
    cat("Variance explained (%):", g$var, "\n")

    # 3 groups on PC1-PC2 (main cloud, PC1 group, PC2 tail); the smallest
    # groups are the ones of interest
    set.seed(1)
    km <- kmeans(scale(d[, c("PC1", "PC2")]), centers = 3, nstart = 50)
    lab <- order(table(km$cluster), decreasing = TRUE)
    d$group <- factor(match(km$cluster, lab), labels = c("main", "groupB", "groupC"))

    cat("\n--- Group composition by population ---\n")
    print(table(d$group, d$pop))

    if (length(args) >= 5 && file.exists(args[5])) {
        m <- read.delim(args[5], header = FALSE, col.names = c("sample", "mapped_reads"))
        d <- merge(d, m, by = "sample", all.x = TRUE)
        cat("\n--- Mapped reads (millions) by group: if one group differs clearly, suspect a technical batch ---\n")
        print(aggregate(mapped_reads / 1e6 ~ group, d, function(x) round(c(mean = mean(x), sd = sd(x)), 1)))
    }
    if (length(args) >= 6 && file.exists(args[6])) {
        f <- read.delim(args[6])
        cat("\nF_individual.tsv columns:", paste(names(f), collapse = ", "), "\n")
        idc <- names(f)[1]; fc <- names(f)[ncol(f)]
        f$sample <- sample_id(f[[idc]])
        d <- merge(d, f[, c("sample", fc)], by = "sample", all.x = TRUE)
        cat("--- Inbreeding F by group ---\n")
        print(aggregate(d[[fc]] ~ d$group, FUN = function(x) round(mean(x, na.rm = TRUE), 3)))
    }
    write.table(d, file.path(out, "rankin_pca_groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    if (requireNamespace("ggplot2", quietly = TRUE)) {
        library(ggplot2)
        p <- ggplot(d, aes(PC1, PC2, colour = pop, shape = group)) + geom_point(size = 2) +
            labs(x = sprintf("PC1 (%.2f%%)", g$var[1]), y = sprintf("PC2 (%.2f%%)", g$var[2])) +
            theme_bw()
        ggsave(file.path(out, "rankin_pca_groups.png"), p, width = 6, height = 4.5, dpi = 200)
    }
    cat("\nWrote", file.path(out, "rankin_pca_groups.tsv"), "\n")

} else if (mode == "chroms") {
    cov_f <- args[2]; cdir <- args[3]; out <- args[4]
    G <- pcs_from_cov(cov_f, 2)$pc
    fs <- list.files(cdir, pattern = "\\.cov$", full.names = TRUE)
    res <- do.call(rbind, lapply(fs, function(f) {
        L <- pcs_from_cov(f, 2)
        r <- abs(cor(G, L$pc))       # 2 x 2: genome-wide PC vs chromosome PC
        data.frame(chrom = sub("\\.cov$", "", basename(f)),
                   n_snps = as.integer(tryCatch(readLines(sub("\\.cov$", ".nsnps", f)), error = function(e) NA)),
                   r_PC1 = round(max(r[1, ]), 3), r_PC2 = round(max(r[2, ]), 3),
                   var_PC1 = L$var[1])
    }))
    res <- res[order(-pmax(res$r_PC1, res$r_PC2)), ]
    write.table(res, file.path(out, "rankin_pca_by_chrom.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("\n--- |r| between genome-wide PCs and each chromosome's PCs (top 10) ---\n")
    print(head(res, 10), row.names = FALSE)
    cat("\nMedian over chromosomes: r_PC1 =", median(res$r_PC1), " r_PC2 =", median(res$r_PC2), "\n")
    top <- res[1, ]; rest <- res[-1, ]
    if (max(top$r_PC1, top$r_PC2) > 0.7 && max(top$r_PC1, top$r_PC2) > 3 * median(pmax(rest$r_PC1, rest$r_PC2))) {
        cat("=> Signal concentrated on", top$chrom, ": consistent with an inversion / haploblock.\n")
    } else {
        cat("=> Signal spread over many chromosomes: genome-wide (relatives, another stock or batch).\n")
    }
}
