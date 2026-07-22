# =============================================================================
# deviant_masking.R -- Create BED masks around paralogous SNPs
# Usage: Rscript deviant_masking.R <mask_length> <ngsparalog_dir> <info_dir>
# =============================================================================

argv <- commandArgs(TRUE)
MASK_LENGTH <- as.numeric(argv[1])
NGSPARALOG_DIR <- argv[2]
INFO_DIR <- argv[3]

DIST <- MASK_LENGTH / 2

suppressPackageStartupMessages({
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos = "https://cran.r-project.org")
        BiocManager::install("GenomicRanges")
    }
    library(GenomicRanges)
})

cat("Mask length:", MASK_LENGTH, "bp (+/-", DIST, "bp)\n")

fai <- read.table(file.path(INFO_DIR, "genome.fasta.fai"))
colnames(fai) <- c("chr", "length", "offset", "linebases", "linewidth")

region_num <- scan(file.path(INFO_DIR, "regions_number.txt"),
                   what = "character", quiet = TRUE)

mask_dir <- file.path(INFO_DIR, "mask_by_chr")
dir.create(mask_dir, showWarnings = FALSE, recursive = TRUE)

total_dev <- 0
total_masked <- 0

for (NUM in region_num) {
    ngsP_file <- file.path(NGSPARALOG_DIR,
        paste0("all_maf0.05_pctind0.50_maxdepth8_chr", NUM, ".ngsparalog"))

    if (!file.exists(ngsP_file)) {
        cat("  Skip chr", NUM, ": file missing\n")
        next
    }

    LR <- read.table(ngsP_file, header = FALSE)
    if (nrow(LR) == 0) next

    LR$p.value <- p.adjust(
        0.5 * pchisq(LR$V5, df = 1, lower.tail = FALSE), method = "BH")
    LR$deviant <- LR$p.value < 0.001

    n_dev <- sum(LR$deviant)
    total_dev <- total_dev + n_dev

    if (n_dev == 0) {
        write.table(data.frame(), file.path(mask_dir, paste0("mask_deviant_chr", NUM, ".bed")),
                    col.names = FALSE, row.names = FALSE)
        next
    }

    dev <- data.frame(
        mid   = LR$V2[LR$deviant],
        chr   = LR$V1[1]
    )
    dev$start <- pmax(dev$mid - DIST, 1)
    chr_len <- fai$length[fai$chr == dev$chr[1]]
    dev$end <- pmin(dev$mid + DIST, chr_len)

    gr <- makeGRangesFromDataFrame(dev, seqnames.field = "chr",
                                    start.field = "start", end.field = "end")
    gr_merged <- reduce(gr)
    gr_filt <- gr_merged[width(gr_merged) >= MASK_LENGTH]
    if (length(gr_filt) == 0) gr_filt <- gr_merged

    df <- as.data.frame(gr_filt)
    masked_bp <- sum(df$width)
    total_masked <- total_masked + masked_bp

    mask <- df[, c("seqnames", "start", "end")]
    write.table(mask, file.path(mask_dir, paste0("mask_deviant_chr", NUM, ".bed")),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

    cat("  chr", NUM, ":", n_dev, "deviants,", masked_bp, "bp masked\n")
}

total_genome <- sum(fai$length)
pct <- round(total_masked * 100 / total_genome, 2)
cat("\n=== SUMMARY ===\n")
cat("Total deviants:", total_dev, "\n")
cat("Total masked:", total_masked, "bp (", pct, "% of genome)\n")