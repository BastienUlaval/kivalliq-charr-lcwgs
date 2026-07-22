# =============================================================================
# classify_paralog.R -- Separate canonical from deviant (paralogous) sites
# Usage: Rscript classify_paralog.R <ngsP_file> <deviant_out> <keep_out> <threshold>
# =============================================================================
argv <- commandArgs(TRUE)
if (length(argv) < 4) {
    stop("Usage: Rscript classify_paralog.R <ngsP_file> <deviant_out> <keep_out> <threshold>")
}
ngsP_file    <- argv[1]
deviant_out  <- argv[2]
keep_out     <- argv[3]
threshold    <- as.numeric(argv[4])

if (is.na(threshold)) stop("Threshold (argv[4]) is not numeric: ", argv[4])

ngsP <- read.table(ngsP_file)
ngsP$pval     <- 0.5 * pchisq(ngsP$V5, df = 1, lower.tail = FALSE)
ngsP$pval.adj <- p.adjust(ngsP$pval, method = "BH")

canonical <- ngsP[ngsP$pval.adj > threshold, ]
deviant   <- ngsP[ngsP$pval.adj <= threshold, ]

write.table(canonical[, 1:2], keep_out,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(deviant[, 1:2], deviant_out,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

cat("Canonical (kept):", nrow(canonical),
    "| Deviant (excluded):", nrow(deviant),
    "(threshold:", threshold, ")\n")
