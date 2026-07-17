# =============================================================================
# sum_sfs.R — Sum multi-line SFS (from realSFS -nSites) into single line
# Usage: Rscript sum_sfs.R <sfs_file>
# =============================================================================

argv <- commandArgs(TRUE)
sfs_file <- argv[1]

sfs <- read.table(sfs_file)
sfs_sum <- colSums(sfs)
write.table(rbind(sfs_sum), paste0(sfs_file, ".dsfs"),
            quote = FALSE, col.names = FALSE, row.names = FALSE)

cat("SFS summed:", length(sfs_sum), "categories\n")
