# =============================================================================
# sum_2dsfs.R -- Sum multi-line 2D-SFS into single line
# Usage: Rscript sum_2dsfs.R <2dsfs_file>
# =============================================================================

argv <- commandArgs(TRUE)
sfs_file <- argv[1]

sfs <- read.table(sfs_file)
sfs_sum <- colSums(sfs)
write.table(rbind(sfs_sum), paste0(sfs_file, ".summed"),
            quote = FALSE, col.names = FALSE, row.names = FALSE)
