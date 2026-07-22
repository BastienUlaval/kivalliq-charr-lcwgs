# =============================================================================
# make_sites_list.R -- Extract chr/pos/major/minor from ANGSD .mafs file
# Usage: Rscript make_sites_list.R <infile.mafs> <outfile_sites>
# =============================================================================

argv <- commandArgs(TRUE)
INFILE <- argv[1]
OUTFILE <- argv[2]

maf <- read.table(INFILE, header = TRUE)
sites <- maf[, 1:4]
sites <- sites[order(sites$chromo), ]
write.table(sites, OUTFILE, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

cat("Sites extracted:", nrow(sites), "\n")
