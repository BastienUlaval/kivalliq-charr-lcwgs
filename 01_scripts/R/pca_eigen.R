# =============================================================================
# pca_eigen.R -- Eigendecomposition of PCAngsd covariance matrix
# Usage: Rscript pca_eigen.R <cov_matrix_file> <bam_filelist>
# =============================================================================

argv <- commandArgs(TRUE)
INPUT <- argv[1]
BAM   <- argv[2]

cov_mat <- as.matrix(read.table(INPUT, header = FALSE))
pca <- eigen(cov_mat)

# Scale eigenvectors by sqrt(eigenvalue)
pca_mat <- pca$vectors %*% diag(sqrt(pmax(pca$values, 0)))

# Column names
nPC <- ncol(pca_mat)
colnames(pca_mat) <- paste0("PC", seq_len(nPC))

# Row names from BAM list
bam_names <- read.table(BAM, header = FALSE)
rownames(pca_mat) <- bam_names$V1

# Variance explained
var_total <- sum(pca$values[pca$values >= 0])
var_explained <- round(pca$values[1:min(10, nPC)] * 100 / var_total, 2)

# Save
write.table(pca_mat[, 1:min(10, nPC)], paste0(INPUT, ".pca"), quote = FALSE)
write.table(var_explained, paste0(INPUT, ".eig"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

cat("PCA done:", nrow(pca_mat), "samples,", nPC, "PCs\n")
cat("Variance explained (PC1-4):", var_explained[1:min(4, length(var_explained))], "\n")
