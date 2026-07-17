# =============================================================================
# subsample_equal_n.R — Downsample populations to smallest N for FST
# Usage: Rscript subsample_equal_n.R <bamlist_dir> <output_dir> <pop_file>
# =============================================================================

argv <- commandArgs(TRUE)
bamlist_dir <- argv[1]
output_dir  <- argv[2]
pop_file    <- argv[3]

pops <- scan(pop_file, what = "character", quiet = TRUE)

# Find minimum N
sizes <- sapply(pops, function(p) {
    f <- file.path(bamlist_dir, paste0(p, ".bamlist"))
    if (file.exists(f)) length(readLines(f)) else 0
})

cat("Population sizes:\n")
print(data.frame(Pop = pops, N = sizes))

min_n <- min(sizes[sizes > 0])
cat("Subsampling to:", min_n, "individuals per population\n")

set.seed(42)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (p in pops) {
    f <- file.path(bamlist_dir, paste0(p, ".bamlist"))
    if (!file.exists(f)) next
    bams <- readLines(f)
    idx <- sample(seq_along(bams), min(min_n, length(bams)))
    out_f <- file.path(output_dir, paste0(p, "_subset.bamlist"))
    writeLines(bams[sort(idx)], out_f)
    cat("  ", p, ":", length(idx), "samples\n")
}
