#!/usr/bin/env Rscript
# =============================================================================
# 19_plot_admix_global_K.R — global NGSadmix barplots for one or several K,
# stacked for comparison (same style as Fig. 2b). Best replicate for each K.
#
# Usage (from the pipeline root, no sbatch needed, < 1 min):
#   $RSCRIPT 01_scripts/R/19_plot_admix_global_K.R $PWD 7 8
# Output: 99_figures/heredity/Admix_global_K7_K8.{pdf,png,tiff}
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 19_plot_admix_global_K.R <base_dir> <K> [<K> ...]")
BASE_DIR <- args[1]
Ks <- as.integer(args[-1])
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "fig_common.R"))

ADMIX_DIR <- file.path(BASE_DIR, "05_admixture", "global")
BAMLIST   <- file.path(BASE_DIR, "02_info", "bam.filelist")
FIG_DIR   <- file.path(BASE_DIR, "99_figures", "heredity")

plots <- lapply(Ks, function(K) {
    plot_admix(best_qopt(ADMIX_DIR, "global", K), BAMLIST) +
        labs(subtitle = sprintf("K = %d", K)) +
        theme(plot.subtitle = element_text(size = 8, hjust = 0))
})
fig <- wrap_plots(plots, ncol = 1)
save_fig(fig, FIG_DIR, paste0("Admix_global_K", paste(Ks, collapse = "_K")),
         WIDTH_2COL, WIDTH_2COL * 0.33 * length(Ks))
