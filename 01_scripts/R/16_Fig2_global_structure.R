#!/usr/bin/env Rscript
# =============================================================================
# 16_Fig2_global_structure.R — Figure 2 and Figure S1 (Heredity)
#
#   Fig. 2   (a) global PCA, PC1 vs PC2, 95% ellipses per zone
#            (b) global NGSadmix at K_MAIN (default 8), best replicate
#   Fig. S1  global NGSadmix, the two K values stacked, smaller K on top
#            (default K = 7 top, K = 8 bottom)
#
# Colours are matched across K: each cluster at K_SUPP takes the colour of the
# K_MAIN cluster it correlates with best. So Fig. S1 uses the same colours as
# Fig. 2b, and at K = 7 PAM shows as a mix of the colours of other clusters.
#
# Usage (from the pipeline root, no sbatch needed, < 1 min):
#   $RSCRIPT 01_scripts/R/16_Fig2_global_structure.R $PWD          # main K = 8, supp K = 7
#   $RSCRIPT 01_scripts/R/16_Fig2_global_structure.R $PWD 8 7      # explicit
# Output: 99_figures/heredity/Fig2_global_structure.{pdf,png,tiff}
#         99_figures/heredity/FigS1_admix_global_K7_K8.{pdf,png,tiff}
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 16_Fig2_global_structure.R <base_dir> [K_main] [K_supp]")
BASE_DIR <- args[1]
K_MAIN <- if (length(args) >= 2) as.integer(args[2]) else 8L
K_SUPP <- if (length(args) >= 3) as.integer(args[3]) else 7L
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "fig_common.R"))

SUFFIX    <- "maf0.05_pctind0.50_maxdepth8_prunednosex"
PCA_DIR   <- file.path(BASE_DIR, "04_pca")
ADMIX_DIR <- file.path(BASE_DIR, "05_admixture", "global")
BAMLIST   <- file.path(BASE_DIR, "02_info", "bam.filelist")
FIG_DIR   <- file.path(BASE_DIR, "99_figures", "heredity")

# --- Fig. 2 -------------------------------------------------------------------
pca <- load_pca(file.path(PCA_DIR, paste0("global_", SUFFIX)), BAMLIST)
p_a <- plot_pca(pca, ellipse = "region")

q_main <- best_qopt(ADMIX_DIR, "global", K_MAIN)
p_b <- plot_admix(q_main, BAMLIST)

fig2 <- (p_a / p_b) + plot_layout(heights = c(1.6, 1)) +
    plot_annotation(tag_levels = "a")
save_fig(fig2, FIG_DIR, "Fig2_global_structure", WIDTH_2COL, WIDTH_2COL * 0.85)

# --- Fig. S1: both K, smaller K on top, colours matched to Fig. 2b ---------
q_supp <- best_qopt(ADMIX_DIR, "global", K_SUPP)
ref    <- admix_colour_ref(q_main, BAMLIST)
lab <- function(K) labs(subtitle = sprintf("K = %d", K))
p_main <- p_b + lab(K_MAIN)
p_supp <- plot_admix(q_supp, BAMLIST, ref = ref) + lab(K_SUPP)
s4 <- (if (K_SUPP < K_MAIN) p_supp / p_main else p_main / p_supp) &
    theme(plot.subtitle = element_text(size = 8, hjust = 0))
Ks <- sort(c(K_MAIN, K_SUPP))
save_fig(s4, FIG_DIR, sprintf("FigS1_admix_global_K%d_K%d", Ks[1], Ks[2]),
         WIDTH_2COL, WIDTH_2COL * 0.66)
