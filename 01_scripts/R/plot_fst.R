# =============================================================================
# plot_fst.R — FST heatmap (red gradient) + regional side-by-side + tables
# Usage: Rscript plot_fst.R <fst_results.tsv> <pop.txt> <fig_dir> <table_dir>
# =============================================================================
argv <- commandArgs(TRUE)
fst_file  <- argv[1]
pop_file  <- argv[2]
fig_dir   <- argv[3]
table_dir <- argv[4]

suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
    library(dplyr)
    library(patchwork)
})

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# ─── Read data ───────────────────────────────────────────────────────────────
fst <- read.delim(fst_file, header = TRUE)
pops <- scan(pop_file, what = "character", quiet = TRUE)
n <- length(pops)

# ─── Build FST matrix ───────────────────────────────────────────────────────
mat <- matrix(0, n, n, dimnames = list(pops, pops))
for (i in seq_len(nrow(fst))) {
    p1 <- fst$Pop1[i]; p2 <- fst$Pop2[i]; val <- fst$Fst_weighted[i]
    if (p1 %in% pops & p2 %in% pops) {
        mat[p1, p2] <- val; mat[p2, p1] <- val
    }
}

# Pop ordering by region
pop_order <- c("AKL","AUL","CRB","DIA","MEL","HOR","ITI","KGJ","NOP","PAM","SUP","TIN","WHI")
pop_order <- intersect(pop_order, pops)
mat <- mat[pop_order, pop_order]

# ─── Region definitions ─────────────────────────────────────────────────────
rankin_pops  <- intersect(c("AKL","AUL","CRB","DIA","MEL"), pop_order)
naujaat_pops <- intersect(c("ITI","KGJ","NOP","PAM","SUP","TIN","WHI"), pop_order)

# ─── Helper: build heatmap from a FST matrix ────────────────────────────────
make_heatmap <- function(fst_mat, pop_ord, title_text, fst_limits) {
    melted <- melt(fst_mat)
    colnames(melted) <- c("Pop1", "Pop2", "FST")
    melted$Pop1 <- factor(melted$Pop1, levels = pop_ord)
    melted$Pop2 <- factor(melted$Pop2, levels = rev(pop_ord))

    ggplot(melted, aes(Pop1, Pop2, fill = FST)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = ifelse(FST > 0, sprintf("%.3f", FST), "")),
                  size = 2.8) +
        scale_fill_gradient(low = "white", high = "darkred",
                            name = expression(F[ST]),
                            limits = fst_limits) +
        labs(title = title_text) +
        theme_minimal(base_size = 11) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(face = "bold"),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5)
        ) +
        coord_fixed()
}

# ─── 1. GLOBAL HEATMAP (all 13 pops) ────────────────────────────────────────
global_max <- max(mat[mat > 0])

p_global <- make_heatmap(mat, pop_order, "Global", c(0, global_max))

ggsave(file.path(fig_dir, "Fig_FST_heatmap.png"), p_global,
       width = 8, height = 7, dpi = 300)
cat("Global FST heatmap saved.\n")

# ─── 2. REGIONAL SIDE-BY-SIDE (Rankin | Naujaat, same scale) ────────────────
mat_rankin  <- mat[rankin_pops, rankin_pops]
mat_naujaat <- mat[naujaat_pops, naujaat_pops]

# Shared scale across both regional panels
regional_max <- max(c(mat_rankin[mat_rankin > 0], mat_naujaat[mat_naujaat > 0]))
regional_limits <- c(0, regional_max)

p_rankin <- make_heatmap(mat_rankin, rankin_pops, "Rankin Inlet", regional_limits)
p_naujaat <- make_heatmap(mat_naujaat, naujaat_pops, "Naujaat", regional_limits)

# Combine side by side with shared legend
p_regional <- p_rankin + p_naujaat +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")

ggsave(file.path(fig_dir, "Fig_FST_heatmap_regional.png"), p_regional,
       width = 14, height = 6, dpi = 300)
cat("Regional FST heatmap saved.\n")

# ─── 3. TABLES ──────────────────────────────────────────────────────────────
write.csv(mat, file.path(table_dir, "Table_FST_matrix.csv"))

fst$Fst_linearized <- fst$Fst_weighted / (1 - fst$Fst_weighted)
write.csv(fst, file.path(table_dir, "Table_FST_pairwise.csv"), row.names = FALSE)

cat("FST tables saved.\n")