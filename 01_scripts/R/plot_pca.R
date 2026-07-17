# =============================================================================
# plot_pca.R — Publication-quality PCA figures
# Handles saal* prefix in BAM names
# Produces: composite + individual global, rankin, naujaat figures
# Usage: Rscript plot_pca.R <pca_dir> <info_file> <fig_dir>
# =============================================================================

argv <- commandArgs(TRUE)
pca_dir  <- argv[1]
info_file <- argv[2]
fig_dir  <- argv[3]

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
})

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ─── Population color palette (colorblind-friendly) ─────────────────────────
pop_colors <- c(
    AKL = "#443983", AUL = "#31B8BD", CRB = "#5DC863",
    DIA = "#3B528B", MEL = "#21918C",
    HOR = "#E69F00",
    ITI = "#2CA02C", KGJ = "#7CAE00", NOP = "#F0E442",
    PAM = "#CD2626", SUP = "#F5A623", TIN = "#E68000", WHI = "#E8601C"
)

region_shapes <- c(Rankin = 16, Baker = 17, Naujaat = 15)

pop_region <- c(
    AKL = "Rankin", AUL = "Rankin", CRB = "Rankin", DIA = "Rankin", MEL = "Rankin",
    HOR = "Baker",
    ITI = "Naujaat", KGJ = "Naujaat", NOP = "Naujaat", PAM = "Naujaat",
    SUP = "Naujaat", TIN = "Naujaat", WHI = "Naujaat"
)

# ─── Extract population code from BAM path ───────────────────────────────────
extract_pop <- function(bam_path) {
    # Get filename without path and extension
    fname <- gsub(".*/", "", bam_path)
    fname <- gsub("\\.sorted\\.bam$|\\.bam$", "", fname)
    # Remove saal prefix if present
    fname <- gsub("^saal", "", fname)
    # Extract population code (uppercase letters before s_ or _ or digits)
    pop <- gsub("[s_].*$", "", fname)
    # Handle cases like TIN_25-20 (no s_ separator)
    pop <- gsub("_.*$", "", pop)
    return(pop)
}

# ─── Helper: load PCA + make plot ────────────────────────────────────────────
make_pca_plot <- function(prefix, label, pc_x = 1, pc_y = 2, point_size = 2) {
    pca_file <- paste0(prefix, ".cov.pca")
    eig_file <- paste0(prefix, ".cov.eig")

    if (!file.exists(pca_file)) {
        warning("Missing: ", pca_file)
        return(NULL)
    }

    pca <- read.table(pca_file, header = TRUE)
    eig <- scan(eig_file, quiet = TRUE)

    pca$bam <- rownames(pca)
    pca$pop <- sapply(pca$bam, extract_pop)
    pca$region <- pop_region[pca$pop]

    # Check for unmatched populations
    unmatched <- unique(pca$pop[is.na(pca$region)])
    if (length(unmatched) > 0) {
        cat("WARNING: Unmatched populations:", paste(unmatched, collapse = ", "), "\n")
    }

    pcx <- paste0("PC", pc_x)
    pcy <- paste0("PC", pc_y)
    varx <- eig[pc_x]
    vary <- eig[pc_y]

    p <- ggplot(pca, aes(x = .data[[pcx]], y = .data[[pcy]],
                    color = pop, shape = region)) +
        geom_point(size = point_size, alpha = 0.8) +
        scale_color_manual(values = pop_colors, name = "Population") +
        scale_shape_manual(values = region_shapes, name = "Region") +
        labs(
            x = paste0(pcx, " (", varx, "%)"),
            y = paste0(pcy, " (", vary, "%)"),
            title = label
        ) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            legend.position = "right",
            plot.title = element_text(face = "bold", size = 12)
        )

    return(p)
}

# ─── Build file prefix ──────────────────────────────────────────────────────
suffix <- "maf0.05_pctind0.50_maxdepth8"

global_prefix <- file.path(pca_dir, paste0("global_", suffix, "_prunednosex"))
rankin_prefix <- file.path(pca_dir, paste0("rankin_", suffix, "_prunednosex"))
naujaat_prefix <- file.path(pca_dir, paste0("naujaat_", suffix, "_prunednosex"))

# ─── Composite figure (4 panels) ────────────────────────────────────────────
cat("Generating composite PCA figure...\n")

p_global_12 <- make_pca_plot(global_prefix, "A) All populations", 1, 2)
p_global_34 <- make_pca_plot(global_prefix, "B) All populations", 3, 4)
p_rankin    <- make_pca_plot(rankin_prefix, "C) Rankin Inlet", 1, 2)
p_naujaat   <- make_pca_plot(naujaat_prefix, "D) Naujaat", 1, 2)

if (!is.null(p_global_12) & !is.null(p_rankin)) {
    combined <- (p_global_12 | p_global_34) / (p_rankin | p_naujaat) +
        plot_layout(guides = "collect") &
        theme(legend.position = "right")

    ggsave(file.path(fig_dir, "Fig_PCA_composite.pdf"),
           combined, width = 14, height = 10, dpi = 300)
    ggsave(file.path(fig_dir, "Fig_PCA_composite.png"),
           combined, width = 14, height = 10, dpi = 300)
    cat("  Composite saved.\n")
}

# ─── Individual figures ──────────────────────────────────────────────────────
cat("Generating individual PCA figures...\n")

# Global PC1-2
p <- make_pca_plot(global_prefix, "All populations — PC1 vs PC2", 1, 2, point_size = 2.5)
if (!is.null(p)) {
    ggsave(file.path(fig_dir, "Fig_PCA_global_PC12.pdf"), p, width = 9, height = 7, dpi = 300)
    ggsave(file.path(fig_dir, "Fig_PCA_global_PC12.png"), p, width = 9, height = 7, dpi = 300)
    cat("  Global PC1-2 saved.\n")
}

# Global PC3-4
p <- make_pca_plot(global_prefix, "All populations — PC3 vs PC4", 3, 4, point_size = 2.5)
if (!is.null(p)) {
    ggsave(file.path(fig_dir, "Fig_PCA_global_PC34.pdf"), p, width = 9, height = 7, dpi = 300)
    ggsave(file.path(fig_dir, "Fig_PCA_global_PC34.png"), p, width = 9, height = 7, dpi = 300)
    cat("  Global PC3-4 saved.\n")
}

# Rankin PC1-2
p <- make_pca_plot(rankin_prefix, "Rankin Inlet — PC1 vs PC2", 1, 2, point_size = 3)
if (!is.null(p)) {
    ggsave(file.path(fig_dir, "Fig_PCA_rankin_PC12.pdf"), p, width = 9, height = 7, dpi = 300)
    ggsave(file.path(fig_dir, "Fig_PCA_rankin_PC12.png"), p, width = 9, height = 7, dpi = 300)
    cat("  Rankin PC1-2 saved.\n")
}

# Rankin PC3-4
p <- make_pca_plot(rankin_prefix, "Rankin Inlet — PC3 vs PC4", 3, 4, point_size = 3)
if (!is.null(p)) {
    ggsave(file.path(fig_dir, "Fig_PCA_rankin_PC34.pdf"), p, width = 9, height = 7, dpi = 300)
    ggsave(file.path(fig_dir, "Fig_PCA_rankin_PC34.png"), p, width = 9, height = 7, dpi = 300)
    cat("  Rankin PC3-4 saved.\n")
}

# Naujaat PC1-2
p <- make_pca_plot(naujaat_prefix, "Naujaat — PC1 vs PC2", 1, 2, point_size = 3)
if (!is.null(p)) {
    ggsave(file.path(fig_dir, "Fig_PCA_naujaat_PC12.pdf"), p, width = 9, height = 7, dpi = 300)
    ggsave(file.path(fig_dir, "Fig_PCA_naujaat_PC12.png"), p, width = 9, height = 7, dpi = 300)
    cat("  Naujaat PC1-2 saved.\n")
}

# Naujaat PC3-4
p <- make_pca_plot(naujaat_prefix, "Naujaat — PC3 vs PC4", 3, 4, point_size = 3)
if (!is.null(p)) {
    ggsave(file.path(fig_dir, "Fig_PCA_naujaat_PC34.pdf"), p, width = 9, height = 7, dpi = 300)
    ggsave(file.path(fig_dir, "Fig_PCA_naujaat_PC34.png"), p, width = 9, height = 7, dpi = 300)
    cat("  Naujaat PC3-4 saved.\n")
}


# ─── Versions with 95% ellipses ─────────────────────────────────────────────
cat("Generating PCA figures with ellipses...\n")

add_ellipse <- function(p) {
    p + stat_ellipse(aes(group = pop, fill = pop), level = 0.95, geom = "polygon", alpha = 0.15,
                     linewidth = 0.5, show.legend = FALSE) + scale_fill_manual(values = pop_colors, guide = "none")
}

for (info in list(
    list(prefix = global_prefix, name = "global_PC12", label = "All populations - PC1 vs PC2", pcx = 1, pcy = 2),
    list(prefix = global_prefix, name = "global_PC34", label = "All populations - PC3 vs PC4", pcx = 3, pcy = 4),
    list(prefix = rankin_prefix, name = "rankin_PC12", label = "Rankin Inlet - PC1 vs PC2", pcx = 1, pcy = 2),
    list(prefix = rankin_prefix, name = "rankin_PC34", label = "Rankin Inlet - PC3 vs PC4", pcx = 3, pcy = 4),
    list(prefix = naujaat_prefix, name = "naujaat_PC12", label = "Naujaat - PC1 vs PC2", pcx = 1, pcy = 2),
    list(prefix = naujaat_prefix, name = "naujaat_PC34", label = "Naujaat - PC3 vs PC4", pcx = 3, pcy = 4)
)) {
    p <- make_pca_plot(info$prefix, info$label, info$pcx, info$pcy, point_size = 2.5)
    if (!is.null(p)) {
        pe <- add_ellipse(p)
        ggsave(file.path(fig_dir, paste0("Fig_PCA_", info$name, "_ellipse.pdf")),
               pe, width = 9, height = 7, dpi = 300)
        ggsave(file.path(fig_dir, paste0("Fig_PCA_", info$name, "_ellipse.png")),
               pe, width = 9, height = 7, dpi = 300)
    }
}

cat("Ellipse figures saved.\n")
EOF



cat("PCA figures saved to:", fig_dir, "\n")