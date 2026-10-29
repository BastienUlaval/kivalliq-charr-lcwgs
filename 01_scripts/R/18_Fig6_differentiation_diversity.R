#!/usr/bin/env Rscript
# =============================================================================
# 18_Fig6_differentiation_diversity.R — Figure 6 (Heredity)
#   (a) Sampling design (zones)
#   (b) Pairwise weighted FST within and between zones
#   (c) Nucleotide diversity (pi) per population
#   (d) Isolation by distance: global, Naujaat, Rankin Inlet
#
# Earlier versions of this figure also showed D-statistics and management units.
# Removed: old panel D (D-statistics, redundant with Fig. 5), old panel F
# (management-unit cards), the figure title, and the partial-Mantel block.
# Four panels -> counts as one display item at Heredity.
# Regression line drawn only where the Mantel test is significant (p < 0.05).
#
# Inputs:
#   <base_dir>/99_tables/main/Table_Diversity.csv   Pop, Region, pi
#   <base_dir>/14_ibd/Table_IBD_pairwise.csv        Pop1, Pop2, Fst_weighted,
#                                                   Fst_lin, Hydro_km, Comparison
#   <base_dir>/14_ibd/Table_IBD_Mantel.csv          Group, Distance, Mantel_r, Mantel_p
#
# Usage:
#   Rscript 18_Fig6_differentiation_diversity.R <base_dir> <shp_dir> [<r_lib_path>]
#   Output: <base_dir>/99_figures/heredity/Fig6_differentiation_diversity.{pdf,png,tiff}
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 18_Fig6_differentiation_diversity.R <base_dir> <shp_dir> [<r_lib_path>]")
BASE_DIR <- args[1]; SHP_DIR <- args[2]
if (length(args) >= 3) .libPaths(args[3])
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "fig_common.R"))
suppressPackageStartupMessages(library(sf))

TABLE_DIR <- file.path(BASE_DIR, "99_tables", "main")
IBD_DIR   <- file.path(BASE_DIR, "14_ibd")
FIG_DIR   <- file.path(BASE_DIR, "99_figures", "heredity")
ALPHA     <- 0.05

norm_region <- function(x) {
    x <- as.character(x)
    x[grepl("^Rankin", x)] <- "Rankin Inlet"
    x[grepl("^Baker|^HOR$", x)] <- "Baker Lake"
    factor(x, levels = REGION_ORDER)
}

# --- (a) sampling design -----------------------------------------------------
XLIM <- c(-99, -82); YLIM <- c(61.5, 68.5)
bm <- load_basemap(SHP_DIR, XLIM, YLIM)
p_a <- ggplot() +
    geom_sf(data = bm$land,  fill = "grey90", colour = "grey60", linewidth = 0.12) +
    geom_sf(data = bm$lakes, fill = "white",  colour = "grey70", linewidth = 0.08) +
    geom_point(data = SITES, aes(lon, lat, colour = region, shape = region), size = 1.4) +
    # Rankin Inlet label left-aligned east of AKL/AUL (-91.3) so it does not cover the sites
    annotate("text", x = c(-95.6, -90.6, -84.6), y = c(64.25, 62.70, 67.15),
             hjust = c(0.5, 0, 0.5),
             label = REGION_ORDER, colour = REGION_COLORS[REGION_ORDER],
             fontface = "bold", size = 2.5) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_shape_manual(values = REGION_SHAPES, guide = "none") +
    coord_sf(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_hered() +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.3),
          axis.text = element_text(size = 5.5))
p_a <- add_scalebar(p_a, location = "bl")

# --- (b) pairwise FST --------------------------------------------------------
pw <- read.csv(file.path(IBD_DIR, "Table_IBD_pairwise.csv"))
pw$cmp <- factor(ifelse(pw$Comparison == "Rankin Inlet", "Within\nRankin Inlet",
                 ifelse(pw$Comparison == "Naujaat", "Within\nNaujaat", "Between\nzones")),
                 levels = c("Within\nRankin Inlet", "Within\nNaujaat", "Between\nzones"))
cmp_cols <- c("Within\nRankin Inlet" = unname(REGION_COLORS["Rankin Inlet"]),
              "Within\nNaujaat"      = unname(REGION_COLORS["Naujaat"]),
              "Between\nzones"       = "grey35")
p_b <- ggplot(pw[is.finite(pw$Fst_weighted), ], aes(cmp, Fst_weighted, colour = cmp)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, fill = NA) +
    geom_jitter(width = 0.12, height = 0, size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = cmp_cols, guide = "none") +
    labs(x = NULL, y = expression(italic(F)[ST]~"(weighted)")) +
    theme_hered()

# --- (c) nucleotide diversity ------------------------------------------------
dv <- read.csv(file.path(TABLE_DIR, "Table_Diversity.csv"))
dv$Pop    <- factor(dv$Pop, levels = POP_ORDER_GEO)
dv$Region <- norm_region(dv$Region)
p_c <- ggplot(dv, aes(Pop, pi * 1e3, fill = Region)) +
    geom_col(width = 0.75) +
    facet_grid(. ~ Region, scales = "free_x", space = "free_x",
               labeller = labeller(Region = ZONE_SHORT)) +
    scale_fill_manual(values = REGION_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = expression(pi~"("%*%10^-3*")")) +
    theme_hered() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          panel.spacing.x = unit(1.5, "mm"),
          strip.text = element_blank())   # zones shown by colour, as in panel a

# --- (d) isolation by distance -----------------------------------------------
mt <- read.csv(file.path(IBD_DIR, "Table_IBD_Mantel.csv"))
mt <- mt[mt$Distance == "Hydro_km", ]
if (!nrow(mt)) stop("No Hydro_km rows in Table_IBD_Mantel.csv")
pw <- pw[is.finite(pw$Fst_lin) & is.finite(pw$Hydro_km), ]
scopes <- data.frame(group = c("Global", "Naujaat", "Rankin Inlet"),
                     label = c("Global (13 populations)", "Naujaat (7)", "Rankin Inlet (5)"),
                     col   = c("grey25", REGION_COLORS["Naujaat"], REGION_COLORS["Rankin Inlet"]))
ibd <- do.call(rbind, lapply(seq_len(nrow(scopes)), function(i) {
    g <- scopes$group[i]
    d <- if (g == "Global") pw else pw[pw$Comparison == g, ]
    if (!nrow(d)) return(NULL)
    d$scope <- scopes$label[i]; d
}))
ibd$scope <- factor(ibd$scope, levels = scopes$label)
stats <- merge(scopes, mt, by.x = "group", by.y = "Group")
stats$scope <- factor(stats$label, levels = scopes$label)
stats$txt <- sprintf("r = %.2f\n%s", stats$Mantel_r,
                     ifelse(stats$Mantel_p < 0.001, "p < 0.001",
                            sprintf("p = %.2f", stats$Mantel_p)))
sig <- as.character(stats$scope[stats$Mantel_p < ALPHA])
col_map <- setNames(scopes$col, scopes$label)
p_d <- ggplot(ibd, aes(Hydro_km, Fst_lin, colour = scope)) +
    geom_point(size = 0.9, alpha = 0.75) +
    geom_smooth(data = ibd[ibd$scope %in% sig, ], method = "lm", formula = y ~ x,
                se = TRUE, linewidth = 0.5, alpha = 0.15,
                aes(fill = scope)) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
              hjust = -0.1, vjust = 1.2, size = 2.3, lineheight = 0.9) +
    facet_wrap(~ scope, nrow = 1, scales = "free") +
    scale_colour_manual(values = col_map, guide = "none") +
    scale_fill_manual(values = col_map, guide = "none") +
    labs(x = "Hydrographic distance (km)",
         y = expression(italic(F)[ST]/(1 - italic(F)[ST]))) +
    theme_hered() +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.3))

# --- assemble ----------------------------------------------------------------
fig6 <- (p_a | p_b | p_c) / p_d +
    plot_layout(heights = c(1, 0.9)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig(fig6, FIG_DIR, "Fig6_differentiation_diversity", WIDTH_2COL, WIDTH_2COL * 0.62)
