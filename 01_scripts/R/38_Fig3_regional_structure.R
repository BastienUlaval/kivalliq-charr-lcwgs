#!/usr/bin/env Rscript
# =============================================================================
# 38_Fig3_regional_structure.R — Figure 3 (Heredity): structure within zones
#   (a) Rankin Inlet sites        (b) Naujaat sites
#   (c) Rankin Inlet PCA          (d) Naujaat PCA      (95% ellipses per pop)
#   (e) Rankin Inlet NGSadmix K=2 (f) Naujaat NGSadmix K=5
#   Panel f bar order: SUP, NOP, PAM, KGJ, ITI, TIN, WHI (ZONES$naujaat$admix_order)
#
# Same conventions as Fig. 2:
# geographic order, population codes under the bars, colour-blind-safe
# cluster palette, thin separators, larger bold site labels on the maps.
#
# Usage:
#   Rscript 38_Fig3_regional_structure.R <base_dir> <shp_dir> [<r_lib_path>]
#   Output: <base_dir>/99_figures/heredity/Fig3_regional_structure.{pdf,png,tiff}
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 38_Fig3_regional_structure.R <base_dir> <shp_dir> [<r_lib_path>]")
BASE_DIR <- args[1]; SHP_DIR <- args[2]
if (length(args) >= 3) .libPaths(args[3])
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "fig_common.R"))
suppressPackageStartupMessages({ library(sf); library(ggrepel) })

SUFFIX  <- "maf0.05_pctind0.50_maxdepth8_prunednosex"
PCA_DIR <- file.path(BASE_DIR, "04_pca")
INFO    <- file.path(BASE_DIR, "02_info", "bamlists")
FIG_DIR <- file.path(BASE_DIR, "99_figures", "heredity")
ZONES <- list(
    rankin  = list(name = "Rankin Inlet", K = 2, pal = PALETTE_RANKIN, bam = file.path(INFO, "rankin.bamlist"),
                   xlim = c(-92.65, -91.05), ylim = c(62.40, 63.00), xbr = c(-92.5, -92, -91.5),
                   town = c(-92.0852853, 62.808375)),
    naujaat = list(name = "Naujaat", K = 5, pal = PALETTE_NAUJAAT, bam = file.path(INFO, "naujaat.bamlist"),
                   xlim = c(-87.5, -84.0), ylim = c(65.85, 67.15), xbr = c(-87, -86, -85, -84),
                   town = c(-86.244657, 66.528295),
                   # bar order of the admixture plot (panel f) only; maps and
                   # PCA keep the geographic order of fig_common.R
                   admix_order = c("SUP", "NOP", "PAM", "KGJ", "ITI", "TIN", "WHI")))

bm <- load_basemap(SHP_DIR, c(-93, -83), c(62, 68), tol = 0.002, rivers = TRUE)

zone_map <- function(z) {
    s <- droplevels(SITES[SITES$region == z$name, ])
    town <- data.frame(lon = z$town[1], lat = z$town[2], name = z$name)
    p <- ggplot() +
        geom_sf(data = bm$land,  fill = "grey90", colour = "grey60", linewidth = 0.15) +
        geom_sf(data = bm$lakes, fill = "white",  colour = "grey70", linewidth = 0.1)
    if (!is.null(bm$rivers)) p <- p + geom_sf(data = bm$rivers, colour = "grey55", linewidth = 0.15)
    p <- p +
        geom_point(data = town, aes(lon, lat), shape = 18, size = 2.5) +
        geom_text_repel(data = town, aes(lon, lat, label = name), size = 2.6,
                        fontface = "italic", colour = "grey20", seed = 2,
                        box.padding = 0.5) +
        geom_point(data = s, aes(lon, lat, colour = pop),
                   shape = REGION_SHAPES[z$name], size = 2.6) +
        geom_label_repel(data = s, aes(lon, lat, label = pop), size = 3.0,
                         fontface = "bold", fill = alpha("white", 0.9),
                         label.size = 0.15, label.padding = unit(0.12, "lines"),
                         segment.size = 0.25, min.segment.length = 0,
                         box.padding = 0.4, max.overlaps = Inf, seed = 1) +
        scale_colour_manual(values = POP_COLORS, guide = "none") +
        scale_x_continuous(breaks = z$xbr) +
        coord_sf(xlim = z$xlim, ylim = z$ylim, expand = FALSE) +
        labs(x = NULL, y = NULL, title = z$name) +
        theme_hered() +
        theme(panel.border = element_rect(fill = NA, linewidth = 0.4),
              plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
              axis.text = element_text(size = 6))
    add_scalebar(p)
}

panels <- lapply(names(ZONES), function(k) {
    z <- ZONES[[k]]
    pca <- load_pca(file.path(PCA_DIR, paste0(k, "_", SUFFIX)), z$bam)
    qopt <- best_qopt(file.path(BASE_DIR, "05_admixture", k), k, z$K)
    # admixture bars: optional zone-specific population order
    geo <- POP_ORDER_GEO
    if (!is.null(z$admix_order)) POP_ORDER_GEO <<- c(z$admix_order, setdiff(geo, z$admix_order))
    adm <- plot_admix(qopt, z$bam, show_regions = FALSE, palette = z$pal) +
        labs(subtitle = sprintf("K = %d", z$K)) +
        theme(plot.subtitle = element_text(size = 7, hjust = 0.5))
    POP_ORDER_GEO <<- geo
    list(map = zone_map(z),
         pca = plot_pca(pca, ellipse = "pop"),
         adm = adm)
})
R <- panels[[1]]; N <- panels[[2]]

fig3 <- (R$map | N$map) / (R$pca | N$pca) / (R$adm | N$adm) +
    plot_layout(heights = c(1.1, 1, 0.55)) +
    plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e", "f")))
save_fig(fig3, FIG_DIR, "Fig3_regional_structure", WIDTH_2COL, WIDTH_2COL * 1.05)
