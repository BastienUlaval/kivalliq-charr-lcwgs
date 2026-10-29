#!/usr/bin/env Rscript
# =============================================================================
# 15_Fig1_map.R — Figure 1 (Heredity): study area and sampling design
#
# Map panel of Fig. 1 (study area).
# Change requested by JS (c145): population codes larger and in bold.
# Codes are now drawn in black on a white label, next to points coloured by
# population and shaped by zone; zone names in their zone colour.
#
# Usage:
#   Rscript 15_Fig1_map.R <base_dir> <shp_dir> [<r_lib_path>]
#   Output: <base_dir>/99_figures/heredity/Fig1_map.{pdf,png,tiff}
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 15_Fig1_map.R <base_dir> <shp_dir> [<r_lib_path>]")
BASE_DIR <- args[1]; SHP_DIR <- args[2]
if (length(args) >= 3) .libPaths(args[3])
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "fig_common.R"))
suppressPackageStartupMessages({ library(sf); library(ggrepel) })
FIG_DIR <- file.path(BASE_DIR, "99_figures", "heredity")

XLIM <- c(-99, -82); YLIM <- c(61.5, 68.5)
bm <- load_basemap(SHP_DIR, XLIM, YLIM)
sites_sf <- st_as_sf(SITES, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

zone_lab <- data.frame(region = factor(REGION_ORDER, levels = REGION_ORDER),
                       lon = c(-95.6, -90.0, -84.4), lat = c(64.35, 62.05, 67.05))

p_main <- ggplot() +
    geom_sf(data = bm$land,  fill = "grey90", colour = "grey60", linewidth = 0.15) +
    geom_sf(data = bm$lakes, fill = "white",  colour = "grey70", linewidth = 0.1) +
    geom_point(data = SITES, aes(lon, lat, colour = pop, shape = region), size = 2.4) +
    geom_label_repel(data = SITES, aes(lon, lat, label = pop),
                     size = 3.2, fontface = "bold", colour = "black",
                     fill = alpha("white", 0.9), label.size = 0.15,
                     label.padding = unit(0.12, "lines"),
                     segment.size = 0.25, segment.colour = "grey40",
                     min.segment.length = 0, box.padding = 0.35,
                     max.overlaps = Inf, seed = 1) +
    geom_label(data = zone_lab, aes(lon, lat, label = region, colour = region),
               fontface = "bold", size = 3.4, fill = alpha("white", 0.92),
               label.size = 0.3) +
    scale_colour_manual(values = c(POP_COLORS, REGION_COLORS), guide = "none") +
    scale_shape_manual(values = REGION_SHAPES, name = "Zone") +
    coord_sf(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_hered() +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.4),
          panel.grid.major = element_line(colour = "grey85", linewidth = 0.15,
                                          linetype = "dotted"),
          legend.position = c(0.86, 0.12),
          legend.background = element_rect(fill = alpha("white", 0.9), colour = NA))
p_main <- add_scalebar(p_main)

# Inset: Hudson Bay with the study area
bm_in <- load_basemap(SHP_DIR, c(-100, -74), c(51, 69), tol = 0.03)
p_inset <- ggplot() +
    geom_sf(data = bm_in$land, fill = "grey75", colour = NA) +
    annotate("rect", xmin = XLIM[1], xmax = XLIM[2], ymin = YLIM[1], ymax = YLIM[2],
             fill = NA, colour = "red", linewidth = 0.5) +
    coord_sf(xlim = c(-100, -74), ylim = c(51, 69), expand = FALSE) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", colour = "grey30", linewidth = 0.4))

fig1 <- p_main + inset_element(p_inset, left = 0.01, bottom = 0.68, right = 0.30,
                               top = 0.99, align_to = "panel")
save_fig(fig1, FIG_DIR, "Fig1_map", WIDTH_2COL, WIDTH_2COL * 0.72)
