#!/usr/bin/env Rscript
# =============================================================================
# plot_panel_global_structure.R  -- Figure 2 (CJFAS Kivalliq Arctic charr)
#
# 2-panel composite:
#   a) Global PCA (13 pops, PC1 vs PC2, 95% ellipses per region)
#   b) Global NGSadmix barplot (best-K replicate)
#
# K=7 chosen as the best-supported consensus across Puechmaille (2016)
# estimators (MedMedK, MedMeaK, MaxMedK, MaxMeaK), computed via the
# StructureSelector web tool (Li & Liu 2018) from the NGSadmix per-replicate
# ancestry matrices and log-likelihoods -- not by a script in this repo.
# Evanno deltaK alone was judged unreliable here (bias under unbalanced
# sampling; see manuscript Methods).
#
# Palette = Okabe-Ito (plot_pca.R convention) for populations.
# Admixture cluster colours = viridis.
# Shapes: Rankin = circle (16), Naujaat = square (15), Baker = triangle (17).
#
# NOTE: basemap fills harmonised with plot_panel_regional_structure.R
#   (continent = grey92, lakes = white, sea/background = white) so that the
#   population points stand out as on the Regional figure.
#
# Usage:
#   Rscript plot_panel_global_structure.R <base_dir> <shp_dir> [<r_lib_path>]
#     base_dir    = pipeline root (contains 04_pca/, 05_admixture/, 02_info/,
#                   99_figures/)
#     shp_dir     = directory with the 3 basemap shapefiles (Statistics Canada
#                   boundary files lpr_000b21a_e.shp + NRCan CanVec
#                   waterbody_2.shp / watercourse_1.shp -- not included in this
#                   repo, see README)
#     r_lib_path  = optional, only needed if packages are not on the default
#                   R library path (passed to .libPaths())
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript plot_panel_global_structure.R <base_dir> <shp_dir> [<r_lib_path>]")
}
BASE_DIR <- args[1]
SHP_DIR  <- args[2]
if (length(args) >= 3) .libPaths(args[3])

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(viridisLite)
    library(sf)
    library(ggspatial)
    library(ggrepel)
})
sf_use_s2(FALSE)

# --- CONFIG (derived from base_dir/shp_dir arguments) ------------------------
PCA_DIR   <- file.path(BASE_DIR, "04_pca")
ADMIX_DIR <- file.path(BASE_DIR, "05_admixture", "global")
LOG_DIR   <- ADMIX_DIR                       # NGSadmix .log files live alongside .qopt
INFO_DIR  <- file.path(BASE_DIR, "02_info")
FIG_DIR   <- file.path(BASE_DIR, "99_figures")
SUFFIX    <- "maf0.05_pctind0.50_maxdepth8_prunednosex"

# Bamlist (individual order, matching .cov.pca and .qopt row order)
GLOBAL_BAMLIST <- file.path(INFO_DIR, "bam.filelist")

# Admixture prefix and the K to display for the global dataset
GLOBAL_ADMIX_PREFIX <- "global"
GLOBAL_K            <- 7

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- PALETTE (Okabe-Ito, from plot_pca_v2.R) ---------------------------------
POP_COLORS <- c(
    "HOR" = "#009E73",
    "PAM" = "#0072B2", "NOP" = "#56B4E9", "TIN" = "#4E79A7",
    "SUP" = "#9966CC", "KGJ" = "#332288", "WHI" = "#117733", "ITI" = "#44AA99",
    "AUL" = "#E69F00", "MEL" = "#D55E00", "DIA" = "#CC6677",
    "AKL" = "#882255", "CRB" = "#AA4499"
)
REGION_SHAPES         <- c("Rankin" = 16, "Naujaat" = 15, "Baker" = 17)
REGION_ELLIPSE_COLORS <- c("Rankin" = "#D62728", "Naujaat" = "#0072B2", "Baker" = "#009E73")
POP_REGION <- c(
    AKL = "Rankin", AUL = "Rankin", CRB = "Rankin", DIA = "Rankin", MEL = "Rankin",
    HOR = "Baker",
    ITI = "Naujaat", KGJ = "Naujaat", NOP = "Naujaat", PAM = "Naujaat",
    SUP = "Naujaat", TIN = "Naujaat", WHI = "Naujaat"
)
POP_ORDER <- c("AKL","AUL","CRB","DIA","MEL","HOR",
               "ITI","KGJ","NOP","PAM","SUP","TIN","WHI")

# --- HELPERS (match plot_pca.R conventions) ---------------------------------
extract_pop <- function(bam_path) {
    # Extracts the 3-letter population code from a BAM filename: ".*saal(XXX).*" -> XXX
    toupper(sub(".*saal([A-Za-z]{3}).*", "\\1", basename(bam_path)))
}

load_pca <- function(prefix) {
    pca_file <- paste0(prefix, ".cov.pca")
    eig_file <- paste0(prefix, ".cov.eig")
    if (!file.exists(pca_file)) stop("Missing PCA file: ", pca_file)
    pca <- read.table(pca_file, header = TRUE)
    bams <- scan(GLOBAL_BAMLIST, what = "character", quiet = TRUE)
    if (nrow(pca) != length(bams))
        warning("PCA rows (", nrow(pca), ") != bamlist length (", length(bams), ")")
    pca$pop    <- sapply(bams, extract_pop)
    pca$region <- POP_REGION[pca$pop]
    present_pops <- intersect(POP_ORDER, unique(pca$pop))
    pca$pop    <- factor(pca$pop, levels = present_pops)
    pca$region <- factor(pca$region,
                         levels = intersect(c("Baker","Naujaat","Rankin"),
                                            unique(pca$region)))
    list(data = pca, eig = scan(eig_file, quiet = TRUE))
}

# NGSadmix: parse all .log lnL, pick best replicate for a given K
parse_all_lnL <- function(log_dir, prefix) {
    lnL <- data.frame(K = integer(), rep = integer(), lnL = numeric())
    all_logs <- list.files(log_dir,
                           pattern = paste0("^", prefix, "_K.*\\.log$"),
                           full.names = TRUE)
    for (lf in all_logs) {
        lines   <- readLines(lf, warn = FALSE)
        ll_line <- grep("best like=", lines, value = TRUE)
        if (length(ll_line) > 0) {
            ll <- as.numeric(sub(".*best like=([^ ]+).*", "\\1", ll_line[1]))
            k  <- as.numeric(sub(".*_K(\\d+)_.*", "\\1", basename(lf)))
            r  <- as.numeric(sub(".*rep(\\d+).*", "\\1", basename(lf)))
            lnL <- rbind(lnL, data.frame(K = k, rep = r, lnL = ll))
        }
    }
    lnL
}
best_rep_for_K <- function(lnL, k) {
    sub <- lnL %>% filter(K == k) %>% arrange(desc(lnL)) %>% slice(1)
    if (nrow(sub) > 0) sub$rep else NA
}

# --- SAMPLING SITES (real coords + town markers) -----------------------------
sites <- tibble::tribble(
    ~pop,  ~region,    ~lon,         ~lat,
    "AKL", "Rankin",   -91.3071148,   62.8372049,
    "AUL", "Rankin",   -91.15114147,  62.905473,
    "CRB", "Rankin",   -92.333333,    62.466667,
    "DIA", "Rankin",   -92.44444,     62.833333,
    "MEL", "Rankin",   -92.116667,    62.866667,
    "HOR", "Baker",    -97.933631,    64.723358,
    "ITI", "Naujaat",  -84.4652392,   66.311737,
    "KGJ", "Naujaat",  -85.254424,    66.459357,
    "NOP", "Naujaat",  -86.883333,    66.633333,
    "PAM", "Naujaat",  -87.055555,    67.083333,
    "TIN", "Naujaat",  -84.4398478,   66.3826439,
    "SUP", "Naujaat",  -86.716667,    66.383333,
    "WHI", "Naujaat",  -85.0426395,   65.96663
)
sites$region <- factor(sites$region, levels = c("Rankin","Naujaat","Baker"))
sites$pop    <- factor(sites$pop, levels = intersect(POP_ORDER, sites$pop))
sites_sf <- st_as_sf(sites, coords = c("lon","lat"), crs = 4326)
towns_sf <- st_as_sf(tibble::tribble(
    ~name,           ~lon,         ~lat,
    "Rankin Inlet",  -92.0852853,   62.808375,
    "Naujaat",       -86.244657,    66.528295,
    "Baker Lake",    -96.077,       64.318
), coords = c("lon","lat"), crs = 4326)

# --- PANEL A: Global map (MU-style: region labels only, no town markers) -----
cat("Loading shapefiles from ", SHP_DIR, "...\n", sep = "")
bbox_A   <- st_bbox(c(xmin = -99, xmax = -82, ymin = 61.5, ymax = 68.5), crs = 4326)
bbox_sfA <- st_as_sfc(bbox_A)
canadaA  <- read_sf(file.path(SHP_DIR, "lpr_000b21a_e.shp")) %>%
    st_transform(4326) %>% st_crop(bbox_sfA) %>% st_simplify(dTolerance = 0.005)
lakesA   <- read_sf(file.path(SHP_DIR, "waterbody_2.shp")) %>%
    st_transform(4326) %>% st_crop(bbox_sfA) %>% st_simplify(dTolerance = 0.005)

cat("Building Panel A (global map)...\n")
pres_pops_map <- intersect(POP_ORDER, unique(as.character(sites_sf$pop)))
pres_regs_map <- intersect(c("Baker","Naujaat","Rankin"),
                           unique(as.character(sites_sf$region)))

# Region label colours (Rankin=red, Naujaat=blue, Baker=green; explicitly named
# to avoid collisions with the per-population palette used by the points).
REGION_LABEL_COLORS <- c("Rankin"  = "#D62728",
                         "Naujaat" = "#0072B2",
                         "Baker"   = "#009E73")
REGION_LABEL_TEXT <- c("Rankin"  = "Rankin Inlet",
                       "Naujaat" = "Naujaat",
                       "Baker"   = "Baker Lake")

# Hand-picked label positions, placed in free ocean / land areas so that the
# community label never sits on top of a population point.
region_labels_df <- data.frame(
    region   = c("Rankin",  "Naujaat", "Baker"),
    lon_lab  = c(-90.0,     -84.5,     -95.5),   # circles drawn by Bastien
    lat_lab  = c(62.0,      67.0,      64.7),
    stringsAsFactors = FALSE)
region_labels_df$region <- factor(region_labels_df$region,
                                  levels = pres_regs_map)
region_labels_df$text   <- REGION_LABEL_TEXT[as.character(region_labels_df$region)]
region_labels_df$colour <- REGION_LABEL_COLORS[as.character(region_labels_df$region)]

build_global_map <- function(legend_position = "none",
                             point_size = 2.0,
                             label_size = 4.0,
                             pop_label_size = 2.8) {
    p <- ggplot() +
        geom_sf(data = canadaA, fill = "grey92", color = "grey55", linewidth = 0.25) +
        geom_sf(data = lakesA,  fill = "white",  color = "grey75", linewidth = 0.15) +
        geom_sf(data = sites_sf,
                aes(colour = pop, shape = region),
                size = point_size, stroke = 0.5) +
        ggrepel::geom_text_repel(
            data = sites_sf,
            aes(label = pop, geometry = geometry, colour = pop),
            stat = "sf_coordinates",
            size = pop_label_size, fontface = "bold",
            segment.size = 0.25, segment.colour = "grey55",
            min.segment.length = 0,
            box.padding = 0.35, point.padding = 0.20,
            max.overlaps = 50, show.legend = FALSE) +
        geom_label(data = region_labels_df,
                   aes(x = lon_lab, y = lat_lab, label = text),
                   colour = region_labels_df$colour,
                   fontface = "bold", size = label_size,
                   label.padding = unit(0.25, "lines"),
                   label.r = unit(0.15, "lines"), label.size = 0.4,
                   fill = alpha("white", 0.93)) +
        scale_colour_manual(values = POP_COLORS[pres_pops_map],
                            breaks = pres_pops_map, name = "Population",
                            guide = "none") +
        scale_shape_manual(values = REGION_SHAPES[pres_regs_map],
                           breaks = pres_regs_map, name = "Region") +
        coord_sf(xlim = c(-99, -82), ylim = c(61.5, 68.5), expand = FALSE) +
        annotation_scale(location = "br", line_width = 0.4,
                         height = unit(0.18, "cm"), text_cex = 0.65) +
        theme_bw(base_size = 11) +
        theme(panel.background = element_rect(fill = "white", color = NA),
              panel.grid.major  = element_line(color = "grey85",
                                               linewidth = 0.2, linetype = "dotted"),
              panel.border = element_rect(color = "black", fill = NA),
              plot.title = element_text(face = "bold", size = 12),
              legend.position = legend_position,
              axis.text = element_text(size = 7.5))
    if (legend_position != "none") {
        p <- p + guides(
            shape = guide_legend(override.aes = list(size = 4, colour = "grey25"),
                                 order = 1))
    }
    p
}

# Map for the composite figure (no legend, no title - legend comes from PCA)
p_map <- build_global_map(legend_position = "none",
                          point_size = 2.0, label_size = 3.6) +
    labs(title = "a) Sampling sites across the Kivalliq region",
         x = NULL, y = NULL)

# Standalone map for Figure 1 of the manuscript (full legend, no panel-letter)
p_map_standalone <- build_global_map(legend_position = "right",
                                     point_size = 2.5, label_size = 4.3) +
    labs(title = NULL, x = NULL, y = NULL)

# --- Inset: whole Hudson Bay overview with a red box on the study region ----
# Low-detail locator map (coarse simplify) showing where the study area sits
# in the NW corner of Hudson Bay. Added to the STANDALONE Fig1 only.
cat("Building Hudson Bay inset...\n")
bbox_inset    <- st_bbox(c(xmin = -100, xmax = -74, ymin = 51, ymax = 69), crs = 4326)
bbox_sf_inset <- st_as_sfc(bbox_inset)
canada_inset  <- read_sf(file.path(SHP_DIR, "lpr_000b21a_e.shp")) %>%
    st_transform(4326) %>% st_crop(bbox_sf_inset) %>% st_simplify(dTolerance = 0.03)

# Red rectangle = extent of the main map (the study region)
study_box <- data.frame(xmin = -99, xmax = -82, ymin = 61.5, ymax = 68.5)

p_inset <- ggplot() +
    geom_sf(data = canada_inset, fill = "grey80", colour = "grey55", linewidth = 0.1) +
    geom_rect(data = study_box,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, colour = "red", linewidth = 0.7) +
    coord_sf(xlim = c(-100, -74), ylim = c(51, 69), expand = FALSE) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", colour = "grey30",
                                          linewidth = 0.5),
          plot.margin = margin(1, 1, 1, 1))

# Drop the inset into the top-left corner of the standalone map panel
p_map_standalone <- p_map_standalone +
    inset_element(p_inset,
                  left = 0.01, bottom = 0.66, right = 0.33, top = 0.99,
                  align_to = "panel")

cat("Saving Fig1_Overview_map (standalone) ...\n")
ggsave(file.path(FIG_DIR, "Fig1_Overview_map.pdf"), p_map_standalone,
       width = 9, height = 7, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig1_Overview_map.png"), p_map_standalone,
       width = 9, height = 7, dpi = 300)

# --- PANEL B: Global PCA -----------------------------------------------------
global_prefix <- file.path(PCA_DIR, paste0("global_", SUFFIX))
cat("Loading global PCA:", global_prefix, "\n")
pca_global <- load_pca(global_prefix)

make_global_pca <- function(pca_obj) {
    pca <- pca_obj$data; eig <- pca_obj$eig
    pres_pops <- levels(pca$pop); pres_regs <- levels(pca$region)
    p <- ggplot(pca, aes(x = PC1, y = PC2, colour = pop, shape = region))
    for (reg in pres_regs) {
        p <- p + stat_ellipse(data = pca[pca$region == reg, ],
                              aes(x = PC1, y = PC2),
                              colour = REGION_ELLIPSE_COLORS[reg],
                              linetype = "solid", linewidth = 0.6,
                              level = 0.95, inherit.aes = FALSE,
                              show.legend = FALSE)
    }
    p + geom_point(size = 2.7, alpha = 0.85, stroke = 0.4) +
        scale_colour_manual(values = POP_COLORS[pres_pops],
                            breaks = pres_pops, name = "Population") +
        scale_shape_manual(values = REGION_SHAPES[pres_regs],
                           breaks = pres_regs, name = "Region") +
        labs(x = paste0("PC1 (", eig[1], "%)"),
             y = paste0("PC2 (", eig[2], "%)"),
             title = "b) Global PCA (all 13 populations)") +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = 12),
              legend.position = "right") +
        guides(colour = guide_legend(override.aes = list(size = 3.5, shape = 16),
                                     order = 1, ncol = 1),
               shape  = guide_legend(override.aes = list(size = 3.5, colour = "grey25"),
                                     order = 2))
}
p_pca <- make_global_pca(pca_global)

# --- PANEL C: Global NGSadmix barplot ----------------------------------------
cat("Selecting best NGSadmix replicate for K =", GLOBAL_K, "\n")
lnL  <- parse_all_lnL(LOG_DIR, GLOBAL_ADMIX_PREFIX)
brep <- best_rep_for_K(lnL, GLOBAL_K)
if (is.na(brep)) stop("No NGSadmix .log found for prefix '", GLOBAL_ADMIX_PREFIX,
                      "' K=", GLOBAL_K, " in ", LOG_DIR)
qopt_file <- file.path(ADMIX_DIR,
                       sprintf("%s_K%d_rep%d.qopt", GLOBAL_ADMIX_PREFIX, GLOBAL_K, brep))
cat("  Best replicate:", brep, "->", qopt_file, "\n")

make_global_admix <- function(qopt_file, bamlist, K) {
    Q    <- read.table(qopt_file, header = FALSE)
    bams <- scan(bamlist, what = "character", quiet = TRUE)
    pops <- sapply(bams, extract_pop)
    df <- data.frame(ind = seq_len(nrow(Q)), pop = pops)
    for (j in seq_len(ncol(Q))) df[[paste0("K", j)]] <- Q[[j]]
    df$region_lab <- c(AKL="Rankin Inlet",AUL="Rankin Inlet",CRB="Rankin Inlet",
                       DIA="Rankin Inlet",MEL="Rankin Inlet",HOR="Baker Lake",
                       ITI="Naujaat",KGJ="Naujaat",NOP="Naujaat",PAM="Naujaat",
                       SUP="Naujaat",TIN="Naujaat",WHI="Naujaat")[df$pop]
    df$pop_f    <- factor(df$pop, levels = POP_ORDER)
    df$region_f <- factor(df$region_lab,
                          levels = c("Rankin Inlet","Baker Lake","Naujaat"))
    df$max_Q <- apply(Q, 1, max)
    df <- df %>% arrange(region_f, pop_f, desc(max_Q))
    df$x <- seq_len(nrow(df))
    long <- df %>% pivot_longer(cols = starts_with("K"),
                                names_to = "cluster", values_to = "Q")
    region_info <- df %>% group_by(region_f) %>%
        summarise(xmin = min(x), xmax = max(x), .groups = "drop") %>%
        mutate(xmid = (xmin + xmax) / 2)
    separators <- head(region_info$xmax, -1) + 0.5
    cols <- viridisLite::viridis(K)
    ggplot(long, aes(x = x, y = Q, fill = cluster)) +
        geom_bar(stat = "identity", width = 1) +
        geom_vline(xintercept = separators, color = "white", linewidth = 2.5) +
        scale_fill_manual(values = cols) +
        scale_x_continuous(breaks = region_info$xmid,
                           labels = as.character(region_info$region_f),
                           expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(title = sprintf("c) Global NGSadmix (K=%d)", K),
             y = "Ancestral proportion", x = "") +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(size = 9, face = "bold",
                                         angle = 45, hjust = 1),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              legend.position = "none",
              plot.title = element_text(face = "bold", size = 12),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))
}
p_admix <- make_global_admix(qopt_file, GLOBAL_BAMLIST, GLOBAL_K)

# --- ASSEMBLE (2 columns: map left, PCA + admixture stacked right) -----------
right_col <- (p_pca / p_admix) + plot_layout(heights = c(1.5, 1.0))
final <- (p_map | right_col) + plot_layout(widths = c(1.0, 1.3))
ggsave(file.path(FIG_DIR, "Fig2_Global_structure.pdf"), final,
       width = 16, height = 10, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig2_Global_structure.png"), final,
       width = 16, height = 10, dpi = 300)
cat("\n=== DONE === Fig2_Global_structure + Fig1_Overview_map saved to", FIG_DIR, "\n")