#!/usr/bin/env Rscript
# =============================================================================
# plot_panel_regional_structure.R  — Figure 3 (CJFAS Kivalliq Arctic charr)
#
# 3-row x 2-column composite:
#   Row 1 (maps)        : a) Rankin Inlet zoom   b) Naujaat zoom
#   Row 2 (PCA)         : c) Rankin PCA          d) Naujaat PCA
#   Row 3 (NGSadmix)    : e) Rankin K=2 barplot  f) Naujaat K=7 barplot
#
# K values updated 2026-07-15 following StructureSelector validation
# (Puechmaille 2016 + fastSTRUCTURE Choose-K method):
#   - Rankin  : K=5 -> K=2  (Puechmaille MedMedK/MedMeaK flat at 1 across K2-7;
#               MaxMedK/MaxMeaK and fastSTRUCTURE converge on K=2. The K=5 peak
#               from the internal Evanno deltaK was an artifact of uneven
#               sampling across the 5 Rankin populations.)
#   - Naujaat : K=5 (unchanged — Evanno, Puechmaille, and fastSTRUCTURE all
#               converge on K=5)
#
# Same palette/shape convention as plot_panel_global_structure.R:
#   Okabe-Ito POP_COLORS, Rankin=circle / Naujaat=square / Baker=triangle.
#
# Usage:
#   Rscript plot_panel_regional_structure.R <base_dir> <shp_dir> [<r_lib_path>]
#     base_dir    = pipeline root (contains 04_pca/, 05_admixture/, 02_info/,
#                   99_figures/)
#     shp_dir     = directory with the basemap shapefiles (Statistics Canada
#                   boundary files + NRCan CanVec lakes/rivers — not included
#                   in this repo, see README)
#     r_lib_path  = optional, only needed if packages are not on the default
#                   R library path
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript plot_panel_regional_structure.R <base_dir> <shp_dir> [<r_lib_path>]")
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

# ─── CONFIG (derived from base_dir/shp_dir arguments) ────────────────────────
PCA_DIR         <- file.path(BASE_DIR, "04_pca")
ADMIX_RANKIN    <- file.path(BASE_DIR, "05_admixture", "rankin")
ADMIX_NAUJAAT   <- file.path(BASE_DIR, "05_admixture", "naujaat")
INFO_DIR        <- file.path(BASE_DIR, "02_info")
FIG_DIR         <- file.path(BASE_DIR, "99_figures")
SUFFIX          <- "maf0.05_pctind0.50_maxdepth8_prunednosex"

RANKIN_BAMLIST  <- file.path(INFO_DIR, "bamlists", "rankin.bamlist")
NAUJAAT_BAMLIST <- file.path(INFO_DIR, "bamlists", "naujaat.bamlist")
RANKIN_K        <- 2
NAUJAAT_K       <- 5

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ─── PALETTE (Okabe-Ito) ─────────────────────────────────────────────────────
POP_COLORS <- c(
    "HOR" = "#009E73",
    "PAM" = "#0072B2", "NOP" = "#56B4E9", "TIN" = "#4E79A7",
    "SUP" = "#9966CC", "KGJ" = "#332288", "WHI" = "#117733", "ITI" = "#44AA99",
    "AUL" = "#E69F00", "MEL" = "#D55E00", "DIA" = "#CC6677",
    "AKL" = "#882255", "CRB" = "#AA4499"
)
REGION_SHAPES <- c("Rankin" = 16, "Naujaat" = 15, "Baker" = 17)
POP_REGION <- c(
    AKL = "Rankin", AUL = "Rankin", CRB = "Rankin", DIA = "Rankin", MEL = "Rankin",
    HOR = "Baker",
    ITI = "Naujaat", KGJ = "Naujaat", NOP = "Naujaat", PAM = "Naujaat",
    SUP = "Naujaat", TIN = "Naujaat", WHI = "Naujaat"
)
POP_ORDER <- c("AKL","AUL","CRB","DIA","MEL","HOR",
               "ITI","KGJ","NOP","PAM","SUP","TIN","WHI")

# Sampling-site coordinates (same as carte_rankin / carte_naujaat scripts)
sites <- tibble::tribble(
    ~pop,  ~region,    ~lon,         ~lat,
    "AKL", "Rankin",   -91.3071148,   62.8372049,
    "AUL", "Rankin",   -91.15114147,  62.905473,
    "CRB", "Rankin",   -92.333333,    62.466667,
    "DIA", "Rankin",   -92.44444,     62.833333,
    "MEL", "Rankin",   -92.116667,    62.866667,
    "ITI", "Naujaat",  -84.4652392,   66.311737,
    "KGJ", "Naujaat",  -85.254424,    66.459357,
    "NOP", "Naujaat",  -86.883333,    66.633333,
    "PAM", "Naujaat",  -87.055555,    67.083333,
    "TIN", "Naujaat",  -84.4398478,   66.3826439,
    "SUP", "Naujaat",  -86.716667,    66.383333,
    "WHI", "Naujaat",  -85.0426395,   65.96663
)
sites$region <- factor(sites$region, levels = c("Rankin","Naujaat","Baker"))
sites_sf <- st_as_sf(sites, coords = c("lon","lat"), crs = 4326)
towns_sf <- st_as_sf(tibble::tribble(
    ~name,           ~lon,         ~lat,
    "Rankin Inlet",  -92.0852853,   62.808375,
    "Naujaat",       -86.244657,    66.528295
), coords = c("lon","lat"), crs = 4326)

# ─── HELPERS ─────────────────────────────────────────────────────────────────
extract_pop <- function(bam_path) {
    toupper(sub(".*saal([A-Za-z]{3}).*", "\\1", basename(bam_path)))
}
load_pca <- function(prefix, bamlist) {
    pca_file <- paste0(prefix, ".cov.pca")
    eig_file <- paste0(prefix, ".cov.eig")
    if (!file.exists(pca_file)) stop("Missing PCA file: ", pca_file)
    pca <- read.table(pca_file, header = TRUE)
    bams <- scan(bamlist, what = "character", quiet = TRUE)
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

# ─── LOAD MAP LAYERS (SHP) ───────────────────────────────────────────────────
bbox_all <- st_bbox(c(xmin = -93, xmax = -83, ymin = 62, ymax = 68), crs = 4326)
bbox_sf  <- st_as_sfc(bbox_all)

cat("Loading shapefiles from ", SHP_DIR, "...\n", sep = "")
canada <- read_sf(file.path(SHP_DIR, "lpr_000b21a_e.shp")) %>%
    st_transform(4326) %>% st_crop(bbox_sf) %>% st_simplify(dTolerance = 0.002)
lakes  <- read_sf(file.path(SHP_DIR, "waterbody_2.shp")) %>%
    st_transform(4326) %>% st_crop(bbox_sf) %>% st_simplify(dTolerance = 0.002)
rivers <- read_sf(file.path(SHP_DIR, "watercourse_1.shp")) %>%
    st_transform(4326) %>% st_crop(bbox_sf) %>% st_simplify(dTolerance = 0.002)
nunavut <- canada[canada$PRNAME == "Nunavut", ]

# ─── PANELS A + B: regional zoom maps ────────────────────────────────────────
make_zoom_map <- function(sites_sf_sub, town_sub, xlim, ylim, label) {
    pop_levels    <- intersect(POP_ORDER, unique(as.character(sites_sf_sub$pop)))
    region_levels <- unique(as.character(sites_sf_sub$region))
    sites_sf_sub$pop    <- factor(sites_sf_sub$pop,    levels = pop_levels)
    sites_sf_sub$region <- factor(sites_sf_sub$region, levels = region_levels)
    ggplot() +
        geom_sf(data = nunavut, fill = "grey92", color = "grey55", linewidth = 0.25) +
        geom_sf(data = lakes,   fill = "white",  color = "grey75", linewidth = 0.15) +
        geom_sf(data = rivers,  color = "grey50", linewidth = 0.25) +
        geom_sf(data = sites_sf_sub,
                aes(colour = pop, shape = region),
                size = 3.5, stroke = 0.7) +
        geom_sf(data = town_sub, shape = 18, size = 4, colour = "black") +
        ggrepel::geom_label_repel(
            data = town_sub,
            aes(label = name, geometry = geometry),
            stat = "sf_coordinates",
            size = 3.0, fontface = "italic",
            colour = "grey25", fill = alpha("white", 0.85),
            label.size = 0.2,
            segment.size = 0.3, segment.colour = "grey50",
            min.segment.length = 0,
            box.padding = 0.5, point.padding = 0.3,
            max.overlaps = 30) +
        ggrepel::geom_label_repel(
            data = sites_sf_sub,
            aes(label = pop, geometry = geometry, colour = pop),
            stat = "sf_coordinates",
            size = 3.0, fontface = "bold",
            label.padding = unit(0.18, "lines"),
            label.size = 0.3,
            segment.size = 0.3, segment.colour = "grey50",
            min.segment.length = 0,
            box.padding = 0.45, point.padding = 0.25,
            max.overlaps = 30, show.legend = FALSE) +
        scale_colour_manual(values = POP_COLORS[pop_levels],
                            breaks = pop_levels, name = "Population") +
        scale_shape_manual(values = REGION_SHAPES[region_levels],
                           breaks = region_levels, name = "Region") +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        annotation_scale(location = "br", line_width = 0.4,
                         height = unit(0.15,"cm"), text_cex = 0.65) +
        labs(title = label) +
        theme_bw(base_size = 10) +
        theme(panel.grid = element_blank(),
              panel.border = element_rect(color = "black", fill = NA),
              plot.title = element_text(face = "bold", size = 11),
              legend.position = "none",
              axis.text = element_text(size = 7))
}
cat("Building regional maps...\n")
p_map_rankin  <- make_zoom_map(
    subset(sites_sf, region == "Rankin"),
    subset(towns_sf, name   == "Rankin Inlet"),
    xlim = c(-92.7, -90.8), ylim = c(62.40, 63.00),
    label = "a) Rankin Inlet sampling sites")
p_map_naujaat <- make_zoom_map(
    subset(sites_sf, region == "Naujaat"),
    subset(towns_sf, name   == "Naujaat"),
    xlim = c(-87.5, -84.0), ylim = c(65.85, 67.15),
    label = "b) Naujaat sampling sites")

# ─── PANELS C + D: regional PCAs (PC1 vs PC2, per-pop 95% ellipses) ──────────
make_regional_pca <- function(pca_obj, label) {
    pca <- pca_obj$data; eig <- pca_obj$eig
    pres_pops <- levels(pca$pop); pres_regs <- levels(pca$region)
    keep_pops <- names(table(pca$pop))[table(pca$pop) >= 4]
    pca_ell   <- pca[pca$pop %in% keep_pops, ]
    p <- ggplot(pca, aes(x = PC1, y = PC2, colour = pop, shape = region))
    if (nrow(pca_ell) > 0) {
        p <- p + stat_ellipse(data = pca_ell,
                              mapping = aes(x = PC1, y = PC2,
                                            group = pop, colour = pop),
                              linetype = "solid", linewidth = 0.45,
                              level = 0.95, inherit.aes = FALSE,
                              show.legend = FALSE)
    }
    p + geom_point(size = 2.4, alpha = 0.85, stroke = 0.4) +
        scale_colour_manual(values = POP_COLORS[pres_pops],
                            breaks = pres_pops, name = "Population") +
        scale_shape_manual(values = REGION_SHAPES[pres_regs],
                           breaks = pres_regs, name = "Region") +
        labs(x = paste0("PC1 (", eig[1], "%)"),
             y = paste0("PC2 (", eig[2], "%)"),
             title = label) +
        theme_bw(base_size = 10) +
        theme(panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = 11),
              legend.position = "right") +
        guides(colour = guide_legend(override.aes = list(size = 3, shape = 16),
                                     order = 1, ncol = 1),
               shape  = "none")
}
cat("Loading regional PCA data...\n")
pca_rankin  <- load_pca(file.path(PCA_DIR, paste0("rankin_",  SUFFIX)),
                        RANKIN_BAMLIST)
pca_naujaat <- load_pca(file.path(PCA_DIR, paste0("naujaat_", SUFFIX)),
                        NAUJAAT_BAMLIST)
p_pca_rankin  <- make_regional_pca(pca_rankin,  "c) Rankin Inlet PCA")
p_pca_naujaat <- make_regional_pca(pca_naujaat, "d) Naujaat PCA")

# ─── PANELS E + F: regional NGSadmix barplots ────────────────────────────────
make_regional_admix <- function(qopt_file, bamlist, K, label, subset_pops) {
    Q    <- read.table(qopt_file, header = FALSE)
    bams <- scan(bamlist, what = "character", quiet = TRUE)
    pops <- sapply(bams, extract_pop)
    df <- data.frame(ind = seq_len(nrow(Q)), pop = pops)
    for (j in seq_len(ncol(Q))) df[[paste0("K", j)]] <- Q[[j]]
    sub_order  <- POP_ORDER[POP_ORDER %in% subset_pops]
    df$pop_f   <- factor(df$pop, levels = sub_order)
    df$max_Q   <- apply(Q, 1, max)
    df <- df %>% arrange(pop_f, desc(max_Q))
    df$x <- seq_len(nrow(df))
    long <- df %>% pivot_longer(cols = starts_with("K"),
                                names_to = "cluster", values_to = "Q")
    pop_info <- df %>% group_by(pop_f) %>%
        summarise(xmin = min(x), xmax = max(x), .groups = "drop") %>%
        mutate(xmid = (xmin + xmax) / 2)
    separators <- head(pop_info$xmax, -1) + 0.5
    cols <- viridisLite::viridis(K)
    ggplot(long, aes(x = x, y = Q, fill = cluster)) +
        geom_bar(stat = "identity", width = 1) +
        geom_vline(xintercept = separators, color = "white", linewidth = 1.8) +
        scale_fill_manual(values = cols) +
        scale_x_continuous(breaks = pop_info$xmid,
                           labels = as.character(pop_info$pop_f),
                           expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(title = label, y = "Ancestral proportion", x = "") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(size = 9, face = "bold",
                                         angle = 45, hjust = 1),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              legend.position = "none",
              plot.title = element_text(face = "bold", size = 11),
              panel.border = element_rect(color = "black", fill = NA,
                                          linewidth = 0.4))
}
cat("Loading NGSadmix logs and selecting best replicates...\n")
lnL_rankin  <- parse_all_lnL(ADMIX_RANKIN,  "rankin")
lnL_naujaat <- parse_all_lnL(ADMIX_NAUJAAT, "naujaat")
brep_rankin  <- best_rep_for_K(lnL_rankin,  RANKIN_K)
brep_naujaat <- best_rep_for_K(lnL_naujaat, NAUJAAT_K)
if (is.na(brep_rankin))  stop("No rankin_K", RANKIN_K, " .log files in ", ADMIX_RANKIN)
if (is.na(brep_naujaat)) stop("No naujaat_K", NAUJAAT_K, " .log files in ", ADMIX_NAUJAAT)
qopt_rankin  <- file.path(ADMIX_RANKIN,
                          sprintf("rankin_K%d_rep%d.qopt",  RANKIN_K,  brep_rankin))
qopt_naujaat <- file.path(ADMIX_NAUJAAT,
                          sprintf("naujaat_K%d_rep%d.qopt", NAUJAAT_K, brep_naujaat))
cat("  rankin  K=", RANKIN_K,  " best_rep=", brep_rankin,  "\n", sep = "")
cat("  naujaat K=", NAUJAAT_K, " best_rep=", brep_naujaat, "\n", sep = "")

p_admix_rankin  <- make_regional_admix(qopt_rankin,  RANKIN_BAMLIST,  RANKIN_K,
                                       sprintf("e) Rankin Inlet NGSadmix (K=%d)", RANKIN_K),
                                       c("AKL","AUL","CRB","DIA","MEL"))
p_admix_naujaat <- make_regional_admix(qopt_naujaat, NAUJAAT_BAMLIST, NAUJAAT_K,
                                       sprintf("f) Naujaat NGSadmix (K=%d)", NAUJAAT_K),
                                       c("ITI","KGJ","NOP","PAM","SUP","TIN","WHI"))

# ─── ASSEMBLE 3 rows x 2 cols ────────────────────────────────────────────────
cat("Assembling Figure 3...\n")
final <- (p_map_rankin   | p_map_naujaat) /
         (p_pca_rankin   | p_pca_naujaat) /
         (p_admix_rankin | p_admix_naujaat) +
         plot_layout(heights = c(1.15, 1.0, 0.55))
ggsave(file.path(FIG_DIR, "Fig3_Regional_structure.pdf"), final,
       width = 13, height = 12, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig3_Regional_structure.png"), final,
       width = 13, height = 12, dpi = 300)
cat("\n=== DONE === Fig3_Regional_structure saved to", FIG_DIR, "\n")