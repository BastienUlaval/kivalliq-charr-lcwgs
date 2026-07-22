#!/usr/bin/env Rscript
# =============================================================================
# plot_management_units.R  -- Figure 7 (CJFAS Kivalliq Arctic charr)
#
# Six-panel synthesis figure - renamed from plot_conservation_units.R to use
# "Management Units" (MU) at the Kivalliq regional scale, while explicitly
# nesting these MUs within the Hudson Bay Conservation Unit defined at the
# pan-Canadian scale by Dallaire et al. 2026 (Evol. Appl. 19:e70259).
#
# Panels:
#   A. Sampling design map (shapefile-based, same SHP as carte_*.R)
#   B. Pairwise FST boxplots (within Rankin / within Naujaat / between regions)
#   C. Nucleotide diversity (pi) per population
#   D. Atlantic ancestry (D-statistic) per population
#   E. Isolation-by-distance scatter panels (Global / Naujaat / Rankin)
#   F. Recommended management framework (tiered MU cards)
#
# Inputs:
#   ${TABLE_DIR}/Table_Diversity.csv      - Pop, Region, pi
#   ${DSTATS_DIR}/Table_Dstats_JAY.tsv    - D-statistic per population
#   ${IBD_DIR}/Table_IBD_pairwise.csv     - Fst_lin, Hydro_km, Comparison
#   ${IBD_DIR}/Table_IBD_Mantel.csv       - Group, Distance, Mantel_r, Mantel_p
#
# Usage:
#   Rscript plot_management_units.R <base_dir> <shp_dir> [<r_lib_path>]
#     base_dir    = pipeline root (contains 99_tables/, figures/dstats/,
#                   14_ibd/, 99_figures/)
#     shp_dir     = directory with the basemap shapefiles (Statistics Canada
#                   boundary files + NRCan CanVec lakes -- not included in
#                   this repo, see README)
#     r_lib_path  = optional, only needed if packages are not on the default
#                   R library path
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript plot_management_units.R <base_dir> <shp_dir> [<r_lib_path>]")
}
BASE_DIR <- args[1]
SHP_DIR  <- args[2]
if (length(args) >= 3) .libPaths(args[3])

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(scales)
    library(grid)
    library(sf)
    library(ggspatial)
    library(ggrepel)
})
sf_use_s2(FALSE)

# --- CONFIG (derived from base_dir/shp_dir arguments) ------------------------
TABLE_DIR  <- file.path(BASE_DIR, "99_tables", "main")
DSTATS_DIR <- file.path(BASE_DIR, "figures", "dstats")
FIG_DIR    <- file.path(BASE_DIR, "99_figures")

# Set FALSE if HOR being a freshwater outlier biases the global IBD regression
INCLUDE_HOR_IN_GLOBAL <- TRUE

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- PALETTE (regions; vivid colours for the synthesis figure) ---------------
region_colors <- c("Rankin"      = "#D62728",  # rouge vif
                   "Naujaat"     = "#0072B2",  # bleu vif
                   "Baker (HOR)" = "#009E73")  # vert vif
region_levels <- c("Rankin", "Naujaat", "Baker (HOR)")

clean_region <- function(x) {
    x <- as.character(x)
    x <- gsub("Rankin Inlet|^Rankin$",    "Rankin",      x)
    x <- gsub("Baker Lake|^Baker$|^HOR$", "Baker (HOR)", x)
    x <- gsub("^Naujaat$",                 "Naujaat",     x)
    factor(x, levels = region_levels)
}

theme_pub <- theme_classic(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0),
          plot.subtitle = element_text(size = 9, color = "grey35", hjust = 0),
          axis.title    = element_text(size = 10),
          axis.text     = element_text(size = 9),
          legend.position = "none",
          plot.margin   = margin(8, 8, 8, 8))
theme_ibd <- theme_classic(base_size = 9) +
    theme(plot.title    = element_text(face = "bold", size = 10, hjust = 0.5),
          axis.title    = element_text(size = 8.5),
          axis.text     = element_text(size = 8),
          legend.position = "none",
          plot.margin   = margin(4, 6, 4, 6),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.25))

# --- SAMPLING SITES (with real coords) ---------------------------------------
sites <- tibble::tribble(
    ~pop,  ~name,                ~region,       ~lon,         ~lat,
    "AKL", "Akaalik",            "Rankin",      -91.3071148,   62.8372049,
    "AUL", "Aulatsivik",         "Rankin",      -91.15114147,  62.905473,
    "CRB", "Corbett Inlet",      "Rankin",      -92.333333,    62.466667,
    "DIA", "Diana River",        "Rankin",      -92.44444,     62.833333,
    "MEL", "Melania River",      "Rankin",      -92.116667,    62.866667,
    "HOR", "HorseShoe Island",   "Baker (HOR)", -97.933631,    64.723358,
    "ITI", "Itirjuk",            "Naujaat",     -84.4652392,   66.311737,
    "KGJ", "Kuugaarjuk",         "Naujaat",     -85.254424,    66.459357,
    "NOP", "North Pole Lake",    "Naujaat",     -86.883333,    66.633333,
    "PAM", "Pamiurluk Lake",     "Naujaat",     -87.055555,    67.083333,
    "TIN", "Tinujivik",          "Naujaat",     -84.4398478,   66.3826439,
    "SUP", "Sipujjaqtuq River",  "Naujaat",     -86.716667,    66.383333,
    "WHI", "White Island",       "Naujaat",     -85.0426395,   65.96663
)
sites$region <- factor(sites$region, levels = region_levels)
sites_sf <- st_as_sf(sites, coords = c("lon","lat"), crs = 4326)

# --- PANEL A: Sampling map (real SHP) ----------------------------------------
cat("Loading shapefiles from ", SHP_DIR, "...\n", sep = "")
bbox_A  <- st_bbox(c(xmin = -99, xmax = -82, ymin = 61.5, ymax = 68.5), crs = 4326)
bbox_sfA <- st_as_sfc(bbox_A)
canadaA <- read_sf(file.path(SHP_DIR, "lpr_000b21a_e.shp")) %>%
    st_transform(4326) %>% st_crop(bbox_sfA) %>% st_simplify(dTolerance = 0.005)
lakesA  <- read_sf(file.path(SHP_DIR, "waterbody_2.shp")) %>%
    st_transform(4326) %>% st_crop(bbox_sfA) %>% st_simplify(dTolerance = 0.005)

labels_df <- sites %>% group_by(region) %>%
    summarise(lon = mean(lon), lat = mean(lat), .groups = "drop") %>%
    mutate(lon_lab = case_when(
        region == "Rankin"      ~ -96.0,
        region == "Naujaat"     ~ -84.5,
        region == "Baker (HOR)" ~ -95.5,
        TRUE                    ~ lon),
        lat_lab = case_when(
            region == "Rankin"      ~ 62.5,
            region == "Naujaat"     ~ 67.7,
            region == "Baker (HOR)" ~ 64.7,
            TRUE                    ~ lat))

panel_A <- ggplot() +
    geom_sf(data = canadaA, fill = "grey92", color = "grey55", linewidth = 0.25) +
    geom_sf(data = lakesA,  fill = "white",  color = "grey75", linewidth = 0.15) +
    geom_sf(data = sites_sf, aes(color = region), size = 1.3, alpha = 0.95) +
    ggrepel::geom_text_repel(
        data = sites_sf,
        aes(label = pop, geometry = geometry, colour = region),
        stat = "sf_coordinates",
        size = 2.2, fontface = "bold",
        segment.size = 0.2, segment.colour = "grey60",
        min.segment.length = 0,
        box.padding = 0.30, point.padding = 0.15,
        max.overlaps = 50, show.legend = FALSE) +
    geom_label(data = labels_df,
               aes(x = lon_lab, y = lat_lab,
                   label = ifelse(region == "Baker (HOR)", "Baker",
                                  as.character(region)),
                   color = region),
               fontface = "bold", size = 2.6,
               label.padding = unit(0.18, "lines"),
               label.r = unit(0.12, "lines"), label.size = 0.3,
               fill = alpha("white", 0.92), show.legend = FALSE) +
    scale_color_manual(values = region_colors) +
    coord_sf(xlim = c(-99, -82), ylim = c(61.5, 68.5), expand = FALSE) +
    annotation_scale(location = "br", width_hint = 0.25, line_width = 0.4,
                     height = unit(0.18, "cm"), text_cex = 0.6) +
    labs(title = "A   Sampling design",
         subtitle = "13 populations, 3 sub-regions, ~600 km",
         x = NULL, y = NULL) +
    theme_pub +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey85", linewidth = 0.2,
                                          linetype = "dotted"),
          axis.text = element_text(size = 7))

# --- PANEL B: Pairwise FST ---------------------------------------------------
cat("Loading data tables...\n")
fst <- read.csv(file.path(TABLE_DIR, "Table_IBD_pairwise.csv"))
fst$comparison_label <- factor(case_when(
    fst$Comparison == "Rankin Inlet" ~ "Within Rankin",
    fst$Comparison == "Naujaat"      ~ "Within Naujaat",
    fst$Comparison == "Inter-region" ~ "Between regions",
    TRUE                              ~ as.character(fst$Comparison)),
    levels = c("Within Rankin", "Within Naujaat", "Between regions"))
fst <- fst[is.finite(fst$Fst_weighted) & !is.na(fst$comparison_label), ]
comp_colors <- c("Within Rankin"   = "#D62728",
                 "Within Naujaat"  = "#0072B2",
                 "Between regions" = "#525252")
panel_B <- ggplot(fst, aes(x = comparison_label, y = Fst_weighted,
                           fill = comparison_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
    geom_jitter(width = 0.15, size = 1.0, alpha = 0.6,
                aes(color = comparison_label)) +
    scale_fill_manual(values = comp_colors) +
    scale_color_manual(values = comp_colors) +
    labs(title = expression(paste("B   Genetic differentiation (",
                                   italic(F)[ST], ")")),
         subtitle = "Within-region vs between-region pairs",
         x = NULL, y = expression(italic(F)[ST]~"(weighted)")) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8))

# --- PANEL C: nucleotide diversity -------------------------------------------
div <- read.csv(file.path(TABLE_DIR, "Table_Diversity.csv"))
div$Region <- clean_region(div$Region)
div <- div %>% arrange(Region, pi)
div$Pop <- factor(div$Pop, levels = div$Pop)
region_strip <- div %>% group_by(Region) %>%
    summarise(x_min = min(as.numeric(Pop)) - 0.5,
              x_max = max(as.numeric(Pop)) + 0.5,
              x_mid = mean(as.numeric(Pop)), .groups = "drop") %>%
    mutate(strip_label = ifelse(Region == "Baker (HOR)", "Baker",
                                 as.character(Region)))
y_max    <- max(div$pi * 1e3)
strip_y0 <- y_max * 1.06; strip_y1 <- y_max * 1.16; label_y <- y_max * 1.22
panel_C <- ggplot(div, aes(x = Pop, y = pi * 1e3, fill = Region)) +
    geom_col(alpha = 0.9, width = 0.8) +
    geom_rect(data = region_strip, inherit.aes = FALSE,
              aes(xmin = x_min, xmax = x_max, ymin = strip_y0, ymax = strip_y1,
                  fill = Region), alpha = 0.9) +
    geom_text(data = region_strip, inherit.aes = FALSE,
              aes(x = x_mid, y = label_y, label = strip_label),
              fontface = "bold", size = 3.0, color = "grey20") +
    geom_vline(data = region_strip[-1, ], aes(xintercept = x_min),
               color = "grey75", linewidth = 0.3, linetype = "dotted") +
    scale_fill_manual(values = region_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits = c(0, label_y * 1.05)) +
    labs(title = expression(paste("C   Nucleotide diversity (", pi, ")")),
         subtitle = "Per population, ordered within region",
         x = NULL, y = expression(pi%*%10^-3)) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# --- PANEL D: Atlantic ancestry (D-statistic, JAY config) --------------------
dstats <- read.delim(file.path(DSTATS_DIR, "Table_Dstats_JAY.tsv"))
dstats_obs <- dstats[dstats$Mode == "Observed", ]
dstats_pop <- dstats_obs[, c("H2","Region","D","SD")]
names(dstats_pop) <- c("Pop","Region","D","SD")
dstats_pop$Region <- clean_region(dstats_pop$Region)
dstats_pop <- dstats_pop[!is.na(dstats_pop$Region), ]
dstats_pop <- dstats_pop %>% arrange(Region, D)
dstats_pop$Pop <- factor(dstats_pop$Pop, levels = dstats_pop$Pop)
reg_strip_d <- dstats_pop %>% group_by(Region) %>%
    summarise(x_min = min(as.numeric(Pop)) - 0.5,
              x_max = max(as.numeric(Pop)) + 0.5,
              x_mid = mean(as.numeric(Pop)), .groups = "drop") %>%
    mutate(strip_label = ifelse(Region == "Baker (HOR)", "Baker",
                                 as.character(Region)))
y_max_d   <- max(dstats_pop$D + dstats_pop$SD) * 1.05
strip_y0d <- y_max_d * 1.06
strip_y1d <- y_max_d * 1.16
label_yd  <- y_max_d * 1.22
panel_D <- ggplot(dstats_pop, aes(x = Pop, y = D, fill = Region)) +
    geom_col(alpha = 0.9, width = 0.8) +
    geom_errorbar(aes(ymin = D - SD, ymax = D + SD),
                  width = 0.25, linewidth = 0.4, color = "grey30") +
    geom_rect(data = reg_strip_d, inherit.aes = FALSE,
              aes(xmin = x_min, xmax = x_max, ymin = strip_y0d, ymax = strip_y1d,
                  fill = Region), alpha = 0.9) +
    geom_text(data = reg_strip_d, inherit.aes = FALSE,
              aes(x = x_mid, y = label_yd, label = strip_label),
              fontface = "bold", size = 3.0, color = "grey20") +
    geom_vline(data = reg_strip_d[-1, ], aes(xintercept = x_min),
               color = "grey75", linewidth = 0.3, linetype = "dotted") +
    scale_fill_manual(values = region_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits = c(0, label_yd * 1.05)) +
    labs(title = "D   Atlantic ancestry (D-statistic)",
         subtitle = "Per population, ordered within region",
         x = NULL, y = "D-statistic (JAY, pop, LLS, DV)") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# --- PANEL E: Isolation-by-distance scatter panels (hydrographic) ------------
IBD_DIR <- file.path(BASE_DIR, "14_ibd")
ibd_pw <- read.csv(file.path(IBD_DIR, "Table_IBD_pairwise.csv"))
ibd_st <- read.csv(file.path(IBD_DIR, "Table_IBD_Mantel.csv"))
ibd_pw <- ibd_pw[is.finite(ibd_pw$Fst_lin) & is.finite(ibd_pw$Hydro_km), ]
make_ibd_data <- function(scope) {
    if (scope == "Global") {
        df <- if (!INCLUDE_HOR_IN_GLOBAL)
                  ibd_pw[ibd_pw$Pop1 != "HOR" & ibd_pw$Pop2 != "HOR", ]
              else ibd_pw
    } else if (scope == "Naujaat") {
        df <- ibd_pw[ibd_pw$Comparison == "Naujaat", ]
    } else if (scope == "Rankin") {
        df <- ibd_pw[ibd_pw$Comparison == "Rankin Inlet", ]
    }
    df$scope <- scope; df
}
get_mantel <- function(scope) {
    grp_lkp <- c("Global" = "Global", "Naujaat" = "Naujaat",
                 "Rankin" = "Rankin Inlet")
    row <- ibd_st[ibd_st$Group == grp_lkp[scope] & ibd_st$Distance == "Hydro_km", ]
    if (nrow(row) == 0) return(list(r = NA, p = NA))
    list(r = row$Mantel_r[1], p = row$Mantel_p[1])
}
fmt_p <- function(p) {
    if (is.na(p)) return("p = NA")
    if (p < 0.001) return("p < 0.001")
    sprintf("p = %.3f", p)
}
build_ibd_panel <- function(scope, color, title_txt) {
    df <- make_ibd_data(scope); mt <- get_mantel(scope)
    sig <- if (!is.na(mt$p)) {
        if (mt$p < 0.001) "***" else if (mt$p < 0.01) "**"
        else if (mt$p < 0.05) "*" else "n.s."
    } else "n.s."
    ann <- sprintf("Mantel r = %.2f\n%s   %s", mt$r, fmt_p(mt$p), sig)
    show_line <- !is.na(mt$p) && mt$p < 0.10
    p <- ggplot(df, aes(x = Hydro_km, y = Fst_lin)) +
        geom_point(color = color, alpha = 0.7, size = 1.6) +
        labs(title = title_txt,
             x = "Hydrographic distance (km)",
             y = expression(italic(F)[ST]~"/"~"("*1*"-"*italic(F)[ST]*")")) +
        theme_ibd
    if (show_line) p <- p + geom_smooth(method = "lm", se = TRUE,
                                        color = color, fill = color,
                                        linewidth = 0.6, alpha = 0.18)
    else p <- p + geom_hline(yintercept = mean(df$Fst_lin, na.rm = TRUE),
                             color = "grey60", linetype = "dashed",
                             linewidth = 0.4)
    p + annotate("text", x = -Inf, y = Inf, label = ann,
                 hjust = -0.08, vjust = 1.3,
                 size = 2.9, fontface = "bold", color = "grey15")
}
p_glob <- build_ibd_panel("Global",  "#525252", "Global (n = 13 pops)")
p_nauj <- build_ibd_panel("Naujaat", "#0072B2", "Naujaat (n = 7)")
p_rank <- build_ibd_panel("Rankin",  "#D62728", "Rankin (n = 5)")
panel_E_titled <- wrap_plots(p_glob, p_nauj, p_rank, nrow = 1) +
    plot_annotation(
        title    = "E   Isolation-by-distance",
        subtitle = "Linearized FST vs hydrographic distance, by scale",
        theme    = theme(plot.title    = element_text(face = "bold", size = 11,
                                                     hjust = 0),
                         plot.subtitle = element_text(size = 9, color = "grey35",
                                                     hjust = 0)))
panel_E_grob <- wrap_elements(full = panel_E_titled)

# --- PARTIAL MANTEL: ancestry vs distance (Naujaat, Global) -------------------
# Tests whether linearised FST tracks heterogeneity in Atlantic ancestry rather
# than hydrographic distance. Ancestry per pop is taken from the D-statistic
# (panel D source); a median hybrid index would be preferable if available.
suppressPackageStartupMessages(library(vegan))

anc_by_pop <- setNames(dstats_pop$D, as.character(dstats_pop$Pop))

.build_mat <- function(df, value_col, pops) {
    m <- matrix(NA_real_, length(pops), length(pops),
                dimnames = list(pops, pops))
    diag(m) <- 0
    for (i in seq_len(nrow(df))) {
        p1 <- df$Pop1[i]; p2 <- df$Pop2[i]
        if (p1 %in% pops && p2 %in% pops) {
            m[p1, p2] <- df[[value_col]][i]; m[p2, p1] <- df[[value_col]][i]
        }
    }
    m
}

run_partial <- function(scope, pops) {
    pops <- pops[pops %in% names(anc_by_pop)]
    if (length(pops) < 4) return(NULL)
    sub   <- ibd_pw[ibd_pw$Pop1 %in% pops & ibd_pw$Pop2 %in% pops, ]
    fst_m <- .build_mat(sub, "Fst_lin",  pops)
    hyd_m <- .build_mat(sub, "Hydro_km", pops)
    anc_m <- as.matrix(dist(anc_by_pop[pops]))          # |D_i - D_j|
    m_anc <- mantel.partial(as.dist(fst_m), as.dist(anc_m), as.dist(hyd_m),
                            permutations = 9999, na.rm = TRUE)
    m_dst <- mantel.partial(as.dist(fst_m), as.dist(hyd_m), as.dist(anc_m),
                            permutations = 9999, na.rm = TRUE)
    data.frame(
        Scope    = scope,
        N_pops   = length(pops),
        Test     = c("FST ~ ancestry | distance", "FST ~ distance | ancestry"),
        Mantel_r = round(c(m_anc$statistic, m_dst$statistic), 4),
        Mantel_p = round(c(m_anc$signif,   m_dst$signif),     4))
}

naujaat_pops_pm <- c("ITI","KGJ","NOP","PAM","SUP","TIN","WHI")
all_pops_pm     <- unique(c(ibd_pw$Pop1, ibd_pw$Pop2))
partial_res <- rbind(run_partial("Naujaat", naujaat_pops_pm),
                     run_partial("Global",  all_pops_pm))
cat("\n=== PARTIAL MANTEL (ancestry vs hydrographic distance) ===\n")
print(partial_res, row.names = FALSE)
write.csv(partial_res, file.path(IBD_DIR, "Table_partialMantel.csv"),
          row.names = FALSE)

# --- PANEL F: Management framework cards -------------------------------------
verdict <- tibble::tribble(
    ~Region,        ~Tier,                              ~Action,
    "Rankin",       "Management Unit (MU)",             "Single, demographically connected stock",
    "Naujaat",      "Management Unit (MU) + sub-MUs",   "Within-region heterogeneity by drainage",
    "Baker (HOR)",  "Management Unit (MU)",             "Distinct interior lineage, low diversity"
)
verdict$Region <- factor(verdict$Region, levels = region_levels)
verdict$row    <- 1:nrow(verdict)
panel_F <- ggplot(verdict) +
    geom_rect(aes(xmin = row - 0.42, xmax = row + 0.42,
                  ymin = 0.05, ymax = 0.95, fill = Region),
              alpha = 0.18, color = NA) +
    geom_rect(aes(xmin = row - 0.42, xmax = row - 0.36,
                  ymin = 0.05, ymax = 0.95, fill = Region),
              alpha = 1, color = NA) +
    geom_text(aes(x = row + 0.03, y = 0.78, label = Region, color = Region),
              fontface = "bold", size = 4.2, hjust = 0.5) +
    geom_text(aes(x = row + 0.03, y = 0.50, label = Tier),
              fontface = "bold", size = 3.3, color = "grey15", hjust = 0.5) +
    geom_text(aes(x = row + 0.03, y = 0.22, label = Action),
              fontface = "italic", size = 2.85, color = "grey35", hjust = 0.5) +
    scale_fill_manual(values = region_colors) +
    scale_color_manual(values = region_colors) +
    coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0, 1), clip = "off") +
    labs(title = "F   Recommended management framework",
         subtitle = "Convergent evidence supports tiered management units",
         x = NULL, y = NULL) +
    theme_pub +
    theme(axis.text = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank())

# --- ASSEMBLE ----------------------------------------------------------------
design <- "
ABC
DEE
FFF
"
combined <- panel_A + panel_B + panel_C +
            panel_D + panel_E_grob +
            panel_F +
    plot_layout(design = design, heights = c(1, 1, 0.55)) +
    plot_annotation(
        title    = "Candidate management units in the Kivalliq region of Nunavut",
        subtitle = "Convergent evidence from population genomics supports tiered management",
        theme    = theme(plot.title    = element_text(face = "bold", size = 14,
                                                     hjust = 0),
                         plot.subtitle = element_text(size = 11, color = "grey30",
                                                     hjust = 0)))

ggsave(file.path(FIG_DIR, "Fig_Management_Units.pdf"), combined,
       width = 13, height = 8.8, units = "in", dpi = 300)
ggsave(file.path(FIG_DIR, "Fig_Management_Units.png"), combined,
       width = 13, height = 8.8, units = "in", dpi = 300)

cat("\nSaved Fig_Management_Units to:", FIG_DIR, "\n")