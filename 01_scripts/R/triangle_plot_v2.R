#!/usr/bin/env Rscript
# =============================================================================
# triangle_plot_v3.R - lcWGS triangle plot (HI vs IH) - publication version
#
# CHANGE FROM v2 (2026-07-15):
#   The full triangle and the zoomed scatter are no longer separate figures.
#   Fig_Triangle_combined is now a SINGLE composite:
#       - main panel : the ZOOMED HI vs IH scatter (where the biology is)
#       - inset      : the full triangle with the F1/F2/BC envelope, small,
#                      top-left of the main panel. Its only job is to show that
#                      no individual sits anywhere near F1/F2 - i.e. no
#                      early-generation hybrids.
#       - bottom     : the hybrid-index boxplot by population.
#   The standalone Fig_Triangle_full / Fig_Triangle_HI_IH_scatter are still
#   written (useful for supplementary / talks) but the combined figure is the
#   manuscript deliverable.
#
# OUTPUT FIGURES:
#   Fig_Triangle_combined.pdf/png        - publication figure: zoom + inset + boxplot
#   Fig_Triangle_full.pdf/png            - standalone full triangle (supplementary)
#   Fig_Triangle_HI_IH_scatter.pdf/png   - standalone zoomed scatter (supplementary)
#   Fig_Triangle_by_pop.pdf/png          - standalone boxplot (supplementary)
#
# OUTPUT TABLES:
#   Table_HybridIndex.csv                - per-individual HI, IH
#   Table_HybridIndex_by_pop.csv         - per-population summary
#   Table_HybridIndex_by_region.csv      - per-region summary
#
# NOTES:
#   - "IH" is the proportion of diagnostic sites where an individual is
#     heterozygous (one Arctic + one Atlantic allele). Following Wang &
#     Schumer (2023) and current charr literature, we use the term
#     "interancestry heterozygosity" rather than "interspecific" since
#     the Arctic and Atlantic lineages are not formally described as
#     distinct species.
#   - Unified region palette + shapes at top of file - edit to swap globally.
#   - Region = shape + colour (works in B&W and for colour-blind readers).
#   - Diagnostic SNPs are restricted to the ngsParalog-canonical (non-deviant)
#     set via CANON_SITES (see section 2). Set CANON_SITES="" to disable.
#
# Usage: Rscript triangle_plot_v3.R <triangle_dir> <popmap_file> <fig_dir>
#   Optional env var: CANON_SITES=/path/to/sites_..._canonical
#
# popmap_file format: <bam>\t<pop>\t<region>
# region values expected: Rankin, Naujaat, Baker, Arctic_anchor, Atlantic_anchor
# =============================================================================

argv <- commandArgs(TRUE)
if (length(argv) < 3) {
    stop("Usage: Rscript triangle_plot_v3.R <triangle_dir> <popmap_file> <fig_dir>")
}
triangle_dir <- argv[1]
popmap_file  <- argv[2]
fig_dir      <- argv[3]
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
})

# =============================================================================
# 1. CONFIGURATION (edit here to change palette / shapes / sizes)
# =============================================================================

# --- Colour palette: hybrid colour-blind-friendly --------------------------- #
# Edit these 5 hex codes to swap the whole figure-set's palette.
REGION_COLORS <- c(
    "Rankin"           = "#F00000",   # orange  (Rankin Inlet, Kivalliq)
    "Naujaat"          = "#0072B2",   # blue    (Naujaat,      Kivalliq)
    "Baker"            = "#009E73",   # green   (Baker Lake = HOR, Kivalliq)
    "Arctic_anchor"    = "#000000",   # black   (JAY, Canadian Arctic outgroup)
    "Atlantic_anchor"  = "#D55E00"    # red     (LLS, Sweden Atlantic outgroup)
)

# --- Shapes: filled for Kivalliq, open/distinctive for anchors -------------- #
REGION_SHAPES <- c(
    "Rankin"           = 16,   # filled circle
    "Naujaat"          = 15,   # filled square
    "Baker"            = 17,   # filled triangle
    "Arctic_anchor"    = 4,    # cross (X) - clearly different from Kivalliq
    "Atlantic_anchor"  = 8     # asterisk-star - clearly different
)

REGION_SIZES <- c(
    "Rankin"           = 1.6,
    "Naujaat"          = 1.6,
    "Baker"            = 1.6,
    "Arctic_anchor"    = 2.4,
    "Atlantic_anchor"  = 2.4
)

REGION_LABELS <- c(
    "Rankin"           = "Rankin Inlet",
    "Naujaat"          = "Naujaat",
    "Baker"            = "Baker Lake (HOR)",
    "Arctic_anchor"    = "JAY (Arctic anchor)",
    "Atlantic_anchor"  = "LLS (Atlantic anchor)"
)

# Plot draw order (last = on top). Anchors on top to be visible.
REGION_ORDER <- c("Naujaat", "Rankin", "Baker", "Arctic_anchor", "Atlantic_anchor")

# X-axis ordering for the by-population panel (b): group by region in this
# left-to-right order, then sort populations by median HI WITHIN each region
# (so Naujaat and Rankin are no longer interleaved). Anchors sit at the
# extremes: JAY (Arctic, HI~0) far left, LLS (Atlantic, HI~1) far right.
REGION_X_ORDER <- c("Arctic_anchor", "Baker", "Naujaat", "Rankin", "Atlantic_anchor")

# --- Inset geometry (full triangle placed inside the zoomed panel) ---------- #
# Fractions of the main panel, measured from its bottom-left corner.
INSET_LEFT   <- 0.015
INSET_BOTTOM <- 0.58
INSET_RIGHT  <- 0.375
INSET_TOP    <- 0.995

# Strict-diagnostic thresholds (overridable via env vars)
FIXED_HI <- as.numeric(Sys.getenv("FIXED_HI", "0.85"))
FIXED_LO <- as.numeric(Sys.getenv("FIXED_LO", "0.15"))

# =============================================================================
# 2. LOAD DATA AND COMPUTE HI / IH
# =============================================================================

beagle_file <- file.path(triangle_dir, "triangle.beagle.gz")
diag_file   <- file.path(triangle_dir, "diagnostic_snps.tsv")

cat("=== Triangle plot v3 ===\n")
cat("  Beagle      :", beagle_file, "\n")
cat("  Diagnostic  :", diag_file,   "\n")
cat("  Popmap      :", popmap_file, "\n")
cat("  Output dir  :", fig_dir,     "\n")
cat("  Strict diag : max(freq) >", FIXED_HI, "AND min(freq) <", FIXED_LO, "\n\n")

# --- Load diagnostic SNP table (Atlantic vs Arctic anchor frequencies) ------ #
# Expected cols: chromo (or chr), position (or pos), major, minor,
#                freq_LLS (MAF), freq_JAY (MAF), [deltaAF or dAF]
diag <- fread(diag_file)
# Normalise column names
if ("chr" %in% names(diag) && !"chromo" %in% names(diag))      setnames(diag, "chr", "chromo")
if ("pos" %in% names(diag) && !"position" %in% names(diag))    setnames(diag, "pos", "position")
if ("dAF" %in% names(diag) && !"deltaAF" %in% names(diag))     setnames(diag, "dAF", "deltaAF")
required <- c("chromo", "position", "major", "minor", "freq_LLS", "freq_JAY")
missing  <- setdiff(required, names(diag))
if (length(missing) > 0) {
    stop("diagnostic_snps.tsv missing columns: ", paste(missing, collapse = ", "))
}
diag[, marker := paste(chromo, position, sep = "_")]
cat("  Loaded", nrow(diag), "candidate diagnostic SNPs\n")

# Strict filter: one anchor near-fixed major, the other near-fixed minor
diag <- diag[ pmax(freq_LLS, freq_JAY) > FIXED_HI &
              pmin(freq_LLS, freq_JAY) < FIXED_LO ]
cat("  After strict filter:", nrow(diag), "SNPs retained\n")
if (nrow(diag) < 100) {
    warning("Fewer than 100 strict-diagnostic SNPs; consider relaxing FIXED_HI/FIXED_LO")
}

# --- Deviant (paralog) masking: keep only canonical (non-deviant) sites ----- #
# Restricts diagnostic SNPs to the ngsParalog-canonical set so paralogous loci
# (which segregate near 0.5 and inflate interancestry heterozygosity) are
# excluded, consistent with the rest of the pipeline. The beagle is matched to
# diag$marker below, so this automatically subsets the plotted SNPs (no ANGSD
# rerun needed). Override the path with:
#   CANON_SITES=/path/to/sites_..._canonical Rscript triangle_plot_v2.R ...
# If the file is absent (or CANON_SITES=""), masking is skipped with a warning
# and the figure still builds on the strict-diagnostic set.
canon_file <- Sys.getenv("CANON_SITES", "")
deviant_masked <- nzchar(canon_file) && file.exists(canon_file)
if (deviant_masked) {
    canon <- fread(canon_file, header = FALSE, select = 1:2,
                   col.names = c("chromo", "position"))
    canon_markers <- paste(canon$chromo, canon$position, sep = "_")
    n_before <- nrow(diag)
    diag <- diag[ marker %in% canon_markers ]
    cat(sprintf("  After deviant masking: %d SNPs retained (%d removed as deviant)\n",
                nrow(diag), n_before - nrow(diag)))
    if (nrow(diag) < 100) {
        warning("Fewer than 100 SNPs after deviant masking; check CANON_SITES path.")
    }
} else {
    warning("CANON_SITES not found (", canon_file,
            "); skipping deviant masking. Set CANON_SITES to enable.")
    cat("  Deviant masking SKIPPED (canonical file not found)\n")
}

# Label used in figure subtitles to reflect whether masking was applied
diag_label <- if (deviant_masked) "strict-diagnostic, non-deviant" else "strict-diagnostic"
cat("\n")

# Polarize: "Atlantic allele" = minor when freq_LLS > freq_JAY (since MAF given)
diag[, atl_is_minor := freq_LLS > freq_JAY]
cat("  SNPs where Atlantic = minor allele:", sum(diag$atl_is_minor), "\n")
cat("  SNPs where Atlantic = major allele:", sum(!diag$atl_is_minor), "\n\n")

# --- Load beagle genotype likelihoods -------------------------------------- #
cat("Loading beagle (this may take a minute)...\n")
beagle <- fread(beagle_file)
orig_cols <- names(beagle)
n_ind  <- (ncol(beagle) - 3L) %/% 3L
cat("  Beagle: ", nrow(beagle), "sites,", n_ind, "individuals\n")

# Match beagle rows to filtered diagnostic SNPs.
# Beagle's first column is "marker" (chr_pos format from ANGSD).
beagle_markers <- beagle[[1]]
match_idx <- match(beagle_markers, diag$marker)
keep <- !is.na(match_idx)
cat("  Beagle sites matched to diagnostic:", sum(keep), "\n\n")
if (sum(keep) < 100) stop("Too few matched diagnostic sites - check chr name formats in beagle vs diagnostic_snps.tsv")

# atl_is_min_vec is aligned to the rows of beagle (after subsetting to keep)
atl_is_min_vec <- diag$atl_is_minor[match_idx[keep]]
beagle <- beagle[keep, ]
# Column order is preserved; first 3 cols are (marker, allele1, allele2) header,
# then GL trios per individual
stopifnot(identical(names(beagle), orig_cols))
gc()

# --- Load popmap ------------------------------------------------------------ #
popmap <- fread(popmap_file, header = FALSE)
if (ncol(popmap) == 3) {
    setnames(popmap, c("bam", "pop", "region"))
    popmap[, ind_id := sub("\\.bam.*$", "", basename(bam))]
} else if (ncol(popmap) == 4) {
    setnames(popmap, c("bam", "ind_id", "pop", "region"))
} else {
    stop("Popmap must have 3 (bam/pop/region) or 4 (bam/ind/pop/region) columns")
}
if (nrow(popmap) != n_ind) {
    stop(sprintf("Popmap n=%d but beagle has %d individuals", nrow(popmap), n_ind))
}
cat("Popmap:\n")
print(popmap[, .N, by = region])
cat("\n")

# --- Compute HI and IH per individual --------------------------------------- #
HI <- numeric(n_ind); IH <- numeric(n_ind); N <- integer(n_ind)
EQUAL_TOL <- 1e-3

cat("Computing HI / IH for", n_ind, "individuals x", nrow(beagle), "SNPs...\n")
chunk <- 200L
for (start in seq(1L, nrow(beagle), by = chunk)) {
    end  <- min(start + chunk - 1L, nrow(beagle))
    rows <- start:end
    atl_is_min <- atl_is_min_vec[rows]
    for (i in seq_len(n_ind)) {
        c0 <- 4L + 3L * (i - 1L)   # P(MajMaj)
        c1 <- c0 + 1L              # P(MajMin)
        c2 <- c0 + 2L              # P(MinMin)
        pmm <- beagle[[c0]][rows]
        pmh <- beagle[[c1]][rows]
        phh <- beagle[[c2]][rows]
        # informative if any GL deviates from uniform 1/3
        inf <- abs(pmm - 1/3) > EQUAL_TOL | abs(pmh - 1/3) > EQUAL_TOL |
               abs(phh - 1/3) > EQUAL_TOL
        # Expected Atlantic allele dosage
        atl_dos <- ifelse(atl_is_min, 2 * phh + pmh, 2 * pmm + pmh)
        het_p   <- pmh
        atl_dos[!inf] <- 0
        het_p[!inf]   <- 0
        HI[i] <- HI[i] + sum(atl_dos)
        IH[i] <- IH[i] + sum(het_p)
        N[i]  <- N[i]  + sum(inf)
    }
}

ind_df <- data.frame(
    ind_id = popmap$ind_id,
    bam    = popmap$bam,
    pop    = popmap$pop,
    region = popmap$region,
    n_diag = N,
    HI     = ifelse(N > 0, HI / (2 * N), NA_real_),
    IH     = ifelse(N > 0, IH / N,       NA_real_)
)

# Drop individuals with < 10% of diagnostic SNPs informative
min_snp <- max(50L, round(0.10 * nrow(beagle)))
n_drop  <- sum(ind_df$n_diag < min_snp, na.rm = TRUE)
cat("  Dropping", n_drop, "individuals with <", min_snp, "informative SNPs\n")
ind_df <- ind_df[!is.na(ind_df$n_diag) & ind_df$n_diag >= min_snp, ]
cat("  Retained", nrow(ind_df), "individuals\n\n")

# Make region a factor in the chosen draw order
ind_df$region <- factor(ind_df$region, levels = REGION_ORDER)

# =============================================================================
# 3. SAVE TABLES
# =============================================================================

write.csv(ind_df, file.path(fig_dir, "Table_HybridIndex.csv"), row.names = FALSE)

pop_df <- ind_df %>%
    group_by(region, pop) %>%
    summarise(
        n = n(),
        HI_mean = mean(HI), HI_sd = sd(HI), HI_median = median(HI),
        HI_min = min(HI),   HI_max = max(HI),
        IH_mean = mean(IH), IH_sd = sd(IH),
        .groups = "drop"
    ) %>%
    arrange(HI_mean)
write.csv(pop_df, file.path(fig_dir, "Table_HybridIndex_by_pop.csv"), row.names = FALSE)

region_df <- ind_df %>%
    group_by(region) %>%
    summarise(
        n = n(),
        HI_mean = mean(HI), HI_sd = sd(HI), HI_median = median(HI),
        HI_min = min(HI),   HI_max = max(HI),
        IH_mean = mean(IH), IH_sd = sd(IH),
        .groups = "drop"
    )
write.csv(region_df, file.path(fig_dir, "Table_HybridIndex_by_region.csv"), row.names = FALSE)

cat("Per-region HI summary:\n"); print(region_df); cat("\n")
cat("Per-population HI summary (sorted by mean HI):\n"); print(pop_df, n = Inf); cat("\n")

# =============================================================================
# 4. FULL TRIANGLE (Fitzpatrick 2012 style, with F1/F2/BC anchors)
#    Built by a helper so the same geometry can be rendered twice:
#      - "standalone" : full-size supplementary figure, with legend + titles
#      - "inset"      : stripped-down miniature for the combined figure
# =============================================================================

# Theoretical anchor points for the hybrid envelope:
#   Pure parental Arctic    : (0,   0)
#   Pure parental Atlantic  : (1,   0)
#   F1 hybrid               : (0.5, 1)
#   F2 hybrid               : (0.5, 0.5)
#   BC (backcross) Arctic   : (0.25, 0.5)
#   BC Atlantic             : (0.75, 0.5)
# Envelope: IH <= 2*HI  AND  IH <= 2*(1-HI)
theory_pts <- data.frame(
    x = c(0, 1, 0.5, 0.5, 0.25, 0.75),
    y = c(0, 0, 1.0, 0.5, 0.5,  0.5),
    label = c("Arctic", "P_Atlan", "F1", "F2", "BC_Arctic", "BC_Atlantic")
)
# Envelope vertices (Arctic - F1 - Atlantic - Arctic)
envelope <- data.frame(
    x = c(0,   0.5, 1,   0),
    y = c(0,   1.0, 0,   0)
)

# Zoom window on the main panel, computed here so the inset can draw the same
# rectangle and show the reader exactly which sliver of the triangle is expanded.
plot_df <- ind_df[ind_df$region != "Atlantic_anchor", ]
plot_df$region <- droplevels(plot_df$region)
PAD_X <- 0.02; PAD_Y <- 0.02
x_min <- max(0, min(plot_df$HI, na.rm = TRUE) - PAD_X)
x_max <-        max(plot_df$HI, na.rm = TRUE) + PAD_X
y_min <- max(0, min(plot_df$IH, na.rm = TRUE) - PAD_Y)
y_max <-        max(plot_df$IH, na.rm = TRUE) + PAD_Y
zoom_box <- data.frame(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max)

make_full_triangle <- function(mode = c("standalone", "inset")) {
    mode <- match.arg(mode)
    is_inset <- identical(mode, "inset")

    # Inset: small points, no theoretical-anchor text labels except F1/F2,
    # no axis titles, no legend. The inset exists ONLY to show that the cloud
    # sits nowhere near F1/F2 - it is evidence of absence, not a data display.
    pt_scale   <- if (is_inset) 0.55 else 1.6
    anchor_lab <- if (is_inset) theory_pts[theory_pts$label %in% c("F1", "F2"), ] else theory_pts
    base_size  <- if (is_inset) 8 else 12

    p <- ggplot(ind_df, aes(HI, IH, colour = region, shape = region, size = region)) +
        # Triangle envelope
        geom_path(data = envelope, aes(x, y), inherit.aes = FALSE,
                  linetype = "dashed", colour = "grey55", linewidth = 0.4) +
        # Horizontal Arctic-Atlantic baseline
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), inherit.aes = FALSE,
                     linetype = "dashed", colour = "grey55", linewidth = 0.4) +
        # Theoretical anchor crosses
        geom_point(data = theory_pts, aes(x, y), inherit.aes = FALSE,
                   shape = 4, size = if (is_inset) 1.6 else 2.5,
                   colour = "grey55", stroke = 0.5) +
        geom_text(data = anchor_lab, aes(x, y, label = label), inherit.aes = FALSE,
                  size = if (is_inset) 2.5 else 3, colour = "grey45",
                  nudge_y = 0.045, fontface = "italic") +
        geom_point(alpha = 0.78, stroke = 0.4) +
        scale_colour_manual(values = REGION_COLORS, labels = REGION_LABELS,
                            name = "Region / Lineage") +
        scale_shape_manual( values = REGION_SHAPES, labels = REGION_LABELS,
                            name = "Region / Lineage") +
        scale_size_manual(  values = REGION_SIZES * pt_scale, guide = "none")

    if (is_inset) {
        # Red rectangle marking the region expanded in the main panel
        p <- p +
            geom_rect(data = zoom_box,
                      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                      inherit.aes = FALSE, fill = NA, colour = "#D62728",
                      linewidth = 0.5) +
            scale_x_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1),
                               labels = c("0", "0.5", "1"), expand = expansion(0)) +
            scale_y_continuous(limits = c(-0.02, 1.10), breaks = c(0, 0.5, 1),
                               expand = expansion(0)) +
            labs(x = "Hybrid Index", y = "Interancestry het.", title = NULL,
                 subtitle = NULL) +
            theme_classic(base_size = base_size) +
            theme(
                legend.position  = "none",
                panel.grid       = element_blank(),
                plot.background  = element_rect(fill = alpha("white", 0.92),
                                                colour = "grey30", linewidth = 0.4),
                plot.margin      = margin(3, 4, 2, 2),
                axis.title       = element_text(size = 6.5, colour = "grey30"),
                axis.text        = element_text(size = 6, colour = "grey30"),
                axis.ticks       = element_line(linewidth = 0.25),
                axis.line        = element_line(linewidth = 0.3)
            )
    } else {
        p <- p +
            scale_x_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.25),
                               labels = c("0\n(Arctic)", "0.25", "0.50", "0.75",
                                          "1\n(Atlantic)"),
                               expand = expansion(0)) +
            scale_y_continuous(limits = c(-0.02, 1.08), breaks = seq(0, 1, 0.25),
                               expand = expansion(0)) +
            labs(
                x = "Hybrid Index (proportion of Atlantic alleles)",
                y = "Interancestry heterozygosity",
                title    = "Hybrid index vs interancestry heterozygosity",
                subtitle = sprintf(
                    "%d individuals, %d %s SNPs   \u00b7   Dashed: theoretical hybrid envelope (Fitzpatrick 2012)",
                    nrow(ind_df), nrow(beagle), diag_label)
            ) +
            theme_classic(base_size = base_size) +
            theme(
                plot.title    = element_text(face = "bold"),
                plot.subtitle = element_text(colour = "grey40", size = 9.5),
                panel.grid = element_blank(),
                legend.position = "right"
            ) +
            guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
                   shape  = guide_legend(override.aes = list(size = 3, alpha = 1)))
    }
    p
}

cat("Generating Fig_Triangle_full (standalone, supplementary)...\n")
p_triangle_full  <- make_full_triangle("standalone")
p_triangle_inset <- make_full_triangle("inset")

ggsave(file.path(fig_dir, "Fig_Triangle_full.pdf"), p_triangle_full,
       width = 8.5, height = 6.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig_Triangle_full.png"), p_triangle_full,
       width = 8.5, height = 6.5, dpi = 300)
cat("  Saved.\n")

# =============================================================================
# 5. ZOOMED HI vs IH SCATTER (LLS excluded so the axes are not compressed)
#    This is the MAIN panel of the combined figure.
# =============================================================================

cat("Generating Fig_Triangle_HI_IH_scatter (zoomed)...\n")

# Subset palette/shapes/sizes/labels for this plot
keep_regions <- levels(plot_df$region)
col_sub  <- REGION_COLORS[keep_regions]
shp_sub  <- REGION_SHAPES[keep_regions]
sz_sub   <- REGION_SIZES[keep_regions]
lab_sub  <- REGION_LABELS[keep_regions]

make_zoom_scatter <- function(with_inset_space = FALSE) {
    # When the inset is overlaid we widen the y-range slightly so the inset
    # does not sit on top of any data point.
    y_hi <- if (with_inset_space) y_max + (y_max - y_min) * 0.42 else y_max

    ggplot(plot_df, aes(HI, IH, colour = region, shape = region, size = region)) +
        geom_point(alpha = 0.75, stroke = 0.6) +
        scale_colour_manual(values = col_sub, labels = lab_sub, name = "Region / Lineage") +
        scale_shape_manual( values = shp_sub, labels = lab_sub, name = "Region / Lineage") +
        scale_size_manual(  values = sz_sub * 1.6, guide = "none") +
        scale_x_continuous(limits = c(x_min, x_max),
                           breaks = pretty(c(x_min, x_max), n = 5)) +
        scale_y_continuous(limits = c(y_min, y_hi),
                           breaks = pretty(c(y_min, y_max), n = 5)) +
        labs(
            x = "Hybrid Index (proportion of Atlantic alleles)",
            y = "Interancestry heterozygosity",
            title    = "Hybrid index vs interancestry heterozygosity",
            subtitle = sprintf(
                "%d individuals; LLS anchor excluded (HI~1) for axis readability  \u00b7  %d %s SNPs",
                nrow(plot_df), nrow(beagle), diag_label)
        ) +
        theme_classic(base_size = 12) +
        theme(
            plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(colour = "grey40", size = 9.5),
            panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
            legend.position = "right"
        ) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
               shape  = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

p_scatter <- make_zoom_scatter(with_inset_space = FALSE)

ggsave(file.path(fig_dir, "Fig_Triangle_HI_IH_scatter.pdf"), p_scatter,
       width = 8, height = 5.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig_Triangle_HI_IH_scatter.png"), p_scatter,
       width = 8, height = 5.5, dpi = 300)
cat("  Saved.\n")

# =============================================================================
# 6. HI BY POPULATION (bottom panel of the combined figure)
# =============================================================================

cat("Generating Fig_Triangle_by_pop (boxplot grouped by region, then median HI)...\n")

# Order populations BY REGION FIRST (REGION_X_ORDER), then by median HI within
# each region, so Naujaat and Rankin populations are not interleaved.
pop_order <- ind_df %>%
    group_by(pop, region) %>%
    summarise(med_HI = median(HI), .groups = "drop") %>%
    mutate(region_rank = match(as.character(region), REGION_X_ORDER)) %>%
    arrange(region_rank, med_HI) %>%
    pull(pop)
ind_df$pop <- factor(ind_df$pop, levels = pop_order)

# Build region "bands" for vertical shading behind the boxplot. With the
# region-first ordering above, each region now occupies a contiguous block.
band_df <- ind_df %>%
    group_by(pop, region) %>%
    summarise(.groups = "drop") %>%
    arrange(pop) %>%
    mutate(x_int = as.integer(pop)) %>%
    group_by(region) %>%
    summarise(xmin = min(x_int) - 0.5,
              xmax = max(x_int) + 0.5,
              x_mid = (min(x_int) + max(x_int)) / 2,
              n_pops = n_distinct(pop[!is.na(pop)]),
              .groups = "drop") %>%
    mutate(label = REGION_LABELS[as.character(region)])

make_by_pop <- function(show_title = TRUE) {
    ggplot() +
        # Region bands behind the boxplot (numeric x positions).
        geom_tile(data = band_df,
                  aes(x = x_mid, y = 0.55,
                      width = xmax - xmin, height = 1.20,
                      fill = region),
                  alpha = 0.12) +
        scale_fill_manual(values = REGION_COLORS, guide = "none") +
        # Boxplots and jitter use discrete pop on x; ggplot promotes the integer
        # x of the tile layer to a discrete scale that aligns with these.
        geom_boxplot(data = ind_df,
                     aes(x = pop, y = HI, colour = region),
                     outlier.shape = NA, fill = "white",
                     linewidth = 0.4, width = 0.65, alpha = 0.95) +
        geom_jitter(data = ind_df,
                    aes(x = pop, y = HI, colour = region, shape = region),
                    width = 0.16, size = 1.5, alpha = 0.75, stroke = 0.4) +
        scale_colour_manual(values = REGION_COLORS, labels = REGION_LABELS,
                            name = "Region / Lineage") +
        scale_shape_manual( values = REGION_SHAPES, labels = REGION_LABELS,
                            name = "Region / Lineage") +
        # Force a discrete x scale based on pop level order
        scale_x_discrete(limits = levels(ind_df$pop)) +
        scale_y_continuous(limits = c(-0.02, 1.12), breaks = seq(0, 1, 0.25),
                           labels = c("0\n(Arctic)", "0.25", "0.50", "0.75",
                                      "1\n(Atlantic)"),
                           expand = expansion(0)) +
        labs(x = NULL,
             y = "Hybrid Index (proportion of Atlantic alleles)",
             title = if (show_title) "Hybrid index distribution by population" else NULL,
             subtitle = if (show_title)
                 "Populations grouped by region, then ordered by median HI within region"
                 else NULL) +
        theme_classic(base_size = 12) +
        theme(
            plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(colour = "grey40", size = 9.5),
            axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
            panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
            legend.position = "right"
        ) +
        guides(colour = guide_legend(override.aes = list(size = 3)),
               shape  = guide_legend(override.aes = list(size = 3)))
}

p_pop <- make_by_pop(show_title = TRUE)

ggsave(file.path(fig_dir, "Fig_Triangle_by_pop.pdf"), p_pop,
       width = 9.5, height = 5.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig_Triangle_by_pop.png"), p_pop,
       width = 9.5, height = 5.5, dpi = 300)
cat("  Saved.\n")

# =============================================================================
# 7. COMBINED PUBLICATION FIGURE
#    panel a = zoomed HI vs IH scatter, with the full triangle as a top-left
#              inset (red box on the inset marks the zoomed window)
#    panel b = hybrid index by population, below
# =============================================================================

cat("Generating combined publication panel (zoom + inset triangle + by-pop)...\n")

# Main panel: extra headroom so the inset never overlaps a data point
p_zoom_main <- make_zoom_scatter(with_inset_space = TRUE) +
    labs(subtitle = sprintf(
        "%d individuals; LLS anchor excluded (HI~1) for axis readability  \u00b7  %d %s SNPs  \u00b7  Inset: full triangle, red box = zoomed window",
        nrow(plot_df), nrow(beagle), diag_label))

# Drop the miniature full triangle into the top-left of the main panel.
# align_to = "panel" pins it to the plotting area, not the whole ggplot object,
# so the legend on the right does not shift it.
p_zoom_with_inset <- p_zoom_main +
    inset_element(p_triangle_inset,
                  left   = INSET_LEFT,
                  bottom = INSET_BOTTOM,
                  right  = INSET_RIGHT,
                  top    = INSET_TOP,
                  align_to = "panel")

# Bottom panel: no title/subtitle (the composite tag carries the numbering)
p_pop_panel <- make_by_pop(show_title = FALSE)

combined <- p_zoom_with_inset / p_pop_panel +
    plot_layout(heights = c(1.35, 1)) +
    plot_annotation(tag_levels = "a")

ggsave(file.path(fig_dir, "Fig_Triangle_combined.pdf"), combined,
       width = 10, height = 11, device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig_Triangle_combined.png"), combined,
       width = 10, height = 11, dpi = 300)
cat("  Saved.\n")

cat("\n=== triangle_plot_v3.R DONE ===\n")
cat("All figures + tables in:", fig_dir, "\n")