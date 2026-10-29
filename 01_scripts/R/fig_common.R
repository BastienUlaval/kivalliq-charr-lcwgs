# =============================================================================
# fig_common.R — shared settings for the Heredity figures (Kivalliq charr)
#
# Sourced by 15_Fig1_map.R, 16_Fig2_global_structure.R,
# 17_Fig3_regional_structure.R and 18_Fig6_differentiation_diversity.R.
#
# Changes compared with the earlier scripts (JS comments on Draft 6):
#   - one GEOGRAPHIC population order used in every figure (c149)
#   - population codes printed under the admixture bars (c146)
#   - colour-blind-safe qualitative palette for ancestry clusters (c147)
#   - thin dark separators instead of thick white bars (c146: the two
#     "white bars" around Baker Lake were the region separators)
#   - lowercase panel letters, no figure titles (Heredity style)
#   - Heredity sizes: single column 86 mm, double column 178 mm
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
})

# --- Population order --------------------------------------------------------
# Interior first (Baker Lake), then the coast from south to north.
# Rankin Inlet: south-west to north-east. Naujaat: along the coast from
# White Island to the northernmost lake (PAM). Edit here only.
POP_ORDER_GEO <- c("HOR",
                   "CRB", "DIA", "MEL", "AUL", "AKL",
                   "WHI", "ITI", "TIN", "KGJ", "SUP", "NOP", "PAM")

POP_REGION <- c(HOR = "Baker Lake",
                CRB = "Rankin Inlet", DIA = "Rankin Inlet", MEL = "Rankin Inlet",
                AKL = "Rankin Inlet", AUL = "Rankin Inlet",
                WHI = "Naujaat", ITI = "Naujaat", TIN = "Naujaat", KGJ = "Naujaat",
                SUP = "Naujaat", NOP = "Naujaat", PAM = "Naujaat")
REGION_ORDER <- c("Baker Lake", "Rankin Inlet", "Naujaat")
# Short strip labels (Baker Lake has a single, narrow population)
ZONE_SHORT <- c("Baker Lake" = "Baker", "Rankin Inlet" = "Rankin Inlet", "Naujaat" = "Naujaat")

# --- Colours -----------------------------------------------------------------
REGION_COLORS <- c("Baker Lake"   = "#009E73",
                   "Rankin Inlet" = "#D55E00",
                   "Naujaat"      = "#0072B2")
REGION_SHAPES <- c("Baker Lake" = 17, "Rankin Inlet" = 16, "Naujaat" = 15)

# Population colours (Okabe-Ito / Tol based, unchanged from earlier figures)
POP_COLORS <- c(
    "HOR" = "#009E73",
    "PAM" = "#0072B2", "NOP" = "#56B4E9", "TIN" = "#4E79A7",
    "SUP" = "#9966CC", "KGJ" = "#332288", "WHI" = "#117733", "ITI" = "#44AA99",
    "AUL" = "#E69F00", "MEL" = "#D55E00", "DIA" = "#CC6677",
    "AKL" = "#882255", "CRB" = "#AA4499")

# Ancestry clusters: Paul Tol "bright" + "muted" (colour-blind safe, distinct
# in greyscale order). Up to 12 clusters.
CLUSTER_COLORS <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
                    "#AA3377", "#BBBBBB", "#332288", "#DDCC77", "#117733",
                    "#882255", "#44AA99")

# --- Sizes and theme ---------------------------------------------------------
MM <- 1 / 25.4
WIDTH_1COL <- 86 * MM
WIDTH_2COL <- 178 * MM

theme_hered <- function(base_size = 8) {
    theme_classic(base_size = base_size) +
        theme(axis.text  = element_text(colour = "black", size = base_size - 1),
              axis.title = element_text(size = base_size),
              axis.line  = element_line(linewidth = 0.3),
              axis.ticks = element_line(linewidth = 0.3),
              legend.title = element_text(size = base_size),
              legend.text  = element_text(size = base_size - 1),
              legend.key.size = unit(3, "mm"),
              plot.tag = element_text(face = "bold", size = base_size + 2),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold", size = base_size))
}

save_fig <- function(plot, dir, name, width, height) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(dir, paste0(name, ".pdf")), plot,
           width = width, height = height, units = "in", device = cairo_pdf, bg = "white")
    ggsave(file.path(dir, paste0(name, ".png")), plot,
           width = width, height = height, units = "in", dpi = 600, bg = "white")
    ggsave(file.path(dir, paste0(name, ".tiff")), plot,
           width = width, height = height, units = "in", dpi = 600,
           compression = "lzw", bg = "white")
    cat("Saved", name, "(pdf, png, tiff) to", dir, "\n")
}

# --- Sampling sites ----------------------------------------------------------
SITES <- data.frame(
    pop = c("AKL","AUL","CRB","DIA","MEL","HOR",
            "ITI","KGJ","NOP","PAM","TIN","SUP","WHI"),
    lon = c(-91.3071148, -91.3700786, -92.333333, -92.44444, -92.116667,
            -97.933631, -84.4652392, -85.254424, -86.883333, -87.055555,
            -84.4398478, -86.716667, -85.0426395),
    lat = c(62.8372049, 62.8118048, 62.466667, 62.833333, 62.866667,
            64.723358, 66.311737, 66.459357, 66.633333, 67.083333,
            66.3826439, 66.383333, 65.96663))
SITES$region <- factor(POP_REGION[SITES$pop], levels = REGION_ORDER)
SITES$pop    <- factor(SITES$pop, levels = POP_ORDER_GEO)

# --- Helpers -----------------------------------------------------------------
extract_pop <- function(bam_path) {
    toupper(sub(".*saal([A-Za-z]{3}).*", "\\1", basename(bam_path)))
}

parse_all_lnL <- function(log_dir, prefix) {
    logs <- list.files(log_dir, pattern = paste0("^", prefix, "_K.*\\.log$"),
                       full.names = TRUE)
    out <- lapply(logs, function(lf) {
        ll <- grep("best like=", readLines(lf, warn = FALSE), value = TRUE)
        if (!length(ll)) return(NULL)
        data.frame(K   = as.numeric(sub(".*_K(\\d+)_.*", "\\1", basename(lf))),
                   rep = as.numeric(sub(".*rep(\\d+).*", "\\1", basename(lf))),
                   lnL = as.numeric(sub(".*best like=([^ ]+).*", "\\1", ll[1])))
    })
    do.call(rbind, out)
}

best_qopt <- function(admix_dir, prefix, K) {
    lnL <- parse_all_lnL(admix_dir, prefix)
    if (is.null(lnL) || !any(lnL$K == K))
        stop("No NGSadmix logs for ", prefix, " K=", K, " in ", admix_dir)
    rep <- lnL$rep[lnL$K == K][which.max(lnL$lnL[lnL$K == K])]
    f <- file.path(admix_dir, sprintf("%s_K%d_rep%d.qopt", prefix, K, rep))
    cat("  ", prefix, " K=", K, " best replicate = ", rep, "\n", sep = "")
    f
}

load_pca <- function(prefix, bamlist) {
    pca <- read.table(paste0(prefix, ".cov.pca"), header = TRUE)
    bams <- scan(bamlist, what = "character", quiet = TRUE)
    if (nrow(pca) != length(bams))
        stop("PCA rows (", nrow(pca), ") != bamlist length (", length(bams), ")")
    pca$pop    <- factor(sapply(bams, extract_pop), levels = POP_ORDER_GEO)
    pca$region <- factor(POP_REGION[as.character(pca$pop)], levels = REGION_ORDER)
    list(data = droplevels(pca), eig = scan(paste0(prefix, ".cov.eig"), quiet = TRUE))
}

# PCA panel: points coloured by population, shaped by region, 95% ellipses
# either per region ("region") or per population ("pop").
plot_pca <- function(pca_obj, ellipse = c("region", "pop"), pc = c(1, 2)) {
    ellipse <- match.arg(ellipse)
    d <- pca_obj$data; eig <- pca_obj$eig
    d$x <- d[[paste0("PC", pc[1])]]; d$y <- d[[paste0("PC", pc[2])]]
    pops <- levels(d$pop); regs <- levels(d$region)
    p <- ggplot(d, aes(x, y))
    if (ellipse == "region") {
        p <- p + stat_ellipse(aes(group = region, colour = region), level = 0.95,
                              linewidth = 0.4, show.legend = FALSE)
    } else {
        keep <- names(which(table(d$pop) >= 4))
        p <- p + stat_ellipse(data = d[d$pop %in% keep, ],
                              aes(group = pop, colour = pop), level = 0.95,
                              linewidth = 0.3, show.legend = FALSE)
    }
    p + geom_point(aes(fill = pop, shape = region, colour = pop),
                   size = 1.3, alpha = 0.85, stroke = 0.2) +
        scale_colour_manual(values = c(POP_COLORS, REGION_COLORS),
                            breaks = pops, name = "Population") +
        scale_fill_manual(values = POP_COLORS, guide = "none") +
        scale_shape_manual(values = REGION_SHAPES[regs], name = "Zone") +
        labs(x = sprintf("PC%d (%.2f%%)", pc[1], eig[pc[1]]),
             y = sprintf("PC%d (%.2f%%)", pc[2], eig[pc[2]])) +
        guides(colour = guide_legend(override.aes = list(size = 2, shape = 16),
                                     ncol = 1, order = 1),
               shape = if (length(regs) > 1)
                   guide_legend(override.aes = list(size = 2), order = 2) else "none") +
        theme_hered() +
        theme(panel.border = element_rect(fill = NA, linewidth = 0.3))
}

# Palettes for the regional runs (Fig. 3): disjoint, so the two zones never
# share a colour (independent NGSadmix runs, clusters are not comparable).
PALETTE_RANKIN  <- c("#882255", "#DDCC77")
PALETTE_NAUJAAT <- c("#4477AA", "#EE6677", "#228833", "#66CCEE", "#AA3377")

# Prepare admixture data: geographic order, within-population sorting and a
# colour order for the clusters. Without `ref`, clusters are coloured by where
# they peak along the geographic axis. With `ref` (matrix returned by
# admix_colour_ref() for another K), each cluster takes the colour of the
# reference cluster it correlates with best, so colours match across K.
prep_admix <- function(qopt_file, bamlist, palette = CLUSTER_COLORS, ref = NULL) {
    Q <- as.matrix(read.table(qopt_file, header = FALSE))
    bams <- scan(bamlist, what = "character", quiet = TRUE)
    if (nrow(Q) != length(bams)) stop("qopt rows != bamlist length")
    K <- ncol(Q); cl <- paste0("C", seq_len(K)); colnames(Q) <- cl
    if (K > length(palette)) stop("Palette too short for K = ", K)
    df <- data.frame(row = seq_len(nrow(Q)),
                     pop = factor(sapply(bams, extract_pop), levels = POP_ORDER_GEO), Q)
    df$region <- factor(POP_REGION[as.character(df$pop)], levels = REGION_ORDER)
    df <- droplevels(df)
    df$dom  <- max.col(Q, ties.method = "first")
    df$domq <- apply(Q, 1, max)
    df <- df[order(df$pop, df$dom, -df$domq), ]
    df$x <- ave(seq_len(nrow(df)), df$region, FUN = seq_along)
    if (is.null(ref)) {
        peak <- sapply(cl, function(k) { w <- df[[k]]; sum(seq_len(nrow(df)) * w) / sum(w) })
        cl_order <- names(sort(peak))
    } else {
        # greedy best-pair matching on the correlation matrix (works whether
        # this K is larger or smaller than the reference K)
        cc <- cor(Q, ref); slot <- setNames(rep(NA_integer_, K), cl)
        while (any(is.na(slot)) && any(is.finite(cc))) {
            ij <- which(cc == max(cc, na.rm = TRUE), arr.ind = TRUE)[1, ]
            slot[rownames(cc)[ij[1]]] <- ij[2]
            cc[ij[1], ] <- NA; cc[, ij[2]] <- NA
        }
        left <- which(is.na(slot))
        slot[left] <- setdiff(seq_along(palette), slot)[seq_along(left)]
        cl_order <- names(sort(slot))
        cols <- setNames(palette[sort(slot)], cl_order)
        return(list(df = df, K = K, cl = cl, cl_order = cl_order, cols = cols,
                    qref = Q[, cl_order, drop = FALSE]))
    }
    list(df = df, K = K, cl = cl, cl_order = cl_order,
         cols = setNames(palette[seq_len(K)], cl_order),
         qref = Q[, cl_order, drop = FALSE])
}
admix_colour_ref <- function(qopt_file, bamlist, palette = CLUSTER_COLORS, ref = NULL)
    prep_admix(qopt_file, bamlist, palette, ref)$qref

# Admixture barplot ordered geographically, population codes below the bars,
# zones as facet strips on top, thin dark separators between populations.
plot_admix <- function(qopt_file, bamlist, show_regions = TRUE,
                       palette = CLUSTER_COLORS, ref = NULL) {
    a <- prep_admix(qopt_file, bamlist, palette, ref)
    df <- a$df
    long <- pivot_longer(df, all_of(a$cl), names_to = "cluster", values_to = "q")
    long$cluster <- factor(long$cluster, levels = rev(a$cl_order))
    info <- df %>% group_by(region, pop) %>%
        summarise(xmin = min(x) - 0.5, xmax = max(x) + 0.5, .groups = "drop") %>%
        mutate(xmid = (xmin + xmax) / 2)
    seps <- info %>% group_by(region) %>% filter(xmax < max(xmax))
    p <- ggplot(long, aes(x, q, fill = cluster)) +
        geom_col(width = 1, colour = NA) +
        geom_vline(data = seps, aes(xintercept = xmax), colour = "grey15",
                   linewidth = 0.25) +
        scale_fill_manual(values = a$cols, guide = "none") +
        scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
        scale_x_continuous(expand = c(0, 0)) +
        geom_text(data = info, aes(x = xmid, y = -0.04, label = pop),
                  inherit.aes = FALSE, angle = 90, hjust = 1, vjust = 0.5,
                  size = 2.2) +
        coord_cartesian(ylim = c(0, 1), clip = "off") +
        labs(x = NULL, y = "Ancestry proportion") +
        theme_hered() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              panel.border = element_rect(fill = NA, linewidth = 0.3),
              panel.spacing.x = unit(1.2, "mm"),
              plot.margin = margin(2, 2, 20, 2))
    if (show_regions) {
        p <- p + facet_grid(. ~ region, scales = "free_x", space = "free_x",
                            labeller = labeller(region = ZONE_SHORT)) +
            theme(strip.clip = "off")
    }
    p
}

# Base map layers from the Statistics Canada / CanVec shapefiles.
load_basemap <- function(shp_dir, xlim, ylim, tol = 0.005, rivers = FALSE) {
    suppressPackageStartupMessages(library(sf))
    sf_use_s2(FALSE)
    bb <- st_as_sfc(st_bbox(c(xmin = xlim[1], xmax = xlim[2],
                              ymin = ylim[1], ymax = ylim[2]), crs = 4326))
    rd <- function(f) {
        x <- st_read(file.path(shp_dir, f), quiet = TRUE)
        x <- st_transform(x, 4326)
        st_simplify(st_crop(x, bb), dTolerance = tol)
    }
    out <- list(land = rd("lpr_000b21a_e.shp"), lakes = rd("waterbody_2.shp"))
    if (rivers && file.exists(file.path(shp_dir, "watercourse_1.shp")))
        out$rivers <- rd("watercourse_1.shp")
    out
}

add_scalebar <- function(p, location = "br") {
    if (requireNamespace("ggspatial", quietly = TRUE)) {
        p + ggspatial::annotation_scale(location = location, line_width = 0.3,
                                        height = unit(1.5, "mm"), text_cex = 0.6)
    } else p
}
