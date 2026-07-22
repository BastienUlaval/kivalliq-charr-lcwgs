# =============================================================================
# ibd_analysis.R -- Isolation by Distance (4-panel: Global, Global-noHOR,
#   Naujaat, Rankin) with Mantel tests
#   Distance types: Euclidean (Haversine) + Hydrographic (from file or placeholder)
# Usage: Rscript ibd_analysis.R <fst_file> <info_dir> <fig_dir> <table_dir> [hydro_dist_file]
# =============================================================================

argv <- commandArgs(TRUE)
fst_file   <- argv[1]
info_dir   <- argv[2]
fig_dir    <- argv[3]
table_dir  <- argv[4]
hydro_file <- if (length(argv) >= 5) argv[5] else NULL

suppressPackageStartupMessages({
    library(ggplot2)
    library(vegan)
    library(dplyr)
    library(patchwork)
})

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. COORDINATES ---------------------------------------------------------
coords <- data.frame(
    Pop = c("AKL","AUL","CRB","DIA","MEL","HOR","ITI","KGJ","NOP","PAM","TIN","SUP","WHI"),
    Lat = c(62.8372049, 62.905473, 62.466667, 62.833333, 62.866667, 64.723358,
            66.311737, 66.459357, 66.633333, 67.083333, 66.3826439, 66.383333, 65.96663),
    Lon = c(-91.3071148, -91.15114147, -92.333333, -92.44444, -92.116667, -97.933631,
            -84.4652392, -85.254424, -86.883333, -87.055555, -84.4398478, -86.716667, -85.0426395)
)

# --- 2. REGION DEFINITIONS --------------------------------------------------
rankin_pops  <- c("AKL","AUL","CRB","DIA","MEL")
naujaat_pops <- c("ITI","KGJ","NOP","PAM","SUP","TIN","WHI")

get_region <- function(p) {
    if (p %in% rankin_pops) return("Rankin Inlet")
    if (p %in% naujaat_pops) return("Naujaat")
    return("Baker Lake")
}

# --- 3. HAVERSINE DISTANCE (km) ---------------------------------------------
haversine_km <- function(lat1, lon1, lat2, lon2) {
    R <- 6371  # Earth radius in km
    dlat <- (lat2 - lat1) * pi / 180
    dlon <- (lon2 - lon1) * pi / 180
    a <- sin(dlat/2)^2 + cos(lat1 * pi/180) * cos(lat2 * pi/180) * sin(dlon/2)^2
    return(2 * R * asin(sqrt(a)))
}

# --- 4. LOAD FST DATA -------------------------------------------------------
fst <- read.delim(fst_file, header = TRUE)
fst$Fst_lin <- fst$Fst_weighted / (1 - fst$Fst_weighted)

# Compute Euclidean (Haversine) distances
fst$Euclidean_km <- mapply(function(p1, p2) {
    c1 <- coords[coords$Pop == p1, ]
    c2 <- coords[coords$Pop == p2, ]
    if (nrow(c1) == 1 & nrow(c2) == 1) {
        return(haversine_km(c1$Lat, c1$Lon, c2$Lat, c2$Lon))
    }
    return(NA)
}, fst$Pop1, fst$Pop2)

# --- 5. HYDROGRAPHIC DISTANCES ----------------------------------------------
# If a hydro distance file is provided, load it
# Expected format: CSV with columns Pop1, Pop2, Hydro_km
if (!is.null(hydro_file) && file.exists(hydro_file)) {
    cat("Loading hydrographic distances from:", hydro_file, "\n")
    hydro <- read.csv(hydro_file)
    # Merge (both directions)
    fst$Hydro_km <- NA
    for (i in seq_len(nrow(fst))) {
        p1 <- fst$Pop1[i]; p2 <- fst$Pop2[i]
        h <- hydro[(hydro$Pop1 == p1 & hydro$Pop2 == p2) |
                    (hydro$Pop1 == p2 & hydro$Pop2 == p1), "Hydro_km"]
        if (length(h) > 0) fst$Hydro_km[i] <- h[1]
    }
    has_hydro <- sum(!is.na(fst$Hydro_km)) > 0
    cat("  Hydro distances loaded for", sum(!is.na(fst$Hydro_km)), "pairs\n")
} else {
    cat("No hydrographic distance file provided. Using Euclidean only.\n")
    fst$Hydro_km <- NA
    has_hydro <- FALSE
}

# Assign regions
fst$Region1 <- sapply(fst$Pop1, get_region)
fst$Region2 <- sapply(fst$Pop2, get_region)
fst$Comparison <- mapply(function(r1, r2) {
    if (r1 == r2) return(r1)
    return("Inter-region")
}, fst$Region1, fst$Region2)

# --- 6. MANTEL TEST FUNCTION ------------------------------------------------
run_mantel <- function(data_sub, dist_col, label) {
    pops_sub <- unique(c(data_sub$Pop1, data_sub$Pop2))
    n <- length(pops_sub)
    if (n < 3) return(data.frame(Group = label, Distance = dist_col,
                                  N_pops = n, N_pairs = nrow(data_sub),
                                  Mantel_r = NA, Mantel_p = NA, R2 = NA))

    # Build matrices
    fst_mat <- matrix(0, n, n, dimnames = list(pops_sub, pops_sub))
    geo_mat <- matrix(0, n, n, dimnames = list(pops_sub, pops_sub))

    for (i in seq_len(nrow(data_sub))) {
        p1 <- data_sub$Pop1[i]; p2 <- data_sub$Pop2[i]
        fst_mat[p1, p2] <- data_sub$Fst_lin[i]
        fst_mat[p2, p1] <- data_sub$Fst_lin[i]
        geo_mat[p1, p2] <- data_sub[[dist_col]][i]
        geo_mat[p2, p1] <- data_sub[[dist_col]][i]
    }

    m <- mantel(as.dist(fst_mat), as.dist(geo_mat), permutations = 9999)
    r2 <- summary(lm(data_sub$Fst_lin ~ data_sub[[dist_col]]))$r.squared

    return(data.frame(Group = label, Distance = dist_col,
                      N_pops = n, N_pairs = nrow(data_sub),
                      Mantel_r = round(m$statistic, 4),
                      Mantel_p = round(m$signif, 4),
                      R2 = round(r2, 4)))
}

# --- 7. RUN MANTEL TESTS ----------------------------------------------------
fst_complete <- fst[!is.na(fst$Euclidean_km), ]
fst_no_hor <- fst_complete[fst_complete$Pop1 != "HOR" & fst_complete$Pop2 != "HOR", ]
fst_rankin <- fst_complete[fst_complete$Region1 == "Rankin Inlet" & fst_complete$Region2 == "Rankin Inlet", ]
fst_naujaat <- fst_complete[fst_complete$Region1 == "Naujaat" & fst_complete$Region2 == "Naujaat", ]

results <- rbind(
    run_mantel(fst_complete, "Euclidean_km", "Global"),
    run_mantel(fst_no_hor, "Euclidean_km", "Global (sans HOR)"),
    run_mantel(fst_rankin, "Euclidean_km", "Rankin Inlet"),
    run_mantel(fst_naujaat, "Euclidean_km", "Naujaat")
)

# Hydro Mantel tests if available
if (has_hydro) {
    fst_hydro <- fst[!is.na(fst$Hydro_km), ]
    fst_hydro_no_hor <- fst_hydro[fst_hydro$Pop1 != "HOR" & fst_hydro$Pop2 != "HOR", ]
    fst_hydro_rankin <- fst_hydro[fst_hydro$Region1 == "Rankin Inlet" & fst_hydro$Region2 == "Rankin Inlet", ]
    fst_hydro_naujaat <- fst_hydro[fst_hydro$Region1 == "Naujaat" & fst_hydro$Region2 == "Naujaat", ]

    results <- rbind(results,
        run_mantel(fst_hydro, "Hydro_km", "Global"),
        run_mantel(fst_hydro_no_hor, "Hydro_km", "Global (sans HOR)"),
        run_mantel(fst_hydro_rankin, "Hydro_km", "Rankin Inlet"),
        run_mantel(fst_hydro_naujaat, "Hydro_km", "Naujaat")
    )
}

cat("\n=== MANTEL TEST RESULTS ===\n")
print(results)
write.csv(results, file.path(table_dir, "Table_IBD_Mantel.csv"), row.names = FALSE)

# --- 8. PLOT HELPER ---------------------------------------------------------
make_ibd_panel <- function(data_sub, dist_col, title_text, point_color,
                           mantel_r, mantel_p, show_y_lab = TRUE) {

    dist_label <- if (dist_col == "Hydro_km") {
        "Hydrographic distance (km)"
    } else {
        "Euclidean distance (km)"
    }

    # Significance annotation
    sig <- if (is.na(mantel_p)) "NA"
           else if (mantel_p < 0.001) "***"
           else if (mantel_p < 0.01) "**"
           else if (mantel_p < 0.05) "*"
           else "ns"

    r2 <- if (nrow(data_sub) >= 3) {
        round(summary(lm(data_sub$Fst_lin ~ data_sub[[dist_col]]))$r.squared, 3)
    } else NA

    annot_label <- paste0("r = ", round(mantel_r, 3),
                          " (p = ", formatC(mantel_p, format = "f", digits = 4), ") ", sig,
                          "\nR² = ", r2)

    p <- ggplot(data_sub, aes(x = .data[[dist_col]], y = Fst_lin)) +
        geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed",
                    fill = "grey85", linewidth = 0.8) +
        geom_point(size = 2.5, alpha = 0.85, color = point_color) +
        annotate("text",
                 x = min(data_sub[[dist_col]], na.rm = TRUE),
                 y = max(data_sub$Fst_lin, na.rm = TRUE),
                 label = annot_label,
                 hjust = 0, vjust = 1, size = 3.2, fontface = "italic") +
        labs(title = title_text,
             x = dist_label,
             y = if (show_y_lab) expression(F[ST] / (1 - F[ST])) else NULL) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
        )

    return(p)
}

# --- 9. GENERATE FIGURES ----------------------------------------------------

# Colors for comparison types
comparison_colors <- c("Rankin Inlet" = "#31688E", "Naujaat" = "#35B779",
                       "Baker Lake" = "#E69F00", "Inter-region" = "grey55")

# Special panel builder for Global views (color by comparison type)
make_ibd_panel_global <- function(data_sub, dist_col, title_text,
                                   mantel_r, mantel_p, show_y_lab = TRUE) {

    dist_label <- if (dist_col == "Hydro_km") {
        "Hydrographic distance (km)"
    } else {
        "Euclidean distance (km)"
    }

    sig <- if (is.na(mantel_p)) "NA"
           else if (mantel_p < 0.001) "***"
           else if (mantel_p < 0.01) "**"
           else if (mantel_p < 0.05) "*"
           else "ns"

    r2 <- if (nrow(data_sub) >= 3) {
        round(summary(lm(data_sub$Fst_lin ~ data_sub[[dist_col]]))$r.squared, 3)
    } else NA

    annot_label <- paste0("r = ", round(mantel_r, 3),
                          " (p = ", formatC(mantel_p, format = "f", digits = 4), ") ", sig,
                          "\nR² = ", r2)

    p <- ggplot(data_sub, aes(x = .data[[dist_col]], y = Fst_lin, color = Comparison)) +
        geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed",
                    fill = "grey85", linewidth = 0.8, inherit.aes = FALSE,
                    aes(x = .data[[dist_col]], y = Fst_lin)) +
        geom_point(size = 2.5, alpha = 0.85) +
        scale_color_manual(values = comparison_colors, name = "Comparison") +
        annotate("text",
                 x = min(data_sub[[dist_col]], na.rm = TRUE),
                 y = max(data_sub$Fst_lin, na.rm = TRUE),
                 label = annot_label,
                 hjust = 0, vjust = 1, size = 3.2, fontface = "italic") +
        labs(title = title_text,
             x = dist_label,
             y = if (show_y_lab) expression(F[ST] / (1 - F[ST])) else NULL) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 9),
            legend.text = element_text(size = 8)
        )

    return(p)
}

# Helper to extract Mantel stats for a group + distance type
get_stats <- function(group, dist) {
    row <- results[results$Group == group & results$Distance == dist, ]
    if (nrow(row) == 0) return(list(r = NA, p = NA))
    return(list(r = row$Mantel_r, p = row$Mantel_p))
}

# --- Euclidean 4-panel ------------------------------------------------------
s1 <- get_stats("Global", "Euclidean_km")
s2 <- get_stats("Global (sans HOR)", "Euclidean_km")
s3 <- get_stats("Naujaat", "Euclidean_km")
s4 <- get_stats("Rankin Inlet", "Euclidean_km")

p1 <- make_ibd_panel_global(fst_complete, "Euclidean_km", "Global",
                              s1$r, s1$p, TRUE)
p2 <- make_ibd_panel_global(fst_no_hor, "Euclidean_km", "Global (sans HOR)",
                              s2$r, s2$p, FALSE)
p3 <- make_ibd_panel(fst_naujaat, "Euclidean_km", "Naujaat",
                      "#35B779", s3$r, s3$p, TRUE)
p4 <- make_ibd_panel(fst_rankin, "Euclidean_km", "Rankin Inlet",
                      "#31688E", s4$r, s4$p, FALSE)

p_eucl <- (p1 | p2) / (p3 | p4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Isolation by Distance -- Euclidean (Haversine)",
                    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))) &
    theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Fig_IBD_euclidean.png"), p_eucl,
       width = 12, height = 9, dpi = 300)

# Individual panels
ggsave(file.path(fig_dir, "Fig_IBD_euclidean_global.png"), p1,
       width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig_IBD_euclidean_global_noHOR.png"), p2,
       width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig_IBD_euclidean_naujaat.png"), p3,
       width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig_IBD_euclidean_rankin.png"), p4,
       width = 6, height = 5, dpi = 300)

cat("Euclidean IBD figures saved (combined + 4 individual).\n")

# --- Hydrographic 4-panel (if available) -------------------------------------
if (has_hydro) {
    fst_hydro <- fst[!is.na(fst$Hydro_km), ]
    fst_hydro_no_hor <- fst_hydro[fst_hydro$Pop1 != "HOR" & fst_hydro$Pop2 != "HOR", ]
    fst_hydro_rankin <- fst_hydro[fst_hydro$Region1 == "Rankin Inlet" & fst_hydro$Region2 == "Rankin Inlet", ]
    fst_hydro_naujaat <- fst_hydro[fst_hydro$Region1 == "Naujaat" & fst_hydro$Region2 == "Naujaat", ]

    h1 <- get_stats("Global", "Hydro_km")
    h2 <- get_stats("Global (sans HOR)", "Hydro_km")
    h3 <- get_stats("Naujaat", "Hydro_km")
    h4 <- get_stats("Rankin Inlet", "Hydro_km")

    ph1 <- make_ibd_panel_global(fst_hydro, "Hydro_km", "Global",
                                   h1$r, h1$p, TRUE)
    ph2 <- make_ibd_panel_global(fst_hydro_no_hor, "Hydro_km", "Global (sans HOR)",
                                   h2$r, h2$p, FALSE)
    ph3 <- make_ibd_panel(fst_hydro_naujaat, "Hydro_km", "Naujaat",
                           "#35B779", h3$r, h3$p, TRUE)
    ph4 <- make_ibd_panel(fst_hydro_rankin, "Hydro_km", "Rankin Inlet",
                           "#31688E", h4$r, h4$p, FALSE)

    p_hydro <- (ph1 | ph2) / (ph3 | ph4) +
        plot_annotation(title = "Isolation by Distance -- Hydrographic",
                        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

    ggsave(file.path(fig_dir, "Fig_IBD_hydrographic.png"), p_hydro,
           width = 12, height = 9, dpi = 300)
    cat("Hydrographic IBD figure saved.\n")
}

# --- Save full pairwise data ------------------------------------------------
write.csv(fst, file.path(table_dir, "Table_IBD_pairwise.csv"), row.names = FALSE)

cat("=== IBD ANALYSIS COMPLETE ===\n")