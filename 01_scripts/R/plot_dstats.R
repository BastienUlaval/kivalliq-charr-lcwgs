#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# plot_dstats.R (v4) — D-statistics plots, with and without significance stars
#
# Usage:
#   Rscript plot_dstats.R <dstat_dir> <fig_dir> [<h1_label>]
#
# Reads:
#   {POP}_vs_{H1}_H3-{H3}_H4-{H4}_Dstat.Observed.txt
#   {POP}_vs_{H1}_H3-{H3}_H4-{H4}_Dstat.TransRem.txt
#
# Outputs in <fig_dir> (TWO versions of each main figure):
#   Fig_Dstats_<H1>_forest.pdf/png            — without stars
#   Fig_Dstats_<H1>_forest_stars.pdf/png      — with stars
#   Fig_Dstats_<H1>_barplot.pdf/png           — without stars
#   Fig_Dstats_<H1>_barplot_stars.pdf/png     — with stars
#   Fig_Dstats_<H1>_combined.pdf              — forest + bar (no stars)
#   Fig_Dstats_<H1>_combined_stars.pdf        — forest + bar (with stars)
#   Fig_Dstats_<H1>_obs_vs_notrans.pdf/png    — Observed vs No Trans
#   Table_Dstats_<H1>.tsv
#
# Significance based on Z-score:
#   *** p < 0.001  ** p < 0.01  * p < 0.05
#
# Author: Bastien Rubin + Claude (2026)
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(patchwork)
    library(scales)
})

# ──────────────────────────────────────────────────────────────────────────────
# 1. Configuration
# ──────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript plot_dstats.R <dstat_dir> <fig_dir> [<h1_label>]")
}

DSTAT_DIR <- args[1]
FIG_DIR   <- args[2]
H1_LABEL  <- ifelse(length(args) >= 3, args[3], "JAY")

if (!dir.exists(DSTAT_DIR)) stop("DSTAT_DIR does not exist: ", DSTAT_DIR)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Plot D-stats (v4: with and without stars) ===\n")
cat("  DSTAT_DIR:", DSTAT_DIR, "\n")
cat("  FIG_DIR:  ", FIG_DIR, "\n")
cat("  H1 label: ", H1_LABEL, "\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Population metadata
# ──────────────────────────────────────────────────────────────────────────────
POP_REGIONS <- c(
    AKL = "Rankin", AUL = "Rankin", CRB = "Rankin",
    DIA = "Rankin", MEL = "Rankin",
    ITI = "Naujaat", KGJ = "Naujaat", NOP = "Naujaat",
    PAM = "Naujaat", SUP = "Naujaat", TIN = "Naujaat",
    WHI = "Naujaat",
    HOR = "Baker"
)

REGION_COLORS <- c(
    Rankin  = "#C0392B",
    Naujaat = "#2E86C1",
    Baker   = "#229954"
)

REGION_FILL <- c(
    Rankin  = "#E74C3C",
    Naujaat = "#3498DB",
    Baker   = "#27AE60"
)

# Ordre : Baker sera en haut, Rankin en bas
REGION_ORDER <- c("Baker", "Naujaat", "Rankin")

# ──────────────────────────────────────────────────────────────────────────────
# 3. Helpers
# ──────────────────────────────────────────────────────────────────────────────
parse_dstat_file <- function(f, mode_label) {
    df <- tryCatch(
        read.table(f, header = TRUE, sep = "\t", check.names = FALSE,
                   stringsAsFactors = FALSE, comment.char = ""),
        error = function(e) {
            warning("Cannot read ", basename(f), ": ", conditionMessage(e))
            return(NULL)
        }
    )
    if (is.null(df) || nrow(df) == 0 || ncol(df) != 12) return(NULL)
    names(df) <- c("D", "JKD", "VJKD", "Z", "pvalue",
                   "nABBA", "nBABA", "nBlocks", "H1", "H2", "H3", "H4")
    df <- df[1, , drop = FALSE]
    df$Mode <- mode_label
    df$SD <- sqrt(df$VJKD)
    df
}

sig_marker <- function(p) {
    ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < 0.05,  "*", ""))))
}

# Determine sig stars from p-value or fallback to Z if p == 0
sig_from_pz <- function(p, z) {
    out <- sig_marker(p)
    fb <- is.na(out) | out == ""
    out[fb] <- ifelse(abs(z[fb]) >= 3.29, "***",
               ifelse(abs(z[fb]) >= 2.58, "**",
               ifelse(abs(z[fb]) >= 1.96, "*", "")))
    out
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Load all D-stat files
# ──────────────────────────────────────────────────────────────────────────────
obs_pattern <- paste0("^.+_vs_", H1_LABEL, "_H3-.+_H4-.+_Dstat\\.Observed\\.txt$")
obs_files <- list.files(DSTAT_DIR, pattern = obs_pattern, full.names = TRUE)

if (length(obs_files) == 0) {
    stop("No D-stat files found for H1=", H1_LABEL, " in ", DSTAT_DIR,
         "\n  Expected: {POP}_vs_", H1_LABEL, "_H3-{X}_H4-{Y}_Dstat.Observed.txt")
}

cat("Found", length(obs_files), "Observed files for H1 =", H1_LABEL, "\n")

results_list <- list()
for (obs_f in obs_files) {
    tr_f <- sub("_Dstat\\.Observed\\.txt$", "_Dstat.TransRem.txt", obs_f)
    obs <- parse_dstat_file(obs_f, "Observed")
    if (is.null(obs)) next
    if (file.exists(tr_f)) {
        tr <- parse_dstat_file(tr_f, "No Trans")
        if (!is.null(tr)) {
            results_list[[length(results_list) + 1]] <- rbind(obs, tr)
        } else {
            results_list[[length(results_list) + 1]] <- obs
        }
    } else {
        results_list[[length(results_list) + 1]] <- obs
    }
}

if (length(results_list) == 0) stop("No results parsed successfully.")

results <- do.call(rbind, results_list)
results$Region  <- POP_REGIONS[results$H2]

# Drop unknown pops
unknown_pops <- unique(results$H2[is.na(results$Region)])
if (length(unknown_pops) > 0) {
    cat("WARNING: dropping unknown populations: ",
        paste(unknown_pops, collapse = ", "), "\n", sep = "")
    results <- results[!is.na(results$Region), ]
}

results$CI_low      <- results$D - 1.96 * results$SD
results$CI_high     <- results$D + 1.96 * results$SD
results$Significant <- abs(results$Z) > 3
results$Sig         <- sig_from_pz(results$pvalue, results$Z)

cat("Parsed", length(unique(results$H2)), "populations\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 5. Save table
# ──────────────────────────────────────────────────────────────────────────────
table_path <- file.path(FIG_DIR, paste0("Table_Dstats_", H1_LABEL, ".tsv"))
results %>%
    select(H1, H2, H3, H4, Region, Mode, D, SD, CI_low, CI_high, Z, pvalue, Sig,
           nABBA, nBABA, nBlocks) %>%
    arrange(Region, desc(D), Mode) %>%
    write_tsv(table_path)
cat("✓ Table written:", table_path, "\n")

# ──────────────────────────────────────────────────────────────────────────────
# 6. Plotting data
# ──────────────────────────────────────────────────────────────────────────────
# 1. On garde l'ordre des régions (Baker en haut, Rankin en bas)
REGION_ORDER <- c("Baker", "Naujaat", "Rankin")

# 2. Préparation des données
plot_data <- results %>%
    filter(Mode == "Observed") %>%
    # CHANGEMENT ICI : On trie par D croissant (pas de 'desc')
    # Les plus petits D de chaque région arrivent en premier dans le dataframe
    arrange(factor(Region, levels = REGION_ORDER), D)

# 3. On capture cet ordre (plus petit D -> plus grand D)
pop_order <- unique(plot_data$H2)

# 4. On transforme en facteurs pour ggplot
plot_data <- plot_data %>%
    mutate(
        Region = factor(Region, levels = REGION_ORDER),
        # On utilise rev() car ggplot dessine le niveau 1 en bas de l'axe Y.
        # En inversant, le premier élément (plus petit D) se retrouve au sommet.
        H2 = factor(H2, levels = rev(pop_order))
    ) # <--- Parenthèse fermée, le script ne plantera plus ici !
# ──────────────────────────────────────────────────────────────────────────────
# 7. Theme
# ──────────────────────────────────────────────────────────────────────────────
theme_paper <- theme_minimal(base_size = 13, base_family = "sans") +
    theme(
        panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(color = "gray95", linewidth = 0.3),
        axis.line.x        = element_line(color = "gray30", linewidth = 0.5),
        axis.ticks.x       = element_line(color = "gray30"),
        axis.text          = element_text(color = "gray20"),
        axis.title         = element_text(color = "gray10", face = "bold"),
        legend.position    = "top",
        legend.title       = element_text(face = "bold", size = 11),
        legend.text        = element_text(size = 11),
        legend.background  = element_rect(fill = "white", color = NA),
        legend.key         = element_blank(),
        plot.title         = element_text(face = "bold", size = 14, hjust = 0,
                                          margin = margin(b = 8)),
        plot.subtitle      = element_text(size = 11, color = "gray30",
                                          margin = margin(b = 14)),
        plot.caption       = element_text(size = 9, color = "gray50",
                                          hjust = 0, margin = margin(t = 10)),
        plot.margin        = margin(15, 20, 15, 15)
    )

# ──────────────────────────────────────────────────────────────────────────────
# 8. Forest plot — function with stars on/off toggle
# ──────────────────────────────────────────────────────────────────────────────
make_forest <- function(data, with_stars = FALSE) {
    if (with_stars) {
        data$annot <- sprintf("Z=%.1f %s", data$Z, data$Sig)
        cap <- "Bars: 95% CI (D ± 1.96 × SE).  *** p < 0.001    ** p < 0.01    * p < 0.05  (Z-score test)."
    } else {
        data$annot <- sprintf("Z=%.1f", data$Z)
        cap <- "Bars: 95% CI (D ± 1.96 × SE).  Z = jackknife Z-score."
    }

    ggplot(data, aes(x = D, y = H2, color = Region, fill = Region)) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   color = "gray40", linewidth = 0.6) +
        geom_segment(aes(x = CI_low, xend = CI_high, y = H2, yend = H2),
                     linewidth = 1.0, alpha = 0.7) +
        geom_point(shape = 21, size = 4.5, stroke = 1.0, color = "gray20") +
        geom_text(aes(x = CI_high + 0.015, label = annot),
                  hjust = 0, size = 3.2, color = "gray35") +
        scale_color_manual(values = REGION_COLORS, name = "Region") +
        scale_fill_manual(values  = REGION_FILL,   name = "Region") +
        scale_x_continuous(
            limits = c(min(data$CI_low) - 0.05,
                       max(data$CI_high) + ifelse(with_stars, 0.13, 0.10)),
            breaks = seq(0, 0.5, 0.1),
            expand = c(0, 0)
        ) +
        labs(
            title    = bquote(bolditalic(D)~bold("-statistics with")~bolditalic(.(H1_LABEL))~bold("as H1")),
            subtitle = bquote(italic("D")~"("*.(H1_LABEL)*", H2 ; LLS, DV)  —  "*
                              italic("D")~">"~"0 indicates Atlantic introgression in H2"),
            x        = expression(italic("D") * "-statistic (Observed)"),
            y        = NULL,
            caption  = cap
        ) +
        theme_paper +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(color = "gray92", linewidth = 0.3),
            axis.text.y        = element_text(face = "bold", size = 12)
        )
}

# ──────────────────────────────────────────────────────────────────────────────
# 9. Barplot horizontal — function with stars on/off toggle
# ──────────────────────────────────────────────────────────────────────────────
make_barplot <- function(data, with_stars = FALSE) {
    if (with_stars) {
        data$label_bar <- sprintf("%.3f%s", data$D,
                                  ifelse(data$Sig != "", paste0(" ", data$Sig), ""))
        cap <- "Error bars: 95% CI.  *** p < 0.001    ** p < 0.01    * p < 0.05  (Z-score test)."
    } else {
        data$label_bar <- sprintf("%.3f", data$D)
        cap <- "Bars sorted by D-value within region.  Error bars: 95% CI."
    }

    ggplot(data, aes(x = D, y = H2, fill = Region)) +
        geom_vline(xintercept = 0, linetype = "solid",
                   color = "gray40", linewidth = 0.5) +
        geom_col(width = 0.65, color = NA, alpha = 0.92) +
        geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = Region),
                      width = 0.25, linewidth = 0.7, alpha = 0.8) +
        geom_text(aes(x = D + ifelse(D > 0, 0.012, -0.012), label = label_bar),
                  hjust = ifelse(data$D > 0, 0, 1),
                  size = 3.3, color = "gray15", fontface = "bold") +
        scale_fill_manual(values = REGION_FILL, name = "Region") +
        scale_color_manual(values = REGION_COLORS, guide = "none") +
        scale_x_continuous(
            limits = c(min(data$CI_low, 0) - 0.08,
                       max(data$CI_high) + ifelse(with_stars, 0.12, 0.08)),
            breaks = seq(0, 0.5, 0.1),
            expand = c(0, 0)
        ) +
        labs(
            title    = bquote(bolditalic(D)~bold("-statistics with")~bolditalic(.(H1_LABEL))~bold("as H1")),
            subtitle = bquote(italic("D")~"("*.(H1_LABEL)*", H2 ; LLS, DV)"),
            x        = expression(italic("D") * "-statistic (Observed)"),
            y        = NULL,
            caption  = cap
        ) +
        theme_paper +
        theme(
            panel.grid.major.y = element_blank(),
            axis.text.y        = element_text(face = "bold", size = 12)
        )
}

# ──────────────────────────────────────────────────────────────────────────────
# 10. Build both versions
# ──────────────────────────────────────────────────────────────────────────────
forest_no    <- make_forest(plot_data,  with_stars = FALSE)
forest_stars <- make_forest(plot_data,  with_stars = TRUE)
bar_no       <- make_barplot(plot_data, with_stars = FALSE)
bar_stars    <- make_barplot(plot_data, with_stars = TRUE)

combined_no    <- forest_no    + bar_no    + plot_layout(guides = "collect", widths = c(1, 1)) & theme(legend.position = "top")
combined_stars <- forest_stars + bar_stars + plot_layout(guides = "collect", widths = c(1, 1)) & theme(legend.position = "top")

# ──────────────────────────────────────────────────────────────────────────────
# 11. Comparison Observed vs No Trans (only one version, no stars needed)
# ──────────────────────────────────────────────────────────────────────────────
comparison_plot <- NULL
if (any(results$Mode == "No Trans")) {
    both_modes <- results %>%
        filter(Mode %in% c("Observed", "No Trans")) %>%
        mutate(
            Region = factor(Region, levels = REGION_ORDER),
            Mode   = factor(Mode, levels = c("Observed", "No Trans")),
            H2     = factor(H2, levels = pop_order)
        ) %>%
        filter(!is.na(H2))

    if (length(unique(both_modes$Mode)) == 2) {
        cat("✓ Both Observed and No Trans available — making comparison plot\n")
        comparison_plot <- ggplot(both_modes,
                                  aes(x = D, y = H2, color = Region, shape = Mode)) +
            geom_vline(xintercept = 0, linetype = "dashed",
                       color = "gray40", linewidth = 0.6) +
            geom_segment(aes(x = CI_low, xend = CI_high, y = H2, yend = H2),
                         linewidth = 0.7, alpha = 0.5,
                         position = position_dodge(width = 0.5)) +
            geom_point(size = 3.5, stroke = 1.2,
                       position = position_dodge(width = 0.5)) +
            scale_color_manual(values = REGION_COLORS, name = "Region") +
            scale_shape_manual(values = c("Observed" = 16, "No Trans" = 1),
                               name = "Mode") +
            labs(
                title    = bquote(bolditalic(D)~bold("-statistics: Observed vs No Trans")),
                subtitle = bquote(italic("D")~"("*.(H1_LABEL)*", H2 ; LLS, DV)"),
                x        = expression(italic("D") * "-statistic"),
                y        = NULL,
                caption  = "Closed: all sites (Observed). Open: transversions only (No Trans)."
            ) +
            theme_paper +
            theme(
                panel.grid.major.y = element_blank(),
                axis.text.y        = element_text(face = "bold", size = 12)
            )
    }
}

# ──────────────────────────────────────────────────────────────────────────────
# 12. Save plots — both versions
# ──────────────────────────────────────────────────────────────────────────────
n_pops <- length(unique(plot_data$H2))
plot_h <- max(4.5, 0.45 * n_pops + 2.5)

cat("\nSaving plots...\n")

save_plot <- function(p, basename, width, height) {
    pdf_out <- file.path(FIG_DIR, paste0(basename, ".pdf"))
    png_out <- file.path(FIG_DIR, paste0(basename, ".png"))
    ggsave(pdf_out, p, width = width, height = height, device = cairo_pdf)
    cat("  ✓", pdf_out, "\n")
    ggsave(png_out, p, width = width, height = height, dpi = 300)
    cat("  ✓", png_out, "\n")
}

# Forest
save_plot(forest_no,    paste0("Fig_Dstats_", H1_LABEL, "_forest"),       9,  plot_h)
save_plot(forest_stars, paste0("Fig_Dstats_", H1_LABEL, "_forest_stars"), 10, plot_h)

# Bar
save_plot(bar_no,    paste0("Fig_Dstats_", H1_LABEL, "_barplot"),       9,  plot_h)
save_plot(bar_stars, paste0("Fig_Dstats_", H1_LABEL, "_barplot_stars"), 10, plot_h)

# Combined (PDF only — these are large)
comb_pdf_no    <- file.path(FIG_DIR, paste0("Fig_Dstats_", H1_LABEL, "_combined.pdf"))
comb_pdf_stars <- file.path(FIG_DIR, paste0("Fig_Dstats_", H1_LABEL, "_combined_stars.pdf"))
ggsave(comb_pdf_no,    combined_no,    width = 16, height = plot_h, device = cairo_pdf)
ggsave(comb_pdf_stars, combined_stars, width = 17, height = plot_h, device = cairo_pdf)
cat("  ✓", comb_pdf_no,    "\n")
cat("  ✓", comb_pdf_stars, "\n")

# Comparison Observed vs No Trans
if (!is.null(comparison_plot)) {
    save_plot(comparison_plot,
              paste0("Fig_Dstats_", H1_LABEL, "_obs_vs_notrans"),
              10, plot_h)
}

# ──────────────────────────────────────────────────────────────────────────────
# 13. Summary
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== Summary (Observed mode) ===\n")
summary_tbl <- plot_data %>%
    select(Region, H2, D, SD, Z, Sig) %>%
    mutate(across(c(D, SD), ~ round(., 4)),
           Z = round(Z, 1)) %>%
    arrange(Region, desc(D))
print(summary_tbl)
cat("\n✓ Done.\n")