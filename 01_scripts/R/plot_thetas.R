# =============================================================================
# plot_thetas.R — Publication-quality diversity summary
#   Dot plots with 95% CI from sliding windows, axes in ×10⁻³
#   Reads both per-scaffold .pestPG and per-window .thetaswindow.pestPG
# Usage: Rscript plot_thetas.R <theta_dir> <pop_file> <fig_dir> <table_dir>
# =============================================================================

argv <- commandArgs(TRUE)
theta_dir <- argv[1]
pop_file  <- argv[2]
fig_dir   <- argv[3]
table_dir <- argv[4]

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
})

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

pops <- scan(pop_file, what = "character", quiet = TRUE)

pop_order <- c("AKL","AUL","CRB","DIA","MEL","HOR","ITI","KGJ","NOP","PAM","SUP","TIN","WHI")
pop_order <- intersect(pop_order, pops)

region_map <- c(
    AKL="Rankin", AUL="Rankin", CRB="Rankin", DIA="Rankin", MEL="Rankin",
    HOR="Baker",
    ITI="Naujaat", KGJ="Naujaat", NOP="Naujaat", PAM="Naujaat",
    SUP="Naujaat", TIN="Naujaat", WHI="Naujaat")

region_colors <- c(Rankin = "#31688E", Baker = "#E69F00", Naujaat = "#35B779")

# ─── 1. Read per-scaffold stats (genome-wide means) ─────────────────────────
results <- data.frame()

for (pop in pops) {
    stat_file <- list.files(theta_dir,
        pattern = paste0("^", pop, ".*\\.thetas\\.idx\\.pestPG$"),
        full.names = TRUE)
    if (length(stat_file) == 0) next

    st <- tryCatch(
        read.table(stat_file[1], header = TRUE, comment.char = ""),
        error = function(e) {
            cat("WARNING: Could not read", stat_file[1], ":", conditionMessage(e), "\n")
            return(NULL)
        }
    )
    if (is.null(st) || nrow(st) == 0) next

    cat("  ", pop, ": read", nrow(st), "scaffolds\n")

    # Genome-wide summaries (weighted by nSites)
    total_sites <- sum(st$nSites)
    tw <- sum(st$tW) / total_sites
    tp <- sum(st$tP) / total_sites
    td <- mean(st$Tajima, na.rm = TRUE)

    results <- rbind(results, data.frame(
        Pop = pop, Region = region_map[pop],
        theta_W = tw, pi = tp, TajimaD = td, nSites = total_sites
    ))
}

if (nrow(results) == 0) {
    cat("ERROR: No valid .pestPG files found. Exiting.\n")
    quit(save = "no", status = 1)
}

# ─── 2. Read sliding window stats (for CI) ──────────────────────────────────
window_data <- data.frame()

for (pop in pops) {
    win_file <- list.files(theta_dir,
        pattern = paste0("^", pop, ".*\\.thetaswindow\\.pestPG$"),
        full.names = TRUE)
    if (length(win_file) == 0) next

    wn <- tryCatch(
        read.table(win_file[1], header = TRUE, comment.char = ""),
        error = function(e) {
            cat("WARNING: Could not read window file for", pop, "\n")
            return(NULL)
        }
    )
    if (is.null(wn) || nrow(wn) == 0) next

    # Per-window per-site thetas
    wn <- wn[wn$nSites > 0, ]
    wn$pi_site <- wn$tP / wn$nSites
    wn$tw_site <- wn$tW / wn$nSites
    wn$Pop <- pop

    window_data <- rbind(window_data,
        wn[, c("Pop", "pi_site", "tw_site", "Tajima", "nSites")])
}

cat("Window data:", nrow(window_data), "windows across", length(unique(window_data$Pop)), "pops\n")

# ─── 3. Compute CI from windows ─────────────────────────────────────────────
if (nrow(window_data) > 0) {
    ci_summary <- window_data %>%
        group_by(Pop) %>%
        summarise(
            pi_se  = sd(pi_site, na.rm = TRUE) / sqrt(n()),
            tw_se  = sd(tw_site, na.rm = TRUE) / sqrt(n()),
            td_se  = sd(Tajima, na.rm = TRUE)  / sqrt(n()),
            .groups = "drop"
        )

    # Merge SE into results — CI centered on genome-wide estimate
    results <- left_join(results, ci_summary, by = "Pop")
    results <- results %>%
        mutate(
            pi_lo = pi - 1.96 * pi_se,
            pi_hi = pi + 1.96 * pi_se,
            tw_lo = theta_W - 1.96 * tw_se,
            tw_hi = theta_W + 1.96 * tw_se,
            td_lo = TajimaD - 1.96 * td_se,
            td_hi = TajimaD + 1.96 * td_se
        )
} else {
    # No window data — CI will be NA
    results$pi_lo <- results$pi_hi <- NA
    results$tw_lo <- results$tw_hi <- NA
    results$td_lo <- results$td_hi <- NA
}

results$Pop <- factor(results$Pop, levels = pop_order)
results$Region <- factor(results$Region, levels = c("Rankin", "Baker", "Naujaat"))

# ─── 4. Summary table ───────────────────────────────────────────────────────
cat("\n=== Diversity summary ===\n")
print(results[, c("Pop", "Region", "pi", "theta_W", "TajimaD", "nSites")], row.names = FALSE)

write.csv(results[, c("Pop", "Region", "pi", "theta_W", "TajimaD", "nSites")],
          file.path(table_dir, "Table_Diversity.csv"), row.names = FALSE)
cat("Table written with", nrow(results), "populations.\n")

# ─── 5. Publication-quality dot plots ────────────────────────────────────────
theme_pub <- theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "none",
        plot.margin = margin(5, 10, 5, 5)
    )

# A) Nucleotide diversity — x10^-3 scale
p_pi <- ggplot(results, aes(x = Pop, y = pi * 1000, color = Region)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = pi_lo * 1000, ymax = pi_hi * 1000),
                  width = 0.3, linewidth = 0.5) +
    scale_color_manual(values = region_colors) +
    labs(y = expression(pi ~ "(×10"^{-3}*")"),
         title = expression(bold("A)") ~ "Nucleotide diversity (" * pi * ")")) +
    theme_pub

# B) Watterson's theta — x10^-3 scale
p_tw <- ggplot(results, aes(x = Pop, y = theta_W * 1000, color = Region)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = tw_lo * 1000, ymax = tw_hi * 1000),
                  width = 0.3, linewidth = 0.5) +
    scale_color_manual(values = region_colors) +
    labs(y = expression(theta[W] ~ "(×10"^{-3}*")"),
         title = expression(bold("B)") ~ "Watterson's" ~ theta[W])) +
    theme_pub

# C) Tajima's D — natural scale, zero-line
p_td <- ggplot(results, aes(x = Pop, y = TajimaD, color = Region)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = td_lo, ymax = td_hi),
                  width = 0.3, linewidth = 0.5) +
    scale_color_manual(values = region_colors, name = "Region") +
    labs(y = "Tajima's D",
         title = expression(bold("C)") ~ "Tajima's D")) +
    theme_pub +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold"))

# Combine
combined <- p_pi / p_tw / p_td + plot_layout(heights = c(1, 1, 1.15))

ggsave(file.path(fig_dir, "Fig_Diversity.pdf"), combined,
       width = 7, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "Fig_Diversity.png"), combined,
       width = 7, height = 9, dpi = 300)

cat("Diversity figures saved.\n")