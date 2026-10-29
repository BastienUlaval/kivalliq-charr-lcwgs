#!/usr/bin/env Rscript
# =============================================================================
# 14B_variance_partitioning.R - Latitude vs marine connectivity as drivers of
#   Atlantic introgression across the Kivalliq
#
# Variance partitioning is the primary analysis here, not model selection.
# AIC comparison forces a winner; the point of this analysis is that both axes
# contribute independently and that their relative weight depends on whether
# HOR (the only interior population) is included.
#
# Produces
#   99_tables/Table_VarPart_Connectivity.csv     fractions, all subsets
#   99_tables/Table_VarPart_PartialTests.csv     permutation tests
#   99_tables/Table_Connectivity_Sensitivity.csv univariate sensitivity
#   99_figures/Fig_VarPart_Connectivity.png      Venn + scatter panels
#
# Inputs
#   02_info/connectivity_covariates.csv          from 14A_marine_distances.R
#   99_figures/main/Table_Dstats_JAY.tsv         from 11B_dstats_figures.sh
#                                                (other location: DSTAT_TABLE=<path>)
#
# Run from the pipeline root (BASE_DIR):
#   $RSCRIPT 01_scripts/R/14B_variance_partitioning.R
# =============================================================================

if (nzchar(Sys.getenv("RLIB"))) .libPaths(c(Sys.getenv("RLIB"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ggplot2); library(patchwork)
})
HAS_REPEL <- requireNamespace("ggrepel", quietly = TRUE)
lab_layer <- function() {
  if (HAS_REPEL) ggrepel::geom_text_repel(size = 2.6, show.legend = FALSE,
        max.overlaps = 30, segment.size = 0.2, min.segment.length = 0,
        seed = 1, box.padding = 0.35)
  else geom_text(vjust = -0.9, size = 2.6, show.legend = FALSE)
}
suppressPackageStartupMessages({
  invisible(NULL)
})

BASE     <- Sys.getenv("BASE_DIR", getwd())
INFO     <- file.path(BASE, "02_info")
DSTAT_TABLE <- Sys.getenv("DSTAT_TABLE", file.path(BASE, "99_figures/main/Table_Dstats_JAY.tsv"))
TABDIR   <- file.path(BASE, "99_tables")
FIGDIR   <- file.path(BASE, "99_figures")
dir.create(TABDIR, showWarnings = FALSE)
dir.create(FIGDIR, showWarnings = FALSE)

N_PERM <- 9999
set.seed(20260825)

REGION_COLORS <- c(Rankin = "#2C7FB8", Naujaat = "#D95F02", Baker = "#1B9E77")

# ---------------------------------------------------------------------------
cov <- fread(file.path(INFO, "connectivity_covariates.csv"))
dst <- fread(DSTAT_TABLE)
setnames(dst, tolower(names(dst)))
dst <- dst[mode == "Observed"]     # one row per population

x <- merge(cov, dst[, .(pop = h2, Dstat = d, SD = sd, Z = z)], by = "pop")
setorder(x, region, pop)
cat(sprintf("Populations: %d\n", nrow(x)))
if (nrow(x) < 10) stop("Too few populations matched - check pop names")

print(x[, .(pop, region, lat, corridor = round(marine_corridor_km),
            d_FH = round(d_FH), D = round(Dstat, 3))])

# ---------------------------------------------------------------------------
# 1. Variance partitioning across subsets
#
# HOR is the only interior population and PAM is a distance outlier (its
# drainage flows north, so its route to Hudson Strait is genuinely long).
# Both are legitimate data points, but each is influential at n = 13, so the
# partition is reported with and without them rather than for one subset.
# ---------------------------------------------------------------------------
subsets <- list(
  "All (n=13)"        = x,
  "No HOR (n=12)"     = x[pop != "HOR"],
  "No PAM (n=12)"     = x[pop != "PAM"],
  "No HOR/PAM (n=11)" = x[!pop %in% c("HOR", "PAM")])

vp_rows <- rbindlist(lapply(names(subsets), function(nm) {
  s <- subsets[[nm]]
  vp <- varpart(s$Dstat, ~ lat, ~ marine_corridor_km, data = s)
  f  <- vp$part$indfract$Adj.R.square
  data.table(subset = nm, n = nrow(s),
             lat_only      = f[1],
             corridor_only = f[2],
             shared        = f[3],
             residual      = f[4],
             total         = sum(f[1:3]))
}))
fwrite(vp_rows, file.path(TABDIR, "Table_VarPart_Connectivity.csv"))

cat("\n=== Variance partitioning (adjusted R^2) ===\n")
print(vp_rows[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])

# ---------------------------------------------------------------------------
# 2. Permutation tests of the unique fractions
#
# varpart does not test fractions. Partial RDA does: each test asks whether one
# predictor explains variance the other cannot account for.
# ---------------------------------------------------------------------------
partial_test <- function(s, target, condition) {
  s <- as.data.frame(s); y <- matrix(s$Dstat, ncol = 1)
  f <- as.formula(sprintf("y ~ %s + Condition(%s)", target, condition))
  a <- anova(rda(f, data = as.data.frame(s)), permutations = N_PERM)
  data.table(target = target, conditioned_on = condition,
             Df = a$Df[1], Variance = a$Variance[1],
             Fstat = a$F[1], p = a$`Pr(>F)`[1])
}

pt <- rbindlist(lapply(names(subsets), function(nm) {
  s <- subsets[[nm]]
  r <- rbind(partial_test(s, "marine_corridor_km", "lat"),
             partial_test(s, "lat", "marine_corridor_km"))
  cbind(subset = nm, n = nrow(s), r)
}))
pt[, sig := fifelse(p < 0.001, "***", fifelse(p < 0.01, "**",
             fifelse(p < 0.05, "*", fifelse(p < 0.1, ".", "ns"))))]
fwrite(pt, file.path(TABDIR, "Table_VarPart_PartialTests.csv"))

cat("\n=== Partial RDA - unique fractions ===\n")
print(pt[, .(subset, target, Fstat = round(Fstat, 2), p = round(p, 4), sig)])

# ---------------------------------------------------------------------------
# 3. Univariate sensitivity
#
# Reported so readers can see that the ranking of single predictors flips with
# HOR and PAM, while the partitioned fractions are more stable. Suppressing
# this would misrepresent how leverage-dependent the univariate models are.
# ---------------------------------------------------------------------------
sens <- rbindlist(lapply(names(subsets), function(nm) {
  s <- subsets[[nm]]
  vars <- c("lat", "marine_corridor_km", "d_FH")
  vars <- vars[vapply(vars, function(v) sd(s[[v]], na.rm = TRUE) > 0, logical(1))]
  rbindlist(lapply(vars, function(v) {
    m <- summary(lm(as.formula(paste("Dstat ~", v)), data = s))
    data.table(subset = nm, n = nrow(s), predictor = v,
               adj_R2 = m$adj.r.squared, slope = coef(m)[2, 1], p = coef(m)[2, 4])
  }))
}))
fwrite(sens, file.path(TABDIR, "Table_Connectivity_Sensitivity.csv"))

cat("\n=== Univariate sensitivity ===\n")
print(sens[, .(subset, predictor, adj_R2 = round(adj_R2, 3), p = round(p, 4))])

# ---------------------------------------------------------------------------
# 4. Figure
# ---------------------------------------------------------------------------
# plot.varpart draws Xnames just outside each circle, so long names run off
# the device. Short names plus generous margins keep them inside.
png(file.path(FIGDIR, "Fig_VarPart_Venn.png"), width = 6.5, height = 5.8,
    units = "in", res = 300)
op <- par(mar = c(3, 4, 4, 4), xpd = NA)
vp_all <- varpart(x$Dstat, ~ lat, ~ marine_corridor_km, data = x)
plot(vp_all, Xnames = c("Latitude", "Corridor"),
     bg = c("#9ECAE1", "#FDBE85"), alpha = 90,
     digits = 2, cex = 1.15, id.size = 1.05)
title(main = sprintf("Atlantic introgression (D), n = %d", nrow(x)), line = 1.8)
gp <- function(tg, col) pt[subset == "All (n=13)" & target == tg][[col]]
mtext(sprintf("corridor | latitude: F = %.2f, p = %.3f      latitude | corridor: F = %.2f, p = %.3f",
              gp("marine_corridor_km", "Fstat"), gp("marine_corridor_km", "p"),
              gp("lat", "Fstat"), gp("lat", "p")),
      side = 1, line = 1, cex = 0.85)
par(op)
dev.off()

p_lat <- ggplot(x, aes(lat, Dstat, colour = region, label = pop)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              fill = "grey88", linewidth = 0.5) +
  geom_point(size = 2.8) +
  lab_layer() +
  scale_colour_manual(values = REGION_COLORS) +
  labs(x = "Latitude (deg N)", y = expression(italic(D))) +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p_cor <- ggplot(x, aes(marine_corridor_km, Dstat, colour = region, label = pop)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              fill = "grey88", linewidth = 0.5) +
  geom_point(size = 2.8) +
  lab_layer() +
  scale_colour_manual(values = REGION_COLORS) +
  labs(x = "Along-coast distance to Hudson Strait (km)", y = NULL) +
  theme_bw(base_size = 10)

frac <- melt(vp_rows[, .(subset, Latitude = lat_only,
                         Corridor = corridor_only, Shared = shared)],
             id.vars = "subset", variable.name = "fraction", value.name = "adjR2")
p_frac <- ggplot(frac, aes(subset, adjR2, fill = fraction)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(Latitude = "#2171B5", Corridor = "#E6550D",
                               Shared = "grey65")) +
  labs(x = NULL, y = expression("Adjusted"~R^2),
       title = "Unique and shared fractions across subsets") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(FIGDIR, "Fig_VarPart_Connectivity.png"),
       (p_lat | p_cor) / p_frac + plot_annotation(tag_levels = "A"),
       width = 9, height = 7.5, dpi = 300)

cat("\nWritten:\n")
cat("  ", file.path(FIGDIR, "Fig_VarPart_Connectivity.png"), "\n")
cat("  ", file.path(FIGDIR, "Fig_VarPart_Venn.png"), "\n")
cat("  3 tables in ", TABDIR, "\n", sep = "")