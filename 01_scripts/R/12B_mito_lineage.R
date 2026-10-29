# =============================================================================
# 12B_mito_lineage.R — Arctic or Atlantic mitochondrial lineage of each fish
#
# Input : 24_mito/mito_all.fasta (16,659 bp per fish, already aligned; 12A)
#         02_info/bamlists/<POP>.bamlist   (population of each sample)
#         04_figures/triangle/Table_HybridIndex.csv (nuclear hybrid index)
# Kept  : 13 Kivalliq populations + JAY + LLS (Greenland removed; Dolly Varden
#         used only to root the tree)
# Steps : 1) diagnostic sites = positions where the Atlantic reference (LLS)
#            and the Arctic reference (JAY) are each nearly fixed (>= 90%)
#            for a different base
#         2) each fish: share of diagnostic sites carrying the LLS base
#            (mtDNA does not recombine, so values should be ~0 or ~1)
#            -> lineage "Atlantic" if > 0.5, else "Arctic"
#         3) table per population, comparison with the nuclear hybrid index,
#            neighbour-joining tree rooted on Dolly Varden
#
# Run from angsd_pipeline2/:
#   $RSCRIPT 01_scripts/R/12B_mito_lineage.R
# =============================================================================
suppressPackageStartupMessages({ library(ape); library(data.table); library(ggplot2) })

OUT  <- "24_mito"
FA   <- file.path(OUT, "mito_all.fasta")
BLD  <- "02_info/bamlists"
HIF  <- "04_figures/triangle/Table_HybridIndex.csv"
FIX  <- 0.90       # within-reference frequency for a site to be diagnostic

# same ID rule as 12A_mito_extract.sh
sample_id <- function(x) {
    x <- basename(x)
    x <- sub("\\.(trimmed|sorted|dedup|no_overlap|merged|markdup|rmdup|realigned|filtered|bam).*$", "", x)
    sub("^.*\\.", "", x)
}

# ---- 1. alignment -----------------------------------------------------------
aln <- toupper(as.character(read.dna(FA, format = "fasta", as.matrix = TRUE)))
aln[!aln %in% c("A", "C", "G", "T")] <- NA          # N, gaps, ambiguity -> missing
cat("Alignment:", nrow(aln), "fish x", ncol(aln), "sites\n")

# ---- 2. population of each fish ---------------------------------------------
skip <- c("rankin", "naujaat", "triangle", "pipe30_kiv_LLS_JAY")
pm <- rbindlist(lapply(list.files(BLD, "\\.bamlist$", full.names = TRUE), function(f) {
    p <- sub("\\.bamlist$", "", basename(f))
    if (p %in% skip) return(NULL)
    data.table(sample = sample_id(readLines(f)), pop = p)
}))
pm <- unique(pm, by = "sample")
pop <- pm$pop[match(rownames(aln), pm$sample)]
# Greenland samples (not in a per-population bamlist): code in the sample name
gl <- is.na(pop) & grepl("sa..(NUK|QAQ|UUM|UMM)s", rownames(aln))
pop[gl] <- sub("^.*sa..(NUK|QAQ|UUM|UMM)s.*$", "\\1", rownames(aln)[gl])
pop[pop == "UUM"] <- "UMM"
pop[is.na(pop)] <- "unknown"
REGION <- c(AKL = "Rankin Inlet", AUL = "Rankin Inlet", CRB = "Rankin Inlet", DIA = "Rankin Inlet", MEL = "Rankin Inlet",
            ITI = "Naujaat", KGJ = "Naujaat", NOP = "Naujaat", PAM = "Naujaat", SUP = "Naujaat", TIN = "Naujaat", WHI = "Naujaat",
            HOR = "Baker Lake", JAY = "Arctic ref. (JAY)", LLS = "Atlantic ref. (LLS)", DV = "Dolly Varden",
            NUK = "Greenland", QAQ = "Greenland", UMM = "Greenland", HRN = "HRN")
region <- unname(REGION[pop]); region[is.na(region)] <- "unknown"
if (any(pop == "unknown")) {
    cat("\nFish without a population (not found in any <POP>.bamlist):\n")
    print(rownames(aln)[pop == "unknown"])
}
# Analysis restricted to the Kivalliq populations and the two references.
# Dolly Varden is kept only to root the tree; Greenland is removed.
KEEP <- c("AKL", "AUL", "CRB", "DIA", "MEL", "ITI", "KGJ", "NOP", "PAM", "SUP", "TIN", "WHI", "HOR", "JAY", "LLS")
dv_aln <- aln[pop == "DV", , drop = FALSE]
k <- pop %in% KEEP
aln <- aln[k, , drop = FALSE]; pop <- pop[k]; region <- region[k]
cat("\nFish per population (analysed):\n"); print(table(pop))

# ---- 3. diagnostic sites ----------------------------------------------------
major <- function(rows) apply(aln[rows, , drop = FALSE], 2, function(x) {
    x <- x[!is.na(x)]; if (!length(x)) return(c(NA, 0))
    t <- sort(table(x), decreasing = TRUE); c(names(t)[1], t[[1]] / length(x)) })
mL <- major(pop == "LLS"); mJ <- major(pop == "JAY")
diag <- which(!is.na(mL[1, ]) & !is.na(mJ[1, ]) & mL[1, ] != mJ[1, ] &
              as.numeric(mL[2, ]) >= FIX & as.numeric(mJ[2, ]) >= FIX)
cat("\nDiagnostic sites (LLS vs JAY):", length(diag), "\n")
if (length(diag) < 20) stop("Too few diagnostic sites: check that LLS and JAY are in the alignment")
atl <- mL[1, diag]; arc <- mJ[1, diag]

# ---- 4. lineage of each fish -----------------------------------------------
sub_aln <- aln[, diag, drop = FALSE]
n_atl <- rowSums(sweep(sub_aln, 2, atl, "=="), na.rm = TRUE)
n_arc <- rowSums(sweep(sub_aln, 2, arc, "=="), na.rm = TRUE)
ind <- data.table(sample = rownames(aln), pop = pop, region = region,
                  n_diag = n_atl + n_arc, prop_atlantic = round(n_atl / (n_atl + n_arc), 3))
ind[, lineage := fifelse(prop_atlantic > 0.5, "Atlantic", "Arctic")]
cat("\nFish with an intermediate value (0.1-0.9; check them):",
    ind[prop_atlantic > 0.1 & prop_atlantic < 0.9, .N], "\n")

# nuclear hybrid index
if (file.exists(HIF)) {
    hi <- fread(HIF)
    hi[, sample := sample_id(bam)]
    ind <- merge(ind, hi[, .(sample, HI)], by = "sample", all.x = TRUE)
}
fwrite(ind, file.path(OUT, "mito_lineage_individuals.tsv"), sep = "\t")

# ---- 5. per population ------------------------------------------------------
ORDER <- c("HOR", "CRB", "DIA", "MEL", "AUL", "AKL", "WHI", "ITI", "TIN", "KGJ", "SUP", "NOP", "PAM",
           "JAY", "LLS")
# haplotypes counted on the sites called in every fish (missing data would inflate the count)
full <- which(colSums(is.na(aln)) == 0)
hap <- apply(aln[, full, drop = FALSE], 1, paste, collapse = "")
cat("\nSites called in all fish (used to count haplotypes):", length(full), "\n")
tab <- ind[, .(n = .N, n_atlantic = sum(lineage == "Atlantic"),
               prop_atlantic = round(mean(lineage == "Atlantic"), 3),
               mean_HI = if ("HI" %in% names(ind)) round(mean(HI, na.rm = TRUE), 3) else NA_real_),
           by = .(region, pop)]
tab[, n_haplotypes := sapply(pop, function(p) length(unique(hap[ind$sample[ind$pop == p]])))]
tab <- tab[order(match(pop, ORDER))]
cat("\n--- Mitochondrial lineage per population ---\n"); print(tab)
fwrite(tab, file.path(OUT, "Table_mito_lineage_by_pop.tsv"), sep = "\t")

# ---- 6. figures ---------------------------------------------------------------
# (a) neighbour-joining tree, rooted on Dolly Varden
aln_t <- rbind(aln, dv_aln)
d  <- dist.dna(as.DNAbin(ifelse(is.na(aln_t), "n", tolower(aln_t))), model = "raw", pairwise.deletion = TRUE)
tr <- nj(d)
if (nrow(dv_aln) > 0) {
    tr <- root(tr, outgroup = rownames(dv_aln), resolve.root = TRUE)
    tr <- drop.tip(tr, rownames(dv_aln))           # DV used only for rooting
}
COLS <- c("Rankin Inlet" = "#D55E00", "Naujaat" = "#0072B2", "Baker Lake" = "#009E73",
          "Arctic ref. (JAY)" = "grey40", "Atlantic ref. (LLS)" = "#CC79A7")
pdf(file.path(OUT, "Fig_mito_NJ_tree.pdf"), width = 7, height = 9)
plot(tr, type = "fan", show.tip.label = FALSE, no.margin = TRUE)
tiplabels(pch = 16, cex = 0.5, col = COLS[region[match(tr$tip.label, rownames(aln))]])
legend("topleft", legend = names(COLS)[names(COLS) %in% region], col = COLS[names(COLS) %in% region],
       pch = 16, bty = "n", cex = 0.7)
invisible(dev.off())

# (b) proportion of Atlantic haplotypes per population
tab[, pop := factor(pop, levels = rev(intersect(ORDER, pop)))]
p <- ggplot(tab, aes(prop_atlantic, pop, fill = region)) + geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%d/%d", n_atlantic, n)), hjust = -0.1, size = 2.6) +
    scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    scale_fill_manual(values = COLS, breaks = intersect(names(COLS), tab$region)) +
    labs(x = "Proportion of Atlantic mtDNA haplotypes", y = NULL, fill = NULL) +
    theme_bw(base_size = 9)
ggsave(file.path(OUT, "Fig_mito_atlantic_by_pop.png"), p, width = 15, height = 11, units = "cm", dpi = 300)
cat("\nWritten in", OUT, ": mito_lineage_individuals.tsv, Table_mito_lineage_by_pop.tsv,",
    "Fig_mito_NJ_tree.pdf, Fig_mito_atlantic_by_pop.png\n")
