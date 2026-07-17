#!/bin/bash
# =============================================================================
# ANGSD Pipeline 2 — Master Configuration
# Arctic char (Salvelinus alpinus) lcWGS population genomics
# 13 Kivalliq populations (405 individuals) + JAY/LLS/DV reference populations
# =============================================================================
# Bastien Rubin — Université Laval, IBIS / Bernatchez Lab
# =============================================================================
#
# EDIT THE "USER-SPECIFIC PATHS" SECTION BELOW before running anything.
# Everything after that section should work as-is once those paths are set.
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# USER-SPECIFIC PATHS — EDIT THESE
# ─────────────────────────────────────────────────────────────────────────────
# Pipeline root: where 00_setup.sh will create the working directory tree
# (02_info/, 04_pca/, 05_admixture/, etc.). Scripts are normally submitted
# from this directory (or with SLURM_SUBMIT_DIR pointing here).
export BASE_DIR="/path/to/your/angsd_pipeline2"

# Reference genome (FASTA) and a version with deviant/paralogous regions
# masked (see 08_mask_deviants.sh). GENOME_MASKED is only needed for
# 07_thetas.sh; if you don't have one yet, 07_thetas.sh will fall back to
# the unmasked GENOME with a warning.
export GENOME="/path/to/your/reference/genome.fasta"
export GENOME_MASKED="/path/to/your/reference/genome_masked_deviants.fasta"

# Existing BAM files + sample metadata to seed the pipeline (see 00_setup.sh):
#   ${SOURCE_INFO_DIR}/bam.filelist   one BAM path per line, all individuals
#   ${SOURCE_INFO_DIR}/pop.txt        one population code per line
#   ${SOURCE_INFO_DIR}/info.txt       sample metadata (used by several R scripts)
#   ${SOURCE_INFO_DIR}/regions.txt    one chromosome/scaffold name per line
export SOURCE_INFO_DIR="/path/to/your/existing/bam_and_metadata"

# R install to use for ALL Rscript calls in this pipeline. Any system R with
# nlme + vegan installed together (e.g. a dedicated conda env) works; do NOT
# rely on a system-wide R module without checking for package conflicts first.
export RSCRIPT="/path/to/your/conda/envs/your_r_env/bin/Rscript"
export RLIB="/path/to/your/R/library"          # optional, only if .libPaths() is needed

# Program binaries / scripts
export NGSPARALOG_PATH="/path/to/programs/ngsParalog"
export NGSLD_PATH="/path/to/programs/ngsLD/ngsLD"
export NGSADMIX_PATH="/path/to/programs/ngsAdmix/ngsAdmix"
export PCANGSD_ENV="/path/to/conda/envs/pcangsd_py2"
export PCANGSD_PATH="/path/to/programs/pcangsd"
export REALSFS_PATH="/path/to/programs/angsd/misc/realSFS"
export THETASTAT_PATH="/path/to/programs/angsd/misc/thetaStat"
export PRUNE_SCRIPT="/path/to/programs/ngsLD/scripts/prune_ngsLD.py"
export FIT_LDDECAY="/path/to/programs/ngsLD/scripts/fit_LDdecay.R"
export ANGSD="/path/to/programs/angsd/angsd"
# Root of the ANGSD install (used to locate R/estAvgError.R for 11_dstats.sh).
# Leave unset to auto-detect from `command -v angsd`.
export ANGSD_DIR="/path/to/programs/angsd"

# =============================================================================
# Nothing below this line should normally need editing.
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# FILTERING PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
export MIN_MAF=0.05          # Minimum allele frequency
export PERCENT_IND=0.50      # Minimum proportion of individuals with data
export MIN_DEPTH=1           # Minimum read depth per individual
export MAX_DEPTH_FACTOR=8    # Max depth = N_IND * factor (removes paralogs/repeats)
export PVAL_THRESHOLD=0.001  # ngsParalog deviant p-value (BH-corrected); also
                              # reused as the SNP_pval for the triangle-plot
                              # diagnostic SNP panel (10B_triangle_call_sites.sh)

# ─────────────────────────────────────────────────────────────────────────────
# LD PRUNING PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
export LD_MAX_KB=200         # Maximum distance for LD calculation (kb)
export LD_RND_SAMPLE=0.5     # Random sampling fraction for ngsLD
export LD_PRUNE_DIST=50000   # Max distance for pruning (bp)
export LD_PRUNE_WEIGHT=0.1   # Minimum r² weight for pruning

# ─────────────────────────────────────────────────────────────────────────────
# SLIDING WINDOW (FST, thetas)
# ─────────────────────────────────────────────────────────────────────────────
export WINDOW=25000
export WINDOW_STEP=5000

# ─────────────────────────────────────────────────────────────────────────────
# ADMIXTURE (NGSadmix)
# ─────────────────────────────────────────────────────────────────────────────
export K_MIN=2
export K_MAX_GLOBAL=15
export K_MAX_RANKIN=7
export K_MAX_NAUJAAT=9
export NGSADMIX_NREP=10      # Number of independent replicates per K

# ─────────────────────────────────────────────────────────────────────────────
# POPULATIONS
# ─────────────────────────────────────────────────────────────────────────────
# Rankin Inlet populations
export RANKIN_POPS="AKL AUL CRB DIA MEL"
# Naujaat populations
export NAUJAAT_POPS="ITI KGJ NOP PAM SUP TIN WHI"
# Baker Lake population
export BAKER_POPS="HOR"
# All 13 Kivalliq populations
export ALL_POPS="AKL AUL CRB DIA MEL HOR ITI KGJ NOP PAM SUP TIN WHI"

# Reference populations for the triangle plot / D-statistics (not part of
# ALL_POPS): JAY (pure Arctic, Kitikmeot), LLS (pure Atlantic, Sweden),
# DV (Dolly Varden outgroup, Babbage River). Their bamlists must exist at
# ${INFO_DIR}/bamlists/{JAY,LLS,DV}.bamlist — see README for provenance.

# Sex-linked chromosomes to exclude (Arctic char sex determination regions)
export SEX_CHRS="NC_036841.1 NC_036842.1 NC_036843.1 NC_036851.1 NC_036862.1"

# ─────────────────────────────────────────────────────────────────────────────
# DIRECTORY STRUCTURE
# ─────────────────────────────────────────────────────────────────────────────
export INFO_DIR="${BASE_DIR}/02_info"
export NGSPARALOG_DIR="${BASE_DIR}/03A_ngsparalog"
export GL_DIR="${BASE_DIR}/03B_genotype_likelihoods"
export LD_DIR="${BASE_DIR}/03C_ld_pruning"
export PCA_DIR="${BASE_DIR}/04_pca"
export ADMIX_DIR="${BASE_DIR}/05_admixture"
export FST_DIR="${BASE_DIR}/07_fst"
export THETA_DIR="${BASE_DIR}/08_thetas"
export IBD_DIR="${BASE_DIR}/14_ibd"
export FIG_DIR="${BASE_DIR}/99_figures"
export TABLE_DIR="${BASE_DIR}/99_tables"
export LOG_DIR="${BASE_DIR}/logs"

# ─────────────────────────────────────────────────────────────────────────────
# FILE NAMING CONVENTION
# ─────────────────────────────────────────────────────────────────────────────
export SUFFIX="maf${MIN_MAF}_pctind${PERCENT_IND}_maxdepth${MAX_DEPTH_FACTOR}"

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

# Compute filtering thresholds from a BAM list
compute_filters() {
    local bamlist="$1"
    local n_ind
    n_ind=$(wc -l < "$bamlist")
    local min_ind
    min_ind=$(echo "($n_ind * $PERCENT_IND)" | bc -l)
    min_ind=${min_ind%.*}
    local max_depth
    max_depth=$(echo "($n_ind * $MAX_DEPTH_FACTOR)" | bc -l)
    max_depth=${max_depth%.*}
    echo "$n_ind $min_ind $max_depth"
}

# Log a timestamped message
log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Check that a required file exists; exit if not
require_file() {
    local f="$1"
    local desc="${2:-$f}"
    if [[ ! -f "$f" ]]; then
        log_msg "ERROR: Required file not found — $desc ($f)"
        exit 1
    fi
}

# Check that a required directory exists
require_dir() {
    local d="$1"
    if [[ ! -d "$d" ]]; then
        log_msg "ERROR: Required directory not found — $d"
        exit 1
    fi
}

log_msg "Configuration loaded: ${BASE_DIR}"
