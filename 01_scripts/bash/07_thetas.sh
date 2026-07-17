#!/bin/bash
# =============================================================================
# 07_thetas.sh — Population diversity statistics (Watterson's theta, pi, Tajima's D)
# FIX: Uses $REALSFS_PATH saf2theta instead of angsd -doThetas (no BAM re-reading)
# Features: skips existing SAF/SFS/thetas, -rf consistency
# Submit: sbatch 01_scripts/bash/07_thetas.sh
# =============================================================================
#SBATCH -J "07_thetas"
#SBATCH -o logs/07_thetas_%j.out
#SBATCH -e logs/07_thetas_%j.err
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --time=7-00:00
#SBATCH --mem=200G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

module load angsd
ulimit -S -n 2048
NB_CPU=10
NSITES=5000000

# Use masked ancestral genome for thetas
ANC="${INFO_DIR}/genome_masked.fasta"
if [[ ! -f "$ANC" ]]; then
    ANC="${INFO_DIR}/genome.fasta"
    log_msg "WARNING: Masked genome unavailable, using reference as ancestral"
fi

POPS_ARRAY=($ALL_POPS)
mkdir -p "${THETA_DIR}"

for POP in "${POPS_ARRAY[@]}"; do
    log_msg "=== THETAS: ${POP} ==="

    BAMLIST="${INFO_DIR}/bamlists/${POP}.bamlist"
    require_file "$BAMLIST" "BAM list for ${POP}"

    read -r N_IND MIN_IND MAX_DEPTH <<< "$(compute_filters "$BAMLIST")"

    PREFIX="${THETA_DIR}/${POP}_pctind${PERCENT_IND}_maxdepth${MAX_DEPTH_FACTOR}"
    SAF_GZ="${PREFIX}.saf.gz"
    SAF_IDX="${PREFIX}.saf.idx"

    # ─── SAF (skip if already complete) ──────────────────────────────────
    if [[ -s "$SAF_GZ" ]] && [[ -s "$SAF_IDX" ]]; then
        log_msg "  SAF already complete for ${POP} ($(du -h "$SAF_GZ" | cut -f1)), skipping SAF..."
    else
        log_msg "  Computing SAF (all sites, no MAF filter)..."

        angsd -P $NB_CPU -underFlowProtect 1 \
          -doSaf 1 -GL 2 -doMajorMinor 5 -doCounts 1 \
          -anc "$ANC" \
          -rf "${INFO_DIR}/regions_nosex.txt" \
          -remove_bads 1 -minMapQ 30 -minQ 20 \
          -minInd "$MIN_IND" -setMaxDepth "$MAX_DEPTH" -setMinDepthInd "$MIN_DEPTH" \
          -b "$BAMLIST" \
          -out "$PREFIX"

        if [[ ! -f "$SAF_IDX" ]]; then
            log_msg "  ERROR: SAF failed for ${POP} (.saf.idx missing). Skipping."
            continue
        fi

        if [[ ! -s "$SAF_GZ" ]]; then
            log_msg "  ERROR: SAF file is empty for ${POP}. Skipping."
            continue
        fi
    fi

    # ─── 1D-SFS (skip if already done) ───────────────────────────────────
    SFS_RAW="${PREFIX}.${NSITES}"
    SFS_DSFS="${SFS_RAW}.dsfs"

    if [[ -s "$SFS_DSFS" ]]; then
        log_msg "  SFS already complete for ${POP}, skipping..."
    else
        log_msg "  Estimating 1D-SFS..."

        $REALSFS_PATH "$SAF_IDX" \
            -P $NB_CPU -nSites $NSITES -maxIter 50 \
            > "$SFS_RAW" 2>> logs/08_thetas_realSFS.log

        "$RSCRIPT" 01_scripts/R/sum_sfs.R "$SFS_RAW"

        if [[ ! -f "$SFS_DSFS" ]]; then
            log_msg "  ERROR: SFS folding failed for ${POP}. Skipping."
            continue
        fi
    fi

    # ─── Thetas via $REALSFS_PATH saf2theta (skip if already done) ─────────────
    THETA_IDX="${PREFIX}.thetas.idx"

    if [[ -s "$THETA_IDX" ]]; then
        log_msg "  Thetas already complete for ${POP}, skipping..."
    else
        log_msg "  Computing thetas with $REALSFS_PATH saf2theta..."

        # saf2theta: converts SAF + SFS prior into per-site thetas
        # This is MUCH faster than angsd -doThetas because it doesn't re-read BAMs
        $REALSFS_PATH saf2theta "$SAF_IDX" \
            -sfs "$SFS_DSFS" \
            -outname "$PREFIX" \
            -P $NB_CPU

        if [[ ! -f "$THETA_IDX" ]] || [[ ! -s "$THETA_IDX" ]]; then
            log_msg "  ERROR: saf2theta failed for ${POP}. Skipping."
            continue
        fi
    fi

    # ─── Per-scaffold + sliding window statistics ────────────────────────
    PESTPG="${PREFIX}.thetas.idx.pestPG"

    if [[ -s "$PESTPG" ]]; then
        log_msg "  thetaStat already complete for ${POP}, skipping..."
    else
        log_msg "  Computing per-scaffold stats..."
        "$THETASTAT_PATH" do_stat "$THETA_IDX"

        log_msg "  Computing sliding window stats (${WINDOW}bp / ${WINDOW_STEP}bp step)..."
        "$THETASTAT_PATH" do_stat "$THETA_IDX" \
            -win $WINDOW -step $WINDOW_STEP \
            -outnames "${PREFIX}.thetaswindow"
    fi

    log_msg "  ${POP} complete."
done

# ─── Summary table + figures ─────────────────────────────────────────────────
log_msg "Generating diversity summary table and figures..."

"$RSCRIPT" 01_scripts/R/plot_thetas.R \
    "${THETA_DIR}" \
    "${INFO_DIR}/pop.txt" \
    "${FIG_DIR}/main" \
    "${TABLE_DIR}/main" 2>&1 | tee -a logs/08_thetas_R.log

log_msg "=== THETAS COMPLETE ==="