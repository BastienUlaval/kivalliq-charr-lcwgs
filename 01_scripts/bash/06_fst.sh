#!/bin/bash
# =============================================================================
# 06_fst.sh — Pairwise FST for all 78 population pairs
# Includes global + sliding window FST
# Submit: sbatch 01_scripts/bash/06_fst.sh
# =============================================================================
#SBATCH -J "06_fst"
#SBATCH -o logs/06_fst_%j.out
#SBATCH -e logs/06_fst_%j.err
#SBATCH -c 6
#SBATCH -p large
#SBATCH --time=21-00:00
#SBATCH --mem=100G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

module load angsd/0.931
ulimit -S -n 2048
NB_CPU=6
NSITES=500000

REALSFS="${REALSFS_PATH}"
POPS_ARRAY=($ALL_POPS)
NUM_POPS=${#POPS_ARRAY[@]}

log_msg "=== PAIRWISE FST (${NUM_POPS} populations, $((NUM_POPS * (NUM_POPS - 1) / 2)) pairs) ==="

# ─── Subsample to equal N for fairness ───────────────────────────────────────
log_msg "Subsampling populations to equal size..."
"$RSCRIPT" 01_scripts/R/subsample_equal_n.R "${INFO_DIR}/bamlists" "${FST_DIR}" "${INFO_DIR}/pop.txt"

# ─── Compute SAF per population (using subsampled BAMs) ─────────────────────
log_msg "Computing SAF per population..."

SITES="${INFO_DIR}/sites_all_${SUFFIX}_canonical_nosex"

for POP in "${POPS_ARRAY[@]}"; do
    SUBSET_BAM="${FST_DIR}/${POP}_subset.bamlist"
    if [[ ! -f "$SUBSET_BAM" ]]; then
        log_msg "  WARNING: Subset BAM list missing for ${POP}, using full list"
        SUBSET_BAM="${INFO_DIR}/bamlists/${POP}.bamlist"
    fi

    read -r N_IND MIN_IND MAX_DEPTH <<< "$(compute_filters "$SUBSET_BAM")"

    log_msg "  SAF: ${POP} (N=${N_IND})"

    angsd -P $NB_CPU \
      -doSaf 1 -GL 2 -doMajorMinor 1 \
      -anc "${INFO_DIR}/genome.fasta" \
      -rf "${INFO_DIR}/regions_nosex.txt" \
      -sites "$SITES" \
      -remove_bads 1 -minMapQ 30 -minQ 20 \
      -minInd "$MIN_IND" -setMinDepthInd "$MIN_DEPTH" \
      -b "$SUBSET_BAM" \
      -out "${FST_DIR}/${POP}_${SUFFIX}"
done

# ─── Pairwise FST ───────────────────────────────────────────────────────────
log_msg "Computing pairwise FST..."

RESULTS="${FST_DIR}/fst_results.tsv"
echo -e "Pop1\tPop2\tFst_unweighted\tFst_weighted" > "$RESULTS"

for ((i = 0; i < NUM_POPS; i++)); do
    for ((j = i + 1; j < NUM_POPS; j++)); do
        POP1="${POPS_ARRAY[$i]}"
        POP2="${POPS_ARRAY[$j]}"
        PREFIX="${FST_DIR}/${POP1}_${POP2}"

        log_msg "  FST: ${POP1} vs ${POP2}"

        # 2D-SFS prior
        $REALSFS \
            "${FST_DIR}/${POP1}_${SUFFIX}.saf.idx" \
            "${FST_DIR}/${POP2}_${SUFFIX}.saf.idx" \
            -P $NB_CPU -maxIter 30 -nSites $NSITES \
            > "${PREFIX}.2dsfs"

        "$RSCRIPT" 01_scripts/R/sum_2dsfs.R "${PREFIX}.2dsfs"

        # FST index
        $REALSFS fst index \
            "${FST_DIR}/${POP1}_${SUFFIX}.saf.idx" \
            "${FST_DIR}/${POP2}_${SUFFIX}.saf.idx" \
            -sfs "${PREFIX}.2dsfs.summed" \
            -P $NB_CPU \
            -fstout "$PREFIX"

        # Global FST
        FST_LINE=$($REALSFS fst stats "${PREFIX}.fst.idx" -P $NB_CPU)
        FST_UW=$(echo "$FST_LINE" | awk '{print $1}')
        FST_W=$(echo "$FST_LINE" | awk '{print $2}')
        echo -e "${POP1}\t${POP2}\t${FST_UW}\t${FST_W}" >> "$RESULTS"

        log_msg "    FST(${POP1},${POP2}) = ${FST_W}"

        # Sliding window FST
        $REALSFS fst stats2 "${PREFIX}.fst.idx" \
            -win $WINDOW -step $WINDOW_STEP -P $NB_CPU \
            > "${PREFIX}.slidingwindow"

        # Per-site FST
        $REALSFS fst print "${PREFIX}.fst.idx" -P $NB_CPU \
            > "${PREFIX}.bypos.fst"
    done
done

# ─── Publication figures ─────────────────────────────────────────────────────
log_msg "Generating FST heatmaps and IBD plots..."

"$RSCRIPT" 01_scripts/R/plot_fst.R \
    "$RESULTS" \
    "${INFO_DIR}/pop.txt" \
    "${FIG_DIR}/main" \
    "${TABLE_DIR}/main"

log_msg "=== FST COMPLETE ==="
