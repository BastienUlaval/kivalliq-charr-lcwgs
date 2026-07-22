#!/bin/bash
#SBATCH --job-name=tri_sites
#SBATCH --partition=medium
#SBATCH --time=7-00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output=99_logs/10B_%A_%a.out
#SBATCH --error=99_logs/10B_%A_%a.err
#SBATCH --array=0-99

# =============================================================================
# 10B_triangle_call_sites.sh (array) -- Call SNPs per chromosome to avoid OOM
#
# Each array task processes ONE chromosome (read from autosomes.chrs).
# Outputs go to 26_triangle/sites/per_chr/triangle_sites_<chr>.{mafs.gz, arg}
#
# After all tasks finish, run 10C_triangle_merge_sites.sh to consolidate.
#
# Strategy:
#   - Build autosome list ONCE before submitting (or inline if missing)
#   - Each task: only ~20G RAM, 4 CPUs (small footprint)
#   - --array=0-99 covers any reasonable number of autosomes; tasks beyond
#     the actual count exit cleanly
# =============================================================================

set -euo pipefail
source config/00_config.sh

OUTDIR="26_triangle/sites/per_chr"
SITES_DIR="26_triangle/sites"
mkdir -p "$OUTDIR" "$SITES_DIR" 99_logs

BAMLIST="02_info/bamlists/triangle.bamlist"
[[ -s "$BAMLIST" ]] || { echo "Missing $BAMLIST. Run 10A_triangle_bamlist.sh first." >&2; exit 1; }

N_IND=$(wc -l < "$BAMLIST")
MIN_IND=$(awk -v n="$N_IND" -v p="$PERCENT_IND" 'BEGIN{printf "%d", n*p+0.999}')

# --- Build autosome list (NC_ only, excluding sex chromosomes) ---------------
AUTO_RF="$SITES_DIR/autosomes.chrs"
FAI="${GENOME}.fai"
[[ -s "$FAI" ]] || { echo "Missing $FAI" >&2; exit 1; }

if [[ ! -s "$AUTO_RF" ]]; then
    if [[ -n "${SEX_CHRS:-}" ]]; then
        cut -f1 "$FAI" | grep "^NC_" | \
            grep -v -F -w -f <(echo "$SEX_CHRS" | tr ' ' '\n') > "$AUTO_RF"
    else
        cut -f1 "$FAI" | grep "^NC_" > "$AUTO_RF"
    fi
fi

N_AUTO=$(wc -l < "$AUTO_RF")
echo "Autosomes available: $N_AUTO"

# --- Skip if task index exceeds chromosome count -----------------------------
if [[ "$SLURM_ARRAY_TASK_ID" -ge "$N_AUTO" ]]; then
    echo "Task $SLURM_ARRAY_TASK_ID >= $N_AUTO autosomes -- exiting cleanly."
    exit 0
fi

# Get the chromosome for THIS task (1-indexed in awk; 0-indexed in SLURM)
CHR=$(awk -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR == i' "$AUTO_RF")
[[ -n "$CHR" ]] || { echo "Could not get chr for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

OUT_PREFIX="$OUTDIR/triangle_sites_${CHR}"

# --- Skip if already done (re-run safe) --------------------------------------
if [[ -s "${OUT_PREFIX}.mafs.gz" ]]; then
    N_LINES=$(zcat "${OUT_PREFIX}.mafs.gz" | wc -l)
    if [[ "$N_LINES" -gt 1 ]]; then
        echo "Skipping $CHR -- already done (${N_LINES} lines)"
        exit 0
    fi
fi

echo "=== Calling SNPs on $CHR (task $SLURM_ARRAY_TASK_ID) ==="
echo "  N individuals: $N_IND"
echo "  Min individuals: $MIN_IND  (PERCENT_IND=$PERCENT_IND)"
echo "  Min MAF: $MIN_MAF"
echo "  SNP p-value: $PVAL_THRESHOLD"
echo "  Output: $OUT_PREFIX"
echo ""

${ANGSD:-angsd} \
    -bam "$BAMLIST" \
    -ref "$GENOME" \
    -out "$OUT_PREFIX" \
    -r "${CHR}:" \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -SNP_pval "$PVAL_THRESHOLD" \
    -minMaf "$MIN_MAF" \
    -minInd "$MIN_IND" \
    -minMapQ 30 -minQ 20 \
    -doCounts 1 \
    -setMaxDepth $((N_IND * MAX_DEPTH_FACTOR)) \
    -nThreads "$SLURM_CPUS_PER_TASK"

N_SNP=$(zcat "${OUT_PREFIX}.mafs.gz" | tail -n +2 | wc -l)
echo ""
echo "Done $CHR: $N_SNP SNPs"