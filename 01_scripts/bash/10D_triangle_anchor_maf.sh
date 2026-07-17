#!/bin/bash
#SBATCH --job-name=tri_anchor_maf
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=8
#SBATCH --output=99_logs/10D_%A_%a.out
#SBATCH --error=99_logs/10D_%A_%a.err
#SBATCH --array=0-1

# =============================================================================
# 10D_triangle_anchor_maf.sh — Compute per-pop MAF for LLS and JAY anchors
# Strategy:
#   1. Call SNPs once on combined Kivalliq+LLS+JAY dataset (-SNP_pval)
#      to obtain a shared site list (sites + chrs file)
#   2. Run ANGSD per anchor pop on the SAME sites with -doMajorMinor 3
#      so major/minor alleles are CONSISTENT across pops
#   3. Output per-pop .mafs.gz files used to compute delta-AF
#
# This array job processes both anchors in parallel:
#   array index 0 -> LLS,  1 -> JAY
# Site list is built ONCE in a separate job (10B_triangle_call_sites.sh).
# =============================================================================

set -euo pipefail

source config/00_config.sh

ANCHORS=(LLS JAY)
ANCHOR=${ANCHORS[$SLURM_ARRAY_TASK_ID]}

OUTDIR="26_triangle/maf_anchors"
SITES_PREFIX="26_triangle/sites/triangle_sites"
mkdir -p "$OUTDIR" 99_logs

BAMLIST="02_info/bamlists/${ANCHOR}.bamlist"
N_IND=$(wc -l < "$BAMLIST")
MIN_IND=$(awk -v n="$N_IND" -v p="$PERCENT_IND" 'BEGIN{printf "%d", n*p+0.999}')

echo "=== Computing MAF for anchor: $ANCHOR ==="
echo "  N individuals: $N_IND"
echo "  Min individuals (PERCENT_IND=$PERCENT_IND): $MIN_IND"
echo "  Sites file: ${SITES_PREFIX}.sites"
echo ""

if [[ ! -s "${SITES_PREFIX}.sites" ]]; then
    echo "ERROR: Site list not built yet. Run 10B_triangle_call_sites.sh first." >&2
    exit 1
fi

${ANGSD:-angsd} \
    -bam "$BAMLIST" \
    -ref "$GENOME" \
    -out "$OUTDIR/${ANCHOR}" \
    -sites "${SITES_PREFIX}.sites" \
    -rf    "${SITES_PREFIX}.chrs" \
    -GL 1 \
    -doMajorMinor 3 \
    -doMaf 1 \
    -doCounts 1 \
    -minMapQ 30 -minQ 20 \
    -minInd "$MIN_IND" \
    -nThreads "$SLURM_CPUS_PER_TASK"

echo ""
echo "Done: $OUTDIR/${ANCHOR}.mafs.gz"
zcat "$OUTDIR/${ANCHOR}.mafs.gz" | wc -l
