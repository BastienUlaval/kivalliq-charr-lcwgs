#!/bin/bash
#SBATCH --job-name=tri_merge
#SBATCH --partition=medium
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=2
#SBATCH --output=99_logs/10C_%j.out
#SBATCH --error=99_logs/10C_%j.err

# =============================================================================
# 10C_triangle_merge_sites.sh -- Consolidate per-chromosome MAFs into a shared site list
#
# Input : 26_triangle/sites/per_chr/triangle_sites_<chr>.mafs.gz  (one per chr)
# Output: 26_triangle/sites/triangle_sites.sites   (chr  pos  major  minor)
#         26_triangle/sites/triangle_sites.chrs    (unique chr list)
#
# These files are used by 10D_triangle_anchor_maf.sh (-sites / -rf flags in ANGSD).
# =============================================================================

set -euo pipefail
source config/00_config.sh

PERCHDIR="26_triangle/sites/per_chr"
SITES_DIR="26_triangle/sites"
SITES_OUT="${SITES_DIR}/triangle_sites.sites"
CHRS_OUT="${SITES_DIR}/triangle_sites.chrs"

mkdir -p "$SITES_DIR"

echo "=== Merging per-chromosome SNP site lists ==="

N_CHR=0
N_SNP=0

# Build sites file: chr  pos  major  minor  (tab-separated, no header)
# ANGSD mafs.gz columns: chromo  position  major  minor  ref  knownEM  nInd
> "$SITES_OUT"

for f in "${PERCHDIR}"/triangle_sites_*.mafs.gz; do
    [[ -s "$f" ]] || { echo "  WARN: empty/missing $f, skipping"; continue; }
    CHR=$(basename "$f" | sed 's/triangle_sites_//; s/\.mafs\.gz//')
    N=$(zcat "$f" | tail -n +2 | wc -l)
    if [[ "$N" -eq 0 ]]; then
        echo "  WARN: $CHR has 0 SNPs, skipping"
        continue
    fi
    echo "  $CHR: $N SNPs"
    zcat "$f" | tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4}' >> "$SITES_OUT"
    N_CHR=$((N_CHR + 1))
    N_SNP=$((N_SNP + N))
done

echo ""
echo "Total chromosomes with SNPs : $N_CHR"
echo "Total SNPs in site list     : $N_SNP"

# Build chr list (unique, sorted)
cut -f1 "$SITES_OUT" | sort -u > "$CHRS_OUT"
N_CHRS_UNIQUE=$(wc -l < "$CHRS_OUT")
echo "Unique chromosomes in chr list : $N_CHRS_UNIQUE"

# Index sites for ANGSD
echo ""
echo "Indexing sites with ANGSD..."
${ANGSD:-angsd} sites index "$SITES_OUT"

echo ""
echo "=== Merge complete ==="
echo "  Sites file : $SITES_OUT"
echo "  Chrs file  : $CHRS_OUT"
echo "  ANGSD index: ${SITES_OUT}.bin / ${SITES_OUT}.idx"