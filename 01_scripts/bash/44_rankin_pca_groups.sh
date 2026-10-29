#!/bin/bash
# =============================================================================
# 44_rankin_pca_groups.sh — Check what drives the two groups in the Rankin PCA
#   1) groups: who is in each group (population, mapped reads, inbreeding F)
#   2) chroms: PCAngsd per chromosome on the Rankin pruned Beagle, compared with
#      the genome-wide PCA (inversion vs genome-wide signal)
# Uses only existing files: 04_pca/rankin_*.cov, rankin.bamlist,
# 03C_ld_pruning/rankin_*.beagle.gz, 18_inbreeding/F_individual.tsv, the BAMs.
# Submit from angsd_pipeline2/:  sbatch 01_scripts/bash/44_rankin_pca_groups.sh
# =============================================================================
#SBATCH -J "44_rankin_groups"
#SBATCH -p small
#SBATCH -o 98_log_files/44_rankin_groups_%j.out
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH --time=0-12:00

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh
module load samtools 2>/dev/null || true

OUT="${BASE_DIR}/04_pca/rankin_groups_check"
CDIR="${OUT}/by_chrom"
mkdir -p "$CDIR"
COV="${PCA_DIR}/rankin_${SUFFIX}_prunednosex.cov"
BAMLIST="${INFO_DIR}/bamlists/rankin.bamlist"
BEAGLE="${LD_DIR}/rankin_${SUFFIX}.prunednosex.beagle.gz"
FIND="${BASE_DIR}/18_inbreeding/F_individual.tsv"
R_SCRIPT=01_scripts/R/44_rankin_pca_groups.R
for f in "$COV" "$BAMLIST" "$BEAGLE" "$R_SCRIPT"; do require_file "$f"; done

# ---- 1. mapped reads per individual (fast: read from the BAM index) ----------
log_msg "Mapped reads per BAM (samtools idxstats)"
MAPPED="${OUT}/mapped_reads.tsv"; : > "$MAPPED"
while read -r bam; do
    id=$(basename "$bam" | sed 's/\..*$//')
    n=$(samtools idxstats "$bam" 2>/dev/null | awk '{s+=$3} END {print s+0}')
    echo -e "${id}\t${n}" >> "$MAPPED"
done < "$BAMLIST"

log_msg "Step 1: groups"
"$RSCRIPT" "$R_SCRIPT" groups "$COV" "$BAMLIST" "$OUT" "$MAPPED" "$FIND"

# ---- 2. PCAngsd per chromosome ----------------------------------------------
log_msg "Step 2: PCAngsd per chromosome"
PY="${PCANGSD_ENV}/bin/python"
# split the Beagle once by chromosome (marker = <chrom>_<pos>; chrom names may contain "_")
zcat "$BEAGLE" | awk -F'\t' -v d="$CDIR" 'NR==1 {h=$0; next}
    {m=$1; sub(/_[0-9]+$/, "", m); f=d"/"m".beagle";
     if (!(m in seen)) {print h > f; seen[m]=1}
     print > f}'
for sub in "${CDIR}"/*.beagle; do
    chr=$(basename "$sub" .beagle)
    n=$(( $(wc -l < "$sub") - 1 ))
    if (( n < 1000 )); then rm -f "$sub"; continue; fi
    echo "$n" > "${CDIR}/${chr}.nsnps"
    gzip -f "$sub"
    log_msg "  ${chr}: ${n} SNPs"
    "$PY" "${PCANGSD_PATH}/pcangsd.py" -beagle "${sub}.gz" -o "${CDIR}/${chr}" -threads 4 \
        -iter 200 -tole 1e-5 > /dev/null 2>&1 || log_msg "  PCAngsd failed on ${chr}"
    rm -f "${sub}.gz"
done

"$RSCRIPT" "$R_SCRIPT" chroms "$COV" "$CDIR" "$OUT"
log_msg "Done. Results in ${OUT}"
